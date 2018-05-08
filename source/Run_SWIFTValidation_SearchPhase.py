import sys
import os
import subprocess
import re

## User controlled ##


# 2015 data
#fitFile = "inputs/TLASearchAndLimitSetting/mjj_fullDataset_06.root"
#fitHistogram = "mjj"
#lowEstimate =440
#highEstimate = 1236
#highFit = 7500
#workTag = "runSWIFT2015"
#
## 2016 data J75
#fitFile = "inputs/UnblindingTLA2016/NominalUnblinding/lookInsideTheBox.root"
#fitFile = "inputs/approvalTLA2016/lookInsideTheBoxAllSRs.root"
#fitHistogram = "Nominal/DSJ75yStar03_TriggerJets_J75_yStar03_mjj_2016binning_TLArange_data"
##fitHistogram = "Nominal/DSJ75yStar06_TriggerJets_J75_yStar06_mjj_2016binning_TLArange_data"
#lowEstimate = 400
#highEstimate = 2079
#workTag = "runSWIFT2016_J75yStar03"
#configInName = "configurations/Step1_SearchPhase_Swift_J75.config"
#windowSize = "13"
#
# Yvonne 2017 MC

#fitFile = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21Rebuild/input/feb2018/trijet_HLT_j380_inclusive.root"
#fitFile = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21Rebuild/input/btagged/reweighted_hist-background_ABCD_trijet.root"
fitFile = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21Rebuild/input/feb2018/trijet_HLT_j380_nbtag2.root"
#fitFile = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21Rebuild/input/btagged/trijet_HLT_j380_nbtag2.root"
#fitFile= "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts/../results2/signalplusbackground/signalplusbackground.Gauss_width3.0p1.ifb.mjj_Gauss_sig__smooth.4.par.0.seed.default.root "
#fitFile="/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts/Fluctuated_SwiftFittrijet_HLT_j380_inclusive.root"
#fitFile="/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts/../results2/signalplusbackground/signalplusbackground.Gauss_width3.0p1.ifb.mjj_Gauss_sig__smooth.4.par.0.seed.default.root"
#fitFile="/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21Rebuild/input/2018-4B-taggedNewWorkingPoint/trijet_mjj_nbtag2.root"

#fitFile = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21Rebuild/input/MC.root"
#fitHistogram = "Zprime_mjj_var"
fitHistogram = "background_mjj_var"
#fitHistogram = "basicBkgFrom4ParamFit_fluctuated"
#fitHistogram = "mjj_Gauss_sig_500_smoothinjectedToBkg"
#fitHistogram = "dijetgamma_g85_2j65/Zprime_mjj_var"
#fitHistogram = "Nominal/DSJ75yStar06_TriggerJets_J75_yStar06_mjj_2016binning_TLArange_data"
lowFit=300
highFit=1200
lowEstimate = 300
highEstimate = 3000
#workTag ="DijetISRMC-TrijetinclusiveNoUseScaled"
workTag ="DijetISRMC-Trijet2btaggedKateSwiftIssueReproduce-feb2018"
folderextension="data2017"
#configInName = "configurations/Step1_SearchPhase_Swift_dijetISR.config"
configInName = "configurations/Step1_SearchPhase_Swift_dijetISR.config"
#configInName = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts/submitConfigs/sensitivity_mjj_Gauss_sig_500_smoothinjectedToBkg_500_ww12/Step1_SearchPhase_case4Param_window12.config"

windowSize = "12"
lumi=35.5
Ecm=13000

'''
##2016 data J100
#fitFile = "inputs/UnblindingTLA2016/NominalUnblinding/lookInsideTheBox.root"
#fitHistogram = "Nominal/DSJ100yStar06_TriggerJets_J100_yStar06_mjj_2016binning_TLArange_data"
#lowEstimate = 530
#highEstimate = 2081
#workTag = "runSWIFT2016_J100"
#configInName = "configurations/Step1_SearchPhase_Swift_J100.config"
#windowSize = "9"
'''

scriptArchive = "submitConfigs/%s/"%workTag
outDir = "results/data2017/%s/"%workTag

### Common

for directory in [outDir,scriptArchive]:
  if not os.path.exists(directory): os.makedirs(directory)

useBatch = False
templatescript = "scripts/batchScript_template.sh"

saveOutput = True

doAlternate = False


## Automatic ##
# stdout_redirect = "1>> results/data2017/validateSWIFT/output_{4}"

#no useScaled
#commandTemplate = "SearchPhase --config {0} --file {1} --histName {2} --noDE --saveEx --outputfile {3}/SearchResultData_{4}.root \n".format("{0}",fitFile,fitHistogram,outDir,"{1}")
commandTemplate = "SearchPhase --config {0} --file {1} --histName {2} --noDE --outputfile {3}/SearchResultData_{4}.root \n".format("{0}",fitFile,fitHistogram,outDir,"{1}")
if saveOutput :
  commandTemplate = commandTemplate.replace("2>/dev/null","")

# Get current working directory
statspath = os.getcwd() # path used in outputFileName in config
headdir = statspath.split("/Bayesian")[0] # directory for whole package

def batchSubmit(command,stringForNaming) :

  # Perform setLimitsOneMassPoint on batch
  batchcommand = command.split("|&")[0]
  print batchcommand

  # Open batch script as fbatchin
  fbatchin = open(templatescript, 'r')
  fbatchindata = fbatchin.read()
  fbatchin.close()

  # open modified batch script (fbatchout) for writing
  batchtempname = '{0}/Step1_BatchScript_Template_SearchPhase_{1}.sh'.format(scriptArchive,stringForNaming)
  fbatchout = open(batchtempname,'w')
  fbatchoutdata = fbatchindata.replace("YYY",headdir) # In batch script replace YYY for path for whole package
  fbatchoutdata = fbatchoutdata.replace("ZZZ",batchcommand) # In batch script replace ZZZ for submit command
  fbatchout.write(fbatchoutdata)
  modcommand = 'chmod 744 {0}'.format(batchtempname)
  subprocess.call(modcommand, shell=True)
  fbatchout.close()
  submitcommand = "qsub {0}".format(batchtempname)
  print submitcommand
  subprocess.call(submitcommand, shell=True)

for case in ["D"]: #["A","B","C","D"] :
 for windowWidth in [windowSize]:#[7,10,13,17] :

   for doSwift in [True] : #True] :

    configOutName = "submitConfigs/"+workTag+"/Step1_SearchPhase_case{0}_window{1}_{2}.config".format(case,windowWidth,("doSwift" if doSwift else "noSwift" ))
    configOut = open(configOutName,'w')
    with open(configInName) as configInData :
      for line in configInData :
        if "doSwift" in line:
          if (doSwift) :
            line = "doSwift true\n"
          else :
            line = "doSwift false\n"
        if "minXForFit" in line :
          line = "minXForFit   {0}\n".format(lowFit)
        if "maxXForFit" in line :
          line = "maxXForFit   {0}\n".format(highFit)
        if "swift_minXAvailable"  in line :
          line = "swift_minXAvailable   {0}\n".format(lowEstimate)
        if "swift_maxXAvailable"  in line :
          #YEDIT
          #line = "swift_maxXAvailable   {0}\n".format(highFit)
          line = "swift_maxXAvailable   {0}\n".format(highEstimate)
#        if "inputHistDir"  in line and fitHistogram == "mjj" :
#          continue

        if doAlternate and "doAlternateFunction" in line :
          line = "doAlternateFunction  true\n"

        if "swift_nBinsLeft" in line :
          line = "swift_nBinsLeft  {0}\n".format(windowWidth)
        elif "swift_nBinsRight" in line :
          line = "swift_nBinsRight  {0}\n".format(windowWidth)
        if "A" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  false\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  false\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       true\n"
        if "B" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  false\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  true\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       true\n"
        if "C" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  true\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  false\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       true\n"
        if "D" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  true\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  true\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       false\n"
        configOut.write(line)
    configOut.close()

    modcommand = "chmod 744 {0}".format(configOutName)
    subprocess.call(modcommand,shell=True)

    extraString = "case{0}_window{1}".format(case,windowWidth)
    if (doSwift) :  extraString = extraString + "_doSwift"
    else : extraString = extraString + "_noSwift"
    thisCommand = commandTemplate.format(configOutName,extraString)

    if useBatch :
      batchSubmit(thisCommand,extraString)
    else :
      print "About to call",thisCommand
      print(thisCommand)
      subprocess.call(thisCommand,shell=True)

