import sys
import os
import subprocess

## User controlled ##

fitFile = "/cluster/warehouse/kpachal/StatisticalAnalysis/samples/TLA2017_SWIFTValidation/outfile.root"
fitHistogram = "J75yStar06_TriggerJets_J75_yStar06_mjj_2016binning_TLArange_data"
inputHistDir = "Nominal/"
outDir = "results/data2017/"

scriptArchive = "submitConfigs/test/"

useBatch = False
templatescript = "scripts/batchScript_template.sh"

saveOutput = True

## Automatic ##

commandTemplate = "SearchPhase --config {0} --file {1} --histName {2} --noDE --outputfile {3}/SearchResultData{4}.root 2>> /dev/null\n".format("{0}",fitFile,fitHistogram,outDir,"{1}")
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

for doSwift in [True] : #False,True] :

    nameExt = ""
    if doSwift :
      nameExt = "_Swift"

    configInName = "configurations/Step1_SearchPhase{0}.config".format(nameExt)
    configOutName = "submitConfigs/Step1_SearchPhase{0}.config".format(nameExt)
    configOut = open(configOutName,'w')
    with open(configInName) as configInData :
      for line in configInData :
        configOut.write(line)
    configOut.close()

    modcommand = "chmod 744 {0}".format(configOutName)
    subprocess.call(modcommand,shell=True)
 
    thisCommand = commandTemplate.format(configOutName,nameExt)

    if useBatch :
      batchSubmit(thisCommand,nameExt)
    else :
      print "About to call",thisCommand
      subprocess.call(thisCommand,shell=True)
