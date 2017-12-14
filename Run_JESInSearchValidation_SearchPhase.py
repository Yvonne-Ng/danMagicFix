import sys
import os
import subprocess

## User controlled ##

fitFile = "/cluster/warehouse/kpachal/StatisticalAnalysis/samples/TLA2017_BumpHunterDevelopment/systInputHist_HistogramsOfAsMuchJ75AsWeCouldRunOnTheGrid.root.root"
fitHistogram = "mjj_Nominal"
inputHistDir = ""
outDir = "results/data2017/"

scriptArchive = "submitConfigs/test/"

useBatch = False
templatescript = "scripts/batchScript_template.sh"

saveOutput = False

## Automatic ##

commandTemplate = "SearchPhase_JESSysts --config {0} --file {1} --histName {2} --noDE --saveToyFits --useScaled --outputfile {3}/SearchResultData_JESTests{4}.root 2>/dev/null\n".format("{0}",fitFile,fitHistogram,outDir,"{1}")
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

for doSwift in [False] :

 for doFit in [True,False] :

  for doJES in [True,False] :

    nameExt = "_test"
    if doFit :
      nameExt = nameExt+"_doFit"
    else :
      nameExt = nameExt+"_noFit"
    if doJES :
      nameExt = nameExt+"_doJES"
    else :
      nameExt = nameExt+"_noJES"

    if doSwift :
      nameExt = nameExt+"_doSwift"

    configInName = "configurations/Step1_SearchPhase_withJES.config"
    configOutName = "submitConfigs/Step1_SearchPhase_withJES{0}.config".format(nameExt)
    configOut = open(configOutName,'w')
    with open(configInName) as configInData :
      for line in configInData :
        if "inputHistDir" in line :
          line = "inputHistDir  {0}\n".format(inputHistDir)
        if "doSwift" in line :
          if doSwift :
            line = "doSwift  true\n"
          else :
            line = "doSwift  false\n"
        if "fitToys" in line :
          if doFit :
            line = "fitToys   true\n"
          else :
            line = "fitToys   false\n"
        if "doJES" in line :
          if doJES :
            line = "doJES     true\n"
          else :
            line = "doJES     false\n"
        if "nPseudoExp" in line : 
          line = "nPseudoExp  200\n"
        if "nPseudoExpFit" in line :
          line = "nPseudoExpFit 100\n"
        configOut.write(line)
    configOut.close()

    modcommand = "chmod 744 {0}".format(configOutName)
    subprocess.call(modcommand,shell=True)
 
    thisCommand = commandTemplate.format(configOutName,nameExt)

    print "this command is", thisCommand

    if useBatch :
      batchSubmit(thisCommand,nameExt)
    else :
#      print 'placeholder'
      subprocess.call(thisCommand,shell=True)
