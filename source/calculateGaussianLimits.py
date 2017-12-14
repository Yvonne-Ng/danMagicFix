#!/usr/bin/env python

# Basic imports
import time
import os
import subprocess

# ROOT functionality
import ROOT

config = 'configurations/GenericGaussians.config'
starttime = time.time()

lumis = ['37p4'] #15p7

#ratios = [-1,0.03,0.05,0.07,0.10,0.15,0.20] # -1 is resolution width!
#ratios = [0,0.01,0.03,0.05,0.07,0.10,0.15,0.20]
ratios = [0,0.03,0.07,0.10,0.15]
ratios = [0,0.03,0.07]

subranges = [[1100,1300],[1300,1500],[1500,1700],[1700,1900],[1900,2200],[2200,2500],[2500,3000],[3000,4000],[4000,5000],[5000,6000],[6000,7000]]

#----------------------
# Preliminary steps
#----------------------

# For storing
plotextension = "Data_Unblinded2016" #Data_ICHEP
useBatch = True

# for folding
doSfold = True
FoldingMatrixName = 'QCDPythia8.finerT5' #QCDPythia8 QCDPythia8.finerT5 QCDPythia8.matching Gaussian Gaussian.finerT5 CI_minus

fold_command = ""
Fname = "no_folding"
if doSfold :
    fold_command = "--doSfold"
    plotextension = plotextension + "_" + FoldingMatrixName + "_Sfold"
    Fname = './FoldingMatrix/TM.%s.root'%FoldingMatrixName #/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/transferMatrix

# Get current working directory
statspath = os.getcwd() # path used in outputFileName in config
headdir = statspath.split("/Bayesian")[0] # directory for whole package
logsdir = "/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/StatisticalAnalysis/Bayesian/LogFiles/" #"/cluster/warehouse/kpachal/DijetsSummer2016/LogFiles/" # Log Files Sent Here
batchdir = logsdir+"/SvnStatisticalAnalysis" #headdir # Folder to be copyied to batch, by default headdir unless otherwise specified

# FIXME
searchOutputDir = "/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/StatisticalAnalysis/Bayesian/results/Step1_SearchPhase/Data_Unblinded2016/" # "/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/StatisticalAnalysis/Bayesian/results/Step1_SearchPhase/Data_Unblinded2016/" #"/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/StatisticalAnalysis/Bayesian/results/Step1_SearchPhase/Data_ICHEP/" # Data_ICHEP Data_LHCP #"/cluster/warehouse/kpachal/DijetsSummer2016/Samples/SearchPhaseResults/"
templatefile = searchOutputDir+"Step1_SearchPhase_mjj_Data_2016_37p4fb.root" # "Step1_SearchPhase_mjj_Data_2016_15p7fb.root" "Step1_SearchPhase_mjj_Data_PeriodD_0p08fb.root"

resultsdir = "/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/StatisticalAnalysis/Bayesian/Gaussians/{0}/".format(plotextension)

#config = config.format(plotextension.split("_")[1])
print "Config",config

if not useBatch :
  resultsdir = resultsdir + "interactive/"

templatescript = './scripts/Step2_BatchScript_Template.sh'

# Make directories to store outputs if they don't exist already!
ConfigArchive = "{0}/LogFiles/{1}/Gaussians/ConfigArchive".format(logsdir,plotextension)
ScriptArchive = "{0}/LogFiles/{1}/Gaussians/ScriptArchive".format(logsdir,plotextension)
CodeOutput = "{0}/LogFiles/{1}/Gaussians/CodeOutput".format(logsdir,plotextension)

directories = [resultsdir,ConfigArchive,ScriptArchive,CodeOutput]#,BATplotDirectory]
for directory in directories:
  if "eos" in directory:
    command = "eos mkdir {0}".format(directory)
    print command
  else :
    if not os.path.exists(directory):
      os.makedirs(directory)

# Functions
def batchSubmit(command, ratio=-1) :

  # Perform setLimitsOneMassPoint on batch
  batchcommand = command
  print batchcommand

  # Open batch script as fbatchin
  fbatchin = open(templatescript, 'r')
  fbatchindata = fbatchin.read()
  fbatchin.close()

  # open modified batch script (fbatchout) for writing
  if ratio < 0 :
    batchtempname = '{0}/Gaussian_BatchScript_Template_{1}fb_resolution_{2}.sh'.format(ScriptArchive,lumi,subrange[0])
  else :
    batchtempname = '{0}/Gaussian_BatchScript_Template_{1}fb_r{2}_{3}.sh'.format(ScriptArchive,lumi,int(100*ratio),subrange[0])
  fbatchout = open(batchtempname,'w')
  fbatchoutdata = fbatchindata.replace("ZZZ",batchcommand) # In batch script replace ZZZ for submit command
  fbatchout.write(fbatchoutdata)
  modcommand = 'chmod 744 {0}'.format(batchtempname)
  subprocess.call(modcommand, shell=True)
  fbatchout.close()
  submitcommand = "bsub -q 1nd {0}".format(batchtempname) #change qsub to bsub
  print submitcommand
  subprocess.call(submitcommand, shell=True)


#-------------------------------------
# Performing Step 2: Limit setting for each model, mass, lumi combo using setLimitsOneMassPoint.cxx
#-------------------------------------

for lumi in lumis :

  # open modified config file (fout) for writing
  fout = open('{0}/Gaussians_{1}fb.config'.format(ConfigArchive,lumi), 'w')

  # read in config file as fin and replace relevat fields with user inout specified at top of this file
  with open(config, 'r') as fin:
    for line in fin:
      if line.startswith("inputFileName"):
        thefile = templatefile.format(plotextension,lumi)
        line = "inputFileName {0}\n".format(thefile)
        fout.write(line)
      elif line.startswith("outputFileName"):
        theoutfile = resultsdir+"GenericGaussians_{0}\n".format(lumi)
        line = "outputFileName {0}\n".format(theoutfile)
        fout.write(line)
      elif line.startswith("plotDirectory"):
        line = "plotDirectory {0}/\n".format(BATplotDirectory)
        fout.write(line)
      else:
        fout.write(line)
  fin.close()
  fout.close()

  for subrange in subranges :
    for ratio in ratios :

#      if not( (subranges.index(subrange) == 1 and ratios.index(ratio) == 2) ) :
#          (subranges.index(subrange) == 1 and ratios.index(ratio) == 2 ) or \
#          (subranges.index(subrange) == 2 and ratios.index(ratio) == 1) or \
#          (subranges.index(subrange) == 2 and ratios.index(ratio) == 2)) :
#          (subranges.index(subrange) == 2 and ratios.index(ratio) == 1 ) ) :
#        continue

      # Setting command to be submitted (use tee to direc output to screen and to log file)
      if ratio >= 0.0:
        command = "doGaussianLimits --config {0}/Gaussians_{1}fb.config --ratio {2} --rangelow {3} --rangehigh {4} {5} --Fname {6} 2>/dev/null".format(ConfigArchive,lumi,ratio,subrange[0],subrange[1],fold_command,Fname)
#        command = "doGaussianLimits --config {0}/Gaussians_{1}fb.config --ratio {2} --rangelow {3} --rangehigh {4}".format(ConfigArchive,lumi,ratio,subrange[0],subrange[1])
      else:
        command = "doGaussianLimits --config {0}/Gaussians_{1}fb.config --ratio {2} --useresolutionwidth --rangelow {3} --rangehigh {4} {5} --Fname {6} 2>/dev/null".format(ConfigArchive,lumi,ratio,subrange[0],subrange[1],fold_command,Fname)
#        command = "doGaussianLimits --config {0}/Gaussians_{1}fb.config --ratio {2} --useresolutionwidth --rangelow {3} --rangehigh {4}".format(ConfigArchive,lumi,ratio,subrange[0],subrange[1])
        #command = "doGaussianLimits --config {0}/Gaussians_{1}fb.config --ratio -1 --useresolutionwidth --rangelow {2} --rangehigh {3} |& tee {5}/Gaussians_{1}fb.txt".format(ConfigArchive,lumi,subrange[0],subrange[1],CodeOutput)

      # Perform setLimitsOneMassPoint locally
      if not useBatch:
        subprocess.call(command, shell=True)

      # Use batch i.e. perform setLimitsOneMassPoint on the batch
      if useBatch:

        # Perform setLimitsOneMassPoint on batch
        batchSubmit(command,ratio)

#    if not(\
#         (subranges.index(subrange) == 3) or \
#         subranges.index(subrange) == 4) :

print "Done."
