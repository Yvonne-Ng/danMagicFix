#!/usr/bin/env python

#=============================================================================
#          Steering script for running the generic Gaussian Limits           #
# 
# Instructions -- after setting up the Bayesian package and running the
# Search Phase, set configurables in this script and run with:
#   python calculateGaussianLimits_SWIFTValidation.py
# Some parameters still need to be set in the configuration file:
#    * fit functions & parameters
#    * resolution function & parameters
#=============================================================================

# Basic imports
import time
import os
import subprocess

# ROOT functionality
import ROOT

starttime = time.time()

#----------------------
# BEGIN user config
#----------------------

# Search phase input
SearchPhaseresults = "results/data2017/runSWIFT2016_J75yStar03/SearchResultData_caseD_window13_doSwift.root"
SearchPhaseconfig = "submitConfigs/runSWIFT2016_J75yStar03/Step1_SearchPhase_caseD_window13_doSwift.config"  
lumis = ['3p57']
config = 'configurations/GenericGaussians_Swift_J75.config'
workTag = "runSWIFT2016_J75yStar03"
subranges = [[400,500],[500,700],[700,900],[900,1300],[1300,2000]]

'''
SearchPhaseresults = "results/data2017/runSWIFT2016_J75yStar06/SearchResultData_caseD_window13_doSwift.root"
SearchPhaseconfig = "submitConfigs/runSWIFT2016_J75yStar06/Step1_SearchPhase_caseD_window13_doSwift.config"  
lumis = ['3p57']
config = 'configurations/GenericGaussians_Swift_J75.config'
workTag = "runSWIFT2016_J75yStar06"
subranges = [[400,500],[500,700],[700,900],[900,1300],[1300,2000]]


#J100
SearchPhaseresults = "results/data2017/runSWIFT2016_J100/SearchResultData_caseD_window9_doSwift.root"
SearchPhaseconfig = "submitConfigs/runSWIFT2016_J100/Step1_SearchPhase_caseD_window9_doSwift.config" 
lumis = ['29p30']
config = 'configurations/GenericGaussians_Swift_J100.config'
workTag = "runSWIFT2016_J100yStar06"
#subranges = [[700,800],[800,1000],[1000,1200],[1200,1500],[1500,1900]]
subranges = [[700,750],[750,800],[800,850],[850,900],[900,950],[950,1000],[1000,1200],[1200,1500],[1500,1900]]
'''
ratios = [0,0.05,0.07,0.10]  # 0 is resolution width!
# Run configuration
useBatch = True # use the batch system or not
# which batch template script to use
templatescript = "scripts/batchScript_template_lunarc.sh"

doSwift = True # do fit with sliding windows or not

# expected limits
doExpected = True
nPETotal = 500 # Total number of pseudo-experiments to run for expected limit bands
nSplits = 10 # Number of divisions to split PEs into
nPEForExpected = nPETotal/nSplits # Calculating how many PEs per division (split)
seedOffset = 500



# output directories

outDir = "/projects/hep/fs4/scratch/ATLAS/TLA/limits/" + "results/data2017/%s/"%workTag
scriptArchive = "submitConfigs/%s/"%workTag

#----------------------
# END user config
#----------------------

headdir = (os.getcwd()).split("/Bayesian")[0] # directory for whole package

directories = [outDir,scriptArchive]#,BATplotDirectory]
for directory in directories:
  if "eos" in directory:
    command = "eos mkdir {0}".format(directory)
    print command
  else :
    if not os.path.exists(directory):
      os.makedirs(directory)

swiftStr = "_doSwift" if doSwift else "_noSwift"

# Functions
def batchSubmit(batchcommand, ratio=-1, seed = -1) :

  print batchcommand

  # Open batch script as fbatchin
  fbatchin = open(templatescript, 'r')
  fbatchindata = fbatchin.read()
  fbatchin.close()
  
  extraStr = swiftStr
  if seed > 0:
    batchcommand += "--seed " + str(seed)
    extraStr += "_seed" + str(seed)

  if ratio < 0 :
    batchtempname = '{0}/Gaussian_BatchScript_Template_{1}fb_resolution_{2}{3}.sh'.format(scriptArchive,lumi,subrange[0],extraStr)
  else :
    batchtempname = '{0}/Gaussian_BatchScript_Template_{1}fb_r{2}_{3}{4}.sh'.format(scriptArchive,lumi,int(100*ratio),subrange[0],extraStr)
  fbatchout = open(batchtempname,'w')
  fbatchoutdata = fbatchindata.replace("YYY",headdir) # In batch script replace YYY for path for whole package
  fbatchoutdata = fbatchoutdata.replace("ZZZ",batchcommand) # In batch script replace ZZZ for submit command
  fbatchout.write(fbatchoutdata)
  modcommand = 'chmod 744 {0}'.format(batchtempname)
  subprocess.call(modcommand, shell=True)
  fbatchout.close()
  
  if 'lunarc' in templatescript:
    submitcommand = "sbatch {0}".format(batchtempname)
  else:
    submitcommand = "bsub -q 1nd {0}".format(batchtempname) #change qsub to bsub

  print submitcommand
  subprocess.call(submitcommand, shell=True)

#-------------------------------------
# Performing Step 2: Limit setting for each model, mass, lumi combo using setLimitsOneMassPoint.cxx
#-------------------------------------

args_to_copy = ["minXForFit","maxXForFit","swift_minXAvailable", "swift_maxXAvailable","swift_nBinsLeft","swift_nBinsRight","swift_fixLow","swift_fixHigh","swift_truncateHigh"]
lines_to_transfer = {}
with open(SearchPhaseconfig, 'r') as fin:
  for line in fin:
    for a in args_to_copy:
      if line.startswith(a):
        lines_to_transfer[a] = line
        break

for lumi in lumis :

  # open modified config file (fout) for writing
  thisConfigFileName = '{0}/Gaussians_{1}fb{2}.config'.format(scriptArchive,lumi,swiftStr)
  fout = open(thisConfigFileName, 'w')

  # read in config file as fin and replace relevat fields with user inout specified at top of this file
  with open(config, 'r') as fin:
    for line in fin:
      if line.startswith("inputFileName"):
        line = "inputFileName {0}\n".format(SearchPhaseresults)
      elif line.startswith("outputFileName"):
        theoutfile = outDir+"/GenericGaussians_{0}{1}\n".format(lumi,swiftStr)
        line = "outputFileName {0}\n".format(theoutfile)
      elif line.startswith("plotDirectory"):
        line = "plotDirectory {0}/\n".format(outDir)
      elif line.startswith("doSwift") :
        line = "doSwift true\n" if doSwift else "doSwift false\n"
      elif line.startswith("doExpected") :
        line = "doExpected true\n" if doExpected else "doExpected false\n"
      elif line.startswith("nPEForExpected") :
        line = "nPEForExpected {0}\n".format(nPEForExpected)
      else:
        for a in args_to_copy:
          if line.startswith(a):
            line = lines_to_transfer[a]
            break
      fout.write(line)
  fin.close()
  fout.close()

  for subrange in subranges :
    for ratio in ratios :

      # Setting command to be submitted (use tee to direc output to screen and to log file)
      if ratio > 0.0:
        command = "doGaussianLimits --config {0} --ratio {1} --rangelow {2} --rangehigh {3} ".format(thisConfigFileName, ratio,subrange[0],subrange[1]) #2>/dev/null
      else:
        command = "doGaussianLimits --config {0} --ratio {1} --useresolutionwidth --rangelow {2} --rangehigh {3} ".format(thisConfigFileName,ratio,subrange[0],subrange[1])

      # Run locally
      if not useBatch:
        print command
        if doExpected:
          for p in range(nSplits):
            subprocess.call(command + "--seed " + str(seedOffset+p+1), shell=True)
        else: subprocess.call(command, shell=True)

      # Run on batch
      if useBatch:
        if doExpected:
          for p in range(nSplits): 
            batchSubmit(command,ratio, seedOffset+p+1)
        else: batchSubmit(command,ratio)

print "Done."
