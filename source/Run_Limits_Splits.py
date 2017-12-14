#!/usr/bin/env python
import subprocess # So can use shell scripting in python
import os
import ROOT
from ROOT import *
import re

# *****************
# Lydia Beresford
# April 2015

# Script to run: Step 2 and Step 3
# - setLimitsOneMassPoint.cxx 
# - LimitSettingPhase.cxx 
# - plotLimitSetting.py to plot results of LimitSettingPhase.cxx
#
# Note:To change fit parameters, mjj cut (minXForFit) and other things, do this in Step2_setLimitsOneMassPoint.config!
# All log files stored 1 directory up in LogFiles directory, unless otherwise specified

# *****************

statspath = os.getcwd() # path used in outputFileName in config 
headdir = statspath.split("/Bayesian")[0] # directory for whole package 
logsdir =  headdir#"/data/atlas/atlasdata/beresford/StatisticalAnalysis" # Log Files Sent Here
batchdir = headdir#logsdir+"/SvnStatisticalAnalysis" #headdir # Folder to be copyied to batch, by default headdir unless otherwise specified 
Signals = {}

#---------------------------
# ***** User specifies *****
#---------------------------

# *NOTE* in path/file names use %d for mass, MMM for model !!!!!!!!!!!

#---------------------------
# Files and directories

dosetLimitsOneMassPoint = True # Set to True to run setLimitOneMassPoint.cxx 

doLimitSettingPhase = False # ONLY set to True when Step 2 setLimitsOneMassPoint has finished running on the batch (or locally)!!!
                            # Set to True to run LimitSettingPhase.cxx
doPlotting = False # ONLY set to True when Step 3 LimitSettingPhase.cxx has/is being run!!!
                   # Set to True to run plotLimitSetting.py
# FIXME
SearchPhaseresults = "%s/Bayesian/results/Step1_SearchPhase/Data_3p57invfb_IBLOff/Step1_SearchPhase_mjj_Data_2015_LLLfb.root"%headdir # Path/file of SearchPhase results
#SearchPhaseresults = "%s/Bayesian/results/Step1_SearchPhase/Data_3p34invfb_Partial/Step1_SearchPhase_mjj_Data_2015_LLLfb.root"%headdir # Path/file of SearchPhase results
                                                                                                                               # Can replace lumi in name with LLL and will be filled by Lumis value (must do this if running on multiple lumis)

##---------------------------
# Analysis quantities

Ecm = 13000.0 # Centre of mass energy in GeV

doPDFAccErr = True # Set to True to use PDF Acceptance Error 
PDFErrSize = 0.01
# FIXME
Lumis = ["3p57"] # Luminosities to scale limits to in fb e.g. "10" or 0p1, should have corresponding search phase ran for this lumi 
#Lumis = ["3p34"] # Luminosities to scale limits to in fb e.g. "10" or 0p1, should have corresponding search phase ran for this lumi 

##---------------------------
# Run controls 
 
# For Step 2 running setLimitsOneMassPoint.cxx  
setLimitsOneMassPointconfig = "./configurations/Step2_setLimitsOneMassPoint_MMM.config"   
useBatch = True # Set to True to run setLimitsOneMassPoint.cxx on the batch, or set to False to run locally. runs code in batchdir 
atOx = False # Set to True to use Oxford batch rather than lxbatch for running!

# Splits
nPETotal = 1000 # Total number of pseudo-experiments to run for expected limit bands
nSplits = 5 # Number of divisions to split PEs into
nPEForExpected = nPETotal/nSplits # Calculating how many PEs per division (split) performiFIXME move

# For Step 3 running LimitSettingPhase.cxx
LimitSettingPhaseconfig = "./configurations/Step3_LimitSettingPhase_MMM.config" # Path/file of LimitSettingPhase config


#---------------------------
# Signal information

plotextension = "Data_3p57invfb_IBLOff_Thousand"  # folder name to store all plots, results, configs and outputs consistently
#plotextension = "Data_3p34invfb_Partial"#_Thousand"  # folder name to store all plots, results, configs and outputs consistently
                                  # MUST keep name same when doing step 2 and 3
                                  # FIXME not doing lumi sub this time as want more accurate lumi SearchPhaseresults = "%s/Bayesian/results/Step1_SearchPhase/Data_fGRL_Hopeful_20150828/Step1_SearchPhase_mjj_Data_PeriodD_LLLfb.root"%headdir # Path/file of SearchPhase results

# Uncomment only 1 at a time depending which model you're running on

# QStar 
signalFileName = "inputs/MC15_20151017/QStar/dataLikeHists_v1/StatisticalHists/1fb/MMM%d_1fb.root" # Path/file of MC signal input file. Can always use 1 inv fb file, as signal normalisation floating in limit setting
Signals["QStar"]=["1000","2000","2500","3000","3500","4000","4500","5000","5500","6000","6500"] # Signals you're using e.g. QStar, followed by mass points for the signal 
outName = "QStar"

# QBH BlackMax  
#signalFileName = "inputs/MC15_20151017/BlackMax/dataLikeHists_v1/StatisticalHists/1fb/MMM%d_1fb.root" 
#Signals["BlackMax"]=["4000","5000","5500","6000","6500","7000","7500","8000","8500","9000","9500"]
#outName = "BlackMax"

# QBHRS0 
#signalFileName = "inputs/MC15_20151017/QBHRS0/dataLikeHists_v1/StatisticalHists/1fb/MMM%d_1fb.root" 
#Signals["QBHRS0"]=["4500","5000","5500","6000","6500","7000","7500","8000","8500"]
#outName = "QBHRS"

# WPrime
#signalFileName = "inputs/MC15_20151017/WPrime/dataLikeHists_v1/StatisticalHists/1fb/MMM%d_1fb.root" 
#Signals["WPrime"]=["1200","1500","1700","2000","2500","3000","3500","4000","4500","5000","5500","6000","6500"] 
#outName = "WPrime"

# ZPrime 
# Two new files: gSM0p10 2000 and gSM0p30 3500
#Coupling = "gSM0p10" 
#Signals["ZPrimemR"]=["1000","1500","2000","2500","3000","3500"] 
#outName = "ZPrime0p10"

#Coupling = "gSM0p20" 
#Signals["ZPrimemR"]=["1000","1500","2000","2500","3000","3500"] 
#outName = "ZPrime0p20"

#Coupling = "gSM0p30"  
#Signals["ZPrimemR"]=["1500","2000","2500","3000","3500"] 
#outName = "ZPrime0p30"

#Coupling = "gSM0p40"  
#Signals["ZPrimemR"]=["2000","2500","3000","3500"] 
#outName = "ZPrime0p40"

#Coupling = "gSM0p50"  
#Signals["ZPrimemR"]=["2500","3000","3500"] 
#outName = "ZPrime0p50"

# SKIPPED 
#Coupling = "gSM0p60" 
#Signals["ZPrimemR"]=["2500","3000","3500"] 
#outName = "ZPrime0p60"

#Coupling = "gSM0p70" 
#Signals["ZPrimemR"]=["3000","3500"]
#outName = "ZPrime0p70"

#Coupling = "gSM0p80" 
#Signals["ZPrimemR"]=["3000","3500"] 
#outName = "ZPrime0p80"

#Coupling = "gSM0p90" 
#Signals["ZPrimemR"]=["3500"] 
#outName = "ZPrime0p90"

#Coupling = "gSM1p00" 
#Signals["ZPrimemR"]=["3500"] 
#outName = "ZPrime1p00"

#LimitSettingPhaseconfig = "./configurations/Step3_LimitSettingPhase_MMM%s.config"%Coupling # Path/file of LimitSettingPhase config
#signalFileName = "inputs/MC15_20151017/ZPrime/dataLikeHists_v1/StatisticalHists/1fb/{0}/MMM%d{0}_1fb.root".format(Coupling)
## TEMPORARY? signalFileName = "inputs/MC15_20151017/ZPrime/dataLikeHists_v1/StatisticalHists/1fb/{0}/Temp/MMM%d{0}_1fb.root".format(Coupling)

# QBH
#signalFileName = "inputs/MC15_20151017/QBH/dataLikeHists_v1/StatisticalHists/1fb/MMM%d_1fb.root" 
#Signals["QBH0"]=["4500","5000","5500","6000","6500","7000","7500","8000","8500","9000","9500"]
#outName = "QBH"


#----------------------------------
# ***** End of User specifies *****
#----------------------------------

#----------------------
# Preliminary steps
#----------------------

# Make directories to store outputs if they don't exist already!
if dosetLimitsOneMassPoint: # For Step2 setLimitsOneMassPoint
  Step2_ConfigArchive = "%s/LogFiles/%s/Step2_setLimitsOneMassPoint/ConfigArchive"%(logsdir,plotextension)
  Step2_ScriptArchive = "%s/LogFiles/%s/Step2_setLimitsOneMassPoint/ScriptArchive"%(logsdir,plotextension)
  Step2_CodeOutput = "%s/LogFiles/%s/Step2_setLimitsOneMassPoint/CodeOutput"%(logsdir,plotextension)
  BATplotDirectory = "%s/LogFiles/%s/Step2_setLimitsOneMassPoint/BATPlots/"%(logsdir,plotextension)
  directories = ["./results/Step2_setLimitsOneMassPoint/%s"%plotextension,Step2_ConfigArchive,Step2_ScriptArchive,Step2_CodeOutput,BATplotDirectory]
  for directory in directories:
    if not os.path.exists(directory):
      os.makedirs(directory)

if doLimitSettingPhase: # For Step3 LimitSettingPhase
  Step3_ConfigArchive = "%s/LogFiles/%s/Step3_LimitSettingPhase/ConfigArchive"%(logsdir,plotextension)
  Step3_CodeOutput = "%s/LogFiles/%s/Step3_LimitSettingPhase/CodeOutput"%(logsdir,plotextension)
  directories = ["./results/Step3_LimitSettingPhase/%s"%plotextension, Step3_ConfigArchive, Step3_CodeOutput]
  for directory in directories:
    if not os.path.exists(directory):
      os.makedirs(directory)


#-------------------------------------
# Performing Step 2: Limit setting for each model, mass, lumi combo using setLimitsOneMassPoint.cxx
#-------------------------------------

# Loop over combination of models, masses and lumis
if dosetLimitsOneMassPoint:
  print "Performing setLimitsOneMassPoint"
  for Model in Signals.keys() :
    for Mass in sorted(Signals[Model]):
      for Lumi in Lumis:
        # open modified config file (fout) for writing
        fout = open('%s/Step2_%s%s_%sfb.config'%(Step2_ConfigArchive,outName,Mass,Lumi), 'w')

        # read in config file as fin and replace relevat fields with user inout specified at top of this file
        with open('%s'%setLimitsOneMassPointconfig.replace("MMM",Model), 'r') as fin:
          for line in fin:
            if (line.startswith("dataFileName") or line.startswith("signalFileName") or line.startswith("nominalSignalHist") or line.startswith("outputFileName") or line.startswith("plotDirectory") or line.startswith("plotNameExtension") or line.startswith("signame") or line.startswith("nominalTemplateJES") or line.startswith("nameTemp") or line.startswith("nPEForExpected")): 
              if line.startswith("dataFileName"):
                line = "dataFileName %s\n"%SearchPhaseresults.replace("LLL",Lumi)
                fout.write(line)
              if line.startswith("signalFileName"):
                line = "signalFileName %s\n"%signalFileName
                line = line.replace("MMM",Model)
                fout.write(line)
              if line.startswith("nominalSignalHist"):
                if Model == "ZPrimemR":
                  line = "nominalSignalHist mjj_MMM%dCCC_1fb_Nominal\n"
                  line = line.replace("MMM",Model)
                  line = line.replace("CCC",Coupling)
                else:
                  line = "nominalSignalHist mjj_MMM%d_1fb_Nominal\n"
                  line = line.replace("MMM",Model)
                fout.write(line)
              if line.startswith("outputFileName"):
                line = "outputFileName %s/results/Step2_setLimitsOneMassPoint/%s/Step2_setLimitsOneMassPoint_%s%s_%sfb.root\n"%(statspath,plotextension,outName,Mass,Lumi)
                fout.write(line)
              if line.startswith("plotDirectory"):
                line = "plotDirectory %s/\n"%BATplotDirectory
                fout.write(line)
              if line.startswith("plotNameExtension"):
                line = "plotNameExtension %s/\n"%outName
                fout.write(line)
              if line.startswith("signame"):
                line = "signame %s/\n"%outName
                fout.write(line)
              if line.startswith("nPEForExpected"):
                line = "nPEForExpected %s\n"%nPEForExpected
                fout.write(line)
              if line.startswith("nominalTemplateJES"):
                if Model == "ZPrimemR":
                  line = "nominalTemplateJES mjj_MMM%dCCC_1fb_Nominal\n"
                  line = line.replace("MMM",Model)
                  line = line.replace("CCC",Coupling)
                  fout.write(line)
                else:
                  fout.write(line)
              if line.startswith("nameTemp1"):
                if Model == "ZPrimemR":
                  line = "nameTemp1 mjj_MMM%dCCC_1fb_JET_GroupedNP_1\n"
                  line = line.replace("MMM",Model)
                  line = line.replace("CCC",Coupling)
                  fout.write(line)
                else:
                  fout.write(line)
              if line.startswith("nameTemp2"):
                if Model == "ZPrimemR":
                  line = "nameTemp2 mjj_MMM%dCCC_1fb_JET_GroupedNP_2\n"
                  line = line.replace("MMM",Model)
                  line = line.replace("CCC",Coupling)
                  fout.write(line)
                else:
                  fout.write(line)
              if line.startswith("nameTemp3"):
                if Model == "ZPrimemR":
                  line = "nameTemp3 mjj_MMM%dCCC_1fb_JET_GroupedNP_3\n"
                  line = line.replace("MMM",Model)
                  line = line.replace("CCC",Coupling)
                  fout.write(line)
                else:
                  fout.write(line)
            else:
              fout.write(line)
        fin.close()
        fout.close()

        # Setting command to be submitted (use tee to direc output to screen and to log file)
        submitcommand = ""
        if doPDFAccErr: # do PDF acceptance error
          command = "setLimitsOneMassPoint --config %s/Step2_%s%s_%sfb.config --mass %s --PDFAccErr %f |& tee %s/Step2_%s%s_%sfb_withPDFAccErr.txt"%(Step2_ConfigArchive,outName,Mass,Lumi,Mass,PDFErrSize,Step2_CodeOutput,outName,Mass,Lumi)
        else: # do not do PDF acceptance error
          command = "setLimitsOneMassPoint --config %s/Step2_%s%s_%sfb.config --mass %s |& tee %s/Step2_%s%s_%sfb.txt"%(Step2_ConfigArchive,outName,Mass,Lumi,Mass,Step2_CodeOutput,outName,Mass,Lumi)

        print "Submitting %i jobs for mass %s"%(nSplits,Mass)
        for split in range (0,nSplits):
          print split
          # Set different seed and output file name for each split
          Seed = eval(Mass)+split
          print "SEED"
          print Seed
          print type(Seed)
          outfile = "%s/results/Step2_setLimitsOneMassPoint/%s/Step2_setLimitsOneMassPoint_%s%s_%sfb_%i.root\n"%(statspath,plotextension,outName,Mass,Lumi,split)
          # Setting command to be submitted (use tee to direc output to screen and to log file)
          submitcommand = ""
          if doPDFAccErr: # do PDF acceptance error
            command = "setLimitsOneMassPoint --config %s/Step2_%s%s_%sfb.config --mass %s --PDFAccErr %f --seed %i --outfile %s |& tee %s/Step2_%s%s_%sfb_%i_withPDFAccErr.txt"%(Step2_ConfigArchive,outName,Mass,Lumi,Mass,PDFErrSize,Seed,outfile,Step2_CodeOutput,outName,Mass,Lumi,split)
          else: # do not do PDF acceptance error
            command = "setLimitsOneMassPoint --config %s/Step2_%s%s_%sfb.config --mass %s -- seed %i --outfile |& tee %s/Step2_%s%s_%sfb_%i.txt"%(Step2_ConfigArchive,outName,Mass,Lumi,Mass,Seed,outfile,Step2_CodeOutput,outName,Mass,Lumi,split)
	         
          # Perform setLimitsOneMassPoint locally
          if not useBatch:  
            subprocess.call(command, shell=True)

          # Use batch i.e. perform setLimitsOneMassPoint on the batch   
          if useBatch:
            if atOx:
              # Perform setLimitsOneMassPoint on batch
              print "Ox Batch!!"
              batchcommand = command.split("|&")[0]
              CodeOutputName = (command.split("|& tee ")[1]).split(".txt")[0] # Name of files for code output to be stored as
              print batchcommand
        
              # Open batch script as fbatchin
              fbatchin = open('./scripts/OxfordBatch/Step2_BatchScript_Template_Oxford.sh', 'r') 
              fbatchindata = fbatchin.read()
              fbatchin.close()
        
              # open modified batch script (fbatchout) for writing
              fbatchout = open('%s/Step2_BatchScript_Template_%s%s_%sfb_%i.sh'%(Step2_ScriptArchive,outName,Mass,Lumi,split),'w')
              fbatchoutdata = fbatchindata.replace("YYY",batchdir) # In batch script replace YYY for path for whole package
              fbatchoutdata = fbatchoutdata.replace("ZZZ",batchcommand) # In batch script replace ZZZ for submit command
              fbatchoutdata = fbatchoutdata.replace("OOO",CodeOutputName) # In batch script replace OOO (i.e. std output stream) to CodeOutput directory
              fbatchoutdata = fbatchoutdata.replace("EEE",CodeOutputName) # In batch script replace EEE (i.e. output error stream) to CodeOutput directory
              fbatchout.write(fbatchoutdata)    
 
              fbatchout.close()
              subprocess.call("qsub < %s/Step2_BatchScript_Template_%s%s_%sfb_%i.sh"%(Step2_ScriptArchive,outName,Mass,Lumi,split), shell=True)
            else:
              # Perform setLimitsOneMassPoint on batch
              print "Batch!!"
              batchcommand = command.split("|&")[0]
              CodeOutputName = (command.split("|& tee ")[1]).split(".txt")[0] # Name of files for code output to be stored as
              print batchcommand
          
              # Open batch script as fbatchin
              fbatchin = open('./scripts/Step2_BatchScript_Template.sh', 'r') 
              fbatchindata = fbatchin.read()
              fbatchin.close()
          
              # open modified batch script (fbatchout) for writing
              fbatchout = open('%s/Step2_BatchScript_Template_%s%s_%sfb_%s.sh'%(Step2_ScriptArchive,outName,Mass,Lumi,split),'w')
              fbatchoutdata = fbatchindata.replace("YYY",batchdir) # In batch script replace YYY for path for whole package
              fbatchoutdata = fbatchoutdata.replace("ZZZ",batchcommand) # In batch script replace ZZZ for submit command
              fbatchoutdata = fbatchoutdata.replace("OOO",CodeOutputName) # In batch script replace OOO (i.e. std output stream) to CodeOutput directory
              fbatchoutdata = fbatchoutdata.replace("EEE",CodeOutputName) # In batch script replace EEE (i.e. output error stream) to CodeOutput directory
              fbatchout.write(fbatchoutdata)    

              modcommand = 'chmod 744 %s/Step2_BatchScript_Template_%s%s_%sfb_%i.sh'%(Step2_ScriptArchive,outName,Mass,Lumi,split)
              print modcommand
              subprocess.call(modcommand, shell=True)
              subprocess.call("ls -l {0}".format(Step2_ScriptArchive), shell=True)
 
              fbatchout.close()
              command = "bsub -q 1nd %s/Step2_BatchScript_Template_%s%s_%sfb_%i.sh"%(Step2_ScriptArchive,outName,Mass,Lumi,split)
              print command
              subprocess.call(command, shell=True)


#-------------------------------------
# Performing Step 3: Combining mass points for Limit setting for each model and lumi combo using LimitSettingPhase.cxx
#-------------------------------------

if doLimitSettingPhase or doPlotting:
  for Model in Signals.keys() :
    for Lumi in Lumis:
      if doLimitSettingPhase: 
        print "Performing LimitSettingPhase"
        # open modified config file (fout) for writing
        fout = open('%s/Step3_%s_%sfb.config'%(Step3_ConfigArchive,outName,Lumi), 'w')

        # read in config file as fin and replace relevant fields
        with open('%s'%LimitSettingPhaseconfig.replace("MMM",Model), 'r') as fin:
          for line in fin:
            if (line.startswith("inputFileForm") or line.startswith("outputFileName") or line.startswith("luminosity") or line.startswith("nSplits1")): 
              if line.startswith("inputFileForm"): # Used %%d in line below as %d is special character in python so need extra % to allow %d
                line = "inputFileForm ./results/Step2_setLimitsOneMassPoint/%s/Step2_setLimitsOneMassPoint_%s%%d_%sfb_%%d.root\n"%(plotextension,outName,Lumi) 
                fout.write(line)
              if line.startswith("outputFileName"): 
                line = "outputFileName ./results/Step3_LimitSettingPhase/%s/Step3_LimitSettingPhase_%s_%sfb.root\n"%(plotextension,outName,Lumi)
                fout.write(line)
    
              if line.startswith("luminosity"):
                if "p" in Lumi: # converting to pb
                  lumi = Lumi.replace("p",".")
                  lumi = str(float(lumi)*1000)
                else: # converting to pb
                  lumi = str(float(Lumi)*1000)
                line = "luminosity %s\n"%lumi
                fout.write(line)
              if line.startswith("nSplits1"): 
                line = "nSplits1 %i\n"%(nSplits)
                fout.write(line)
            else:
              fout.write(line)
        fin.close()
        fout.close()
        # Perform LimitSettingPhase locally
        command = "LimitSettingPhase --config %s/Step3_%s_%sfb.config |& tee %s/Step3_%s%sfb.txt"%(Step3_ConfigArchive,outName,Lumi,Step3_CodeOutput,outName,Lumi)
        print command
        subprocess.call(command, shell=True)
 
      #-------------------------------------
      # Plotting results of running LimitSettingPhase.cxx using plotLimitSetting.py
      #-------------------------------------
    
      if doPlotting:
        if "p" in Lumi: # converting to pb
          lumi = Lumi.replace("p",".")
          lumi = float(lumi)*1000
        else: # converting to pb
          lumi = float(Lumi)*1000

        # open modified plotLimitSetting.py (fout) for writing
        fout = open('plotting/LimitSettingPhase/plotLimitSetting_%s_%sfb.py'%(outName,Lumi), 'w')

        # read in plotLimitSetting as fin
        with open('./plotting/LimitSettingPhase/plotLimitSetting.py', 'r') as fin:
          for line in fin:
            if (line.startswith("searchInputFile") or line.startswith("folderextension") or line.startswith("limitFileNameTemplate") or line.startswith("luminosity") or line.startswith("Ecm") or line.startswith("individualLimitFiles")): 
              
              if line.startswith("searchInputFile"):
                line = "searchInputFile = ROOT.TFile('%s')\n"%(SearchPhaseresults.replace("LLL",Lumi))
                fout.write(line)
              if line.startswith("folderextension"): 
                line = "folderextension = './plotting/LimitSettingPhase/plots/%s/%s/'\n"%(plotextension,Lumi)
                fout.write(line)
              if line.startswith("limitFileNameTemplate"): 
                line = "limitFileNameTemplate = './results/Step3_LimitSettingPhase/%s/Step3_LimitSettingPhase_{0}_%sfb.root'\n"%(plotextension,Lumi)
                fout.write(line)
              if line.startswith("individualLimitFiles"): 
                line = "individualLimitFiles = './results/Step2_setLimitsOneMassPoint/%s/Step2_setLimitsOneMassPoint_{0}{1}_%sfb_0.root'\n"%(plotextension,Lumi)
                fout.write(line)
              if line.startswith("luminosity"):
                line = "luminosity = %s\n"%str(lumi)
                fout.write(line)
              if "Ecm" in line:
                line = "Ecm = %d\n"%(Ecm/1000) 
                fout.write(line)
            else:
              fout.write(line)  
        fin.close()
        fout.close()
        # do plotting locally
        subprocess.call("python plotting/LimitSettingPhase/plotLimitSetting_%s_%sfb.py -b"%(outName,Lumi), shell=True)
        os.remove("./plotting/LimitSettingPhase/plotLimitSetting_%s_%sfb.py"%(outName,Lumi)) # Remove modified plotLimitSetting after plotting
