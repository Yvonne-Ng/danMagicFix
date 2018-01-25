#!/usr/bin/env python
import subprocess # So can use shell scripting in python
import os
import ROOT
from ROOT import *
import re

# *****************


#---------------------------
# Files and directories

# J75
'''
SearchPhaseresults = "results/data2017/runSWIFT2016_J75/SearchResultData_caseD_window13_doSwift.root"
SearchPhaseconfig = "submitConfigs/runSWIFT2016_J75/Step1_SearchPhase_caseD_window13_doSwift.config" 
workTag = "runSWIFT2016_J75_fixed" 
ZPrime_masses = ["350","450","550","650","750","850"]
signalFileNameTemplate = "inputs/TLA2016_MorphedSignalsFixed/J75_ystar03/dataLikeHistograms.%s.root"
Lumi = "3p57" # Luminosities to scale limits to in fb e.g. "10" or 0p1, should have corresponding search phase ran for this lumi 
'''

SearchPhaseresults = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/r21StatisticalAnalysis/source/results/data2017/DijetISRMC/SearchResultData_caseA_window13_doSwift.root"
#SearchPhaseconfig = "submitConfigs/runSWIFT2016_J75yStar03/Step1_SearchPhase_caseD_window13_doSwift.config"  
SearchPhaseconfig ="configurations/Step1_SearchPhase_Swift_dijetISR.config" 
Lumi = "3p57"
workTag = "DijetISRMC"
#ZPrime_masses = ["700","725","750","800","850","900","950","1000","1050","1100","1150","1200","1250","1300","1350","1400","1450","1500","1600","1700","1800"]
ZPrime_masses = ["500","550","600","650","700","725","750","800","850"] 
ZPrime_masses +=["900","950","1000","1050","1100","1150","1200","1250","1300","1350","1400","1450","1500","1600","1700","1800"]

#ZPrime_masses = /lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21Rebuild/input/Chopped_Ph100_ZPrimemR1500_gSM0p3.root

#signalFileNameTemplate = "inputs/TLA2016_MorphedSignalsFixed/dataLikeHists_J7503_1GeVBins_fixed_morphed/dataLikeHistograms.%s.root"

signalFileNameTemplate = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/r21StatisticalAnalysis/source/inputs/TLASearchAndLimitSetting/dataLikeHists_yStar03/dataLikeHistograms.%s.root"
#signalFileNameTemplate = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/r21StatisticalAnalysis/source/inputs/TLA2015_Signals/ZPrimemR%s.root"
#ZPrime_masses = ["450"]
#signalFileNameTemplate = "inputs/TLA2016_Signals/J75_ystar03/dataLikeHistograms.%s.root"

'''
SearchPhaseresults = "results/data2017/runSWIFT2016_J75yStar06/SearchResultData_caseD_window13_doSwift.root"
SearchPhaseconfig = "submitConfigs/runSWIFT2016_J75yStar06/Step1_SearchPhase_caseD_window13_doSwift.config"  
Lumi = "3p57"
workTag = "runSWIFT2016_J75yStar06"
#ZPrime_masses = ["700","725","750","800","850","900","950","1000","1050","1100","1150","1200","1250","1300","1350","1400","1450","1500","1600","1700","1800"]
#ZPrime_masses = ["700","725","750","800","850"] 
#ZPrime_masses = ["900","950","1000","1050","1100","1150","1200","1250","1300","1350","1400","1450","1500","1600","1700","1800"]
signalFileNameTemplate = "inputs/TLA2016_MorphedSignalsFixed/dataLikeHists_J7506_1GeVBins_fixed_morphed/dataLikeHistograms.%s.root"


#J100
SearchPhaseresults = "results/data2017/runSWIFT2016_J100/SearchResultData_caseD_window9_doSwift.root"
SearchPhaseconfig  = "submitConfigs/runSWIFT2016_J100/Step1_SearchPhase_caseD_window9_doSwift.config" 
workTag = "runSWIFT2016_J100_exoticsapproval"  #"runSWIFT2016_J100_fixed" 
ZPrime_masses = ["650","750","850","950","1050","1500","1700","2000"]
ZPrime_masses = ["700","725","750","800","850","900","950","1000","1050","1100","1150","1200","1250","1300","1350","1400","1450","1500","1600","1700","1800"]
#signalFileNameTemplate = "inputs/TLA2016_Signals/J100_ystar06/dataLikeHistograms.%s.root" 
signalFileNameTemplate = "inputs/TLA2016_MorphedSignalsFixed/dataLikeHists_J10006_1GeVBins_fixed_morphed/dataLikeHistograms.%s.root" 
Lumi = "29p30" 
'''

statspath = os.getcwd() # path used in outputFileName in config 

outDir = "results/data2017/{0}/".format(workTag)
#scriptArchive = "/cluster/warehouse/kpachal/StatisticalAnalysis/Bayesian/submitConfigs/{0}/".format(workTag)
scriptArchive = "submitConfigs/%s/"%workTag
headdir = (os.getcwd()).split("/Bayesian")[0] # directory for whole package
Signals = {}

directories = [outDir,scriptArchive]#,BATplotDirectory]
for directory in directories:
  if "eos" in directory:
    command = "eos mkdir {0}".format(directory)
    print command
  else :
    if not os.path.exists(directory):
      os.makedirs(directory)

# For easy cross checking
#doSwift = False
suppressPrintout = False

# In case I need to remake output file
doSearch = False

#---------------------------
# ***** User specifies *****
#---------------------------

# *NOTE* in path/file names use %d for mass, MMM for model !!!!!!!!!!!


if (doSearch) :
  commandTemplate = "SearchPhase --config {0} --file {1} --histName {2} --noDE --outputfile {3}/SearchResultData_caseD_window13_Kate.root \n".format(SearchPhaseconfig,SearchPhaseresults,"basicData",outDir)
  print commandTemplate
  subprocess.call(commandTemplate, shell=True)

##---------------------------
# Analysis quantities

Ecm = 13000.0 # Centre of mass energy in GeV
doPDFAccErr = True # Set to True to use PDF Acceptance Error 
PDFErrSize = 0.01
# FIXME


##---------------------------
# Run controls 
 
# For Step 2 running setLimitsOneMassPoint.cxx  
setLimitsOneMassPointconfig = "./configurations/Step2_setLimitsOneMassPoint_MMM.config"   
useBatch = True # Set to True to run setLimitsOneMassPoint.cxx on the batch, or set to False to run locally. runs code in batchdir 
templatescript = "scripts/batchScript_template_lunarc.sh"

# expected limits
doExpected = True
nPETotal = 200 # Total number of pseudo-experiments to run for expected limit bands
nSplits = 2 # Number of divisions to split PEs into
nPEForExpected = nPETotal/nSplits # Calculating how many PEs per division (split)
seedOffset = 400

#---------------------------
# Signal information

Signals["ZPrime0p05"]= ZPrime_masses
Signals["ZPrime0p10"]= ZPrime_masses
Signals["ZPrime0p20"]= ZPrime_masses
Signals["ZPrime0p30"]= ZPrime_masses
Signals["ZPrime0p40"]= ZPrime_masses

#Signals["ZPrime0p20"]= ["800"]

#----------------------------------
# ***** End of User specifies *****
#----------------------------------

def batchSubmit(command,stringForNaming, seed = -1) :

  # Perform setLimitsOneMassPoint on batch
  batchcommand = command.split("|&")[0]
  
  if seed > 0:
    batchcommand += " --seed " + str(seed)
    stringForNaming += "_seed" + str(seed)
  
  print batchcommand

  # Open batch script as fbatchin
  fbatchin = open(templatescript, 'r')
  fbatchindata = fbatchin.read()
  fbatchin.close()

  # open modified batch script (fbatchout) for writing
  batchtempname = '{0}/Step2_BatchScript_Template_Limits_{1}.sh'.format(scriptArchive,stringForNaming)
  fbatchout = open(batchtempname,'w')
  fbatchoutdata = fbatchindata.replace("YYY",headdir) # In batch script replace YYY for path for whole package
  fbatchoutdata = fbatchoutdata.replace("ZZZ",batchcommand) # In batch script replace ZZZ for submit command
  fbatchout.write(fbatchoutdata)
  modcommand = 'chmod 744 {0}'.format(batchtempname)
  subprocess.call(modcommand, shell=True)
  fbatchout.close()
  #submitcommand = "qsub {0}".format(batchtempname)
  #submitcommand = "bsub -q 1nd {0}".format(batchtempname)
  submitcommand = "sbatch {0}".format(batchtempname)
  print submitcommand
  subprocess.call(submitcommand, shell=True)

limitsString = '''


##########################################
# Background
##########################################

doFitError    true
nFitsInBkgError   100
doFitFunctionChoiceError true
doExtendedRange false

##########################################
# Lumi
##########################################

doLumiError       true
# percent lumi error
luminosityErr   0.05

##########################################
# Beam energy systematic
##########################################

doBeam     false
BeamFile          ./inputs/BeamUncertainty/AbsoluteBEAMUncertaintiesForPlotting.root

'''

JESString = '''

##--------------------------------------##
# nJES is number of extensions +1 
nJES      13
extension1        __3down
extension2        __2down5
extension3        __2down
extension4        __1down5
extension5        __1down
extension6        __0down5
extension7        __0up5
extension8        __1up
extension9        __1up5
extension10       __2up
extension11       __2up5
extension12       __3up
'''

#-------------------------------------
# Performing Step 2: Limit setting for each model, mass, lumi combo using setLimitsOneMassPoint.cxx
#-------------------------------------

# Loop over combination of models, masses and lumis
print "Performing setLimitsOneMassPoint"
for Model in Signals.keys():
  Coupling = Model.replace("ZPrime","")
  print "Using coupling",Coupling

  #signalFileName = "/cluster/warehouse/kpachal/StatisticalAnalysis/samples/TLA2017_validateWith2015/signals/ZPrimemR{0}gSM"+Coupling+"_1fb.root"
  #signalHist = "mjj_Scaled_mR{0}_mDM10_gSM"+Coupling+"_1fb_Nominal"
  #signalFileName = "inputs/TLA2015_Signals/ZPrimemR{0}gSM"+Coupling+"_1fb.root"
  #signalHist = "mjj_Scaled_mR{0}_mDM10_gSM"+Coupling+"_1fb_Nominal"
 

  for Mass in sorted(Signals[Model]):
  
    # try first naming convention
    special_mass_str = Mass[:-1][::-1][2:][::-1] + '.' + Mass[:-1][::-1][:2][::-1]
    if special_mass_str.startswith('.'): special_mass_str = '0' + special_mass_str
    
    signalName = "m{0}_g{1}".format(special_mass_str,Coupling.replace('p','.'))
    signalFileName = signalFileNameTemplate%signalName
    
    if not os.path.exists(signalFileName):
      mass_float = float(Mass)
      if mass_float%100 == 0:
        special_mass_str = "%1.1f"%(mass_float/1000.)
      elif mass_float%10 == 0:
        special_mass_str = "%1.2f"%(mass_float/1000.)
      else: special_mass_str = "%1.3f"%(mass_float/1000.)

      signalName = "m{0}_g{1}".format(special_mass_str,Coupling.replace('p','.'))
      signalFileName = signalFileNameTemplate%signalName

 #     if not os.path.exists(signalFileName):
 #       print "File for Mass = ", Mass, " Coupling = ", Coupling, "not found"
 #       continue

    signalHist = "mjj_Scaled_%s_1fb_Nominal"%signalName
    outName = "ZPrime"+Coupling+"_mZ{0}"
    for doSwift in [True] :
      # open modified config file (fout) for writing

      myString = "doSwift"
      if not doSwift : myString = "noSwift"

      configName = scriptArchive + '/Step2_{0}_{1}fb.config'.format(outName.format(Mass),Lumi)
      print "Creating config",configName
      fout = open(configName, 'w')

      # read in config file as fin and replace relevat fields with user inout specified at top of this file
      with open(SearchPhaseconfig, 'r') as fin:
        doWriteSection = True
        for line in fin:
          if line.startswith("dataHist"):
            line = "dataHist  basicData"
          if line.startswith("outputFileName"):
            MassString = "0p" + Mass[:-1]
            fout.write("##########################################\n# input/output for limit setting \n##########################################\n")
            fout.write("dataFileName {0}\n".format(SearchPhaseresults))
            fout.write("signalFileName {0}\n".format(signalFileName.format(Mass)))
            fout.write("nominalSignalHist {0}\n".format(signalHist.format(MassString)))
            outname = "outputFileName "+outDir+"Step2_setLimitsOneMassPoint_"+outName.format(int(Mass))+"_{0}.root\n"
            outname = outname.format(myString)
            fout.write(outname)
            fout.write("plotDirectory {0}\n".format(outDir))
            fout.write("plotNameExtension ZPrime\n")
            fout.write("signame     ZPrime\n")
            continue
          if line.startswith("doSwift") :
            if (doSwift) : line = "doSwift true\n"
            else : line = "doSwift false\n"
          fout.write(line)
        fout.write("##########################################\n# for limits\n##########################################\n")
        fout.write("nSigmas     3.\n")
        fout.write("doExpected    {0}\n".format( "true" if doExpected else "false"))
        fout.write("nPEForExpected    {0}\n".format(nPEForExpected))
        fout.write(limitsString)
        fout.write("##########################################\n# JES\n##########################################\n")
        fout.write("doJES true\nuseMatrices false\nuseTemplates true\n")
        sigHistNom = signalHist.format(MassString)
        fout.write("nominalTemplateJES {0}\n".format(sigHistNom))
        fout.write("nComponentsTemp 5\n")
        fout.write("nameTemp1 {0}\n".format(sigHistNom.replace("Nominal","JET_GroupedNP_1")))
        fout.write("nameTemp2 {0}\n".format(sigHistNom.replace("Nominal","JET_GroupedNP_2")))
        fout.write("nameTemp3 {0}\n".format(sigHistNom.replace("Nominal","JET_GroupedNP_3")))
        fout.write("nameTemp4 {0}\n".format(sigHistNom.replace("Nominal","JET_EtaIntercalibration_NonClosure")))
        fout.write("nameTemp5 {0}\n".format(sigHistNom.replace("Nominal","JET_TLA_ScaleFactor")))
        #
        fout.write(JESString)

      fin.close()
      fout.close()

      # Setting command to be submitted (use tee to direc output to screen and to log file)
      command = ""
      if doPDFAccErr: # do PDF acceptance error
        #command = "setLimitsOneMassPoint --config {0} --mass {1} --PDFAccErr {2} 2>/dev/null 1>output_{3}.txt".format(configName, Mass,PDFErrSize,myString)
        command = "setLimitsOneMassPoint --config {0} --mass {1} --PDFAccErr {2}".format(configName, Mass,PDFErrSize)
      else: # do not do PDF acceptance error
        command = "setLimitsOneMassPoint --config {0} --mass {1} ".format(configName,Mass) #2>/dev/null 1>output_{2}.txt

      if useBatch :
        if doExpected:
          for p in range(nSplits): batchSubmit(command,"{0}_{1}fb".format(outName.format(Mass),Lumi), seedOffset+p+1)
        else: batchSubmit(command,"{0}_{1}fb".format(outName.format(Mass),Lumi))
        
      # Perform setLimitsOneMassPoint locally
      else:  
        print command
        if doExpected:
          for p in range(nSplits):
            subprocess.call(command + " --seed " + str(seedOffset+p+1), shell=True)
        else: subprocess.call(command, shell=True)
    print("command: ", command)
