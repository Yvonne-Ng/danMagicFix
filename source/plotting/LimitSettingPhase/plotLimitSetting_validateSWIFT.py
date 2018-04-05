#!/usr/bin/env python

import ROOT
import os,sys
#import pprint
from art.morisot import Morisot
from array import array
from pprint import pprint
from decimal import Decimal
import SignalDictionaries
#import MorphedDictionaries

SignalTitles = {  "ZPrime0p05": "Z' (0.05)"
                , "ZPrime0p10": "Z' (0.10)"
                , "ZPrime0p20": "Z' (0.20)"
                , "ZPrime0p30": "Z' (0.30)"
                , "ZPrime0p40": "Z' (0.40)"
                , "ZPrime0p50": "Z' (0.50)"
                , "ZPrime0p60": "Z' (0.60)"
                }
SignalCouplings ={"ZPrime0p05": "0.05"
                , "ZPrime0p10": "0.10"
                , "ZPrime0p20": "0.20"
                , "ZPrime0p30": "0.30"
                , "ZPrime0p40": "0.40"
                , "ZPrime0p50": "0.50"
                , "ZPrime0p60": "0.60"
                }

ptCutList=["50","100"]
couplingList=["0p1","0p2", "0p3", "0p4"]
masses=[300, 350,400, 450,500, 550,750, 950,1500]

def rreplace(s, old, new, occurrence):
  li = s.rsplit(old, occurrence) 
  return new.join(li)

def makeBand(graph1, graph2):
  points = []
  for i in range(graph1.GetN()):
    points += [(i,graph1.GetX()[i],graph1.GetY()[i])]
  for i in range(graph2.GetN()-1,-1,-1):
    points += [(i,graph2.GetX()[i],graph2.GetY()[i])]
  graph_band = ROOT.TGraph();
  for i in range (len(points)): graph_band.SetPoint(i,points[i][1],points[i][2])
  return graph_band

def GetCenterAndSigmaDeviations(inputs) :
  inputs = sorted(inputs)
  statVals = []
  quantiles = [0.02275,0.1587,0.5,0.8413,0.9772]
  for q in quantiles:
    wantEvents = len(inputs)*q
    statVals.append(inputs[int(wantEvents)])
  return statVals

def getModelLimits(ptCut, coupling, these_massvals,individualLimitFiles, sig_dict, limitsDictOut, \
  cutstring = '', xname = "M_{Z'} [GeV]", yname = "#sigma #times #it{A} #times BR [pb]" ):
  print ptCut,"  ", coupling

  # Initialize painter
  myPainter = Morisot()
  myPainter.setColourPalette("ATLAS")
  myPainter.cutstring = cutstring
  #myPainter.setEPS(True)
  myPainter.setLabelType(2) # Sets label type i.e. Internal, Work in progress etc.
                            # See below for label explanation
  # 0 Just ATLAS    
  # 1 "Preliminary"
  # 2 "Internal"
  # 3 "Simulation Preliminary"
  # 4 "Simulation Internal"
  # 5 "Simulation"
  # 6 "Work in Progress"
  
  thisobserved = ROOT.TGraph()
  thisexpected = ROOT.TGraph()
  thisexpected_plus1  = ROOT.TGraph()
  thisexpected_minus1 = ROOT.TGraph()
  thisexpected_plus2  = ROOT.TGraph()
  thisexpected_minus2 = ROOT.TGraph()
  thistheory = ROOT.TGraph()
  #for mass in these_massvals:
  for  mass in masses:
    import glob
    #file_list = glob.glob(individualLimitFiles.format(signal,mass))
    individualLimitFiles=individualLimitFiles.replace("PPP", ptCut)
    individualLimitFiles=individualLimitFiles.replace("CCC", coupling)
    individualLimitFiles=individualLimitFiles.replace("MMM", str(mass))
    file_list = glob.glob(individualLimitFiles.format(ptCut,coupling, mass))
    print "individualLimitFiles", individualLimitFiles
    print("file_list: ", file_list)
    if len(file_list) == 0: continue
    allCLs = []
    PE_CLs = []
    for f in file_list:
      file = ROOT.TFile.Open(f)
      if not file or not file.Get("CLOfRealLikelihood"): continue
      CL = file.Get("CLOfRealLikelihood")[0]
      PE_tree = file.Get("ensemble_test")
      
      if not PE_tree or not CL: continue
      
      for event in PE_tree:
          PE_CLs.append( event.GetBranch("95quantile_marginalized_2").GetListOfLeaves().At(0).GetValue() )
      allCLs.append(CL)
    if len(allCLs) == 0: continue
    expCLs = GetCenterAndSigmaDeviations(PE_CLs)
    print mass, allCLs[0], expCLs[2], len(PE_CLs)
    m = float(mass)/1000.
    obsCL = allCLs[0]/luminosity
    expCLs = [e/luminosity for e in expCLs]
    thisobserved.SetPoint(       thisobserved.GetN(),m,obsCL)
    thisexpected_minus2.SetPoint(thisexpected_minus2.GetN(),m,expCLs[0])
    thisexpected_minus1.SetPoint(thisexpected_minus1.GetN(),m,expCLs[1])
    thisexpected.SetPoint(       thisexpected.GetN(),m,expCLs[2])
    thisexpected_plus1.SetPoint( thisexpected_plus1.GetN(),m,expCLs[3])
    thisexpected_plus2.SetPoint( thisexpected_plus2.GetN(),m,expCLs[4])

    #c = SignalCouplings[signal]
    #print sig_dict[c]
    #signal_info        = sig_dict[c]['%1.2f'%m]
    #signal_acc         = signal_info['acc']
    #signal_thxsec      = signal_info['theory']
    #signal_info['exp'] = expCLs[2]
    #signal_info['obs'] = obsCL
    #signal_info['exp+1'] = expCLs[3]    
    #signal_info['exp+2'] = expCLs[4] 
    #signal_info['exp-1'] = expCLs[1]  
    #signal_info['exp-2'] = expCLs[0] 
    #if c not in limitsDictOut: limitsDictOut[c] = {}
    limitsDictOut = {}
    #limitsDictOut['%1.2f'%m] = signal_info
    #t#histheory.SetPoint(thistheory.GetN(),m,signal_acc*signal_thxsec)
    
    #if not c in ZPrimeLimits: ZPrimeLimits[c] = {}
    #ZPrimeLimits[c][m] = {'obs':obsCL,'exp':expCLs[2],'th':signal_acc*signal_thxsec}

  if thisobserved.GetN() == 0:
    print "No limits found for couping: ",coupling, "ptCut: ",ptCut
    return limitsDictOut

  thisexpected1 = makeBand(thisexpected_minus1,thisexpected_plus1)
  thisexpected2 = makeBand(thisexpected_minus2,thisexpected_plus2)
  outputName = folderextension+"Limits_pH"+ptCut+'_gSM'+coupling+"_"+dataset+plotextension

  xlow  = 'automatic'# (int(masses[signal][0]) - 100)/1000.
  xhigh = 'automatic'#(int(masses[signal][-1]) + 100)/1000.


  #myPainter.drawLimitSettingPlotObservedExpected(thisobserved,thisexpected,thisexpected1, thisexpected2, thistheory,SignalTitles[signal],\
  #   outputName, xname,yname,luminosity,Ecm,xlow,xhigh,2E-4,100,False)
  myPainter.drawLimitSettingPlotObservedExpected(thisobserved,thisexpected,thisexpected1, thisexpected2,"",\
     "",outputName, xname,yname,luminosity,Ecm,xlow,xhigh,2E-4,100,False)
  return limitsDictOut

def writeLimitsDict(folderextension, dataset,limitsDictOut):

  limitsDictFileName = folderextension+'LimitsDict_'+dataset+'.py'
  print 'Writing signal limit info to', limitsDictFileName
  limitsDictFile=open(limitsDictFileName, 'w')
  import pprint
  pp = pprint.PrettyPrinter(indent=4,stream=limitsDictFile)
  pp.pprint(limitsDictOut)
  limitsDictFile.close()
  
if __name__ == "__main__":
  #====================================================================================
  # User Options
  #====================================================================================


  # Options
  folderextension = "./plots/"
  plotextension = ""
  # make plots folder i.e. make folder extension
  if not os.path.exists(folderextension):
      os.makedirs(folderextension)

  # Define necessary quantities.
  Ecm = 13
  signalInputFileTypes = ["ZPrime0p05","ZPrime0p10","ZPrime0p20","ZPrime0p30","ZPrime0p40"]
  masses_Zprime = ["450","500","550","600","650","700","725","750","800","850","900","950","1000","1050","1100","1150","1200","1250","1300","1350","1400","1450","1500","1600","1700","1800"]
  #signalInputFileTypes = ["ZPrime0p40"]

  #====================================================================================
  '''
  limitsDictOut = {}
  (dataset, luminosity, cutstring, sig_dict) = ("J100yStar06", 29.3*1000, "J100 |y*| < 0.6",SignalDictionaries.sig_J100y06)
  individualLimitFiles = "results/data2017/runSWIFT2016_J100_fixed/Step2_setLimitsOneMassPoint_{0}_mZ{1}_*seed*.root"
  #individualLimitFiles = "results/data2017/runSWIFT2016_J100_morphed/Step2_setLimitsOneMassPoint_{0}_mZ{1}_*seed1.root"
  individualLimitFiles = "results/data2017/runSWIFT2016_J100_exoticsapproval/Step2_setLimitsOneMassPoint_{0}_mZ{1}_*seed*.root"

  sig_dict = MorphedDictionaries.J10006_Dict
  

  # Loop over signals in signalInputFileTypes
  for signal in signalInputFileTypes :
    limitsDictOut = getModelLimits(signal, masses_Zprime,individualLimitFiles, sig_dict, limitsDictOut, cutstring)

  writeLimitsDict(folderextension, dataset,limitsDictOut)

  #====================================================================================

  limitsDictOut = {}
  (dataset, luminosity, cutstring, sig_dict) = ("J75yStar03", 3.57*1000, "J75 |y*| < 0.3",MorphedDictionaries.J7503_Dict)
  individualLimitFiles = "results/data2017/runSWIFT2016_J75yStar03/Step2_setLimitsOneMassPoint_{0}_mZ{1}_*.root"
  

  # Loop over signals in signalInputFileTypes
  for signal in signalInputFileTypes :
    limitsDictOut = getModelLimits(signal, masses_Zprime,individualLimitFiles, sig_dict, limitsDictOut, cutstring)

  writeLimitsDict(folderextension, dataset,limitsDictOut)
  
  '''
  limitsDictOut = {}
  #(dataset, luminosity, cutstring, sig_dict) = ("J75yStar06", 3.57*1000, "J75 |y*| < 0.6")#,MorphedDictionaries.J7506_Dict)
  (dataset, luminosity, cutstring) = ("J75yStar06", 3.57*1000, "J75 |y*| < 0.6")#,MorphedDictionaries.J7506_Dict)
  #individualLimitFiles = "results/data2017/runSWIFT2016_J75yStar06/Step2_setLimitsOneMassPoint_{0}_mZ{1}_*.root"
  individualLimitFiles = "/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/r21StatisticalAnalysis/source/results/Step2_setLimitsOneMassPoint/test_dijet_g150_2j25/Step2_setLimitsOneMassPoint_JDMPhPPP_ZprimeCCCMMM_35p45fb_0_seedMMM.root"

  # Loop over signals in signalInputFileTypes
#  for signal in signalInputFileTypes :
  for ptCut in ptCutList:
    for coupling in couplingList:
        limitsDictOut = getModelLimits(ptCut, coupling, masses_Zprime,individualLimitFiles, {}, limitsDictOut, cutstring)

  writeLimitsDict(folderextension, dataset,limitsDictOut)


  #====================================================================================


