
#/usr/bin/env python

import os
import sys
import ROOT
from art.morisot import Morisot
from itertools import repeat


names = { '0.00':  '#sigma_{G}/m_{G} = Res.',
          '0.05': '#sigma_{G}/m_{G} = 0.05',
          '0.07': '#sigma_{G}/m_{G} = 0.07',
          '0.10': '#sigma_{G}/m_{G} = 0.10',
          '0.15': '#sigma_{G}/m_{G} = 0.15'}

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
  
def getGaussianLimits( inputfileform, ratios, dataset, luminosity, cutstring, makePlot = True, outfolder = "./plots/"):

  basicInputFiles = {}
  for r in ratios: basicInputFiles[r] = []

  import glob
  file_list = glob.glob(inputfileform)
  for f in file_list:
    ratio_str = f.split('.')[-2].split('_')[-1]
    if ratio_str == 'resolutionwidth': ratio = 0.0
    else: ratio = float(ratio_str)/1000.
    if ratio in basicInputFiles: basicInputFiles[ratio].append(f)

  # Initialize painter
  myPainter = Morisot()
  # Internal
  myPainter.setLabelType(2)
  myPainter.cutstring = cutstring
  #myPainter.setLabelType(1)

  minMassVal = {}
  values = {}
  mass_list = []
  allobserved = []
  allexpected = []
  allexpected1Sigma = []
  allexpected2Sigma = []
  
  results = {}
  
  for r in ratios:
    values[r] = {}
    for f in basicInputFiles[r]:
      file = ROOT.TFile.Open(f)
      CLs = file.Get("CLsPerMass_widthToMass%d"%(r*1000))
      masses = file.Get("massesUsed")

      if masses == None: continue

      for i,mass in enumerate(masses) :
      
          #if dataset == "J100" and mass < 700: continue
          if "J75" in dataset and not makePlot and mass > 700: continue
          if mass > 1850: continue
          
          mass_list += [mass]
          PE_tree = file.Get("ensemble_tree_%d_%d"%(mass,r*1000))
          PE_CLs = []
          for event in PE_tree:
              PE_CLs.append( event.GetBranch("95quantile_marginalized_1").GetListOfLeaves().At(0).GetValue() )
          expCLs = GetCenterAndSigmaDeviations(PE_CLs)

          #print mass, CLs[i]/luminosity, [e/luminosity for e in expCLs]
          if mass not in values[r]:
            values[r][mass] = {'obs': [], 'exp': [], 'PEs': [] }
          values[r][mass]['obs'].append(CLs[i]/luminosity)
          values[r][mass]['exp'].append(expCLs[2]/luminosity)
          values[r][mass]['PEs'] += [e/luminosity for e in PE_CLs]        

    mass_list = sorted(list(set(mass_list)))
    thisobserved = ROOT.TGraph()
    thisexpected = ROOT.TGraph()
    thisexpected_plus1  = ROOT.TGraph()
    thisexpected_minus1 = ROOT.TGraph()
    thisexpected_plus2  = ROOT.TGraph()
    thisexpected_minus2 = ROOT.TGraph()
    for m in mass_list :
      if m not in values[r]: continue
      expCLs = GetCenterAndSigmaDeviations(values[r][m]['PEs'])
      print r, m, values[r][m]['obs'][0], values[r][m]['exp'][0], len(values[r][m]['PEs'])
      thisobserved.SetPoint(       thisobserved.GetN(),m,values[r][m]['obs'][0])
      thisexpected_minus2.SetPoint(thisexpected_minus2.GetN(),m,expCLs[0])
      thisexpected_minus1.SetPoint(thisexpected_minus1.GetN(),m,expCLs[1])
      thisexpected.SetPoint(       thisexpected.GetN(),m,expCLs[2])
      thisexpected_plus1.SetPoint( thisexpected_plus1.GetN(),m,expCLs[3])
      thisexpected_plus2.SetPoint( thisexpected_plus2.GetN(),m,expCLs[4])

    allobserved.append(thisobserved)
    allexpected.append(thisexpected)
    
    thisexpected1Sigma = makeBand(thisexpected_minus1,thisexpected_plus1)
    thisexpected2Sigma = makeBand(thisexpected_minus2,thisexpected_plus2)
    
    if r == 0:
    	allexpected1Sigma.append(thisexpected1Sigma)
    	allexpected2Sigma.append(thisexpected2Sigma)
    else:
        allexpected1Sigma.append(ROOT.TGraph())
        allexpected2Sigma.append(ROOT.TGraph())
    results[r] = {'obs':thisobserved, 'exp': thisexpected,'exp1':  thisexpected1Sigma,'exp2': thisexpected2Sigma}
    
  if makePlot:
    #print ratios
    #print [names['%1.2f'%r] for r in ratios]
    myPainter.drawSeveralObservedExpectedLimits(allobserved,allexpected,allexpected1Sigma,allexpected2Sigma,[names['%1.2f'%r] for r in ratios],outfolder+"GenericGaussians_"+dataset,"m_{G} [GeV]",\
     "#sigma #times #it{A} #times BR [pb]",luminosity,13,400,2000,1E-2,50,[],ATLASLabelLocation="BottomR",cutLocation="Left")
  
  return results

 
if __name__ == "__main__":

  #==========================
  #   User configurables
  #==========================

  makeCombinedPlot = False

  ratios_J100 = [0.0]#[0.0, 0.05, 0.07]
  ratios_J75  = [0.0]#[0.0, 0.05, 0.07, 0.10]

  indir = "/projects/hep/fs4/scratch/ATLAS/TLA/limits/"
  J100inputfileform = indir+"results/data2017/runSWIFT2016_J100yStar06/GenericGaussians_29p30_doSwift_*seed*"
  J75inputfileform  = indir+"results/data2017/runSWIFT2016_J75yStar03/GenericGaussians_3p57_doSwift_*"

  #==========================
  
  (dataset, luminosity, cutstring) = ("J100yStar06", 29.3*1000, "J100 |y*| < 0.6" )
  J100_plots = getGaussianLimits( J100inputfileform, ratios_J100, dataset, luminosity, cutstring, makePlot = not makeCombinedPlot)
   
  (dataset, luminosity, cutstring) = ("J75yStar03",  3.57*1000, "J75 |y*| < 0.3")
  J75_plots = getGaussianLimits( inputfileform, ratios_J75, dataset, luminosity, cutstring, makePlot = not makeCombinedPlot)
  
  '''
  if not makeCombinedPlot:
    (dataset, luminosity, cutstring) = ("J75yStar06", 3.57*1000, "J75 |y*| < 0.6")
    ratios_J75 = [0.0, 0.05, 0.07, 0.10]
    inputfileform = indir+"results/data2017/runSWIFT2016_J75yStar06/GenericGaussians_3p57_doSwift_*"
    
    J75_plots = getGaussianLimits( inputfileform, ratios_J75, dataset, luminosity, cutstring, makePlot = not makeCombinedPlot)
  '''

  if makeCombinedPlot:  
    allobserved = []
    allexpected = []
    allexpected1Sigma = []
    allexpected2Sigma = []
    
    all_ratios = list(set().union(ratios_J100,ratios_J75))
    all_ratios.sort()
    for r in all_ratios:
        SR_dicts = [J100_plots,J75_plots]
        allobserved.append([d[r]['obs'] if r in d else None for d in SR_dicts])
        allexpected.append([d[r]['exp'] if r in d else None for d in SR_dicts])
        allexpected1Sigma.append([d[r]['exp1'] if r in d else None for d in SR_dicts if r == 0])
        allexpected2Sigma.append([d[r]['exp2'] if r in d else None for d in SR_dicts if r == 0])
      
    # Initialize painter
    myPainter = Morisot()
    # Internal
    myPainter.setLabelType(2)
    myPainter.cutstring = 'J75 |y*| < 0.3, 3.57 fb^{-1}    J100 |y*| < 0.6, 29.3 fb^{-1}'
    outfolder = "./plots/"
    
    myPainter.drawSeveralObservedExpectedLimits(allobserved,allexpected,allexpected1Sigma,allexpected2Sigma,[names['%1.2f'%r] for r in all_ratios],outfolder+"GenericGaussians_Combined","m_{G} [GeV]",\
       "#sigma #times #it{A} #times BR [pb]",0,13,400,2000,1E-2,50,[],ATLASLabelLocation="byLegend",cutLocation="BottomL",labels=["J75 |y*| < 0.3, 3.57 fb^{-1}","J100 |y*| < 0.6, 29.3 fb^{-1} "])
