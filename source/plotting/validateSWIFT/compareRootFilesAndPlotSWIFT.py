#!/usr/bin/env python

import ROOT
from art.morisot import Morisot
from array import array
import sys

#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
#import matplotlib as mpl

cases = ["D"]#,"B","C","D"]
windowWidths = ["13"]#[7,10,13,17]

testRootTemplate = "/afs/cern.ch/user/e/etolley/work/DijetTLA/DijetTLAStatisticalAnalysis_test/Bayesian/results/data2017/validateSWIFT/SearchResultData_{0}.root"
testPrintTemplate = "/afs/cern.ch/user/e/etolley/work/DijetTLA/DijetTLAStatisticalAnalysis_test/Bayesian/results/data2017/validateSWIFT/output_{0}"

dataName = "basicData"
testFitName = "basicBkgFrom4ParamFit"
testParGraphNames = "evolution_parameter{0}"
testWindows = "swiftBinsUsed"

validateTemplate = "/afs/cern.ch/user/e/etolley/work/DijetTLA/SWiFt/output/results/TLA_bkg_J75_yStar06_0_1x_{0}.root"

validateFitName = "BkgDefault"
#validateFitName = "Background_smooth0"
validateParGraphNames = "BG1func_p{0}"
validateWindows_low = "WindowLowEdges" # TVectorT<double>
validateWindows_high = "WindowHighEdges" # TVectorT<double>
validateWindows_centers = "WindowCenters"

folder = "plots/"

luminosity = "0.0"

#fig = plt.figure(figsize=(10,8))
#ax = fig.add_subplot(111)

plotstart = 0
plotend = 4200

                  # greens,  reds,      blues,     # yellows
fillcolours    = ['#CAFF70', '#FF6A6A', '#B0E2FF', '#FFF68F']
binedgecolours = ['#B3EE3A', '#CD5C5C', '#87CEFA', '#EEE685']
rangecolours   = ['#6E8B3D', '#FF3030', '#4F94CD', '#EEEE00']
centrecolours  = ['#556B2F', '#FF0000', '#36648B', '#8B7500']


# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Teals")
myPainter.setLabelType(2)
myPainter.setEPS(True)

# Useful
firstXBinPlot = 4
lastXBinPlot = 50

for case in cases :
  for windowWidth in windowWidths : 

    if "D" in case and windowWidth == 17 : continue

    mystring = "case{0}_window{1}".format(case,windowWidth)

    inProbeFile = ROOT.TFile.Open(testRootTemplate.format(mystring),"READ")
    dataHist = inProbeFile.Get(dataName)
    dataHist.SetDirectory(0)
    testHist = inProbeFile.Get(testFitName)
    testHist.SetDirectory(0)
    parGraphs = []
    for par in range(4) :
      testGraph = inProbeFile.Get(testParGraphNames.format(par))
      parGraphs.append(testGraph)
    testWindowsPlot = inProbeFile.Get(testWindows)
    testWindowsPlot.SetDirectory(0)
    inProbeFile.Close()

    validateFile = ROOT.TFile.Open(validateTemplate.format(mystring),"READ")
    validateHist = validateFile.Get(validateFitName)
    validateHist.SetDirectory(0)
    validateGraphs = []
    for par in range(4) :
      validateGraph = validateFile.Get(validateParGraphNames.format(par))
      validateGraphs.append(validateGraph)
    vec_validateWindows_low = validateFile.Get(validateWindows_low)
    vec_validateWindows_high = validateFile.Get(validateWindows_high)
    vec_validateWindows_center = validateFile.Get(validateWindows_centers)
    validateFile.Close()

    ratio = testHist.Clone("ratio_{0}".format(mystring))
    diff = testHist.Clone("diff_{0}".format(mystring))
    for bin in range(1,ratio.GetNbinsX()+1) :
      if validateHist.GetBinContent(bin) != 0 :
        ratioval = (testHist.GetBinContent(bin)-validateHist.GetBinContent(bin))/validateHist.GetBinContent(bin)
      else :
        ratioval = 0
      ratio.SetBinContent(bin,ratioval)
      diff.SetBinContent(bin,testHist.GetBinContent(bin) - validateHist.GetBinContent(bin))

    # Find range
    firstBin = 0
    lastBin = testHist.GetNbinsX()
    while (testHist.GetBinContent(firstBin)==0 and firstBin < testHist.GetNbinsX()) :
      firstBin+=1
    while (testHist.GetBinContent(lastBin)==0 and lastBin > 0) :
      lastBin-=1
    if (firstBin > lastBin) :
      firstBin=1
      lastBin = testHist.GetNbinsX()
    print "First bin = ",firstBin,": lower edge at",testHist.GetBinLowEdge(firstBin)
    print "Last bin = ",lastBin,": upper edge at",testHist.GetBinLowEdge(lastBin+1)

    myPainter.drawMultipleFitsAndResiduals(dataHist,[testHist,validateHist],[ratio],["New fit","SWIFT fit"],"m_{jj} [TeV]","Events",["(new-SWIFT)/SWIFT"],"{0}/compareOnData_{1}".format(folder,mystring),luminosity,13,firstBin,lastBin)

    myPainter.drawSignificanceHistAlone(ratio,"m_{jj} [TeV]","(new-SWIFT)/SWIFT","{0}/fitsRatio_{1}".format(folder,mystring),firstBin=firstBin,lastBin=lastBin)

    myPainter.drawSignificanceHistAlone(diff,"m_{jj} [TeV]","(new-SWIFT)","{0}/fitsDifference_{1}".format(folder,mystring),firstBin=firstBin,lastBin=lastBin)

    plotxlow = dataHist.GetBinLowEdge(firstXBinPlot)
    plotxhigh = dataHist.GetBinLowEdge(lastXBinPlot+1)

    # Compare parameter values
    for par in range(4) :
      test = parGraphs[par]
      val = validateGraphs[par]
      myPainter.drawSeveralObservedLimits([test,val],["Ours","SWIFT"],"{0}/parameterEvolution_par{1}_{2}".format(folder,par,mystring),"Bin estimated","Parameter value",luminosity,13,plotxlow,plotxhigh,"automatic","automatic",doLogY=False,doLogX=False,doRectangular=False,doLegendLocation="Right",ATLASLabelLocation="BottomL",isTomBeingDumb=False,addHorizontalLines=[],pairNeighbouringLines=False,cutLocation="Right")

    # Make plots of bins used in estimation
    theirs = testWindowsPlot.Clone("validateWindowRanges")
    theirs.Reset()
    for ybin in range(1,len(vec_validateWindows_center)+1) :
      xLow = vec_validateWindows_low[ybin-1]
      xHigh = vec_validateWindows_high[ybin-1]
      xEstimate = vec_validateWindows_center[ybin-1]
      print "Estimate for",xEstimate,"comes from",xLow,",",xHigh
      for xbin in range(1,theirs.GetNbinsX()+1) :
        if theirs.GetXaxis().GetBinCenter(xbin) < xLow or theirs.GetXaxis().GetBinCenter(xbin) > xHigh :
          theirs.SetBinContent(xbin,theirs.GetYaxis().FindBin(xEstimate),0.0)
        else :
          theirs.SetBinContent(xbin,theirs.GetYaxis().FindBin(xEstimate),2.0)

    theirs.SetContour(2)
    myPainter.drawOverlaid2DPlots(testWindowsPlot,[theirs],"{0}/windowsUsed_{1}".format(folder,mystring),"Window range",plotxlow,plotxhigh,"Bin estimated",plotxlow,plotxhigh,"",luminosity=-1,CME=-1,doRectangular=False)


