#!/usr/bin/env python

import ROOT
from art.morisot import Morisot
from array import array
import sys

def GetZVal (p, excess) :
  #the function normal_quantile converts a p-value into a significance,
  #i.e. the number of standard deviations corresponding to the right-tail of 
  #a Gaussian
  if excess :
    zval = ROOT.Math.normal_quantile(1-p,1);
  else :
    zval = ROOT.Math.normal_quantile(p,1);

  return zval

def MakeHistoFromStats(statistics) :

  nentries = len(statistics)
  nBins = int(float(nentries)/10.0)

  maxVal = max(statistics)
  minVal = min(statistics)
  axisrange = maxVal - minVal;

  thismin = minVal-0.05*axisrange;
  thismax = maxVal+0.05*axisrange;

  statPlot = ROOT.TH1D("statPlot","",nBins,thismin,thismax)
  for val in range(len(statistics)) :
    statPlot.Fill(statistics[val])

  return statPlot

# Get input
filename = "/cluster/warehouse/kpachal/TLA2016/StatisticalAnalysis/Bayesian/TestResults_dressRehearsal/noSignal_permitWindow/SearchResultNoSig_our3Par_seed1_noSig_from496_permitWindow.root"

folder = "TestAnomalousPoints"
ext = "_our3Par_496"

doAlternate = False
doStatSearch = False

luminosity = 3500

# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Teals")
myPainter.setLabelType(2)
myPainter.setEPS(True)

searchInputFile = ROOT.TFile.Open(filename,"READ")

# Retrieve search phase inputs
basicData = searchInputFile.Get("basicData")
normalizedData = searchInputFile.Get("normalizedData")
basicBkgFrom4ParamFit = searchInputFile.Get("basicBkgFrom4ParamFit")
normalizedBkgFrom4ParamFit = searchInputFile.Get("normalizedBkgFrom4ParamFit")
residualHist = searchInputFile.Get("residualHist")
relativeDiffHist = searchInputFile.Get("relativeDiffHist")
sigOfDiffHist = searchInputFile.Get("sigOfDiffHist")
logLikelihoodPseudoStatHist = searchInputFile.Get("logLikelihoodStatHistNullCase")
chi2PseudoStatHist = searchInputFile.Get("chi2StatHistNullCase")
bumpHunterStatHist = searchInputFile.Get("bumpHunterStatHistNullCase")
theFitFunction = searchInputFile.Get('theFitFunction')
bumpHunterTomographyPlot = searchInputFile.Get('bumpHunterTomographyFromPseudoexperiments')
bumpHunterStatOfFitToData = searchInputFile.Get('bumpHunterStatOfFitToData')

nomPlus1 = searchInputFile.Get("nominalBkgFromFit_plus1Sigma")
nomMinus1 = searchInputFile.Get("nominalBkgFromFit_minus1Sigma")

if doAlternate :
  alternateBkg = searchInputFile.Get("alternateFitOnRealData")
  nomWithNewFuncErrSymm = searchInputFile.Get("nomOnDataWithSymmetricRMSScaleFuncChoiceErr")
  valueNewFuncErrDirected = searchInputFile.Get("nomOnDataWithDirectedRMSScaleFuncChoiceErr")

logLOfFitToDataVec = searchInputFile.Get('logLOfFitToData')
chi2OfFitToDataVec = searchInputFile.Get('chi2OfFitToData')
statOfFitToData = searchInputFile.Get('bumpHunterPLowHigh')
logLOfFitToData = logLOfFitToDataVec[0]
logLPVal = logLOfFitToDataVec[1]
chi2OfFitToData = chi2OfFitToDataVec[0]
chi2PVal = chi2OfFitToDataVec[1]
bumpHunterStatFitToData = statOfFitToData[0]
bumpHunterPVal = bumpHunterStatOfFitToData[1]
bumpLowEdge = statOfFitToData[1]
bumpHighEdge = statOfFitToData[2]

excludeWindowNums = searchInputFile.Get('excludeWindowNums')
excludeWindow = int(excludeWindowNums[0]+0.5)
bottomWindowEdge = excludeWindowNums[1]
topWindowEdge = excludeWindowNums[2]

NDF = searchInputFile.Get('NDF')[0]

fitparams = searchInputFile.Get('fittedParameters')

print "logL of fit to data is",logLOfFitToData
print "logL pvalue is",logLPVal
print "chi2 of fit to data is",chi2OfFitToData
print "NDF is",NDF
print "chi2/NDF is",chi2OfFitToData/NDF
print "chi2 pvalue is",chi2PVal
print "bump hunter stat of fit to data is",bumpHunterStatFitToData
print "bumpLowEdge, bumpHighEdge are",bumpLowEdge,bumpHighEdge
print "BumpHunter pvalue is",bumpHunterPVal
print "which is Z value of",GetZVal(bumpHunterPVal,True)

print "Fitted parameters were:",fitparams

# Find range
firstBin = 0
lastBin = basicData.GetNbinsX()
while (basicData.GetBinContent(firstBin)==0 and firstBin < basicData.GetNbinsX()) :
  firstBin+=1
while (basicData.GetBinContent(lastBin-1)==0 and lastBin > 0) :
  lastBin-=1
if (firstBin > lastBin) :
  firstBin=1
  lastBin = basicData.GetNbinsX()
print "First bin = ",firstBin,": lower edge at",basicData.GetBinLowEdge(firstBin)

# Calculate from fit range
fitRange = searchInputFile.Get("FitRange")
firstBin = basicData.FindBin(fitRange[0])-1
lastBin = basicData.FindBin(fitRange[1])+2

# Convert plots into desired final form
nbins = basicData.GetNbinsX()
newbins = [basicData.GetBinLowEdge(1)] # []
for np in range(nbins+1) :
  newbins.append((basicData.GetBinLowEdge(np)+basicData.GetBinWidth(np))/1000)

print "New bins are",newbins

# Make never versions of old plots
newbasicdata = ROOT.TH1D("basicData_TeV","basicData_TeV",len(newbins)-1,array('d',newbins))
newbasicBkgFrom4ParamFit = ROOT.TH1D("basicBkgFrom4ParamFit_TeV","basicBkgFrom4ParamFit_TeV",len(newbins)-1,array('d',newbins))
newresidualHist = ROOT.TH1D("residualHist_TeV","residualHist_TeV",len(newbins)-1,array('d',newbins))
newrelativeDiffHist = ROOT.TH1D("relativeDiffHist_TeV","relativeDiffHist_TeV",len(newbins)-1,array('d',newbins))
newsigOfDiffHist = ROOT.TH1D("sigOfDiffHist_TeV","sigOfDiffHist_TeV",len(newbins)-1,array('d',newbins))

newAlternateBkg = ROOT.TH1D("alternateBkg_TeV","alternateBkg_TeV",len(newbins)-1,array('d',newbins))
newNomPlus1= ROOT.TH1D("nomPlus1_TeV","nomPlus1_TeV",len(newbins)-1,array('d',newbins))
newNomMinus1= ROOT.TH1D("nomMinus1_TeV","nomMinus1_TeV",len(newbins)-1,array('d',newbins))
newnomWithNewFuncErrSymm = ROOT.TH1D("nomWithNewFuncErrSymm_TeV","nomWithNewFuncErrSymm_TeV",len(newbins)-1,array('d',newbins))
newValueNewFuncErrDirected= ROOT.TH1D("nomWithNewFuncErrDirected_TeV","nomWithNewFuncErrDirected_TeV",len(newbins)-1,array('d',newbins))

for histnew,histold in [[newbasicdata,basicData],[newbasicBkgFrom4ParamFit,basicBkgFrom4ParamFit],\
        [newresidualHist,residualHist],[newrelativeDiffHist,relativeDiffHist],[newsigOfDiffHist,sigOfDiffHist]] :
  for bin in range(histnew.GetNbinsX()+2) :
    valToUse = histold.GetBinCenter(bin)/1000.0
#    histnew.SetBinContent(bin,histold.GetBinContent(bin))
#    histnew.SetBinError(bin,histold.GetBinError(bin))
    histnew.SetBinContent(histnew.FindBin(valToUse),histold.GetBinContent(bin))
    histnew.SetBinError(histnew.FindBin(valToUse),histold.GetBinError(bin))
 
if doAlternate :
  for histnew,histold in [[newAlternateBkg,alternateBkg],[newNomPlus1,nomPlus1],[newNomMinus1,nomMinus1],\
        [newnomWithNewFuncErrSymm,nomWithNewFuncErrSymm],[newValueNewFuncErrDirected,valueNewFuncErrDirected]] :
    for bin in range(histnew.GetNbinsX()+2) :
      histnew.SetBinContent(bin,histold.GetBinContent(bin))
      histnew.SetBinError(bin,histold.GetBinError(bin))

# Search phase plots
print "Writing window limits at",bottomWindowEdge,topWindowEdge
myPainter.drawDataAndFitOverSignificanceHist(basicData,basicBkgFrom4ParamFit,residualHist,\
        'm_{jj} [TeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
        luminosity,13,fitRange[0],fitRange[1],firstBin,lastBin,True,bumpLowEdge,bumpHighEdge,doWindowLimits=excludeWindow,windowLow=bottomWindowEdge,windowHigh=topWindowEdge)
myPainter.drawPseudoExperimentsWithObservedStat(logLikelihoodPseudoStatHist,float(logLOfFitToData),logLPVal,0,luminosity,13,\
        'logL statistic','Pseudo-exeperiments',"{0}/logLStatPlot".format(folder)+ext)
myPainter.drawPseudoExperimentsWithObservedStat(chi2PseudoStatHist,float(chi2OfFitToData),chi2PVal,0,luminosity,13,\
        "#chi^{2}",'Pseudo-exeperiments',"{0}/chi2StatPlot".format(folder)+ext)
myPainter.drawPseudoExperimentsWithObservedStat(bumpHunterStatHist,float(bumpHunterStatFitToData),bumpHunterPVal,0,luminosity,13,\
        'BumpHunter','Pseudo-exeperiments',"{0}/bumpHunterStatPlot".format(folder)+ext)
myPainter.drawBumpHunterTomographyPlot(bumpHunterTomographyPlot,"{0}/bumpHunterTomographyPlot".format(folder)+ext)

# Various significance plots
myPainter.drawSignificanceHistAlone(newrelativeDiffHist,"m_{jj} [TeV]","(D - B)/B","{0}/significanceonlyplot".format(folder)+ext)
myPainter.drawSignificanceHistAlone(newsigOfDiffHist,"m_{jj} [TeV]","(D - B)/#sqrt{Derr^{2}+Berr^{2}}","{0}/sigofdiffonlyplot".format(folder)+ext)

if excludeWindow :
  pValuesRemainderOfWindow = searchInputFile.Get("BHLogLAndChi2OfRemainderAfterWindow")
  BHOfRemainder = pValuesRemainderOfWindow[0]
  LogLOfRemainder = pValuesRemainderOfWindow[1]
  Chi2OfRemainder = pValuesRemainderOfWindow[2]
  print "Window was excluded!"
  print "After window exclusion, p-values of remaining spectrum were:"
  print "\tChi2 p-value:\t",Chi2OfRemainder
  print "\tBH p-value:\t",BHOfRemainder
  print "\tLogL p-value:\t",LogLOfRemainder,"\n"

print "Final values were:"
print "\tChi2 p-value:\t",chi2PVal
print "\tBH p-value:\t",bumpHunterPVal
print "\tLogL p-value:\t",logLPVal

if doAlternate :

  # Now to make the one comparing uncertainties
  placeHolderNom = newbasicBkgFrom4ParamFit.Clone()
  placeHolderNom.SetName("placeHolderNom")
  nomPlusSymmFuncErr = newbasicBkgFrom4ParamFit.Clone()
  nomPlusSymmFuncErr.SetName("nomPlusNewFuncErr")
  nomMinusSymmFuncErr = newbasicBkgFrom4ParamFit.Clone()
  nomMinusSymmFuncErr.SetName("nomMinusNewFuncErr")
  for bin in range(nomPlusSymmFuncErr.GetNbinsX()+2) :
    nomPlusSymmFuncErr.SetBinContent(bin,newnomWithNewFuncErrSymm.GetBinContent(bin) + newnomWithNewFuncErrSymm.GetBinError(bin))
    nomMinusSymmFuncErr.SetBinContent(bin,newnomWithNewFuncErrSymm.GetBinContent(bin) - newnomWithNewFuncErrSymm.GetBinError(bin))

  myPainter.drawDataWithFitAsHistogram(newbasicdata,newbasicBkgFrom4ParamFit,luminosity,13,"m_{jj} [TeV]","Events",["Data","Fit","Fit uncertainty","Function choice"],"{0}/compareFitQualityAndAlternateFit".format(folder)+ext,True,[[newNomPlus1,newNomMinus1],[placeHolderNom,newAlternateBkg]],firstBin,lastBin+2,True,True,True,False,False)
  myPainter.drawDataWithFitAsHistogram(newbasicdata,newbasicBkgFrom4ParamFit,luminosity,13,"m_{jj} [TeV]","Events",["Data","Fit","Function choice","Alternate function"],"{0}/compareFitChoiceAndAlternateFit".format(folder)+ext,True,[[nomPlusSymmFuncErr,nomMinusSymmFuncErr],[placeHolderNom,newAlternateBkg]],firstBin,lastBin+2,True,True,True,False,False)
  myPainter.drawDataWithFitAsHistogram(newbasicdata,newbasicBkgFrom4ParamFit,luminosity,13,"m_{jj} [TeV]","Events",["Data","Fit","Fit uncertainty","Function choice"],"{0}/compareFitQualtiyAndFitChoice".format(folder)+ext,True,[[newNomPlus1,newNomMinus1],[nomPlusSymmFuncErr,nomMinusSymmFuncErr]],firstBin,lastBin+2,True,True,True,False,False)
  myPainter.drawDataWithFitAsHistogram(newbasicdata,newbasicBkgFrom4ParamFit,luminosity,13,"m_{jj} [TeV]","Events",["Data","Fit","Fit uncertainty","RMS times sign"],"{0}/compareFitQualityAndFitChoice_Asymm".format(folder)+ext,True,[[newNomPlus1,newNomMinus1],[placeHolderNom,newValueNewFuncErrDirected]],firstBin,lastBin+2,True,True,True,False,False)

  # Make a ratio histogram for Sasha's plot.
  altFitRatio = ROOT.TH1D("altFitRatio","altFitRatio",len(newbins)-1,array('d',newbins))
  for bin in range(0,altFitRatio.GetNbinsX()+1) :
    if newbasicBkgFrom4ParamFit.GetBinContent(bin) == 0 :
      altFitRatio.SetBinContent(bin,0)
    else :
      altFitRatio.SetBinContent(bin,(valueNewFuncErrDirected.GetBinContent(bin)-newbasicBkgFrom4ParamFit.GetBinContent(bin))/newbasicBkgFrom4ParamFit.GetBinContent(bin))
  # Make a ratio histogram for paper plot.
  PlusNomRatio = ROOT.TH1D("PlusNomRatio","PlusNomRatio",len(newbins)-1,array('d',newbins))
  for bin in range(0,PlusNomRatio.GetNbinsX()+1) :
    if newbasicBkgFrom4ParamFit.GetBinContent(bin) == 0 :
      PlusNomRatio.SetBinContent(bin,0)
    else :
      PlusNomRatio.SetBinContent(bin,(newNomPlus1.GetBinContent(bin)-newbasicBkgFrom4ParamFit.GetBinContent(bin))/newbasicBkgFrom4ParamFit.GetBinContent(bin))
  MinusNomRatio = ROOT.TH1D("MinusNomRatio","MinusNomRatio",len(newbins)-1,array('d',newbins))
  for bin in range(0,MinusNomRatio.GetNbinsX()+1) :
    if newbasicBkgFrom4ParamFit.GetBinContent(bin) == 0 :
      MinusNomRatio.SetBinContent(bin,0)
    else :
      MinusNomRatio.SetBinContent(bin,(newNomMinus1.GetBinContent(bin)-newbasicBkgFrom4ParamFit.GetBinContent(bin))/newbasicBkgFrom4ParamFit.GetBinContent(bin))

  myPainter.drawMultipleFitsAndResiduals(newbasicdata,[newbasicBkgFrom4ParamFit,newValueNewFuncErrDirected],[altFitRatio],["Nominal fit","Func choice unc"],"m_{jj} [TeV]","Events",["(alt-nom)/nom"],"{0}/directedFuncChoiceVersusNominal_withRatio".format(folder)+ext,luminosity,13,firstBin,lastBin+2)

  myPainter.drawDataWithFitAsHistogramAndResidual(newbasicdata,newbasicBkgFrom4ParamFit,luminosity,13,"m_{jj} [TeV]","Events",["Data","Fit","Statistical uncertainty on fit","Function choice"],"{0}/compareFitQualityAndFitChoice_Asymm_WithRatio".format(folder),True,[[newNomPlus1,newNomMinus1],[placeHolderNom,newValueNewFuncErrDirected]],[altFitRatio,MinusNomRatio,PlusNomRatio],firstBin,lastBin+2,True,True,True,False,False) # changed from lastBin+2 to lastBin+18 to match FancyFigure

  myPainter.drawDataWithFitAsHistogramAndResidualPaper(newbasicdata,newbasicBkgFrom4ParamFit,luminosity,13,"m_{jj} [TeV]","Events",["Data","Fit","Statistical uncertainty on fit","Function choice"],"{0}/compareFitQualityAndFitChoice_Asymm_WithRatioPaper".format(folder),True,[[newNomPlus1,newNomMinus1],[placeHolderNom,newValueNewFuncErrDirected]],[altFitRatio,MinusNomRatio,PlusNomRatio],firstBin,lastBin+2,True,True,True,False,False) # changed from lastBin+2 to lastBin+18 to match FancyFigure

if doStatSearch :
  TomographyPlotWithStats = searchInputFile.Get("TomographyPlotWithStats")
  bumpHunterStatsWSyst = searchInputFile.Get("bumpHunterStatsWSyst")
  bumpHunterPValVec = searchInputFile.Get("bumpHunterPValErrWithStats")
  bumpStat = bumpHunterPValVec[0]
  bumpPVal = bumpHunterPValVec[1]
  myPainter.drawPseudoExperimentsWithObservedStat(bumpHunterStatsWSyst,bumpStat,bumpPVal,0,luminosity,13,\
          'BumpHunter','Pseudo-experiments',"{0}/bumpHunterStatPlot_withUncertainties".format(folder)+ext)
  myPainter.drawBumpHunterTomographyPlot(TomographyPlotWithStats,"{0}/bumpHunterTomographyPlot_withStats".format(folder)+ext)

  print "BumpHunter pvalue is",bumpPVal

searchInputFile.Close()

print "Done."
