#!/usr/bin/env python
import os
import ROOT
from art.morisot import Morisot

# get bash environment variables into python
STATS12PATH=str(os.environ["STATS12PATH"])

# Get input (the rootfile name should match that specified in the SearchPhase.config file)
searchInputFile = ROOT.TFile(STATS12PATH+'/results/SearchPhase_results.root')
folderextension = ''

# Initialize painter
myPainter = Morisot()

# Define necessary quantities.
luminosity = 20300

# Retrieve search phase inputs
basicData = searchInputFile.Get("basicData")
basicData.Print("all")
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

logLOfFitToData = searchInputFile.Get('logLOfFitToData')[0]
chi2OfFitToData = searchInputFile.Get('chi2OfFitToData')[0]
statOfFitToData = searchInputFile.Get('bumpHunterPLowHigh')
bumpHunterPvalFitToData = statOfFitToData[0]
bumpLowEdge = statOfFitToData[1]
bumpHighEdge = statOfFitToData[2]

fitparams = searchInputFile.Get('fittedParameters')

print "logL of fit to data is",logLOfFitToData
print "chi2 of fit to data is",chi2OfFitToData
print "bump hunter stat of fit to data is",bumpHunterPvalFitToData
print "bumpLowEdge, bumpHighEdge are",bumpLowEdge,bumpHighEdge

# Search phase plots
myPainter.drawDataAndFitOverSignificanceHist(basicData,basicBkgFrom4ParamFit,residualHist,\
                              'Mass[GeV]','Events','Significance',folderextension+'figure1',luminosity,8,250,6000,True,bumpLowEdge,bumpHighEdge)

myPainter.drawDataAndFitOverSignificanceHist(basicData,basicBkgFrom4ParamFit,residualHist,\
                              'Mass[GeV]','Events','Significance',folderextension+'figure1_nobump',luminosity,8,200,5000,False,bumpLowEdge,bumpHighEdge)

myPainter.drawPseudoExperimentsWithObservedStat(logLikelihoodPseudoStatHist,float(logLOfFitToData),luminosity,8,\
                              'logL statistic','Pseudo-exeperiments',folderextension+"logLStatPlot")
myPainter.drawPseudoExperimentsWithObservedStat(chi2PseudoStatHist,float(chi2OfFitToData),luminosity,8,\
                              "#chi^{2}",'Pseudo-exeperiments',folderextension+"chi2StatPlot")
myPainter.drawPseudoExperimentsWithObservedStat(bumpHunterStatHist,float(bumpHunterPvalFitToData),luminosity,8,\
                              'BumpHunter','Pseudo-exeperiments',folderextension+"bumpHunterStatPlot")

myPainter.drawBumpHunterTomographyPlot(bumpHunterTomographyPlot,folderextension+"bumpHunterTomographyPlot")

