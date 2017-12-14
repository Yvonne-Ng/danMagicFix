#!/bin/python

#******************************************
#NOTE some of these functions are taken from CommonFunctions.py

#******************************************
#import
import ROOT
import re, sys, os, math

#******************************************
def printHist (hist):
    """print histogram content
    """
    for bin in xrange(hist.GetNbinsX()):
        print "bin %5s = %10.4f +/- %10.4f"%(bin, hist.GetBinContent(bin), hist.GetBinError(bin))

#******************************************
def checkHistPoissonianErrors (hist):
    """check if the errors of the histogram are Poissonian
    """
    for bin in xrange(hist.GetNbinsX()):
        if hist.GetBinContent(bin) != 0:
            if hist.GetBinError(bin) != math.sqrt(hist.GetBinContent(bin)):
                print "bin %5s = %10.4f +/- %10.4f\t->\t%10.4f"%(bin, hist.GetBinContent(bin), hist.GetBinError(bin), math.sqrt(hist.GetBinContent(bin)))
            else:
                print "bin %5s is ok"%bin
        else:
            print "bin %5s is empty"%bin

#******************************************
def setHistPoissonianErrors (hist, debug=False):
    """set Poissonian errors
    """
    for ii in xrange(hist.GetNbinsX()):
        hist.SetBinError(ii, math.sqrt(hist.GetBinContent(ii)))
        if debug:
            print "bin %5s = %10.4f +/- %10.4f"%(ii, hist.GetBinContent(ii), hist.GetBinError(ii))
    return hist

#******************************************
#NOTE this method does not take into account the sum of the weights
#NOTE use carefully
def simpleRoundHistogram (hist):
    """round the histogram bin content to integer numbers; set the bin error as Poissonian; NOTE this deos not take into account the sum of weights
    """
    newHist = ROOT.TH1D("newHist","newHist",hist.GetXaxis().GetNbins(),hist.GetXaxis().GetXbins().GetArray())
    newHist.SetDefaultSumw2(ROOT.kTRUE)
    for ii in xrange(newHist.GetNbinsX()):
        newHist.SetBinContent(ii, int(round(hist.GetBinContent(ii))))
        newHist.SetBinError(ii, math.sqrt(newHist.GetBinContent(ii)))
    return newHist

#******************************************
#NEW #NOTE this the new version of roundHistogram(); for the older version see simpleRoundHist()
def roundHistogram (hist):
    """round the histogram bin content to integer numbers; NOTE fill histogram one entry at a time
    """
    newHist = ROOT.TH1D("newHist","newHist",hist.GetXaxis().GetNbins(),hist.GetXaxis().GetXbins().GetArray())
    newHist.SetDefaultSumw2(ROOT.kTRUE)
    for ii in xrange(hist.GetNbinsX()):
        bincontent = int(round(hist.GetBinContent(ii)))
        for jj in xrange(bincontent):
            newHist.Fill(newHist.GetBinCenter(ii))
    return newHist

#******************************************
#function taken from CommonFunctions.py
#previously named getSample(name)
def getSignalMass(name):
    if "ExcitedQ" in name:
        regEx = re.compile("ExcitedQ(.*)Lambda")
    elif "BlackMax" in name:
        regEx = re.compile("MthMD(.*)\.e2303")
    else:
        regEx = re.compile("_14TeV(.*)_AU2")

    match = regEx.search(name)
    sample = match.group(1)
    return sample

#******************************************
#function taken from CommonFunctions.py
def getQStarXSec(mass):
    #cross sections are in nb
    sampleDict = {}
    sampleDict["2000"]  = 3.680e-02
    sampleDict["3000"]  = 2.985e-03
    sampleDict["4000"]  = 3.534e-04
    sampleDict["5000"]  = 4.833e-05
    sampleDict["6000"]  = 7.047e-06
    sampleDict["7000"]  = 1.039e-06
    sampleDict["8000"]  = 1.626e-07
    sampleDict["9000"]  = 3.186e-08
    sampleDict["10000"] = 9.205e-09
    sampleDict["11000"] = 3.653e-09
    sampleDict["12000"] = 1.684e-09
    sampleDict["13000"] = 8.456e-10

    #return cross section in nb
    #print 'QStar [%s GeV] cross section: %s nb'%(mass, str(sampleDict[mass]/1e3))
    return sampleDict[mass]

#******************************************
#function taken from CommonFunctions.py
def getQBHXSec(mass):
	#cross sections are in pb
    sampleDict = {}
    sampleDict["2000"] = 2.570e+03
    sampleDict["3000"] = 2.194e+02
    sampleDict["4000"] = 2.742e+01
    sampleDict["5000"] = 3.850e+00
    sampleDict["6000"] = 5.507e-01
    sampleDict["7000"] = 7.374e-02
    sampleDict["8000"] = 8.615e-03
    sampleDict["9000"] = 8.318e-04
    sampleDict["10000"] = 5.754e-05
    sampleDict["11000"] = 2.431e-06
    sampleDict["12000"] = 4.606e-08

    #return cross section in nb
    #print 'QBH [%s GeV] cross section: %s nb'%(mass, str(sampleDict[mass]/1e3))
    return sampleDict[mass]/1e3

#******************************************
#function taken from CommonFunctions.py
def getScaleFactor(model, mass, lumi, nEvents):
    #print "lumi = %s nb^-1"%lumi
    if model == "QStar":
        return (getQStarXSec(mass) * lumi) / nEvents
    elif model == "QBH":
        return (getQBHXSec(mass) * lumi) / nEvents
    else:
        raise SystemExit("\n***ERROR*** unknown signal sample: %s"%name)

#******************************************
#function taken from CommonFunctions.py
def getHistogramContainingPercentage(hist, percentage, lowBin=-999, highBin=-999):

    if lowBin == -999 : lowBin = 1
    if highBin == -999 : highBin = hist.GetNbinsX()

    nbins = hist.GetNbinsX()
    if (nbins==0) : return None

    integral = hist.Integral()
    if (integral==0) : return None

    firstBin = 1
    while (hist.GetBinContent(firstBin)==0) : firstBin = firstBin+1

    lastBin = nbins
    while (hist.GetBinContent(lastBin)==0) : lastBin = lastBin-1

    binsToScan = lastBin - firstBin
    if (binsToScan==0) : return None

    scanStep = binsToScan/100
    if (scanStep<1) : scanStep=1

    thisPercentage=0
    thisInterval=0

    remember1=1
    remember2=1

    #a really big interval to start with.
    smallestInterval= hist.GetBinLowEdge(lastBin) + hist.GetBinWidth(lastBin) - hist.GetBinLowEdge(firstBin) + 1e12 

    rememberThisPercentage=0

    for bin1 in xrange(firstBin, lastBin+1, scanStep) :
        for bin2 in xrange(bin1, lastBin+1, scanStep) :

            #print "bins ", bin1, " to ", bin2
            thisPercentage = (hist.Integral(bin1,bin2))/integral
            thisInterval = hist.GetBinLowEdge(bin2) +  hist.GetBinWidth(bin2) - hist.GetBinLowEdge(bin1)
            #print " -. contain ", thisPercentage

            if thisPercentage >= percentage :
                if thisInterval < smallestInterval :
                    remember1=bin1
                    remember2=bin2
                    smallestInterval=thisInterval
                    rememberThisPercentage=thisPercentage
                break
 
    #print "best interval are bins : ", remember1, " , ", remember2
    #print " which span ", smallestInterval, "GeV"
    lowBin=remember1
    highBin=remember2

    #now make the histogram
    choppedHist = hist.Clone()
    choppedHist.Reset("ICE")
    for iBin in xrange(0, hist.GetNbinsX()+1) :
        if iBin < lowBin : continue
        if iBin > highBin : continue
        choppedHist.SetBinContent(iBin, hist.GetBinContent(iBin))
        choppedHist.SetBinError(iBin, hist.GetBinError(iBin))

    #print "bins [" , lowBin, "," , highBin , "] contain " , rememberThisPercentage*100, "percent of the signal"
    return choppedHist

#******************************************
#function taken from CommonFunctions.py
def getEffectiveEntriesHistogram(hist):
    hee = ROOT.TH1D("hee","hee",hist.GetXaxis().GetNbins(),hist.GetXaxis().GetXbins().GetArray())
    for ii in xrange(hist.GetNbinsX()):
        if hist.GetBinError(ii) != 0.:
            nee = pow(hist.GetBinContent(ii), 2) /  pow(hist.GetBinError(ii), 2)
        else:
            nee = 0.
        hee.SetBinContent(ii, nee)
        hee.SetBinError(ii, math.sqrt(nee))
    return hee

#******************************************
def createPseudoDataHist(inHist):

    #nEvents = inHist.GetSumOfWeights()
    nEvents = inHist.Integral()
    #nEvents = inHist.GetEntries() * lumi

    print "number of events: %s"%nEvents
    random = ROOT.TRandom3(1986)
    gRandom = random

    normalizedMassDist = inHist.Clone("normalized_mjj")
    normalizedMassDist.Scale(1./normalizedMassDist.Integral())
    nEventsFluctuated = random.Poisson(floor(nEvents))
    pseudoData = inHist.Clone("pseudo_data")
    pseudoData.Reset()
    pseudoData.FillRandom(normalizedMassDist, int(nEventsFluctuated))

    return pseudoData
