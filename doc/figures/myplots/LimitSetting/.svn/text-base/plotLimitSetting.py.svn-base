#!/usr/bin/env python

import ROOT
import os,sys
from art.morisot import Morisot
from array import array

# get bash environment variables into python
STATS12PATH=str(os.environ["STATS12PATH"])
DATA=str(os.environ["DATA"])

luminosity = 20300

# Will run for each signal type you include here.
signalInputFileTypes = ["QStar"]
#signalInputFileTypes = ["QStar","S8","BlackMax","WPrime","QBH","WStarSinX0"]#,"WStar"]
#signalInputFileTypes = ["QStar","S8","BlackMax","WPrime","QBH"]

SignalTitles = {"QStar": "#it{q}*"
                , "WPrime": "#it{W}'"
                , "S8":"#it{s}8"
                , "WStarSinX0" : "Leptophobic #it{W}* (sin#phi_{X}=0)"
                , "WStarSinX1": "Leptophilic #it{W}* (sin#phi_{X}=1)"
                , "QBH": "QBH, QBH gen."
                , "BlackMax": "QBH, BlackMax gen."
                }

SignalAxes = {"QStar": {"X" : "m_{#it{q}*} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "BlackMax": {"X" : "m_{th} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "QBH": {"X" : "m_{th} [GeV]", "Y": "#sigma #times #it{A}  [pb]"},
              "S8": {"X" : "m_{#it{s}8} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "WPrime": {"X" : "m_{#it{W}'} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "WStarSinX0": {"X" : "m_{#it{W}*} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "WStarSinX1": {"X" : "m_{#it{W}*}} [GeV]", "Y": "#sigma #times #it{A} [pb]"}
             }

yranges = {}
yranges['BlackMax'] = [1E-4,100] # was 10
yranges['QBH'] = [1E-4,100] # was 10
yranges['QStar'] = [5E-4,5E3] # was 1000
yranges['S8'] = [5E-3,5E3] # was 1000
yranges['WPrime'] = [5E-3,1E3]
yranges['WStarSinX0'] = [1E-3,100]
yranges['WStarSinX1'] = [1E-3,100]

# Get input
searchInputFile = ROOT.TFile(STATS12PATH+'/results/SearchPhase_results.root')
#limitInputFile = ROOT.TFile('/home/pachal/oxfordDijets/statisticsCode/StatisticalAnalysis2012/results/LimitSettingPhase_results.root')
xseceffInputFile = ROOT.TFile(DATA+'/inputs/xsecandacceptance/NormalizeTemplates/CrossSectionsForPlotting.root')
templateInputFile = ROOT.TFile(DATA+'/inputs/xsecandacceptance/NormalizeTemplates/TemplatesForPlotting.root')

folderextension = ""

# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("ATLAS")
myPainter.setLabelType(0)

for signal in signalInputFileTypes :
  # My files
#  limitFileName = "/home/pachal/oxfordDijets/statisticsCode/StatisticalAnalysis2012/results/{0}/LimitSettingPhase_results.root".format(signal)
  # Caterina's files
#  limitFileName = "/home/pachal/oxfordDijets/statisticsCode/StatisticalAnalysis2012/results/LimitSettingPhase{0}_Caterina_results.root".format(signal)
  # Oliver's files
#  limitFileName = "/data/atlas/atlasdata/pachal/fromOliver/LimitSettingPhase_results_{0}.root".format(signal)
  # "Final" files
  limitFileName = STATS12PATH+"/results/{0}/LimitSettingPhase_results.root".format(signal)
  limitInputFile = ROOT.TFile(limitFileName)

  axisnames = SignalAxes[signal]
  xname = axisnames["X"]
  yname = axisnames["Y"]

  # Define necessary quantities.
  signalName = signal+"_graph"

  # Limit-setting plots
  outputName = folderextension+"Limits_"+signal
  observedGraph = limitInputFile.Get('observedXSecAccVersusMass')
  expectedGraph1Sigma = limitInputFile.Get('expectedXSecAccVersusMass_oneSigma')
  expectedGraph2Sigma = limitInputFile.Get('expectedXSecAccVersusMass_twoSigma')
  signalGraph = xseceffInputFile.Get(signalName)
  extrasignalGraph = 0

  if signal=="BlackMax" :
    extrasignalGraph = xseceffInputFile.Get("QBH_graph")
  elif signal=="WStarSinX0" :
    extrasignalGraph = xseceffInputFile.Get("WStarSinX1_graph")

  # Make everything shifted by 1000 for TeV plots
  shiftedgraphs = []
  d1, d2 = ROOT.Double(0), ROOT.Double(0)
  for graph in [observedGraph,expectedGraph1Sigma,expectedGraph2Sigma,signalGraph,extrasignalGraph] :
    if graph==0 :
      continue
    newgraph = graph.Clone()
    newgraph.SetName(graph.GetName()+"_scaled")
    for np in range(newgraph.GetN()) :
      newgraph.GetPoint(np,d1,d2)
      newgraph.SetPoint(np,d1/1000.0,d2)
    shiftedgraphs.append(newgraph)


  thisx, thisy = ROOT.Double(0), ROOT.Double(0)

  thisyrange = yranges[signal]

  observedGraph.GetPoint(0,thisx,thisy)
  xlow = thisx - 100
  observedGraph.GetPoint(observedGraph.GetN()-1,thisx, thisy)
  xhigh = thisx + 100

  # make regular plots
  thisname = outputName+"_GeV"
  [obs,exp] = myPainter.drawLimitSettingPlot2Sigma(observedGraph,expectedGraph1Sigma,expectedGraph2Sigma,signalGraph,SignalTitles[signal],\
     thisname,xname,yname,luminosity,8,xlow,xhigh,thisyrange[0],thisyrange[1],False)

  # make TeV plots
  newxname = xname.replace("G","T")
  [obs,exp] = myPainter.drawLimitSettingPlot2Sigma(shiftedgraphs[0],shiftedgraphs[1],shiftedgraphs[2],shiftedgraphs[3],SignalTitles[signal],\
     outputName,newxname,yname,luminosity,8,xlow/1000,xhigh/1000,thisyrange[0],thisyrange[1],False)

  if signal=="BlackMax" :
    outputName = outputName+"_withQBH"
    thisname = thisname+"_withQBH"
    [[obs1,exp1],[obs2,exp2]] = myPainter.drawLimitSettingPlot2Sigma(observedGraph,expectedGraph1Sigma,expectedGraph2Sigma,\
       [signalGraph,extrasignalGraph],[SignalTitles[signal],SignalTitles["QBH"]],thisname,xname,yname,luminosity,8,\
       xlow,xhigh,thisyrange[0],thisyrange[1],False)
    [[obs1,exp1],[obs2,exp2]] = myPainter.drawLimitSettingPlot2Sigma(shiftedgraphs[0],shiftedgraphs[1],shiftedgraphs[2],\
       [shiftedgraphs[3],shiftedgraphs[4]],[SignalTitles[signal],SignalTitles["QBH"]],outputName,xname,yname,luminosity,8,\
       xlow/1000,xhigh/1000,thisyrange[0],thisyrange[1],False)

  if signal=="WStarSinX0" :
    outputName = outputName+"_withLeptophilic"
    thisname = thisname+"_withLeptophilic"
    [[obs1,exp1],[obs2,exp2]] = myPainter.drawLimitSettingPlot2Sigma(observedGraph,expectedGraph1Sigma,expectedGraph2Sigma,\
       [signalGraph,extrasignalGraph],[SignalTitles[signal],SignalTitles["WStarSinX1"]],thisname,newxname,yname,luminosity,8,\
       xlow,xhigh,thisyrange[0],thisyrange[1],False)
    [[obs1,exp1],[obs2,exp2]] = myPainter.drawLimitSettingPlot2Sigma(shiftedgraphs[0],shiftedgraphs[1],shiftedgraphs[2],\
       [shiftedgraphs[3],shiftedgraphs[4]],[SignalTitles[signal],SignalTitles["WStarSinX1"]],outputName,newxname,yname,luminosity,8,\
       xlow/1000,xhigh/1000,thisyrange[0],thisyrange[1],False)

  print "-"*50
  print "Signal",signal
  print "Observed limit at 95% CL:",obs
  print "Expected limit at 95% CL:",exp,"\n"

  if signal=="BlackMax" or signal=="WStarSinX0" :
    print "Second signal curve: "
    print "Observed limit from this plot:",obs2
    print "Expected limit from this plot:",exp2

# Draw signal overlay plots
SignalMasses = {"QStar": [600, 2000, 3500],
                "WPrime": [600, 1600, 2600],
                "S8":[600, 2000, 3500],
                "WStarSinX0" : [500, 1000, 2400], 
                "WStarSinX1": [500, 1000, 2500], 
                "QBH": [1000, 2000, 3500, 5000], 
                "BlackMax": [1000, 2000, 3500, 5000], 
                }
SignalScalingFactors = {"QStar": 3,
                "WPrime": 100,
                "S8":10,
                "WStarSinX0" : 100,
                "WStarSinX1": 100, 
                "QBH": 1, 
                "BlackMax": 1}

# Get basic hists
datahist = searchInputFile.Get("basicData")
datahist.SetDirectory(0)
fithist = searchInputFile.Get("basicBkgFrom4ParamFit")
fithist.SetDirectory(0)
basicSignificancePlot = searchInputFile.Get("relativeDiffHist")
basicSignificancePlot.SetDirectory(0)
residual = searchInputFile.Get("residualHist")
residual.SetDirectory(0)

# Make scaled versions
standardbins = datahist.GetXaxis().GetXbins()
newbins = []#ROOT.TArrayD(standardbins.GetSize())
for np in range(standardbins.GetSize()) :
  newbins.append(standardbins[np]/1000)
newdatahist = ROOT.TH1D("basicData_TeV","basicData_TeV",len(newbins)-1,array('d',newbins))
newfithist = ROOT.TH1D("basicBkgFrom4ParamFit_TeV","basicBkgFrom4ParamFit_TeV",len(newbins)-1,array('d',newbins))
newbasicSignificancePlot = ROOT.TH1D("relativeDiffHist_TeV","relativeDiffHist_TeV",len(newbins)-1,array('d',newbins))
newresidual = ROOT.TH1D("residualHist_TeV","residualHist_TeV",len(newbins)-1,array('d',newbins))
for histnew,histold in [[newdatahist,datahist],[newfithist,fithist],[newbasicSignificancePlot,basicSignificancePlot],[newresidual,residual]] :
  for bin in range(histnew.GetNbinsX()+2) :
    histnew.SetBinContent(bin,histold.GetBinContent(bin))
    histnew.SetBinError(bin,histold.GetBinError(bin))
newsigtemplate = ROOT.TH1D("basicsignal_TeV","basicsignal_TeV",len(newbins)-1,array('d',newbins))

# Find range
firstBin =0
while (fithist.GetBinContent(firstBin)<1 and firstBin < fithist.GetNbinsX()) :
  firstBin+=1
lastBin = fithist.GetNbinsX()+1
while (fithist.GetBinContent(lastBin-1)==0 and lastBin > 0) :
  lastBin-=1
if (firstBin > lastBin) :
  firstBin=1
  lastBin = fithist.GetNbinsX()

for signal in signalInputFileTypes :
  print "in signal",signal

  signalMasses = SignalMasses[signal][:]
  signalMassesTeV = SignalMasses[signal][:]
  for index in range(len(SignalMasses[signal])) :
    signalMassesTeV[index] = signalMasses[index]/1000.0
  signalPlots = []
  signalPlotsTeV = []
  sigratioPlots = []
  sigratioPlotsTeV = []
  legendlist = []
  legendlistTeV = []
  for mass in signalMasses :
    sigplot = templateInputFile.Get(signal+"_%d" % mass)
    sigplot.SetDirectory(0)

    sigplottev = newsigtemplate.Clone()
    sigplottev.SetName("sigplot_{0}_{1}_TeV".format(signal,mass))
    for bin in range(sigplottev.GetNbinsX()+2) :
      sigplottev.SetBinContent(bin,sigplot.GetBinContent(bin))
      sigplottev.SetBinError(bin,sigplot.GetBinError(bin))

    index = 0
    for thissigplot,thissuffix in [[sigplot,""],[sigplottev,"_TeV"]] :

      # Normalise to correct amount
      thissigplot.Scale(luminosity)

      sigplotforfitplusbkg = thissigplot.Clone()
      sigplotforfitplusbkg.SetName(thissigplot.GetName()+"_forfitplusbkg"+thissuffix)
      sigplotforfitplusbkg.Scale(SignalScalingFactors[signal])

      sigplotforratio = thissigplot.Clone()
      sigplotforratio.SetName(thissigplot.GetName()+"_forratio"+thissuffix)

      if index==0 :
        sigplotforratio.Divide(fithist)
        signalPlots.append(sigplotforfitplusbkg)
        sigratioPlots.append(sigplotforratio)
        thistitle = SignalTitles[signal] + ", m = %d GeV" % mass
        legendlist.append(thistitle)

      else :
        sigplotforratio.Divide(newfithist)
        signalPlotsTeV.append(sigplotforfitplusbkg)
        sigratioPlotsTeV.append(sigplotforratio)
        thistitle = SignalTitles[signal] + ", m = {0} TeV".format(mass/1000.0)
        legendlistTeV.append(thistitle)

      index = index+1

  # Do regular plots

  outputName = folderextension+"SignalsOnSignificancePlot_"+signal+"_GeV"
  myPainter.drawSignalOverlaidOnBkgPlot(basicSignificancePlot,sigratioPlots,signalMasses,legendlist,luminosity,8,"[data - fit]/fit",outputName,firstBin,lastBin+2)

  outputName = folderextension+"SignalsOnFitPlusBkg_"+signal+"_GeV"
  myPainter.drawSignalOverlaidOnDataAndFit(datahist,fithist,signalPlots,signalMasses,legendlist,luminosity,8,"Events",outputName,firstBin,lastBin+2,True,True)

  outputName = folderextension+"SignalsOnFitPlusBkg_nolog_"+signal+"_GeV"
  myPainter.drawSignalOverlaidOnDataAndFit(datahist,fithist,signalPlots,signalMasses,legendlist,luminosity,8,"Events",outputName,firstBin,lastBin+2,False,True)

  outputName = folderextension+"FancyFigure1_"+signal+"_GeV"
  myPainter.drawDataAndFitWithSignalsOverSignificances(datahist,fithist,basicSignificancePlot,residual,\
       signalPlots,sigratioPlots,signalMasses,legendlist,"Reconstructed m_{jj} [GeV]","Prescale-weighted events","[data-fit]/fit",\
       "Signif.    ",outputName,luminosity,8,firstBin,lastBin+2,False,-1,-1,True,False)

  # Do TeV plots

  outputName = folderextension+"SignalsOnSignificancePlot_"+signal
  myPainter.drawSignalOverlaidOnBkgPlot(newbasicSignificancePlot,sigratioPlotsTeV,signalMassesTeV,legendlistTeV,\
            luminosity,8,"[data - fit]/fit",outputName,firstBin,lastBin+2)

  outputName = folderextension+"SignalsOnFitPlusBkg_"+signal
  myPainter.drawSignalOverlaidOnDataAndFit(newdatahist,newfithist,signalPlotsTeV,signalMassesTeV,legendlistTeV,\
            luminosity,8,"Events",outputName,firstBin,lastBin+2,True,True)

  outputName = folderextension+"SignalsOnFitPlusBkg_nolog_"+signal
  myPainter.drawSignalOverlaidOnDataAndFit(newdatahist,newfithist,signalPlotsTeV,signalMassesTeV,legendlistTeV,\
            luminosity,8,"Events",outputName,firstBin,lastBin+2,False,True)

  outputName = folderextension+"FancyFigure1_"+signal
  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"Reconstructed m_{jj} [TeV]","Prescale-weighted events","[data-fit]/fit","Signif.    ",\
     outputName,luminosity,8,firstBin,lastBin+2,False,-1,-1,True,False)

  outputName = folderextension+"FancyFigure1_nologx_"+signal
  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"Reconstructed m_{jj} [TeV]","Prescale-weighted events","[data-fit]/fit","Signif.    ",\
     outputName,luminosity,8,firstBin,lastBin+2,False,-1,-1,False,False)

