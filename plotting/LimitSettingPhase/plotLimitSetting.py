#!/usr/bin/env python

import ROOT
import os,sys
from art.morisot import Morisot
from array import array
from pprint import pprint
from decimal import Decimal

def rreplace(s, old, new, occurrence):
  li = s.rsplit(old, occurrence) 
  return new.join(li)

# Get input
xseceffInputFile = ROOT.TFile('./inputs/xsecandacceptance/CrossSectionsFor13Plotting.root')
templateInputFile = ROOT.TFile('./inputs/xsecandacceptance/MC15_20151017_1fb.root')
mcfile = ROOT.TFile("./inputs/MCForComparison/mjj_mc15_13TeV_361023_Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W_total_final.root")		
mchist = mcfile.Get("mjj_mc15_13TeV_361023_Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W_total_final")		
mchist.SetDirectory(0)		
mcfile.Close()

mcUPfile = ROOT.TFile("./inputs/MCForComparison/mjj_mc15_13TeV_361023_Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W_total_final_UP.root")		
mcUPhist = mcUPfile.Get("mjj_mc15_13TeV_361023_Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W_total_final_UP")		
mcUPhist.SetDirectory(0)		
mcUPfile.Close()

mcDOWNfile = ROOT.TFile("./inputs/MCForComparison/mjj_mc15_13TeV_361023_Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W_total_final_DOWN.root")		
mcDOWNhist = mcDOWNfile.Get("mjj_mc15_13TeV_361023_Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W_total_final_DOWN")	
mcDOWNhist.SetDirectory(0)		
mcDOWNfile.Close()

# Get results files
searchInputFile = ROOT.TFile('./results/Step1_SearchPhase/Step1_SearchPhase.root')
limitFileNameTemplate = "./results/Step3_LimitSettingPhase/Step3_LimitSettingPhase_{0}.root"
individualLimitFiles = "./results/Step2_setLimitsOneMassPoint/Step2_setLimitsOneMassPoint_{0}{1}_1p04fb.root"

# Options
doSyst = False
doQuadrants = True
folderextension = "./plots/"
plotextension = ""

# Define necessary quantities.
Ecm = 13
luminosity = 1000

# make plots folder i.e. make folder extension
if not os.path.exists(folderextension):
    os.makedirs(folderextension)

# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("ATLAS")
#myPainter.setEPS(True)
myPainter.setLabelType(0) # Sets label type i.e. Internal, Work in progress etc.
                          # See below for label explanation
# 0 Just ATLAS    
# 1 "Preliminary"
# 2 "Internal"
# 3 "Simulation Preliminary"
# 4 "Simulation Internal"
# 5 "Simulation"
# 6 "Work in Progress"

# Will run for each signal type you include here.
#signalInputFileTypes = ["QStar","BlackMax","ZPrime0p10","ZPrime0p20","ZPrime0p30","ZPrime0p60","ZPrime0p70","ZPrime0p80""ZPrime0p90","WPrime","QBHRS","QBHRS_FullSim"]
#signalInputFileTypes = ["QStar","BlackMax","ZPrime0p10","ZPrime0p20","ZPrime0p30","ZPrime0p40","ZPrime0p50","ZPrime0p60","WPrime","QBHRS","QBH"]
signalInputFileTypes = ["QStar","BlackMax","ZPrime0p10","ZPrime0p20","ZPrime0p30","ZPrime0p40","ZPrime0p50","WPrime","QBHRS","QBH"]

# Setup signal info

masses = {} # Define mass points for signals
masses["QStar"]         = ["2000","2500","3000","3500","4000","4500","5000","5500","6000","6500"]
masses["BlackMax"]      = ["4000","5000","5500","6000","6500","7000","7500","8000","8500","9000","9500"]
masses["QBHRS"]         = ["4500","5000","5500","6000","6500","7000","7500","8000","8500"]
masses["WPrime"]        = ["1200","1500","1700","2000","2500","3000","3500","4000","4500","5000","5500","6000","6500"]
masses["ZPrime0p10"]    = ["1500","2000","2500","3000","3500"]
masses["ZPrime0p20"]    = ["1500","2000","2500","3000","3500"]
masses["ZPrime0p30"]    = ["1500","2000","2500","3000","3500"]
masses["ZPrime0p40"]    = ["2000","2500","3000","3500"]
masses["ZPrime0p50"]    = ["2500","3000","3500"]
masses["ZPrime0p60"]    = ["2500","3000","3500"]
masses["QBH"]           = ["04500","05000","05500","06000","06500","07000","07500","08000","08500","09000"] # removed 09500, 10000. totally blank
masses["QBH_FullSim"]   = ["04500"] #,"10000"] last point currently totally blank
masses["QBHRS_FullSim"] = ["4500","8500"]

SignalTitles = {"QStar": "#it{q}*"  
                , "BlackMax": "QBH  (BM)"
                , "QBHRS" : "QBH  (RS)"
                , "WPrime" : "W'"
                , "ZPrime0p10": "Z' (0.10)"
                , "ZPrime0p20": "Z' (0.20)"
                , "ZPrime0p30": "Z' (0.30)"
                , "ZPrime0p40": "Z' (0.40)"
                , "ZPrime0p50": "Z' (0.50)"
                , "ZPrime0p60": "Z' (0.60)"
                , "QBH": "QBH  (QBH)"
                , "QBH_FullSim": "QBH  (QBH)"
                , "QBHRS_FullSim" : "QBH  (RS)"
                }

SignalAxes = {"QStar": {"X" : "M_{#it{q}*} [GeV]", "Y": "#sigma #times #it{A} #times BR [pb]"},
              "BlackMax": {"X" : "M_{th} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "QBHRS": {"X" : "M_{th} [GeV]", "Y":"#sigma #times #it{A} [pb]"},
              "WPrime": {"X" : "M_{W'} [GeV]", "Y":"#sigma #times #it{A} #times BR [pb]"},
              "ZPrime0p10": {"X" : "M_{Z'} [GeV]", "Y":"#sigma #times #it{A} #times BR [pb]"},
              "ZPrime0p20": {"X" : "M_{Z'} [GeV]", "Y":"#sigma #times #it{A} #times BR [pb]"},
              "ZPrime0p30": {"X" : "M_{Z'} [GeV]", "Y":"#sigma #times #it{A} #times BR [pb]"},
              "ZPrime0p40": {"X" : "M_{Z'} [GeV]", "Y":"#sigma #times #it{A} #times BR [pb]"},
              "ZPrime0p50": {"X" : "M_{Z'} [GeV]", "Y":"#sigma #times #it{A} #times BR [pb]"},
              "ZPrime0p60": {"X" : "M_{Z'} [GeV]", "Y":"#sigma #times #it{A} #times BR [pb]"},
              "QBH": {"X" : "M_{th} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "QBH_FullSim": {"X" : "M_{th} [GeV]", "Y": "#sigma #times #it{A} [pb]"},
              "QBHRS_FullSim": {"X" : "M_{th} [GeV]", "Y":"#sigma #times #it{A} [pb]"},
             }

yranges = {}
yranges['QStar'] = [1E-4,100] # was 1000, then was 5E3
yranges['BlackMax'] = [2E-7,2] # was 10
yranges['QBHRS'] = [3E-6,0.1] # was 10
yranges['WPrime'] = [5E-5,5] # was 10
yranges['ZPrime0p10'] = [2E-4,5]
yranges['ZPrime0p20'] = [2E-4,5]
yranges['ZPrime0p30'] = [2E-4,5]
yranges['ZPrime0p40'] = [2E-4,5]
yranges['ZPrime0p50'] = [2E-4,5]
yranges['ZPrime0p60'] = [2E-4,5]
yranges['QBH'] = [5E-5,1] # was 10
yranges['QBH_FullSim'] = [1E-3,10] # was 10
yranges['QBHRS_FullSim'] = [1E-4,1] # was 10

xranges = {}
AllPlotMaterials = {}

ZPrimeLimits = {} # For 2D plot		
ZPrimexsec = { '0.10' : {},		
	      '0.20' : {},		
	      '0.30' : {},		
	      '0.40' : {},		
	      '0.50' : {},		
	      '0.60' : {},		
	      '0.70' : {},		
	      '0.80' : {},		
	      '0.90' : {}		
	     }

# Loop over signals in signalInputFileTypes

for signal in signalInputFileTypes :

  limitInputFile = ROOT.TFile(limitFileNameTemplate.format(signal))

  axisnames = SignalAxes[signal]
  xname = axisnames["X"]
  yname = axisnames["Y"]

  # Define necessary quantities.
  signalName = signal+"_graph"

  # Limit-setting plots
  outputName = folderextension+"Limits_"+signal+plotextension
  observedGraph = limitInputFile.Get('observedXSecAccVersusMass')

  # X-check: Scaling observed graph by square root of the ratio of lumis to make equiv to 100 inv pb
  # to check if limits roughly agree with PUB note result for 1000 inv pb (1 inv fb)
  # Should be commented out, except for when doing x-checks!
  #for i in range (0,observedGraph.GetN()):
  # scale = pow(100.0/1000.0,0.5)
  # scale = float(scale) # needed by python to make sure get correct value
  # print "SCALING by "+str(scale)
  # observedGraph.GetY()[i] *= scale
  #for i in range (0,expectedGraph1Sigma.GetN()):
  # expectedGraph1Sigma.GetY()[i] *= scale
  #for i in range (0,expectedGraph2Sigma.GetN()):
  # expectedGraph2Sigma.GetY()[i] *= scale
  #### Scaled MC hereprint "SCALING by "+str(pow(1000/100,0.5))
  ####observedGraph.GetY()[i] *= pow(1000/100,0.5)

  expectedGraph1Sigma = limitInputFile.Get('expectedXSecAccVersusMass_oneSigma')
  expectedGraph2Sigma = limitInputFile.Get('expectedXSecAccVersusMass_twoSigma')

  if "ZPrime" in signal :
    code = signal.replace("0p","0")
    code = rreplace(code,"0","",1)
    signalGraph = xseceffInputFile.Get("{0}_graph".format(code))
  elif "QBHRS" in signal :
    signalGraph = xseceffInputFile.Get("RSQBH_graph")
  else :
    signalGraph = xseceffInputFile.Get("{0}_graph".format(signal))
    if "WPrime" in signal :
      for point in range(signalGraph.GetN()) :
        x = ROOT.Double(0) 
        y = ROOT.Double(0)
        if point < signalGraph.GetN()-1 :
          signalGraph.GetPoint(point+1,x,y)
          signalGraph.SetPoint(point,x,y)
        else :
          signalGraph.Set(signalGraph.GetN()-1)

  if signalGraph == None :
    signalGraph = ROOT.TGraph()
    signalGraph.SetPoint(0,-1,1)
    signalGraph.SetPoint(1,0,1)
  extrasignalGraph = 0
  extraextrasignalGraph = 0

  print "-"*50
  print "SIGNAL: "+signal
  print signalName
  print "-"*50
  print signal

  # Two theory lines
  if signal=="BlackMax" :
    extrasignalGraph = xseceffInputFile.Get("QBH_graph")
    extraextrasignalGraph = xseceffInputFile.Get("RSQBH_graph")
  elif signal=="QBH" :
    extrasignalGraph = xseceffInputFile.Get("BlackMax_graph")
    extraextrasignalGraph = xseceffInputFile.Get("RSQBH_graph")

  # Make everything shifted by 1000 for TeV plots
  shiftedgraphs = []
  d1, d2 = ROOT.Double(0), ROOT.Double(0)
  for graph in [observedGraph,expectedGraph1Sigma,expectedGraph2Sigma,signalGraph,extrasignalGraph,extraextrasignalGraph] :
    if graph==0 :
      continue
    print graph
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

  xranges[signal] = [float(xlow/1000),float(xhigh/1000)]
 
  if "ZPrime" in signal :
    localZPrimeDict = {}
    coupling = signal[-4:].replace("p",".")
    print signal, coupling
    print "mass\tobs\texp\ttheory"
    for masspoint in range(shiftedgraphs[0].GetN()) :
      b1, b2, b3, b4,b5 = ROOT.Double(0), ROOT.Double(0),ROOT.Double(0), ROOT.Double(0), ROOT.Double(0)
      shiftedgraphs[0].GetPoint(masspoint,b1,b2)
      thisdict = {}
      thisdict["obs"] = b2
      print b1,"\t",b2,"\t",
      shiftedgraphs[1].GetPoint(masspoint,b1,b3)
      thisdict["exp"] = b3
      print b3,"\t",
      shiftedgraphs[3].GetPoint(masspoint,b4,b5)
      print b5
      localZPrimeDict[b1] = thisdict
      ZPrimexsec[coupling][b4] = b5
    ZPrimeLimits[coupling] = localZPrimeDict

  # make TeV plots
  newxname = xname.replace("G","T")
  [obs,exp] = myPainter.drawLimitSettingPlot2Sigma(shiftedgraphs[0],shiftedgraphs[1],shiftedgraphs[2],shiftedgraphs[3],SignalTitles[signal],\
     outputName,newxname,yname,luminosity,Ecm,xlow/1000,xhigh/1000,thisyrange[0],thisyrange[1],False)

  # Two theory lines
  if signal=="BlackMax" :
    outputName = outputName+"_withQBH"
    #thisname = thisname+"_withBlackMax"
    #[[obs1,exp1],[obs2,exp2]] = myPainter.drawLimitSettingPlot2Sigma(observedGraph,expectedGraph1Sigma,expectedGraph2Sigma,\
    #   [signalGraph,extrasignalGraph],[SignalTitles[signal],SignalTitles["QBH0"]],thisname,xname,yname,luminosity,Ecm,\
    #   xlow,xhigh,thisyrange[0],thisyrange[1],False)
    [[obs1,exp1],[obs2,exp2],[obs3,exp3]] = myPainter.drawLimitSettingPlot2Sigma(shiftedgraphs[0],shiftedgraphs[1],shiftedgraphs[2],\
       [shiftedgraphs[3],shiftedgraphs[4],shiftedgraphs[5]],[SignalTitles[signal],SignalTitles["QBH"],SignalTitles["QBHRS"]],outputName,newxname,yname,luminosity,Ecm,\
       xlow/1000,xhigh/1000,thisyrange[0],thisyrange[1],False)

  if signal=="QBH" :
    outputName = outputName+"_withBlackMax"
    #[[obs1,exp1],[obs2,exp2]] = myPainter.drawLimitSettingPlot2Sigma(observedGraph,expectedGraph1Sigma,expectedGraph2Sigma,\
    #   [signalGraph,extrasignalGraph],[SignalTitles[signal],SignalTitles["QBH0"]],thisname,xname,yname,luminosity,Ecm,\
    #   xlow,xhigh,thisyrange[0],thisyrange[1],False)
    [[obs1,exp1],[obs2,exp2],[obs3,exp3]] = myPainter.drawLimitSettingPlot2Sigma(shiftedgraphs[0],shiftedgraphs[1],shiftedgraphs[2],\
       [shiftedgraphs[3],shiftedgraphs[4],shiftedgraphs[5]],[SignalTitles[signal],SignalTitles["BlackMax"],SignalTitles["QBHRS"]],outputName,newxname,yname,luminosity,Ecm,\
       xlow/1000,xhigh/1000,thisyrange[0],thisyrange[1],False)

  # Save stuff in dict for the four-plot plot
  thisPlotMaterials = {}
  thisPlotMaterials["observed"] = shiftedgraphs[0]
  thisPlotMaterials["expected_1sigma"] = shiftedgraphs[1]
  thisPlotMaterials["expected_2sigma"] = shiftedgraphs[2]
  thisPlotMaterials["signal"] = shiftedgraphs[3]
  if signal=="BlackMax" :
    thisPlotMaterials["extrasignal"] = shiftedgraphs[4]
    thisPlotMaterials["extraextrasignal"] = shiftedgraphs[5]
  thisPlotMaterials["signalLabel"] = SignalTitles[signal]
  if signal=="BlackMax" :
    thisPlotMaterials["extrasignalLabel"] = SignalTitles["QBH"]
    thisPlotMaterials["extraextrasignalLabel"] = SignalTitles["QBHRS"]
  thisPlotMaterials["xAxisName"] = newxname
  AllPlotMaterials[signal] = thisPlotMaterials

  # Print out limit values!
  print "-"*50
  print "File",limitFileNameTemplate
  print "Signal",signal
  print "Observed limit at 95% CL:",obs
  print "Expected limit at 95% CL:",exp,"\n"

  if signal=="BlackMax" or signal=="WStarSinX0" :
    print "Second signal curve: QBH"
    print "Observed limit from this plot:",obs2
    print "Expected limit from this plot:",exp2
    print "Third signal curve: QBHRS"
    print "Observed limit from this plot:",obs3
    print "Expected limit from this plot:",exp3
  print "-"*50

# Make 2D Z' limit plot
TwoDPlotMaterials = {}
aZPrime = False
for signal in signalInputFileTypes :
  if "ZPrime" in signal : aZPrime = True
if aZPrime :
  print "ZPrime cross sections times acceptances are"
  pprint(ZPrimexsec)
  print "And ZPrime observed limits are"
  pprint(ZPrimeLimits)

  couplings = [ '0.10', '0.20', '0.30' , '0.40', '0.50']#, '0.60','0.70','0.80','0.90' ]
  masses = [1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]

  h = ROOT.TH2F("mulimit","",20,0.25,10.25,10,0.,1.)

  for g in couplings:
    for m in masses:
        if not g in ZPrimexsec :
            continue
        if not m in ZPrimeLimits[g].keys():
            continue

        print "doing g, m",g,m

        print "Calculation is",ZPrimeLimits[g][m]["obs"],"/",(ZPrimexsec[g][m]),"=",ZPrimeLimits[g][m]["obs"]/(ZPrimexsec[g][m])


        thisval = ZPrimeLimits[g][m]["obs"]/(ZPrimexsec[g][m])
        writeval = float("%.2g" %thisval) #round(thisval,1)#'{0}'.format(float('%.3g' % thisval))
        #writeval = Decimal(thisval)
        h.Fill( m , g , writeval )

  TwoDPlotMaterials["hist"] = h
  TwoDPlotMaterials["xAxisName"] = "M_{Z'} [GeV]"
  TwoDPlotMaterials["yAxisName"] = "g_{q}"
  TwoDPlotMaterials["zAxisName"] = "#sigma_{limit}/#sigma_{theory}"
  TwoDPlotMaterials["signalLabel"] = "Z'"

  # Do the 2D Z' plot
  outputName = folderextension+"ZPrime_2D"
  myPainter.draw2DLimit(h,outputName,"M_{Z'} [TeV]",1.25,3.75,"g_{q}",0.0,0.50,"#sigma_{limit}/#sigma_{theory}",luminosity,Ecm)
 

# Do the new big quadranty plot for Tom
if doQuadrants:
  outputName = folderextension+"threePlusOnePlotQuadrant"
  myPainter.drawThreePlusOneLimitPlots_Grid(AllPlotMaterials["QStar"],AllPlotMaterials["WPrime"],AllPlotMaterials["BlackMax"],TwoDPlotMaterials,outputName,"Mass [TeV]","#sigma #times #it{A} [pb]",luminosity,Ecm,2.0,9.0,1.2,6.5,2E-3,10,1E-3,13)

  outputName = folderextension+"fourPlotQuadrant"
  myPainter.drawFourLimitPlots_Grid(AllPlotMaterials["BlackMax"],AllPlotMaterials["ZPrime0p30"],AllPlotMaterials["QStar"],AllPlotMaterials["WPrime"],outputName,"Mass [TeV]","#sigma #times #it{A} [pb]",luminosity,Ecm,2.0,9.0,1.2,6.5,5E-4,10,1E-3,13)

  outputName = folderextension+"fourPlotLine"
  myPainter.drawFourLimitPlots_Line(AllPlotMaterials["QStar"],AllPlotMaterials["BlackMax"],AllPlotMaterials["WPrime"],AllPlotMaterials["ZPrime0p30"],outputName,"#sigma #times #it{A} [pb]",luminosity,Ecm,xranges["QStar"][0],xranges["QStar"][1],xranges["BlackMax"][0],xranges["BlackMax"][1],xranges["WPrime"][0],xranges["WPrime"][1],xranges["ZPrime0p20"][0],xranges["ZPrime0p20"][1],5E-4,15)

# Draw signal overlay plots
SignalMasses = {"QStar": [4000, 5000], #[2000,5000] SWITCH
                "BlackMax": [5000,6500],#,6000],#[5000,5500] SWITCH
                "QBHRS" : [4500,5500],
                "WPrime" : [1500,2000],
                "ZPrime0p10" : [1500,2000],
                "ZPrime0p20" : [1500,2000],
                "ZPrime0p30" : [1500,2000],
                "ZPrime0p40" : [2000,2500],
                "ZPrime0p50" : [2500,3000],
                "ZPrime0p60" : [2500,3000],
                "QBH": [5500, 6500], 
                "QBHRS_FullSim" : [4500]
                }

SignalScalingFactors = {"QStar": 3, # 2 SWITCH
                        "BlackMax": 1,
                        "QBHRS": 3,
                        "WPrime": 100, 
                        "ZPrime0p10": 100,
                        "ZPrime0p20": 100,
                        "ZPrime0p30": 100,
                        "ZPrime0p40": 100,  
                        "ZPrime0p50": 100,
                        "ZPrime0p60": 100,
                        "QBH": 1, 
                        "QBH_FullSim": 1, 
                        "QBHRS_FullSim": 1
		       }

# Get basic hists
datahist = searchInputFile.Get("basicData")
datahist.SetDirectory(0)
fithist = searchInputFile.Get("basicBkgFrom4ParamFit")
fithist.SetDirectory(0)
basicSignificancePlot = searchInputFile.Get("relativeDiffHist")
basicSignificancePlot.SetDirectory(0)
residual = searchInputFile.Get("residualHist")
residual.SetDirectory(0)

# get bump info
statOfFitToData = searchInputFile.Get('bumpHunterPLowHigh')
bumpHunterStatFitToData = statOfFitToData[0]
bumpLowEdge = statOfFitToData[1]
bumpHighEdge = statOfFitToData[2]
bumpHunterStatOfFitToData = searchInputFile.Get('bumpHunterStatOfFitToData')
bumpHunterPVal = bumpHunterStatOfFitToData[1]

# and fit info
fitRange = searchInputFile.Get("FitRange")

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

mixedSignalCollection = []
mixedTitleCollection = []
mixedUserScaleTextCollection = ""
mixedSignalMassCollection = []
mixedExtLastBin = 0

for signal in signalInputFileTypes :
  print "in signal",signal

  if "FullSim" in signal :
    continue

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
  sigToGet = signal

  if signal == "QBH" or signal == "QBHRS":
    sigToGet = signal+"0"
  elif signal == "QBHRS_FullSim" :
    sigToGet = "QBHRS0"
  elif "ZPrime" in signal :
    sigToGet = "ZPrime"
  print "sigToGet is",sigToGet,"for signal",signal
  for mass in signalMasses :

    if "ZPrime" in signal :
      sigplotname = signal.replace("Prime","PrimemR{0}gSM".format(mass))
      sigplot = templateInputFile.Get("mjj_"+sigplotname+"_1fb_Nominal")
    else :
      sigplot = templateInputFile.Get("mjj_"+sigToGet+"{0}_1fb_Nominal".format(mass))
    sigplot.SetDirectory(0)
    # TEMPORARY
    newsigGeV = datahist.Clone()
    newsigGeV.SetName("sigplot_{0}_{1}_GeV".format(signal,mass))
    for bin in range(newsigGeV.GetNbinsX()+2) :
      if bin < sigplot.GetNbinsX()+1 :
        newsigGeV.SetBinContent(bin,sigplot.GetBinContent(bin))
        newsigGeV.SetBinError(bin,sigplot.GetBinError(bin))
      else :
        newsigGeV.SetBinContent(bin,0.0)
        newsigGeV.SetBinError(bin,0.0)

    sigplottev = newdatahist.Clone()
    sigplottev.SetName("sigplot_{0}_{1}_TeV".format(signal,mass))
    # TEMPORARY: REBIN TO WHAT WE NEED
    # Also find new last bin.
    for bin in range(sigplottev.GetNbinsX()+2) :
      sigplottev.SetBinContent(bin,newsigGeV.GetBinContent(bin))
      sigplottev.SetBinError(bin,newsigGeV.GetBinError(bin))

    index = 0
    extLastBin = lastBin
    for thissigplot,thissuffix in [[newsigGeV,""],[sigplottev,"_TeV"]] :

      # Normalise to correct amount
      # Lydia Updated from thissigplot.Scale(luminosity) as new input templates are 1 inv fb not 1 inv pb
      luminosity = float(luminosity) # cast as float otherwise if divide e.g. 7 by 1000 get 0
      thissigplot.Scale(luminosity/1000)
      sigplotforfitplusbkg = thissigplot.Clone()
      sigplotforfitplusbkg.SetDirectory(0)
      sigplotforfitplusbkg.SetName(thissigplot.GetName()+"_forfitplusbkg"+thissuffix)
      sigplotforfitplusbkg.Scale(SignalScalingFactors[signal])

      sigplotforratio = thissigplot.Clone()
      sigplotforratio.SetDirectory(0)
      sigplotforratio.SetName(thissigplot.GetName()+"_forratio"+thissuffix)

      if index==0 :
        sigplotforratio.Divide(fithist)
        signalPlots.append(sigplotforfitplusbkg)
        sigratioPlots.append(sigplotforratio)
        #thistitle = sigToGet + ", m = %d GeV" % mass
        thistitle = sigToGet + ", %s= %d GeV" % (SignalAxes[signal]["X"].split("[GeV]")[0].replace("M","m"),mass)
        legendlist.append(thistitle)

      else :
        sigplotforratio.Divide(newfithist)
        signalPlotsTeV.append(sigplotforfitplusbkg)
        sigratioPlotsTeV.append(sigplotforratio)
        #thistitle = SignalTitles[signal] + ", m = {0} TeV".format(mass/1000.0)
        thistitle = SignalTitles[signal] + ", {0}= {1} TeV".format(SignalAxes[signal]["X"].split("[GeV]")[0].replace("M","m"),mass/1000.0)
        legendlistTeV.append(thistitle)

        for bin in range(sigplotforfitplusbkg.GetNbinsX()+2) :
          if bin > extLastBin and sigplotforfitplusbkg.GetBinContent(bin) > 0.01 :
            extLastBin = bin
          if sigplotforfitplusbkg.GetBinLowEdge(bin) > 1.3*mass/1000.0 :
            continue
        if extLastBin < lastBin :
          extLastBin = lastBin 

        # For making mixed signals plot
        if (signal == "QStar" and signalMasses.index(mass) == 0) or (signal == "BlackMax" and signalMasses.index(mass) == 1) :
          mixedSignalCollection.append(sigplotforfitplusbkg)
          mixedTitleCollection.append(thistitle)
          if SignalScalingFactors[signal] == 1 :
            mixedUserScaleTextCollection += "{"+SignalTitles[signal]+"}"
          else :
            mixedUserScaleTextCollection+="{"+SignalTitles[signal]+",  #sigma #times "+str(SignalScalingFactors[signal])+"}"
          if mixedExtLastBin < extLastBin :
            mixedExtLastBin = extLastBin
          mixedSignalMassCollection.append(signalMassesTeV[signalMasses.index(mass)])

      index = index+1

  UserScaleText = SignalTitles[signal]
  if SignalScalingFactors[signal] == 1 :
    UserScaleText = SignalTitles[signal]
  else :
    UserScaleText = UserScaleText+",  #sigma #times "+str(SignalScalingFactors[signal])

  # Do TeV plots
  outputName = folderextension+"SignalsOnSignificancePlot_"+signal+plotextension
  myPainter.drawSignalOverlaidOnBkgPlot(newbasicSignificancePlot,sigratioPlotsTeV,signalMassesTeV,legendlistTeV,\
            luminosity,Ecm,"[data - fit]/fit",outputName,firstBin,lastBin+2)

  outputName = folderextension+"SignalsOnFitPlusBkg_"+signal+plotextension
  myPainter.drawSignalOverlaidOnDataAndFit(newdatahist,newfithist,signalPlotsTeV,signalMassesTeV,legendlistTeV,\
            luminosity,Ecm,"Events",outputName,firstBin,lastBin+2,True,True)

  doRightLeg = False

  if "QStar" in signal or "WPrime" in signal or "ZPrime" in signal:
    doRightLeg = True

  outputName = folderextension+"FancyFigure1_"+signal+plotextension
  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
     outputName,luminosity,Ecm,firstBin,extLastBin+5,True,bumpLowEdge/1000,bumpHighEdge/1000,True,False,doRightLeg,UserScaleText,True,bumpHunterPVal)#Lydia changed Prescale-weighted events to Events AND added label to say how much signal scaled by

  outputName = folderextension+"FancyFigure1_nobump_"+signal+plotextension
  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
     outputName,luminosity,Ecm,firstBin,extLastBin+5,False,bumpLowEdge/1000,bumpHighEdge/1000,True,False,doRightLeg,UserScaleText,True, bumpHunterPVal)#Lydia changed Prescale-weighted events to Events AND added label to say how much signal scaled by

  outputName = folderextension+"FancyFigure1WithFitLabels_"+signal+plotextension
  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
     outputName,luminosity,Ecm,firstBin,extLastBin+5,True,bumpLowEdge/1000,bumpHighEdge/1000,True,False,doRightLeg,UserScaleText,True,bumpHunterPVal,True,fitRange[0],fitRange[1])#Lydia changed Prescale-weighted events to Events AND added label to say how much signal scaled by

  outputName = folderextension+"FancyFigure1WithFitLabels_nobump_"+signal+plotextension
  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
     outputName,luminosity,Ecm,firstBin,extLastBin+5,False,bumpLowEdge/1000,bumpHighEdge/1000,True,False,doRightLeg,UserScaleText,True,bumpHunterPVal,True,fitRange[0],fitRange[1])#Lydia changed Prescale-weighted events to Events AND added label to say how much signal scaled by

  outputName = folderextension+"FancyFigure1WithMCComparison_"+signal+plotextension		
  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
     outputName,luminosity,Ecm,firstBin,extLastBin+15,True,bumpLowEdge/1000,bumpHighEdge/1000,True,False,doRightLeg,UserScaleText,True,bumpHunterPVal,True,fitRange[0],fitRange[1],mchist)

#  outputName = folderextension+"SignalsOnFitPlusBkg_nolog_"+signal+plotextension
#  myPainter.drawSignalOverlaidOnDataAndFit(newdatahist,newfithist,signalPlotsTeV,signalMassesTeV,legendlistTeV,\
#            luminosity,Ecm,"Events",outputName,firstBin,lastBin+2,False,True)

#  outputName = folderextension+"FancyFigure1_"+signal+plotextension
#  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
#     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"m_{jj} [TeV]","Events","[data-fit]/fit","Signif.    ",\
#     outputName,luminosity,Ecm,firstBin,lastBin+2,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False,UserScaleText,True,bumpHunterPVal)

#  outputName = folderextension+"FancyFigure1_nologx_"+signal+plotextension
#  myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,signalPlotsTeV,\
#     sigratioPlotsTeV,signalMassesTeV,legendlistTeV,"m_{jj} [TeV]","Events","[data-fit]/fit","Signif.    ",\
#     outputName,luminosity,Ecm,firstBin,lastBin+2,False,-1,-1,False,False,False,UserScaleText)#Lydia changed Prescale-weighted events to Events

# Systematics x-check plots

  if (doSyst) :
    # Cross-check plots
    systDict = {"BKG_normalisation": "Fit quality", "FUNCCHOICE": "Function choice", "LUMI":"Luminosity", "PDFAcc":"PDF acceptance","mjj_{0}{1}_1fb_JET_GroupedNP_1":"Grouped NP 1","mjj_{0}{1}_1fb_JET_GroupedNP_2":"Grouped NP 2","mjj_{0}{1}_1fb_JET_GroupedNP_3":"Grouped NP 3"}

    for mass in masses[signal]:

      # Cross-check plots: compare residual from search phase to that
      # from comparing the best fit in all nuisance parameters, background, and signal to data
      filename = individualLimitFiles.format(signal, mass)

      inputfile = ROOT.TFile(filename)
      print "in input file",filename
      fit2 = inputfile.Get("BestFullPrediction")
      fit2.SetDirectory(0)
      residual2 = inputfile.Get("residual_bestfitToData")
      residual2.SetDirectory(0)

      newfit2 = ROOT.TH1D("fit2_TeV","fit2_TeV",len(newbins)-1,array('d',newbins))
      newresidual2 = ROOT.TH1D("residual2_TeV","residual2_TeV",len(newbins)-1,array('d',newbins))
      for histnew,histold in [[newfit2,fit2], [newresidual2,residual2]] :
        for bin in range(histnew.GetNbinsX()+2) :
          histnew.SetBinContent(bin,histold.GetBinContent(bin))
          histnew.SetBinError(bin,histold.GetBinError(bin))

      outputName = folderextension+"NominalAndBestFits_{0}_m{1}TeV".format(signal,mass)+plotextension
      myPainter.drawMultipleFitsAndResiduals(newdatahist,[newfithist,newfit2],[newresidual,newresidual2],["Nominal fit","Best fit in all #theta"],\
	"m_{jj} [TeV]","Events",["Nominal","Best Fit  "],outputName,luminosity,Ecm,\
	firstBin,lastBin+2)

#      outputName = folderextension+"NominalAndBestFits_{0}_m{1}TeV_megazoom".format(signal,mass)+plotextension
#      myPainter.drawMultipleFitsAndResiduals(newdatahist,[newfithist,newfit2],[newresidual,newresidual2],["Nominal fit","Best fit in all #theta"],\
#	"m_{jj} [TeV]","Events",["Nominal","Best Fit  "],outputName,luminosity,Ecm,\
#	newdatahist.FindBin(1.55),newdatahist.FindBin(1.74),False,0,0,True,False,notLogY=True)

      # Posteriors with CLs
      signalpost = inputfile.Get("likelihoodFunction")
      signalpost.SetDirectory(0)
      signalCL = inputfile.Get("CLOfRealLikelihood")[0]

      outputName = folderextension+"SignalPosteriorWith95CL_{0}_m{1}TeV".format(signal,mass)
      inlist = [signalpost,signalCL]
      myPainter.drawPosteriorsWithCLs([inlist],["Signal posterior"],luminosity,Ecm,outputName,2,True,False,True)

      for syst in systDict.keys() :
        gaus = ROOT.TF1("template_f", "gaus", 3.5,-3.5)
        gaus.SetParameters(1,0,1.0)
        if 'GroupedNP' in syst:
          print "Getting",syst.format(sigToGet,int(mass))
          gaushist = inputfile.Get(syst.format(sigToGet,int(mass))).Clone()
        else:
          gaushist = inputfile.Get(syst).Clone()
        gaushist.Reset("ICE")
        for bin in range(1,gaushist.GetNbinsX()+1) :
          cent = gaushist.GetBinCenter(bin)
          width = gaushist.GetBinWidth(bin)
          gaushist.SetBinContent(bin,gaus.Eval(cent)*width)
          gaushist.SetBinError(bin,0)

        gaushist.Scale(1.0/gaushist.Integral())
        gaushist.SetDirectory(0)

        if 'GroupedNP' in syst:
          systhist = inputfile.Get(syst.format(sigToGet,int(mass))).Clone(syst+"_mine")
          systhist.SetDirectory(0)
        else:
          systhist = inputfile.Get(syst.format(mass)).Clone(syst+"_mine")
          systhist.SetDirectory(0)

        cl = -100
        systhist.Scale(1.0/systhist.Integral())
        pair = [systhist,[cl]]
        pairs = [pair]
        if 'GroupedNP' in syst:
          outputname = folderextension+"posterior_"+syst.format(signal,mass)+"_{0}_{1}".format(signal,mass)+plotextension
        else:
          outputname = folderextension+"posterior_"+syst+"_{0}_{1}".format(signal,mass)+plotextension
        shortnames = [systDict[syst]]
        if 'FUNC' in syst :
          myPainter.drawPosteriorsWithCLs(pairs,shortnames,luminosity,Ecm,outputname,0,True,False,False,False,[gaushist],False,"Nuisance parameter [0 ,1]")
        else :
          myPainter.drawPosteriorsWithCLs(pairs,shortnames,luminosity,Ecm,outputname,0,True,False,False,False,[gaushist],False,"Nuisance parameter, #sigma")

      inputfile.Close()

outputName = folderextension+"FancyFigure1_BothSignals"+plotextension
myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,mixedSignalCollection,\
  [],mixedSignalMassCollection,mixedTitleCollection,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
  outputName,luminosity,Ecm,firstBin,mixedExtLastBin,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False)

outputName = folderextension+"FancyFigure1WithFitLabels_BothSignals"+plotextension
myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,mixedSignalCollection,\
  [],mixedSignalMassCollection,mixedTitleCollection,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
  outputName,luminosity,Ecm,firstBin,mixedExtLastBin,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False,mixedUserScaleTextCollection,True,bumpHunterPVal,True,fitRange[0],fitRange[1])

outputName = folderextension+"FancyFigure1WithFitLabels_nobump_BothSignals"+plotextension
myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,mixedSignalCollection,\
  [],mixedSignalMassCollection,mixedTitleCollection,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
  outputName,luminosity,Ecm,firstBin,mixedExtLastBin,False,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False,mixedUserScaleTextCollection,True,bumpHunterPVal,True,fitRange[0],fitRange[1])

outputName = folderextension+"FancyFigure1WithMCComparison_BothSignals"+plotextension		
myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,mixedSignalCollection,\
  [],mixedSignalMassCollection,mixedTitleCollection,"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
  outputName,luminosity,Ecm,firstBin,mixedExtLastBin,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False,mixedUserScaleTextCollection,True,bumpHunterPVal,True,fitRange[0],fitRange[1],mchist)

outputName = folderextension+"FancyFigure1_NoSignals"+plotextension
myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,[],\
  [],[],[],"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
  outputName,luminosity,Ecm,firstBin,lastBin+2,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False)

outputName = folderextension+"FancyFigure1WithFitLabels_NoSignals"+plotextension
myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,[],\
  [],[],[],"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
  outputName,luminosity,Ecm,firstBin,lastBin+2,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False,"",True,bumpHunterPVal,True,fitRange[0],fitRange[1])

outputName = folderextension+"FancyFigure1WithFitLabels_nobump_NoSignals"+plotextension
myPainter.drawDataAndFitWithSignalsOverSignificances(newdatahist,newfithist,newbasicSignificancePlot,newresidual,[],\
  [],[],[],"m_{jj} [TeV]","Events","[data-fit]/fit","Significance ",\
  outputName,luminosity,Ecm,firstBin,lastBin+2,False,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False,"",True,bumpHunterPVal,True,fitRange[0],fitRange[1])

# ---------------------------- 

# FIXME temporary solution to chop plot at 8 TeV, where JES bands end 
# SWITCH 
mixedExtLastBin = 130

# Scaling not necessary as we have updated inputs now!
#mchist.Scale(3.57/3.27)
#mcUPhist.Scale(3.57/3.27)
#mcDOWNhist.Scale(3.57/3.27)

tmpRatioHist = newdatahist.Clone()
tmpRatioHist.SetMarkerColor(ROOT.kBlack)
tmpRatioHist.Add(mchist,-1)
tmpRatioHist.Divide(mchist)

UpDownHists = []

print "UP"
print mcUPhist.GetEntries()
if mcUPhist.GetEntries() >= 0:
  print "UP!"
  tmpJESRatioHistup = mcUPhist
  tmpJESRatioHistup.Add( mchist, -1. )
  tmpJESRatioHistup.Divide( mchist )
  tmpJESRatioHistup.SetMarkerColorAlpha( ROOT.kBlue,0.15)
  tmpJESRatioHistup.SetLineColorAlpha( ROOT.kBlue,0.15)
  tmpJESRatioHistup.SetFillColorAlpha( ROOT.kBlue, 0.15)
  tmpJESRatioHistup.SetFillStyle(1001)
  UpDownHists.append(tmpJESRatioHistup) 

print "DOWN"
print mcDOWNhist.GetEntries()
if mcDOWNhist.GetEntries() >= 0:
  print "DOWN!"
  tmpJESRatioHistdown = mcDOWNhist.Clone()
  tmpJESRatioHistdown.Add( mchist, -1. )
  tmpJESRatioHistdown.Divide( mchist )
  tmpJESRatioHistdown.SetMarkerColorAlpha( ROOT.kBlue,0.15)
  tmpJESRatioHistdown.SetLineColorAlpha( ROOT.kBlue,0.15)
  tmpJESRatioHistdown.SetFillColorAlpha( ROOT.kBlue, 0.15)
  tmpJESRatioHistdown.SetFillStyle(1001)
  UpDownHists.append(tmpJESRatioHistdown) 

## If data is 0 then there should be no ratio drawn
for iBin in range(1, tmpRatioHist.GetNbinsX()+1):
  if newdatahist.GetBinContent(iBin) == 0:
    tmpRatioHist.SetBinContent(iBin, 0)
    tmpRatioHist.SetBinError(iBin, 0)

outputName = folderextension+"FancyFigure1WithFitLabels_WithMCRatio_BothSignals"+plotextension		
myPainter.test(newdatahist,newfithist,newbasicSignificancePlot,newresidual,mixedSignalCollection,\
  [],mixedSignalMassCollection,mixedTitleCollection,"m_{jj} [TeV]","Events","#frac{Data-MC}{MC}","#splitline{Significance}{data - fit}",\
  outputName,luminosity,Ecm,firstBin,mixedExtLastBin,True,bumpLowEdge/1000.0,bumpHighEdge/1000.0,True,False,False,mixedUserScaleTextCollection,True,bumpHunterPVal,True,fitRange[0],fitRange[1],mchist,tmpRatioHist,UpDownHists[0],UpDownHists[1])

# -----------------------------
