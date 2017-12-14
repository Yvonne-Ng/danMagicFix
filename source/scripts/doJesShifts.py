#!/usr/bin/env python
import os
import ROOT
from ROOT import *
import numpy as np

#from art.options import Options # Include it to apply ATLAS style settings
#from art.morisot import Morisot

# define input path, extension
inDir = '../inputs/MC15_20151017/QStar/dataLikeHists_v1/StatisticalHists/1fb'
nameofpdf="QStarJESShifts1invfbJES1.pdf"
nameofpdf1="QStarJESShifts1invfbJES1Component1.pdf"
nameofpdf2="QStarJESShifts1invfbJES1Component2.pdf"
nameofpdf3="QStarJESShifts1invfbJES1Component3.pdf"
nameofroot3="QStarJESShifts1invfbJES1Component3"

markerstyle=20
markersize=0.8

gStyle.SetOptStat(0)

Signals = {}
MassList=["2000","2500","3000","3500","4000","4500","5000","5500","6000","6500"] 
nBins = len(MassList)
InputFileList = []

canvas=ROOT.TCanvas("Beginning","Beginning",0,0,800,600)
canvas.Print(nameofpdf+"(")

for Mass in MassList:
  myfile = inDir+'/QStar%s_1fb.root'%(Mass)
  InputFileList.append(ROOT.TFile(myfile))

CompDown1Diff = ROOT.TH1D("CompDown1Diff","",nBins,0,nBins)
CompUp1Diff = ROOT.TH1D("CompUp1Diff","",nBins,0,nBins)
CompDown2Diff = ROOT.TH1D("CompDown2Diff","",nBins,0,nBins)
CompUp2Diff = ROOT.TH1D("CompUp2Diff","",nBins,0,nBins)
CompDown3Diff = ROOT.TH1D("CompDown3Diff","",nBins,0,nBins)
CompUp3Diff = ROOT.TH1D("CompUp3Diff","",nBins,0,nBins)

canvas1=ROOT.TCanvas("Canvas1","Canvas1",0,0,800,600)
canvas2=ROOT.TCanvas("Canvas2","Canvas2",0,0,800,600)
canvas3=ROOT.TCanvas("Canvas3","Canvas3",0,0,800,600)

lXmin = 0.60
lXmax = 0.80
lYmin = 0.25
lYmax = 0.45
leg1 = ROOT.TLegend(lXmin, lYmin, lXmax, lYmax)
leg2 = ROOT.TLegend(lXmin, lYmin, lXmax, lYmax)
leg3 = ROOT.TLegend(lXmin, lYmin, lXmax, lYmax)

for iBin,(Mass,InputFile) in enumerate(zip(MassList,InputFileList)):
  #print iBin
  #print Mass
  nominalHisto =InputFile.Get("mjj_QStar%s_1fb_Nominal"%Mass)
  # 1 2 3 are components
  Comp1DownHist =InputFile.Get("mjj_QStar%s_1fb_JET_GroupedNP_1__3down"%Mass)
  Comp1UpHist =InputFile.Get("mjj_QStar%s_1fb_JET_GroupedNP_1__3up"%Mass)
  Comp2DownHist =InputFile.Get("mjj_QStar%s_1fb_JET_GroupedNP_2__3down"%Mass)
  Comp2UpHist =InputFile.Get("mjj_QStar%s_1fb_JET_GroupedNP_2__3up"%Mass)
  Comp3DownHist =InputFile.Get("mjj_QStar%s_1fb_JET_GroupedNP_3__3down"%Mass)
  Comp3UpHist =InputFile.Get("mjj_QStar%s_1fb_JET_GroupedNP_3__3up"%Mass)

  Comp1maxDownValue=abs(Comp1DownHist.GetMean()-nominalHisto.GetMean())/nominalHisto.GetMean()
  Comp1maxUpValue=abs(Comp1UpHist.GetMean()-nominalHisto.GetMean())/nominalHisto.GetMean()
  Comp2maxDownValue=abs(Comp2DownHist.GetMean()-nominalHisto.GetMean())/nominalHisto.GetMean()
  Comp2maxUpValue=abs(Comp2UpHist.GetMean()-nominalHisto.GetMean())/nominalHisto.GetMean()
  Comp3maxDownValue=abs(Comp3DownHist.GetMean()-nominalHisto.GetMean())/nominalHisto.GetMean()
  Comp3maxUpValue=abs(Comp3UpHist.GetMean()-nominalHisto.GetMean())/nominalHisto.GetMean()
  
  CompDown1Diff.SetBinContent(iBin+1,Comp1maxDownValue)
  CompDown1Diff.SetBinError(iBin+1,0.0)
  CompDown1Diff.GetXaxis().SetBinLabel(iBin+1, Mass) # +1 as Bins start at 1

  CompUp1Diff.SetLineColor(kRed)
  CompUp1Diff.SetBinContent(iBin+1,Comp1maxUpValue)
  CompUp1Diff.SetBinError(iBin+1,0.0)
  CompUp1Diff.GetXaxis().SetBinLabel(iBin+1, Mass) # +1 as Bins start at 1

  CompDown2Diff.SetLineColor(kViolet)
  CompDown2Diff.SetBinContent(iBin+1,Comp2maxDownValue)
  CompDown2Diff.SetBinError(iBin+1,0.0)
  CompDown2Diff.GetXaxis().SetBinLabel(iBin+1, Mass) # +1 as Bins start at 1

  CompUp2Diff.SetLineColor(kOrange+1)
  CompUp2Diff.SetBinContent(iBin+1,Comp2maxUpValue)
  CompUp2Diff.SetBinError(iBin+1,0.0)
  CompUp2Diff.GetXaxis().SetBinLabel(iBin+1, Mass) # +1 as Bins start at 1

  CompDown3Diff.SetLineColor(kGreen+1)
  CompDown3Diff.SetBinContent(iBin+1,Comp3maxDownValue)
  CompDown3Diff.SetBinError(iBin+1,0.0)
  CompDown3Diff.GetXaxis().SetBinLabel(iBin+1, Mass) # +1 as Bins start at 1

  CompUp3Diff.SetLineColor(kYellow+3)
  CompUp3Diff.SetBinContent(iBin+1,Comp3maxUpValue)
  CompUp3Diff.SetBinError(iBin+1,0.0)
  CompUp3Diff.GetXaxis().SetBinLabel(iBin+1, Mass) # +1 as Bins start at 1

  if iBin == 1:   
    leg1.AddEntry(CompDown1Diff, "Component 1, down ", "l")
    leg1.AddEntry(CompUp1Diff, "Component 1, up", "l")
    leg2.AddEntry(CompDown2Diff, "Component 2, down ", "l")
    leg2.AddEntry(CompUp2Diff, "Component 2, up ", "l")
    leg3.AddEntry(CompDown3Diff, "Component 3, down ", "l")
    leg3.AddEntry(CompUp3Diff, "Component 3, up ", "l")
  #n = len(y)
  #print relativeUncHisto.GetMean()
  #print nominalHisto.GetMean()
  #print "Relative difference: " , abs(relativeUncHisto.GetMean()-nominalHisto.GetMean())/nominalHisto.GetMean()
  ##print "Error on difference, assuming uncorrelated " , sqrt(pow(relativeUncHisto.GetMeanError(),2)+pow(nominalHisto.GetMeanError(),2))

  #canvas.SetLogy(True)
XName = "Mass Point [GeV]"
YName = "Relative uncertainty"
Offset = 1.5
canvas1.cd()
CompDown1Diff.GetXaxis().SetTitle(XName)
CompDown1Diff.GetYaxis().SetTitle(YName)
CompDown1Diff.GetYaxis().SetTitleOffset(Offset)
CompDown1Diff.SetMinimum(0.015)
CompDown1Diff.Draw("HIST")

CompUp1Diff.GetXaxis().SetTitle(XName)
CompUp1Diff.GetYaxis().SetTitle(YName)
CompUp1Diff.GetYaxis().SetTitleOffset(Offset)
CompUp1Diff.SetMinimum(0.015)
CompUp1Diff.Draw("HIST same")
leg1.Draw() 
canvas1.Print(nameofpdf)
canvas1.Print(nameofpdf1)

canvas2.cd()
CompDown2Diff.GetXaxis().SetTitle(XName)
CompDown2Diff.GetYaxis().SetTitle(YName)
CompDown2Diff.GetYaxis().SetTitleOffset(Offset)
CompDown2Diff.SetMinimum(0.023)
CompDown2Diff.Draw("HIST")

CompUp2Diff.GetXaxis().SetTitle(XName)
CompUp2Diff.GetYaxis().SetTitle(YName)
CompUp2Diff.GetYaxis().SetTitleOffset(Offset)
CompUp2Diff.SetMinimum(0.023)
CompUp2Diff.Draw("HIST same")
leg2.Draw() 
canvas2.Print(nameofpdf)
canvas2.Print(nameofpdf2)

canvas3.cd()
CompDown3Diff.GetXaxis().SetTitle(XName)
CompDown3Diff.GetYaxis().SetTitle(YName)
CompDown3Diff.GetYaxis().SetTitleOffset(Offset)
CompDown3Diff.Draw("HIST")
MyDownFile = ROOT.TFile(nameofroot3+"Down.root","RECREATE")
CompDown3Diff.Write("CompDown3")#, ROOT.TObject.kOverwrite)
MyDownFile.Close()

CompUp3Diff.GetXaxis().SetTitle(XName)
CompUp3Diff.GetYaxis().SetTitle(YName)
CompUp3Diff.GetYaxis().SetTitleOffset(Offset)
CompUp3Diff.Draw("HIST same")
leg3.Draw() 
canvas3.Print(nameofpdf)
canvas3.Print(nameofpdf3)
MyUpFile = ROOT.TFile(nameofroot3+"Up.root","RECREATE")
CompUp3Diff.Write("CompUp3")#, ROOT.TObject.kOverwrite)
MyUpFile.Close()

canvas=ROOT.TCanvas("Ending","Ending",0,0,800,600)
canvas.Print(nameofpdf+")")


