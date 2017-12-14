
#/usr/bin/env python

import os
import sys
import ROOT
from art.morisot import Morisot
from itertools import repeat

luminosity = 3.5

# Get input: basic gaussian
#ratios = [-1,0.07,0.10,0.15]
#names = ['#sigma_{G}/m_{G} = 0.15','#sigma_{G}/m_{G} = 0.10','#sigma_{G}/m_{G} = 0.07','#sigma_{G}/m_{G} = Res.']
ratios = [1E-10]
names = ['Truth: #sigma_{G}/m_{G} = Res.']

#names = ['#sigma_{G}/m_{G} = 0.10','#sigma_{G}/m_{G} = 0.07','#sigma_{G}/m_{G} = Res.']

subranges = [[1100,1300],[1300,1500],[1500,1700],[1700,1900],[1900,2200],[2200,2500],[2500,3000],[3000,4000],[4000,5000],[6000,7000]]
#[5000,6000],
folderextension = 'Gaussians/'

extensions = ['_LHCP_fold']

#plotextension = '_yStar0p6'
#plotextension = '_yStar0p3'

basicInputFileTemplates = [
"/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/StatisticalAnalysis/Bayesian/Gaussians/Dijet_LHCP_fold/GenericGaussians_0p08_low{0}_high{1}_{2}.root"
]
#basicInputFileTemplates = [
#"/cluster/warehouse/kpachal/DijetsSummer2016//Results/Gaussians/Dijet2016p2015/GenericGaussians_6p1_low{0}_high{1}_{2}.root"
#]

inputsForFailedPoints = '' #'/cluster/warehouse/kpachal/ResultFiles_DijetStat13/Data15_C_fGRL_20150718/GaussianLimits/GenericGaussians_mass{0}_{1}.root'
failedPoints = []

# Initialize painter
myPainter = Morisot()
# Internal
myPainter.setLabelType(2)
#myPainter.setLabelType(1)

minMassVal = {}

for basicInputFileTemplate in basicInputFileTemplates :

  plotextension = extensions[basicInputFileTemplates.index(basicInputFileTemplate)]

  graphs = []

  results = {}

  print "Beginning",plotextension,"minMassVal is",minMassVal

  # Retrieve search phase inputs
  for width in sorted(ratios,reverse=True) :

    if "0p6" in plotextension:
      minMassVal[width] = -1

    print "--------------------------------------"
    print "Beginning width",width,":"

    thisobserved = ROOT.TGraph()
    massAndWidths = {}

    if width>0 :
      widthfornames = int(1000*width)
      internalwidth = widthfornames
    else :
      widthfornames = 'resolutionwidth'
      internalwidth = -1000

    myrange = subranges

    outindex = 0
    for thisrange in sorted(myrange) :

      print "Using subrange2",thisrange
      outindex = outindex + 1

      # TEMP
      myrangelow = thisrange[0]
      myrangehigh = thisrange[1]
      filename = basicInputFileTemplate.format(thisrange[0],thisrange[1],widthfornames)
      print "Getting",filename

      if not os.path.isfile(filename) :
        continue
      file = ROOT.TFile.Open(filename)
      vectorName = "CLsPerMass_widthToMass{0}".format(internalwidth)
      cls = file.Get(vectorName)

      masses = file.Get("massesUsed")#massPoints")

      index = 0

      for i in range(len(masses)) :
        if cls[i]<0 :
          continue
        mass = masses[i]
        if mass > 6700 :
          continue
        print "Using mass",mass
        if "0p6" in plotextension and minMassVal[width] < 0 :
          minMassVal[width] = mass
        if "0p3" in plotextension and not mass < minMassVal[width] :
          print "Comparing",mass,"to",minMassVal[width]
          continue
        print mass,",",cls[i]
        if mass < myrangelow :
          continue

        # Replace with a fix-file if this mass point initially failed
        print int(mass),width,failedPoints
        if [int(mass),width] in failedPoints :
          otherfile = ROOT.TFile.Open(inputsForFailedPoints.format(int(mass),widthfornames))
          print "Getting new val out of",inputsForFailedPoints.format(int(mass),widthfornames)
          cl = otherfile.Get("CLOfRealLikelihood")[0]
          massAndWidths[mass] = cl/luminosity

        else :
          massAndWidths[mass] = cls[i]/luminosity

    index = 0
    for mass in sorted(massAndWidths.keys()) :
      thisobserved.SetPoint(index,mass,massAndWidths[mass])
      index = index+1

    print "setting results[",width,"] = ",massAndWidths
    results[width] = massAndWidths
    thisobserved.SetName("width_"+str(width))
    graphs.append(thisobserved)

  # Make everything shifted by 1000 for TeV plots
  shiftedgraphs = []
  d1, d2 = ROOT.Double(0), ROOT.Double(0)
  findMaxRange = 0
  for graph in graphs :
    newgraph = graph.Clone()
    newgraph.SetName(graph.GetName()+"_scaled")
    for np in range(newgraph.GetN()) :
      newgraph.GetPoint(np,d1,d2)
      newgraph.SetPoint(np,d1/1000.0,d2)
      if d1/1000 > findMaxRange :
        findMaxRange = d1/1000
    shiftedgraphs.append(newgraph)

  trueMaxRange = round(findMaxRange * 2) / 2 + 0.5

  myPainter.cutstring = "|y*| < 0.6"
  myPainter.drawSeveralObservedLimits(shiftedgraphs,names,folderextension+"GenericGaussians"+plotextension,"m_{G} [TeV]",\
     "#sigma #times #it{A} #times BR [pb]",luminosity,13,0.9,7.2,1E-2,50,[],ATLASLabelLocation="BottomR",cutLocation="Left")

  print "For table in note:"
#  mostMasses = results[sorted(results.keys())[0]]
  print results.keys()[0]
  groupMasses = []
  for key, value in results.iteritems() :#
    groupMasses = groupMasses + value.keys()
  groupMasses = list(set(groupMasses))
  print groupMasses
  mostMasses = sorted(groupMasses)
  for mass in sorted(mostMasses) :

     sys.stdout.write("{0}".format(float('%.3g' % (mass/1000.0))))
     sys.stdout.write(" & ")
     for width in sorted(results.keys()) :
       if mass in results[width].keys() :
         #sys.stdout.write('{0}'.format('%.2E' % results[width][mass]))
         sys.stdout.write('{0}'.format(float('%.3g' % (1000.0*results[width][mass]))))
         #sys.stdout.write("{0}".format(results[width][mass]))
       else :
         sys.stdout.write("-")
       if sorted(results.keys())[-1]!=width :
         sys.stdout.write(" & ")
       else :
         sys.stdout.write(" \\\\ ")

     sys.stdout.write("  \n")

  #froot = ROOT.TFile.Open(folderextension+"GenericGaussians"+plotextension+'.root','UPDATE')
  #for graph in shiftedgraphs :
      #graph.SetMarkerColor(ROOT.kBlack)
      #graph.SetMarkerStyle(20)
      #graph.Write(graph.GetName(),ROOT.TObject.kOverwrite)
