#!/bin/python

#******************************************
#perform sensitivity scan by looping over all masses on a range of luminosity values
#EXAMPLE python -u doSensitivityScan.py --config <config file> --tag <tag> --batch --debug

#******************************************
#import stuff
import sys, os, math, argparse, ROOT
import sensitivityTools
import numpy as np
import os.path
from fileNamingTool import *


def removeOldLabelledFile(spFileDir, label, mass, window, gauss=7):
    """test to see if old labelled file of the mass point /window width and gaussian width is still arround, if so delete them"""
    OldTaggedFileMayExist=True
    i=0
    while(OldTaggedFileMayExist and i<20):
        try:
            testFile=findLabelledFileName(spFileDir, label, mass , window, gauss=15)
            #print("TaggedFileToBe ", nosignalFN)
            #os.rename(testFile,noSignalFN)
            #if the above returns something instead of raising an error,
            i=i+1
        except ValueError:
            OldTaggedFileMayExist=False
            print("no other duplicatedly labelled file from previous runs ")

        else :
            print ("removing this file: ", testFile)
            os.remove(testFile)

#******************************************
def doSensitivityScan(args):
    spFileDir="/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts/../results2/searchphase/"
    print("# of justabove files afte before eveyrhing : ", countFilesInDirWithKeyword(spFileDir, "JUSTABOVE"))

    print '\n******************************************'
    print '\n******************************************'
    print '\n******************************************'
    print '\n******************************************'
    print 'sensitivity scan starting'

    #------------------------------------------
    #input parameters
    print '\ninput parameters:'
    argsdict = vars(args)
    for ii in xrange(len(argsdict)):
        print '  %s = %s'%(argsdict.keys()[ii], argsdict.values()[ii],)

    #------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #get settings
    print '\nconfig settings:'
    settings = ROOT.TEnv()
    if settings.ReadFile(args.configFileName,ROOT.EEnvLevel(0)) != 0:
        raise SystemExit('***ERROR*** could not find sensitivity scan config file: %s'%args.configFileName)

    model = settings.GetValue('signalModel','')
    print '  signal model = %s'%model

    modelLabel = settings.GetValue('signalModelLabel','').replace('"','')
    print '  signal model label = %s'%modelLabel

    massValuesConfig = settings.GetValue('signalMasses','2000,3000,4000').split(',')
    massValuesConfig = [float(m) for m in massValuesConfig]
    print '  signal masses [GeV] = %s'%massValuesConfig

    lumiMin = float(settings.GetValue('luminosityMin','0.1'))
    #if lumiMin < 0.1: #fb^-1
    #    lumiMin = 0.1 #fb^-1
    print '  minimum luminosity = %s'%lumiMin

    lumiMax = float(settings.GetValue('luminosityMax','10.'))
    print '  maximum luminosity = %s'%lumiMax
    if lumiMax > 10000.: #fb^-1
        lumiMax = 10000. #fb^-1

    #QCDFileName = settings.GetValue('QCDFile','../inputs/QCD/histograms.mc.dijet.1p0.ifb.root')
    #Yvonne edit for dijetISR
    QCDFileName = settings.GetValue('QCDFile','../inputs/QCD/histograms.mc.QCD.1p0.ifb.root')
    print '  QCD input file = %s'%QCDFileName

    histBaseName = settings.GetValue('histBaseName','mjj')
    print '  hist base name = %s'%histBaseName

    bTaggingWP = settings.GetValue('bTaggingWP','') #fix_8585
    print '  b-tagging WP = %s'%bTaggingWP

    axisLabel = settings.GetValue('axisLabel','m [GeV]')
    print '  hist x-axis label = %s'%axisLabel

    nPar = int(settings.GetValue('nFitParameters','3'))
    print '  n fit parameters = %s'%nPar

    nPseudoExps_withSig = settings.GetValue('nPseudoExperimentsWithSig','1')
    print '  number of pseudo-experiments = %s'%nPseudoExps_withSig

    nPseudoExps_bkg = settings.GetValue('nPseudoExperiments_bkg','1000')
    print '  number of pseudo-experiments = %s'%nPseudoExps_bkg

    thresholdMass = float(settings.GetValue('thresholdMass','1100.'))
    print '  threshold mass = %s'%thresholdMass

    seed = float(settings.GetValue('randomSeed','0'))
    print '  random seed = %s'%seed

    notes = settings.GetValue('notes','').split(',')
    print '  notes = %s'%notes

    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    #------------------------------------------
    #set variables
    histName = histBaseName
    if bTaggingWP != '':
        histName+='_'+bTaggingWP
    if args.debug:
        print '\nhist name = %s'%histName

    #------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #check input file
    if not os.path.isfile(QCDFileName):
        raise SystemExit('\n***ERROR*** couldn\'t find QCD file: %s'%QCDFileName)

    #------------------------------------------

    #for windowWidth in [13, 12, 10, 9, 8]:#add 9 later
    #for windowWidth in [12, 10, 8, 7, 6]:#add 9 later
    for windowWidth in [12]:#add 9 later


      print '\n******************************************'
      print "***************** Window Width "+str(windowWidth)+" ********************"

      #get list of available mass points
      massValuesAvailable = []

      fileList = os.listdir(localdir+'/../inputs2/'+model+'/')
      print "mass directory: ", localdir+'/../inputs2/'+model+'/'
      for sigFileName in sorted(fileList):
          if not '.root' in sigFileName:
              continue
          massValuesAvailable.append( float( sensitivityTools.getSignalMass(sigFileName)))
      massValuesAvailable.sort(key=float)
      #massValues = [m for m in massValues if m >= 10000.0] #DEBUG
      #massValues = [3000.0] #DEBUG
      #print '\navailable mass values [GeV]: %s'%massValuesAvailable

      massValues = list( set(massValuesConfig) & set(massValuesAvailable) )
      massValues.sort(key=float)
      massValues.reverse()
      print 'using mass values [GeV]: %s'%massValues

      #------------------------------------------
      #TEST
      #raise SystemExit('\n***TEST*** exit')
      #------------------------------------------

      #------------------------------------------
      #initial luminosity value
      #lumi = lumiMin

      #------------------------------------------
      #arrays for sensitivity scan graphs
      gmass = np.array(massValues)
      glumi = np.zeros_like(gmass)

      #lumiSteps =[0.0, 0.1,0.2,0.3,0.5,0.7]+range(1,10)+range(10,20,1)+range(20,30,2)+range(30,50,3)+range(50,100,5)+range(100,200,10)+range(200,1000,50))
      gausswidthStartStep={"15": 29}
      lumiSteps=[0.0, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 25,  30, 31, 32, 33, 34, 35, 36, 40, 42, 44, 45, 47, 48, 49, 50, 52, 54, 55, 56, 57]+range(100, 1000, 20)
      #lumiSteps = [0.1]
      print "lumiSteps", lumiSteps

      #reset lumi to lowest value
      lumiStep = 1
      lumiIncrement = 1.
      #------------------------------------------
      #loop over mass values
      for mass in massValues : #and lumi <= lumiMax:

        print '\n******************************************'
        print '\n******************************************'

        setInitialLumi = False
        #doing 1000/fb doesn't make much sense, may as well give up?
        #Yvonne:^^ HUHHH??
        #while lumi<(lumiSteps[-1]) :
        for lumi in lumiSteps:

          ##add algorithm here...

          #check it makes sense to run this mass/lumi combination
          #####TLA2016 thresholds: don't test if above a given lumi as too high

          #ww9
          #mass = [ 1850.  1750.  1050.   750.   650.]
          #lumi = [  6.2942564   6.2942564  14.95725    45.          0.       ]
          #ww11
          #mass = [ 1850.  1750.  1050.   750.   650.]
          #lumi = [ 5.46587075  5.46587075  0.          0.          0.        ]

          #just hardwire because i'm so tired
          #Yvonne: ^^LOL

        #  if mass == 650 and windowWidth == 9 and setInitialLumi==False:
        #    lumi = 74 #136 for 7*
        #    setInitialLumi = True
        #  elif mass == 650 and windowWidth == 10 and setInitialLumi==False:
        #    lumi = 72
        #    setInitialLumi = True
        #  elif mass == 650 and windowWidth == 11 and setInitialLumi==False:
        #    lumi = 70
        #    setInitialLumi = True
        #  elif mass == 650 and windowWidth == 12 and setInitialLumi==False:
        #    lumi = 55
        #    setInitialLumi = True
        #  if mass == 750 and windowWidth == 9 and setInitialLumi==False:
        #    lumi = 35
        #    setInitialLumi = True
        #  elif mass == 750 and windowWidth == 10 and setInitialLumi==False:
        #    lumi = 30
        #    setInitialLumi = True
        #  elif mass == 750 and windowWidth == 11 and setInitialLumi==False:
        #    lumi = 26
        #    setInitialLumi = True
        #  elif mass == 750 and windowWidth == 12 and setInitialLumi==False:
        #    lumi = 15
        #    setInitialLumi = True
        #  if mass == 1050 and windowWidth == 9 and setInitialLumi==False:
        #    lumi = 13
        #    setInitialLumi = True
        #  elif mass == 1050 and windowWidth == 10 and setInitialLumi==False:
        #    lumi = 12
        #    setInitialLumi = True
        #  elif mass == 1050 and windowWidth == 11 and setInitialLumi==False:
        #    lumi = 11
        #    setInitialLumi = True
        #  elif mass == 1050 and windowWidth == 12 and setInitialLumi==False:
        #    lumi = 10
        #    setInitialLumi = True
        #  if mass == 1450 and windowWidth == 9 and setInitialLumi==False:
        #    lumi = 5
        #    setInitialLumi = True
        #  elif mass == 1450 and windowWidth == 10 and setInitialLumi==False:
        #    lumi = 4
        #    setInitialLumi = True
        #  elif mass == 1450 and windowWidth == 11 and setInitialLumi==False:
        #    lumi = 3
        #    setInitialLumi = True
        #  elif mass == 1450 and windowWidth == 12 and setInitialLumi==False:
        #    lumi = 2
        #    setInitialLumi = True
        #  if mass == 1850 and windowWidth == 9 and setInitialLumi==False:
        #    lumi = 5
        #    lumiIncrement = 1
        #    setInitialLumi = True
        #  elif mass == 1850 and windowWidth == 10 and setInitialLumi==False:
        #    lumi = 4
        #    lumiIncrement = 1
        #    setInitialLumi = True
        #  elif mass == 1850 and windowWidth == 11 and setInitialLumi==False:
        #    lumi = 3
        #    lumiIncrement = 1
        #    setInitialLumi = True
        #  elif mass == 1850 and windowWidth == 12 and setInitialLumi==False:
        #    lumi = 2
        #    lumiIncrement = 1
        #    setInitialLumi = True

          #------------------------------------------
          #print '\n\n******************************************'
          #print '******************************************'
          #print 'testing:'
          #print 'luminosity = %s ^fb-1'%lumi
          #print 'mass values [GeV]: %s'%massValues
          #print '\Testing mass:', mass
          #print '******************************************'
          #print '******************************************\n'

          #luminosity string
          slumi = ('%.1f'% float( str(lumi).replace('p','.'))).replace('.','p')

          print "what's the lumi", lumi

          #------------------------------------------
          #STEP 01 - get data-like QCD for the given luminosity
          #os.system('python -u step01.getDataLikeQCD.py --config %s --lumi %.1f --tag %s --patch --plot --batch --debug -b'%(args.configFileName, lumi, args.tag)) #HANNO: Commented out to skip

          #------------------------------------------
          #TEST
          #raise SystemExit('\n***TEST*** exit')
          #------------------------------------------

          #------------------------------------------
          #STEP 02 - inject signal (fast, do it anyway)

          #dataLikeQCDFileName = localdir+'/../results/datalikeQCD/datalikeQCD.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root' #HANNO: Fixed Lumi in line below!!!
          debugString = ""
          if args.debug :
            debugString = "--debug"
          #dataLikeQCDFileName = localdir+'/../results/datalikeQCD/NLO_yStar_06_29p7ifb_forSig.root' #HANNO: FIXED NAME
          #yvonne hack
          dataLikeQCDFileName = QCDFileName #HANNO: FIXED NAME
          os.system('python -u step02.injectDataLikeSignal.py --config %s --QCDFile %s --lumi %.1f --tag %s --plot --batch -b %s'%(args.configFileName, dataLikeQCDFileName, lumi, args.tag,debugString))
          print('python -u step02.injectDataLikeSignal.py --config %s --QCDFile %s --lumi %.1f --tag %s --plot --batch -b --debug'%(args.configFileName, dataLikeQCDFileName, lumi, args.tag))

          #------------------------------------------
          #TEST
          #raise SystemExit('\n***TEST*** exit')
          #------------------------------------------

          #signal mass
          print '\n******************************************'
          print 'Running SearchPhase'
          print '\n%s mass = %s GeV'%(model, mass)
          print 'lumi = %s fb^-1'%lumi
          print 'WindowWidth for SWiFt = %s'%windowWidth
          print '******************************************'
          #print("# of justabove files afte before eveyrhing : ", countFilesInDirWithKeyword(spFileDir, "JUSTABOVE"))

          #------------------------------------------
          #step 03 - search fordijet mass resonances (slow, only do it if needed)
          #signalPlusBackgroundFileName = localdir+'/../results2/signalplusbackground/signalplusbackground.'+model+'.'+slumi+'.ifb.'+"mjj_Gauss_sig__smooth"+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
          signalPlusBackgroundFileName=signalPlusBkgFileName(localdir, model, slumi, nPar, seed, tag=args.tag, signalName="mjj_Gauss_sig__smooth")
          print("YvonnesignalPlusBkgFileName: ", signalPlusBackgroundFileName)
          #signalPlusBackgroundFileName = localdir+'/../results/signalplusbackground/signalplusbackground.'+model+'.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'

          if lumi==0.0:
              nPseudoExps=nPseudoExps_bkg
          else:
              nPseudoExps=nPseudoExps_withSig

          os.system('python -u step03.searchPhase.py --config %s --file %s --mass %s --lumi %s --window %s --functionParam %s --nPseudoExps %s --tag %s --batch -b --debug'%(args.configFileName, signalPlusBackgroundFileName, int(mass), lumi, windowWidth, nPar, nPseudoExps, args.tag)) #Hanno: Removed --plot for now
          print('python -u step03.searchPhase.py --config %s --file %s --mass %s --lumi %s --window %s --functionParam %s --nPseudoExps %s --tag %s --batch -b --debug'%(args.configFileName, signalPlusBackgroundFileName, int(mass), lumi, windowWidth, nPar, nPseudoExps, args.tag)) #Hanno: Removed --plot for now

          #------------------------------------------
          #TEST
          #raise SystemExit('\n***TEST*** exit')
          #------------------------------------------

          #------------------------------------------
          #check SearchPhase results

          #this needs to be the same name as the file we made in the search phase - leave the 4param junk in for now
          #searchphase.Gauss_width15.620.GeV.500p0.ifb.mjj_Gauss.4.par.102.seed.NLO2015_29p7_sensitivityScan_ww13_case5Param
          ## these are commented out by yvonne
          #spfilename = localdir+'/../results/searchphase/searchphase.'+model+'.%i'%int(mass)+'.gev.'+slumi+'.ifb.'+histname+'.%i'%npar+'.par.%i'%seed+'.seed.'+args.tag+'_ww'+str(windowwidth)+'.root'
          #I believe this doesn't depend on mass
          spFileName = searchPhaseResultName(localdir, model, mass, slumi, histName, nPar, seed, args.tag, windowWidth)

          if not os.path.isfile(spFileName):
              raise SystemExit('\n***ERROR*** couldn\'t find SearchPhase output file for %s mass %s GeV: %s'%(model, int(mass), spFileName))

          spFile = ROOT.TFile(spFileName,'READ')
          spSignificance = spFile.Get('residualHist')#this is not quite the significance graph but it works the same for bin edges
          spSignificance.SetAxisRange( spSignificance.GetBinLowEdge( spSignificance.FindBin(2000.) ), 2e4, "X")

          #------------------------------------------
          #fill sensitivity scan graph and remove discovered signal mass values from the list
          bumpHunterStatOfFitToData = None
          try :
            bumpHunterStatOfFitToData = spFile.Get("bumpHunterStatOfFitToData")
          except :
            print "FIT FAILED!!! Try previous lumi with more pseudoexperiments"
            nPseudoExps = nPseudoExps+100
            lumi = lumi-1
            continue

          bumpHunterStatValue = bumpHunterStatOfFitToData[0]
          bumpHunterPValue    = bumpHunterStatOfFitToData[1]
          bumpHunterPValueErr = bumpHunterStatOfFitToData[2]

          bumpHunterPLowHigh = spFile.Get('bumpHunterPLowHigh')
          #bumpHunterStatValue = bumpHunterPLowHigh[0]
          bumpLowEdge         = bumpHunterPLowHigh[1]
          bumpHighEdge        = bumpHunterPLowHigh[2]

          excludeWindowVector = spFile.Get('excludeWindowNums')
          excludedWindow = excludeWindowVector[0]
          excludedWindowLow = excludeWindowVector[1]
          excludedWindowHigh = excludeWindowVector[2]

          print '\n******************************************'
          print 'SearchPhase results summary'
          print '\n%s mass = %s GeV'%(model, mass)
          print 'lumi = %s fb^-1'%lumi
          print "bump range: %s GeV - %s GeV"%(bumpLowEdge,bumpHighEdge)
          print "BumpHunter stat = %s"%bumpHunterStatValue
          print "BumpHunter p-value = %s +/- %s"%(bumpHunterPValue, bumpHunterPValueErr)
          print "excluded window (1=yes, 0=no): ", excludedWindow
          print "window low edge = %s"%excludedWindowLow
          print "window high edge = %s"%excludedWindowHigh

          bumpHunterSigmas = ROOT.Math.normal_quantile(1.-bumpHunterPValue, 1.)
          print "BumpHunter sigmas = %s"%bumpHunterSigmas
#--------removing old labelled files
          spFileDir=findDirFromFilePath(spFileName)
          print("spFileDir in code: ", spFileDir)
          #spFileDir="/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts/../results2/searchphase/"
          #print("spFileDir hardCoded: ", spFileDir)
          print ("mass", mass)
          print ("windowWidth", windowWidth)
          removeOldLabelledFile(spFileDir, "JUSTABOVE", mass, windowWidth)
          removeOldLabelledFile(spFileDir, "JUSTBELOW", mass, windowWidth)
          removeOldLabelledFile(spFileDir, "NOSIGNAL", mass, windowWidth)

          print("# of justabove files afte remove: ", countFilesInDirWithKeyword(spFileDir, "JUSTABOVE"))
          if bumpHunterPValue < 0.01 and int(excludedWindow)>0:
              massIndex = np.where(gmass == mass)
              glumi[massIndex] = lumi
              print '******************************************'
              #removeMassValues.append(mass)
              print "Discovery, with window removal"
              print "just above"
              print "SPFile: ", spFileName
              print "does GPFile exist: ", os.path.isfile(spFileName)
              #---finding the justabove file
              justAboveFN=justAboveFileName(spFileName)
              print("justAboveFileName: ", justAboveFN)
              os.rename(spFileName, justAboveFN)

            #setting the SP fileName previous to be one lumi step below exclusion window kicks in
              #---finding the justbelow file
              justBelowFN=justBelowFileName(spFileNamePrevious)
              os.rename(spFileNamePrevious, justBelowFN)
              print("justBelowFileName: ", justBelowFN)
              #---finding the nosignal file
              zeroLumiFN=zeroLumiFileName(spFileName)
              noSignalFN=noSignalFileName(zeroLumiFN)
              print("no signal file name: ", noSignalFN)
              print("zeoroLumi name:", zeroLumiFN)
              os.rename(zeroLumiFN,noSignalFN)

              print("# of justabove files:(after setting names) ", countFilesInDirWithKeyword(spFileDir, "JUSTABOVE"))


              break # get out of the lumi loop, go to next mass point
          #else :#increment lumi
          elif bumpHunterPValue < 0.01 and int(excludedWindow)==0:
              print '******************************************'
              print "Discovery, without window removal - increment lumi"
              lumi = lumi+lumiIncrement
          else:
            if lumi*bumpHunterPValue > 1:
              lumi = (lumi+2*lumi*bumpHunterPValue) #not even close, take five lumi steps, was 0.5 before
              print "Not too close. try next time with lumi: ", lumi
            else:
              lumi = lumi+lumiIncrement
              print "Close enough. try next time with lumi: ", lumi

            #setting a lumi file that is one step lower
          spFileNamePrevious=spFileName

            #------------------------------------------

            #------------------------------------------
            #TEST
            #raise SystemExit('\n***TEST*** exit')
            #------------------------------------------


        #------------------------------------------
        #print sensitivity scan results
          print '\n******************************************'
          print 'sensitivity scan results'
          print '******************************************\n'
          print 'mass = %s'%gmass
          print 'lumi = %s'%glumi

          print("# of justabove files close to end: ", countFilesInDirWithKeyword(spFileDir, "JUSTABOVE"))


      #------------------------------------------
      #NOTE the section for plotting the sensitivity scan results has been removed from here
      #     results can be plotted using the dedicated macro: plotSensitivityScan.py
      #     this macro needs some further work
#******************************************
if __name__ == '__main__':

    spFileDir="/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts/../results2/searchphase/"
    #------------------------------------------
    #parse input arguments
    parser = argparse.ArgumentParser(description='%prog [options]')
    parser.add_argument('--config', dest='configFileName', default='../configs/sensitivityScan.Test.config', required=True, help='sensitivity scan config file')
    parser.add_argument('--tag', dest='tag', default='default', help='tag for output files')
    parser.add_argument('-b', '--batch', dest='batch', action='store_true', default=False, help='batch mode')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    args = parser.parse_args()

    #------------------------------------------
    #do sensitivity scan
    print '******************************************'
    print 'sensitivity scan starting'
    doSensitivityScan(args)
    print '******************************************'
    print 'sensitivity scan done'
