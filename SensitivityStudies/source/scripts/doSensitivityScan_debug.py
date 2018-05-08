#!/bin/python

#******************************************
#perform sensitivity scan by looping over all masses on a range of luminosity values
#EXAMPLE python -u doSensitivityScan.py --config <config file> --tag <tag> --batch --debug

#******************************************
#import stuff
import sys, os, math, argparse, ROOT
import sensitivityTools
import numpy as np
import subprocess
#******************************************

def doSensitivityScan(args):
    
    fout.write('\n******************************************')
    fout.write('\n******************************************')

    #------------------------------------------
    #input parameters
    fout.write('\nSetting input parameters:')
    argsdict = vars(args)
    for ii in xrange(len(argsdict)):
        print '  %s = %s'%(argsdict.keys()[ii], argsdict.values()[ii],)
    
    #------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))
    #------------------------------------------
    #get settings
    fout.write('\n config settings:')
    settings = ROOT.TEnv()
    if settings.ReadFile(args.configFileName,ROOT.EEnvLevel(0)) != 0:
        raise SystemExit('***ERROR*** could not find sensitivity scan config file: %s'%args.configFileName)

    model = settings.GetValue('signalModel','')
    fout.write('\n  signal model = %s'%model)

    modelLabel = settings.GetValue('signalModelLabel','').replace('"','')
    fout.write('\n  signal model label = %s'%modelLabel)

    massValuesConfig = settings.GetValue('signalMasses','2000,3000,4000').split(',')
    massValuesConfig = [float(m) for m in massValuesConfig]
    fout.write('\n  signal masses [GeV] = %s'%massValuesConfig)
    
    lumiMin = float(settings.GetValue('luminosityMin','0.1'))
    #if lumiMin < 0.1: #fb^-1
    #    lumiMin = 0.1 #fb^-1
    fout.write('\n  minimum luminosity = %s'%lumiMin)

    lumiMax = float(settings.GetValue('luminosityMax','10.'))
    fout.write('\n  maximum luminosity = %s'%lumiMax)
    if lumiMax > 10000.: #fb^-1
        lumiMax = 10000. #fb^-1

    QCDFileName = settings.GetValue('QCDFile','../inputs/QCD/histograms.mc.dijet.1p0.ifb.root')
    fout.write('\n  QCD input file = %s'%QCDFileName)
    
    histBaseName = settings.GetValue('histBaseName','mjj')
    fout.write('\n  hist base name = %s'%histBaseName)

    bTaggingWP = settings.GetValue('bTaggingWP','') #fix_8585
    fout.write('\n  b-tagging WP = %s'%bTaggingWP)

    axisLabel = settings.GetValue('axisLabel','m [GeV]')
    fout.write('\n  hist x-axis label = %s'%axisLabel)

    nPar = int(settings.GetValue('nFitParameters','3'))
    fout.write('\n  n fit parameters = %s'%nPar)

    nPseudoExps = settings.GetValue('nPseudoExperiments','1000')
    fout.write('\n  number of pseudo-experiments = %s'%nPseudoExps)

    thresholdMass = float(settings.GetValue('thresholdMass','1100.'))
    fout.write('\n  threshold mass = %s'%thresholdMass)

    seed = float(settings.GetValue('randomSeed','0'))
    fout.write('\n  random seed = %s'%seed)

    notes = settings.GetValue('notes','').split(',')
    fout.write('\n  notes = %s'%notes+"\n")

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

    fout.write('\nDone setting initial parameters and checking inputs\n')
    fout.write('******************************************\n')
    fout.write('******************************************\n')


    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    for windowWidth in [9]:#add 9 later
    
      fout.write("******* Testing Window Width "+str(windowWidth)+" ******\n")
      
      #get list of available mass points
      massValuesAvailable = []
      fileList = os.listdir(localdir+'/../inputs/'+model+'/')
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
      fout.write('******* Using list of mass values [GeV]: %s'%massValues+"\n")
      #------------------------------------------
      #TEST
      #raise SystemExit('\n***TEST*** exit')
      #------------------------------------------

      #initial luminosity value, coming from config file
      lumi = lumiMin

      #------------------------------------------
      #arrays for sensitivity scan graphs
      gmass = np.array(massValues)
      glumi = np.zeros_like(gmass)

      #not used now, just upper bound
      lumiSteps = [0.1,0.2,0.3,0.5,0.7]+range(1,10)+range(10,20,1)+range(20,30,2)+range(30,50,3)+range(50,100,5)+range(100,200,10)+range(200,1000,50)
      #reset lumi to lowest value
      lumiStep = 0

      lumiIncrement = 1#this gets decided mass point by mass point
      #loop over mass values
      for mass in massValues : #and lumi <= lumiMax:
          
        setInitialLumi = False
        #print "******* Testing Mass "+str(mass)+" ******"

        while lumi<(lumiSteps[-1]) :

          if mass == 650 and windowWidth == 9 and setInitialLumi==False:
            lumi = 60
            setInitialLumi = True
          if mass == 675 and windowWidth == 9 and setInitialLumi==False:
              lumi = 40
              setInitialLumi = True
          if mass == 700 and windowWidth == 9 and setInitialLumi==False:
              lumi = 26
              setInitialLumi = True
          if mass == 725 and windowWidth == 9 and setInitialLumi==False:
              lumi = 20
              setInitialLumi = True
          if mass == 750 and windowWidth == 9 and setInitialLumi==False:
            lumi = 10
            setInitialLumi = True
          if mass == 1050 and windowWidth == 9 and setInitialLumi==False:
            lumi = 12
            setInitialLumi = True
          if mass == 1450 and windowWidth == 9 and setInitialLumi==False:
            lumi = 10
            setInitialLumi = True
          if mass == 1750 and windowWidth == 9 and setInitialLumi==False:
            lumi = 5
            setInitialLumi = True
          if mass == 1850 and windowWidth == 9 and setInitialLumi==False:
            lumi = 3
            setInitialLumi = True

          #luminosity string
          slumi = ('%.1f'% float( str(lumi).replace('p','.'))).replace('.','p')
          fout.write("******* Testing Mass, Lumi "+str(mass)+" GeV, "+str(lumi)+"/fb ******\n")

          #------------------------------------------
          #STEP 01 - get data-like QCD for the given luminosity
          #os.system('python -u step01.getDataLikeQCD.py --config %s --lumi %.1f --tag %s --patch --plot --batch --debug -b'%(args.configFileName, lumi, args.tag)) #HANNO: Commented out to skip

          #------------------------------------------
          #TEST
          #raise SystemExit('\n***TEST*** exit')
          #------------------------------------------
          
          #------------------------------------------
          fout.write('******************************************\n')
          fout.write('******************************************\n')
          fout.write('Injecting signal\n')

          #STEP 02 - inject signal (fast, do it anyway)
          #dataLikeQCDFileName = localdir+'/../results/datalikeQCD/datalikeQCD.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root' #HANNO: Fixed Lumi in line below!!!
          debugString = ""
          if args.debug :
            debugString = "--debug"
          #dataLikeQCDFileName = localdir+'/../results/datalikeQCD/NLO_yStar_06_29p7ifb_forSig.root' #HANNO: FIXED NAME
          dataLikeQCDFileName = localdir+'/../results/datalikeQCD/Pseudodata_from_DSJ100yStar06_TriggerJets_J100_yStar06_mjj_2016binning_TLArange_data_4param_G_upTo4000.29p7.ifb.root' #CD: FIXED NAME
          #why calling python within python....but not doing so as it would require major rewrite
          print('python -u step02.injectDataLikeSignal.py --config %s --QCDFile %s --lumi %.1f --tag %s --plot --batch -b --debug'%(args.configFileName, dataLikeQCDFileName, lumi, args.tag))
          os.system('python -u step02.injectDataLikeSignal.py --config %s --QCDFile %s --lumi %.1f --tag %s --plot --batch -b %s'%(args.configFileName, dataLikeQCDFileName, lumi, args.tag,debugString))
          #processOutput = ""
          #try :
          #  processOutput = subprocess.check_output('python -u step02.injectDataLikeSignal.py --config %s --QCDFile %s --lumi %.1f --tag %s --plot --batch -b %s'%(args.configFileName, dataLikeQCDFileName, lumi, args.tag,debugString), shell=True)
          #except :
          #  print "Signal-injected files existed already!"
          #print processOutput

          fout.write('\nDone injecting signal\n')
          fout.write('******************************************\n')
          fout.write('******************************************\n')

          #lumi = lumiSteps[-1]
          #------------------------------------------
          #TEST
          #raise SystemExit('\n***TEST*** exit')
          #------------------------------------------

          #run search phase for this signal mass
          fout.write('******************************************\n')
          fout.write('******************************************\n')
          fout.write('Running SearchPhase\n')
          fout.write('\n%s mass = %s GeV'%(model, mass)+'\n')
          fout.write('lumi = %s fb^-1'%lumi+'\n')
          fout.write('WindowWidth for SWiFt = %s'%windowWidth+'\n')
          fout.write('******************************************\n')
          fout.write('******************************************\n')

          #------------------------------------------
          #step 03 - search for dijet mass resonances (slow, only do it if needed)
          signalPlusBackgroundFileName = localdir+'/../results/signalplusbackground/signalplusbackground.'+model+'.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
          print('python -u step03.searchPhase.py --config %s --file %s --mass %s --lumi %s --window %s --functionParam %s --nPseudoExps %s --tag %s --batch -b --debug'%(args.configFileName, signalPlusBackgroundFileName, int(mass), lumi, windowWidth, nPar, nPseudoExps, args.tag)) #Hanno: Removed --plot for now
          #try :
          #  process = subprocess.check_output('python -u step03.searchPhase.py --config %s --file %s --mass %s --lumi %s --window %s --functionParam %s --nPseudoExps %s --tag %s --batch -b --debug'%(args.configFileName, signalPlusBackgroundFileName, int(mass), lumi, windowWidth, nPar, nPseudoExps, args.tag), shell=True)
          #except :
          #  raise SystemExit('\n***Something went wrong with the search phase***')

          os.system('python -u step03.searchPhase.py --config %s --file %s --mass %s --lumi %s --window %s --functionParam %s --nPseudoExps %s --tag %s --batch -b --debug'%(args.configFileName, signalPlusBackgroundFileName, int(mass), lumi, windowWidth, nPar, nPseudoExps, args.tag)) #Hanno: Removed --plot for now
          #------------------------------------------
          #check SearchPhase results
          
          #this needs to be the same name as the file we made in the search phase - leave the 4param junk in for now
          #searchphase.Gauss_width15.620.GeV.500p0.ifb.mjj_Gauss.4.par.102.seed.NLO2015_29p7_sensitivityScan_ww13_case5Param
          spFileName = localdir+'/../results/searchphase/searchphase.'+model+'.%i'%int(mass)+'.GeV.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'_ww'+str(windowWidth)+'.root'
          
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

          fout.write( '\n******************************************')
          fout.write( 'SearchPhase results summary\n')
          fout.write( '\n%s mass = %s GeV'%(model, mass)+'\n')
          fout.write( 'lumi = %s fb^-1'%lumi+'\n')
          fout.write( "bump range: %s GeV - %s GeV"%(bumpLowEdge,bumpHighEdge)+'\n')
          fout.write( "BumpHunter stat = %s"%bumpHunterStatValue+'\n')
          fout.write( "BumpHunter p-value = %s +/- %s"%(bumpHunterPValue, bumpHunterPValueErr)+'\n')
          fout.write( "excluded window (1=yes, 0=no): "+str(excludedWindow)+'\n')
          fout.write( "window low edge = %s"%excludedWindowLow+'\n')
          fout.write("window high edge = %s"%excludedWindowHigh+'\n')

          bumpHunterSigmas = ROOT.Math.normal_quantile(1.-bumpHunterPValue, 1.)
          fout.write("BumpHunter sigmas = %s"%bumpHunterSigmas+'\n')
          
          #if bumpHunterSigmas > 5.:
          if bumpHunterPValue < 0.01 and int(excludedWindow)>0:
              massIndex = np.where(gmass == mass)
              glumi[massIndex] = lumi
              fout.write("******************************************"+'\n')
              #removeMassValues.append(mass)
              fout.write("Discovery, with window removal"+'\n')
              break # get out of the lumi loop, go to next mass point
          #else :#increment lumi
          elif bumpHunterPValue < 0.01 and int(excludedWindow)==0:
              #massIndex = np.where(gmass == mass)
              #glumi[massIndex] = lumi
              fout.write("******************************************"+'\n')
              #removeMassValues.append(mass)
              fout.write("Discovery, without window removal - increment lumi"+'\n')
              lumi = lumi+lumiIncrement
          else:
            if lumi*bumpHunterPValue > 1:
              lumi = (lumi+1*lumi*bumpHunterPValue) #not even close, take five lumi steps, was 0.5 before
              fout.write( "Not too close. try next time with lumi: "+str(lumi)+'\n')
            else:
              lumi = lumi+lumiIncrement
              fout.write( "Close enough. try next time with lumi: "+str(lumi)+'\n')

        #------------------------------------------
        #print sensitivity scan results
          fout.write( '\n******************************************\n')
          fout.write( 'sensitivity scan results for all points so far, discovery happened at:\n')
          fout.write( '******************************************\n')
          fout.write( 'mass = %s'%gmass )
          fout.write( 'lumi = %s'%glumi )

        #lumi = lumiSteps[-1]

#******************************************
if __name__ == '__main__':

    #------------------------------------------
    #parse input arguments
    parser = argparse.ArgumentParser(description='%prog [options]')
    parser.add_argument('--config', dest='configFileName', default='', required=True, help='sensitivity scan config file')
    parser.add_argument('--tag', dest='tag', default='default', help='tag for output files')
    parser.add_argument('-b', '--batch', dest='batch', action='store_true', default=False, help='batch mode')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    args = parser.parse_args()
    fout = open(args.configFileName+"_log.txt","w")


    #------------------------------------------
    #do sensitivity scan
    print '******************************************'
    print 'sensitivity scan starting'
    doSensitivityScan(args)
    print '******************************************'
    print 'sensitivity scan done'

