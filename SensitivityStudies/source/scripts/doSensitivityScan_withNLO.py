#!/bin/python

#******************************************
#perform sensitivity scan by looping over all masses on a range of luminosity values
#EXAMPLE python -u doSensitivityScan.py --config <config file> --tag <tag> --batch --debug

#******************************************
#import stuff
import sys, os, math, argparse, ROOT
import sensitivityTools
import numpy as np

#******************************************
def doSensitivityScan(args):

    print '\n******************************************'
    print 'sensitivity scan'
    
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
    if lumiMin < 0.1: #fb^-1
        lumiMin = 0.1 #fb^-1
    print '  minimum luminosity = %s'%lumiMin

    lumiMax = float(settings.GetValue('luminosityMax','10.'))
    print '  maximum luminosity = %s'%lumiMax
    if lumiMax > 1000.: #fb^-1
        lumiMax = 1000. #fb^-1

    QCDFileName = settings.GetValue('QCDFile','../inputs/QCD/histograms.mc.dijet.1p0.ifb.root')
    print '  QCD input file = %s'%QCDFileName
    
    histBaseName = settings.GetValue('histBaseName','mjj')
    print '  hist base name = %s'%histBaseName

    bTaggingWP = settings.GetValue('bTaggingWP','') #fix_8585
    print '  b-tagging WP = %s'%bTaggingWP

    axisLabel = settings.GetValue('axisLabel','m [GeV]')
    print '  hist x-axis label = %s'%axisLabel

    nPar = int(settings.GetValue('nFitParameters','3'))
    print '  n fit parameters = %s'%nPar

    nPseudoExps = settings.GetValue('nPseudoExperiments','1000')
    print '  number of pseudo-experiments = %s'%nPseudoExps

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
    print '\navailable mass values [GeV]: %s'%massValuesAvailable

    massValues = list( set(massValuesConfig) & set(massValuesAvailable) )
    massValues.sort(key=float)
    print 'using mass values [GeV]: %s'%massValues
    
    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    #------------------------------------------
    #initial luminosity value
    lumi = lumiMin

    #------------------------------------------
    #arrays for sensitivity scan graphs
    gmass = np.array(massValues)
    glumi = np.zeros_like(gmass)
    
    #------------------------------------------
    #loop over mass values and luminosity
    while len(massValues) > 0 and lumi <= lumiMax:

        #------------------------------------------
        print '\n\n******************************************'
        print '******************************************'
        print 'luminosity = %s ^fb-1'%lumi
        print 'mass values [GeV]: %s'%massValues
        print '******************************************'
        print '******************************************\n'

        #luminosity string
        slumi = ('%.1f'% float( str(lumi).replace('p','.'))).replace('.','p')
        
        #------------------------------------------
        #STEP 01 - get data-like QCD for the given luminosity
        #os.system('python -u step01.getDataLikeQCD.py --config %s --lumi %.1f --patch --tag %s --plot'%(args.configFileName, lumi, args.tag))

        #------------------------------------------
        #TEST
        #raise SystemExit('\n***TEST*** exit')
        #------------------------------------------
        
        #------------------------------------------
        #STEP 02 - inject signal
        #dataLikeQCDFileName = localdir+'/../results/datalikeQCD/datalikeQCD.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
        dataLikeQCDFileName = localdir+'/../results/datalikeQCD/NLO_yStar_06_29p7ifb_forSig.root'
        os.system('python -u step02.injectDataLikeSignal.py --config %s --QCDFile %s --lumi %.1f --tag %s --plot'%(args.configFileName, dataLikeQCDFileName, lumi, args.tag))

        #------------------------------------------
        #TEST
        #raise SystemExit('\n***TEST*** exit')
        #------------------------------------------
        
        #------------------------------------------
        #loop over signal mass values
        removeMassValues = []
        for mass in massValues:

            #------------------------------------------
            #signal mass
            print '\n******************************************'
            print 'SearchPhase'
            print '\n%s mass = %s GeV'%(model, mass)
            print 'lumi = %s fb^-1'%lumi
            print '******************************************'

            #------------------------------------------
            #step 03 - search for dijet mass resonances
            signalPlusBackgroundFileName = localdir+'/../results/signalplusbackground/signalplusbackground.'+model+'.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
            os.system('python -u step03.searchPhase.py --config %s --file %s --mass %s --lumi %s --tag %s --plot'%(args.configFileName, signalPlusBackgroundFileName, int(mass), lumi, args.tag))

            #------------------------------------------
            #TEST
            #raise SystemExit('\n***TEST*** exit')
            #------------------------------------------

            #------------------------------------------
            #check SearchPhase results
            print '\n******************************************'
            print 'SearchPhase results summary'
            print '\n%s mass = %s GeV'%(model, mass)
            print 'lumi = %s fb^-1'%lumi

            spFileName = localdir+'/../results/searchphase/searchphase.'+model+'.%i'%int(mass)+'.GeV.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
                
            if not os.path.isfile(spFileName):
                raise SystemExit('\n***ERROR*** couldn\'t find SearchPhase output file for %s mass %s GeV: %s'%(model, int(mass), spFileName))

            spFile = ROOT.TFile(spFileName,'READ')
            spSignificance = spFile.Get('residualHist')
            spSignificance.SetAxisRange( spSignificance.GetBinLowEdge( spSignificance.FindBin(2000.) ), 2e4, "X")

            #------------------------------------------
            #fill sensitivity scan graph and remove discovered signal mass values from the list
            bumpHunterStatOfFitToData = spFile.Get("bumpHunterStatOfFitToData")
            bumpHunterStatValue = bumpHunterStatOfFitToData[0]
            bumpHunterPValue    = bumpHunterStatOfFitToData[1]
            bumpHunterPValueErr = bumpHunterStatOfFitToData[2]

            bumpHunterPLowHigh = spFile.Get('bumpHunterPLowHigh')
            #bumpHunterStatValue = bumpHunterPLowHigh[0]
            bumpLowEdge         = bumpHunterPLowHigh[1]
            bumpHighEdge        = bumpHunterPLowHigh[2]

            print "bump range: %s GeV - %s GeV"%(bumpLowEdge,bumpHighEdge)
            print "BumpHunter stat = %s"%bumpHunterStatValue
            print "BumpHunter p-value = %s +/- %s"%(bumpHunterPValue, bumpHunterPValueErr)

            bumpHunterSigmas = ROOT.Math.normal_quantile(1.-bumpHunterPValue, 1.)
            print "BumpHunter sigmas = %s"%bumpHunterSigmas
            
            if bumpHunterSigmas > 5.:
                massIndex = np.where(gmass == mass)
                glumi[massIndex] = lumi
                removeMassValues.append(mass)

            print '******************************************'
            #------------------------------------------

            #------------------------------------------
            #TEST
            #raise SystemExit('\n***TEST*** exit')
            #------------------------------------------

        #------------------------------------------
        #remove mass points discovered
        for removeMass in removeMassValues:
            massValues.remove(removeMass)
        print '\n******************************************'
        print 'available mass values [GeV]: %s'%massValues
        print '******************************************'

        #------------------------------------------
        #increase luminosity
        lumiSteps = [0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0, 200.0, 300.0, 500.0, 1000.0]
        #lumiSteps = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0] #SHORT
        if lumi < lumiSteps[0]:
            lumi = lumiSteps[0]
        elif lumi >= lumiSteps[-1]:
            lumi*=2.
        else:
            for ii in xrange( len(lumiSteps)-1):
                if lumi >= lumiSteps[ii] and lumi < lumiSteps[ii+1]:
                    lumi = lumiSteps[ii+1]
                    break

    #------------------------------------------
    #print sensitivity scan results
    print '\n******************************************'
    print 'sensitivity scan results'
    print '******************************************\n'
    print 'mass = %s'%gmass
    print 'lumi = %s'%glumi

    #------------------------------------------
    #NOTE the section for plotting the sensitivity scan results has been removed from here
    #     results can be plotted using the dedicated macro: plotSensitivityScan.py
    #     this macro needs some further work


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

    #------------------------------------------
    #do sensitivity scan    
    doSensitivityScan(args)
    print '\n******************************************'
    print 'sensitivity scan done'
