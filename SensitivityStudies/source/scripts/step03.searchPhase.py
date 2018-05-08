#!/bin/python

#******************************************
#perform search on given input file and signal mass value
#EXAMPLE python -u step03.searchPhase.py --config <config file> --file <file> --mass <mass> --lumi <luminosity [fb^-1]> --fit --tag <tag> --window <lumi> --plot --batch --debug

#******************************************
#write output file to run
runFile = open("run.sh", "w")

#******************************************
#import stuff
import sys, os, math, argparse, ROOT
import searchPhaseTools
from subprocess import call
from fileNamingTool import *

#******************************************
def runSearchPhase(args):

    print '\n******************************************'
    print 'run search phase'

    #******************************************
    #set ATLAS style
    if os.path.isfile(os.path.expanduser('~/RootUtils/AtlasStyle.C')):
        ROOT.gROOT.LoadMacro('~/RootUtils/AtlasStyle.C')
        ROOT.SetAtlasStyle()
        #ROOT.set_color_env()
    else:
        print '\n***WARNING*** couldn\'t find ATLAS Style'
        #import AtlasStyle
        #AtlasStyle.SetAtlasStyle()

    #------------------------------------------
    #set error sum and overflow
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows()
    ROOT.TH2.SetDefaultSumw2()
    ROOT.TH2.StatOverflows()

    #------------------------------------------
    #input parameters
    print '\ninput parameters:'
    argsdict = vars(args)
    for ii in xrange(len(argsdict)):
        print '  %s = %s'%(argsdict.keys()[ii], argsdict.values()[ii],)

    slumi = ('%.1f'% float( str(args.lumi).replace('p','.'))).replace('.','p')

    #------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #get settings
    print '\nconfig settings:'
    settings = ROOT.TEnv()
    if settings.ReadFile(args.configFileName,ROOT.EEnvLevel(0)) != 0:
        raise IOError('could not find sensitivity scan config file: %s'%args.configFileName)

    model = settings.GetValue('signalModel','')
    print '  signal model = %s'%model

    modelLabel = settings.GetValue('signalModelLabel','').replace('"','')
    print '  signal model label = %s'%modelLabel

    histBaseName = settings.GetValue('histBaseName','mjj')
    print '  hist base name = %s'%histBaseName

    bTaggingWP = settings.GetValue('bTaggingWP','') #fix_8585
    print '  b-tagging WP = %s'%bTaggingWP

    axisLabel = settings.GetValue('axisLabel','m [GeV]')
    print '  hist x-axis label = %s'%axisLabel

    nPar = int(settings.GetValue('nFitParameters','3'))
    print '  n fit parameters = %s'%nPar

    thresholdMass = float(settings.GetValue('thresholdMass','1100.'))
    print '  threshold mass = %s'%thresholdMass

    seed = float(settings.GetValue('randomSeed','0'))
    print '  random seed = %s'%seed

    configNotes = settings.GetValue('notes','').split(',')
    print '  notes = %s'%configNotes

	#------------------------------------------
    #set variables
    slumi = ('%.1f'% float( str(args.lumi).replace('p','.'))).replace('.','p')

    histName = histBaseName
    if bTaggingWP != '':
        histName+='_'+bTaggingWP
    if args.debug:
        print '\nhist name = %s'%histName

    textSize=20

	#------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #check data-like QCD file
    if not os.path.isfile(args.inputFileName):
        raise SystemExit('\n***ERROR*** couldn\'t find input file: %s'%args.inputFileName)

    #------------------------------------------
    #check luminosity
    if not slumi+'.ifb.' in args.inputFileName:
        raise SystemExit('\n***ERROR*** is the lumi value right?')

    #------------------------------------------
    #check model
    if not model in args.inputFileName:
        raise SystemExit('\n***ERROR*** is the model name right?')

    #------------------------------------------
    #check output file
    #outFileName = localdir+'/../results2/searchphase/searchphase.'+model+'.%s'%int(args.mass)+'.GeV.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
    outFileName=searchPhaseResultName(localdir, model, args.mass, slumi, histName, nPar, seed, args.tag, args.window)
    print outFileName
#    if os.path.isfile(outFileName):
#        raise SystemExit('\n***WARNING*** output file exists already: %s'%outFileName)
#
    #------------------------------------------
    #get fit initial parameters
    #print 'get fit initial parameters'

    pars = []

    if args.fit:
        #do a quick fit to get initial fit parameters
        pars = searchPhaseTools.simpleFit(localdir+'/'+args.inputFileName, '', histName+'_%s'%int(args.mass), thresholdMass, 13e3, nPar, False)
        print '\ninitial parameters from quick fit: %s'%pars
        #------------------------------------------
        #TEST
        #raise SystemExit('\n***TEST*** exit')
        #------------------------------------------

    else:
        #predefined fit parameters
        #NOTE the available values are limited *FOR NOW* to those obtained from the fits of 1 fb^-1 mjj spectrum

       # #5 parameters fit
       # if nPar == 5:
       #     if args.lumi <= 1.: #fb^-1
       #         print '\n1 fb^-1 mjj fit 5 parameters'
       #         pars = [0.242016, 10.627, 1.40793, -1.46081, -0.197971]
       #     else:
       #         print '\n1 fb^-1 mjj fit 5 parameters'
       #         pars = [0.242016, 10.627, 1.40793, -1.46081, -0.197971]

       # #4 parameters fit
       # if nPar == 4:
       #     if args.lumi <= 1.: #fb^-1
       #         print '\n1 fb^-1 mjj fit 4 parameters'
       #         pars = [0.00280214, 7.2046, 5.63296, 0.0571367]
       #     else:
       #         print '\n1 fb^-1 mjj fit 4 parameters'
       #         pars = [0.00280214, 7.2046, 5.63296, 0.0571367]

       # #3 parameters fit
       # else:
       #     if args.lumi <= 1.: #fb^-1
       #         print '\n1 fb^-1 mjj fit 3 parameters'
       #         pars = [0.00458656, 7.75728, -5.31253]
       #     else:
       #         print '\n1 fb^-1 mjj fit 3 parameters'
       #         pars = [0.00458656, 7.75728, -5.31253]

       # print 'initial parameters: %s'%pars
       #Yvonne Frustratingly hardcoded this
       pars=[2114090.3221715163, -0.10249143346059597, 110.96753410590006, -41.90626047681503]



    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------
    # 2016 data
    fitFile = args.inputFileName
    print("YvonneFitFileName: ", fitFile)
    #Yvonne hardcoding this.....
    print "slumi", slumi
    if slumi=="0p0":
        histName=histBaseName+"_QCD"
    else:
        histName="mjj_Gauss_sig_{0}_smoothinjectedToBkg".format(str(args.mass))
    #histName="basicBkgFrom4ParamFit_fluctuated_QCD"
    fitHistogram = histName
    lowEstimate = 300
    highEstimate = 1500
    highFit = 7500
    doAlternate = False
    permitWindow = False
    case = args.functionParam+"Param" #hardwired for now, can get from args
    #outFileName = outFileName[0:len(outFileName)-5]
    #outFileName = outFileName+"_ww"+str(args.window)+".root"
    #not sure this is needed?
    statspath = os.getcwd() # path used in outputFileName in config
    headdir = statspath.split("/Bayesian")[0] # directory for whole package

    workTag = "sensitivity_"+histName+'_'+str(int(args.mass))+"_ww"+args.window
    outDir = "results/%s/"%workTag

    commandTemplate = "SearchPhase --config {0} --file {1} --histName {2} --noDE --outputfile {3}".format("{}",fitFile,fitHistogram,outFileName)
    scriptArchive = "submitConfigs/%s/"%workTag
    for directory in [outDir,scriptArchive]:
      if not os.path.exists(directory): os.makedirs(directory)

    configInName = "configurations/Step1_SearchPhase_Swift.config"
    configOutName = "submitConfigs/"+workTag+"/Step1_SearchPhase_case{0}_window{1}.config".format(case,args.window)
    configOut = open(configOutName,'w')
    with open(configInName) as configInData :
      for line in configInData :
        if "doSwift" in line: #always do Swift
          line = "doSwift true\n"
        if "minXForFit" in line :
          line = "minXForFit   {0}\n".format(lowEstimate)
        if "maxXForFit" in line :
          line = "maxXForFit   {0}\n".format(highEstimate)
        if "swift_minXAvailable"  in line :
          line = "swift_minXAvailable   {0}\n".format(lowEstimate)
        if "swift_maxXAvailable"  in line :
          line = "swift_maxXAvailable   {0}\n".format(highFit)
#        if "inputHistDir"  in line and fitHistogram == "mjj" :
#          continue

        if doAlternate and "doAlternateFunction" in line :
          line = "doAlternateFunction  true\n"
        if permitWindow and "permitWindow" in line :
          line = "permitWindow  false\n"
        if "nPseudoExp" in line :
          line = "nPseudoExp  {0}\n".format(args.nPseudoExps)

        if "swift_nBinsLeft" in line :
          line = "swift_nBinsLeft  {0}\n".format(args.window)
        elif "swift_nBinsRight" in line :
          line = "swift_nBinsRight  {0}\n".format(args.window)
        if "4Param" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  true\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  true\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       false\n"
          if "functionCode" in line :
            line = "functionCode  1\n" #seriously why are we hardcoding everything
          if "nParameters" in line :
            line = "nParameters  4\n"
          if "alternateFunctionCode" in line :
            line = "alternateFunctionCode  4\n"
          if "alternateNParameters" in line :
            line = "alternateNParameters  4\n"

        if "UA2" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  true\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  true\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       false\n"
          if "functionCode" in line :
            line = "functionCode  1\n"
          if "nParameters" in line :
            line = "nParameters  4\n"
          if "alternateFunctionCode" in line :
            line = "alternateFunctionCode  1\n"
          if "alternateNParameters" in line :
            line = "alternateNParameters  4\n"

        if "5Param" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  true\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  true\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       false\n"
          if "functionCode" in line :
            line = "functionCode  7\n"
          if "nParameters" in line :
            line = "nParameters  5\n"
          if "alternateFunctionCode" in line :
            line = "alternateFunctionCode  7\n"
          if "alternateNParameters" in line :
            line = "alternateNParameters  5\n"

        if "5ParamUA2Log" in case :
          if "swift_fixLow" in line :
            line = "swift_fixLow  true\n"
          if "swift_fixHigh" in line :
            line = "swift_fixHigh  true\n"
          if "swift_truncateHigh" in line :
            line = "swift_truncateHigh       false\n"
          if "functionCode" in line :
            line = "functionCode  24\n"
          if "nParameters" in line :
            line = "nParameters  5\n"
          if "alternateFunctionCode" in line :
            line = "alternateFunctionCode  24\n"
          if "alternateNParameters" in line :
            line = "alternateNParameters  5\n"

        configOut.write(line)
    configOut.close()
    #???
    modcommand = "chmod 744 {0}".format(configOutName)
    print("YvonneConfigOutName: ", configOutName)
    #call(modcommand,shell=True)

    os.system(modcommand)
    thisCommand = commandTemplate.replace("{}",configOutName)

#    #------------------------------------------
#    #create config file
#    configFileName    = localdir+'/../data/searchPhase.config'
#    newConfigFileName = localdir+'/../configs/searchphase.'+model+'.%s'%int(args.mass)+'.GeV.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.config'
#    print 'config file = %s'%newConfigFileName
#
#    with open(configFileName,'r') as configFile:
#        with open(newConfigFileName,'w') as newConfigFile:
#            for line in configFile:
#                newLine=line
#                newLine = newLine.replace('dummyInputFileName', args.inputFileName)
#                if len(args.inputHistDir)>0:
#                    newLine = newLine.replace('#inputHistDir', 'inputHistDir\t\t'+args.inputHistDir)
#                newLine = newLine.replace('dummyHistName', histName+'_'+str(int(args.mass)))
#                newLine = newLine.replace('dummyOutputFileName', outFileName)
#                newLine = newLine.replace('dummyMinX', str(thresholdMass))
#
#                #3 parameters
#                if nPar == 3:
#                    newLine = newLine.replace('dummyFuncCode', str(9))
#                    newLine = newLine.replace('dummyNPar', str(3))
#                    newLine = newLine.replace('dummyP1', str(pars[0]))
#                    newLine = newLine.replace('dummyP2', str(pars[1]))
#                    newLine = newLine.replace('dummyP3', str(pars[2]))
#
#                #4 parameters
#                elif nPar == 4:
#                    newLine = newLine.replace('dummyFuncCode', str(4))
#                    newLine = newLine.replace('dummyNPar', str(4))
#                    newLine = newLine.replace('dummyP1', str(pars[0]))
#                    newLine = newLine.replace('dummyP2', str(pars[1]))
#                    newLine = newLine.replace('dummyP3', str(pars[2]))
#                    newLine = newLine.replace('#parameter4', 'parameter4\t\t'+str(pars[3]))
#
#                #5 parameters
#                elif nPar == 5:
#                    newLine = newLine.replace('dummyFuncCode', str(7))
#                    newLine = newLine.replace('dummyNPar', str(5))
#                    newLine = newLine.replace('dummyP1', str(pars[0]))
#                    newLine = newLine.replace('dummyP2', str(pars[1]))
#                    newLine = newLine.replace('dummyP3', str(pars[2]))
#                    newLine = newLine.replace('#parameter4', 'parameter4\t\t'+str(pars[3]))
#                    newLine = newLine.replace('#parameter5', 'parameter5\t\t'+str(pars[4]))
#
#                newConfigFile.write(newLine)

    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    #------------------------------------------
    #run Bump Hunter
    #------------------------------------------
    print '\n******************************************'
    print 'setting up to run SearchPhase\n'
    runFile.write('SearchPhase --config %s --noDE'%configOutName)
    #print "SearchPhase --config %s --noDE"%configOutName+"| tee "+outDir+"/"+workTag+"_log.txt \n"
    #call(["SearchPhase --config %s --noDE"%configOutName+"| tee "+outDir+"/"+workTag+"_log.txt \n"], shell=True)
    print "Search Phase Run command", thisCommand+"| tee "+outDir+"/"+workTag+"_log.txt \n"
    #call([thisCommand+"| tee "+outDir+"/"+workTag+"_log.txt \n"], shell=True)
    os.system(thisCommand+"| tee "+outDir+"/"+workTag+"_log.txt \n")

    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    #plot search phase results
    snotes = '" "'.join(configNotes)
    print "To plot:"
    plotCommand = 'python plotSearchPhase.py --bump --file %s --lumi %s --notes "m_{%s} = %s GeV" "%s" --xlabel "%s"'%(outFileName, args.lumi, modelLabel.replace('\\',''), int(args.mass), snotes, axisLabel)
    print plotCommand
    if args.wait:
        plotCommand+=' --wait'
    if args.plot:
        print plotCommand #DEBUG
        os.system(plotCommand)
    runFile.write(plotCommand)


#******************************************
if __name__ == '__main__':

    #------------------------------------------
    #parse input arguments
    parser = argparse.ArgumentParser(description='%prog [options]')
    parser.add_argument('--config', dest='configFileName', default='', required=True, help='sensitivity scan config file')
    parser.add_argument('--file', dest='inputFileName', default='', required=True, help='input file')
    parser.add_argument('--dir', dest='inputHistDir', default='', help='input hist directory')
    parser.add_argument('--mass', dest='mass', default=0, required=True, help='signal mass value')
    parser.add_argument('--lumi', dest='lumi', default=0., required=True, help='luminosity [fb^-1]')
    parser.add_argument('--window', dest='window', required=True, help='SWiFt window')
    parser.add_argument('--nPseudoExps', dest='nPseudoExps', required=True, help='Number of PE')
    parser.add_argument('--functionParam', dest='functionParam', required=True, help='Number of function parameters')
    parser.add_argument('--fit', dest='fit', action='store_true', default=False, help='run a quick fit to get initial parameters?')
    parser.add_argument('--tag', dest='tag', default='default', help='tag for output files')
    parser.add_argument('--wait', dest='wait', action='store_true', default=False, help='wait?')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true', default=False, help='plot histograms')
    parser.add_argument('-b', '--batch', dest='batch', action='store_true', default=False, help='batch mode')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    args = parser.parse_args()

    #------------------------------------------
    #run search phase
    runSearchPhase(args)
    print '\nran search phase'
    print '******************************************'
