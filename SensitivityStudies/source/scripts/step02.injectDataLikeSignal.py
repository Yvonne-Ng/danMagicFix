#!/bin/python

#******************************************
#script to produce sig+bkg MC distributions ready for SearchPhase
#NOTE all distributions are normalized to the desired luminosity
#NOTE the signal is data-like, if not enough events are vailable, the signal is not injected
#EXAMPLE python -u step02.injectDataLikeSignal.py --config <config file> --QCDFile <QCD file> --lumi <luminosity [fb^-1]> --tag <tag> --fixCS --wait --plot --batch --debug

#******************************************
#import stuff
import sys, os, math, argparse, ROOT
import plotTools, sensitivityTools

#******************************************
def injectDataLikeSignal(args):

    print '\n******************************************'
    print 'inject data-like signal'

    #******************************************
    #set ATLAS style
    #Yvonne hacking this
    #using the python version
#just import the python one later
    if os.path.isfile(os.path.expanduser("/lustre/SCRATCH/atlas/ywng/WorkSpace/signalInjection/20171122_SensitivityScan/PlotSensitivity/RootStyle/AtlasStyle.C")):
        ROOT.gROOT.LoadMacro('/lustre/SCRATCH/atlas/ywng/WorkSpace/signalInjection/20171122_SensitivityScan/PlotSensitivity/RootStyle/AtlasStyle.C')
        #ROOT.set_color_env()

    #if os.path.isfile(os.path.expanduser('~/RootUtils/AtlasStyle.C')):
    #    ROOT.gROOT.LoadMacro('~/RootUtils/AtlasStyle.C')
    #    ROOT.SetAtlasStyle()
    #    #ROOT.set_color_env()
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

    massValuesConfig = settings.GetValue('signalMasses','2000,3000,4000').split(',')
    massValuesConfig = [float(m) for m in massValuesConfig]
    print '  signal masses [GeV] = %s'%massValuesConfig

    histBaseNameBkg = settings.GetValue('histBaseNameBkg','mjj')
    print '  hist base name = %s'%histBaseNameBkg

    histBaseNameSig = settings.GetValue('histBaseNameSig','mjj')
    print '  hist base name = %s'%histBaseNameSig

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

    histNameSig = histBaseNameSig

    histNameBkg = histBaseNameBkg
    print("Step2 histNameBkg: ", histNameBkg)

    if bTaggingWP != '':
        histNameBkg+='_'+bTaggingWP
    if args.debug:
        print '\nhist name = %s'%histNameBkg

    textSize=20

	#------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #check data-like QCD file
    if not os.path.isfile(args.dataLikeQCDFileName):
        raise KeyboardInterrupt('\n***ERROR*** couldn\'t find data-like QCD file: %s'%args.dataLikeQCDFileName)

    #------------------------------------------
    #check luminosity
    #if not slumi+'.ifb.' in args.dataLikeQCDFileName: #HANNO: commented out
    #    raise SystemExit('\n***ERROR*** is the lumi value right?')

    #------------------------------------------
    #check output file
    #outFileName = localdir+'/../results/signalplusbackground/signalplusbackground.'+model+'.'+slumi+'.ifb.'+histNameSig+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
    histNameFromInput=histNameSig.format("")
    outFileName = localdir+'/../results2/signalplusbackground/signalplusbackground.'+model+"."+slumi+'.ifb.'+histNameFromInput+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
    #if os.path.isfile(outFileName):
    #    raise SystemExit('\n***WARNING*** output file exists already: %s'%outFileName)
    outFile = ROOT.TFile(outFileName, 'RECREATE')

    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    #------------------------------------------
    #data-like QCD hist
    dataLikeQCDFile = ROOT.TFile(args.dataLikeQCDFileName, 'READ')
    if not dataLikeQCDFile:
        raise SystemExit('\n***ERROR*** couldn\'t open data-like QCD input file: %s'%args.dataLikeQCDFileName)

    #QCDHist = dataLikeQCDFile.Get(histName+'_DL')
    QCDHist = dataLikeQCDFile.Get(histNameBkg)
    if not QCDHist:
        raise SystemExit('\n***ERROR*** couldn\'t find data-like QCD input histogram: %s'%histNameBkg)

    QCDHist.SetName(histNameBkg+'_QCD')
    print ("histNameBkg: ",histNameBkg)
    QCDHist.SetTitle(histNameBkg+'_QCD')
    outFile.cd()
    QCDHist.Write()

    #------------------------------------------
    #define canvas before entering loop over input signals
    if args.plot:
        c1 = ROOT.TCanvas('c1', 'c1', 100, 50, 800, 600)
        c1.SetLogx(1)
        c1.SetLogy(1)

    #------------------------------------------
    #get signal samples
    #fileList = os.listdir(localdir+'/../inputs/'+model+'/')
    fileList = os.listdir(localdir+'/../inputs2/'+model+'/')
    print("sigFiles: ", fileList)
    for sigFileName in sorted(fileList):

        #------------------------------------------
        #check that it is a valid signal file
        if not '.root' in sigFileName:
            continue

        #------------------------------------------
        #get signal mass value from file name
        mass = sensitivityTools.getSignalMass(sigFileName)
        print '\n%s: %s GeV'%(model,mass)

        #check that the signal mass value is contained in the input list of signal masses (config file)
        #if float(mass) not in massValuesConfig:
        #    print '  skip'
        #    continue

        #------------------------------------------
        #TEST
        #if float(mass) != 3000.:
        #    continue

        #------------------------------------------
        #get signal hist
        sigFile = ROOT.TFile(localdir+'/../inputs2/'+model+'/'+sigFileName)
        print("path: ", localdir+'/../inputs2/'+model+'/'+sigFileName)
        print("histNameSig: ", histNameSig)

        histNameSigUpdated=histNameSig.format(str(mass))

        print("histNameSigUpdated: ", histNameSigUpdated)

        sigHist = sigFile.Get(histNameSigUpdated)

        print("sigHistIntegral: ", sigHist.Integral())

        print(sigHist, sigHist)

        print("histNameSig: ", histNameSigUpdated)

        sigHist.SetName(histNameSigUpdated+'_sig')

        sigHist.SetTitle(histNameSigUpdated+'_sig')

        if sigHist is not None and args.debug:
            print '\n  signal events:           ', sigHist.GetEntries()
            print '  effective signal events: ', sigHist.GetEffectiveEntries()
            print '  sum of signal weights:   ', sigHist.GetSumOfWeights()

        #------------------------------------------
        #scaled signal hist
        scaleFactor = float(args.lumi)*1000
        #NOTE cross sections are included already
        '''
        if args.fixCS:
            ZprimebbCrossSection = sensitivityTools.getZprimebbCrossSection(mass)
            print '\n  Z\'->bb cross section applied: %s fb'%ZprimebbCrossSection
            scaleFactor*=ZprimebbCrossSection
        '''
        if args.debug:
            print '\n  scale factor: %s'%scaleFactor

        #scale histogram by factor
        scaledSigHist = sigHist.Clone()
        scaledSigHist.Scale(scaleFactor)
        scaledSigHist.SetName(histNameSig+'_scaledSig')

        if scaledSigHist is not None and args.debug:
            print '\n  events after scaling:            ', scaledSigHist.GetEntries()
            print '  effective entries after scaling: ', scaledSigHist.GetEffectiveEntries()
            print '  sum of weights after scaling:    ', scaledSigHist.GetSumOfWeights()
            print '  first bin above 0 after scaline: ', scaledSigHist.FindFirstBinAbove(0,1)

        #------------------------------------------
        #effective entries signal hist
        effEntSigHist = plotTools.getEffectiveEntriesHistogram(sigHist, histNameSig+'_eff')
        if args.debug:
            print '\n  entries in effective entries:           ', effEntSigHist.GetEntries()
            print '  effective entries in effective entries: ', effEntSigHist.GetEffectiveEntries()
            print '  sum of weights in effective entries:    ', effEntSigHist.GetSumOfWeights()

        #------------------------------------------
        #data-like signal hist
        histNameSig2=histNameSig.format(str(int(mass)))
        dataLikeSigHist = plotTools.getDataLikeHist(effEntSigHist, scaledSigHist, histNameSig2+'_sig', seed, thresholdMass)
        firstBin = -1
        for ii in xrange(dataLikeSigHist.GetNbinsX()):
            if dataLikeSigHist.GetBinContent(ii)>0.:
                firstBin=ii
                break

        if firstBin > 0.:
            print '  data-like signal histogram starts at %s GeV'%dataLikeSigHist.GetBinLowEdge(firstBin)
        else:
            print '  data-like signal histogram is empty: n. entries = %s'%dataLikeSigHist.GetEntries()

        #------------------------------------------
        #smooth signal hist
        histNameSigSmooth=histNameSig.format(str(int(mass)))
        smoothSigHist = plotTools.getSmoothHistogram(scaledSigHist, histNameSigSmooth+'_sig_smooth')

        #------------------------------------------
        #remove any signal entry if there are no QCD entries (low mass region)
        if dataLikeSigHist is not None:
            for ii in xrange(QCDHist.GetNbinsX()):
                if QCDHist.GetBinContent(ii) == 0:
                    dataLikeSigHist.SetBinContent(ii,0)
                    dataLikeSigHist.SetBinError(ii,0)
                    smoothSigHist.SetBinContent(ii,0)
                    smoothSigHist.SetBinError(ii,0)
                else:
                    break

        #------------------------------------------
        #signal+QCD hist
        if QCDHist is not None and args.debug:
            print '\n  QCD events:           ', QCDHist.GetEntries()
            print '  QCD effective events: ', QCDHist.GetEffectiveEntries()
            print '  QCD sum of weights:   ', QCDHist.GetSumOfWeights()

        totHist = QCDHist.Clone()

        #..........................................
        #ORIGINAL
        #if dataLikeSigHist is not None:
        #    totHist.Add(dataLikeSigHist)

        #..........................................
        #NEW
        if dataLikeSigHist.GetBinContent(ii)>0.:
            print '  injecting data-like signal'
            totHist.Add(dataLikeSigHist)
        else:
            print '  not enough effective entries: injecting smoooth signal'
            totHist.Add(smoothSigHist)
        #END NEW
        #..........................................
        histNameToHist=histNameSig.format(str( int(mass)))

        histNameToHist=histNameToHist+"injectedToBkg"
        totHist.SetName(histNameToHist)
        totHist.SetTitle(histNameToHist)

        if totHist is not None and args.debug:
            print '\n  tot events:           ', totHist.GetEntries()
            print '  tot effective events: ', totHist.GetEffectiveEntries()
            print '  tot sum of weights:   ', totHist.GetSumOfWeights()

        #------------------------------------------
        #write to output file
        outFile.cd()
        dataLikeSigHist.Write()
        totHist.Write()

        #------------------------------------------
        #plot
        if args.plot:

            QCDHist.SetMarkerStyle(24)

            totHist.SetMarkerStyle(20)

            effEntSigHist.SetMarkerStyle(20)
            effEntSigHist.SetMarkerColor(ROOT.kGreen+1)
            effEntSigHist.SetLineColor(ROOT.kGreen+1)

            sigHist.SetMarkerStyle(25)
            sigHist.SetMarkerColor(ROOT.kAzure+1)
            sigHist.SetLineColor(ROOT.kAzure+1)

            scaledSigHist.SetMarkerStyle(26)
            scaledSigHist.SetMarkerColor(ROOT.kOrange+1)
            scaledSigHist.SetLineColor(ROOT.kOrange+1)

            dataLikeSigHist.SetMarkerStyle(24)
            dataLikeSigHist.SetMarkerColor(ROOT.kRed+1)
            dataLikeSigHist.SetLineColor(ROOT.kRed+1)

            smoothSigHist.SetMarkerStyle(32)
            smoothSigHist.SetMarkerColor(ROOT.kMagenta+1)
            smoothSigHist.SetLineColor(ROOT.kMagenta+1)

            hs = ROOT.THStack('hs','hs')
            hs.Add(QCDHist)
            hs.Add(totHist)
            hs.Add(effEntSigHist)
            hs.Add(sigHist)
            hs.Add(scaledSigHist)
            hs.Add(dataLikeSigHist)
            hs.Add(smoothSigHist)
            hs.Draw('nostack')

            hs.GetXaxis().SetTitle(axisLabel)
            hs.GetXaxis().SetTitleFont(43)
            hs.GetXaxis().SetTitleSize(textSize)
            #hs.GetXaxis().SetTitleOffset(1.5)
            hs.GetXaxis().SetLabelFont(43)
            hs.GetXaxis().SetLabelSize(textSize)

            hs.GetYaxis().SetTitle('entries')
            hs.GetYaxis().SetTitleFont(43)
            hs.GetYaxis().SetTitleSize(textSize)
            #hs.GetYaxis().SetTitleOffset(1.5)
            hs.GetYaxis().SetLabelFont(43)
            hs.GetYaxis().SetLabelSize(textSize)

            #------------------------------------------
            #labels and legends
            ax = 0.65
            ay = 0.88 #0.85
            a = plotTools.getATLAS()
            p = plotTools.getInternal()
            n = plotTools.getNote(textSize)
            l = plotTools.getLegend(ax,ay,textSize)

            #ATLAS internal
            a.DrawLatex(ax,ay,'ATLAS')
            p.DrawLatex(ax+0.13,ay,'internal')

            #notes
            notes=[]
            notes.append('#sqrt{s} = 13 TeV')
            if float(args.lumi) < -0.1:
                notes.append('L_{int} = %.0f pb^{-1}'%float(args.lumi)*1e3)
            else:
                notes.append('L_{int} = %.1f fb^{-1}'%float(args.lumi))
            if modelLabel != '':
                notes.append('m_{%s} = %0.f GeV'%(modelLabel,float(mass)))
            else:
                notes.append('m_{%s} = %0.f GeV'%(model,float(mass)))
            #notes.append('%s par. fit func.'%nPar)
            notes+=configNotes

            for ii, note in enumerate(notes):
                n.DrawLatex(ax,ay-0.04*(ii+1),note)

            #legend
            l.Clear()
            l.SetTextSize(textSize)
            l.AddEntry(QCDHist,"QCD background","pl")
            l.AddEntry(totHist,"signal plus backgournd","pl")
            l.AddEntry(sigHist,"non-scaled signal","pl")
            l.AddEntry(scaledSigHist,"scaled signal","pl")
            l.AddEntry(effEntSigHist,"effective signal","pl")
            l.AddEntry(dataLikeSigHist,"data-like signal","pl")
            l.AddEntry(smoothSigHist,"smooth signal","pl")
            l.SetX1(ax)
            l.SetY1(ay - 0.04*(len(notes)+1) - 0.04*l.GetNRows())
            l.SetX2(ax+0.15)
            l.SetY2(ay - 0.04*(len(notes)+1))
            l.Draw()

            c1.Update()
            if args.wait:
                c1.WaitPrimitive()
            c1.SaveAs('../figures/signalplusbackground.'+model+'.'+mass+'.GeV.'+slumi+'.ifb.'+histNameSig+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.pdf')

    outFile.Write()
    outFile.Close()

#******************************************
if __name__ == '__main__':

    #------------------------------------------
    #parse input arguments
    parser = argparse.ArgumentParser(description='%prog [options]')
    parser.add_argument('--config', dest='configFileName', default='', required=True, help='sensitivity scan config file')
    parser.add_argument('--QCDFile', dest='dataLikeQCDFileName', default='', required=True, help='input data-like QCD file')
    parser.add_argument('--lumi', dest='lumi', default='0.', required=True, help='luminosity [fb^-1]')
    parser.add_argument('--tag', dest='tag', default='default', help='tag for output files')
    #parser.add_argument('--fixCS', dest='fixCS', action='store_true', default=False, help='fix Z\'->bb cross section')
    parser.add_argument('--wait', dest='wait', action='store_true', default=False, help='wait?')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true', default=False, help='plot histograms')
    parser.add_argument('-b', '--batch', dest='batch', action='store_true', default=False, help='batch mode')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    args = parser.parse_args()

    #------------------------------------------
    #inject data-like signal
    injectDataLikeSignal(args)
    print "done with step2"
    print '\ndata-like signal injected'
    print '******************************************'
