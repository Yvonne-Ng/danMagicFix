#!/bin/python

#******************************************
#script to get data-like distribution from MC
#EXAMPLE python -u step01.getDataLikeQCD.py --config <config file> --lumi <luminosity [fb^-1]> --tag <tag> --patch --wait --plot --batch --debug

#******************************************
#import stuff
import sys, os, math, argparse, ROOT
import plotTools

#******************************************
def getDataLikeQCD(args):

    print '\n******************************************'
    print 'get data-like QCD'

    #******************************************
    #set ATLAS style
    if os.path.isfile(os.path.expanduser('~/RootUtils/AtlasStyle.C')):
        ROOT.gROOT.LoadMacro('~/RootUtils/AtlasStyle.C')
        ROOT.SetAtlasStyle()
        ROOT.set_color_env()
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

    #------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #get settings
    print '\nconfig settings:'
    settings = ROOT.TEnv()
    if settings.ReadFile(args.configFileName,ROOT.EEnvLevel(0)) != 0:
        raise IOError('could not find sensitivity scan config file: %s'%args.configFileName)

    QCDFileName = settings.GetValue('QCDFile','../inputs/QCD/histograms.mc.dijet.1p0.ifb.root')
    print '  QCD input file = %s'%QCDFileName

    patchFileName = settings.GetValue('patchFile','')
    print '  patch file = %s'%patchFileName

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
    #check QCD file
    if not os.path.isfile(QCDFileName):
        raise SystemExit('\n***ERROR*** couldn\'t find QCD file: %s'%QCDFileName)

	#------------------------------------------
    #check output file
    outFileName = localdir+'/../results/datalikeQCD/datalikeQCD.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.root'
    if not (args.overwrite):
        if os.path.isfile(outFileName):
            raise SystemExit('\n***WARNING*** output file exists already: %s'%outFileName)

	#------------------------------------------
    #QCD hist
    QCDFile = ROOT.TFile(QCDFileName, 'READ')
    if not QCDFile:
        raise SystemExit('\n***ERROR*** couldn\'t open QCD input file: %s'%QCDFilename)

    QCDHist = QCDFile.Get(histName)
    if not QCDHist:
        raise SystemExit('\n***ERROR*** couldn\'t find QCD input histogram: %s'%histName)

    QCDHist.SetName(histName+'_input')

	#------------------------------------------
	#scaled QCD hist
    scaledHist = QCDHist.Clone()
    scaledHist.SetName(histName+'_scaled')
    scaledHist.SetTitle(histName+'_scaled')
    scaleFactor = float(args.lumi)
    if args.debug:
        print 'scale factor = %s'%scaleFactor
    scaledHist.Scale(scaleFactor)

    #------------------------------------------
    #effective entries QCD hist
    effEntHist = plotTools.getEffectiveEntriesHistogram(QCDHist, histName+'_eff')

    #------------------------------------------
    #data-like QCD hist
    pureDataLikeHist = plotTools.getDataLikeHist(effEntHist, scaledHist, histName+'_pureDL', seed, thresholdMass)
    firstBin = -1
    for ii in xrange(pureDataLikeHist.GetNbinsX()):
        if pureDataLikeHist.GetBinContent(ii)>0.:
            firstBin=ii
            break

    if firstBin > 0.:
        print '\ndata-like histogram starts at %s GeV'%pureDataLikeHist.GetBinLowEdge(firstBin)
    else:
        raise SystemExit('\n***ERROR*** data-like histogram has no entries')

    #------------------------------------------
    #smooth QCD hist
    smoothHist = plotTools.getSmoothHistogram(scaledHist, histName+'_smooth')

    #------------------------------------------
    #OPTION 1 PATCH: patch from smooth hist
    #NOTE this introduces features in the background, use option 2 instead
    '''
    #patch data-like hist with smooth hist + Poisson noise to extend low mass region when there are not enough effective entries
    dataLikeHist = pureDataLikeHist.Clone()
    dataLikeHist.SetName(histName+'_DL')
    dataLikeHist.SetTitle(histName+'_DL')

    patchHist = ROOT.TH1D(histName+'_patch', histName+'_patch', QCDHist.GetXaxis().GetNbins(), QCDHist.GetXaxis().GetXbins().GetArray())
    if args.patch:
        for ii in xrange( int( round(pureDataLikeHist.GetNbinsX()))):
            if ii > firstBin:
                break
            if pureDataLikeHist.GetXaxis().GetBinUpEdge(ii) > thresholdMass and dataLikeHist.GetBinContent(ii) == 0.:

                #bincontent
                bincontent = int( round(scaledHist.GetBinContent(ii)))

                #random number generator
                #NOTE the seed for each bin must be always the same
                binSeed = int( round( patchHist.GetBinCenter(ii) + seed*1e5))
                rand3 = ROOT.TRandom3(binSeed)

                #add Poisson noise
                bincontent = int( round( rand3.PoissonD(bincontent)))

                #fill
                for jj in xrange(bincontent):
                    patchHist.Fill(patchHist.GetBinCenter(ii))

        #patch data-like histogram
        dataLikeHist.Add(patchHist)
    '''
    #------------------------------------------
    #END OPTION 1 PATCH

    #------------------------------------------
    #OPTION 2 PATCH: patch from smooth fit
    #patch data-like hist with fit of smooth hist; then apply Poisson noise
    dataLikeHist = pureDataLikeHist.Clone()
    dataLikeHist.SetName(histName+'_DL')
    dataLikeHist.SetTitle(histName+'_DL')

    #..........................................
    #mjj or mbj or mbb? 85% or 77%? fixed or flat? 3 or 4 parameter fit?
    #this is setting the output file name TODO Yvonne: fix this
    if len(patchFileName) == 0:
        if 'mbb' in histName:
            if nPar == 3:
                if 'fix_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_fix_8585_smooth.0.seed.3.par.root'
                elif 'flt_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_flt_8585_smooth.0.seed.3.par.root'
                elif 'fix_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_fix_7777_smooth.0.seed.3.par.root'
                elif 'flt_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_flt_7777_smooth.0.seed.3.par.root'
                else:
                    args.patch = False
            elif nPar == 4:
                if 'fix_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_fix_8585_smooth.0.seed.4.par.root'
                elif 'flt_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_flt_8585_smooth.0.seed.4.par.root'
                elif 'fix_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_fix_7777_smooth.0.seed.4.par.root'
                elif 'flt_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbb_flt_7777_smooth.0.seed.4.par.root'
                else:
                    args.patch = False
            else:
                args.patch = False
        elif 'mbj' in histName:
            if nPar == 3:
                if 'fix_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_fix_8585_smooth.0.seed.3.par.root'
                elif 'flt_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_flt_8585_smooth.0.seed.3.par.root'
                elif 'fix_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_fix_7777_smooth.0.seed.3.par.root'
                elif 'flt_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_flt_7777_smooth.0.seed.3.par.root'
                else:
                    args.patch = False
            elif nPar == 4:
                if 'fix_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_fix_8585_smooth.0.seed.4.par.root'
                elif 'flt_8585' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_flt_8585_smooth.0.seed.4.par.root'
                elif 'fix_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_fix_7777_smooth.0.seed.4.par.root'
                elif 'flt_7777' in histName:
                    patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mbj_flt_7777_smooth.0.seed.4.par.root'
                else:
                    args.patch = False
            else:
                args.patch = False
        else:
            if nPar == 3:
                patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mjj_smooth.0.seed.3.par.root'
            elif nPar == 4:
                patchFileName = localdir+'/../data/searchPhase.patch.QCD.20p0.ifb.mjj_smooth.0.seed.4.par.root'
            else:
                args.patch = False

    #..........................................
    #open patch file
    if not os.path.isfile(patchFileName):
        raise SystemExit('\n***ERROR*** couldn\'t find patch file: %s'%patchFileName)
    patchFile = ROOT.TFile(patchFileName,'READ')
    if not patchFile:
        raise SystemExit('\n***ERROR*** couldn\'t open patch file: %s'%patchFileName)

    #..........................................
    #get patch hist
    fullPatchHist = patchFile.Get("basicBkgFrom4ParamFit")

    #check the binninng of the patch and the nominal hist are the same
    #if QCDHist.GetXaxis().GetXbins().GetArray() != fullPatchHist.GetXaxis().GetXbins().GetArray():
    if QCDHist.GetNbinsX() != fullPatchHist.GetNbinsX() or QCDHist.GetBinLowEdge(1) != fullPatchHist.GetBinLowEdge(1) or QCDHist.GetBinLowEdge( QCDHist.GetNbinsX() ) != fullPatchHist.GetBinLowEdge( fullPatchHist.GetNbinsX() ):
        raise SystemExit('\n***ERROR*** QCD hist and patch hist have different binning')

    #..........................................
    #scale patch to right luminosity
    #NOTE the patch are obtained from 20/fb fits
    patchScaleFactor = float(args.lumi)/20.
    print '\npatch scale factor = %s'%patchScaleFactor
    fullPatchHist.Scale(patchScaleFactor)

    #..........................................
    #apply patch
    patchHist = ROOT.TH1D(histName+'_patch', histName+'_patch', QCDHist.GetXaxis().GetNbins(), QCDHist.GetXaxis().GetXbins().GetArray())
    if args.patch:
        for ii in xrange( int( round(pureDataLikeHist.GetNbinsX()))):
            if ii > firstBin:
                break
            if pureDataLikeHist.GetXaxis().GetBinUpEdge(ii) > thresholdMass and dataLikeHist.GetBinContent(ii) == 0.:

                #bincontent
                bincontent = int( round(fullPatchHist.GetBinContent(ii)))

                #random number generator
                #NOTE the seed for each bin must be always the same
                binSeed = int( round( patchHist.GetBinCenter(ii) + seed*1e5))
                rand3 = ROOT.TRandom3(binSeed)

                #add Poisson noise
                bincontent = int( round( rand3.PoissonD(bincontent)))

                #fill
                for jj in xrange(bincontent):
                    patchHist.Fill(patchHist.GetBinCenter(ii))

        #patch data-like histogram
        dataLikeHist.Add(patchHist)
    else:
        print 'patch not applied'
    #END OPTION 2 PATCH
    #------------------------------------------

    #------------------------------------------
    #check entries and weights of the data-like QCD background
    if args.debug:
        print '\nevents:		   ', pureDataLikeHist.GetEntries()
        print 'effective events: ', pureDataLikeHist.GetEffectiveEntries()
        print 'sum of weights:   ', pureDataLikeHist.GetSumOfWeights()

    #------------------------------------------
    #write data-like hist to output file
    outFile = ROOT.TFile(outFileName, 'RECREATE')
    outFile.cd()
    scaledHist.Write()
    effEntHist.Write()
    pureDataLikeHist.Write()
    smoothHist.Write()
    patchHist.Write()
    dataLikeHist.Write()
    outFile.Write()

	#------------------------------------------
	#plot
    if args.plot:
        c1 = ROOT.TCanvas('c1', 'c1', 100, 50, 800, 600)
        c1.SetLogx(not args.linx)
        c1.SetLogy(1)

        effEntHist.SetMarkerColor(ROOT.kGreen+1)
        effEntHist.SetMarkerStyle(20)
        effEntHist.SetLineColor(ROOT.kGreen+1)

        QCDHist.SetMarkerColor(ROOT.kAzure+1)
        QCDHist.SetMarkerStyle(25)
        QCDHist.SetLineColor(ROOT.kAzure+1)

        scaledHist.SetMarkerColor(ROOT.kOrange+1)
        scaledHist.SetMarkerStyle(26)
        scaledHist.SetLineColor(ROOT.kOrange+1)

        patchHist.SetMarkerColor(ROOT.kMagenta+1)
        patchHist.SetMarkerStyle(32)
        patchHist.SetLineColor(ROOT.kMagenta+1)

        pureDataLikeHist.SetMarkerColor(ROOT.kRed+1)
        pureDataLikeHist.SetMarkerStyle(24)
        pureDataLikeHist.SetLineColor(ROOT.kRed+1)

        hs = ROOT.THStack('hs','hs')
        hs.Add(effEntHist)
        hs.Add(QCDHist)
        hs.Add(scaledHist)
        hs.Add(patchHist)
        hs.Add(pureDataLikeHist)
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

        #use new axis in TeV instead of GeV
        #NOTE this is over-complicated and will be used only if requested
        #NOTE when the x-axis is logarithmic, GetUxm* will return the exponent only
        '''
        #remove the current axis
        hs.GetXaxis().SetLabelOffset(999)
        hs.GetXaxis().SetTickLength(0)

        #map the axis values to -x
        func = ROOT.TF1('func', 'x/1000.', ROOT.gPad.GetUxmin()/1000., ROOT.gPad.GetUxmax()/1000.)
        newaxis = ROOT.TGaxis(  ROOT.gPad.GetUxmin(),
                                ROOT.gPad.GetUymin(),
                                ROOT.gPad.GetUxmax(),
                                ROOT.gPad.GetUymin(),
                                'func',510,'+')

        newaxis.SetLabelOffset(0.005)
        newaxis.SetLabelFont(43)
        newaxis.SetLabelSize(textSize)
        newaxis.Draw()
        '''

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
        if float(args.lumi) < 0.1:
            notes.append('L_{int} = %.0f pb^{-1}'%float(args.lumi)*1e3)
        else:
            notes.append('L_{int} = %.1f fb^{-1}'%float(args.lumi))
        #notes.append('%s par. fit func.'%nPar)
        notes+=configNotes

        for ii, note in enumerate(notes):
            n.DrawLatex(ax,ay-0.04*(ii+1),note)

		#legend
        l.Clear()
        l.SetTextSize(textSize)
        l.AddEntry(QCDHist,'non-scaled','pl')
        l.AddEntry(scaledHist,'scaled','pl')
        l.AddEntry(effEntHist,'effective','pl')
        l.AddEntry(pureDataLikeHist,'data-like','pl')
        l.AddEntry(patchHist,'patch','pl')
        l.SetX1(ax)
        l.SetY1(ay - 0.04*(len(notes)+1) - 0.04*l.GetNRows())
        l.SetX2(ax+0.15)
        l.SetY2(ay - 0.04*(len(notes)+1))
        l.Draw()

        c1.Update()
        if args.wait:
            c1.WaitPrimitive()
        c1.SaveAs('../figures/datalikeQCD.'+slumi+'.ifb.'+histName+'.%i'%nPar+'.par.%i'%seed+'.seed.'+args.tag+'.pdf')

#******************************************
if __name__ == '__main__':

    #------------------------------------------
    #parse input arguments
    parser = argparse.ArgumentParser(description='%prog [options]')
    parser.add_argument('--config', dest='configFileName', default='', required=True, help='sensitivity scan config file')
    parser.add_argument('--lumi', dest='lumi', default='0.', required=True, help='luminosity [fb^-1]')
    parser.add_argument('--tag', dest='tag', default='default', help='tag for output files')
    parser.add_argument('--patch', dest='patch', action='store_true', default=False, help='patch data-like hist with smooth hist to extend low mass region?')
    parser.add_argument('--wait', dest='wait', action='store_true', default=False, help='wait?')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true', default=False, help='plot histograms')
    parser.add_argument('--linx', dest='linx', action='store_true', default=False, help='x axis linear scale?')
    parser.add_argument('-b', '--batch', dest='batch', action='store_true', default=False, help='batch mode')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    parser.add_argument('-ow', '--overwrite', dest='overwrite', action='store_true', default=False, help='overwrite output file')
    args = parser.parse_args()

    #------------------------------------------
    #get data-like QCD
    getDataLikeQCD(args)
    print '\ngot data-like QCD'
    print '******************************************'
