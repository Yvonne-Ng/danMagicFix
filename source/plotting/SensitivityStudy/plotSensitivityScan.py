#!/bin/python

#******************************************
#plot sensitivity scan by looping over all masses on a range of luminosity values
#EXAMPLE python -u plotSensitivityScan.py ./results/searchphase/ QStar 100. 3000e3 3

#******************************************
#import stuff
import ROOT
import re, sys, os, math, shutil, fileinput
import tools
import numpy as np

#******************************************
#set ATLAS style
from art.morisot import Morisot
ROOT.SetAtlasStyle()

#******************************************
#main
#******************************************
if __name__ == "__main__":

    #------------------------------------------
    #check input parameters
    if len(sys.argv) != 6:
        raise SystemExit(
            '\n***ERROR*** wrong input parameters (%s/%s) \
            \nNOTE: plot sensitivity scan by looping over all masses on a range of luminosity values \
            \nHOW TO: python doSensitivityScan.py model luminosity-min[pb^-1] luminosity-max[pb^-1] par' \
            %(len(sys.argv),6) )

    #------------------------------------------
    #get input parameters
    path = sys.argv[1].strip()
    model = sys.argv[2].strip()
    lumimin = float(sys.argv[3].strip())
    if lumimin <=100. : #pb^-1
        lumimin =100. #pb^-1
    lumimax = float(sys.argv[4].strip())
    par = sys.argv[5].strip()

    #------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #check number of parameters for fit
    if par is str(3) or par is str(4) or par is str(5):
        parstag = par+'par'
    else:
        raise SystemExit('\n***ERROR*** Number of parameters used in search phase fit must 3, 4 or 5')

    #------------------------------------------
    #get list of available mass points
    massValues = []
    fileList=os.listdir(path) 
    for sigFileName in sorted(fileList):
        if (model not in sigFileName or '.root' not in sigFileName):
            continue
        mymass=sigFileName.split('QStar')[1].split('_')[0]
        if mymass not in massValues: # dia ok to do this?
          massValues.append(mymass)
    massValues.sort(key=int)
    #massValues.remove("5000") # Lydia remove these when issue fixed FIXME
    #massValues.remove("5500")
    #massValues.remove("6000")

    # FIXME do like before: massValues.append( float( tools.getSignalMass(sigFileName) ) )
    #massValues = [m for m in massValues if m <= 7000.0] #TEST
    #massValues = [4000.0, 5000.0, 6000.0] #TEST
    print '\navailable mass values [GeV]: %s'%massValues
    
    print "LYDIA ORGANISE OUTPUT DIRS"
    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    #------------------------------------------
    #initial luminosity value
    lumi = lumimin

    #------------------------------------------
    #arrays for sensitivity scan graphs
    gmass = np.array(massValues)
    glumi = np.zeros_like(gmass)
    
    #------------------------------------------
    #loop over mass values and luminosity
    while len(massValues) > 0 and lumi <= lumimax:

        #------------------------------------------
        print '\n\n******************************************'
        print '******************************************'
        print 'luminosity = %s ^pb-1'%lumi
        print 'available mass values [GeV]: %s'%massValues
        print '******************************************'
        print '******************************************\n'

        #------------------------------------------
        #loop over signal mass values
        removeMassValues = []
        for mass in massValues:

            #------------------------------------------
            #signal mass
            print '\n******************************************\nSearchPhase\n'
            print '%s mass = %s GeV'%(model, mass)
            print 'lumi = %s pb^-1'%lumi
            print '******************************************'

            #------------------------------------------
            #check SearchPhase results
            print '******************************************\nSearchPhase results\n'
            print '%s mass = %s GeV'%(model, mass)
            print 'lumi = %s pb^-1'%lumi

            if float(lumi) <1000: # i.e. less than 1fb 
              filenamelumi = float(lumi)/1000
              filenamelumi=str(round(filenamelumi,1))
              filenamelumi = filenamelumi.replace(".","-")
            else:
              filenamelumi = float(lumi)/1000
              filenamelumi = str(int(filenamelumi))
            print "FNL!!!!!!!!!!!!!"+filenamelumi
            spFileName = path+'Step1_SearchPhase_DiJet_'+model+mass+'_'+filenamelumi+'fb.root' #ORIGINAL
            # FIXME do like before: spFileName = localdir+'/searchphase.'+model+'.%.0f.GeV.%.0f.ipb.'%(mass, lumi)+parstag+'.root' # FIXME so organised '../results/searchphase/searchphase.'+model+'.%.0f.GeV.%.0f.ipb.'%(mass, lumi)+parstag+'.root' #ORIGINAL
            #spFileName = '../results/searchphase/backups/20141124.new.injection.AixLesBains/searchphase.'+model+'.%.0f.GeV.%.0f.ipb.'%(mass, lumi)+parstag+'.root' #TEST 1
            #spFileName = '../results/searchphase/backups/20141112.4.and.5.par.scan/searchphase.'+model+'.%.0f.GeV.%.0f.ipb.'%(mass, lumi)+parstag+'.root' #TEST 2
            #spFileName = '../results/searchphase/backups/20141203.double.fluctuation/searchphase.'+model+'.%.0f.GeV.%.0f.ipb.'%(mass, lumi)+parstag+'.root' #TEST 3
            #spFileName = '../results/searchphase/backups/20150108.note.with.additional.points.20150102/searchphase.'+model+'.%.0f.GeV.%.0f.ipb.'%(mass, lumi)+parstag+'.root' #TEST 4
            #spFileName = '../results/searchphase/backups/20150127.100.ipb.seed.0/searchphase.'+model+'.%.0f.GeV.%.0f.ipb.'%(mass, lumi)+parstag+'.seed.%i'%jobSeed+'.root' #TEST 5
            if os.path.isfile(spFileName):
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

                    #TEST remove by hand fake discovery points
                    #if float(mass) >= 8000.:
                    #    glumi[massIndex] = 0.
                    #else:
                    glumi[massIndex] = lumi #NOTE always keep this uncommented
                    
                    removeMassValues.append(mass)

                #------------------------------------------

            else:
                print 'BumpHunter results not available for this %s mass value'%model

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
        # current lumi binning
        if lumi < 500.:
            lumi = 500.
        elif lumi < 1.e3:
            lumi = 1.e3
        elif lumi < 3.e3:
            lumi = 3.e3
        elif lumi < 5.e3:
            lumi = 5.e3
        elif lumi < 7.e3:
            lumi = 7.e3
        elif lumi < 10.e3:
            lumi = 10.e3
        elif lumi < 15.e3:
            lumi = 15.e3
        elif lumi < 30.e3:
            lumi = 30.e3

        #tight lumi binning
        ####if lumi < 200.:
        ####    lumi = 200.
        ####elif lumi < 300.:
        ####    lumi = 300.
        ####elif lumi < 500.:
        ####    lumi = 500.
        ####elif lumi < 1.e3:
        ####    lumi = 1.e3
        ####elif lumi < 2.e3:
        ####    lumi = 2.e3
        ####elif lumi < 3.e3:
        ####    lumi = 3.e3
        ####elif lumi < 5.e3:
        ####    lumi = 5.e3
        ####elif lumi < 10.e3:
        ####    lumi = 10.e3
        ####elif lumi < 15.e3:
        ####    lumi = 15.e3
        ####elif lumi < 20.e3:
        ####    lumi = 20.e3
        ####elif lumi < 25.e3:
        ####    lumi = 25.e3
        ####elif lumi < 50.e3:
        ####    lumi = 50.e3
        #NOTE intermediate steps
        #elif lumi < 60.e3:
        #    lumi = 60.e3
        #elif lumi < 70.e3:
        #    lumi = 70.e3
        #elif lumi < 80.e3:
        #    lumi = 80.e3
        #elif lumi < 90.e3:
        #    lumi = 90.e3
        #END intermediate steps
        ####elif lumi < 100.e3:
        ####    lumi = 100.e3
        ####elif lumi < 200.e3:
        ####    lumi = 200.e3
        ####elif lumi < 300.e3:
        ####    lumi = 300.e3
        ####elif lumi < 500.e3:
        ####    lumi = 500.e3
        ####elif lumi < 1000.e3:
        ####    lumi = 1000.e3
        ####elif lumi < 2000.e3:
        ####    lumi = 2000.e3
        ####elif lumi < 3000.e3:
        ####    lumi = 3000.e3
        ####elif lumi < 10000.e3:
        ####    lumi = 10000.e3
        ####else:
        ####    lumi += 10000.e3
        
        #coarse lumi binning
        '''
        if lumi < 5.e3:
            lumi = 5.e3
        elif lumi < 10.e3:
            lumi = 10.e3
        elif lumi < 25.e3:
            lumi = 25.e3
        elif lumi < 100.e3:
            lumi = 100.e3
        elif lumi < 300.e3:
            lumi = 300.e3
        elif lumi < 1000.e3:
            lumi = 1000.e3
        elif lumi < 3000.e3:
            lumi = 3000.e3
        elif lumi < 10000.e3:
            lumi = 10000.e3
        else:
            lumi += 10000.e3
        '''
        
    #------------------------------------------
    #print sensitivity scan results
    print '\n******************************************'
    print 'sensitivity scan results'
    print '******************************************\n'
    print 'mass = %s'%gmass
    print 'lumi = %s'%glumi

    #------------------------------------------
    #TEST
    #raise SystemExit('\n***TEST*** exit')
    #------------------------------------------

    #canvas
    c1 = ROOT.TCanvas('c1', 'c1', 400, 50, 800, 600)
    c1.SetLogy(1)

    #model name
    if model == 'QStar':
        modelName = 'q*'
    else:
        modelName = model

    #graph
    g = ROOT.TGraph(len(gmass), gmass, glumi)
    g.SetTitle('sensitivity studies')
    g.GetXaxis().SetTitle('m_{'+modelName+'} [GeV]')
    g.GetYaxis().SetTitle('discovery luminosity [pb^{-1}]')
    g.SetMinimum(90) #HERE before Draw() Lydia OK Changed??
    g.SetMaximum(31e3)#HERE before Draw()
    g.Draw("ap")

    #lines
    l1000 = ROOT.TLine(g.GetXaxis().GetXmin(), 1.e3, g.GetXaxis().GetXmax(), 1.e3)
    l1000.SetLineColor(ROOT.kRed)
    l1000.SetLineWidth(2)
    l1000.Draw()
    
    l25000 = ROOT.TLine(g.GetXaxis().GetXmin(), 25.e3, g.GetXaxis().GetXmax(), 25.e3)
    l25000.SetLineColor(ROOT.kRed)
    l25000.SetLineWidth(2)
    l25000.Draw()

    l300000 = ROOT.TLine(g.GetXaxis().GetXmin(), 300.e3, g.GetXaxis().GetXmax(), 300.e3)
    l300000.SetLineColor(ROOT.kRed)
    l300000.SetLineWidth(2)
    l300000.Draw()

    l3000000 = ROOT.TLine(g.GetXaxis().GetXmin(), 3000.e3, g.GetXaxis().GetXmax(), 3000.e3)
    l3000000.SetLineColor(ROOT.kRed)
    l3000000.SetLineWidth(2)
    l3000000.Draw()

    #draw graph (again but) on top of the lines
    g.Draw("same p")
    
    #position
    lXmin=0.70; lXmax=0.90; lYmin=0.61; lYmax=0.74;
    
    #ATLAS
    a = ROOT.TLatex()
    a.SetNDC()
    a.SetTextFont(72)
    a.SetTextColor(1)
    a.SetTextSize(0.04)
    a.DrawLatex(lXmin,lYmax+0.12,'ATLAS')
        
    #internal
    p = ROOT.TLatex()
    p.SetNDC()
    p.SetTextFont(42)
    p.SetTextColor(1)
    p.SetTextSize(0.04)
    p.DrawLatex(lXmin+0.1,lYmax+0.12,'internal')

    #notes
    n = ROOT.TLatex()
    n.SetNDC()
    n.SetTextFont(42)
    n.SetTextColor(1)
    n.SetTextSize(0.03)
    n.DrawLatex(lXmin,lYmax+0.08,'#sqrt{s} = 13 TeV')
    n.DrawLatex(lXmin,lYmax+0.04,modelName+' sensitivity scan')
    if par is str(3):
        n.DrawLatex(lXmin,lYmax+0.0,'3 parameter fit')
    if par is str(4):
        n.DrawLatex(lXmin,lYmax+0.0,'4 parameter fit')
    if par is str(5):
        n.DrawLatex(lXmin,lYmax+0.0,'5 parameter fit')

    c1.Update()
    #c1.WaitPrimitive() Lydia
    c1.SaveAs(localdir+'sensitivityscan.'+model+'.'+str(int(round(lumimin)))+'.ipb.'+str(int(round(lumimax)))+'.ipb.'+parstag+'.pdf') # FIXME so organised'../figures/sensitivityscan.'+model+'.'+str(int(round(lumimin)))+'.ipb.'+str(int(round(lumimax)))+'.ipb.'+parstag+'.pdf') #ORIGINAL
