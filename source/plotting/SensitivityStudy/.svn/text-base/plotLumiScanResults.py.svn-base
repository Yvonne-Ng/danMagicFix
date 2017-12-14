#!/bin/python
# FIXME remove!! Step1_SearchPhase_QCDDiJetWithQStar2000_0p1fb.root
# FIXME FIX MASS LABELLING ON PLOTS FOR E.G. 4.5 writes 4????

#******************************************
#plot luminosity scan results for each mass point separately
# EXAMPLE python -u plotting/SensitivityStudy/plotLumiScanResults_MC15Inputs.py ./results/Step1_SearchPhase/Test/ QStar 3 Test
# FIXME remove!! Step1_SearchPhase_QCDDiJetWithQStar2000_0p1fb.root
# FIXME tell off if e.g. no / in argument paths?? check?
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
             \nNOTE: plot luminosity scan results for each mass point separately \
             \nHOW TO: python -u plotting/SensitivityStudy/plotLumiScanResults_MC15Inputs.py path/ model par folderextension \
             \npath/ is path to directory containing all the search phase results to be included in the sesitivity study \
             \n model is the model you want to perform the sensitivity study for e.g. QStar \
             \n par is the number of parameters used in the fit for the search phase e.g. 3 \
             \n folderextension is the name of the folder you want to store your outputs in \
             \n usingScaled set to 1 if using scaled data, or 0 if otherwise (i.e. if using data-like statistics)' \
            %(len(sys.argv),6) )
            
    #------------------------------------------
    #get input parameters
    path = sys.argv[1].strip()
    model = sys.argv[2].strip()
    par = sys.argv[3].strip()
    folderextension = sys.argv[4].strip() # Chosen by user to put plots in common folder 
    useScaled = sys.argv[5].strip() 

    #------------------------------------------
    #check input file
    if not os.path.isdir(path):
        raise SystemExit('\n***ERROR*** not a valid path')

    #------------------------------------------
    #check number of parameters for fit
    if par is str(3) or par is str(4) or par is str(5):
        parstag = par+'par'
    else:
        raise SystemExit('\n***ERROR*** Number of parameters used in search phase fit must 3, 4 or 5')
    
    #------------------------------------------
    # make plots folder i.e. make folder extension
    if not os.path.exists("plotting/SensitivityStudy/plots/%s"%folderextension):
        os.makedirs("plotting/SensitivityStudy/plots/%s"%folderextension)


    #------------------------------------------
    #get directory of this script
    localdir = os.path.dirname(os.path.realpath(__file__))

    #------------------------------------------
    #get list of mass points
    massValues = []
    searchphaseFileList=os.listdir(path)
    for searchphaseFileName in sorted(searchphaseFileList):
        if not '%s'%model in searchphaseFileName:
            continue
        massValues.append(searchphaseFileName.split(model)[1].split('_')[0])
    massValuesSet = set(massValues)
    massValues = list(massValuesSet)
    massValues.sort(key=int)
    print '\nmass values [GeV]: %s'%massValues


    #------------------------------------------
    # get list of luminosity values
    # and in addition store dictionary of lumi:filename lumi
    # i.e. filename lumi in inverse fb (femtobarns) with less than 1 inverse fb in the form 0pXfb e.g. 10fb, 0p5fb
    lumi = 0
    lumiValues = []
    filenamelumiDict = {}
    searchphaseFileList=os.listdir(path)
    for searchphaseFileName in sorted(searchphaseFileList):
        if not 'fb' in searchphaseFileName:
            continue
        # Use regex to find lumi of hist
        print searchphaseFileName #[0-9]+
        lumi = 0
        if (re.search('_[0-9]+fb',searchphaseFileName) is not None):
          lumi = re.search('_[0-9]+fb',searchphaseFileName).group()
          print "LUM"
          print lumi
          lumi = lumi.strip("_")
          filenamelumi = lumi
          lumi = lumi.strip("fb") 
          lumi = float(lumi)*1000
          print lumi
          lumi = str(int(lumi))
          filenamelumiDict[lumi] = filenamelumi

        if (re.search('_[0-9]+p[0-9]+fb',searchphaseFileName) is not None):
          lumi = re.search('_[0-9]+p[0-9]+fb',searchphaseFileName).group()
          print "LUM"
          print lumi
          lumi = lumi.strip("_")
          filenamelumi = lumi
          lumi = lumi.replace("p",".")
          lumi = lumi.strip("fb") 
          lumi = float(lumi)*1000
          print lumi
          lumi = str(int(lumi))
          filenamelumiDict[lumi] = filenamelumi
#        if (re.search('_\dp\dfb',searchphaseFileName) is not None):
 #           lumi = re.search('_[0-9]+p[0-9]+fb',searchphaseFileName).group()
 #           print "LUM"
 #           print lumi
 #           lumi = lumi.strip("_")
 #           filenamelumi = lumi
 #           lumi = lumi.replace("p",".")
 #           lumi = lumi.strip("fb") 
 #           lumi = float(lumi)*1000
 #           print lumi
 #           lumi = str(int(lumi))
 #           filenamelumiDict[lumi] = filenamelumi
        if lumi == 0: raise SystemExit('\n***Zero lumi*** regex issue')
        lumiValues.append(lumi)
    lumiValuesSet = set(lumiValues)
    lumiValues = list(lumiValuesSet)
    lumiValues.sort(key=int)
    print 'luminosity values [pb^-1]: %s'%lumiValues

    #------------------------------------------
    #loop over mass points
    for mass in massValues:

        #------------------------------------------
        #TEST just one mass value
        #if float(mass) < 12000:
        #    continue
        #------------------------------------------

        print '\n******************************************'
        print '%s m = %s GeV'%(model, mass)
        print '******************************************'
        glumi = []
        gpval = []
        gpvalerr = []
        # Lydia
        gpvalsigmas = []
        gpvalerrsigmas = []
        #
        glowedge = []
        ghighedge = []

        #------------------------------------------
        #p-values histogram #TEST
        #hpval = ROOT.TH1D('pvalues','pvalues',20,0.,1.)
        #hpval.SetDefaultSumw2(ROOT.kTRUE)
        #hpval.SetDirectory(0)
        
        #------------------------------------------
        #loop over luminosity values
        highestlumi = 0
        for lumi in lumiValues:
            print '\nluminosity = '+lumi+' pb^-1'

            discoverlumi = ''
            
            #get SearchPhase results
            if (useScaled):
              spFileName = path+'Step1_SearchPhase_mjj_Scaled_QCDDiJetWith'+model+mass+'_'+filenamelumiDict[lumi]+'.root' #ORIGINAL
            else:
              spFileName = path+'Step1_SearchPhase_mjj_DataLike_QCDDiJetWith'+model+mass+'_'+filenamelumiDict[lumi]+'.root' #ORIGINAL
            print spFileName
            if os.path.isfile(spFileName):
                spFile = ROOT.TFile(spFileName,'READ')

                bumpHunterStatOfFitToData = spFile.Get("bumpHunterStatOfFitToData")
                bumpHunterStatValue = bumpHunterStatOfFitToData[0]
                bumpHunterPValue    = bumpHunterStatOfFitToData[1]
                bumpHunterPValueErr = bumpHunterStatOfFitToData[2]

                bumpHunterPLowHigh = spFile.Get('bumpHunterPLowHigh')
                bumpLowEdge         = bumpHunterPLowHigh[1]
                bumpHighEdge        = bumpHunterPLowHigh[2]

                print "bump range: %s GeV - %s GeV"%(bumpLowEdge,bumpHighEdge)
                print "BumpHunter stat = %s"%bumpHunterStatValue
                print "BumpHunter p-value = %s +/- %s"%(bumpHunterPValue, bumpHunterPValueErr)

                bumpHunterSigmas = ROOT.Math.normal_quantile(1.-bumpHunterPValue, 1.)
                bumpHunterErrSigmas = ROOT.Math.normal_quantile(1.-bumpHunterPValueErr, 1.)
                print "BumpHunter sigmas = %s"%bumpHunterSigmas
                                
                gpval.append(bumpHunterPValue)
                gpvalerr.append(bumpHunterPValueErr)

                # Lydia
                gpvalsigmas.append(bumpHunterSigmas)
                gpvalerrsigmas.append(bumpHunterErrSigmas)
                #

                glowedge.append(bumpLowEdge)
                ghighedge.append(bumpHighEdge)
                glumi.append(float(lumi))
                highestlumi = lumi

                #hpval.Fill(float(bumpHunterPValue))
        if highestlumi == 0: raise SystemExit('\n***Zero highest lumi*** regex issue')
                
        print '\nhighest luminosity tested = %s pb^-1'%highestlumi
    
        #------------------------------------------
        #TEST
        #raise SystemExit('\n***TEST*** exit')
        #continue
        #------------------------------------------        
        
        #------------------------------------------
        #plot sensitivity scan for given mass value

        glumi = np.array(glumi)
        glumierr = np.zeros_like(glumi)
        gpval = np.array(gpval)
        gpvalerr = np.array(gpvalerr)
        # Lydia
        #gpvalsigmas = np.array(gpvalsigmas)
        #gpvalerrsigmas = np.array(gpvalerrsigmas)
        # 
        glowedge = np.array(glowedge)
        ghighedge = np.array(ghighedge)

        ROOT.gStyle.SetPadTickY(0)
        ROOT.gStyle.SetPadRightMargin(0.15)

        #------------------------------------------
        #canvas and pads
        c1 = ROOT.TCanvas('c1', 'c1', 400, 50, 800, 600)
        yd = 0.5 #bottom pad (pad2) height
        pad1 = ROOT.TPad("pad1","pad1",0, yd, 1, 1)#top
        pad2 = ROOT.TPad("pad2","pad2",0, 0, 1, yd)#bottom

        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.0)
        #pad1.SetBorderMode(0)
        pad1.Draw()

        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(0.2)
        #pad2.SetBorderMode(0)
        pad2.Draw() 


        #------------------------------------------
        #pad 1
        pad1.cd()
        pad1.Clear()
        pad1.SetLogx(0)
        pad1.SetLogy(0)

        #band
        #gr.GetXaxis().SetTitle('luminosity [fb^{-1}]')
        gr = ROOT.TGraphAsymmErrors(len(glumi), glumi, (ghighedge+glowedge)/2, np.zeros_like(glumi), np.zeros_like(glumi), (ghighedge-glowedge)/2, (ghighedge-glowedge)/2) #n,x,y,exl,exh,eyl,eyh
        gr.SetTitle('range')
        gr.GetYaxis().SetTitle('m_{jj} [GeV]')
        gr.GetYaxis().SetLabelFont(43)
        gr.GetYaxis().SetLabelSize(20)
        gr.GetYaxis().SetTitleFont(43)
        gr.GetYaxis().SetTitleSize(20)
 	
        gr.GetYaxis().SetDecimals(True) 

        gr.SetMarkerColor(ROOT.kRed)
        gr.SetLineColor(ROOT.kRed)
        gr.SetFillColor(ROOT.kRed)

        #..........................................
 	#draw 
        gr.SetFillStyle(3003)
        #gr.Draw("same pe") #draw points
        #grb = gr.Clone() 
        #grb.Draw("same 3") #draw band

        gr.Draw("ape") #draw points
        grb = gr.Clone()
        grb.Draw("same3") #draw band
        c1.cd()
        c1.Update()

        #------------------------------------------
 	#pad 2
        pad2.cd()
        pad2.Clear()
        pad2.SetLogx(0)

        pad2.SetLogy(1)
 	
        #p-value graph
        gpv = ROOT.TGraphErrors(int(len(glumi)), glumi, gpval, glumierr, gpvalerr)
        gpv.SetTitle('p-value')
        gpv.GetXaxis().SetTitle('luminosity [pb^{-1}]')
        gpv.GetXaxis().SetLabelFont(43)
 	
        gpv.GetXaxis().SetLabelSize(20)
        gpv.GetXaxis().SetTitleFont(43)
        gpv.GetXaxis().SetTitleSize(20)
        gpv.GetYaxis().SetTitle('BumpHunter p-value')
        gpv.GetYaxis().SetLabelFont(43)
        gpv.GetYaxis().SetLabelSize(20)
        gpv.GetYaxis().SetTitleFont(43)
        gpv.GetYaxis().SetTitleSize(20)
        gpv.GetXaxis().SetTitleOffset(2.5)#this needs to be adjusted depending on yd
        
        gpv.SetFillColor(ROOT.kBlack)
        gpv.SetFillStyle(3002)
        gpv.SetMaximum(2.)
        gpv.SetMinimum(1e-10)
   
        gpv.SetMinimum(1e-7)
        gpv.Draw("apel") #HERE after SetMinimum()
        #gpv.Draw("same 3") #draw band #NOTE clashes with the line
          
	#draw sigma lines
        s = ROOT.TLatex()
        s.SetNDC(False)
        s.SetTextFont(42)
        s.SetTextColor(2)
        s.SetTextSize(0.05)
        lxs = ROOT.TLine()
        lxs.SetLineColor(ROOT.kRed)
        lxs.SetLineWidth(2)
  
        for i in range(6):
            if i<1: continue
            pvalxs = 1. - ROOT.Math.normal_cdf(i)
            lxs.DrawLine(gpv.GetXaxis().GetXmin(), pvalxs, gpv.GetXaxis().GetXmax(), pvalxs)
            s.DrawLatex(gpv.GetXaxis().GetXmax()*1.1,pvalxs,'%0.f#sigma'%i)
       
        # ---------- TEST Finding intersection -------------------- FIXME remove
        # save graph and 5 sigma line so can find 5 sigma intersection later on 
        #outputFileName = './plotting/SensitivityStudy/plots/%s/lumiscan'%(folderextension)+model+''+mass+'GeV'+parstag+'.root'
        #outputFile = ROOT.TFile.Open(outputFileName,'RECREATE')
        #outputFile.cd()
        #gpvswapped = ROOT.TGraphErrors(int(len(gpval)), gpval, glumi, gpvalerr, glumierr)
        #fivesig = 1. - ROOT.Math.normal_cdf(5)
        #print "5SIG"
        #print fivesig
        #print "LUM"
        #print gpvswapped.Eval(fivesig)
        #for i in range(10001):
        #  bumpHunterPValue = gpv.Eval(i)
        #  sigma = ROOT.Math.normal_quantile(1.-bumpHunterPValue, 1.)
        #  if sigma >=5:
        #    print "DISCOVERY"
        #    print i
        #    print sigma
        #    break
        #gpvsigmas = ROOT.TGraphErrors(int(len(glumi)), glumi, gpvalsigmas, glumierr, gpvalerrsigmas)
        #x = np.linspace(gpvsigmas.GetXaxis().GetXmin(),gpvsigmas.GetXaxis().GetXmax(), len(gpval)) 
        #y = np.empty(100)
        #y.fill(5)
        #fivesiggraph = ROOT.TGraph(len(x), x,y)
        #print gpvalsigmas
        #print gpvalerrsigmas 
        #tools = Morisot()
        #print "INTERCEPT"
        #test = float(gpval[0])# is 0.0:#print np.isinf(gpvsigmas)
        #print test
        #print type(test)# is 0.0:#print np.isinf(gpvsigmas)
        #if float(gpval[0]) != 0.0:
        #  print "NONZERO"
        #  print tools.calculateIntersectionOfGraphs(gpvsigmas,fivesiggraph,True) #,True,False) # doLogGraph1=False, doLogGraph2=False) :
        #else:
        #  print "ZERO"
        #print gpvalsigmas[0]
        #print type(gpvalsigmas[0])
        #print tools.calculateIntersectionOfGraphs(fivesiggraph,gpv)#,False,True) #,True,False) # doLogGraph1=False, doLogGraph2=False) :
        # --------------- End of TEST ---------------------------------------- 

        #draw graph (again) on top
        gpv.Draw("same pel") #HERE after SetMinimum()
        #gpv.Draw("same 3") #draw band #NOTE clashes with the line
 	
        #..........................................
        #position
        lXmin=0.20 #0.60
        lXmax=0.30 #0.70
        lYmin=0.10 #0.61
        lYmax=0.24 #0.75
 	
        #ATLAS
        a = ROOT.TLatex()
        a.SetNDC()
        a.SetTextFont(73) #72
        a.SetTextColor(1)
        a.SetTextSize(30) #0.04
        a.DrawLatex(lXmin,lYmax+0.24,'ATLAS')
 
        #internal
        p = ROOT.TLatex()
        p.SetNDC()
        p.SetTextFont(43) #42
        p.SetTextColor(1)
        p.SetTextSize(30) #0.04
        #p.DrawLatex(lXmin+0.13,lYmax+0.24,'internal')
        p.DrawLatex(lXmin+0.13,lYmax+0.24,'simulation internal')

        #notes
        n = ROOT.TLatex()
        n.SetNDC()
        n.SetTextFont(43) #42
        n.SetTextColor(1)
        n.SetTextSize(20) #0.03
        n.DrawLatex(lXmin,lYmax+0.16,'#sqrt{s} = 13 TeV')
        if model == 'QStar':
            n.DrawLatex(lXmin,lYmax+0.08,'m_{q*} = %.1f'%(float(mass)/1e3)+' TeV')
            print 'm_{q*} = %.1f'%(float(mass)/1e3)+' TeV'
            print
        else:
            n.DrawLatex(lXmin,lYmax+0.08,'m_{'+model+'} = %.1f'%(float(mass)/1e3)+' TeV')

        if par is str(3):
            n.DrawLatex(lXmin,lYmax+0.0,'3 parameter fit')
        if par is str(4):
            n.DrawLatex(lXmin,lYmax+0.0,'4 parameter fit')
        if par is str(5):
            n.DrawLatex(lXmin,lYmax+0.0,'5 parameter fit')

        c1.SetLogx(0)
        c1.SetLogy(0)
        c1.Update()
        #c1.WaitPrimitive()
        c1.SaveAs('./plotting/SensitivityStudy/plots/%s/lumiscan'%(folderextension)+model+''+mass+'GeV'+parstag+'.pdf') #ORIGINAL

        #------------------------------------------
        #plot p-values histogram
        #NOTE plot just one (high, low cross section) mass value
        if ( (float(mass) == 13000. and model == 'QStar')  or (float(mass) == 12000. and model == 'QBH0') ):

            c1.Clear()
            c1.SetLogx(0)
            c1.SetLogy(0)
            c1.Update()
            
            #hpval.GetYaxis().SetTitle('entries/bin')
            #hpval.GetXaxis().SetTitle('BumpHunter p-value')
            #hpval.Draw('histo')

            #ATLAS
            a.DrawLatex(lXmin,lYmax+0.12,'ATLAS')
            
            #internal
            p.DrawLatex(lXmin+0.13,lYmax+0.12,'internal')
    
            #notes
            n.DrawLatex(lXmin,lYmax+0.08,'#sqrt{s} = 13 TeV')
            if model == 'QStar':
                n.DrawLatex(lXmin,lYmax+0.04,'q* mass = '+str(float(mass)/1e3)+' TeV') 
            else:
                n.DrawLatex(lXmin,lYmax+0.04,model+' mass = '+str(float(mass)/1e3)+' TeV') 

            if par is 3:
                n.DrawLatex(lXmin,lYmax+0.0,'3 parameter fit')
            if par is 4:
                n.DrawLatex(lXmin,lYmax+0.0,'4 parameter fit')
            if par is 5:
                n.DrawLatex(lXmin,lYmax+0.0,'5 parameter fit')

            #c1.WaitPrimitive()
            c1.SaveAs('./plotting/SensitivityStudy/plots/%s/pvalues.'%(folderextension)+model+'_'+mass+'_GeV_'+parstag+'.pdf')
