import sys, os, math, argparse, ROOT
import plotTools, sensitivityTools
import ROOT as r
import io
import json

#--------Rewrite summary-----------#
#in the loop function, only go to this part of the code if the root file of this luminosity does not already exist. ( or if the mass point looked for doesn't exist
# Remove all the luminosity use SignalEvent instead

#-------End of rewrite Summary-----#

def printConfig(config):
    try:
        print("Starting signal injection for: ", config["SeriesName"])
        print("signalModel: ", config["signalModel"])
        print("signal Masses of: ", config["signalMasses"])
        print("signal hist name config: ", config["histBasedNameSig"])
        print("QCDFile: ", config["QCDFile"])
        print ("bkg histogram name: ", config["histBaseNameBkg"])
        print("signal histogram name: ", config["histBasedNameSig"])
        print("nPseudoExperiment for Bkg: ", config["nPseudoExperiments_bkg"])
        print("nPseudoExperiments_withSig: ", config["nPseudoExperiments_withSig"])
    except:
        print "---Config file may be mission something or is broken--"
        print "----ABORTING----"
        raise RuntimeError

def doesRootFileExist(rootFile):
    if not os.path.isfile(rootFile):
        print ("can't find ", rootFile)
        raise RuntimeError

def makeOutputFileName(config, args):
    localdir = os.path.dirname(os.path.realpath(__file__))
    histNameFromInput=config["histBasedNameSig"].format("")
    outFileName = localdir+config["signalInjectedFileDir"]+"signalplusbackground."+config["SeriesName"]+"."+config["signalModel"]+".SigNum"+str(args.sigScale)+"."+histNameFromInput+".root"
    return outFileName

def makeSignalFileName(inputDir,model, signalMass):
    fileName=inputDir+"/"+model+"/"+"Gauss_mass"+str(signalMass)+"_"+model[6:]+".root"
    return fileName

def makeSignalFileList(config):
    """find the list of signal files to be injected"""
    fileList=[]
    for signalMass in config["signalMasses"]:
        sigFile=makeSignalFileName(config["signalFileDir"], config["signalModel"], signalMass)
        if not os.path.isfile(sigFile):
            print("error this file does not exist:", sigFile)
            raise ValueError
        fileList.append(sigFile)
    return fileList

def injectDataLikeSignal(args):
#---Get dir of the script
    localdir = os.path.dirname(os.path.realpath(__file__))

#---setting the correct error for MC
    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows()
    ROOT.TH2.SetDefaultSumw2()
    ROOT.TH2.StatOverflows()
#---Opening json file
    try:
        json_data = open(args.config)
        config = json.load(json_data)
    except:
        print "Cannot open json config file: ", args.config

        print "---Aborting---"
        raise RuntimeError

    #---debug
    if args.debug:
        printConfig(config)
#--Test if RootFiles exist
    #bkg
    doesRootFileExist(localdir+"/"+config["QCDFileDir"]+config["QCDFile"])
    #signal
    signalFileList=makeSignalFileList(config)
    if args.debug:
        print ("the signal file list being injected is ", signalFileList)
    for signalFile in signalFileList:
        doesRootFileExist(signalFile)
    #TODO test if signal inputfile exist
#create outputFile
    outFileName=makeOutputFileName(config,args)
    outFile=r.TFile(outFileName, "RECREATE")

    #---debug
    if args.debug:
        print("output file name : ", outFile)
#---Get bkg histogram and writing it
# skip if bkg histogram already exist

    bkgFile = ROOT.TFile("/lustre/SCRATCH/atlas/ywng/WorkSpace/r21/r21SwiftNew/SensitivityStudies/source/scripts//../input_dijetISR2018/bkg//Fluctuated_SwiftFittrijet_HLT_j380_inclusiveAprRewdo.root", 'READ')
    bkgHist= bkgFile.Get("basicBkgFrom4ParamFit_fluctuated").Clone()
    #bkgFile = ROOT.TFile(localdir+"/"+config["QCDFileDir"]+"/"+config["QCDFile"], 'READ')
    #print(bkgFile)
    #print(localdir+"/"+config["QCDFileDir"]+"/"+config["QCDFile"])
    #print(config["histBaseNameBkg"])
    #bkgHist = bkgFile.Get(config["histBaseNameBkg"]).Clone()

    bkgHist.SetTitle(config["histBaseNameBkg"]+"_QCD")

    outFile.cd()
    bkgHist.Write()

    histNameFromInput=config["histBasedNameSig"].format("")

    sigHist={}
    for signalFile in signalFileList:
        mass = sensitivityTools.getSignalMass(signalFile)
        sigFile=ROOT.TFile(signalFile)
        print(sigFile)
        # getting the original size histogram
        print("chekc5 in loop")
        histNameSigUpdated=config["histBasedNameSig"].format(str(mass))
        histNameSigUpdated=histNameSigUpdated.encode("ascii")
        sigHist["ori"]=sigFile.Get(histNameSigUpdated)
        print(sigHist)
        sigHist["ori"].SetName(histNameSigUpdated+"_sig")
        sigHist["ori"].SetTitle(histNameSigUpdated+"_sig")
# Scaling the signal histograms up
        sigHist["scaled"]= sigHist["ori"].Clone()
        # the scaling is actually a bit off....
        sigHist["scaled"].Scale(args.sigScale/sigHist["scaled"].Integral())
        sigHist["eff"]=plotTools.getEffectiveEntriesHistogram(sigHist["ori"], config["histBasedNameSig"].format("")+'_eff')
        #fluctuating the scaled histogram
        sigHist["dataLike"] = plotTools.getDataLikeHistYvonneVersion(sigHist["eff"], sigHist["scaled"],histNameSigUpdated+"_sig", 1, config["thresholdMass"])
        # I am skippig the smoothing part
        totHist= bkgHist.Clone()
        totHist.Add(sigHist["dataLike"])
        histNameToHist=histNameSigUpdated+"injectedToBkg"
        totHist.SetName(histNameToHist)
        totHist.SetTitle(histNameToHist)
        if args.debug:
            print("background hist entries: ", bkgHist.Integral())
            print("scale: ",args.sigScale)
            print("signalHistEntries: ", sigHist["dataLike"].Integral())
            print("total his entries: ", totHist.Integral())
        outFile.cd()
        totHist.Write()
        if args.plot:
            totHist.SetLineColor(r.kGreen)
            scaledSigHist=sigHist["dataLike"].Clone()
            scaledSigHist.Scale(totHist.Integral()/sigHist["dataLike"].Integral())
            scaledSigHist.SetLineColor(r.kRed)
            c1=r.TCanvas(1)
            c1.SetLogy()
            bkgHist.Draw()
            scaledSigHist.Draw("same")
            totHist.Draw("same")


            l = plotTools.getLegend(0.65,0.85,2)
            l.AddEntry(bkgHist,"QCDBkgnd", "l")
            l.AddEntry(scaledSigHist,"signal scaled up", "l")
            l.AddEntry(totHist,"total", "l")
            l.Draw()
            c1.SaveAs("../figures/step02Result."+model+"."+mass+"."+args.sigScale+".pdf")

    outFile.Write()
    outFile.Close()




if __name__=="__main__":
    parser = argparse.ArgumentParser(description='%prog [options]')
    parser.add_argument('--config', dest='config', default='', required=True, help='sensitivity scan config file')
    parser.add_argument('--sigScale', dest='sigScale', type=int, required=True, help="signal scale")
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    parser.add_argument('-p', '--plot', dest='plot', action='store_true', default=False, help='plot mode')
    args = parser.parse_args()
    injectDataLikeSignal(args)



