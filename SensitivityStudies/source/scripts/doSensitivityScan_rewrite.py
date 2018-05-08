#!/bin/python

#******************************************
#perform sensitivity scan by looping over all masses on a range of luminosity values
#EXAMPLE python -u doSensitivityScan.py --config <config file> --tag <tag> --batch --debug
#based on code from Hanno
#re-written by yvonne


#******************************************
#import stuff
import sys, os, math, argparse, ROOT
import sensitivityTools
import numpy as np
import os.path
from fileNamingTool import *
import json

#def removeOldLabelledFile(spResultDir, label, mass, window, gauss):
#    """test to see if old labelled file of the mass point /window width and gaussian width is still arround, if so delete them"""
#    OldTaggedFileMayExist=True
#    i=0
#    while(OldTaggedFileMayExist and i<20):
#        try:
#            testFile=findLabelledFileName(spResultDir, label, mass , window, gauss=15)
#            #print("TaggedFileToBe ", nosignalFN)
#            #os.rename(testFile,noSignalFN)
#            #if the above returns something instead of raising an error,
#            i=i+1
#        except ValueError:
#            OldTaggedFileMayExist=False
#            print("no other duplicatedly labelled file from previous runs ")
#
#        else :
#            print ("removing this file: ", testFile)
#            os.remove(testFile)

def createWorkTag(window, sigScale, mass):
    workTag="ww"+Window+"sigScale"+sigScale+"_mass"+mass
    return workTag

def injectionFileExist(localdir,signalScale, config):
    """Check if the singal injected bkg file exist for this signal scale and window"""
    """signalplusbackground.TrijetAprSelection.Gauss_width7.SigNum100.mjj_Gauss_sig__smooth.root"""
    injectionFile="signalplusbackground."+config["SeriesName"]+"."+config["signalModel"]+"."+"SigNum"+str(signalScale)+"."+config["histBasedNameSig"].format("")+".root"
    if not os.path.isdir(localdir+"/"+config["signalInjectedFileDir"]):
        raise RuntimeError
    if injectionFile in os.listdir(localdir+"/"+config["signalInjectedFileDir"]):
        return True
    else:
        return False

def injectionFile(localdir, signalScale, config):
    """returning the step2 output file of bkgnd data file with signal injected"""
    """signalplusbackground.TrijetAprSelection.Gauss_width7.SigNum100.mjj_Gauss_sig__smooth.root"""
    injectionFile=localdir+"/"+config["signalInjectedFileDir"]+"/"+"signalplusbackground."+config["SeriesName"]+"."+config["signalModel"]+"."+"SigNum"+str(signalScale)+"."+config["histBasedNameSig"].format("")+".root"
    return injectionFile

def fluctuatedBkgFile(localdir, config):
    print("fluctuatedBkgFile: ", localdir+"/"+config["QCDFileDir"]+"/"+config["QCDFile"])
    return localdir+"/"+config["QCDFileDir"]+"/"+config["QCDFile"]

def spExcludeWindow(windowVector):
    excludeWindow=windowVector[0]
    excludeWindowLow=windowVector[1]
    excludeWindowHigh=windowVector[2]
    return (excludeWindow, excludeWindowLow, excludeWindowHigh)

def doSensitivityScan(args):
    #---------Setting local directory
    localdir = os.path.dirname(os.path.realpath(__file__))
    #---------Open config file
    json_data = open(args.config)
    config=json.load(json_data)
    try:
        json_data = open(args.config)
        config=json.load(json_data)
    except:
        print("can't open json file. Abort.")
        raise RuntimeError
    #---------fluctuates bkg data file --using a fixed step1 file
    #eliminating step 1
        #** Should fluctuate the bkg data everytime
#        #** Use the same file now for testing
#    bkgFile=TFile(config["QCDFile"])
#    bkgHist=bkgFile.Open(config["histBaseNameBkg"])
#    #bkg File
#    bkgFile=TFile(config["QCDFile"])
#    bkgHist=bkgFile.Open(config["histBaseNameBkg"])
    #signalInjectedFiles
    for window in config["windows"]:
        # do the no signal file
        noSignalFile=FluctuatedBkgFile(localdir, config)
        command3NoSig="python step03_rewrite.py --config {0} --file {1} --mass {2} --window {3}".format(args.config, noSignalFile, "NOSIGNAL", window)
        #setting mass to 0 for no signal files
        os.system(command3NoSig)
        print("nosig search phase done")
        #try:
        #    os.system(command3NoSig)
        #    print(
        #except:
        #    print("noSig Search phase failed")

        previousSPFile=None
        for mass in config["signalMasses"]:
            for signalScale in config["signalScales"]:
            # if the scaled file for this window and signal scale doesn't exsit yet. run step2
                if not injectionFileExist(localdir,signalScale, config):
                    #this is the equilvalant of step02
                    #generateInjectedFile(window, signalScale)
                    command2="python -u step02_rewrite.py --config {0} --sigScale {1} --debug".format(args.config,signalScale)
                    print("command", command2)
                    if args.debug:
                        print ("command for step2: ", command2)
                    try:
                        os.system(command2)
                        print("command step2 ran")
                    except:
                        print("step2 failed. aborting")
                        raise RuntimeError
                    if not injectionFileExist(localdir, signalScale, config):
                        print("didn't produce the correct injection file:", )
                        raise ValueError
                else:
                    print ("file already exist:",  injectionFile(localdir, signalScale, config))
                # Now that we know the injection file exist, run step03

                # step3
                step2File=injectionFile(localdir, signalScale, config)
                command3="python step03_rewrite.py --config {0} --mass {1} --window {2} --file {3} --debug ".format(args.config, mass, window,step2File)
                os.system(command3)
#                try:
#                    print("running step3 for the injection file: ", step2File)
#                    if args.debug:
#                        print ("command for step3: ", command3)
#                    print ("command for step3: ", command3)
#                    os.system(command3)
#                except:
#                    print("step03 failed. aborting")
#                    raise RuntimeError
#

                spFileName=searchPhaseResultName(config["signalModel"],mass, signalScale,window, config["SeriesName"])
                spFile = ROOT.TFile(localdir+"/"+config["spResultDir"]+"/"+spFileName,'READ')
                print(" does the file exsist?")
                print(os.path.isfile(localdir+"/"+config["spResultDir"]+"/"+spFileName))

                print(spFile.Get('excludeWindowNums'))

                try:
                    (excludeWindow, excludeWindowLow, excludeWindowHigh)=spExcludeWindow(spFile.Get('excludeWindowNums'))
                except:
                    print("fit faied for mass {}, window: {}, sigScale: {}")
                    print("skipping mass point")
                    continue
                #-----Tagging files
                #remove old tags
                spBHPValue=spFile.Get("bumpHunterStatOfFitToData")[1]

                removeOldLabelledFile(config["spResultDir"], "JUSTABOVE", mass, window, config["signalModel"], config["SeriesName"])
                removeOldLabelledFile(config["spResultDir"], "JUSTBELOW", mass, window, config["signalModel"], config["SeriesName"])
                removeOldLabelledFile(config["spResultDir"], "NOSIGNAL", mass, window, config["signalModel"], config["SeriesName"])
                #if spBHPValue<0.01 and excludeWindow==1:
                if excludeWindow==1:
                  print "Discovery, with window removal"
                  print "just above SPFile: ", spFileName
                  #---finding the justabove file
                  justAboveFN=justAboveFileName(spFileName)
                  os.rename(localdir+"/"+config["spResultDir"]+"/"+spFileName, localdir+"/"+config["spResultDir"]+"/"+justAboveFN)
                  print("Search phase BH p value: ", spBHPValue)
                  if previousSPFile==None:
                      print("window: {}, mass: {}, signalScale: {} need to start at a lower signal scale".format(window, mass, signalScale))
                      #skipping the mass point
                      break
                  else:# if previousSP file is not none
                  #----finding the just below file name
                      justBelowFN=justBelowFileName(previousSPFile)
                      os.rename(localdir+"/"+config["spResultDir"]+"/"+previousSPFile,localdir+"/"+config["spResultDir"]+"/"+justBelowFN)
                      break #going to the next mass point
                else: #if exclusion of widnow didn't happen
                  previousSPFile=spFileName
                  print("window: {}, mass: {}, signalScale: {} No exclusion yet.".format(window, mass, signalScale))
                  print("Search phase BH p value: ", spBHPValue)



                  #TODO:
                  # fix the signal injection scale increments
                  # write out all the missiong function defined above
                  # run to test plotting results
                  #beautifing:
                  # remove step1 step2 and step3 and turn them into functions


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='%prog [options]')
    parser.add_argument('--config', dest='config', default='../configs/sensitivityScan.Test.config', required=True, help='sensitivity scan config file')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    args = parser.parse_args()
    doSensitivityScan(args)

