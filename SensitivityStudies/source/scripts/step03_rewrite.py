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
import json

def getSigNumFromStep2FileName(fileName, mass):
    #signalplusbackground.TrijetAprSelection.Gauss_width7.SigNum100.mjj_Gauss_sig__smooth.root
    if mass=="NOSIGNAL":
        return 0
    else:
        print("step2 file name : ", fileName)
        pos0=fileName.find("SigNum")+len("SigNum")
        print("pos0", pos0)
        pos1 = fileName[pos0:].find(".")+pos0
        print("pos1", pos1)
        SigNum=fileName[pos0:pos1]
        return SigNum

def makeSearchPhaseConfig(seriesName,configTemplate, window,mass,  nPseudoExps):
    """using the input config template, make new config for search phase that has different window and nPsuedoExp"""
    #----check for file
    if not os.path.isfile(configTemplate):
        print(" the template input file does not exist: ", configTemplate)
        print(" ----------Aborting-------------")
        raise RuntimeError
    workTag = "sensitivity_"+seriesName+'_'+str(mass)+"ww"+window
    configOutDir="submitConfigs/"+workTag
    if not os.path.isdir(configOutDir):
        os.mkdir(configOutDir)
    configOutName = "submitConfigs/"+workTag+"/Step1_SearchPhase_window{0}.config".format(window)
    configOut = open(configOutName,'w')
    with open(configTemplate) as configInData :
      for line in configInData :
        if "nPseudoExp" in line :
          line = "nPseudoExp  {0}\n".format(nPseudoExps)
        if "swift_nBinsLeft" in line :
          line = "swift_nBinsLeft  {0}\n".format(window)
        elif "swift_nBinsRight" in line :
          line = "swift_nBinsRight  {0}\n".format(window)
        configOut.write(line)
    configOut.close()
    modcommand = "chmod 744 {0}".format(configOut)
    os.system(modcommand)
    return configOutName

def makeCommand(fitFile, fitHistogram, outFileName , configOut):
    """make searchphase command"""
    #remake this outFileName
    commandTemplate = "SearchPhase --config {0} --file {1} --histName {2} --noDE --outputfile {3}".format(configOut,fitFile,fitHistogram,outFileName)
    return commandTemplate




#******************************************
def runSearchPhase(args):


#Setting the local directory
    localdir = os.path.dirname(os.path.realpath(__file__))
#pulling the corrrect json file
    try:
        json_data = open(args.config)
        config = json.load(json_data)
    except:
        print "Cannot open json config file: ", args.config
        print "---Aborting---"
        raise RuntimeError
    #checking if input file exist
    if not os.path.isfile(args.inputFileName):
        print("Search phase input file not found :", args.inputFileName)
        raise keyboardInterupt
    #Get Hist Name
    sigNum=getSigNumFromStep2FileName(args.inputFileName, args.mass)
    if sigNum==0:
        histName=config["histBaseNameBkg"]
        pseudoExps=100
    else:
        histName="mjj_Gauss_sig_{0}_smoothinjectedToBkg".format(str(args.mass))
        pseudoExps=1
    #making a specific output config for this run
    outputConfig=makeSearchPhaseConfig(config["SeriesName"],config["searchPhaseConfigTemplate"], args.window,args.mass,  pseudoExps)

    searchPhaseResultFile=config["spResultDir"]+"/"+searchPhaseResultName(config["signalModel"], args.mass, sigNum, args.window, config["SeriesName"])

    print("sigNum: ", sigNum)
    print(searchPhaseResultFile)
    if sigNum==0:
        searchPhaseResultFile[:-5]+"NOSIGNAL"+".root"

    command=makeCommand(args.inputFileName, histName, searchPhaseResultFile , outputConfig)
    if args.debug:
        print ("command: ", command)


    os.system(command+"| tee "+config["outputLogDir"]+"/"+config["SeriesName"]+"ww"+args.window+"_mass"+args.mass+"_log.txt \n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', dest='config', default='', required=True, help='sensitivity scan config file')
    parser.add_argument('--file', dest='inputFileName', default='', required=True, help='input file')
    # need mass to save output name
    parser.add_argument('--mass', dest='mass', default=0, required=True, help='signal mass value')
    #swiftWindow size
    #should grab ths from the config file in doSensitivity scan
    parser.add_argument('--window', dest='window', required=True, help='SWiFt window')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False, help='debug mode')
    args = parser.parse_args()
    runSearchPhase(args)
