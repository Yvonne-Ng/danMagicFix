import os

def searchPhaseResultName(model,mass, sigScale,window, seriesName):
#searchphase.Gauss_width3.400.gev.0p2.ifb.mjj.4.par.0.seed.default_ww11.root
    template="searchphase.{0}.{1}.{2}.gev.SigEvent{3}.mjj._ww{4}.root"
    searchPhaseFile=template.format(seriesName,model, str(mass), sigScale, window)
    return searchPhaseFile

def labelTest(label):
    if label not in ["JUSTABOVE", "JUSTBELOW", "NOSIGNAL"]:
        print("wrong label: ", label)
        print("must be JUSTABOVE, JUSTBELOW, NOSIGNAL")
        raise ValueError
    else:
        return

#def justAboveFileName(model, mass, sigScale, window):
#    template="searchphase.{0}.{1}.gev.SigEvent{2}.mjj._ww{3}JUSTABOVE.root"
#    searchPhaseFile=template.format(model, str(mass), sigScale, window)
#    return searchPhaseFile
#
#def justBelowFileName(model, mass, sigScale, window):
#    template="searchphase.{0}.{1}.gev.SigEvent{2}.mjj._ww{3}JUSTBELOW.root"
#    searchPhaseFile=template.format(model, str(mass), sigScale, window)
#    return searchPhaseFile
#
def getDiscoveryEventNum(fileName):
    """get discovery event num"""
    #searchphase.Gauss_width7.400.gev.SigEvent50.mjj._ww12JUSTABOVE.root
    pos0=fileName.find("SigEvent")+len("SigEvent")
    pos1=fileName[pos0:].find(".")+pos0
    eventNum=int(fileName[pos0:pos1])
    return eventNum



def justAboveFileName(fileName):
    jafile=fileName[:-5]+"JUSTABOVE"+".root"
    return jafile
def justBelowFileName(fileName):
    jafile=fileName[:-5]+"JUSTBELOW"+".root"
    return jafile

def removeOldLabelledFile(fileDir, label, mass, windowWidth, model, seriesName):
    try:
        fileList=os.listdir(fileDir)
    except:
        print("listing fileDir failed: " , fileDir)

    massStr=str(mass)
    windowStr=str(windowWidth)

    for fileName in os.listdir(fileDir):
        if all( x in fileName for x in [massStr, label, windowStr, model, SeriesName]) :
            print("old labelled file removed:",  fileDir+"/"+fileName)
            os.remove(fileDir+"/"+fileName)

def FluctuatedBkgFile(localdir,config):
    return localdir+"/"+config["QCDFileDir"]+"/"+config["QCDFile"]

def modelNameFromGaussianWidth(width):
  return "Gauss_width"+str(width)




#need a different one for signal
def findLabelledFileName(directory,label, model, mass, window, seriesName):
    labelTest(label)
    #find file name
    massString="."+str(mass)+"."
    modelString=model
    windowString="ww"+str(window)
    if label=="NOSIGNAL":
        targetFiles=[f for f in os.listdir(directory) if modelString in f and windowString in f and label in f and seriesName in f]
    else:
        targetFiles=[f for f in os.listdir(directory) if massString in f and modelString in f and windowString in f and label in f and seriesName in f]
    if len(targetFiles)>1:
        print("Too many target files found: ", targetFiles)
        print("mass string: ", massString)
        print("aborting. ")
        raise RuntimeError
    elif len(targetFiles)==0:
        print("no file found for: ", model, " ", mass, " ", "window: ", window, " ", label)
        print("abordting,")
        raise RuntimeError
    elif len(targetFiles)==1:
        return directory+"/"+targetFiles[0]
    else:
        print ("undefined behavior")
        print("aborting")
        raise RuntimeError

def setLabelForFile(directory,f, label):
    labelTest(label)
    pos0=f.find(".root")
    newName=f[:pos0]+"."+label+f[pos0:]
    os.rename(directory+"/"+f, directory+"/"+newName)







