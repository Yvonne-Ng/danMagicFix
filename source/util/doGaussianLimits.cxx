// ***************************************************************
// This file was created using the CreateProject.sh script
// for project MjjProject.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "TH1.h"
#include <iostream>
#include "TEnv.h"
#include "TError.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TVectorD.h"
#include "TVector2.h"
#include "TStopwatch.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TClass.h"

#include "Bayesian/MjjFitter.h"
#include "Bayesian/MjjHistogram.h"
#include "Bayesian/MjjFitFunction.h"
#include "Bayesian/MjjStatisticalTest.h"
#include "Bayesian/MjjChi2Test.h"
#include "Bayesian/MjjPseudoExperimenter.h"
#include "Bayesian/MjjStatisticsBundle.h"
#include "Bayesian/MjjSignificanceTests.h"
#include "Bayesian/MjjBumpHunter.h"
#include "Bayesian/BonusFitFunctions.h"

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>
#include <BAT/BCH1D.h>

#include "Bayesian/MjjBATModel.h"
#include "Bayesian/MjjBATProcess.h"
#include "Bayesian/MjjBATShapeChangingSyst.h"
#include "Bayesian/MjjBATScaleChangingSyst.h"
#include "Bayesian/MjjBATTemplateSyst.h"

#include "Bayesian/FoldingFunctions.h"

void UsagedoGaussianLimits(){
  std::cout << "Usage: " << std::endl;
  std::cout << " --ratio     : width to mass ratio, 0.10 -> 10% width, 0 uses variable resolution" << std::endl;
  std::cout << " --config    : config file to read details from (filenames, JES)" << std::endl;
  std::cout << " --infile    : Gaussian sample input file" << std::endl;
  std::cout << " --useweight : " << std::endl;
  std::cout << "Examples: " << std::endl;
  std::cout << " Fixed Gaussians 10% width: \n"
	    << "  doGaussianLimits --ratio 0.10 --config configurations/GenericGaussians.config"
	    << std::endl;
  std::cout << " Variable resolution width, from sample file, apply PDF and acceptance weights: \n"
	    << "  doGaussianLimits --ratio 0.0 --config configurations/GenericGaussiansAK.config --infile mtestpdf43_8000_w7.root --useweight"
	    << std::endl;
}

double interpretJESFile(TH1D* theHist, double massGeV) {
  int inbin;
  double frac;
  // Hard coded values correspond to high mass analysis beginning with 1 TeV point, summer 2016
  // Updated by Kate on July 4
  if (massGeV < 2000) {inbin = 1; if (massGeV < 1000) frac = 0.0; else frac = (massGeV-1000.0)/1000.0;}
  else if (massGeV >= 2000 && massGeV < 2500) {inbin = 2; frac = (massGeV-2000)/500.0;}
  else if (massGeV >= 2500 && massGeV < 3000) {inbin = 3; frac = (massGeV-2500)/500.0;}
  else if (massGeV >= 3000 && massGeV < 3500) {inbin = 4; frac = (massGeV-3000)/500.0;}
  else if (massGeV >= 3500 && massGeV < 4000) {inbin = 5; frac = (massGeV-3500)/500.0;}
  else if (massGeV >= 4000 && massGeV < 4500) {inbin = 6; frac = (massGeV-4000)/500.0;}
  else if (massGeV >= 4500 && massGeV < 5000) {inbin = 7; frac = (massGeV-4500)/500.0;}
  else if (massGeV >= 5000 && massGeV < 5500) {inbin = 8; frac = (massGeV-5000)/500.0;}
  else if (massGeV >= 5500 && massGeV < 6000) {inbin = 9; frac = (massGeV-5500)/500.0;}
  else if (massGeV >= 6000 && massGeV < 6500) {inbin = 10; frac = (massGeV-6000)/500.0;}
  else if (massGeV >= 6500 && massGeV < 7000) {inbin = 11; frac = (massGeV-6500)/500.0;}
  else if (massGeV >= 7000) {inbin = 12; frac = 1.0;}
  else {
    std::cout << "Unexpected result!" << std::endl;
    return -1.0;
  }
  double val;
  if (inbin==10) val = theHist->GetBinContent(inbin);
  else val = (theHist->GetBinContent(inbin+1) - theHist->GetBinContent(inbin))*frac + theHist->GetBinContent(inbin);
  std::cout << "Assigning val " << val << " for a mass " << massGeV << std::endl;
  return val;
}
//void Fold(TH1D *h_MJJtruth, TH2D *FM, TH1D *h_MJJreco);
//void FillArrayToHisto(TH1D *h_MJJtruth, double *recoFold);
//void fillDirac(TH1D *histo, double shiftMass);

int main(int argc,char **argv)
{
  ////////////////////////////////////////////////////////////
  // Initialisation: this part should be user-modified

  // Start counting time
  TStopwatch totaltime;
  totaltime.Start();

  // Start reading input configuration
  TString configFile;
  double thisRatio = -1;
  bool useResolutionWidth = false;
  bool useSubRange = false;
  bool doSfold = false;
  bool useSeed = false;
  TString Fname;// = "QCDPythia8";
  double rangeLow = 0;
  double rangeHigh = -1;
  int ip=1;
  unsigned int seed = 1234;
  while (ip<argc) {

    if (string(argv[ip]).substr(0,2)=="--") {
      //config file
      if (string(argv[ip])=="--config") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          configFile=argv[ip+1];
          ip+=2;
        } else {std::cout<<"\nno config file inserted"<<std::endl; break;}
      }
      //number of chains to run
      else if (string(argv[ip])=="--ratio") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          sscanf(argv[ip+1],"%lf",&thisRatio);
          ip+=2;
        } else {std::cout<<"\nno ratio specified"<<std::endl; break;}
      }
      //use resolution width
      else if (string(argv[ip])=="--useresolutionwidth") {
          useResolutionWidth = true;
          ip+=1;
      }
      
      // seed for random number
			else if (string(argv[ip]) == "--seed") {
				if (ip+1<argc) {
					sscanf(argv[ip+1],"%u",&seed);
					ip+=2;
					std::cout << "Random seed: " << seed << std::endl;
          useSeed = true;
				} else {std::cout<<"\nno random seed specified"<<std::endl; break;}
			}
      // do only subrange of masses
      else if (string(argv[ip])=="--rangelow") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          sscanf(argv[ip+1],"%lf",&rangeLow);
          useSubRange = true;
          ip+=2;
        } else {std::cout<<"\nno ratio specified"<<std::endl; break;}
      }
      else if (string(argv[ip])=="--rangehigh") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          sscanf(argv[ip+1],"%lf",&rangeHigh);
          useSubRange = true;
          ip+=2;
        } else {std::cout<<"\nno ratio specified"<<std::endl; break;}
      }
      //use folding
      else if (string(argv[ip])=="--doSfold") {
          doSfold = true;
          ip+=1;
      }
      //config file
      else if (string(argv[ip])=="--Fname") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          Fname=argv[ip+1];
          ip+=2;
        }
      }
      //unknown command
      else {
        std::cout<<"\n" << string(argv[0]) <<": command '"<<string(argv[ip])<<"' unknown"<<std::endl;
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") ip+=2;
        else ip+=1;
	UsagedoGaussianLimits();
      } }//end if "--command"

    else { //if command does not start with "--"
      std::cout << "\n" << string(argv[0]) << ": command '"<<string(argv[ip])<<"' unknown"<<std::endl;
      UsagedoGaussianLimits();
      break;
    }//end if "--"

  }//end while loop

  std::cout << "got to end of reading inputs." << std::endl;
  std::cout << "doSfold: "<< doSfold << '\n';
  if (doSfold) std::cout << "Fname: " << Fname << '\n';
  int ratioForNames = int(1000*thisRatio);
  int internalRatio = int(100*thisRatio);

  // Get config file
  TEnv * settings = new TEnv();
  int status = settings->ReadFile(configFile.Data(),EEnvLevel(0));
  if (status!=0) {
    std::cout<<"cannot read config file"<<std::endl;
    std::cout<<"******************************************\n"<<std::endl;
    return 1;
  }

  // Swift fit or regular fit?
  bool doSwift = settings->GetValue("doSwift",false);

  // Specify files
  string inputFileName = settings->GetValue("inputFileName","");
  string outputFileName = settings->GetValue("outputFileName","");

  // Get data
  TString dataMjjHistoName = settings->GetValue("dataHist","");

  // Do we use extended region for fit?
  bool doExtendedRange = settings->GetValue("doExtendedRange",true);

  // Get center of mass energy
  double Ecm = settings->GetValue("Ecm",13000.0);
  if (rangeHigh<0) rangeHigh = Ecm;

  // PDF information
  bool readFromFile = settings->GetValue("useBW",false);
  TString inFileFormat = settings->GetValue("BWFileFormat","");
  string PDF = settings->GetValue("PDF","");

  // Resolution file name
  string resfilename = settings->GetValue("resolutionFile","");
  vector<double> respars;
  bool useFile = resfilename != ""; //(useResolutionWidth && resfilename != "");
  bool useFunc = false;
  if (resfilename == "") //(useResolutionWidth && resfilename == "")
  {
    useFunc = true;
    for (int i = 1; i < 8; i++) {
      string title = Form("respar%i",i);
      double param = settings->GetValue(title.c_str(),1.0);
      respars.push_back(param);
    }
  }

  // JES size file name
  bool doVarJES = settings->GetValue("doVarJES",true);
  TString JESfilename = settings->GetValue("JESFile","");
  bool useJESFile = (doVarJES && JESfilename != "");
  TH1D * JESVariationHist;
  if (useJESFile) {
    TFile * infileJES = TFile::Open(JESfilename, "READ");
    assert(infileJES);
    infileJES->ls();
    JESVariationHist = (TH1D*) infileJES->Get("QStar_JESUp");
    std::cout << "\nThis is the JES variation hist:" << std::endl;
    JESVariationHist->Print("all");
  }

  // Store them for later
  vector<double> JESShiftsUsed;
  // Open files
  TFile * infile = TFile::Open(inputFileName.c_str(), "READ");
  std::cout << inputFileName << std::endl;
  assert(infile);
  TH1::AddDirectory(kFALSE);
  // Range for data fit
  double minX = settings->GetValue("minXForFit",-1);
  double maxX = settings->GetValue("maxXForFit",-1);

  //------------------------------------------
  // SWIFT only parameters!

  // Range usable (will in fit but won't necessarily be in background estimate)
  // Default is same as minX, maxX
  double minXAvailable = settings->GetValue("swift_minXAvailable",minX);
  double maxXAvailable = settings->GetValue("swift_maxXAvailable",maxX);
  
  // Number of bins to left and right of the window
  int nBinsWindowLeft = settings->GetValue("swift_nBinsLeft",13);
  int nBinsWindowRight = settings->GetValue("swift_nBinsRight",13);

  // What does the window do when you get near an edge?

  // If true, total number of bins remains the same at low mass
  // but eventually by the last bin all of them are on the right hand side.
  // If false, by the final bin the window is only equal to nBinsWindowRight
  bool fixWidthAtLowEnd = settings->GetValue("swift_fixLow",true);

  // If true, total number of bins remains the same at high mass
  // but eventually by the last bin all of them are on the left hand side.
  // If false, by the final bin the window is only equal to nBinsWindowLeft
  // Only applies if truncateHigh is true
  bool fixWidthAtHighEnd = settings->GetValue("swift_fixHigh",false);

  // An alternative to fixWidthAtHighEnd: can make this False and
  // simply have the windows move naturally to the right. This is an
  // option when the end of the fit range isn't the end of the data.
  bool truncateHigh = settings->GetValue("swift_truncateHigh",false);

  // Set up numbers.
  double nSigmas = settings->GetValue("nSigmas",3.0);
  // Do we do bkg error?
  bool doFitError = settings->GetValue("doFitError",false);
  int nFitsInBkgError = settings->GetValue("nFitsInBkgError",100);

  // New: control for choosing nominal & alternate functions.
  int functionNumber = settings->GetValue("functionCode",9);
  int nPars = settings->GetValue("nParameters",3);
  int alternateFuncNumber = 0, altNPars = 0;

  // Parameters for fit
  std::cout << "Setting the "<< nPars << " parameters to "; // Lydia
  vector<double> paramDefaults;
  vector<bool> areParamsFixed;
  for (int par = 1; par<nPars+1; par++) {
    string title = Form("parameter%i",par);
    double param = settings->GetValue(title.c_str(),1.0);
    paramDefaults.push_back(param);
    title = Form("fixParameter%i",par);
    bool isFixed = settings->GetValue(title.c_str(),false);
    areParamsFixed.push_back(isFixed);
    std::cout << "   (" << par-1 << ") :" << param;
  }
  std::cout << std::endl;

  //folding: open the folding matrix
  //TFile *f_FM;
  //TH2D *h_FM_nominal;
  if(doSfold)
  {
    LoadFoldingMatrix(Fname);
    //f_FM = new TFile(Form("/afs/cern.ch/work/r/rhankach/workDir/jet_exotic/transferMatrix/TM.%s.root",Fname.Data()));
    //h_FM_nominal = (TH2D*)f_FM->Get("Nominal/TM_normalized_with_eff");
  }
  //h_FM_nominal->Print("all");

  // Parameters for alternate fit, if requested
  // Note it may not have been run in search phase, so
  // do not count on picking up defaults from the input file.
  bool doFitFunctionChoiceError = settings->GetValue("doFitFunctionChoiceError",false);

  vector<double> altParDefaults;
  vector<bool> altAreParsFixed;
  if (doFitFunctionChoiceError) {
    alternateFuncNumber = settings->GetValue("alternateFunctionCode",4);
    altNPars = settings->GetValue("alternateNParameters",4);
    for (int par = 1; par < altNPars+1; par++) {
      string title = Form("altparameter%i",par);
      double param = settings->GetValue(title.c_str(),1.0);
      altParDefaults.push_back(param);
      std::cout<<"altparam "<<param<<std::endl; // Lydia
      title = Form("fixAltParameter%i",par);
      bool isFixed = settings->GetValue(title.c_str(),false);
      altAreParsFixed.push_back(isFixed);
      std::cout << "Setting alternative function param "<< par << " to "<< param << std::endl;
    }
  }

  // Do we do lumi?
  bool doLumiError = settings->GetValue("doLumiError", false);
  double luminosityErr = settings->GetValue("luminosityErr",0.90);

  // Do we do JES?
  bool doJES = settings->GetValue("doJES",true);
  double sigmaJESShift = settings->GetValue("sigmaJESShift",0.10);
  int nJES = settings->GetValue("nJES",25);

  // Do we calculate expected limits?
  bool doExpected = settings->GetValue("doExpected",false);
  // If so, how many?
  int nPseudoExpForExpected = settings->GetValue("nPEForExpected",5);

  // Print various settings
  std::cout.setf(std::ios::boolalpha);
  std::cout << "Doing fit quality uncertainty: " << doFitError << std::endl;
  std::cout << "Doing fit choice uncertainty: " << doFitFunctionChoiceError << std::endl;

  // Get ranges for density of gaussian points
  double startBonusPoints = settings->GetValue("startBonusPoints",-1);
  double startPoints = settings->GetValue("startPoints",500);
  double stopDensePoints = settings->GetValue("stopDensePoints",2000);
  double stopMediumPoints = settings->GetValue("stopMediumPoints",3000);
  double stopPoints = settings->GetValue("stopPoints",5000);
  
  

  // Create output name with ratio if necessary
  if (useSubRange) {
     TString addon = Form("_low%d_high%d",(int)(rangeLow),(int)(rangeHigh));
     outputFileName.append(addon);
  }
  if (useSeed) {
     TString addon = Form("_seed%d",(int)(seed));
     outputFileName.append(addon);
  }
  if (useResolutionWidth) {
    outputFileName.append("_resolutionwidth.root");
  } else if (thisRatio >= 0) {
    outputFileName.append(Form("_%d.root",(int)(1000*thisRatio)));
    std::cout << "Ratio specified. Now using output name " << outputFileName << std::endl;
  } else {
    outputFileName.append(".root");
    std::cout << "No ratio specified. Now using output name " << outputFileName << std::endl;
  }
  
  TFile * outfile = TFile::Open(outputFileName.c_str(), "RECREATE");

  // Start code proper

  double interval = (nSigmas*2)/((double)(nJES-1));
  vector<double> jesSigmas;
  jesSigmas.push_back(0.);
  for (int i=0; i<nJES; i++) {
    double value = -1*nSigmas + i*interval;
    if (!(fabs(value - 0) < 1e-3)) jesSigmas.push_back(value);
  }

  // Get and store data
  TH1D * basicInputHisto = (TH1D*) infile->Get(dataMjjHistoName);
  MjjHistogram theHistogram((TH1D*)basicInputHisto, 111);

  // Set up values for calculating gaussians
  // Fill with big range of values, but only look at those that fall within reasonable
  // bounds of data spectrum.
  vector<double> masses;
  if (startBonusPoints > 0) for (double m=startBonusPoints; m<startPoints; m+= 25.) masses.push_back(m);
  for (double m=startPoints; m<stopDensePoints; m += 25.) masses.push_back(m);
  for (double m=stopDensePoints; m<stopMediumPoints; m += 50.) masses.push_back(m);
  for (double m=stopMediumPoints; m<=stopPoints; m += 100.) masses.push_back(m);

  vector<double> widthToMass;
  if (thisRatio>=0)  widthToMass.push_back(thisRatio);
  else if (!useResolutionWidth){
    widthToMass.push_back(0.07);
    widthToMass.push_back(0.10);
    widthToMass.push_back(0.15);
  } else {
    std::cout << "Adding one value" << std::endl;
    widthToMass.push_back(0);
  }

  // Get resolution tgraph if making Gaussians of resolution width
  TGraphAsymmErrors * resolutiongraph = 0;
  TFormula myfunc("resolutionFunc", "sqrt(([2]*[2])/(x*x) + ([1]*[1])/x + ([0]*[0]))");
  if (useFile) {
    TH1::AddDirectory(kFALSE);
    TFile * mresf = TFile::Open(resfilename.c_str());
    TCanvas * c3 = (TCanvas*) mresf->Get("C1");
//    resolutiongraph = (TGraphAsymmErrors*) ((c3->GetListOfPrimitives())->At(0));
    resolutiongraph = (TGraphAsymmErrors*) (c3->GetPrimitive("Graph"));
    mresf->Close();
    delete mresf;
  } else if (useFunc) {
    //TFormula myfunc("resolutionFunc","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x");//,770,8000); // was TF1
    
    for (int par = 0; par < 3; par++) {
      std::cout << "Setting parameter value " << respars.at(par) << std::endl;
      myfunc.SetParameter(par,respars.at(par));
    }
    resolutiongraph = new TGraphAsymmErrors();
    for (int x = 421; x < 2200; x++) {
      resolutiongraph->SetPoint(x-421,double(x),myfunc.Eval(x));
      std::cout << "resolution at " << x << " = " << myfunc.Eval(x) << std::endl;
    }
  }

  std::cout << "Made function." << std::endl;

  ////////////////////////////////////////////////////////////
  // Set up fit function.

  // Lydia changed this part so picks start and end bins that are fully spanned by range minX-maxX
  // e.g. user picks minX maxX of 2000-3000 only bins fully within this range are fitted
  // rather than fitting e.g. 1950-3500
  // Now matches SearchPhase.cxx

  int firstBin, lastBin;
  if (minX < theHistogram.GetHistogram().GetBinLowEdge(theHistogram.GetFirstBinWithData()) || minX < 0) firstBin = theHistogram.GetFirstBinWithData();
  else firstBin = theHistogram.GetHistogram().FindBin(minX)+1;

  if (maxX > theHistogram.GetHistogram().GetBinLowEdge(theHistogram.GetLastBinWithData()+1) || maxX < 0) lastBin = theHistogram.GetLastBinWithData();
  else lastBin = theHistogram.GetHistogram().FindBin(maxX)-1;
  double minXForFit = theHistogram.GetHistogram().GetBinLowEdge(firstBin);
  double maxXForFit = theHistogram.GetHistogram().GetBinLowEdge(lastBin+1);

  std::cout << "Fitting " << minXForFit << " to " << maxXForFit << " GeV, or bins " << firstBin << " to " << lastBin << std::endl;

  /*int firstBin, lastBin;
  if (theHistogram.GetHistogram().FindBin(minX) < theHistogram.GetFirstBinWithData() || minX < 0) firstBin = theHistogram.GetFirstBinWithData();
  else firstBin = theHistogram.GetHistogram().FindBin(minX)+1;
  if (theHistogram.GetHistogram().FindBin(maxX) > theHistogram.GetLastBinWithData() || maxX < 0) lastBin = theHistogram.GetLastBinWithData();
  else lastBin = theHistogram.GetHistogram().FindBin(maxX)-1;
  double minXForFit = theHistogram.GetHistogram().GetBinLowEdge(firstBin);
  double maxXForFit = theHistogram.GetHistogram().GetBinLowEdge(lastBin+1);*/

  // Create fit function(s)
  MjjFitFunction * theMjjFitFunction;
  MjjFitFunction * theAlternateFunction;
  vector<std::pair<int,MjjFitFunction**> > functionsAndCodes;
  functionsAndCodes.push_back(std::make_pair(functionNumber,&theMjjFitFunction));
  functionsAndCodes.push_back(std::make_pair(alternateFuncNumber,&theAlternateFunction));

  for (unsigned int index = 0; index < functionsAndCodes.size(); index++) {

    int thisFuncCode = functionsAndCodes.at(index).first;
    std::cout << functionsAndCodes.at(index).second << std::endl;
    switch (thisFuncCode) {
      case 1 :
        std::cout << "Creating UA2 fit function." << std::endl;
        *functionsAndCodes.at(index).second = new UA2FitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 2 :
        std::cout << "Creating CDF (1995) fit function." << std::endl;
        *functionsAndCodes.at(index).second = new CDFFitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 3 :
        std::cout << "Creating CDF (1997) fit function." << std::endl;
        *functionsAndCodes.at(index).second = new CDF1997FitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 4 :
        std::cout << "Creating standard dijet function." << std::endl;
        *functionsAndCodes.at(index).second = new FourParamFitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 5 :
        std::cout << "Creating TeV Gravity function." << std::endl;
        *functionsAndCodes.at(index).second = new ThreeParamFitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 6 :
        std::cout << "Creating 5-parameter function, floating sqrt(s)." << std::endl;
        *functionsAndCodes.at(index).second = new FiveParamSqrtsFitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 7 :
        std::cout << "Creating 5-parameter function, log(x)^2 term." << std::endl;
        *functionsAndCodes.at(index).second = new FiveParamLog2FitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 8 :
        std::cout << "Creating 6-parameter function." << std::endl;
        *functionsAndCodes.at(index).second = new SixParamFitFunction(minXForFit,maxXForFit,Ecm);
        break;
      case 9 :
        std::cout << "Creating 3-parameter function for Run II search." << std::endl;
        *functionsAndCodes.at(index).second = new ThreeParam2015FitFunction(minXForFit,maxXForFit,Ecm);
        break;
            case 10 :
                std::cout << "Creating second Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction2(minXForFit,maxXForFit,Ecm);
                break;
            case 11 :
                std::cout << "Creating third Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction3(minXForFit,maxXForFit,Ecm);
                break;
            case 12 :
                std::cout << "Creating fourth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction4(minXForFit,maxXForFit,Ecm);
                break;
            case 13 :
                std::cout << "Creating fifth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction5(minXForFit,maxXForFit,Ecm);
                break;
            case 14 :
                std::cout << "Creating sixth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction6(minXForFit,maxXForFit,Ecm);
                break;
            case 15 :
                std::cout << "Creating seventh Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction7(minXForFit,maxXForFit,Ecm);
                break;
            case 16 :
                std::cout << "Creating eigth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction8(minXForFit,maxXForFit,Ecm);
                break;
            case 17 :
                std::cout << "Creating ninth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction9(minXForFit,maxXForFit,Ecm);
                break;
            case 18 :
                std::cout << "Creating tenth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction10(minXForFit,maxXForFit,Ecm);
                break;
            case 19 :
                std::cout << "Creating gamma gamma with no exponent, 2 par." << std::endl;
                *functionsAndCodes.at(index).second = new TwoParGammaGamma(minXForFit,maxXForFit,Ecm);
                break;
            case 20 :
                std::cout << "Creating gamma gamma with no exponent, 4 par." << std::endl;
                *functionsAndCodes.at(index).second = new FourParGammaGamma(minXForFit,maxXForFit,Ecm);
                break;
            case 21 :
                std::cout << "Creating gamma gamma with 1/3 exponent, 2 par." << std::endl;
                *functionsAndCodes.at(index).second = new TwoParGammaGammaWithThird(minXForFit,maxXForFit,Ecm);
                break;
            case 22 :
                std::cout << "Creating gamma gamma with 1/3 exponent, 3 par." << std::endl;
                *functionsAndCodes.at(index).second = new ThreeParGammaGammaWithThird(minXForFit,maxXForFit,Ecm);
                break;
            case 23 :
                std::cout << "Creating gamma gamma with 1/3 exponent, 4 par." << std::endl;
                *functionsAndCodes.at(index).second = new FourParGammaGammaWithThird(minXForFit,maxXForFit,Ecm);
                break;


     }
    std::cout << functionsAndCodes.at(index).second << std::endl;
  }

  theMjjFitFunction->SetParameterDefaults(paramDefaults);
  theMjjFitFunction->RestoreParameterDefaults();
  for (int par = 0; par<nPars; par++) {
    std::cout << "Fixing parameter " << par << " to " << areParamsFixed.at(par) << std::endl;
    theMjjFitFunction->GetParameter(par)->SetFixParameter(areParamsFixed.at(par));
  }

  // Create alternate function
  // This one really needs you to tell it what the defaults should be.
  if (doFitFunctionChoiceError) {
    theAlternateFunction->SetParameterDefaults(altParDefaults);
    theAlternateFunction->RestoreParameterDefaults();
  }

  std::cout << "Starting initial fit." << std::endl;

  // Create fitter
  MjjFitter theSilentFitter;

  // Set up for Swift if necessary
  if (doSwift) {
    theSilentFitter.SetSwiftParameters(minXForFit,maxXForFit,minXAvailable,maxXAvailable,nBinsWindowLeft,nBinsWindowRight,fixWidthAtLowEnd,fixWidthAtHighEnd,truncateHigh);
  }

  ////////////////////////////////////////////////////////////
  // Loop over all mass points in all ranges.

  vector<TH1D> bestEstimates;
  vector<double> massesUsed;
  vector<vector<double> >  modelCLVectors;
  for (unsigned int model=0; model<widthToMass.size(); model++) {
    std::cout << "In width " << model << std::endl;
    vector<double> CLsForModel;
    for (unsigned int thisMass=0; thisMass<masses.size(); thisMass++) {
      double thisMeanMass = masses.at(thisMass);
      std::cout << "Mass " << thisMeanMass << std::endl;
      if (useSubRange) {
        if (thisMeanMass < rangeLow || rangeHigh <=thisMeanMass) {
//          CLsForModel.push_back(-100);
          continue;
        }
      }

      ////////////////////////////////////////////////////////////
      // Create signal histogram and JES variations

      // Histogram base for setting binning
      TH1D binbase(*(TH1D*) infile->Get(dataMjjHistoName));
      vector<TH1D> SignalTemplates;

      double thisWidthToMass;
      if (useResolutionWidth) thisWidthToMass =  myfunc.Eval(thisMeanMass);//resolutiongraph->Eval(thisMeanMass);
      else thisWidthToMass = widthToMass.at(model);

      std::cout << myfunc.Eval(thisMeanMass) << " " << thisWidthToMass << " " << useResolutionWidth<< std::endl;
      std::cout << "this width is " << thisMeanMass << " * " << thisWidthToMass << " = " << thisMeanMass*thisWidthToMass <<std::endl;

      double ResolutionWidth=0; //when using folding, this value will be added to thisWidthToMass just for the test of distance to edge
      if(doSfold) ResolutionWidth = resolutiongraph->Eval(thisMeanMass);

      // Only compute point if it is more than 2 sigma from an edge
      if ((thisMeanMass - 1.0*thisMeanMass*sqrt(pow(thisWidthToMass,2)+pow(ResolutionWidth,2)) < minXForFit) // changed from 2
         or (thisMeanMass + 1.0*thisMeanMass*sqrt(pow(thisWidthToMass,2)+pow(ResolutionWidth,2)) > maxXForFit)) {
//        CLsForModel.push_back(-100);
        std::cout<<"Mass point too close to edge (i.e. less than 2 sigma away)"<<std::endl;
        continue;
      }

      if((!readFromFile)||(inFileFormat==NULL)){

	// Create JES shifted samples
    double JESShift;
    if (useJESFile) {
      JESShift = interpretJESFile(JESVariationHist,thisMeanMass);
    } else JESShift = sigmaJESShift;
    JESShiftsUsed.push_back(JESShift);

	for (unsigned int i=0; i<jesSigmas.size(); i++) {

	  double shiftMass = thisMeanMass + jesSigmas.at(i)*JESShift*thisMeanMass;
	  double onesigma = shiftMass*thisWidthToMass;
	  double shift3SigDown = shiftMass - 3*onesigma;
	  double shift3SigUp = shiftMass + 3*onesigma;


	  TF1 GenericGaussian("signal",Form("TMath::Gaus(x,%f,%f)",shiftMass,onesigma),
			      shift3SigDown, shift3SigUp);

	  TH1D JESShifted(binbase);
    JESShifted.ClearUnderflowAndOverflow(); JESShifted.Scale(0);
	  JESShifted.SetName(Form("JES_%f",jesSigmas.at(i)));

    TH1D *JESShifted_truth;
    TH1D *SignalFill;
    if(doSfold)
    {
      JESShifted_truth = new TH1D("truth_signal","",h_FM_nominal->GetXaxis()->GetNbins(), h_FM_nominal->GetXaxis()->GetXbins()->GetArray());
      SignalFill=JESShifted_truth;
    }
    else SignalFill=&JESShifted;

	  for (int ibin=1; ibin<SignalFill->GetNbinsX()+1; ibin++) {
	    if (SignalFill->GetBinLowEdge(ibin) >= shift3SigUp ||
		  SignalFill->GetBinLowEdge(ibin+1) <= shift3SigDown) {
	      SignalFill->SetBinContent(ibin,0.);
	      SignalFill->SetBinError(ibin,0.);
	    } else {
	      double a = SignalFill->GetBinLowEdge(ibin);
	      double b = SignalFill->GetBinLowEdge(ibin)+SignalFill->GetBinWidth(ibin);
	      double content = GenericGaussian.Integral(a,b);
	      SignalFill->SetBinContent(ibin,content);
	      SignalFill->SetBinError(ibin,0.);
	    }
	  }

    if(doSfold)
    {
      if(thisWidthToMass<1E-9) fillDirac(JESShifted_truth, shiftMass);
      //double recoFold[JESShifted.GetXaxis()->GetNbins()];
      Fold(JESShifted_truth, h_FM_nominal, &JESShifted);
      delete JESShifted_truth;
      //FillArrayToHisto(&JESShifted,recoFold);
    }//

	  SignalTemplates.push_back(JESShifted);
	}
      }else{

        //exit if no PDF specified
        if (PDF == "") return 0;

        // Open file with this mass and width
        std::cout << "thisMeanMass,ratioForNames are " << thisMeanMass << " " << ratioForNames << std::endl;
        TString infilename = Form(inFileFormat,(int) thisMeanMass,ratioForNames);
	TFile mf(infilename);

        double integral = 0;
        for (unsigned int i=0; i<jesSigmas.size(); i++) {

          TH1D JESShifted(binbase);
          JESShifted.SetName(Form("JES_%f",jesSigmas.at(i)));

          // Get hist from file
          std::cout << "About to get " << Form("finalsignal_PDF%s_width%d_JES%d",PDF.c_str(),internalRatio,i) << std::endl;
          TString histname = Form("finalsignal_PDF%s_width%d_JES%d",PDF.c_str(),internalRatio,i);
          TH1D fullhist(*(TH1D*)mf.Get(histname));

          // Rebin into correct histogram
          for (int bin = 1; bin < JESShifted.GetNbinsX()+1; bin++) {
            double binstart = JESShifted.GetBinLowEdge(bin) + 0.1;
            double binend = JESShifted.GetBinLowEdge(bin) + JESShifted.GetBinWidth(bin) - 0.1;
            JESShifted.SetBinContent(bin,fullhist.Integral(fullhist.FindBin(binstart),fullhist.FindBin(binend)));
            JESShifted.SetBinError(bin,0.0);
          }

          // Normalise to 1
          if (i==0) integral = JESShifted.Integral();
          JESShifted.Scale(1/integral);

          SignalTemplates.push_back(JESShifted);

        }
      }

      int nomIndex=-1, highestJESIndex=-1;
      for (unsigned int i=0; i<jesSigmas.size(); i++) {
        if (jesSigmas.at(i)==0.) nomIndex = i;
        if (jesSigmas.at(i)==nSigmas) highestJESIndex = i;
      }
      TH1D nominal(SignalTemplates.at(nomIndex));
      nominal.SetName("nominal");

      std::cout << "Highest JES index is " << highestJESIndex << std::endl;

      TH1D * nominalFitTemplate = (TH1D*) nominal.Clone();
      nominalFitTemplate->SetName("nominalFitTemplate");
      MjjHistogram thisSigHisto(nominalFitTemplate);

      ///////////////////////////////////////////////////////////////////////////
      // Recalculate range for process using new signal
      int extendedLastBin = lastBin;
      if (doExtendedRange) {
        for (int i=0; i<theHistogram.GetHistogram().GetNbinsX()+2; i++) {
          if (basicInputHisto->GetBinContent(i) > 0 ||
            nominalFitTemplate->GetBinContent(i) > 0 ||
            SignalTemplates.at(highestJESIndex).GetBinContent(i) > 0)
           extendedLastBin = i;
        }
      }
      std::cout << "New extended range ends at " << basicInputHisto->GetBinLowEdge(extendedLastBin+1);

      ////////////////////////////////////////////////////////////
      // Create background template

      theSilentFitter.SetSignalTemplate(nominalFitTemplate);

      MjjHistogram theBackgroundFromSimultaneousFit = theSilentFitter.FitAndGetBkgWithDataErr(*theMjjFitFunction,theHistogram,nFitsInBkgError);
      TH1D bkgTemplate = (TH1D) theBackgroundFromSimultaneousFit.GetHistogram();

      std::cout << "Finished initial fit." << std::endl;

      bkgTemplate.Print("all");

      // What were fitted params?
      vector<double> fittedParameters = theMjjFitFunction->GetCurrentParameterValues();

      std::cout << "After fit, parameters are  ";
      for (unsigned int i=0; i<fittedParameters.size(); i++) {
        std::cout << "   (" << i << "):" << fittedParameters.at(i);
      }
      std::cout << std::endl;

      // Set up variation template for background
      TH1D bkgVariation(bkgTemplate);
      bkgVariation.SetName("BkgVariation");
      for (int i=0; i<bkgTemplate.GetNbinsX()+1; i++) {
        bkgVariation.SetBinContent(i,bkgTemplate.GetBinError(i));
        bkgVariation.SetBinError(i,0.);
      }

      // Set up fit function choice templates for background
      // Old method: Use difference between nominal and alternate in this data.
/*
      MjjHistogram bkgFromFiveParam = theSilentFitter.FitAndGetBkgWithMCErr(*theAlternateFunction,theHistogram);
      TH1D variedfit = bkgFromFiveParam.GetHistogram();
      variedfit.SetName("fitplusonesigmachoice");
      vector<std::pair<double,TH1D*> > fitvariedtemplates;
      fitvariedtemplates.push_back(std::make_pair(0,&bkgTemplate));
      fitvariedtemplates.push_back(std::make_pair(1,&variedfit));
*/
      // New method: Use average difference between the two functions across a range of pseudoexperiments
      // where the pseudoexperiments are thrown from data.
      MjjHistogram bkgFromAltFunc;
      MjjHistogram bkgWithFuncChoiceErr;
      TH1D fitChoicePlus1Sig(bkgTemplate);
      if (doFitFunctionChoiceError) {

        std::cout << "Beginning alternate fit" << std::endl;
        bkgFromAltFunc = theSilentFitter.FitAndGetBkgWithNoErr(*theAlternateFunction,theHistogram);
        vector<double> altfittedParameters = theAlternateFunction->GetCurrentParameterValues();

        // Use results of fit to data to set start parameters for alternate function in fit to PEs
        theAlternateFunction->SetParameterDefaults(altfittedParameters);

        // What were fitted params?
        std::cout << "After alternate fit, parameters are  ";
        for (unsigned int i=0; i<altfittedParameters.size(); i++) {
          std::cout << "   (" << i << "):" << altfittedParameters.at(i);
        }
        std::cout << std::endl;

        // Get histogram with errors equal to RMS of distance between functions
        // across PEs in each bin. This is the uncertainty that will be used for calculating
        // the search phase p-value using systematics.
        std::cout << "Getting fit choice uncertainty" << std::endl;
        bkgWithFuncChoiceErr = theSilentFitter.FitAndGetBkgWithFitDiffErr(*theMjjFitFunction,*theAlternateFunction,theHistogram,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN(),nFitsInBkgError,true);
        std::cout << "Fit function choice uncertainty computed" << std::endl;

        bkgWithFuncChoiceErr.GetHistogram().Print("all");

        // Create histogram of format which will be used as alternate function
        // to define function choice uncertainty when directionality matters:
        // that is, in the limit setting stage.
        fitChoicePlus1Sig.SetName("alternateFitChoiceSyst_setByRMSTimesDirectionality");
        for (int bin = 1; bin < fitChoicePlus1Sig.GetNbinsX()+1; bin++) {
  	  double nomquantity = theBackgroundFromSimultaneousFit.GetHistogram().GetBinContent(bin);
	  double varquantity = bkgFromAltFunc.GetHistogram().GetBinContent(bin);
	  if (varquantity < nomquantity) fitChoicePlus1Sig.SetBinContent(bin,nomquantity - bkgWithFuncChoiceErr.GetHistogram().GetBinError(bin));
	  else fitChoicePlus1Sig.SetBinContent(bin,nomquantity + bkgWithFuncChoiceErr.GetHistogram().GetBinError(bin));
        }
      }
      vector<std::pair<double,TH1D*> > fitvariedtemplates;
      fitvariedtemplates.push_back(std::make_pair(0,&bkgTemplate));
      fitvariedtemplates.push_back(std::make_pair(1,&fitChoicePlus1Sig));

      ///////////////////////////////////////////////////////////////////////////
      // Estimate appropriate signal range
      double maxNSignal = 0;
      double stepSize;
      if (thisMass < 1000) stepSize = 200;
      else if (thisMass < 1500) stepSize = 50;
      else stepSize = 10;  // increase nSignalEvents nStepSize events at a time
      double maxLsoFar = 0;
      double thisL = 1;
      MjjLogLikelihoodTest theLogLTest;
      TH1D sumHist(bkgTemplate);
      sumHist.SetName("HistSigEstimate");
      double integral = nominalFitTemplate->Integral();
      nominalFitTemplate->Scale(1./integral);
      while (maxLsoFar / thisL < 1e5) {

        for (int bin=0; bin<sumHist.GetNbinsX(); bin++) {
          sumHist.SetBinContent(bin,bkgTemplate.GetBinContent(bin) + nominalFitTemplate->GetBinContent(bin)*maxNSignal);
          sumHist.SetBinError(bin,0.);
        }
        MjjHistogram theSignalTemplate(&sumHist);
        TH1D weighthist = theHistogram.GetWeightsHistogram();
        theSignalTemplate.SetEffectiveFromBasicAndWeights(&weighthist);

        double logL = theLogLTest.DoTest(theHistogram,theSignalTemplate,firstBin,extendedLastBin);
        thisL = exp(-logL);

        if (std::isnan(thisL)) break;
        if (thisL > maxLsoFar) maxLsoFar = thisL;

        maxNSignal+=stepSize;
      }

      std::cout << "Estimated maximum signal range = " << maxNSignal << std::endl;

      ////////////////////////////////////////////////////////////
      // Create signal template

      TH1D nominalSignal(nominal);
      nominalSignal.SetName("nominalSignal");
      TH1D variationSignal(nominal);
      variationSignal.SetName("variationSignal");
      double totalBinContent = 0;
      for (int i=0; i<nominalSignal.GetNbinsX()+1; i++) {
        totalBinContent+=nominalSignal.GetBinContent(i);
        nominalSignal.SetBinContent(i,0.);
        nominalSignal.SetBinError(i,0.);
      }
      for (int i=0; i<variationSignal.GetNbinsX()+1; i++) {
        if (variationSignal.GetBinContent(i)==0) continue;
        double cont = variationSignal.GetBinContent(i);
        variationSignal.SetBinContent(i,cont/totalBinContent);
        variationSignal.SetBinError(i,0.);
      }

      ////////////////////////////////////////////////////////////
      // Create histograms for systematics

      // JES histos
      vector<std::pair<double,TH1D*> > signalJESVariations; signalJESVariations.clear();
      for (unsigned int i=0; i<jesSigmas.size(); i++) signalJESVariations.push_back(std::make_pair(jesSigmas.at(i),&SignalTemplates.at(i)));

      ////////////////////////////////////////////////////////////
      // Create limit-setting tools.

      // set nicer style for drawing than the ROOT default
      BCAux::SetStyle();

      // open log file
      TString logfile("log.txt");
      BCLog::OpenLog(logfile);
      BCLog::SetLogLevel(BCLog::detail);

      // create new MjjBATModel object
      MjjBATModel * m = new MjjBATModel();

      m->MCMCSetRandomSeed(seed);
      std::cout<<"Seed of m is set"<<std::endl;

      // create a new summary tool object
      BCSummaryTool * summary = new BCSummaryTool(m);

      // My stuff: adding signal, processes, systematics ...
      // Background and signal normalization are assumed (for now) to be gaussian around given values
      m->SetData(theHistogram,firstBin,extendedLastBin);

      if(doFitError) m->AddProcess("BKG",true,-nSigmas,nSigmas);
      else m->AddProcess("BKG",1,0.0,0.0);
      m->AddProcess("SIGNAL",true,0.0,maxNSignal);

      ////////////////////////////////////////////////////////////
      // Set BAT tools to all the new values.

      m->SetTemplate("BKG",bkgTemplate,bkgVariation);
      vector<string> newparams = m->GetProcess(m->GetProcessIndex("BKG"))->GetParamNames();
      for (unsigned int i=0; i<newparams.size(); i++) {
        m->SetPriorGauss(newparams.at(i).c_str(),0.,1.);
      }

      m->SetTemplate("SIGNAL",nominalSignal,variationSignal);
      vector<string> sigparams = m->GetProcess(m->GetProcessIndex("SIGNAL"))->GetParamNames();
      for (unsigned int i=0; i<sigparams.size(); i++) {
        m->SetPriorConstant(sigparams.at(i).c_str());
      }

      // Create JES uncertainty.
      if (doJES) {
        m->AddSystematic("JES",-nSigmas,nSigmas);
        m->SetPriorGauss("JES",0.,1.);
        MjjBATTemplateSyst * JESOnSignal = new MjjBATTemplateSyst(false);
        JESOnSignal->SetSpectra(signalJESVariations);
        m->SetSystematicVariation("SIGNAL","JES",JESOnSignal);
      }

      if (doFitFunctionChoiceError) {

        m->AddSystematic("FUNCCHOICE",0,1);
        m->SetPriorGauss("FUNCCHOICE",0.,1.);
        MjjBATTemplateSyst * FuncChoiceOnBkg = new MjjBATTemplateSyst(true);
        FuncChoiceOnBkg->SetSpectra(fitvariedtemplates);
        m->SetSystematicVariation("BKG","FUNCCHOICE",FuncChoiceOnBkg);

      }

      // Create luminosity uncertainty. Add variations in processes.
      if (doLumiError) {
        m->AddSystematic("LUMI",-nSigmas,nSigmas);
        m->SetPriorGauss("LUMI",0.,1.);
        MjjBATScaleChangingSyst * LumiOnSignal = new MjjBATScaleChangingSyst(luminosityErr);
        LumiOnSignal->SetScaleFromBinContent();
        m->SetSystematicVariation("SIGNAL","LUMI",LumiOnSignal);
      }

      ////////////////////////////////////////////////////////////
      // Marginalise and get observed limits.

      // run MCMC and marginalize posterior wrt. all parameters
      // and all combinations of two parameters
      m->MarginalizeAll();

      BCParameter * normPar = m->GetParameter(m->GetProcess(1)->GetNormalisationParamIndex());
      BCH1D * marginalizedSig = m->GetMarginalized(normPar);
      double trueCL = marginalizedSig->GetLimit(0.95);
      std::cout << "CL for mass " << thisMeanMass << " is " << trueCL << std::endl;

      ////////////////////////////////////////////////////////////
      // Perform pseudoexperiments
      vector<double> CLPseudo;
      if (doExpected) {
        std::cout << "*******************\n* EXPECTED LIMITS\n*******************" <<std::endl;
        // Create output object
        // For now, just use best fit output parameters
        // However, first, need to fix this to give background-only hypothesis
        int sigparam = m->GetParIndicesProcess(m->GetProcessIndex("SIGNAL")).at(0); // signal has just one param
        std::cout << "Fixing parameter " << m->GetParameter(sigparam)->GetName() << std::endl;

        // Marginalise with signal fixed to zero
        m->GetParameter(sigparam)->Fix(0.0);
        m->MarginalizeAll();

        // Fix signal to zero for generation of pseudoexperiments
        vector<double> useParams;
        std::cout << "Best-fit parameter values for PE generation are: " << std::endl;
        for (unsigned int i=0; i<m->GetNParameters(); i++) {
          double bestval = m->GetBestFitParameter(i);
          std::cout << m->GetParameter(i)->GetName() << ": " << bestval << std::endl;
          useParams.push_back(bestval);
        }

        // Un-fix signal.
        m->GetParameter(sigparam)->Unfix();

        // Make analysis facility
        MjjBATAnalysisFacility * maf = new MjjBATAnalysisFacility(m,seed);
        maf->SetFlagMCMC(true);

        // Current model params still at best-fit with no signal,
        // so use those for pseudoexperiments
        outfile->cd();
        TTree * pseudotree = maf->BuildEnsembles(useParams, 10000 );
        TTree * outputtree = maf->PerformEnsembleTest(pseudotree, nPseudoExpForExpected, firstBin, lastBin );
        
        std::ostringstream ss;
        ss << thisMeanMass << "_" << ratioForNames;
        outputtree->Write(("ensemble_tree_"+ss.str()).c_str());
        
        /*std::ostringstream sss;
        sss << m->GetProcessIndex("SIGNAL");

        double expCL;
        
        outputtree->SetBranchAddress(("95quantile_marginalized_" + sss.str()).c_str(), &expCL);
        for (int i = 0; i<outputtree->GetEntries(); i++) {
          outputtree->GetEvent(i);
          CLPseudo.push_back(expCL);
          std::cout << "individual pseudo CL = " << expCL  << std::endl;
        }*/
      }

      delete m;
      delete summary;

      // close log file
      BCLog::CloseLog();

      CLsForModel.push_back(trueCL);
      massesUsed.push_back(thisMeanMass);
      //PECLsForModel.push_back(CLPseudo);

    }
    modelCLVectors.push_back(CLsForModel);
  }

  ////////////////////////////////////////////////////////////
  // Save everything.

  // Print plots of input variations into file for checking
  std::cout << "Writing output files " << outputFileName << std::endl;

  outfile->cd();

  std::cout << "In output file " << std::endl;

  for (unsigned int i=0; i<widthToMass.size(); i++) {
    vector<double> theseCLs = modelCLVectors.at(i);
    TVectorD CLsForModelMasses(theseCLs.size());
    for (unsigned int j=0; j<theseCLs.size(); j++) {
      CLsForModelMasses[j] = theseCLs.at(j);
    }
    CLsForModelMasses.Write(Form("CLsPerMass_widthToMass%d",ratioForNames));
    
  }

  TVectorD writeMassesUsed(massesUsed.size());
  for (unsigned int i=0; i<massesUsed.size(); i++) {
    writeMassesUsed[i] = massesUsed.at(i);
  }
  writeMassesUsed.Write("massesUsed");

  TVectorD writeJESShiftsUsed(JESShiftsUsed.size());
  for (unsigned int i=0; i<JESShiftsUsed.size(); i++) {
    writeJESShiftsUsed[i] = JESShiftsUsed.at(i);
  }
  writeJESShiftsUsed.Write("JESShifts");

  if (useResolutionWidth) //(useFile || useFunc)
    resolutiongraph->Write("resolution");

  outfile->Close();
  std::cout << "Closed outfile" << std::endl;
  infile->Close();
  std::cout << "Closed infile" << std::endl;

  totaltime.Stop();
  std::cout << "Ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

  delete outfile;
  delete infile;

  return 0;

}
