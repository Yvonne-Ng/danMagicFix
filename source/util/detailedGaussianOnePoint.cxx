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
#include "Bayesian/MjjBATAnalysisFacility.h"

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
  if (massGeV < 2500) return theHist->GetBinContent(1);
  else if (massGeV >= 2500 && massGeV < 3000) return theHist->GetBinContent(2);
  else if (massGeV >= 3000 && massGeV < 4000) return theHist->GetBinContent(3);
  else if (massGeV >= 4000 && massGeV < 4500) return theHist->GetBinContent(4);
  else if (massGeV >= 4500 && massGeV < 5000) return theHist->GetBinContent(5);
  else if (massGeV >= 5000 && massGeV < 5500) return theHist->GetBinContent(6);
  else if (massGeV >= 5500 && massGeV < 6000) return theHist->GetBinContent(7);
  else if (massGeV >= 6000 && massGeV < 6500) return theHist->GetBinContent(8);
  else if (massGeV >= 6500) return theHist->GetBinContent(9);
  else {
    std::cout << "Unexpected result!" << std::endl;
    return -1.0;
  } 
}

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
  int ip=1;
  double thisMeanMass = 0;
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

      // do only subrange of masses
      else if (string(argv[ip])=="--mass") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          sscanf(argv[ip+1],"%lf",&thisMeanMass);
          ip+=2;
        } else {std::cout<<"\nno ratio specified"<<std::endl; break;}
      }

	  // seed for random number
      else if (string(argv[ip]) == "--seed") {
		if (ip+1<argc) {
          sscanf(argv[ip+1],"%u",&seed);
		  ip+=2;
          std::cout << "Random seed: " << seed << std::endl;
		} else {std::cout<<"\nno random seed specified"<<std::endl; break;}
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
  
  std::cout << "Reading config file " <<configFile << std::endl;

  // Specify files
  string inputFileName = settings->GetValue("inputFileName","");
  string outputFileName = settings->GetValue("outputFileName","");
  std::cout << "Intend to generate output file " << outputFileName << std::endl;

  // Get data
  TString dataMjjHistoName = settings->GetValue("dataHist","");

  // Get center of mass energy
  double Ecm = settings->GetValue("Ecm",8000.0);

  // PDF information
  bool readFromFile = settings->GetValue("useBW",false);
  TString inFileFormat = settings->GetValue("BWFileFormat","");
  string PDF = settings->GetValue("PDF","");

  // Resolution file name
  string resfilename = settings->GetValue("resolutionFile","");
  vector<double> respars;
  bool useFile = (useResolutionWidth && resfilename != "");
  bool useFunc = false;
  if (useResolutionWidth && resfilename == "") {
    useFunc = true;
    for (int i = 1; i < 8; i++) {
      string title = Form("respar%i",i);
      double param = settings->GetValue(title.c_str(),1.0);
      respars.push_back(param);
    }
  }

  // JES size file name
  bool doVarJES = settings->GetValue("doVarJES",true);
  TString JESfilename = settings->GetValue("JESFile","./inputs/JESshifts/QStarJESShifts1invfbJES1Component3Up.root");
  bool useJESFile = (doVarJES && JESfilename != "");
  TH1D JESVariationHist;
  if (useJESFile) {
    TFile infileJES(JESfilename);
    JESVariationHist = *(TH1D*)infileJES.Get("CompUp3");
  }
  std::cout << "\nThis is the JES variation hist:" << std::endl;
  JESVariationHist.Print("all");

  // Open files
  TFile * infile = TFile::Open(inputFileName.c_str(), "READ");
  TH1::AddDirectory(kFALSE);

  // Range for data fit
  double minX = settings->GetValue("minXForFit",-1);
  double maxX = settings->GetValue("maxXForFit",-1);

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
  double luminosityErr = settings->GetValue("luminosityErr",0.036);

  // Do we do JES?
  bool doJES = settings->GetValue("doJES",true);
  double sigmaJESShift = settings->GetValue("sigmaJESShift",0.03);
  int nJES = settings->GetValue("nJES",25);

  // Create output name with ratio if necessary
  TString addon = Form("_mass%d",(int)(thisMeanMass));
  outputFileName.append(addon);
  if (useResolutionWidth) {
    outputFileName.append("_resolutionwidth.root");
  } else if (thisRatio >= 0) {
    outputFileName.append(Form("_%d.root",(int)(1000*thisRatio)));
    std::cout << "Ratio specified. Now using output name " << outputFileName << std::endl;
  } else {
    outputFileName.append(".root");
    std::cout << "No ratio specified. Now using output name " << outputFileName << std::endl;
  }

  // Do we calculate expected limits?
  bool doExpected = settings->GetValue("doExpected",false);
  // If so, how many?
  int nPseudoExpForExpected = settings->GetValue("nPEForExpected",5);

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
  MjjHistogram theHistogram((TH1D*)basicInputHisto);

  std::cout << "Contents of the data histogram: " << std::endl;
  basicInputHisto->Print("all");

  // Set up values for calculating gaussians
  // Fill with big range of values, but only look at those that fall within reasonable
  // bounds of data spectrum.

  double widthToMass;
  if (thisRatio>=0)  widthToMass = thisRatio;
  else widthToMass = 0.0;

  // Get resolution tgraph if making Gaussians of resolution width
  TGraphAsymmErrors * resolutiongraph = 0;
  if (useFile) {
    TH1::AddDirectory(kFALSE);
    TFile * mresf = TFile::Open(resfilename.c_str());
    TCanvas * c3 = (TCanvas*) mresf->Get("C1");
    //resolutiongraph = (TGraphAsymmErrors*) ((c3->GetListOfPrimitives())->At(0));
    resolutiongraph = (TGraphAsymmErrors*) (c3->GetPrimitive("Graph"));
    mresf->Close();
    delete mresf;
  } else if (useFunc) {
    TF1 myfunc("resolutionFunc","[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x + [6]*x*x*x*x*x*x;",770,8000);
    for (int par = 0; par < 7; par++) {
      myfunc.SetParameter(par,respars.at(par));
    }
    resolutiongraph = new TGraphAsymmErrors();
    for (int x = 771; x < 8000; x++) {
      resolutiongraph->SetPoint(x-771,double(x),myfunc.Eval(x));
    }
  }

  ////////////////////////////////////////////////////////////
  // Set up fit function.

  int firstBin, lastBin;
  if (minX < theHistogram.GetFirstBinWithData() || minX < 0) firstBin = theHistogram.GetFirstBinWithData();
  else firstBin = theHistogram.GetHistogram().FindBin(minX+1);
  if (maxX > theHistogram.GetLastBinWithData() || maxX < 0) lastBin = theHistogram.GetLastBinWithData();
  else lastBin = theHistogram.GetHistogram().FindBin(maxX-1);
  double minXForFit = theHistogram.GetHistogram().GetBinLowEdge(firstBin);
  double maxXForFit = theHistogram.GetHistogram().GetBinLowEdge(lastBin+1);

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
      case 9 :
        std::cout << "Creating 3-parameter function for Run II search." << std::endl;
        *functionsAndCodes.at(index).second = new ThreeParam2015FitFunction(minXForFit,maxXForFit,Ecm);
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
  theAlternateFunction->SetParameterDefaults(altParDefaults);
  theAlternateFunction->RestoreParameterDefaults();

  std::cout << "Starting initial fit." << std::endl;

  // Create fitter
  MjjFitter theSilentFitter;

  ////////////////////////////////////////////////////////////
  // Loop over all mass points in all ranges.

  std::cout << "Mass " << thisMeanMass << std::endl;

  ////////////////////////////////////////////////////////////
  // Create signal histogram and JES variations
      
  // Histogram base for setting binning
  TH1D binbase(*(TH1D*) infile->Get(dataMjjHistoName));
  vector<TH1D> SignalTemplates;

  double thisWidthToMass;
  if (useResolutionWidth) thisWidthToMass = resolutiongraph->Eval(thisMeanMass);
  else thisWidthToMass = widthToMass;

  std::cout << "this width is " << thisMeanMass << " * " << thisWidthToMass << " = " << thisMeanMass*thisWidthToMass <<std::endl;

  // Only compute point if it is more than 1 sigma from an edge
  if ((thisMeanMass - 1*thisMeanMass*thisWidthToMass < minXForFit)
    or (thisMeanMass + 1*thisMeanMass*thisWidthToMass > maxXForFit)) {
    std::cout << "Calculated " << thisMeanMass << " - " << thisMeanMass*thisWidthToMass << " < " << minXForFit << std::endl;
//        CLsForModel.push_back(-100);
    return 0;
  }

  if((!readFromFile)||(inFileFormat==NULL)){
    // Create JES shifted samples
    for (unsigned int i=0; i<jesSigmas.size(); i++) {
      double JESShift; 
      if (useJESFile) {
        JESShift = interpretJESFile(&JESVariationHist,thisMeanMass);
      } else JESShift = sigmaJESShift;


      std::cout << "About to create JES " << jesSigmas.at(i) << std::endl;
 
      double shiftMass = thisMeanMass + jesSigmas.at(i)*JESShift*thisMeanMass;
      double onesigma = shiftMass*thisWidthToMass;
      double shift3SigDown = shiftMass - 3*onesigma;
      double shift3SigUp = shiftMass + 3*onesigma;
	  	  
      TF1 GenericGaussian("signal",Form("TMath::Gaus(x,%f,%f)",shiftMass,onesigma), 
	     shift3SigDown, shift3SigUp);  

	  
      TH1D JESShifted(binbase);
      JESShifted.SetName(Form("JES_%f",jesSigmas.at(i)));
	  
      for (int ibin=0; ibin<JESShifted.GetNbinsX()+2; ibin++) {
          if (JESShifted.GetBinLowEdge(ibin) >= shift3SigUp ||
              JESShifted.GetBinLowEdge(ibin+1) <= shift3SigDown) {
          JESShifted.SetBinContent(ibin,0.0);
          JESShifted.SetBinError(ibin,0.0);
          } else {
          double a = JESShifted.GetBinLowEdge(ibin);
          double b = JESShifted.GetBinLowEdge(ibin)+JESShifted.GetBinWidth(ibin);
          double content = GenericGaussian.Integral(a,b);
          JESShifted.SetBinContent(ibin,content);
          JESShifted.SetBinError(ibin,0.);
        }
      }
	  
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

  int nomIndex=-1;
  int highestJESIndex = -1;
  for (unsigned int i=0; i<jesSigmas.size(); i++) { 
     if (jesSigmas.at(i)==0.) nomIndex = i;
     if (jesSigmas.at(i)==3.) highestJESIndex = i;
  }
  TH1D nominal(SignalTemplates.at(nomIndex));
  nominal.SetName("nominal");

  TH1D * nominalFitTemplate = (TH1D*) nominal.Clone();
  nominalFitTemplate->SetName("nominalFitTemplate");
  MjjHistogram thisSigHisto(nominalFitTemplate);

  ///////////////////////////////////////////////////////////////////////////
  // Recalculate range for process using new signal
  int extendedLastBin = 0;
  for (int i=0; i<theHistogram.GetHistogram().GetNbinsX()+2; i++) {
    if (basicInputHisto->GetBinContent(i) > 0 ||
        nominalFitTemplate->GetBinContent(i) > 0 ||
        SignalTemplates.at(highestJESIndex).GetBinContent(i) > 0)
      extendedLastBin = i;
  }
  std::cout << "Extended range ends at " << nominalFitTemplate->GetBinLowEdge(extendedLastBin+1) << std::endl;

  ////////////////////////////////////////////////////////////
  // Create background template

  theSilentFitter.SetSignalTemplate(nominalFitTemplate);

  MjjHistogram theBackgroundFromSimultaneousFit = theSilentFitter.FitAndGetBkgWithDataErr(*theMjjFitFunction,theHistogram,nFitsInBkgError);
  TH1D bkgTemplate = (TH1D) theBackgroundFromSimultaneousFit.GetHistogram();

  std::cout << "Finished initial fit." << std::endl;

  // For the random collection of histograms to be saved
  std::vector<TH1D> extraHistograms;

  // Store +/- 1 sigma in nominal fit
  TH1D Bkg_plus1(theBackgroundFromSimultaneousFit.GetHistogram());
  Bkg_plus1.SetName("nominalBkgFromFit_plus1Sigma");
  TH1D Bkg_minus1(theBackgroundFromSimultaneousFit.GetHistogram());
  Bkg_minus1.SetName("nominalBkgFromFit_minus1Sigma");
  for (int bin =0; bin < theBackgroundFromSimultaneousFit.GetHistogram().GetNbinsX()+2; bin++) {
	double bincontent = theBackgroundFromSimultaneousFit.GetHistogram().GetBinContent(bin);
	double binerr = theBackgroundFromSimultaneousFit.GetHistogram().GetBinError(bin);
	Bkg_plus1.SetBinContent(bin, bincontent + binerr);
	Bkg_minus1.SetBinContent(bin, bincontent - binerr);
  }
  extraHistograms.push_back(Bkg_plus1);
  extraHistograms.push_back(Bkg_minus1);

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

    bkgFromAltFunc = theSilentFitter.FitAndGetBkgWithNoErr(*theAlternateFunction,theHistogram);
    vector<double> altfittedParameters = theAlternateFunction->GetCurrentParameterValues();

    // Use results of fit to data to set start parameters for alternate function in fit to PEs
    theAlternateFunction->SetParameterDefaults(altfittedParameters);

    // Get histogram with errors equal to RMS of distance between functions
    // across PEs in each bin. This is the uncertainty that will be used for calculating
    // the search phase p-value using systematics.
    bkgWithFuncChoiceErr = theSilentFitter.FitAndGetBkgWithFitDiffErr(*theMjjFitFunction,*theAlternateFunction,theHistogram,100,true);

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
 
    extraHistograms.push_back(bkgFromAltFunc.GetHistogram());
    extraHistograms.at(extraHistograms.size()-1).SetName("alternateFitOnRealData");
    extraHistograms.push_back(bkgWithFuncChoiceErr.GetHistogram());
    extraHistograms.at(extraHistograms.size()-1).SetName("nomOnDataWithSymmetricRMSScaleFuncChoiceErr");
    extraHistograms.push_back(fitChoicePlus1Sig);
    extraHistograms.at(extraHistograms.size()-1).SetName("nomOnDataWithDirectedRMSScaleFuncChoiceErr");

 }
  vector<std::pair<double,TH1D*> > fitvariedtemplates;
  fitvariedtemplates.push_back(std::make_pair(0,&bkgTemplate));
  fitvariedtemplates.push_back(std::make_pair(1,&fitChoicePlus1Sig));

  ///////////////////////////////////////////////////////////////////////////
  // Estimate appropriate signal range
  double maxNSignal = 0;
  double stepSize;
  if (thisMeanMass < 1000) stepSize = 200;
  else if (thisMeanMass < 1500) stepSize = 50;
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
  BCLog::OpenLog("log.txt");
  BCLog::SetLogLevel(BCLog::detail);

  // create new MjjBATModel object
  MjjBATModel * m = new MjjBATModel();

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
  TH1D * sigLikeVsNumber = (TH1D*) marginalizedSig->GetHistogram()->Clone();
  double trueCL = marginalizedSig->GetLimit(0.95);
  std::cout << "CL for mass " << thisMeanMass << " is " << trueCL << std::endl;

  vector<TH1D*> systematicposteriors;
  for (unsigned int i=0; i<m->GetNParameters(); i++) {
    if ((int)i==m->GetProcess(1)->GetNormalisationParamIndex()) continue;
    BCH1D * margparam = m->GetMarginalized(m->GetParameter(i));
    TH1D * histparam = (TH1D*) margparam->GetHistogram()->Clone();
    string thisname = m->GetParameter(i)->GetName().data();
    std::cout << "This name is " << thisname << std::endl;
    histparam->SetName(thisname.c_str());
    systematicposteriors.push_back(histparam);
  }

 
  // close log file
  BCLog::CloseLog();

  ////////////////////////////////////////////////////////////
  // Save everything.

  // Print plots of input variations into file for checking
  std::cout << "Writing output file " << outputFileName << std::endl;

  TFile * outfile = TFile::Open(outputFileName.c_str(), "RECREATE");
  outfile->cd();

  ////////////////////////////////////////////////////////////
  // Perform pseudoexperiments

  if (doExpected) {
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
    TTree * pseudotree = maf->BuildEnsembles(useParams, 10000 );
    TTree * outputtree = maf->PerformEnsembleTest(pseudotree, nPseudoExpForExpected, firstBin, extendedLastBin );

    outfile->cd();
    outputtree->Write();

  }


  ////////////////////////////////////////////////////////////
  // Create and save best fit prediction

  MjjSignificanceTests theTestMaker;
  MjjChi2Test thechicalc;
  TVectorD chi2vec(1);
  for (int index=0; index<11; index++) {

    vector<double> globalParamValues;
    if (index==0) for (unsigned int i=0; i<m->GetNParameters(); i++) globalParamValues.push_back(m->GetBestFitParameter(i));
    else {
      int sigparam = m->GetParIndicesProcess(m->GetProcessIndex("SIGNAL")).at(0); // signal has just one param
      for (unsigned int i=0; i<m->GetNParameters(); i++) {
        if (int(i)==sigparam) globalParamValues.push_back(double(index)*trueCL/10.0);
	else {globalParamValues.push_back(m->GetBestFitParameter(i));
	  std::cout << "Using signal " << double(index)*trueCL/10.0 << std::endl; 
        }
      }
    }

    TH1D sigPrediction(*basicInputHisto);
    TH1D bkgPrediction(*basicInputHisto);
    TH1D fullPrediction(*basicInputHisto);

    if (index==0) {
      sigPrediction.SetName("BestSignalPrediction");
      bkgPrediction.SetName("BestBackgroundPrediction");
      fullPrediction.SetName("BestFullPrediction");
    } else {
      TString signame = Form("SignalPrediction_%dpercentOfMax",10*index);
      sigPrediction.SetName(signame);
      TString bkgname = Form("BkgPrediction_%dpercentOfMax",10*index);
      bkgPrediction.SetName(bkgname);
      TString fullname = Form("FullPrediction_%dpercentOfMax",10*index);
      fullPrediction.SetName(fullname);
    }

    vector<std::pair<double,double> > sigvalues = m->ProcessExpectation(1,globalParamValues);
    vector<std::pair<double,double> > bkgvalues = m->ProcessExpectation(0,globalParamValues);
    std::cout << "Length of sigvalues, bkgvalues are: " << sigvalues.size() << ", " << bkgvalues.size() << std::endl;
    std::cout << "Number of bins to check is 1 to " << basicInputHisto->GetNbinsX() << std::endl;
    for (int bin = 1; bin < basicInputHisto->GetNbinsX()+1; bin++) {

      sigPrediction.SetBinContent(bin,sigvalues.at(bin).first);
      bkgPrediction.SetBinContent(bin,bkgvalues.at(bin).first);
      fullPrediction.SetBinContent(bin,sigvalues.at(bin).first+bkgvalues.at(bin).first);

    }  

    MjjHistogram bestFull(&fullPrediction);
    TH1D weights = theHistogram.GetWeightsHistogram();
    bestFull.SetEffectiveFromBasicAndWeights(&weights);
    MjjHistogram bestBkg(&bkgPrediction);
    bestBkg.SetEffectiveFromBasicAndWeights(&weights);

    // Want to get this compared to data
    TH1D residualBestData = theTestMaker.GetResidual(theHistogram, bestFull, firstBin, lastBin);
    TH1D residualBkgToData = theTestMaker.GetResidual(theHistogram, bestBkg, firstBin, lastBin);
    double chi2 = thechicalc.DoTest(theHistogram, bestFull, firstBin, extendedLastBin);

    sigPrediction.Write();
    bkgPrediction.Write();
    fullPrediction.Write();

    if (index==0) {
	residualBestData.Write("residual_bestfitToData");
	residualBkgToData.Write("residual_bkgOnly");
	chi2vec[0] = chi2;
	chi2vec.Write("Chi2BestFitToData");
	std::cout << "Chi2 of best prediction to data is " << chi2 << std::endl;
    } else {
	TString residualFullName = Form("residual_%dpercentOfMax",10*index);
	residualBestData.Write(residualFullName);
    }
  }

  ////////////////////////////////////////////////////////////
  // Create and save signal test plots

  std::cout << "About to make test plots" << std::endl;

  vector<double> parameters;
  vector<TH1D> signalplots;
  vector<TH1D> backgroundplots;

  if (doJES) {

  // Create parameters for plots: max and min of jes (m 3000)
  const int nJESSlices = 5; 
  double paramslices [nJESSlices] = {-1.41,-0.21,0.99,-0.75,0.45}; // first 3 are maximum probability locations, last 2 are minima

  for (int jesindex=0; jesindex < nJESSlices; jesindex++) {
    TH1D thissig(nominalSignal);
    string name;
    if (jesindex < 3) name = "maximum";
    else name = "minimum";
    parameters.clear();
    parameters.push_back(0.0); //background
    parameters.push_back(1.0); //signal
    double jesval = paramslices[jesindex];
    parameters.push_back(jesval); // jes component
    name = Form((name+"_jes%d").c_str(),(int)(jesval*10));
    if (doFitFunctionChoiceError) parameters.push_back(0.0); // function choice
    if(doLumiError) parameters.push_back(0.0); //lumi
    thissig.SetName(name.c_str());
    vector<std::pair<double,double> > sigexpectation = m->ProcessExpectation(1,parameters);
    for (int bin=0; bin<thissig.GetNbinsX(); bin++) {
	thissig.SetBinContent(bin,sigexpectation.at(bin).first);
    }
    thissig.Scale(1.0/thissig.GetMaximum());
    signalplots.push_back(thissig);
  }
  const int provescaling = 25;
  double regslices [provescaling] = {-3.0,-2.75,-2.50,-2.25,-2.0,-1.75,-1.50,-1.25,-1.00,-0.75,-0.50,-0.25,0.0,0.25,0.50,0.75,1.0,1.25,1.50,1.75,2.0,2.25,2.50,2.75,3.00};
  for (int jesindex=0; jesindex < provescaling; jesindex++) {
    TH1D thissig(nominalSignal);
    parameters.clear();
    parameters.push_back(0.0); //background
    parameters.push_back(1.0); //signal
    double jesval = regslices[jesindex];
    parameters.push_back(jesval); // jes component
    string name = "showscaling";
    name = Form((name+"_jes%d").c_str(),(int)(jesval*100));
    if (doFitFunctionChoiceError) parameters.push_back(0.0); // function choice
    if(doLumiError) parameters.push_back(0.0); //lumi
    thissig.SetName(name.c_str());
    vector<std::pair<double,double> > sigexpectation = m->ProcessExpectation(1,parameters);
    for (int bin=0; bin<thissig.GetNbinsX(); bin++) {
	thissig.SetBinContent(bin,sigexpectation.at(bin).first);
    }
//    thissig.Scale(1.0/thissig.GetMaximum());
    signalplots.push_back(thissig);
  }

/*  for (unsigned int jesindex=0; jesindex < jesSigmas.size(); jesindex++) {
    TH1D shouldbesig(SignalTemplates.at(jesindex));
    string name2 = "shouldbesig";
    name2 = Form((name2+"_jes%d").c_str(),(int)(jesSigmas.at(jesindex)*10.0));
    shouldbesig.SetName(name2.c_str());
    signalplots.push_back(shouldbesig);
  }*/

  } // end of doJES flag

  ////////////////////////////////////////////////////////////
  // Create and save background test plots 

  std::cout << "About to make background test plots" << std::endl;

  if (doFitError) {

    double bkgslices [7];
    for (int i=0; i<7; i++) {
	bkgslices[i] = -3.0 + i;
    }

    backgroundplots.clear();
    for (int bkgindex=0; bkgindex < 7; bkgindex++) {
	TH1D thissig(nominalSignal);
	string name = Form("bkg_fiterr%d",bkgindex);
	thissig.SetName(name.c_str());
	parameters.clear();
	parameters.push_back(bkgslices[bkgindex]); //background
	parameters.push_back(0.0); //signal
	parameters.push_back(0.0); // jes
	if (doFitFunctionChoiceError) parameters.push_back(0.0); // function choice
	if (doLumiError) parameters.push_back(0.0); //lumi
	vector<std::pair<double,double> > bkgexpectation = m->ProcessExpectation(0,parameters);
	for (int i=0; i<thissig.GetNbinsX()+2; i++) {
		thissig.SetBinContent(i,bkgexpectation.at(i).first);
	}
	backgroundplots.push_back(thissig);
    }
  }

  basicInputHisto->SetName("data");
  basicInputHisto->Write();

  sigLikeVsNumber->SetName("likelihoodFunction");
  sigLikeVsNumber->Write();

  std::cout << "Wrote likelihood function" << std::endl;

  TVectorD CLOfRealLikelihood(1);
  CLOfRealLikelihood[0] = trueCL;
  CLOfRealLikelihood.Write("CLOfRealLikelihood");

  std::cout << "Wrote trueCL" << std::endl;

  // Original signal shape
  TH1D nomSigHisto = SignalTemplates.at(nomIndex);
  nomSigHisto.Write("nominalSignal");

  // Control plots
  if (doJES)
    for (unsigned int i=0; i<signalplots.size(); i++)
	signalplots.at(i).Write();
  if (doFitError)
    for (unsigned int i=0; i<backgroundplots.size(); i++)
	backgroundplots.at(i).Write();

  // Nuisance parameter posteriors
  for (unsigned int i=0; i<systematicposteriors.size(); i++) {
	systematicposteriors.at(i)->Write();
  }

  // JES input plots
  if (doFitFunctionChoiceError)
    for (unsigned int i=0; i<extraHistograms.size(); i++)
	extraHistograms.at(i).Write();

  if (useFile || useFunc)
    resolutiongraph->Write("resolution");

  outfile->Close();
  std::cout << "Closed outfile" << std::endl;
  infile->Close();
  std::cout << "Closed infile" << std::endl;
 
  totaltime.Stop();
  std::cout << "Ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

  delete m;
  delete summary;
 
  delete outfile;
  delete infile;

  return 0;

}

