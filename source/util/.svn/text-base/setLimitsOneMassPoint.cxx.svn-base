// ***************************************************************
// This file was created using the CreateProject.sh script
// for project MjjProject.
// CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://www.mppmu.mpg.de/bat
// ***************************************************************

#include "TH1.h"
#include <iostream>
#include "TError.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TVectorD.h"
#include "TStopwatch.h"
#include "TEnv.h"

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
#include <BAT/BCParameter.h>

#include "Bayesian/MjjBATModel.h"
#include "Bayesian/MjjBATProcess.h"
#include "Bayesian/MjjBATShapeChangingSyst.h"
#include "Bayesian/MjjBATScaleChangingSyst.h"
#include "Bayesian/MjjBATAnalysisFacility.h"


int main(int argc,char **argv)
{

	////////////////////////////////////////////////////////////
	// Initialisation of values from arguments, configuration

	// Start counting time
	TStopwatch totaltime;
	totaltime.Start();

	double trimValue = 0.95;
	double CLpercentage = 0.95;

	// Start reading input configuration
	TString configFile;
	TString outputFileName = "";
	double thisMass = -1;
	unsigned int seed = 1234;
    bool useSeed = false;
	float PDFAccErr = 0.;
	bool doPDFAccError = false;
	int ip=1;
	while (ip<argc) {

		if (string(argv[ip]).substr(0,2)=="--") {

			//config file
			if (string(argv[ip])=="--config") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					configFile=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno config file inserted"<<std::endl; break;}
			}

			//mass
			else if (string(argv[ip])=="--mass") {
				if (ip+1<argc) {
					sscanf(argv[ip+1],"%lf",&thisMass);
					ip+=2;
					std::cout << "Reading mass: " << thisMass << std::endl;
				} else {std::cout<<"\nno mass inserted"<<std::endl; break;}
			}

			// override output name from config file, if doing slices for example
			else if (string(argv[ip])=="--outfile") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					outputFileName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno output file specified"<<std::endl; break;}
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

			// acceptance uncertainty due to PDF
			else if (string(argv[ip]) == "--PDFAccErr") {
				if (ip+1<argc) {
					sscanf(argv[ip+1],"%f",&PDFAccErr);
					ip+=2;
					std::cout << "PDF acceptance uncertainty: " << PDFAccErr << std::endl;
					doPDFAccError = true;
				} else {std::cout<<"\nno PDF acceptance uncertainty specified"<<std::endl; break;}
			}

			//unknown command
			else {
				std::cout<<"\nsetLimitsOneMassPoint: command '"<<string(argv[ip])<<"' unknown"<<std::endl;
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") ip+=2;
				else ip+=1;
			} }//end if "--command"

		else { //if command does not start with "--"
			std::cout << "\nsetLimitsOneMassPoint: command '"<<string(argv[ip])<<"' unknown"<<std::endl;
			break;
		}//end if "--"

	}//end while loop

	// Get config file
	TEnv * settings = new TEnv();
	int status = settings->ReadFile(configFile.Data(),EEnvLevel(0));
	if (status!=0) {
		std::cout<<"cannot read config file"<<std::endl;
		std::cout<<"******************************************\n"<<std::endl;
		return 1;
	}

	if (thisMass < 0 ) {
		std::cout << "Please specify a mass! Use --mass " << std::endl;
		return 1;
	}

	//------------------------------------------
	// Data 
	TString inputFileName = settings->GetValue("dataFileName","");
	TString dataHist = settings->GetValue("dataHist","");

	//------------------------------------------
    // Fitting options

    // Swift fit or regular fit?
    bool doSwift = settings->GetValue("doSwift",false);

	// Range for data fit
	double minX = settings->GetValue("minXForFit",-1.0);
	double maxX = settings->GetValue("maxXForFit",-1.0);

	// Control for choosing nominal & alternate functions.
	int functionNumber = settings->GetValue("functionCode",9);
	int nPars = settings->GetValue("nParameters",3);
	int alternateFuncNumber = 0, altNPars = 0;

	// Check whether any parameters in fit are to be fixed
	vector<bool> areParamsFixed;
	for (int par = 1; par<=nPars; par++) {
		string title = Form("fixParameter%i",par);
		bool isFixed = settings->GetValue(title.c_str(),false);
		areParamsFixed.push_back(isFixed);
	}

	// Parameters for alternate fit, if requested
	// Note it may not have been run in search phase, so
	// do not count on picking up defaults from the input file.
	bool doFitFunctionChoiceError = settings->GetValue("doFitFunctionChoiceError",false);
    int nFitsInFcnError = settings->GetValue("nFitsInFcnError",100);
	double nFitFSigmas = settings->GetValue("nFitFSigmas",1);

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

	// Do we use extended region for fit?
	bool doExtendedRange = settings->GetValue("doExtendedRange",true);

	// Center of mass energy
	double Ecm = settings->GetValue("Ecm",8000.0);

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

	//------------------------------------------
  // Limit setting parameters

	// Signal 
	TString inputSignalFileName = Form(settings->GetValue("signalFileName",""), int(thisMass));
	TString nominalSignalHist = Form(settings->GetValue("nominalSignalHist",""), int(thisMass)); // Lydia changed to now 'form' mass with nominalSignalHist, due to name change convention
	string signame = settings->GetValue("signame","");

	// Output 
	if (outputFileName.CompareTo("")==0) {
		std::cout << "Output filename was: " << outputFileName << std::endl;
		outputFileName = Form(settings->GetValue("outputFileName",""), int(thisMass));
		std::cout << "Now going to be: " << outputFileName << std::endl;
	}

	if (useSeed) {
     TString addon = Form("_seed%d.root",int(seed));
     outputFileName.ReplaceAll(".root",addon);
   }

	string plotDirectory = settings->GetValue("plotDirectory","");
	string plotExtension = settings->GetValue("plotNameExtension","");

	// Do we calculate expected limits?
	bool doExpected = settings->GetValue("doExpected",false);
	// If so, how many?
	int nPseudoExpForExpected = settings->GetValue("nPEForExpected",5);

	// How many standard deviations do we cover?
	double nSigmas = settings->GetValue("nSigmas",0.);

	// Do we do bkg error?
	bool doFitError = settings->GetValue("doFitError",false);
	int nFitsInBkgError = settings->GetValue("nFitsInBkgError",100);

	// Do we do lumi?
	bool doLumiError = settings->GetValue("doLumiError", false);
	double luminosityErr = settings->GetValue("luminosityErr",0.);

	// Do we do JES?
	bool doJES = settings->GetValue("doJES",false);
	int nJES = settings->GetValue("nJES",0);

	// Do we do Beam choice?
	bool doBeam = settings->GetValue("doBeam",false);
	string BeamFileName = settings->GetValue("BeamFile","");

	// Are we using templates or matrices?
	bool useMatrices = settings->GetValue("useMatrices",true);
	bool useTemplates = settings->GetValue("useTemplates",false);
	if (useMatrices && useTemplates) {
		std::cout << "This takes matrices or templates, but not both!" << std::endl;
		return 1;
	} else if (doJES && !(useMatrices || useTemplates)) {
		std::cout << "In order to convolute the JES uncertainty,"
			<< " either templates or matrices are required!" << std::endl;
		return 1;
	}

	// Read and store JES component names
	int nComp = 0;
	string nominalJES;
	vector<string> JESComponentNames;
	vector<vector<string> > JESVariationNames;
	if (useMatrices) {
		nComp = settings->GetValue("nComponents",-1);
		nominalJES = Form(settings->GetValue("nominalJES",""),(int)thisMass);// Lydia added Form for 13 TeV

	} else if (useTemplates) {
		nComp = settings->GetValue("nComponentsTemp",-1); // Lydia added for multi template
		nominalJES = Form(settings->GetValue("nominalTemplateJES",""),(int)thisMass);// Lydia added Form for 13 TeV

	}

	for (int comp=0; comp<nComp; comp++) {
		string getname = Form("name%d",comp+1);
		string gettempname = Form("nameTemp%d",comp+1); //Lydia added for multi templates
		string componentname;
		if (useMatrices) componentname = Form(settings->GetValue(getname.c_str(),""),(int)thisMass); // Lydia added Form for 13 TeV
		else componentname = Form(settings->GetValue(gettempname.c_str(),""),(int)thisMass); // Lydia added Form for 13 TeV
		JESComponentNames.push_back(componentname);
		std::cout<<"COMP "<<componentname<<std::endl;
		vector<string> componentnames;
		for (int i=1; i<nJES; i++) {
			string extension = Form("extension%d",i);
			string JESextension = settings->GetValue(extension.c_str(),"");
                        std::cout << "JES extension is " << JESextension << std::endl;
			// Lydia string variedJESName = "matrix_mjj_"+componentname+"_"+JESextension;
			// Lydia 8 TeV string variedJESName = componentname+"_"+JESextension;
			string variedJESName = componentname+JESextension;
			std::cout<<"JNAME "<<variedJESName<<std::endl;
			componentnames.push_back(variedJESName);
		}
		JESVariationNames.push_back(componentnames);
	}

	// Open files
	std::cout << "Opening files " << inputFileName << ", " << inputSignalFileName << std::endl;
	TFile * infile = TFile::Open(inputFileName,"READ");
	assert(infile);
	TFile * insignalfile = TFile::Open(inputSignalFileName,"READ");
	TH1::AddDirectory(kFALSE);

	// Get and store histograms
	TH1D * basicInputHisto = (TH1D*) infile->Get(dataHist);
	std::cout << "Data hist " << dataHist << " from " << inputFileName << std::endl;
	assert(basicInputHisto);
	MjjHistogram theHistogram((TH1D*)basicInputHisto);
    std::cout << "Looking for hist " << nominalSignalHist << " in file " << inputSignalFileName << std::endl;
	TH1D * thisSigHisto = (TH1D*) insignalfile->Get(nominalSignalHist);
	assert(thisSigHisto);
    std::cout << "Printing hist " << nominalSignalHist << " from file " << inputSignalFileName<< std::endl;
//    thisSigHisto->Print("all");

 	std::vector<TH1D> extraHistograms;

	////////////////////////////////////////////////////////////
	// Any re-sets of values from search phase

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

	TStopwatch fittime;
	fittime.Start();

	vector<double> fittedParametersFromSearchPhase;
	TVectorD *fittedParameters = (TVectorD*) infile->Get("fittedParameters");
	fittedParameters->Print();
	for (int par = 0; par < nPars; par++) {
		fittedParametersFromSearchPhase.push_back((*fittedParameters)[par]);
		std::cout << "Added par "<< par << " with value "<< (*fittedParameters)[par] << std::endl;
	}

	////////////////////////////////////////////////////////////
	// Create fit functions and get background from simultaneous fit

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

	//	ThreeParam2015FitFunction theMjjFitFunction(minXForFit,maxXForFit,Ecm);
	// Use old fitted pars as defaults for this fit.
	theMjjFitFunction->SetParameterDefaults(fittedParametersFromSearchPhase);
	theMjjFitFunction->RestoreParameterDefaults();
	for (int par = 0; par<nPars; par++) {
		std::cout << "Fixing parameter " << par << " to " << areParamsFixed.at(par) << std::endl;
		theMjjFitFunction->GetParameter(par)->SetFixParameter(areParamsFixed.at(par));
	}
	if (doFitFunctionChoiceError) {
	  // Create four-param fit function
          // This one really needs you to tell it what the defaults should be.
	  theAlternateFunction->SetParameterDefaults(altParDefaults);
	  theAlternateFunction->RestoreParameterDefaults();
	}
	std::cout << "Starting initial fit." << std::endl;

  	// This signal histogram isn't used for the limit setting, so we can trim it now. Makes fit more stable...
	double thisPercentage=0;
	double thisInterval=0;
	int keeplowbin=1;
	int keephighbin=1;
	double smallestInterval= thisSigHisto->GetBinLowEdge(thisSigHisto->GetNbinsX()) + 
		thisSigHisto->GetBinWidth(thisSigHisto->GetNbinsX()) - thisSigHisto->GetBinLowEdge(0) + 1e12;
	//double rememberThisPercentage=0;
	for (int bin1=1; bin1<=thisSigHisto->GetNbinsX()+1; bin1++) { 
		for (int bin2=bin1; bin2<=thisSigHisto->GetNbinsX()+1; bin2++) {
			//std::cout << thisSigHisto->GetNbinsX() << " "<< bin1 << " " << bin2 << std::endl;
			thisPercentage = thisSigHisto->Integral(bin1,bin2)/thisSigHisto->Integral();
			thisInterval = thisSigHisto->GetBinLowEdge(bin2) + thisSigHisto->GetBinWidth(bin2) - thisSigHisto->GetBinLowEdge(bin1);
			if (thisPercentage >= trimValue) {
				if (thisInterval < smallestInterval) {
					keeplowbin=bin1;
					keephighbin=bin2;
					smallestInterval=thisInterval;
					//rememberThisPercentage=thisPercentage;
				}
				break;
			}
		} // end of inner loop
	} // end of outer loop */
    std::cout << "keeplowbin, keephighbin are " << keeplowbin << ", " << keephighbin << std::endl;
	for (int bin = 0; bin<thisSigHisto->GetNbinsX()+2; bin++) {
		if ((bin < keeplowbin || bin > keephighbin)&& thisSigHisto->GetBinContent(bin)>0)  {
            //std::cout << "Setting bin " << bin << " to zero. " << std::endl;
			thisSigHisto->SetBinContent(bin,0.0);
			thisSigHisto->SetBinError(bin,0.0);
		}
	}

  	// Fit with this and retrieve background.
	// Note: do in fact need uncertainties here, so go with 100 fits
	MjjFitter theSilentFitter;
  
    // Set up for Swift if necessary
    if (doSwift) {
      theSilentFitter.SetSwiftParameters(minXForFit,maxXForFit,minXAvailable,maxXAvailable,nBinsWindowLeft,nBinsWindowRight,fixWidthAtLowEnd,fixWidthAtHighEnd,truncateHigh);
    }
  
    // Make sure signal used in running.
	theSilentFitter.SetSignalTemplate(thisSigHisto);
 
	int nFits = 1;
  MjjHistogram theBackgroundFromFit;
	if (doFitError) {
      nFits = nFitsInBkgError;
	  std::cout << "Doing " << nFits << " fits for error." << std::endl;
	  theBackgroundFromFit = theSilentFitter.FitAndGetBkgWithDataErr(*theMjjFitFunction,theHistogram,nFitsInBkgError);
    } else {
      theBackgroundFromFit = theSilentFitter.FitAndGetBkgWithNoErr(*theMjjFitFunction,theHistogram);
    }
	TH1D bkgTemplate(theBackgroundFromFit.GetHistogram());
	bkgTemplate.SetName("BkgTemplate");

	fittime.Stop();
	std::cout << "Finished initial fit in " << fittime.CpuTime() << " seconds." << std::endl;

	// Set up vectors of fitted parameters and errors.
	vector<double> fittedPars, errorPars;
	for (int i=0; i<theMjjFitFunction->GetNParams(); i++) {
		std::cout << "par, error: " << theMjjFitFunction->GetFitFunction()->GetParameter(i) 
			<< " " << theMjjFitFunction->GetFitFunction()->GetParError(i) << std::endl;
		fittedPars.push_back(theMjjFitFunction->GetFitFunction()->GetParameter(i));
		errorPars.push_back(theMjjFitFunction->GetFitFunction()->GetParError(i));
	}

	// Set up variation template for background
	TH1D bkgVariation(bkgTemplate);
	bkgVariation.SetName("BkgVariation");
	for (int i=0; i<bkgTemplate.GetNbinsX()+2; i++) {
		bkgVariation.SetBinContent(i,bkgTemplate.GetBinError(i));
		bkgVariation.SetBinError(i,0.);
	}


        // Old method: Use difference between nominal and alternate in this data.
        // Set up fit function choice templates for background
/*	MjjHistogram bkgFromAlternate;
	TH1D variedfit;
	vector<std::pair<double,TH1D*> > fitvariedtemplates;

	if (doFitFunctionChoiceError) {
		bkgFromAlternate = theSilentFitter.FitAndGetBkgWithDataErr(*theAlternateFunction,theHistogram,100);
		//		bkgFromAlternate = theSilentFitter.FitAndGetBkgWithMCErr(*theAlternateFunction,theHistogram);
		for (int i=0; i<theAlternateFunction->GetNParams(); i++) {
			std::cout << "par, error: " << theAlternateFunction->GetFitFunction()->GetParameter(i) 
				<< " " << theAlternateFunction->GetFitFunction()->GetParError(i) << std::endl;
			fittedPars.push_back(theAlternateFunction->GetFitFunction()->GetParameter(i));
			errorPars.push_back(theAlternateFunction->GetFitFunction()->GetParError(i));
		}
		variedfit = bkgFromAlternate.GetHistogram();
		variedfit.SetName("fitplusonesigmachoice");
		fitvariedtemplates.push_back(std::make_pair(0,&bkgTemplate));
		fitvariedtemplates.push_back(std::make_pair(1,&variedfit));
	}
*/

  // New method: Use average difference between the two functions across a range of pseudoexperiments
  // where the pseudoexperiments are thrown from data.
  MjjHistogram bkgFromAlternate;
  MjjHistogram bkgWithFuncChoiceErr;
  TH1D fitChoicePlus1Sig(bkgTemplate);
  if (doFitFunctionChoiceError) {
    std::cout << "******\n* doFitFunctionChoiceError\n******" << std::endl; //debug
    bkgFromAlternate = theSilentFitter.FitAndGetBkgWithNoErr(*theAlternateFunction,theHistogram);
     std::cout << "* First fit complete" << std::endl; //debug
    vector<double> altfittedParameters = theAlternateFunction->GetCurrentParameterValues();

    // Use results of fit to data to set start parameters for alternate function in fit to PEs
    theAlternateFunction->SetParameterDefaults(altfittedParameters);

        	// Get histogram with errors equal to RMS of distance between functions
        	// across PEs in each bin. This is the uncertainty that will be used for calculating
        	// the search phase p-value using systematics.
 
        	bkgWithFuncChoiceErr = theSilentFitter.FitAndGetBkgWithFitDiffErr(*theMjjFitFunction,*theAlternateFunction,theHistogram,minXForFit,maxXForFit, nFitsInFcnError, true);

           // Create histogram of format which will be used as alternate function
           // to define function choice uncertainty when directionality matters:
          // that is, in the limit setting stage.
          fitChoicePlus1Sig.SetName("alternateFitChoiceSyst_setByRMSTimesDirectionality");
          for (int bin = 1; bin < fitChoicePlus1Sig.GetNbinsX()+1; bin++) {
             double nomquantity = theBackgroundFromFit.GetHistogram().GetBinContent(bin);
             double varquantity = bkgFromAlternate.GetHistogram().GetBinContent(bin);
             if (varquantity < nomquantity) fitChoicePlus1Sig.SetBinContent(bin,nomquantity - bkgWithFuncChoiceErr.GetHistogram().GetBinError(bin));
             else fitChoicePlus1Sig.SetBinContent(bin,nomquantity + bkgWithFuncChoiceErr.GetHistogram().GetBinError(bin));
          }
          extraHistograms.push_back(bkgWithFuncChoiceErr.GetHistogram());
          extraHistograms.at(extraHistograms.size()-1).SetName("nomWithSymmetricFuncChoiceErr");
          extraHistograms.push_back(fitChoicePlus1Sig);
          extraHistograms.at(extraHistograms.size()-1).SetName("RMSScaleAlternateFunctionForErr");
  }
      	vector<std::pair<double,TH1D*> > fitvariedtemplates;
      	fitvariedtemplates.push_back(std::make_pair(0,&bkgTemplate));
      	fitvariedtemplates.push_back(std::make_pair(1,&fitChoicePlus1Sig));

	// Set up Beam choice templates
	TString centralhistname, uphistname, downhistname;
	vector<TH1D> saveHists;
	vector<std::pair<double,TH1D*> > Beamvariedtemplates;
	if (doBeam) {
		TFile BeamFile(BeamFileName.c_str(),"READ");
		TH1::AddDirectory(kFALSE);
		TString centralhistname = Form("NOMINAL_%s_%d",signame.c_str(),(int)thisMass);
		TString uphistname = Form("BEAM_UP_%s_%d",signame.c_str(),(int)thisMass);
		TString downhistname = Form("BEAM_DOWN_%s_%d",signame.c_str(),(int)thisMass);
		std::cout << "About to get " << centralhistname << std::endl;
		TH1D centralBeam(*(TH1D*) BeamFile.Get(centralhistname));
		TH1D plus3Beam(*(TH1D*) BeamFile.Get(uphistname));
		plus3Beam.SetName(uphistname);
		TH1D minus3Beam(*(TH1D*) BeamFile.Get(downhistname));
		minus3Beam.SetName(downhistname);
		saveHists.push_back(centralBeam);
		saveHists.push_back(plus3Beam);
		saveHists.push_back(minus3Beam);
		BeamFile.Close();
		Beamvariedtemplates.push_back(std::make_pair(0.0,&saveHists.at(0)));
		Beamvariedtemplates.push_back(std::make_pair(3.0,&saveHists.at(1)));
		Beamvariedtemplates.push_back(std::make_pair(-3.0,&saveHists.at(2)));
	}

	///////////////////////////////////////////////////////////////////////////
	// Create signal template and variation 
	// Note template is blank; all is in the variation

	std::cout << "Creating signal templates" << std::endl;

	TH1D nominalSignal(*(TH1D*)insignalfile->Get(nominalSignalHist));
	nominalSignal.SetName("nominalSignal");
	TH1D variationSignal(*(TH1D*) insignalfile->Get(nominalSignalHist));
	variationSignal.SetName("variationSignal");
	double totalBinContent = 0;
	for (int i=0; i<nominalSignal.GetNbinsX()+2; i++) {

		if (i>0 && i<nominalSignal.GetNbinsX()+1) totalBinContent+=nominalSignal.GetBinContent(i);

		nominalSignal.SetBinContent(i,0.);
		nominalSignal.SetBinError(i,0.);
	}
	for (int i=0; i<variationSignal.GetNbinsX()+2; i++) {

		if (variationSignal.GetBinContent(i)==0) continue;

		double cont = variationSignal.GetBinContent(i);
		variationSignal.SetBinContent(i,cont/totalBinContent);
		variationSignal.SetBinError(i,0.);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recalculate range for process
	int extendedLastBin = lastBin;
	if (doExtendedRange) {
		for (int i=0; i<theHistogram.GetHistogram().GetNbinsX()+2; i++) {
			//      std::cout << "In bin " << i << " basic content, signal content are " << basicInputHisto->GetBinContent(i) << " " << thisSigHisto->GetBinContent(i) << std::endl;
			if (basicInputHisto->GetBinContent(i) > 0 ||
					thisSigHisto->GetBinContent(i) > 0)
				extendedLastBin = i;
		}
	}

	extendedLastBin = lastBin; // DEBUG EMMA


	std::cout << "Comparison in limit setting will be between bins " << firstBin <<
		", " << extendedLastBin << std::endl;
	std::cout << "This covers mjj range " << basicInputHisto->GetBinLowEdge(firstBin) <<
		", " << basicInputHisto->GetBinLowEdge(extendedLastBin+1) << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Estimate appropriate signal range
	double maxNSignal = 0;
	double stepSize = 10;  // increase nSignalEvents nStepSize events at a time
	double maxLsoFar = 0;
	double thisL = 1;
	MjjLogLikelihoodTest theLogLTest;
	TH1D sumHist(bkgTemplate);
	sumHist.SetName("HistSigEstimate");
    TH1D * testSigHisto = (TH1D*) insignalfile->Get(nominalSignalHist);
	double integral = testSigHisto->Integral();
	testSigHisto->Scale(1./integral);
	TH1D weighthist = theHistogram.GetWeightsHistogram();
    std::vector<std::pair<double,double> > probValues;
	while (maxLsoFar / thisL < 1e5) {

		for (int bin=0; bin<sumHist.GetNbinsX()+2; bin++) {
			sumHist.SetBinContent(bin,bkgTemplate.GetBinContent(bin) + testSigHisto->GetBinContent(bin)*maxNSignal);
			sumHist.SetBinError(bin,0.);
		}
		MjjHistogram theSignalTemplate(&sumHist);
		theSignalTemplate.SetEffectiveFromBasicAndWeights(&weighthist);

		double logL = theLogLTest.DoTest(theHistogram,theSignalTemplate,firstBin,extendedLastBin);
		thisL = exp(-logL);

		if (std::isnan(thisL)) break;
		if (thisL > maxLsoFar) maxLsoFar = thisL;
        probValues.push_back(std::make_pair(maxNSignal,thisL));
      
		maxNSignal+=stepSize;
	}
    TH1D byHandLikelihoods("byHandLikelihoods","byHandLikelihoods",probValues.size(),0,stepSize*probValues.size());
    for (unsigned int index=0; index < probValues.size(); index++) {
      byHandLikelihoods.SetBinContent(index,probValues.at(index).second);
    }
    byHandLikelihoods.Scale(1.0/byHandLikelihoods.Integral());
    extraHistograms.push_back(byHandLikelihoods);

	std::cout << "Estimated maximum signal range = " << maxNSignal << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Create tools for limit calculation.

	std::cout << "Creating BAT tools. " << std::endl;

	// set nicer style for drawing than the ROOT default
	BCAux::SetStyle();

	// open log file
 	TString logfile(Form((plotDirectory+"log_"+plotExtension.c_str()+"_%d_%d.txt").c_str(), (int) thisMass, seed));
        const char* logfilestr = logfile.Data();
	BCLog::OpenLog(logfilestr);
	BCLog::SetLogLevel(BCLog::detail);

	// create new MjjBATModel object
	MjjBATModel * m = new MjjBATModel();

	BCLog::OutSummary("Test model created");

	// create a new summary tool object
	BCSummaryTool * summary = new BCSummaryTool(m);

	///////////////////////////////////////////////////////////////////////////
	// Set data

	std::cout << "Setting data histogram. " << std::endl;

	m->SetData(theHistogram,firstBin,extendedLastBin);

	///////////////////////////////////////////////////////////////////////////
	// Set processes: background, signal

	std::cout << "Creating process BKG" << std::endl;
	if (doFitError) {
      std::cout << "Making background with uncertainty" << std::endl;
      m->AddProcess("BKG",true,-nSigmas,nSigmas);
	} else {
      m->AddProcess("BKG",false,0.,0.);
    }
	m->SetTemplate("BKG",bkgTemplate,bkgVariation);
	vector<string> newparams = m->GetProcess(m->GetProcessIndex("BKG"))->GetParamNames();
	for (unsigned int i=0; i<newparams.size(); i++) {
		m->SetPriorGauss(newparams.at(i).c_str(),0.,1.);
	}

	if (doFitFunctionChoiceError) {

		m->AddSystematic("FUNCCHOICE",0,nFitFSigmas);
		m->SetPriorGauss("FUNCCHOICE",0.,1.);
		MjjBATTemplateSyst * FuncChoiceOnBkg = new MjjBATTemplateSyst(true);
		FuncChoiceOnBkg->SetSpectra(fitvariedtemplates);
		m->SetSystematicVariation("BKG","FUNCCHOICE",FuncChoiceOnBkg);

	}

	std::cout <<"Creating process SIGNAL with range 0, "<< maxNSignal << std::endl;
	m->AddProcess("SIGNAL",true,0,maxNSignal);
	m->SetTemplate("SIGNAL",nominalSignal,variationSignal);
	vector<string> sigparams = m->GetProcess(m->GetProcessIndex("SIGNAL"))->GetParamNames();
	for (unsigned int i=0; i<sigparams.size(); i++) {
		m->SetPriorConstant(sigparams.at(i).c_str());
	}
	// NOT USED!!! m->GetProcess(m->GetProcessIndex("SIGNAL"))->SetTrimPercentage(trimValue);

	///////////////////////////////////////////////////////////////////////////
	// Create shape changing nuisance parameter

	std::cout << "Adding shape changing uncertainties " << std::endl;

	// Step through input JES shifted histograms
	// step size uses nSigmas to calculate step manually, if fails then use manual approach, rather than nSigmas
	// Example of manual approach shown commented out above
    vector<double> jesSigmas(nJES,0.0);
	double step = (2*nSigmas)/(nJES-1); // Fencepost problems
    // Nominal is first, but only do if it there is more than one JES parameter (otherwise the vector is size zero)
    if (nJES>0) jesSigmas[0] = 0.;
	// Below and above nominal
	for (int i=0; i<int(nJES/2.); i++) {
		jesSigmas[i+1] = step*i-nSigmas;
		std::cout<<"JESSIG "<< i+1 <<": "<<jesSigmas[i+1]<<std::endl;
		jesSigmas[nJES-i-1] = nSigmas-step*i;
		std::cout<<"JESSIG "<< nJES-i-1 <<": "<<jesSigmas[nJES-i-1]<<std::endl;
	}

	vector<vector<TH2D> > storeAll2DHistograms;
	vector<vector<TH1D> > storeAll1DHistograms;// Lydia
	
	//doJES = false; //HACK by EMMA
	if (doJES) {
		if (useMatrices) {
			for (int component=0; component<nComp; component++) {

				vector<std::pair<double,TH2D*> > signalJESVariations; signalJESVariations.clear();
				std::vector<TH2D> storeHistograms;

				TH2D nominal(*(TH2D*)insignalfile->Get(nominalJES.c_str()));

				string nomname = nominalJES + "_nom";
				nominal.SetName(nomname.c_str());
				for (int i=0; i<basicInputHisto->GetNbinsX()+2; i++) {
					for (int j=0; j<basicInputHisto->GetNbinsX()+2; j++) {
						if (i==j) nominal.SetBinContent(i,j,1.0);
						else nominal.SetBinContent(i,j,0.0);
					}
				}
				storeHistograms.push_back(nominal);
				for (int i=1; i<nJES; i++) {
					string histname = JESVariationNames.at(component).at(i-1);
					std::cout << "Getting " << histname << std::endl;
					TH2D thisVariation(*(TH2D*)insignalfile->Get(histname.c_str()));
					thisVariation.SetName(histname.c_str());
					storeHistograms.push_back(thisVariation);
				}
				storeAll2DHistograms.push_back(storeHistograms);

				for (int i=0; i<nJES; i++) {
					signalJESVariations.push_back(std::make_pair(jesSigmas[i],&storeAll2DHistograms.at(component).at(i)));
				}

				string name = JESComponentNames.at(component);
				if (doJES) {
					std::cout << "Creating systematic " << name << std::endl;
					m->AddSystematic(name.c_str(),-nSigmas,nSigmas);
					m->SetPriorGauss(name.c_str(),0.,1.);
					MjjBATShapeChangingSyst * thisShapeChangingSyst = new MjjBATShapeChangingSyst();
					thisShapeChangingSyst->SetSpectra(signalJESVariations);
					m->SetSystematicVariation("SIGNAL",name.c_str(),thisShapeChangingSyst);
				}
			} // End of loop over JES components

		} else if (useTemplates ) {
			std::cout << "Using template JES"<< std::endl;
			for (int component=0; component<nComp; component++) { // Lydia
				vector<std::pair<double,TH1D*> > signalJESVariations; signalJESVariations.clear();
				std::vector<TH1D> storeHistograms;// Lydia
				TH1D nominal(*(TH1D*)insignalfile->Get(nominalJES.c_str()));
				string nomname = nominalJES + "_nom";
				nominal.SetName(nomname.c_str());
				storeHistograms.push_back(nominal); // Lydia

				for (int i=1; i<nJES; i++) {
					string histname = JESVariationNames.at(component).at(i-1);
          if (!insignalfile->Get(histname.c_str())){
            std::cout<<"skipping "<<histname<<std::endl;
            continue;
          }
          std::cout<<"Found syst hist "<<histname<<std::endl;
					TH1D thisVariation(*(TH1D*)insignalfile->Get(histname.c_str()));
					thisVariation.SetName(histname.c_str());
					storeHistograms.push_back(thisVariation); //Lydia
				}

				storeAll1DHistograms.push_back(storeHistograms); // Lydia

				for (int i=0; i<nJES; i++) {
					signalJESVariations.push_back(std::make_pair(jesSigmas[i],&storeAll1DHistograms.at(component).at(i)));
				}

				string name = JESComponentNames.at(component);
				if (doJES) {
					std::cout << "Creating systematic " << name << std::endl;
					m->AddSystematic(name.c_str(),-nSigmas,nSigmas);
					m->SetPriorGauss(name.c_str(),0.,1.);
					MjjBATTemplateSyst * thisTemplateSyst = new MjjBATTemplateSyst(true);
					thisTemplateSyst->SetSpectra(signalJESVariations);
					m->SetSystematicVariation("SIGNAL",name.c_str(),thisTemplateSyst);
				}

			} // End of loop over JES components Lydia
		}
	}

	///////////////////////////////////////////////////////////////////////////
	// Create scale changing nuisance parameters

	std::cout << "Adding scale changing uncertainties " << std::endl;

	if (doLumiError) {

		// Create luminosity uncertainty. Add variations in processes.
		m->AddSystematic("LUMI",-nSigmas,nSigmas);
		m->SetPriorGauss("LUMI",0.,1.);
		MjjBATScaleChangingSyst * LumiOnSignal = new MjjBATScaleChangingSyst(luminosityErr);
		LumiOnSignal->SetScaleFromBinContent();
		m->SetSystematicVariation("SIGNAL","LUMI",LumiOnSignal);

	}

	if (doPDFAccError) {

		// Create luminosity uncertainty. Add variations in processes.
		m->AddSystematic("PDFAcc",-nSigmas,nSigmas);
		m->SetPriorGauss("PDFAcc",0.,1.);
		MjjBATScaleChangingSyst * PDFAccOnSignal = new MjjBATScaleChangingSyst(PDFAccErr);
		PDFAccOnSignal->SetScaleFromBinContent();
		m->SetSystematicVariation("SIGNAL","PDFAcc",PDFAccOnSignal);

	}

	std::cout << "doBeam is " << doBeam << std::endl;
	if (doBeam) {
		m->AddSystematic("Beam",-nSigmas,nSigmas);
		m->SetPriorGauss("Beam",0.,1.);
		MjjBATTemplateSyst * BeamEffects = new MjjBATTemplateSyst(true);
		BeamEffects->SetSpectra(Beamvariedtemplates);
		m->SetSystematicVariation("SIGNAL","Beam",BeamEffects);

	}

	///////////////////////////////////////////////////////////////////////////
	// Begin limit setting calculations.

	for (unsigned int i=0; i<m->GetNParameters(); i++) {
		std::cout << "Parameter " << i << " is " << m->GetParameter(i)->GetName() << std::endl;
	}

	std::cout << "Beginning calculation." << std::endl;
	TStopwatch marginalisationTime;
	marginalisationTime.Start();

	// run MCMC and marginalize posterior wrt. all parameters
	// and all combinations of two parameters
	//   m->SetMarginalizationMethod(BCIntegrate::kMargDefault);
//   m->SetMarginalizationMethod(BCIntegrate::kMargMetropolis);
//   m->SetMarginalizationMethod(BCIntegrate::kMargMonteCarlo);
//   m->SetMarginalizationMethod(BCIntegrate::kMargGrid);
	//m->PrintTemplates("battemplates_before.pdf");
	m->MarginalizeAll();
	m->PrintAllMarginalized((outputFileName + '_').ReplaceAll(".root","") + "ModelPlots.pdf");
	//m->PrintTemplates("battemplates_after.pdf");

	marginalisationTime.Stop();

	BCParameter * normPar = m->GetParameter(m->GetProcess(1)->GetNormalisationParamIndex());
	BCH1D * marginalizedSig = m->GetMarginalized(normPar);
	TH1D * sigLikeVsNumber = (TH1D*) marginalizedSig->GetHistogram()->Clone();
	double trueCL = marginalizedSig->GetLimit(CLpercentage);
	std::cout << "CL for this mass is " << trueCL << std::endl;

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

	std::cout << "Making plots" << std::endl;
	// draw all marginalized distributions into a PostScript file
	TString marginalizedPlotFilename(Form((plotDirectory+"MjjBATModel_plots_"+plotExtension.c_str()+"_%d_%d.ps").c_str(), (int) thisMass, seed));

	// print all summary plots
	TString parameterPlotFilename(Form((plotDirectory+"MjjBATModel_parameters_"+plotExtension.c_str()+"_%d_%d.eps").c_str(), (int) thisMass, seed));
	summary->PrintParameterPlot(parameterPlotFilename);
	TString correlationPlotFilename(Form((plotDirectory+"MjjBATModel_correlation_"+plotExtension.c_str()+"_%d_%d.eps").c_str(), (int) thisMass, seed));
	summary->PrintCorrelationPlot(correlationPlotFilename);
	if(doJES && (doFitError || doLumiError)) {
		TString updatePlotFilename(Form((plotDirectory+"MjjBATModel_update_"+plotExtension.c_str()+"_%d_%d.ps").c_str(), (int) thisMass, seed));
		summary->PrintKnowledgeUpdatePlots(updatePlotFilename); 
	}

	// print results of the analysis into a text file
	TString resultPlotFilename(Form((plotDirectory+"MjjBATModel_results_"+plotExtension.c_str()+"_%d_%d.txt").c_str(), (int) thisMass, seed));
	m->PrintResults(resultPlotFilename);

	////////////////////////////////////////////////////////////
	// Create output file in which to write pseudoexperiments
 

	std::cout << "Creating output file" << std::endl;
	TFile * outfile = TFile::Open(outputFileName,"RECREATE");

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

	if ((doJES) && (nComp<5)) {

		// Create parameters for plots
		const int nJESSlices = 5;
		double paramslices [nJESSlices] = {-3.0,-1.0,0.0,1.0,3.0};

		/*double signalslices [10];
		  for (unsigned int i=0; i<10; i++) {
		  signalslices[i] = maxNSignal/10.0 * i;
		  }*/

		// All combinations of possible JES parameter values
		int nCombinations = (int)pow(nJESSlices,nComp);
		double ** paramArray = new double*[nCombinations];
		for (int i=0; i<nCombinations; i++) {
			paramArray[i] = new double[nComp];
			for (int j=0; j<nComp; j++) {
				int nInBlock = (int)pow(nJESSlices,(nComp-j-1));
				int val = ceil(i/nInBlock);
				int index = val%nJESSlices;
				paramArray[i][j] = paramslices[index];
			}
		}

		for (int jesindex=0; jesindex < nCombinations; jesindex++) {
			TH1D thissig(nominalSignal);
			string name = "sigsample";
			parameters.clear();
			parameters.push_back(0.0); //background
			if (doFitFunctionChoiceError) parameters.push_back(0.0); // function choice
			parameters.push_back(1.0); //signal
			for (int npar=0; npar<nComp; npar++) {
				double jesval = paramArray[jesindex][npar];
				parameters.push_back(jesval); // jes component
				name = Form((name+"_jes%d").c_str(),(int)jesval);
			}
			if(doLumiError) parameters.push_back(0.0); //lumi
			if(doPDFAccError) parameters.push_back(0.0); //PDF acceptance
			if (doBeam) parameters.push_back(0.0); // Beam
			thissig.SetName(name.c_str());
			vector<std::pair<double,double> > sigexpectation = m->ProcessExpectation(1,parameters);
			for (int bin=0; bin<thissig.GetNbinsX(); bin++) {
				thissig.SetBinContent(bin,sigexpectation.at(bin).first);
			}
			signalplots.push_back(thissig);
		}

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
			if (doFitFunctionChoiceError) parameters.push_back(0.0); // function choice
			parameters.push_back(0.0); //signal
			for (int npar=0; npar<nComp; npar++) parameters.push_back(0.0); // jes
			if (doLumiError) parameters.push_back(0.0); //lumi
			if(doPDFAccError) parameters.push_back(0.0); //PDF acceptance
			if (doBeam) parameters.push_back(0.0); // Beam
			vector<std::pair<double,double> > bkgexpectation = m->ProcessExpectation(0,parameters);
			for (int i=0; i<thissig.GetNbinsX()+2; i++) {
				thissig.SetBinContent(i,bkgexpectation.at(i).first);
			}
			backgroundplots.push_back(thissig);
		}

	}

	////////////////////////////////////////////////////////////
	// Perform pseudoexperiments


	if (doExpected) {

		std::cout << "*******************\n* PSEUDOEXPERIMENTS\n*******************" <<std::endl;

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

		// Make analysis facility
		MjjBATAnalysisFacility * maf = new MjjBATAnalysisFacility(m,seed);
		maf->SetFlagMCMC(true);

		// Current model params still at best-fit with no signal,
		// so use those for pseudoexperiments
		TTree * pseudotree = maf->BuildEnsembles(useParams, 10000 );

		// Un-fix signal.
        m->GetParameter(sigparam)->Unfix();
      
		TTree * outputtree = maf->PerformEnsembleTest(pseudotree, nPseudoExpForExpected, firstBin, extendedLastBin );

		outfile->cd();
		outputtree->Write();

	}

	////////////////////////////////////////////////////////////
	// Save everything.

	// Print plots of input variations into file for checking
	std::cout << "Writing output files " << outputFileName << std::endl;

	outfile->cd();

	std::cout << "In output file " << std::endl;

	sigLikeVsNumber->SetName("likelihoodFunction");
	sigLikeVsNumber->Write();

	std::cout << "Wrote likelihood function" << std::endl;

	TVectorD CLOfRealLikelihood(1);
	CLOfRealLikelihood[0] = trueCL;
	CLOfRealLikelihood.Write("CLOfRealLikelihood");

	std::cout << "Wrote trueCL" << std::endl;

	// Original signal shape
	TH1D * nomSigHisto = (TH1D*) insignalfile->Get(nominalSignalHist);
	nomSigHisto->Write();

	// Transfer matrices
	if (useMatrices) {
		for (unsigned int i=0; i<storeAll2DHistograms.size(); i++) {
			vector<TH2D>  compHistograms = storeAll2DHistograms.at(i);
			for (unsigned int j=0; j<compHistograms.size(); j++)
				compHistograms.at(j).Write();
		}
	}

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

	// This is storage for any extra plots we wanted
    for (unsigned int i=0; i<extraHistograms.size(); i++)
        extraHistograms.at(i).Write();

	outfile->Close();
	std::cout << "Closed outfile" << std::endl;
	infile->Close();
	std::cout << "Closed infile" << std::endl;
	insignalfile->Close();
	std::cout << "Closed insignalfile" << std::endl;

	delete m;
	delete summary;
        delete infile;
        delete insignalfile;

	std::cout << "Deleted m, summary" << std::endl;

	BCLog::OutSummary("Test program ran successfully");
	BCLog::OutSummary("Exiting");

	// close log file
	BCLog::CloseLog();

	std::cout << "Closed log file" << std::endl;

	totaltime.Stop();
	std::cout << "Marginalisation took " << marginalisationTime.CpuTime() << " seconds." << std::endl;
	std::cout << "Entire process ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

	return 0;

}

