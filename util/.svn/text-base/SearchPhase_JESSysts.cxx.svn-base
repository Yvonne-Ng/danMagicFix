#include "TH1.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "Bayesian/MathFunctions.h"
#include "Bayesian/MjjFitter.h"
#include "Bayesian/MjjHistogram.h"
#include "Bayesian/MjjFitFunction.h"
#include "Bayesian/BonusFitFunctions.h"
#include "Bayesian/MjjStatisticalTest.h"
#include "Bayesian/MjjChi2Test.h"
#include "Bayesian/MjjPseudoExperimenter.h"
#include "Bayesian/MjjStatisticsBundle.h"
#include "Bayesian/MjjSignificanceTests.h"
#include "Bayesian/MjjBumpHunter.h"

#include "Bayesian/MjjBATModel.h"
#include "Bayesian/MjjBATProcess.h"
#include "Bayesian/MjjBATShapeChangingSyst.h"
#include "Bayesian/MjjBATScaleChangingSyst.h"
#include "Bayesian/MjjBATAnalysisFacility.h"

#include <BAT/BCAux.h>
#include <BAT/BCSummaryTool.h>

#include "TEnv.h"
#include "TFile.h"
#include "TString.h"
#include "TVector.h"
#include "TStopwatch.h"

int main (int argc,char **argv)
{

	////////////////////////////////////////////////////////////
	// Initialisation: Read from config file

	// Start counting time
	TStopwatch totaltime;
	totaltime.Start();

	bool f_noDataErr = false;
	bool f_useScaled = false;
	bool f_saveExclusion = false;
    bool f_saveToyFits = false;
	float thresholdPVal_FindSignal = 0.01;
	float thresholdPVal_RemoveSignal = 0.01;
	float minBHMass = -1;
	float maxBHMass = -1;
	TString inputFileName;
	TString outputFileName;
	TString inputHistDir = "";
	TString dataMjjHistoName;

	// Start reading input configuration
	TString configFile;
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

			//in file
			else if (string(argv[ip])=="--file") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					inputFileName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno input file name inserted"<<std::endl; break;}
			}

			//output file
			else if (string(argv[ip])=="--outputfile") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					outputFileName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno output file name inserted"<<std::endl; break;}
			}

			//directory of the input histogram
			else if (string(argv[ip])=="--dir") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					inputHistDir=std::stof(string(argv[ip+1]));
					ip+=2;
				} else {std::cout<<"\nno histogram directory given"<<std::endl; break;}
			}

			//histogram name
			else if (string(argv[ip])=="--histName") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					dataMjjHistoName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno histogram name inserted"<<std::endl; break;}
			}

			//threshold to use
			else if (string(argv[ip])=="--threshold") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					thresholdPVal_RemoveSignal=std::stof(string(argv[ip+1]));
					ip+=2;
				} else {std::cout<<"\nno p-value threshold value given "<<std::endl; break;}
			}

			//minimum mass for Bump Hunter
			else if (string(argv[ip])=="--minBH") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					minBHMass=std::stof(string(argv[ip+1]));
					ip+=2;
				} else {std::cout<<"\nNo BumpHunter minimum value given "<<std::endl; break;}
			}

			//maximum mass for Bump Hunter
			else if (string(argv[ip])=="--maxBH") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					maxBHMass=std::stof(string(argv[ip+1]));
					ip+=2;
				} else {std::cout<<"\nNo BumpHunter maximum value given "<<std::endl; break;}
			}

			//Save exclusion results to text file
			else if (string(argv[ip])=="--saveEx") {
				f_saveExclusion = true;
				ip+=1;
			}

			//Don't Run Data Err
			else if (string(argv[ip])=="--noDE") {
				f_noDataErr = true;
				ip+=1;
			}
          
            // Keep every fitted pseudoexperiment. If false keep only 25.
            else if (string(argv[ip])=="--saveToyFits") {
				f_saveToyFits = true;
				ip+=1;
            }

			//Add this flag to run on Scaled MC rather than Data-like MC or Data 
			else if (string(argv[ip])=="--useScaled") {
				f_useScaled = true;
				ip+=1;
			}

			//unknown command
			else {
				std::cout<<"\nSearchPhase: command '"<<string(argv[ip])<<"' unknown"<<std::endl;
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") ip+=2;
				else ip+=1;
			} }//end if "--command"

		else { //if command does not start with "--"
			std::cout << "\nSearchPhase: command '"<<string(argv[ip])<<"' unknown"<<std::endl;
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

	// Specify files
	if( inputFileName.Length() == 0){
		inputFileName = settings->GetValue("inputFileName","");
	}
	std::cout<<"inputFileName: "<<inputFileName<<std::endl;
	if ( outputFileName.Length() == 0) {
		outputFileName = settings->GetValue("outputFileName","");
	}

	//input histogram directory
	if (inputHistDir.Length() == 0){
		inputHistDir = settings->GetValue("inputHistDir","");
	}
	std::cout<<"inputHistDir: "<<inputHistDir<<std::endl;

	// Get center of mass energy
	double Ecm = settings->GetValue("Ecm",13000.0);

	// For certain tests (ie unphysical MC stats) do not want window exclusion possible.
	bool permitWindow = settings->GetValue("permitWindow",true);

	// Do we use the systematic uncertainties in the p-value calculation?
	bool doPValWithSysts = settings->GetValue("doPValWithSysts",false);
	std::cout << "Doing pvalue with systematics? " << doPValWithSysts << "!" << std::endl;

	// Open files
	TFile * infile = TFile::Open(inputFileName,"READ");
	TH1::AddDirectory(kFALSE);

	int nPseudoExpBH = settings->GetValue("nPseudoExp",1e3);
	std::cout << "nPseudoExpBH is "<< nPseudoExpBH << std::endl;

	int nPseudoExpFit = settings->GetValue("nPseudoExpFit",100);
	std::cout << "nPseudoExpFit is "<< nPseudoExpFit << std::endl;

	// Get and store histograms
	if( dataMjjHistoName.Length() == 0){
		dataMjjHistoName = settings->GetValue("dataHist","");
	}
	std::cout<<" Hist Name: "<<dataMjjHistoName<<std::endl;
	TH1D* basicInputHisto = new TH1D();
	if (inputHistDir.Length() == 0) {
		basicInputHisto = (TH1D*) infile->Get(dataMjjHistoName);
	} else {
		basicInputHisto = (TH1D*) infile->GetDirectory(inputHistDir)->Get(dataMjjHistoName);
	}

	// If this is a scaled histogram, the errors need to be correct
	if( f_useScaled == true ){
                std::cout<<"Using Scaled MC"<<std::endl;
		for( int iBin=1; iBin < basicInputHisto->GetNbinsX()+1; ++iBin){
			basicInputHisto->SetBinError(iBin, sqrt(basicInputHisto->GetBinContent(iBin)) );
		}
	}

	MjjHistogram theHistogram(basicInputHisto);

	gErrorIgnoreLevel=kWarning; // Want default to ignore info messages from minimizers

    // Swift fit or regular fit?
    bool doSwift = settings->GetValue("doSwift",false);

	// Range for data fit
	double minX = settings->GetValue("minXForFit",-1.0);
	double maxX = settings->GetValue("maxXForFit",-1.0);

	int functionNumber = settings->GetValue("functionCode",9);
	int nPars = settings->GetValue("nParameters",3);

	bool doPEOnData = settings->GetValue("doPEOnData",false);

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

	bool doAlternate = settings->GetValue("doAlternateFunction",true);
	int alternateFuncNumber = 0, altNPars = 0;
	// Parameters for alternate fit, if requested
	vector<double> altParDefaults;
	vector<bool> altAreParsFixed;
	if (doAlternate) {
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
    // Uncertainty options!

	// Do we do JES?
	bool doJES = settings->GetValue("doJES",false);
	int nJES = settings->GetValue("nJES",0);

    // Do we fit the toys?
    bool fitToys = settings->GetValue("fitToys",false);

	// How many standard deviations do we cover?
	double nSigmas = settings->GetValue("nSigmas",3.0);
  
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
		nominalJES = settings->GetValue("nominalJES","");

	} else if (useTemplates) {
		nComp = settings->GetValue("nComponentsTemp",-1); // Lydia added for multi template
		nominalJES = settings->GetValue("nominalTemplateJES","");

	}

	for (int comp=0; comp<nComp; comp++) {
		string getname = Form("name%d",comp+1);
		string gettempname = Form("nameTemp%d",comp+1); //Lydia added for multi templates
		string componentname;
		if (useMatrices) componentname = settings->GetValue(getname.c_str(),""); // Lydia added Form for 13 TeV
		else componentname = settings->GetValue(gettempname.c_str(),""); // Lydia added Form for 13 TeV
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
  
	//------------------------------------------
	//initial BH values, before window exclusion
	double initialBHPValue;
	double initialBHPValueErr;
	double initialBHValue;
	double initialBHRangeLow;
	double initialBHRangeHigh;

	//refined BH values, after window exclusion  
	double refinedBHPValue = -1;
	double refinedBHPValueErr = -1;
	double refinedBHValue = -1;
	double refinedBHRangeLow = -1;
	double refinedBHRangeHigh = -1;

	//flag for (significant) bump identification
	Bool_t bumpFound = false;

	//------------------------------------------
	// Set up fit ranges

	std::cout << "minX, maxX: " << minX << " " << maxX << std::endl;

	// Lydia changed this part so picks start and end bins that are fully spanned by range minX-maxX
	// e.g. user picks minX maxX of 2000-3000 only bins fully within this range are fitted
	// rather than fitting e.g. 1950-350

	int firstBin, lastBin;
	int firstBinBH, lastBinBH;
	if (minX < theHistogram.GetHistogram().GetBinLowEdge(theHistogram.GetFirstBinWithData()) || minX < 0) firstBin = theHistogram.GetFirstBinWithData();
	else firstBin = theHistogram.GetHistogram().FindBin(minX)+1;

	if (maxX > theHistogram.GetHistogram().GetBinLowEdge(theHistogram.GetLastBinWithData()+1) || maxX < 0) lastBin = theHistogram.GetLastBinWithData();
	else lastBin = theHistogram.GetHistogram().FindBin(maxX)-1;
    // Below can be different from specified minX, maxX if they weren't bin edges. This is a bin edge
	double minXForFit = theHistogram.GetHistogram().GetBinLowEdge(firstBin);
	double maxXForFit = theHistogram.GetHistogram().GetBinLowEdge(lastBin+1);

	if ( minBHMass == -1)
		firstBinBH = firstBin;
	else
		firstBinBH = minBHMass;
	if ( maxBHMass == -1)
		lastBinBH = lastBin;
	else
		lastBinBH = maxBHMass;

    std::cout << "Actual fit range is " << minXForFit << " - " << maxXForFit << std::endl;
    std::cout << "corresponding to including the bins " << firstBin << ", " << lastBin << "." << std::endl;

	////////////////////////////////////////////////////////////
	// Beginning search phase
  
	// Create fit function(s)
	MjjFitFunction * theMjjFitFunction;
	MjjFitFunction * theAlternateFunction;
	vector<std::pair<int,MjjFitFunction**> > functionsAndCodes;
	functionsAndCodes.push_back(std::make_pair(functionNumber,&theMjjFitFunction));
	if (doAlternate) functionsAndCodes.push_back(std::make_pair(alternateFuncNumber,&theAlternateFunction));

	for (unsigned int index = 0; index < functionsAndCodes.size(); index++) {

		int thisFuncCode = functionsAndCodes.at(index).first;
		switch (thisFuncCode) {
			case 1 :
				std::cout << "Creating UA2 fit function." << std::endl;
				*functionsAndCodes.at(index).second = new UA2FitFunction(minXAvailable,maxXAvailable,Ecm);
				break;
			case 2 :
				std::cout << "Creating CDF (1995) fit function." << std::endl;
				*functionsAndCodes.at(index).second = new CDFFitFunction(minXAvailable,maxXAvailable,Ecm);
				break;
			case 3 :
				std::cout << "Creating CDF (1997) fit function." << std::endl;
				*functionsAndCodes.at(index).second = new CDF1997FitFunction(minXAvailable,maxXAvailable,Ecm);
				break;
			case 4 :
				std::cout << "Creating standard dijet function." << std::endl;
				*functionsAndCodes.at(index).second = new FourParamFitFunction(minXAvailable,maxXAvailable,Ecm);
				break;
			case 5 :
				std::cout << "Creating TeV Gravity function." << std::endl;
				*functionsAndCodes.at(index).second = new ThreeParamFitFunction(minXAvailable,maxXAvailable,Ecm);
				break;
			case 6 :
				std::cout << "Creating 5-parameter function, floating sqrt(s)." << std::endl;
				*functionsAndCodes.at(index).second = new FiveParamSqrtsFitFunction(minXAvailable,maxXAvailable,Ecm);
				break;
			case 7 :
				std::cout << "Creating 5-parameter function, log(x)^2 term." << std::endl;
				*functionsAndCodes.at(index).second = new FiveParamLog2FitFunction(minXAvailable,maxXAvailable,Ecm);
				break;
			case 8 :
				std::cout << "Creating 6-parameter function." << std::endl;
				*functionsAndCodes.at(index).second = new SixParamFitFunction(minXAvailable,maxXAvailable,Ecm);
                                break;
			case 9 :
				std::cout << "Creating 3-parameter function for Run II search." << std::endl;
				*functionsAndCodes.at(index).second = new ThreeParam2015FitFunction(minXAvailable,maxXAvailable,Ecm);
                                break;
            case 10 :
                std::cout << "Creating second Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction2(minXAvailable,maxXAvailable,Ecm);
                break;
            case 11 :
                std::cout << "Creating third Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction3(minXAvailable,maxXAvailable,Ecm);
                break;
            case 12 :
                std::cout << "Creating fourth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction4(minXAvailable,maxXAvailable,Ecm);
                break;
            case 13 :
                std::cout << "Creating fifth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction5(minXAvailable,maxXAvailable,Ecm);
                break;
            case 14 :
                std::cout << "Creating sixth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction6(minXAvailable,maxXAvailable,Ecm);
                break;
            case 15 :
                std::cout << "Creating seventh Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction7(minXAvailable,maxXAvailable,Ecm);
                break;
            case 16 :
                std::cout << "Creating eigth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction8(minXAvailable,maxXAvailable,Ecm);
                break;
            case 17 :
                std::cout << "Creating ninth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction9(minXAvailable,maxXAvailable,Ecm);
                break;
            case 18 :
                std::cout << "Creating tenth Multijet function." << std::endl;
 				*functionsAndCodes.at(index).second = new MultijetFitFunction10(minXAvailable,maxXAvailable,Ecm);
                break;
            case 19 :
                std::cout << "Creating gamma gamma with no exponent, 2 par." << std::endl;
                *functionsAndCodes.at(index).second = new TwoParGammaGamma(minXAvailable,maxXAvailable,Ecm);
                break;
            case 20 :
                std::cout << "Creating gamma gamma with no exponent, 4 par." << std::endl;
                *functionsAndCodes.at(index).second = new FourParGammaGamma(minXAvailable,maxXAvailable,Ecm);
                break;
            case 21 :
                std::cout << "Creating gamma gamma with 1/3 exponent, 2 par." << std::endl;
                *functionsAndCodes.at(index).second = new TwoParGammaGammaWithThird(minXAvailable,maxXAvailable,Ecm);
                break;
            case 22 :
                std::cout << "Creating gamma gamma with 1/3 exponent, 3 par." << std::endl;
                *functionsAndCodes.at(index).second = new ThreeParGammaGammaWithThird(minXAvailable,maxXAvailable,Ecm);
                break;
            case 23 :
                std::cout << "Creating gamma gamma with 1/3 exponent, 4 par." << std::endl;
                *functionsAndCodes.at(index).second = new FourParGammaGammaWithThird(minXAvailable,maxXAvailable,Ecm);
                break;

 		}

	}

	std::cout<<"paramDefaultsLength"<< paramDefaults.size()<<std::endl; // Lydia
	theMjjFitFunction->SetParameterDefaults(paramDefaults);
	theMjjFitFunction->RestoreParameterDefaults();
	for (int par = 0; par<nPars; par++)
		theMjjFitFunction->GetParameter(par)->SetFixParameter(areParamsFixed.at(par));

	if (doAlternate) {
		theAlternateFunction->SetParameterDefaults(altParDefaults);
		theAlternateFunction->RestoreParameterDefaults();
		for (int par = 0; par<nPars; par++)
			theAlternateFunction->GetParameter(par)->SetFixParameter(altAreParsFixed.at(par));
	}

	MjjFitter theFitter;

    // Set up for Swift if necessary
    std::cout << "doSwift is " << doSwift << std::endl;
    if (doSwift) {
      theFitter.SetSwiftParameters(minXForFit,maxXForFit,minXAvailable,maxXAvailable,nBinsWindowLeft,nBinsWindowRight-1,fixWidthAtLowEnd,fixWidthAtHighEnd,truncateHigh);
    }

	// Perform fit and retrieve background histogram
	// MjjHistogram backgroundFromFunc = theFitter.FitAndGetBkgWithDataErr(*theMjjFitFunction,theHistogram,100);
	// Lydia changed -- Kate has added a line after this step with full pseudoexperiments and thus accepts the change.
	MjjHistogram backgroundFromFunc = theFitter.FitAndGetBkgWithNoErr(*theMjjFitFunction,theHistogram,minXForFit,maxXForFit);

	// Save result of fit to use later
	vector<double> fittedParameters = theMjjFitFunction->GetCurrentParameterValues();

	// What were fitted params?
	std::cout << "After fit, parameters are  ";
	for (unsigned int i=0; i<fittedParameters.size(); i++) {
		std::cout << "   (" << i << "):" << fittedParameters.at(i);
	}
	std::cout << std::endl;

	///////////////////////////////////////////////////////////////////////////
	// Create tools for systematic variations.

	std::cout << "Creating BAT tools. " << std::endl;

	// set nicer style for drawing than the ROOT default
	BCAux::SetStyle();

	// create new MjjBATModel object
	MjjBATModel * m = new MjjBATModel();

	BCLog::OutSummary("Test model created");


    // Create a process which will apply variations taken from data
    // to a background estimate taken from the fit.
	std::cout << "Creating process BKG" << std::endl;
    TH1D bkgTemplate = backgroundFromFunc.GetHistogram();
    bkgTemplate.SetName("extendedBkgFromFit_forModelBuilding");
    // In the case that SWIFT is used and the fit range
    // extended beyond the background estimation range, fill those
    // bins too.
    // Just take them from the last SWIFT fit.
//    std::cout << "BEFORE:" << std::endl;
//    bkgTemplate.Print("all");
//    int lastBinFitted = -1;
//    if (doSwift && maxXAvailable > maxX) {
//      TH2D usedBins = theFitter.GetSwiftBinsUsed();
//      int lastBinEstimated = backgroundFromFunc.GetLastBinWithData();
//      bool foundLastFitted = false;
//      int lastBinFitted = -1;
//      int bin = 1;
//      while (!(foundLastFitted)) {
//        if (usedBins.GetBinContent(bin,lastBinEstimated) < 1.0 && lastBinFitted > 0) foundLastFitted = true;
//        else {
//          if (usedBins.GetBinContent(bin,lastBinEstimated) > 0.0) lastBinFitted = bin;
//          bin++;
//        }
//      }
//      for (int bin=lastBinEstimated+1; bin < lastBinFitted+1; bin++) {
//        theMjjFitFunction->MakeHistFromFunction(&bkgTemplate,
//          bkgTemplate.GetBinLowEdge(bin),
//          bkgTemplate.GetBinLowEdge(bin+1),
//          false);
//      }
//    } else {
//      lastBinFitted = lastBin;
//    }
//    std::cout << "AFTER:" << std::endl;
//    bkgTemplate.Print("all");
    // Need a data histogram
	std::cout << "Setting data histogram. " << std::endl;
	m->SetData(theHistogram,firstBin,lastBin);
	m->AddProcess("BKG",false);
	m->SetTemplate("BKG",bkgTemplate);

    // Add variation templates due to JES, if desired
	std::vector<double> jesSigmas(nJES,0.0);
	double step = (2*nSigmas)/(nJES-1); // Fencepost problems
	// Nominal is first
	jesSigmas.at(0) = 0.;
	// Below and above nominal
	for (int i=0; i<int(nJES/2.); i++) {
		jesSigmas.at(i+1) = step*i-nSigmas;
		std::cout<<"JESSIG "<< i+1 <<": "<<jesSigmas.at(i+1)<<std::endl;
		jesSigmas.at(nJES-i-1) = nSigmas-step*i;
		std::cout<<"JESSIG "<< nJES-i-1 <<": "<<jesSigmas.at(nJES-i-1)<<std::endl;
	}
  
	vector<vector<TH2D> > storeAll2DHistograms;
	vector<vector<TH1D> > storeAll1DHistograms;
    if (doJES) {
      for (int component=0; component<nComp; component++) { // Lydia
        vector<std::pair<double,TH1D*> > signalJESVariations; signalJESVariations.clear();
        std::vector<TH1D> storeHistograms;// Lydia
        TH1D nominal(*(TH1D*)infile->Get(nominalJES.c_str()));

        string nomname = nominalJES + "_nom";
        nominal.SetName(nomname.c_str());
        storeHistograms.push_back(nominal); // Lydia

        for (int i=1; i<nJES; i++) {
          string histname = JESVariationNames.at(component).at(i-1);
          std::cout<<"HISTNAME "<<histname<<std::endl;
          TH1D thisVariation(*(TH1D*)infile->Get(histname.c_str()));
          thisVariation.SetName(histname.c_str());
          // Because we are using the fit here, and the fit has different normalisation than the data
          storeHistograms.push_back(thisVariation); //Lydia
        }

        storeAll1DHistograms.push_back(storeHistograms); // Lydia

        for (int i=0; i<nJES; i++) {
          signalJESVariations.push_back(std::make_pair(jesSigmas[i],&storeAll1DHistograms.at(component).at(i)));
        }

        string name = JESComponentNames.at(component);

        std::cout << "Creating systematic " << name << std::endl;
        m->AddSystematic(name.c_str(),-nSigmas,nSigmas);
        m->SetPriorGauss(name.c_str(),0.,1.);
        MjjBATTemplateSyst * thisTemplateSyst = new MjjBATTemplateSyst(true);
        thisTemplateSyst->SetSpectra(signalJESVariations);
        m->SetSystematicVariation("BKG",name.c_str(),thisTemplateSyst);
        
      } // End of loop over JES components Lydia
    }


	// Create bump hunter and set up specs.
	MjjBumpHunter theBumpHunter;
	theBumpHunter.SetMinBumpWidth(2);
	theBumpHunter.SetUseSidebands(false);

	// Make pseudoexperimenter.
	MjjPseudoExperimenter thePseudinator(m);

	// Make statistical tests.
	MjjChi2Test theChi2Test;
	MjjLogLikelihoodTest theLikelihoodTest;
  
    // Make a vector to store them.
	vector<MjjStatisticalTest*> theStatsTests;
	theStatsTests.push_back(&theLikelihoodTest);
	theStatsTests.push_back(&theChi2Test);
	theStatsTests.push_back(&theBumpHunter);

	// Obtain estimate of BumpHunter p-value to see whether any window needs to be removed.
    std::cout << "About to call refitting pseudoexps for the first time" << std::endl;
  
	MjjStatisticsBundle initialStats;
    if (fitToys) {
        initialStats = thePseudinator.GetPseudoExperimentStats_Refitting(&theFitter, theMjjFitFunction, theHistogram, &theBumpHunter,firstBinBH,lastBinBH,100); // std::min(100,nPseudoExpBH)
    } else {
        initialStats = thePseudinator.GetPseudoExperimentStatsOnHistogram(backgroundFromFunc,theHistogram,&theBumpHunter,firstBinBH,lastBinBH,100); //std::min(100,nPseudoExpBH)
    }

	std::pair<double,double> initialBHpval = GetFrequentistPValAndError
		(initialStats.statisticsFromPseudoexperiments,initialStats.originalStatistic);
	double initialpval = initialBHpval.first;
	std::cout << "initialpval is " << initialpval << std::endl;

	//------------------------------------------
	//save initial BH p-value
	initialBHPValue = initialBHpval.first;
	initialBHPValueErr = initialBHpval.second;
	initialBHRangeLow = initialStats.originalFurtherInformation.at(0);
	initialBHRangeHigh = initialStats.originalFurtherInformation.at(1);
	initialBHValue = initialStats.originalStatistic;
	//------------------------------------------

	// IF THE P-VALUE IS TOO LOW:
	// Exclude window corresponding to bump and re-fit to get unbiased bkg estimate.
	bool excludeWindow = (initialpval < thresholdPVal_FindSignal);
	(excludeWindow == true) ? printf("Wish to use an Exclusion Window\n") : printf("Not excluding any window\n");

	double lowEdgeOfWindow = initialStats.originalFurtherInformation.at(0);
	double highEdgeOfWindow = initialStats.originalFurtherInformation.at(1);
	int firstBinInWindow = basicInputHisto->FindBin(lowEdgeOfWindow+1);
	int lastBinInWindow = basicInputHisto->FindBin(highEdgeOfWindow-1);
  
    double myRange = basicInputHisto->GetBinLowEdge(lastBinInWindow+1) - basicInputHisto->GetBinLowEdge(firstBinInWindow);
    if (firstBinInWindow - firstBin >= 2) myRange += basicInputHisto->GetBinWidth(firstBinInWindow-1);
    double mywidth = (myRange)/2.0;
    double mymass = basicInputHisto->GetBinLowEdge(firstBinInWindow)+mywidth;
  
    std::cout << "permitWindow is " << permitWindow << std::endl;
  
	while (initialpval < thresholdPVal_RemoveSignal && permitWindow
        && ((mywidth/mymass)<0.2)) {

        std::cout << "Window width is at " << (mywidth/mymass) << std::endl;

		bumpFound = true;

		// We don't permit windows at the end of the fit.
		// If this is happening it means real excess is in middle somewhere
		// and is pulling fit up so much that ends get really biased.
		// But removing bins at the end of the fit is sort of meaningless.
		//if (firstBinInWindow == firstBin || lastBinInWindow == lastBin) 

		// Lydia Updated so only doesn't permit windows at the beginning of the fit, so now bumps are allowed at the end of the spectrum
		if (firstBinInWindow == firstBin) {
			TGraphErrors thisBHTomography = theBumpHunter.GetBumpHunterTomography();
			thisBHTomography.Sort(&TGraph::CompareY);
			double x; double y;
			thisBHTomography.GetPoint(0,x,y);
			double mosterr = thisBHTomography.GetErrorX(0);
			std::cout << "Wanted window between " << x - mosterr << " and " << x+mosterr << std::endl;
			thisBHTomography.GetPoint(1,x,y);
			double seconderr = thisBHTomography.GetErrorX(1);
			firstBinInWindow = basicInputHisto->FindBin(x - seconderr + 1.0);
			lastBinInWindow = basicInputHisto->FindBin(x + seconderr - 1.0);
		}

		std::cout << "First bin and last bin to exclude are " << firstBinInWindow << " " << lastBinInWindow << ", corresponding to " << lowEdgeOfWindow << ", " << highEdgeOfWindow << std::endl;
		printf( "Window width / max allowed width (half-spectrum) is %i / %i \n", (lastBinInWindow - firstBinInWindow+1), (lastBin - firstBin+1)/2 );

		// Stop making the window bigger if it is larger than half the spectrum
		// or it is too close to either endpoint -- this background estimation is as good
		// as it's going to get.
		if (((lastBinInWindow - firstBinInWindow+1) > (lastBin - firstBin+1)/2.0) ||
				(firstBinInWindow - firstBin < 2) || (lastBin - lastBinInWindow < 2)) break;

		lowEdgeOfWindow = basicInputHisto->GetBinLowEdge(firstBinInWindow);
		highEdgeOfWindow = basicInputHisto->GetBinLowEdge(lastBinInWindow) + basicInputHisto->GetBinWidth(lastBinInWindow);

		std::cout << "Trying to refit." << std::endl;
		theMjjFitFunction->SetExclusionWindowFromRange(lowEdgeOfWindow,highEdgeOfWindow);
		theMjjFitFunction->SetDoWindowExclusion(true);

		// backgroundFromFunc = theFitter.FitAndGetBkgWithDataErr(*theMjjFitFunction,theHistogram,100);
		// Lydia changed -- Kate has added a line after this step with full pseudoexperiments and thus accepts the change.
        std::cout.setstate(std::ios_base::failbit);
		backgroundFromFunc = theFitter.FitAndGetBkgWithNoErr(*theMjjFitFunction,theHistogram,minXForFit,maxXForFit);
        std::cout.clear();

		// Obtain estimate of BumpHunter p-value in region outside window to see if we need
		// to keep going.
		theBumpHunter.SetWindowToExclude(firstBinInWindow,lastBinInWindow);

		MjjStatisticsBundle theseStats;
        if (fitToys) initialStats = thePseudinator.GetPseudoExperimentStats_Refitting(&theFitter, theMjjFitFunction,theHistogram,&theBumpHunter,firstBinBH,lastBinBH,nPseudoExpBH);
        else initialStats = thePseudinator.GetPseudoExperimentStatsOnHistogram
			(backgroundFromFunc,theHistogram,&theBumpHunter,firstBinBH,lastBinBH,nPseudoExpBH);
      
		initialBHpval = GetFrequentistPValAndError
			(theseStats.statisticsFromPseudoexperiments,theseStats.originalStatistic);
		initialpval = initialBHpval.first;
		std::cout << "pval of remaining spectrum is " << initialpval << std::endl;

		double biggestRemainingBumpLowEdge = theseStats.originalFurtherInformation.at(0);
		double biggestRemainingBumpHighEdge = theseStats.originalFurtherInformation.at(1);

		// Make new first bin, last bin in window /if/ pval is still low
		// otherwise save final window size with these variables
		if (initialpval < thresholdPVal_RemoveSignal) {
			if (basicInputHisto->FindBin(biggestRemainingBumpLowEdge+1)
					== lastBinInWindow+1) {
				lastBinInWindow = lastBinInWindow+1;
			} else if (basicInputHisto->FindBin(biggestRemainingBumpHighEdge-1)
					== firstBinInWindow-1) {
				firstBinInWindow = firstBinInWindow-1;
			} else {
				lastBinInWindow = lastBinInWindow+1;
				firstBinInWindow = firstBinInWindow-1;
			}
		} // End of control of window region

      myRange = basicInputHisto->GetBinLowEdge(lastBinInWindow+1) - basicInputHisto->GetBinLowEdge(firstBinInWindow);
      if (firstBinInWindow - firstBin > 2) myRange += basicInputHisto->GetBinWidth(firstBinInWindow-1);
      mywidth = (myRange)/2.0;
      mymass = basicInputHisto->GetBinLowEdge(firstBinInWindow)+mywidth;

	} // End of background recalculation

	// Now outside the loop, add 1 extra bin to low edge of exclusion window and re-fit 
	// Only add one to lower end of window if previously excluded a window AND if lower end of window isn't too close to starting point
	//  -- this background estimation is as good as it's going to get.
    double remainingLogLPVal = -1.0, remainingChi2PVal = -1.0;
	if (excludeWindow && permitWindow){ // - 1 off if >= 2  equivalent to break if < 2 seen in code above
    
        if (firstBinInWindow - firstBin > 2) {
    
  		  std::cout<<" Adding 1 extra bin to exclusion window at low mass end"<<std::endl;

		  firstBinInWindow --; // -1 off low mass end of window
		  lowEdgeOfWindow = basicInputHisto->GetBinLowEdge(firstBinInWindow);
		  highEdgeOfWindow = basicInputHisto->GetBinLowEdge(lastBinInWindow) + basicInputHisto->GetBinWidth(lastBinInWindow);

		  std::cout << "First bin and last bin to exclude are " << firstBinInWindow << " " << lastBinInWindow << ", corresponding to " << lowEdgeOfWindow << ", " << highEdgeOfWindow << std::endl;
		  printf( "Window width / max allowed width (half-spectrum) is %i / %i \n", (lastBinInWindow - firstBinInWindow+1), (lastBin - firstBin+1)/2 );

		  std::cout << "Trying to refit." << std::endl;
		  theMjjFitFunction->SetExclusionWindowFromRange(lowEdgeOfWindow,highEdgeOfWindow);
		  theMjjFitFunction->SetDoWindowExclusion(true);

		  // backgroundFromFunc = theFitter.FitAndGetBkgWithDataErr(*theMjjFitFunction,theHistogram,100);
		  // Lydia changed -- Kate has added a line after this step with full pseudoexperiments and thus accepts the change.
          std::cout.setstate(std::ios_base::failbit);
		  backgroundFromFunc = theFitter.FitAndGetBkgWithNoErr(*theMjjFitFunction,theHistogram,minXForFit,maxXForFit);
          std::cout.clear();
        }

		// Obtain estimate of BumpHunter p-value in region outside window to see if we need
		// to keep going.
		theBumpHunter.SetWindowToExclude(firstBinInWindow,lastBinInWindow);
        theChi2Test.SetWindowToExclude(firstBinInWindow,lastBinInWindow);
        theLikelihoodTest.SetWindowToExclude(firstBinInWindow,lastBinInWindow);

		vector<MjjStatisticsBundle> theseStats;
        if (fitToys) theseStats = thePseudinator.GetPseudoExperimentStats_Refitting
		(&theFitter, theMjjFitFunction, theHistogram,theStatsTests,firstBinBH,lastBinBH,nPseudoExpBH);
        else theseStats = thePseudinator.GetPseudoExperimentStatsOnHistogram
			(backgroundFromFunc,theHistogram,theStatsTests,firstBinBH,lastBinBH,nPseudoExpBH);

	    MjjStatisticsBundle theseLogLStats = theseStats.at(0);
	    MjjStatisticsBundle theseChi2Stats = theseStats.at(1);
	    MjjStatisticsBundle theseBHStats = theseStats.at(2);

		initialpval = GetFrequentistPValAndError
			(theseBHStats.statisticsFromPseudoexperiments,theseBHStats.originalStatistic).first;
		std::cout << "pval of remaining spectrum is " << initialpval << std::endl;
        remainingLogLPVal = GetFrequentistPValAndError
            (theseLogLStats.statisticsFromPseudoexperiments, theseLogLStats.originalStatistic).first;
    	remainingChi2PVal = GetFrequentistPValAndError
		(theseChi2Stats.statisticsFromPseudoexperiments,theseChi2Stats.originalStatistic).first;


		//------------------------------------------
		//save refined BH p-value once the signal window is removed from the fit range
		refinedBHPValue = initialBHpval.first;
		refinedBHPValueErr = initialBHpval.second;
		refinedBHRangeLow = theseBHStats.originalFurtherInformation.at(0);
		refinedBHRangeHigh = theseBHStats.originalFurtherInformation.at(1);
		refinedBHValue = theseBHStats.originalStatistic;
		//------------------------------------------
    }

	// Fit is now ok: repeat with all pseudoexperiments to get background error estimate.
	if( f_noDataErr == false )
		backgroundFromFunc = theFitter.FitAndGetBkgWithDataErr(*theMjjFitFunction,theHistogram,nPseudoExpFit,minXForFit,maxXForFit);

	// Get final fit information.
	fittedParameters.clear();
	fittedParameters = theMjjFitFunction->GetCurrentParameterValues();

	// What were fitted params?
	std::cout << "After fit, parameters are ";
	for (unsigned int i=0; i<fittedParameters.size(); i++) {
		std::cout << "   (" << i << "):" << fittedParameters.at(i);
	}
	std::cout << std::endl;

	// Un-set any exclusions and bump hunt the spectrum
	theBumpHunter.SetUseWindowExclusion(false);
	double bumpHunterStat = theBumpHunter.DoTest(theHistogram,backgroundFromFunc,firstBinBH,lastBinBH);
	vector<double> bumpEdges = theBumpHunter.GetFurtherInformation();
	double lowEdgeOfBump = bumpEdges.at(0);
	double highEdgeOfBump = bumpEdges.at(1);
	std::cout << "BumpHunter results: stat = " << bumpHunterStat << std::endl;
	std::cout << "Low edge, high edge of bump: " << lowEdgeOfBump << " " << highEdgeOfBump << std::endl;

	// Added for checking deficits
	MjjBumpHunter theDeficitHunter;
	theDeficitHunter.SetMinBumpWidth(2);
	theDeficitHunter.SetUseSidebands(false);
	theDeficitHunter.AllowDeficit(true);
	double deficitHunterStat = theDeficitHunter.DoTest(theHistogram,backgroundFromFunc,firstBinBH,lastBinBH);
	vector<double> defEdges = theDeficitHunter.GetFurtherInformation();
	double lowEdgeOfDef = defEdges.at(0);
	double highEdgeOfDef = defEdges.at(1);
    theStatsTests.push_back(&theDeficitHunter);
	std::cout << "DeficitHunter results: stat = " << deficitHunterStat << std::endl;
	std::cout << "Low edge, high edge of bump: " << lowEdgeOfDef << " " << highEdgeOfDef << std::endl;

	// Get bump hunter tomography plot
	TGraphErrors bumpHunterTomography = theBumpHunter.GetBumpHunterTomography();

    // For the random collection of histograms to be saved
	std::vector<TH1D> extraHistograms;

    // Store +/- 1 sigma in nominal fit
	TH1D Bkg_plus1(backgroundFromFunc.GetHistogram());
	Bkg_plus1.SetName("nominalBkgFromFit_plus1Sigma");
	TH1D Bkg_minus1(backgroundFromFunc.GetHistogram());
	Bkg_minus1.SetName("nominalBkgFromFit_minus1Sigma");
	for (int bin =0; bin < backgroundFromFunc.GetHistogram().GetNbinsX()+2; bin++) {
		double bincontent = backgroundFromFunc.GetHistogram().GetBinContent(bin);
		double binerr = backgroundFromFunc.GetHistogram().GetBinError(bin);
		Bkg_plus1.SetBinContent(bin, bincontent + binerr);
		Bkg_minus1.SetBinContent(bin, bincontent - binerr);
	}
	extraHistograms.push_back(Bkg_plus1);
	extraHistograms.push_back(Bkg_minus1);
  
    std::cout << "Made +/-1 hists" << std::endl;

    // Collect cross check info from SWIFT
    TH2D swiftBinsUsed;
    std::vector<TGraph> swiftFitPars;
    if (doSwift) {
      std::cout << "In doSwift loop" << std::endl;
      swiftBinsUsed = theFitter.GetSwiftBinsUsed();
      swiftBinsUsed.SetName("swiftBinsUsed");
      swiftBinsUsed.SetTitle("swiftBinsUsed");
      
  
      std::cout << "Made swift binning hists" << std::endl;
  
      swiftFitPars = theFitter.GetSwiftParameterGraphs();
  
      std::cout << "Collected swift fit parameters" << std::endl;
    }

	// Do alternate fit now, if we are doing that.
	MjjHistogram bkgFromAltFunc;
	MjjHistogram bkgWithFuncChoiceErr;
	TH1D fitChoicePlus1Sig(backgroundFromFunc.GetHistogram());
	std::vector<TH1D> binerrs;
//	TH1D * Alt_bkg_plus1 = 0; TH1D * Alt_bkg_minus1 = 0;
	if (doAlternate) {

		// If there was an exclusion above we want to use the same region here.
		if (bumpFound) {
			theAlternateFunction->SetExclusionWindowFromRange(lowEdgeOfWindow,highEdgeOfWindow);
			theAlternateFunction->SetDoWindowExclusion(true);
		}

        std::cout.setstate(std::ios_base::failbit);
		bkgFromAltFunc = theFitter.FitAndGetBkgWithNoErr(*theAlternateFunction,theHistogram,minXForFit,maxXForFit);
        std::cout.clear();

		// Save result of fit to use later
		vector<double> altfittedParameters = theAlternateFunction->GetCurrentParameterValues();

		// What were fitted params?
		std::cout << "After alternate fit, parameters are  ";
		for (unsigned int i=0; i<altfittedParameters.size(); i++) {
			std::cout << "   (" << i << "):" << altfittedParameters.at(i);
		}
		std::cout << std::endl;

		std::cout << "Using doPEOnData="<<doPEOnData << std::endl;

        // Use results of fit to data to set start parameters for alternate function in fit to PEs
		theAlternateFunction->SetParameterDefaults(altfittedParameters);

		// Get histogram with errors equal to RMS of distance between functions
		// across PEs in each bin. This is the uncertainty that will be used for calculating
		// the search phase p-value using systematics.
        std::cout << "Doing fit diff error calculation" << std::endl;
		bkgWithFuncChoiceErr = theFitter.FitAndGetBkgWithFitDiffErr(*theMjjFitFunction,*theAlternateFunction,theHistogram,minXForFit,maxXForFit, 100,doPEOnData);

		// For now, save distribution of distances in each bin for testing.
		std::vector<std::vector<double> > bincontents = theFitter.getBinContentVectors();
		for (unsigned int bin = 0; bin < bincontents.size(); bin++) {
			TH1D thisbinerr = MakeHistoFromStats(bincontents.at(bin));
			binerrs.push_back(thisbinerr);
		}

		// Create histogram of format which will be used as alternate function
		// to define function choice uncertainty when directionality matters:
		// that is, in the limit setting stage.
		fitChoicePlus1Sig.SetName("alternateFitChoiceSyst_setByRMSTimesDirectionality");
		for (int bin = 1; bin < fitChoicePlus1Sig.GetNbinsX()+1; bin++) {
			double nomquantity = backgroundFromFunc.GetHistogram().GetBinContent(bin);
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

	} // End of if doAlternate

	// Make bump hunter for doing version using uncertainties.
	MjjBumpHunter simpleBumpHunter;
	simpleBumpHunter.SetMinBumpWidth(2);
	simpleBumpHunter.SetUseSidebands(false);
    // Stat uncertainties only
	if (doPValWithSysts) theStatsTests.push_back(&simpleBumpHunter);

	// Make a histogram of errors summed in quadrature for use in finding a p-value with systematics.
	TH1D basicBkgFrom4ParamFit = backgroundFromFunc.GetHistogram();
    TH1D ErrHist(basicBkgFrom4ParamFit);
	if (doPValWithSysts) {

		std::cout << "Starting calculation with uncertainties!" << std::endl;

		for (int bin = 0; bin < ErrHist.GetNbinsX()+2; bin++) {
			double thisval = basicBkgFrom4ParamFit.GetBinError(bin);
			if (doAlternate) {
				double fitchoiceval = bkgWithFuncChoiceErr.GetHistogram().GetBinError(bin);
				std::cout << "Two vals to combine are " << thisval << ", " << fitchoiceval << std::endl;
				ErrHist.SetBinContent(bin,sqrt(pow(thisval,2) + pow(fitchoiceval,2)));
			} else ErrHist.SetBinContent(bin,thisval);
		}
		simpleBumpHunter.SetUseError(ErrHist);

		double newbumpHunterStat = simpleBumpHunter.DoTest(theHistogram,backgroundFromFunc,firstBinBH,lastBinBH);
		std::cout << "Found bump stat " << newbumpHunterStat << std::endl;

	}
  
	// Make histos for residual, diff, significance of diff of fit to data
	MjjSignificanceTests theTestMaker;
	TH1D residualHist = theTestMaker.GetResidual(theHistogram, backgroundFromFunc, firstBin, lastBin);
	TH1D relativeDiffHist = theTestMaker.GetRelativeDifference(theHistogram,backgroundFromFunc, firstBin, lastBin);
	TH1D sigOfDiffHist = theTestMaker.GetSignificanceOfDifference(theHistogram,backgroundFromFunc, firstBin, lastBin);
  
    TH1D residualHistWithSyst;
    if (doPValWithSysts) {
      residualHistWithSyst = theTestMaker.GetResidual(theHistogram,backgroundFromFunc,firstBin,lastBin,&ErrHist);
    }
  
    theBumpHunter.SetUseWindowExclusion(false);
	vector<MjjStatisticsBundle> theStatsBundles;
    if (fitToys) theStatsBundles = thePseudinator.GetPseudoExperimentStats_Refitting
		(&theFitter, theMjjFitFunction,theHistogram,theStatsTests,firstBinBH,lastBinBH,nPseudoExpBH);
    else theStatsBundles = thePseudinator.GetPseudoExperimentStatsOnHistogram
		(backgroundFromFunc,theHistogram,theStatsTests,firstBinBH,lastBinBH,nPseudoExpBH);
	MjjStatisticsBundle logLPseudoStatsOfFit = theStatsBundles.at(0);
	MjjStatisticsBundle chi2PseudoStatsOfFit = theStatsBundles.at(1);
	MjjStatisticsBundle bumpHunterPseudoStatsOfFit = theStatsBundles.at(2);
	MjjStatisticsBundle defHunterPseudoStatsOfFit = theStatsBundles.at(3);

    std::vector<TH1D> theToyDists;
    std::vector<TH1D> theToyFits;
    std::vector<std::vector<double> > theNPsForToys;
    if (fitToys) {
      theToyDists = thePseudinator.GetToyHists();
      theToyFits = thePseudinator.GetFitHists();
      theNPsForToys = thePseudinator.GetToyNPs();
    }

	int nbins = lastBin - firstBin + 1;
	int NDF = nbins - 1 - nPars;

	std::cout << "Absolute values of logl, chi2, BH tests are: " << logLPseudoStatsOfFit.originalStatistic <<
		" " << chi2PseudoStatsOfFit.originalStatistic << " " << bumpHunterPseudoStatsOfFit.originalStatistic << std::endl;

	std::cout << "NDF is " << NDF << std::endl;
	std::cout << "chi2/ndf is " << chi2PseudoStatsOfFit.originalStatistic/NDF << std::endl;

	// Calculate p-values for pseudoexperiment statistics.
	std::pair<double,double> logLPValAndErr = GetFrequentistPValAndError
		(logLPseudoStatsOfFit.statisticsFromPseudoexperiments,logLPseudoStatsOfFit.originalStatistic);
	std::pair<double,double> chi2PValAndErr = GetFrequentistPValAndError
		(chi2PseudoStatsOfFit.statisticsFromPseudoexperiments,chi2PseudoStatsOfFit.originalStatistic);
	std::pair<double,double> bumpHunterPValAndErr = GetFrequentistPValAndError
		(bumpHunterPseudoStatsOfFit.statisticsFromPseudoexperiments,bumpHunterPseudoStatsOfFit.originalStatistic);
	std::pair<double,double> defHunterPValAndErr = GetFrequentistPValAndError
		(defHunterPseudoStatsOfFit.statisticsFromPseudoexperiments,defHunterPseudoStatsOfFit.originalStatistic);

	std::cout << "Log likelihood for our fit was " << logLPseudoStatsOfFit.originalStatistic << std::endl;

	std::cout << "Pval of logL, chi2, and bumpHunter are: " << logLPValAndErr.first << " "
		<< chi2PValAndErr.first << " " << bumpHunterPValAndErr.first << std::endl;

	std::cout << "Pval of largest deficit was: "<< defHunterPValAndErr.first << std::endl;

	// Evaluate output from statistical p-value calculation, if any
	MjjStatisticsBundle BHWithStatsPseudoStatsOfFit;
	std::pair<double,double> BHWithStatsPValAndErr;
	TGraphErrors tomographyWithStats;
	if (doPValWithSysts) {
		BHWithStatsPseudoStatsOfFit = theStatsBundles.at(4);
		BHWithStatsPValAndErr = GetFrequentistPValAndError
			(BHWithStatsPseudoStatsOfFit.statisticsFromPseudoexperiments,BHWithStatsPseudoStatsOfFit.originalStatistic);

		std::cout << "PVal when we consider systematics is: " << BHWithStatsPValAndErr.first << std::endl;
		tomographyWithStats = simpleBumpHunter.GetBumpHunterTomography();
	}

	//------------------------------------------
    // Evaluate distribution of residual(s)
    std::vector<std::pair<TString,TH1D> > residuals;
    residuals.push_back(std::make_pair("residual",residualHist));
    if (doPValWithSysts) residuals.push_back(std::make_pair("residualWithSysts",residualHistWithSyst));
    std::vector<TH1D> distributionHists;
    std::vector<double> resultGausMeans;
    std::vector<double> resultGausWidths;
    std::vector<double> resultHistMeans;
    std::vector<double> resultHistRMS;
    std::vector<TF1> fittedHists;
    for (unsigned int i=0; i<residuals.size(); i++) {
      TString title = residuals.at(i).first;
      TH1D thisResHist = residuals.at(i).second;
      TH1D residualResults(title+"_histOfVals",title+"_histOfVals",40,-10,+10);
      for (int bin=firstBin; bin < lastBin; bin++)
        residualResults.Fill(thisResHist.GetBinContent(bin));
      residualResults.Fit("gaus");
      TF1 fittedGauss(*(TF1*)residualResults.GetFunction("gaus"));
      fittedHists.push_back(fittedGauss);
      double gausmean = fittedGauss.GetParameter(1);
      double gauswidth = fittedGauss.GetParameter(2);
      double histmean = residualResults.GetMean();
      double histRMS = residualResults.GetRMS();
      std::cout << "From fit, mean and width are " << gausmean << " " << gauswidth << std::endl;
      std::cout << "Calculating quantities from function itself, mean and RMS are " << histmean << " " << histRMS << std::endl;
      distributionHists.push_back(residualResults);
      resultGausMeans.push_back(gausmean);
      resultGausWidths.push_back(gauswidth);
      resultHistMeans.push_back(histmean);
      resultHistRMS.push_back(histRMS);
    }

	//------------------------------------------
	//print BH p-values before and after windowd exclusion (if any)
	std::cout << "\n******************************************" << std::endl;
	std::cout << "*** initial values " << std::endl;
	std::cout << "*** BH p-value = " << initialBHPValue << " +/- " << initialBHPValueErr << std::endl;
	std::cout << "*** BH value = " << initialBHValue << std::endl;
	std::cout << "*** BH range = " << initialBHRangeLow << " - " << initialBHRangeHigh << std::endl;
	std::cout << "******************************************" << std::endl;
	std::cout << "*** refined values (after window removal)" << std::endl;
	std::cout << "*** BH p-value = " << refinedBHPValue << " +/- " << refinedBHPValueErr << std::endl;
	std::cout << "*** BH value = " << refinedBHValue << std::endl;
	std::cout << "*** BH range = " << refinedBHRangeLow << " - " << refinedBHRangeHigh << std::endl;
	std::cout << "******************************************" << std::endl;
	std::cout << "*** final values " << std::endl;
	std::cout << "*** BH p-value = " << bumpHunterPValAndErr.first << " +/- " << bumpHunterPValAndErr.second << std::endl;
	std::cout << "*** BH value = " << bumpHunterPseudoStatsOfFit.originalStatistic << std::endl;
	std::cout << "*** BH range = " << bumpHunterPseudoStatsOfFit.originalFurtherInformation.at(0) << " - " << bumpHunterPseudoStatsOfFit.originalFurtherInformation.at(1) << std::endl;
	if (doPValWithSysts) {
	std::cout << "******************************************" << std::endl;
	std::cout << "*** values including systematics " << std::endl;
	std::cout << "*** BH p-value = " << BHWithStatsPValAndErr.first << " +/- " << BHWithStatsPValAndErr.second << std::endl;
	std::cout << "*** BH value = " << BHWithStatsPseudoStatsOfFit.originalStatistic << std::endl;
	std::cout << "*** BH range = " << BHWithStatsPseudoStatsOfFit.originalFurtherInformation.at(0) << " - " << BHWithStatsPseudoStatsOfFit.originalFurtherInformation.at(1) << std::endl;
	}
	std::cout << "******************************************\n" << std::endl;
	//------------------------------------------

	////////////////////////////////////////////////////////////
	// Save everything.

	std::cout << "Writing output file " << outputFileName << std::endl; std::cout.flush();

	TFile * outfile = TFile::Open(outputFileName,"RECREATE");
	outfile->cd();

	// Save fit range
	TVectorD FitRange(2);
	FitRange[0] = minXForFit;
	FitRange[1] = maxXForFit;
	FitRange.Write("FitRange");

	// set correct values to preserve stuff below
	theMjjFitFunction->SetCurrentParameterValues(fittedParameters);

	TH1D basicData = theHistogram.GetHistogram();
	basicData.SetName("basicData");
	basicData.Write();
	TH1D normalizedData = theHistogram.GetNormalizedHistogram();
	normalizedData.SetName("normalizedData");
	normalizedData.Write();

	theMjjFitFunction->GetFitFunction()->SetName("theFitFunction");
	theMjjFitFunction->GetFitFunction()->Write();

	basicBkgFrom4ParamFit.SetName("basicBkgFrom4ParamFit");
	basicBkgFrom4ParamFit.Write();
	TH1D normalizedBkgFrom4ParamFit = backgroundFromFunc.GetNormalizedHistogram();
	normalizedBkgFrom4ParamFit.SetName("normalizedBkgFrom4ParamFit");
	normalizedBkgFrom4ParamFit.Write();

	residualHist.SetName("residualHist");
	residualHist.Write();
	relativeDiffHist.SetName("relativeDiffHist");
	relativeDiffHist.Write();
	sigOfDiffHist.SetName("sigOfDiffHist");
	sigOfDiffHist.Write();

    if (doPValWithSysts) {
      residualHistWithSyst.SetName("residualHistWithSysts");
      residualHistWithSyst.Write();
    }

	logLPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.SetName("logLikelihoodStatHistNullCase");
	logLPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.Write();
	TVectorD statPValErrOfFitToData(3);
	statPValErrOfFitToData[0] = logLPseudoStatsOfFit.originalStatistic;
	statPValErrOfFitToData[1] = logLPValAndErr.first;
	statPValErrOfFitToData[2] = logLPValAndErr.second;
	statPValErrOfFitToData.Write("logLOfFitToData");

	chi2PseudoStatsOfFit.statisticsFromPseudoexperimentsHist.SetName("chi2StatHistNullCase");
	chi2PseudoStatsOfFit.statisticsFromPseudoexperimentsHist.Write();
	statPValErrOfFitToData[0] = chi2PseudoStatsOfFit.originalStatistic;
	statPValErrOfFitToData[1] = chi2PValAndErr.first;
	statPValErrOfFitToData[2] = chi2PValAndErr.second;
	statPValErrOfFitToData.Write("chi2OfFitToData");

	bumpHunterPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.SetName("bumpHunterStatHistNullCase");
	bumpHunterPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.Write();
	statPValErrOfFitToData[0] = bumpHunterPseudoStatsOfFit.originalStatistic;
	statPValErrOfFitToData[1] = bumpHunterPValAndErr.first;
	statPValErrOfFitToData[2] = bumpHunterPValAndErr.second;
	statPValErrOfFitToData.Write("bumpHunterStatOfFitToData");
	bumpHunterTomography.SetName("bumpHunterTomographyFromPseudoexperiments");
	bumpHunterTomography.Write();

	defHunterPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.SetName("deficitHunterStatHistNullCase");
	defHunterPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.Write();
	TVectorD statVector5(5);
	statVector5[0] = defHunterPseudoStatsOfFit.originalStatistic;
	statVector5[1] = defHunterPValAndErr.first;
	statVector5[2] = defHunterPValAndErr.second;
        statVector5[3] = lowEdgeOfDef;
        statVector5[4] = highEdgeOfDef;
	statVector5.Write("deficitHunterStatPValErrLowHigh");

	TVectorD ndf(1);
	ndf[0] = NDF;
	ndf.Write("NDF");

	TVectorD bumpHunterStatLowHigh(3);
	bumpHunterStatLowHigh[0] = bumpHunterStat;
	bumpHunterStatLowHigh[1] = lowEdgeOfBump;
	bumpHunterStatLowHigh[2] = highEdgeOfBump;
	bumpHunterStatLowHigh.Write("bumpHunterPLowHigh");

	TVectorD paramsForBasicFit(nPars);
	for (unsigned int i=0; i<fittedParameters.size(); i++) paramsForBasicFit[i] = fittedParameters.at(i);
	paramsForBasicFit.Write("fittedParameters");

	TVectorD excludeWindowNums(3);
	excludeWindowNums[0] = (double) excludeWindow;
	excludeWindowNums[1] = basicData.GetBinLowEdge(firstBinInWindow);
	excludeWindowNums[2] = basicData.GetBinLowEdge(lastBinInWindow)+basicData.GetBinWidth(lastBinInWindow);
    excludeWindowNums.Write("excludeWindowNums");

	TVectorD bumpFoundVector(1);
	bumpFoundVector[0] = bumpFound;
	bumpFoundVector.Write("bumpFound");

	TVectorD statPValErrOfFitToDataInitial(3);
	statPValErrOfFitToDataInitial[0] = initialBHValue;
	statPValErrOfFitToDataInitial[1] = initialBHPValue;
	statPValErrOfFitToDataInitial[2] = initialBHPValueErr;
	statPValErrOfFitToDataInitial.Write("bumpHunterStatOfFitToDataInitial");

	TVectorD bumpHunterStatLowHighInitial(3);
	bumpHunterStatLowHighInitial[0] = initialBHValue;
	bumpHunterStatLowHighInitial[1] = initialBHRangeLow;
	bumpHunterStatLowHighInitial[2] = initialBHRangeHigh;
	bumpHunterStatLowHighInitial.Write("bumpHunterPLowHighInitial");

	TVectorD statPValErrOfFitToDataRefined(3);
	statPValErrOfFitToDataRefined[0] = refinedBHValue;
	statPValErrOfFitToDataRefined[1] = refinedBHPValue;
	statPValErrOfFitToDataRefined[2] = refinedBHPValueErr;
	statPValErrOfFitToDataRefined.Write("bumpHunterStatOfFitToDataRefined");

	TVectorD bumpHunterStatLowHighRefined(3);
	bumpHunterStatLowHighRefined[0] = refinedBHValue;
	bumpHunterStatLowHighRefined[1] = refinedBHRangeLow;
	bumpHunterStatLowHighRefined[2] = refinedBHRangeHigh;
	bumpHunterStatLowHighRefined.Write("bumpHunterPLowHighRefined");
 
      if (permitWindow && excludeWindow) {
        TVectorD pValuesRemainderOfWindow(3);
        pValuesRemainderOfWindow[0] = initialpval;
        pValuesRemainderOfWindow[1] = remainingLogLPVal;
        pValuesRemainderOfWindow[2] = remainingChi2PVal;
        pValuesRemainderOfWindow.Write("BHLogLAndChi2OfRemainderAfterWindow");
      }

        // Write bump hunter with statistics results if applicable
	if (doPValWithSysts) {
		tomographyWithStats.Write("TomographyPlotWithStats");
		TVectorD thesePVals(3);
		BHWithStatsPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.SetName("bumpHunterStatsWSyst");
		BHWithStatsPseudoStatsOfFit.statisticsFromPseudoexperimentsHist.Write();
		thesePVals[0] = BHWithStatsPseudoStatsOfFit.originalStatistic;
		thesePVals[1] = bumpHunterPValAndErr.first;
		thesePVals[2] = bumpHunterPValAndErr.second;
		thesePVals.Write("bumpHunterPValErrWithStats");
	}
 
	// Write additional histograms saved by optional tests
	for (unsigned int i=0; i<extraHistograms.size(); i++)
		extraHistograms.at(i).Write();
        // Save cross-check histograms from function choice uncertainty
        // calculation, if applicable

	if (doAlternate) {
		for (unsigned int bin = 0; bin < binerrs.size(); bin++) {
			TString histname(Form("BinErr_%d",bin));
			binerrs.at(bin).SetName(histname);
			binerrs.at(bin).Write();
		}
	}
  
    // Write results of residual distribution analysis
    for (unsigned int i=0; i<residuals.size(); i++) {
    
      TString name = residuals.at(i).first;
      TH1D distHist = distributionHists.at(i);
      distHist.SetName(name+"_distributionHist");
      distHist.Write();
      TVectorD meansAndWidths(4);
      meansAndWidths[0] = resultGausMeans.at(i);
      meansAndWidths[1] = resultGausWidths.at(i);
      meansAndWidths[2] = resultHistMeans.at(i);
      meansAndWidths[3] = resultHistRMS.at(i);
      meansAndWidths.Write(name+"_gausMeanWidthHistMeanWidth");
      
      std::cout << "doing fittedhists" << std::endl;
      fittedHists.at(i).Write(name+"_fittedHist");
      std::cout << "done" << std::endl;
    }

    // Write tgraphs of parameters
    if (doSwift) {
      for (unsigned int i=0; i< swiftFitPars.size(); i++) {
        string saveas = Form("evolution_parameter%i",i);
        swiftFitPars.at(i).Write(saveas.c_str());
      }
  
      swiftBinsUsed.Write();
    }

    // Sanity check for JES variation implementation: are things moved by the right amount?
    // Reproduce each of the input templates by extracting predictions from model.
    if (doJES) {
          
      // Make a directory for JES cross checks
      TDirectory * myDir = outfile->mkdir("JESchecks");
      myDir->cd();
      
      // Store some cross check plots that ought to be the same as
      for (unsigned int parIndex = 0; parIndex < m->GetNParameters(); parIndex++) {

        std::vector<double> testPoints = {-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0};
        std::vector<std::string> pointNames = {"m3p0","m2p0","m1p0","0","1p0","2p0","3p0"};
        for (unsigned int testPoint = 0; testPoint < testPoints.size(); testPoint++) {
          std::string histName = Form("VerifyJESApplication_%s_%s",m->GetSystematic(parIndex)->GetName().c_str(),pointNames.at(testPoint).c_str());
          TH1D JESHist(basicData);
          JESHist.SetName(histName.c_str());
          JESHist.SetTitle(histName.c_str());
          JESHist.Reset(0);
          vector<double> parVals(m->GetNParameters(),0.0);
          parVals.at(parIndex) = testPoints.at(testPoint);
          vector<double> binContents = m->Expectation(parVals);
          for (unsigned int bin=0; bin < binContents.size(); bin++) {
            if (std::isnan(binContents.at(bin))) JESHist.SetBinContent(bin,0.0);
            else JESHist.SetBinContent(bin,binContents.at(bin));
          }
          JESHist.Write();
        }

      }
      
      // Now store histograms and fits from toys. Compute a residual for each
      int nToKeep = std::min(25,nPseudoExpBH);
      if (f_saveToyFits) nToKeep = theToyDists.size();
      for (int toy=0; toy< nToKeep; toy++) {
        TH1D thisToy = theToyDists.at(toy);
        std::string toyName = Form("toy%i",toy);
        thisToy.Write(toyName.c_str());
        TH1D thisFit = theToyFits.at(toy);
        std::string fitName = Form("fitToToy%i",toy);
        thisFit.Write(fitName.c_str());
        MjjHistogram mjjToy(&thisToy);
        MjjHistogram mjjFit(&thisFit);
        TH1D thisResidual = theTestMaker.GetResidual(mjjToy,mjjFit, firstBin, lastBin);
        std::string residualName = Form("residualToy%i",toy);
        thisResidual.Write(residualName.c_str());
        std::vector<double> theseNPs = theNPsForToys.at(toy);
    	TVectorD toypars(theseNPs.size());
        for (unsigned int par=0; par < theseNPs.size(); par++) toypars[par] = theseNPs.at(par);
        std::string parsName = Form("NPsForToy%i",toy);
	    toypars.Write(parsName.c_str());
      }
      outfile->cd();
      myDir->Write();
    }

	outfile->Close();
	infile->Close();

	if( f_saveExclusion ){
		std::ofstream fileOut;
		fileOut.open( "exclusionList.txt", ios::out | ios::app);
		fileOut << inputFileName << " " << dataMjjHistoName << " " << thresholdPVal_RemoveSignal << "\n";
		fileOut << firstBinInWindow << " " << lastBinInWindow << " " << logLPValAndErr.first << " " << chi2PValAndErr.first << " " << bumpHunterPValAndErr.first << "\n";
		fileOut.close();
	}

	totaltime.Stop();
	std::cout << "Process ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

    delete infile;
    delete outfile;

	return 0;

}
