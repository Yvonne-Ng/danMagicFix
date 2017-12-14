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

	bool f_useScaled = false;
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
	double Ecm = settings->GetValue("Ecm",8000.0);

	// Do we use the systematic uncertainties in the p-value calculation?
	bool doPValWithSysts = settings->GetValue("doPValWithSysts",false);
	std::cout << "Doing pvalue with systematics? " << doPValWithSysts << "!" << std::endl;

	// Open files
	TFile * infile = TFile::Open(inputFileName,"READ");
	TH1::AddDirectory(kFALSE);

	int nPseudoExp = settings->GetValue("nPseudoExp",1e3);
	std::cout << "nPseudoExp is "<< nPseudoExp << std::endl;

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

	// Range for data fit
	double minX = settings->GetValue("minXForFit",-1.0);
	double maxX = settings->GetValue("maxXForFit",-1.0);

	int functionNumber = settings->GetValue("functionCode",9);
	int nPars = settings->GetValue("nParameters",3);

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

	////////////////////////////////////////////////////////////
	// Beginning search phase


	std::cout << "minX, maxX: " << minX << " " << maxX << std::endl;

	// Lydia changed this part so picks start and end bins that are fully spanned by range minX-maxX
	// e.g. user picks minX maxX of 2000-3000 only bins fully within this range are fitted
	// rather than fitting e.g. 1950-350

	int firstBin, lastBin;
	if (minX < theHistogram.GetHistogram().GetBinLowEdge(theHistogram.GetFirstBinWithData()) || minX < 0) firstBin = theHistogram.GetFirstBinWithData();
	else firstBin = theHistogram.GetHistogram().FindBin(minX)+1;


	if (maxX > theHistogram.GetHistogram().GetBinLowEdge(theHistogram.GetLastBinWithData()+1) || maxX < 0) lastBin = theHistogram.GetLastBinWithData();
	else lastBin = theHistogram.GetHistogram().FindBin(maxX)-1;
	double minXForFit = theHistogram.GetHistogram().GetBinLowEdge(firstBin);
	double maxXForFit = theHistogram.GetHistogram().GetBinLowEdge(lastBin+1);

	// Create fit function(s)
	MjjFitFunction * theMjjFitFunction;
	MjjFitFunction * theAlternateFunction;
	vector<std::pair<int,MjjFitFunction**> > functionsAndCodes;
	functionsAndCodes.push_back(std::make_pair(functionNumber,&theMjjFitFunction));
	if (doAlternate) functionsAndCodes.push_back(std::make_pair(alternateFuncNumber,&theAlternateFunction));

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

	TH1D basicData = theHistogram.GetHistogram();
	basicData.SetName("basicData");
	basicData.Write();
	TH1D normalizedData = theHistogram.GetNormalizedHistogram();
	normalizedData.SetName("normalizedData");
	normalizedData.Write();

	theMjjFitFunction->GetFitFunction()->SetName("theFitFunction");
	theMjjFitFunction->GetFitFunction()->Write();


	outfile->Close();
	infile->Close();

	totaltime.Stop();
	std::cout << "Process ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

        delete infile;
        delete outfile;

	return 0;

}
