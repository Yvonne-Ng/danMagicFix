#include "TH1.h"
#include "TF1.h"
#include <iostream>
#include <iomanip>
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

#include "Bayesian/MjjBATBumpHunter.h"
#include <BAT/BCLog.h>

//#include "Bayesian/TBumpHunter.h"

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

  // Start reading input configuration
  TString configFile;
  TString searchphaseFile;
  TString outputFileName;
  int ip=1;
  bool isNominal = true;
  int nPE = 0;
  while (ip<argc) {

    if (string(argv[ip]).substr(0,2)=="--") {

      //config file
      if (string(argv[ip])=="--config") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          configFile=argv[ip+1];
          ip+=2;
        } else {std::cout<<"\nno config file inserted"<<std::endl; break;}
      }

      // Search phase file to use
      else if (string(argv[ip])=="--searchphasefile") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          searchphaseFile=argv[ip+1];
          ip+=2;
        } else {std::cout<<"\nno search phase file inserted"<<std::endl; break;}
      }
      
      // Output file name to use
      else if (string(argv[ip])=="--outputfile") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          outputFileName=argv[ip+1];
          ip+=2;
        } else {std::cout<<"\nno output file name inserted"<<std::endl; break;}
      }

      // Specify which number PE to do.
      else if (string(argv[ip])=="--PEnumber") {
        if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
          isNominal = false;
          sscanf(argv[ip+1],"%i",&nPE);
          ip+=2;
        } else {std::cout<<"\nno PE number specified"<<std::endl; break;}
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

  // Open files
  TFile * infile = TFile::Open(searchphaseFile,"READ");
  TH1::AddDirectory(kFALSE);

  // Range for data fit
  double minX = settings->GetValue("minXForFit",-1.0);
//  double maxX = settings->GetValue("maxXForFit",-1.0);

  bool doAlternate = settings->GetValue("doAlternateFunction",true);
  bool doFunctionErrors = settings->GetValue("doFunctionErrors",true);

  TH1D * basicData = (TH1D*) infile->Get("basicData");
  MjjHistogram theHistogram(basicData, nPE);

  TH1D * basicBkgFrom4ParamFit = (TH1D*) infile->Get("basicBkgFrom4ParamFit");
  MjjHistogram backgroundFromFunc(basicBkgFrom4ParamFit);
  TH1D weightsHist = theHistogram.GetWeightsHistogram();
  backgroundFromFunc.SetEffectiveFromBasicAndWeights(&weightsHist);
  
  TH1D * Bkg_plus1 = (TH1D*) infile->Get("nominalBkgFromFit_plus1Sigma");
  TH1D * variedFitHist = (TH1D*) infile->Get("alternateBkgFromFit_0Sigma");
  
  int firstBin = 0, lastBin = 0;
  if (minX < theHistogram.GetHistogram().GetBinLowEdge(theHistogram.GetFirstBinWithData()) || minX < 0) firstBin = theHistogram.GetFirstBinWithData();
  else firstBin = theHistogram.GetHistogram().FindBin(minX)+1;

  ////////////////////////////////////////////////////////////
  // Obtain p-value estimate from this background model taking
  // systematics into account, if systematics are activated.
  
  // open log file
  std::string logfile = "/log.txt";
  const char* logfilestr = logfile.c_str();
  BCLog::OpenLog(logfilestr);
  BCLog::SetLogLevel(BCLog::detail);

  // Create new MjjBATModel object with only one process: background
  MjjBATModel * m = new MjjBATModel();
  
  TH1D bkgTemplate(*basicBkgFrom4ParamFit);
  bkgTemplate.SetName("BkgTemplate");
  TH1D bkgVariation(*Bkg_plus1);
  bkgVariation.SetName("BkgVariation");
  bkgVariation.Add(&bkgTemplate,-1.0);
  
  if (doFunctionErrors) m->AddProcess("BKG",true,-3.0,3.0);
  else m->AddProcess("BKG",false,0.,0.);
  m->SetTemplate("BKG",bkgTemplate,bkgVariation);
  vector<string> newparams = m->GetProcess(m->GetProcessIndex("BKG"))->GetParamNames();
  for (unsigned int i=0; i<newparams.size(); i++) {
      m->SetPriorGauss(newparams.at(i).c_str(),0.,1.);
  }

  TH1D * variedfit, * variedfitp1,* variedfitm1, *sizeOfDiff;
  vector<std::pair<double,TH1D*> > fitvariedtemplates;
/*  if (doAlternate) {
	variedfit = new TH1D(*variedFitHist);
	variedfit->SetName("fitplusonesigmachoice");
	fitvariedtemplates.push_back(std::make_pair(0,&bkgTemplate));
        fitvariedtemplates.push_back(std::make_pair(1,variedfit));

	m->AddSystematic("FUNCCHOICE",0,1.0);
	m->SetPriorGauss("FUNCCHOICE",0.,1.);
	MjjBATTemplateSyst * FuncChoiceOnBkg = new MjjBATTemplateSyst(true);
	FuncChoiceOnBkg->SetSpectra(fitvariedtemplates);
    m->SetSystematicVariation("BKG","FUNCCHOICE",FuncChoiceOnBkg); */

//  }

  MjjHistogram MjjScaledBkg(&bkgTemplate);
  TH1D ErrHist(*basicBkgFrom4ParamFit);
  for (int bin = 0; bin < ErrHist.GetNbinsX()+2; bin++) {
    double thisval = basicBkgFrom4ParamFit->GetBinError(bin);//pow(sizeOfDiff->GetBinContent(bin),2) + pow(bkgVariation.GetBinContent(bin),2);
    double fitchoiceval = fabs(basicBkgFrom4ParamFit->GetBinContent(bin) - variedFitHist->GetBinContent(bin));
    std::cout << "Two vals to combine are " << thisval << ", " << fitchoiceval << std::endl;
    ErrHist.SetBinContent(bin,sqrt(pow(thisval,2) + pow(fitchoiceval,2)));
  }

  MjjBumpHunter simpleBumpHunter;
  simpleBumpHunter.SetMinBumpWidth(2);
  simpleBumpHunter.SetUseSidebands(false);
  simpleBumpHunter.SetUseError(ErrHist);
  simpleBumpHunter.SetNoisy(true);
  double theTestStatWithSyst = simpleBumpHunter.DoTest(theHistogram,backgroundFromFunc,firstBin,lastBin);
  TGraphErrors newBHTomographyStat = simpleBumpHunter.GetBumpHunterTomography();
  vector<double> bumpLowHigh = simpleBumpHunter.GetFurtherInformation();

  // Now have model set up: feed it into a BumpHunter.
  MjjBATBumpHunter theBHWithSyst(m);
  theBHWithSyst.SetMinBumpWidth(2);
  theBHWithSyst.SetUseSidebands(false);

  double otherTestStatWithSyst;
  if (isNominal) otherTestStatWithSyst = theBHWithSyst.DoTest(theHistogram,MjjScaledBkg,firstBin,lastBin);
  else {
    TH1D pseudo(*basicData);
    TString pseudoname(Form("%s_internal_pseudo",pseudo.GetName()));
    pseudo.SetName(pseudoname);
    pseudo.Clear();
    theHistogram.PoissonFluctuateBinByBin(&pseudo);
    MjjHistogram pseudoHist(&pseudo);
    otherTestStatWithSyst = theBHWithSyst.DoTest(pseudoHist,pseudoHist,firstBin,lastBin);
  }
  
  // Get a tomography plot and the p-value.
  TGraphErrors oldBHTomographyStat = theBHWithSyst.GetBumpHunterTomography();
  TGraphErrors WillStat = theBHWithSyst.GetWillTomography();

  std::cout << "With systematics, calculated a BH test statistic of " << theTestStatWithSyst << std::endl;
  std::cout << "My older version gave " << otherTestStatWithSyst << std::endl;

  //vector<double> bumpLowHigh = theBHWithSyst.GetFurtherInformation();
  
  std::vector<TH1D> theposts = theBHWithSyst.GetPosteriors();
  std::vector<TH1D> thewillposts = theBHWithSyst.GetWillPosteriors();

  ////////////////////////////////////////////////////////////
  // Save everything.

  std::cout << "Writing output file " << outputFileName << std::endl; std::cout.flush();

  TFile * outfile = TFile::Open(outputFileName,"RECREATE");
  outfile->cd();
  newBHTomographyStat.Write("newTomography");
  oldBHTomographyStat.Write("oldTomography");
//  WillStat.Write("WillTomography");
  TVectorD bumpHunterStatLowHigh(3);
  bumpHunterStatLowHigh[0] = theTestStatWithSyst;
  bumpHunterStatLowHigh[1] = bumpLowHigh.at(0);
  bumpHunterStatLowHigh[2] = bumpLowHigh.at(1);
  bumpHunterStatLowHigh.Write("bumpHunterPLowHigh");
/*  for (unsigned int i=0; i< theposts.size(); i++) {
    TString name(Form("posterior_%i", (int) i));
    theposts.at(i).Scale(1.0/theposts.at(i).Integral());
    theposts.at(i).Write(name);
    TString alternatename(Form("will_posterior_%i", (int) i));
    thewillposts.at(i).Write(alternatename);
  } */
  ErrHist.Write("ErrHist");
  outfile->Close();
  infile->Close();

  totaltime.Stop();
  std::cout << "Process ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

  return 0;
  
}

