#include "TH1.h"
#include <iostream>
#include "Bayesian/MjjFitter.h"
#include "Bayesian/MjjHistogram.h"
#include "Bayesian/MjjFitFunction.h"
#include "Bayesian/MjjStatisticalTest.h"
#include "Bayesian/MjjChi2Test.h"
#include "Bayesian/MjjPseudoExperimenter.h"
#include "Bayesian/MjjStatisticsBundle.h"
#include "Bayesian/MjjSignificanceTests.h"
#include "Bayesian/MjjBumpHunter.h"
#include "Bayesian/MathFunctions.h"

#include "TError.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TVectorD.h"
#include "TGraphAsymmErrors.h"
#include "TEnv.h"

int main (int argc,char **argv) {

  ////////////////////////////////////////////////////////////
  // Initialisation

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

  TString inputFileForm = settings->GetValue("inputFileForm","");
  std::string outputFileName = settings->GetValue("outputFileName","");

  int nSplitPatterns = settings->GetValue("nSplitPatterns",0);

  std::cout << "NSplitPatterns is " << nSplitPatterns << std::endl;

  vector<int> nsplits;
  vector<vector<double> > masses;
  for (int split=1; split<nSplitPatterns+1; split++) {
    TString namenummasses = Form("nMasses%d",split);
    int nummasses = settings->GetValue(namenummasses,0);
    std::cout << "Got nMasses " << nummasses << std::endl;
    TString namenumsplits = Form("nSplits%d",split);
    int numsplits = settings->GetValue(namenumsplits,0);
    std::cout << "Each has numsplits " << numsplits << std::endl;
    nsplits.push_back(numsplits);
    vector<double> holdmasses;
    for (int thismass=0; thismass<nummasses; thismass++) {
      TString namemass = "split%dmass%d";
      TString fullname = Form(namemass,split,thismass);
      //holdmasses.push_back(settings->GetValue(Form(namemass,split,thismass),0));
      holdmasses.push_back(settings->GetValue(fullname,0));
      std::cout << "Just added " << holdmasses.at(thismass) << std::endl;
    }
    masses.push_back(holdmasses);
  }

  double luminosity = settings->GetValue("luminosity",0.0);

  // Prevent segfaults all over the place
  TH1::AddDirectory(kFALSE);

  ////////////////////////////////////////////////////////////
  // Creation of limit-setting plots

  // Create graphs that will define limit-setting plot
  TGraphErrors * observed = new TGraphErrors();
  TGraphAsymmErrors * expectedWithOneSigma = new TGraphAsymmErrors();
  TGraphAsymmErrors * expectedWithTwoSigma = new TGraphAsymmErrors();

  int index = 0;

  for (int thissplitpattern=0; thissplitpattern<nSplitPatterns; thissplitpattern++) {

    vector<double> thesesplitmasses = masses.at(thissplitpattern);
    int thisnsplits = nsplits.at(thissplitpattern);

    std::cout << "This nsplits is " << thisnsplits << std::endl;

    for (unsigned int mass = 0; mass<thesesplitmasses.size(); mass++) {

      double CLObserved = -1;
      vector<double> CLPseudo;

      for (int split=0; split<thisnsplits; split++) {

        int thismass = int(thesesplitmasses.at(mass));
        std::cout << "Split and mass are " << split << ", " << thesesplitmasses.at(mass) << std::endl;  
        TString thisFileName;
        if (thisnsplits > 1) 
          thisFileName = Form(inputFileForm,thismass,split);
        else thisFileName = Form(inputFileForm,thismass);
        std::cout << "Opening file " << thisFileName << std::endl;

        // Read in from file containing this likelihood function
        TFile inLikelihoodFile(thisFileName, "READ");
        bool NotOK = inLikelihoodFile.IsZombie();
        if (NotOK) continue; 

        // Should all have same value for real data. Set once.
        if (split==0 ) {
          // Get CL
          TVectorD* CLOfRealLikelihood = (TVectorD*) inLikelihoodFile.Get("CLOfRealLikelihood");
          CLObserved = (*CLOfRealLikelihood)[0]/luminosity;

          std::cout << "True CL is " << (*CLOfRealLikelihood)[0]<<"/" <<luminosity << " = "<<  (*CLOfRealLikelihood)[0]/luminosity << std::endl;

          // Set point in observed graph
          observed->SetPoint(index, thesesplitmasses.at(mass), CLObserved);
          observed->SetPointError(index, 0, 0);

        }

        // Read tree from file
        TTree * tree = (TTree*) inLikelihoodFile.Get("ensemble_test");
        
        if (tree){
          double CL;
          tree->SetBranchAddress("95quantile_marginalized_2", &CL);
          for (int i = 0; i<tree->GetEntries(); i++) {
            tree->GetEvent(i);
            CLPseudo.push_back(CL);
            std::cout << "individual pseudo CL = " << CL << "/" << luminosity << " = " << CL/luminosity << std::endl;
          }
        } else std::cout << "pseudo CL tree does not exist" << std::endl;

      inLikelihoodFile.Close();
 
      } // End of splits

      std::cout << "For mass " << thesesplitmasses.at(mass) << " got real CL " << CLObserved << std::endl;

      if (CLPseudo.size() > 0){
        vector<double> bandLimits = GetCenterAndSigmaDeviations(CLPseudo);
        for(std::vector<double>::size_type p=0; p!=bandLimits.size(); ++p) bandLimits[p] /= luminosity;

        // Set points in expected graphs
        expectedWithOneSigma->SetPoint(index, thesesplitmasses.at(mass), bandLimits.at(2));
        expectedWithOneSigma->SetPointError(index, 0, 0, (bandLimits.at(2)-bandLimits.at(1)),
                                 (bandLimits.at(3)-bandLimits.at(2)));
        expectedWithTwoSigma->SetPoint(index, thesesplitmasses[mass], bandLimits.at(2));
        expectedWithTwoSigma->SetPointError(index, 0, 0, (bandLimits.at(2)-bandLimits.at(0)),
                                 (bandLimits.at(4)-bandLimits.at(2)));
        ///////////std::cout<<"\nbandLimits.at"<<"\n(4)"<<bandLimits.at(4)<<"\n(3)"<<bandLimits.at(3)<<"\n(2)"<<bandLimits.at(2)<<"\n(1)"<<bandLimits.at(1)<<"\n(0)"<<bandLimits.at(0)<<std::endl;

        std::cout << "Center of expected, observed: " << bandLimits.at(2) << " " << CLObserved << std::endl;
      }

      index +=1;

    }
  }

  ////////////////////////////////////////////////////////////
  // Save new TGraphs to the output file

  std::cout << "Writing output file " << outputFileName << std::endl;

  TFile outfile(outputFileName.c_str(),"RECREATE");
  outfile.cd();

  observed->SetName("observedXSecAccVersusMass");
  observed->Write();
  expectedWithOneSigma->SetName("expectedXSecAccVersusMass_oneSigma");
  expectedWithOneSigma->Write();
  expectedWithTwoSigma->SetName("expectedXSecAccVersusMass_twoSigma");
  expectedWithTwoSigma->Write();

  outfile.Close();

  return 0;

}

