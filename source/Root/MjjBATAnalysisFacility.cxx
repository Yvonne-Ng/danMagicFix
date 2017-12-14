/*
 * Copyright (C) 2008-2012, Daniel Kollar and Kevin Kroeninger.
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

// ---------------------------------------------------------

#include "Bayesian/MjjBATAnalysisFacility.h"

// ---------------------------------------------------------
MjjBATAnalysisFacility::MjjBATAnalysisFacility(MjjBATModel * mtf, unsigned int seed)
 : fRandom(new TRandom3(seed))
 , fFlagMCMC(false)
 , fLogLevel(BCLog::nothing)
{
   fMBM = mtf;
   BCLog::OutDetail(Form("Prepared Analysis Facility for MTF model \'%s\'",mtf->GetName().c_str()));
}

// ---------------------------------------------------------
MjjBATAnalysisFacility::~MjjBATAnalysisFacility()
{
}

// ---------------------------------------------------------
TH1D MjjBATAnalysisFacility::BuildEnsemble(const std::vector<double> & parameters)
{

   // create new histogram
   TH1D hist( *(fMBM->GetData()) );

   TH1D * weights = fMBM->GetWeights();

   // get number of bins
   int nbins = hist.GetNbinsX();

   vector<double> expectation = fMBM->Expectation(parameters);

   // loop over all bins
   for (int ibin = 1; ibin <= nbins; ++ibin) {

      double weight = weights->GetBinContent(ibin);

      double observation = fRandom->Poisson(expectation.at(ibin));

      hist.SetBinContent(ibin, observation*weight);
      hist.SetBinError(ibin, sqrt(observation)*weight);
   }

   // return histograms
   return hist;
}

// ---------------------------------------------------------
TTree * MjjBATAnalysisFacility::BuildEnsembles(TTree * tree, int nensembles)
{

   BCLog::OutDetail(Form("MTF Building %d ensambles.",nensembles));

   // get number of parameters
   int nparameters = fMBM->GetNParameters();

   // create tree variables for input tree
   std::vector<double> parameters(nparameters);

   // set branch addresses
   for (int i = 0; i < nparameters; ++i) {
      tree->SetBranchAddress(Form("Parameter%i", i), &(parameters[i]));
   }

   // create tree
   TTree * tree_out = new TTree("ensembles", "ensembles");

   // get number of bins
   int nbins = fMBM->GetData()->GetNbinsX();

   std::cout << "making branches" << std::endl;

   // create vector of bins
   std::vector<double> nbins_matrix(nbins);

   std::vector<double> in_parameters(nparameters);

   // loop over bins
   for (int ibin = 1; ibin <= nbins; ++ibin) {
      // create branches
      tree_out->Branch(Form("bin_%i", ibin),
                          &(nbins_matrix[ibin-1]), "n/D");
   }

   for (int i = 0; i < nparameters; ++i) {
      tree_out->Branch(Form("parameter_%i", i), &in_parameters[i], Form("parameter_%i/D", i));
   }

   std::cout << "at ensemble loop" << std::endl;

   // loop over ensembles
   for (int iensemble = 0; iensemble < nensembles; ++iensemble) {
      // get random event from tree
      int index = (int) fRandom->Uniform(tree->GetEntries());
      tree->GetEntry(index);

      std::cout << "building ensemble " << iensemble << std::endl;
      // create ensembles
      TH1D histogram = BuildEnsemble(parameters);

      // loop over bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
        // fill tree variable
        nbins_matrix[ibin-1] = histogram.GetBinContent(ibin);
      }

      // copy parameter information
      for (int i = 0; i < nparameters; ++i) {
         in_parameters[i] = parameters.at(i);
      }

      tree_out->Fill();
   }

   return tree_out;
}

// ---------------------------------------------------------
TTree * MjjBATAnalysisFacility::BuildEnsembles(const std::vector<double> & parameters, int nensembles)
{

   BCLog::OutDetail(Form("MTF Building %d ensambles.",nensembles));

   // create tree
   TTree * tree = new TTree("ensembles", "ensembles");

   // prepare the tree variables
   // get number of bins
   int nbins = fMBM->GetData()->GetNbinsX();

   // create new matrix row
   std::vector<double> nbins_matrix(nbins);

   // get number of parameters
   int nparameters = fMBM->GetNParameters();

   std::vector<double> in_parameters(nparameters);

   // create branches
   // loop over bins
   for (int ibin = 1; ibin <= nbins; ++ibin) {
      // create branches
      tree->Branch(Form("bin_%i", ibin),
                   &(nbins_matrix)[ibin-1], "n/D");
   }

   for (int i = 0; i < nparameters; ++i) {
      tree->Branch(Form("parameter_%i", i), &in_parameters[i], Form("parameter_%i/D", i));
   }

   // loop over ensembles
   for (int iensemble = 0; iensemble < nensembles; ++iensemble) {
      // create ensembles
      TH1D histogram = BuildEnsemble(parameters);

      // copy information from histograms into tree variables
      // loop over bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {
         // fill tree variable
         (nbins_matrix)[ibin-1] = histogram.GetBinContent(ibin);
      }

      // copy parameter information
      for (int i = 0; i < nparameters; ++i) {
         in_parameters[i] = parameters.at(i);
      }

      tree->Fill();
   }

   return tree;
}

// ---------------------------------------------------------
TTree * MjjBATAnalysisFacility::PerformEnsembleTest(const std::vector<double> & parameters, int nensembles, int firstBinToUse, int lastBinToUse)
{
   // create new tree
   TTree * tree = 0;

   // create ensembles
   tree = BuildEnsembles(parameters, nensembles);

   // perform ensemble test
   return PerformEnsembleTest(tree, nensembles, firstBinToUse, lastBinToUse);
}

// ---------------------------------------------------------
TTree * MjjBATAnalysisFacility::PerformEnsembleTest(TTree * tree, int nensembles, int firstBinToUse, int lastBinToUse, int start)
{
   BCLog::OutSummary("Running ensemble test.");
   if (fFlagMCMC) {
      BCLog::OutSummary("Fit for each ensemble is going to be run using MCMC. It can take a while.");
   }

   // set log level
   // It would be better to set the level for the screen and for the file
   // separately. Perhaps in the future.
   BCLog::LogLevel lls = BCLog::GetLogLevelScreen();
   BCLog::LogLevel llf = BCLog::GetLogLevelFile();
   if(fLogLevel==BCLog::nothing) {
      BCLog::OutSummary("No log messages for the ensemble fits are going to be printed.");
      BCLog::SetLogLevel(fLogLevel);
   }
   else if(fLogLevel!=lls) {
      BCLog::OutSummary(Form("The log level for the ensemble test is set to \'%s\'.",BCLog::ToString(fLogLevel)));
      BCLog::SetLogLevel(fLogLevel);
   }

   // define histogram for the original data set
   MjjHistogram * histogram_data;

   // prepare the tree
   // get number of bins
   int nbins = fMBM->GetData()->GetNbinsX();

   // create new matrix row
   std::vector<double> nbins_matrix(nbins);

   // loop over bins
   for (int ibin = 1; ibin <= nbins; ++ibin) {
      // create branches
      tree->SetBranchAddress(Form("bin_%i", ibin),
                                        &(nbins_matrix)[ibin-1]);
   }

   // get number of parameters
   int nparameters = fMBM->GetNParameters();

   // define tree variables
   std::vector<double> out_parameters(nparameters);

   for (int i = 0; i < nparameters; ++i) {
      tree->SetBranchAddress(Form("parameter_%i", i), &out_parameters[i]);
   }

   // copy the original data set

   // set data pointer
   histogram_data = new MjjHistogram(fMBM->GetMjjHistogram());

   // create output tree
   TTree * tree_out = new TTree("ensemble_test", "ensemble test");

   // define tree variables
   std::vector<double> out_mode_global(nparameters);
   std::vector<double> out_std_global(nparameters);
   std::vector<double> out_mode_marginalized(nparameters);
   std::vector<double> out_mean_marginalized(nparameters);
   std::vector<double> out_median_marginalized(nparameters);
   std::vector<double> out_5quantile_marginalized(nparameters);
   std::vector<double> out_10quantile_marginalized(nparameters);
   std::vector<double> out_16quantile_marginalized(nparameters);
   std::vector<double> out_84quantile_marginalized(nparameters);
   std::vector<double> out_90quantile_marginalized(nparameters);
   std::vector<double> out_95quantile_marginalized(nparameters);
   std::vector<double> out_std_marginalized(nparameters);
   double out_chi2_generated;
   double out_chi2_mode;
   double out_cash_generated;
   double out_cash_mode;
   int out_nevents;

   // create branches
   for (int i = 0; i < nparameters; ++i) {
      tree_out->Branch(Form("parameter_%i", i), &out_parameters[i], Form("parameter %i/D", i));
      tree_out->Branch(Form("mode_global_%i", i), &out_mode_global[i], Form("global mode of par. %i/D", i));
      tree_out->Branch(Form("std_global_%i", i), &out_std_global[i], Form("global std of par. %i/D", i));
      if (fFlagMCMC) {
         tree_out->Branch(Form("mode_marginalized_%i", i), &out_mode_marginalized[i], Form("marginalized mode of par. %i/D", i));
         tree_out->Branch(Form("mean_marginalized_%i", i), &out_mean_marginalized[i], Form("marginalized mean of par. %i/D", i));
         tree_out->Branch(Form("median_marginalized_%i", i), &out_median_marginalized[i], Form("median of par. %i/D", i));
         tree_out->Branch(Form("5quantile_marginalized_%i", i), &out_5quantile_marginalized[i], Form("marginalized 5 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("10quantile_marginalized_%i", i), &out_10quantile_marginalized[i], Form("marginalized 10 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("16quantile_marginalized_%i", i), &out_16quantile_marginalized[i], Form("marginalized 16 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("84quantile_marginalized_%i", i), &out_84quantile_marginalized[i], Form("marginalized 84 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("90quantile_marginalized_%i", i), &out_90quantile_marginalized[i], Form("marginalized 90 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("95quantile_marginalized_%i", i), &out_95quantile_marginalized[i], Form("marginalized 95 per cent quantile of par. %i/D", i));
         tree_out->Branch(Form("std_marginalized_%i", i), &out_std_marginalized[i], Form("marginalized std of par. %i/D", i));
      }
   }
   tree_out->Branch("chi2_generated", &out_chi2_generated, "chi2 (generated par.)/D");
   tree_out->Branch("chi2_mode", &out_chi2_mode, "chi2 (mode of par.)/D");
   tree_out->Branch("cash_generated", &out_cash_generated, "cash statistic (generated par.)/D");
   tree_out->Branch("cash_mode", &out_cash_mode, "cash statistic (mode of par.)/D");
   tree_out->Branch("nevents", &out_nevents, "total number of events/I");

   // loop over ensembles
   for (int iensemble = 0; iensemble < nensembles; ++iensemble) {
      // print status
      if ((iensemble+1)%100 == 0 && iensemble > 0) {
         BCLog::SetLogLevel(lls,llf);
         int frac = double(iensemble+1) / double(nensembles) * 100.;
         BCLog::OutDetail(Form("Fraction of ensembles analyzed: %i%%",frac));
         BCLog::SetLogLevel(fLogLevel);
      }

      // get next (commented out: random) event from tree
      //      int index = (int) fRandom->Uniform(tree->GetEntries());
      int index = iensemble + start;
      tree->GetEntry(index);

      // transform matrix into histograms
      TH1D histogram = MatrixToHistograms(nbins_matrix);

      // Make MjjHistogram out of it
      MjjHistogram pseudodata(&histogram);

      // set data histogram
      fMBM->SetData(pseudodata,firstBinToUse,lastBinToUse);

      // check if MCMC should be run and perform analysis
      if (fFlagMCMC) {
         BCLog::SetLogLevel(lls,llf);
         BCLog::OutDetail(Form("Running MCMC for ensemble %i",iensemble));
         BCLog::SetLogLevel(fLogLevel);

         // work-around: force initialization
         fMBM->ResetResults();

         // run mcmc
         fMBM->MarginalizeAll();

         // find mode
         fMBM->FindMode( fMBM->GetBestFitParameters() );
      }
      else {
         // find mode
         fMBM->FindMode();
      }

      // fill tree variables
      out_mode_global = fMBM->GetBestFitParameters();
      out_std_global = fMBM->GetBestFitParameterErrors();

      for (int i = 0; i < nparameters; ++i) {
         if (fFlagMCMC) {
            BCH1D * hist = fMBM->GetMarginalized( fMBM->GetParameter(i) );
            out_mode_marginalized[i] = hist->GetMode();
            out_mean_marginalized[i] = hist->GetMean();
            out_median_marginalized[i] = hist->GetMedian();
            out_5quantile_marginalized[i] = hist->GetQuantile(0.05);
            out_10quantile_marginalized[i] = hist->GetQuantile(0.10);
            out_16quantile_marginalized[i] = hist->GetQuantile(0.16);
            out_84quantile_marginalized[i] = hist->GetQuantile(0.84);
            out_90quantile_marginalized[i] = hist->GetQuantile(0.90);
            out_95quantile_marginalized[i] = hist->GetQuantile(0.95);
            out_std_marginalized[i]=hist->GetRMS();
            std::cout << "Calculated 95th quantile is " << hist->GetQuantile(0.95) << std::endl;
         }
      }

      // Print out results of PE marginalisation for monitoring
      if (fFlagMCMC) {
        for (int i = 0; i < nparameters; ++i) {
          std::cout << "Selected value of parameter " << fMBM->GetParameter(i)->GetName() << " is " << fMBM->GetBestFitParameter(i) << std::endl;
        }
        int sigparam = fMBM->GetParIndicesProcess(fMBM->GetProcessIndex("SIGNAL")).at(0); // signal has just one param
        std::cout << "95th percentile of signal is " << out_95quantile_marginalized[sigparam] << std::endl;
      }

      // get number of events
      out_nevents = (int) fMBM->GetData()->Integral();

      // calculate chi2
      out_chi2_generated = fMBM->CalculateChi2( out_parameters );
      out_chi2_mode = fMBM->CalculateChi2( fMBM->GetBestFitParameters() );

      // calculate cash statistic
      out_cash_generated = fMBM->CalculateCash( out_parameters );
      out_cash_mode = fMBM->CalculateCash( fMBM->GetBestFitParameters() );

      // fill tree
      tree_out->Fill();
   }

   // put the original data back in place
   fMBM->SetData(*histogram_data,firstBinToUse,lastBinToUse);

   // reset log level
   BCLog::SetLogLevel(lls,llf);

   // work-around: force initialization
   fMBM->ResetResults();

   BCLog::OutSummary("Ensemble test ran successfully.");

   // return output tree
   return tree_out;
}

// ---------------------------------------------------------
TH1D MjjBATAnalysisFacility::MatrixToHistograms(std::vector<double> & matrix)
{

   // create new histogram
   TH1D hist( *(fMBM->GetData()) );

   // get number of bins
   int nbins = hist.GetNbinsX();

   // fill bin content
   for (int ibin = 1; ibin <= nbins; ++ibin) {
      hist.SetBinContent(ibin, matrix.at(ibin-1));
   }

   return hist;

}

// ---------------------------------------------------------
int MjjBATAnalysisFacility::PerformSingleSystematicAnalyses(std::string dirname, std::string options)
{
   BCLog::OutSummary(Form("Running systematic analysis in directory \'%s\'.",dirname.c_str()));

   // ---- create new directory ---- //

   mkdir(dirname.c_str(), 0777);
   chdir(dirname.c_str());

   // ---- check options ---- //

   bool flag_mcmc = true;

   if (std::string(options).find("mcmc") < std::string(options).size())
      flag_mcmc = true;

   // get number of systematics
   int nsystematics = fMBM->GetNSystematics();

   // container of flags for systematics
   std::vector<bool> flag_systematic(nsystematics);

   // create new container of comparison tools
   std::vector<MjjBATComparisonTool *> ctc;

   // get number of parameters
   int nparameters = fMBM->GetNParameters();

   // ---- add one comparison tool for each systematic ---- //

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // create new comparison tool
      MjjBATComparisonTool * ct = new MjjBATComparisonTool(fMBM->GetParameter(i)->GetName().c_str());

      // add comparison tool to container
      ctc.push_back(ct);
   }

   // ---- switch on all systematics ---- //

   for (int isystematic = 0; isystematic < nsystematics; ++ isystematic) {
      // get systematic
      MjjBATSystematic * systematic = fMBM->GetSystematic(isystematic);

      // remember old setting
      flag_systematic[isystematic] = systematic->GetFlagSystematicActive();

      // switch on
      systematic->SetFlagSystematicActive(true);
   }

   if (flag_mcmc) {
      // work-around: force initialization
      fMBM->ResetResults();

      // run mcmc
      fMBM->MarginalizeAll();

      // find mode
      fMBM->FindMode( fMBM->GetBestFitParameters() );
   }
   else {
      // find mode
      fMBM->FindMode();
   }

   // print results
   if (flag_mcmc)
      fMBM->PrintAllMarginalized("marginalized_all.ps");
   fMBM->PrintResults("results_all.txt");

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // get comparison tool
      MjjBATComparisonTool * ct = ctc.at(i);

      ct->AddContribution("all systematics",
                                    fMBM->GetBestFitParameters().at(i),
                                    fMBM->GetBestFitParameterErrors().at(i));
   }

   // ---- switch off all systematics ---- //

   // loop over all systematics
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      MjjBATSystematic * systematic = fMBM->GetSystematic(isystematic);

      // switch off
      systematic->SetFlagSystematicActive(false);
   }

   // ---- perform analysis with all systematics separately ---- //

   // loop over all systematics
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      MjjBATSystematic * systematic = fMBM->GetSystematic(isystematic);

      // switch on systematic
      systematic->SetFlagSystematicActive(true);

      // perform analysis
      if (flag_mcmc) {
         // work-around: force initialization
         fMBM->ResetResults();

         // run mcmc
         fMBM->MarginalizeAll();

         // find mode
         fMBM->FindMode( fMBM->GetBestFitParameters() );
      }
      else {
         // find mode
         fMBM->FindMode();
      }

      // print results
      if (flag_mcmc)
         fMBM->PrintAllMarginalized(Form("marginalized_systematic_%i.ps", isystematic));
      fMBM->PrintResults(Form("results_systematic_%i.txt", isystematic));

      // ---- update comparison tools ---- //

      // loop over all parameters
      for (int i = 0; i < nparameters; ++ i) {
         // get comparison tool
         MjjBATComparisonTool * ct = ctc.at(i);

         ct->AddContribution(systematic->GetName().c_str(),
                                       fMBM->GetBestFitParameters().at(i),
                                       fMBM->GetBestFitParameterErrors().at(i));
      }

      // switch off systematic
      systematic->SetFlagSystematicActive(false);
   }

   // ---- analysis without any systematic uncertainty ---- //

   if (flag_mcmc) {
      // work-around: force initialization
      fMBM->ResetResults();

      // run mcmc
      fMBM->MarginalizeAll();

      // find mode
      fMBM->FindMode( fMBM->GetBestFitParameters() );
   }
   else {
      // find mode
      fMBM->FindMode();
   }

   // print results
   if (flag_mcmc)
      fMBM->PrintAllMarginalized("marginalized_none.ps");
   fMBM->PrintResults("results_none.txt");

   // loop over all parameters
   for (int i = 0; i < nparameters; ++ i) {
      // get comparison tool
      MjjBATComparisonTool * ct = ctc.at(i);

      ct->AddContribution("no systematics",
                                    fMBM->GetBestFitParameters().at(i),
                                    fMBM->GetBestFitParameterErrors().at(i));
   }

   // ---- reset all systematics ---- //

   // loop over all systematics
   for (int isystematic = 0; isystematic < nsystematics; ++isystematic) {
      // get systematic
      MjjBATSystematic * systematic = fMBM->GetSystematic(isystematic);

      // switch off systematic
      if (flag_systematic[isystematic])
         systematic->SetFlagSystematicActive(flag_systematic[isystematic]);
   }

   // ---- workaround: reset MCMC ---- //
   fMBM->ResetResults();

   // ---- print everything ---- //
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // draw first one
   MjjBATComparisonTool * ct =  ctc.at(0);
   ct->DrawOverview();
   //   c1->Print((std::string(filename)+std::string("(")).c_str());
   c1->Print((std::string("overview.ps")+std::string("(")).c_str());

   // loop over all parameters
   for (int i = 1; i < nparameters-1; ++i) {
      // create new comparison tool
      MjjBATComparisonTool * ct = ctc.at(i);

      ct->DrawOverview();
      c1->Print("overview.ps");
   }

   // draw last one
   ct =  ctc.at(nparameters-1);
   ct->DrawOverview();
   c1->Print((std::string("overview.ps")+std::string(")")).c_str());

   // ---- free memory ---- //
   for (int i = 0; i < nparameters; ++i) {
      MjjBATComparisonTool * ct = ctc[i];
      delete ct;
   }
   ctc.clear();

   delete c1;

   // ---- change directory ---- //

   chdir("../");

   BCLog::OutSummary("Systematic analysis ran successfully");

   // no error
   return 1;
}

// ---------------------------------------------------------
int MjjBATAnalysisFacility::PerformCalibrationAnalysis(std::string dirname, const std::vector<double> & default_parameters, int index, const std::vector<double> & parametervalues, int nensembles)
{
   BCLog::OutSummary(Form("Running calibration analysis in directory \'%s\'.",dirname.c_str()));

   // ---- create new directory ---- //

   mkdir(dirname.c_str(), 0777);
   chdir(dirname.c_str());

   // ---- loop over parameter values and perform analysis  ---- //

   int nvalues = int(parametervalues.size());
   for (int ivalue = 0; ivalue < nvalues; ++ivalue) {

      // open file
      TFile * file = new TFile(Form("ensemble_%i.root", ivalue), "RECREATE");
      file->cd();

      // set parameters
      std::vector<double> parameters = default_parameters;
      parameters[index] = parametervalues.at(ivalue);

      // create ensemble
      TTree * tree = PerformEnsembleTest(parameters, nensembles);

      // write tree
      tree->Write();

      // close file
      file->Close();

      // free memory
      delete file;
   }

   // ---- change directory ---- //

   chdir("../");

   BCLog::OutSummary("Calibration analysis ran successfully");

   // no error
   return 1;
}

// ---------------------------------------------------------
