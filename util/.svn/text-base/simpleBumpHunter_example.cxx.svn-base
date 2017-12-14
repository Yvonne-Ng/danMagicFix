#include "TH1.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include "Bayesian/MjjHistogram.h"
#include "Bayesian/MjjStatisticalTest.h"
#include "Bayesian/MjjPseudoExperimenter.h"
#include "Bayesian/MjjStatisticsBundle.h"
#include "Bayesian/MjjSignificanceTests.h"
#include "Bayesian/MjjBumpHunter.h"

#include "TEnv.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVector.h"
#include "TStopwatch.h"
#include "TLine.h"
#include "TLegend.h"

int main (int argc,char **argv)
{

	//------------------------------------------

	// Start counting time
	TStopwatch totaltime;
	totaltime.Start();

	float minBHMass = -1;
	float maxBHMass = -1;
	TString inputFileName;
	TString outputPlotName = "BumpHunterResult";
	TString dataHistoName;
    TString bkgHistoName;
    int nPseudoExpBH = 1e3;
  
  	//------------------------------------------

	// Start reading input configuration
	int ip=1;
	while (ip<argc) {

		if (string(argv[ip]).substr(0,2)=="--") {

			//input file
			if (string(argv[ip])=="--file") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					inputFileName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno input file name inserted"<<std::endl; break;}
			}

			//output file
			else if (string(argv[ip])=="--output") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					outputPlotName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno name of plot to be created not included"<<std::endl; break;}
			}

			//data histogram name
			else if (string(argv[ip])=="--dataHist") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					dataHistoName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno histogram name inserted"<<std::endl; break;}
			}

			//background prediction histogram name
			else if (string(argv[ip])=="--bkgHist") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					bkgHistoName=argv[ip+1];
					ip+=2;
				} else {std::cout<<"\nno histogram name inserted"<<std::endl; break;}
			}
          
			//minimum mass for Bump Hunter, if you don't want to examine all bins
			else if (string(argv[ip])=="--minBH") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					minBHMass=std::stof(string(argv[ip+1]));
					ip+=2;
				} else {std::cout<<"\nNo BumpHunter minimum value given "<<std::endl; break;}
			}

			//maximum mass for Bump Hunter, if you don't want to examine all bins
			else if (string(argv[ip])=="--maxBH") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					maxBHMass=std::stof(string(argv[ip+1]));
					ip+=2;
				} else {std::cout<<"\nNo BumpHunter maximum value given "<<std::endl; break;}
			}

			//number of pseudoexperiments to use in BH (default 1000)
			else if (string(argv[ip])=="--nPseudoExpBH") {
				if (ip+1<argc && string(argv[ip+1]).substr(0,2)!="--") {
					nPseudoExpBH=std::stoi(string(argv[ip+1]));
					ip+=2;
				} else {std::cout<<"\nNumber of pseudoexperiments not specified "<<std::endl; break;}
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

	//------------------------------------------

    // Retrieve histograms and set up for use

	// Open files
	TFile * infile = TFile::Open(inputFileName,"READ");
	TH1::AddDirectory(kFALSE);

	// Get and store histograms
	TH1D * rawDataHisto = (TH1D*) infile->Get(dataHistoName);
    TH1D * rawBkgHisto = (TH1D*) infile->Get(bkgHistoName);

    // Make histogram wrapper class that the code requires.
    // The "false" in the background histogram constructor just
    // tells the code not to assume its effective number of events
    // is related to the error bars (sorry for inconvenience)
    MjjHistogram dataHistogram(rawDataHisto);
    MjjHistogram backgroundHistogram(rawBkgHisto,false);

  
	//------------------------------------------

    // Set up range to use. You may not always want to bump-hunt
    // the entire set of filled bins (say if your background is poorly modeled
    // below some point.
    // This script will use the first and last bins with data to define the range
    // unless instructed otherwise.

	int firstBin = dataHistogram.GetFirstBinWithData();
    int lastBin = dataHistogram.GetLastBinWithData();
    int firstBinBH; int lastBinBH;

	if ( minBHMass == -1) firstBinBH = firstBin;
	else firstBinBH = rawDataHisto->FindBin(minBHMass);
	if ( maxBHMass == -1) lastBinBH = lastBin;
	else lastBinBH = rawDataHisto->FindBin(maxBHMass);

    std::cout << "Will bump hunt the spectrum in bins [" << firstBinBH << " - " << lastBinBH << "]\n"
          << "\tcorresponding to a range [" << rawDataHisto->GetBinLowEdge(firstBinBH) << " - "
          << rawDataHisto->GetBinLowEdge(lastBinBH)+rawDataHisto->GetBinWidth(lastBinBH) << "]" << std::endl;

	//------------------------------------------

    // Set up tools for calculation.

	// Create bump hunter.
	MjjBumpHunter theBumpHunter;
  
    // We will not consider one-bin bumps. Change this if
    // you want tigher or stricter minimum widths, but 2 is usually fine.
	theBumpHunter.SetMinBumpWidth(2);
  
    // The default maximum bin width is 1/2 the number of bins in the spectrum.
    // If you want another maximum, uncomment this and set it like so.
    // theBumpHunter.SetMaxBumpWidth(15)
  
    // The sideband option is rarely used in ATLAS but has been
    // included in BH implementations for a long time.
    // The default in the current code is "false" but you can
    // control it here
	theBumpHunter.SetUseSidebands(false);
  
    // The standard operating mode of the BumpHunter is only interested in
    // excesses, since this is how we expect new physics to appear. If you
    // would like to consider both excesses and deficits when finding the
    // region of greatest discrepancy, activate this:
   	// theBumpHunter.AllowDeficit(true);

  
    // The BumpHunter class contains a function you can call like so to
    // compute the test statistic comparing any data histogram to a background estimate:
    // theBumpHunter.DoTest(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse=-1, int lastBinToUse=-1)
    // However, the test statistic is not the important quantity for the BumpHunter.
    // We need to compute a test statistic for each of a range of pseudoexperiments
    // and find the proportion with a more extreme value than that found in our data.

	// Make a pseudoexperimenter. I keep a toy-generating engine
    // as a separate class because I use it for a range of statistical tests
	MjjPseudoExperimenter thePseudinator;

	// Obtain estimate of BumpHunter p-value.
    // The MjjStatisticsBundle is just a struct containing the test statistic in data
    // and the test statistic in all of the histograms (they will be extracted below).
    // The higher the number of pseudoexperiments used, the more precisely known
    // the final p-value will be.
    // Sorry about putting the prediction first and observation second: it's counterintuitive
    // but historical at this point.
    MjjStatisticsBundle myStats = thePseudinator.GetPseudoExperimentStatsOnHistogram
		(backgroundHistogram,dataHistogram,&theBumpHunter,firstBinBH,lastBinBH,nPseudoExpBH);

    // This is the test statistic in data
	double bumpHunterStat = myStats.originalStatistic;
  
  	// This uses the test statistics in pseudoexperiments and the value in data to
    // compute a p-value (with stat uncertainties)
	std::pair<double,double> bumpHunterPValAndErr = GetFrequentistPValAndError
		(myStats.statisticsFromPseudoexperiments,myStats.originalStatistic);
  
    // The BumpHunter stored the edges of the bins corresponding to the region of
    // greatest excess. If you ran only DoTest on a single spectrum, you
    // could retrieve it like so:
    // vector<double> bumpEdges = theBumpHunter.GetFurtherInformation();
    // double lowEdgeOfBump = bumpEdges.at(0);
    // double highEdgeOfBump = bumpEdges.at(1);
    // However, this contains the information corresponding to the last spectrum
    // examined (in this case, one of the pseudoexperiments).
    // The pseudoexperimenter saved the results in data, so we can access it like this:
    double lowEdgeOfBump = myStats.originalFurtherInformation.at(0);
	double highEdgeOfBump = myStats.originalFurtherInformation.at(1);
  
	std::cout << "BumpHunter results: stat = " << bumpHunterStat << std::endl;
  
	// Make histo of residual: we'll put this in the plot because it's helpful
    // for visualisation
    MjjSignificanceTests theTestMaker;
	TH1D residualHist = theTestMaker.GetResidual(dataHistogram, backgroundHistogram, firstBinBH, lastBinBH);

	//------------------------------------------
	// Print BH quantities for inspection/checking
	std::cout << "******************************************" << std::endl;
	std::cout << "*** Final values " << std::endl;
	std::cout << "*** BH p-value = " << bumpHunterPValAndErr.first << std::endl;
	std::cout << "*** BH test statistic value = " << myStats.originalStatistic << std::endl;
	std::cout << "*** Selected most discrepant range = " << lowEdgeOfBump << " - " << highEdgeOfBump << std::endl;
	std::cout << "******************************************" << std::endl;
	//------------------------------------------

    // Make a plot! This uses an approximation of the usual dijet style, since plotting
    // is handled in the analyses by a standalone Python library.

    //outputname = name+epsorpdf
    TCanvas c("plot_canvas","plot_canvas",0,0,600,600);
    c.SetLogx(1);
    c.SetLogy(true);
    c.SetGridx(0);
    c.SetGridy(0);

    // Dimensions: xlow, ylow, xup, yup
    TPad outpad("extpad","extpad",0,0,1,1); // For marking outermost dimensions
    TPad pad1("pad1","pad1",0,0.27,1,1); // For main histo
    TPad pad2("pad2","pad2",0,0,1,0.27); // For residuals histo

    // Set up to draw in right orientations
    outpad.SetFillStyle(4000); // transparent
    pad1.SetBottomMargin(0.00001);
    pad1.SetTopMargin(0.02);
    pad1.SetRightMargin(0.02);
    pad1.SetBorderMode(0);
    pad1.SetLogy(true);
    pad1.SetLogx(true);
    pad2.SetTopMargin(0.00001);
    pad2.SetBottomMargin(0.3);
    pad2.SetRightMargin(0.02);
    pad2.SetBorderMode(0);
    pad2.SetLogx(true);
    pad1.Draw();
    pad2.Draw();
    outpad.Draw();

    // Set up to draw data and background histograms
    pad1.cd();

    rawDataHisto->GetYaxis()->SetTitleSize(0.05);
    rawDataHisto->GetYaxis()->SetTitleOffset(1.0);
    rawDataHisto->GetYaxis()->SetLabelSize(0.05);
  
    // Draw data histogram
    rawDataHisto->SetMarkerColor(kBlack);
    rawDataHisto->SetLineColor(kBlack);
    rawDataHisto->SetTitle("");
    rawDataHisto->GetXaxis()->SetTitle("Mass");
    rawDataHisto->GetYaxis()->SetTitle("Entries");
    rawDataHisto->GetXaxis()->SetRange(firstBin,lastBin);
    rawDataHisto->SetMarkerSize(0.75);
    rawDataHisto->GetYaxis()->SetRangeUser(0.3,rawDataHisto->GetBinContent(rawDataHisto->GetMaximumBin())*5);
    rawDataHisto->GetYaxis()->SetNdivisions(805,kTRUE);
    rawDataHisto->GetXaxis()->SetMoreLogLabels(kTRUE);
    rawDataHisto->Draw("E");

    // Draw background histogram
    rawBkgHisto->SetLineColor(kRed);
    rawBkgHisto->SetLineStyle(1);
    rawBkgHisto->SetFillStyle(0);
    rawBkgHisto->SetLineWidth(2);
    rawBkgHisto->SetTitle("");
    //rawBkgHisto->GetXaxis()->SetTitle(xname);
    //rawBkgHisto->GetYaxis()->SetTitle(yname);
    rawBkgHisto->GetXaxis()->SetRange(firstBin,lastBin);
    rawBkgHisto->GetXaxis()->SetNdivisions(805,kTRUE);
    rawBkgHisto->Draw("HIST L ][ SAME");

    // Draw residual histogram
    pad2.cd();
    residualHist.GetYaxis()->SetTitleSize(0.1);
    residualHist.GetYaxis()->SetTitleOffset(0.42); // 1.2 = 20% larger
    residualHist.GetYaxis()->SetLabelSize(0.1);
    residualHist.GetXaxis()->SetLabelSize(0.1);
    residualHist.GetXaxis()->SetTitleSize(0.1);
    residualHist.GetXaxis()->SetTitleOffset(1.2);
    residualHist.SetLineColor(kBlack);
    residualHist.SetLineWidth(2);
    residualHist.SetFillColor(kRed);
    residualHist.SetFillStyle(1001);
    residualHist.SetTitle("");
    residualHist.GetXaxis()->SetRange(firstBin,lastBin);
    residualHist.GetYaxis()->SetRangeUser(-3.7,3.7);
    residualHist.GetYaxis()->SetTickLength(0.055);
    residualHist.GetXaxis()->SetNdivisions(805,kTRUE);
    residualHist.GetYaxis()->SetTitle("Significance");
    residualHist.GetXaxis()->SetTitle("Mass");
    residualHist.Draw("HIST");

    c.Update();

    // Draw blue lines

    double heightLowEdge = 0;
    double heightHighEdge = 0;
    double minYvalue = rawDataHisto->GetMinimum();
    for (int i = 0; i < rawDataHisto->GetNbinsX()+2; i++) {
        double locationOfTallEdge = rawDataHisto->GetBinLowEdge(i);
        if (fabs(locationOfTallEdge - lowEdgeOfBump) < 2*std::numeric_limits<double>::epsilon())
          heightLowEdge = rawDataHisto->GetBinContent(i);
        if (fabs(locationOfTallEdge - highEdgeOfBump) < 2*std::numeric_limits<double>::epsilon())
          heightHighEdge = rawDataHisto->GetBinContent(i-1);
    }

    //double lowYVal = significance.GetMinimum();
    //double highYVal = significance.GetMaximum();

    TLine line1(lowEdgeOfBump, minYvalue, lowEdgeOfBump, heightLowEdge);
    TLine line2(highEdgeOfBump, minYvalue, highEdgeOfBump, heightHighEdge);
    TLine line3(lowEdgeOfBump, -3.7, lowEdgeOfBump, 3.7);
    TLine line4(highEdgeOfBump, -3.7, highEdgeOfBump, 3.7);
  
    pad1.cd();
    line1.SetLineColor(kBlue);
    line1.Draw();
    line2.SetLineColor(kBlue);
    line2.Draw();
    pad2.cd();
    line3.SetLineColor(kBlue);
    line3.Draw();
    line4.SetLineColor(kBlue);
    line4.Draw();
    c.Update();

    // Draw legend
    outpad.cd();
    TLegend legend(0.48,0.75,0.9,0.9);
    legend.SetTextFont(42);
    legend.SetTextSize(0.04);
    legend.SetBorderSize(0);
    legend.SetLineColor(0);
    legend.SetLineStyle(1);
    legend.SetLineWidth(1);
    legend.SetFillColor(0);
    legend.SetFillStyle(0);
  
    outpad.cd();
    legend.AddEntry(rawDataHisto,"Data","LFP");
    legend.AddEntry(rawBkgHisto,"Background fit","LF");
    legend.AddEntry(&line4,"BumpHunter interval","L");
    legend.Draw();
    c.Update();
    
    pad1.RedrawAxis();
    pad2.RedrawAxis();
    c.RedrawAxis();
    c.Update();
    TString outname = outputPlotName+".eps";
    c.SaveAs(outname);

  
  	//------------------------------------------

    // Cleanup

	infile->Close();

	totaltime.Stop();
	std::cout << "Process ran in " << totaltime.CpuTime() << " seconds. " << std::endl;

    delete infile;

	return 0;

}
