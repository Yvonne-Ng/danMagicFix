//still not implemented in the code
#ifndef __FoldingFunctions__
#define  __FoldingFunctions__

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

TFile *f_FM;
TH2D *h_FM_nominal;

void LoadFoldingMatrix(TString FileName)
{
  f_FM = new TFile(FileName);
  h_FM_nominal = (TH2D*)f_FM->Get("Nominal/TM_normalized_with_eff");
}

void Fold(TH1D *h_MJJtruth, TH2D *FM, TH1D *h_MJJreco)
{
  int Nbins_Mjj_truth = h_MJJtruth->GetXaxis()->GetNbins();
  int Nbins_Mjj_reco = h_MJJreco->GetXaxis()->GetNbins();
  int Nbins_FM_X = FM->GetXaxis()->GetNbins();
  int Nbins_FM_Y = FM->GetYaxis()->GetNbins();
  int diff = FM->GetYaxis()->FindBin(5000)-h_MJJreco->GetXaxis()->FindBin(5000);
  std::cout<<"\ndifference between the begin of binning of Mjj and FM = "<<diff<<"\n\n";
  double recoFold;
  for (int iMjj=1; iMjj< Nbins_Mjj_reco+1; iMjj++)
  {
    recoFold = 0;
    for (int iTruth=1; iTruth<Nbins_FM_X+1; iTruth++)
      recoFold += h_MJJtruth->GetBinContent(iTruth)*FM->GetBinContent(iTruth,iMjj+diff);

    if(recoFold<1E-5) h_MJJreco->SetBinContent(iMjj, 0);
    else h_MJJreco->SetBinContent(iMjj, recoFold);
    h_MJJreco->SetBinError(iMjj,0.);
  }
  /*std::cout << "MJJtruth" << '\n';
  h_MJJtruth->Print("all");
  std::cout << "MJJreco" << '\n';
  h_MJJreco->Print("all");*/
}

/*void Fold(TH1D *h_MJJtruth, TH2D *FM, double *recoFold)
{
  int Nbins_Mjj = h_MJJtruth->GetXaxis()->GetNbins();
  int Nbins_FM = FM->GetXaxis()->GetNbins();
  int diff = FM->GetXaxis()->FindBin(5000)-h_MJJtruth->GetXaxis()->FindBin(5000);
  std::cout<<"\ndifference between binning of Mjj and FM = "<<diff<<"\n\n";
  for (int iMjj=1; iMjj< Nbins_Mjj+1; iMjj++)
  {
    recoFold[iMjj-1] = 0;
    for (int iTruth=1; iTruth<Nbins_FM+1; iTruth++)
      recoFold[iMjj-1] += h_MJJtruth->GetBinContent(iTruth)*FM->GetBinContent(iTruth,iMjj+diff);
  }
}*/

void FillArrayToHisto(TH1D *histo, double *array)
{
  int Nbins = histo->GetXaxis()->GetNbins();
  for (int iBin=1; iBin< Nbins+1; iBin++)
  {
    if(array[iBin-1]<1E-5) histo->SetBinContent(iBin, 0);
    else histo->SetBinContent(iBin, array[iBin-1]);
    histo->SetBinError(iBin,0.);
  }
}

void fillDirac(TH1D *histo, double shiftMass)
{
  int Nbins = histo->GetXaxis()->GetNbins();
  for (int iBin=1; iBin <= Nbins; iBin++)
  {
    histo->SetBinContent(iBin, 0);
  }
  int diracBin = histo->GetXaxis()->FindBin(shiftMass);
  if(shiftMass == histo->GetXaxis()->GetBinLowEdge(diracBin))
  {
    histo->SetBinContent(diracBin-1,0.5);
    histo->SetBinContent(diracBin,0.5);
  }
  else histo->SetBinContent(diracBin,1);
}

#endif //__FoldingFunctions__
