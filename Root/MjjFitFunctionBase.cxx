// ---------------------------------------------------------

#include "Bayesian/MjjFitFunctionBase.h"

// ---------------------------------------------------------
vector<double> MjjFitFunction::GetCurrentParameterValues()
{

  vector<double> theParameters;
  theParameters.clear();
  for (int i=0; i<fNParameters; i++) {
    double par = fFitFunction->GetParameter(i);
    theParameters.push_back(par);
  }
  return theParameters;

}

// ---------------------------------------------------------
void MjjFitFunction::SetParameterDefaults(vector<double> newParameterDefaults) 
{

  if (newParameterDefaults.size() != fNParameters) {
    fFitFunction->Fatal("MjjFitFunction::setParameterDefaults",
         "Length of vector does not match number of function parameters!");
  }
  for (unsigned int i=0; i<fNParameters; i++) {
    fFitParameters.at(i)->SetParamDefault(newParameterDefaults.at(i));
  }

  return;

}

// ---------------------------------------------------------
void MjjFitFunction::RestoreParameterDefaults() 
{

  for (int i=0; i<fNParameters; i++) {
    double thisParamDefault = fFitParameters.at(i)->GetParamDefault();
    fFitFunction->SetParameter(i,thisParamDefault);
  }

  return;

}

// ---------------------------------------------------------
void MjjFitFunction::SetCurrentParameterValues(vector<double> desiredParameterValues) 
{

  if (desiredParameterValues.size() != fNParameters) {
    fFitFunction->Fatal("MjjFitFunction::setCurrentParameterValues",
         "Length of vector does not match number of function parameters!");
  }
  for (unsigned int i=0; i<fNParameters; i++) {
    fFitFunction->SetParameter(i,desiredParameterValues.at(i));
  }

  return;

}

// ---------------------------------------------------------
void MjjFitFunction::SetParameterLimitLow(int paramIndex, double limit) 
{

  fFitParameters.at(paramIndex)->SetParamLimitLow(limit);
  return;

}

// ---------------------------------------------------------
void MjjFitFunction::SetParameterLimitHigh(int paramIndex, double limit) 
{

  fFitParameters.at(paramIndex)->SetParamLimitHigh(limit);
  return;

}

// ---------------------------------------------------------
void MjjFitFunction::SetDoWindowExclusion(bool yesOrNo) 
{

  fPersonalisedFunction->SetDoWindowExclusion(yesOrNo);
  return;

}

// ---------------------------------------------------------
void MjjFitFunction::SetExclusionWindowFromRange(double xlow, double xhigh) 
{

  fPersonalisedFunction->SetExclusionWindowFromRange(xlow,xhigh);
  return;

}

// ---------------------------------------------------------
void MjjFitFunction::SetExclusionWindowFromHisto(MjjHistogram & signalTemplate, double percentage) 
{

  double lowLimit, highLimit;
  // Fit is always performed on basic histo, so that's the interval we need.
  std::pair<int, int> edges = signalTemplate.GetIntervalContainingPercentage(percentage);
  TH1D theHist = (TH1D) signalTemplate.GetHistogram();
  lowLimit = theHist.GetBinLowEdge(edges.first);
  highLimit = theHist.GetBinLowEdge(edges.second) + theHist.GetBinWidth(edges.second);
  SetExclusionWindowFromRange(lowLimit, highLimit);

  return;

}

// ---------------------------------------------------------
void MjjFitFunction::MakeHistFromFunction(TH1D* blankHist, double xmin, double xmax, bool overwrite)
{

  bool saveExc = fPersonalisedFunction->GetDoWindowExclusion();
  fPersonalisedFunction->SetDoWindowExclusion(false);

  if (overwrite) {
    for (int i=0; i<=blankHist->GetNbinsX(); i++) {
      blankHist->SetBinContent(i,0);
      blankHist->SetBinError(i,0);
    }
  }

  int firstBinToFill, lastBinToFill;
  if (xmin==-1 && xmax==-1) {
    xmin = blankHist->GetXaxis()->GetXmin();
    xmax = blankHist->GetXaxis()->GetXmax();
    if (xmin < fMjjLow) {
      xmin = fMjjLow;
      firstBinToFill = blankHist->FindBin(xmin+1);
    } else {firstBinToFill = 1;}
    if (xmax > fMjjHigh) {
      xmax = fMjjHigh;
      lastBinToFill = blankHist->FindBin(xmax-1);
    } else {lastBinToFill = blankHist->GetNbinsX();}
  }
  else {
    assert (xmin <= xmax);
    firstBinToFill = blankHist->FindBin(xmin+1);
    lastBinToFill = blankHist->FindBin(xmax-1);
  }
  for (int bin=firstBinToFill; bin<lastBinToFill+1; bin++) {
    double x1 = blankHist->GetBinLowEdge(bin);
    double x2 = x1 + blankHist->GetBinWidth(bin);
//    double eventsGenerated = fFitFunction->Eval(blankHist->GetBinCenter(bin))*blankHist->GetBinWidth(bin);
    double eventsGenerated = fFitFunction->Integral(x1,x2);
    if (isnan(eventsGenerated)) eventsGenerated = 0;
    blankHist->SetBinContent(bin,eventsGenerated);
    blankHist->SetBinError(bin,sqrt(eventsGenerated));
  }

  fPersonalisedFunction->SetDoWindowExclusion(saveExc);


  return;

}

// ---------------------------------------------------------

