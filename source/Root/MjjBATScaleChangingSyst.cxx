// ---------------------------------------------------------

#include "Bayesian/MjjBATScaleChangingSyst.h"

// ---------------------------------------------------------
MjjBATScaleChangingSyst::MjjBATScaleChangingSyst(double scale) {

  fBinObserved = 0;
  fScale = scale;

  // default setting uses bin content
  fUseBinContent = true;
  fUseBinError = false;

}

// ---------------------------------------------------------
MjjBATScaleChangingSyst::~MjjBATScaleChangingSyst() {}

// ---------------------------------------------------------
void MjjBATScaleChangingSyst::SetScaleFromBinContent() {

  fUseBinContent = true;
  fUseBinError = false;

}

// ---------------------------------------------------------
void MjjBATScaleChangingSyst::SetScaleFromBinError() {

  fUseBinError = true;
  fUseBinContent = false;

}

// ---------------------------------------------------------
TVectorD MjjBATScaleChangingSyst::GetAdjusted(TVectorD bins, double param) {

  double lumiUncertainty = param * fScale;
  if (fUseBinContent) {
    TVectorD result(bins);
    result*=(1+lumiUncertainty);
    return result;
  } else {
    TVectorD result(bins);
    result.Abs();
    result.Sqrt();
    result*=(1+lumiUncertainty);
    return result;
  }

}

// ---------------------------------------------------------
double MjjBATScaleChangingSyst::GetAdjustedBin(double bincontent, double param) {

  if (fUseBinContent) return bincontent*(1+fScale*param);
  else return sqrt(fabs(bincontent))*(1+fScale*param);

}

// ---------------------------------------------------------

