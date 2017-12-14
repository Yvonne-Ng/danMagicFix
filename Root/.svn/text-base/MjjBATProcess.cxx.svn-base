// ---------------------------------------------------------

#include "Bayesian/MjjBATProcess.h"

// ---------------------------------------------------------
MjjBATProcess::MjjBATProcess(std::string name, bool doNormalisationUnc) //, bool doStatUnc)
 : fNominal(0)
 , fVariation(0)
 , fFunction(0)
 , fNParams(0)
{
   fName = name;
   fLuminosity = 1;
   templateVariation = NULL;
   fTrimPercentage = 0.95;
   fTrimProcess = false;
//   fUseStatUnc = doStatUnc;
   fUseNormUnc = doNormalisationUnc;
   fNormalisationParameter = -1;
}

// ---------------------------------------------------------
MjjBATProcess::~MjjBATProcess()
{
   delete templateVariation;
}

// ---------------------------------------------------------
double MjjBATProcess::GetBinValue(int bin, vector<double> params)
{

   if (!fUseNormUnc) return fNominalBins[bin];
   else {
   if (fFunction == 0) {
      return fNominalBins[bin] + params.at(fParamIndices.at(0)) * fVariationBins[bin];
   } else {
         double x = fNominal->GetBinCenter(bin);
         vector<double> funcpars;
         for (int i=0; i<fNParams; i++) {
            int index = fParamIndices.at(i);
            funcpars.push_back(fFuncParamsNominal.at(i) + fFuncParamsError.at(i) * params.at(index));
         }
         fFunction->SetCurrentParameterValues(funcpars);
         return fFunction->GetFitFunction()->Eval(x) * fNominal->GetBinWidth(bin);
      }
   }

}

// ---------------------------------------------------------
TVectorD MjjBATProcess::GetBins(vector<double> params) {

   if (!fUseNormUnc) {
      TVectorD bins(fNominalBins);
      return bins;
   } else {
      if (fFunction == 0) {
  
         TVectorD bins(fVariationBins);
         bins*=params.at(fParamIndices.at(0));
         bins+=fNominalBins;
         return bins;

      } else {

         TVectorD bins(fNBinsInSpectrum);
         vector<double> funcpars;
         for (int i=0; i<fNParams; i++) {
            int index = fParamIndices.at(i);
            funcpars.push_back(fFuncParamsNominal.at(i) + fFuncParamsError.at(i) * params.at(index));
         }
         fFunction->SetCurrentParameterValues(funcpars);
         for (int i=0; i<fNBinsInSpectrum; i++) {
            bins[i] = fFunction->GetFitFunction()->Eval(fNominal->GetBinCenter(i)) * fNominal->GetBinWidth(i);
         }
         return bins;
      }
   }
}

// ---------------------------------------------------------
double MjjBATProcess::GetNEvents(vector<double> params) {

   double nEvents=0;
   // Including under and overflow bins would be meaningless here.
   for (int bin=1; bin<fNominal->GetNbinsX()+1; bin++) {
     nEvents += this -> GetBinValue(bin,params);
   }
  
   return nEvents;

}

// ---------------------------------------------------------
void MjjBATProcess::SetTemplateFromHistograms(TH1D * centralTemplate, TH1D * variationTemplate)
{

   // Store histograms
   fNominal = centralTemplate;
   fNBinsInSpectrum = centralTemplate->GetNbinsX()+2;

   if (variationTemplate == 0) {
     if (fUseNormUnc) {
        std::cout << "You need to supply a variation template if you want to use a normalisation uncertainty!" << std::endl;
        exit(EXIT_FAILURE);
     }
   } else {
      fVariation = variationTemplate;
      fVariationBins.ResizeTo(fNBinsInSpectrum);
      for (int bin=0; bin<fNBinsInSpectrum; bin++) {
        fVariationBins[bin] = fVariation->GetBinContent(bin);
      }
   }

   // Store bins as TVectorD's for increased speed elsewhere
   fNominalBins.ResizeTo(fNBinsInSpectrum);
   for (int bin=0; bin<fNBinsInSpectrum; bin++) {
     fNominalBins[bin] = fNominal->GetBinContent(bin);
   }

   // Store number of parameters
   fNParams = 1;
   return;
}

// ---------------------------------------------------------
void MjjBATProcess::SetTemplateFromFunction(MjjFitFunction * functionTemplate, TH1D * blankHist, vector<double> nominalFuncParams, vector<double> errorFuncParams, double luminosity)
{

   // set up function
   fFunction = new MjjFitFunction(*functionTemplate);
   fFunction->SetCurrentParameterValues(nominalFuncParams);

   // store other information
   fNParams = functionTemplate->GetNParams();
   fFuncParamsNominal = nominalFuncParams;
   fFuncParamsError = errorFuncParams;
   fLuminosity = luminosity;

   assert (nominalFuncParams.size()==(unsigned int)fNParams);
   assert (errorFuncParams.size()==(unsigned int)fNParams);

   // set up nominal from function and nominal params
   fNominal = new TH1D(*blankHist);
   fNominal->SetName(Form("NominalHistFromFunc_%s",fName.c_str()));
   fFunction->MakeHistFromFunction(fNominal);

   return;
}

// ---------------------------------------------------------
int MjjBATProcess::SetSystematicShapeVariation(std::string systematic, MjjBATShapeChangingSyst * systvar) {

  // Guard against having both template and matrices
  if (templateVariation!=0) {
    std::cout << "You already set a SystematicTemplateVariation! You cannot use both." << std::endl;
    exit (EXIT_FAILURE);
  }

  systematics.push_back(systematic);
  shapeVariations.push_back(systvar);
  return 0;
}

// ---------------------------------------------------------
void MjjBATProcess::SetSystematicScaleVariation(std::string systematic, MjjBATScaleChangingSyst * systvar) {

  systematics.push_back(systematic);
  scaleVariations.push_back(systvar);
  return;
}

// ---------------------------------------------------------
int MjjBATProcess::SetSystematicTemplateVariation(std::string systematic, MjjBATTemplateSyst * systvar) {

  bool isScale = systvar->GetIsScale();

  if (!isScale) {

    // Guard against having both template and matrices
    if (shapeVariations.size()>0) {
      std::cout << "You already set a ShapeChangingVariation! You cannot use both." << std::endl;
      exit (EXIT_FAILURE);
    }
    if (templateVariation) {
      std::cout << "You already set a template uncertainty! You cannot use another one." << std::endl;
      exit (EXIT_FAILURE);
    }
    
    templateVariation = systvar;
  } else {
    templateScaleVariations.push_back(systvar);
  }

  systematics.push_back(systematic);
  return 0;
}

// ---------------------------------------------------------
