// ---------------------------------------------------------

#include "Bayesian/MjjBATModel.h"

// ---------------------------------------------------------
MjjBATModel::MjjBATModel()
 : BCModel("Multi-template Fitter")
 , fNProcesses(0)
 , fNSystematics(0)
 , fData(0)
 , fEffectiveData(0)
 , fWeights(0)
 , fBinStructure(0)
 , fHistUncertaintyBandExpectation(0)
 , fHistUncertaintyBandPoisson(0)
{ 
  fLikelihood = DEFAULT;
  fDebug = false;
}

// ---------------------------------------------------------
MjjBATModel::MjjBATModel(std::string name)
 : BCModel(name.c_str())
 , fNProcesses(0)
 , fNSystematics(0)
 , fData(0)
 , fEffectiveData(0)
 , fWeights(0)
 , fBinStructure(0)
 , fHistUncertaintyBandExpectation(0)
 , fHistUncertaintyBandPoisson(0)
{
  fLikelihood = DEFAULT;
  fDebug = false;
}

// ---------------------------------------------------------
MjjBATModel::~MjjBATModel()
   // default destructor
{ 

}

// ---------------------------------------------------------
int MjjBATModel::GetProcessIndex(std::string name)
{
   // loop over all processs and compare names
   for (int i = 0; i < fNProcesses; ++i) {
      // get process
      MjjBATProcess * process = GetProcess(i);

      // compare names
      if (!process->GetName().compare(name))
         return i;
   }

   // if process does not exist, return -1
   return -1;
}

// ---------------------------------------------------------
int MjjBATModel::GetSystematicIndex(std::string name)
{
   // loop over all systematics and compare names
   for (int i = 0; i < fNSystematics; ++i) {
      // get systematic
      MjjBATSystematic * systematic = GetSystematic(i);

      // compare names
      if (!systematic->GetName().compare(name))
         return i;
   }

   // if systematic does not exist, return -1
   return -1;
}

// ---------------------------------------------------------
bool MjjBATModel::DoBinsMatch(const TArrayD * firstBins, const TArrayD * secondBins) {

  if (firstBins->GetSize() != secondBins->GetSize()) return false;
  for (int i=0; i<firstBins->GetSize(); i++) {
    if (firstBins->At(i) != secondBins->At(i)) return false;
  }
  return true;

}

// ---------------------------------------------------------
bool MjjBATModel::DoBinsMatchStored(const TArrayD * firstBins) {

  const TArrayD * compareBins = new TArrayD(fBinStructure);
  bool result = DoBinsMatch(firstBins,compareBins);
  delete compareBins;
  return result;

}

// ---------------------------------------------------------
int MjjBATModel::SetTemplate(std::string processname, TH1D centralhist)
{

   // get process index
   int processindex = GetProcessIndex(processname);
   MjjBATProcess * process = GetProcess(processindex);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("MjjBATModel::SetTemplate() : Process does not exist.");
         return -1;
   }
  
   // check if process expects a variation histogram
   if (process->GetDoNormalisationUnc()) {
        std::cout << "MjjBATModel::SetTemplate(): This process requires a variation histogram!" << std::endl;
        exit(EXIT_FAILURE);
   }
  
   // check that binning of histograms match and match data, if set
   const TArrayD * centralbinstructure = centralhist.GetXaxis()->GetXbins();
   if (fBinStructure.GetSize()==0) {
      fBinStructure = *centralbinstructure;
   } else if (!(DoBinsMatch(centralbinstructure,&fBinStructure))) {
     fData->Fatal("MjjBATModel::SetTemplate()","Binning of template does not match some previously set histogram.");
     return -1;
   }

   // remove statistics box
   centralhist.SetStats(kFALSE);

   // set color and fill style
   centralhist.SetFillColor(2 + processindex);
   centralhist.SetFillStyle(1001);

   // create pointers to histograms
   TH1D * tempcentral = new TH1D(centralhist);

   // set histograms
   process->SetTemplateFromHistograms(tempcentral);

   // no error
   return 1;


}

// ---------------------------------------------------------
int MjjBATModel::SetTemplate(std::string processname, TH1D centralhist, TH1D variationhist)
{

   // get process index
   int processindex = GetProcessIndex(processname);
   MjjBATProcess * process = GetProcess(processindex);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("MjjBATModel::SetTemplate() : Process does not exist.");
         return -1;
   }

   // check that binning of histograms match and match data, if set
   const TArrayD * centralbinstructure = centralhist.GetXaxis()->GetXbins();
   const TArrayD * variationbinstructure = variationhist.GetXaxis()->GetXbins();
   if ((DoBinsMatch(centralbinstructure,variationbinstructure)) && (fBinStructure.GetSize()==0)) 
      fBinStructure = *centralbinstructure;
   else if (!(DoBinsMatch(centralbinstructure,variationbinstructure)) && (fBinStructure.GetSize()==0)) {
     BCLog::OutWarning("MjjBATModel::SetTemplate() : Binning of central template does not match that of variation template.");
     return -1;
   } else if (!(DoBinsMatch(centralbinstructure,&fBinStructure))) {
     fData->Fatal("MjjBATModel::SetTemplate()","Binning of template does not match some previously set histogram.");
     return -1;
   }

   // remove statistics box
   centralhist.SetStats(kFALSE);
   variationhist.SetStats(kFALSE);

   // set color and fill style
   centralhist.SetFillColor(2 + processindex);
   centralhist.SetFillStyle(1001);
   variationhist.SetFillStyle(1001);

   // create pointers to histograms
   TH1D * tempcentral = new TH1D(centralhist);
   TH1D * tempvar = new TH1D(variationhist);

   // set histograms
   process->SetTemplateFromHistograms(tempcentral,tempvar);

   // no error
   return 1;
}

// ---------------------------------------------------------
int MjjBATModel::SetTemplate(std::string processname, MjjFitFunction & function, TH1D templatehist, vector<double> central, vector<double> error, double luminosity)
{

   // get process index
   int processindex = GetProcessIndex(processname);
   MjjBATProcess * process = GetProcess(processindex);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("MjjBATModel::SetTemplate() : Process does not exist.");
         return -1;
   }

   // check that binning matches any previous histograms
   const TArrayD * binstructure = templatehist.GetXaxis()->GetXbins();
   if (fBinStructure.GetSize()==0) fBinStructure = *binstructure;
   else if (!(DoBinsMatchStored(binstructure))) {
     fData->Fatal("MjjBATModel::SetTemplate()","Binning of template does not match some previously set histogram.");
     return -1;
   }

   // remove statistics box
   templatehist.SetStats(kFALSE);

   // set color and fill style
   templatehist.SetFillColor(2 + processindex);
   templatehist.SetFillStyle(1001);

   // create new histogram
   TH1D * temp = new TH1D(templatehist);

   // create new function
   MjjFitFunction * fitfunction = new MjjFitFunction(function);

   // set histograms
   process->SetTemplateFromFunction(fitfunction,temp,central,error,luminosity);

   // no error
   return 1;
}


// ---------------------------------------------------------
int MjjBATModel::SetData(MjjHistogram mjjhist, int firstBinToUse, int lastBinToUse, double minimum, double maximum)
{

   TH1D hist = mjjhist.GetHistogram();
   hist.SetName("BATData");

   // If any signal histograms set yet, make sure they have the same bin structure.
   const TArrayD * binstructure = hist.GetXaxis()->GetXbins();
   if (fBinStructure.GetSize()==0) fBinStructure = *binstructure;
   else if (!(DoBinsMatchStored(binstructure))) {
     fData->Fatal("MjjBATModel::SetTemplate()","Binning of template does not match some previously set histogram.");
     return -1;
   }

   fMjjHist = mjjhist;

   // remove statistics box
   hist.SetStats(kFALSE);

   // set marker
   hist.SetMarkerStyle(20);
   hist.SetMarkerSize(1.1);

   // set divisions
   hist.SetNdivisions(509);

   // remove old data set if it exists
   if (fData) {
      delete fData;
      fData=0;
   }

   // TESTING
   if (fEffectiveData) {
      delete fEffectiveData;
      fEffectiveData=0;
   }
  if (fWeights) {
      delete fWeights;
      fWeights=0;
   }

   // remove old uncertainty histograms if they exist
   if (fHistUncertaintyBandExpectation) {
      delete fHistUncertaintyBandExpectation;
      fHistUncertaintyBandExpectation=0;
   }
   if (fHistUncertaintyBandPoisson) {
      delete fHistUncertaintyBandPoisson;
      fHistUncertaintyBandPoisson=0;
   }

   // create new histograms for uncertainty bands
   //	 double minimum = floor(TMath::Max(0., hist.GetMinimum() - 7.*sqrt(hist.GetMinimum())));
   if (minimum==-1)
      minimum = 0;
   if (maximum==-1)
      maximum = ceil(hist.GetMaximum() + 5.*sqrt(hist.GetMaximum()));
   int nbinsy = int(maximum-minimum);
   while (nbinsy>1000) nbinsy/=10;

   std::vector<double> a(hist.GetNbinsX()+1);
   for (int i = 0; i < hist.GetNbinsX()+1; ++i) {
      a[i] = hist.GetXaxis()->GetBinLowEdge(i+1);
   }

   TH2D* hist_uncbandexp = new TH2D(TString::Format("UncertaintyBandExpectation_%i", BCLog::GetHIndex()), "",											hist.GetNbinsX(), &a[0], 1000, minimum, maximum);
   hist_uncbandexp->SetStats(kFALSE);

   TH2D* hist_uncbandpoisson = new TH2D(TString::Format("UncertaintyBandPoisson_%i", BCLog::GetHIndex()), "",											hist.GetNbinsX(), &a[0], nbinsy, minimum, maximum);
   hist_uncbandpoisson->SetStats(kFALSE);

   // set histograms
   fData = new TH1D(hist);
   fEffectiveData = new TH1D(mjjhist.GetEffectiveHistogram());
   fWeights = new TH1D(mjjhist.GetWeightsHistogram());
   fHistUncertaintyBandExpectation = hist_uncbandexp; 
   fHistUncertaintyBandPoisson = hist_uncbandpoisson; 

   // set first and last bins to use from specified range or from hist
   fFirstBin = 1;
   while (fData->GetBinContent(fFirstBin)==0 && fFirstBin < fData->GetNbinsX()) {fFirstBin++; }
   fLastBin = fData->GetNbinsX();
   while (fData->GetBinContent(fLastBin)==0 && fLastBin > 1) {fLastBin--; }

   // safeguard against empty histograms
   if (fFirstBin==fData->GetNbinsX() && fLastBin == 1) {
      std::cout << "No data in histogram! Resetting limits to first and last bin." << std::endl;
      fFirstBin = 1; fLastBin = fData->GetNbinsX();
   }

   if (firstBinToUse>0 && firstBinToUse<fData->GetNbinsX()+1) fFirstBin = firstBinToUse;
   if(lastBinToUse>0 && /*lastBinToUse<fData->GetNbinsX()+1 && */ lastBinToUse > firstBinToUse) {
     if (lastBinToUse==(fData->GetNbinsX()+1)) 
       fLastBin = lastBinToUse-1; // overflow bin
     else fLastBin = lastBinToUse;
   }

   // set y-range for printing
   SetRangeY(minimum, maximum);

   // no error
   return 1;
}

// ---------------------------------------------------------
int MjjBATModel::AddProcess(std::string name, vector<double> min, vector<double> max)
{
   // check if process exists
   for (int i = 0; i < fNProcesses; ++i) {
      // compare names
      if (GetProcessIndex(name) >= 0) {
      BCLog::OutWarning("MjjBATModel::AddProcess() : Process with this name exists already.");
         return -1;
      }
   }

   // check if parameter limits are appropriate
   if (min.size() != max.size()) {
      BCLog::OutWarning("MjjBATModel::AddProcess() : Minimum and maximum sigma vectors must have same length.");
         return -1;
   }

   // create new process
   // this comes from a function so it should have a normalisation
   // uncertainty but no uncorrelated uncertainty.
   MjjBATProcess * process = new MjjBATProcess(name,true);

   // add process
   fProcessContainer.push_back(process);

   // increase number of processes
   fNProcesses++;

   // add all parameters
   for (unsigned int i=0; i<min.size(); i++) {

      // awkward string gymnastics because C++ is ridiculous
      std::ostringstream s;
      s << name << "_" << i;
      std::string paramname(s.str());

      // add parameter index to container
      process->AddProcessParameterIndex(GetNParameters());
      process->AddProcessParameterName(paramname.c_str());

      // add parameter
      AddParameter(paramname.c_str(), min.at(i), max.at(i));

   }

   // no error
   return 1;
}

// ---------------------------------------------------------
int MjjBATModel::AddProcess(std::string name, bool doUnc, double nmin, double nmax)
{

     // check if process exists
   for (int i = 0; i < fNProcesses; ++i) {
      // compare names
      if (GetProcessIndex(name) >= 0) {
      BCLog::OutWarning("MjjBATModel::AddProcess() : Process with this name exists already.");
         return -1;
      }
   }

   // create new process
   // this comes from a function so it should have a normalisation
   // uncertainty but no uncorrelated uncertainty.
   MjjBATProcess * process = new MjjBATProcess(name,doUnc);

   // add process
   fProcessContainer.push_back(process);

   // increase number of processes
   fNProcesses++;

   // add a parameter for the normalisation, if desired.
   if (doUnc) {

      // awkward string gymnastics because C++ is ridiculous
      std::ostringstream s;
      s << name << "_normalisation";
      std::string paramname(s.str());

      // add parameter index to container
      process->AddProcessParameterIndex(GetNParameters());
      process->AddProcessParameterName(paramname.c_str());
      process->AddNormalisationParIndex(GetNParameters());

      // add parameter
      AddParameter(paramname.c_str(), nmin, nmax);

   }

   // no error
   return 1;

}

// ---------------------------------------------------------
int MjjBATModel::AddSystematic(std::string name, double min, double max)
{
   // check if systematic exists
   for (int i = 0; i < fNSystematics; ++i) {
      // compare names
      if (GetSystematicIndex(name) >= 0) {
      BCLog::OutWarning("MjjBATModel::AddSystematic() : Systematic with this name exists already.");
         return -1;
      }
   }

   // create new systematic
   MjjBATSystematic * systematic = new MjjBATSystematic(name);

   // add systematic
   fSystematicContainer.push_back(systematic);

   // increase number of systematices
   fNSystematics++;

   // add parameter index to container
   systematic->SetSystematicParameterIndex(GetNParameters());

   // add parameter
   AddParameter(systematic->GetName().c_str(),min,max); 

   // no error
   return 1;
}

// ---------------------------------------------------------
int MjjBATModel::SetSystematicVariation(std::string processname,  std::string systematicname, MjjBATSystematicVariation * systVariation)
{

   // get process index
   int processindex = GetProcessIndex(processname);

   // check if process exists
   if (processindex < 0) {
      BCLog::OutWarning("MjjBATModel::SetSystematicVaration() : Process does not exist.");
         return -1;
   }

   // get systematic index
   int systematicindex = GetSystematicIndex(systematicname);

   // check if systematic exists
   if (systematicindex < 0) {
      BCLog::OutWarning("MjjBATModel::SetSystematicVariation() : Systematic does not exist.");
         return -1;
   }

   // get process
   MjjBATProcess * process = GetProcess(processindex);

   systVariation->SetParentSystematic(systematicname);

   MjjBATShapeChangingSyst * thisshapevar = dynamic_cast<MjjBATShapeChangingSyst*> (systVariation);
   MjjBATScaleChangingSyst * thisscalevar = dynamic_cast<MjjBATScaleChangingSyst*> (systVariation);
   MjjBATTemplateSyst * thistemplatevar = dynamic_cast<MjjBATTemplateSyst*> (systVariation);

   // Check binning where appropriate
   if (thisshapevar!=NULL || thistemplatevar!=NULL) {
     const TArrayD * compareBins = systVariation->GetBinStructure();
     if (!(DoBinsMatchStored(compareBins))) fData->Fatal("MjjBATModel::SetSystematicVariation",
                                      "Bins in systematic template/matrix do not match data or process");
   }

   if (thisshapevar!=NULL) {
     process->SetSystematicShapeVariation(systematicname,thisshapevar);
   } else if (thisscalevar!=NULL) {
     process->SetSystematicScaleVariation(systematicname,thisscalevar);
   } else {
     process->SetSystematicTemplateVariation(systematicname,thistemplatevar);
   }

   // no error
   return 1;
}

int MjjBATModel::SetInterestingBins(int firstBin, int lastBin)
{

  // Check to make sure these are not outside the scope of permitted bins
  if ((firstBin >= lastBin) || (firstBin < fFirstBin) || (lastBin > fLastBin)) {
    std::cout << "This bin range is not allowed!" << std::endl;
    exit (EXIT_FAILURE);

  }
  fFirstInterestingBin = firstBin;
  fLastInterestingBin = lastBin;
  
  return 1;

}

// ---------------------------------------------------------
int MjjBATModel::PrintSummary(std::string filename)
{
   // open file
   std::ofstream ofi(filename);
   ofi.precision(3);

   // check if file is open
   if(!ofi.is_open()) {
      BCLog::OutWarning(Form("MjjBATModel::PrintSummary() : Could not open file %s", filename.c_str()));
      return 0;
   }

   ofi
      << " Multi template fitter summary " << std::endl
      << " ----------------------------- " << std::endl
      << std::endl
      << " Number of processes     : " << fNProcesses << std::endl
      << " Number of systematics   : " << fNSystematics << std::endl
      << std::endl;

   ofi
      << " Processes :" << std::endl;
      for (int i = 0; i < GetNProcesses(); ++i) {
         ofi
            << " " << i
            << " : \"" << GetProcess(i)->GetName().c_str()  << "\""
            << std::endl;
      }
      ofi
         << std::endl;

   ofi
      << " Systematics :" << std::endl;
      for (int i = 0; i < GetNSystematics(); ++i) {
         ofi
            << " " << i
            << " : \"" << GetSystematic(i)->GetName().c_str()  << "\""
            << " (par index " << GetSystematic(i)->GetParamIndex() << ")"
            << std::endl;
      }
      ofi
         << std::endl;
      if (GetNSystematics() == 0)
         ofi
            << " - none - " << std::endl;

      ofi
         << " Goodness-of-fit: " 
         << CalculateChi2( GetBestFitParameters() )
         << std::endl;
      
      std::cout << std::endl;


   // close file
   ofi.close();

   // no error
   return 1;
}

// ---------------------------------------------------------
vector<double> MjjBATModel::Expectation(const std::vector<double> & parameters)
{

   int nbins = fData->GetNbinsX()+2;
   vector<vector<std::pair<double,double> > > expectation;

   // loop over all processes
   for (int i = 0; i < fNProcesses; ++i) {

     expectation.push_back(this->ProcessExpectation(i,parameters));

   }

   // check if expectation is positive
   double weight, binexpectation;
   vector<double> result;
   for (int bin=0; bin<nbins; bin++) {
     binexpectation = 0;
     weight = fWeights->GetBinContent(bin);
     for (unsigned int i=0; i<expectation.size(); i++) binexpectation+=expectation.at(i).at(bin).first;
     if (binexpectation < 0) result.push_back(0);
     else result.push_back(binexpectation/weight);
   }

   return result;
}

// ---------------------------------------------------------
std::vector<std::pair<double,double> > MjjBATModel::ProcessExpectation(int processindex, const std::vector<double> & parameters)
{
   int nbins = fData->GetNbinsX()+2;
   vector<std::pair<double,double> > expectation;


   // get expected nominal shape of process
   fProcessExpectation.ResizeTo(nbins);
   fProcessExpectation = GetProcess(processindex)->GetBins(parameters);
//   TH1D * nominalHist = GetProcess(processindex)->GetNominalHistogram();

   // get template-based scale-changing systematic variations for process
   vector<MjjBATTemplateSyst*> tempscalevars = GetProcess(processindex)->GetTempScaleChangingSysts();

   // loop over all scale-changing systematic variations
   for (unsigned int j=0; j<tempscalevars.size(); ++j) {

     MjjBATTemplateSyst * thisshapevar = tempscalevars.at(j);
     MjjBATSystematic * parentsyst = fSystematicContainer.at(GetSystematicIndex(thisshapevar->GetParentSystematic()));

     if (!(parentsyst->GetFlagSystematicActive()))
        continue;

     // get parameter index
     int parindex = parentsyst->GetParamIndex();

//     std::cout << "Got scale changing systematic " << j << std::endl;

     for (int bin = 0; bin<nbins; bin++) {
       double shift = thisshapevar->GetBinAdjustment(bin, parameters.at(parindex));
       double norm = GetProcess(processindex)->GetBinValue(bin,parameters);
       fProcessExpectation[bin] = (fProcessExpectation[bin] + norm*shift);
//       std::cout << "In bin " << bin << " will return " << fProcessExpectation[bin] << " + " << norm << " * " << shift << " = " << fProcessExpectation[bin] << std::endl;
     }
   }

   // get shape-changing systematic variations for process
   vector<MjjBATShapeChangingSyst*> shapevars = GetProcess(processindex)->GetShapeChangingSysts();

   // get template-based systematic variations, if one exists
   MjjBATTemplateSyst * templateVar = GetProcess(processindex)->GetTemplateSyst();

   // do template-changing systematic if one exists
   if (templateVar!=0 && shapevars.size()==0) {

//       std::cout << "Doing the one template var" << std::endl;

       MjjBATSystematic * parentsyst = fSystematicContainer.at(GetSystematicIndex(templateVar->GetParentSystematic()));
       if ((parentsyst->GetFlagSystematicActive())) {

         // get parameter index
         int parindex = parentsyst->GetParamIndex();

         for (int bin = 0; bin<nbins; bin++) {
           double shift = templateVar->GetBinAdjustment(bin, parameters.at(parindex));
           double norm = GetProcess(processindex)->GetBinValue(bin,parameters);
           fProcessExpectation[bin] = (fProcessExpectation[bin] + norm*shift);
         }
      }
   }

   // loop over all shape-changing systematic variations
   for (unsigned int j=0; j<shapevars.size(); ++j) {

//      std::cout << "Doing the shape vars. On " << j << std::endl;

     MjjBATShapeChangingSyst * thisshapevar = shapevars.at(j);
     MjjBATSystematic * parentsyst = fSystematicContainer.at(GetSystematicIndex(thisshapevar->GetParentSystematic()));

     if (!(parentsyst->GetFlagSystematicActive()))
        continue;

     // get parameter index
     int parindex = parentsyst->GetParamIndex();

     fProcessExpectation*=thisshapevar->GetMatrix(parameters.at(parindex));

   }

   // get scale-changing systematic variations for process
   vector<MjjBATScaleChangingSyst*> scalevars = GetProcess(processindex)->GetScaleChangingSysts();

   // For finding central region
   double thisPercentage=0;
   double thisInterval=0;
   int lowbin=1;
   int highbin=1;
   double smallestInterval= fData->GetBinLowEdge(fLastBin) + fData->GetBinWidth(fLastBin) - fData->GetBinLowEdge(fFirstBin) + 1e12;
   bool go = true;
   if (GetProcess(processindex)->GetTrimProcess()) {
     for (int bin1=0; bin1<=fProcessExpectation.GetNoElements(); bin1++) { // Lydia, changed fFirstBin to 0 to solve error band issue
       for (int bin2=fProcessExpectation.GetNoElements()-1; bin2>=bin1; bin2--) {
         thisPercentage = fProcessExpectation.GetSub(bin1,bin2).Sum()/fProcessExpectation.Sum();
         if (thisPercentage < GetProcess(processindex)->GetTrimPercentage()) {
           if (bin2==fProcessExpectation.GetNoElements()-1) go=false;
           break;
         }
         thisInterval = fData->GetBinLowEdge(bin2) + fData->GetBinWidth(bin2) - fData->GetBinLowEdge(bin1);
         if (thisInterval < smallestInterval) {
           lowbin=bin1;
           highbin=bin2;
           smallestInterval=thisInterval;
         }
       } // end of inner loop
       if (!go) break;
     } // end of outer loop
   
   } // end of if statement

   // loop over all bins and get adjustments
   for (int bin = 0; bin<nbins; bin++) {

     if (GetProcess(processindex)->GetTrimProcess() && (bin < lowbin || bin > highbin)) {
       expectation.push_back(std::make_pair(0.0,0.0));
     } else {

       double adjusted = fProcessExpectation[bin];

       // loop over all scale-changing systematic variations
       for (unsigned int j=0; j<scalevars.size(); ++j) {
         MjjBATScaleChangingSyst * thisscalevar = scalevars.at(j);
         MjjBATSystematic * parentsyst = fSystematicContainer.at(GetSystematicIndex(thisscalevar->GetParentSystematic()));

         if (!(parentsyst->GetFlagSystematicActive()))
            continue;

         int parindex = parentsyst->GetParamIndex();
         adjusted = scalevars.at(j)->GetAdjustedBin(adjusted,parameters.at(parindex));

       }
       
       // If this is a process with a correlated uncertainty on the normalization, such as the fit function:
       // the variance will be the size of the uncertainty recorded when defining the process.
       // If this process has a statistical uncertainty bin by bin,
       // the variance comes from the number of MC events that were used, I suppose.
       // Does an uncertainty on the normalization of the sample come into play here as well?
       // For now, assume no. That will simply be what the process parameter corresponds to,
       // but will not change its variance except by scaling.
       // Return sqrt(variance) not variance because we want it in this form most often anyway.

       double var;
        if (GetProcess(processindex)->GetDoNormalisationUnc())
          var = GetProcess(processindex)->GetVariationHistogram()->GetBinContent(bin);
        else
          var = 0.0;
//       if (GetProcess(processindex)->GetDoStatUnc()) {
//          // Variance will be sqrt(unweighted number of events) *
//          // current normalization / unweighted normalization.
//          // What is in histogram we put in: in each bin, true n * some weight and sqrt(true n) * some weight
//          // What I have now: some other normalization H. So:
//          // E = e*w ; B = e*e*w
//          // e = E/w = sqrt(B/w)
//          // Desired: e * H/(w*e*e) = H/(w*e) = H/E = adjusted/error I grab off of the nominal histogram?
//          var = adjusted / nominalHist->GetBinError(bin);
//       } else if (GetProcess(processindex)->GetDoNormalisationUnc() && !(GetProcess(processindex)->GetDoStatUnc())){
//          var = GetProcess(processindex)->GetVariationHistogram()->GetBinContent(bin);
//       }
       expectation.push_back(std::make_pair(adjusted,var));
     }
   }

   // If we have anything out farther than we have decided to compare so far, extend range now
   if (highbin > fLastBin) fLastBin = highbin;

   return expectation;

}

// ---------------------------------------------------------
int MjjBATModel::PrintStack(const std::vector<double> & parameters, std::string filename, std::string options)
{
	// todo:
	// - remove x-error on data points
	// - use hatched fill for error band

   // check if parameters are filled
   if (!parameters.size())
      return -1;

   // check options
   bool flag_logx   = false; // plot x-axis in log-scale
   bool flag_logy   = false; // plot y-axis in log-scale
   bool flag_bw     = false; // plot in black and white

   bool flag_sum    = false; // plot sum of all templates
   bool flag_stack  = false; // plot stack of templates

   bool flag_e0     = false; // do not draw error bars on data
   bool flag_e1     = false; // draw sqrt(N) error bars on data

   bool flag_b0     = false; // draw an error band on the expectation
   bool flag_b1     = false; // draw an error band on the number of events

   if (std::string(options).find("logx") < std::string(options).size())
      flag_logx = true;

   if (std::string(options).find("logy") < std::string(options).size())
      flag_logy = true;

   if (std::string(options).find("bw") < std::string(options).size())
      flag_bw = true;

   if (std::string(options).find("sum") < std::string(options).size())
      flag_sum = true;

   if (std::string(options).find("stack") < std::string(options).size())
      flag_stack = true;

   if (std::string(options).find("e0") < std::string(options).size())
      flag_e0 = true;

   if (std::string(options).find("e1") < std::string(options).size())
      flag_e1 = true;

   if (std::string(options).find("b0") < std::string(options).size())
      flag_b0 = true;

   if (std::string(options).find("b1") < std::string(options).size())
      flag_b1 = true;

	 if (!flag_e0)
		 flag_e1=true;

   // create canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // set log or linear scale
   if (flag_logx)
      c1->SetLogx();

   if (flag_logy)
      c1->SetLogy();
	 
   // get number of bins
   int nbins = fData->GetNbinsX();

   // define sum of templates
   TH1D* hist_sum = new TH1D(*fData);
   hist_sum->SetLineColor(kBlack);
   for (int i = 1; i <= nbins; ++i)
      hist_sum->SetBinContent(i, 0);

   // define error band
   TH1D* hist_error_band = new TH1D(*fData);
   hist_error_band->SetFillColor(kBlack);
   hist_error_band->SetFillStyle(3005);
   hist_error_band->SetLineWidth(1);
   hist_error_band->SetStats(kFALSE);
   hist_error_band->SetMarkerSize(0);

   TGraphAsymmErrors * graph_error_exp = new TGraphAsymmErrors(nbins);
   graph_error_exp->SetMarkerStyle(0);
   graph_error_exp->SetFillColor(kBlack);
   graph_error_exp->SetFillStyle(3005);

   // fill error band
   if (flag_b0) {
      for (int i = 1; i <= nbins; ++i) {
	 TH1D * proj = fHistUncertaintyBandExpectation->ProjectionY("_py", i, i);
	 if (proj->Integral() > 0)
		 proj->Scale(1.0 / proj->Integral());
	 double quantiles[3];
	 double sums[3] = {0.16, 0.5, 0.84};
	 proj->GetQuantiles(3, quantiles, sums);
	 graph_error_exp->SetPoint(i-1, fData->GetBinCenter(i), quantiles[1]);
	 graph_error_exp->SetPointError(i-1, 0.0, 0.0, quantiles[1] - quantiles[0], quantiles[2]-quantiles[1]);
	 hist_error_band->SetBinContent(i, 0.5*(quantiles[2]+quantiles[0]));
	 hist_error_band->SetBinError(i, 0, 0.5*(quantiles[2]-quantiles[0]));
	 delete proj;
      }
   }
	 
   // create stack
   THStack * stack = new THStack("", "");

   // create a container of temporary histograms
   std::vector<TH1D *> histcontainer;

   // get number of templates
   unsigned int ntemplates = GetNProcesses();

   // loop over all templates
   for (unsigned int i = 0; i < ntemplates; ++i) {

      // get histogram
      TH1D * temphist = GetProcess(i)->GetNominalHistogram();

      // create new histogram
      TH1D * hist(0);

      if (temphist)
         hist = new TH1D( *(temphist) );
      else
         continue;

      if (flag_bw) {
         hist->SetFillColor(0);
         hist->SetLineWidth(1);
         hist->SetLineStyle(int(1+i));
      }
      else {
         hist->SetFillColor(2 + i);
         hist->SetFillStyle(1001);
      }

      vector<double> expectation = Expectation(parameters);

      // scale histogram
      for (int ibin = 1; ibin <= nbins; ++ibin) {

         // set bin content
         hist->SetBinContent(ibin, expectation.at(ibin));

	 // add bin content
	 hist_sum->SetBinContent(ibin, hist_sum->GetBinContent(ibin) + expectation.at(ibin));
      }

      // add histogram to container (for memory management)
      histcontainer.push_back(hist);

      // add histogram to stack
      stack->Add(hist);
   }

   //draw data
   fData->Draw("P0");

   // set range user
   fData->GetYaxis()->SetRangeUser(fRangeYMin, fRangeYMax);

   // draw stack
   if (flag_stack)
         stack->Draw("SAMEHIST");

   // draw error band on number of observed events
   if (flag_b1) {
      CalculateUncertaintyBandPoisson(0.001, 0.999, kRed)->Draw("SAMEE2");
      CalculateUncertaintyBandPoisson(0.023, 0.977, kOrange)->Draw("SAMEE2");
      CalculateUncertaintyBandPoisson(0.159, 0.841, kGreen)->Draw("SAMEE2");
   }

   // draw error band on expectation
   if (flag_b0) {
      hist_error_band->Draw("SAMEE2");
   }

   if (flag_sum)
      	hist_sum->Draw("SAME");

   //draw data again
   if (flag_e0)
	 fData->Draw("SAMEP0");

   if (flag_e1)
	 fData->Draw("SAMEP0E");

   // redraw the axes
   gPad->RedrawAxis();

   // print
   c1->Print(filename.c_str());

   // free memory
   for (unsigned int i = 0; i < histcontainer.size(); ++i) {
      TH1D * hist = histcontainer.at(i);
      delete hist;
   }
   delete stack;
   delete c1;
	 delete graph_error_exp;
	 delete hist_error_band;
	 delete hist_sum;

   // no error
   return 1;
}

// ---------------------------------------------------------
double MjjBATModel::CalculateChi2(const std::vector<double> & parameters)
{
   if (parameters.size() == 0)
      return -1;

   double chi2 = 0;

   // check if histogram exists
   if (fData) {
      // get number of bins in data
      int nbins = fData->GetNbinsX();

      vector<double> expectation = Expectation(parameters);

      // loop over all bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {

         // get observation
         double observation = fEffectiveData->GetBinContent(ibin);

         // add Poisson term
         chi2 += (expectation.at(ibin) - observation) * (expectation.at(ibin) - observation) / expectation.at(ibin);
      }
   }

   return chi2;
}

// ---------------------------------------------------------
double MjjBATModel::CalculateCash(const std::vector<double> & parameters)
{
   if (parameters.size() == 0)
      return -1;

   double cash = 0;

   // check if histogram exists
   if (fData) {
      // get number of bins in data
      int nbins = fData->GetNbinsX();

      vector<double> expectation = Expectation(parameters);

      // loop over all bins
      for (int ibin = 1; ibin <= nbins; ++ibin) {

         // get observation
         double observation = fEffectiveData->GetBinContent(ibin);

         // calculate Cash statistic
         cash += 2. * (expectation.at(ibin) - observation);

         // check negative log
         if (observation > 0)
            cash += 2. * observation * log (observation/expectation.at(ibin));
      }
   }

   // return cash;
   return cash;

}

// ---------------------------------------------------------
double MjjBATModel::LogLikelihood(const std::vector<double> & parameters)
{

   if (fLikelihood == BUMPHUNTER) return BHLogLikelihood(parameters);

   double logprob = 0.;

   if (!fData)
      return 0;

   vector<double> expectation = Expectation(parameters);

   // loop over all bins
   for (int ibin = fFirstBin; ibin < fLastBin+1; ++ibin) {

      // get observation
      double observation = fEffectiveData->GetBinContent(ibin);

      // add Poisson term
      logprob += SimpleLogPoisson(observation,expectation.at(ibin));
//      std::cout << "In bin " << ibin << " adding L(" << observation << ","<< expectation.at(ibin) << " = " << SimpleLogPoisson(observation,expectation.at(ibin)) << std::endl;

   }

//   std::cout << "With sig param " << parameters.at(1) << " log l is " << logprob << std::endl;

   return logprob;
}

// ---------------------------------------------------------
double MjjBATModel::BHLogLikelihood(const std::vector<double> & parameters)
{

  std::cout << "YOU SHOULD NEVER BE SEEING THIS!"  << std::endl;
  return 0.0;

  // Need to implement the model for the BumpHunter now.
  // This implementation is only valid for a single background contribution with
  // an assumed Gamma prior. Check that before calling this.
  if (fProcessContainer.size() > 1) {
    std::cout << "This likelihood definition only works for a single background process!" << std::endl;
    return -1;
  }
  vector<std::pair<double,double> > processExp = ProcessExpectation(0, parameters);
  // For the BumpHunter we want the sum of the relevant bins to be treated as a single bin.
  // This is not the same as fFirstBin, fLastBin, which we don't want to change due to
  // using them for limiting the overall range.
  // Add new "interesting bins" and let them be set between minimizations from some external location.
  
  double total=0, totalErr = 0;
  for (int bin = fFirstInterestingBin; bin < fLastInterestingBin+1; bin++) {
    total += processExp.at(bin).first/fWeights->GetBinContent(bin);
    totalErr += processExp.at(bin).second/fWeights->GetBinContent(bin);
  }
//  double n = std::max(parameters.at(GetParameter("dummyN")->GetIndex()),0.0);
  double n = std::max(parameters.at(GetParameters().Index("dummyN")),0.0);
  n = std::floor(n);
  double beta = (total)/(totalErr*totalErr);
  double alpha = total*beta;

  if (!(GetProcess(0)->GetDoNormalisationUnc()))
      return SimpleLogPoisson(n,total);
  
  // This is the expression:
  //double myval = pow(beta,alpha)/[n! * pow((1+beta),(n+alpha))] * TMath::Gamma(n+alpha)/TMath::Gamma(alpha);
  // This is the log of the expression:
  // ln((a*b)/(c*d)) = ln(a) + ln(b) - ln(c) - ln(d)
  // and ln(x^y) = y ln x
  // Note also that ln(n!) = n*ln(n) - n + 1 (Stirling's approximation).
  // So log(myval) is:


  if (fDebug) std::cout << "total, totalErr are " << total << ", " << totalErr << std::endl;
  if (fDebug) std::cout << "alpha = " << alpha << ", beta = " << beta << ", alpha+n = " << alpha+n << ", n = " << n << std::endl;

  double myval;
  if ((alpha) > 120.0 && n > 120.0) {  if (fDebug) std::cout << "one" << std::endl;
    myval = alpha*log(beta) - (n+alpha)*log(1.0+beta) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha) - TMath::LnGamma(n+1);
  } else if ((n+alpha-1) > 120.0 && n < 120.0 && alpha < 120.0) { if (fDebug) std::cout << "two" << std::endl;
    myval = log(pow(beta,alpha)/(pow((1.0+beta),(n+alpha)) * TMath::Gamma(alpha) * TMath::Factorial(n))) + TMath::LnGamma(n+alpha);
 } else if (alpha+n > 120.0 && n <= 120.0) { if (fDebug) std::cout << "three" << std::endl;
    myval = log(pow(beta,alpha)/(pow(1.0+beta,n+alpha) * TMath::Factorial(n))) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha);
  } else if (alpha+n > 120.0 && alpha <= 120.0) { if (fDebug) std::cout << "four" << std::endl;
    myval = log(pow(beta,alpha)/(pow((1.0+beta),(n+alpha)) * TMath::Gamma(alpha))) + TMath::LnGamma(n+alpha) - TMath::LnGamma(n+1);
  } else { if (fDebug) std::cout << "five" << std::endl;
    myval = log(pow(beta,alpha) * TMath::Gamma(n+alpha)/(TMath::Factorial(n) * pow((1.0+beta),(n+alpha)) * TMath::Gamma(alpha)));
  }
  if (!(std::isfinite(myval))) { if (fDebug) std::cout << "replacing " << myval << " with "; myval = alpha*log(beta) - (n+alpha)*log(1.0+beta) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha) - TMath::LnGamma(n+1); if (fDebug) std::cout << myval << std::endl;}

  if (fDebug) std::cout << "Returning " << myval << std::endl;
  return myval;

  // Really detailed version -- accurate but very slow
  /*
  // Simplest case (small numbers): we can do everything as a log.
  if (AreIdentical(n,0.0)) return alpha * log(beta/(1.0+beta));
  double holdval = pow(beta,alpha) * TMath::Gamma(n+alpha)/(TMath::Factorial(n) * pow((1.0+beta),(n+alpha)) * TMath::Gamma(alpha));
  if (std::isfinite(log(holdval))) return log(holdval);
  // Next largest term will be Gamma(n+alpha). Try removing.
  holdval = pow(beta,alpha)/(TMath::Factorial(n) * pow((1.0+beta),(n+alpha)) * TMath::Gamma(alpha));
  if (std::isfinite(log(holdval))) return log(holdval) + TMath::LnGamma(n+alpha);

  // Still a problem: remove either Gamma(alpha) or Factorial(n) depending on which is bigger.
  if (n > alpha - 1) {
    // Remove n term.
    holdval = pow(beta,alpha)/(pow((1.0+beta),(n+alpha)) * TMath::Gamma(alpha));
    if (std::isfinite(log(holdval)))
         return log(holdval) + TMath::LnGamma(n+alpha) - TMath::LnGamma(n+1); 
    // Alpha is a problem too: take it out.
    holdval = pow(beta,alpha)/pow((1.0+beta),(n+alpha));
    if (std::isfinite(log(holdval))) 
         return log(holdval) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha) - TMath::LnGamma(n+1);
    // Still broken? Log it all.
    return alpha*log(beta) - (n+alpha)*log(1.0+beta) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha) - TMath::LnGamma(n+1);

  // Then it was the alpha that was a problem first.
  } else {
    // Remove alpha term.
    holdval = pow(beta,alpha)/(TMath::Factorial(n) * pow((1.0+beta),(n+alpha)));
    if (std::isfinite(log(holdval))) return log(holdval) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha);
    // Still not good? Remove n term too.
    holdval = pow(beta,alpha)/pow((1.0+beta),(n+alpha));
    if (std::isfinite(log(holdval))) 
         return log(holdval) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha) - TMath::LnGamma(n+1); 
  }
  // Still broken? Log it all.
  return n > 5 ? alpha*log(beta) - (n+alpha)*log(1.0+beta) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha) - n*log(n) + n - 1.0 :
                 alpha*log(beta) - (n+alpha)*log(1.0+beta) + TMath::LnGamma(n+alpha) - TMath::LnGamma(alpha) - log(TMath::Factorial(n));
  */

}

// ---------------------------------------------------------
double MjjBATModel::SimpleLogPoisson(double x, double par) {

  int ix = Int_t(x);

  if (ix<0)
    // Poisson does return 0;
    return -1*std::numeric_limits<double>::max();
  else if (ix == 0.0)
    // Poisson does return 1./Exp(par);
    return -1*par;
  else {
    Double_t lnpoisson = ix*log(par)-par-TMath::LnGamma(ix+1.);
    // Poisson does return Exp(lnpoisson);
    return lnpoisson;
  }

}

// ---------------------------------------------------------
void MjjBATModel::MCMCUserIterationInterface()
{

   // check if histogram exists
   if (!fData)
      return;

   // check if histogram exists
   if (!fHistUncertaintyBandExpectation)
      return;

   // get number of bins in data
   int nbins = fData->GetNbinsX();

   vector<double> expectation = Expectation(fMCMCx);

   // loop over all bins
   for (int ibin = 1; ibin <= nbins; ++ibin) {

      // fill uncertainty band on expectation
      fHistUncertaintyBandExpectation->Fill(fData->GetBinCenter(ibin), expectation.at(ibin));
   }

}

// ---------------------------------------------------------
void MjjBATModel::PrintTemplates(std::string filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->Divide(2, 2);

   std::string f(filename);

   // get number of templates
   unsigned int ntemplates = fProcessContainer.size();

   // calculate number of pages
   unsigned int npages =  ntemplates / 4;
   if (ntemplates % 4 > 0)
      npages++;

   // loop over pages
   for (unsigned int i = 0; i < npages; ++i) {
      // loop over pads
      for (unsigned int j = 0; j < 4; ++j) {
         // calculate template index
         unsigned int templateindex = 4 * i + j;

         if (templateindex < ntemplates && GetProcess(templateindex)->GetNominalHistogram()) {
            // cd into pad
            c1->cd(j+1);

            // get histogram
            TH1D * hist = new TH1D( *( (TH1D *) GetProcess(templateindex)->GetNominalHistogram()->Clone() ) );

            // draw histogram
            hist->Draw("HIST");
         }
         else {
            // clear pad
            c1->cd(j+1)->Clear();
         }
      }

      if (npages == 1)
         c1->Print(f.c_str());
      else if (i == 0) {
         c1->Print( (f+std::string("[")).c_str() );
         c1->Print(f.c_str());
      }
      else if (i < npages-1) {
         c1->Print(f.c_str());
      }
      else {
         c1->Print(f.c_str());
         c1->Print( (f+std::string("]")).c_str() );
      }
   }

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void MjjBATModel::PrintHistUncertaintyBandExpectation(std::string filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // draw histogram
   fHistUncertaintyBandExpectation->Draw("COLZ");
   c1->Draw();

   // print
   c1->Print(filename.c_str());

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void MjjBATModel::CalculateHistUncertaintyBandPoisson()
{
   // calculate histogram
   int nbinsy_exp = fHistUncertaintyBandExpectation->GetNbinsY();
   int nbinsx_poisson = fHistUncertaintyBandPoisson->GetNbinsX();
   int nbinsy_poisson = fHistUncertaintyBandPoisson->GetNbinsY();

   // loop over x-axis of observation
   for (int ix = 1; ix <= nbinsx_poisson; ++ix) {
      double sum_w = 0;
      // loop over y-axis of expectation and calculate sum of weights
      for (int iy = 1; iy <= nbinsy_exp; ++iy) {
         double w = fHistUncertaintyBandExpectation->GetBinContent(ix, iy);
         sum_w += w;
      }
      // loop over y-axis of expectation
      for (int iy = 1; iy <= nbinsy_exp; ++iy) {
         double w = fHistUncertaintyBandExpectation->GetBinContent(ix, iy)/sum_w;
         double expectation = fHistUncertaintyBandExpectation->GetYaxis()->GetBinCenter(iy);
         // loop over y-axis of observation
         for (int jbin = 1; jbin <= nbinsy_poisson; ++jbin) {
            double p = TMath::PoissonI(double(jbin-1), expectation);
            double bincontent = 0;
            if (iy>1)
               bincontent=fHistUncertaintyBandPoisson->GetBinContent(ix, jbin);
            fHistUncertaintyBandPoisson->SetBinContent(ix, jbin, bincontent+p*w);
         }
      }
   }
}

// ---------------------------------------------------------
TH1D* MjjBATModel::CalculateUncertaintyBandPoisson(double minimum, double maximum, int color)
{
   TH1D* hist = new TH1D(*fData);
   hist->SetMarkerSize(0);
   hist->SetFillColor(color);
   hist->SetFillStyle(1001);

   int nbinsx_poisson = fHistUncertaintyBandPoisson->GetNbinsX();
   int nbinsy_poisson = fHistUncertaintyBandPoisson->GetNbinsY();

   // loop over x-axis of observation
   for (int ix = 1; ix <= nbinsx_poisson; ++ix) {
      double sum_p = 0;  // sum of all probabilities inside the interval
      int limit_min = 0;
      int limit_max = nbinsx_poisson-1;

      // loop over y-axis of observation
      for (int jbin = 1; jbin <= nbinsy_poisson; ++jbin) {
         double p = fHistUncertaintyBandPoisson->GetBinContent(ix, jbin);
         sum_p+=p;
         if (sum_p < minimum)
            limit_min=jbin;
         if (sum_p > maximum && (sum_p - p) < maximum )
            limit_max=jbin-1;
      }
      hist->SetBinContent(ix, 0.5*double(limit_min+limit_max));
      hist->SetBinError(ix, 0.5*double(limit_max-limit_min));
   }

   return hist;
}

// ---------------------------------------------------------
void MjjBATModel::PrintHistCumulativeUncertaintyBandPoisson(std::string filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // calculate error band
   CalculateHistUncertaintyBandPoisson();

   TH2D hist(*fHistUncertaintyBandPoisson);

   int nbinsx_poisson = hist.GetNbinsX();
   int nbinsy_poisson = hist.GetNbinsY();

   // loop over x-axis of observation
   for (int ix = 1; ix <= nbinsx_poisson; ++ix) {
      double sum_p = 0;  // sum of all probabilities inside the interval

      // loop over y-axis of observation
      for (int jbin = 1; jbin <= nbinsy_poisson; ++jbin) {
         double p = hist.GetBinContent(ix, jbin);
         sum_p+=p;
         hist.SetBinContent(ix, jbin, sum_p);
      }
   }

   // draw histogram
   hist.Draw("COLZ");
   c1->Draw();

   // print
   c1->Print(filename.c_str());

   // free memory
   delete c1;
}

// ---------------------------------------------------------
void MjjBATModel::PrintHistUncertaintyBandPoisson(std::string filename)
{
   // create new canvas
   TCanvas * c1 = new TCanvas();
   c1->cd();

   // calculate error band
   CalculateHistUncertaintyBandPoisson();

   // draw histogram
   fHistUncertaintyBandPoisson->Draw("COLZ");
   c1->Draw();

   // print
   c1->Print(filename.c_str());

   // free memory
   delete c1;
}

// ---------------------------------------------------------
vector<double> MjjBATModel::GetProcessNorms(const std::vector<double> params) {

   vector<double> norms; norms.clear();

   // loop over all processes
   for (int i = 0; i < fNProcesses; ++i) {

      // get normalisation of process
      double nEvents = GetProcess(i)->GetNEvents(params);
      norms.push_back(nEvents);

   }

   return norms;

}

// ---------------------------------------------------------

