// ---------------------------------------------------------

#include "Bayesian/MjjFitter.h"

// ---------------------------------------------------------
MjjFitter::MjjFitter() :
usedFitWindows("null","null",1,0.0,1.0,1,0.0,1.0),residualSwiftBins("null","null",1,0.0,1.0,1,0.0,1.0)
{
    
    // create minimizer and function based on LikeMin
    // Minimizer uses specified minimization type
    // kCombined is combined Migrad + Simplex minimization
    // From mnusersguide.pdf: Migrad is generally supposed to be best. Thus combined should be okay.
    fMinuitMinimizer = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kCombined);
    
    // Use non-chatty default print level.
    fPrintLevel = -1; // WAS -1
    fMinuitMinimizer->SetPrintLevel(fPrintLevel);
    
    // Set various options for minimizer
    //ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("AdaptiveSingular");
    //ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("GaussLegendre");
    
    // Accuracy options
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);
    fMinuitMinimizer->SetMaxFunctionCalls(500000);
    fMinuitMinimizer->SetMaxIterations(100000);
    fMinuitMinimizer->SetTolerance(0.001); // changed from 0.0001
    
    // Indicates fitter can take longer to ensure best fit
    fMinuitMinimizer->SetStrategy(2);
    
    // See user manual. Value depends on type of minimization used.
    // Chi^2, negative log likelihood have diff values.
    // This is for -log(Likelihood). Sets FCN::up in Minuit.
    fMinuitMinimizer->SetErrorDef(0.5);
    
    // Parameter step size
    fStepSize = 0.0001;
    
    // Adjust printout levels for errors
    // to suppress massages like Info in <Minuit2>:[...]
    gErrorIgnoreLevel=kWarning; //1001;
    
    // By default, do not use signal.
    fFitWithSignal = false;
    fSignalParameterIndex = -1;
    fFixNormalisation = -2.0;
    fDoExtraPreliminaryFits = false;
    
    // Initialize retry value
    fRetry = 0;
    
    // Set up SWIFT defaults
    fDoSwift = false;
    fNBinsLeft = 13;
    fNBinsRight = 13;
    fLowXEstimate = -1;
    fHighXEstimate = -1;
    fFixWidthAtLowEnd = true;
    fFixWidthAtHighEnd = false;
    fTruncateHigh = false;
    fMinXForFit = std::numeric_limits<double>::quiet_NaN();
    fMaxXForFit = std::numeric_limits<double>::quiet_NaN();
    
    fSignalTemplate = NULL;
    
    latestFitPars.clear();
    
}

// ---------------------------------------------------------
MjjFitter::~MjjFitter()
{
    
}

// ---------------------------------------------------------
double MjjFitter::LikeMin(const double *params)
{
    
    for (int i=0; i<fFitFunction->GetNpar(); i++) {
        
        fFitFunction->SetParameter(i,params[i]);
    }
    
    double nSignal = 0;
    if (fFitWithSignal) nSignal = params[fSignalParameterIndex];
    
    double logL = 0;
    double data, weight, fitEquivalentOfNEvents;
    //std::cout << "Using bins " << "[" << fStartBin << ", " << fStopBin << "]"<< std::endl;
    for (int bin = fStartBin; bin < fStopBin+1; bin++) {
        
        data = fEffectiveHistogramToFit.GetBinContent(bin);
        weight = fWeightsHistogram.GetBinContent(bin);
        fitEquivalentOfNEvents = fFitFunction->Integral(fBasicHistogramToFit.GetBinLowEdge(bin),fBasicHistogramToFit.GetBinLowEdge(bin+1));
        
        if (fMjjFitFunction.GetDoWindowExclusion()) {
            std::pair <double,double> boundaries = fMjjFitFunction.GetWindowBoundaries();
            
            // Exclude all bins that overlap with the bump
            // if top edge of bin touches any region of exclusion: continue
            if ((fBasicHistogramToFit.GetBinLowEdge(bin+1) >= boundaries.first &&
                 fBasicHistogramToFit.GetBinLowEdge(bin+1) <= boundaries.second) ||
                (fBasicHistogramToFit.GetBinLowEdge(bin) >= boundaries.first &&
                 fBasicHistogramToFit.GetBinLowEdge(bin) <= boundaries.second) ) continue;
        }
        
        double sigEquivalentOfNEvents = 0;
        if (fFitWithSignal) sigEquivalentOfNEvents = fSignalTemplate->GetBinContent(bin)*nSignal;
        
        if (fitEquivalentOfNEvents<0) fitEquivalentOfNEvents = 0;
        
        logL += -SimpleLogPoisson( data , (fitEquivalentOfNEvents + sigEquivalentOfNEvents)/weight);
        
        //    std::cout << "In bin " << bin << " adding L(" << data << ", " << (fitEquivalentOfNEvents +  sigEquivalentOfNEvents/weight) << ") = " << -SimpleLogPoisson(data , (fitEquivalentOfNEvents + sigEquivalentOfNEvents)/weight)<< std::endl;
        
    }
    
    //  std::cout << "With sig parameter " << params[4] << " log l is " << logL << std::endl;
    
    return logL;
    
}

// ---------------------------------------------------------
double MjjFitter::SimpleLogPoisson(double x, double par)
{
    
    // Replacement for call to log(TMath::PoissonI(x,par))
    
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
void MjjFitter::SetSignalTemplate(TH1D * sigTemplate, double sigMass)
{
    
    // Create new TH1D with unique name
    fSignalTemplate = (TH1D*) sigTemplate->Clone();
    std::string newName = "sigtemp_";
    newName.append(sigTemplate->GetName());
    fSignalTemplate->SetName(newName.c_str());
    
    // Normalize it so that amount by which it is scaled
    // will equal number of events
    double integral = fSignalTemplate->Integral();
    for (int i=0; i<fSignalTemplate->GetNbinsX()+2; i++) {
        double oldContent = fSignalTemplate->GetBinContent(i);
        fSignalTemplate->SetBinContent(i, oldContent/integral);
    }
    
    // Switch fitting default to use signal.
    fFitWithSignal = true;
    fSignalMass = sigMass;
    
    return;
    
}

// ---------------------------------------------------------
bool MjjFitter::MinuitMLFit ()
{
    
    // Erase parameters from previous fit, if any
    fMinuitMinimizer->Clear();
    
    // Number of params
    int numberOfParameters = (int) fFitFunction->GetNpar();
    
    // initializes all function parameters and sets limits if applicable
    vector <double> initialParamVals;
    for (int iLoop=0;iLoop<numberOfParameters;iLoop++){
        
        initialParamVals.push_back(fFitFunction->GetParameter(iLoop));
        
        
        bool isParamFixed = fMjjFitFunction.GetParameter(iLoop)->GetIsFixed();
        std::pair<bool,bool> hasLimits = fMjjFitFunction.GetParameter(iLoop)->GetHasParamLimits();
        double minLimit = fMjjFitFunction.GetParameter(iLoop)->GetParamLimitLow();
        double maxLimit = fMjjFitFunction.GetParameter(iLoop)->GetParamLimitHigh();
        
        // Is this a retry? If so, deviate slightly from initial start value.
        double adjustment = (fRetry%2)*(-0.1*fRetry)*initialParamVals[iLoop];
        
        if (isParamFixed) {
            fMinuitMinimizer->SetFixedVariable(iLoop,Form("a%d",iLoop),initialParamVals[iLoop]);
        } else if (hasLimits.first==true && hasLimits.second==true) {
            fMinuitMinimizer->SetLimitedVariable(iLoop,Form("a%d",iLoop),
                                                 initialParamVals[iLoop]+adjustment, fStepSize, minLimit, maxLimit);
        } else if (hasLimits.first==true && hasLimits.second==false) {
            fMinuitMinimizer->SetLowerLimitedVariable(iLoop,Form("a%d",iLoop),
                                                      initialParamVals[iLoop]+adjustment, fStepSize, minLimit);
        } else if (hasLimits.first==false && hasLimits.second==true) {
            fMinuitMinimizer->SetUpperLimitedVariable(iLoop,Form("a%d",iLoop),
                                                      initialParamVals[iLoop]+adjustment, fStepSize, maxLimit);
        } else {
            fMinuitMinimizer->SetVariable(iLoop,Form("a%d",iLoop),
                                          initialParamVals[iLoop]+adjustment, fStepSize);
        }
    }
    
    if (fFitWithSignal) {
        int sigParam = numberOfParameters;
        if (fDoSwift && fFixNormalisation > -1.0)
            fMinuitMinimizer->SetFixedVariable(sigParam,Form("a%d",sigParam),fFixNormalisation);
        else
            fMinuitMinimizer->SetLowerLimitedVariable(sigParam,Form("a%d",sigParam),0.0, 1,0.);
        initialParamVals.push_back(0);
        fSignalParameterIndex = sigParam;
        numberOfParameters++;
    }
    
    // This is just a wrapper for the function we defined in LikeMin, letting it be used by Minuit
    ROOT::Math::Functor functor = ROOT::Math::Functor(this, &MjjFitter::LikeMin, numberOfParameters);
    fMinuitMinimizer->SetFunction(functor);
    
    // Perform minimization.
    fMinuitMinimizer->Minimize();
    fFitFunction->SetParErrors(fMinuitMinimizer->Errors());
    
    
    // Did it work?
    bool didWork = false;
    if (fMinuitMinimizer->Status() < 2) didWork = true;
    
    //   ///DEBUG:
    if (fPrintLevel != -1) {
        cout << "Printing Minuit results: " << endl;
        fMinuitMinimizer->PrintResults();
        const double  *par_result = fMinuitMinimizer->X();
        std::cout << "Parameter results were: " << par_result[0] << " " << par_result[1] << " " << par_result[2] << " " << par_result[3] << " " << par_result[4] << std::endl;
        if (fPrintLevel != -1) cout <<  "-logL: " << LikeMin(par_result) << endl;
    }
    
    return didWork;
    
}

// ---------------------------------------------------------
// This helper function checks whether the function has limits on the parameters
bool isParamFixed (TF1* function, int i)
{
    double al, bl;
    function -> GetParLimits(i,al,bl);
    return (al*bl != 0 && al >= bl);
}


// ---------------------------------------------------------
// This helper function randomizes the parameters
void fudge(TF1 *function, TRandom3 random, int numTotParam = 0)
{
    
    //cout << "Parameters that are fudged: " << endl;
    //for (int i = 0; i < numTotParam; i++)
    //  cout << function -> GetParameter(i) << endl;
    
    bool isSB = 0;
    //    if (numTotParam > 8) isSB=1;
    
    double start = -999;
    double end = -999;
    
    for (int j = 0; j < numTotParam; j++)
    {
        // skip random twiddle if param is fixed
        if ( isParamFixed (function,j) ) continue;
        
        // uniform random numbers
        function -> GetParLimits(j, start, end);
        double rand = random.Uniform(start, end);
        function -> SetParameter( j, rand);
        
        //        if (isSB) function -> SetParameter(0, config.norm);
    }
    
    //cout << "Fudged parameters: " << endl;
    //for (int i = 0; i < numTotParam; i++)
    //  cout << function -> GetParameter(i) << endl;
    //cout << endl;
}

// ---------------------------------------------------------
bool MjjFitter::Fit(MjjHistogram & mjjHistogram, MjjFitFunction & mjjFitFunction, TH1D* blankhist)
{
    
    fMjjFitFunction = mjjFitFunction;
    
    //get fit function
    fFitFunction = mjjFitFunction.GetFitFunction();
    
    // Get histogram for preliminary fit
    TH1D normalizedHistogramToFit = (TH1D) mjjHistogram.GetNormalizedHistogram();
    
    // get histograms for main fit
    fBasicHistogramToFit = (TH1D) mjjHistogram.GetHistogram();
    fEffectiveHistogramToFit = (TH1D) mjjHistogram.GetEffectiveHistogram();
    fWeightsHistogram = (TH1D) mjjHistogram.GetWeightsHistogram();
    
    // Fit within range specified by function
    fMinXForFit = mjjFitFunction.GetMinXVal();
    fMaxXForFit = mjjFitFunction.GetMaxXVal();
    
    // Get start and stop bin values
    fStartBin = fBasicHistogramToFit.FindBin(fMinXForFit+1);
    if (fBasicHistogramToFit.FindBin(fMaxXForFit-1) > fBasicHistogramToFit.GetXaxis()->GetNbins())
        fStopBin = fBasicHistogramToFit.GetXaxis()->GetNbins();
    else fStopBin = fBasicHistogramToFit.FindBin(fMaxXForFit-1);
    
    //  std::cout << "Going to scan bin range " << fStartBin << " " << fStopBin << std::endl;
    
    mjjFitFunction.RestoreParameterDefaults();
    
    // Preliminary ROOT fit to get parameters near expected values
    for (int i=0; i<5; i++)  normalizedHistogramToFit.Fit(fFitFunction,"R0q","",fMinXForFit,fMaxXForFit);
    double * ROOTFitParams = fFitFunction->GetParameters();
    
    if (fDoExtraPreliminaryFits) {
        //do a few more preliminary ROOT fits to help convergence
        
        // Preliminary ROOT fit to get parameters near expected values
        int nGoodPreliminaryFits = 0;
        int nRequiredGoodPreliminaryFits = 5;
        //this will be a vector of parameters, one per function, index-parallel
        vector < vector <double> > allPrelimParamVals;
        int bestPreliminaryFitIndex = -999;
        double bestPreliminaryFitChi2 = 999;
        
        while (nGoodPreliminaryFits < nRequiredGoodPreliminaryFits) {
            std::cout << "Preliminary fits" << std::endl;
            //fit
            TFitResultPtr prelimFitPtr = normalizedHistogramToFit.Fit(fFitFunction,"R0Sq");
            //is the fit good? if not, randomize the input values and do it again
            if (!prelimFitPtr->IsValid()) {
                fudge(fFitFunction, fRnd, fFitFunction->GetNpar());
                continue;
            }
            else {
                nGoodPreliminaryFits++;
                //save parameters and chi2/NDF
                vector<double> thisPrelimParamVals;
                for (int iLoop=0;iLoop<fFitFunction->GetNpar();iLoop++){
                    thisPrelimParamVals.push_back(fFitFunction->GetParameter(iLoop));
                }
                allPrelimParamVals.push_back(thisPrelimParamVals);
                //check whether this fit is best
                double thisFitChi2 = prelimFitPtr->Chi2()/prelimFitPtr->Ndf();
                if (bestPreliminaryFitChi2 > thisFitChi2) {
                    std::cout << "This fit is best!" << std::endl;
                    std::cout << "Chi2 before/after" << bestPreliminaryFitChi2 << ", " << thisFitChi2 << std::endl;
                    bestPreliminaryFitChi2 = prelimFitPtr->Chi2()/prelimFitPtr->Ndf();
                    bestPreliminaryFitIndex = nGoodPreliminaryFits-1;//subtract 1 because we increase this counter at the start of the else
                }
            }
        }
    }
    
    // ML fit using Minuit minimizer
    fLatestFitStatus = MinuitMLFit();
    
    // If it didn't work, keep trying up to
    // 5 times until it does.
    fRetry = 0;
    while (!fLatestFitStatus && fRetry<5) {
        fRetry++;
        std::cout << "Retry number " << fRetry << std::endl;
        fFitFunction->SetParameters(ROOTFitParams);
        fLatestFitStatus = MinuitMLFit();
    }
    
    if (fLatestFitStatus == false) std::cout << "Fit failed to converge after five retries." << std::endl;
    
    fMjjFitFunction.MakeHistFromFunction(blankhist);
    return fLatestFitStatus;
    
}

void MjjFitter::setSwiftFitRange(const int binToEstimate, const int lowBinOfFitRange, const int highBinOfFitRange, const int highBinOfEstimate){                     
  // Getting low window edge. Default: either edge of range on left or current
  // position minus window width on left
  int lowBinToInclude = max(binToEstimate - fNBinsLeft, lowBinOfFitRange);
  // If we are truncating on the high edge and using the fixWidthAtHighEnd option,
  // this will actually be farther to the left than anticipated.
  if (fTruncateHigh && fFixWidthAtHighEnd) {
      if (highBinOfFitRange - binToEstimate < fNBinsRight) {
          lowBinToInclude = max(highBinOfFitRange - fNBinsLeft - fNBinsRight, lowBinOfFitRange);
          std::cout << "Truncating at high mjj and fixing width: reset lowBinToInclude to " << lowBinToInclude << std::endl;
      }
  }
  
  // Getting high window edge. Default: either edge of range on right or current
  // position plus window width on right
  int highBinToInclude = min(binToEstimate + fNBinsRight, highBinOfFitRange);
  // If we are truncating on the high edge, this is instead the min with
  // respect to highBinOfEstimate
  if (fTruncateHigh) highBinToInclude = min(highBinToInclude,highBinOfEstimate);
  // If we are near the low edge and using the fixWidthAtLowEnd option,
  // this will be farther to the right than anticipated.
  if (fFixWidthAtLowEnd) {
      if (binToEstimate - fNBinsLeft < lowBinOfFitRange)
          highBinToInclude = lowBinOfFitRange + fNBinsLeft + fNBinsRight;
  }
  
  // Now we know what bins to fit between: lowBinToInclude, highBinToInclude.
  // Set the fit to use these.
  fStartBin = lowBinToInclude;
  fStopBin = highBinToInclude;
}

// ---------------------------------------------------------
bool MjjFitter::SwiftFit(MjjHistogram & mjjHistogram, MjjFitFunction & mjjFitFunction, TH1D* blankhist)
{
    //fPrintLevel = 1;
    
    fMjjFitFunction = mjjFitFunction;
    
    //get fit function
    fFitFunction = mjjFitFunction.GetFitFunction();
    
    // Get histogram for preliminary fit
    TH1D normalizedHistogramToFit = (TH1D) mjjHistogram.GetNormalizedHistogram();
    
    // get histograms for main fit
    fBasicHistogramToFit = (TH1D) mjjHistogram.GetHistogram();
    fEffectiveHistogramToFit = (TH1D) mjjHistogram.GetEffectiveHistogram();
    fWeightsHistogram = (TH1D) mjjHistogram.GetWeightsHistogram();
    
    // Get start and stop bin values allowed for fitting as well as
    // intended for background estimate. These may have been
    // overwritten by SetSwiftParameters(): in that case, leave them.
    double minXForFit, maxXForFit;
    if (std::isnan(fMinXForFit)) minXForFit = mjjFitFunction.GetMinXVal();
    else minXForFit = fMinXForFit;
    if (std::isnan(fMaxXForFit)) maxXForFit = mjjFitFunction.GetMaxXVal();
    else maxXForFit = fMaxXForFit;
    
    // Get initial start and stop bin values
    fStartBin = fBasicHistogramToFit.FindBin(minXForFit+1);
    if ((fBasicHistogramToFit.FindBin(maxXForFit-1) > fBasicHistogramToFit.GetXaxis()->GetNbins())
        || maxXForFit == 0)
        fStopBin = fBasicHistogramToFit.GetXaxis()->GetNbins();
    else fStopBin = fBasicHistogramToFit.FindBin(maxXForFit-1);
    
    // If SWIFT wasn't turned on, we actually just want to do a global fit.
    // Make sure these parameters don't alter behaviour of the global fit.
    // When SWIFT is on, then we can have different full fit range than the
    // range in which we want an estimate. We only allow this on the high end though!
    if (!(fDoSwift)) fHighXEstimate = maxXForFit;
    
    if(fPrintLevel > 0) std::cout << "maxXforFit, fHighXEstimate: " << maxXForFit << " " << fHighXEstimate << std::endl;
    
    // Now set these.
    int lowBinOfEstimate = fStartBin;
    int lowBinOfFitRange = fStartBin;
    int highBinOfEstimate = fStopBin;
    int highBinOfFitRange;
    if (fBasicHistogramToFit.FindBin(fHighXEstimate-1) > fBasicHistogramToFit.GetXaxis()->GetNbins())
        highBinOfFitRange = fBasicHistogramToFit.GetXaxis()->GetNbins();
    else highBinOfFitRange = fBasicHistogramToFit.FindBin(fHighXEstimate-1);
    
    mjjFitFunction.RestoreParameterDefaults();
    
    //-----------------------------------------------
    // When SWIFT is turned on and a signal is used, step 1
    // is to do a fit centered where the signal sample is
    // We take a suitable estimate by looking at the peak bin.
    // This should work for sensible resonance shapes
    // It can be overridden in SetSignalTemplate if you have
    // something with a crazy low mass tail and so forth
    if (fDoSwift && fFitWithSignal) {
        
        // find mass to use unless overwritten
        double signalMass;
        int signalBinCenter = 0;
        if (fSignalMass > 0.0) {
            signalMass = fSignalMass;
            for (int b = 1; b < fSignalTemplate->GetNbinsX(); b++){
                double lowbinedge = fSignalTemplate->GetBinLowEdge(b);
                double hibinedge  = fSignalTemplate->GetBinLowEdge(b) + fSignalTemplate->GetBinWidth(b);
                //GetBinWidth
                cout << lowbinedge << " " << hibinedge<< endl;
                if (lowbinedge < signalMass && hibinedge> signalMass){
                    signalBinCenter = b;
                    break;
                }
            }
        } else {
            signalMass = fSignalTemplate->GetBinCenter(fSignalTemplate->GetMaximumBin());
            signalBinCenter = fSignalTemplate->GetMaximumBin();
        }
        if(fPrintLevel > 0) std::cout << "Using signal mass " << signalMass << " in preliminary fit from bin " << fStartBin << " to " << fStopBin << std::endl;

        ////////////// Determine bins to fit //////////////
        setSwiftFitRange(signalBinCenter, lowBinOfFitRange,highBinOfFitRange, highBinOfEstimate);
        
        //debug
        if (fPrintLevel > 0) {
            std::cout << "In signal fit set sliding window range: " << fStartBin << " " << fStopBin << std::endl;
            std::cout << "Which is " << fBasicHistogramToFit.GetBinLowEdge(fStartBin) << " " << fBasicHistogramToFit.GetBinLowEdge(fStopBin+1) << std::endl;
        }
        
        ////////////////////////////
        //fStartBin = max(fStartBin,signalBinCenter - fNBinsLeft );
        //fStopBin = min(fStopBin, signalBinCenter + fNBinsRight );
        
        // Now do the fit with signal template and extract the normalisation it chooses.
        double lowPoint = normalizedHistogramToFit.GetBinLowEdge(fStartBin);
        double highPoint = normalizedHistogramToFit.GetBinLowEdge(fStopBin+1);
        
        //%%%%preliminary fits...
        
        // Preliminary ROOT fit to get parameters near expected values
        for (int i=0; i<5; i++)  normalizedHistogramToFit.Fit(fFitFunction,"R0q","",lowPoint,highPoint);
        double * ROOTFitParams = fFitFunction->GetParameters();
        
        //do a few more preliminary ROOT fits to help convergence
        if(fPrintLevel > 0) std::cout << "Preliminary fits " << fDoExtraPreliminaryFits << std::endl;
        if (fDoExtraPreliminaryFits) {
            
            // Preliminary ROOT fit to get parameters near expected values
            int nGoodPreliminaryFits = 0;
            int nRequiredGoodPreliminaryFits = 5;
            //this will be a vector of parameters, one per function, index-parallel
            vector < vector <double> > allPrelimParamVals;
            int bestPreliminaryFitIndex = -999;
            double bestPreliminaryFitChi2 = 999;
            
            while (nGoodPreliminaryFits < nRequiredGoodPreliminaryFits) {
                std::cout << "Preliminary fits" << std::endl;
                //fit
                TFitResultPtr prelimFitPtr = normalizedHistogramToFit.Fit(fFitFunction,"R0Sq","",lowPoint,highPoint);
                //is the fit good? if not, randomize the input values and do it again
                if (!prelimFitPtr->IsValid()) {
                    fudge(fFitFunction, fRnd, fFitFunction->GetNpar());
                    continue;
                }
                else {
                    nGoodPreliminaryFits++;
                    //save parameters and chi2/NDF
                    vector<double> thisPrelimParamVals;
                    for (int iLoop=0;iLoop<fFitFunction->GetNpar();iLoop++){
                        thisPrelimParamVals.push_back(fFitFunction->GetParameter(iLoop));
                    }
                    allPrelimParamVals.push_back(thisPrelimParamVals);
                    //check whether this fit is best
                    double thisFitChi2 = prelimFitPtr->Chi2()/prelimFitPtr->Ndf();
                    if (bestPreliminaryFitChi2 > thisFitChi2) {
                        std::cout << "This fit is best!" << std::endl;
                        std::cout << "Chi2 before/after" << bestPreliminaryFitChi2 << ", " << thisFitChi2 << std::endl;
                        bestPreliminaryFitChi2 = prelimFitPtr->Chi2()/prelimFitPtr->Ndf();
                        bestPreliminaryFitIndex = nGoodPreliminaryFits-1;//subtract 1 because we increase this counter at the start of the else
                    }
                }//end else on good fit
            }//end number of required preliminary fits achieved
            //pick up the best of the fits
            for (int iLoop=0;iLoop<fFitFunction->GetNpar();iLoop++){
                fFitFunction->SetParameter(iLoop, allPrelimParamVals[bestPreliminaryFitIndex][iLoop]);
            }//end loop on parameters
        }//end loop on do preliminary fits
        
        fFixNormalisation = -2.0;
        MinuitMLFit();
        fFixNormalisation = fMinuitMinimizer->X()[fSignalParameterIndex];
        
        if(fPrintLevel > 0) std::cout << "After fit with signal, signal normalization set to " << fFixNormalisation << endl;
        //fFitFunction->FixParameter(fSignalParameterIndex,fFixNormalisation);
        
    } // end of SWIFT plus signal
    
    //-----------------------------------------------
    // Back to the SWIFT fit itself!
    if(fPrintLevel > 0) std::cout << "About to fit!" << std::endl;
    if(fPrintLevel > 0) std::cout << "Going to scan bin range " << lowBinOfEstimate << " " << highBinOfEstimate << std::endl;
    if(fPrintLevel > 0)  std::cout << "Corresponding to " << fBasicHistogramToFit.GetBinLowEdge(lowBinOfEstimate) << " @ " << fBasicHistogramToFit.GetBinLowEdge(highBinOfEstimate+1) << std::endl;
    if(fPrintLevel > 0)  std::cout << "Using permitted fit range " << lowBinOfFitRange << " " << highBinOfFitRange << std::endl;
    if(fPrintLevel > 0)  std::cout << "Corresponding to " << fBasicHistogramToFit.GetBinLowEdge(lowBinOfFitRange) << " @ " << fBasicHistogramToFit.GetBinLowEdge(highBinOfFitRange+1) << std::endl;
    
    // Make a graph to store my parameter values
    latestFitPars.clear();
    int npar = fFitFunction->GetNpar();
    for (int par = 0; par < npar; par++) {
        TGraph fitPars(highBinOfEstimate-lowBinOfEstimate+1);
        latestFitPars.push_back(fitPars);
    }
    
    // Make something to store the bins included in the fit
    int nBinsInHist = fBasicHistogramToFit.GetNbinsX();
    const TArrayD * ourbins = fBasicHistogramToFit.GetXaxis()->GetXbins();
    if (ourbins) {
        TH2D thisHist("usedSwiftBins","usedSwiftBins",nBinsInHist,ourbins->GetArray(),nBinsInHist,ourbins->GetArray());
        usedFitWindows = thisHist;
        TH2D thatHist("residualsSwiftBins","residualsSwiftBins",nBinsInHist,ourbins->GetArray(),nBinsInHist,ourbins->GetArray());
        residualSwiftBins = thatHist;

    } else {
        TH2D thisHist("usedSwiftBins","usedSwiftBins",nBinsInHist,fBasicHistogramToFit.GetBinLowEdge(1),fBasicHistogramToFit.GetBinLowEdge(fBasicHistogramToFit.GetNbinsX()+1),nBinsInHist,fBasicHistogramToFit.GetBinLowEdge(1),fBasicHistogramToFit.GetBinLowEdge(fBasicHistogramToFit.GetNbinsX()+1));
        usedFitWindows = thisHist;
        TH2D thatHist("residualsSwiftBins","residualsSwiftBins",nBinsInHist,fBasicHistogramToFit.GetBinLowEdge(1),fBasicHistogramToFit.GetBinLowEdge(fBasicHistogramToFit.GetNbinsX()+1),nBinsInHist,fBasicHistogramToFit.GetBinLowEdge(1),fBasicHistogramToFit.GetBinLowEdge(fBasicHistogramToFit.GetNbinsX()+1));
        residualSwiftBins = thatHist;
    }
    
    //  std::cout << "High bin of estimate is " << highBinOfEstimate << std::endl;
    
    fLatestFitStatus = false;
    int lastStartBin = -1; int lastStopBin = -1;
    int testID = -1;
    for (int binToEstimate = lowBinOfEstimate; binToEstimate < highBinOfEstimate+1; binToEstimate++) {
        
        testID++;
        
        // If SWIFT wasn't turned on, we actually just want to do a global fit.
        // We already set the right parameters outside.
        // Otherwise, calculate correct range for this fit.
        if (fDoSwift) {
            
            // No need to restore defaults -- use ending parameter values of previous fit for this one.
            
            // Define fit range to use. Set fStartBin and fStopBin in each case,
            // as these specify what range to use in actual fitting process later.
            
            // Getting low window edge. Default: either edge of range on left or current
            // position minus window width on left
            int lowBinToInclude = max(binToEstimate - fNBinsLeft, lowBinOfFitRange);
            // If we are truncating on the high edge and using the fixWidthAtHighEnd option,
            // this will actually be farther to the left than anticipated.
            if (fTruncateHigh && fFixWidthAtHighEnd) {
                if (highBinOfFitRange - binToEstimate < fNBinsRight) {
                    lowBinToInclude = max(highBinOfFitRange - fNBinsLeft - fNBinsRight, lowBinOfFitRange);
                    std::cout << "Truncating at high mjj and fixing width: reset lowBinToInclude to " << lowBinToInclude << std::endl;
                }
            }
            
            // Getting high window edge. Default: either edge of range on right or current
            // position plus window width on right
            int highBinToInclude = min(binToEstimate + fNBinsRight, highBinOfFitRange);
            // If we are truncating on the high edge, this is instead the min with
            // respect to highBinOfEstimate
            if (fTruncateHigh) highBinToInclude = min(highBinToInclude,highBinOfEstimate);
            // If we are near the low edge and using the fixWidthAtLowEnd option,
            // this will be farther to the right than anticipated.
            if (fFixWidthAtLowEnd) {
                if (binToEstimate - fNBinsLeft < lowBinOfFitRange)
                    highBinToInclude = lowBinOfFitRange + fNBinsLeft + fNBinsRight;
            }
            
            // Now we know what bins to fit between: lowBinToInclude, highBinToInclude.
            // Set the fit to use these.
            fStartBin = lowBinToInclude;
            fStopBin = highBinToInclude;
            
            // Store which bins were used
            for (int ibin=1; ibin < nBinsInHist+1; ibin++) {
                if ((ibin < fStartBin) || (ibin > fStopBin)) usedFitWindows.SetBinContent(ibin,binToEstimate,0.0);
                else usedFitWindows.SetBinContent(ibin,binToEstimate,1.0);
            }
            
            //debug
            if (fPrintLevel > 0) {
                std::cout << "Set sliding window range: " << fStartBin << " " << fStopBin << std::endl;
                std::cout << "Which is " << fBasicHistogramToFit.GetBinLowEdge(fStartBin) << " " << fBasicHistogramToFit.GetBinLowEdge(fStopBin+1) << std::endl;
            }
        }
        
        // If range is identical to last window (i.e. special edge cases)
        // or is a global fit, then we don't need to refit --
        // save time and take the last fit result
        if (!((fStartBin==lastStartBin) && (fStopBin==lastStopBin))) {
            
            //debug
            if (fPrintLevel > 0) {
                std::cout << "Start preliminary ROOT fits: " << std::endl;
            }
            
            // Preliminary ROOT fit to get parameters near expected values
            double lowPoint = normalizedHistogramToFit.GetBinLowEdge(fStartBin);
            double highPoint = normalizedHistogramToFit.GetBinLowEdge(fStopBin+1);
            
            // Preliminary ROOT fit to get parameters near expected values
            for (int i=0; i<5; i++)  normalizedHistogramToFit.Fit(fFitFunction,"R0q","",lowPoint,highPoint);
            double * ROOTFitParams = fFitFunction->GetParameters();
            
            //do a few more preliminary ROOT fits to help convergence
            if (fDoExtraPreliminaryFits) {
                
                // Preliminary ROOT fit to get parameters near expected values
                int nGoodPreliminaryFits = 0;
                int nRequiredGoodPreliminaryFits = 5;
                //this will be a vector of parameters, one per function, index-parallel
                vector < vector <double> > allPrelimParamVals;
                int bestPreliminaryFitIndex = -999;
                double bestPreliminaryFitChi2 = 999;
                
                while (nGoodPreliminaryFits < nRequiredGoodPreliminaryFits) {
                    //std::cout << "Preliminary fits" << std::endl;
                    //fit
                    TFitResultPtr prelimFitPtr = normalizedHistogramToFit.Fit(fFitFunction,"R0S","", lowPoint, highPoint);//add q
                    //is the fit good? if not, randomize the input values and do it again
                    if (!prelimFitPtr->IsValid()) {
                        std::cout << "fudging" << std::endl;
                        fudge(fFitFunction, fRnd, fFitFunction->GetNpar());
                        continue;
                    }
                    else {
                        nGoodPreliminaryFits++;
                        //save parameters and chi2/NDF
                        vector<double> thisPrelimParamVals;
                        for (int iLoop=0;iLoop<fFitFunction->GetNpar();iLoop++){
                            thisPrelimParamVals.push_back(fFitFunction->GetParameter(iLoop));
                        }
                        allPrelimParamVals.push_back(thisPrelimParamVals);
                        //check whether this fit is best
                        double thisFitChi2 = prelimFitPtr->Chi2()/prelimFitPtr->Ndf();
                        if (bestPreliminaryFitChi2 > thisFitChi2) {
                            //std::cout << "This fit is best!" << std::endl;
                            //std::cout << "Chi2 before/after" << bestPreliminaryFitChi2 << ", " << thisFitChi2 << std::endl;
                            bestPreliminaryFitChi2 = prelimFitPtr->Chi2()/prelimFitPtr->Ndf();
                            bestPreliminaryFitIndex = nGoodPreliminaryFits-1;//subtract 1 because we increase this counter at the start of the else
                        }
                    }//end else
                    
                }//end while
                //pick up the best of the fits
                for (int iLoop=0;iLoop<fFitFunction->GetNpar();iLoop++){
                    std::cout << "Setting parameter: " << allPrelimParamVals[bestPreliminaryFitIndex][iLoop] << std::endl;
                    fFitFunction->SetParameter(iLoop, allPrelimParamVals[bestPreliminaryFitIndex][iLoop]);
                }
                
            }//end doExtraPreliminaryFits
            
            if (fPrintLevel > -1) {
                std::cout << "In SwiftFit without signal, setting function parameters to:" << std::endl;
                fFitFunction->Print("all"); // for debugging purposes
                for (int iLoop=0;iLoop<fFitFunction->GetNpar();iLoop++){
                    std::cout << fFitFunction->GetParameter(iLoop) << std::endl;
                }
            }
            
            // Use preliminary fit parameters from function
            // ML fit using Minuit minimizer
            fLatestFitStatus = MinuitMLFit();

            // If it didn't work, keep trying up to
            // 5 times until it does.
            fRetry = 0;
            while (!fLatestFitStatus && fRetry<5) {
                fRetry++;
                std::cout << "Retry number " << fRetry << std::endl;
                fFitFunction->SetParameters(ROOTFitParams);
                bool status = MinuitMLFit();
                if (!(status)) fLatestFitStatus = false;
            }
            
            if (fLatestFitStatus == false) std::cout << "Fit failed to converge after five retries." << std::endl;
            
            if(fPrintLevel > 0) std::cout << "Signal normalization at " << fFixNormalisation << endl;
            
            //      std::cout << "Taking new fit result for bin: " << binToEstimate << std::endl;
            
            // Mark these edges to use next time.
            lastStartBin = fStartBin;
            lastStopBin = fStopBin;
        } else {
            
            //          std::cout << "Taking last fit result for this window too: " << binToEstimate << std::endl;
        }
        
        // Now fill appropriate bin of histogram.
        //    std::cout << "Filling blankhist in between bins " << binToEstimate <<  " " << binToEstimate+1 << std::endl;
        fMjjFitFunction.MakeHistFromFunction(blankhist,
                                             fBasicHistogramToFit.GetBinLowEdge(binToEstimate),
                                             fBasicHistogramToFit.GetBinLowEdge(binToEstimate+1),
                                             false);
        
        if (fDoSwiftChecks) {
            // Fill all histogram bins in 2D with the residual data-function/function for this fit
            //loop on all the bins in the histogram
            for (int ibin=1; ibin < nBinsInHist+1; ibin++) {
                //check whether they are used in swift, if not set to zero
                //binToEstimate is the bin where we are taking the swift estimation (y axis), while ibin is the bin of the original histogram to fit (x axis)
                if ((ibin < fStartBin) || (ibin > fStopBin)) residualSwiftBins.SetBinContent(ibin,binToEstimate,0.0);
                // if not, fill with data minus the integral of the function, divided by sqrt(data)
                else {
                    double x1 = blankhist->GetBinLowEdge(ibin);
                    double x2 = x1 + blankhist->GetBinWidth(ibin);
                    double eventsGenerated = fFitFunction->Integral(x1,x2);
                    if (std::isnan(eventsGenerated)) eventsGenerated = 0;
                    double data = fBasicHistogramToFit.GetBinContent(ibin);
                    if (data < 1) residualSwiftBins.SetBinContent(ibin,binToEstimate, 0, 0); //assume we're dealing with data, unweighted...
                    else residualSwiftBins.SetBinContent(ibin,binToEstimate,(data-eventsGenerated)/sqrt(data));
                    
                    //std::cout << "binToEstimate:" << binToEstimate << std::endl;
                    //std::cout << "ibin:" << ibin << std::endl;
                    //std::cout << "x1:" << x1 << std::endl;
                    //std::cout << "x2:" << x2 << std::endl;
                    //std::cout << "data:" << data << std::endl;
                    //std::cout << "bkg:" << eventsGenerated <<  std::endl;
                    //std::cout << "residual:" << (data-eventsGenerated)/sqrt(data) << std::endl;
                }
            }//end loop on bins for filling residual histogram
        }
        // Store parameters used
        for (int ipar=0; ipar < fFitFunction->GetNpar(); ipar++) {
            latestFitPars.at(ipar).SetPoint(testID,fBasicHistogramToFit.GetBinLowEdge(binToEstimate),fFitFunction->GetParameter(ipar));
        }
        //store various fit outcomes
        g_chi2NDF.SetPoint(testID,fBasicHistogramToFit.GetBinLowEdge(binToEstimate),fFitFunction->GetChisquare()/fFitFunction->GetNDF());
        
        // TEST
        // Extend filled bins to the end of the range used in the fit. This turns out
        // to be important for toy throwing later. Just use last SWIFT parameters.
        if (fDoSwift && highBinOfFitRange > highBinOfEstimate) {
            for (int bin=highBinOfEstimate; bin < highBinOfFitRange+1; bin++) {
                fMjjFitFunction.MakeHistFromFunction(blankhist,
                                                     fBasicHistogramToFit.GetBinLowEdge(bin),
                                                     fBasicHistogramToFit.GetBinLowEdge(bin+1),
                                                     false);
            }
        }
        
        
    }
    
    return fLatestFitStatus;
    
}

// ---------------------------------------------------------
MjjHistogram MjjFitter::FitAndGetBkgWithNoErr(MjjFitFunction & fitFunction, MjjHistogram & histToFit, float minXForFit, float maxXForFit)
{
    
    TH1D backgroundFromFunc((TH1D) histToFit.GetHistogram());
    TString bkgname(Form("%s_internal_bkg",backgroundFromFunc.GetName()));
    backgroundFromFunc.SetName(bkgname);
    backgroundFromFunc.Reset(0);
    TH1D weights((TH1D) histToFit.GetWeightsHistogram());
    TString weightname(Form("%s_internal_weights",weights.GetName()));
    weights.SetName(weightname);
    
    // If minXForFit, maxXForFit are set, they can overwrite previously
    // set values of fMinXForFit, fMaxXForFit.
    // Otherwise, take existing values (even if NaN)
    if (!(std::isnan(minXForFit))) fMinXForFit = minXForFit;
    if (!(std::isnan(maxXForFit))) fMaxXForFit = maxXForFit;
    
    fitFunction.RestoreParameterDefaults();
    SwiftFit(histToFit,fitFunction,&backgroundFromFunc);
    
    //  fitFunction.MakeHistFromFunction(&backgroundFromFunc);
    MjjHistogram backgroundNoErr(&backgroundFromFunc);
    backgroundNoErr.SetEffectiveFromBasicAndWeights(&weights);
    
    fFitSuccessRate = (double) fLatestFitStatus;
    
    if (fFitWithSignal) std::cout << "Finally selected nSignal: " << fMinuitMinimizer->X()[fSignalParameterIndex] << std::endl;
    
    return backgroundNoErr;
    
}

// ---------------------------------------------------------
MjjHistogram MjjFitter::FitAndGetBkgWithMCErr(MjjFitFunction & fitFunction, MjjHistogram & histToFit,float minXForFit, float maxXForFit, bool resetDefaults)
{
    
    TH1D backgroundFromFunc((TH1D) histToFit.GetHistogram());
    TString bkgname(Form("%s_internal_bkg",backgroundFromFunc.GetName()));
    backgroundFromFunc.SetName(bkgname);
    backgroundFromFunc.Reset(0);
    TH1D weights((TH1D) histToFit.GetWeightsHistogram());
    TString weightname(Form("%s_internal_weights",weights.GetName()));
    weights.SetName(weightname);
    
    fMinXForFit = minXForFit;
    fMaxXForFit = maxXForFit;
    
    if (resetDefaults) fitFunction.RestoreParameterDefaults();
    SwiftFit(histToFit,fitFunction,&backgroundFromFunc);
    
    //  fitFunction.MakeHistFromFunction(&backgroundFromFunc);
    for (int i=1; i<=backgroundFromFunc.GetNbinsX(); i++) {
        double thisBinErr = sqrt(backgroundFromFunc.GetBinContent(i));
        backgroundFromFunc.SetBinError(i,thisBinErr);
    }
    MjjHistogram backgroundMCErr(&backgroundFromFunc);
    backgroundMCErr.SetEffectiveFromBasicAndWeights(&weights);
    
    fFitSuccessRate = (double) fLatestFitStatus;
    return backgroundMCErr;
    
}

// ---------------------------------------------------------
MjjHistogram MjjFitter::FitAndGetBkgWithDataErr(MjjFitFunction & fitFunction, MjjHistogram & histToFit, int nFluctuations,float minXForFit, float maxXForFit)
{
    
    if (fPrintLevel > 0 ) {
        std::cout << "Entering FitAndGetBkgWithDataErr" << std::endl;
    }
    
    //int successes = nFluctuations+1;
    //int total = nFluctuations+1;
    int successes = nFluctuations;
    int total = nFluctuations;


    
    TH1D background((TH1D) histToFit.GetHistogram());
    TString bkgname(Form("%s_internal_bkg",background.GetName()));
    background.SetName(bkgname);
    background.Reset(0);
    TH1D weights((TH1D) histToFit.GetWeightsHistogram());
    TString weightname(Form("%s_internal_weights",weights.GetName()));
    weights.SetName(weightname);
    
    fMinXForFit = minXForFit;
    fMaxXForFit = maxXForFit;

    std::cout << "debug "  << " " << fMinXForFit << " " << fMaxXForFit <<std::endl;
    
    fitFunction.RestoreParameterDefaults();
    bool didImportantFitWork;
    didImportantFitWork = SwiftFit(histToFit,fitFunction,&background);
    vector<double> storeImportantCovarianceMatrix = fLatestCovarianceMatrix;
    vector<double> correctFitParameters = fitFunction.GetCurrentParameterValues();
    
    if (!didImportantFitWork) background.Fatal("MjjFitter::FitAndGetBkgWithDataErr","Main fit failed!");

    std::cout << "debug 2" << std::endl;
    
    //  fitFunction.MakeHistFromFunction(&background);
    
    // Use MjjHistogram to get proper fluctuations.
    MjjHistogram bkgTemplate(&background,11);
    bkgTemplate.PrintSeed("bkg error");
    bkgTemplate.SetEffectiveFromBasicAndWeights(&weights);
    
    // Beginning of calculation of error bars on background.
    ///////////////////////////////////////////////////////////////////////////////
    
    // This construction seems complicated but it's a good way to rearrange data
    // from ordering by fluctuation to ordering by bin
    fStoreBinContentVectors.clear();
    for (int i=1; i<=background.GetNbinsX(); i++) {
        vector<double> veccontent;
        veccontent.clear();
        fStoreBinContentVectors.push_back(veccontent);
    }
    
    
    // Create histograms for pseudoexperiments
    TH1D fluctData(background);
    TString fluctdataname(Form("%s_internal_pseudodata",fluctData.GetName()));
    fluctData.SetName(fluctdataname);
    
    fluctData.Clear();
    
    for (int i = 0; i<nFluctuations; i++) {
        
        bkgTemplate.PoissonFluctuateBinByBin(&fluctData);
        MjjHistogram fluctDataHist(&fluctData);
        MjjHistogram fluctBkgHist = FitAndGetBkgWithNoErr(fitFunction,fluctDataHist);
        
        TH1D fluctBkg((TH1D)fluctBkgHist.GetHistogram());
        TString fluctbkgname(Form("%s_internal_pseudobkg",fluctBkg.GetName()));
        fluctBkg.SetName(fluctbkgname);
        std::cout<<"PE #: "<<i<<std::endl;
        // If fit fails don't count this.
        if (fLatestFitStatus==false) {
            nFluctuations++;
            successes--;
            continue;
        }
        else{
            if (fluctBkgHist.GetHistogram().Integral() <= 0.){
                nFluctuations++;
                successes--;
                continue;
            }
        }
        // Lydia
        std::cout<<"In DataErr"<<std::endl;
        //
        // Store bin contents in vector<double> corresponding to this bin
        for (int j=1; j<=fluctBkg.GetNbinsX(); j++) {
            fStoreBinContentVectors.at(j-1).push_back(fluctBkg.GetBinContent(j));
        }
        
        
    }
    
    double RMSOfThisBin;
    for (int i=1; i<=background.GetNbinsX(); i++) {
        if (fStoreBinContentVectors.at(i-1).size() == 0) background.SetBinError(i,sqrt(background.GetBinContent(i)));
        else {
            RMSOfThisBin = GetRMS(fStoreBinContentVectors.at(i-1));
            background.SetBinError(i,RMSOfThisBin);
        }
    }
    ///////////////////////////////////////////////////////////////////////////////
    // End of calculation of error bars on background.
    
    // Make an MjjHistogram
    MjjHistogram theResult(&background);
    theResult.SetEffectiveFromBasicAndWeights(&weights);
    
    // Restore values from primary fit
    fitFunction.SetCurrentParameterValues(correctFitParameters);
    fLatestFitStatus = didImportantFitWork;  // Keep result of main fit
    fLatestCovarianceMatrix = storeImportantCovarianceMatrix; // Keep result of main fit
    fFitSuccessRate = (double)successes/(double)total;
    
    if (fFitWithSignal) std::cout << "Finally selected nSignal: " << fMinuitMinimizer->X()[fSignalParameterIndex] << std::endl;
    
    return theResult;
    
}

// ---------------------------------------------------------
MjjHistogram MjjFitter::FitAndGetBkgWithFitDiffErr(MjjFitFunction & nominal, MjjFitFunction & alternate, MjjHistogram & histToFit,float minXForFit, float maxXForFit, int nFluctuations, bool throwFromData)
{
    int successes = nFluctuations+1;
    int total = nFluctuations+1;
    
    TH1D background((TH1D) histToFit.GetHistogram());
    TString bkgname(Form("%s_internal_bkg",background.GetName()));
    background.SetName(bkgname);
    background.Reset(0);
    TH1D weights((TH1D) histToFit.GetWeightsHistogram());
    TString weightname(Form("%s_internal_weights",weights.GetName()));
    weights.SetName(weightname);
    
    fMinXForFit = minXForFit;
    fMaxXForFit = maxXForFit;
    
    nominal.RestoreParameterDefaults(); alternate.RestoreParameterDefaults();
    bool didImportantFitWork;
    didImportantFitWork = SwiftFit(histToFit,nominal,&background);
    vector<double> storeImportantCovarianceMatrix = fLatestCovarianceMatrix;
    vector<double> correctFitParameters = nominal.GetCurrentParameterValues();
    
    if (!didImportantFitWork) background.Fatal("MjjFitter::FitAndGetBkgWithDiffErr","Main fit failed!");
    
    //  nominal.MakeHistFromFunction(&background);
    
    // Use MjjHistogram to get proper fluctuations.
    MjjHistogram bkgTemplate(&background, 22);
    bkgTemplate.SetEffectiveFromBasicAndWeights(&weights);
    
    // Beginning of calculation of error bars on background.
    ///////////////////////////////////////////////////////////////////////////////
    
    // This construction seems complicated but it's a good way to rearrange data
    // from ordering by fluctuation to ordering by bin
    fStoreBinContentVectors.clear();
    for (int i=1; i<=background.GetNbinsX(); i++) {
        vector<double> veccontent;
        veccontent.clear();
        fStoreBinContentVectors.push_back(veccontent);
    }
    
    //print seed
    if (throwFromData) histToFit.PrintSeed("function choice error (histToFit)");
    else bkgTemplate.PrintSeed("function choice error (bkgTemplate)");
    
    // Create histograms for pseudoexperiments
    TH1D fluctData(background);
    TString fluctdataname(Form("%s_internal_pseudodata",fluctData.GetName()));
    fluctData.SetName(fluctdataname);
    
    fluctData.Clear();
    for (int i = 0; i<nFluctuations; i++) {
        
        std::cout << "Doing new error on PE #" << i << std::endl;
        //    std::cout << "Using throwFromData="<<throwFromData << std::endl;
        std::cout << "At this point, fMaxXForFit is " << fMaxXForFit << std::endl;
        
        if (throwFromData) histToFit.PoissonFluctuateBinByBin(&fluctData);
        else bkgTemplate.PoissonFluctuateBinByBin(&fluctData);
        MjjHistogram fluctDataHist(&fluctData);
        MjjHistogram fluctBkgHistNom = FitAndGetBkgWithNoErr(nominal,fluctDataHist);
        bool fit1 = fLatestFitStatus;
        MjjHistogram fluctBkgHistVar = FitAndGetBkgWithNoErr(alternate,fluctDataHist);
        // If either fit fails don't count this.
        if (fit1==false || fLatestFitStatus==false) {
            nFluctuations++;
            successes--;
            continue;
        }
        else{
            if (fluctBkgHistNom.GetHistogram().Integral() <= 0. || fluctBkgHistVar.GetHistogram().Integral() <= 0.){
                nFluctuations++;
                successes--;
                continue;
            }
        }
        
        // Store difference between two fits in vector<double> corresponding to this bin
        for (int j=1; j<=fluctBkgHistNom.GetHistogram().GetNbinsX(); j++)
            fStoreBinContentVectors.at(j-1).push_back(fluctBkgHistNom.GetHistogram().GetBinContent(j) - fluctBkgHistVar.GetHistogram().GetBinContent(j));
        
        
    }
    
    double RMSOfThisBin;
    for (int i=1; i<=background.GetNbinsX(); i++) {
        if (fStoreBinContentVectors.at(i-1).size() == 0) background.SetBinError(i,sqrt(background.GetBinContent(i)));
        else {
            RMSOfThisBin = GetRMS(fStoreBinContentVectors.at(i-1));
            background.SetBinError(i,RMSOfThisBin);
        }
    }
    ///////////////////////////////////////////////////////////////////////////////
    // End of calculation of error bars on background.
    
    // Make an MjjHistogram
    MjjHistogram theResult(&background);
    theResult.SetEffectiveFromBasicAndWeights(&weights);
    
    // Restore values from primary fit
    nominal.SetCurrentParameterValues(correctFitParameters);
    fLatestFitStatus = didImportantFitWork;  // Keep result of main fit
    fLatestCovarianceMatrix = storeImportantCovarianceMatrix; // Keep result of main fit
    fFitSuccessRate = (double)successes/(double)total;
    
    return theResult;
    
}

// ---------------------------------------------------------
MjjHistogram MjjFitter::FitAndGetBkgWithConfidenceBand(MjjFitFunction & fitFunction, MjjHistogram & histToFit,float minXForFit, float maxXForFit)
{
    
    TH1D data((TH1D) histToFit.GetHistogram());
    TString dataname(Form("%s_internal_bkg",data.GetName()));
    data.SetName(dataname);
    TH1D weights((TH1D) histToFit.GetWeightsHistogram());
    TString weightname(Form("%s_internal_weights",weights.GetName()));
    weights.SetName(weightname);
    
    fMinXForFit = minXForFit;
    fMaxXForFit = maxXForFit;
    
    fitFunction.RestoreParameterDefaults();
    SwiftFit(histToFit,fitFunction);
    
    TF1 * fittedFunction = fitFunction.GetFitFunction();
    
    double * covarianceMatrix = new double [fittedFunction->GetNpar()*fittedFunction->GetNpar()];
    fMinuitMinimizer->GetCovMatrix(covarianceMatrix);
    
    std::cout << "Using xmin, xmax: " << fitFunction.GetMinXVal() << " " << fitFunction.GetMaxXVal() << std::endl;
    TH1D backgroundWithErrors(*GetConfidenceBand(&data, fittedFunction, covarianceMatrix, fitFunction.GetMinXVal(), fitFunction.GetMaxXVal(), 0.683));
    MjjHistogram result(&backgroundWithErrors);
    result.SetEffectiveFromBasicAndWeights(&weights);
    
    fFitSuccessRate = (double) fLatestFitStatus;
    return result;
    
}

// ---------------------------------------------------------
void MjjFitter::StoreLatestFitMatrices(int dimensionOfMatrix)
{
    
    fLatestCovarianceMatrix.clear();
    fLatestHessianMatrix.clear();
    
    double hessElem [dimensionOfMatrix*dimensionOfMatrix];
    fMinuitMinimizer->GetHessianMatrix(hessElem);
    
    // Store items of covariance matrix in the order:
    // Left to right, starting on the top row, then second, etc.
    for (int i=0; i<dimensionOfMatrix; i++) {
        for (int j=0; j<dimensionOfMatrix; j++) {
            double thisElement = fMinuitMinimizer->CovMatrix(i,j);
            fLatestCovarianceMatrix.push_back(thisElement);
            fLatestHessianMatrix.push_back(hessElem[i*dimensionOfMatrix + j]);
        }
    }
    
    return;
    
}

// ---------------------------------------------------------
TH1D * MjjFitter::GetConfidenceBand(TH1D * data, TF1 * fittedFunction,  double * covarianceMatrix, double xmin, double xmax, double CL)
{
    
    int npar = fittedFunction->GetNpar();
    double * grad = new double [npar];
    double x[1];
    int binfirst, binlast;
    if (xmin < 0) binfirst = data->GetXaxis()->GetFirst();
    else binfirst = data->FindBin(xmin+1);
    if (xmax < 0) binlast = data->GetXaxis()->GetLast();
    else binlast = data->FindBin(xmax-1);
    double t = TMath::StudentQuantile(0.5 + CL/2, fittedFunction->GetNDF());
    double chidf = TMath::Sqrt(fittedFunction->GetChisquare()/fittedFunction->GetNDF());
    
    TH1D * answer = (TH1D *) data->Clone();
    
    for (int bin=1; bin < answer->GetNbinsX()+1; bin++) {
        if (bin>binfirst-1 && bin<binlast+1) {
            x[0] = data->GetXaxis()->GetBinCenter(bin);
            fittedFunction->GradientPar(x,grad);
            double sum = 0;
            for (int i=0; i<npar; i++) {
                for (int j=0; j<npar; j++) {
                    sum += covarianceMatrix[i*npar+j]*grad[i]*grad[j];
                }
            }
            double c = sqrt(sum);
            answer->SetBinContent(bin,fittedFunction->Eval(x[0])*answer->GetBinWidth(bin));
            answer->SetBinError(bin,c*t*chidf*answer->GetBinWidth(bin));
        } else {
            answer->SetBinContent(bin,0);
            answer->SetBinError(bin,0);
        }
    }
    delete [] grad;
    return answer;
}

// ---------------------------------------------------------
