// ---------------------------------------------------------

#include "Bayesian/MjjBATBumpHunter.h"

// ---------------------------------------------------------
MjjBATBumpHunter::MjjBATBumpHunter(MjjBATModel * fitter) :
  fSpectrumTomography() {

  // Set defaults as specified in header

  fAllowDeficit = false;
  fUseSidebands = true;
  fMinBinsInBump = 1;
  fMaxBinsInBump = 1e5;
  fNBinsInSideband = 1;
  fNBinsInShift = 1;
  fLowEdge = fHighEdge = 0;
  fLowEdgesAllBumps.clear();
  fHighEdgesAllBumps.clear();
  fProbAllBumps.clear();
  fWillProbBumps.clear();
  differences.clear();
  fFitter = fitter;
  fFitter->SetLikelihoodChoice(BUMPHUNTER);
  
  // If the background process is supposed to use both statistical and normalisation
  // errors, we need to leave the normalisation floating.
  // If on the other hand it is supposed to use normalisation alone, we can treat this
  // the same way we do the stats, because of the one-bin interpretation of the
  // BumpHunter. So we freeze the normalisation parameter to save time.
//  if (!(fFitter->GetProcess(0)->GetDoStatUnc()) && fFitter->GetProcess(0)->GetDoNormalisationUnc()) {
 if (fFitter->GetProcess(0)->GetDoNormalisationUnc()) {
    int index = fFitter->GetProcess(0)->GetNormalisationParamIndex();
//    fFitter->SetParameterRange(index,0.0,0.0);
    fFitter->GetParameter(index)->Fix(0.0);
  }
  
  fFitter->AddParameter("dummyN", 0,1e5);
  fFitter->SetPriorConstant("dummyN");

}

// ---------------------------------------------------------
MjjBATBumpHunter::~MjjBATBumpHunter()
   { }

// ---------------------------------------------------------
double MjjBATBumpHunter::DoTest(MjjHistogram & dataHist, MjjHistogram & bkgHist, int firstBinToUse, int lastBinToUse){

  TH1D h_data(dataHist.GetEffectiveHistogram());
  TString dataname(Form("%s_bh_data",h_data.GetName()));
  h_data.SetName(dataname); 

  TH1D h_bkg(bkgHist.GetEffectiveHistogram());
  TString bkgname(Form("%s_bh_bkg",h_bkg.GetName()));
  h_bkg.SetName(bkgname);

  int nBinsData = h_data.GetNbinsX();
  int nBinsBkg = h_bkg.GetNbinsX();
  assert(nBinsData==nBinsBkg);

  // Find first and last bins with data
  // If reasonable, overwrite with user's choice
  int firstBin = dataHist.GetFirstBinWithData();
  int lastBin = dataHist.GetLastBinWithData();
  if (firstBinToUse>0 && firstBinToUse > firstBin && firstBinToUse < lastBin) firstBin = firstBinToUse;
  if (lastBinToUse > firstBinToUse && lastBinToUse>0 && lastBinToUse > firstBin && lastBinToUse < lastBin) lastBin = lastBinToUse;
  
  fFitter->SetData(dataHist, firstBin, lastBin);

  vector<std::pair<int,int> > regionsdef;
  if (fExcludeWindow) {
    regionsdef.push_back(std::make_pair(firstBin,fFirstBinToExclude-1));
    regionsdef.push_back(std::make_pair(fLastBinToExclude+1,lastBin));
  } else {
    regionsdef.push_back(std::make_pair(firstBin,lastBin));
  }

  // To hold final p-value, low edge, and high edge of most interesting bump  
  fMostInterestingP = 1;
  fLowEdge = fHighEdge = 0;

  // Store low edges, high edges, p-values of all tested bumps
  fLowEdgesAllBumps.clear();
  fHighEdgesAllBumps.clear();
  fProbAllBumps.clear();
  
  fWillProbBumps.clear();

  for (unsigned int index = 0; index < regionsdef.size(); index++) {

    int thisfirstBin = regionsdef.at(index).first;
    int thislastBin = regionsdef.at(index).second;

    // Set number of bins in spectrum, max number of bins in bump
    int nBins = thislastBin-thisfirstBin+1;
    int minWidth = std::max(fMinBinsInBump,1);
    int maxWidth = (int) std::min(fMaxBinsInBump,nBins/2);

    DoTestCore(h_data, h_bkg, minWidth,maxWidth,thisfirstBin,thislastBin);

  }

  // Store fProbAllBumps contents into tomography plot
  TGraphErrors tomography(fProbAllBumps.size());
  TGraphErrors willtomography(fWillProbBumps.size());
  for (unsigned int i=0; i<fProbAllBumps.size(); i++) {
    tomography.SetPoint(i,(fLowEdgesAllBumps.at(i)+fHighEdgesAllBumps.at(i))/2., fProbAllBumps.at(i));
    tomography.SetPointError(i,(fHighEdgesAllBumps.at(i)-fLowEdgesAllBumps.at(i))/2., 0);
    willtomography.SetPoint(i,(fLowEdgesAllBumps.at(i)+fHighEdgesAllBumps.at(i))/2., fWillProbBumps.at(i));
    willtomography.SetPointError(i,(fHighEdgesAllBumps.at(i)-fLowEdgesAllBumps.at(i))/2., 0);
  }
  fSpectrumTomography = tomography;
  fWillTomography = willtomography;

  if (fMostInterestingP==0) FindBumpInCaseOfIncalculable(h_data, h_bkg, firstBin,lastBin);

  return -log(fMostInterestingP);

}

void MjjBATBumpHunter::DoTestCore(TH1D h_data, TH1D h_bkg, int minWidth, int maxWidth, int firstBin, int lastBin) {

  // Declare variables to hold everything during bump hunt:
  double sidebandWidth = 0;
  double smallestPforWidth, lowEdgeForWidth, highEdgeForWidth;
  // Left hand side of leftmost bin in bump at minimum of its range
  int minBinL = 0;
  // Left hand side of leftmost bin in bump at maximum of its range
  int maxBinL = 0;
  // Number of bins by which we shift window between tests
  int nbinsinstep = 0;
  // Left of bump window, right of bump window, left of left sideband, right of right sideband
  int binL , binR , binLL, binRR;
  // Background, background error in center, left sideband, right sideband
  double bC, bL, bR, deltaBC, deltaBL, deltaBR;
  // Data in center, left sideband, right sideband
  double dC, dL, dR;
  // Newly added: data error
  double deltaDC, deltaDL, deltaDR;
  // p-value in left sideband, center, right sideband
  double probL, probC, probR;

  // Initialize all.
  dC = dL = dR = 0.;
  deltaBC = deltaBL = deltaBR = 0.;
  bC = bL = bR = deltaBC = deltaBL = deltaBR = 0.;
  binL = binR = binLL = binRR = 0;
  smallestPforWidth = lowEdgeForWidth = highEdgeForWidth = 0.;
  probL = probC = probR = 1.;

  // Test each bin width in range
  int counter = 0;
  for (int width = minWidth; width < maxWidth+1; width+=1) {
  
    // Use user specification for sideband width if reasonable
    // Else use default of width/2 to a minimum of 1 bin
    if (fNBinsInSideband>=1) sidebandWidth = fNBinsInSideband;
    else sidebandWidth = std::max(1,(int)width/2);

    // Set number of bins by which to shift window in checks
    if (!(fNBinsInShift<1)) nbinsinstep = fNBinsInShift;
    else nbinsinstep = std::max(1,(int)width/2);

    // Reset values before looping over bins
    smallestPforWidth = 1;
    lowEdgeForWidth = 0;
    highEdgeForWidth = 0;
    if (fUseSidebands) {
       minBinL = firstBin + sidebandWidth;
       maxBinL = lastBin - width - sidebandWidth + 1;
    } else {
       minBinL = firstBin;
       maxBinL = lastBin - width + 1;
    }

//    std::cout << "In width " << width << " minBinL is " << minBinL << " where firstBin, lastBin were " << firstBin << ", " << lastBin << std::endl;

    // Loop while leftmost edge of bump is below or equal to leftmost allowed limit
    for (binL = minBinL; binL <= maxBinL; binL += nbinsinstep) {
      
      // Find upper edge of bump, sideband edges
      binR = binL+width-1;
      binLL = binL-sidebandWidth;
      binRR = binR+sidebandWidth;

      // Get probabilities.
      //std::cout << "Doing number " << counter << ": bins " << binL << " - " << binR << std::endl;
      counter++;

      // TEST
      //if (counter == 93); // fFitter->SetDebug(true);
      //else { probC = 0.0; continue; }

//std::cout << "Doing number " << counter << ": bins " << binL << " - " << binR << std::endl;

      probC = probL = probR = 0.0;


      // Get effective data and background content and uncertainty in bump window and sidebands
      // from histogram and errors
      // Use the same functions as before to count nominal background and data.
      GetEffectiveBandContentsWithError(h_data, h_bkg, binL, binR, dC, deltaDC, bC, deltaBC);
      if (fUseSidebands) {
        GetEffectiveBandContentsWithError(h_data, h_bkg, binLL, binL-1, dL, deltaDL, bL, deltaBL);
        GetEffectiveBandContentsWithError(h_data, h_bkg, binR+1, binRR, dR, deltaDR, bR, deltaBR);
      }

      // I can use the nominal background value to determine if this is an excess or not.
      // Determine if it is an excess. Only keep if flagged to do so.
      if (!(fAllowDeficit)) { if (dC<=bC) continue; }

      // TEST
//      fFitter->SetDebug(true);
/*      std::vector<double> pars;
      for (int i=0; i< 150; i++)  {
        pars.clear();
        pars.push_back(0.0); pars.push_back(0.0); pars.push_back(float(i));
        std::cout << "For " << i << " is " << posteriors.at(posteriors.size()-1).GetBinContent(i)/posteriors.at(posteriors.size()-1).Integral() << ". Should be " << exp(fFitter->BHLogLikelihood(pars)) << std::endl;
      } 

      probC = GetBATProbability(binL, binR, dC, bC, deltaBC);
      if (fUseSidebands) {
        probL = GetBATProbability(binLL, binL-1, dL, bL, deltaBC);
        probR = GetBATProbability(binR+1,binRR, dR, bR, deltaBC);
      }

       // TEST
      std::vector<double> pars;
      for (int i=100; i< 105; i++)  {
        pars.clear();
        pars.push_back(0.0); pars.push_back(0.0); pars.push_back(float(i));
        std::cout << "For " << i << " is " << posteriors.at(posteriors.size()-1).GetBinContent(i)/posteriors.at(posteriors.size()-1).Integral() << ". Should be " << exp(fFitter->BHLogLikelihood(pars)) << std::endl;
      }*/


      fFitter->SetDebug(true);
      double willprob = GetPoissonConvGammaPValue(dC, bC, deltaBC);
      fFitter->SetDebug(false);    
      probC = willprob;
 
//      std::cout << "Me = " << probC << "; Will = " << willprob << std::endl;
      if (fabs(probC - willprob) > 0.03) std::cout << "Found a bigger difference!!!!" << std::endl;
      if (fabs(probC - willprob) > 0.1) {std::cout << "AAAAAGHGHGHGH" << std::endl; return;}
      differences.push_back(probC - willprob);
 
      // Question: if we want sidebands, do marginalizations have to be done over all region together?
      // Can't do that right now.

      // Ignore cases where a significant discrepancy is observed in sidebands
      if (fUseSidebands) {
        if (probL <= 1e-3 || probR <= 1e-3)
          continue;
      }

      // Save information to use later in tomography plot
      fLowEdgesAllBumps.push_back(h_data.GetBinLowEdge(binL));
      fHighEdgesAllBumps.push_back(h_data.GetBinLowEdge(binR) + h_data.GetBinWidth(binR));
      fProbAllBumps.push_back(probC);
      
      // Save for Will-tomography plot
      fWillProbBumps.push_back(willprob);
      
      if (probC < smallestPforWidth) {
        smallestPforWidth = probC;
        lowEdgeForWidth = h_bkg.GetBinLowEdge(binL);
        highEdgeForWidth = h_bkg.GetBinLowEdge(binR) + h_bkg.GetBinWidth(binR);
      }
    } // next center
    if (smallestPforWidth < fMostInterestingP) {
      fMostInterestingP = smallestPforWidth;
      fLowEdge = lowEdgeForWidth;
      fHighEdge = highEdgeForWidth;
    }
  } // next width

}

double MjjBATBumpHunter::GetBATProbability(int firstBin, int lastBin, double nData, double nBNom, double nBErr) {

  fFitter->SetInterestingBins(firstBin, lastBin);
  
  double scaleby = 6;
  
  int maxN = int(std::max(nBNom + scaleby*nBErr,80.0));
  int minN = int(std::max(nBNom - scaleby*nBErr,0.0));
 
  // Change parameter range and precision to get a good look at this.
//  fFitter->SetParameterRange(fFitter->GetParameter("dummyN")->GetIndex(),minN,maxN);
  fFitter->GetParameter("dummyN")->Unfix();
  fFitter->GetParameter("dummyN")->SetLimits(minN,maxN);
  fFitter->SetNbins(maxN-minN);
  fFitter->ResetResults();
  fFitter->MarginalizeAll();
  BCH1D * marginalizedSig = fFitter->GetMarginalized("dummyN");
  
  TH1D * sigLikeVsNumber = (TH1D*) marginalizedSig->GetHistogram()->Clone();
  TH1D * willPosterior = ((TH1D*) sigLikeVsNumber->Clone());
  willPosterior->Reset();
  posteriors.push_back(*sigLikeVsNumber);
  willposteriors.push_back(*willPosterior);

  // Determine what sum over n we need to give a valuable result.
  int limit = int(nData);
  double sum = 0;
  if (nData < nBNom) { //limit--;
    sum = sigLikeVsNumber->Integral(sigLikeVsNumber->FindBin(minN),sigLikeVsNumber->FindBin(limit));
  } else {
    sum = sigLikeVsNumber->Integral(sigLikeVsNumber->FindBin(limit),sigLikeVsNumber->FindBin(maxN));
  }
  sum /= sigLikeVsNumber->Integral();

  return sum;

}

void MjjBATBumpHunter::FindBumpInCaseOfIncalculable(TH1D h_data, TH1D h_bkg, int firstBin, int lastBin) {

  vector<int> singlebinsinf;
  bool lastBinWasInf = false;
  bool allInfsConsecutive = true;
  for (int bin = firstBin; bin < lastBin+1; bin ++) {
    double D = h_data.GetBinContent(bin);
    double B = h_bkg.GetBinContent(bin);
    double thisbinpval = PoissonPval(D,B);
    if (thisbinpval==0 && D>B) {
      if (singlebinsinf.size()>0 && lastBinWasInf==false) allInfsConsecutive = false;
      singlebinsinf.push_back(bin);
      lastBinWasInf = true;
    } else {
      lastBinWasInf = false;
    }
  }
  
   if (singlebinsinf.size() > 0 && allInfsConsecutive) {
    fMostInterestingP = 0;
    fLowEdge = h_data.GetBinLowEdge(singlebinsinf.at(0));
    fHighEdge = h_data.GetBinLowEdge(singlebinsinf.at(singlebinsinf.size()-1))
                + h_data.GetBinWidth(singlebinsinf.at(singlebinsinf.size()-1));
  }

}

// ---------------------------------------------------------
void MjjBATBumpHunter::GetEffectiveBandContentsWithError(TH1D histData, TH1D histBkg, int firstBin, int lastBin, double& dataIntegral, double& dataError, double& bkgIntegral, double& bkgError) {

  dataIntegral = dataError = bkgIntegral = bkgError = 0;
  for (int bin=firstBin; bin<lastBin+1; bin++) {
    dataIntegral+=histData.GetBinContent(bin);
    dataError+=histData.GetBinError(bin);
    bkgIntegral+=fFitter->GetProcess(0)->GetNominalHistogram()->GetBinContent(bin);
    bkgError+=fFitter->GetProcess(0)->GetVariationHistogram()->GetBinContent(bin);
  }
  return; 
}

// ---------------------------------------------------------
vector<double> MjjBATBumpHunter::GetFurtherInformation() {

  vector<double> lowHighLimits;
  lowHighLimits.clear();
  lowHighLimits.push_back(fLowEdge);
  lowHighLimits.push_back(fHighEdge);
  return lowHighLimits;

}

//return the p-value (and associated zValue = number of sigma) of seeing nObs "events"
//under the model a poisson distribution with the value of the mean parameter being itself uncertain
//. Assumption is that the mean parameter is distributed as a gamma density, with the expectation (mean) of the gamma function = E
// and the variance of the gamma density equal to err^2
double MjjBATBumpHunter::GetPoissonConvGammaPValue(double nObs, double E, double err) {
   //parameters (a,b) of gamma density can be written in terms of the expectation and variance:
   //Gammma(x, a,b);   - note this is not equal to gamma(x) or gamma(a,b), which are different functions
   double b = E/(err*err); // = E/V
   double a = E*b; // = E^2/V

  //std::cout << "N, B, Berr are " << nObs << ", " << E << ", " << err << std::endl;
  //std::cout << "So alpha, beta are " << a << ", " << b << std::endl;

   double pval = 0.;
   //decide if we should ignore systematics or not 
   if (a>100*nObs) {
      //stat error is big enough to ignore the syst error
      //considering only stat error means the p-value is given by:
      // (nData>nMC): pval = sum(x = nData->inf, Poisson(x,nMC)) = 1 - sum(x = 0->nData-1,Poisson(x,nMC))
      // (nData<=nMC): pval = sum(x = 0->nData, Poisson(x,nMC))
      // But sum(x = 0->nData-1,Poisson(x,nMC)) = gamma(nData,nMC)/gamma(nData); <---- see maths websites
      // so we have:
      // (nData>nMC): pval = 1 - gamma(nData,nMC)/gamma(nData);
      // (nData<=nMC): pval = gamma(nData+1,nMC)/gamma(nData);
      // .....And ROOT provides gamma(a,b)/gamma(a) = ROOT::Math::inc_gamma_c(a,b)
      //pval = (nObs<=E) ? ROOT::Math::inc_gamma_c(nObs+1,E) : (1. - ROOT::Math::inc_gamma_c(nObs,E));
      pval = GetPoissonPValue(nObs,E);
   } else {
      //use recursive formula to solve:
      // (nData>nMC): pval = 1 - sum(x=0->nData-1, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
      // (nData<=nMC): pval = sum(x=0->nData, Integral(y=0->inf, Poisson(x,y)Gamma(y,a,b) dy))
      //i.e. integrating out the unknown parameter y
      // Recursive formula: P(n;A,B) = P(n-1;A,B) (A+n-1)/(n*(1+B))
      unsigned stop=nObs;
      if (nObs>E) --stop;
      double sum = 0;
      vector<double> pars;
      if(a>100) {
         /// NB: must work in log-scale otherwise troubles!
         double logProb = a*log(b/(1+b));
         sum=exp(logProb); // P(n=0)
         for (unsigned u=1; u<stop+1; ++u) {
            logProb += log((a+u-1)/(u*(1+b)));
            sum += exp(logProb);
            //std::cout << "in bin " << u << " will adds " << exp(logProb) << std::endl;
            /*pars.clear();
            pars.push_back(0.0); pars.push_back(0.0); pars.push_back(float(u));
            std::cout << ". I add " << exp(fFitter->BHLogLikelihood(pars)) << ". Actually have " << posteriors.at(posteriors.size()-1).GetBinContent(u)/posteriors.at(posteriors.size()-1).Integral() << std::endl; */
            //willposteriors.at(willposteriors.size()-1).SetBinContent(u,exp(logProb));
         }
      } else {
         double p0 = pow(b/(1+b),a); // P(0;A,B)
         double pLast = p0;
         sum = p0;
         for (unsigned k=1; k<stop+1; ++k) {
            double p = pLast * (a+k-1) / (k*(1+b));
            sum += p;pLast = p;
            //std::cout << "in bin " << k << " will adds " << p << std::endl;
            /*pars.clear();
            pars.push_back(0.0); pars.push_back(0.0); pars.push_back(float(k));
            std::cout << ". I add " << exp(fFitter->BHLogLikelihood(pars)) << ". Actually have " << posteriors.at(posteriors.size()-1).GetBinContent(k)/posteriors.at(posteriors.size()-1).Integral() <<  std::endl;*/
            //willposteriors.at(willposteriors.size()-1).SetBinContent(k,p);
        }
      }
      pval = (nObs>E) ?  1-sum : sum; 
   } 

   //bool overflow(false);
   //zValue = GetZValue(pValue,overflow) //large z-values correspond to small p-values.... i.e. significant diff
   //zValue = (nObs<E) ? -1.*zValue : zValue; //flip the z-values of deficits

   return pval;
}

double MjjBATBumpHunter::GetPoissonPValue(double nObs, double E) {
   if(nObs>E) { //excess
      double p = ROOT::Math::inc_gamma_c(nObs,E);
      if(p==1. && nObs>100.) { //excess very excessive, and have a large nObs so no big problem summing one extra term
         return ROOT::Math::inc_gamma(nObs,E);
      } else {
         return 1-p;
      }
   } 
   return ROOT::Math::inc_gamma_c(nObs+1,E); //deficit
}



// ---------------------------------------------------------

