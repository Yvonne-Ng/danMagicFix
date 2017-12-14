// ---------------------------------------------------------

#include "Bayesian/MjjBATShapeChangingSyst.h"

// ---------------------------------------------------------
MjjBATShapeChangingSyst::MjjBATShapeChangingSyst() 
{

}

// ---------------------------------------------------------
MjjBATShapeChangingSyst::~MjjBATShapeChangingSyst() 
{

}

// ---------------------------------------------------------
TMatrixDSparse MjjBATShapeChangingSyst::GetMatrix(double param) 
{

  int index=0;
  for (unsigned int i=0; i<(fCorrespondingSigmas.size()-1); i++) {
    if (param > fCorrespondingSigmas[i+1]) continue;
    index=i;
    break;
  }

  double change = (param-fCorrespondingSigmas.at(index)) / fSigmaDifferences.at(index);
  TMatrixDSparse matrix(fMatrixDifferences.at(index));
  matrix*=change;
  matrix+=fTransferMatrices.at(index); 

  return matrix;

}

// ------------------------------------------
void MjjBATShapeChangingSyst::SetSpectra(vector<std::pair<double,TH2D*> > sigmasAndTMs) 
{

  if (fTransferMatrices.size() > 0 || fCorrespondingSigmas.size() > 0 ) 
          std::cout << "Caution: overwriting old matrices" << std::endl;
  fTransferMatrices.clear();
  fCorrespondingSigmas.clear();

  std::sort(sigmasAndTMs.begin(), sigmasAndTMs.end(), sort_pairs_TH2D());

  // determine structure of histograms from a sample
  int nBinsInSpectra = sigmasAndTMs.at(0).second->GetNbinsX()+2;
  fBinStructure = *sigmasAndTMs.at(0).second->GetXaxis()->GetXbins();

  // create transfer matrices in form TMatrixDSparse
  fStoreRows.clear();
  fStoreColumns.clear();
  fStoreData.clear(); 
  for (unsigned int i=0; i<sigmasAndTMs.size(); i++) {
    TArrayI rows(nBinsInSpectra*nBinsInSpectra);
    TArrayI cols(nBinsInSpectra*nBinsInSpectra);
    TArrayD data(nBinsInSpectra*nBinsInSpectra);
    TMatrixDSparse thisMatrix(nBinsInSpectra,nBinsInSpectra);
    int nNonZero = 0;
    for (int binx=0; binx<nBinsInSpectra; binx++) {
      for (int biny=0; biny<nBinsInSpectra; biny++) {
        double val = sigmasAndTMs.at(i).second->GetBinContent(binx,biny);
        if (val==0) continue;
        rows[nNonZero] = binx;
        cols[nNonZero] = biny;
        data[nNonZero] = val;
        nNonZero++;
      }
    }
    fStoreRows.push_back(rows);
    fStoreColumns.push_back(cols);
    fStoreData.push_back(data);
    thisMatrix.SetMatrixArray(nNonZero,fStoreRows.at(i).GetArray(),
                     fStoreColumns.at(i).GetArray(),fStoreData.at(i).GetArray());
    fTransferMatrices.push_back(thisMatrix);
    fCorrespondingSigmas.push_back(sigmasAndTMs.at(i).first);
  }

  // store differences: cutting off one computation may help.
  fSigmaDifferences.clear();
  fMatrixDifferences.clear();
  for (unsigned int i=0; i<sigmasAndTMs.size()-1; i++) {
    fSigmaDifferences.push_back(fCorrespondingSigmas.at(i+1)-fCorrespondingSigmas.at(i));
    TMatrixDSparse diff(fTransferMatrices.at(i+1));
    diff-=fTransferMatrices.at(i);
    fMatrixDifferences.push_back(diff);
  }

  return;

}

// ------------------------------------------

