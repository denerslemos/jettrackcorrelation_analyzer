#include "call_libraries.h"  // call libraries from ROOT and C++

// inputs are 2D histograms: X --> reco pT and Y --> ref pT
// as well as the pT of reco jet (pT)
// and also multiplicity (to be added)

double getUnfolding(TH2* unf_histo, int mult, double pT){

  int bin_number = unf_histo->GetXaxis()->FindBin(pT);
  TH1D* histo1D = (TH1D*) unf_histo->ProjectionY('histo',bin_number,bin_number);
  double factor = histo1D->GetMean();
  return factor;

}
