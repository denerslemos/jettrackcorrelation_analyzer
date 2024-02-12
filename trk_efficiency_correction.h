#include "call_libraries.h"  // call libraries from ROOT and C++

// Obtain weight from tracking efficiency correction - default used by CMS
// it returns a efficiency correction to be applied in histograms as -> Fill(variable,this_weight)

bool checkBounds(double pt, double eta){
  if( TMath::Abs(eta) > 2.4 ){return false;}
  if( pt < 0 || pt > 500 ){return false;}
  return true;
}

// inputs are 2D histograms: reff2D for efficiency, rfak2D for fakes, rsec2D for secondary (decays), rmul2D for multiple reconstruction
// pT and eta are the transverse momentum and pseudorapidity of the track (considering a 2D histogram where X is eta axis and Y pT axis)
double getTrkCorrWeight(TFile *trkeff_file, bool use_centrality, string system, int year, int energy, int mult, double pT, double eta, double phi){

  if( !checkBounds(pT, eta) ) return 0;
  double factor = 1.0;

  // add it for each system here

  if(system == "pPb" && year == 2016 && energy == 8160){
    TH2 *eff_factor = nullptr; 
    trkeff_file->GetObject("rTotalEff3D_0", eff_factor);  // eff / (1 - fake rate)
    double eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(eta),eff_factor->GetYaxis()->FindBin(pT) );
    factor = (1. / eff); //only efficiency
  }

  return factor;

}
