#include "call_libraries.h"  // call libraries from ROOT and C++

// inputs are 2D histograms: reff2D for efficiency, rfak2D for fakes, rsec2D for secondary (decays), rmul2D for multiple reconstruction
// pT and eta are the transverse momentum and pseudorapidity of the track (considering a 2D histogram where X is eta axis and Y pT axis)
double getUnfCorrWeight(TFile *unffile, float leadpt, float subleadpt, int multiplicity, float etaregion){

  double factor = 1.0;
  double eff = 1.0;
  TH2 *eff_factor = nullptr; 
  if(multiplicity >= 10 && multiplicity < 60){
  	  if(etaregion == 0.5){
		  unffile->GetObject("ratios/reco_ratio2D_10-60_mid-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 1.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_mid-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 2.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_mid-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 3.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_fwd-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 4.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_bkw-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 5.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_fwd-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 6.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_fwd-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 7.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_bkw-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 8.5){
 		  unffile->GetObject("ratios/reco_ratio2D_10-60_bkw-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
	  } else {
          unffile->GetObject("ratios/reco_ratio2D_10-60_other", eff_factor);  // data / MC
          eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );	  
	  }
  } else if(multiplicity >= 60 && multiplicity < 120){
  	  if(etaregion == 0.5){
		  unffile->GetObject("ratios/reco_ratio2D_60-120_mid-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 1.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_mid-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 2.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_mid-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 3.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_fwd-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 4.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_bkw-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 5.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_fwd-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 6.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_fwd-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 7.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_bkw-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 8.5){
 		  unffile->GetObject("ratios/reco_ratio2D_60-120_bkw-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
	  } else {
          unffile->GetObject("ratios/reco_ratio2D_60-120_other", eff_factor);  // data / MC
          eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );	  
	  }
  } else if(multiplicity >= 120.0 && multiplicity < 185.0){
  	  if(etaregion == 0.5){
		  unffile->GetObject("ratios/reco_ratio2D_120-185_mid-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 1.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_mid-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 2.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_mid-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 3.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_fwd-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 4.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_bkw-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 5.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_fwd-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 6.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_fwd-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 7.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_bkw-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 8.5){
 		  unffile->GetObject("ratios/reco_ratio2D_120-185_bkw-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
	  } else {
          unffile->GetObject("ratios/reco_ratio2D_120-185_other", eff_factor);  // data / MC
          eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );	  
	  }
  } else if(multiplicity >= 185.0){
  	  if(etaregion == 0.5){
		  unffile->GetObject("ratios/reco_ratio2D_185-250_mid-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 1.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_mid-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 2.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_mid-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 3.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_fwd-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 4.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_bkw-mid", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 5.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_fwd-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 6.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_fwd-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 7.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_bkw-fwd", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
  	  } else if(etaregion == 8.5){
 		  unffile->GetObject("ratios/reco_ratio2D_185-250_bkw-bkw", eff_factor);  // data / MC
  		  eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );
	  } else {
          unffile->GetObject("ratios/reco_ratio2D_185-250_other", eff_factor);  // data / MC
          eff = eff_factor->GetBinContent( eff_factor->GetXaxis()->FindBin(leadpt),eff_factor->GetYaxis()->FindBin(subleadpt) );	  
	  }
  }

  if(eff > 1.0e-20) { factor = eff; } else {factor = 1.0;}

  return factor;

}
