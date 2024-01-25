#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // call inputs

const int trkbinsize = (int) trk_pt_bins.size(); // track bins for jet-track correlation
const int multbinsize = (int) multiplicity_centrality_bins.size();// multiplicity or centrality bins for jet-track correlation
const int extrabinsize = (int) extra_bins.size();// any additional dependency you wanna add (be carefull about memory)

// binning definition
const double binnerShift = 0.0; // shift if starts at 0, Log(0) -> error
// Needed to define log binning
// Xp and XPb -> see sumw2 function
const int nXBins = 40; // number of bins
const double minX = 3e-04;  // minimum
const double maxX = 1.0; 	   // maximum
double XlogBinWidth = (TMath::Log(maxX+binnerShift) - TMath::Log(minX+binnerShift)) / nXBins; // binwidth
//double XBins[nXBins+1] = {minX, 0.000125893, 0.000158489, 0.000199526, 0.000251189, 0.000316228, 0.000398107, 0.000501187, 0.000630957, 0.000794328, 0.001, 0.00125893, 0.00158489, 0.00199526, 0.00251189, 0.00316228, 0.00398107, 0.00501187, 0.00630957, 0.00794328, 0.01, 0.0125893, 0.0158489, 0.0199526, 0.0251189, 0.0316228, 0.0398107, 0.0501187, 0.0630957, 0.0794328, 0.1, 0.125893, 0.158489, 0.199526, 0.251189, 0.316228, 0.398107, 0.501187, 0.630957, 0.794328, maxX};
// M12 -> see sumw2 function
const int nMBins = 15; // number of bins
const double minM = 70.0;  // minimum
const double maxM = 3000.0; // maximum
double MlogBinWidth = (TMath::Log(maxM+binnerShift) - TMath::Log(minM+binnerShift)) / nMBins; // binwidth
// PT average -> see sumw2 function
const int nPtaveBins = 25; // number of bins
const double minPtave = 50.0;  // minimum
const double maxPtave = 2000.0; // maximum
double PtavelogBinWidth = (TMath::Log(maxPtave+binnerShift) - TMath::Log(minPtave+binnerShift)) / nPtaveBins; // binwidth
// PT leading/subleading -> see sumw2 function
const int nPtLSLBins = 50; // number of bins
const double minPtLSL = 20.0;  // minimum
const double maxPtLSL = 1000.0; // maximum
double PtLSLlogBinWidth = (TMath::Log(maxPtLSL+binnerShift) - TMath::Log(minPtLSL+binnerShift)) / nPtLSLBins; // binwidth
// Xj and Aj bins
const int nXjAjBins = 21; // number of bins
double XjBins[nXjAjBins+1] = {0.0,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5};
double AjBins[nXjAjBins+1] = {0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
// less bins
//const int nDphiBins = 18; // number of bins
//double DphiBins[nDphiBins+1] = {0.0, TMath::Pi()/5. ,TMath::Pi()/3., (3./7.)*TMath::Pi(), TMath::Pi()/2., (4./7.)*TMath::Pi(), (3./5.)*TMath::Pi(), 1.98967535,  2.0943951 ,  2.19911486,  2.30383461, 2.40855437,  2.51327412,  2.61799388,  2.72271363,  2.82743339, 2.93215314,  3.0368729, TMath::Pi()};
// More bins
const int nDphiBins = 30; // number of bins
double DphiBins[nDphiBins+1] = {0.0, TMath::Pi()/5. ,TMath::Pi()/3., (3./7.)*TMath::Pi(), TMath::Pi()/2., (4./7.)*TMath::Pi(), (3./5.)*TMath::Pi(), 1.93731547,  1.98967535,  2.04203522,  2.0943951 , 2.14675498,  2.19911486,  2.25147474,  2.30383461,  2.35619449, 2.40855437,  2.46091425,  2.51327412,  2.565634,  2.61799388, 2.67035376,  2.72271363,  2.77507351,  2.82743339,  2.87979327, 2.93215314,  2.98451302,  3.0368729 ,  3.08923278,  TMath::Pi()};

// for 3rd jet studies
const int nXjBins3rd = 17; // number of bins
double XjBins3rd[nXjBins3rd+1] = {-1.0,0.0,0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
const int nDphiBins3rd = 31; // number of bins
double DphiBins3rd[nDphiBins3rd+1] = {-1.0, 0.0, TMath::Pi()/5. ,TMath::Pi()/3., (3./7.)*TMath::Pi(), TMath::Pi()/2., (4./7.)*TMath::Pi(), (3./5.)*TMath::Pi(), 1.93731547,  1.98967535,  2.04203522,  2.0943951 , 2.14675498,  2.19911486,  2.25147474,  2.30383461,  2.35619449, 2.40855437,  2.46091425,  2.51327412,  2.565634,  2.61799388, 2.67035376,  2.72271363,  2.77507351,  2.82743339,  2.87979327, 2.93215314,  2.98451302,  3.0368729 ,  3.08923278,  TMath::Pi()};

// leading and subleading jet pTs
const int nPtLSLBins2 = 38; // number of bins
double PtLSLBins2[nPtLSLBins2+1] = {30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 350.0, 400.0, 500.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 8160.0};


// ============================ Event quantities ============================

// Number of events
TH1I *Nevents = new TH1I("Nevents", "Nevents", 10, 0, 10);
TH1I *Nev_recoreco = new TH1I("Nev_recoreco", "Nev_recoreco", 1, 0, 1);
TH1I *Nev_recoreco_lead = new TH1I("Nev_recoreco_lead", "Nev_recoreco_lead", 1, 0, 1);
TH1I *Nev_recoreco_subl = new TH1I("Nev_recoreco_subl", "Nev_recoreco_subl", 1, 0, 1);
TH1I *Nev_recogen = new TH1I("Nev_recogen", "Nev_recogen", 1, 0, 1);
TH1I *Nev_recogen_lead = new TH1I("Nev_recogen_lead", "Nev_recogen_lead", 1, 0, 1);
TH1I *Nev_recogen_subl = new TH1I("Nev_recogen_subl", "Nev_recogen_subl", 1, 0, 1);
TH1I *Nev_genreco = new TH1I("Nev_genreco", "Nev_genreco", 1, 0, 1);
TH1I *Nev_genreco_lead = new TH1I("Nev_genreco_lead", "Nev_genreco_lead", 1, 0, 1);
TH1I *Nev_genreco_subl = new TH1I("Nev_genreco_subl", "Nev_genreco_subl", 1, 0, 1);
TH1I *Nev_gengen = new TH1I("Nev_gengen", "Nev_gengen", 1, 0, 1);
TH1I *Nev_gengen_lead = new TH1I("Nev_gengen_lead", "Nev_gengen_lead", 1, 0, 1);
TH1I *Nev_gengen_subl = new TH1I("Nev_gengen_subl", "Nev_gengen_subl", 1, 0, 1);

// Multiplicities
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 80, 0.0, 400.0);
TH1D *multiplicity_weighted = new TH1D("multiplicity_weighted", "multiplicity_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet_weighted = new TH1D("multiplicity_withonejet_weighted", "multiplicity_withonejet_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_withdijets_weighted = new TH1D("multiplicity_withdijets_weighted", "multiplicity_withdijets_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_midmid_weighted = new TH1D("multiplicity_midmid_weighted", "multiplicity_midmid_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_midfwd_weighted = new TH1D("multiplicity_midfwd_weighted", "multiplicity_midfwd_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_midbkw_weighted = new TH1D("multiplicity_midbkw_weighted", "multiplicity_midbkw_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_fwdmid_weighted = new TH1D("multiplicity_fwdmid_weighted", "multiplicity_fwdmid_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_fwdfwd_weighted = new TH1D("multiplicity_fwdfwd_weighted", "multiplicity_fwdfwd_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_fwdbkw_weighted = new TH1D("multiplicity_fwdbkw_weighted", "multiplicity_fwdbkw_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_bkwmid_weighted = new TH1D("multiplicity_bkwmid_weighted", "multiplicity_bkwmid_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_bkwfwd_weighted = new TH1D("multiplicity_bkwfwd_weighted", "multiplicity_bkwfwd_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_bkwbkw_weighted = new TH1D("multiplicity_bkwbkw_weighted", "multiplicity_bkwbkw_weighted", 80, 0.0, 400.0);

TH1D *multiplicity_weighted_at1 = new TH1D("multiplicity_weighted_at1", "multiplicity_weighted_at1", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet_weighted_at1 = new TH1D("multiplicity_withonejet_weighted_at1", "multiplicity_withonejet_weighted_at1", 80, 0.0, 400.0);
TH1D *multiplicity_withdijets_weighted_at1 = new TH1D("multiplicity_withdijets_weighted_at1", "multiplicity_withdijets_weighted_at1", 80, 0.0, 400.0);

TH1D *reco_mult = new TH1D("reco_mult", "reco_mult", 80, 0.0, 400.0);
TH1D *reco_mult_weighted = new TH1D("reco_mult_weighted", "reco_mult_weighted", 80, 0.0, 400.0);
TH1D *reco_mult_withonejet_weighted = new TH1D("reco_mult_withonejet_weighted", "reco_mult_withonejet_weighted", 80, 0.0, 400.0);
TH1D *reco_mult_withdijets_weighted = new TH1D("reco_mult_withdijets_weighted", "reco_mult_withdijets_weighted", 80, 0.0, 400.0);
TH1D *gen_mult = new TH1D("gen_mult", "gen_mult", 80, 0.0, 400.0);
TH1D *gen_mult_weighted = new TH1D("gen_mult_weighted", "gen_mult_weighted", 80, 0.0, 400.0);
TH1D *gen_mult_withonejet_weighted = new TH1D("gen_mult_withonejet_weighted", "gen_mult_withonejet_weighted", 80, 0.0, 400.0);
TH1D *gen_mult_withdijets_weighted = new TH1D("gen_mult_withdijets_weighted", "gen_mult_withdijets_weighted", 80, 0.0, 400.0);

// Hadron Forward (HF) Calorimeter information
// Axis : 0 -> HF+, 1 -> HF-, 2 -> multbin
int	bins_HF[3]   =      { 200  ,  200 ,  80};
double xmin_HF[3]   =   { 0.0  ,  0.0 ,  0.0};
double xmax_HF[3]   =   { 400  ,  400 ,  400};
THnSparseD *hfhist = new THnSparseD("hfhist", "hfhist", 3, bins_HF, xmin_HF, xmax_HF);
THnSparseD *hfhist_weighted = new THnSparseD("hfhist_weighted", "hfhist_weighted", 3, bins_HF, xmin_HF, xmax_HF);
THnSparseD *hfhistEta4 = new THnSparseD("hfhistEta4", "hfhistEta4", 3, bins_HF, xmin_HF, xmax_HF);
THnSparseD *hfhistEta4_weighted = new THnSparseD("hfhistEta4_weighted", "hfhistEta4_weighted", 3, bins_HF, xmin_HF, xmax_HF);
THnSparseD *hfhist_onejet_weighted = new THnSparseD("hfhist_onejet_weighted", "hfhist_onejet_weighted", 3, bins_HF, xmin_HF, xmax_HF);
THnSparseD *hfhistEta4_onejet_weighted = new THnSparseD("hfhistEta4_onejet_weighted", "hfhistEta4_onejet_weighted", 3, bins_HF, xmin_HF, xmax_HF);
THnSparseD *hfhist_dijet_weighted = new THnSparseD("hfhist_dijet_weighted", "hfhist_dijet_weighted", 3, bins_HF, xmin_HF, xmax_HF);
THnSparseD *hfhistEta4_dijet_weighted = new THnSparseD("hfhistEta4_dijet_weighted", "hfhistEta4_dijet_weighted", 3, bins_HF, xmin_HF, xmax_HF);

// Zero Degree Calorimeter (ZDC) information
// Axis : 0 -> ZDC+, 1 -> ZDC-, 2 -> multbin
int	bins_ZDC[3]      =   {  200   ,   200   ,  80};
double xmin_ZDC[3]   =   { -10000 ,  -10000 ,  0.0};
double xmax_ZDC[3]   =   {  10000 ,   10000 ,  400};
THnSparseD *zdchist = new THnSparseD("zdchist", "zdchist", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);
THnSparseD *zdchist_weighted = new THnSparseD("zdchist_weighted", "zdchist_weighted", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);
THnSparseD *zdchist_onejet_weighted = new THnSparseD("zdchist_onejet_weighted", "zdchist_onejet_weighted", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);
THnSparseD *zdchist_dijet_weighted = new THnSparseD("zdchist_dijet_weighted", "zdchist_dijet_weighted", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);

// HFSum
int	bins_HFSum[2]   =      { 500  ,  80};
double xmin_HFSum[2]   =   { 0.0  ,  0.0};
double xmax_HFSum[2]   =   { 500  ,  400};
THnSparseD *hfhistSum_weighted = new THnSparseD("hfhistSum_weighted", "hfhistSum_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistEta4Sum_weighted = new THnSparseD("hfhistEta4Sum_weighted", "hfhistEta4Sum_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistSum_onejet_weighted = new THnSparseD("hfhistSum_onejet_weighted", "hfhistSum_onejet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistEta4Sum_onejet_weighted = new THnSparseD("hfhistEta4Sum_onejet_weighted", "hfhistEta4Sum_onejet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistSum_dijet_weighted = new THnSparseD("hfhistSum_dijet_weighted", "hfhistSum_dijet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistEta4Sum_dijet_weighted = new THnSparseD("hfhistEta4Sum_dijet_weighted", "hfhistEta4Sum_dijet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);

// Vertex
// Axis : 0 -> Vz, 1 -> event multiplicity, 2 -> extra dependency
int	bins_vz[3]   	=   {  60   ,   multbinsize-1    		 					  ,  extrabinsize-1};
double xmin_vz[3]   =   { -15.5 ,   multiplicity_centrality_bins[0]  			  ,  extra_bins[0]};
double xmax_vz[3]   =   {  15.5 ,   multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *vzhist = new THnSparseD("vzhist", "vzhist", 3, bins_vz, xmin_vz, xmax_vz);
THnSparseD *vzhist_weighted = new THnSparseD("vzhist_weighted", "vzhist_weighted", 3, bins_vz, xmin_vz, xmax_vz);
THnSparseD *vzhist_jet_weighted = new THnSparseD("vzhist_jet_weighted", "vzhist_jet_weighted", 3, bins_vz, xmin_vz, xmax_vz);
THnSparseD *vzhist_dijet_weighted = new THnSparseD("vzhist_dijet_weighted", "vzhist_dijet_weighted", 3, bins_vz, xmin_vz, xmax_vz);
// Axis : 0 -> Vx, 1 -> Vy, 2 -> event multiplicity, 3 -> extra dependency
int	bins_vxy[4]   	 =   {  100 , 100  , multbinsize-1    		 					  ,  extrabinsize-1};
double xmin_vxy[4]   =   { -1.0 , -1.0 , multiplicity_centrality_bins[0]  			  ,  extra_bins[0]};
double xmax_vxy[4]   =   {  1.0 ,  1.0 , multiplicity_centrality_bins[multbinsize-1]  ,  extra_bins[extrabinsize-1]};
THnSparseD *vxyhist = new THnSparseD("vxyhist", "vxyhist", 4, bins_vxy, xmin_vxy, xmax_vxy);
THnSparseD *vxyhist_weighted = new THnSparseD("vxyhist_weighted", "vxyhist_weighted", 4, bins_vxy, xmin_vxy, xmax_vxy);

// Pthat
// Axis : 0 -> Pthat, 1 -> event multiplicity, 2 -> extra dependency
int	bins_pthat[3]      =   {  100 	 ,   multbinsize-1    		  					   ,  extrabinsize-1};
double xmin_pthat[3]   =   {  0.0 	 ,   multiplicity_centrality_bins[0]  			   ,  extra_bins[0]};
double xmax_pthat[3]   =   {  1000.0 ,   multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *pthathist = new THnSparseD("pthathist", "pthathist", 3, bins_pthat, xmin_pthat, xmax_pthat);
THnSparseD *pthathist_weighted = new THnSparseD("pthathist_weighted", "pthathist_weighted", 3, bins_pthat, xmin_pthat, xmax_pthat);

// event plane histograms
// Axis : 0 -> EP multiplicity, 1 -> qvector, 2 -> PsiEP, 3 -> event multiplicity, 4 -> extra dependency
int	bins_EP[5]   	=   { 40	,  200 ,   32		     , multbinsize-1								,  extrabinsize-1};
double xmin_EP[5]   =   { 0	 	,  0   ,   -TMath::Pi()  , multiplicity_centrality_bins[0]  			,  extra_bins[0]};
double xmax_EP[5]   =   { 400   ,  100 ,   TMath::Pi()   , multiplicity_centrality_bins[multbinsize-1]	,  extra_bins[extrabinsize-1]};
THnSparseD *EP2_plus_flat = new THnSparseD("EP2_plus_flat", "EP2_plus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP2_minus_flat = new THnSparseD("EP2_minus_flat", "EP2_minus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP3_plus_flat = new THnSparseD("EP3_plus_flat", "EP3_plus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP3_minus_flat = new THnSparseD("EP3_minus_flat", "EP3_minus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP4_plus_flat = new THnSparseD("EP4_plus_flat", "EP4_plus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP4_minus_flat = new THnSparseD("EP4_minus_flat", "EP4_minus_flat", 5, bins_EP, xmin_EP, xmax_EP);

// Track/Particle histograms
// Axis : 0 -> track pT, 1 -> trk eta, 2 -> trk phi, 3 -> multiplicity bin, 4 -> extra dependency
int	bins_trk[5]      =   { 100   ,  30  ,   32		     , multbinsize-1								,  extrabinsize-1};
double xmin_trk[5]   =   { 0.0   , -3.0 ,   -TMath::Pi() , multiplicity_centrality_bins[0]  			,  extra_bins[0]};
double xmax_trk[5]   =   { 50.0  ,  3.0 ,   TMath::Pi()  , multiplicity_centrality_bins[multbinsize-1]	,  extra_bins[extrabinsize-1]};
// --> Reco
THnSparseD *hist_reco_trk = new THnSparseD("hist_reco_trk", "hist_reco_trk", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_reco_trk_corr = new THnSparseD("hist_reco_trk_corr", "hist_reco_trk_corr", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_reco_trk_weighted = new THnSparseD("hist_reco_trk_weighted", "hist_reco_trk_weighted", 5, bins_trk, xmin_trk, xmax_trk);
// --> Gen
THnSparseD *hist_gen_trk = new THnSparseD("hist_gen_trk", "hist_gen_trk", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_gen_trk_weighted = new THnSparseD("hist_gen_trk_weighted", "hist_gen_trk_weighted", 5, bins_trk, xmin_trk, xmax_trk);

// Tracks from Correlations
// Inclusive
THnSparseD *hist_trk_from_reco_reco_sig = new THnSparseD("hist_trk_from_reco_reco_sig", "hist_trk_from_reco_reco_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_trk_from_reco_gen_sig = new THnSparseD("hist_trk_from_reco_gen_sig", "hist_trk_from_reco_gen_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_trk_from_gen_reco_sig = new THnSparseD("hist_trk_from_gen_reco_sig", "hist_trk_from_gen_reco_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_trk_from_gen_gen_sig = new THnSparseD("hist_trk_from_gen_gen_sig", "hist_trk_from_gen_gen_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_trk_from_reco_reco_mix = new THnSparseD("hist_trk_from_reco_reco_mix", "hist_trk_from_reco_reco_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_trk_from_reco_gen_mix = new THnSparseD("hist_trk_from_reco_gen_mix", "hist_trk_from_reco_gen_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_trk_from_gen_reco_mix = new THnSparseD("hist_trk_from_gen_reco_mix", "hist_trk_from_gen_reco_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_trk_from_gen_gen_mix = new THnSparseD("hist_trk_from_gen_gen_mix", "hist_trk_from_gen_gen_mix", 5, bins_trk, xmin_trk, xmax_trk);
// Leading
THnSparseD *hist_LJ_trk_from_reco_reco_sig = new THnSparseD("hist_LJ_trk_from_reco_reco_sig", "hist_LJ_trk_from_reco_reco_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_LJ_trk_from_reco_gen_sig = new THnSparseD("hist_LJ_trk_from_reco_gen_sig", "hist_LJ_trk_from_reco_gen_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_LJ_trk_from_gen_reco_sig = new THnSparseD("hist_LJ_trk_from_gen_reco_sig", "hist_LJ_trk_from_gen_reco_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_LJ_trk_from_gen_gen_sig = new THnSparseD("hist_LJ_trk_from_gen_gen_sig", "hist_LJ_trk_from_gen_gen_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_LJ_trk_from_reco_reco_mix = new THnSparseD("hist_LJ_trk_from_reco_reco_mix", "hist_LJ_trk_from_reco_reco_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_LJ_trk_from_reco_gen_mix = new THnSparseD("hist_LJ_trk_from_reco_gen_mix", "hist_LJ_trk_from_reco_gen_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_LJ_trk_from_gen_reco_mix = new THnSparseD("hist_LJ_trk_from_gen_reco_mix", "hist_LJ_trk_from_gen_reco_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_LJ_trk_from_gen_gen_mix = new THnSparseD("hist_LJ_trk_from_gen_gen_mix", "hist_LJ_trk_from_gen_gen_mix", 5, bins_trk, xmin_trk, xmax_trk);
// Subleading
THnSparseD *hist_SLJ_trk_from_reco_reco_sig = new THnSparseD("hist_SLJ_trk_from_reco_reco_sig", "hist_SLJ_trk_from_reco_reco_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_SLJ_trk_from_reco_gen_sig = new THnSparseD("hist_SLJ_trk_from_reco_gen_sig", "hist_SLJ_trk_from_reco_gen_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_SLJ_trk_from_gen_reco_sig = new THnSparseD("hist_SLJ_trk_from_gen_reco_sig", "hist_SLJ_trk_from_gen_reco_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_SLJ_trk_from_gen_gen_sig = new THnSparseD("hist_SLJ_trk_from_gen_gen_sig", "hist_SLJ_trk_from_gen_gen_sig", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_SLJ_trk_from_reco_reco_mix = new THnSparseD("hist_SLJ_trk_from_reco_reco_mix", "hist_SLJ_trk_from_reco_reco_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_SLJ_trk_from_reco_gen_mix = new THnSparseD("hist_SLJ_trk_from_reco_gen_mix", "hist_SLJ_trk_from_reco_gen_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_SLJ_trk_from_gen_reco_mix = new THnSparseD("hist_SLJ_trk_from_gen_reco_mix", "hist_SLJ_trk_from_gen_reco_mix", 5, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_SLJ_trk_from_gen_gen_mix = new THnSparseD("hist_SLJ_trk_from_gen_gen_mix", "hist_SLJ_trk_from_gen_gen_mix", 5, bins_trk, xmin_trk, xmax_trk);

// Tracks for event plane
// Axis : 0 -> delta phi between track and EP, 1 -> trkbin, 2 -> multbin, 3 -> extra variable
int	bins_TRKEP[4]      =   { 40				   		,  trkbinsize-1 			 ,  multbinsize-1								 ,  extrabinsize-1};
double xmin_TRKEP[4]   =   { -TMath::Pi()/2.0	    ,  trk_pt_bins[0] 			 ,  multiplicity_centrality_bins[0]  			 ,  extra_bins[0]};
double xmax_TRKEP[4]   =   { 3.0*TMath::Pi()/2.0 	,  trk_pt_bins[trkbinsize-1] ,  multiplicity_centrality_bins[multbinsize-1]  ,  extra_bins[extrabinsize-1]};
THnSparseD *Dphi_EP2_flat_trk_minus = new THnSparseD("Dphi_EP2_flat_trk_minus", "Dphi_EP2_flat_trk_minus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_EP2_flat_trk_plus = new THnSparseD("Dphi_EP2_flat_trk_plus", "Dphi_EP2_flat_trk_plus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_EP3_flat_trk_minus = new THnSparseD("Dphi_EP3_flat_trk_minus", "Dphi_EP3_flat_trk_minus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_EP3_flat_trk_plus = new THnSparseD("Dphi_EP3_flat_trk_plus", "Dphi_EP3_flat_trk_plus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_EP4_flat_trk_minus = new THnSparseD("Dphi_EP4_flat_trk_minus", "Dphi_EP4_flat_trk_minus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_EP4_flat_trk_plus = new THnSparseD("Dphi_EP4_flat_trk_plus", "Dphi_EP4_flat_trk_plus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_GEN_EP2_flat_trk_minus = new THnSparseD("Dphi_GEN_EP2_flat_trk_minus", "Dphi_GEN_EP2_flat_trk_minus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_GEN_EP2_flat_trk_plus = new THnSparseD("Dphi_GEN_EP2_flat_trk_plus", "Dphi_GEN_EP2_flat_trk_plus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_GEN_EP3_flat_trk_minus = new THnSparseD("Dphi_GEN_EP3_flat_trk_minus", "Dphi_GEN_EP3_flat_trk_minus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_GEN_EP3_flat_trk_plus = new THnSparseD("Dphi_GEN_EP3_flat_trk_plus", "Dphi_GEN_EP3_flat_trk_plus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_GEN_EP4_flat_trk_minus = new THnSparseD("Dphi_GEN_EP4_flat_trk_minus", "Dphi_GEN_EP4_flat_trk_minus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);
THnSparseD *Dphi_GEN_EP4_flat_trk_plus = new THnSparseD("Dphi_GEN_EP4_flat_trk_plus", "Dphi_GEN_EP4_flat_trk_plus", 4, bins_TRKEP, xmin_TRKEP, xmax_TRKEP);

// until here it is fine
// Jet histograms
// Number of jets per event
int	bins_NJETS[3]      =   {  30   ,   multbinsize-1    		 					 ,  extrabinsize-1};
double xmin_NJETS[3]   =   {  0.0  ,   multiplicity_centrality_bins[0]  			 ,  extra_bins[0]};
double xmax_NJETS[3]   =   {  30.0 ,   multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *NJets = new THnSparseD("NJets", "NJets", 3, bins_NJETS, xmin_NJETS, xmax_NJETS);

// trackmax histogram
// Axis : 0 -> max track pt in a jet, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_trkmax[3]      =   {  200   ,   multbinsize-1    							    ,  extrabinsize-1};
double xmin_trkmax[3]   =   {  0.0 	  ,   multiplicity_centrality_bins[0]  			    ,  extra_bins[0]};
double xmax_trkmax[3]   =   {  100.0 ,   multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *trackmaxptinjethisto = new THnSparseD("trackmaxptinjethisto", "trackmaxptinjethisto", 3, bins_trkmax, xmin_trkmax, xmax_trkmax);

// trackmax/rawjet histogram
// Axis : 0 -> max track pt in a jet, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_trkmaxjet[3]      =   {  100  ,   multbinsize-1    							    ,  extrabinsize-1};
double xmin_trkmaxjet[3]   =   {  0.0   ,   multiplicity_centrality_bins[0]  			    ,  extra_bins[0]};
double xmax_trkmaxjet[3]   =   {  1.0   ,   multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *jettrackmaxptinjethisto = new THnSparseD("jettrackmaxptinjethisto", "jettrackmaxptinjethisto", 3, bins_trkmaxjet, xmin_trkmaxjet, xmax_trkmaxjet);

// UE histogram
// Axis : 0 -> UE, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_UE[3]      =   {  1000  	,   multbinsize-1    							    ,  extrabinsize-1};
double xmin_UE[3]   =   {  -1000.0   	,   multiplicity_centrality_bins[0]  			    ,  extra_bins[0]};
double xmax_UE[3]   =   {  1000.0   ,   multiplicity_centrality_bins[multbinsize-1]   	,  extra_bins[extrabinsize-1]};
THnSparseD *histo_jetUE 			 = new THnSparseD("histo_jetUE", "histo_jetUE", 3, bins_UE, xmin_UE, xmax_UE);
THnSparseD *histo_jetAverageRho 	 = new THnSparseD("histo_jetAverageRho", "histo_jetAverageRho", 3, bins_UE, xmin_UE, xmax_UE);
THnSparseD *histo_jetcheckcorrection = new THnSparseD("histo_jetcheckcorrection", "histo_jetcheckcorrection", 3, bins_UE, xmin_UE, xmax_UE);

//correlations to EP
// Axis : 0 -> delta phi between jet and EP, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_JETEP[3]      =   { 200			   		,  multbinsize-1							 	 ,  extrabinsize-1};
double xmin_JETEP[3]   =   { -TMath::Pi()/2.0	    ,  multiplicity_centrality_bins[0]  			 ,  extra_bins[0]};
double xmax_JETEP[3]   =   { 3.0*TMath::Pi()/2.0 	,  multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};

THnSparseD *Dphi_flat_EP2_inclusive_minus = new THnSparseD("Dphi_flat_EP2_inclusive_minus", "Dphi_flat_EP2_inclusive_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP2_leading_minus = new THnSparseD("Dphi_flat_EP2_leading_minus", "Dphi_flat_EP2_leading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP2_subleading_minus = new THnSparseD("Dphi_flat_EP2_subleading_minus", "Dphi_flat_EP2_subleading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP2_inclusive_plus = new THnSparseD("Dphi_flat_EP2_inclusive_plus", "Dphi_flat_EP2_inclusive_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP2_leading_plus = new THnSparseD("Dphi_flat_EP2_leading_plus", "Dphi_flat_EP2_leading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP2_subleading_plus = new THnSparseD("Dphi_flat_EP2_subleading_plus", "Dphi_flat_EP2_subleading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP3_inclusive_minus = new THnSparseD("Dphi_flat_EP3_inclusive_minus", "Dphi_flat_EP3_inclusive_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP3_leading_minus = new THnSparseD("Dphi_flat_EP3_leading_minus", "Dphi_flat_EP3_leading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP3_subleading_minus = new THnSparseD("Dphi_flat_EP3_subleading_minus", "Dphi_flat_EP3_subleading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP3_inclusive_plus = new THnSparseD("Dphi_flat_EP3_inclusive_plus", "Dphi_flat_EP3_inclusive_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP3_leading_plus = new THnSparseD("Dphi_flat_EP3_leading_plus", "Dphi_flat_EP3_leading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP3_subleading_plus = new THnSparseD("Dphi_flat_EP3_subleading_plus", "Dphi_flat_EP3_subleading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP4_inclusive_minus = new THnSparseD("Dphi_flat_EP4_inclusive_minus", "Dphi_flat_EP4_inclusive_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP4_leading_minus = new THnSparseD("Dphi_flat_EP4_leading_minus", "Dphi_flat_EP4_leading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP4_subleading_minus = new THnSparseD("Dphi_flat_EP4_subleading_minus", "Dphi_flat_EP4_subleading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP4_inclusive_plus = new THnSparseD("Dphi_flat_EP4_inclusive_plus", "Dphi_flat_EP4_inclusive_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP4_leading_plus = new THnSparseD("Dphi_flat_EP4_leading_plus", "Dphi_flat_EP4_leading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_flat_EP4_subleading_plus = new THnSparseD("Dphi_flat_EP4_subleading_plus", "Dphi_flat_EP4_subleading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);

THnSparseD *Dphi_GEN_flat_EP2_inclusive_minus = new THnSparseD("Dphi_GEN_flat_EP2_inclusive_minus", "Dphi_GEN_flat_EP2_inclusive_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP2_leading_minus = new THnSparseD("Dphi_GEN_flat_EP2_leading_minus", "Dphi_GEN_flat_EP2_leading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP2_subleading_minus = new THnSparseD("Dphi_GEN_flat_EP2_subleading_minus", "Dphi_GEN_flat_EP2_subleading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP2_inclusive_plus = new THnSparseD("Dphi_GEN_flat_EP2_inclusive_plus", "Dphi_GEN_flat_EP2_inclusive_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP2_leading_plus = new THnSparseD("Dphi_GEN_flat_EP2_leading_plus", "Dphi_GEN_flat_EP2_leading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP2_subleading_plus = new THnSparseD("Dphi_GEN_flat_EP2_subleading_plus", "Dphi_GEN_flat_EP2_subleading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP3_inclusive_minus = new THnSparseD("Dphi_GEN_flat_EP3_inclusive_minus", "Dphi_GEN_flat_EP3_inclusive_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP3_leading_minus = new THnSparseD("Dphi_GEN_flat_EP3_leading_minus", "Dphi_GEN_flat_EP3_leading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP3_subleading_minus = new THnSparseD("Dphi_GEN_flat_EP3_subleading_minus", "Dphi_GEN_flat_EP3_subleading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP3_inclusive_plus = new THnSparseD("Dphi_GEN_flat_EP3_inclusive_plus", "Dphi_GEN_flat_EP3_inclusive_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP3_leading_plus = new THnSparseD("Dphi_GEN_flat_EP3_leading_plus", "Dphi_GEN_flat_EP3_leading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP3_subleading_plus = new THnSparseD("Dphi_GEN_flat_EP3_subleading_plus", "Dphi_GEN_flat_EP3_subleading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP4_inclusive_minus = new THnSparseD("Dphi_GEN_flat_EP4_inclusive_minus", "Dphi_GEN_flat_EP4_inclusive_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP4_leading_minus = new THnSparseD("Dphi_GEN_flat_EP4_leading_minus", "Dphi_GEN_flat_EP4_leading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP4_subleading_minus = new THnSparseD("Dphi_GEN_flat_EP4_subleading_minus", "Dphi_GEN_flat_EP4_subleading_minus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP4_inclusive_plus = new THnSparseD("Dphi_GEN_flat_EP4_inclusive_plus", "Dphi_GEN_flat_EP4_inclusive_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP4_leading_plus = new THnSparseD("Dphi_GEN_flat_EP4_leading_plus", "Dphi_GEN_flat_EP4_leading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);
THnSparseD *Dphi_GEN_flat_EP4_subleading_plus = new THnSparseD("Dphi_GEN_flat_EP4_subleading_plus", "Dphi_GEN_flat_EP4_subleading_plus", 3, bins_JETEP, xmin_JETEP, xmax_JETEP);

// Jet quantities
int	bins_jet[5]      =   { 100	  ,  40  ,   32		      , multbinsize-1          						, extrabinsize-1};
double xmin_jet[5]   =   { 0.0	  , -4.0 ,   -TMath::Pi() , multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jet[5]   =   { 1000.0  ,  4.0 ,   TMath::Pi() , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
// --> Reco
THnSparseD *hist_reco_jet = new THnSparseD("hist_reco_jet", "hist_reco_jet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_jet_corr = new THnSparseD("hist_reco_jet_corr", "hist_reco_jet_corr", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_jet_weighted = new THnSparseD("hist_reco_jet_weighted", "hist_reco_jet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_jet_corr_weighted = new THnSparseD("hist_reco_jet_corr_weighted", "hist_reco_jet_corr_weighted", 5, bins_jet, xmin_jet, xmax_jet);
// --> Gen
THnSparseD *hist_gen_jet = new THnSparseD("hist_gen_jet", "hist_gen_jet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_jet_weighted = new THnSparseD("hist_gen_jet_weighted", "hist_gen_jet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
// Leading and Subleading Jets
// --> Reco
THnSparseD *hist_reco_leadjet = new THnSparseD("hist_reco_leadjet", "hist_reco_leadjet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_leadjet_weighted = new THnSparseD("hist_reco_leadjet_weighted", "hist_reco_leadjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_subljet = new THnSparseD("hist_reco_subljet", "hist_reco_subljet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_subljet_weighted = new THnSparseD("hist_reco_subljet_weighted", "hist_reco_subljet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_thrdjet_weighted = new THnSparseD("hist_reco_thrdjet_weighted", "hist_reco_thrdjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
// --> Gen
THnSparseD *hist_gen_leadjet = new THnSparseD("hist_gen_leadjet", "hist_gen_leadjet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_leadjet_weighted = new THnSparseD("hist_gen_leadjet_weighted", "hist_gen_leadjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_subljet = new THnSparseD("hist_gen_subljet", "hist_gen_subljet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_subljet_weighted = new THnSparseD("hist_gen_subljet_weighted", "hist_gen_subljet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_thrdjet_weighted = new THnSparseD("hist_gen_thrdjet_weighted", "hist_gen_thrdjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
// --> Ref
THnSparseD *hist_ref_jet_weighted = new THnSparseD("hist_ref_jet_weighted", "hist_ref_jet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_leadjet_weighted = new THnSparseD("hist_ref_leadjet_weighted", "hist_ref_leadjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_subljet_weighted = new THnSparseD("hist_ref_subljet_weighted", "hist_ref_subljet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_thrdjet_weighted = new THnSparseD("hist_ref_thrdjet_weighted", "hist_ref_thrdjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);

// Jets from Correlations
// Inclusive
THnSparseD *hist_jet_from_reco_reco_sig = new THnSparseD("hist_jet_from_reco_reco_sig", "hist_jet_from_reco_reco_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_jet_from_reco_gen_sig = new THnSparseD("hist_jet_from_reco_gen_sig", "hist_jet_from_reco_gen_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_jet_from_gen_reco_sig = new THnSparseD("hist_jet_from_gen_reco_sig", "hist_jet_from_gen_reco_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_jet_from_gen_gen_sig = new THnSparseD("hist_jet_from_gen_gen_sig", "hist_jet_from_gen_gen_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_jet_from_reco_reco_mix = new THnSparseD("hist_jet_from_reco_reco_mix", "hist_jet_from_reco_reco_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_jet_from_reco_gen_mix = new THnSparseD("hist_jet_from_reco_gen_mix", "hist_jet_from_reco_gen_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_jet_from_gen_reco_mix = new THnSparseD("hist_jet_from_gen_reco_mix", "hist_jet_from_gen_reco_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_jet_from_gen_gen_mix = new THnSparseD("hist_jet_from_gen_gen_mix", "hist_jet_from_gen_gen_mix", 5, bins_jet, xmin_jet, xmax_jet);
// Leading
THnSparseD *hist_lead_jet_from_reco_reco_sig = new THnSparseD("hist_lead_jet_from_reco_reco_sig", "hist_lead_jet_from_reco_reco_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_lead_jet_from_reco_gen_sig = new THnSparseD("hist_lead_jet_from_reco_gen_sig", "hist_lead_jet_from_reco_gen_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_lead_jet_from_gen_reco_sig = new THnSparseD("hist_lead_jet_from_gen_reco_sig", "hist_lead_jet_from_gen_reco_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_lead_jet_from_gen_gen_sig = new THnSparseD("hist_lead_jet_from_gen_gen_sig", "hist_lead_jet_from_gen_gen_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_lead_jet_from_reco_reco_mix = new THnSparseD("hist_lead_jet_from_reco_reco_mix", "hist_lead_jet_from_reco_reco_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_lead_jet_from_reco_gen_mix = new THnSparseD("hist_lead_jet_from_reco_gen_mix", "hist_lead_jet_from_reco_gen_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_lead_jet_from_gen_reco_mix = new THnSparseD("hist_lead_jet_from_gen_reco_mix", "hist_lead_jet_from_gen_reco_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_lead_jet_from_gen_gen_mix = new THnSparseD("hist_lead_jet_from_gen_gen_mix", "hist_lead_jet_from_gen_gen_mix", 5, bins_jet, xmin_jet, xmax_jet);
// Subleading
THnSparseD *hist_subl_jet_from_reco_reco_sig = new THnSparseD("hist_subl_jet_from_reco_reco_sig", "hist_subl_jet_from_reco_reco_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_subl_jet_from_reco_gen_sig = new THnSparseD("hist_subl_jet_from_reco_gen_sig", "hist_subl_jet_from_reco_gen_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_subl_jet_from_gen_reco_sig = new THnSparseD("hist_subl_jet_from_gen_reco_sig", "hist_subl_jet_from_gen_reco_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_subl_jet_from_gen_gen_sig = new THnSparseD("hist_subl_jet_from_gen_gen_sig", "hist_subl_jet_from_gen_gen_sig", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_subl_jet_from_reco_reco_mix = new THnSparseD("hist_subl_jet_from_reco_reco_mix", "hist_subl_jet_from_reco_reco_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_subl_jet_from_reco_gen_mix = new THnSparseD("hist_subl_jet_from_reco_gen_mix", "hist_subl_jet_from_reco_gen_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_subl_jet_from_gen_reco_mix = new THnSparseD("hist_subl_jet_from_gen_reco_mix", "hist_subl_jet_from_gen_reco_mix", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_subl_jet_from_gen_gen_mix = new THnSparseD("hist_subl_jet_from_gen_gen_mix", "hist_subl_jet_from_gen_gen_mix", 5, bins_jet, xmin_jet, xmax_jet);

// Jet Energy Scale (JES) and Jet Energy Resolution (JER)
int	bins_jes[6]   =      { 200  ,  100  ,  40  ,  8, multbinsize-1          					 , extrabinsize-1};
double xmin_jes[6]   =   { 0.0  ,  0    , -4.0 ,  0, multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jes[6]   =   { 5.0 ,  1000 ,  4.0 ,  8, multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_jes_reco_weighted = new THnSparseD("hist_jes_reco_weighted", "hist_jes_reco_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_jes_reco_fromB_weighted = new THnSparseD("hist_jes_reco_fromB_weighted", "hist_jes_reco_fromB_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_leadjes_reco_weighted = new THnSparseD("hist_leadjes_reco_weighted", "hist_leadjes_reco_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_leadjes_reco_fromB_weighted = new THnSparseD("hist_leadjes_reco_fromB_weighted", "hist_leadjes_reco_fromB_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_subleadjes_reco_weighted = new THnSparseD("hist_subleadjes_reco_weighted", "hist_subleadjes_reco_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_subleadjes_reco_fromB_weighted = new THnSparseD("hist_subleadjes_reco_fromB_weighted", "hist_subleadjes_reco_fromB_weighted", 6, bins_jes, xmin_jes, xmax_jes);

// for jet closure
int    bins_jetptclos[6]   =   { 100  ,  100  , 40  , 40  , multbinsize-1                               , extrabinsize-1};
double xmin_jetptclos[6]   =   { 0    ,  0    , -4.0, -4.0, multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jetptclos[6]   =   { 1000 ,  1000 , 4.0 , 4.0 , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_jetptclos_weighted = new THnSparseD("hist_jetptclos_weighted", "hist_jetptclos_weighted", 6, bins_jetptclos, xmin_jetptclos, xmax_jetptclos);
THnSparseD *hist_leadjetptclos_weighted = new THnSparseD("hist_leadjetptclos_weighted", "hist_leadjetptclos_weighted", 6, bins_jetptclos, xmin_jetptclos, xmax_jetptclos);
THnSparseD *hist_subljetptclos_weighted = new THnSparseD("hist_subljetptclos_weighted", "hist_subljetptclos_weighted", 6, bins_jetptclos, xmin_jetptclos, xmax_jetptclos);
THnSparseD *hist_leadjetptclosremovesome_weighed = new THnSparseD("hist_leadjetptclosremovesome_weighed", "hist_leadjetptclosremovesome_weighed", 6, bins_jetptclos, xmin_jetptclos, xmax_jetptclos);
THnSparseD *hist_subljetptclosremovesome_weighed = new THnSparseD("hist_subljetptclosremovesome_weighed", "hist_subljetptclosremovesome_weighed", 6, bins_jetptclos, xmin_jetptclos, xmax_jetptclos);

int    bins_jetptclosave[8]   =   { 50  ,  50  , 80  , 80  , 80  , 80  , multbinsize-1                               , extrabinsize-1};
double xmin_jetptclosave[8]   =   { 0    ,  0   , -4.0, -4.0, -4.0, -4.0, multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jetptclosave[8]   =   { 1000 ,  1000 , 4.0 , 4.0 , 4.0 , 4.0 , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_averjetptclos_weighted = new THnSparseD("hist_averjetptclos_weighted", "hist_averjetptclos_weighted", 8, bins_jetptclosave, xmin_jetptclosave, xmax_jetptclosave);
THnSparseD *hist_averjetptclosremovesome_weighed = new THnSparseD("hist_averjetptclosremovesome_weighed", "hist_averjetptclosremovesome_weighed", 8, bins_jetptclosave, xmin_jetptclosave, xmax_jetptclosave);

int    bins_xjclos[8]   =   { 10  ,  10  , 80  , 80  , 80  , 80  , multbinsize-1                               , extrabinsize-1};
double xmin_xjclos[8]   =   { 0   ,  0   , -4.0, -4.0, -4.0, -4.0, multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_xjclos[8]   =   { 1.0 ,  1.0 , 4.0 , 4.0 , 4.0 , 4.0 , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_xjclos_weighted = new THnSparseD("hist_xjclos_weighted", "hist_xjclos_weighted", 8, bins_xjclos, xmin_xjclos, xmax_xjclos);
THnSparseD *hist_xjclos_removesome_weighted = new THnSparseD("hist_xjclos_removesome_weighted", "hist_xjclos_removesome_weighted", 8, bins_xjclos, xmin_xjclos, xmax_xjclos);

// for jet unfolding
int    bins_jetunf[4]   =   { nPtLSLBins2    , nPtLSLBins2   , multbinsize-1                               , extrabinsize-1};
double xmin_jetunf[4]   =   { PtLSLBins2[0]  , PtLSLBins2[0] , multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jetunf[4]   =   { 8160.0		 , 8160.0 		 , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_jetunf_weighted = new THnSparseD("hist_jetunf_weighted", "hist_jetunf_weighted", 4, bins_jetunf, xmin_jetunf, xmax_jetunf);
THnSparseD *hist_leadjetunf_weighted = new THnSparseD("hist_leadjetunf_weighted", "hist_leadjetunf_weighted", 4, bins_jetunf, xmin_jetunf, xmax_jetunf);
THnSparseD *hist_subljetunf_weighted = new THnSparseD("hist_subljetunf_weighted", "hist_subljetunf_weighted", 4, bins_jetunf, xmin_jetunf, xmax_jetunf);
THnSparseD *hist_leadjetunf_match_weighted = new THnSparseD("hist_leadjetunf_match_weighted", "hist_leadjetunf_match_weighted", 4, bins_jetunf, xmin_jetunf, xmax_jetunf);
THnSparseD *hist_subljetunf_match_weighted = new THnSparseD("hist_subljetunf_match_weighted", "hist_subljetunf_match_weighted", 4, bins_jetunf, xmin_jetunf, xmax_jetunf);
THnSparseD *hist_leadjetunf_swap_weighted = new THnSparseD("hist_leadjetunf_swap_weighted", "hist_leadjetunf_swap_weighted", 4, bins_jetunf, xmin_jetunf, xmax_jetunf);
THnSparseD *hist_subljetunf_swap_weighted = new THnSparseD("hist_subljetunf_swap_weighted", "hist_subljetunf_swap_weighted", 4, bins_jetunf, xmin_jetunf, xmax_jetunf);

int    bins_jetunf4D[6]   =   { nPtLSLBins2  	,  nPtLSLBins2   , nPtLSLBins2  	,  nPtLSLBins2   , multbinsize-1                               , extrabinsize-1};
double xmin_jetunf4D[6]   =   { PtLSLBins2[0] ,  PtLSLBins2[0] , PtLSLBins2[0] 	,  PtLSLBins2[0] , multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jetunf4D[6]   =   { 8160.0 		,  8160.0 		 , 8160.0 			,  8160.0 		 , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_jetunf_weighted_4D = new THnSparseD("hist_jetunf_weighted_4D", "hist_jetunf_weighted_4D", 6, bins_jetunf4D, xmin_jetunf4D, xmax_jetunf4D);
THnSparseD *hist_jetunf_match_weighted_4D = new THnSparseD("hist_jetunf_match_weighted_4D", "hist_jetunf_match_weighted_4D", 6, bins_jetunf4D, xmin_jetunf4D, xmax_jetunf4D);
THnSparseD *hist_jetunf_swap_weighted_4D = new THnSparseD("hist_jetunf_swap_weighted_4D", "hist_jetunf_swap_weighted_4D", 6, bins_jetunf4D, xmin_jetunf4D, xmax_jetunf4D);
THnSparseD *hist_jetunf_match_weighted_sym_4D = new THnSparseD("hist_jetunf_match_weighted_sym_4D", "hist_jetunf_match_weighted_sym_4D", 6, bins_jetunf4D, xmin_jetunf4D, xmax_jetunf4D);
THnSparseD *hist_jetunf_match_weighted_symhalf_4D = new THnSparseD("hist_jetunf_match_weighted_symhalf_4D", "hist_jetunf_match_weighted_symhalf_4D", 6, bins_jetunf4D, xmin_jetunf4D, xmax_jetunf4D);

int    bins_xjunf[4]   =   { nXjAjBins  ,  nXjAjBins  , multbinsize-1                               , extrabinsize-1};
double xmin_xjunf[4]   =   { 0    		,  0    	  , multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_xjunf[4]   =   { 1.5 		,  1.5 		  , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_xjunf_weighted = new THnSparseD("hist_xjunf_weighted", "hist_xjunf_weighted", 4, bins_xjunf, xmin_xjunf, xmax_xjunf);
THnSparseD *hist_xjunf_match_weighted = new THnSparseD("hist_xjunf_match_weighted", "hist_xjunf_match_weighted", 4, bins_xjunf, xmin_xjunf, xmax_xjunf);
THnSparseD *hist_xjunf_swap_weighted = new THnSparseD("hist_xjunf_swap_weighted", "hist_xjunf_swap_weighted", 4, bins_xjunf, xmin_xjunf, xmax_xjunf);

int    bins_xjunfptave[6]   =   { nXjAjBins  ,  nXjAjBins , nPtLSLBins2  , nPtLSLBins2  , multbinsize-1                               , extrabinsize-1};
double xmin_xjunfptave[6]   =   { 0    		,  0    	  , PtLSLBins2[0], PtLSLBins2[0], multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_xjunfptave[6]   =   { 1.5 		,  1.5 		  , 8160.0       , 8160.0       ,  multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_xjunfptave_weighted = new THnSparseD("hist_xjunfptave_weighted", "hist_xjunfptave_weighted", 6, bins_xjunfptave, xmin_xjunfptave, xmax_xjunfptave);


// cross-check it
int    bins_jetunf_smear[3]       =   { nPtLSLBins2     , multbinsize-1                               , extrabinsize-1};
double xmin_jetunf_smear[3]   	  =   { PtLSLBins2[0]   , multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jetunf_smear[3]   	  =   { 8160.0 			, multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_leadjetunf_gensmear = new THnSparseD("hist_leadjetunf_gensmear", "hist_leadjetunf_gensmear", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_leadjetunf_recosmear = new THnSparseD("hist_leadjetunf_recosmear", "hist_leadjetunf_recosmear", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_leadjetunf_gensmear_4D = new THnSparseD("hist_leadjetunf_gensmear_4D", "hist_leadjetunf_gensmear_4D", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_leadjetunf_recosmear_4D = new THnSparseD("hist_leadjetunf_recosmear_4D", "hist_leadjetunf_recosmear_4D", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_subljetunf_gensmear = new THnSparseD("hist_subljetunf_gensmear", "hist_subljetunf_gensmear", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_subljetunf_recosmear = new THnSparseD("hist_subljetunf_recosmear", "hist_subljetunf_recosmear", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_subljetunf_gensmear_4D = new THnSparseD("hist_subljetunf_gensmear_4D", "hist_subljetunf_gensmear_4D", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_subljetunf_recosmear_4D = new THnSparseD("hist_subljetunf_recosmear_4D", "hist_subljetunf_recosmear_4D", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);

THnSparseD *hist_leadjetunf_gensmear_fromInclJet = new THnSparseD("hist_leadjetunf_gensmear_fromInclJet", "hist_leadjetunf_gensmear_fromInclJet", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_leadjetunf_recosmear_fromInclJet = new THnSparseD("hist_leadjetunf_recosmear_fromInclJet", "hist_leadjetunf_recosmear_fromInclJet", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_subljetunf_gensmear_fromInclJet = new THnSparseD("hist_subljetunf_gensmear_fromInclJet", "hist_subljetunf_gensmear_fromInclJet", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);
THnSparseD *hist_subljetunf_recosmear_fromInclJet = new THnSparseD("hist_subljetunf_recosmear_fromInclJet", "hist_subljetunf_recosmear_fromInclJet", 3, bins_jetunf_smear, xmin_jetunf_smear, xmax_jetunf_smear);

int    bins_xjunf_smear[3]   =   { nXjAjBins  , multbinsize-1                               , extrabinsize-1};
double xmin_xjunf_smear[3]   =   { 0    	  , multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_xjunf_smear[3]   =   { 1.5 		  , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_xjunf_gensmear = new THnSparseD("hist_xjunf_gensmear", "hist_xjunf_gensmear", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);
THnSparseD *hist_xjunf_recosmear = new THnSparseD("hist_xjunf_recosmear", "hist_xjunf_recosmear", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);
THnSparseD *hist_xjunf_gensmear_fromLSL = new THnSparseD("hist_xjunf_gensmear_fromLSL", "hist_xjunf_gensmear_fromLSL", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);
THnSparseD *hist_xjunf_recosmear_fromLSL = new THnSparseD("hist_xjunf_recosmear_fromLSL", "hist_xjunf_recosmear_fromLSL", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);
THnSparseD *hist_xjunf_gensmear_fromLSL_4D= new THnSparseD("hist_xjunf_gensmear_fromLSL_4D", "hist_xjunf_gensmear_fromLSL_4D", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);
THnSparseD *hist_xjunf_recosmear_fromLSL_4D = new THnSparseD("hist_xjunf_recosmear_fromLSL_4D", "hist_xjunf_recosmear_fromLSL_4D", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);
THnSparseD *hist_xjunf_gensmear_fromLSL_InclJet = new THnSparseD("hist_xjunf_gensmear_fromLSL_InclJet", "hist_xjunf_gensmear_fromLSL_InclJet", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);
THnSparseD *hist_xjunf_recosmear_fromLSL_InclJet = new THnSparseD("hist_xjunf_recosmear_fromLSL_InclJet", "hist_xjunf_recosmear_fromLSL_InclJet", 3, bins_xjunf_smear, xmin_xjunf_smear, xmax_xjunf_smear);

// --------------------------------------------------------------------------------------------------------
// 3rd jet studies
// Axis : 0 -> XJ12, 1 -> DPhi12, 2 -> XJ13, 3 ->DPhi13, 4 -> XJ23, 5 ->DPhi23 , 6 -> multiplicity, 7 -> extra dependency
int	bins_3rdjet[8]   	=     { nXjBins3rd   , nDphiBins3rd		, nXjBins3rd	, nDphiBins3rd	 , nXjBins3rd	, nDphiBins3rd   , multbinsize-1		  	 					 ,  extrabinsize-1};
double xmin_3rdjet[8]   =  	  { XjBins3rd[0] , DphiBins3rd[0]	, XjBins3rd[0]	, DphiBins3rd[0] , XjBins3rd[0]	, DphiBins3rd[0] , multiplicity_centrality_bins[0]		   	     ,  extra_bins[0]};
double xmax_3rdjet[8]   =     { 1.0  		 , TMath::Pi()  	, 1.0 			, TMath::Pi() 	 , 1.0 			, TMath::Pi() 	 , multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *hist_reco_3rdjet = new THnSparseD("hist_reco_3rdjet", "hist_reco_3rdjet", 8, bins_3rdjet, xmin_3rdjet, xmax_3rdjet);
THnSparseD *hist_ref_3rdjet = new THnSparseD("hist_ref_3rdjet", "hist_ref_3rdjet", 8, bins_3rdjet, xmin_3rdjet, xmax_3rdjet);
THnSparseD *hist_gen_3rdjet = new THnSparseD("hist_gen_3rdjet", "hist_gen_3rdjet", 8, bins_3rdjet, xmin_3rdjet, xmax_3rdjet);

// Quenching studies
// Axis : 0 -> Xj, 1 -> Aj, 2 -> delta phi, 3 -> multiplicity , 4 -> jet pT average, 5 -> extra dependency, 6 -> pT leading jet, 7 -> pT subleading jet
int	bins_quenc[8]   =      { nXjAjBins   , nXjAjBins	  , nDphiBins		, multbinsize-1		  	 					  ,	 nPtaveBins  ,  extrabinsize-1 				, nPtLSLBins, nPtLSLBins};
double xmin_quenc[8]   =   { 0.0 		 , 0.0  		  , 0.0		    	, multiplicity_centrality_bins[0]		   	  ,	 minPtave	 ,	extra_bins[0]				, minPtLSL	, minPtLSL};
double xmax_quenc[8]   =   { 1.5  		 , 1.5   		  , TMath::Pi() 	, multiplicity_centrality_bins[multbinsize-1] ,  maxPtave    ,  extra_bins[extrabinsize-1]  , maxPtLSL	, maxPtLSL};
THnSparseD *hist_reco_lead_reco_subl_quench_mid_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_mid", "hist_reco_lead_reco_subl_quench_mid_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_mid_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_fwd", "hist_reco_lead_reco_subl_quench_mid_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_mid_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_bkw", "hist_reco_lead_reco_subl_quench_mid_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_mid", "hist_reco_lead_reco_subl_quench_fwd_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_mid", "hist_reco_lead_reco_subl_quench_bkw_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_fwd", "hist_reco_lead_reco_subl_quench_fwd_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_bkw", "hist_reco_lead_reco_subl_quench_fwd_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_fwd", "hist_reco_lead_reco_subl_quench_bkw_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_bkw", "hist_reco_lead_reco_subl_quench_bkw_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_mid", "hist_gen_lead_gen_subl_quench_mid_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_fwd", "hist_gen_lead_gen_subl_quench_mid_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_bkw", "hist_gen_lead_gen_subl_quench_mid_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_mid", "hist_gen_lead_gen_subl_quench_fwd_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_mid", "hist_gen_lead_gen_subl_quench_bkw_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_fwd", "hist_gen_lead_gen_subl_quench_fwd_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_bkw", "hist_gen_lead_gen_subl_quench_fwd_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_fwd", "hist_gen_lead_gen_subl_quench_bkw_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_bkw", "hist_gen_lead_gen_subl_quench_bkw_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_mid", "hist_ref_lead_ref_subl_quench_mid_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_fwd", "hist_ref_lead_ref_subl_quench_mid_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_bkw", "hist_ref_lead_ref_subl_quench_mid_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_mid", "hist_ref_lead_ref_subl_quench_fwd_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_mid", "hist_ref_lead_ref_subl_quench_bkw_mid", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_fwd", "hist_ref_lead_ref_subl_quench_fwd_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_bkw", "hist_ref_lead_ref_subl_quench_fwd_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_fwd", "hist_ref_lead_ref_subl_quench_bkw_fwd", 8, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_bkw", "hist_ref_lead_ref_subl_quench_bkw_bkw", 8, bins_quenc, xmin_quenc, xmax_quenc);

// EP dependency
// Axis : 0 -> Xj, 1 -> delta phi, 2 -> 2*|Psi2 - jetphi|, 3 -> 3*|Psi3 - jetphi|, 4 -> 4*|Psi4 - jetphi|, 5 -> multiplicity , 6 -> jet pT average
int	bins_quencEP[7]   =      { nXjAjBins   , nDphiBins		, 16 		 , 16 			, 16 			, multbinsize-1		  	 					  ,	 nPtaveBins};
double xmin_quencEP[7]   =   { 0.0 		   , 0.0		  	, 0.0		 , 0.0			, 0.0			, multiplicity_centrality_bins[0]		   	  ,	 minPtave};
double xmax_quencEP[7]   =   { 1.5  	   , TMath::Pi() 	, TMath::Pi(), TMath::Pi()	, TMath::Pi()	, multiplicity_centrality_bins[multbinsize-1] ,  maxPtave};
THnSparseD *hist_reco_leadEP_quench_plus_mid_mid = new THnSparseD("hist_reco_leadEP_quench_plus_mid_mid", "hist_reco_leadEP_quench_plus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_mid_fwd = new THnSparseD("hist_reco_leadEP_quench_plus_mid_fwd", "hist_reco_leadEP_quench_plus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_mid_bkw = new THnSparseD("hist_reco_leadEP_quench_plus_mid_bkw", "hist_reco_leadEP_quench_plus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_fwd_mid = new THnSparseD("hist_reco_leadEP_quench_plus_fwd_mid", "hist_reco_leadEP_quench_plus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_bkw_mid = new THnSparseD("hist_reco_leadEP_quench_plus_bkw_mid", "hist_reco_leadEP_quench_plus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_fwd_fwd = new THnSparseD("hist_reco_leadEP_quench_plus_fwd_fwd", "hist_reco_leadEP_quench_plus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_fwd_bkw = new THnSparseD("hist_reco_leadEP_quench_plus_fwd_bkw", "hist_reco_leadEP_quench_plus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_bkw_fwd = new THnSparseD("hist_reco_leadEP_quench_plus_bkw_fwd", "hist_reco_leadEP_quench_plus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_plus_bkw_bkw = new THnSparseD("hist_reco_leadEP_quench_plus_bkw_bkw", "hist_reco_leadEP_quench_plus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_mid_mid = new THnSparseD("hist_reco_sublEP_quench_plus_mid_mid", "hist_reco_sublEP_quench_plus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_mid_fwd = new THnSparseD("hist_reco_sublEP_quench_plus_mid_fwd", "hist_reco_sublEP_quench_plus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_mid_bkw = new THnSparseD("hist_reco_sublEP_quench_plus_mid_bkw", "hist_reco_sublEP_quench_plus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_fwd_mid = new THnSparseD("hist_reco_sublEP_quench_plus_fwd_mid", "hist_reco_sublEP_quench_plus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_bkw_mid = new THnSparseD("hist_reco_sublEP_quench_plus_bkw_mid", "hist_reco_sublEP_quench_plus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_fwd_fwd = new THnSparseD("hist_reco_sublEP_quench_plus_fwd_fwd", "hist_reco_sublEP_quench_plus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_fwd_bkw = new THnSparseD("hist_reco_sublEP_quench_plus_fwd_bkw", "hist_reco_sublEP_quench_plus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_bkw_fwd = new THnSparseD("hist_reco_sublEP_quench_plus_bkw_fwd", "hist_reco_sublEP_quench_plus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus_bkw_bkw = new THnSparseD("hist_reco_sublEP_quench_plus_bkw_bkw", "hist_reco_sublEP_quench_plus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_mid_mid = new THnSparseD("hist_ref_leadEP_quench_plus_mid_mid", "hist_ref_leadEP_quench_plus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_mid_fwd = new THnSparseD("hist_ref_leadEP_quench_plus_mid_fwd", "hist_ref_leadEP_quench_plus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_mid_bkw = new THnSparseD("hist_ref_leadEP_quench_plus_mid_bkw", "hist_ref_leadEP_quench_plus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_fwd_mid = new THnSparseD("hist_ref_leadEP_quench_plus_fwd_mid", "hist_ref_leadEP_quench_plus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_bkw_mid = new THnSparseD("hist_ref_leadEP_quench_plus_bkw_mid", "hist_ref_leadEP_quench_plus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_fwd_fwd = new THnSparseD("hist_ref_leadEP_quench_plus_fwd_fwd", "hist_ref_leadEP_quench_plus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_fwd_bkw = new THnSparseD("hist_ref_leadEP_quench_plus_fwd_bkw", "hist_ref_leadEP_quench_plus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_bkw_fwd = new THnSparseD("hist_ref_leadEP_quench_plus_bkw_fwd", "hist_ref_leadEP_quench_plus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus_bkw_bkw = new THnSparseD("hist_ref_leadEP_quench_plus_bkw_bkw", "hist_ref_leadEP_quench_plus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_mid_mid = new THnSparseD("hist_ref_sublEP_quench_plus_mid_mid", "hist_ref_sublEP_quench_plus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_mid_fwd = new THnSparseD("hist_ref_sublEP_quench_plus_mid_fwd", "hist_ref_sublEP_quench_plus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_mid_bkw = new THnSparseD("hist_ref_sublEP_quench_plus_mid_bkw", "hist_ref_sublEP_quench_plus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_fwd_mid = new THnSparseD("hist_ref_sublEP_quench_plus_fwd_mid", "hist_ref_sublEP_quench_plus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_bkw_mid = new THnSparseD("hist_ref_sublEP_quench_plus_bkw_mid", "hist_ref_sublEP_quench_plus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_fwd_fwd = new THnSparseD("hist_ref_sublEP_quench_plus_fwd_fwd", "hist_ref_sublEP_quench_plus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_fwd_bkw = new THnSparseD("hist_ref_sublEP_quench_plus_fwd_bkw", "hist_ref_sublEP_quench_plus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_bkw_fwd = new THnSparseD("hist_ref_sublEP_quench_plus_bkw_fwd", "hist_ref_sublEP_quench_plus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus_bkw_bkw = new THnSparseD("hist_ref_sublEP_quench_plus_bkw_bkw", "hist_ref_sublEP_quench_plus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_mid_mid = new THnSparseD("hist_gen_leadEP_quench_plus_mid_mid", "hist_gen_leadEP_quench_plus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_mid_fwd = new THnSparseD("hist_gen_leadEP_quench_plus_mid_fwd", "hist_gen_leadEP_quench_plus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_mid_bkw = new THnSparseD("hist_gen_leadEP_quench_plus_mid_bkw", "hist_gen_leadEP_quench_plus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_fwd_mid = new THnSparseD("hist_gen_leadEP_quench_plus_fwd_mid", "hist_gen_leadEP_quench_plus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_bkw_mid = new THnSparseD("hist_gen_leadEP_quench_plus_bkw_mid", "hist_gen_leadEP_quench_plus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_fwd_fwd = new THnSparseD("hist_gen_leadEP_quench_plus_fwd_fwd", "hist_gen_leadEP_quench_plus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_fwd_bkw = new THnSparseD("hist_gen_leadEP_quench_plus_fwd_bkw", "hist_gen_leadEP_quench_plus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_bkw_fwd = new THnSparseD("hist_gen_leadEP_quench_plus_bkw_fwd", "hist_gen_leadEP_quench_plus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus_bkw_bkw = new THnSparseD("hist_gen_leadEP_quench_plus_bkw_bkw", "hist_gen_leadEP_quench_plus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_mid_mid = new THnSparseD("hist_gen_sublEP_quench_plus_mid_mid", "hist_gen_sublEP_quench_plus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_mid_fwd = new THnSparseD("hist_gen_sublEP_quench_plus_mid_fwd", "hist_gen_sublEP_quench_plus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_mid_bkw = new THnSparseD("hist_gen_sublEP_quench_plus_mid_bkw", "hist_gen_sublEP_quench_plus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_fwd_mid = new THnSparseD("hist_gen_sublEP_quench_plus_fwd_mid", "hist_gen_sublEP_quench_plus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_bkw_mid = new THnSparseD("hist_gen_sublEP_quench_plus_bkw_mid", "hist_gen_sublEP_quench_plus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_fwd_fwd = new THnSparseD("hist_gen_sublEP_quench_plus_fwd_fwd", "hist_gen_sublEP_quench_plus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_fwd_bkw = new THnSparseD("hist_gen_sublEP_quench_plus_fwd_bkw", "hist_gen_sublEP_quench_plus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_bkw_fwd = new THnSparseD("hist_gen_sublEP_quench_plus_bkw_fwd", "hist_gen_sublEP_quench_plus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus_bkw_bkw = new THnSparseD("hist_gen_sublEP_quench_plus_bkw_bkw", "hist_gen_sublEP_quench_plus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_mid_mid = new THnSparseD("hist_reco_leadEP_quench_minus_mid_mid", "hist_reco_leadEP_quench_minus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_mid_fwd = new THnSparseD("hist_reco_leadEP_quench_minus_mid_fwd", "hist_reco_leadEP_quench_minus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_mid_bkw = new THnSparseD("hist_reco_leadEP_quench_minus_mid_bkw", "hist_reco_leadEP_quench_minus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_fwd_mid = new THnSparseD("hist_reco_leadEP_quench_minus_fwd_mid", "hist_reco_leadEP_quench_minus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_bkw_mid = new THnSparseD("hist_reco_leadEP_quench_minus_bkw_mid", "hist_reco_leadEP_quench_minus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_fwd_fwd = new THnSparseD("hist_reco_leadEP_quench_minus_fwd_fwd", "hist_reco_leadEP_quench_minus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_fwd_bkw = new THnSparseD("hist_reco_leadEP_quench_minus_fwd_bkw", "hist_reco_leadEP_quench_minus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_bkw_fwd = new THnSparseD("hist_reco_leadEP_quench_minus_bkw_fwd", "hist_reco_leadEP_quench_minus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus_bkw_bkw = new THnSparseD("hist_reco_leadEP_quench_minus_bkw_bkw", "hist_reco_leadEP_quench_minus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_mid_mid = new THnSparseD("hist_reco_sublEP_quench_minus_mid_mid", "hist_reco_sublEP_quench_minus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_mid_fwd = new THnSparseD("hist_reco_sublEP_quench_minus_mid_fwd", "hist_reco_sublEP_quench_minus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_mid_bkw = new THnSparseD("hist_reco_sublEP_quench_minus_mid_bkw", "hist_reco_sublEP_quench_minus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_fwd_mid = new THnSparseD("hist_reco_sublEP_quench_minus_fwd_mid", "hist_reco_sublEP_quench_minus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_bkw_mid = new THnSparseD("hist_reco_sublEP_quench_minus_bkw_mid", "hist_reco_sublEP_quench_minus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_fwd_fwd = new THnSparseD("hist_reco_sublEP_quench_minus_fwd_fwd", "hist_reco_sublEP_quench_minus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_fwd_bkw = new THnSparseD("hist_reco_sublEP_quench_minus_fwd_bkw", "hist_reco_sublEP_quench_minus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_bkw_fwd = new THnSparseD("hist_reco_sublEP_quench_minus_bkw_fwd", "hist_reco_sublEP_quench_minus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus_bkw_bkw = new THnSparseD("hist_reco_sublEP_quench_minus_bkw_bkw", "hist_reco_sublEP_quench_minus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_mid_mid = new THnSparseD("hist_ref_leadEP_quench_minus_mid_mid", "hist_ref_leadEP_quench_minus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_mid_fwd = new THnSparseD("hist_ref_leadEP_quench_minus_mid_fwd", "hist_ref_leadEP_quench_minus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_mid_bkw = new THnSparseD("hist_ref_leadEP_quench_minus_mid_bkw", "hist_ref_leadEP_quench_minus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_fwd_mid = new THnSparseD("hist_ref_leadEP_quench_minus_fwd_mid", "hist_ref_leadEP_quench_minus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_bkw_mid = new THnSparseD("hist_ref_leadEP_quench_minus_bkw_mid", "hist_ref_leadEP_quench_minus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_fwd_fwd = new THnSparseD("hist_ref_leadEP_quench_minus_fwd_fwd", "hist_ref_leadEP_quench_minus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_fwd_bkw = new THnSparseD("hist_ref_leadEP_quench_minus_fwd_bkw", "hist_ref_leadEP_quench_minus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_bkw_fwd = new THnSparseD("hist_ref_leadEP_quench_minus_bkw_fwd", "hist_ref_leadEP_quench_minus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus_bkw_bkw = new THnSparseD("hist_ref_leadEP_quench_minus_bkw_bkw", "hist_ref_leadEP_quench_minus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_mid_mid = new THnSparseD("hist_ref_sublEP_quench_minus_mid_mid", "hist_ref_sublEP_quench_minus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_mid_fwd = new THnSparseD("hist_ref_sublEP_quench_minus_mid_fwd", "hist_ref_sublEP_quench_minus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_mid_bkw = new THnSparseD("hist_ref_sublEP_quench_minus_mid_bkw", "hist_ref_sublEP_quench_minus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_fwd_mid = new THnSparseD("hist_ref_sublEP_quench_minus_fwd_mid", "hist_ref_sublEP_quench_minus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_bkw_mid = new THnSparseD("hist_ref_sublEP_quench_minus_bkw_mid", "hist_ref_sublEP_quench_minus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_fwd_fwd = new THnSparseD("hist_ref_sublEP_quench_minus_fwd_fwd", "hist_ref_sublEP_quench_minus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_fwd_bkw = new THnSparseD("hist_ref_sublEP_quench_minus_fwd_bkw", "hist_ref_sublEP_quench_minus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_bkw_fwd = new THnSparseD("hist_ref_sublEP_quench_minus_bkw_fwd", "hist_ref_sublEP_quench_minus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus_bkw_bkw = new THnSparseD("hist_ref_sublEP_quench_minus_bkw_bkw", "hist_ref_sublEP_quench_minus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_mid_mid = new THnSparseD("hist_gen_leadEP_quench_minus_mid_mid", "hist_gen_leadEP_quench_minus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_mid_fwd = new THnSparseD("hist_gen_leadEP_quench_minus_mid_fwd", "hist_gen_leadEP_quench_minus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_mid_bkw = new THnSparseD("hist_gen_leadEP_quench_minus_mid_bkw", "hist_gen_leadEP_quench_minus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_fwd_mid = new THnSparseD("hist_gen_leadEP_quench_minus_fwd_mid", "hist_gen_leadEP_quench_minus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_bkw_mid = new THnSparseD("hist_gen_leadEP_quench_minus_bkw_mid", "hist_gen_leadEP_quench_minus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_fwd_fwd = new THnSparseD("hist_gen_leadEP_quench_minus_fwd_fwd", "hist_gen_leadEP_quench_minus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_fwd_bkw = new THnSparseD("hist_gen_leadEP_quench_minus_fwd_bkw", "hist_gen_leadEP_quench_minus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_bkw_fwd = new THnSparseD("hist_gen_leadEP_quench_minus_bkw_fwd", "hist_gen_leadEP_quench_minus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus_bkw_bkw = new THnSparseD("hist_gen_leadEP_quench_minus_bkw_bkw", "hist_gen_leadEP_quench_minus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_mid_mid = new THnSparseD("hist_gen_sublEP_quench_minus_mid_mid", "hist_gen_sublEP_quench_minus_mid_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_mid_fwd = new THnSparseD("hist_gen_sublEP_quench_minus_mid_fwd", "hist_gen_sublEP_quench_minus_mid_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_mid_bkw = new THnSparseD("hist_gen_sublEP_quench_minus_mid_bkw", "hist_gen_sublEP_quench_minus_mid_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_fwd_mid = new THnSparseD("hist_gen_sublEP_quench_minus_fwd_mid", "hist_gen_sublEP_quench_minus_fwd_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_bkw_mid = new THnSparseD("hist_gen_sublEP_quench_minus_bkw_mid", "hist_gen_sublEP_quench_minus_bkw_mid", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_fwd_fwd = new THnSparseD("hist_gen_sublEP_quench_minus_fwd_fwd", "hist_gen_sublEP_quench_minus_fwd_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_fwd_bkw = new THnSparseD("hist_gen_sublEP_quench_minus_fwd_bkw", "hist_gen_sublEP_quench_minus_fwd_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_bkw_fwd = new THnSparseD("hist_gen_sublEP_quench_minus_bkw_fwd", "hist_gen_sublEP_quench_minus_bkw_fwd", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus_bkw_bkw = new THnSparseD("hist_gen_sublEP_quench_minus_bkw_bkw", "hist_gen_sublEP_quench_minus_bkw_bkw", 7, bins_quencEP, xmin_quencEP, xmax_quencEP);
// Assymetry studies
// Axis : 0 -> etaDijet, 1 -> delta eta / 2, 2 -> Xj, 3 -> Aj, 4 -> delta phi, 5 -> x_p, 6 -> x_Pb, 7-> M12, 8 -> multiplicity, 9 -> jet pT average, 10 -> extra dependency
int	bins_etaDijet[11]      =   {  40   ,  16  , nXjAjBins	  , nXjAjBins , nDphiBins	    , nXBins ,  nXBins  ,  nMBins,	 multbinsize-1		  	 					  ,	 nPtaveBins	 ,  extrabinsize-1};
double xmin_etaDijet[11]   =   { -4.0  , -4.0 ,  0.0   		  , 0.0		  , 0.0 	 	    , minX   ,  minX    ,  minM, 	 multiplicity_centrality_bins[0]		   	  ,	 minPtave	 ,  extra_bins[0]};
double xmax_etaDijet[11]   =   {  4.0  ,  4.0 ,  1.5   		  , 1.5		  , TMath::Pi() 	, maxX   ,  maxX    ,  maxM, 	 multiplicity_centrality_bins[multbinsize-1]  ,  maxPtave    ,  extra_bins[extrabinsize-1]};
THnSparseD *hist_etaDijet_reco = new THnSparseD("hist_etaDijet_reco", "hist_etaDijet_reco", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_reco = new THnSparseD("hist_etaDijet_CM_reco", "hist_etaDijet_CM_reco", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_ref = new THnSparseD("hist_etaDijet_ref", "hist_etaDijet_ref", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_ref = new THnSparseD("hist_etaDijet_CM_ref", "hist_etaDijet_CM_ref", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_gen = new THnSparseD("hist_etaDijet_gen", "hist_etaDijet_gen", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_gen = new THnSparseD("hist_etaDijet_CM_gen", "hist_etaDijet_CM_gen", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
// rapidity dependency
THnSparseD *hist_yDijet_CM_reco = new THnSparseD("hist_yDijet_CM_reco", "hist_yDijet_CM_reco", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_yDijet_CM_ref = new THnSparseD("hist_yDijet_CM_ref", "hist_yDijet_CM_ref", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_yDijet_CM_gen = new THnSparseD("hist_yDijet_CM_gen", "hist_yDijet_CM_gen", 11, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);

// Axis : 0 -> in-jet multiplicity, 1 -> multiplicity, 2 -> extra dimension
int	bins_injettrk[3]   	  =   { 100 , multbinsize-1			  						  ,  extrabinsize-1};
double xmin_injettrk[3]   =   { 0.0 , multiplicity_centrality_bins[0]				  ,	 extra_bins[0]};
double xmax_injettrk[3]   =   { 100 , multiplicity_centrality_bins[multbinsize-1]	  ,  extra_bins[extrabinsize-1]};
THnSparseD *hist_injet_reco_track_reco = new THnSparseD("hist_injet_reco_track_reco","hist_injet_reco_track_reco",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_injet_reco_track_gen = new THnSparseD("hist_injet_reco_track_gen","hist_injet_reco_track_gen",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_injet_gen_track_reco = new THnSparseD("hist_injet_gen_track_reco","hist_injet_gen_track_reco",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_injet_gen_track_gen = new THnSparseD("hist_injet_gen_track_gen","hist_injet_gen_track_gen",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inLeadjet_reco_track_reco = new THnSparseD("hist_inLeadjet_reco_track_reco","hist_inLeadjet_reco_track_reco",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inLeadjet_reco_track_gen = new THnSparseD("hist_inLeadjet_reco_track_gen","hist_inLeadjet_reco_track_gen",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inLeadjet_gen_track_reco = new THnSparseD("hist_inLeadjet_gen_track_reco","hist_inLeadjet_gen_track_reco",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inLeadjet_gen_track_gen = new THnSparseD("hist_inLeadjet_gen_track_gen","hist_inLeadjet_gen_track_gen",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inSubljet_reco_track_reco = new THnSparseD("hist_inSubljet_reco_track_reco","hist_inSubljet_reco_track_reco",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inSubljet_reco_track_gen = new THnSparseD("hist_inSubljet_reco_track_gen","hist_inSubljet_reco_track_gen",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inSubljet_gen_track_reco = new THnSparseD("hist_inSubljet_gen_track_reco","hist_inSubljet_gen_track_reco",3,bins_injettrk,xmin_injettrk,xmax_injettrk);
THnSparseD *hist_inSubljet_gen_track_gen = new THnSparseD("hist_inSubljet_gen_track_gen","hist_inSubljet_gen_track_gen",3,bins_injettrk,xmin_injettrk,xmax_injettrk);

// Correlation studies
// Axis : 0 -> delta phi, 1 -> delta eta, 2 -> track pT, 3 -> multiplicity, 4 -> extra dimension
int	bins_jettrk[5]      =   { 200				    , 400  ,   trkbinsize-1			  	 	 	 , multbinsize-1			  						  ,  extrabinsize-1};
double xmin_jettrk[5]   =   { -TMath::Pi()/2.0		, -4.0 ,   trk_pt_bins[0]					 , multiplicity_centrality_bins[0]					  ,	 extra_bins[0]};
double xmax_jettrk[5]   =   { 3.0*TMath::Pi()/2.0 	, 4.0  ,   trk_pt_bins[trkbinsize-1] 		 , multiplicity_centrality_bins[multbinsize-1]  	  ,  extra_bins[extrabinsize-1]};

// Correlation: Reco Jet + Reco Track
THnSparseD *hist_correlation_signal_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_jet_reco_track_reco","hist_correlation_signal_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_jet_reco_track_reco = new THnSparseD("hist_correlation_rotation_jet_reco_track_reco","hist_correlation_rotation_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_jet_reco_track_reco = new THnSparseD("hist_correlation_mixing_jet_reco_track_reco","hist_correlation_mixing_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_lead_jet_reco_track_reco","hist_correlation_signal_lead_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_rotation_lead_jet_reco_track_reco","hist_correlation_rotation_lead_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_mixing_lead_jet_reco_track_reco","hist_correlation_mixing_lead_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subl_jet_reco_track_reco","hist_correlation_signal_subl_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_rotation_subl_jet_reco_track_reco","hist_correlation_rotation_subl_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_mixing_subl_jet_reco_track_reco","hist_correlation_mixing_subl_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);

// Correlation: Reco Jet + Gen Track
THnSparseD *hist_correlation_signal_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_jet_reco_track_gen","hist_correlation_signal_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_jet_reco_track_gen = new THnSparseD("hist_correlation_rotation_jet_reco_track_gen","hist_correlation_rotation_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_jet_reco_track_gen = new THnSparseD("hist_correlation_mixing_jet_reco_track_gen","hist_correlation_mixing_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_lead_jet_reco_track_gen","hist_correlation_signal_lead_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_rotation_lead_jet_reco_track_gen","hist_correlation_rotation_lead_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_mixing_lead_jet_reco_track_gen","hist_correlation_mixing_lead_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subl_jet_reco_track_gen","hist_correlation_signal_subl_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_rotation_subl_jet_reco_track_gen","hist_correlation_rotation_subl_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_mixing_subl_jet_reco_track_gen","hist_correlation_mixing_subl_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);

// Correlation: Gen Jet + Reco Track
THnSparseD *hist_correlation_signal_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_jet_gen_track_reco","hist_correlation_signal_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_jet_gen_track_reco = new THnSparseD("hist_correlation_rotation_jet_gen_track_reco","hist_correlation_rotation_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_jet_gen_track_reco = new THnSparseD("hist_correlation_mixing_jet_gen_track_reco","hist_correlation_mixing_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_lead_jet_gen_track_reco","hist_correlation_signal_lead_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_rotation_lead_jet_gen_track_reco","hist_correlation_rotation_lead_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_mixing_lead_jet_gen_track_reco","hist_correlation_mixing_lead_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subl_jet_gen_track_reco","hist_correlation_signal_subl_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_rotation_subl_jet_gen_track_reco","hist_correlation_rotation_subl_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_mixing_subl_jet_gen_track_reco","hist_correlation_mixing_subl_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);

// Correlation: Gen Jet + Gen Track
THnSparseD *hist_correlation_signal_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_jet_gen_track_gen","hist_correlation_signal_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_jet_gen_track_gen = new THnSparseD("hist_correlation_rotation_jet_gen_track_gen","hist_correlation_rotation_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_jet_gen_track_gen = new THnSparseD("hist_correlation_mixing_jet_gen_track_gen","hist_correlation_mixing_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_lead_jet_gen_track_gen","hist_correlation_signal_lead_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_rotation_lead_jet_gen_track_gen","hist_correlation_rotation_lead_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_mixing_lead_jet_gen_track_gen","hist_correlation_mixing_lead_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subl_jet_gen_track_gen","hist_correlation_signal_subl_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_rotation_subl_jet_gen_track_gen","hist_correlation_rotation_subl_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_mixing_subl_jet_gen_track_gen","hist_correlation_mixing_subl_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);

// Correlation: Include sube > 0
THnSparseD *hist_correlation_signal_subg0_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subg0_jet_reco_track_reco","hist_correlation_signal_subg0_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subg0_jet_reco_track_gen","hist_correlation_signal_subg0_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subg0_jet_gen_track_reco","hist_correlation_signal_subg0_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subg0_jet_gen_track_gen","hist_correlation_signal_subg0_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);

THnSparseD *hist_correlation_signal_subg0_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subg0_lead_jet_reco_track_reco","hist_correlation_signal_subg0_lead_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subg0_lead_jet_reco_track_gen","hist_correlation_signal_subg0_lead_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subg0_lead_jet_gen_track_reco","hist_correlation_signal_subg0_lead_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subg0_lead_jet_gen_track_gen","hist_correlation_signal_subg0_lead_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);

THnSparseD *hist_correlation_signal_subg0_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subg0_subl_jet_reco_track_reco","hist_correlation_signal_subg0_subl_jet_reco_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subg0_subl_jet_reco_track_gen","hist_correlation_signal_subg0_subl_jet_reco_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subg0_subl_jet_gen_track_reco","hist_correlation_signal_subg0_subl_jet_gen_track_reco",5,bins_jettrk,xmin_jettrk,xmax_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subg0_subl_jet_gen_track_gen","hist_correlation_signal_subg0_subl_jet_gen_track_gen",5,bins_jettrk,xmin_jettrk,xmax_jettrk);

// Two particle correlation studies
// Axis : 0 -> delta phi, 1 -> delta eta, 2 -> track pT, 3 -> multiplicity, 4 -> extra dimension
int	bins_2pc[5]      =   { 40				     ,  80   ,   trkbinsize-1		    			, multbinsize-1			  					   ,  extrabinsize-1};
double xmin_2pc[5]   =   { -TMath::Pi()/2.0	     , -4.0  ,   trk_pt_bins[0]					    , multiplicity_centrality_bins[0]			   ,  extra_bins[0]};
double xmax_2pc[5]   =   { 3.0*TMath::Pi()/2.0   ,  4.0  ,   trk_pt_bins[trkbinsize-1]  		, multiplicity_centrality_bins[multbinsize-1]  ,  extra_bins[extrabinsize-1]};
// 2 particle correlations for flow analysis
THnSparseD *hist_reco_reco_2pcorrelation_signal = new THnSparseD("hist_reco_reco_2pcorrelation_signal","hist_reco_reco_2pcorrelation_signal",5,bins_2pc,xmin_2pc,xmax_2pc);
THnSparseD *hist_reco_reco_2pcorrelation_signal_subg0 = new THnSparseD("hist_reco_reco_2pcorrelation_signal_subg0","hist_reco_reco_2pcorrelation_signal_subg0",5,bins_2pc,xmin_2pc,xmax_2pc);
THnSparseD *hist_reco_reco_2pcorrelation_signal_subcross = new THnSparseD("hist_reco_reco_2pcorrelation_signal_subcross","hist_reco_reco_2pcorrelation_signal_subcross",5,bins_2pc,xmin_2pc,xmax_2pc);
THnSparseD *hist_reco_reco_2pcorrelation_mixing = new THnSparseD("hist_reco_reco_2pcorrelation_mixing","hist_reco_reco_2pcorrelation_mixing",5,bins_2pc,xmin_2pc,xmax_2pc);
THnSparseD *hist_gen_gen_2pcorrelation_signal = new THnSparseD("hist_gen_gen_2pcorrelation_signal","hist_gen_gen_2pcorrelation_signal",5,bins_2pc,xmin_2pc,xmax_2pc);
THnSparseD *hist_gen_gen_2pcorrelation_signal_subg0 = new THnSparseD("hist_gen_gen_2pcorrelation_signal_subg0","hist_gen_gen_2pcorrelation_signal_subg0",5,bins_2pc,xmin_2pc,xmax_2pc);
THnSparseD *hist_gen_gen_2pcorrelation_signal_subcross = new THnSparseD("hist_gen_gen_2pcorrelation_signal_subcross","hist_gen_gen_2pcorrelation_signal_subcross",5,bins_2pc,xmin_2pc,xmax_2pc);
THnSparseD *hist_gen_gen_2pcorrelation_mixing = new THnSparseD("hist_gen_gen_2pcorrelation_mixing","hist_gen_gen_2pcorrelation_mixing",5,bins_2pc,xmin_2pc,xmax_2pc);

// Evaluate uncertainties correctly at ROOT
void sw2(){

double XBins[nXBins+1];
// Xp and XpPb dependency
for(int a = 0; a <= nXBins; a++){XBins[a] = (minX+binnerShift)*TMath::Exp(a*XlogBinWidth)-binnerShift;}
//adjust the bins for Xp and XPb
hist_etaDijet_reco->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_reco->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_CM_reco->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_CM_reco->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_ref->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_ref->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_CM_ref->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_CM_ref->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_gen->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_gen->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_CM_gen->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_CM_gen->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_yDijet_CM_reco->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_yDijet_CM_reco->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_yDijet_CM_ref->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_yDijet_CM_ref->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_yDijet_CM_gen->GetAxis(5)->Set(bins_etaDijet[5],XBins);
double MBins[nMBins+1];
for(int a = 0; a <= nMBins; a++){MBins[a] = (minM+binnerShift)*TMath::Exp(a*MlogBinWidth)-binnerShift;} // add bins by hand (bellow), because cannot loop here
hist_etaDijet_reco->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_etaDijet_CM_reco->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_etaDijet_ref->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_etaDijet_CM_ref->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_etaDijet_gen->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_etaDijet_CM_gen->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_yDijet_CM_reco->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_yDijet_CM_ref->GetAxis(7)->Set(bins_etaDijet[7],MBins);
hist_yDijet_CM_gen->GetAxis(7)->Set(bins_etaDijet[7],MBins);
// Xj, Aj and dphi
hist_etaDijet_reco->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_CM_reco->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_ref->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_CM_ref->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_gen->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_CM_gen->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_yDijet_CM_reco->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_yDijet_CM_ref->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_yDijet_CM_gen->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_reco_lead_reco_subl_quench_mid_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_mid_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_mid_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_fwd_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_bkw_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_fwd_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_fwd_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_bkw_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_reco_subl_quench_bkw_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_mid_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_mid_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_mid_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_fwd_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_bkw_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_fwd_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_fwd_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_bkw_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench_bkw_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_mid_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_mid_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_mid_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_fwd_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_bkw_mid->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_fwd_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_fwd_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_bkw_fwd->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench_bkw_bkw->GetAxis(0)->Set(bins_quenc[0],XjBins);

hist_reco_leadEP_quench_plus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_plus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_mid_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_mid_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_mid_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_fwd_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_bkw_mid->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_fwd_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_fwd_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_bkw_fwd->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus_bkw_bkw->GetAxis(0)->Set(bins_quencEP[0],XjBins);

hist_etaDijet_reco->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_CM_reco->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_ref->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_CM_ref->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_gen->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_CM_gen->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_yDijet_CM_reco->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_yDijet_CM_ref->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_yDijet_CM_gen->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_reco_lead_reco_subl_quench_mid_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_mid_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_mid_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_fwd_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_bkw_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_fwd_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_fwd_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_bkw_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_reco_subl_quench_bkw_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_mid_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_mid_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_mid_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_fwd_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_bkw_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_fwd_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_fwd_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_bkw_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench_bkw_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_mid_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_mid_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_mid_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_fwd_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_bkw_mid->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_fwd_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_fwd_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_bkw_fwd->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench_bkw_bkw->GetAxis(1)->Set(bins_quenc[1],AjBins);

hist_etaDijet_reco->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_CM_reco->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_ref->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_CM_ref->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_gen->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_CM_gen->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_yDijet_CM_reco->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_yDijet_CM_ref->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_yDijet_CM_gen->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_reco_lead_reco_subl_quench_mid_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_mid_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_mid_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_fwd_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_bkw_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_fwd_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_fwd_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_bkw_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_reco_subl_quench_bkw_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_mid_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_mid_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_mid_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_fwd_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_bkw_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_fwd_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_fwd_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_bkw_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench_bkw_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_mid_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_mid_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_mid_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_fwd_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_bkw_mid->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_fwd_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_fwd_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_bkw_fwd->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench_bkw_bkw->GetAxis(2)->Set(bins_quenc[2],DphiBins);

hist_reco_leadEP_quench_plus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_plus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);

hist_reco_leadEP_quench_minus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_mid_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_mid_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_mid_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_fwd_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_bkw_mid->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_fwd_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_fwd_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_bkw_fwd->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus_bkw_bkw->GetAxis(1)->Set(bins_quencEP[1],DphiBins);


// pT average
double PtaveBins[nPtaveBins+1];
// pT average
for(int a = 0; a <= nPtaveBins; a++){PtaveBins[a] = (minPtave+binnerShift)*TMath::Exp(a*PtavelogBinWidth)-binnerShift;} // add bins by hand (bellow), because cannot loop here
//adjust the bins for pT average
hist_etaDijet_reco->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_etaDijet_CM_reco->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_etaDijet_ref->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_etaDijet_CM_ref->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_etaDijet_gen->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_etaDijet_CM_gen->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_yDijet_CM_reco->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_yDijet_CM_ref->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
hist_yDijet_CM_gen->GetAxis(9)->Set(bins_etaDijet[9],PtaveBins);
//adjust the bins for pT average
hist_reco_lead_reco_subl_quench_mid_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_mid_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_mid_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_fwd_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_bkw_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_fwd_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_fwd_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_bkw_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_reco_lead_reco_subl_quench_bkw_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_mid_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_mid_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_mid_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_fwd_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_bkw_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_fwd_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_fwd_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_bkw_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_gen_lead_gen_subl_quench_bkw_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_mid_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_mid_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_mid_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_fwd_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_bkw_mid->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_fwd_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_fwd_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_bkw_fwd->GetAxis(4)->Set(bins_quenc[4],PtaveBins);
hist_ref_lead_ref_subl_quench_bkw_bkw->GetAxis(4)->Set(bins_quenc[4],PtaveBins);

hist_reco_leadEP_quench_plus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_plus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_plus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_plus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_plus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_plus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_plus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);

hist_reco_leadEP_quench_minus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_leadEP_quench_minus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_reco_sublEP_quench_minus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_leadEP_quench_minus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_ref_sublEP_quench_minus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_leadEP_quench_minus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_mid_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_mid_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_mid_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_fwd_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_bkw_mid->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_fwd_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_fwd_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_bkw_fwd->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);
hist_gen_sublEP_quench_minus_bkw_bkw->GetAxis(6)->Set(bins_quencEP[6],PtaveBins);

//Trk pt binning
double TrkPtbins[trkbinsize-1];
for(int a = 0; a<trk_pt_bins.size();a++){TrkPtbins[a] = trk_pt_bins[a];}
hist_correlation_signal_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_lead_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_lead_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_lead_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subl_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_subl_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_subl_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_lead_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_lead_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_lead_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subl_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_subl_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_subl_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_lead_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_lead_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_lead_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subl_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_subl_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_subl_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_lead_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_lead_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_lead_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subl_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_rotation_subl_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_mixing_subl_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_lead_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_lead_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_lead_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_lead_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_subl_jet_reco_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_subl_jet_reco_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_subl_jet_gen_track_reco->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_correlation_signal_subg0_subl_jet_gen_track_gen->GetAxis(2)->Set(bins_jettrk[2],TrkPtbins);
hist_reco_reco_2pcorrelation_signal->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_reco_reco_2pcorrelation_signal_subg0->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_reco_reco_2pcorrelation_signal_subcross->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_reco_reco_2pcorrelation_mixing->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_signal->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_signal_subg0->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_signal_subcross->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_mixing->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
Dphi_EP2_flat_trk_minus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_EP2_flat_trk_plus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_EP3_flat_trk_minus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_EP3_flat_trk_plus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_EP4_flat_trk_minus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_EP4_flat_trk_plus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_GEN_EP2_flat_trk_minus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_GEN_EP2_flat_trk_plus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_GEN_EP3_flat_trk_minus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_GEN_EP3_flat_trk_plus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_GEN_EP4_flat_trk_minus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
Dphi_GEN_EP4_flat_trk_plus->GetAxis(1)->Set(bins_TRKEP[1],TrkPtbins);
//Multiplicity/Centrality binning
double MultCentbins[multbinsize-1];
for(int a = 0; a<multiplicity_centrality_bins.size();a++){MultCentbins[a] = multiplicity_centrality_bins[a];}
vzhist->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vzhist_weighted->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vzhist_jet_weighted->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vzhist_dijet_weighted->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vxyhist->GetAxis(2)->Set(bins_vxy[2],MultCentbins);
vxyhist_weighted->GetAxis(2)->Set(bins_vxy[2],MultCentbins);
pthathist->GetAxis(1)->Set(bins_pthat[1],MultCentbins);
pthathist_weighted->GetAxis(1)->Set(bins_pthat[1],MultCentbins);
EP2_plus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP2_minus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP3_plus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP3_minus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP4_plus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP4_minus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
hist_reco_trk->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_reco_trk_corr->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_reco_trk_weighted->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_gen_trk->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_gen_trk_weighted->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_reco_reco_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_reco_gen_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_gen_reco_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_gen_gen_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_reco_reco_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_reco_gen_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_gen_reco_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_trk_from_gen_gen_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_reco_reco_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_reco_gen_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_gen_reco_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_gen_gen_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_reco_reco_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_reco_gen_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_gen_reco_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_LJ_trk_from_gen_gen_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_reco_reco_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_reco_gen_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_gen_reco_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_gen_gen_sig->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_reco_reco_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_reco_gen_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_gen_reco_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_SLJ_trk_from_gen_gen_mix->GetAxis(3)->Set(bins_trk[3],MultCentbins);
Dphi_EP2_flat_trk_minus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_EP2_flat_trk_plus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_EP3_flat_trk_minus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_EP3_flat_trk_plus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_EP4_flat_trk_minus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_EP4_flat_trk_plus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_GEN_EP2_flat_trk_minus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_GEN_EP2_flat_trk_plus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_GEN_EP3_flat_trk_minus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_GEN_EP3_flat_trk_plus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_GEN_EP4_flat_trk_minus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
Dphi_GEN_EP4_flat_trk_plus->GetAxis(2)->Set(bins_TRKEP[2],MultCentbins);
NJets->GetAxis(1)->Set(bins_NJETS[1],MultCentbins);
trackmaxptinjethisto->GetAxis(1)->Set(bins_trkmax[1],MultCentbins);
jettrackmaxptinjethisto->GetAxis(1)->Set(bins_trkmaxjet[1],MultCentbins);
histo_jetUE->GetAxis(1)->Set(bins_UE[1],MultCentbins);
histo_jetAverageRho->GetAxis(1)->Set(bins_UE[1],MultCentbins);
histo_jetcheckcorrection->GetAxis(1)->Set(bins_UE[1],MultCentbins);
Dphi_flat_EP2_inclusive_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP2_leading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP2_subleading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP2_inclusive_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP2_leading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP2_subleading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP3_inclusive_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP3_leading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP3_subleading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP3_inclusive_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP3_leading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP3_subleading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP4_inclusive_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP4_leading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP4_subleading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP4_inclusive_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP4_leading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_flat_EP4_subleading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP2_inclusive_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP2_leading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP2_subleading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP2_inclusive_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP2_leading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP2_subleading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP3_inclusive_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP3_leading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins); 
Dphi_GEN_flat_EP3_subleading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP3_inclusive_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP3_leading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP3_subleading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP4_inclusive_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP4_leading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP4_subleading_minus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP4_inclusive_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP4_leading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
Dphi_GEN_flat_EP4_subleading_plus->GetAxis(1)->Set(bins_JETEP[1],MultCentbins);
hist_reco_jet->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_jet_corr->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_jet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_jet_corr_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_jet->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_jet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_leadjet->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_leadjet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_subljet->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_subljet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_leadjet->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_leadjet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_subljet->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_subljet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_ref_jet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_ref_leadjet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_ref_subljet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_thrdjet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_thrdjet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_ref_thrdjet_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_reco_reco_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_reco_gen_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_gen_reco_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_gen_gen_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_reco_reco_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_reco_gen_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_gen_reco_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jet_from_gen_gen_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_reco_reco_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_reco_gen_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_gen_reco_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_gen_gen_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_reco_reco_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_reco_gen_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_gen_reco_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_lead_jet_from_gen_gen_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_reco_reco_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_reco_gen_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_gen_reco_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_gen_gen_sig->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_reco_reco_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_reco_gen_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_gen_reco_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_subl_jet_from_gen_gen_mix->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_jes_reco_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_jes_reco_fromB_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_leadjes_reco_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_leadjes_reco_fromB_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_subleadjes_reco_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_subleadjes_reco_fromB_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_jetptclos_weighted->GetAxis(4)->Set(bins_jetptclos[4],MultCentbins);
hist_leadjetptclos_weighted->GetAxis(4)->Set(bins_jetptclos[4],MultCentbins);
hist_subljetptclos_weighted->GetAxis(4)->Set(bins_jetptclos[4],MultCentbins);
hist_leadjetptclosremovesome_weighed->GetAxis(4)->Set(bins_jetptclos[4],MultCentbins);
hist_subljetptclosremovesome_weighed->GetAxis(4)->Set(bins_jetptclos[4],MultCentbins);
hist_averjetptclos_weighted->GetAxis(6)->Set(bins_jetptclosave[6],MultCentbins);
hist_averjetptclosremovesome_weighed->GetAxis(6)->Set(bins_jetptclosave[6],MultCentbins);
hist_xjclos_weighted->GetAxis(6)->Set(bins_xjclos[6],MultCentbins);
hist_xjclos_removesome_weighted->GetAxis(6)->Set(bins_xjclos[6],MultCentbins);
hist_reco_lead_reco_subl_quench_mid_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_mid_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_mid_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_fwd_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_bkw_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_fwd_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_fwd_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_bkw_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_reco_subl_quench_bkw_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_mid_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_mid_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_mid_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_fwd_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_bkw_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_fwd_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_fwd_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_bkw_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench_bkw_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_mid_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_mid_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_mid_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_fwd_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_bkw_mid->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_fwd_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_fwd_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_bkw_fwd->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench_bkw_bkw->GetAxis(3)->Set(bins_quenc[3],MultCentbins);

hist_reco_leadEP_quench_plus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_plus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);

hist_reco_leadEP_quench_minus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_mid_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_mid_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_mid_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_fwd_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_bkw_mid->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_fwd_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_fwd_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_bkw_fwd->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus_bkw_bkw->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);

hist_etaDijet_reco->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_etaDijet_CM_reco->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_etaDijet_ref->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_etaDijet_CM_ref->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_etaDijet_gen->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_etaDijet_CM_gen->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_yDijet_CM_reco->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_yDijet_CM_ref->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_yDijet_CM_gen->GetAxis(8)->Set(bins_etaDijet[8],MultCentbins);
hist_injet_reco_track_reco->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_injet_reco_track_gen->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_injet_gen_track_reco->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_injet_gen_track_gen->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inLeadjet_reco_track_reco->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inLeadjet_reco_track_gen->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inLeadjet_gen_track_reco->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inLeadjet_gen_track_gen->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inSubljet_reco_track_reco->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inSubljet_reco_track_gen->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inSubljet_gen_track_reco->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_inSubljet_gen_track_gen->GetAxis(1)->Set(bins_injettrk[1],MultCentbins);
hist_correlation_signal_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_lead_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_lead_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_lead_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subl_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_subl_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_subl_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_lead_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_lead_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_lead_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subl_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_subl_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_subl_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_lead_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_lead_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_lead_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subl_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_subl_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_subl_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_lead_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_lead_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_lead_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subl_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_rotation_subl_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_mixing_subl_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_lead_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_lead_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_lead_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_lead_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_subl_jet_reco_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_subl_jet_reco_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_subl_jet_gen_track_reco->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_correlation_signal_subg0_subl_jet_gen_track_gen->GetAxis(3)->Set(bins_jettrk[3],MultCentbins);
hist_reco_reco_2pcorrelation_signal->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_reco_reco_2pcorrelation_signal_subg0->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_reco_reco_2pcorrelation_signal_subcross->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_reco_reco_2pcorrelation_mixing->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_signal->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_signal_subg0->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_signal_subcross->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_mixing->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
//Extra binning
double Extrabins[extrabinsize-1];
for(int a = 0; a<extra_bins.size();a++){Extrabins[a] = extra_bins[a];}
vzhist->GetAxis(2)->Set(bins_vz[2],Extrabins);
vzhist_weighted->GetAxis(2)->Set(bins_vz[2],Extrabins);
vzhist_jet_weighted->GetAxis(2)->Set(bins_vz[2],Extrabins);
vzhist_dijet_weighted->GetAxis(2)->Set(bins_vz[2],Extrabins);
vxyhist->GetAxis(3)->Set(bins_vxy[3],Extrabins);
vxyhist_weighted->GetAxis(3)->Set(bins_vxy[3],Extrabins);
pthathist->GetAxis(2)->Set(bins_pthat[2],Extrabins);
pthathist_weighted->GetAxis(2)->Set(bins_pthat[2],Extrabins);
EP2_plus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP2_minus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP3_plus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP3_minus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP4_plus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP4_minus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
hist_reco_trk->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_reco_trk_corr->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_reco_trk_weighted->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_gen_trk->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_gen_trk_weighted->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_reco_reco_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_reco_gen_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_gen_reco_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_gen_gen_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_reco_reco_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_reco_gen_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_gen_reco_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_trk_from_gen_gen_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_reco_reco_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_reco_gen_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_gen_reco_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_gen_gen_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_reco_reco_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_reco_gen_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_gen_reco_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_LJ_trk_from_gen_gen_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_reco_reco_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_reco_gen_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_gen_reco_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_gen_gen_sig->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_reco_reco_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_reco_gen_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_gen_reco_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_SLJ_trk_from_gen_gen_mix->GetAxis(4)->Set(bins_trk[4],Extrabins);
Dphi_EP2_flat_trk_minus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_EP2_flat_trk_plus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_EP3_flat_trk_minus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_EP3_flat_trk_plus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_EP4_flat_trk_minus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_EP4_flat_trk_plus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_GEN_EP2_flat_trk_minus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_GEN_EP2_flat_trk_plus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_GEN_EP3_flat_trk_minus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_GEN_EP3_flat_trk_plus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_GEN_EP4_flat_trk_minus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
Dphi_GEN_EP4_flat_trk_plus->GetAxis(3)->Set(bins_TRKEP[3],Extrabins);
NJets->GetAxis(2)->Set(bins_NJETS[2],Extrabins);
trackmaxptinjethisto->GetAxis(2)->Set(bins_trkmax[2],Extrabins);
jettrackmaxptinjethisto->GetAxis(2)->Set(bins_trkmaxjet[2],Extrabins);
histo_jetUE->GetAxis(2)->Set(bins_UE[2],Extrabins);
histo_jetAverageRho->GetAxis(2)->Set(bins_UE[2],Extrabins);
histo_jetcheckcorrection->GetAxis(2)->Set(bins_UE[2],Extrabins);
Dphi_flat_EP2_inclusive_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP2_leading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP2_subleading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP2_inclusive_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP2_leading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP2_subleading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP3_inclusive_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP3_leading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP3_subleading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP3_inclusive_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP3_leading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP3_subleading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP4_inclusive_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP4_leading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP4_subleading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP4_inclusive_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP4_leading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_flat_EP4_subleading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP2_inclusive_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP2_leading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP2_subleading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP2_inclusive_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP2_leading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP2_subleading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP3_inclusive_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP3_leading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins); 
Dphi_GEN_flat_EP3_subleading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP3_inclusive_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP3_leading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP3_subleading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP4_inclusive_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP4_leading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP4_subleading_minus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP4_inclusive_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP4_leading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
Dphi_GEN_flat_EP4_subleading_plus->GetAxis(2)->Set(bins_JETEP[2],Extrabins);
hist_reco_jet->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_jet_corr->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_jet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_jet_corr_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_jet->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_jet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_leadjet->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_leadjet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_subljet->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_subljet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_leadjet->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_leadjet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_subljet->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_subljet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_ref_jet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_ref_leadjet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_ref_subljet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_thrdjet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_thrdjet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_ref_thrdjet_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_reco_reco_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_reco_gen_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_gen_reco_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_gen_gen_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_reco_reco_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_reco_gen_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_gen_reco_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jet_from_gen_gen_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_reco_reco_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_reco_gen_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_gen_reco_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_gen_gen_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_reco_reco_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_reco_gen_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_gen_reco_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_lead_jet_from_gen_gen_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_reco_reco_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_reco_gen_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_gen_reco_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_gen_gen_sig->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_reco_reco_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_reco_gen_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_gen_reco_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_subl_jet_from_gen_gen_mix->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_jes_reco_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_jes_reco_fromB_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_leadjes_reco_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_leadjes_reco_fromB_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_subleadjes_reco_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_subleadjes_reco_fromB_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_jetptclos_weighted->GetAxis(5)->Set(bins_jetptclos[5],Extrabins);
hist_leadjetptclos_weighted->GetAxis(5)->Set(bins_jetptclos[5],Extrabins);
hist_subljetptclos_weighted->GetAxis(5)->Set(bins_jetptclos[5],Extrabins);
hist_leadjetptclosremovesome_weighed->GetAxis(5)->Set(bins_jetptclos[5],Extrabins);
hist_subljetptclosremovesome_weighed->GetAxis(5)->Set(bins_jetptclos[5],Extrabins);
hist_averjetptclos_weighted->GetAxis(7)->Set(bins_jetptclosave[7],Extrabins);
hist_averjetptclosremovesome_weighed->GetAxis(7)->Set(bins_jetptclosave[7],Extrabins);
hist_xjclos_weighted->GetAxis(7)->Set(bins_xjclos[7],Extrabins);
hist_xjclos_removesome_weighted->GetAxis(7)->Set(bins_xjclos[7],Extrabins);
hist_reco_lead_reco_subl_quench_mid_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_mid_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_mid_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_fwd_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_bkw_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_fwd_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_fwd_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_bkw_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_reco_subl_quench_bkw_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_mid_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_mid_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_mid_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_fwd_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_bkw_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_fwd_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_fwd_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_bkw_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench_bkw_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_mid_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_mid_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_mid_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_fwd_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_bkw_mid->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_fwd_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_fwd_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_bkw_fwd->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench_bkw_bkw->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_etaDijet_reco->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_etaDijet_CM_reco->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_etaDijet_ref->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_etaDijet_CM_ref->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_etaDijet_gen->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_etaDijet_CM_gen->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_yDijet_CM_reco->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_yDijet_CM_ref->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_yDijet_CM_gen->GetAxis(10)->Set(bins_etaDijet[10],Extrabins);
hist_injet_reco_track_reco->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_injet_reco_track_gen->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_injet_gen_track_reco->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_injet_gen_track_gen->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inLeadjet_reco_track_reco->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inLeadjet_reco_track_gen->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inLeadjet_gen_track_reco->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inLeadjet_gen_track_gen->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inSubljet_reco_track_reco->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inSubljet_reco_track_gen->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inSubljet_gen_track_reco->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_inSubljet_gen_track_gen->GetAxis(2)->Set(bins_injettrk[2],Extrabins);
hist_correlation_signal_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_lead_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_lead_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_lead_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subl_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_subl_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_subl_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_lead_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_lead_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_lead_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subl_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_subl_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_subl_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_lead_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_lead_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_lead_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subl_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_subl_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_subl_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_lead_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_lead_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_lead_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subl_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_rotation_subl_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_mixing_subl_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_lead_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_lead_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_lead_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_lead_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_subl_jet_reco_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_subl_jet_reco_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_subl_jet_gen_track_reco->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_correlation_signal_subg0_subl_jet_gen_track_gen->GetAxis(4)->Set(bins_jettrk[4],Extrabins);
hist_reco_reco_2pcorrelation_signal->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_reco_reco_2pcorrelation_signal_subg0->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_reco_reco_2pcorrelation_signal_subcross->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_reco_reco_2pcorrelation_mixing->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_signal->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_signal_subg0->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_signal_subcross->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_mixing->GetAxis(4)->Set(bins_2pc[4],Extrabins);

// 3rd jet studies
hist_reco_3rdjet->GetAxis(0)->Set(bins_3rdjet[0],XjBins3rd);
hist_ref_3rdjet->GetAxis(0)->Set(bins_3rdjet[0],XjBins3rd);
hist_gen_3rdjet->GetAxis(0)->Set(bins_3rdjet[0],XjBins3rd);
hist_reco_3rdjet->GetAxis(1)->Set(bins_3rdjet[1],DphiBins3rd);
hist_ref_3rdjet->GetAxis(1)->Set(bins_3rdjet[1],DphiBins3rd);
hist_gen_3rdjet->GetAxis(1)->Set(bins_3rdjet[1],DphiBins3rd);
hist_reco_3rdjet->GetAxis(2)->Set(bins_3rdjet[2],XjBins3rd);
hist_ref_3rdjet->GetAxis(2)->Set(bins_3rdjet[2],XjBins3rd);
hist_gen_3rdjet->GetAxis(2)->Set(bins_3rdjet[2],XjBins3rd);
hist_reco_3rdjet->GetAxis(3)->Set(bins_3rdjet[3],DphiBins3rd);
hist_ref_3rdjet->GetAxis(3)->Set(bins_3rdjet[3],DphiBins3rd);
hist_gen_3rdjet->GetAxis(3)->Set(bins_3rdjet[3],DphiBins3rd);
hist_reco_3rdjet->GetAxis(4)->Set(bins_3rdjet[4],XjBins3rd);
hist_ref_3rdjet->GetAxis(4)->Set(bins_3rdjet[4],XjBins3rd);
hist_gen_3rdjet->GetAxis(4)->Set(bins_3rdjet[4],XjBins3rd);
hist_reco_3rdjet->GetAxis(5)->Set(bins_3rdjet[5],DphiBins3rd);
hist_ref_3rdjet->GetAxis(5)->Set(bins_3rdjet[5],DphiBins3rd);
hist_gen_3rdjet->GetAxis(5)->Set(bins_3rdjet[5],DphiBins3rd);
hist_reco_3rdjet->GetAxis(6)->Set(bins_3rdjet[6],MultCentbins);
hist_ref_3rdjet->GetAxis(6)->Set(bins_3rdjet[6],MultCentbins);
hist_gen_3rdjet->GetAxis(6)->Set(bins_3rdjet[6],MultCentbins);
hist_reco_3rdjet->GetAxis(7)->Set(bins_3rdjet[7],Extrabins);
hist_ref_3rdjet->GetAxis(7)->Set(bins_3rdjet[7],Extrabins);
hist_gen_3rdjet->GetAxis(7)->Set(bins_3rdjet[7],Extrabins);

hist_jetunf_weighted->GetAxis(0)->Set(bins_jetunf[0],PtLSLBins2);
hist_jetunf_weighted->GetAxis(1)->Set(bins_jetunf[1],PtLSLBins2);
hist_jetunf_weighted->GetAxis(2)->Set(bins_jetunf[2],MultCentbins);
hist_jetunf_weighted->GetAxis(3)->Set(bins_jetunf[3],Extrabins);
hist_leadjetunf_weighted->GetAxis(0)->Set(bins_jetunf[0],PtLSLBins2);
hist_leadjetunf_weighted->GetAxis(1)->Set(bins_jetunf[1],PtLSLBins2);
hist_leadjetunf_weighted->GetAxis(2)->Set(bins_jetunf[2],MultCentbins);
hist_leadjetunf_weighted->GetAxis(3)->Set(bins_jetunf[3],Extrabins);
hist_leadjetunf_match_weighted->GetAxis(0)->Set(bins_jetunf[0],PtLSLBins2);
hist_leadjetunf_match_weighted->GetAxis(1)->Set(bins_jetunf[1],PtLSLBins2);
hist_leadjetunf_match_weighted->GetAxis(2)->Set(bins_jetunf[2],MultCentbins);
hist_leadjetunf_match_weighted->GetAxis(3)->Set(bins_jetunf[3],Extrabins);
hist_leadjetunf_swap_weighted->GetAxis(0)->Set(bins_jetunf[0],PtLSLBins2);
hist_leadjetunf_swap_weighted->GetAxis(1)->Set(bins_jetunf[1],PtLSLBins2);
hist_leadjetunf_swap_weighted->GetAxis(2)->Set(bins_jetunf[2],MultCentbins);
hist_leadjetunf_swap_weighted->GetAxis(3)->Set(bins_jetunf[3],Extrabins);
hist_subljetunf_weighted->GetAxis(0)->Set(bins_jetunf[0],PtLSLBins2);
hist_subljetunf_weighted->GetAxis(1)->Set(bins_jetunf[1],PtLSLBins2);
hist_subljetunf_weighted->GetAxis(2)->Set(bins_jetunf[2],MultCentbins);
hist_subljetunf_weighted->GetAxis(3)->Set(bins_jetunf[3],Extrabins);
hist_subljetunf_match_weighted->GetAxis(0)->Set(bins_jetunf[0],PtLSLBins2);
hist_subljetunf_match_weighted->GetAxis(1)->Set(bins_jetunf[1],PtLSLBins2);
hist_subljetunf_match_weighted->GetAxis(2)->Set(bins_jetunf[2],MultCentbins);
hist_subljetunf_match_weighted->GetAxis(3)->Set(bins_jetunf[3],Extrabins);
hist_subljetunf_swap_weighted->GetAxis(0)->Set(bins_jetunf[0],PtLSLBins2);
hist_subljetunf_swap_weighted->GetAxis(1)->Set(bins_jetunf[1],PtLSLBins2);
hist_subljetunf_swap_weighted->GetAxis(2)->Set(bins_jetunf[2],MultCentbins);
hist_subljetunf_swap_weighted->GetAxis(3)->Set(bins_jetunf[3],Extrabins);
hist_xjunf_weighted->GetAxis(0)->Set(bins_xjunf[0],XjBins);
hist_xjunf_weighted->GetAxis(1)->Set(bins_xjunf[1],XjBins);
hist_xjunf_weighted->GetAxis(2)->Set(bins_xjunf[2],MultCentbins);
hist_xjunf_weighted->GetAxis(3)->Set(bins_xjunf[3],Extrabins);
hist_xjunf_match_weighted->GetAxis(0)->Set(bins_xjunf[0],XjBins);
hist_xjunf_match_weighted->GetAxis(1)->Set(bins_xjunf[1],XjBins);
hist_xjunf_match_weighted->GetAxis(2)->Set(bins_xjunf[2],MultCentbins);
hist_xjunf_match_weighted->GetAxis(3)->Set(bins_xjunf[3],Extrabins);
hist_xjunf_swap_weighted->GetAxis(0)->Set(bins_xjunf[0],XjBins);
hist_xjunf_swap_weighted->GetAxis(1)->Set(bins_xjunf[1],XjBins);
hist_xjunf_swap_weighted->GetAxis(2)->Set(bins_xjunf[2],MultCentbins);
hist_xjunf_swap_weighted->GetAxis(3)->Set(bins_xjunf[3],Extrabins);

hist_jetunf_weighted_4D->GetAxis(0)->Set(bins_jetunf4D[0],PtLSLBins2);
hist_jetunf_weighted_4D->GetAxis(1)->Set(bins_jetunf4D[1],PtLSLBins2);
hist_jetunf_weighted_4D->GetAxis(2)->Set(bins_jetunf4D[2],PtLSLBins2);
hist_jetunf_weighted_4D->GetAxis(3)->Set(bins_jetunf4D[3],PtLSLBins2);
hist_jetunf_weighted_4D->GetAxis(4)->Set(bins_jetunf4D[4],MultCentbins);
hist_jetunf_weighted_4D->GetAxis(5)->Set(bins_jetunf4D[5],Extrabins);
hist_jetunf_match_weighted_4D->GetAxis(0)->Set(bins_jetunf4D[0],PtLSLBins2);
hist_jetunf_match_weighted_4D->GetAxis(1)->Set(bins_jetunf4D[1],PtLSLBins2);
hist_jetunf_match_weighted_4D->GetAxis(2)->Set(bins_jetunf4D[2],PtLSLBins2);
hist_jetunf_match_weighted_4D->GetAxis(3)->Set(bins_jetunf4D[3],PtLSLBins2);
hist_jetunf_match_weighted_4D->GetAxis(4)->Set(bins_jetunf4D[4],MultCentbins);
hist_jetunf_match_weighted_4D->GetAxis(5)->Set(bins_jetunf4D[5],Extrabins);
hist_jetunf_swap_weighted_4D->GetAxis(0)->Set(bins_jetunf4D[0],PtLSLBins2);
hist_jetunf_swap_weighted_4D->GetAxis(1)->Set(bins_jetunf4D[1],PtLSLBins2);
hist_jetunf_swap_weighted_4D->GetAxis(2)->Set(bins_jetunf4D[2],PtLSLBins2);
hist_jetunf_swap_weighted_4D->GetAxis(3)->Set(bins_jetunf4D[3],PtLSLBins2);
hist_jetunf_swap_weighted_4D->GetAxis(4)->Set(bins_jetunf4D[4],MultCentbins);
hist_jetunf_swap_weighted_4D->GetAxis(5)->Set(bins_jetunf4D[5],Extrabins);
hist_jetunf_match_weighted_sym_4D->GetAxis(0)->Set(bins_jetunf4D[0],PtLSLBins2);
hist_jetunf_match_weighted_sym_4D->GetAxis(1)->Set(bins_jetunf4D[1],PtLSLBins2);
hist_jetunf_match_weighted_sym_4D->GetAxis(2)->Set(bins_jetunf4D[2],PtLSLBins2);
hist_jetunf_match_weighted_sym_4D->GetAxis(3)->Set(bins_jetunf4D[3],PtLSLBins2);
hist_jetunf_match_weighted_sym_4D->GetAxis(4)->Set(bins_jetunf4D[4],MultCentbins);
hist_jetunf_match_weighted_sym_4D->GetAxis(5)->Set(bins_jetunf4D[5],Extrabins);
hist_jetunf_match_weighted_symhalf_4D->GetAxis(0)->Set(bins_jetunf4D[0],PtLSLBins2);
hist_jetunf_match_weighted_symhalf_4D->GetAxis(1)->Set(bins_jetunf4D[1],PtLSLBins2);
hist_jetunf_match_weighted_symhalf_4D->GetAxis(2)->Set(bins_jetunf4D[2],PtLSLBins2);
hist_jetunf_match_weighted_symhalf_4D->GetAxis(3)->Set(bins_jetunf4D[3],PtLSLBins2);
hist_jetunf_match_weighted_symhalf_4D->GetAxis(4)->Set(bins_jetunf4D[4],MultCentbins);
hist_jetunf_match_weighted_symhalf_4D->GetAxis(5)->Set(bins_jetunf4D[5],Extrabins);

hist_leadjetunf_gensmear->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_subljetunf_gensmear->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_leadjetunf_gensmear->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_subljetunf_gensmear->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_leadjetunf_gensmear->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_subljetunf_gensmear->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_leadjetunf_recosmear->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_subljetunf_recosmear->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_leadjetunf_recosmear->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_subljetunf_recosmear->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_leadjetunf_recosmear->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_subljetunf_recosmear->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);


hist_leadjetunf_gensmear_fromInclJet->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_leadjetunf_gensmear_fromInclJet->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_subljetunf_gensmear_fromInclJet->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_subljetunf_gensmear_fromInclJet->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_leadjetunf_recosmear_fromInclJet->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_leadjetunf_recosmear_fromInclJet->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_subljetunf_recosmear_fromInclJet->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_subljetunf_recosmear_fromInclJet->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);


hist_leadjetunf_gensmear_4D->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_subljetunf_gensmear_4D->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_leadjetunf_gensmear_4D->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_subljetunf_gensmear_4D->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_leadjetunf_gensmear_4D->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_subljetunf_gensmear_4D->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_leadjetunf_recosmear_4D->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_subljetunf_recosmear_4D->GetAxis(0)->Set(bins_jetunf_smear[0],PtLSLBins2);
hist_leadjetunf_recosmear_4D->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_subljetunf_recosmear_4D->GetAxis(1)->Set(bins_jetunf_smear[1],MultCentbins);
hist_leadjetunf_recosmear_4D->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);
hist_subljetunf_recosmear_4D->GetAxis(2)->Set(bins_jetunf_smear[2],Extrabins);

hist_xjunf_gensmear->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_gensmear->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_gensmear->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);
hist_xjunf_recosmear->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_recosmear->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_recosmear->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);
hist_xjunf_gensmear_fromLSL->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_gensmear_fromLSL->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_gensmear_fromLSL->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);
hist_xjunf_recosmear_fromLSL->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_recosmear_fromLSL->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_recosmear_fromLSL->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);

hist_xjunf_gensmear_fromLSL_4D->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_gensmear_fromLSL_4D->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_gensmear_fromLSL_4D->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);
hist_xjunf_recosmear_fromLSL_4D->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_recosmear_fromLSL_4D->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_recosmear_fromLSL_4D->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);

hist_xjunf_gensmear_fromLSL_InclJet->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_gensmear_fromLSL_InclJet->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_gensmear_fromLSL_InclJet->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);
hist_xjunf_recosmear_fromLSL_InclJet->GetAxis(0)->Set(bins_xjunf_smear[0],XjBins);
hist_xjunf_recosmear_fromLSL_InclJet->GetAxis(1)->Set(bins_xjunf_smear[1],MultCentbins);
hist_xjunf_recosmear_fromLSL_InclJet->GetAxis(2)->Set(bins_xjunf_smear[2],Extrabins);


//Sumw2 starts here!
Nevents->Sumw2();
Nev_recoreco->Sumw2();
Nev_recoreco_lead->Sumw2();
Nev_recoreco_subl->Sumw2();
Nev_recogen->Sumw2();
Nev_recogen_lead->Sumw2();
Nev_recogen_subl->Sumw2();
Nev_genreco->Sumw2();
Nev_genreco_lead->Sumw2();
Nev_genreco_subl->Sumw2();
Nev_gengen->Sumw2();
Nev_gengen_lead->Sumw2();
Nev_gengen_subl->Sumw2();
multiplicity->Sumw2();
multiplicity_weighted->Sumw2();
multiplicity_withonejet_weighted->Sumw2();
multiplicity_withdijets_weighted->Sumw2();
multiplicity_midmid_weighted->Sumw2();
multiplicity_midfwd_weighted->Sumw2();
multiplicity_midbkw_weighted->Sumw2();
multiplicity_fwdmid_weighted->Sumw2();
multiplicity_fwdfwd_weighted->Sumw2();
multiplicity_fwdbkw_weighted->Sumw2();
multiplicity_bkwmid_weighted->Sumw2();
multiplicity_bkwfwd_weighted->Sumw2();
multiplicity_bkwbkw_weighted->Sumw2();
multiplicity_weighted_at1->Sumw2();
multiplicity_withonejet_weighted_at1->Sumw2();
multiplicity_withdijets_weighted_at1->Sumw2();
reco_mult->Sumw2();
reco_mult_weighted->Sumw2();
reco_mult_withonejet_weighted->Sumw2();
reco_mult_withdijets_weighted->Sumw2();
gen_mult->Sumw2();
gen_mult_weighted->Sumw2();
gen_mult_withonejet_weighted->Sumw2();
gen_mult_withdijets_weighted->Sumw2();
hfhist->Sumw2();
hfhist_weighted->Sumw2();
hfhistEta4->Sumw2();
hfhistEta4_weighted->Sumw2();
hfhist_onejet_weighted->Sumw2();
hfhistEta4_onejet_weighted->Sumw2();
hfhist_dijet_weighted->Sumw2();
hfhistEta4_dijet_weighted->Sumw2();
zdchist->Sumw2();
zdchist_weighted->Sumw2();
zdchist_onejet_weighted->Sumw2();
zdchist_dijet_weighted->Sumw2();
hfhistSum_weighted->Sumw2();
hfhistEta4Sum_weighted->Sumw2();
hfhistSum_onejet_weighted->Sumw2();
hfhistEta4Sum_onejet_weighted->Sumw2();
hfhistSum_dijet_weighted->Sumw2();
hfhistEta4Sum_dijet_weighted->Sumw2();
vzhist->Sumw2();
vzhist_weighted->Sumw2();
vzhist_jet_weighted->Sumw2();
vzhist_dijet_weighted->Sumw2();
vxyhist->Sumw2();
vxyhist_weighted->Sumw2();
pthathist->Sumw2();
pthathist_weighted->Sumw2();
EP2_plus_flat->Sumw2();
EP2_minus_flat->Sumw2();
EP3_plus_flat->Sumw2();
EP3_minus_flat->Sumw2();
EP4_plus_flat->Sumw2();
EP4_minus_flat->Sumw2();
hist_reco_trk->Sumw2();
hist_reco_trk_corr->Sumw2();
hist_reco_trk_weighted->Sumw2();
hist_gen_trk->Sumw2();
hist_gen_trk_weighted->Sumw2();
hist_trk_from_reco_reco_sig->Sumw2();
hist_trk_from_reco_gen_sig->Sumw2();
hist_trk_from_gen_reco_sig->Sumw2();
hist_trk_from_gen_gen_sig->Sumw2();
hist_trk_from_reco_reco_mix->Sumw2();
hist_trk_from_reco_gen_mix->Sumw2();
hist_trk_from_gen_reco_mix->Sumw2();
hist_trk_from_gen_gen_mix->Sumw2();
hist_LJ_trk_from_reco_reco_sig->Sumw2();
hist_LJ_trk_from_reco_gen_sig->Sumw2();
hist_LJ_trk_from_gen_reco_sig->Sumw2();
hist_LJ_trk_from_gen_gen_sig->Sumw2();
hist_LJ_trk_from_reco_reco_mix->Sumw2();
hist_LJ_trk_from_reco_gen_mix->Sumw2();
hist_LJ_trk_from_gen_reco_mix->Sumw2();
hist_LJ_trk_from_gen_gen_mix->Sumw2();
hist_SLJ_trk_from_reco_reco_sig->Sumw2();
hist_SLJ_trk_from_reco_gen_sig->Sumw2();
hist_SLJ_trk_from_gen_reco_sig->Sumw2();
hist_SLJ_trk_from_gen_gen_sig->Sumw2();
hist_SLJ_trk_from_reco_reco_mix->Sumw2();
hist_SLJ_trk_from_reco_gen_mix->Sumw2();
hist_SLJ_trk_from_gen_reco_mix->Sumw2();
hist_SLJ_trk_from_gen_gen_mix->Sumw2();
Dphi_EP2_flat_trk_minus->Sumw2();
Dphi_EP2_flat_trk_plus->Sumw2();
Dphi_EP3_flat_trk_minus->Sumw2();
Dphi_EP3_flat_trk_plus->Sumw2();
Dphi_EP4_flat_trk_minus->Sumw2();
Dphi_EP4_flat_trk_plus->Sumw2();
Dphi_GEN_EP2_flat_trk_minus->Sumw2();
Dphi_GEN_EP2_flat_trk_plus->Sumw2();
Dphi_GEN_EP3_flat_trk_minus->Sumw2();
Dphi_GEN_EP3_flat_trk_plus->Sumw2();
Dphi_GEN_EP4_flat_trk_minus->Sumw2();
Dphi_GEN_EP4_flat_trk_plus->Sumw2();
NJets->Sumw2();
trackmaxptinjethisto->Sumw2();
jettrackmaxptinjethisto->Sumw2();
histo_jetUE->Sumw2();
histo_jetAverageRho->Sumw2();
histo_jetcheckcorrection->Sumw2();
Dphi_flat_EP2_inclusive_minus->Sumw2();
Dphi_flat_EP2_leading_minus->Sumw2();
Dphi_flat_EP2_subleading_minus->Sumw2();
Dphi_flat_EP2_inclusive_plus->Sumw2();
Dphi_flat_EP2_leading_plus->Sumw2();
Dphi_flat_EP2_subleading_plus->Sumw2();
Dphi_flat_EP3_inclusive_minus->Sumw2();
Dphi_flat_EP3_leading_minus->Sumw2();
Dphi_flat_EP3_subleading_minus->Sumw2();
Dphi_flat_EP3_inclusive_plus->Sumw2();
Dphi_flat_EP3_leading_plus->Sumw2();
Dphi_flat_EP3_subleading_plus->Sumw2();
Dphi_flat_EP4_inclusive_minus->Sumw2();
Dphi_flat_EP4_leading_minus->Sumw2();
Dphi_flat_EP4_subleading_minus->Sumw2();
Dphi_flat_EP4_inclusive_plus->Sumw2();
Dphi_flat_EP4_leading_plus->Sumw2();
Dphi_flat_EP4_subleading_plus->Sumw2();
Dphi_GEN_flat_EP2_inclusive_minus->Sumw2();
Dphi_GEN_flat_EP2_leading_minus->Sumw2();
Dphi_GEN_flat_EP2_subleading_minus->Sumw2();
Dphi_GEN_flat_EP2_inclusive_plus->Sumw2();
Dphi_GEN_flat_EP2_leading_plus->Sumw2();
Dphi_GEN_flat_EP2_subleading_plus->Sumw2();
Dphi_GEN_flat_EP3_inclusive_minus->Sumw2();
Dphi_GEN_flat_EP3_leading_minus->Sumw2(); 
Dphi_GEN_flat_EP3_subleading_minus->Sumw2();
Dphi_GEN_flat_EP3_inclusive_plus->Sumw2();
Dphi_GEN_flat_EP3_leading_plus->Sumw2();
Dphi_GEN_flat_EP3_subleading_plus->Sumw2();
Dphi_GEN_flat_EP4_inclusive_minus->Sumw2();
Dphi_GEN_flat_EP4_leading_minus->Sumw2();
Dphi_GEN_flat_EP4_subleading_minus->Sumw2();
Dphi_GEN_flat_EP4_inclusive_plus->Sumw2();
Dphi_GEN_flat_EP4_leading_plus->Sumw2();
Dphi_GEN_flat_EP4_subleading_plus->Sumw2();
hist_reco_jet->Sumw2();
hist_reco_jet_corr->Sumw2();
hist_reco_jet_weighted->Sumw2();
hist_reco_jet_corr_weighted->Sumw2();
hist_gen_jet->Sumw2();
hist_gen_jet_weighted->Sumw2();
hist_reco_leadjet->Sumw2();
hist_reco_leadjet_weighted->Sumw2();
hist_reco_subljet->Sumw2();
hist_reco_subljet_weighted->Sumw2();
hist_gen_leadjet->Sumw2();
hist_gen_leadjet_weighted->Sumw2();
hist_gen_subljet->Sumw2();
hist_gen_subljet_weighted->Sumw2();
hist_ref_jet_weighted->Sumw2();
hist_ref_leadjet_weighted->Sumw2();
hist_ref_subljet_weighted->Sumw2();
hist_reco_thrdjet_weighted->Sumw2();
hist_gen_thrdjet_weighted->Sumw2();
hist_ref_thrdjet_weighted->Sumw2();
hist_jes_reco_weighted->Sumw2();
hist_jes_reco_fromB_weighted->Sumw2();
hist_leadjes_reco_weighted->Sumw2();
hist_leadjes_reco_fromB_weighted->Sumw2();
hist_subleadjes_reco_weighted->Sumw2();
hist_subleadjes_reco_fromB_weighted->Sumw2();
hist_jetptclos_weighted->Sumw2();
hist_leadjetptclos_weighted->Sumw2();
hist_subljetptclos_weighted->Sumw2();
hist_averjetptclos_weighted->Sumw2();
hist_leadjetptclosremovesome_weighed->Sumw2();
hist_subljetptclosremovesome_weighed->Sumw2();
hist_averjetptclosremovesome_weighed->Sumw2();
hist_xjclos_weighted->Sumw2();
hist_xjclos_removesome_weighted->Sumw2();
hist_jetunf_weighted->Sumw2();
hist_leadjetunf_weighted->Sumw2();
hist_subljetunf_weighted->Sumw2();
hist_xjunf_weighted->Sumw2();
hist_leadjetunf_match_weighted->Sumw2();
hist_subljetunf_match_weighted->Sumw2();
hist_leadjetunf_swap_weighted->Sumw2();
hist_subljetunf_swap_weighted->Sumw2();
hist_xjunf_match_weighted->Sumw2();
hist_xjunf_swap_weighted->Sumw2();
hist_leadjetunf_gensmear->Sumw2();
hist_leadjetunf_recosmear->Sumw2();
hist_subljetunf_gensmear->Sumw2();
hist_subljetunf_recosmear->Sumw2();
hist_xjunf_gensmear->Sumw2();
hist_xjunf_recosmear->Sumw2();
hist_xjunf_gensmear_fromLSL->Sumw2();
hist_xjunf_recosmear_fromLSL->Sumw2();

hist_leadjetunf_gensmear_fromInclJet->Sumw2();
hist_leadjetunf_recosmear_fromInclJet->Sumw2();
hist_subljetunf_gensmear_fromInclJet->Sumw2();
hist_subljetunf_recosmear_fromInclJet->Sumw2();
hist_xjunf_gensmear_fromLSL_InclJet->Sumw2();
hist_xjunf_recosmear_fromLSL_InclJet->Sumw2();

hist_xjunf_gensmear_fromLSL_4D->Sumw2();
hist_xjunf_recosmear_fromLSL_4D->Sumw2();
hist_leadjetunf_gensmear_4D->Sumw2();
hist_subljetunf_gensmear_4D->Sumw2();
hist_leadjetunf_recosmear_4D->Sumw2();
hist_subljetunf_recosmear_4D->Sumw2();
hist_jetunf_weighted_4D->Sumw2();
hist_jetunf_match_weighted_4D->Sumw2();
hist_jetunf_match_weighted_sym_4D->Sumw2();
hist_jetunf_match_weighted_symhalf_4D->Sumw2();
hist_jetunf_swap_weighted_4D->Sumw2();

hist_jet_from_reco_reco_sig->Sumw2();
hist_jet_from_reco_gen_sig->Sumw2();
hist_jet_from_gen_reco_sig->Sumw2();
hist_jet_from_gen_gen_sig->Sumw2();
hist_jet_from_reco_reco_mix->Sumw2();
hist_jet_from_reco_gen_mix->Sumw2();
hist_jet_from_gen_reco_mix->Sumw2();
hist_jet_from_gen_gen_mix->Sumw2();
hist_lead_jet_from_reco_reco_sig->Sumw2();
hist_lead_jet_from_reco_gen_sig->Sumw2();
hist_lead_jet_from_gen_reco_sig->Sumw2();
hist_lead_jet_from_gen_gen_sig->Sumw2();
hist_lead_jet_from_reco_reco_mix->Sumw2();
hist_lead_jet_from_reco_gen_mix->Sumw2();
hist_lead_jet_from_gen_reco_mix->Sumw2();
hist_lead_jet_from_gen_gen_mix->Sumw2();
hist_subl_jet_from_reco_reco_sig->Sumw2();
hist_subl_jet_from_reco_gen_sig->Sumw2();
hist_subl_jet_from_gen_reco_sig->Sumw2();
hist_subl_jet_from_gen_gen_sig->Sumw2();
hist_subl_jet_from_reco_reco_mix->Sumw2();
hist_subl_jet_from_reco_gen_mix->Sumw2();
hist_subl_jet_from_gen_reco_mix->Sumw2();
hist_subl_jet_from_gen_gen_mix->Sumw2();
hist_reco_lead_reco_subl_quench_mid_mid->Sumw2();
hist_reco_lead_reco_subl_quench_mid_fwd->Sumw2();
hist_reco_lead_reco_subl_quench_mid_bkw->Sumw2();
hist_reco_lead_reco_subl_quench_fwd_mid->Sumw2();
hist_reco_lead_reco_subl_quench_bkw_mid->Sumw2();
hist_reco_lead_reco_subl_quench_fwd_fwd->Sumw2();
hist_reco_lead_reco_subl_quench_fwd_bkw->Sumw2();
hist_reco_lead_reco_subl_quench_bkw_fwd->Sumw2();
hist_reco_lead_reco_subl_quench_bkw_bkw->Sumw2();
hist_gen_lead_gen_subl_quench_mid_mid->Sumw2();
hist_gen_lead_gen_subl_quench_mid_fwd->Sumw2();
hist_gen_lead_gen_subl_quench_mid_bkw->Sumw2();
hist_gen_lead_gen_subl_quench_fwd_mid->Sumw2();
hist_gen_lead_gen_subl_quench_bkw_mid->Sumw2();
hist_gen_lead_gen_subl_quench_fwd_fwd->Sumw2();
hist_gen_lead_gen_subl_quench_fwd_bkw->Sumw2();
hist_gen_lead_gen_subl_quench_bkw_fwd->Sumw2();
hist_gen_lead_gen_subl_quench_bkw_bkw->Sumw2();
hist_ref_lead_ref_subl_quench_mid_mid->Sumw2();
hist_ref_lead_ref_subl_quench_mid_fwd->Sumw2();
hist_ref_lead_ref_subl_quench_mid_bkw->Sumw2();
hist_ref_lead_ref_subl_quench_fwd_mid->Sumw2();
hist_ref_lead_ref_subl_quench_bkw_mid->Sumw2();
hist_ref_lead_ref_subl_quench_fwd_fwd->Sumw2();
hist_ref_lead_ref_subl_quench_fwd_bkw->Sumw2();
hist_ref_lead_ref_subl_quench_bkw_fwd->Sumw2();
hist_ref_lead_ref_subl_quench_bkw_bkw->Sumw2();
hist_reco_leadEP_quench_plus_mid_mid->Sumw2();
hist_reco_leadEP_quench_plus_mid_fwd->Sumw2();
hist_reco_leadEP_quench_plus_mid_bkw->Sumw2();
hist_reco_leadEP_quench_plus_fwd_mid->Sumw2();
hist_reco_leadEP_quench_plus_bkw_mid->Sumw2();
hist_reco_leadEP_quench_plus_fwd_fwd->Sumw2();
hist_reco_leadEP_quench_plus_fwd_bkw->Sumw2();
hist_reco_leadEP_quench_plus_bkw_fwd->Sumw2();
hist_reco_leadEP_quench_plus_bkw_bkw->Sumw2();
hist_reco_sublEP_quench_plus_mid_mid->Sumw2();
hist_reco_sublEP_quench_plus_mid_fwd->Sumw2();
hist_reco_sublEP_quench_plus_mid_bkw->Sumw2();
hist_reco_sublEP_quench_plus_fwd_mid->Sumw2();
hist_reco_sublEP_quench_plus_bkw_mid->Sumw2();
hist_reco_sublEP_quench_plus_fwd_fwd->Sumw2();
hist_reco_sublEP_quench_plus_fwd_bkw->Sumw2();
hist_reco_sublEP_quench_plus_bkw_fwd->Sumw2();
hist_reco_sublEP_quench_plus_bkw_bkw->Sumw2();
hist_ref_leadEP_quench_plus_mid_mid->Sumw2();
hist_ref_leadEP_quench_plus_mid_fwd->Sumw2();
hist_ref_leadEP_quench_plus_mid_bkw->Sumw2();
hist_ref_leadEP_quench_plus_fwd_mid->Sumw2();
hist_ref_leadEP_quench_plus_bkw_mid->Sumw2();
hist_ref_leadEP_quench_plus_fwd_fwd->Sumw2();
hist_ref_leadEP_quench_plus_fwd_bkw->Sumw2();
hist_ref_leadEP_quench_plus_bkw_fwd->Sumw2();
hist_ref_leadEP_quench_plus_bkw_bkw->Sumw2();
hist_ref_sublEP_quench_plus_mid_mid->Sumw2();
hist_ref_sublEP_quench_plus_mid_fwd->Sumw2();
hist_ref_sublEP_quench_plus_mid_bkw->Sumw2();
hist_ref_sublEP_quench_plus_fwd_mid->Sumw2();
hist_ref_sublEP_quench_plus_bkw_mid->Sumw2();
hist_ref_sublEP_quench_plus_fwd_fwd->Sumw2();
hist_ref_sublEP_quench_plus_fwd_bkw->Sumw2();
hist_ref_sublEP_quench_plus_bkw_fwd->Sumw2();
hist_ref_sublEP_quench_plus_bkw_bkw->Sumw2();
hist_gen_leadEP_quench_plus_mid_mid->Sumw2();
hist_gen_leadEP_quench_plus_mid_fwd->Sumw2();
hist_gen_leadEP_quench_plus_mid_bkw->Sumw2();
hist_gen_leadEP_quench_plus_fwd_mid->Sumw2();
hist_gen_leadEP_quench_plus_bkw_mid->Sumw2();
hist_gen_leadEP_quench_plus_fwd_fwd->Sumw2();
hist_gen_leadEP_quench_plus_fwd_bkw->Sumw2();
hist_gen_leadEP_quench_plus_bkw_fwd->Sumw2();
hist_gen_leadEP_quench_plus_bkw_bkw->Sumw2();
hist_gen_sublEP_quench_plus_mid_mid->Sumw2();
hist_gen_sublEP_quench_plus_mid_fwd->Sumw2();
hist_gen_sublEP_quench_plus_mid_bkw->Sumw2();
hist_gen_sublEP_quench_plus_fwd_mid->Sumw2();
hist_gen_sublEP_quench_plus_bkw_mid->Sumw2();
hist_gen_sublEP_quench_plus_fwd_fwd->Sumw2();
hist_gen_sublEP_quench_plus_fwd_bkw->Sumw2();
hist_gen_sublEP_quench_plus_bkw_fwd->Sumw2();
hist_gen_sublEP_quench_plus_bkw_bkw->Sumw2();
hist_reco_leadEP_quench_minus_mid_mid->Sumw2();
hist_reco_leadEP_quench_minus_mid_fwd->Sumw2();
hist_reco_leadEP_quench_minus_mid_bkw->Sumw2();
hist_reco_leadEP_quench_minus_fwd_mid->Sumw2();
hist_reco_leadEP_quench_minus_bkw_mid->Sumw2();
hist_reco_leadEP_quench_minus_fwd_fwd->Sumw2();
hist_reco_leadEP_quench_minus_fwd_bkw->Sumw2();
hist_reco_leadEP_quench_minus_bkw_fwd->Sumw2();
hist_reco_leadEP_quench_minus_bkw_bkw->Sumw2();
hist_reco_sublEP_quench_minus_mid_mid->Sumw2();
hist_reco_sublEP_quench_minus_mid_fwd->Sumw2();
hist_reco_sublEP_quench_minus_mid_bkw->Sumw2();
hist_reco_sublEP_quench_minus_fwd_mid->Sumw2();
hist_reco_sublEP_quench_minus_bkw_mid->Sumw2();
hist_reco_sublEP_quench_minus_fwd_fwd->Sumw2();
hist_reco_sublEP_quench_minus_fwd_bkw->Sumw2();
hist_reco_sublEP_quench_minus_bkw_fwd->Sumw2();
hist_reco_sublEP_quench_minus_bkw_bkw->Sumw2();
hist_ref_leadEP_quench_minus_mid_mid->Sumw2();
hist_ref_leadEP_quench_minus_mid_fwd->Sumw2();
hist_ref_leadEP_quench_minus_mid_bkw->Sumw2();
hist_ref_leadEP_quench_minus_fwd_mid->Sumw2();
hist_ref_leadEP_quench_minus_bkw_mid->Sumw2();
hist_ref_leadEP_quench_minus_fwd_fwd->Sumw2();
hist_ref_leadEP_quench_minus_fwd_bkw->Sumw2();
hist_ref_leadEP_quench_minus_bkw_fwd->Sumw2();
hist_ref_leadEP_quench_minus_bkw_bkw->Sumw2();
hist_ref_sublEP_quench_minus_mid_mid->Sumw2();
hist_ref_sublEP_quench_minus_mid_fwd->Sumw2();
hist_ref_sublEP_quench_minus_mid_bkw->Sumw2();
hist_ref_sublEP_quench_minus_fwd_mid->Sumw2();
hist_ref_sublEP_quench_minus_bkw_mid->Sumw2();
hist_ref_sublEP_quench_minus_fwd_fwd->Sumw2();
hist_ref_sublEP_quench_minus_fwd_bkw->Sumw2();
hist_ref_sublEP_quench_minus_bkw_fwd->Sumw2();
hist_ref_sublEP_quench_minus_bkw_bkw->Sumw2();
hist_gen_leadEP_quench_minus_mid_mid->Sumw2();
hist_gen_leadEP_quench_minus_mid_fwd->Sumw2();
hist_gen_leadEP_quench_minus_mid_bkw->Sumw2();
hist_gen_leadEP_quench_minus_fwd_mid->Sumw2();
hist_gen_leadEP_quench_minus_bkw_mid->Sumw2();
hist_gen_leadEP_quench_minus_fwd_fwd->Sumw2();
hist_gen_leadEP_quench_minus_fwd_bkw->Sumw2();
hist_gen_leadEP_quench_minus_bkw_fwd->Sumw2();
hist_gen_leadEP_quench_minus_bkw_bkw->Sumw2();
hist_gen_sublEP_quench_minus_mid_mid->Sumw2();
hist_gen_sublEP_quench_minus_mid_fwd->Sumw2();
hist_gen_sublEP_quench_minus_mid_bkw->Sumw2();
hist_gen_sublEP_quench_minus_fwd_mid->Sumw2();
hist_gen_sublEP_quench_minus_bkw_mid->Sumw2();
hist_gen_sublEP_quench_minus_fwd_fwd->Sumw2();
hist_gen_sublEP_quench_minus_fwd_bkw->Sumw2();
hist_gen_sublEP_quench_minus_bkw_fwd->Sumw2();
hist_gen_sublEP_quench_minus_bkw_bkw->Sumw2();
hist_etaDijet_reco->Sumw2();
hist_etaDijet_CM_reco->Sumw2();
hist_etaDijet_ref->Sumw2();
hist_etaDijet_CM_ref->Sumw2();
hist_etaDijet_gen->Sumw2();
hist_etaDijet_CM_gen->Sumw2();
hist_yDijet_CM_reco->Sumw2();
hist_yDijet_CM_ref->Sumw2();
hist_yDijet_CM_gen->Sumw2();
hist_injet_reco_track_reco->Sumw2();
hist_injet_reco_track_gen->Sumw2();
hist_injet_gen_track_reco->Sumw2();
hist_injet_gen_track_gen->Sumw2();
hist_inLeadjet_reco_track_reco->Sumw2();
hist_inLeadjet_reco_track_gen->Sumw2();
hist_inLeadjet_gen_track_reco->Sumw2();
hist_inLeadjet_gen_track_gen->Sumw2();
hist_inSubljet_reco_track_reco->Sumw2();
hist_inSubljet_reco_track_gen->Sumw2();
hist_inSubljet_gen_track_reco->Sumw2();
hist_inSubljet_gen_track_gen->Sumw2();
hist_correlation_signal_jet_reco_track_reco->Sumw2();
hist_correlation_rotation_jet_reco_track_reco->Sumw2();
hist_correlation_mixing_jet_reco_track_reco->Sumw2();
hist_correlation_signal_lead_jet_reco_track_reco->Sumw2();
hist_correlation_rotation_lead_jet_reco_track_reco->Sumw2();
hist_correlation_mixing_lead_jet_reco_track_reco->Sumw2();
hist_correlation_signal_subl_jet_reco_track_reco->Sumw2();
hist_correlation_rotation_subl_jet_reco_track_reco->Sumw2();
hist_correlation_mixing_subl_jet_reco_track_reco->Sumw2();
hist_correlation_signal_jet_reco_track_gen->Sumw2();
hist_correlation_rotation_jet_reco_track_gen->Sumw2();
hist_correlation_mixing_jet_reco_track_gen->Sumw2();
hist_correlation_signal_lead_jet_reco_track_gen->Sumw2();
hist_correlation_rotation_lead_jet_reco_track_gen->Sumw2();
hist_correlation_mixing_lead_jet_reco_track_gen->Sumw2();
hist_correlation_signal_subl_jet_reco_track_gen->Sumw2();
hist_correlation_rotation_subl_jet_reco_track_gen->Sumw2();
hist_correlation_mixing_subl_jet_reco_track_gen->Sumw2();
hist_correlation_signal_jet_gen_track_reco->Sumw2();
hist_correlation_rotation_jet_gen_track_reco->Sumw2();
hist_correlation_mixing_jet_gen_track_reco->Sumw2();
hist_correlation_signal_lead_jet_gen_track_reco->Sumw2();
hist_correlation_rotation_lead_jet_gen_track_reco->Sumw2();
hist_correlation_mixing_lead_jet_gen_track_reco->Sumw2();
hist_correlation_signal_subl_jet_gen_track_reco->Sumw2();
hist_correlation_rotation_subl_jet_gen_track_reco->Sumw2();
hist_correlation_mixing_subl_jet_gen_track_reco->Sumw2();
hist_correlation_signal_jet_gen_track_gen->Sumw2();
hist_correlation_rotation_jet_gen_track_gen->Sumw2();
hist_correlation_mixing_jet_gen_track_gen->Sumw2();
hist_correlation_signal_lead_jet_gen_track_gen->Sumw2();
hist_correlation_rotation_lead_jet_gen_track_gen->Sumw2();
hist_correlation_mixing_lead_jet_gen_track_gen->Sumw2();
hist_correlation_signal_subl_jet_gen_track_gen->Sumw2();
hist_correlation_rotation_subl_jet_gen_track_gen->Sumw2();
hist_correlation_mixing_subl_jet_gen_track_gen->Sumw2();
hist_correlation_signal_subg0_jet_reco_track_reco->Sumw2();
hist_correlation_signal_subg0_jet_reco_track_gen->Sumw2();
hist_correlation_signal_subg0_jet_gen_track_reco->Sumw2();
hist_correlation_signal_subg0_jet_gen_track_gen->Sumw2();
hist_correlation_signal_subg0_lead_jet_reco_track_reco->Sumw2();
hist_correlation_signal_subg0_lead_jet_reco_track_gen->Sumw2();
hist_correlation_signal_subg0_lead_jet_gen_track_reco->Sumw2();
hist_correlation_signal_subg0_lead_jet_gen_track_gen->Sumw2();
hist_correlation_signal_subg0_subl_jet_reco_track_reco->Sumw2();
hist_correlation_signal_subg0_subl_jet_reco_track_gen->Sumw2();
hist_correlation_signal_subg0_subl_jet_gen_track_reco->Sumw2();
hist_correlation_signal_subg0_subl_jet_gen_track_gen->Sumw2();
hist_reco_reco_2pcorrelation_signal->Sumw2();
hist_reco_reco_2pcorrelation_signal_subg0->Sumw2();
hist_reco_reco_2pcorrelation_signal_subcross->Sumw2();
hist_reco_reco_2pcorrelation_mixing->Sumw2();
hist_gen_gen_2pcorrelation_signal->Sumw2();
hist_gen_gen_2pcorrelation_signal_subg0->Sumw2();
hist_gen_gen_2pcorrelation_signal_subcross->Sumw2();
hist_gen_gen_2pcorrelation_mixing->Sumw2();

hist_reco_3rdjet->Sumw2();
hist_ref_3rdjet->Sumw2();
hist_gen_3rdjet->Sumw2();

}

// write QA histograms
/*
--> Arguments
isMC: true for MC and false for Data
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_QA_hist(bool isMC){
	Nevents->Write();
	Nev_recoreco->Write();
	Nev_recoreco_lead->Write();
	Nev_recoreco_subl->Write();
	if(isMC){ 
		Nev_recogen->Write();
		Nev_genreco->Write();
		Nev_gengen->Write();
		Nev_recogen_lead->Write();
		Nev_genreco_lead->Write();
		Nev_gengen_lead->Write();
		Nev_recogen_subl->Write();
		Nev_genreco_subl->Write();
		Nev_gengen_subl->Write();
		gen_mult->Write();
		gen_mult_weighted->Write();
		gen_mult_withonejet_weighted->Write();
		gen_mult_withdijets_weighted->Write();
		pthathist->Write(); 
		pthathist_weighted->Write();
	}
	reco_mult->Write();
	reco_mult_weighted->Write();
	reco_mult_withonejet_weighted->Write();
	reco_mult_withdijets_weighted->Write();
 	multiplicity->Write();
	multiplicity_weighted->Write();
	multiplicity_withonejet_weighted->Write();
	multiplicity_withdijets_weighted->Write();
	multiplicity_midmid_weighted->Write();
	multiplicity_midfwd_weighted->Write();
	multiplicity_midbkw_weighted->Write();
	multiplicity_fwdmid_weighted->Write();
	multiplicity_fwdfwd_weighted->Write();
	multiplicity_fwdbkw_weighted->Write();
	multiplicity_bkwmid_weighted->Write();
	multiplicity_bkwfwd_weighted->Write();
	multiplicity_bkwbkw_weighted->Write();
	multiplicity_weighted_at1->Write();
	multiplicity_withonejet_weighted_at1->Write();
	multiplicity_withdijets_weighted_at1->Write();
 	vzhist->Write();
 	vzhist_weighted->Write();
	vzhist_jet_weighted->Write();
	vzhist_dijet_weighted->Write();
	vxyhist->Write();
	vxyhist_weighted->Write();
	hfhist->Write();
	hfhist_weighted->Write();
	hfhist_onejet_weighted->Write();
	hfhist_dijet_weighted->Write();
	hfhistEta4->Write();
	hfhistEta4_weighted->Write();
	hfhistEta4_onejet_weighted->Write();
	hfhistEta4_dijet_weighted->Write();
	zdchist->Write();
	zdchist_weighted->Write();
	zdchist_onejet_weighted->Write();
	zdchist_dijet_weighted->Write();
	hfhistSum_weighted->Write();
	hfhistEta4Sum_weighted->Write();
	hfhistSum_onejet_weighted->Write();
	hfhistEta4Sum_onejet_weighted->Write();
	hfhistSum_dijet_weighted->Write();
	hfhistEta4Sum_dijet_weighted->Write();
	//tracks 
	//reco
	hist_reco_trk->Write();
	hist_reco_trk_corr->Write();
	hist_reco_trk_weighted->Write();
	//gen
	if(isMC){hist_gen_trk->Write(); hist_gen_trk_weighted->Write();}
	//jets 
	NJets->Write();
	trackmaxptinjethisto->Write();
	jettrackmaxptinjethisto->Write();
	histo_jetUE->Write();
	histo_jetAverageRho->Write();
	histo_jetcheckcorrection->Write();
	//reco
	hist_reco_jet->Write();
	hist_reco_jet_corr->Write();
	hist_reco_jet_weighted->Write();
	hist_reco_jet_corr_weighted->Write();
 	hist_reco_leadjet->Write();
	hist_reco_leadjet_weighted->Write();
	hist_reco_subljet->Write();
	hist_reco_subljet_weighted->Write();
	hist_reco_thrdjet_weighted->Write();
	if(isMC){
		hist_jes_reco_weighted->Write();
		hist_jes_reco_fromB_weighted->Write();
		hist_leadjes_reco_weighted->Write();
		hist_leadjes_reco_fromB_weighted->Write();
		hist_subleadjes_reco_weighted->Write();
		hist_subleadjes_reco_fromB_weighted->Write();
		hist_jetptclos_weighted->Write();
		hist_leadjetptclos_weighted->Write();
		hist_subljetptclos_weighted->Write();
		hist_averjetptclos_weighted->Write();
		hist_leadjetptclosremovesome_weighed->Write();
		hist_subljetptclosremovesome_weighed->Write();
		hist_averjetptclosremovesome_weighed->Write();
		hist_xjclos_weighted->Write();
		hist_xjclos_removesome_weighted->Write();
		hist_gen_jet->Write();
		hist_gen_jet_weighted->Write();
		hist_gen_leadjet->Write();
		hist_gen_leadjet_weighted->Write();
		hist_gen_subljet->Write();
		hist_gen_subljet_weighted->Write();
		hist_gen_thrdjet_weighted->Write();
		hist_ref_jet_weighted->Write();
		hist_ref_leadjet_weighted->Write();
		hist_ref_subljet_weighted->Write();
		hist_ref_thrdjet_weighted->Write();
	}
}

// Reco-Reco correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_recoreco_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl){

	if(doinclusive){
		hist_correlation_signal_jet_reco_track_reco->Write();
		hist_correlation_signal_subg0_jet_reco_track_reco->Write();
		hist_jet_from_reco_reco_sig->Write();
		hist_trk_from_reco_reco_sig->Write();
		hist_injet_reco_track_reco->Write();
		if(rotation) hist_correlation_rotation_jet_reco_track_reco->Write();
		if(mixing){
			hist_correlation_mixing_jet_reco_track_reco->Write();
			hist_jet_from_reco_reco_mix->Write();
			hist_trk_from_reco_reco_mix->Write();
		}
	}
	if(doleadsubl){
		hist_correlation_signal_lead_jet_reco_track_reco->Write();
		hist_correlation_signal_subg0_lead_jet_reco_track_reco->Write();
		hist_lead_jet_from_reco_reco_sig->Write();
		hist_LJ_trk_from_reco_reco_sig->Write();
		hist_inLeadjet_reco_track_reco->Write();
		hist_inSubljet_reco_track_reco->Write();
		if(rotation) hist_correlation_rotation_lead_jet_reco_track_reco->Write();
		if(mixing){
			hist_correlation_mixing_lead_jet_reco_track_reco->Write();
			hist_lead_jet_from_reco_reco_mix->Write();
			hist_LJ_trk_from_reco_reco_mix->Write();
		}
		hist_correlation_signal_subl_jet_reco_track_reco->Write();
		hist_correlation_signal_subg0_subl_jet_reco_track_reco->Write();
		hist_subl_jet_from_reco_reco_sig->Write();
		hist_SLJ_trk_from_reco_reco_sig->Write();
		if(rotation) hist_correlation_rotation_subl_jet_reco_track_reco->Write();
		if(mixing){
			hist_correlation_mixing_subl_jet_reco_track_reco->Write();
			hist_subl_jet_from_reco_reco_mix->Write();
			hist_SLJ_trk_from_reco_reco_mix->Write();
		}
	}	
}

// Reco-Gen correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_recogen_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl){

	if(doinclusive){
		hist_correlation_signal_jet_reco_track_gen->Write();
		hist_correlation_signal_subg0_jet_reco_track_gen->Write();
		hist_jet_from_reco_gen_sig->Write();
		hist_trk_from_reco_gen_sig->Write();
		hist_injet_reco_track_gen->Write();
		if(rotation) hist_correlation_rotation_jet_reco_track_gen->Write();
		if(mixing){
			hist_correlation_mixing_jet_reco_track_gen->Write();
			hist_jet_from_reco_gen_mix->Write();
			hist_trk_from_reco_gen_mix->Write();
		}
	}

	if(doleadsubl){
		hist_correlation_signal_lead_jet_reco_track_gen->Write();
		hist_correlation_signal_subg0_lead_jet_reco_track_gen->Write();
		hist_lead_jet_from_reco_gen_sig->Write();
		hist_LJ_trk_from_reco_gen_sig->Write();
		hist_inLeadjet_reco_track_gen->Write();
		hist_inSubljet_reco_track_gen->Write();
		if(rotation) hist_correlation_rotation_lead_jet_reco_track_gen->Write();
		if(mixing){
			hist_correlation_mixing_lead_jet_reco_track_gen->Write();
			hist_lead_jet_from_reco_gen_mix->Write();
			hist_LJ_trk_from_reco_gen_mix->Write();
		}
		hist_correlation_signal_subl_jet_reco_track_gen->Write();
		hist_correlation_signal_subg0_subl_jet_reco_track_gen->Write();
		hist_subl_jet_from_reco_gen_sig->Write();
		hist_SLJ_trk_from_reco_gen_sig->Write();
		if(rotation) hist_correlation_rotation_subl_jet_reco_track_gen->Write();
		if(mixing){
			hist_correlation_mixing_subl_jet_reco_track_gen->Write();
			hist_subl_jet_from_reco_gen_mix->Write();
			hist_SLJ_trk_from_reco_gen_mix->Write();
		}
	}	
}

// Gen-Reco correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_genreco_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl){

	if(doinclusive){
		hist_correlation_signal_jet_gen_track_reco->Write();
		hist_correlation_signal_subg0_jet_gen_track_reco->Write();
		hist_jet_from_gen_reco_sig->Write();
		hist_trk_from_gen_reco_sig->Write();
		hist_injet_gen_track_reco->Write();
		if(rotation) hist_correlation_rotation_jet_gen_track_reco->Write();
		if(mixing){
			hist_correlation_mixing_jet_gen_track_reco->Write();
			hist_jet_from_gen_reco_mix->Write();
			hist_trk_from_gen_reco_mix->Write();
		}
	}

	if(doleadsubl){
		hist_correlation_signal_lead_jet_gen_track_reco->Write();
		hist_correlation_signal_subg0_lead_jet_gen_track_reco->Write();
		hist_lead_jet_from_gen_reco_sig->Write();
		hist_LJ_trk_from_gen_reco_sig->Write();
		hist_inLeadjet_gen_track_reco->Write();
		hist_inSubljet_gen_track_reco->Write();
		if(rotation) hist_correlation_rotation_lead_jet_gen_track_reco->Write();
		if(mixing){
			hist_correlation_mixing_lead_jet_gen_track_reco->Write();
			hist_lead_jet_from_gen_reco_mix->Write();
			hist_LJ_trk_from_gen_reco_mix->Write();
		}
		hist_correlation_signal_subl_jet_gen_track_reco->Write();
		hist_correlation_signal_subg0_subl_jet_gen_track_reco->Write();
		hist_subl_jet_from_gen_reco_sig->Write();
		hist_SLJ_trk_from_gen_reco_sig->Write();
		if(rotation) hist_correlation_rotation_subl_jet_gen_track_reco->Write();
		if(mixing){
			hist_correlation_mixing_subl_jet_gen_track_reco->Write();
			hist_subl_jet_from_gen_reco_mix->Write();
			hist_SLJ_trk_from_gen_reco_mix->Write();
		}
	}	
}

// Gen-Gen correlation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_gengen_hist(bool mixing, bool rotation, bool doinclusive, bool doleadsubl){

	if(doinclusive){
		hist_correlation_signal_jet_gen_track_gen->Write();
		hist_correlation_signal_subg0_jet_gen_track_gen->Write();
		hist_jet_from_gen_gen_sig->Write();
		hist_trk_from_gen_gen_sig->Write();
		hist_injet_gen_track_gen->Write();
		if(rotation) hist_correlation_rotation_jet_gen_track_gen->Write();
		if(mixing){
			hist_correlation_mixing_jet_gen_track_gen->Write();
			hist_jet_from_gen_gen_mix->Write();
			hist_trk_from_gen_gen_mix->Write();
		}
	}

	if(doleadsubl){
		hist_correlation_signal_lead_jet_gen_track_gen->Write();
		hist_correlation_signal_subg0_lead_jet_gen_track_gen->Write();
		hist_lead_jet_from_gen_gen_sig->Write();
		hist_LJ_trk_from_gen_gen_sig->Write();
		hist_inLeadjet_gen_track_gen->Write();
		hist_inSubljet_gen_track_gen->Write();
		if(rotation) hist_correlation_rotation_lead_jet_gen_track_gen->Write();
		if(mixing){
			hist_correlation_mixing_lead_jet_gen_track_gen->Write();
			hist_lead_jet_from_gen_gen_mix->Write();
			hist_LJ_trk_from_gen_gen_mix->Write();
		}
		hist_correlation_signal_subl_jet_gen_track_gen->Write();
		hist_correlation_signal_subg0_subl_jet_gen_track_gen->Write();
		hist_subl_jet_from_gen_gen_sig->Write();
		hist_SLJ_trk_from_gen_gen_sig->Write();
		if(rotation) hist_correlation_rotation_subl_jet_gen_track_gen->Write();
		if(mixing){
			hist_correlation_mixing_subl_jet_gen_track_gen->Write();
			hist_subl_jet_from_gen_gen_mix->Write();
			hist_SLJ_trk_from_gen_gen_mix->Write();
		}
	}	
}

// Jet quenching histograms
/*
--> Arguments
isMC: true for MC and false for Data
*/
void w_dijet_hist(bool isMC){
	hist_reco_lead_reco_subl_quench_mid_mid->Write();
	hist_reco_lead_reco_subl_quench_mid_fwd->Write();
	hist_reco_lead_reco_subl_quench_mid_bkw->Write();
	hist_reco_lead_reco_subl_quench_fwd_mid->Write();
	hist_reco_lead_reco_subl_quench_bkw_mid->Write();
	hist_reco_lead_reco_subl_quench_fwd_fwd->Write();
	hist_reco_lead_reco_subl_quench_fwd_bkw->Write();
	hist_reco_lead_reco_subl_quench_bkw_fwd->Write();
	hist_reco_lead_reco_subl_quench_bkw_bkw->Write();
	hist_reco_leadEP_quench_plus_mid_mid->Write();
	hist_reco_leadEP_quench_plus_mid_fwd->Write();
	hist_reco_leadEP_quench_plus_mid_bkw->Write();
	hist_reco_leadEP_quench_plus_fwd_mid->Write();
	hist_reco_leadEP_quench_plus_bkw_mid->Write();
	hist_reco_leadEP_quench_plus_fwd_fwd->Write();
	hist_reco_leadEP_quench_plus_fwd_bkw->Write();
	hist_reco_leadEP_quench_plus_bkw_fwd->Write();
	hist_reco_leadEP_quench_plus_bkw_bkw->Write();
	hist_reco_sublEP_quench_plus_mid_mid->Write();
	hist_reco_sublEP_quench_plus_mid_fwd->Write();
	hist_reco_sublEP_quench_plus_mid_bkw->Write();
	hist_reco_sublEP_quench_plus_fwd_mid->Write();
	hist_reco_sublEP_quench_plus_bkw_mid->Write();
	hist_reco_sublEP_quench_plus_fwd_fwd->Write();
	hist_reco_sublEP_quench_plus_fwd_bkw->Write();
	hist_reco_sublEP_quench_plus_bkw_fwd->Write();
	hist_reco_sublEP_quench_plus_bkw_bkw->Write();	
	hist_reco_leadEP_quench_minus_mid_mid->Write();
	hist_reco_leadEP_quench_minus_mid_fwd->Write();
	hist_reco_leadEP_quench_minus_mid_bkw->Write();
	hist_reco_leadEP_quench_minus_fwd_mid->Write();
	hist_reco_leadEP_quench_minus_bkw_mid->Write();
	hist_reco_leadEP_quench_minus_fwd_fwd->Write();
	hist_reco_leadEP_quench_minus_fwd_bkw->Write();
	hist_reco_leadEP_quench_minus_bkw_fwd->Write();
	hist_reco_leadEP_quench_minus_bkw_bkw->Write();
	hist_reco_sublEP_quench_minus_mid_mid->Write();
	hist_reco_sublEP_quench_minus_mid_fwd->Write();
	hist_reco_sublEP_quench_minus_mid_bkw->Write();
	hist_reco_sublEP_quench_minus_fwd_mid->Write();
	hist_reco_sublEP_quench_minus_bkw_mid->Write();
	hist_reco_sublEP_quench_minus_fwd_fwd->Write();
	hist_reco_sublEP_quench_minus_fwd_bkw->Write();
	hist_reco_sublEP_quench_minus_bkw_fwd->Write();
	hist_reco_sublEP_quench_minus_bkw_bkw->Write();	
	hist_etaDijet_reco->Write();
	hist_etaDijet_CM_reco->Write();
	hist_yDijet_CM_reco->Write();
	hist_reco_3rdjet->Write();
	if(isMC){
		hist_gen_lead_gen_subl_quench_mid_mid->Write();
		hist_gen_lead_gen_subl_quench_mid_fwd->Write();
		hist_gen_lead_gen_subl_quench_mid_bkw->Write();
		hist_gen_lead_gen_subl_quench_fwd_mid->Write();
		hist_gen_lead_gen_subl_quench_bkw_mid->Write();
		hist_gen_lead_gen_subl_quench_fwd_fwd->Write();
		hist_gen_lead_gen_subl_quench_fwd_bkw->Write();
		hist_gen_lead_gen_subl_quench_bkw_fwd->Write();
		hist_gen_lead_gen_subl_quench_bkw_bkw->Write();
		hist_ref_lead_ref_subl_quench_mid_mid->Write();
		hist_ref_lead_ref_subl_quench_mid_fwd->Write();
		hist_ref_lead_ref_subl_quench_mid_bkw->Write();
		hist_ref_lead_ref_subl_quench_fwd_mid->Write();
		hist_ref_lead_ref_subl_quench_bkw_mid->Write();
		hist_ref_lead_ref_subl_quench_fwd_fwd->Write();
		hist_ref_lead_ref_subl_quench_fwd_bkw->Write();
		hist_ref_lead_ref_subl_quench_bkw_fwd->Write();
		hist_ref_lead_ref_subl_quench_bkw_bkw->Write();
		hist_ref_leadEP_quench_plus_mid_mid->Write();
		hist_ref_leadEP_quench_plus_mid_fwd->Write();
		hist_ref_leadEP_quench_plus_mid_bkw->Write();
		hist_ref_leadEP_quench_plus_fwd_mid->Write();
		hist_ref_leadEP_quench_plus_bkw_mid->Write();
		hist_ref_leadEP_quench_plus_fwd_fwd->Write();
		hist_ref_leadEP_quench_plus_fwd_bkw->Write();
		hist_ref_leadEP_quench_plus_bkw_fwd->Write();
		hist_ref_leadEP_quench_plus_bkw_bkw->Write();
		hist_ref_sublEP_quench_plus_mid_mid->Write();
		hist_ref_sublEP_quench_plus_mid_fwd->Write();
		hist_ref_sublEP_quench_plus_mid_bkw->Write();
		hist_ref_sublEP_quench_plus_fwd_mid->Write();
		hist_ref_sublEP_quench_plus_bkw_mid->Write();
		hist_ref_sublEP_quench_plus_fwd_fwd->Write();
		hist_ref_sublEP_quench_plus_fwd_bkw->Write();
		hist_ref_sublEP_quench_plus_bkw_fwd->Write();
		hist_ref_sublEP_quench_plus_bkw_bkw->Write();
		hist_gen_leadEP_quench_plus_mid_mid->Write();
		hist_gen_leadEP_quench_plus_mid_fwd->Write();
		hist_gen_leadEP_quench_plus_mid_bkw->Write();
		hist_gen_leadEP_quench_plus_fwd_mid->Write();
		hist_gen_leadEP_quench_plus_bkw_mid->Write();
		hist_gen_leadEP_quench_plus_fwd_fwd->Write();
		hist_gen_leadEP_quench_plus_fwd_bkw->Write();
		hist_gen_leadEP_quench_plus_bkw_fwd->Write();
		hist_gen_leadEP_quench_plus_bkw_bkw->Write();
		hist_gen_sublEP_quench_plus_mid_mid->Write();
		hist_gen_sublEP_quench_plus_mid_fwd->Write();
		hist_gen_sublEP_quench_plus_mid_bkw->Write();
		hist_gen_sublEP_quench_plus_fwd_mid->Write();
		hist_gen_sublEP_quench_plus_bkw_mid->Write();
		hist_gen_sublEP_quench_plus_fwd_fwd->Write();
		hist_gen_sublEP_quench_plus_fwd_bkw->Write();
		hist_gen_sublEP_quench_plus_bkw_fwd->Write();
		hist_gen_sublEP_quench_plus_bkw_bkw->Write();
		hist_ref_leadEP_quench_minus_mid_mid->Write();
		hist_ref_leadEP_quench_minus_mid_fwd->Write();
		hist_ref_leadEP_quench_minus_mid_bkw->Write();
		hist_ref_leadEP_quench_minus_fwd_mid->Write();
		hist_ref_leadEP_quench_minus_bkw_mid->Write();
		hist_ref_leadEP_quench_minus_fwd_fwd->Write();
		hist_ref_leadEP_quench_minus_fwd_bkw->Write();
		hist_ref_leadEP_quench_minus_bkw_fwd->Write();
		hist_ref_leadEP_quench_minus_bkw_bkw->Write();
		hist_ref_sublEP_quench_minus_mid_mid->Write();
		hist_ref_sublEP_quench_minus_mid_fwd->Write();
		hist_ref_sublEP_quench_minus_mid_bkw->Write();
		hist_ref_sublEP_quench_minus_fwd_mid->Write();
		hist_ref_sublEP_quench_minus_bkw_mid->Write();
		hist_ref_sublEP_quench_minus_fwd_fwd->Write();
		hist_ref_sublEP_quench_minus_fwd_bkw->Write();
		hist_ref_sublEP_quench_minus_bkw_fwd->Write();
		hist_ref_sublEP_quench_minus_bkw_bkw->Write();
		hist_gen_leadEP_quench_minus_mid_mid->Write();
		hist_gen_leadEP_quench_minus_mid_fwd->Write();
		hist_gen_leadEP_quench_minus_mid_bkw->Write();
		hist_gen_leadEP_quench_minus_fwd_mid->Write();
		hist_gen_leadEP_quench_minus_bkw_mid->Write();
		hist_gen_leadEP_quench_minus_fwd_fwd->Write();
		hist_gen_leadEP_quench_minus_fwd_bkw->Write();
		hist_gen_leadEP_quench_minus_bkw_fwd->Write();
		hist_gen_leadEP_quench_minus_bkw_bkw->Write();
		hist_gen_sublEP_quench_minus_mid_mid->Write();
		hist_gen_sublEP_quench_minus_mid_fwd->Write();
		hist_gen_sublEP_quench_minus_mid_bkw->Write();
		hist_gen_sublEP_quench_minus_fwd_mid->Write();
		hist_gen_sublEP_quench_minus_bkw_mid->Write();
		hist_gen_sublEP_quench_minus_fwd_fwd->Write();
		hist_gen_sublEP_quench_minus_fwd_bkw->Write();
		hist_gen_sublEP_quench_minus_bkw_fwd->Write();
		hist_gen_sublEP_quench_minus_bkw_bkw->Write();
		hist_etaDijet_ref->Write();
		hist_etaDijet_CM_ref->Write();
		hist_yDijet_CM_ref->Write();
		hist_etaDijet_gen->Write();
		hist_etaDijet_CM_gen->Write();
		hist_yDijet_CM_gen->Write();
		hist_ref_3rdjet->Write();
		hist_gen_3rdjet->Write();
	}
}

// 2PCorrelation histograms
/*
--> Arguments
mixing: true when doing mixing reference sample otherwise is false
rotation: true when doing rotation reference sample otherwise is false
doinclusive: in the case of measuring inclusive jet correlation
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_2pc_hist(bool isMC, bool mixing){
	hist_reco_reco_2pcorrelation_signal->Write();
	hist_reco_reco_2pcorrelation_signal_subg0->Write();
	hist_reco_reco_2pcorrelation_signal_subcross->Write();
	if(mixing) hist_reco_reco_2pcorrelation_mixing->Write();
	if(isMC){
		hist_gen_gen_2pcorrelation_signal->Write();
		hist_gen_gen_2pcorrelation_signal_subg0->Write();
		hist_gen_gen_2pcorrelation_signal_subcross->Write();
		if(mixing) hist_gen_gen_2pcorrelation_mixing->Write();
	}
}

// event plane histograms
void w_ep_hist(bool isMC){ 

	EP2_plus_flat->Write();
	EP2_minus_flat->Write();
	EP3_plus_flat->Write();
	EP3_minus_flat->Write();
	EP4_plus_flat->Write();
	EP4_minus_flat->Write();

	Dphi_flat_EP2_inclusive_minus->Write();
	Dphi_flat_EP2_leading_minus->Write();
	Dphi_flat_EP2_subleading_minus->Write();
	Dphi_flat_EP2_inclusive_plus->Write();
	Dphi_flat_EP2_leading_plus->Write();
	Dphi_flat_EP2_subleading_plus->Write();

	Dphi_flat_EP3_inclusive_minus->Write();
	Dphi_flat_EP3_leading_minus->Write();
	Dphi_flat_EP3_subleading_minus->Write();
	Dphi_flat_EP3_inclusive_plus->Write();
	Dphi_flat_EP3_leading_plus->Write();
	Dphi_flat_EP3_subleading_plus->Write();

	Dphi_flat_EP4_inclusive_minus->Write();
	Dphi_flat_EP4_leading_minus->Write();
	Dphi_flat_EP4_subleading_minus->Write();
	Dphi_flat_EP4_inclusive_plus->Write();
	Dphi_flat_EP4_leading_plus->Write();
	Dphi_flat_EP4_subleading_plus->Write();
	
	if(isMC){

		Dphi_GEN_flat_EP2_inclusive_minus->Write();
		Dphi_GEN_flat_EP2_leading_minus->Write();
		Dphi_GEN_flat_EP2_subleading_minus->Write();
		Dphi_GEN_flat_EP2_inclusive_plus->Write();
		Dphi_GEN_flat_EP2_leading_plus->Write();
		Dphi_GEN_flat_EP2_subleading_plus->Write();

		Dphi_GEN_flat_EP3_inclusive_minus->Write();
		Dphi_GEN_flat_EP3_leading_minus->Write();
		Dphi_GEN_flat_EP3_subleading_minus->Write();
		Dphi_GEN_flat_EP3_inclusive_plus->Write();
		Dphi_GEN_flat_EP3_leading_plus->Write();
		Dphi_GEN_flat_EP3_subleading_plus->Write();

		Dphi_GEN_flat_EP4_inclusive_minus->Write();
		Dphi_GEN_flat_EP4_leading_minus->Write();
		Dphi_GEN_flat_EP4_subleading_minus->Write();
		Dphi_GEN_flat_EP4_inclusive_plus->Write();
		Dphi_GEN_flat_EP4_leading_plus->Write();
		Dphi_GEN_flat_EP4_subleading_plus->Write();
		
	}
	
	Dphi_EP2_flat_trk_minus->Write();
	Dphi_EP2_flat_trk_plus->Write();
	Dphi_EP3_flat_trk_minus->Write();
	Dphi_EP3_flat_trk_plus->Write();
	Dphi_EP4_flat_trk_minus->Write();
	Dphi_EP4_flat_trk_plus->Write();
	
	if(isMC){
		Dphi_GEN_EP2_flat_trk_minus->Write();
		Dphi_GEN_EP2_flat_trk_plus->Write();
		Dphi_GEN_EP3_flat_trk_minus->Write();
		Dphi_GEN_EP3_flat_trk_plus->Write();
		Dphi_GEN_EP4_flat_trk_minus->Write();
		Dphi_GEN_EP4_flat_trk_plus->Write();
	}
	
}


void w_unf_hist(){ 
	hist_jetunf_weighted->Write();
	hist_leadjetunf_weighted->Write();
	hist_subljetunf_weighted->Write();
	hist_xjunf_weighted->Write();
	hist_leadjetunf_match_weighted->Write();
	hist_subljetunf_match_weighted->Write();
	hist_leadjetunf_swap_weighted->Write();
	hist_subljetunf_swap_weighted->Write();
	hist_xjunf_match_weighted->Write();
	hist_xjunf_swap_weighted->Write();
	hist_leadjetunf_gensmear->Write();
	hist_leadjetunf_recosmear->Write();
	hist_subljetunf_gensmear->Write();
	hist_subljetunf_recosmear->Write();
	hist_xjunf_gensmear->Write();
	hist_xjunf_recosmear->Write();
	hist_xjunf_gensmear_fromLSL->Write();
	hist_xjunf_recosmear_fromLSL->Write();
	hist_xjunf_gensmear_fromLSL_4D->Write();
	hist_xjunf_recosmear_fromLSL_4D->Write();
	hist_leadjetunf_gensmear_fromInclJet->Write();
	hist_leadjetunf_recosmear_fromInclJet->Write();
	hist_subljetunf_gensmear_fromInclJet->Write();
	hist_subljetunf_recosmear_fromInclJet->Write();
	hist_xjunf_gensmear_fromLSL_InclJet->Write();
	hist_xjunf_recosmear_fromLSL_InclJet->Write();
	hist_leadjetunf_gensmear_4D->Write();
	hist_subljetunf_gensmear_4D->Write();
	hist_leadjetunf_recosmear_4D->Write();
	hist_subljetunf_recosmear_4D->Write();
	hist_jetunf_weighted_4D->Write();
	hist_jetunf_match_weighted_4D->Write();
	hist_jetunf_match_weighted_sym_4D->Write();
	hist_jetunf_match_weighted_symhalf_4D->Write();
	hist_jetunf_swap_weighted_4D->Write();
}
