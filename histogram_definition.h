#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // call inputs

int trkbinsize = (int) trk_pt_bins.size(); // track bins for jet-track correlation
int multbinsize = (int) multiplicity_centrality_bins.size();// multiplicity or centrality bins for jet-track correlation
int ptavebinsize = (int) pt_ave_bins.size();// average pT between leading and subleading jets
int extrabinsize = (int) extra_bins.size();// any additional dependency you wanna add (be carefull about memory)

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

// Vertex
// Axis : 0 -> Vz, 1 -> event multiplicity, 2 -> extra dependency
int	bins_vz[3]   	=   {  60   ,   multbinsize-1    		 ,  extrabinsize-1};
double xmin_vz[3]   =   { -15.5 ,   0.0  					 ,  0.0};
double xmax_vz[3]   =   {  15.5 ,   (double) multbinsize-1   ,  (double)extrabinsize-1};
THnSparseD *vzhist = new THnSparseD("vzhist", "vzhist", 3, bins_vz, xmin_vz, xmax_vz);
THnSparseD *vzhist_weighted = new THnSparseD("vzhist_weighted", "vzhist_weighted", 3, bins_vz, xmin_vz, xmax_vz);
THnSparseD *vzhist_jet_weighted = new THnSparseD("vzhist_jet_weighted", "vzhist_jet_weighted", 3, bins_vz, xmin_vz, xmax_vz);
THnSparseD *vzhist_dijet_weighted = new THnSparseD("vzhist_dijet_weighted", "vzhist_dijet_weighted", 3, bins_vz, xmin_vz, xmax_vz);
// Axis : 0 -> Vx, 1 -> Vy, 1 -> event multiplicity, 2 -> extra dependency
int	bins_vxy[4]   	 =   {  100 , 100  , multbinsize-1    		 ,  extrabinsize-1};
double xmin_vxy[4]   =   { -1.0 , -1.0 , 0.0  					 ,  0.0};
double xmax_vxy[4]   =   {  1.0 ,  1.0 , (double) multbinsize-1  ,  (double)extrabinsize-1};
THnSparseD *vxyhist = new THnSparseD("vxyhist", "vxyhist", 3, bins_vxy, xmin_vxy, xmax_vxy);
THnSparseD *vxyhist_weighted = new THnSparseD("vxyhist_weighted", "vxyhist_weighted", 4, bins_vxy, xmin_vxy, xmax_vxy);

// Pthat
// Axis : 0 -> Pthat, 1 -> event multiplicity, 2 -> extra dependency
int	bins_pthat[3]      =   {  100 	 ,   multbinsize-1    		  ,  extrabinsize-1};
double xmin_pthat[3]   =   {  0.0 	 ,   0.0  					  ,  0.0};
double xmax_pthat[3]   =   {  1000.0 ,   (double) multbinsize-1   ,  (double)extrabinsize-1};
THnSparseD *pthathist = new THnSparseD("pthathist", "pthathist", 3, bins_pthat, xmin_pthat, xmax_pthat);
THnSparseD *pthathist_weighted = new THnSparseD("pthathist_weighted", "pthathist_weighted", 3, bins_pthat, xmin_pthat, xmax_pthat);

// event plane histograms
// Axis : 0 -> EP multiplicity, 1 -> qvector, 2 -> PsiEP, 3 -> event multiplicity, 4 -> extra dependency
int	bins_EP[5]   	=   { 40	,  200 ,   32		     , multbinsize-1			,  extrabinsize-1};
double xmin_EP[5]   =   { 0	 	,  0   ,   -TMath::Pi()  , 0.0  	    			,  0.0};
double xmax_EP[5]   =   { 400   ,  100 ,   TMath::Pi()   , (double) multbinsize-1	,  (double)extrabinsize-1};
THnSparseD *EP2_plus_flat = new THnSparseD("EP2_plus_flat", "EP2_plus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP2_minus_flat = new THnSparseD("EP2_minus_flat", "EP2_minus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP3_plus_flat = new THnSparseD("EP3_plus_flat", "EP3_plus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP3_minus_flat = new THnSparseD("EP3_minus_flat", "EP3_minus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP4_plus_flat = new THnSparseD("EP4_plus_flat", "EP4_plus_flat", 5, bins_EP, xmin_EP, xmax_EP);
THnSparseD *EP4_minus_flat = new THnSparseD("EP4_minus_flat", "EP4_minus_flat", 5, bins_EP, xmin_EP, xmax_EP);

// Track/Particle histograms
// Axis : 0 -> track pT, 1 -> trk eta, 2 -> trk phi, 3 -> multiplicity bin, 4 -> extra dependency
int	bins_trk[5]      =   { 100   ,  30  ,   32		     , multbinsize-1			,  extrabinsize-1};
double xmin_trk[5]   =   { 0.0   , -3.0 ,   -TMath::Pi() , 0.0  	    			,  0.0};
double xmax_trk[5]   =   { 50.0  ,  3.0 ,   TMath::Pi()  , (double) multbinsize-1	,  (double)extrabinsize-1};
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
// Axis : 0 -> delta phi between track and EP, 1 -> trkbin, 2 -> multbin
int	bins_TRKEP[4]      =   { 40				   		,  trkbinsize-1 		,  multbinsize-1			,  extrabinsize-1};
double xmin_TRKEP[4]   =   { -TMath::Pi()/2.0	    ,  0.0 					,  0.0  	    			,  0.0};
double xmax_TRKEP[4]   =   { 3.0*TMath::Pi()/2.0 	,  (double)trkbinsize-1 ,  (double) multbinsize-1   ,  (double)extrabinsize-1};
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
int	bins_NJETS[3]      =   {  30   ,   multbinsize-1    		 ,  extrabinsize-1};
double xmin_NJETS[3]   =   {  0.0  ,   0.0  					 ,  0.0};
double xmax_NJETS[3]   =   {  30.0 ,   (double) multbinsize-1   ,  (double)extrabinsize-1};
THnSparseD *NJets = new THnSparseD("NJets", "NJets", 3, bins_NJETS, xmin_NJETS, xmax_NJETS);

// trackmax histogram
// Axis : 0 -> delta phi between jet and EP, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_trkmax[3]      =   {  1000   ,   multbinsize-1    		 ,  extrabinsize-1};
double xmin_trkmax[3]   =   {  0.0 	  ,   0.0  					 ,  0.0};
double xmax_trkmax[3]   =   {  1000.0 ,   (double) multbinsize-1   ,  (double)extrabinsize-1};
THnSparseD *trackmaxptinjethisto = new THnSparseD("trackmaxptinjethisto", "trackmaxptinjethisto", 3, bins_trkmax, xmin_trkmax, xmax_trkmax);

//correlations to EP
// Axis : 0 -> delta phi between jet and EP, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_JETEP[3]      =   { 200			   		,  multbinsize-1			,  extrabinsize-1};
double xmin_JETEP[3]   =   { -TMath::Pi()/2.0	    ,  0.0  	    			,  0.0};
double xmax_JETEP[3]   =   { 3.0*TMath::Pi()/2.0 	,  (double) multbinsize-1   ,  (double)extrabinsize-1};

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
int	bins_jet[5]      =   { 100	  ,  40  ,   32		      , multbinsize-1          , extrabinsize-1};
double xmin_jet[5]   =   { 0.0	  , -4.0 ,   -TMath::Pi() , 0.0                    , 0.0};
double xmax_jet[5]   =   { 500.0  ,  4.0 ,   TMath::Pi()  , (double) multbinsize-1 , (double)extrabinsize-1};
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
// --> Gen
THnSparseD *hist_gen_leadjet = new THnSparseD("hist_gen_leadjet", "hist_gen_leadjet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_leadjet_weighted = new THnSparseD("hist_gen_leadjet_weighted", "hist_gen_leadjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_subljet = new THnSparseD("hist_gen_subljet", "hist_gen_subljet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_subljet_weighted = new THnSparseD("hist_gen_subljet_weighted", "hist_gen_subljet_weighted", 5, bins_jet, xmin_jet, xmax_jet);

// Jet Energy Scale (JES) and Jet Energy Resolution (JER)
int	bins_jes[6]   =      { 200  ,  100  ,  80  ,  8, multbinsize-1          , extrabinsize-1};
double xmin_jes[6]   =   { 0.0  ,  0    , -4.0 ,  0, 0                      , 0.0};
double xmax_jes[6]   =   { 10.0 ,  1000 ,  4.0 ,  8, (double) multbinsize-1 , (double)extrabinsize-1};
THnSparseD *hist_jes_reco_weighted = new THnSparseD("hist_jes_reco_weighted", "hist_jes_reco_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_jes_reco_fromB_weighted = new THnSparseD("hist_jes_reco_fromB_weighted", "hist_jes_reco_fromB_weighted", 6, bins_jes, xmin_jes, xmax_jes);

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

// --------------------------------------------------------------------------------------------------------
// Quenching studies
// Axis : 0 -> Xj, 1 -> Aj, 2 -> delta phi, 3 -> multiplicity , 4 -> jet pT average, 5 -> extra dependency
int	bins_quenc[6]   =      { 20   , 20	  , 30		    , multbinsize-1		  	 ,	 ptavebinsize-1			  ,  extrabinsize-1};
double xmin_quenc[6]   =   { 0.0  , 0.0   , 0.0		    , 0.0		   			 ,	 0.0					  ,	 0};
double xmax_quenc[6]   =   { 1.0  , 1.0   , TMath::Pi() , (double) multbinsize-1 ,   (double) ptavebinsize-1  ,  (double) extrabinsize-1};
THnSparseD *hist_reco_lead_reco_subl_quench_mid_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_mid", "hist_reco_lead_reco_subl_quench_mid_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_mid_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_fwd", "hist_reco_lead_reco_subl_quench_mid_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_mid_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_bkw", "hist_reco_lead_reco_subl_quench_mid_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_mid", "hist_reco_lead_reco_subl_quench_fwd_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_mid", "hist_reco_lead_reco_subl_quench_bkw_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_fwd", "hist_reco_lead_reco_subl_quench_fwd_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_bkw", "hist_reco_lead_reco_subl_quench_fwd_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_fwd", "hist_reco_lead_reco_subl_quench_bkw_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_bkw", "hist_reco_lead_reco_subl_quench_bkw_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_mid", "hist_gen_lead_gen_subl_quench_mid_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_fwd", "hist_gen_lead_gen_subl_quench_mid_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_bkw", "hist_gen_lead_gen_subl_quench_mid_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_mid", "hist_gen_lead_gen_subl_quench_fwd_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_mid", "hist_gen_lead_gen_subl_quench_bkw_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_fwd", "hist_gen_lead_gen_subl_quench_fwd_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_bkw", "hist_gen_lead_gen_subl_quench_fwd_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_fwd", "hist_gen_lead_gen_subl_quench_bkw_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_bkw", "hist_gen_lead_gen_subl_quench_bkw_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_mid", "hist_ref_lead_ref_subl_quench_mid_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_fwd", "hist_ref_lead_ref_subl_quench_mid_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_bkw", "hist_ref_lead_ref_subl_quench_mid_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_mid", "hist_ref_lead_ref_subl_quench_fwd_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_mid", "hist_ref_lead_ref_subl_quench_bkw_mid", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_fwd", "hist_ref_lead_ref_subl_quench_fwd_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_bkw", "hist_ref_lead_ref_subl_quench_fwd_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_fwd", "hist_ref_lead_ref_subl_quench_bkw_fwd", 6, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_bkw", "hist_ref_lead_ref_subl_quench_bkw_bkw", 6, bins_quenc, xmin_quenc, xmax_quenc);

// Assymetry studies
// Axis : 0 -> etaDijet, 1 -> delta eta / 2, 2 -> Xj, 3 -> Aj, 4 -> delta phi, 5 -> x_p, 6 -> x_Pb, 7 -> multiplicity, 8 -> jet pT average, 9 -> extra dependency
int	bins_etaDijet[10]      =   {  40   ,  16  , 20	  , 20 , 6		     , 30    ,  30     ,	 multbinsize-1		  	 ,	 ptavebinsize-1			  ,  extrabinsize-1};
double xmin_etaDijet[10]   =   { -4.0  , -4.0,  0.0   , 0.0, 0.0 	     , 0.001 ,  0.001  ,     0.0		   			 ,	 0.0					  ,	 0.0};
double xmax_etaDijet[10]   =   {  4.0  ,  4.0,  1.0   , 1.0, TMath::Pi() , 1.0   ,  1.0    ,     (double) multbinsize-1  ,   (double) ptavebinsize-1  ,  (double) extrabinsize-1};
THnSparseD *hist_etaDijet_reco = new THnSparseD("hist_etaDijet_reco", "hist_etaDijet_reco", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_reco = new THnSparseD("hist_etaDijet_CM_reco", "hist_etaDijet_CM_reco", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_ref = new THnSparseD("hist_etaDijet_ref", "hist_etaDijet_ref", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_ref = new THnSparseD("hist_etaDijet_CM_ref", "hist_etaDijet_CM_ref", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_gen = new THnSparseD("hist_etaDijet_gen", "hist_etaDijet_gen", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_gen = new THnSparseD("hist_etaDijet_CM_gen", "hist_etaDijet_CM_gen", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);

// Axis : 0 -> in-jet multiplicity, 1 -> multiplicity, 2 -> extra dimension
int	bins_injettrk[3]   	  =   { 100 ,  multbinsize-1			  ,  extrabinsize-1};
double xmin_injettrk[3]   =   { 0.0 , 0.0					  	  ,	 0.0};
double xmax_injettrk[3]   =   { 100 , (double) multbinsize-1	  ,  (double) extrabinsize-1};
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
int	bins_jettrk[5]      =   { 200				    , 400  ,   trkbinsize-1		 	 , multbinsize-1			  ,  extrabinsize-1};
double xmin_jettrk[5]   =   { -TMath::Pi()/2.0		, -4.0 ,   0					 , 0.0					  	  ,	 0.0};
double xmax_jettrk[5]   =   { 3.0*TMath::Pi()/2.0 	, 4.0  ,   (double) trkbinsize-1 , (double) multbinsize-1  	  ,  (double) extrabinsize-1};

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
int	bins_2pc[5]      =   { 40				     ,  80   ,   trkbinsize-1		    , multbinsize-1			  ,  extrabinsize-1};
double xmin_2pc[5]   =   { -TMath::Pi()/2.0	     , -4.0  ,   0					    , 0.0					  ,	 0.0};
double xmax_2pc[5]   =   { 3.0*TMath::Pi()/2.0   ,  4.0  ,   (double) trkbinsize-1  , (double) multbinsize-1  ,  (double) extrabinsize-1};

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
hist_jes_reco_weighted->Sumw2();
hist_jes_reco_fromB_weighted->Sumw2();
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
hist_etaDijet_reco->Sumw2();
hist_etaDijet_CM_reco->Sumw2();
hist_etaDijet_ref->Sumw2();
hist_etaDijet_CM_ref->Sumw2();
hist_etaDijet_gen->Sumw2();
hist_etaDijet_CM_gen->Sumw2();
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

}

// write QA histograms
/*
--> Arguments
isMC: true for MC and false for Data
doleadsubl: in the case of measuring leading and subleading jet correlation
*/
void w_QA_hist(bool isMC, bool doleadsubl){
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
	//reco
	hist_reco_jet->Write();
	hist_reco_jet_corr->Write();
	hist_reco_jet_weighted->Write();
	hist_reco_jet_corr_weighted->Write();
 	hist_reco_leadjet->Write();
	hist_reco_leadjet_weighted->Write();
	hist_reco_subljet->Write();
	hist_reco_subljet_weighted->Write();
	if(isMC){
		hist_jes_reco_weighted->Write();
		hist_jes_reco_fromB_weighted->Write();
		hist_gen_jet->Write();
		hist_gen_jet_weighted->Write();
		hist_gen_leadjet->Write();
		hist_gen_leadjet_weighted->Write();
		hist_gen_subljet->Write();
		hist_gen_subljet_weighted->Write();
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
	hist_etaDijet_reco->Write();
	hist_etaDijet_CM_reco->Write();
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
		hist_etaDijet_ref->Write();
		hist_etaDijet_CM_ref->Write();
		hist_etaDijet_gen->Write();
		hist_etaDijet_CM_gen->Write();
	}
}

// JES and JER histograms
void w_jes_jer_hist(){ 
	hist_jes_reco_weighted->Write();
	hist_jes_reco_fromB_weighted->Write();
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