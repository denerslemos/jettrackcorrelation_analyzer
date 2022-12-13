#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // call inputs

int trkbinsize = (int) trk_pt_bins.size(); // track bins for jet-track correlation
int multbinsize = (int) multiplicity_centrality_bins.size();// multiplicity or centrality bins for jet-track correlation

// -------------------------------- QA plots --------------------------------
// Event quantities
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
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 500, 0.0, 500.0);
TH1D *multiplicity_weighted = new TH1D("multiplicity_weighted", "multiplicity_weighted", 500, 0.0, 500.0);
TH1D *vzhist = new TH1D("vzhist", "vzhist", 80, -20, 20);
TH1D *vzhist_weighted = new TH1D("vzhist_weighted", "vzhist_weighted", 80, -20, 20);
TH1D *pthathist = new TH1D("pthathist", "pthathist", 230, 0, 460);
TH1D *pthathist_weighted = new TH1D("pthathist_weighted", "pthathist_weighted", 230, 0, 460);

// Track/Particle histograms
int    bins4D_trk[4]   =   { 500   ,  50  ,   64           , multbinsize-1};
double xmin4D_trk[4]   =   { 0.0   , -2.5 ,   -TMath::Pi() , 0};
double xmax4D_trk[4]   =   { 20.0  ,  2.5 ,   TMath::Pi()  , (double) multbinsize-1};

// --> Reco
THnSparseD *hist_reco_trk = new THnSparseD("hist_reco_trk", "hist_reco_trk", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_reco_trk_corr = new THnSparseD("hist_reco_trk_corr", "hist_reco_trk_corr", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_reco_trk_weighted = new THnSparseD("hist_reco_trk_weighted", "hist_reco_trk_weighted", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// --> Gen
THnSparseD *hist_gen_trk = new THnSparseD("hist_gen_trk", "hist_gen_trk", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_gen_trk_weighted = new THnSparseD("hist_gen_trk_weighted", "hist_gen_trk_weighted", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// Tracks from Correlations
// Inclusive
THnSparseD *hist_trk_from_reco_reco_sig = new THnSparseD("hist_trk_from_reco_reco_sig", "hist_trk_from_reco_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_reco_gen_sig = new THnSparseD("hist_trk_from_reco_gen_sig", "hist_trk_from_reco_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_reco_sig = new THnSparseD("hist_trk_from_gen_reco_sig", "hist_trk_from_gen_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_gen_sig = new THnSparseD("hist_trk_from_gen_gen_sig", "hist_trk_from_gen_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_reco_reco_mix = new THnSparseD("hist_trk_from_reco_reco_mix", "hist_trk_from_reco_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_reco_gen_mix = new THnSparseD("hist_trk_from_reco_gen_mix", "hist_trk_from_reco_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_reco_mix = new THnSparseD("hist_trk_from_gen_reco_mix", "hist_trk_from_gen_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_trk_from_gen_gen_mix = new THnSparseD("hist_trk_from_gen_gen_mix", "hist_trk_from_gen_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
// Leading
THnSparseD *hist_LJ_trk_from_reco_reco_sig = new THnSparseD("hist_LJ_trk_from_reco_reco_sig", "hist_LJ_trk_from_reco_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_reco_gen_sig = new THnSparseD("hist_LJ_trk_from_reco_gen_sig", "hist_LJ_trk_from_reco_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_reco_sig = new THnSparseD("hist_LJ_trk_from_gen_reco_sig", "hist_LJ_trk_from_gen_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_gen_sig = new THnSparseD("hist_LJ_trk_from_gen_gen_sig", "hist_LJ_trk_from_gen_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_reco_reco_mix = new THnSparseD("hist_LJ_trk_from_reco_reco_mix", "hist_LJ_trk_from_reco_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_reco_gen_mix = new THnSparseD("hist_LJ_trk_from_reco_gen_mix", "hist_LJ_trk_from_reco_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_reco_mix = new THnSparseD("hist_LJ_trk_from_gen_reco_mix", "hist_LJ_trk_from_gen_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_LJ_trk_from_gen_gen_mix = new THnSparseD("hist_LJ_trk_from_gen_gen_mix", "hist_LJ_trk_from_gen_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
// Subleading
THnSparseD *hist_SLJ_trk_from_reco_reco_sig = new THnSparseD("hist_SLJ_trk_from_reco_reco_sig", "hist_SLJ_trk_from_reco_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_reco_gen_sig = new THnSparseD("hist_SLJ_trk_from_reco_gen_sig", "hist_SLJ_trk_from_reco_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_reco_sig = new THnSparseD("hist_SLJ_trk_from_gen_reco_sig", "hist_SLJ_trk_from_gen_reco_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_gen_sig = new THnSparseD("hist_SLJ_trk_from_gen_gen_sig", "hist_SLJ_trk_from_gen_gen_sig", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_reco_reco_mix = new THnSparseD("hist_SLJ_trk_from_reco_reco_mix", "hist_SLJ_trk_from_reco_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_reco_gen_mix = new THnSparseD("hist_SLJ_trk_from_reco_gen_mix", "hist_SLJ_trk_from_reco_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_reco_mix = new THnSparseD("hist_SLJ_trk_from_gen_reco_mix", "hist_SLJ_trk_from_gen_reco_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);
THnSparseD *hist_SLJ_trk_from_gen_gen_mix = new THnSparseD("hist_SLJ_trk_from_gen_gen_mix", "hist_SLJ_trk_from_gen_gen_mix", 4, bins4D_trk, xmin4D_trk, xmax4D_trk);

// Jet histograms
int    bins4D_jet[4]   =   { 200   ,  50  ,   64           , multbinsize-1};
double xmin4D_jet[4]   =   { 0.0   , -2.5 ,   -TMath::Pi() , 0};
double xmax4D_jet[4]   =   { 1000.0  ,  2.5 ,   TMath::Pi() , (double) multbinsize-1};

// --> Reco
TH1D *hist_reco_jet_weighted_nocut = new TH1D("hist_reco_jet_weighted_nocut", "hist_reco_jet_weighted_nocut", 100, 0.0, 500.0);
THnSparseD *hist_reco_jet = new THnSparseD("hist_reco_jet", "hist_reco_jet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_jet_corr = new THnSparseD("hist_reco_jet_corr", "hist_reco_jet_corr", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_jet_weighted = new THnSparseD("hist_reco_jet_weighted", "hist_reco_jet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_jet_corr_weighted = new THnSparseD("hist_reco_jet_corr_weighted", "hist_reco_jet_corr_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
// --> Gen
TH1D *hist_gen_jet_weighted_nocut = new TH1D("hist_gen_jet_weighted_nocut", "hist_gen_jet_weighted_nocut", 100, 0.0, 500.0);
THnSparseD *hist_gen_jet = new THnSparseD("hist_gen_jet", "hist_gen_jet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_jet_weighted = new THnSparseD("hist_gen_jet_weighted", "hist_gen_jet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// Leading and Subleading Jets
// --> Reco
TH1D *hist_reco_leadjet_pt_nocut = new TH1D("hist_reco_leadjet_pt_nocut", "hist_reco_leadjet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_reco_leadjet_pt_nocut_weighted  = new TH1D("hist_reco_leadjet_pt_nocut_weighted", "hist_reco_leadjet_pt_nocut_weighted", 100, 0.0, 500.0);
TH1D *hist_reco_subljet_pt_nocut = new TH1D("hist_reco_subljet_pt_nocut", "hist_reco_subljet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_reco_subljet_pt_nocut_weighted = new TH1D("hist_reco_subljet_pt_nocut_weighted", "hist_reco_subljet_pt_nocut_weighted", 100, 0.0, 500.0);

THnSparseD *hist_reco_leadjet = new THnSparseD("hist_reco_leadjet", "hist_reco_leadjet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_leadjet_weighted = new THnSparseD("hist_reco_leadjet_weighted", "hist_reco_leadjet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_subljet = new THnSparseD("hist_reco_subljet", "hist_reco_subljet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_reco_subljet_weighted = new THnSparseD("hist_reco_subljet_weighted", "hist_reco_subljet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// --> Gen
TH1D *hist_gen_leadjet_pt_nocut = new TH1D("hist_gen_leadjet_pt_nocut", "hist_gen_leadjet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_gen_leadjet_pt_nocut_weighted = new TH1D("hist_gen_leadjet_pt_nocut_weighted", "hist_gen_leadjet_pt_nocut_weighted", 100, 0.0, 500.0);
TH1D *hist_gen_subljet_pt_nocut = new TH1D("hist_gen_subljet_pt_nocut", "hist_gen_subljet_pt_nocut", 100, 0.0, 500.0);
TH1D *hist_gen_subljet_pt_nocut_weighted = new TH1D("hist_gen_subljet_pt_nocut_weighted", "hist_gen_subljet_pt_nocut_weighted", 100, 0.0, 500.0);

THnSparseD *hist_gen_leadjet = new THnSparseD("hist_gen_leadjet", "hist_gen_leadjet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_leadjet_weighted = new THnSparseD("hist_gen_leadjet_weighted", "hist_gen_leadjet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_subljet = new THnSparseD("hist_gen_subljet", "hist_gen_subljet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_gen_subljet_weighted = new THnSparseD("hist_gen_subljet_weighted", "hist_gen_subljet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// Jets from Correlations
// Inclusive
THnSparseD *hist_jet_from_reco_reco_sig = new THnSparseD("hist_jet_from_reco_reco_sig", "hist_jet_from_reco_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_reco_gen_sig = new THnSparseD("hist_jet_from_reco_gen_sig", "hist_jet_from_reco_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_reco_sig = new THnSparseD("hist_jet_from_gen_reco_sig", "hist_jet_from_gen_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_gen_sig = new THnSparseD("hist_jet_from_gen_gen_sig", "hist_jet_from_gen_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_reco_reco_mix = new THnSparseD("hist_jet_from_reco_reco_mix", "hist_jet_from_reco_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_reco_gen_mix = new THnSparseD("hist_jet_from_reco_gen_mix", "hist_jet_from_reco_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_reco_mix = new THnSparseD("hist_jet_from_gen_reco_mix", "hist_jet_from_gen_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_jet_from_gen_gen_mix = new THnSparseD("hist_jet_from_gen_gen_mix", "hist_jet_from_gen_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
// Leading
THnSparseD *hist_lead_jet_from_reco_reco_sig = new THnSparseD("hist_lead_jet_from_reco_reco_sig", "hist_lead_jet_from_reco_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_reco_gen_sig = new THnSparseD("hist_lead_jet_from_reco_gen_sig", "hist_lead_jet_from_reco_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_reco_sig = new THnSparseD("hist_lead_jet_from_gen_reco_sig", "hist_lead_jet_from_gen_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_gen_sig = new THnSparseD("hist_lead_jet_from_gen_gen_sig", "hist_lead_jet_from_gen_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_reco_reco_mix = new THnSparseD("hist_lead_jet_from_reco_reco_mix", "hist_lead_jet_from_reco_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_reco_gen_mix = new THnSparseD("hist_lead_jet_from_reco_gen_mix", "hist_lead_jet_from_reco_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_reco_mix = new THnSparseD("hist_lead_jet_from_gen_reco_mix", "hist_lead_jet_from_gen_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_lead_jet_from_gen_gen_mix = new THnSparseD("hist_lead_jet_from_gen_gen_mix", "hist_lead_jet_from_gen_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
// Subleading
THnSparseD *hist_subl_jet_from_reco_reco_sig = new THnSparseD("hist_subl_jet_from_reco_reco_sig", "hist_subl_jet_from_reco_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_reco_gen_sig = new THnSparseD("hist_subl_jet_from_reco_gen_sig", "hist_subl_jet_from_reco_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_reco_sig = new THnSparseD("hist_subl_jet_from_gen_reco_sig", "hist_subl_jet_from_gen_reco_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_gen_sig = new THnSparseD("hist_subl_jet_from_gen_gen_sig", "hist_subl_jet_from_gen_gen_sig", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_reco_reco_mix = new THnSparseD("hist_subl_jet_from_reco_reco_mix", "hist_subl_jet_from_reco_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_reco_gen_mix = new THnSparseD("hist_subl_jet_from_reco_gen_mix", "hist_subl_jet_from_reco_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_reco_mix = new THnSparseD("hist_subl_jet_from_gen_reco_mix", "hist_subl_jet_from_gen_reco_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_subl_jet_from_gen_gen_mix = new THnSparseD("hist_subl_jet_from_gen_gen_mix", "hist_subl_jet_from_gen_gen_mix", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

// --------------------------------------------------------------------------------------------------------
// Quenching studies
// Axis : 0 -> Aj, 1 -> Xj, 2 -> delta phi, 3 -> multiplicity
int    bins4D_quenc[4]   =   { 20   , 20    , 40          , multbinsize-1          };
double xmin4D_quenc[4]   =   { 0.0  , 0.0   , 0.0         , 0                      };
double xmax4D_quenc[4]   =   { 1.0  , 1.0   , TMath::Pi() , (double) multbinsize-1 };
THnSparseD *hist_reco_lead_reco_subl_quench = new THnSparseD("hist_reco_lead_reco_subl_quench", "hist_reco_lead_reco_subl_quench", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench = new THnSparseD("hist_gen_lead_gen_subl_quench", "hist_gen_lead_gen_subl_quench", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench = new THnSparseD("hist_ref_lead_ref_subl_quench", "hist_ref_lead_ref_subl_quench", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);

// Axis : 0 -> Aj, 1 -> Xj, 2 -> delta phi 2PC, 3 -> multiplicity
int    bins4D_quenc_2pc[4]   =   { 20   , 20    , 40                  , multbinsize-1          };
double xmin4D_quenc_2pc[4]   =   { 0.0  , 0.0   , -TMath::Pi()/2.0    , 0                      };
double xmax4D_quenc_2pc[4]   =   { 1.0  , 1.0   , 3.0*TMath::Pi()/2.0 , (double) multbinsize-1 };
THnSparseD *hist_reco_lead_reco_subl_quench2pc = new THnSparseD("hist_reco_lead_reco_subl_quench2pc", "hist_reco_lead_reco_subl_quench2pc", 4, bins4D_quenc_2pc, xmin4D_quenc_2pc, xmax4D_quenc_2pc);
THnSparseD *hist_gen_lead_gen_subl_quench2pc = new THnSparseD("hist_gen_lead_gen_subl_quench2pc", "hist_gen_lead_gen_subl_quench2pc", 4, bins4D_quenc_2pc, xmin4D_quenc_2pc, xmax4D_quenc_2pc);
THnSparseD *hist_ref_lead_ref_subl_quench2pc = new THnSparseD("hist_ref_lead_ref_subl_quench2pc", "hist_ref_lead_ref_subl_quench2pc", 4, bins4D_quenc_2pc, xmin4D_quenc_2pc, xmax4D_quenc_2pc);

// Correlation studies
// Axis : 0 -> delta phi, 1 -> delta eta, 2 -> track pT, 3 -> multiplicity
int    bins4D_jettrk[4]   =   { 200                 , 500  ,   trkbinsize-1          , multbinsize-1};
double xmin4D_jettrk[4]   =   { -TMath::Pi()/2.0    , -5.0 ,   0                     , 0};
double xmax4D_jettrk[4]   =   { 3.0*TMath::Pi()/2.0 , 5.0  ,   (double) trkbinsize-1 , (double) multbinsize-1};

// Correlation: Reco Jet + Reco Track
THnSparseD *hist_correlation_signal_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_jet_reco_track_reco","hist_correlation_signal_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_reco_track_reco = new THnSparseD("hist_correlation_rotation_jet_reco_track_reco","hist_correlation_rotation_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_reco_track_reco = new THnSparseD("hist_correlation_mixing_jet_reco_track_reco","hist_correlation_mixing_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_lead_jet_reco_track_reco","hist_correlation_signal_lead_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_rotation_lead_jet_reco_track_reco","hist_correlation_rotation_lead_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_mixing_lead_jet_reco_track_reco","hist_correlation_mixing_lead_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subl_jet_reco_track_reco","hist_correlation_signal_subl_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_rotation_subl_jet_reco_track_reco","hist_correlation_rotation_subl_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_mixing_subl_jet_reco_track_reco","hist_correlation_mixing_subl_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// Correlation: Reco Jet + Gen Track
THnSparseD *hist_correlation_signal_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_jet_reco_track_gen","hist_correlation_signal_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_reco_track_gen = new THnSparseD("hist_correlation_rotation_jet_reco_track_gen","hist_correlation_rotation_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_reco_track_gen = new THnSparseD("hist_correlation_mixing_jet_reco_track_gen","hist_correlation_mixing_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_lead_jet_reco_track_gen","hist_correlation_signal_lead_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_rotation_lead_jet_reco_track_gen","hist_correlation_rotation_lead_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_mixing_lead_jet_reco_track_gen","hist_correlation_mixing_lead_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subl_jet_reco_track_gen","hist_correlation_signal_subl_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_rotation_subl_jet_reco_track_gen","hist_correlation_rotation_subl_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_mixing_subl_jet_reco_track_gen","hist_correlation_mixing_subl_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// Correlation: Gen Jet + Reco Track
THnSparseD *hist_correlation_signal_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_jet_gen_track_reco","hist_correlation_signal_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_gen_track_reco = new THnSparseD("hist_correlation_rotation_jet_gen_track_reco","hist_correlation_rotation_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_gen_track_reco = new THnSparseD("hist_correlation_mixing_jet_gen_track_reco","hist_correlation_mixing_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_lead_jet_gen_track_reco","hist_correlation_signal_lead_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_rotation_lead_jet_gen_track_reco","hist_correlation_rotation_lead_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_mixing_lead_jet_gen_track_reco","hist_correlation_mixing_lead_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subl_jet_gen_track_reco","hist_correlation_signal_subl_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_rotation_subl_jet_gen_track_reco","hist_correlation_rotation_subl_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_mixing_subl_jet_gen_track_reco","hist_correlation_mixing_subl_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// Correlation: Gen Jet + Gen Track
THnSparseD *hist_correlation_signal_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_jet_gen_track_gen","hist_correlation_signal_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_gen_track_gen = new THnSparseD("hist_correlation_rotation_jet_gen_track_gen","hist_correlation_rotation_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_gen_track_gen = new THnSparseD("hist_correlation_mixing_jet_gen_track_gen","hist_correlation_mixing_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_lead_jet_gen_track_gen","hist_correlation_signal_lead_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_rotation_lead_jet_gen_track_gen","hist_correlation_rotation_lead_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_mixing_lead_jet_gen_track_gen","hist_correlation_mixing_lead_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subl_jet_gen_track_gen","hist_correlation_signal_subl_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_rotation_subl_jet_gen_track_gen","hist_correlation_rotation_subl_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_mixing_subl_jet_gen_track_gen","hist_correlation_mixing_subl_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// Correlation: Ref Jet + Ref Track
THnSparseD *hist_correlation_signal_jet_ref_track_ref = new THnSparseD("hist_correlation_signal_jet_ref_track_ref","hist_correlation_signal_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_ref_track_ref = new THnSparseD("hist_correlation_rotation_jet_ref_track_ref","hist_correlation_rotation_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_ref_track_ref = new THnSparseD("hist_correlation_mixing_jet_ref_track_ref","hist_correlation_mixing_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_ref_track_ref = new THnSparseD("hist_correlation_signal_lead_jet_ref_track_ref","hist_correlation_signal_lead_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_ref_track_ref = new THnSparseD("hist_correlation_rotation_lead_jet_ref_track_ref","hist_correlation_rotation_lead_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_ref_track_ref = new THnSparseD("hist_correlation_mixing_lead_jet_ref_track_ref","hist_correlation_mixing_lead_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_ref_track_ref = new THnSparseD("hist_correlation_signal_subl_jet_ref_track_ref","hist_correlation_signal_subl_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_ref_track_ref = new THnSparseD("hist_correlation_rotation_subl_jet_ref_track_ref","hist_correlation_rotation_subl_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_ref_track_ref = new THnSparseD("hist_correlation_mixing_subl_jet_ref_track_ref","hist_correlation_mixing_subl_jet_ref_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// Correlation: Ref Jet + Gen Track
THnSparseD *hist_correlation_signal_jet_ref_track_gen = new THnSparseD("hist_correlation_signal_jet_ref_track_gen","hist_correlation_signal_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_ref_track_gen = new THnSparseD("hist_correlation_rotation_jet_ref_track_gen","hist_correlation_rotation_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_ref_track_gen = new THnSparseD("hist_correlation_mixing_jet_ref_track_gen","hist_correlation_mixing_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_ref_track_gen = new THnSparseD("hist_correlation_signal_lead_jet_ref_track_gen","hist_correlation_signal_lead_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_ref_track_gen = new THnSparseD("hist_correlation_rotation_lead_jet_ref_track_gen","hist_correlation_rotation_lead_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_ref_track_gen = new THnSparseD("hist_correlation_mixing_lead_jet_ref_track_gen","hist_correlation_mixing_lead_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_ref_track_gen = new THnSparseD("hist_correlation_signal_subl_jet_ref_track_gen","hist_correlation_signal_subl_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_ref_track_gen = new THnSparseD("hist_correlation_rotation_subl_jet_ref_track_gen","hist_correlation_rotation_subl_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_ref_track_gen = new THnSparseD("hist_correlation_mixing_subl_jet_ref_track_gen","hist_correlation_mixing_subl_jet_ref_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// Correlation: Gen Jet + Ref Track
THnSparseD *hist_correlation_signal_jet_gen_track_ref = new THnSparseD("hist_correlation_signal_jet_gen_track_ref","hist_correlation_signal_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_jet_gen_track_ref = new THnSparseD("hist_correlation_rotation_jet_gen_track_ref","hist_correlation_rotation_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_jet_gen_track_ref = new THnSparseD("hist_correlation_mixing_jet_gen_track_ref","hist_correlation_mixing_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_lead_jet_gen_track_ref = new THnSparseD("hist_correlation_signal_lead_jet_gen_track_ref","hist_correlation_signal_lead_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_lead_jet_gen_track_ref = new THnSparseD("hist_correlation_rotation_lead_jet_gen_track_ref","hist_correlation_rotation_lead_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_lead_jet_gen_track_ref = new THnSparseD("hist_correlation_mixing_lead_jet_gen_track_ref","hist_correlation_mixing_lead_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subl_jet_gen_track_ref = new THnSparseD("hist_correlation_signal_subl_jet_gen_track_ref","hist_correlation_signal_subl_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_rotation_subl_jet_gen_track_ref = new THnSparseD("hist_correlation_rotation_subl_jet_gen_track_ref","hist_correlation_rotation_subl_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_mixing_subl_jet_gen_track_ref = new THnSparseD("hist_correlation_mixing_subl_jet_gen_track_ref","hist_correlation_mixing_subl_jet_gen_track_ref",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// Correlation: Include sube > 0
THnSparseD *hist_correlation_signal_subg0_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subg0_jet_reco_track_reco","hist_correlation_signal_subg0_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subg0_jet_reco_track_gen","hist_correlation_signal_subg0_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subg0_jet_gen_track_reco","hist_correlation_signal_subg0_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subg0_jet_gen_track_gen","hist_correlation_signal_subg0_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

THnSparseD *hist_correlation_signal_subg0_lead_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subg0_lead_jet_reco_track_reco","hist_correlation_signal_subg0_lead_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subg0_lead_jet_reco_track_gen","hist_correlation_signal_subg0_lead_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subg0_lead_jet_gen_track_reco","hist_correlation_signal_subg0_lead_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_lead_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subg0_lead_jet_gen_track_gen","hist_correlation_signal_subg0_lead_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

THnSparseD *hist_correlation_signal_subg0_subl_jet_reco_track_reco = new THnSparseD("hist_correlation_signal_subg0_subl_jet_reco_track_reco","hist_correlation_signal_subg0_subl_jet_reco_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_reco_track_gen = new THnSparseD("hist_correlation_signal_subg0_subl_jet_reco_track_gen","hist_correlation_signal_subg0_subl_jet_reco_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_gen_track_reco = new THnSparseD("hist_correlation_signal_subg0_subl_jet_gen_track_reco","hist_correlation_signal_subg0_subl_jet_gen_track_reco",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);
THnSparseD *hist_correlation_signal_subg0_subl_jet_gen_track_gen = new THnSparseD("hist_correlation_signal_subg0_subl_jet_gen_track_gen","hist_correlation_signal_subg0_subl_jet_gen_track_gen",4,bins4D_jettrk,xmin4D_jettrk,xmax4D_jettrk);

// histograms for matched jets and parton flavor studies
TH1D *hist_matched_jet_weighted_nocut = new TH1D("hist_matched_jet_weighted_nocut", "hist_matched_jet_weighted_nocut", 100, 0.0, 500.0);
THnSparseD *hist_matched_jet = new THnSparseD("hist_matched_jet", "hist_matched_jet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_matched_jet_weighted = new THnSparseD("hist_matched_jet_weighted", "hist_matched_jet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

int    bins4D_jetflavor[4]   =   { 200  ,  50  , 8, multbinsize-1};
double xmin4D_jetflavor[4]   =   { 0    , -2.5 , 0, 0};
double xmax4D_jetflavor[4]   =   { 1000 ,  2.5 , 8, (double) multbinsize-1};
THnSparseD *hist_matched_jet_parton = new THnSparseD("hist_matched_jet_parton", "hist_matched_jet_parton", 4, bins4D_jetflavor, xmin4D_jetflavor, xmax4D_jetflavor);
THnSparseD *hist_matched_jet_parton_fromB = new THnSparseD("hist_matched_jet_parton_fromB", "hist_matched_jet_parton_fromB", 4, bins4D_jetflavor, xmin4D_jetflavor, xmax4D_jetflavor);

//histograms for Jet Energy Scale (JES)
//double jetptbin[18] = {30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,140.0,160.0,180.0,210.0,250.0,300.0,380.0,500.0}; //just to remember
int    bins4D_jes[4]   =   { 500  ,  50  , 8, multbinsize-1};
double xmin4D_jes[4]   =   { 0.0  ,  0   , 0, 0};
double xmax4D_jes[4]   =   { 5.0  ,  1000, 8, (double) multbinsize-1};
THnSparseD *hist_jes_raw = new THnSparseD("hist_jes_raw", "hist_jes_raw", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD *hist_jes_raw_weighted = new THnSparseD("hist_jes_raw_weighted", "hist_jes_raw_weighted", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD *hist_jes_reco = new THnSparseD("hist_jes_reco", "hist_jes_reco", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD *hist_jes_reco_weighted = new THnSparseD("hist_jes_reco_weighted", "hist_jes_reco_weighted", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD *hist_jes_raw_fromB = new THnSparseD("hist_jes_raw_fromB", "hist_jes_raw_fromB", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD *hist_jes_raw_fromB_weighted = new THnSparseD("hist_jes_raw_fromB_weighted", "hist_jes_raw_fromB_weighted", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD *hist_jes_reco_fromB = new THnSparseD("hist_jes_reco_fromB", "hist_jes_reco_fromB", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);
THnSparseD *hist_jes_reco_fromB_weighted = new THnSparseD("hist_jes_reco_fromB_weighted", "hist_jes_reco_fromB_weighted", 4, bins4D_jes, xmin4D_jes, xmax4D_jes);

//histograms for Jet Energy Resolution (JER)
//double jetptbin[18] = {30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,140.0,160.0,180.0,210.0,250.0,300.0,380.0,500.0}; //just to remember
int    bins4D_jer[4]   =   { 400  ,  50  , 8, multbinsize-1};
double xmin4D_jer[4]   =   { -2.0 ,  0   , 0, 0};
double xmax4D_jer[4]   =   { 2.0  ,  1000, 8, (double) multbinsize-1};
THnSparseD *hist_jer_raw = new THnSparseD("hist_jer_raw", "hist_jer_raw", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD *hist_jer_raw_weighted = new THnSparseD("hist_jer_raw_weighted", "hist_jer_raw_weighted", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD *hist_jer_reco = new THnSparseD("hist_jer_reco", "hist_jer_reco", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD *hist_jer_reco_weighted = new THnSparseD("hist_jer_reco_weighted", "hist_jer_reco_weighted", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD *hist_jer_raw_fromB = new THnSparseD("hist_jer_raw_fromB", "hist_jer_raw_fromB", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD *hist_jer_raw_fromB_weighted = new THnSparseD("hist_jer_raw_fromB_weighted", "hist_jer_raw_fromB_weighted", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD *hist_jer_reco_fromB = new THnSparseD("hist_jer_reco_fromB", "hist_jer_reco_fromB", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);
THnSparseD *hist_jer_reco_fromB_weighted = new THnSparseD("hist_jer_reco_fromB_weighted", "hist_jer_reco_fromB_weighted", 4, bins4D_jer, xmin4D_jer, xmax4D_jer);

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
    vzhist->Sumw2();
   	vzhist_weighted->Sumw2();
    pthathist->Sumw2();
    pthathist_weighted->Sumw2();
	hist_reco_trk->Sumw2();
	hist_reco_trk_corr->Sumw2();
	hist_reco_trk_weighted->Sumw2();
	hist_gen_trk->Sumw2();
	hist_gen_trk_weighted->Sumw2();
	hist_reco_jet_weighted_nocut->Sumw2();
    hist_reco_jet->Sumw2();
    hist_reco_jet_corr->Sumw2();
    hist_reco_jet_weighted->Sumw2();
    hist_reco_jet_corr_weighted->Sumw2();
    hist_reco_leadjet_pt_nocut->Sumw2();
    hist_reco_leadjet_pt_nocut_weighted->Sumw2();
    hist_reco_subljet_pt_nocut->Sumw2();
    hist_reco_subljet_pt_nocut_weighted->Sumw2();
    hist_reco_leadjet->Sumw2();
    hist_reco_leadjet_weighted->Sumw2();
    hist_reco_subljet->Sumw2();
    hist_reco_subljet_weighted->Sumw2();
	hist_gen_jet_weighted_nocut->Sumw2();
	hist_gen_jet->Sumw2();
	hist_gen_jet_weighted->Sumw2();
	hist_gen_leadjet_pt_nocut->Sumw2();
	hist_gen_leadjet_pt_nocut_weighted->Sumw2();
	hist_gen_subljet_pt_nocut->Sumw2();
	hist_gen_subljet_pt_nocut_weighted->Sumw2();
	hist_gen_leadjet->Sumw2();
	hist_gen_leadjet_weighted->Sumw2();
	hist_gen_subljet->Sumw2();
	hist_gen_subljet_weighted->Sumw2();
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
    hist_reco_lead_reco_subl_quench->Sumw2(); 
    hist_reco_lead_reco_subl_quench2pc->Sumw2();
   	hist_gen_lead_gen_subl_quench->Sumw2(); 
   	hist_gen_lead_gen_subl_quench2pc->Sumw2();
	hist_ref_lead_ref_subl_quench->Sumw2();
	hist_ref_lead_ref_subl_quench2pc->Sumw2();
	hist_matched_jet_weighted_nocut->Sumw2();
	hist_matched_jet->Sumw2();
	hist_matched_jet_weighted->Sumw2();
	hist_matched_jet_parton->Sumw2();
	hist_matched_jet_parton_fromB->Sumw2();
	hist_jes_raw->Sumw2();
	hist_jes_raw_weighted->Sumw2();
	hist_jes_reco->Sumw2();
	hist_jes_reco_weighted->Sumw2();
	hist_jer_raw->Sumw2();
	hist_jer_raw_weighted->Sumw2();
	hist_jer_reco->Sumw2();
	hist_jer_reco_weighted->Sumw2();
	hist_jes_raw_fromB->Sumw2();
	hist_jes_raw_fromB_weighted->Sumw2();
	hist_jes_reco_fromB->Sumw2();
	hist_jes_reco_fromB_weighted->Sumw2();
	hist_jer_raw_fromB->Sumw2();
	hist_jer_raw_fromB_weighted->Sumw2();
	hist_jer_reco_fromB->Sumw2();
	hist_jer_reco_fromB_weighted->Sumw2();

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
	if(doleadsubl){
		Nev_recoreco_lead->Write();
		Nev_recoreco_subl->Write();
	}
	if(isMC){
		Nev_recogen->Write();
		Nev_genreco->Write();
		Nev_gengen->Write();
		if(doleadsubl){
			Nev_recogen_lead->Write();
			Nev_genreco_lead->Write();
			Nev_gengen_lead->Write();
			Nev_recogen_subl->Write();
			Nev_genreco_subl->Write();
			Nev_gengen_subl->Write();
		}
	}
    multiplicity->Write();
   	multiplicity_weighted->Write();
    vzhist->Write();
   	vzhist_weighted->Write();
    if(isMC){
        pthathist->Write();
        pthathist_weighted->Write();
    }
    //tracks 
    //reco
	hist_reco_trk->Write();
	hist_reco_trk_corr->Write();
	hist_reco_trk_weighted->Write();

    //gen
	if(isMC){
		hist_gen_trk->Write();
		hist_gen_trk_weighted->Write();
	}
    //jets 
    //reco
	hist_reco_jet_weighted_nocut->Write();
    hist_reco_jet->Write();
    hist_reco_jet_corr->Write();
    hist_reco_jet_weighted->Write();
    hist_reco_jet_corr_weighted->Write();

    if(doleadsubl){
    	hist_reco_leadjet_pt_nocut->Write();
    	hist_reco_leadjet_pt_nocut_weighted->Write();
    	hist_reco_subljet_pt_nocut->Write();
    	hist_reco_subljet_pt_nocut_weighted->Write();
    	hist_reco_leadjet->Write();
    	hist_reco_leadjet_weighted->Write();
    	hist_reco_subljet->Write();
    	hist_reco_subljet_weighted->Write();
	}
	if(isMC){
		hist_matched_jet_weighted_nocut->Write();
		hist_matched_jet->Write();
		hist_matched_jet_weighted->Write();
		hist_matched_jet_parton->Write();
		hist_matched_jet_parton_fromB->Write();
		if(doleadsubl){

		}
		hist_gen_jet_weighted_nocut->Write();
		hist_gen_jet->Write();
		hist_gen_jet_weighted->Write();
		if(doleadsubl){
			hist_gen_leadjet_pt_nocut->Write();
			hist_gen_leadjet_pt_nocut_weighted->Write();
			hist_gen_subljet_pt_nocut->Write();
			hist_gen_subljet_pt_nocut_weighted->Write();
			hist_gen_leadjet->Write();
			hist_gen_leadjet_weighted->Write();
			hist_gen_subljet->Write();
			hist_gen_subljet_weighted->Write();
		}
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
void w_jetquenching_hist(bool isMC){
    hist_reco_lead_reco_subl_quench->Write(); 
    hist_reco_lead_reco_subl_quench2pc->Write();
    if(isMC){
    	hist_gen_lead_gen_subl_quench->Write(); 
    	hist_gen_lead_gen_subl_quench2pc->Write();
    	hist_ref_lead_ref_subl_quench->Write();
		hist_ref_lead_ref_subl_quench2pc->Write();
    }
}

// JES and JER histograms
void w_jes_jer_hist(){ 
	hist_jes_raw->Write();
	hist_jes_raw_weighted->Write();
	hist_jes_reco->Write();
	hist_jes_reco_weighted->Write();
	hist_jer_raw->Write();
	hist_jer_raw_weighted->Write();
	hist_jer_reco->Write();
	hist_jer_reco_weighted->Write();
	hist_jes_raw_fromB->Write();
	hist_jes_raw_fromB_weighted->Write();
	hist_jes_reco_fromB->Write();
	hist_jes_reco_fromB_weighted->Write();
	hist_jer_raw_fromB->Write();
	hist_jer_raw_fromB_weighted->Write();
	hist_jer_reco_fromB->Write();
	hist_jer_reco_fromB_weighted->Write();
}


// Unfolding histograms
void w_unfold_hist(){ 

}
