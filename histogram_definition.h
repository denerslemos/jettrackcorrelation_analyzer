#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // call inputs

int trkbinsize = (int) trk_pt_bins.size(); // track bins for jet-track correlation
int multbinsize = (int) multiplicity_centrality_bins.size();// multiplicity or centrality bins for jet-track correlation
int ptavebinsize = (int) pt_ave_bins.size();// multiplicity or centrality bins for jet-track correlation


// -------------------------------- QA plots --------------------------------
// Event quantities
//number of events
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


TH1I *NJets = new TH1I("NJets", "NJets", 20, 0, 20);
TH1I *NJetsSub = new TH1I("NJetsSub", "NJetsSub", 20, 0, 20);
TH1I *NJetsLead = new TH1I("NJetsLead", "NJetsLead", 20, 0, 20);
TH1I *NJetsLJSLJ = new TH1I("NJetsLJSLJ", "NJetsLJSLJ", 20, 0, 20);

// multiplicity
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 80, 0.0, 400.0);
TH1D *multiplicity_weighted = new TH1D("multiplicity_weighted", "multiplicity_weighted", 80, 0.0, 400.0);

TH1D *reco_mult = new TH1D("reco_mult", "reco_mult", 80, 0.0, 400.0);
TH1D *reco_mult_weighted = new TH1D("reco_mult_weighted", "reco_mult_weighted", 80, 0.0, 400.0);
TH1D *gen_mult = new TH1D("gen_mult", "gen_mult", 80, 0.0, 400.0);
TH1D *gen_mult_weighted = new TH1D("gen_mult_weighted", "gen_mult_weighted", 80, 0.0, 400.0);

// Z vertex
TH2D *vzhist = new TH2D("vzhist", "vzhist", 60, -15.5, 15.5, 80, 0.0, 400.0);
TH2D *vzhist_weighted = new TH2D("vzhist_weighted", "vzhist_weighted", 60, -15.5, 15.5, 80, 0.0, 400.0);
TH2D *vzhist_jet_weighted = new TH2D("vzhist_jet_weighted", "vzhist_jet_weighted", 60, -15.5, 15.5, 80, 0.0, 400.0);
TH2D *vzhist_dijet_weighted = new TH2D("vzhist_dijet_weighted", "vzhist_dijet_weighted", 60, -15.5, 15.5, 80, 0.0, 400.0);

// pthat
TH2D *pthathist = new TH2D("pthathist", "pthathist", 100, 0, 1000, 80, 0.0, 400.0);
TH2D *pthathist_weighted = new TH2D("pthathist_weighted", "pthathist_weighted", 100, 0, 1000, 80, 0.0, 400.0);

// trackmax histogram
TH2D *trackmaxptinjethisto = new TH2D("trackmaxptinjethisto", "trackmaxptinjethisto", 1000, 0, 1000, 80, 0.0, 400.0);

// Axis : 0 -> HF+, 1 -> HF-, 2 -> multbin
int	bins3D_HF[3]   =   { 200  ,  200 ,  80};
double xmin3D_HF[3]   =   { 0.0  ,  0.0 ,  0.0};
double xmax3D_HF[3]   =   { 400  ,  400 ,  400};

THnSparseD *hfhist = new THnSparseD("hfhist", "hfhist", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);
THnSparseD *hfhist_weighted = new THnSparseD("hfhist_weighted", "hfhist_weighted", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);
THnSparseD *hfhistEta4 = new THnSparseD("hfhistEta4", "hfhistEta4", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);
THnSparseD *hfhistEta4_weighted = new THnSparseD("hfhistEta4_weighted", "hfhistEta4_weighted", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);
THnSparseD *hfhist_onejet_weighted = new THnSparseD("hfhist_onejet_weighted", "hfhist_onejet_weighted", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);
THnSparseD *hfhistEta4_onejet_weighted = new THnSparseD("hfhistEta4_onejet_weighted", "hfhistEta4_onejet_weighted", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);
THnSparseD *hfhist_dijet_weighted = new THnSparseD("hfhist_dijet_weighted", "hfhist_dijet_weighted", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);
THnSparseD *hfhistEta4_dijet_weighted = new THnSparseD("hfhistEta4_dijet_weighted", "hfhistEta4_dijet_weighted", 3, bins3D_HF, xmin3D_HF, xmax3D_HF);

// Axis : 0 -> ZDC+, 1 -> ZDC-, 2 -> multbin
int	bins3D_ZDC[3]   =   {  200   ,   200   ,  80};
double xmin3D_ZDC[3]   =   { -10000 ,  -10000 ,  0.0};
double xmax3D_ZDC[3]   =   {  10000 ,   10000 ,  400};

THnSparseD *zdchist = new THnSparseD("zdchist", "zdchist", 3, bins3D_ZDC, xmin3D_ZDC, xmax3D_ZDC);
THnSparseD *zdchist_weighted = new THnSparseD("zdchist_weighted", "zdchist_weighted", 3, bins3D_ZDC, xmin3D_ZDC, xmax3D_ZDC);
THnSparseD *zdchist_onejet_weighted = new THnSparseD("zdchist_onejet_weighted", "zdchist_onejet_weighted", 3, bins3D_ZDC, xmin3D_ZDC, xmax3D_ZDC);
THnSparseD *zdchist_dijet_weighted = new THnSparseD("zdchist_dijet_weighted", "zdchist_dijet_weighted", 3, bins3D_ZDC, xmin3D_ZDC, xmax3D_ZDC);


//quantities with at least 1 jet
TH1D *multiplicity_withonejet40 = new TH1D("multiplicity_withonejet40", "multiplicity_withonejet40", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet50 = new TH1D("multiplicity_withonejet50", "multiplicity_withonejet50", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet60 = new TH1D("multiplicity_withonejet60", "multiplicity_withonejet60", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet80 = new TH1D("multiplicity_withonejet80", "multiplicity_withonejet80", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet100 = new TH1D("multiplicity_withonejet100", "multiplicity_withonejet100", 80, 0.0, 400.0);

//quantities with at least 1 jet pT > jet min pT in GeV
TH1D *multiplicity_withonejet = new TH1D("multiplicity_withonejet", "multiplicity_withonejet", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet_weighted = new TH1D("multiplicity_withonejet_weighted", "multiplicity_withonejet_weighted", 80, 0.0, 400.0);
TH1D *reco_mult_withonejet = new TH1D("reco_mult_withonejet", "reco_mult_withonejet", 80, 0.0, 400.0);
TH1D *reco_mult_withonejet_weighted = new TH1D("reco_mult_withonejet_weighted", "reco_mult_withonejet_weighted", 80, 0.0, 400.0);
TH1D *gen_mult_withonejet = new TH1D("gen_mult_withonejet", "gen_mult_withonejet", 80, 0.0, 400.0);
TH1D *gen_mult_withonejet_weighted = new TH1D("gen_mult_withonejet_weighted", "gen_mult_withonejet_weighted", 80, 0.0, 400.0);

//quantities with dijets
// --> multiplicity
TH1D *multiplicity_withdijets = new TH1D("multiplicity_withdijets", "multiplicity_withdijets", 80, 0.0, 400.0);
TH1D *multiplicity_withdijets_weighted = new TH1D("multiplicity_withdijets_weighted", "multiplicity_withdijets_weighted", 80, 0.0, 400.0);
TH1D *reco_mult_withdijets = new TH1D("reco_mult_withdijets", "reco_mult_withdijets", 80, 0.0, 400.0);
TH1D *reco_mult_withdijets_weighted = new TH1D("reco_mult_withdijets_weighted", "reco_mult_withdijets_weighted", 80, 0.0, 400.0);
TH1D *gen_mult_withdijets = new TH1D("gen_mult_withdijets", "gen_mult_withdijets", 80, 0.0, 400.0);
TH1D *gen_mult_withdijets_weighted = new TH1D("gen_mult_withdijets_weighted", "gen_mult_withdijets_weighted", 80, 0.0, 400.0);

// Axis : 0 -> delta R, 1 -> gen jet pT, 2 -> flavour, 3 -> event multiplicity --> for future analysis
/*
int	bins4D_jetaxis[4]   =   { 50  ,  100  , 8, multbinsize-1};
double xmin4D_jetaxis[4]   =   { 0.0 ,   0   , 0, 0};
double xmax4D_jetaxis[4]   =   { 1.0  ,  1000, 8, (double) multbinsize-1};
THnSparseD *genjetaxischeck = new THnSparseD("genjetaxischeck", "genjetaxischeck", 4, bins4D_jetaxis, xmin4D_jetaxis, xmax4D_jetaxis);
*/

// event plane histograms

// Axis : 0 -> EP multiplicity, 1 -> qvector, 2 -> PsiEP, 3 -> event multiplicity
int	bins4D_EP[4]   =   { 50	,  200 ,   64		   , 100};
double xmin4D_EP[4]   =   { 0	 ,	0 ,   -TMath::Pi() , 0};
double xmax4D_EP[4]   =   { 500   ,  100 ,   TMath::Pi()  , 500};

//after flattening
THnSparseD *EP2_plus_flat = new THnSparseD("EP2_plus_flat", "EP2_plus_flat", 4, bins4D_EP, xmin4D_EP, xmax4D_EP);
THnSparseD *EP2_minus_flat = new THnSparseD("EP2_minus_flat", "EP2_minus_flat", 4, bins4D_EP, xmin4D_EP, xmax4D_EP);
THnSparseD *EP3_plus_flat = new THnSparseD("EP3_plus_flat", "EP3_plus_flat", 4, bins4D_EP, xmin4D_EP, xmax4D_EP);
THnSparseD *EP3_minus_flat = new THnSparseD("EP3_minus_flat", "EP3_minus_flat", 4, bins4D_EP, xmin4D_EP, xmax4D_EP);
THnSparseD *EP4_plus_flat = new THnSparseD("EP4_plus_flat", "EP4_plus_flat", 4, bins4D_EP, xmin4D_EP, xmax4D_EP);
THnSparseD *EP4_minus_flat = new THnSparseD("EP4_minus_flat", "EP4_minus_flat", 4, bins4D_EP, xmin4D_EP, xmax4D_EP);

// Axis : 0 -> Dphi track-EP, 1 -> trkbin, 2 -> multbin
int	bins3D_TRKEP[3]   =   { 100				,  trkbinsize-1 		,  multbinsize-1};
double xmin3D_TRKEP[3]   =   { -TMath::Pi()/2.0	 ,  0 					, 0};
double xmax3D_TRKEP[3]   =   { 3.0*TMath::Pi()/2.0  ,  (double)trkbinsize-1 , (double) multbinsize-1};

THnSparseD *Dphi_EP2_flat_trk_minus = new THnSparseD("Dphi_EP2_flat_trk_minus", "Dphi_EP2_flat_trk_minus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_EP2_flat_trk_plus = new THnSparseD("Dphi_EP2_flat_trk_plus", "Dphi_EP2_flat_trk_plus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_EP3_flat_trk_minus = new THnSparseD("Dphi_EP3_flat_trk_minus", "Dphi_EP3_flat_trk_minus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_EP3_flat_trk_plus = new THnSparseD("Dphi_EP3_flat_trk_plus", "Dphi_EP3_flat_trk_plus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_EP4_flat_trk_minus = new THnSparseD("Dphi_EP4_flat_trk_minus", "Dphi_EP4_flat_trk_minus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_EP4_flat_trk_plus = new THnSparseD("Dphi_EP4_flat_trk_plus", "Dphi_EP4_flat_trk_plus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);

THnSparseD *Dphi_GEN_EP2_flat_trk_minus = new THnSparseD("Dphi_GEN_EP2_flat_trk_minus", "Dphi_GEN_EP2_flat_trk_minus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_GEN_EP2_flat_trk_plus = new THnSparseD("Dphi_GEN_EP2_flat_trk_plus", "Dphi_GEN_EP2_flat_trk_plus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_GEN_EP3_flat_trk_minus = new THnSparseD("Dphi_GEN_EP3_flat_trk_minus", "Dphi_GEN_EP3_flat_trk_minus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_GEN_EP3_flat_trk_plus = new THnSparseD("Dphi_GEN_EP3_flat_trk_plus", "Dphi_GEN_EP3_flat_trk_plus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_GEN_EP4_flat_trk_minus = new THnSparseD("Dphi_GEN_EP4_flat_trk_minus", "Dphi_GEN_EP4_flat_trk_minus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);
THnSparseD *Dphi_GEN_EP4_flat_trk_plus = new THnSparseD("Dphi_GEN_EP4_flat_trk_plus", "Dphi_GEN_EP4_flat_trk_plus", 3, bins3D_TRKEP, xmin3D_TRKEP, xmax3D_TRKEP);

//correlations to EP
// Axis : X -> delta phi, Y -> multiplicity bins

TH2D *Dphi_flat_EP2_inclusive_minus = new TH2D("Dphi_flat_EP2_inclusive_minus", "Dphi_flat_EP2_inclusive_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP2_leading_minus = new TH2D("Dphi_flat_EP2_leading_minus", "Dphi_flat_EP2_leading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP2_subleading_minus = new TH2D("Dphi_flat_EP2_subleading_minus", "Dphi_flat_EP2_subleading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP2_inclusive_plus = new TH2D("Dphi_flat_EP2_inclusive_plus", "Dphi_flat_EP2_inclusive_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP2_leading_plus = new TH2D("Dphi_flat_EP2_leading_plus", "Dphi_flat_EP2_leading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP2_subleading_plus = new TH2D("Dphi_flat_EP2_subleading_plus", "Dphi_flat_EP2_subleading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);

TH2D *Dphi_flat_EP3_inclusive_minus = new TH2D("Dphi_flat_EP3_inclusive_minus", "Dphi_flat_EP3_inclusive_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP3_leading_minus = new TH2D("Dphi_flat_EP3_leading_minus", "Dphi_flat_EP3_leading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP3_subleading_minus = new TH2D("Dphi_flat_EP3_subleading_minus", "Dphi_flat_EP3_subleading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP3_inclusive_plus = new TH2D("Dphi_flat_EP3_inclusive_plus", "Dphi_flat_EP3_inclusive_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP3_leading_plus = new TH2D("Dphi_flat_EP3_leading_plus", "Dphi_flat_EP3_leading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP3_subleading_plus = new TH2D("Dphi_flat_EP3_subleading_plus", "Dphi_flat_EP3_subleading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);

TH2D *Dphi_flat_EP4_inclusive_minus = new TH2D("Dphi_flat_EP4_inclusive_minus", "Dphi_flat_EP4_inclusive_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP4_leading_minus = new TH2D("Dphi_flat_EP4_leading_minus", "Dphi_flat_EP4_leading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP4_subleading_minus = new TH2D("Dphi_flat_EP4_subleading_minus", "Dphi_flat_EP4_subleading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP4_inclusive_plus = new TH2D("Dphi_flat_EP4_inclusive_plus", "Dphi_flat_EP4_inclusive_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP4_leading_plus = new TH2D("Dphi_flat_EP4_leading_plus", "Dphi_flat_EP4_leading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_flat_EP4_subleading_plus = new TH2D("Dphi_flat_EP4_subleading_plus", "Dphi_flat_EP4_subleading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);

TH2D *Dphi_GEN_flat_EP2_inclusive_minus = new TH2D("Dphi_GEN_flat_EP2_inclusive_minus", "Dphi_GEN_flat_EP2_inclusive_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP2_leading_minus = new TH2D("Dphi_GEN_flat_EP2_leading_minus", "Dphi_GEN_flat_EP2_leading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP2_subleading_minus = new TH2D("Dphi_GEN_flat_EP2_subleading_minus", "Dphi_GEN_flat_EP2_subleading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP2_inclusive_plus = new TH2D("Dphi_GEN_flat_EP2_inclusive_plus", "Dphi_GEN_flat_EP2_inclusive_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP2_leading_plus = new TH2D("Dphi_GEN_flat_EP2_leading_plus", "Dphi_GEN_flat_EP2_leading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP2_subleading_plus = new TH2D("Dphi_GEN_flat_EP2_subleading_plus", "Dphi_GEN_flat_EP2_subleading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);

TH2D *Dphi_GEN_flat_EP3_inclusive_minus = new TH2D("Dphi_GEN_flat_EP3_inclusive_minus", "Dphi_GEN_flat_EP3_inclusive_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP3_leading_minus = new TH2D("Dphi_GEN_flat_EP3_leading_minus", "Dphi_GEN_flat_EP3_leading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP3_subleading_minus = new TH2D("Dphi_GEN_flat_EP3_subleading_minus", "Dphi_GEN_flat_EP3_subleading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP3_inclusive_plus = new TH2D("Dphi_GEN_flat_EP3_inclusive_plus", "Dphi_GEN_flat_EP3_inclusive_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP3_leading_plus = new TH2D("Dphi_GEN_flat_EP3_leading_plus", "Dphi_GEN_flat_EP3_leading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP3_subleading_plus = new TH2D("Dphi_GEN_flat_EP3_subleading_plus", "Dphi_GEN_flat_EP3_subleading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);

TH2D *Dphi_GEN_flat_EP4_inclusive_minus = new TH2D("Dphi_GEN_flat_EP4_inclusive_minus", "Dphi_GEN_flat_EP4_inclusive_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP4_leading_minus = new TH2D("Dphi_GEN_flat_EP4_leading_minus", "Dphi_GEN_flat_EP4_leading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP4_subleading_minus = new TH2D("Dphi_GEN_flat_EP4_subleading_minus", "Dphi_GEN_flat_EP4_subleading_minus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP4_inclusive_plus = new TH2D("Dphi_GEN_flat_EP4_inclusive_plus", "Dphi_GEN_flat_EP4_inclusive_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP4_leading_plus = new TH2D("Dphi_GEN_flat_EP4_leading_plus", "Dphi_GEN_flat_EP4_leading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);
TH2D *Dphi_GEN_flat_EP4_subleading_plus = new TH2D("Dphi_GEN_flat_EP4_subleading_plus", "Dphi_GEN_flat_EP4_subleading_plus", 200, -TMath::Pi()/2.0, 3.0*TMath::Pi()/2.0, multbinsize-1, 0, (double) multbinsize-1);


// Track/Particle histograms
int	bins4D_trk[4]   =   { 100   ,  60  ,   64		   , multbinsize-1};
double xmin4D_trk[4]   =   { 0.0   , -3.0 ,   -TMath::Pi() , 0};
double xmax4D_trk[4]   =   { 50.0  ,  3.0 ,   TMath::Pi()  , (double) multbinsize-1};

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
int	bins4D_jet[4]   =   { 200	 ,  80  ,   64		   , multbinsize-1};
double xmin4D_jet[4]   =   { 0.0	 , -4.0 ,   -TMath::Pi() , 0};
double xmax4D_jet[4]   =   { 1000.0  ,  4.0 ,   TMath::Pi() , (double) multbinsize-1};

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


// extra jet histograms for ET4
// Axis : 0 -> jet pT, 1 -> jet eta, 2 -> jet phi, 3 -> ET4Sum
int	bins4D_jetET4[4]   =   { 100	 ,  40  ,   32		   , 80};
double xmin4D_jetET4[4]   =   { 0.0	 , -4.0 ,   -TMath::Pi() , 0};
double xmax4D_jetET4[4]   =   { 1000.0  ,  4.0 ,   TMath::Pi() ,  400.0};
THnSparseD *hist_reco_jetET = new THnSparseD("hist_reco_jetET", "hist_reco_jetET", 4, bins4D_jetET4, xmin4D_jetET4, xmax4D_jetET4);
THnSparseD *hist_ref_jetET = new THnSparseD("hist_ref_jetET", "hist_ref_jetET", 4, bins4D_jetET4, xmin4D_jetET4, xmax4D_jetET4);
THnSparseD *hist_reco_LjetET = new THnSparseD("hist_reco_LjetET", "hist_reco_LjetET", 4, bins4D_jetET4, xmin4D_jetET4, xmax4D_jetET4);
THnSparseD *hist_reco_SLjetET = new THnSparseD("hist_reco_SLjetET", "hist_reco_SLjetET", 4, bins4D_jetET4, xmin4D_jetET4, xmax4D_jetET4);
THnSparseD *hist_gen_jetET = new THnSparseD("hist_gen_jetET", "hist_gen_jetET", 4, bins4D_jetET4, xmin4D_jetET4, xmax4D_jetET4);
THnSparseD *hist_gen_LjetET = new THnSparseD("hist_gen_LjetET", "hist_gen_LjetET", 4, bins4D_jetET4, xmin4D_jetET4, xmax4D_jetET4);
THnSparseD *hist_gen_SLjetET = new THnSparseD("hist_gen_SLjetET", "hist_gen_SLjetET", 4, bins4D_jetET4, xmin4D_jetET4, xmax4D_jetET4);

// Axis : 0 -> eta dijet, 1 -> HF+, 2 -> HF-, 3 -> multiplicity
int	bins4D_dijetEpEPb[4]   =   {  60	,  80	,   80	, multbinsize-1};
double xmin4D_dijetEpEPb[4]   =   { -6.0   ,  0.0   ,   0.0   , 0};
double xmax4D_dijetEpEPb[4]   =   {  6.0   ,  400.0 ,   400.0 , (double) multbinsize-1};
THnSparseD *hist_reco_jet_dijetEpEPb = new THnSparseD("hist_reco_jet_dijetEpEPb", "hist_reco_jet_dijetEpEPb", 4, bins4D_dijetEpEPb, xmin4D_dijetEpEPb, xmax4D_dijetEpEPb);
THnSparseD *hist_gen_jet_dijetEpEPb = new THnSparseD("hist_gen_jet_dijetEpEPb", "hist_gen_jet_dijetEpEPb", 4, bins4D_dijetEpEPb, xmin4D_dijetEpEPb, xmax4D_dijetEpEPb);

// Axis : 0 -> pTave, 1 -> delta phi, 2 -> ET or ETEta4, 3 -> multiplicity
int	bins4D_avejetpt[4] =   { 100	 ,  30		  ,   80	,	 multbinsize-1};
double xmin4D_avejetpt[4] =   { 0.0	 ,  0.0		 ,   0	 ,	 0};
double xmax4D_avejetpt[4] =   { 1000.0  ,  TMath::Pi() ,   400.0 ,	 (double) multbinsize-1};
THnSparseD *hist_reco_jetavept = new THnSparseD("hist_reco_jetavept", "hist_reco_jetavept", 4, bins4D_avejetpt, xmin4D_avejetpt, xmax4D_avejetpt);
THnSparseD *hist_ref_jetavept = new THnSparseD("hist_ref_jetavept", "hist_ref_jetavept", 4, bins4D_avejetpt, xmin4D_avejetpt, xmax4D_avejetpt);
THnSparseD *hist_gen_jetavept = new THnSparseD("hist_gen_jetavept", "hist_gen_jetavept", 4, bins4D_avejetpt, xmin4D_avejetpt, xmax4D_avejetpt);


// --------------------------------------------------------------------------------------------------------
// Quenching studies
// Axis : 0 -> Aj, 1 -> Xj, 2 -> delta phi, 3 -> multiplicity
int	bins4D_quenc[4]   =   { 20   , 20	, 30		  , multbinsize-1		  };
double xmin4D_quenc[4]   =   { 0.0  , 0.0   , 0.0		 , 0					  };
double xmax4D_quenc[4]   =   { 1.0  , 1.0   , TMath::Pi() , (double) multbinsize-1 };
THnSparseD *hist_reco_lead_reco_subl_quench_mid_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_mid", "hist_reco_lead_reco_subl_quench_mid_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_mid_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_fwd", "hist_reco_lead_reco_subl_quench_mid_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_mid_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_mid_bkw", "hist_reco_lead_reco_subl_quench_mid_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_mid", "hist_reco_lead_reco_subl_quench_fwd_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_mid = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_mid", "hist_reco_lead_reco_subl_quench_bkw_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_fwd", "hist_reco_lead_reco_subl_quench_fwd_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_fwd_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_fwd_bkw", "hist_reco_lead_reco_subl_quench_fwd_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_fwd = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_fwd", "hist_reco_lead_reco_subl_quench_bkw_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_reco_lead_reco_subl_quench_bkw_bkw = new THnSparseD("hist_reco_lead_reco_subl_quench_bkw_bkw", "hist_reco_lead_reco_subl_quench_bkw_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_mid", "hist_gen_lead_gen_subl_quench_mid_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_fwd", "hist_gen_lead_gen_subl_quench_mid_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_mid_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_mid_bkw", "hist_gen_lead_gen_subl_quench_mid_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_mid", "hist_gen_lead_gen_subl_quench_fwd_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_mid = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_mid", "hist_gen_lead_gen_subl_quench_bkw_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_fwd", "hist_gen_lead_gen_subl_quench_fwd_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_fwd_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_fwd_bkw", "hist_gen_lead_gen_subl_quench_fwd_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_fwd = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_fwd", "hist_gen_lead_gen_subl_quench_bkw_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench_bkw_bkw = new THnSparseD("hist_gen_lead_gen_subl_quench_bkw_bkw", "hist_gen_lead_gen_subl_quench_bkw_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_mid", "hist_ref_lead_ref_subl_quench_mid_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_fwd", "hist_ref_lead_ref_subl_quench_mid_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_mid_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_mid_bkw", "hist_ref_lead_ref_subl_quench_mid_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_mid", "hist_ref_lead_ref_subl_quench_fwd_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_mid = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_mid", "hist_ref_lead_ref_subl_quench_bkw_mid", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_fwd", "hist_ref_lead_ref_subl_quench_fwd_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_fwd_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_fwd_bkw", "hist_ref_lead_ref_subl_quench_fwd_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_fwd = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_fwd", "hist_ref_lead_ref_subl_quench_bkw_fwd", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench_bkw_bkw = new THnSparseD("hist_ref_lead_ref_subl_quench_bkw_bkw", "hist_ref_lead_ref_subl_quench_bkw_bkw", 4, bins4D_quenc, xmin4D_quenc, xmax4D_quenc);

//Assymetry studies
// Axis : 0 -> etaDijet, 1 -> Xj, 2 -> delta phi, 3 -> ET or ETEta4, 4 -> multiplicity, 5 -> jet pT average
int	bins6D_etaDijet[6]   =   {  60   , 20	, 30		  , 80	, multbinsize-1		  , ptavebinsize-1};
double xmin6D_etaDijet[6]   =   { -6.0  , 0.0   , 0.0 	  , 0.0   , 0					  , 0};
double xmax6D_etaDijet[6]   =   {  6.0  , 1.0   , TMath::Pi() , 400.0 , (double) multbinsize-1 , (double) ptavebinsize-1};
THnSparseD *hist_etaDijet_ETEta4_reco = new THnSparseD("hist_etaDijet_ETEta4_reco", "hist_etaDijet_ETEta4_reco", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDijet_ETEta4_ref = new THnSparseD("hist_etaDijet_ETEta4_ref", "hist_etaDijet_ETEta4_ref", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDijet_ETEta4_gen = new THnSparseD("hist_etaDijet_ETEta4_gen", "hist_etaDijet_ETEta4_gen", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDijet_AJ_ETEta4_reco = new THnSparseD("hist_etaDijet_AJ_ETEta4_reco", "hist_etaDijet_AJ_ETEta4_reco", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDijet_AJ_ETEta4_ref = new THnSparseD("hist_etaDijet_AJ_ETEta4_ref", "hist_etaDijet_AJ_ETEta4_ref", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDijet_AJ_ETEta4_gen = new THnSparseD("hist_etaDijet_AJ_ETEta4_gen", "hist_etaDijet_AJ_ETEta4_gen", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDiff_ETEta4_reco = new THnSparseD("hist_etaDiff_ETEta4_reco", "hist_etaDiff_ETEta4_reco", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDiff_ETEta4_ref = new THnSparseD("hist_etaDiff_ETEta4_ref", "hist_etaDiff_ETEta4_ref", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDiff_ETEta4_gen = new THnSparseD("hist_etaDiff_ETEta4_gen", "hist_etaDiff_ETEta4_gen", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDiff_AJ_ETEta4_reco = new THnSparseD("hist_etaDiff_AJ_ETEta4_reco", "hist_etaDiff_AJ_ETEta4_reco", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDiff_AJ_ETEta4_ref = new THnSparseD("hist_etaDiff_AJ_ETEta4_ref", "hist_etaDiff_AJ_ETEta4_ref", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);
THnSparseD *hist_etaDiff_AJ_ETEta4_gen = new THnSparseD("hist_etaDiff_AJ_ETEta4_gen", "hist_etaDiff_AJ_ETEta4_gen", 6, bins6D_etaDijet, xmin6D_etaDijet, xmax6D_etaDijet);

// Axis : 0 -> in-jet multiplicity, 2 -> multiplicity
int	bins2D_injettrk[2]   =   { 100 ,  multbinsize-1			};
double xmin2D_injettrk[2]   =   { 0   , 0						};
double xmax2D_injettrk[2]   =   { 100 , (double) multbinsize-1   };
THnSparseD *hist_injet_reco_track_reco = new THnSparseD("hist_injet_reco_track_reco","hist_injet_reco_track_reco",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_injet_reco_track_gen = new THnSparseD("hist_injet_reco_track_gen","hist_injet_reco_track_gen",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_injet_gen_track_reco = new THnSparseD("hist_injet_gen_track_reco","hist_injet_gen_track_reco",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_injet_gen_track_gen = new THnSparseD("hist_injet_gen_track_gen","hist_injet_gen_track_gen",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inLeadjet_reco_track_reco = new THnSparseD("hist_inLeadjet_reco_track_reco","hist_inLeadjet_reco_track_reco",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inLeadjet_reco_track_gen = new THnSparseD("hist_inLeadjet_reco_track_gen","hist_inLeadjet_reco_track_gen",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inLeadjet_gen_track_reco = new THnSparseD("hist_inLeadjet_gen_track_reco","hist_inLeadjet_gen_track_reco",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inLeadjet_gen_track_gen = new THnSparseD("hist_inLeadjet_gen_track_gen","hist_inLeadjet_gen_track_gen",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inSubljet_reco_track_reco = new THnSparseD("hist_inSubljet_reco_track_reco","hist_inSubljet_reco_track_reco",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inSubljet_reco_track_gen = new THnSparseD("hist_inSubljet_reco_track_gen","hist_inSubljet_reco_track_gen",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inSubljet_gen_track_reco = new THnSparseD("hist_inSubljet_gen_track_reco","hist_inSubljet_gen_track_reco",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);
THnSparseD *hist_inSubljet_gen_track_gen = new THnSparseD("hist_inSubljet_gen_track_gen","hist_inSubljet_gen_track_gen",2,bins2D_injettrk,xmin2D_injettrk,xmax2D_injettrk);


// Correlation studies
// Axis : 0 -> delta phi, 1 -> delta eta, 2 -> track pT, 3 -> multiplicity
int	bins4D_jettrk[4]   =   { 200				 , 400  ,   trkbinsize-1		  , multbinsize-1			};
double xmin4D_jettrk[4]   =   { -TMath::Pi()/2.0	, -4.0 ,   0					 , 0						};
double xmax4D_jettrk[4]   =   { 3.0*TMath::Pi()/2.0 , 4.0  ,   (double) trkbinsize-1 , (double) multbinsize-1   };


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

// Two particle correlation studies
// Axis : 0 -> delta phi, 1 -> delta eta, 2 -> track pT, 3 -> multiplicity
int	bins4D_2pc[4]   =   { 40				  , 100  ,   trkbinsize-1		  , multbinsize-1};
double xmin4D_2pc[4]   =   { -TMath::Pi()/2.0	, -4.0 ,   0					 , 0};
double xmax4D_2pc[4]   =   { 3.0*TMath::Pi()/2.0 , 4.0  ,   (double) trkbinsize-1 , (double) multbinsize-1};

// 2 particle correlations for flow analysis
THnSparseD *hist_reco_reco_2pcorrelation_signal = new THnSparseD("hist_reco_reco_2pcorrelation_signal","hist_reco_reco_2pcorrelation_signal",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);
THnSparseD *hist_reco_reco_2pcorrelation_signal_subg0 = new THnSparseD("hist_reco_reco_2pcorrelation_signal_subg0","hist_reco_reco_2pcorrelation_signal_subg0",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);
THnSparseD *hist_reco_reco_2pcorrelation_signal_subcross = new THnSparseD("hist_reco_reco_2pcorrelation_signal_subcross","hist_reco_reco_2pcorrelation_signal_subcross",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);
THnSparseD *hist_reco_reco_2pcorrelation_mixing = new THnSparseD("hist_reco_reco_2pcorrelation_mixing","hist_reco_reco_2pcorrelation_mixing",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);

THnSparseD *hist_gen_gen_2pcorrelation_signal = new THnSparseD("hist_gen_gen_2pcorrelation_signal","hist_gen_gen_2pcorrelation_signal",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);
THnSparseD *hist_gen_gen_2pcorrelation_signal_subg0 = new THnSparseD("hist_gen_gen_2pcorrelation_signal_subg0","hist_gen_gen_2pcorrelation_signal_subg0",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);
THnSparseD *hist_gen_gen_2pcorrelation_signal_subcross = new THnSparseD("hist_gen_gen_2pcorrelation_signal_subcross","hist_gen_gen_2pcorrelation_signal_subcross",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);
THnSparseD *hist_gen_gen_2pcorrelation_mixing = new THnSparseD("hist_gen_gen_2pcorrelation_mixing","hist_gen_gen_2pcorrelation_mixing",4,bins4D_2pc,xmin4D_2pc,xmax4D_2pc);


// histograms for matched jets and parton flavor studies
TH1D *hist_matched_jet_weighted_nocut = new TH1D("hist_matched_jet_weighted_nocut", "hist_matched_jet_weighted_nocut", 100, 0.0, 500.0);
THnSparseD *hist_matched_jet = new THnSparseD("hist_matched_jet", "hist_matched_jet", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);
THnSparseD *hist_matched_jet_weighted = new THnSparseD("hist_matched_jet_weighted", "hist_matched_jet_weighted", 4, bins4D_jet, xmin4D_jet, xmax4D_jet);

int	bins4D_jetflavor[4]   =   { 200  ,  80  , 8, multbinsize-1};
double xmin4D_jetflavor[4]   =   { 0	, -4.0 , 0, 0};
double xmax4D_jetflavor[4]   =   { 1000 ,  4.0 , 8, (double) multbinsize-1};
THnSparseD *hist_matched_jet_parton = new THnSparseD("hist_matched_jet_parton", "hist_matched_jet_parton", 4, bins4D_jetflavor, xmin4D_jetflavor, xmax4D_jetflavor);
THnSparseD *hist_matched_jet_parton_fromB = new THnSparseD("hist_matched_jet_parton_fromB", "hist_matched_jet_parton_fromB", 4, bins4D_jetflavor, xmin4D_jetflavor, xmax4D_jetflavor);

//histograms for Jet Energy Scale (JES)
//double jetptbin[18] = {30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0,110.0,120.0,140.0,160.0,180.0,210.0,250.0,300.0,380.0,500.0}; //just to remember
int	bins5D_jes[5]   =   { 75  ,  100  , 80,  8, multbinsize-1};
double xmin5D_jes[5]   =   { 0.0  ,  0   , -4.0, 0, 0};
double xmax5D_jes[5]   =   { 3.0  ,  1000, 4.0, 8, (double) multbinsize-1};
THnSparseD *hist_jes_reco_weighted = new THnSparseD("hist_jes_reco_weighted", "hist_jes_reco_weighted", 5, bins5D_jes, xmin5D_jes, xmax5D_jes);
THnSparseD *hist_jes_reco_fromB_weighted = new THnSparseD("hist_jes_reco_fromB_weighted", "hist_jes_reco_fromB_weighted", 5, bins5D_jes, xmin5D_jes, xmax5D_jes);

//histograms for Jet Energy Resolution (JER)
int	bins5D_jer[5]   =   { 50  ,  100  , 80,   8, multbinsize-1};
double xmin5D_jer[5]   =   { -1.0 ,  0   , -4.0, 0, 0};
double xmax5D_jer[5]   =   { 1.0  ,  1000, 4.0,  8, (double) multbinsize-1};
THnSparseD *hist_jer_reco_weighted = new THnSparseD("hist_jer_reco_weighted", "hist_jer_reco_weighted", 5, bins5D_jer, xmin5D_jer, xmax5D_jer);
THnSparseD *hist_jer_reco_fromB_weighted = new THnSparseD("hist_jer_reco_fromB_weighted", "hist_jer_reco_fromB_weighted", 5, bins5D_jer, xmin5D_jer, xmax5D_jer);

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
	NJets->Sumw2();
	NJetsSub->Sumw2();
	NJetsLead->Sumw2();
	NJetsLJSLJ->Sumw2();
	gen_mult->Sumw2();
	gen_mult_weighted->Sumw2();
	reco_mult->Sumw2();
	reco_mult_weighted->Sumw2();
	multiplicity->Sumw2();
	multiplicity_weighted->Sumw2();
	vzhist->Sumw2();
	vzhist_weighted->Sumw2();
	vzhist_jet_weighted->Sumw2();
	vzhist_dijet_weighted->Sumw2();
	pthathist->Sumw2();
	pthathist_weighted->Sumw2();
	hfhist->Sumw2();
	hfhist_weighted->Sumw2();
	hfhistEta4->Sumw2();
	hfhistEta4_weighted->Sumw2();
	zdchist->Sumw2();
	zdchist_weighted->Sumw2();
	multiplicity_withonejet40->Sumw2();
	multiplicity_withonejet50->Sumw2();
	multiplicity_withonejet60->Sumw2();
	multiplicity_withonejet80->Sumw2();
	multiplicity_withonejet100->Sumw2();
	multiplicity_withonejet->Sumw2();
	multiplicity_withonejet_weighted->Sumw2();
	reco_mult_withonejet->Sumw2();
	reco_mult_withonejet_weighted->Sumw2();
	gen_mult_withonejet->Sumw2();
	gen_mult_withonejet_weighted->Sumw2();
	hfhist_onejet_weighted->Sumw2();
	hfhistEta4_onejet_weighted->Sumw2();
	zdchist_onejet_weighted->Sumw2();
	multiplicity_withdijets->Sumw2();
	multiplicity_withdijets_weighted->Sumw2();
	reco_mult_withdijets->Sumw2();
	reco_mult_withdijets_weighted->Sumw2();
	gen_mult_withdijets->Sumw2();
	gen_mult_withdijets_weighted->Sumw2();
	hfhist_dijet_weighted->Sumw2();
	hfhistEta4_dijet_weighted->Sumw2();
	zdchist_dijet_weighted->Sumw2();
	trackmaxptinjethisto->Sumw2();
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
	hist_matched_jet_weighted_nocut->Sumw2();
	hist_matched_jet->Sumw2();
	hist_matched_jet_weighted->Sumw2();
	hist_matched_jet_parton->Sumw2();
	hist_matched_jet_parton_fromB->Sumw2();
	hist_jes_reco_weighted->Sumw2();
	hist_jer_reco_weighted->Sumw2();
	hist_jes_reco_fromB_weighted->Sumw2();
	hist_jer_reco_fromB_weighted->Sumw2();
	hist_reco_reco_2pcorrelation_signal->Sumw2();
	hist_reco_reco_2pcorrelation_signal_subg0->Sumw2();
	hist_reco_reco_2pcorrelation_signal_subcross->Sumw2();
	hist_reco_reco_2pcorrelation_mixing->Sumw2();
	hist_gen_gen_2pcorrelation_signal->Sumw2();
	hist_gen_gen_2pcorrelation_signal_subg0->Sumw2();
	hist_gen_gen_2pcorrelation_signal_subcross->Sumw2();
	hist_gen_gen_2pcorrelation_mixing->Sumw2();
	EP2_plus_flat->Sumw2();
	EP2_minus_flat->Sumw2();
	EP3_plus_flat->Sumw2();
	EP3_minus_flat->Sumw2();
	EP4_plus_flat->Sumw2();
	EP4_minus_flat->Sumw2();
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
	hist_etaDijet_ETEta4_reco->Sumw2();
	hist_etaDijet_ETEta4_ref->Sumw2();
	hist_etaDijet_ETEta4_gen->Sumw2();
	hist_etaDijet_AJ_ETEta4_reco->Sumw2();
	hist_etaDijet_AJ_ETEta4_ref->Sumw2();
	hist_etaDijet_AJ_ETEta4_gen->Sumw2();
	hist_etaDiff_ETEta4_reco->Sumw2();
	hist_etaDiff_ETEta4_ref->Sumw2();
	hist_etaDiff_ETEta4_gen->Sumw2();
	hist_etaDiff_AJ_ETEta4_reco->Sumw2();
	hist_etaDiff_AJ_ETEta4_ref->Sumw2();
	hist_etaDiff_AJ_ETEta4_gen->Sumw2();
	hist_reco_jetET->Sumw2();
	hist_ref_jetET->Sumw2();
	hist_gen_jetET->Sumw2();
	hist_reco_LjetET->Sumw2();
	hist_gen_LjetET->Sumw2();
	hist_reco_SLjetET->Sumw2();
	hist_gen_SLjetET->Sumw2();
	hist_reco_jet_dijetEpEPb->Sumw2();
	hist_gen_jet_dijetEpEPb->Sumw2();
	hist_reco_jetavept->Sumw2();
	hist_ref_jetavept->Sumw2();
	hist_gen_jetavept->Sumw2();
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
	NJets->Write();
	NJetsSub->Write();
	NJetsLead->Write();
	NJetsLJSLJ->Write();
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
	gen_mult->Write();
	gen_mult_weighted->Write();
	}
	reco_mult->Write();
	reco_mult_weighted->Write();
 	multiplicity->Write();
	multiplicity_weighted->Write();
 	vzhist->Write();
 	vzhist_weighted->Write();
	vzhist_jet_weighted->Write();
	vzhist_dijet_weighted->Write();
	if(isMC){
		pthathist->Write();
		pthathist_weighted->Write();
	}
	trackmaxptinjethisto->Write();
	hfhist->Write();
	hfhist_weighted->Write();
	hfhistEta4->Write();
	hfhistEta4_weighted->Write();
	zdchist->Write();
	zdchist_weighted->Write();
	multiplicity_withonejet40->Write();
	multiplicity_withonejet50->Write();
	multiplicity_withonejet60->Write();
	multiplicity_withonejet80->Write();
	multiplicity_withonejet100->Write();
	multiplicity_withonejet->Write();
	multiplicity_withonejet_weighted->Write();
	reco_mult_withonejet->Write();
	reco_mult_withonejet_weighted->Write();
	if(isMC){
		gen_mult_withonejet->Write();
		gen_mult_withonejet_weighted->Write();
	}
	hfhist_onejet_weighted->Write();
	hfhistEta4_onejet_weighted->Write();
	zdchist_onejet_weighted->Write();
	multiplicity_withdijets->Write();
	multiplicity_withdijets_weighted->Write();
	reco_mult_withdijets->Write();
	reco_mult_withdijets_weighted->Write();
	if(isMC){
		gen_mult_withdijets->Write();
		gen_mult_withdijets_weighted->Write();
	}
	hfhist_dijet_weighted->Write();
	hfhistEta4_dijet_weighted->Write();
	zdchist_dijet_weighted->Write();
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
	hist_reco_jetET->Write();
	if(doleadsubl){
 	hist_reco_leadjet_pt_nocut->Write();
 	hist_reco_leadjet_pt_nocut_weighted->Write();
 	hist_reco_subljet_pt_nocut->Write();
 	hist_reco_subljet_pt_nocut_weighted->Write();
 	hist_reco_leadjet->Write();
 	hist_reco_leadjet_weighted->Write();
 	hist_reco_subljet->Write();
 	hist_reco_subljet_weighted->Write();
		hist_reco_LjetET->Write();
		hist_reco_SLjetET->Write();
	}
	if(isMC){
		//genjetaxischeck->Write();
		hist_matched_jet_weighted_nocut->Write();
		hist_matched_jet->Write();
		hist_matched_jet_weighted->Write();
		hist_matched_jet_parton->Write();
		hist_matched_jet_parton_fromB->Write();
		hist_ref_jetET->Write();
		hist_gen_jet_weighted_nocut->Write();
		hist_gen_jet->Write();
		hist_gen_jet_weighted->Write();
		hist_gen_jetET->Write();
		if(doleadsubl){
			hist_gen_leadjet_pt_nocut->Write();
			hist_gen_leadjet_pt_nocut_weighted->Write();
			hist_gen_subljet_pt_nocut->Write();
			hist_gen_subljet_pt_nocut_weighted->Write();
			hist_gen_leadjet->Write();
			hist_gen_leadjet_weighted->Write();
			hist_gen_subljet->Write();
			hist_gen_subljet_weighted->Write();
			hist_gen_LjetET->Write();
			hist_gen_SLjetET->Write();
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
void w_jetquenching_hist(bool isMC){
	hist_reco_lead_reco_subl_quench_mid_mid->Write();
	hist_reco_lead_reco_subl_quench_mid_fwd->Write();
	hist_reco_lead_reco_subl_quench_mid_bkw->Write();
	hist_reco_lead_reco_subl_quench_fwd_mid->Write();
	hist_reco_lead_reco_subl_quench_bkw_mid->Write();
	hist_reco_lead_reco_subl_quench_fwd_fwd->Write();
	hist_reco_lead_reco_subl_quench_fwd_bkw->Write();
	hist_reco_lead_reco_subl_quench_bkw_fwd->Write();
	hist_reco_lead_reco_subl_quench_bkw_bkw->Write();
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
	}
}

// JES and JER histograms
void w_jes_jer_hist(){ 
	hist_jes_reco_weighted->Write();
	hist_jer_reco_weighted->Write();
	hist_jes_reco_fromB_weighted->Write();
	hist_jer_reco_fromB_weighted->Write();
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

// eta assymetry histograms
void w_etassym_hist(bool isMC){ 
	hist_etaDijet_ETEta4_reco->Write();
	hist_etaDijet_AJ_ETEta4_reco->Write();
	hist_etaDiff_ETEta4_reco->Write();
	hist_etaDiff_AJ_ETEta4_reco->Write();
	hist_reco_jet_dijetEpEPb->Write();
	hist_reco_jetavept->Write();
	if(isMC){
		hist_etaDijet_ETEta4_ref->Write();
		hist_etaDijet_AJ_ETEta4_ref->Write();
		hist_etaDijet_ETEta4_gen->Write();
		hist_etaDijet_AJ_ETEta4_gen->Write();
		hist_etaDiff_ETEta4_ref->Write();
		hist_etaDiff_AJ_ETEta4_ref->Write();
		hist_etaDiff_ETEta4_gen->Write();
		hist_etaDiff_AJ_ETEta4_gen->Write();
		hist_ref_jetavept->Write();
		hist_gen_jet_dijetEpEPb->Write();
		hist_gen_jetavept->Write();
	}
}