#include "call_libraries.h"  // call libraries from ROOT and C++
#include "input_variables.h" // call inputs

// ------------------------------------------------ //
// ----------------- Bin Definition --------------- //
// ------------------------------------------------ //

const int trkbinsize = (int) trk_pt_bins.size(); // track bins for jet-track correlation
const int multbinsize = (int) multiplicity_centrality_bins.size();// multiplicity or centrality bins for jet-track correlation
const int extrabinsize = (int) extra_bins.size();// any additional dependency you wanna add (be carefull about memory)

double minpthist = (double) ((jet_pt_min_cut < subleading_pT_min) ? jet_pt_min_cut : subleading_pT_min);
double maxpthist = (double) jet_pt_max_cut;
double minxjhist = (double) xjmin;
double maxxjhist = (double) xjmax;
double mindphihist = 0.0;
double maxdphihist = TMath::Pi();
int Netadijethist = 60;
double minetadijethist = -6.0;
double maxetadijethist = 6.0;


// binning definition
const double binnerShift = 0.0; // shift if starts at 0, Log(0) -> error
// Needed to define log binning
// Xp and XPb -> see sumw2 function
const int nXBins = 40; // number of bins
const double minX = 3e-04;  // minimum
const double maxX = 1.0; 	   // maximum
double XlogBinWidth = (TMath::Log(maxX+binnerShift) - TMath::Log(minX+binnerShift)) / nXBins; // binwidth

// Xj and Aj bins
const int nXjAjBins = 20; // number of bins
double XjBins[nXjAjBins+1] = {minxjhist,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,maxxjhist};
double AjBins[nXjAjBins+1] = {minxjhist,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,maxxjhist};
//Dphibins
const int nDphiBins = 30; // number of bins
double DphiBins[nDphiBins+1] = {0.0, TMath::Pi()/5. ,TMath::Pi()/3., (3./7.)*TMath::Pi(), TMath::Pi()/2., (4./7.)*TMath::Pi(), (3./5.)*TMath::Pi(), 1.93731547,  1.98967535,  2.04203522,  2.0943951 , 2.14675498,  2.19911486,  2.25147474,  2.30383461,  2.35619449, 2.40855437,  2.46091425,  2.51327412,  2.565634,  2.61799388, 2.67035376,  2.72271363,  2.77507351,  2.82743339,  2.87979327, 2.93215314,  2.98451302,  3.0368729 ,  3.08923278,  TMath::Pi()};

// leading and subleading jet pTs
const int nPtLSLBins = 33; // number of bins
const double minPtLSL = minpthist-10.0;  // minimum
const double maxPtLSL = maxpthist; // maximum
double PtLSLBins[nPtLSLBins+1] = {minPtLSL, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 140.0, 160.0, 180.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 350.0, 400.0, 450.0, 500.0, 600.0, maxPtLSL};

// -------------------------------------------------------------------------- //
// ============================= Unfolding stuff ============================ //
// -------------------------------------------------------------------------- //

// Binning for the two-dimensional unfolding
// --> xj vs pt average
const int nUnfoldingBins_xjptave = nXjAjBins*nPtLSLBins;
const double minUnfoldingBin_xjptave = 0;
const double maxUnfoldingBin_xjptave = nPtLSLBins*maxxjhist; //nJetPtBinsEEC*maxDeltaREEC;
double fullUnfoldingBinning_xjptave[nUnfoldingBins_xjptave+1];
// histograms for unfolding
// Axis : 0 -> reco, 1 -> gen, 2 -> delta phi reco, 3 -> delta phi gen, 4 -> mid-fwd-bkw ranges reco, 5 -> mid-fwd-bkw ranges gen,  6 -> multiplicity
int	bins_unfxjptave[7] 		=   { nUnfoldingBins_xjptave   , nUnfoldingBins_xjptave  , nDphiBins   , nDphiBins   , 11  , 11  , multbinsize-1};
double xmin_unfxjptave[7]   =   { minUnfoldingBin_xjptave  , minUnfoldingBin_xjptave , mindphihist , mindphihist	     , 0.0 , 0.0 , multiplicity_centrality_bins[0]};
double xmax_unfxjptave[7]   =   { maxUnfoldingBin_xjptave  , maxUnfoldingBin_xjptave , maxdphihist , maxdphihist , 11.0, 11.0, multiplicity_centrality_bins[multbinsize-1]};
THnSparseD *fhUnfoldingResponse_xjptave = new THnSparseD("fhUnfoldingResponse_xjptave", "fhUnfoldingResponse_xjptave", 7, bins_unfxjptave, xmin_unfxjptave, xmax_unfxjptave);
// Axis : 0 -> Measured or Truth, 1 -> delta phi measured or truth, 2 ->  mid-fwd-bkw ranges, 3 -> multiplicity
int	bins_unfxjptaveMT[4] 	  =   { nUnfoldingBins_xjptave   , nDphiBins   , 11  , multbinsize-1};
double xmin_unfxjptaveMT[4]   =   { minUnfoldingBin_xjptave  , mindphihist	       , 0.0 , multiplicity_centrality_bins[0]};
double xmax_unfxjptaveMT[4]   =   { maxUnfoldingBin_xjptave  , maxdphihist , 11.0, multiplicity_centrality_bins[multbinsize-1]};
THnSparseD *fhUnfoldingMeasu_xjptave = new THnSparseD("fhUnfoldingMeasu_xjptave", "fhUnfoldingMeasu_xjptave", 4, bins_unfxjptaveMT, xmin_unfxjptaveMT, xmax_unfxjptaveMT);
THnSparseD *fhUnfoldingTruthRef_xjptave = new THnSparseD("fhUnfoldingTruthRef_xjptave", "fhUnfoldingTruthRef_xjptave", 4, bins_unfxjptaveMT, xmin_unfxjptaveMT, xmax_unfxjptaveMT);
THnSparseD *fhUnfoldingTruthGen_xjptave = new THnSparseD("fhUnfoldingTruthGen_xjptave", "fhUnfoldingTruthGen_xjptave", 4, bins_unfxjptaveMT, xmin_unfxjptaveMT, xmax_unfxjptaveMT);

// --> pt1 vs pt2
const int nUnfoldingBins_pt1pt2 = nPtLSLBins*nPtLSLBins;
const double minUnfoldingBin_pt1pt2 = 0;
const double maxUnfoldingBin_pt1pt2 = nPtLSLBins*maxpthist;
double fullUnfoldingBinning_pt1pt2[nUnfoldingBins_pt1pt2+1];
// histograms for unfolding
// Axis : 0 -> reco, 1 -> gen, 2 -> delta phi reco, 3 -> delta phi gen, 4 -> mid-fwd-bkw ranges reco, 5 -> mid-fwd-bkw ranges gen,  6 -> multiplicity
int	bins_unfpt1pt2[7] 	   =   { nUnfoldingBins_pt1pt2   , nUnfoldingBins_pt1pt2  , nDphiBins   , nDphiBins   , 11  , 11   , multbinsize-1};
double xmin_unfpt1pt2[7]   =   { minUnfoldingBin_pt1pt2  , minUnfoldingBin_pt1pt2 , mindphihist , mindphihist , 0.0 , 0.0  , multiplicity_centrality_bins[0]};
double xmax_unfpt1pt2[7]   =   { maxUnfoldingBin_pt1pt2  , maxUnfoldingBin_pt1pt2 , maxdphihist , maxdphihist , 11.0, 11.0 , multiplicity_centrality_bins[multbinsize-1]};
THnSparseD *fhUnfoldingResponse_pt1pt2 = new THnSparseD("fhUnfoldingResponse_pt1pt2", "fhUnfoldingResponse_pt1pt2", 7, bins_unfpt1pt2, xmin_unfpt1pt2, xmax_unfpt1pt2);
// Axis : 0 -> Measured or Truth, 1 -> delta phi measured or truth, 2 ->  mid-fwd-bkw ranges, 3 -> multiplicity
int	bins_unfpt1pt2MT[4] 	 =   { nUnfoldingBins_pt1pt2   , nDphiBins   , 11  , multbinsize-1};
double xmin_unfpt1pt2MT[4]   =   { minUnfoldingBin_pt1pt2  , mindphihist , 0.0 , multiplicity_centrality_bins[0]};
double xmax_unfpt1pt2MT[4]   =   { maxUnfoldingBin_pt1pt2  , maxdphihist , 11.0, multiplicity_centrality_bins[multbinsize-1]};
THnSparseD *fhUnfoldingMeasu_pt1pt2 = new THnSparseD("fhUnfoldingMeasu_pt1pt2", "fhUnfoldingMeasu_pt1pt2", 4, bins_unfpt1pt2MT, xmin_unfpt1pt2MT, xmax_unfpt1pt2MT);
THnSparseD *fhUnfoldingTruthRef_pt1pt2 = new THnSparseD("fhUnfoldingTruthRef_pt1pt2", "fhUnfoldingTruthRef_pt1pt2", 4, bins_unfpt1pt2MT, xmin_unfpt1pt2MT, xmax_unfpt1pt2MT);
THnSparseD *fhUnfoldingTruthGen_pt1pt2 = new THnSparseD("fhUnfoldingTruthGen_pt1pt2", "fhUnfoldingTruthGen_pt1pt2", 4, bins_unfpt1pt2MT, xmin_unfpt1pt2MT, xmax_unfpt1pt2MT);

// xj vs pt2 in bins of pt1
const int nUnfoldingBins_xjpt2 = nXjAjBins*nPtLSLBins;
const double minUnfoldingBin_xjpt2 = 0;
const double maxUnfoldingBin_xjpt2 = nPtLSLBins*maxxjhist; //nJetPtBinsEEC*maxDeltaREEC;
double fullUnfoldingBinning_xjpt2[nUnfoldingBins_xjpt2+1];
// histograms for unfolding
// Axis : 0 -> reco, 1 -> gen, 2 -> delta phi reco, 3 -> delta phi gen, 4 -> mid-fwd-bkw ranges reco, 5 -> mid-fwd-bkw ranges gen, 6 -> multiplicity, 7 -> Leading jet pT reco, 8 -> Leading jet pT gen
int	bins_unfxjpt2[9] 	  =   { nUnfoldingBins_xjpt2   , nUnfoldingBins_xjpt2  , nDphiBins   , nDphiBins   , 11  , 11  , multbinsize-1								, nPtLSLBins , nPtLSLBins};
double xmin_unfxjpt2[9]   =   { minUnfoldingBin_xjpt2  , minUnfoldingBin_xjpt2 , mindphihist , mindphihist , 0.0 , 0.0 , multiplicity_centrality_bins[0]			, minPtLSL	 , minPtLSL};
double xmax_unfxjpt2[9]   =   { maxUnfoldingBin_xjpt2  , maxUnfoldingBin_xjpt2 , maxdphihist , maxdphihist , 11.0, 11.0, multiplicity_centrality_bins[multbinsize-1], minPtLSL	 , minPtLSL};
THnSparseD *fhUnfoldingResponse_xjpt2 = new THnSparseD("fhUnfoldingResponse_xjpt2", "fhUnfoldingResponse_xjpt2", 9, bins_unfxjpt2, xmin_unfxjpt2, xmax_unfxjpt2);
// Axis : 0 -> Measured or Truth, 1 -> delta phi measured or truth, 2 ->  mid-fwd-bkw ranges, 3 -> multiplicity, 4 -> leading jet pT measured or truth
int	bins_unfxjpt2MT[5] 	    =   { nUnfoldingBins_xjpt2   , nDphiBins   , 11  , multbinsize-1								, nPtLSLBins};
double xmin_unfxjpt2MT[5]   =   { minUnfoldingBin_xjpt2  , mindphihist , 0.0 , multiplicity_centrality_bins[0]  			, minPtLSL};
double xmax_unfxjpt2MT[5]   =   { maxUnfoldingBin_xjpt2  , maxdphihist , 11.0, multiplicity_centrality_bins[multbinsize-1]  , maxPtLSL};
THnSparseD *fhUnfoldingMeasu_xjpt2 = new THnSparseD("fhUnfoldingMeasu_xjpt2", "fhUnfoldingMeasu_xjpt2", 5, bins_unfxjpt2MT, xmin_unfxjpt2MT, xmax_unfxjpt2MT);
THnSparseD *fhUnfoldingTruthRef_xjpt2 = new THnSparseD("fhUnfoldingTruthRef_xjpt2", "fhUnfoldingTruthRef_xjpt2", 5, bins_unfxjpt2MT, xmin_unfxjpt2MT, xmax_unfxjpt2MT);
THnSparseD *fhUnfoldingTruthGen_xjpt2 = new THnSparseD("fhUnfoldingTruthGen_xjpt2", "fhUnfoldingTruthGen_xjpt2", 5, bins_unfxjpt2MT, xmin_unfxjpt2MT, xmax_unfxjpt2MT);

// -------------------------------------------------------------------------- //
// ============================ Event quantities ============================ //
// -------------------------------------------------------------------------- //

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

TH1D *Nev_alljetfromalltrk = new TH1D("Nev_alljetfromalltrk", "Nev_alljetfromalltrk", 80, 0, 400);
TH1D *Nev_jetwithlowpttrk = new TH1D("Nev_jetwithlowpttrk", "Nev_jetwithlowpttrk", 80, 0, 400);
TH1D *Nev_jetfromonetrk = new TH1D("Nev_jetfromonetrk", "Nev_jetfromonetrk", 80, 0, 400);
TH1D *Nev_jetsfrombothlowpttrkandonetrk = new TH1D("Nev_jetsfrombothlowpttrkandonetrk", "Nev_jetsfrombothlowpttrkandonetrk", 80, 0, 400);
TH1D *Nev_jetwithlowpttrk_lead = new TH1D("Nev_jetwithlowpttrk_lead", "Nev_jetwithlowpttrk_lead", 80, 0, 400);
TH1D *Nev_jetwithlowpttrk_sublead = new TH1D("Nev_jetwithlowpttrk_sublead", "Nev_jetwithlowpttrk_sublead", 80, 0, 400);
TH1D *Nev_jetfromonetrk_lead = new TH1D("Nev_jetfromonetrk_lead", "Nev_jetfromonetrk_lead", 80, 0, 400);
TH1D *Nev_jetfromonetrk_sublead = new TH1D("Nev_jetfromonetrk_sublead", "Nev_jetfromonetrk_sublead", 80, 0, 400);

// Multiplicities
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 80, 0.0, 400.0);
TH1D *multiplicity_weighted = new TH1D("multiplicity_weighted", "multiplicity_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_withonejet_weighted = new TH1D("multiplicity_withonejet_weighted", "multiplicity_withonejet_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_withdijets_weighted = new TH1D("multiplicity_withdijets_weighted", "multiplicity_withdijets_weighted", 80, 0.0, 400.0);
TH1D *multiplicity_nocut = new TH1D("multiplicity_nocut", "multiplicity_nocut", 500, 0.0, 500.0);
TH1D *multiplicity_corrected = new TH1D("multiplicity_corrected", "multiplicity_corrected", 500, 0.0, 500.0);
TH1D *multiplicity_tight = new TH1D("multiplicity_tight", "multiplicity_tight", 500, 0.0, 500.0);
TH1D *multiplicity_loose = new TH1D("multiplicity_loose", "multiplicity_loose", 500, 0.0, 500.0);
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
TH2D *multiplicity_etadep_weighted = new TH2D("multiplicity_etadep_weighted", "multiplicity_etadep_weighted", 80, 0.0, 400.0, 11, 0.0, 11.0);
TH2D *multiplicity2D = new TH2D("multiplicity2D", "multiplicity2D", 500, 0.0, 500.0, 500, 0.0, 500.0);

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
int	bins_ZDC[3]      =   {  100   ,   200   ,  80};
double xmin_ZDC[3]   =   { -10000 ,  -10000 ,  0.0};
double xmax_ZDC[3]   =   {  10000 ,   10000 ,  400};
THnSparseD *zdchist = new THnSparseD("zdchist", "zdchist", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);
THnSparseD *zdchist_weighted = new THnSparseD("zdchist_weighted", "zdchist_weighted", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);
THnSparseD *zdchist_onejet_weighted = new THnSparseD("zdchist_onejet_weighted", "zdchist_onejet_weighted", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);
THnSparseD *zdchist_dijet_weighted = new THnSparseD("zdchist_dijet_weighted", "zdchist_dijet_weighted", 3, bins_ZDC, xmin_ZDC, xmax_ZDC);

// HFSum
int	bins_HFSum[2]   =      { 250  ,  80};
double xmin_HFSum[2]   =   { 0.0  ,  0.0};
double xmax_HFSum[2]   =   { 500  ,  400};
THnSparseD *hfhistSum_weighted = new THnSparseD("hfhistSum_weighted", "hfhistSum_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistEta4Sum_weighted = new THnSparseD("hfhistEta4Sum_weighted", "hfhistEta4Sum_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistSum_onejet_weighted = new THnSparseD("hfhistSum_onejet_weighted", "hfhistSum_onejet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistEta4Sum_onejet_weighted = new THnSparseD("hfhistEta4Sum_onejet_weighted", "hfhistEta4Sum_onejet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistSum_dijet_weighted = new THnSparseD("hfhistSum_dijet_weighted", "hfhistSum_dijet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
THnSparseD *hfhistEta4Sum_dijet_weighted = new THnSparseD("hfhistEta4Sum_dijet_weighted", "hfhistEta4Sum_dijet_weighted", 2, bins_HFSum, xmin_HFSum, xmax_HFSum);
TH1D *hfhistEta4Sum_nocut = new TH1D("hfhistEta4Sum_nocut", "hfhistEta4Sum_nocut",500,0.0,500.0);
TH1D *hfhistEta4Plus_nocut = new TH1D("hfhistEta4Plus_nocut", "hfhistEta4Plus_nocut",500,0.0,500.0);
TH1D *hfhistEta4Minus_nocut = new TH1D("hfhistEta4Minus_nocut", "hfhistEta4Minus_nocut",500,0.0,500.0);


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

// -------------------------------------------------------------------------- //
// ============================ Track quantities ============================ //
// -------------------------------------------------------------------------- //
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

// Jet histograms
// Number of jets per event
int	bins_NJETS[3]      =   {  30   ,   multbinsize-1    		 					 ,  extrabinsize-1};
double xmin_NJETS[3]   =   {  0.0  ,   multiplicity_centrality_bins[0]  			 ,  extra_bins[0]};
double xmax_NJETS[3]   =   {  30.0 ,   multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *NJets = new THnSparseD("NJets", "NJets", 3, bins_NJETS, xmin_NJETS, xmax_NJETS);

// trackmax histogram
// Axis : 0 -> max track pt in a jet, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_trkmax[3]      =   {  100   ,   multbinsize-1    							    ,  extrabinsize-1};
double xmin_trkmax[3]   =   {  0.0 	  ,   multiplicity_centrality_bins[0]  			    ,  extra_bins[0]};
double xmax_trkmax[3]   =   {  100.0 ,   multiplicity_centrality_bins[multbinsize-1]   ,  extra_bins[extrabinsize-1]};
THnSparseD *trackmaxptinjethisto = new THnSparseD("trackmaxptinjethisto", "trackmaxptinjethisto", 3, bins_trkmax, xmin_trkmax, xmax_trkmax);

// trackmax/rawjet histogram
// Axis : 0 -> max track pt in a jet over raw pT, 1 -> raw pT, 2 -> multiplicity bins
int	bins_trkmaxjet[3]      =   {  100   , 100  , multbinsize-1};
double xmin_trkmaxjet[3]   =   {  0.0   , 0    , multiplicity_centrality_bins[0]};
double xmax_trkmaxjet[3]   =   {  1.0   , 1000 , multiplicity_centrality_bins[multbinsize-1]};
THnSparseD *jettrackmaxptinjethisto = new THnSparseD("jettrackmaxptinjethisto", "jettrackmaxptinjethisto", 3, bins_trkmaxjet, xmin_trkmaxjet, xmax_trkmaxjet);
THnSparseD *jettrackmaxptinjethisto_no0 = new THnSparseD("jettrackmaxptinjethisto_no0", "jettrackmaxptinjethisto_no0", 3, bins_trkmaxjet, xmin_trkmaxjet, xmax_trkmaxjet);
THnSparseD *jettrackmaxptinjethisto_ref = new THnSparseD("jettrackmaxptinjethisto_ref", "jettrackmaxptinjethisto_ref", 3, bins_trkmaxjet, xmin_trkmaxjet, xmax_trkmaxjet);

// UE histogram
// Axis : 0 -> UE, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_UE[3]      =   {  400  	,   multbinsize-1    							    ,  extrabinsize-1};
double xmin_UE[3]   =   {  -1000.0   	,   multiplicity_centrality_bins[0]  			    ,  extra_bins[0]};
double xmax_UE[3]   =   {  1000.0   ,   multiplicity_centrality_bins[multbinsize-1]   	,  extra_bins[extrabinsize-1]};
THnSparseD *histo_jetUE 			 = new THnSparseD("histo_jetUE", "histo_jetUE", 3, bins_UE, xmin_UE, xmax_UE);
THnSparseD *histo_jetAverageRho 	 = new THnSparseD("histo_jetAverageRho", "histo_jetAverageRho", 3, bins_UE, xmin_UE, xmax_UE);
THnSparseD *histo_jetcheckcorrection = new THnSparseD("histo_jetcheckcorrection", "histo_jetcheckcorrection", 3, bins_UE, xmin_UE, xmax_UE);

// Jet correlations to EP
// Axis : 0 -> delta phi between jet and EP, 1 -> multiplicity bins, 2 -> extra dependence bins
int	bins_JETEP[3]      =   { 40			   		,  multbinsize-1							 	 ,  extrabinsize-1};
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
int	bins_jet[5]      =   { 80	  ,  60   ,   32		      , multbinsize-1          						, extrabinsize-1};
double xmin_jet[5]   =   { 0.0	  , -6.0  ,   -TMath::Pi()    , multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jet[5]   =   { 800.0  ,  6.0 ,   TMath::Pi()     , multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
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
THnSparseD *hist_reco_thrdproj_weighted = new THnSparseD("hist_reco_thrdproj_weighted", "hist_reco_thrdproj_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_reco_thrdprojdiff_weighted = new THnSparseD("hist_reco_thrdprojdiff_weighted", "hist_reco_thrdprojdiff_weighted", 5, bins_jet, xmin_jet, xmax_jet);
// --> Gen
THnSparseD *hist_gen_leadjet = new THnSparseD("hist_gen_leadjet", "hist_gen_leadjet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_leadjet_weighted = new THnSparseD("hist_gen_leadjet_weighted", "hist_gen_leadjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_subljet = new THnSparseD("hist_gen_subljet", "hist_gen_subljet", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_subljet_weighted = new THnSparseD("hist_gen_subljet_weighted", "hist_gen_subljet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_thrdjet_weighted = new THnSparseD("hist_gen_thrdjet_weighted", "hist_gen_thrdjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_thrdproj_weighted = new THnSparseD("hist_gen_thrdproj_weighted", "hist_gen_thrdproj_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_gen_thrdprojdiff_weighted = new THnSparseD("hist_gen_thrdprojdiff_weighted", "hist_gen_thrdprojdiff_weighted", 5, bins_jet, xmin_jet, xmax_jet);
// --> Ref
THnSparseD *hist_ref_jet_weighted = new THnSparseD("hist_ref_jet_weighted", "hist_ref_jet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_leadjet_weighted = new THnSparseD("hist_ref_leadjet_weighted", "hist_ref_leadjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_subljet_weighted = new THnSparseD("hist_ref_subljet_weighted", "hist_ref_subljet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_thrdjet_weighted = new THnSparseD("hist_ref_thrdjet_weighted", "hist_ref_thrdjet_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_thrdproj_weighted = new THnSparseD("hist_ref_thrdproj_weighted", "hist_ref_thrdproj_weighted", 5, bins_jet, xmin_jet, xmax_jet);
THnSparseD *hist_ref_thrdprojdiff_weighted = new THnSparseD("hist_ref_thrdprojdiff_weighted", "hist_ref_thrdprojdiff_weighted", 5, bins_jet, xmin_jet, xmax_jet);

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
int	bins_jes[6]   =      { 200  ,  80    ,  60  ,  8, multbinsize-1          					 , extrabinsize-1};
double xmin_jes[6]   =   { 0.0  ,  0     , -6.0 ,  0, multiplicity_centrality_bins[0]             , extra_bins[0]};
double xmax_jes[6]   =   { 5.0  ,  800.0 ,  6.0 ,  8, multiplicity_centrality_bins[multbinsize-1] , extra_bins[extrabinsize-1]};
THnSparseD *hist_jes_reco_weighted = new THnSparseD("hist_jes_reco_weighted", "hist_jes_reco_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_jes_reco_fromB_weighted = new THnSparseD("hist_jes_reco_fromB_weighted", "hist_jes_reco_fromB_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_leadjes_reco_weighted = new THnSparseD("hist_leadjes_reco_weighted", "hist_leadjes_reco_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_leadjes_reco_fromB_weighted = new THnSparseD("hist_leadjes_reco_fromB_weighted", "hist_leadjes_reco_fromB_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_subleadjes_reco_weighted = new THnSparseD("hist_subleadjes_reco_weighted", "hist_subleadjes_reco_weighted", 6, bins_jes, xmin_jes, xmax_jes);
THnSparseD *hist_subleadjes_reco_fromB_weighted = new THnSparseD("hist_subleadjes_reco_fromB_weighted", "hist_subleadjes_reco_fromB_weighted", 6, bins_jes, xmin_jes, xmax_jes);

// Quenching studies
// Axis : 0 -> Xj, 1 -> Aj, 2 -> delta phi, 3 -> multiplicity , 4 -> jet pT average, 5 -> extra dependency, 6 -> pT leading jet, 7 -> pT subleading jet, 8 -> mid-fwd-bkw combinations
int	bins_quenc[9]   =      { nXjAjBins   , nXjAjBins	  , nDphiBins		, multbinsize-1		  	 					  ,	 nPtLSLBins  ,  extrabinsize-1 				, nPtLSLBins, nPtLSLBins, 11};
double xmin_quenc[9]   =   { minxjhist 	 , minxjhist	  , mindphihist    	, multiplicity_centrality_bins[0]		   	  ,	 minPtLSL	 ,	extra_bins[0]				, minPtLSL	, minPtLSL  , 0.0};
double xmax_quenc[9]   =   { maxxjhist   , maxxjhist	  , maxdphihist 	, multiplicity_centrality_bins[multbinsize-1] ,  maxPtLSL    ,  extra_bins[extrabinsize-1]  , maxPtLSL	, maxPtLSL  , 11.0};
THnSparseD *hist_reco_lead_reco_subl_quench = new THnSparseD("hist_reco_lead_reco_subl_quench", "hist_reco_lead_reco_subl_quench", 9, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_fake_lead_reco_subl_quench = new THnSparseD("hist_fake_lead_reco_subl_quench", "hist_fake_lead_reco_subl_quench", 9, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_reco_lead_fake_subl_quench = new THnSparseD("hist_reco_lead_fake_subl_quench", "hist_reco_lead_fake_subl_quench", 9, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_fake_lead_fake_subl_quench = new THnSparseD("hist_fake_lead_fake_subl_quench", "hist_fake_lead_fake_subl_quench", 9, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_ref_lead_ref_subl_quench = new THnSparseD("hist_ref_lead_ref_subl_quench", "hist_ref_lead_ref_subl_quench", 9, bins_quenc, xmin_quenc, xmax_quenc);
THnSparseD *hist_gen_lead_gen_subl_quench = new THnSparseD("hist_gen_lead_gen_subl_quench", "hist_gen_lead_gen_subl_quench", 9, bins_quenc, xmin_quenc, xmax_quenc);

// EP dependency
// Axis : 0 -> Xj, 1 -> delta phi, 2 -> 2*|Psi2 - jetphi|, 3 -> 3*|Psi3 - jetphi|, 4 -> 4*|Psi4 - jetphi|, 5 -> multiplicity , 6 -> jet pT average, 7 -> extra dependency , 8 ->  mid-fwd-bkw combinations
int	bins_quencEP[9]   =      { nXjAjBins  , nDphiBins		, 16 		 , 16 			, 16 			, multbinsize-1		  	 					  ,	 nPtLSLBins,  extrabinsize-1 			, 11};
double xmin_quencEP[9]   =   { minxjhist  , mindphihist	    , 0.0		 , 0.0			, 0.0			, multiplicity_centrality_bins[0]		   	  ,	 minPtLSL  ,  extra_bins[0]   			, 0.0};
double xmax_quencEP[9]   =   { maxxjhist  , maxdphihist 	, TMath::Pi(), TMath::Pi()	, TMath::Pi()	, multiplicity_centrality_bins[multbinsize-1] ,  maxPtLSL  ,  extra_bins[extrabinsize-1], 11.0};
THnSparseD *hist_reco_leadEP_quench_plus = new THnSparseD("hist_reco_leadEP_quench_plus", "hist_reco_leadEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_leadEP_quench_minus = new THnSparseD("hist_reco_leadEP_quench_minus", "hist_reco_leadEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_plus = new THnSparseD("hist_reco_sublEP_quench_plus", "hist_reco_sublEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_reco_sublEP_quench_minus = new THnSparseD("hist_reco_sublEP_quench_minus", "hist_reco_sublEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_plus = new THnSparseD("hist_ref_sublEP_quench_plus", "hist_ref_sublEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_sublEP_quench_minus = new THnSparseD("hist_ref_sublEP_quench_minus", "hist_ref_sublEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_plus = new THnSparseD("hist_ref_leadEP_quench_plus", "hist_ref_leadEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_ref_leadEP_quench_minus = new THnSparseD("hist_ref_leadEP_quench_minus", "hist_ref_leadEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_plus = new THnSparseD("hist_gen_leadEP_quench_plus", "hist_gen_leadEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_leadEP_quench_minus = new THnSparseD("hist_gen_leadEP_quench_minus", "hist_gen_leadEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_plus = new THnSparseD("hist_gen_sublEP_quench_plus", "hist_gen_sublEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_gen_sublEP_quench_minus = new THnSparseD("hist_gen_sublEP_quench_minus", "hist_gen_sublEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_fake_leadEP_quench_plus = new THnSparseD("hist_fake_leadEP_quench_plus", "hist_fake_leadEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_fake_leadEP_quench_minus = new THnSparseD("hist_fake_leadEP_quench_minus", "hist_fake_leadEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_fake_sublEP_quench_plus = new THnSparseD("hist_fake_sublEP_quench_plus", "hist_fake_sublEP_quench_plus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);
THnSparseD *hist_fake_sublEP_quench_minus = new THnSparseD("hist_fake_sublEP_quench_minus", "hist_fake_sublEP_quench_minus", 9, bins_quencEP, xmin_quencEP, xmax_quencEP);

// eta dijet studies
// Axis : 0 -> etaDijet, 1 -> delta eta / 2, 2 -> Xj, 3 -> Aj, 4 -> delta phi, 5 -> x_p, 6 -> x_Pb, 7 -> multiplicity, 8 -> jet pT average, 9 -> extra dependency
int	bins_etaDijet[10]      =   {  Netadijethist    ,  24  , nXjAjBins	  , nXjAjBins , nDphiBins	    , nXBins ,  nXBins  ,	 multbinsize-1		  	 					  ,	 nPtLSLBins	 ,  extrabinsize-1};
double xmin_etaDijet[10]   =   {  minetadijethist  , -6.0 , minxjhist	  , minxjhist , mindphihist     , minX   ,  minX    , 	 multiplicity_centrality_bins[0]		   	  ,	 minPtLSL	 ,  extra_bins[0]};
double xmax_etaDijet[10]   =   {  maxetadijethist  ,  6.0 , maxxjhist	  , maxxjhist , maxdphihist 	, maxX   ,  maxX    , 	 multiplicity_centrality_bins[multbinsize-1]  ,  maxPtLSL    ,  extra_bins[extrabinsize-1]};
THnSparseD *hist_etaDijet_reco = new THnSparseD("hist_etaDijet_reco", "hist_etaDijet_reco", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_reco = new THnSparseD("hist_etaDijet_CM_reco", "hist_etaDijet_CM_reco", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_ref = new THnSparseD("hist_etaDijet_ref", "hist_etaDijet_ref", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_ref = new THnSparseD("hist_etaDijet_CM_ref", "hist_etaDijet_CM_ref", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_gen = new THnSparseD("hist_etaDijet_gen", "hist_etaDijet_gen", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_etaDijet_CM_gen = new THnSparseD("hist_etaDijet_CM_gen", "hist_etaDijet_CM_gen", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_yDijet_CM_reco = new THnSparseD("hist_yDijet_CM_reco", "hist_yDijet_CM_reco", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_yDijet_CM_ref = new THnSparseD("hist_yDijet_CM_ref", "hist_yDijet_CM_ref", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);
THnSparseD *hist_yDijet_CM_gen = new THnSparseD("hist_yDijet_CM_gen", "hist_yDijet_CM_gen", 10, bins_etaDijet, xmin_etaDijet, xmax_etaDijet);

// Axis : 0 -> in-jet multiplicity, 1 -> multiplicity, 2 -> extra dimension
int	bins_injettrk[3]   	  =   { 50 , multbinsize-1			  						  ,  extrabinsize-1};
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
int	bins_jettrk[5]      =   { 200				    , 500  ,   trkbinsize-1			  	 	 	 , multbinsize-1			  						  ,  extrabinsize-1};
double xmin_jettrk[5]   =   { -TMath::Pi()/2.0		, -5.0 ,   trk_pt_bins[0]					 , multiplicity_centrality_bins[0]					  ,	 extra_bins[0]};
double xmax_jettrk[5]   =   { 3.0*TMath::Pi()/2.0 	, 5.0  ,   trk_pt_bins[trkbinsize-1] 		 , multiplicity_centrality_bins[multbinsize-1]  	  ,  extra_bins[extrabinsize-1]};

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

// Jet Delta's
// Axis : 0 -> delta phi, 1 -> delta eta, 2 -> multiplicity
int	bins_jetjet[3]      =   { 100 	     , 250  , 80};
double xmin_jetjet[3]   =   { 0.0		 , -5.0 , 0.0};
double xmax_jetjet[3]   =   { TMath::Pi(), 5.0  , 400.0};
THnSparseD *jetjet_Dphi_Deta_reco = new THnSparseD("jetjet_Dphi_Deta_reco", "jetjet_Dphi_Deta_reco", 3, bins_jetjet, xmin_jetjet, xmax_jetjet);
THnSparseD *jetjet_Dphi_Deta_ref = new THnSparseD("jetjet_Dphi_Deta_ref", "jetjet_Dphi_Deta_ref", 3, bins_jetjet, xmin_jetjet, xmax_jetjet);
THnSparseD *jetjet_Dphi_Deta_gen = new THnSparseD("jetjet_Dphi_Deta_gen", "jetjet_Dphi_Deta_gen", 3, bins_jetjet, xmin_jetjet, xmax_jetjet);
THnSparseD *jetjet_Dphi_Deta_reco_diff = new THnSparseD("jetjet_Dphi_Deta_reco_diff", "jetjet_Dphi_Deta_reco_diff", 3, bins_jetjet, xmin_jetjet, xmax_jetjet);
THnSparseD *jetjet_Dphi_Deta_ref_diff = new THnSparseD("jetjet_Dphi_Deta_ref_diff", "jetjet_Dphi_Deta_ref_diff", 3, bins_jetjet, xmin_jetjet, xmax_jetjet);
THnSparseD *jetjet_Dphi_Deta_gen_diff = new THnSparseD("jetjet_Dphi_Deta_gen_diff", "jetjet_Dphi_Deta_gen_diff", 3, bins_jetjet, xmin_jetjet, xmax_jetjet);


// Evaluate uncertainties correctly at ROOT
void sw2(){

// ----------------------------------------------------------
//loops to setup binning
// --> Xp and XpPb 
double XBins[nXBins+1];
for(int a = 0; a <= nXBins; a++){XBins[a] = (minX+binnerShift)*TMath::Exp(a*XlogBinWidth)-binnerShift;}

// --> Trk pT
double TrkPtbins[trkbinsize-1];
for(int a = 0; a<trk_pt_bins.size();a++){TrkPtbins[a] = trk_pt_bins[a];}

// --> Multiplicity binning
double MultCentbins[multbinsize-1];
for(int a = 0; a<multiplicity_centrality_bins.size();a++){MultCentbins[a] = multiplicity_centrality_bins[a];}

// --> Extra binning
double Extrabins[extrabinsize-1];
for(int a = 0; a<extra_bins.size();a++){Extrabins[a] = extra_bins[a];}

// --> Unfolding
for(int ixj = 0; ixj < nXjAjBins; ixj++){
  for(int iJetPtAve = 0; iJetPtAve < nPtLSLBins; iJetPtAve++){
    fullUnfoldingBinning_xjptave[ixj+iJetPtAve*nXjAjBins] = XjBins[ixj]+maxxjhist*iJetPtAve;
  }
}
fullUnfoldingBinning_xjptave[nUnfoldingBins_xjptave] = maxUnfoldingBin_xjptave;

for(int ipt1 = 0; ipt1 < nPtLSLBins; ipt1++){
  for(int ipt2 = 0; ipt2 < nPtLSLBins; ipt2++){
    fullUnfoldingBinning_pt1pt2[ipt1+ipt2*nPtLSLBins] = PtLSLBins[ipt1]+maxpthist*ipt2;
  }
}
fullUnfoldingBinning_pt1pt2[nUnfoldingBins_pt1pt2] = maxUnfoldingBin_pt1pt2;

for(int ixj = 0; ixj < nXjAjBins; ixj++){
  for(int ipt2 = 0; ipt2 < nPtLSLBins; ipt2++){
    fullUnfoldingBinning_xjpt2[ixj+ipt2*nXjAjBins] = XjBins[ixj]+maxxjhist*ipt2;
  }
}
fullUnfoldingBinning_xjpt2[nUnfoldingBins_xjpt2] = maxUnfoldingBin_xjpt2;

// --> Event information histos
// -> sumw2 of the histograms that not needed to be adjusted
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
Nev_alljetfromalltrk->Sumw2();
Nev_jetwithlowpttrk->Sumw2();
Nev_jetfromonetrk->Sumw2();
Nev_jetsfrombothlowpttrkandonetrk->Sumw2();
Nev_jetwithlowpttrk_lead->Sumw2();
Nev_jetwithlowpttrk_sublead->Sumw2();
Nev_jetfromonetrk_lead->Sumw2();
Nev_jetfromonetrk_sublead->Sumw2();
multiplicity->Sumw2();
multiplicity_weighted->Sumw2();
multiplicity_withonejet_weighted->Sumw2();
multiplicity_withdijets_weighted->Sumw2();
multiplicity_nocut->Sumw2();
multiplicity_corrected->Sumw2();
multiplicity_tight->Sumw2();
multiplicity_loose->Sumw2();
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
multiplicity_etadep_weighted->Sumw2();
multiplicity2D->Sumw2();
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
hfhistEta4Sum_nocut->Sumw2();
hfhistEta4Plus_nocut->Sumw2();
hfhistEta4Minus_nocut->Sumw2();
jetjet_Dphi_Deta_reco->Sumw2();
jetjet_Dphi_Deta_ref->Sumw2();
jetjet_Dphi_Deta_gen->Sumw2();
jetjet_Dphi_Deta_reco_diff->Sumw2();
jetjet_Dphi_Deta_ref_diff->Sumw2();
jetjet_Dphi_Deta_gen_diff->Sumw2();

// Multiplicity binning
vzhist->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vzhist_weighted->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vzhist_jet_weighted->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vzhist_dijet_weighted->GetAxis(1)->Set(bins_vz[1],MultCentbins);
vxyhist->GetAxis(2)->Set(bins_vxy[2],MultCentbins);
vxyhist_weighted->GetAxis(2)->Set(bins_vxy[2],MultCentbins);
pthathist->GetAxis(1)->Set(bins_pthat[1],MultCentbins);
pthathist_weighted->GetAxis(1)->Set(bins_pthat[1],MultCentbins);
// Extra binning
vzhist->GetAxis(2)->Set(bins_vz[2],Extrabins);
vzhist_weighted->GetAxis(2)->Set(bins_vz[2],Extrabins);
vzhist_jet_weighted->GetAxis(2)->Set(bins_vz[2],Extrabins);
vzhist_dijet_weighted->GetAxis(2)->Set(bins_vz[2],Extrabins);
vxyhist->GetAxis(3)->Set(bins_vxy[3],Extrabins);
vxyhist_weighted->GetAxis(3)->Set(bins_vxy[3],Extrabins);
pthathist->GetAxis(2)->Set(bins_pthat[2],Extrabins);
pthathist_weighted->GetAxis(2)->Set(bins_pthat[2],Extrabins);
// Sumw2
vzhist->Sumw2();
vzhist_weighted->Sumw2();
vzhist_jet_weighted->Sumw2();
vzhist_dijet_weighted->Sumw2();
vxyhist->Sumw2();
vxyhist_weighted->Sumw2();
pthathist->Sumw2();
pthathist_weighted->Sumw2();

// --> Event plane information histos
// Multiplicity binning
EP2_plus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP2_minus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP3_plus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP3_minus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP4_plus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
EP4_minus_flat->GetAxis(3)->Set(bins_EP[3],MultCentbins);
// Extra binning
EP2_plus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP2_minus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP3_plus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP3_minus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP4_plus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
EP4_minus_flat->GetAxis(4)->Set(bins_EP[4],Extrabins);
// Sumw2
EP2_plus_flat->Sumw2();
EP2_minus_flat->Sumw2();
EP3_plus_flat->Sumw2();
EP3_minus_flat->Sumw2();
EP4_plus_flat->Sumw2();
EP4_minus_flat->Sumw2();

// --> Track QA information histos
// Multiplicity binning
hist_reco_trk->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_reco_trk_corr->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_reco_trk_weighted->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_gen_trk->GetAxis(3)->Set(bins_trk[3],MultCentbins);
hist_gen_trk_weighted->GetAxis(3)->Set(bins_trk[3],MultCentbins);
// Extra binning
hist_reco_trk->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_reco_trk_corr->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_reco_trk_weighted->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_gen_trk->GetAxis(4)->Set(bins_trk[4],Extrabins);
hist_gen_trk_weighted->GetAxis(4)->Set(bins_trk[4],Extrabins);
// Sumw2
hist_reco_trk->Sumw2();
hist_reco_trk_corr->Sumw2();
hist_reco_trk_weighted->Sumw2();
hist_gen_trk->Sumw2();
hist_gen_trk_weighted->Sumw2();

// --> Track information used in the correlations histos
// Multiplicity binning
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
// Extra binning
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
// Sumw2
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

// --> Track correlations to EP histos
// Track Pt binning
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
// Multiplicity binning
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
// Extra binning
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
// Sumw2
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

// --> TrackMax/RawPt and Njets histos
// Multiplicity binning
NJets->GetAxis(1)->Set(bins_NJETS[1],MultCentbins);
trackmaxptinjethisto->GetAxis(1)->Set(bins_trkmax[1],MultCentbins);
jettrackmaxptinjethisto->GetAxis(2)->Set(bins_trkmaxjet[2],MultCentbins);
jettrackmaxptinjethisto_no0->GetAxis(2)->Set(bins_trkmaxjet[2],MultCentbins);
jettrackmaxptinjethisto_ref->GetAxis(2)->Set(bins_trkmaxjet[2],MultCentbins);
// Extra binning
NJets->GetAxis(2)->Set(bins_NJETS[2],Extrabins);
trackmaxptinjethisto->GetAxis(2)->Set(bins_trkmax[2],Extrabins);
//Sumw2
NJets->Sumw2();
trackmaxptinjethisto->Sumw2();
jettrackmaxptinjethisto->Sumw2();
jettrackmaxptinjethisto_no0->Sumw2();
jettrackmaxptinjethisto_ref->Sumw2();

// --> Rho underlying events subtraction histos
// Multiplicity binning
histo_jetUE->GetAxis(1)->Set(bins_UE[1],MultCentbins);
histo_jetAverageRho->GetAxis(1)->Set(bins_UE[1],MultCentbins);
histo_jetcheckcorrection->GetAxis(1)->Set(bins_UE[1],MultCentbins);
// Extra binning
histo_jetUE->GetAxis(2)->Set(bins_UE[2],Extrabins);
histo_jetAverageRho->GetAxis(2)->Set(bins_UE[2],Extrabins);
histo_jetcheckcorrection->GetAxis(2)->Set(bins_UE[2],Extrabins);
//Sumw2
histo_jetUE->Sumw2();
histo_jetAverageRho->Sumw2();
histo_jetcheckcorrection->Sumw2();

// --> Jet correlations to EP histos
// Multiplicity binning
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
// Extra binning
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
// Sumw2
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

// --> Jet information: inclusive, leading, subleading and third jets histos
// Multiplicity binning
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
hist_reco_thrdproj_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_thrdproj_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_ref_thrdproj_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_reco_thrdprojdiff_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_gen_thrdprojdiff_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
hist_ref_thrdprojdiff_weighted->GetAxis(3)->Set(bins_jet[3],MultCentbins);
// Extra binning
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
hist_reco_thrdproj_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_thrdproj_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_ref_thrdproj_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_reco_thrdprojdiff_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_gen_thrdprojdiff_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
hist_ref_thrdprojdiff_weighted->GetAxis(4)->Set(bins_jet[4],Extrabins);
// Sumw2
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
hist_reco_thrdproj_weighted->Sumw2();
hist_gen_thrdproj_weighted->Sumw2();
hist_ref_thrdproj_weighted->Sumw2();
hist_reco_thrdprojdiff_weighted->Sumw2();
hist_gen_thrdprojdiff_weighted->Sumw2();
hist_ref_thrdprojdiff_weighted->Sumw2();

// --> Jet information used in the correlations histos
// Multiplicity binning
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
// Extra binning
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
//Sumw2
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

// --> JES studies histos
// Multiplicity binning
hist_jes_reco_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_jes_reco_fromB_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_leadjes_reco_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_leadjes_reco_fromB_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_subleadjes_reco_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
hist_subleadjes_reco_fromB_weighted->GetAxis(4)->Set(bins_jes[4],MultCentbins);
// Extra binning
hist_jes_reco_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_jes_reco_fromB_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_leadjes_reco_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_leadjes_reco_fromB_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_subleadjes_reco_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
hist_subleadjes_reco_fromB_weighted->GetAxis(5)->Set(bins_jes[5],Extrabins);
// Sumw2
hist_jes_reco_weighted->Sumw2();
hist_jes_reco_fromB_weighted->Sumw2();
hist_leadjes_reco_weighted->Sumw2();
hist_leadjes_reco_fromB_weighted->Sumw2();
hist_subleadjes_reco_weighted->Sumw2();
hist_subleadjes_reco_fromB_weighted->Sumw2();

// --> Xj search for quenching histos
// Xj binning
hist_reco_lead_reco_subl_quench->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_fake_lead_reco_subl_quench->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_reco_lead_fake_subl_quench->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_fake_lead_fake_subl_quench->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_gen_lead_gen_subl_quench->GetAxis(0)->Set(bins_quenc[0],XjBins);
hist_ref_lead_ref_subl_quench->GetAxis(0)->Set(bins_quenc[0],XjBins);
// Aj binning
hist_reco_lead_reco_subl_quench->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_fake_lead_reco_subl_quench->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_reco_lead_fake_subl_quench->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_fake_lead_fake_subl_quench->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_gen_lead_gen_subl_quench->GetAxis(1)->Set(bins_quenc[1],AjBins);
hist_ref_lead_ref_subl_quench->GetAxis(1)->Set(bins_quenc[1],AjBins);
// DPhi binning
hist_reco_lead_reco_subl_quench->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_fake_lead_reco_subl_quench->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_reco_lead_fake_subl_quench->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_fake_lead_fake_subl_quench->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_gen_lead_gen_subl_quench->GetAxis(2)->Set(bins_quenc[2],DphiBins);
hist_ref_lead_ref_subl_quench->GetAxis(2)->Set(bins_quenc[2],DphiBins);
// Multiplicity binning
hist_reco_lead_reco_subl_quench->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_fake_lead_reco_subl_quench->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_reco_lead_fake_subl_quench->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_fake_lead_fake_subl_quench->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_gen_lead_gen_subl_quench->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
hist_ref_lead_ref_subl_quench->GetAxis(3)->Set(bins_quenc[3],MultCentbins);
// Pt average binning
hist_reco_lead_reco_subl_quench->GetAxis(4)->Set(bins_quenc[4],PtLSLBins);
hist_fake_lead_reco_subl_quench->GetAxis(4)->Set(bins_quenc[4],PtLSLBins);
hist_reco_lead_fake_subl_quench->GetAxis(4)->Set(bins_quenc[4],PtLSLBins);
hist_fake_lead_fake_subl_quench->GetAxis(4)->Set(bins_quenc[4],PtLSLBins);
hist_gen_lead_gen_subl_quench->GetAxis(4)->Set(bins_quenc[4],PtLSLBins);
hist_ref_lead_ref_subl_quench->GetAxis(4)->Set(bins_quenc[4],PtLSLBins);
// Extra binning
hist_reco_lead_reco_subl_quench->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_fake_lead_reco_subl_quench->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_reco_lead_fake_subl_quench->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_fake_lead_fake_subl_quench->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_gen_lead_gen_subl_quench->GetAxis(5)->Set(bins_quenc[5],Extrabins);
hist_ref_lead_ref_subl_quench->GetAxis(5)->Set(bins_quenc[5],Extrabins);
// Leading jet pT binning
hist_reco_lead_reco_subl_quench->GetAxis(6)->Set(bins_quenc[6],PtLSLBins);
hist_fake_lead_reco_subl_quench->GetAxis(6)->Set(bins_quenc[6],PtLSLBins);
hist_reco_lead_fake_subl_quench->GetAxis(6)->Set(bins_quenc[6],PtLSLBins);
hist_fake_lead_fake_subl_quench->GetAxis(6)->Set(bins_quenc[6],PtLSLBins);
hist_gen_lead_gen_subl_quench->GetAxis(6)->Set(bins_quenc[6],PtLSLBins);
hist_ref_lead_ref_subl_quench->GetAxis(6)->Set(bins_quenc[6],PtLSLBins);
// Sub Leading jet pT binning
hist_reco_lead_reco_subl_quench->GetAxis(7)->Set(bins_quenc[7],PtLSLBins);
hist_fake_lead_reco_subl_quench->GetAxis(7)->Set(bins_quenc[7],PtLSLBins);
hist_reco_lead_fake_subl_quench->GetAxis(7)->Set(bins_quenc[7],PtLSLBins);
hist_fake_lead_fake_subl_quench->GetAxis(7)->Set(bins_quenc[7],PtLSLBins);
hist_gen_lead_gen_subl_quench->GetAxis(7)->Set(bins_quenc[7],PtLSLBins);
hist_ref_lead_ref_subl_quench->GetAxis(7)->Set(bins_quenc[7],PtLSLBins);
//Sumw2
hist_reco_lead_reco_subl_quench->Sumw2();
hist_fake_lead_reco_subl_quench->Sumw2();
hist_reco_lead_fake_subl_quench->Sumw2();
hist_fake_lead_fake_subl_quench->Sumw2();
hist_gen_lead_gen_subl_quench->Sumw2();
hist_ref_lead_ref_subl_quench->Sumw2();

// --> Leading and Subleading for EP studies histos
// Xj binning
hist_reco_leadEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_leadEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_reco_sublEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_sublEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_ref_leadEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_leadEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_gen_sublEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_fake_leadEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_fake_leadEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_fake_sublEP_quench_plus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
hist_fake_sublEP_quench_minus->GetAxis(0)->Set(bins_quencEP[0],XjBins);
// Dphi binning
hist_reco_leadEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_leadEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_reco_sublEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_sublEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_ref_leadEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_leadEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_gen_sublEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_fake_leadEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_fake_leadEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_fake_sublEP_quench_plus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
hist_fake_sublEP_quench_minus->GetAxis(1)->Set(bins_quencEP[1],DphiBins);
// Multiplicity binning
hist_reco_leadEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_leadEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_reco_sublEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_sublEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_ref_leadEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_leadEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_gen_sublEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_fake_leadEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_fake_leadEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_fake_sublEP_quench_plus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
hist_fake_sublEP_quench_minus->GetAxis(5)->Set(bins_quencEP[5],MultCentbins);
// Pt average binning
hist_reco_leadEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_reco_leadEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_reco_sublEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_reco_sublEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_ref_sublEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_ref_sublEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_ref_leadEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_ref_leadEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_gen_leadEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_gen_leadEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_gen_sublEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_gen_sublEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_fake_leadEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_fake_leadEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_fake_sublEP_quench_plus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
hist_fake_sublEP_quench_minus->GetAxis(6)->Set(bins_quencEP[6],PtLSLBins);
// Extra binning
hist_reco_leadEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_reco_leadEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_reco_sublEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_reco_sublEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_ref_sublEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_ref_sublEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_ref_leadEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_ref_leadEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_gen_leadEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_gen_leadEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_gen_sublEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_gen_sublEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_fake_leadEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_fake_leadEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_fake_sublEP_quench_plus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
hist_fake_sublEP_quench_minus->GetAxis(7)->Set(bins_quencEP[7],Extrabins);
// Sumw2
hist_reco_leadEP_quench_plus->Sumw2();
hist_reco_leadEP_quench_minus->Sumw2();
hist_reco_sublEP_quench_plus->Sumw2();
hist_reco_sublEP_quench_minus->Sumw2();
hist_ref_sublEP_quench_plus->Sumw2();
hist_ref_sublEP_quench_minus->Sumw2();
hist_ref_leadEP_quench_plus->Sumw2();
hist_ref_leadEP_quench_minus->Sumw2();
hist_gen_leadEP_quench_plus->Sumw2();
hist_gen_leadEP_quench_minus->Sumw2();
hist_gen_sublEP_quench_plus->Sumw2();
hist_gen_sublEP_quench_minus->Sumw2();
hist_fake_leadEP_quench_plus->Sumw2();
hist_fake_leadEP_quench_minus->Sumw2();
hist_fake_sublEP_quench_plus->Sumw2();
hist_fake_sublEP_quench_minus->Sumw2();

// --> Eta dijet histos
// Xj binning
hist_etaDijet_reco->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_CM_reco->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_ref->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_CM_ref->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_gen->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_etaDijet_CM_gen->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_yDijet_CM_reco->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_yDijet_CM_ref->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
hist_yDijet_CM_gen->GetAxis(2)->Set(bins_etaDijet[2],XjBins);
// Aj binning
hist_etaDijet_reco->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_CM_reco->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_ref->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_CM_ref->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_gen->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_etaDijet_CM_gen->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_yDijet_CM_reco->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_yDijet_CM_ref->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
hist_yDijet_CM_gen->GetAxis(3)->Set(bins_etaDijet[3],AjBins);
// DPhi binning
hist_etaDijet_reco->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_CM_reco->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_ref->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_CM_ref->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_gen->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_etaDijet_CM_gen->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_yDijet_CM_reco->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_yDijet_CM_ref->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
hist_yDijet_CM_gen->GetAxis(4)->Set(bins_etaDijet[4],DphiBins);
/*
// XPb binning
hist_etaDijet_reco->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_CM_reco->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_ref->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_CM_ref->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_gen->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_etaDijet_CM_gen->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_yDijet_CM_reco->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_yDijet_CM_ref->GetAxis(5)->Set(bins_etaDijet[5],XBins);
hist_yDijet_CM_gen->GetAxis(5)->Set(bins_etaDijet[5],XBins);
// Xp binning
hist_etaDijet_reco->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_CM_reco->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_ref->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_CM_ref->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_gen->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_etaDijet_CM_gen->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_yDijet_CM_reco->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_yDijet_CM_ref->GetAxis(6)->Set(bins_etaDijet[6],XBins);
hist_yDijet_CM_gen->GetAxis(6)->Set(bins_etaDijet[6],XBins);
*/
// Multiplicity binning
hist_etaDijet_reco->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_etaDijet_CM_reco->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_etaDijet_ref->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_etaDijet_CM_ref->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_etaDijet_gen->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_etaDijet_CM_gen->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_yDijet_CM_reco->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_yDijet_CM_ref->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
hist_yDijet_CM_gen->GetAxis(7)->Set(bins_etaDijet[7],MultCentbins);
// Pt average binning
hist_etaDijet_reco->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_etaDijet_CM_reco->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_etaDijet_ref->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_etaDijet_CM_ref->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_etaDijet_gen->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_etaDijet_CM_gen->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_yDijet_CM_reco->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_yDijet_CM_ref->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
hist_yDijet_CM_gen->GetAxis(8)->Set(bins_etaDijet[8],PtLSLBins);
// Extra binning
hist_etaDijet_reco->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_etaDijet_CM_reco->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_etaDijet_ref->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_etaDijet_CM_ref->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_etaDijet_gen->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_etaDijet_CM_gen->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_yDijet_CM_reco->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_yDijet_CM_ref->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
hist_yDijet_CM_gen->GetAxis(9)->Set(bins_etaDijet[9],Extrabins);
// Sumw2
hist_etaDijet_reco->Sumw2();
hist_etaDijet_CM_reco->Sumw2();
hist_etaDijet_ref->Sumw2();
hist_etaDijet_CM_ref->Sumw2();
hist_etaDijet_gen->Sumw2();
hist_etaDijet_CM_gen->Sumw2();
hist_yDijet_CM_reco->Sumw2();
hist_yDijet_CM_ref->Sumw2();
hist_yDijet_CM_gen->Sumw2();

// --> In jet histos
// Multiplicity binning
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
// Extra binning
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
// Sumw2
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

// --> Jet-Track correlation histos
// Trk pt binning
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
// Multiplicity binning
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
// Extra binning
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
// Sumw2
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

// --> Two particle correlation histos
//Trk pT binning
hist_reco_reco_2pcorrelation_signal->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_reco_reco_2pcorrelation_signal_subg0->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_reco_reco_2pcorrelation_signal_subcross->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_reco_reco_2pcorrelation_mixing->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_signal->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_signal_subg0->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_signal_subcross->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
hist_gen_gen_2pcorrelation_mixing->GetAxis(2)->Set(bins_2pc[2],TrkPtbins);
// Multiplicity binning
hist_reco_reco_2pcorrelation_signal->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_reco_reco_2pcorrelation_signal_subg0->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_reco_reco_2pcorrelation_signal_subcross->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_reco_reco_2pcorrelation_mixing->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_signal->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_signal_subg0->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_signal_subcross->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
hist_gen_gen_2pcorrelation_mixing->GetAxis(3)->Set(bins_2pc[3],MultCentbins);
// Extra binning
hist_reco_reco_2pcorrelation_signal->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_reco_reco_2pcorrelation_signal_subg0->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_reco_reco_2pcorrelation_signal_subcross->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_reco_reco_2pcorrelation_mixing->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_signal->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_signal_subg0->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_signal_subcross->GetAxis(4)->Set(bins_2pc[4],Extrabins);
hist_gen_gen_2pcorrelation_mixing->GetAxis(4)->Set(bins_2pc[4],Extrabins);
// Sumw2
hist_reco_reco_2pcorrelation_signal->Sumw2();
hist_reco_reco_2pcorrelation_signal_subg0->Sumw2();
hist_reco_reco_2pcorrelation_signal_subcross->Sumw2();
hist_reco_reco_2pcorrelation_mixing->Sumw2();
hist_gen_gen_2pcorrelation_signal->Sumw2();
hist_gen_gen_2pcorrelation_signal_subg0->Sumw2();
hist_gen_gen_2pcorrelation_signal_subcross->Sumw2();
hist_gen_gen_2pcorrelation_mixing->Sumw2();


// --> Unfolding histos
// -> Xj vs Pt average
// Unfold binning
fhUnfoldingResponse_xjptave->SetBinEdges(0, fullUnfoldingBinning_xjptave);
fhUnfoldingResponse_xjptave->SetBinEdges(1, fullUnfoldingBinning_xjptave);
fhUnfoldingMeasu_xjptave->SetBinEdges(0, fullUnfoldingBinning_xjptave);
fhUnfoldingTruthRef_xjptave->SetBinEdges(0, fullUnfoldingBinning_xjptave);
fhUnfoldingTruthGen_xjptave->SetBinEdges(0, fullUnfoldingBinning_xjptave);
// DPhi binning
fhUnfoldingResponse_xjptave->GetAxis(2)->Set(bins_unfxjptave[2],DphiBins);
fhUnfoldingResponse_xjptave->GetAxis(3)->Set(bins_unfxjptave[3],DphiBins);
fhUnfoldingMeasu_xjptave->GetAxis(1)->Set(bins_unfxjptaveMT[1],DphiBins);
fhUnfoldingTruthRef_xjptave->GetAxis(1)->Set(bins_unfxjptaveMT[1],DphiBins);
fhUnfoldingTruthGen_xjptave->GetAxis(1)->Set(bins_unfxjptaveMT[1],DphiBins);
// Multiplicity binning
fhUnfoldingResponse_xjptave->GetAxis(6)->Set(bins_unfxjptave[6],MultCentbins);
fhUnfoldingMeasu_xjptave->GetAxis(3)->Set(bins_unfxjptaveMT[3],MultCentbins);
fhUnfoldingTruthRef_xjptave->GetAxis(3)->Set(bins_unfxjptaveMT[3],MultCentbins);
fhUnfoldingTruthGen_xjptave->GetAxis(3)->Set(bins_unfxjptaveMT[3],MultCentbins);
// Sumw2
fhUnfoldingResponse_xjptave->Sumw2();
fhUnfoldingMeasu_xjptave->Sumw2();
fhUnfoldingTruthRef_xjptave->Sumw2();
fhUnfoldingTruthGen_xjptave->Sumw2();

// -> Pt1 vs Pt2
// Unfold binning
fhUnfoldingResponse_pt1pt2->SetBinEdges(0, fullUnfoldingBinning_pt1pt2);
fhUnfoldingResponse_pt1pt2->SetBinEdges(1, fullUnfoldingBinning_pt1pt2);
fhUnfoldingMeasu_pt1pt2->SetBinEdges(0, fullUnfoldingBinning_pt1pt2);
fhUnfoldingTruthRef_pt1pt2->SetBinEdges(0, fullUnfoldingBinning_pt1pt2);
fhUnfoldingTruthGen_pt1pt2->SetBinEdges(0, fullUnfoldingBinning_pt1pt2);
// DPhi binning
fhUnfoldingResponse_pt1pt2->GetAxis(2)->Set(bins_unfpt1pt2[2],DphiBins);
fhUnfoldingResponse_pt1pt2->GetAxis(3)->Set(bins_unfpt1pt2[3],DphiBins);
fhUnfoldingMeasu_pt1pt2->GetAxis(1)->Set(bins_unfpt1pt2MT[1],DphiBins);
fhUnfoldingTruthRef_pt1pt2->GetAxis(1)->Set(bins_unfpt1pt2MT[1],DphiBins);
fhUnfoldingTruthGen_pt1pt2->GetAxis(1)->Set(bins_unfpt1pt2MT[1],DphiBins);
// Multiplicity binning
fhUnfoldingResponse_pt1pt2->GetAxis(6)->Set(bins_unfpt1pt2[6],MultCentbins);
fhUnfoldingMeasu_pt1pt2->GetAxis(3)->Set(bins_unfpt1pt2MT[3],MultCentbins);
fhUnfoldingTruthRef_pt1pt2->GetAxis(3)->Set(bins_unfpt1pt2MT[3],MultCentbins);
fhUnfoldingTruthGen_pt1pt2->GetAxis(3)->Set(bins_unfpt1pt2MT[3],MultCentbins);
// Sumw2
fhUnfoldingResponse_pt1pt2->Sumw2();
fhUnfoldingMeasu_pt1pt2->Sumw2();
fhUnfoldingTruthRef_pt1pt2->Sumw2();
fhUnfoldingTruthGen_pt1pt2->Sumw2();

// -> xj vs Pt2 in bins of Pt1
// Unfold binning
fhUnfoldingResponse_xjpt2->SetBinEdges(0, fullUnfoldingBinning_xjpt2);
fhUnfoldingResponse_xjpt2->SetBinEdges(1, fullUnfoldingBinning_xjpt2);
fhUnfoldingMeasu_xjpt2->SetBinEdges(0, fullUnfoldingBinning_xjpt2);
fhUnfoldingTruthRef_xjpt2->SetBinEdges(0, fullUnfoldingBinning_xjpt2);
fhUnfoldingTruthGen_xjpt2->SetBinEdges(0, fullUnfoldingBinning_xjpt2);
// DPhi binning
fhUnfoldingResponse_xjpt2->GetAxis(2)->Set(bins_unfxjpt2[2],DphiBins);
fhUnfoldingResponse_xjpt2->GetAxis(3)->Set(bins_unfxjpt2[3],DphiBins);
fhUnfoldingMeasu_xjpt2->GetAxis(1)->Set(bins_unfxjpt2MT[1],DphiBins);
fhUnfoldingTruthRef_xjpt2->GetAxis(1)->Set(bins_unfxjpt2MT[1],DphiBins);
fhUnfoldingTruthGen_xjpt2->GetAxis(1)->Set(bins_unfxjpt2MT[1],DphiBins);
// Multiplicity binning
fhUnfoldingResponse_xjpt2->GetAxis(6)->Set(bins_unfxjpt2[6],MultCentbins);
fhUnfoldingMeasu_xjpt2->GetAxis(3)->Set(bins_unfxjpt2MT[3],MultCentbins);
fhUnfoldingTruthRef_xjpt2->GetAxis(3)->Set(bins_unfxjpt2MT[3],MultCentbins);
fhUnfoldingTruthGen_xjpt2->GetAxis(3)->Set(bins_unfxjpt2MT[3],MultCentbins);
// Pt1 binning
fhUnfoldingResponse_xjpt2->GetAxis(7)->Set(bins_unfxjpt2[7],PtLSLBins);
fhUnfoldingResponse_xjpt2->GetAxis(8)->Set(bins_unfxjpt2[8],PtLSLBins);
fhUnfoldingMeasu_xjpt2->GetAxis(4)->Set(bins_unfxjpt2MT[4],PtLSLBins);
fhUnfoldingTruthRef_xjpt2->GetAxis(4)->Set(bins_unfxjpt2MT[4],PtLSLBins);
fhUnfoldingTruthGen_xjpt2->GetAxis(4)->Set(bins_unfxjpt2MT[4],PtLSLBins);
// Sumw2
fhUnfoldingResponse_xjpt2->Sumw2();
fhUnfoldingMeasu_xjpt2->Sumw2();
fhUnfoldingTruthRef_xjpt2->Sumw2();
fhUnfoldingTruthGen_xjpt2->Sumw2();


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
	Nev_alljetfromalltrk->Write();
	Nev_jetwithlowpttrk->Write();
	Nev_jetfromonetrk->Write();
	Nev_jetsfrombothlowpttrkandonetrk->Write();
	Nev_jetwithlowpttrk_lead->Write();
	Nev_jetwithlowpttrk_sublead->Write();
	Nev_jetfromonetrk_lead->Write();
	Nev_jetfromonetrk_sublead->Write();

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
	multiplicity_etadep_weighted->Write();
	multiplicity_nocut->Write();
	multiplicity_corrected->Write();
	multiplicity_tight->Write();
	multiplicity_loose->Write();
	multiplicity2D->Write();
 	multiplicity->Write();
	multiplicity_weighted->Write();
	multiplicity_withonejet_weighted->Write();
	multiplicity_withdijets_weighted->Write();
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
	hfhistEta4Sum_nocut->Write();
	hfhistEta4Plus_nocut->Write();
	hfhistEta4Minus_nocut->Write();
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
	jettrackmaxptinjethisto_no0->Write();
	if(isMC)jettrackmaxptinjethisto_ref->Write();
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
	hist_reco_thrdproj_weighted->Write();
	hist_reco_thrdprojdiff_weighted->Write();
	if(isMC){
		hist_jes_reco_weighted->Write();
		hist_jes_reco_fromB_weighted->Write();
		hist_leadjes_reco_weighted->Write();
		hist_leadjes_reco_fromB_weighted->Write();
		hist_subleadjes_reco_weighted->Write();
		hist_subleadjes_reco_fromB_weighted->Write();
		hist_gen_jet->Write();
		hist_gen_jet_weighted->Write();
		hist_gen_leadjet->Write();
		hist_gen_leadjet_weighted->Write();
		hist_gen_subljet->Write();
		hist_gen_subljet_weighted->Write();
		hist_gen_thrdjet_weighted->Write();
		hist_gen_thrdproj_weighted->Write();
		hist_gen_thrdprojdiff_weighted->Write();
		hist_ref_jet_weighted->Write();
		hist_ref_leadjet_weighted->Write();
		hist_ref_subljet_weighted->Write();
		hist_ref_thrdjet_weighted->Write();
		hist_ref_thrdproj_weighted->Write();
		hist_ref_thrdprojdiff_weighted->Write();
	}
	
	jetjet_Dphi_Deta_reco->Write();	if(isMC) { jetjet_Dphi_Deta_ref->Write(); jetjet_Dphi_Deta_gen->Write();}
	jetjet_Dphi_Deta_reco_diff->Write();	if(isMC) { jetjet_Dphi_Deta_ref_diff->Write(); jetjet_Dphi_Deta_gen_diff->Write();}
	
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
	hist_reco_lead_reco_subl_quench->Write();
	hist_reco_leadEP_quench_plus->Write();
	hist_reco_leadEP_quench_minus->Write();
	hist_etaDijet_reco->Write();
	hist_etaDijet_CM_reco->Write();
	hist_yDijet_CM_reco->Write();
	if(isMC){
		hist_fake_lead_reco_subl_quench->Write();
		hist_reco_lead_fake_subl_quench->Write();
		hist_fake_lead_fake_subl_quench->Write();
		hist_gen_lead_gen_subl_quench->Write();
		hist_fake_leadEP_quench_plus->Write();
		hist_fake_leadEP_quench_minus->Write();
		hist_fake_sublEP_quench_plus->Write();
		hist_fake_sublEP_quench_minus->Write();
		hist_gen_leadEP_quench_plus->Write();
		hist_gen_leadEP_quench_minus->Write();
		hist_ref_lead_ref_subl_quench->Write();
		hist_ref_leadEP_quench_plus->Write();
		hist_ref_leadEP_quench_minus->Write();
		hist_etaDijet_ref->Write();
		hist_etaDijet_CM_ref->Write();
		hist_yDijet_CM_ref->Write();
		hist_etaDijet_gen->Write();
		hist_etaDijet_CM_gen->Write();
		hist_yDijet_CM_gen->Write();
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

	fhUnfoldingResponse_xjptave->Write();
	fhUnfoldingMeasu_xjptave->Write();
	fhUnfoldingTruthRef_xjptave->Write();
	fhUnfoldingTruthGen_xjptave->Write();

	fhUnfoldingResponse_pt1pt2->Write();
	fhUnfoldingMeasu_pt1pt2->Write();
	fhUnfoldingTruthRef_pt1pt2->Write();
	fhUnfoldingTruthGen_pt1pt2->Write();

	fhUnfoldingResponse_xjpt2->Write();
	fhUnfoldingMeasu_xjpt2->Write();
	fhUnfoldingTruthRef_xjpt2->Write();
	fhUnfoldingTruthGen_xjpt2->Write();

}
