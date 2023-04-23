#include "call_libraries.h"  // call libraries from ROOT and C++
#include "particleid.h"  // call for particle id

// Input quatities used the codes

TString colliding_system = "pPb"; // use one of this options = "pp", "pPb", "XeXe" and "PbPb" (OO and pO in the future)
int sNN_energy_GeV = 8160; //center of mass colliding energy (GeV)
int year_of_datataking = 2016;

bool do_CM_pPb = false; // do center-of-mass correction in pPb?
bool is_pgoing = false; // is p-going direction?

bool use_centrality = false; // only true for: "XeXe" and "PbPb" (but also can be set as false to see evolution with multiplicity)

float vz_cut_min = -10.0; //vz acceptance
float vz_cut_max = 10.0; //vz acceptance

const std::vector<double> multiplicity_centrality_bins{10.0, 400.0}; //multiplicity range
//event filters
std::vector<int> event_filter_bool; // event filter booleans
//std::vector<TString> event_filter_str{"pBeamScrapingFilter", "pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose"}; // event filters to be applied (pp ref - 2017)
std::vector<TString> event_filter_str{"pBeamScrapingFilter", "pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "phfCoincFilter", "pVertexFilterCutdz1p0"}; // event filters to be applied (pPb - 2016)
//std::vector<TString> event_filter_str{"pprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "collisionEventSelectionAOD", "phfCoincFilter2Th4", "pclusterCompatibilityFilter"}; // event filters to be applied (PbPb - 2018)
//std::vector<TString> event_filter_str{"pBeamScrapingFilter", "pPAprimaryVertexFilter", "HBHENoiseFilterResultRun2Loose", "phfCoincFilter", "pVertexFilterCutdz1p0"}; // event filters to be applied (XeXe - 2017)

// by default the code will calculated QA plots for jets and tracks, you can turn on or off the flags bellow based on your studies
// be carefull about memory usage (e.g. do not run incluse + fragmentation at same time)
bool do_inclusejettrack_correlation = true; // Inclusive jets + track correlation
bool do_leading_subleading_jettrack_correlation = true; // Leading jets + track correlation and Sub-Leading jets + track correlation
bool do_flow = false; // if true if makes correlation for Jet-Vn flow if false it multiply by trk pT to get jet shapes

//=========================================================

//============= Jet information =========================== 

TString jet_collection = "ak4PFJetAnalyzer"; // jet collection in forest
bool dojettrigger = false; // apply jet trigger
TString jet_trigger = "HLT_PAAK4PFJet80_Eta5p1_v3"; // jet trigger in forest 

float jet_pt_min_cut = 100.0; // jet min pT cut 
float jet_pt_max_cut = 8160.0; // jet max pT cut 
float jet_eta_min_cut = -1.6; // jet min eta cut 
float jet_eta_max_cut = 1.6; // jet min eta cut 

TString JEC_file = "Autumn16_HI_pPb_pgoing_Unembedded_MC_L2Relative_AK4PF.txt"; //JEC file
TString JEU_file = "Autumn16_HI_pPb_pgoing_Unembedded_MC_L2Relative_AK4PF.txt"; //JEU file (future)
bool use_WTA = true; // use WTA or E-Scheme axis 
float leading_subleading_deltaphi_min = (5./6.)*TMath::Pi(); //used for jet leading and subleading correlation and jet quenching analysis; use 0 for all;
float leading_pT_min = jet_pt_min_cut; //used for jet leading and subleading correlation and jet quenching analysis
float subleading_pT_min = 50.0; //used for jet leading and subleading correlation and jet quenching analysis

// if we want to make Xj or Aj selections [0,1] are full inclusive
bool do_Xj_or_Ajcut = false;
float xjmin = 0.0;//xj minimum
float xjmax = 1.0;//xj maximum
float Ajmin = 0.0;//Aj minimum
float Ajmax = 1.0;//Aj maximum

// for ejt smearing
bool do_jet_smearing = false; 

//=========================================================

//============= Track information =========================

const std::vector<double> trk_pt_bins{0.7, 1.0, 2.0, 3.0, 4.0, 8.0, 12.0, 300.0}; //trk pT bin range for correlations
float trk_eta_cut = 2.4; // trk +/- eta range
float trk_pt_resolution_cut = 0.1; // trk pt resolution cut
float trk_dca_xy_cut = 3.0; // trk XY DCA cut
float trk_dca_z_cut = 3.0; // trk Z DCA cut
float chi2_ndf_nlayer_cut = 0.18;  // trk chi2/ndf/nlayer cut
float calo_matching = 0.5; // trk calo matching cut 
int nhits = 11; // trk Nhits cut

float trk_pt_min_cut = trk_pt_bins[0]; // min track pT
TString trk_eff_file = "eff_table_p-going_HIJING.root"; //track efficiency table

//=========================================================

//============= Reference samples =========================

// use just one ref sample due memory issues

//--> Mixing ref. samples quantities
bool do_mixing = true; // use mixing method?
bool similar_events = true; // if true we consider only tracks coming for similar events (onl if jet requirement is satisfied), if false all tracks are used
int N_ev_mix = 20; // number of events to mix
int Mult_or_Cent_range = 100; // multiplicity or centrality interval allowed between event and mixed event
float DVz_range = 0.5;  // Vertex Z interval allowed between event and mixed event

//--> rotation ref. samples quantities
bool do_rotation = false; // use rotation method?
int N_of_rot = N_ev_mix; // setup number of rotations

//=========================================================

//For MC only
bool double_weight_mix = false; // double weighting in the mixing
bool do_pid = false; // apply PID? // choose the value between [] based on particleid.h
int particlepid = pid[Pion];   
TString particles = pid_str[Pion];

/*
Print out the inputs
--> Arguments
data_or_mc: MC for Monte Carlo and Data from data
fileeff: efficiency file
coll_system: colliding system
*/
void print_input(TString data_or_mc, TFile *fileeff, TString coll_system, float pthatmin, float pthatmax){
	cout << "From input:" << endl;
	cout << endl;
	cout << "Running over " << data_or_mc.Data() << endl;
	cout << "Colliding system: " << colliding_system.Data() << endl;
	cout << "Colliding energy: " << sNN_energy_GeV/1000. << " TeV"<< endl;
	cout << "Data taking in  " << year_of_datataking << endl;
	cout << "Event filters applied: {"; for(int a=0;a<event_filter_str.size();a++){if(a==event_filter_str.size()-1){cout << "" << event_filter_str[a] << "}";}else{cout << "" << event_filter_str[a] << ",";}} cout << endl;
	cout << "Vz acceptance: " << vz_cut_min << " < Vz < " << vz_cut_max << " cm" << endl; 
	if(!use_centrality){cout << "Multiplicity bins: {"; for(int a=0;a<multiplicity_centrality_bins.size();a++){if(a==multiplicity_centrality_bins.size()-1){cout << "" << multiplicity_centrality_bins[a] << "}";}else{cout << "" << multiplicity_centrality_bins[a] << ",";}} cout << endl; }
	if(use_centrality){cout << "Centrality bins: {"; for(int a=0;a<multiplicity_centrality_bins.size();a++){if(a==multiplicity_centrality_bins.size()-1){cout << "" << multiplicity_centrality_bins[a] << "}";}else{cout << "" << multiplicity_centrality_bins[a] << ",";}} cout << endl; }
	if(do_CM_pPb) cout << "Center-of-mass correction for side: " << endl;
	cout << endl;
	cout << "=========== Jets ===========" << endl;
	cout << endl;
	cout << "Jet collection: " << jet_collection.Data() << endl;
	if(dojettrigger){cout << "Jet trigger: " << jet_trigger.Data() << endl;}else{cout << "No jet trigger applied!" << endl;}
	if(use_WTA){cout << "Using WTA jet axis" << endl;}else{cout << "Using E-Scheme jet axis" << endl;}
	cout << "Jet eta range: [" << jet_eta_min_cut << "," << jet_eta_max_cut << "]" << endl;
	cout << "Jet pT min: " << jet_pt_min_cut << " GeV"<< endl;
	cout << "Jet pT max: " << jet_pt_max_cut << " GeV"<< endl;
	cout << "JEC file: " << JEC_file.Data() << endl;
	cout << "JEU file: " << JEU_file.Data() << endl;
	if(do_jet_smearing) cout << "Smearing applied" << endl;
	if(do_Xj_or_Ajcut){
		cout << "Jet Xj min: " << xjmin << endl;
		cout << "Jet Xj max: " << xjmax << endl;
		cout << "Jet Aj min: " << Ajmin << endl;
		cout << "Jet Aj max: " << Ajmax << endl;
	}
	if(do_leading_subleading_jettrack_correlation){
		cout << "Leading jet pT min: " << leading_pT_min  << " GeV"<< endl;
		cout << "Sub-Leading jet pT min: " << subleading_pT_min  << " GeV"<< endl;
		cout << "Delta phi (leading,subleading) jets > " << leading_subleading_deltaphi_min << endl;
	}
	if(do_pthatcut){
		cout << "pThat min: " << pthatmin  << " GeV"<< endl;
		cout << "pThat max: " << pthatmax  << " GeV"<< endl;
	}
	cout << endl;
	cout << "=========== Tracks/Particles ===========" << endl;
	cout << endl;
	cout << "Track eta range: [" << -trk_eta_cut << "," << trk_eta_cut << "]" << endl;
	cout << "Reco track pT resolution < " << trk_pt_resolution_cut*100 << "%" << endl;
	cout << "Reco track XY DCA significance < " << trk_dca_xy_cut << endl;
	cout << "Reco track Z DCA significance < " << trk_dca_z_cut << endl;
	if(coll_system=="PbPb" || coll_system=="XeXe"){
		cout << "Calo calo_matching: for pT > 20 GeV, ET/pT > " << calo_matching << endl;
		cout << "Reco track chi2/ndf/nlayer < " << chi2_ndf_nlayer_cut << endl;
		cout << "Reco track number of hits >= " << nhits << endl;
		if(coll_system=="PbPb" && sNN_energy_GeV==5020 && year_of_datataking==2018)	cout << "For reco track algorithm 6 remove MVA values less than 0.98" << endl;
	}
	cout << "Track pt bins for correlations: {"; for(int a=0;a<trk_pt_bins.size();a++){if(a==trk_pt_bins.size()-1){cout << "" << trk_pt_bins[a] << "} GeV";}else{cout << "" << trk_pt_bins[a] << ",";}} cout << endl;
	if(do_pid){
		cout << "Using PDG ID:" << particlepid << " --> " << particles.Data() << endl;
	}
	// track or particle efficiency file --> adapt as needed
	cout << endl;
	if(!fileeff->IsOpen()){cout << "Cannot find the track efficiency file. Force quit!" << endl; return;}else{cout << "Track/particle efficiency file: " << trk_eff_file.Data() << endl;}
	cout << endl;
	cout << "=========== Reference Samples ===========" << endl;
	if(do_mixing){
		cout << endl;
		cout << "Reference sample: mixing" << endl;
		cout << "Number of events to mix = " << N_ev_mix << endl;
   		cout << "Multiplicity or Centrality range = " << Mult_or_Cent_range << endl;
  		cout << "Delta Vz range = " << DVz_range << endl;
  		if(similar_events) cout << "Using similar events!" << endl;
  	}
	if(do_rotation){
		cout << endl;
		cout << "Reference sample: rotation" << endl;
		cout << "Number of rotations = " << N_of_rot << endl;
	}
	if(!do_mixing && !do_rotation){
		cout << endl;
		cout << "No reference sample used!" << endl;
	}
	cout << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "+          This code will return the histograms for            +" << endl;
	cout << "+ --> Jets and tracks QA quantities                            +" << endl;
	cout << "+ --> Jet quenching observables                                +" << endl;
	if(do_inclusejettrack_correlation) cout << "+ --> Inclusive Jets + tracks (or particle) correlations       +" << endl;
	if(do_leading_subleading_jettrack_correlation) cout << "+ --> Leading Jets + tracks (or particle) correlations         +" << endl;
	if(do_leading_subleading_jettrack_correlation) cout << "+ --> Sub-Leading Jets + tracks (or particle) correlations     +" << endl;
	if(do_flow) cout << "+                          For Flow                            +" << endl;
	if(!do_flow) cout << "+                       For Jet Shapes                        +" << endl;
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}


