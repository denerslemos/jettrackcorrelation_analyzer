#include "call_libraries.h"  // call libraries from ROOT and C++
#include "read_tree.h" // read the TChains
#include "vector_definition.h"  // define the vectors for mixing
#include "histogram_definition.h" // define histograms
#include "random_mixing.h" // random mixing function
#include "uiclogo.h" // print UIC jets and start/stop time
#include "JetCorrector.h" // reader for JEC
#include "JetUncertainty.h" // reader for JEU
const double pPbRapidityBoost = 0.4654094531;

/*
Main code to run Jet+Track correlation

Written by Dener Lemos (dener.lemos@cern.ch)

--> Arguments
input_file: text file with a list of root input files: Forest or Skims
ouputfilename: just a text to run on Condor (can include the path for output here)
MCSim: 0 for data and > 0 for MC
pthatmin: pthat min cut for MC only
pthatmax: pthat max cut for MC only
*/
void jettrackcorrelation_analyzer(TString input_file, TString ouputfilename, int MCSim, float pthatmin, float pthatmax){

	TApplication *a = new TApplication("a", 0, 0); // avoid issues with corrupted files

	clock_t sec_start, sec_end, sec_start_mix, sec_end_mix; // For timing 
	sec_start = clock(); // start timing measurement
	TDatime* date = new TDatime(); // to add date in the output file
	printwelcome(true); // welcome message
	print_start(); // start timing print
	bool is_MC; if(MCSim==0){is_MC = false;}else{is_MC = true;} // boolean for MC or data
	bool do_pthatcut = true; // always true for MC
	if(!is_MC) do_pthatcut = false; // MC only
	if(!is_MC) do_pid = false; // MC only
	if(!do_pid) particles = "";
	if(colliding_system!="pPb") do_CM_pPb = false; // Only do center-of-mass for pPb (or future asymmetric systems like pO)

	//print important informations in the output file
	TString data_or_mc;
	if(!is_MC){data_or_mc="Data";}else{data_or_mc="MC";}
	if(colliding_system == "pPb" && is_pgoing && invert_pgoing){data_or_mc+="_invertside";}
	TString simev; if(similar_events){simev = "simevs";}else{simev = "";}
	TString ref_sample = "norefsample"; if(do_mixing && !do_rotation){ref_sample = Form("mix%ievsMult%iDVz%.1f%s",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data());}else if(!do_mixing && do_rotation){ref_sample = Form("rot%ievs",N_of_rot);}else if(do_mixing && do_rotation){ref_sample = Form("mix%ievsMult%iDVz%.1f%s_rot%ievs",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data(),N_of_rot);}
	TString jet_axis; if(use_WTA){jet_axis = "WTA";}else{jet_axis = "ESC";}
	TString smear; if(do_jet_smearing){smear = "_smearing_";}else{smear = "";}
	TString jet_type; if(do_inclusejettrack_correlation) jet_type += "Incl"; if(do_leading_subleading_jettrack_correlation) jet_type += Form("LeadSubl_%s",fwdbkw_jettrk_option.Data()); if(do_dijetstudies) jet_type += "Quench";
	TString XjAj; if(do_Xj_or_Ajcut){XjAj = Form("_Ajmin_%.1f_Ajmax_%.1f_Xjmin_%.1f_Xjmax_%.1f",Ajmin,Ajmax,xjmin,xjmax);}else{XjAj = "";}
	TString isflow;	if(do_flow){isflow="flow";}else{isflow="jetshape";}

	// In case of wrong input, printout error message and kill the job
	if(year_of_datataking!=2012 && year_of_datataking!=2016 && year_of_datataking!=2017 && year_of_datataking!=2018){cout << "Data and MC not supported: choose 2012 for pp at 8 TeV, 2016 for pPb at 8.16 TeV, 2017 for pp at 5.02 TeV or XeXe at 5.44 TeV and 2018 for PbPb at 5.02 TeV" << endl; return;}
	if(colliding_system!="pp" && colliding_system!="pPb" && colliding_system!="XeXe" && colliding_system!="PbPb"){cout << "Data and MC not supported: choose pp for proton-proton, pPb for proton-lead, PbPb for lead-lead and XeXe for xenon-xenon" << endl; return;}
	if(sNN_energy_GeV!=5020 && sNN_energy_GeV!=5440 && sNN_energy_GeV!=8000 && sNN_energy_GeV!=8160 && sNN_energy_GeV!=13000){cout << "Data and MC not supported: 5020 for pp 2017 or PbPb 2018, 5440 for XeXe, 8000 for pp 2018, 8160 for pPb 2016" << endl; return;}
	float sqrts = (float)sNN_energy_GeV;

	// Read Jet Energy Correction file
	vector<string> Files;
	Files.push_back(Form("aux_files/%s_%i/JEC/%s",colliding_system.Data(),sNN_energy_GeV,JEC_file.Data()));
	JetCorrector JEC(Files);

	// Track or particle efficiency file
	TFile *fileeff = TFile::Open(Form("aux_files/%s_%i/trk_eff_table/%s",colliding_system.Data(),sNN_energy_GeV,trk_eff_file.Data()));
	cout << endl;
	
	// Print the input in the screen/log 
	print_input(data_or_mc,fileeff,colliding_system,pthatmin,pthatmax);
	cout << endl;

	// Read the list of input file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(file_chain.c_str());}
	inputfile.close();

	// Read the trees to be added in the Chain
	TChain *hlt_tree = new TChain("hltanalysis/HltTree"); // for HLT trigger
	TChain *jet_tree = new TChain(Form("%s/t",jet_collection.Data())); // for jet collection
	TChain *trk_tree = new TChain("ppTrack/trackTree"); // for tracking
	TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree"); // event quantities
	TChain *gen_tree;
	if(is_MC){gen_tree = new TChain("HiGenParticleAna/hi");} // MC gen particles
	TChain *ski_tree = new TChain("skimanalysis/HltTree"); // event filters
	TChain *ep_tree;
	if(colliding_system=="pPb" && year_of_datataking==2016){ep_tree = new TChain("checkflattening/tree");} // event plane for 2016 pPb data
	
	// loop to add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		TFile *testfile = TFile::Open(*listIterator);
		if(!testfile || testfile->IsZombie() || testfile->TestBit(TFile::kRecovered)){cout << "File: " << *listIterator << " failed to open" << endl; continue;} // safety against corrupted files
		cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
		hlt_tree->Add(*listIterator);
		trk_tree->Add(*listIterator);
		hea_tree->Add(*listIterator);
		jet_tree->Add(*listIterator);
		ski_tree->Add(*listIterator);
		if(is_MC){gen_tree->Add(*listIterator);}
		if(colliding_system=="pPb" && year_of_datataking==2016){ep_tree->Add(*listIterator);}
	}
    file_name_vector.clear();
	
	// Connect all chains
	hlt_tree->AddFriend(trk_tree);
	hlt_tree->AddFriend(hea_tree);
	hlt_tree->AddFriend(jet_tree);
	hlt_tree->AddFriend(ski_tree);
	if(is_MC){hlt_tree->AddFriend(gen_tree);}
	if(colliding_system=="pPb" && year_of_datataking==2016){hlt_tree->AddFriend(ep_tree);}

    // Read the desired branchs in the trees
	read_tree(hlt_tree, is_MC, use_WTA, jet_trigger.Data(), colliding_system.Data(), sNN_energy_GeV, year_of_datataking, event_filter_str, event_filter_bool); // access the tree informations
	if(!dojettrigger) jet_trigger="nojettrig"; // just for output name
	
    // Use sumw2() to make sure about histogram uncertainties in ROOT
	sw2(); 

	int nevents = hlt_tree->GetEntries(); // number of events
	cout << "Total number of events in those files: "<< nevents << endl;
	cout << endl;
	cout << "-------------------------------------------------" << endl;
	
    //boosting
    double boost = 0;
    if(colliding_system=="pPb" && do_CM_pPb && is_pgoing){boost = pPbRapidityBoost;}else if(colliding_system=="pPb" && do_CM_pPb && !is_pgoing){boost = -pPbRapidityBoost;}
    // if(isMC) boost = -boost; // pPb MC has proton in + direction, pPb data has it in minus, and both are reversed for the 'Pbp definition'

	// Start loop over events
	double nev = (double)nevents;
	for(int i = 0; i < nevents; i++){

		hlt_tree->GetEntry(i);
		if(i != 0 && (i % 10000) == 0){double alpha = (double)i; cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;} // % processed
		if(do_quicktest) if(i != 0 && i % 20000 == 0 ) break; // just for quick tests
		Nevents->Fill(0); // filled after each event cut

		// Booleans to remove events which does not pass the Aj or Xj selection
		bool pass_Aj_or_Xj_reco_cut = true;
		bool pass_Aj_or_Xj_gen_cut = true;
		if(do_Xj_or_Ajcut){pass_Aj_or_Xj_reco_cut = false; pass_Aj_or_Xj_gen_cut = false;}

		// Apply trigger
		if(dojettrigger){if(jet_trigger_bit != 1) continue;}
		Nevents->Fill(1);

		// Apply event filters
		for(int ii = 0; ii < event_filter_bool.size(); ii++) if(event_filter_bool[ii] != 1) continue;
		Nevents->Fill(2);

		// Vectors used for objects
		// reco jets and tracks
		std::vector<TVector3> tracks_reco;
		std::vector<int> sube_tracks_reco;
		std::vector<TVector3> jets_reco;
		std::vector<TVector3> lead_jets_reco;
		std::vector<TVector3> subl_jets_reco;
		std::vector<double> track_w_reco;
		std::vector<double> jet_w_reco;
		std::vector<double> lead_jet_w_reco;
		std::vector<double> subl_jet_w_reco;
		// ref jets --> gen matched quantities
		std::vector<TVector3> jets_ref;
		std::vector<TVector3> lead_jets_ref;
		std::vector<TVector3> subl_jets_ref;
		std::vector<double> jet_w_ref;
		std::vector<double> lead_jet_w_ref;
		std::vector<double> subl_jet_w_ref;
		// gen jets and tracks
		std::vector<TVector3> tracks_gen;
		std::vector<int> sube_tracks_gen;
		std::vector<TVector3> jets_gen;
		std::vector<TVector3> lead_jets_gen;
		std::vector<TVector3> subl_jets_gen;
		std::vector<double> track_w_gen;
		std::vector<double> jet_w_gen;
		std::vector<double> lead_jet_w_gen;
		std::vector<double> subl_jet_w_gen;

		//apply event selections
		//Vz
		if(colliding_system == "pPb" && is_pgoing && invert_pgoing) { vertexz = -vertexz; } // in pPb we should invert all z quantities
		if(vertexz < vz_cut_min || vertexz > vz_cut_max) continue;
		Nevents->Fill(3);

		//pthat (MC only)
		if(do_pthatcut){if(pthat <= pthatmin || pthat > pthatmax) continue;} //pthat ranges
		if(is_MC){if(gen_jtpt[0] > 3.*pthat || rawpt[0] > 3.*pthat) continue;} //safety to remove some high-pT ejts from low pthat samples 
		Nevents->Fill(4);

		//multiplicity or centrality
		int trksize = (int)ntrk;
		int mult; // centrality or Ntroff
		if(use_centrality){mult = hiBin;}else{mult = get_Ntrkoff(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);}
		int recomult = get_simple_mult_reco(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge); // reco multiplicity (only pt and eta cuts)
		int genmult;
		if(is_MC) genmult = get_simple_mult_gen(colliding_system, sNN_energy_GeV, year_of_datataking, (int) gen_trkpt->size(), gen_trketa, gen_trkpt, gen_trkchg); // gen multiplicity (only pt and eta cuts)
		if(mult < multiplicity_centrality_bins[0] || mult > multiplicity_centrality_bins[multiplicity_centrality_bins.size()-1])continue; //centrality of multiplicity range
		//int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float) mult);
		double multcentbin = (double) mult;
		Nevents->Fill(5);
		
		// invert the HF/ZDC sides for pPb
		if(colliding_system=="pPb" && is_pgoing && invert_pgoing){
			float hfplus_temp = hfplus;
			float hfminus_temp = hfminus;
			float hfplusE4_temp = hfplusEta4;
			float hfminusE4_temp = hfminusEta4;
			float zdcplus_temp = zdcplus;
			float zdcminus_temp = zdcminus;
			hfplus = hfminus_temp; hfminus = hfplus_temp;
			hfplusEta4 = hfminusE4_temp; hfminusEta4 = hfplusE4_temp;
			zdcplus = zdcminus_temp; zdcminus = zdcplus_temp;
		} 

		//  add extra dependency
		float extra_variable = hfminusEta4;
		//int extrabin = (int) find_my_bin(extra_bins, (double) extra_variable);
		double extrabin = (double) extra_variable;

		// event weight(s), this must be applied in all histograms
		double event_weight = get_event_weight(nevents,is_MC, use_centrality, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, vertexz, mult, weight, pthat, extra_variable); // get the event weight

		// Fill vertex, pthat and multiplicity/centrality histograms
		multiplicity->Fill(mult);
		multiplicity_weighted->Fill(mult, event_weight);	
		reco_mult->Fill(recomult);
		reco_mult_weighted->Fill(recomult, event_weight);
		gen_mult->Fill(genmult);
		gen_mult_weighted->Fill(genmult, event_weight);
	
		double x_vz[3]={(double) vertexz, (double) multcentbin, (double) extrabin}; vzhist->Fill(x_vz); vzhist_weighted->Fill(x_vz,event_weight);
		double x_vxy[4]={(double) vertexx, (double) vertexy, (double) multcentbin, (double) extrabin}; vxyhist->Fill(x_vxy); vxyhist_weighted->Fill(x_vxy,event_weight);
		double x_pthat[3]={(double) pthat, (double) multcentbin, (double) extrabin}; pthathist->Fill(x_pthat); pthathist_weighted->Fill(x_pthat,event_weight);

		if(colliding_system=="pPb" && year_of_datataking==2016){
			// HF and ZDC histograms (fill histograms)		
			double x3D_hiHF[3]={hfplus,hfminus,(double) mult}; hfhist->Fill(x3D_hiHF); hfhist_weighted->Fill(x3D_hiHF,event_weight);
			double x3D_hiHFEta4[3]={hfplusEta4,hfminusEta4,(double) mult}; hfhistEta4->Fill(x3D_hiHFEta4); hfhistEta4_weighted->Fill(x3D_hiHFEta4,event_weight);
			double x3D_hiZDC[3]={zdcplus,zdcminus,(double) mult}; zdchist->Fill(x3D_hiZDC); zdchist_weighted->Fill(x3D_hiZDC,event_weight);
			double x2D_hiHFSum[2]={hfplus+hfminus,(double) mult}; hfhistSum_weighted->Fill(x2D_hiHFSum,event_weight);
			double x2D_hiHFEta4Sum[2]={hfplusEta4+hfminusEta4,(double) mult}; hfhistEta4Sum_weighted->Fill(x2D_hiHFEta4Sum,event_weight);			
			
			if(do_flow){
				//event plane information and filling of histograms
				// Psi 2
				double mult2_plus = (double)EP_Mult2_plus;
				double mult2_minus = (double)EP_Mult2_minus;
				double q2vec_plus = (double) sqrt(EP_Qx2_plus*EP_Qx2_plus + EP_Qy2_plus*EP_Qy2_plus);
				double q2vec_minus = (double)sqrt(EP_Qx2_minus*EP_Qx2_minus + EP_Qy2_minus*EP_Qy2_minus);
				double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;
				double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
				double x_reco_plus_ep2_flat[5]={mult2_plus,q2vec_plus,Psi2_EP_flat_plus,(double) multcentbin,(double) extrabin}; 
				EP2_plus_flat->Fill(x_reco_plus_ep2_flat,event_weight);
				double x_reco_minus_ep2_flat[5]={mult2_minus,q2vec_minus,Psi2_EP_flat_minus,(double) multcentbin,(double) extrabin}; 
				EP2_minus_flat->Fill(x_reco_minus_ep2_flat,event_weight);

				// Psi 3
				double mult3_plus = (double)EP_Mult3_plus;
				double mult3_minus = (double)EP_Mult3_minus;
				double q3vec_plus = (double) sqrt(EP_Qx3_plus*EP_Qx3_plus + EP_Qy3_plus*EP_Qy3_plus);
				double q3vec_minus = (double)sqrt(EP_Qx3_minus*EP_Qx3_minus + EP_Qy3_minus*EP_Qy3_minus);
				double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;
				double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
				double x_reco_plus_ep3_flat[5]={mult3_plus,q3vec_plus,Psi3_EP_flat_plus,(double) multcentbin,(double) extrabin}; 
				EP3_plus_flat->Fill(x_reco_plus_ep3_flat,event_weight);
				double x_reco_minus_ep3_flat[5]={mult3_minus,q3vec_minus,Psi3_EP_flat_minus,(double) multcentbin,(double) extrabin}; 
				EP3_minus_flat->Fill(x_reco_minus_ep3_flat,event_weight);

				// Psi 4
				double mult4_plus = (double)EP_Mult4_plus;
				double mult4_minus = (double)EP_Mult4_minus;
				double q4vec_plus = (double) sqrt(EP_Qx4_plus*EP_Qx4_plus + EP_Qy4_plus*EP_Qy4_plus);
				double q4vec_minus = (double)sqrt(EP_Qx4_minus*EP_Qx4_minus + EP_Qy4_minus*EP_Qy4_minus);
				double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;
				double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
				double x_reco_plus_ep4_flat[5]={mult4_plus,q4vec_plus,Psi4_EP_flat_plus,(double) multcentbin,(double) extrabin}; 
				EP4_plus_flat->Fill(x_reco_plus_ep4_flat,event_weight);
				double x_reco_minus_ep4_flat[5]={mult4_minus,q4vec_minus,Psi4_EP_flat_minus,(double) multcentbin,(double) extrabin}; 
				EP4_minus_flat->Fill(x_reco_minus_ep4_flat,event_weight);
			}
		}

		// ------------------- Reconstruction level (Data and MC) ----------------------------
		// Start loop over reco tracks (trksize is number of reco tracks)
		for (int j = 0; j < trksize; j++){ 

			// Define track/particle kinematics
			float trk_pt = trkpt[j];
			float trk_eta = trketa[j];
			float trk_phi = trkphi[j];
			
			int NDF = (int) trkndof[j];
			int NLayer = (int) trknlayer[j];
			int NHits = (int) trknhits[j];
			int Npixelhit = (int) trkpixhits[j];

			// Apply track selection (see read_tree.h to see what each variable means)
			if(fabs(trk_eta) > trk_eta_cut) continue;
			if(trkpt[j] <= trk_pt_min_cut) continue;
			if(highpur[j] == false) continue;
			if(fabs(trkpterr[j]/trkpt[j]) >= trk_pt_resolution_cut) continue;
			if(fabs(trkdcaxy[j]/trkdcaxyerr[j]) >= trk_dca_xy_cut) continue;
			if(fabs(trkdcaz[j]/trkdcazerr[j]) >= trk_dca_z_cut) continue;
			double calomatching = ((pfEcal[j]+pfHcal[j])/cosh(trketa[j]))/trkpt[j];
			if(colliding_system == "PbPb" || colliding_system == "XeXe"){
				if((trkchi2[j]/NDF)/NLayer >= chi2_ndf_nlayer_cut) continue;
				if(NHits < nhits) continue;
				if(trkpt[j] > 20.0 && fabs(calomatching) <= calo_matching) continue;
			}
			if(colliding_system=="PbPb" && sNN_energy_GeV==5020 && year_of_datataking==2018){if(trkalgo[j] == 6 && trkmva[j] < 0.98) continue;}

			// Track efficiency correction
			double trk_weight = 1.0;
			trk_weight = trk_weight*getTrkCorrWeight(fileeff, use_centrality, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, mult, trk_pt, trk_eta, trk_phi);

			trk_eta = trk_eta + boost; // In pPb case, for the center-of-mass correction if needed
			if(colliding_system == "pPb" && is_pgoing && invert_pgoing) trk_eta = -trk_eta; // in case of merging pgoing and Pbgoing use this

			// Track QA histogram filling
			double x_reco_trk[5]={trk_pt,trk_eta,trk_phi,(double) multcentbin,(double) extrabin}; 
			hist_reco_trk->Fill(x_reco_trk);
			hist_reco_trk_corr->Fill(x_reco_trk,trk_weight);
			hist_reco_trk_weighted->Fill(x_reco_trk,trk_weight*event_weight);

			// Track vector filling
			TVector3 GoodTracks;
			GoodTracks.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
			tracks_reco.push_back(GoodTracks);
			sube_tracks_reco.push_back(0); // set == 0 because this is only valid for gen
			double trk_etamix_weight = get_trketamix_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, trk_eta, true); // weight to deal with Seagull (test)
			track_w_reco.push_back(trk_weight*trk_etamix_weight); // save weight to apply in the mixing

			//int trackbin = (int) find_my_bin(trk_pt_bins, (float) trk_pt);
			double trackbin = (double) trk_pt;
			if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
				//event plane information
				// Psi 2
				double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;
				double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
				double x_reco_plus_dphi2_flat[4]={deltaphi2PC(trk_phi, Psi2_EP_flat_plus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
				Dphi_EP2_flat_trk_plus->Fill(x_reco_plus_dphi2_flat,trk_weight*event_weight);
				double x_reco_minus_dphi2_flat[4]={deltaphi2PC(trk_phi, Psi2_EP_flat_minus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
				Dphi_EP2_flat_trk_minus->Fill(x_reco_minus_dphi2_flat,trk_weight*event_weight);
			
				// Psi 3
				double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;
				double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
				double x_reco_plus_dphi3_flat[4]={deltaphi2PC(trk_phi, Psi3_EP_flat_plus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
				Dphi_EP3_flat_trk_plus->Fill(x_reco_plus_dphi3_flat,trk_weight*event_weight);
				double x_reco_minus_dphi3_flat[4]={deltaphi2PC(trk_phi, Psi3_EP_flat_minus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
				Dphi_EP3_flat_trk_minus->Fill(x_reco_minus_dphi3_flat,trk_weight*event_weight);
				
				// Psi 4
				double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;
				double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
				double x_reco_plus_dphi4_flat[4]={deltaphi2PC(trk_phi, Psi4_EP_flat_plus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
				Dphi_EP4_flat_trk_plus->Fill(x_reco_plus_dphi4_flat,trk_weight*event_weight);
				double x_reco_minus_dphi4_flat[4]={deltaphi2PC(trk_phi, Psi4_EP_flat_minus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
				Dphi_EP4_flat_trk_minus->Fill(x_reco_minus_dphi4_flat,trk_weight*event_weight);
			}
		} // End loop over tracks

		// to find leading, subleading and third jets
		float leadrecojet_pt=-999, leadrecojet_eta=-999, leadrecojet_phi=-999, leadrecojet_mass=-999; // leading jet quantities
		float sublrecojet_pt=-999, sublrecojet_eta=-999, sublrecojet_phi=-999, sublrecojet_mass=-999; // subleading jet quantities
		float thirdrecojet_pt=-999, thirdrecojet_eta=-999, thirdrecojet_phi=-999, thirdrecojet_mass=-999; // third jet quantities
		float leadrefjet_pt=-999, leadrefjet_eta=-999, leadrefjet_phi=-999, leadrefjet_mass=-999;; // leading jet ref quantities
		float sublrefjet_pt=-999, sublrefjet_eta=-999, sublrefjet_phi=-999, sublrefjet_mass=-999;; // subleading jet ref quantities
		float thirdrefjet_pt=-999, thirdrefjet_eta=-999, thirdrefjet_phi=-999, thirdrefjet_mass=-999; // third jet quantities
		
		bool isjetincluded = false;
		int njets = 0;
		
		int jetsize = (int)nref; // number of jets in an event
		// Start loop over jets
		for (int j = 0; j < jetsize; j++){

			if(trackMax[j]/rawpt[j] < 0.01)continue; // Cut for jets with only very low pT particles
			if(trackMax[j]/rawpt[j] > 0.98)continue; // Cut for jets where all the pT is taken by one track
			if(jteta[j] < -5.0 || jteta[j] > 5.0) continue; // no accept jets with |eta| > 5
			if(trackMax[j] < trackmaxpt) continue; // Can be use to remove jets from low pT tracks

			// Define jet kinematics
			float jet_rawpt = rawpt[j];
			float jet_eta = jteta[j];
			float jet_phi = jtphi[j];
			float jet_mass = jtmass[j];

			// Apply JEC
			JEC.SetJetPT(jet_rawpt); 
			JEC.SetJetEta(jet_eta); 
			JEC.SetJetPhi(jet_phi);
			float jet_pt_corr = JEC.GetCorrectedPT();

			//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
			double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr); // Jet weight (specially for MC)
			if(do_jet_smearing) jet_weight = jet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr, do_jet_smearing, 0.663); // Jet smearing (For systematics)
			
			//leading and subleading
			find_leading_subleading_third(jet_pt_corr,jet_eta,jet_phi,jet_mass,leadrecojet_pt,leadrecojet_eta,leadrecojet_phi,leadrecojet_mass,sublrecojet_pt,sublrecojet_eta,sublrecojet_phi,sublrecojet_mass,thirdrecojet_pt,thirdrecojet_eta,thirdrecojet_phi,thirdrecojet_mass); // Find leading and subleading jets

			jet_eta = jet_eta + boost; // In pPb case, for the center-of-mass correction if needed
			if(colliding_system == "pPb" && is_pgoing && invert_pgoing)jet_eta = -jet_eta;

			// Fill reco jet QA histograms
			double x_reco_jet[5]={jet_rawpt,jet_eta,jet_phi,(double) multcentbin,(double) extrabin}; 
			hist_reco_jet->Fill(x_reco_jet);
			hist_reco_jet_weighted->Fill(x_reco_jet,event_weight*jet_weight);
			double x_reco_jet_corr[5]={jet_pt_corr,jet_eta,jet_phi,(double) multcentbin,(double) extrabin}; 
			hist_reco_jet_corr->Fill(x_reco_jet_corr);
			hist_reco_jet_corr_weighted->Fill(x_reco_jet_corr,event_weight*jet_weight);

			if((jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut) && (jet_eta > jet_eta_min_cut && jet_eta < jet_eta_max_cut)){ // Jet pT and eta cut		

				njets=njets+1;
				isjetincluded = true;
				double x_trkmax[3]={trackMax[j],(double) multcentbin,(double) extrabin}; 
				trackmaxptinjethisto->Fill(x_trkmax,event_weight*jet_weight);
				// Fill reco jet vectors
				TVector3 GoodJets;
				GoodJets.SetPtEtaPhi(jet_pt_corr, jet_eta, jet_phi);
				jets_reco.push_back(GoodJets);
				jet_w_reco.push_back(jet_weight);
					
				if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
					// Psi 2
					double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;
					double x_jetep2_plus[3]={deltaphi2PC(jet_phi, Psi2_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP2_inclusive_plus->Fill(x_jetep2_plus,event_weight*jet_weight);
					double Psi2_EP_flat_minus = (double) EP_Psi2_minus_flat;
					double x_jetep2_minus[3]={deltaphi2PC(jet_phi, Psi2_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP2_inclusive_minus->Fill(x_jetep2_minus,event_weight*jet_weight);
					// Psi 3
					double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;
					double x_jetep3_plus[3]={deltaphi2PC(jet_phi, Psi3_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP3_inclusive_plus->Fill(x_jetep3_plus,event_weight*jet_weight);
					double Psi3_EP_flat_minus = (double) EP_Psi3_minus_flat;
					double x_jetep3_minus[3]={deltaphi2PC(jet_phi, Psi3_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP3_inclusive_minus->Fill(x_jetep3_minus,event_weight*jet_weight);
					// Psi 4
					double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;
					double x_jetep4_plus[3]={deltaphi2PC(jet_phi, Psi4_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP4_inclusive_plus->Fill(x_jetep4_plus,event_weight*jet_weight);
					double Psi4_EP_flat_minus = (double) EP_Psi4_minus_flat;
					double x_jetep4_minus[3]={deltaphi2PC(jet_phi, Psi4_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP4_inclusive_minus->Fill(x_jetep4_minus,event_weight*jet_weight);
				}
				
			}

			if(is_MC){ 

				float ref_pt = refpt[j];
				float ref_eta = refeta[j];
				float ref_phi = refphi[j];
				float ref_mass = refmass[j];
				
				if(ref_pt <= 0) continue; // just remove non-matched
				if(jet_rawpt <= 0) continue;
				if(jet_pt_corr <= 0) continue;
				if(ref_eta < -5.0 || ref_eta > 5.0) continue; // max jet eta

				double refjet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, ref_pt); // Jet weight (specially for MC)

				find_leading_subleading_third(ref_pt,ref_eta,ref_phi,ref_mass,leadrefjet_pt,leadrefjet_eta,leadrefjet_phi,leadrefjet_mass,sublrefjet_pt,sublrefjet_eta,sublrefjet_phi,sublrefjet_mass,thirdrefjet_pt,thirdrefjet_eta,thirdrefjet_phi,thirdrefjet_mass); // Find leading and subleading ref jets

				ref_eta = ref_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				if(colliding_system == "pPb" && is_pgoing && invert_pgoing)ref_eta = -ref_eta;

				double JES_ratio_reco_vs_ref = jet_pt_corr/ref_pt;
				int refparton; if(fabs(refparton_flavor[j]) >= 1 && fabs(refparton_flavor[j]) <= 6){refparton = fabs(refparton_flavor[j]);}else if(fabs(refparton_flavor[j]) == 21){refparton = 7;}else{refparton = 0;}
				double x_JES_ratio_reco_vs_ref[6]={JES_ratio_reco_vs_ref,ref_pt,ref_eta,(double)refparton,(double)multcentbin,(double) extrabin}; 
				hist_jes_reco_weighted->Fill(x_JES_ratio_reco_vs_ref,event_weight*refjet_weight*jet_weight);
				int refpartonfromB; if(fabs(refparton_flavorForB[j]) >= 1 && fabs(refparton_flavorForB[j]) <= 6){refpartonfromB = fabs(refparton_flavorForB[j]);}else if(fabs(refparton_flavorForB[j]) == 21){refpartonfromB = 7;}else{refpartonfromB = 0;}
				double x_JES_ratio_reco_vs_reffromB[6]={JES_ratio_reco_vs_ref,ref_pt,ref_eta,(double)refpartonfromB,(double)multcentbin,(double) extrabin}; 
				hist_jes_reco_fromB_weighted->Fill(x_JES_ratio_reco_vs_reffromB,event_weight*refjet_weight*jet_weight);
				
			} // End if over ref (MC)

		} // End loop over jets

		if(isjetincluded){
			multiplicity_withonejet_weighted->Fill(mult,event_weight);
			reco_mult_withonejet_weighted->Fill(recomult, event_weight);
			vzhist_jet_weighted->Fill(x_vz, event_weight);
			if(colliding_system=="pPb" && year_of_datataking==2016){
				double x3D_hiHF_onejet[3]={hfplus,hfminus,(double) mult}; hfhist_onejet_weighted->Fill(x3D_hiHF_onejet,event_weight);
				double x3D_hiHFEta4_onejet[3]={hfplusEta4,hfminusEta4,(double) mult}; hfhistEta4_onejet_weighted->Fill(x3D_hiHFEta4_onejet,event_weight);
				double x3D_hiZDC_onejet[3]={zdcplus,zdcminus,(double) mult}; zdchist_onejet_weighted->Fill(x3D_hiZDC_onejet,event_weight);
				double x2D_hiHFSum_onejet[2]={hfplus+hfminus,(double) mult}; hfhistSum_onejet_weighted->Fill(x2D_hiHFSum_onejet,event_weight);
				double x2D_hiHFEta4Sum_onejet[2]={hfplusEta4+hfminusEta4,(double) mult}; hfhistEta4Sum_onejet_weighted->Fill(x2D_hiHFEta4Sum_onejet,event_weight);			
			}
		}

		double xnjets[3]={(double) njets, (double) multcentbin, (double) extrabin}; NJets->Fill(xnjets,event_weight);
		bool isdijet = false;

		bool removethirdjet = false;
		if(do_thirdjet_removal){if(thirdrecojet_pt > 0.5*sublrecojet_pt) removethirdjet = true;}

		//dijets
		if(jetsize > 1 && !removethirdjet){

			Nevents->Fill(6);
			double ljet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrecojet_pt);  // Jet weight (specially for MC)
			//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
			if(do_jet_smearing) ljet_weight = ljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrecojet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

			double sljet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrecojet_pt);  // Jet weight (specially for MC)
			//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
			if(do_jet_smearing) sljet_weight = sljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrecojet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

			//leading/subleading pT cuts
			if(leadrecojet_pt > leading_pT_min && sublrecojet_pt > subleading_pT_min){
				
				Nevents->Fill(7);
				double leadrecojet_eta_lab = leadrecojet_eta;  // before boost for eta dijet
				double sublrecojet_eta_lab = sublrecojet_eta;  // before boost for eta dijet
				leadrecojet_eta = leadrecojet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				sublrecojet_eta = sublrecojet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				
				if(colliding_system == "pPb" && is_pgoing && invert_pgoing){
					leadrecojet_eta = -leadrecojet_eta; 
					sublrecojet_eta = -sublrecojet_eta;
					leadrecojet_eta_lab = -leadrecojet_eta_lab; 
					sublrecojet_eta_lab = -sublrecojet_eta_lab;
				}
				
				Nevents->Fill(8);
				
				// Fill leading/subleading jet quenching quantities
				double delta_phi_reco = fabs(deltaphi(leadrecojet_phi, sublrecojet_phi));
				double Aj_reco = asymmetry(leadrecojet_pt,sublrecojet_pt);
				double Xj_reco = xjvar(leadrecojet_pt,sublrecojet_pt);
				float ptdijet = 0.5*(leadrecojet_pt + sublrecojet_pt);
				//int ptdijetbin = (int) find_my_bin(pt_ave_bins, (float) ptdijet);
				double ptdijetbin = (double) ptdijet;
				double x_reco[6]={Xj_reco,Aj_reco,delta_phi_reco,(double)multcentbin,(double)ptdijetbin,(double)extrabin}; 
				// combinations of midrapidity, forward and backward
				bool leadmidrap = (leadrecojet_eta > jet_eta_min_cut && leadrecojet_eta < jet_eta_max_cut);
				bool sublmidrap = (sublrecojet_eta > jet_eta_min_cut && sublrecojet_eta < jet_eta_max_cut);
				bool leadfwdrap = (leadrecojet_eta > jet_fwd_eta_min_cut && leadrecojet_eta < jet_fwd_eta_max_cut);
				bool sublfwdrap = (sublrecojet_eta > jet_fwd_eta_min_cut && sublrecojet_eta < jet_fwd_eta_max_cut);
				bool leadbkwrap = (leadrecojet_eta > jet_bkw_eta_min_cut && leadrecojet_eta < jet_bkw_eta_max_cut);
				bool sublbkwrap = (sublrecojet_eta > jet_bkw_eta_min_cut && sublrecojet_eta < jet_bkw_eta_max_cut);
				if(do_dijetstudies){
					if(leadmidrap && sublmidrap) hist_reco_lead_reco_subl_quench_mid_mid->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadmidrap && sublfwdrap) hist_reco_lead_reco_subl_quench_mid_fwd->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadmidrap && sublbkwrap) hist_reco_lead_reco_subl_quench_mid_bkw->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadfwdrap && sublmidrap) hist_reco_lead_reco_subl_quench_fwd_mid->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadbkwrap && sublmidrap) hist_reco_lead_reco_subl_quench_bkw_mid->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadfwdrap && sublfwdrap) hist_reco_lead_reco_subl_quench_fwd_fwd->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadfwdrap && sublbkwrap) hist_reco_lead_reco_subl_quench_fwd_bkw->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadbkwrap && sublfwdrap) hist_reco_lead_reco_subl_quench_bkw_fwd->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
					if(leadbkwrap && sublbkwrap) hist_reco_lead_reco_subl_quench_bkw_bkw->Fill(x_reco,event_weight*ljet_weight*sljet_weight);
				}
				
				if(colliding_system=="pPb" && year_of_datataking==2016 && do_dijetstudies){

					double etadijet = 0.5*(leadrecojet_eta_lab + sublrecojet_eta_lab);
					double etadiff = 0.5*deltaeta(leadrecojet_eta_lab,sublrecojet_eta_lab);
					double xp_reco = 2.0*(ptdijet*exp(etadijet)*cosh(etadiff))/sqrts;
					double xPb_reco = 2.0*(ptdijet*exp(-etadijet)*cosh(etadiff))/sqrts;
					double m12_reco = 2.0*ptdijet*cosh(etadiff); // for future, if needed 		
					double x_dijet_reco[11] = {etadijet, etadiff, Xj_reco, Aj_reco, delta_phi_reco, xp_reco, xPb_reco, m12_reco, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

					double etadijet_CM = 0.5*(leadrecojet_eta + sublrecojet_eta);
					double etadiff_CM = 0.5*deltaeta(leadrecojet_eta,sublrecojet_eta);
					double xp_reco_CM = 2.0*(ptdijet*exp(etadijet_CM)*cosh(etadiff_CM))/sqrts;
					double xPb_reco_CM = 2.0*(ptdijet*exp(-etadijet_CM)*cosh(etadiff_CM))/sqrts;
					double m12_reco_CM = 2.0*ptdijet*cosh(etadiff_CM); // for future, if needed 		
					double x_dijet_reco_CM[11] = {etadijet_CM, etadiff_CM, Xj_reco, Aj_reco, delta_phi_reco, xp_reco_CM, xPb_reco_CM, m12_reco_CM, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

					ROOT::Math::PtEtaPhiMVector Lead_reco_jet(leadrecojet_pt,leadrecojet_eta,leadrecojet_phi,leadrecojet_mass);
					ROOT::Math::PtEtaPhiMVector Subl_reco_jet(sublrecojet_pt,sublrecojet_eta,sublrecojet_phi,sublrecojet_mass);			
					double etadijet_CM_y = 0.5*(Lead_reco_jet.Rapidity() + Subl_reco_jet.Rapidity());
					double etadiff_CM_y = 0.5*deltaeta(Lead_reco_jet.Rapidity(),Subl_reco_jet.Rapidity());
					double xp_reco_CM_y = 2.0*(ptdijet*exp(etadijet_CM_y)*cosh(etadiff_CM_y))/sqrts;
					double xPb_reco_CM_y = 2.0*(ptdijet*exp(-etadijet_CM_y)*cosh(etadiff_CM_y))/sqrts;
					double m12_reco_CM_y = 2.0*ptdijet*cosh(etadiff_CM_y); // for future, if needed 					
					double x_dijet_reco_CM_y[11] = {etadijet_CM_y, etadiff_CM_y, Xj_reco, Aj_reco, delta_phi_reco, xp_reco_CM_y, xPb_reco_CM_y, m12_reco_CM_y, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

					if(fabs(leadrecojet_eta_lab) < dijetetamax && fabs(sublrecojet_eta_lab) < dijetetamax){
						hist_etaDijet_reco->Fill(x_dijet_reco,event_weight*ljet_weight*sljet_weight); 
						hist_etaDijet_CM_reco->Fill(x_dijet_reco_CM,event_weight*ljet_weight*sljet_weight);
						hist_yDijet_CM_reco->Fill(x_dijet_reco_CM_y,event_weight*ljet_weight*sljet_weight);
					}
					
				}


				// leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
				if(delta_phi_reco > leading_subleading_deltaphi_min){

					if((Xj_reco >= xjmin && Xj_reco <= xjmax) && (Aj_reco >= Ajmin && Aj_reco <= Ajmax)){

						Nevents->Fill(9);
						pass_Aj_or_Xj_reco_cut = true; // if we apply Xj or Aj cuts
						isdijet = true;

						// Fill leading and subleading jet QA histograms
						double x_lead[5]={leadrecojet_pt,leadrecojet_eta,leadrecojet_phi,(double) multcentbin,(double)extrabin}; 
						hist_reco_leadjet->Fill(x_lead);
						hist_reco_leadjet_weighted->Fill(x_lead,event_weight*ljet_weight);
						
						double x_sublead[5]={sublrecojet_pt,sublrecojet_eta,sublrecojet_phi,(double) multcentbin,(double)extrabin}; 
						hist_reco_subljet->Fill(x_sublead);
						hist_reco_subljet_weighted->Fill(x_sublead,event_weight*sljet_weight);

						bool usethisjets = false;
						if(fwdbkw_jettrk_option == "mid_mid"){ if(leadmidrap && sublmidrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "mid_fwd"){ if(leadmidrap && sublfwdrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "mid_bkw"){ if(leadmidrap && sublbkwrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "fwd_mid"){ if(leadfwdrap && sublmidrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "bkw_mid"){ if(leadbkwrap && sublmidrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "fwd_fwd"){ if(leadfwdrap && sublfwdrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "fwd_bkw"){ if(leadfwdrap && sublbkwrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "bkw_fwd"){ if(leadbkwrap && sublfwdrap) usethisjets = true; }
						if(fwdbkw_jettrk_option == "bkw_bkw"){ if(leadbkwrap && sublbkwrap) usethisjets = true; }

						// Fill leading and subleading jet vectors
						TVector3 GoodLeadingJets_reco;
						TVector3 GoodSubLeadingJets_reco;
						if(usethisjets){
							GoodLeadingJets_reco.SetPtEtaPhi(leadrecojet_pt, leadrecojet_eta, leadrecojet_phi);
							lead_jets_reco.push_back(GoodLeadingJets_reco);
							lead_jet_w_reco.push_back(ljet_weight);
							GoodSubLeadingJets_reco.SetPtEtaPhi(sublrecojet_pt, sublrecojet_eta, sublrecojet_phi);
							subl_jets_reco.push_back(GoodSubLeadingJets_reco);
							subl_jet_w_reco.push_back(sljet_weight);
						}

						if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
							// Psi 2
							double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;
							double Psi2_EP_flat_minus = (double) EP_Psi2_minus_flat;
							double x_ljetep2_plus[3]={deltaphi2PC(leadrecojet_phi, Psi2_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP2_leading_plus->Fill(x_ljetep2_plus,event_weight*ljet_weight);
							double x_ljetep2_minus[3]={deltaphi2PC(leadrecojet_phi, Psi2_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP2_leading_minus->Fill(x_ljetep2_minus,event_weight*ljet_weight);
							double x_sljetep2_plus[3]={deltaphi2PC(sublrecojet_phi, Psi2_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP2_subleading_plus->Fill(x_sljetep2_plus,event_weight*sljet_weight);
							double x_sljetep2_minus[3]={deltaphi2PC(sublrecojet_phi, Psi2_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP2_subleading_minus->Fill(x_sljetep2_minus,event_weight*sljet_weight);
							// Psi 3
							double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;
							double Psi3_EP_flat_minus = (double) EP_Psi3_minus_flat;
							double x_ljetep3_plus[3]={deltaphi2PC(leadrecojet_phi, Psi3_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP3_leading_plus->Fill(x_ljetep3_plus,event_weight*ljet_weight);
							double x_ljetep3_minus[3]={deltaphi2PC(leadrecojet_phi, Psi3_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP3_leading_minus->Fill(x_ljetep3_minus,event_weight*ljet_weight);
							double x_sljetep3_plus[3]={deltaphi2PC(sublrecojet_phi, Psi3_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP3_subleading_plus->Fill(x_sljetep3_plus,event_weight*sljet_weight);
							double x_sljetep3_minus[3]={deltaphi2PC(sublrecojet_phi, Psi3_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP3_subleading_minus->Fill(x_sljetep3_minus,event_weight*sljet_weight);
							// Psi 4
							double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;
							double Psi4_EP_flat_minus = (double) EP_Psi4_minus_flat;
							double x_ljetep4_plus[3]={deltaphi2PC(leadrecojet_phi, Psi4_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP4_leading_plus->Fill(x_ljetep4_plus,event_weight*ljet_weight);
							double x_ljetep4_minus[3]={deltaphi2PC(leadrecojet_phi, Psi4_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP4_leading_minus->Fill(x_ljetep4_minus,event_weight*ljet_weight);
							double x_sljetep4_plus[3]={deltaphi2PC(sublrecojet_phi, Psi4_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP4_subleading_plus->Fill(x_sljetep4_plus,event_weight*sljet_weight);
							double x_sljetep4_minus[3]={deltaphi2PC(sublrecojet_phi, Psi4_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_flat_EP4_subleading_minus->Fill(x_sljetep4_minus,event_weight*sljet_weight);
						}
					}
				}
			}
		}

		if(isdijet){
			multiplicity_withdijets_weighted->Fill(mult,event_weight);
			reco_mult_withdijets_weighted->Fill(recomult, event_weight);
			vzhist_dijet_weighted->Fill(x_vz, event_weight);
			if(colliding_system=="pPb" && year_of_datataking==2016){
				double x3D_hiHF_dijet[3]={hfplus,hfminus,(double) mult}; hfhist_dijet_weighted->Fill(x3D_hiHF_dijet,event_weight);
				double x3D_hiHFEta4_dijet[3]={hfplusEta4,hfminusEta4,(double) mult}; hfhistEta4_dijet_weighted->Fill(x3D_hiHFEta4_dijet,event_weight);
				double x3D_hiZDC_dijet[3]={zdcplus,zdcminus,(double) mult}; zdchist_dijet_weighted->Fill(x3D_hiZDC_dijet,event_weight);
				double x2D_hiHFSum_dijet[2]={hfplus+hfminus,(double) mult}; hfhist_dijet_weighted->Fill(x2D_hiHFSum_dijet,event_weight);
				double x2D_hiHFEta4Sum_dijet[2]={hfplusEta4+hfminusEta4,(double) mult}; hfhistEta4_dijet_weighted->Fill(x2D_hiHFEta4Sum_dijet,event_weight);
			}
		}

		bool removethirdjet_ref = false;
		if(do_thirdjet_removal){if(thirdrefjet_pt > 0.5*sublrefjet_pt) removethirdjet_ref = true;}
		
		if(jetsize > 1 && !removethirdjet_ref){
			//leading/subleading pT cuts
			if(is_MC && leadrefjet_pt > leading_pT_min && sublrefjet_pt > subleading_pT_min){

				double lrefjet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrefjet_pt);  // Jet weight (specially for MC)
				double slrefjet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrefjet_pt);  // Jet weight (specially for MC)
			
				double leadrefjet_eta_lab = leadrefjet_eta; // before boost for eta dijet
				double sublrefjet_eta_lab = sublrefjet_eta; // before boost for eta dijet
				leadrefjet_eta = leadrefjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				sublrefjet_eta = sublrefjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				if(colliding_system == "pPb" && is_pgoing && invert_pgoing){
					leadrefjet_eta = -leadrefjet_eta; 
					sublrefjet_eta = -sublrefjet_eta;
					leadrefjet_eta_lab = -leadrefjet_eta_lab; 
					sublrefjet_eta_lab = -sublrefjet_eta_lab;
				}

				// Fill leading/subleading jet quenching quantities
				double delta_phi_ref = fabs(deltaphi(leadrefjet_phi, sublrefjet_phi));
				double Aj_ref = asymmetry(leadrefjet_pt,sublrefjet_pt);
				double Xj_ref = xjvar(leadrefjet_pt,sublrefjet_pt);
				float ptdijet = 0.5*(leadrefjet_pt + sublrefjet_pt);
				//int ptdijetbin = (int) find_my_bin(pt_ave_bins, (float) ptdijet);
				double ptdijetbin = (double) ptdijet;
				double x_ref[6]={Xj_ref,Aj_ref,delta_phi_ref,(double)multcentbin,(double)ptdijetbin,(double)extrabin};
				// combinations of midrapidity, forward and backward
				bool leadmidrap = (leadrefjet_eta > jet_eta_min_cut && leadrefjet_eta < jet_eta_max_cut);
				bool sublmidrap = (sublrefjet_eta > jet_eta_min_cut && sublrefjet_eta < jet_eta_max_cut);
				bool leadfwdrap = (leadrefjet_eta > jet_fwd_eta_min_cut && leadrefjet_eta < jet_fwd_eta_max_cut);
				bool sublfwdrap = (sublrefjet_eta > jet_fwd_eta_min_cut && sublrefjet_eta < jet_fwd_eta_max_cut);
				bool leadbkwrap = (leadrefjet_eta > jet_bkw_eta_min_cut && leadrefjet_eta < jet_bkw_eta_max_cut);
				bool sublbkwrap = (sublrefjet_eta > jet_bkw_eta_min_cut && sublrefjet_eta < jet_bkw_eta_max_cut);
				if(do_dijetstudies){
					if(leadmidrap && sublmidrap) hist_ref_lead_ref_subl_quench_mid_mid->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadmidrap && sublfwdrap) hist_ref_lead_ref_subl_quench_mid_fwd->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadmidrap && sublbkwrap) hist_ref_lead_ref_subl_quench_mid_bkw->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadfwdrap && sublmidrap) hist_ref_lead_ref_subl_quench_fwd_mid->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadbkwrap && sublmidrap) hist_ref_lead_ref_subl_quench_bkw_mid->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadfwdrap && sublfwdrap) hist_ref_lead_ref_subl_quench_fwd_fwd->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadfwdrap && sublbkwrap) hist_ref_lead_ref_subl_quench_fwd_bkw->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadbkwrap && sublfwdrap) hist_ref_lead_ref_subl_quench_bkw_fwd->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
					if(leadbkwrap && sublbkwrap) hist_ref_lead_ref_subl_quench_bkw_bkw->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);
				}

				if(colliding_system=="pPb" && year_of_datataking==2016 && do_dijetstudies){

					double etadijet = 0.5*(leadrefjet_eta_lab + sublrefjet_eta_lab);
					double etadiff = 0.5*deltaeta(leadrefjet_eta_lab,sublrefjet_eta_lab);
					double xp_ref = 2.0*(ptdijet*exp(etadijet)*cosh(etadiff))/sqrts;
					double xPb_ref = 2.0*(ptdijet*exp(-etadijet)*cosh(etadiff))/sqrts;
					double m12_ref = 2.0*ptdijet*cosh(etadiff); // for future, if needed 	
					double x_dijet_ref[11] = {etadijet, etadiff, Xj_ref, Aj_ref, delta_phi_ref, xp_ref, xPb_ref, m12_ref, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

					double etadijet_CM = 0.5*(leadrefjet_eta + sublrefjet_eta);
					double etadiff_CM = 0.5*deltaeta(leadrefjet_eta,sublrefjet_eta);
					double xp_ref_CM = 2.0*(ptdijet*exp(etadijet_CM)*cosh(etadiff_CM))/sqrts;
					double xPb_ref_CM = 2.0*(ptdijet*exp(-etadijet_CM)*cosh(etadiff_CM))/sqrts;
					double m12_ref_CM = 2.0*ptdijet*cosh(etadiff_CM); // for future, if needed 	
					double x_dijet_ref_CM[11] = {etadijet_CM, etadiff_CM, Xj_ref, Aj_ref, delta_phi_ref, xp_ref_CM, xPb_ref_CM, m12_ref_CM, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

					ROOT::Math::PtEtaPhiMVector Lead_ref_jet(leadrefjet_pt,leadrefjet_eta,leadrefjet_phi,leadrefjet_mass);
					ROOT::Math::PtEtaPhiMVector Subl_ref_jet(sublrefjet_pt,sublrefjet_eta,sublrefjet_phi,sublrefjet_mass);			
					double etadijet_CM_y = 0.5*(Lead_ref_jet.Rapidity() + Subl_ref_jet.Rapidity());
					double etadiff_CM_y = 0.5*deltaeta(Lead_ref_jet.Rapidity(),Subl_ref_jet.Rapidity());
					double xp_ref_CM_y = 2.0*(ptdijet*exp(etadijet_CM_y)*cosh(etadiff_CM_y))/sqrts;
					double xPb_ref_CM_y = 2.0*(ptdijet*exp(-etadijet_CM_y)*cosh(etadiff_CM_y))/sqrts;
					double m12_ref_CM_y = 2.0*ptdijet*cosh(etadiff_CM_y); // for future, if needed 					
					double x_dijet_ref_CM_y[11] = {etadijet_CM_y, etadiff_CM_y, Xj_ref, Aj_ref, delta_phi_ref, xp_ref_CM_y, xPb_ref_CM_y, m12_ref_CM_y, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

					if(fabs(leadrefjet_eta_lab) < dijetetamax && fabs(sublrefjet_eta_lab) < dijetetamax){
						hist_etaDijet_ref->Fill(x_dijet_ref,event_weight*lrefjet_weight*slrefjet_weight); 
						hist_etaDijet_CM_ref->Fill(x_dijet_ref_CM,event_weight*lrefjet_weight*slrefjet_weight);
						hist_yDijet_CM_ref->Fill(x_dijet_ref_CM_y,event_weight*lrefjet_weight*slrefjet_weight);
					}
				}
			}
		}
		
		// Measure correlations and filling mixing vectors
		// Reco-Reco
		// Inclusive jets
		if(do_inclusejettrack_correlation && pass_Aj_or_Xj_reco_cut && !removethirdjet_ref){
			correlation(jets_reco, jet_w_reco, tracks_reco, track_w_reco, hist_correlation_signal_jet_reco_track_reco, hist_jet_from_reco_reco_sig, hist_trk_from_reco_reco_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_jet_reco_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_jet_reco_track_reco,JetR,hist_injet_reco_track_reco,do_flow); // calculate correlations
			fillvectors(similar_events, Nev_recoreco, jets_reco, jet_w_reco, tracks_reco, track_w_reco, mult, vertexz, event_weight, extra_variable, ev_jet_vector_reco_reco, jet_weights_reco_reco, ev_track_vector_reco_reco, trk_weights_reco_reco, multvec_reco_reco, vzvec_reco_reco, weights_reco_reco,extravec_reco_reco);	// for mixing --> store vectors for mixing
		}

		// Leading/SubLeading jets
		if(do_leading_subleading_jettrack_correlation && pass_Aj_or_Xj_reco_cut && !removethirdjet_ref){
			correlation(lead_jets_reco, lead_jet_w_reco, tracks_reco, track_w_reco, hist_correlation_signal_lead_jet_reco_track_reco, hist_lead_jet_from_reco_reco_sig, hist_LJ_trk_from_reco_reco_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_reco_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_lead_jet_reco_track_reco,JetR,hist_inLeadjet_reco_track_reco,do_flow); // calculate correlations
			correlation(subl_jets_reco, subl_jet_w_reco, tracks_reco, track_w_reco, hist_correlation_signal_subl_jet_reco_track_reco, hist_subl_jet_from_reco_reco_sig, hist_SLJ_trk_from_reco_reco_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_reco_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_subl_jet_reco_track_reco,JetR,hist_inSubljet_reco_track_reco,do_flow); // calculate correlations
			fillvectors(similar_events, Nev_recoreco_lead, lead_jets_reco, lead_jet_w_reco, tracks_reco, track_w_reco, mult, vertexz, event_weight, extra_variable, ev_jet_vector_leadjet_reco_reco, jet_weights_leadjet_reco_reco, ev_track_vector_leadjet_reco_reco, trk_weights_leadjet_reco_reco, multvec_leadjet_reco_reco, vzvec_leadjet_reco_reco, weights_leadjet_reco_reco, extravec_leadjet_reco_reco);	// for mixing --> store vectors for mixing
			fillvectors(similar_events, Nev_recoreco_subl, subl_jets_reco, subl_jet_w_reco, tracks_reco, track_w_reco, mult, vertexz, event_weight, extra_variable, ev_jet_vector_subleadjet_reco_reco, jet_weights_subleadjet_reco_reco, ev_track_vector_subleadjet_reco_reco, trk_weights_subleadjet_reco_reco, multvec_subleadjet_reco_reco, vzvec_subleadjet_reco_reco, weights_subleadjet_reco_reco, extravec_subleadjet_reco_reco); // for mixing --> store vectors for mixing
		}

		if(do_flow){twoparticlecorrelation(tracks_reco, track_w_reco, hist_reco_reco_2pcorrelation_signal, event_weight, mult, extra_variable, sube_tracks_reco, hist_reco_reco_2pcorrelation_signal_subg0, hist_reco_reco_2pcorrelation_signal_subcross);} // calculate 2 particle correlations

		// Generator level --> MC only
		int gentrksize; // number of gen tracks/particles
		if(is_MC) gentrksize = (int)gen_trkpt->size(); 
		int gen_jetsize; // number of gen jets
		if(is_MC) gen_jetsize = (int)ngen;

		if(is_MC){
			// Start loop over gen particles
			for(int j = 0; j < gentrksize; j++){ 

				// Define track/particle kinematics
				float gtrk_pt = gen_trkpt->at(j);
				float gtrk_eta = gen_trketa->at(j);
				float gtrk_phi = gen_trkphi->at(j);

				// Kinematic and charge cuts
				if(fabs(gtrk_eta) > trk_eta_cut) continue;
				if(gen_trkpt->at(j) <= trk_pt_min_cut)continue;
				if(!do_pid) if(gen_trkchg->at(j) == 0) continue;
				if(do_pid){if(fabs(gen_trkpid->at(j)) != particlepid) continue;}

				gtrk_eta = gtrk_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				if(colliding_system == "pPb" && is_pgoing && invert_pgoing){gtrk_eta = -gtrk_eta;}

				// Track/particle QA histogram filling
				double x_gen_trk[5]={gtrk_pt,gtrk_eta,gtrk_phi,(double) multcentbin,(double) extrabin}; 
				hist_gen_trk->Fill(x_gen_trk);
				hist_gen_trk_weighted->Fill(x_gen_trk,event_weight);

				// Track/particle vector filling
				TVector3 GoodTracks_gen;
				GoodTracks_gen.SetPtEtaPhi(gtrk_pt, gtrk_eta, gtrk_phi);
				tracks_gen.push_back(GoodTracks_gen);
				sube_tracks_gen.push_back(gen_trksube->at(j)); // get sube from the tree
				double trk_etamix_weight = get_trketamix_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gtrk_eta, false); // weight to deal with Seagull (test)
				track_w_gen.push_back(trk_etamix_weight); // save weight to apply in the mixing
				
				//int trackbin = (int) find_my_bin(trk_pt_bins, (float) gtrk_pt);
				double trackbin = (double) gtrk_pt;
				if(do_flow){
					//event plane information
					// Psi 2
					double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;
					double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
					double x_gen_plus_dphi2_flat[4]={deltaphi2PC(gtrk_phi, Psi2_EP_flat_plus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
					Dphi_GEN_EP2_flat_trk_plus->Fill(x_gen_plus_dphi2_flat,event_weight);
					double x_gen_minus_dphi2_flat[4]={deltaphi2PC(gtrk_phi, Psi2_EP_flat_minus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
					Dphi_GEN_EP2_flat_trk_minus->Fill(x_gen_minus_dphi2_flat,event_weight);

					// Psi 3
					double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;
					double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
					double x_gen_plus_dphi3_flat[4]={deltaphi2PC(gtrk_phi, Psi3_EP_flat_plus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
					Dphi_GEN_EP3_flat_trk_plus->Fill(x_gen_plus_dphi3_flat,event_weight);
					double x_gen_minus_dphi3_flat[4]={deltaphi2PC(gtrk_phi, Psi3_EP_flat_minus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
					Dphi_GEN_EP3_flat_trk_minus->Fill(x_gen_minus_dphi3_flat,event_weight);

					// Psi 4
					double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;
					double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
					double x_gen_plus_dphi4_flat[4]={deltaphi2PC(gtrk_phi, Psi4_EP_flat_plus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
					Dphi_GEN_EP4_flat_trk_plus->Fill(x_gen_plus_dphi4_flat,event_weight);
					double x_gen_minus_dphi4_flat[4]={deltaphi2PC(gtrk_phi, Psi4_EP_flat_minus), (double) trackbin, (double) multcentbin,(double) extrabin}; 
					Dphi_GEN_EP4_flat_trk_minus->Fill(x_gen_minus_dphi4_flat,event_weight);
					
				}
			}

			// Start loop over gen jets
			float leadgenjet_pt=-999, leadgenjet_eta=-999, leadgenjet_phi=-999, leadgenjet_mass=-999; // leading jet quantities
			float sublgenjet_pt=-999, sublgenjet_eta=-999, sublgenjet_phi=-999, sublgenjet_mass=-999;; // subleading jet quantities
			float thirdgenjet_pt=-999, thirdgenjet_eta=-999, thirdgenjet_phi=-999, thirdgenjet_mass=-999;; // third jet quantities

			bool isgjetincluded = false;

			for(int j = 0; j < gen_jetsize; j++){
				
				// Define jet kinematics
				float gjet_pt = gen_jtpt[j];
				float gjet_eta = gen_jteta[j];
				float gjet_phi = gen_jtphi[j];
				float gjet_mass = gen_jtmass[j];
				
				double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gjet_pt); // Jet weight (specially for MC)
				//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
				jet_weight = jet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gjet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)
		
				if(gjet_eta < -5.0 || gjet_eta > 5.0) continue; // no accept jets with |eta| > 4

				find_leading_subleading_third(gjet_pt,gjet_eta,gjet_phi,gjet_mass,leadgenjet_pt,leadgenjet_eta,leadgenjet_phi,leadgenjet_mass,sublgenjet_pt,sublgenjet_eta,sublgenjet_phi,sublgenjet_mass,thirdgenjet_pt,thirdgenjet_eta,thirdgenjet_phi,thirdgenjet_mass); // Find leading and subleading jets

				gjet_eta = gjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				if(colliding_system == "pPb" && is_pgoing && invert_pgoing){gjet_eta = -gjet_eta;}

				// Fill gen jet QA histograms
				double x_gen_jet[5]={gjet_pt,gjet_eta,gjet_phi,(double) multcentbin, (double) extrabin}; 
				hist_gen_jet->Fill(x_gen_jet);
				hist_gen_jet_weighted->Fill(x_gen_jet,event_weight*jet_weight);

				if((gjet_pt > jet_pt_min_cut && gjet_pt < jet_pt_max_cut) && (gjet_eta > jet_eta_min_cut && gjet_eta < jet_eta_max_cut)){  // Gen jet pT and eta cut

					isgjetincluded=true;
					// Fill gen jet vectors
					TVector3 GoodJets_gen;
					GoodJets_gen.SetPtEtaPhi(gjet_pt, gjet_eta, gjet_phi);
					jets_gen.push_back(GoodJets_gen);
					jet_w_gen.push_back(jet_weight);

					if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
						// Psi 2
						double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;
						double x_jetep2_plus[3]={deltaphi2PC(gjet_phi, Psi2_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP2_inclusive_plus->Fill(x_jetep2_plus,event_weight*jet_weight);
						double Psi2_EP_flat_minus = (double) EP_Psi2_minus_flat;
						double x_jetep2_minus[3]={deltaphi2PC(gjet_phi, Psi2_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP2_inclusive_minus->Fill(x_jetep2_minus,event_weight*jet_weight);
						// Psi 3
						double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;
						double x_jetep3_plus[3]={deltaphi2PC(gjet_phi, Psi3_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP3_inclusive_plus->Fill(x_jetep3_plus,event_weight*jet_weight);
						double Psi3_EP_flat_minus = (double) EP_Psi3_minus_flat;
						double x_jetep3_minus[3]={deltaphi2PC(gjet_phi, Psi3_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP3_inclusive_minus->Fill(x_jetep3_minus,event_weight*jet_weight);
						// Psi 4
						double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;
						double x_jetep4_plus[3]={deltaphi2PC(gjet_phi, Psi4_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP4_inclusive_plus->Fill(x_jetep4_plus,event_weight*jet_weight);
						double Psi4_EP_flat_minus = (double) EP_Psi4_minus_flat;
						double x_jetep4_minus[3]={deltaphi2PC(gjet_phi, Psi4_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP4_inclusive_minus->Fill(x_jetep4_minus,event_weight*jet_weight);
					}
				}
			}
			
			if(isgjetincluded){gen_mult_withonejet_weighted->Fill(genmult, event_weight);}
			bool isgdijet = false;

			bool removethirdjet_gen = false;
			if(do_thirdjet_removal){if(thirdgenjet_pt > 0.5*sublgenjet_pt) removethirdjet_gen = true;}

			
			//leading/subleading jets
			if(gen_jetsize > 1 && !removethirdjet_gen){

				double ljet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadgenjet_pt); // Jet weight (specially for MC)
				//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
				if(do_jet_smearing) ljet_weight = ljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadgenjet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

				double sljet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublgenjet_pt); // Jet weight (specially for MC)
				//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
				if(do_jet_smearing) sljet_weight = sljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublgenjet_pt, do_jet_smearing, 0.663);  // Jet smearing (For systematics)

				//leading/subleading pT cuts
				if(leadgenjet_pt > leading_pT_min && sublgenjet_pt > subleading_pT_min){ 

					double leadgenjet_eta_lab = leadgenjet_eta; // before boost for eta dijet
					double sublgenjet_eta_lab = sublgenjet_eta; // before boost for eta dijet
					leadgenjet_eta = leadgenjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
					sublgenjet_eta = sublgenjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
					if(colliding_system == "pPb" && is_pgoing && invert_pgoing){
						leadgenjet_eta = -leadgenjet_eta; 
						sublgenjet_eta = -sublgenjet_eta;
						leadgenjet_eta_lab = -leadgenjet_eta_lab; 
						sublgenjet_eta_lab = -sublgenjet_eta_lab;
					}

					// Fill leading/subleading jet quenching quantities
					double delta_phi_gen = fabs(deltaphi(leadgenjet_phi, sublgenjet_phi));
					double Aj_gen = asymmetry(leadgenjet_pt,sublgenjet_pt);
					double Xj_gen = xjvar(leadgenjet_pt,sublgenjet_pt);
					float ptdijet = 0.5*(leadrecojet_pt + sublrecojet_pt);
					//int ptdijetbin = (int) find_my_bin(pt_ave_bins, (float) ptdijet);
					double ptdijetbin = (double) ptdijet;
					double x_gen[6]={Xj_gen,Aj_gen,delta_phi_gen,(double)multcentbin,(double)ptdijetbin,(double)extrabin}; 
					// combinations of midrapidity, forward and backward
					bool leadmidrap = (leadgenjet_eta > jet_eta_min_cut && leadgenjet_eta < jet_eta_max_cut);
					bool sublmidrap = (sublgenjet_eta > jet_eta_min_cut && sublgenjet_eta < jet_eta_max_cut);
					bool leadfwdrap = (leadgenjet_eta > jet_fwd_eta_min_cut && leadgenjet_eta < jet_fwd_eta_max_cut);
					bool sublfwdrap = (sublgenjet_eta > jet_fwd_eta_min_cut && sublgenjet_eta < jet_fwd_eta_max_cut);
					bool leadbkwrap = (leadgenjet_eta > jet_bkw_eta_min_cut && leadgenjet_eta < jet_bkw_eta_max_cut);
					bool sublbkwrap = (sublgenjet_eta > jet_bkw_eta_min_cut && sublgenjet_eta < jet_bkw_eta_max_cut);
					if(do_dijetstudies){
						if(leadmidrap && sublmidrap) hist_gen_lead_gen_subl_quench_mid_mid->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadmidrap && sublfwdrap) hist_gen_lead_gen_subl_quench_mid_fwd->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadmidrap && sublbkwrap) hist_gen_lead_gen_subl_quench_mid_bkw->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadfwdrap && sublmidrap) hist_gen_lead_gen_subl_quench_fwd_mid->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadbkwrap && sublmidrap) hist_gen_lead_gen_subl_quench_bkw_mid->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadfwdrap && sublfwdrap) hist_gen_lead_gen_subl_quench_fwd_fwd->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadfwdrap && sublbkwrap) hist_gen_lead_gen_subl_quench_fwd_bkw->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadbkwrap && sublfwdrap) hist_gen_lead_gen_subl_quench_bkw_fwd->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
						if(leadbkwrap && sublbkwrap) hist_gen_lead_gen_subl_quench_bkw_bkw->Fill(x_gen,event_weight*ljet_weight*sljet_weight);
					}

					if(colliding_system=="pPb" && year_of_datataking==2016 && do_dijetstudies){

						double etadijet = 0.5*(leadgenjet_eta_lab + sublgenjet_eta_lab);
						double etadiff = 0.5*deltaeta(leadgenjet_eta_lab,sublgenjet_eta_lab);
						double xp_gen = 2.0*(ptdijet*exp(etadijet)*cosh(etadiff))/sqrts;
						double xPb_gen = 2.0*(ptdijet*exp(-etadijet)*cosh(etadiff))/sqrts;
						double m12_gen = 2.0*ptdijet*cosh(etadiff); // for future, if needed 	
						double x_dijet_gen[11] = {etadijet, etadiff, Xj_gen, Aj_gen, delta_phi_gen, xp_gen, xPb_gen, m12_gen, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

						double etadijet_CM = 0.5*(leadgenjet_eta + sublgenjet_eta);
						double etadiff_CM = 0.5*deltaeta(leadgenjet_eta,sublgenjet_eta);
						double xp_gen_CM = 2.0*(ptdijet*exp(etadijet_CM)*cosh(etadiff_CM))/sqrts;
						double xPb_gen_CM = 2.0*(ptdijet*exp(-etadijet_CM)*cosh(etadiff_CM))/sqrts;
						double m12_gen_CM = 2.0*ptdijet*cosh(etadiff_CM); // for future, if needed 	
						double x_dijet_gen_CM[11] = {etadijet_CM, etadiff_CM, Xj_gen, Aj_gen, delta_phi_gen, xp_gen_CM, xPb_gen_CM, m12_gen_CM, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

						ROOT::Math::PtEtaPhiMVector Lead_gen_jet(leadgenjet_pt,leadgenjet_eta,leadgenjet_phi,leadgenjet_mass);
						ROOT::Math::PtEtaPhiMVector Subl_gen_jet(sublgenjet_pt,sublgenjet_eta,sublgenjet_phi,sublgenjet_mass);			
						double etadijet_CM_y = 0.5*(Lead_gen_jet.Rapidity() + Subl_gen_jet.Rapidity());
						double etadiff_CM_y = 0.5*deltaeta(Lead_gen_jet.Rapidity(),Subl_gen_jet.Rapidity());
						double xp_gen_CM_y = 2.0*(ptdijet*exp(etadijet_CM_y)*cosh(etadiff_CM_y))/sqrts;
						double xPb_gen_CM_y = 2.0*(ptdijet*exp(-etadijet_CM_y)*cosh(etadiff_CM_y))/sqrts;
						double m12_gen_CM_y = 2.0*ptdijet*cosh(etadiff_CM_y); // for future, if needed 					
						double x_dijet_gen_CM_y[11] = {etadijet_CM_y, etadiff_CM_y, Xj_gen, Aj_gen, delta_phi_gen, xp_gen_CM_y, xPb_gen_CM_y, m12_gen_CM_y, (double)multcentbin, (double)ptdijetbin, (double)extrabin};

						if(fabs(leadgenjet_eta_lab) < dijetetamax && fabs(sublgenjet_eta_lab) < dijetetamax){
							hist_etaDijet_gen->Fill(x_dijet_gen,event_weight*ljet_weight*sljet_weight); 
							hist_etaDijet_CM_gen->Fill(x_dijet_gen_CM,event_weight*ljet_weight*sljet_weight);
							hist_yDijet_CM_gen->Fill(x_dijet_gen_CM_y,event_weight*ljet_weight*sljet_weight);
							}
					}

					// leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
					if(delta_phi_gen > leading_subleading_deltaphi_min){

						if((Xj_gen >= xjmin && Xj_gen < xjmax) && (Aj_gen >= Ajmin && Aj_gen < Ajmax)){

							// Fill leading and subleading jet QA histograms
							double x_lead[5]={leadgenjet_pt,leadgenjet_eta,leadgenjet_phi,(double) multcentbin,(double)extrabin}; 
							hist_gen_leadjet->Fill(x_lead);
							hist_gen_leadjet_weighted->Fill(x_lead,event_weight*ljet_weight);
							double x_sublead[5]={sublgenjet_pt,sublgenjet_eta,sublgenjet_phi,(double) multcentbin,(double)extrabin}; 
							hist_gen_subljet->Fill(x_sublead);
							hist_gen_subljet_weighted->Fill(x_sublead,event_weight*sljet_weight);
							
							isgdijet = true;
							pass_Aj_or_Xj_gen_cut = true; // if we apply Xj or Aj cuts

							bool usethisjets = false;
							if(fwdbkw_jettrk_option == "mid_mid"){ if(leadmidrap && sublmidrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "mid_fwd"){ if(leadmidrap && sublfwdrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "mid_bkw"){ if(leadmidrap && sublbkwrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "fwd_mid"){ if(leadfwdrap && sublmidrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "bkw_mid"){ if(leadbkwrap && sublmidrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "fwd_fwd"){ if(leadfwdrap && sublfwdrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "fwd_bkw"){ if(leadfwdrap && sublbkwrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "bkw_fwd"){ if(leadbkwrap && sublfwdrap) usethisjets = true; }
							if(fwdbkw_jettrk_option == "bkw_bkw"){ if(leadbkwrap && sublbkwrap) usethisjets = true; }

							// Fill leading and subleading jet vectors
							TVector3 GoodLeadingJets_gen;
							TVector3 GoodSubLeadingJets_gen;

							if(usethisjets){
								GoodLeadingJets_gen.SetPtEtaPhi(leadgenjet_pt, leadgenjet_eta, leadgenjet_phi);
								lead_jets_gen.push_back(GoodLeadingJets_gen);
								lead_jet_w_gen.push_back(ljet_weight);
								GoodSubLeadingJets_gen.SetPtEtaPhi(sublgenjet_pt, sublgenjet_eta, sublgenjet_phi);
								subl_jets_gen.push_back(GoodSubLeadingJets_gen);
								subl_jet_w_gen.push_back(sljet_weight);
							}
						
							if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
								// Psi 2
								double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;
								double Psi2_EP_flat_minus = (double) EP_Psi2_minus_flat;
								double x_ljetep2_plus[3]={deltaphi2PC(leadgenjet_pt, Psi2_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP2_leading_plus->Fill(x_ljetep2_plus,event_weight*ljet_weight);
								double x_ljetep2_minus[3]={deltaphi2PC(leadgenjet_pt, Psi2_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP2_leading_minus->Fill(x_ljetep2_minus,event_weight*ljet_weight);
								double x_sljetep2_plus[3]={deltaphi2PC(sublgenjet_phi, Psi2_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP2_subleading_plus->Fill(x_sljetep2_plus,event_weight*sljet_weight);
								double x_sljetep2_minus[3]={deltaphi2PC(sublgenjet_phi, Psi2_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP2_subleading_minus->Fill(x_sljetep2_minus,event_weight*sljet_weight);
								// Psi 3
								double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;
								double Psi3_EP_flat_minus = (double) EP_Psi3_minus_flat;
								double x_ljetep3_plus[3]={deltaphi2PC(leadgenjet_pt, Psi3_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP3_leading_plus->Fill(x_ljetep3_plus,event_weight*ljet_weight);
								double x_ljetep3_minus[3]={deltaphi2PC(leadgenjet_pt, Psi3_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP3_leading_minus->Fill(x_ljetep3_minus,event_weight*ljet_weight);
								double x_sljetep3_plus[3]={deltaphi2PC(sublgenjet_phi, Psi3_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP3_subleading_plus->Fill(x_sljetep3_plus,event_weight*sljet_weight);
								double x_sljetep3_minus[3]={deltaphi2PC(sublgenjet_phi, Psi3_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP3_subleading_minus->Fill(x_sljetep3_minus,event_weight*sljet_weight);
								// Psi 4
								double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;
								double Psi4_EP_flat_minus = (double) EP_Psi4_minus_flat;
								double x_ljetep4_plus[3]={deltaphi2PC(leadgenjet_pt, Psi4_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP4_leading_plus->Fill(x_ljetep4_plus,event_weight*ljet_weight);
								double x_ljetep4_minus[3]={deltaphi2PC(leadgenjet_pt, Psi4_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP4_leading_minus->Fill(x_ljetep4_minus,event_weight*ljet_weight);
								double x_sljetep4_plus[3]={deltaphi2PC(sublgenjet_phi, Psi4_EP_flat_plus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP4_subleading_plus->Fill(x_sljetep4_plus,event_weight*sljet_weight);
								double x_sljetep4_minus[3]={deltaphi2PC(sublgenjet_phi, Psi4_EP_flat_minus),(double) multcentbin,(double) extrabin}; Dphi_GEN_flat_EP4_subleading_minus->Fill(x_sljetep4_minus,event_weight*sljet_weight);
							}
						}
					}
				}
			}
			
			if(isgdijet){gen_mult_withdijets_weighted->Fill(genmult, event_weight);}

			// Measure correlations and fill mixing vectors for inclusive jet+track correlations
			if(do_inclusejettrack_correlation && !removethirdjet_gen){
				// Reco-Gen
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(jets_reco, jet_w_reco, tracks_gen, track_w_gen, hist_correlation_signal_jet_reco_track_gen, hist_jet_from_reco_gen_sig, hist_trk_from_reco_gen_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_jet_reco_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_jet_reco_track_gen,JetR,hist_injet_reco_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_recogen, jets_reco, jet_w_reco, tracks_gen, track_w_gen, mult, vertexz, event_weight, extra_variable, ev_jet_vector_reco_gen, jet_weights_reco_gen, ev_track_vector_reco_gen, trk_weights_reco_gen, multvec_reco_gen, vzvec_reco_gen, weights_reco_gen, extravec_reco_gen);	// for mixing --> store vectors for mixing
				// Gen-Reco
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(jets_gen, jet_w_gen, tracks_reco, track_w_reco, hist_correlation_signal_jet_gen_track_reco, hist_jet_from_gen_reco_sig, hist_trk_from_gen_reco_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_jet_gen_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_jet_gen_track_reco,JetR,hist_injet_gen_track_reco,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_genreco, jets_gen, jet_w_gen, tracks_reco, track_w_reco, mult, vertexz, event_weight, extra_variable, ev_jet_vector_gen_reco, jet_weights_gen_reco, ev_track_vector_gen_reco, trk_weights_gen_reco, multvec_gen_reco, vzvec_gen_reco, weights_gen_reco, extravec_gen_reco);	// for mixing --> store vectors for mixing
				// Gen-Gen
				if(pass_Aj_or_Xj_gen_cut) correlation(jets_gen, jet_w_gen, tracks_gen, track_w_gen, hist_correlation_signal_jet_gen_track_gen, hist_jet_from_gen_gen_sig, hist_trk_from_gen_gen_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_jet_gen_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_jet_gen_track_gen,JetR,hist_injet_gen_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_gengen, jets_gen, jet_w_gen, tracks_gen, track_w_gen, mult, vertexz, event_weight, extra_variable, ev_jet_vector_gen_gen, jet_weights_gen_gen, ev_track_vector_gen_gen, trk_weights_gen_gen, multvec_gen_gen, vzvec_gen_gen, weights_gen_gen, extravec_gen_gen);	// for mixing --> store vectors for mixing
			}
			// Measure correlations and fill mixing vectors for (leading/subleading) jet+track correlations
			if(do_leading_subleading_jettrack_correlation && !removethirdjet_gen){
				// Reco-Gen
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(lead_jets_reco, lead_jet_w_reco, tracks_gen, track_w_gen, hist_correlation_signal_lead_jet_reco_track_gen, hist_lead_jet_from_reco_gen_sig, hist_LJ_trk_from_reco_gen_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_reco_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_lead_jet_reco_track_gen,JetR,hist_inLeadjet_reco_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(subl_jets_reco, subl_jet_w_reco, tracks_gen, track_w_gen, hist_correlation_signal_subl_jet_reco_track_gen, hist_subl_jet_from_reco_gen_sig, hist_SLJ_trk_from_reco_gen_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_reco_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_subl_jet_reco_track_gen,JetR,hist_inSubljet_reco_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_recogen_lead, lead_jets_reco, lead_jet_w_reco, tracks_gen, track_w_gen, mult, vertexz, event_weight, extra_variable, ev_jet_vector_leadjet_reco_gen, jet_weights_leadjet_reco_gen, ev_track_vector_leadjet_reco_gen, trk_weights_leadjet_reco_gen, multvec_leadjet_reco_gen, vzvec_leadjet_reco_gen, weights_leadjet_reco_gen, extravec_leadjet_reco_gen);	// for mixing --> store vectors for mixing
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_recogen_subl, subl_jets_reco, subl_jet_w_reco, tracks_gen, track_w_gen, mult, vertexz, event_weight, extra_variable, ev_jet_vector_subleadjet_reco_gen, jet_weights_subleadjet_reco_gen, ev_track_vector_subleadjet_reco_gen, trk_weights_subleadjet_reco_gen, multvec_subleadjet_reco_gen, vzvec_subleadjet_reco_gen, weights_subleadjet_reco_gen, extravec_subleadjet_reco_gen);	// for mixing --> store vectors for mixing
				// Gen-Reco	
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(lead_jets_gen, lead_jet_w_gen, tracks_reco, track_w_reco,hist_correlation_signal_lead_jet_gen_track_reco, hist_lead_jet_from_gen_reco_sig, hist_LJ_trk_from_gen_reco_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_gen_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_lead_jet_gen_track_reco,JetR,hist_inLeadjet_gen_track_reco,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(subl_jets_gen, subl_jet_w_gen, tracks_reco, track_w_reco,hist_correlation_signal_subl_jet_gen_track_reco, hist_subl_jet_from_gen_reco_sig, hist_SLJ_trk_from_gen_reco_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_gen_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_subl_jet_gen_track_reco,JetR,hist_inSubljet_gen_track_reco,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_genreco_lead, lead_jets_gen, lead_jet_w_gen, tracks_reco, track_w_reco, mult, vertexz, event_weight, extra_variable, ev_jet_vector_leadjet_gen_reco, jet_weights_leadjet_gen_reco, ev_track_vector_leadjet_gen_reco, trk_weights_leadjet_gen_reco, multvec_leadjet_gen_reco, vzvec_leadjet_gen_reco, weights_leadjet_gen_reco, extravec_leadjet_gen_reco);	// for mixing --> store vectors for mixing
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_genreco_subl, subl_jets_gen, subl_jet_w_gen, tracks_reco, track_w_reco, mult, vertexz, event_weight, extra_variable, ev_jet_vector_subleadjet_gen_reco, jet_weights_subleadjet_gen_reco, ev_track_vector_subleadjet_gen_reco, trk_weights_subleadjet_gen_reco, multvec_subleadjet_gen_reco, vzvec_subleadjet_gen_reco, weights_subleadjet_gen_reco, extravec_subleadjet_gen_reco);	// for mixing --> store vectors for mixing
				// Gen-Gen
				if(pass_Aj_or_Xj_gen_cut) correlation(lead_jets_gen, lead_jet_w_gen, tracks_gen, track_w_gen, hist_correlation_signal_lead_jet_gen_track_gen, hist_lead_jet_from_gen_gen_sig, hist_LJ_trk_from_gen_gen_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_gen_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_lead_jet_gen_track_gen,JetR,hist_inLeadjet_gen_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_gen_cut) correlation(subl_jets_gen, subl_jet_w_gen, tracks_gen, track_w_gen, hist_correlation_signal_subl_jet_gen_track_gen, hist_subl_jet_from_gen_gen_sig, hist_SLJ_trk_from_gen_gen_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_gen_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_subl_jet_gen_track_gen,JetR,hist_inSubljet_gen_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_gengen_lead, lead_jets_gen, lead_jet_w_gen, tracks_gen, track_w_gen, mult, vertexz, event_weight, extra_variable, ev_jet_vector_leadjet_gen_gen, jet_weights_leadjet_gen_gen, ev_track_vector_leadjet_gen_gen, trk_weights_leadjet_gen_gen, multvec_leadjet_gen_gen, vzvec_leadjet_gen_gen, weights_leadjet_gen_gen, extravec_leadjet_gen_gen);	// for mixing --> store vectors for mixing
				if(pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_gengen_subl, subl_jets_gen, subl_jet_w_gen, tracks_gen, track_w_gen, mult, vertexz, event_weight, extra_variable, ev_jet_vector_subleadjet_gen_gen, jet_weights_subleadjet_gen_gen, ev_track_vector_subleadjet_gen_gen, trk_weights_subleadjet_gen_gen, multvec_subleadjet_gen_gen, vzvec_subleadjet_gen_gen, weights_subleadjet_gen_gen, extravec_subleadjet_gen_gen);	// for mixing --> store vectors for mixing
			}
			
			if(do_flow){twoparticlecorrelation(tracks_gen, track_w_gen, hist_gen_gen_2pcorrelation_signal, event_weight, mult, extra_variable, sube_tracks_gen, hist_gen_gen_2pcorrelation_signal_subg0, hist_gen_gen_2pcorrelation_signal_subcross);} // calculate 2 particle correlations
			
		}

	} // End loop over events

	// Event mixing 
	if(do_mixing){
		cout << endl;
		cout << "============= It is mixing time! =============" << endl;
		cout << endl;
		sec_start_mix = clock(); // Start timing measurement
		cout << "Running --> RECO-RECO" << endl;
		// RECO-RECO
		if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_reco_reco, multiplicity_centrality_bins, extravec_reco_reco, extra_bins, vzvec_reco_reco, DVz_range, ev_jet_vector_reco_reco, jet_weights_reco_reco, ev_track_vector_reco_reco, trk_weights_reco_reco, hist_correlation_mixing_jet_reco_track_reco, trk_pt_bins, weights_reco_reco, hist_jet_from_reco_reco_mix, hist_trk_from_reco_reco_mix, double_weight_mix, do_flow);
		if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_reco_reco, multiplicity_centrality_bins, extravec_reco_reco, extra_bins, vzvec_leadjet_reco_reco, DVz_range, ev_jet_vector_leadjet_reco_reco, jet_weights_leadjet_reco_reco, ev_track_vector_leadjet_reco_reco, trk_weights_leadjet_reco_reco, hist_correlation_mixing_lead_jet_reco_track_reco, trk_pt_bins, weights_leadjet_reco_reco, hist_lead_jet_from_reco_reco_mix, hist_LJ_trk_from_reco_reco_mix, double_weight_mix, do_flow);
		if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_reco_reco, multiplicity_centrality_bins, extravec_reco_reco, extra_bins, vzvec_subleadjet_reco_reco, DVz_range, ev_jet_vector_subleadjet_reco_reco, jet_weights_subleadjet_reco_reco, ev_track_vector_subleadjet_reco_reco, trk_weights_subleadjet_reco_reco, hist_correlation_mixing_subl_jet_reco_track_reco, trk_pt_bins, weights_subleadjet_reco_reco, hist_subl_jet_from_reco_reco_mix, hist_SLJ_trk_from_reco_reco_mix, double_weight_mix, do_flow);
		if(do_flow){cout << "Running 2PC --> GEN-GEN" << endl; call_mix_random_2pc(N_ev_mix, Mult_or_Cent_range, multvec_reco_reco, multiplicity_centrality_bins, extravec_reco_reco, extra_bins, vzvec_reco_reco, DVz_range, ev_track_vector_reco_reco, trk_weights_reco_reco, hist_reco_reco_2pcorrelation_mixing, trk_pt_bins, weights_reco_reco, double_weight_mix);}
		if(is_MC){
			cout << "Running --> RECO-GEN" << endl;
			// RECO-GEN
			if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_reco_gen, multiplicity_centrality_bins, extravec_reco_gen, extra_bins, vzvec_reco_gen, DVz_range, ev_jet_vector_reco_gen, jet_weights_reco_gen, ev_track_vector_reco_gen, trk_weights_reco_gen, hist_correlation_mixing_jet_reco_track_gen, trk_pt_bins, weights_reco_gen, hist_jet_from_reco_gen_mix, hist_trk_from_reco_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_reco_gen, multiplicity_centrality_bins, extravec_reco_gen, extra_bins, vzvec_leadjet_reco_gen, DVz_range, ev_jet_vector_leadjet_reco_gen, jet_weights_leadjet_reco_gen, ev_track_vector_leadjet_reco_gen, trk_weights_leadjet_reco_gen, hist_correlation_mixing_lead_jet_reco_track_gen, trk_pt_bins, weights_leadjet_reco_gen, hist_lead_jet_from_reco_gen_mix, hist_LJ_trk_from_reco_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_reco_gen, multiplicity_centrality_bins, extravec_reco_gen, extra_bins, vzvec_subleadjet_reco_gen, DVz_range, ev_jet_vector_subleadjet_reco_gen, jet_weights_subleadjet_reco_gen, ev_track_vector_subleadjet_reco_gen, trk_weights_subleadjet_reco_gen, hist_correlation_mixing_subl_jet_reco_track_gen, trk_pt_bins, weights_subleadjet_reco_gen, hist_subl_jet_from_reco_gen_mix, hist_SLJ_trk_from_reco_gen_mix, double_weight_mix, do_flow);
			cout << "Running --> GEN-RECO" << endl;
			// GEN-RECO
			if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_gen_reco, multiplicity_centrality_bins, extravec_gen_reco, extra_bins, vzvec_gen_reco, DVz_range, ev_jet_vector_gen_reco, jet_weights_gen_reco, ev_track_vector_gen_reco, trk_weights_gen_reco, hist_correlation_mixing_jet_gen_track_reco, trk_pt_bins, weights_gen_reco, hist_jet_from_gen_reco_mix, hist_trk_from_gen_reco_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_gen_reco, multiplicity_centrality_bins, extravec_gen_reco, extra_bins, vzvec_leadjet_gen_reco, DVz_range, ev_jet_vector_leadjet_gen_reco, jet_weights_leadjet_gen_reco, ev_track_vector_leadjet_gen_reco, trk_weights_leadjet_gen_reco, hist_correlation_mixing_lead_jet_gen_track_reco, trk_pt_bins, weights_leadjet_gen_reco, hist_lead_jet_from_gen_reco_mix, hist_LJ_trk_from_gen_reco_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_gen_reco, multiplicity_centrality_bins, extravec_gen_reco, extra_bins, vzvec_subleadjet_gen_reco, DVz_range, ev_jet_vector_subleadjet_gen_reco, jet_weights_subleadjet_gen_reco, ev_track_vector_subleadjet_gen_reco, trk_weights_subleadjet_gen_reco, hist_correlation_mixing_subl_jet_gen_track_reco, trk_pt_bins, weights_subleadjet_gen_reco, hist_subl_jet_from_gen_reco_mix, hist_SLJ_trk_from_gen_reco_mix, double_weight_mix, do_flow);
			cout << "Running --> GEN-GEN" << endl;
			// GEN-GEN
			if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_gen_gen, multiplicity_centrality_bins, extravec_gen_gen, extra_bins, vzvec_gen_gen, DVz_range, ev_jet_vector_gen_gen, jet_weights_gen_gen, ev_track_vector_gen_gen, trk_weights_gen_gen, hist_correlation_mixing_jet_gen_track_gen, trk_pt_bins, weights_gen_gen, hist_jet_from_gen_gen_mix, hist_trk_from_gen_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_gen_gen, multiplicity_centrality_bins, extravec_gen_gen, extra_bins, vzvec_leadjet_gen_gen, DVz_range, ev_jet_vector_leadjet_gen_gen, jet_weights_leadjet_gen_gen, ev_track_vector_leadjet_gen_gen, trk_weights_leadjet_gen_gen, hist_correlation_mixing_lead_jet_gen_track_gen, trk_pt_bins, weights_leadjet_gen_gen, hist_lead_jet_from_gen_gen_mix, hist_LJ_trk_from_gen_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_gen_gen, multiplicity_centrality_bins, extravec_gen_gen, extra_bins, vzvec_subleadjet_gen_gen, DVz_range, ev_jet_vector_subleadjet_gen_gen, jet_weights_subleadjet_gen_gen, ev_track_vector_subleadjet_gen_gen, trk_weights_subleadjet_gen_gen, hist_correlation_mixing_subl_jet_gen_track_gen, trk_pt_bins, weights_subleadjet_gen_gen, hist_subl_jet_from_gen_gen_mix, hist_SLJ_trk_from_gen_gen_mix, double_weight_mix, do_flow);
			if(do_flow){ cout << "Running 2PC --> GEN-GEN" << endl; call_mix_random_2pc(N_ev_mix, Mult_or_Cent_range, multvec_gen_gen, multiplicity_centrality_bins, extravec_gen_gen, extra_bins, vzvec_gen_gen, DVz_range, ev_track_vector_gen_gen, trk_weights_gen_gen, hist_gen_gen_2pcorrelation_mixing, trk_pt_bins, weights_gen_gen, double_weight_mix);}
			cout << " --> Mixing DONE! " << endl;
		}
		sec_end_mix = clock(); // Stop time counting
		cout << endl;

		cout << "========================================" << endl;
		cout << "Mixing time: " << (double)(sec_end_mix - sec_start_mix) / (CLOCKS_PER_SEC) << " [s]" << endl;
		cout << "========================================" << endl;
	}

	// Output file name
	cout << endl;
	cout << "Writing histograms on ";
	cout << endl;

	// Make an output file
	string file_output = Form("%s_%s_%s_%iGeV_%s_%s_%s_Jptmin_%.1f_Jptmax_%.1f_Jetamin_%.1f_Jetamax_%.1f_%s%s%s_%s_%s_%s_%i",ouputfilename.Data(),colliding_system.Data(),data_or_mc.Data(),sNN_energy_GeV,jet_type.Data(),jet_collection.Data(),jet_trigger.Data(),jet_pt_min_cut,jet_pt_max_cut,jet_eta_min_cut,jet_eta_max_cut,jet_axis.Data(),smear.Data(),XjAj.Data(),ref_sample.Data(),particles.Data(),isflow.Data(),date->GetDate()); // output file
	std::replace(file_output.begin(), file_output.end(), '.', 'p'); // replace . to p
	std::replace(file_output.begin(), file_output.end(), '-', 'N'); // replace - to N for negative

	// Open, write and close the output file
	TFile *MyFile = new TFile(Form("%s.root", file_output.c_str()), "RECREATE");
	if(MyFile->IsOpen()) cout << "output file: " << file_output.c_str() << ".root" << endl;
	MyFile->cd(); 

	// Write in different folders (see histogram_definition.h)
	// Control plots
	MyFile->mkdir("QA_histograms"); 
	MyFile->cd("QA_histograms"); 
	w_QA_hist(is_MC); 

	// Jet+Track correlations
	if(do_inclusejettrack_correlation || do_leading_subleading_jettrack_correlation){
		MyFile->mkdir("correlation_reco_reco_histograms"); 
		MyFile->cd("correlation_reco_reco_histograms");
		w_recoreco_hist(do_mixing,do_rotation,do_inclusejettrack_correlation,do_leading_subleading_jettrack_correlation);
		if(is_MC){
			MyFile->mkdir("correlation_reco_gen_histograms"); 
			MyFile->cd("correlation_reco_gen_histograms");
			w_recogen_hist(do_mixing,do_rotation,do_inclusejettrack_correlation,do_leading_subleading_jettrack_correlation);
			MyFile->mkdir("correlation_gen_reco_histograms"); 
			MyFile->cd("correlation_gen_reco_histograms");
			w_genreco_hist(do_mixing,do_rotation,do_inclusejettrack_correlation,do_leading_subleading_jettrack_correlation);
			MyFile->mkdir("correlation_gen_gen_histograms"); 
			MyFile->cd("correlation_gen_gen_histograms");
			w_gengen_hist(do_mixing,do_rotation,do_inclusejettrack_correlation,do_leading_subleading_jettrack_correlation);		
		}
	}

	// Dijet histograms	
	if(do_dijetstudies){
		MyFile->mkdir("dijet_histograms"); 
		MyFile->cd("dijet_histograms");  
		w_dijet_hist(is_MC);
	}

	// 2PC histograms		
	if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
		MyFile->mkdir("eventplane_histograms"); 
		MyFile->cd("eventplane_histograms");  
		w_ep_hist(is_MC);
		MyFile->mkdir("TwoPC_histograms"); 
		MyFile->cd("TwoPC_histograms");  
		w_2pc_hist(is_MC, do_mixing);
	}

	MyFile->Close();

	cout << endl;
	cout << "------------------------------------- DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

	print_stop(); // Print time, date and hour when it stops
	return 0;
}
