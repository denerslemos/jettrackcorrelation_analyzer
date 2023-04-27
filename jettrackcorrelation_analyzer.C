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

	clock_t sec_start, sec_end, sec_start_mix, sec_end_mix;
	sec_start = clock(); // start timing measurement

	TDatime* date = new TDatime();

	printwelcome(true); // welcome message

	print_start(); // start timing print
	bool is_MC;
	if(MCSim==0){is_MC = false;}else{is_MC = true;}

	if(!do_inclusejettrack_correlation && do_leading_subleading_jettrack_correlation) do_Xj_or_Ajcut = true; // if only leading and subleading jet+track correlation are measured we make sure that Xj and Aj are true

	bool do_pthatcut = true;
   	if(!is_MC) do_pthatcut = false; // MC only
	if(!is_MC) do_pid = false; // MC only

   	if(colliding_system!="pPb") do_CM_pPb = false; // Only do center-of-mass for pPb

	//print important informations in the output file
	TString data_or_mc;
	if(!is_MC){data_or_mc="Data";}else{data_or_mc="MC";}
	TString simev;
	if(similar_events){simev = "simevs";}else{simev = "";}
	TString ref_sample = "norefsample";
	if(do_mixing && !do_rotation){ref_sample = Form("mix%ievsMult%iDVz%.1f%s",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data());}else if(!do_mixing && do_rotation){ref_sample = Form("rot%ievs",N_of_rot);}else if(do_mixing && do_rotation){ref_sample = Form("mix%ievsMult%iDVz%.1f%s_rot%ievs",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data(),N_of_rot);}
	TString jet_axis;
	if(use_WTA){jet_axis = "WTA";}else{jet_axis = "ESC";}
	TString smear;	
	if(do_jet_smearing){smear = "_smearing_";}else{smear = "";}
	if(!do_pid) particles = "CH";
	TString jet_type;
	if(do_inclusejettrack_correlation && !do_leading_subleading_jettrack_correlation) jet_type = "Incl";
	if(!do_inclusejettrack_correlation && do_leading_subleading_jettrack_correlation) jet_type = "LeadSubl";
	if(do_inclusejettrack_correlation && do_leading_subleading_jettrack_correlation) jet_type = "InclLeadSubl";
	TString XjAj;
	if(do_Xj_or_Ajcut){XjAj = Form("_Ajmin_%.1f_Ajmax_%.1f_Xjmin_%.1f_Xjmax_%.1f",Ajmin,Ajmax,xjmin,xjmax);}else{XjAj = "";}
	TString isflow;
	if(do_flow){isflow="flow";}else{isflow="jetshape";}

	
	// In case of wrong input, printout error message and kill the job
	if(year_of_datataking!=2012 && year_of_datataking!=2016 && year_of_datataking!=2017 && year_of_datataking!=2018){cout << "Data and MC not supported: choose 2012 for pp at 8 TeV, 2016 for pPb at 8.16 TeV, 2017 for pp at 5.02 TeV or XeXe at 5.44 TeV and 2018 for PbPb at 5.02 TeV" << endl; return;}
	if(colliding_system!="pp" && colliding_system!="pPb" && colliding_system!="XeXe" && colliding_system!="PbPb"){cout << "Data and MC not supported: choose pp for proton-proton, pPb for proton-lead, PbPb for lead-lead and XeXe for xenon-xenon" << endl; return;}
	if(sNN_energy_GeV!=5020 && sNN_energy_GeV!=5440 && sNN_energy_GeV!=8000 && sNN_energy_GeV!=8160 && sNN_energy_GeV!=13000){cout << "Data and MC not supported: 5020 for pp 2017 or PbPb 2018, 5440 for XeXe, 8000 for pp 2018, 8160 for pPb 2016" << endl; return;}

	// Read JEC file
	vector<string> Files;
	Files.push_back(Form("aux_files/%s_%i/JEC/%s",colliding_system.Data(),sNN_energy_GeV,JEC_file.Data()));
	JetCorrector JEC(Files);

	// Track or particle efficiency file
	TFile *fileeff = TFile::Open(Form("aux_files/%s_%i/trk_eff_table/%s",colliding_system.Data(),sNN_energy_GeV,trk_eff_file.Data()));
	cout << endl;
	// Print the input in the screen/log 
	print_input(data_or_mc,fileeff,colliding_system,pthatmin,pthatmax);
	cout << endl;

	// Read the input file(s)
	fstream inputfile;
	inputfile.open(Form("%s",input_file.Data()), ios::in);
	if(!inputfile.is_open()){cout << "List of input files not founded!" << endl; return;}{cout << "List of input files founded! --> " << input_file.Data() << endl;}

	// Make a chain and a vector of file names
	std::vector<TString> file_name_vector;
	string file_chain;
	while(getline(inputfile, file_chain)){file_name_vector.push_back(file_chain.c_str());}
	inputfile.close();

	// Read the trees to be added in the Chain
	TChain *hlt_tree = new TChain("hltanalysis/HltTree");
	TChain *jet_tree = new TChain(Form("%s/t",jet_collection.Data()));
	TChain *trk_tree = new TChain("ppTrack/trackTree");
	TChain *hea_tree = new TChain("hiEvtAnalyzer/HiTree");
	TChain *gen_tree;
	if(is_MC){gen_tree = new TChain("HiGenParticleAna/hi");}
	TChain *ski_tree = new TChain("skimanalysis/HltTree");
	TChain *ep_tree;
	if(colliding_system=="pPb" && year_of_datataking==2016){ep_tree = new TChain("checkflattening/tree");}
	
	// add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		TFile testfile(*listIterator, "READ"); 
		if(testfile.IsZombie() || testfile.TestBit(TFile::kRecovered)) continue;
		cout << "Adding file " << *listIterator << " to the chains" << endl;
		hlt_tree->Add(*listIterator);
		trk_tree->Add(*listIterator);
		hea_tree->Add(*listIterator);
		jet_tree->Add(*listIterator);
		ski_tree->Add(*listIterator);
		if(is_MC){gen_tree->Add(*listIterator);}
		if(colliding_system=="pPb" && year_of_datataking==2016){ep_tree->Add(*listIterator);}
	}

	// Connect all chains
	hlt_tree->AddFriend(trk_tree);
	hlt_tree->AddFriend(hea_tree);
	hlt_tree->AddFriend(jet_tree);
	hlt_tree->AddFriend(ski_tree);	
	if(is_MC){hlt_tree->AddFriend(gen_tree);}
	if(colliding_system=="pPb" && year_of_datataking==2016){hlt_tree->AddFriend(ep_tree);}

    // Read the desired branchs in the trees
	read_tree(hlt_tree, is_MC, use_WTA, jet_trigger.Data(), colliding_system.Data(), sNN_energy_GeV, year_of_datataking, event_filter_str, event_filter_bool); // access the tree informations

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
		
		if(i != 0 && (i % 10000) == 0){double alpha = (double)i; cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;}

		//if(i != 0 && i % 10000 == 0 ) break; // just for tests (need to remove)

		Nevents->Fill(0); // filled after each event cut

		// Booleans to remove events which does not pass the Aj or Xj selection
		bool pass_Aj_or_Xj_reco_cut = true;
		bool pass_Aj_or_Xj_gen_cut = true;
		if(do_Xj_or_Ajcut){pass_Aj_or_Xj_reco_cut = false; pass_Aj_or_Xj_gen_cut = false;}

		// Apply trigger
		if(dojettrigger){
			if(jet_trigger_bit != 1) continue;
		}else{jet_trigger="nojettrig";}
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

		// ref jets
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
		if(vertexz <= vz_cut_min || vertexz >= vz_cut_max) continue;
		Nevents->Fill(3);
		//pthat (MC only)
		if(do_pthatcut){if(pthat <= pthatmin || pthat > pthatmax) continue;}
		Nevents->Fill(4);
		//multiplicity or centrality
		int trksize = (int)ntrk;
		int mult;
		if(use_centrality){mult = hiBin;}else{mult = get_Ntrkoff(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);}
		int recomult = get_simple_mult_reco(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge);
		int genmult;
		if(is_MC)genmult = get_simple_mult_gen(colliding_system, sNN_energy_GeV, year_of_datataking, (int) gen_trkpt->size(), gen_trketa, gen_trkpt, gen_trkchg);
		if(mult < multiplicity_centrality_bins[0] || mult > multiplicity_centrality_bins[multiplicity_centrality_bins.size()-1])continue; //centrality of multiplicity range	
		int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float) mult);
		Nevents->Fill(5);

		// event weight(s), this must be applied in all histograms
		double event_weight = get_event_weight(is_MC, use_centrality, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, vertexz, mult, weight, pthat); // get the event weight

		// Fill vertex, pthat and multiplicity/centrality histograms
		multiplicity->Fill(mult);
		multiplicity_weighted->Fill(mult, event_weight);
		
		reco_mult->Fill(recomult);
		reco_mult_weighted->Fill(recomult, event_weight);
		gen_mult->Fill(genmult);
		gen_mult_weighted->Fill(genmult, event_weight);
		
		vzhist->Fill(vertexz, (double) mult);
		vzhist_weighted->Fill(vertexz, (double) mult, event_weight);

		pthathist->Fill(pthat, (double) mult);
		pthathist_weighted->Fill(pthat, (double) mult, event_weight);

		if(colliding_system=="pPb" && year_of_datataking==2016){

			double x3D_hiHF[3]={hfplus,hfminus,(double) mult}; hfhist->Fill(x3D_hiHF); hfhist_weighted->Fill(x3D_hiHF,event_weight);
			double x3D_hiHFEta4[3]={hfplusEta4,hfminusEta4,(double) mult};hfhistEta4->Fill(x3D_hiHFEta4); hfhistEta4_weighted->Fill(x3D_hiHFEta4,event_weight);
			double x3D_hiZDC[3]={zdcplus,zdcminus,(double) mult}; zdchist->Fill(x3D_hiZDC); zdchist_weighted->Fill(x3D_hiZDC,event_weight);
			
			if(do_flow){
				//event plane information
				// Psi 2
				double mult2_plus = (double)EP_Mult2_plus;
				double mult2_minus = (double)EP_Mult2_minus;
				double q2vec_plus = (double) sqrt(EP_Qx2_plus*EP_Qx2_plus + EP_Qy2_plus*EP_Qy2_plus);
				double q2vec_minus = (double)sqrt(EP_Qx2_minus*EP_Qx2_minus + EP_Qy2_minus*EP_Qy2_minus);
				double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;	
				double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
				double x4D_reco_plus_ep2_flat[4]={mult2_plus,q2vec_plus,Psi2_EP_flat_plus,(double) mult}; 
				EP2_plus_flat->Fill(x4D_reco_plus_ep2_flat,event_weight);
				double x4D_reco_minus_ep2_flat[4]={mult2_minus,q2vec_minus,Psi2_EP_flat_minus,(double) mult}; 
				EP2_minus_flat->Fill(x4D_reco_minus_ep2_flat,event_weight);

				// Psi 3
				double mult3_plus = (double)EP_Mult3_plus;
				double mult3_minus = (double)EP_Mult3_minus;
				double q3vec_plus = (double) sqrt(EP_Qx3_plus*EP_Qx3_plus + EP_Qy3_plus*EP_Qy3_plus);
				double q3vec_minus = (double)sqrt(EP_Qx3_minus*EP_Qx3_minus + EP_Qy3_minus*EP_Qy3_minus);
				double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;	
				double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
				double x4D_reco_plus_ep3_flat[4]={mult3_plus,q3vec_plus,Psi3_EP_flat_plus,(double) mult}; 
				EP3_plus_flat->Fill(x4D_reco_plus_ep3_flat,event_weight);
				double x4D_reco_minus_ep3_flat[4]={mult3_minus,q3vec_minus,Psi3_EP_flat_minus,(double) mult}; 
				EP3_minus_flat->Fill(x4D_reco_minus_ep3_flat,event_weight);

				// Psi 4
				double mult4_plus = (double)EP_Mult4_plus;
				double mult4_minus = (double)EP_Mult4_minus;
				double q4vec_plus = (double) sqrt(EP_Qx4_plus*EP_Qx4_plus + EP_Qy4_plus*EP_Qy4_plus);
				double q4vec_minus = (double)sqrt(EP_Qx4_minus*EP_Qx4_minus + EP_Qy4_minus*EP_Qy4_minus);
				double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;	
				double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
				double x4D_reco_plus_ep4_flat[4]={mult4_plus,q4vec_plus,Psi4_EP_flat_plus,(double) mult}; 
				EP4_plus_flat->Fill(x4D_reco_plus_ep4_flat,event_weight);
				double x4D_reco_minus_ep4_flat[4]={mult4_minus,q4vec_minus,Psi4_EP_flat_minus,(double) mult}; 
				EP4_minus_flat->Fill(x4D_reco_minus_ep4_flat,event_weight);
			}
		}

		// Reconstruction level	(Data and MC)
		// Start loop over reco tracks (trksize is number of reco tracks)
		for (int j = 0; j < trksize; j++){ 

			// Define track/particle kinematics
		 	float trk_pt = trkpt[j];
		 	float trk_eta = trketa[j];
	 		float trk_phi = trkphi[j];

			// Apply track selection (see read_tree.h to see what each variable means)
			if(fabs(trk_eta) > trk_eta_cut) continue;
		 	if(trkpt[j] <= trk_pt_min_cut) continue;
		 	if(highpur[j] == false) continue;
		 	if(fabs(trkpterr[j]/trkpt[j]) >= trk_pt_resolution_cut) continue;
		 	if(fabs(trkdcaxy[j]/trkdcaxyerr[j]) >= trk_dca_xy_cut) continue;
		 	if(fabs(trkdcaz[j]/trkdcazerr[j]) >= trk_dca_z_cut) continue;
		 	double calomatching = ((pfEcal[j]+pfHcal[j])/cosh(trketa[j]))/trkpt[j];
		 	if(colliding_system == "PbPb" || colliding_system == "XeXe"){
		 		if((trkchi2[j]/trkndof[j])/trknlayer[j] >= chi2_ndf_nlayer_cut) continue;
		 		if(trknhits[j] < nhits) continue;
		 		if(trkpt[j] > 20.0 && fabs(calomatching) <= calo_matching) continue;
		 	}
		 	if(colliding_system=="PbPb" && sNN_energy_GeV==5020 && year_of_datataking==2018){if(trkalgo[j] == 6 && trkmva[j] < 0.98) continue;}

			// Track efficiency correction
		 	double trk_weight = 1.0;
		 	trk_weight = trk_weight*getTrkCorrWeight(fileeff, use_centrality, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, mult, trk_pt, trk_eta);

			trk_eta = trk_eta + boost; // In pPb case, for the center-of-mass correction if needed

			// Track QA histogram filling
			double x4D_reco_trk[4]={trk_pt,trk_eta,trk_phi,(double) multcentbin}; 
			hist_reco_trk->Fill(x4D_reco_trk);
			hist_reco_trk_corr->Fill(x4D_reco_trk,trk_weight);
			hist_reco_trk_weighted->Fill(x4D_reco_trk,trk_weight*event_weight);

			// Track vector filling
		 	TVector3 GoodTracks;
		 	GoodTracks.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
		 	tracks_reco.push_back(GoodTracks);
		 	sube_tracks_reco.push_back(0); // set == 0 because this is only valid for gen
		 	double trk_etamix_weight = get_trketamix_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, trk_eta, true); // weight to deal with Seagull (test)
		 	track_w_reco.push_back(trk_weight*trk_etamix_weight); // save weight to apply in the mixing

		 	int trackbin = (int) find_my_bin(trk_pt_bins, (float) trk_pt);
		 	
		 	if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
				//event plane information
				// Psi 2
				double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;	
				double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
				double x3D_reco_plus_dphi2_flat[3]={deltaphi2PC(trk_phi, Psi2_EP_flat_plus), (double) trackbin, (double) multcentbin}; 
				Dphi_EP2_flat_trk_plus->Fill(x3D_reco_plus_dphi2_flat,trk_weight*event_weight);
				double x3D_reco_minus_dphi2_flat[3]={deltaphi2PC(trk_phi, Psi2_EP_flat_minus), (double) trackbin, (double) multcentbin}; 
				Dphi_EP2_flat_trk_minus->Fill(x3D_reco_minus_dphi2_flat,trk_weight*event_weight);
			
				// Psi 3
				double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;	
				double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
				double x3D_reco_plus_dphi3_flat[3]={deltaphi2PC(trk_phi, Psi3_EP_flat_plus), (double) trackbin, (double) multcentbin}; 
				Dphi_EP3_flat_trk_plus->Fill(x3D_reco_plus_dphi3_flat,trk_weight*event_weight);
				double x3D_reco_minus_dphi3_flat[3]={deltaphi2PC(trk_phi, Psi3_EP_flat_minus), (double) trackbin, (double) multcentbin}; 
				Dphi_EP3_flat_trk_minus->Fill(x3D_reco_minus_dphi3_flat,trk_weight*event_weight);
				
				// Psi 4
				double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;	
				double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
				double x3D_reco_plus_dphi4_flat[3]={deltaphi2PC(trk_phi, Psi4_EP_flat_plus), (double) trackbin, (double) multcentbin}; 
				Dphi_EP4_flat_trk_plus->Fill(x3D_reco_plus_dphi4_flat,trk_weight*event_weight);
				double x3D_reco_minus_dphi4_flat[3]={deltaphi2PC(trk_phi, Psi4_EP_flat_minus), (double) trackbin, (double) multcentbin}; 
				Dphi_EP4_flat_trk_minus->Fill(x3D_reco_minus_dphi4_flat,trk_weight*event_weight);
			}
		} // End loop over tracks

		// Start loop over jets
		float leadrecojet_pt=-999, leadrecojet_eta=-999, leadrecojet_phi=-999; // leading jet quantities
		float sublrecojet_pt=-999, sublrecojet_eta=-999, sublrecojet_phi=-999; // subleading jet quantities
		float leadrefjet_pt=-999, leadrefjet_eta=-999, leadrefjet_phi=-999; // leading jet ref quantities
		float sublrefjet_pt=-999, sublrefjet_eta=-999, sublrefjet_phi=-999; // subleading jet ref quantities

		bool isjet40included = false;
		bool isjet50included = false;
		bool isjet60included = false;
		bool isjet80included = false;
		bool isjet100included = false;
		bool isjetincluded = false;
		bool outsideetarange = false;
		bool outsidegenetarange = false;
		
		int njets=0;
		int njetssub=0;
		int njetslead=0;
		
		int jetsize = (int)nref; // number of jets in an event
		for (int j = 0; j < jetsize; j++){

	        if(trackMax[j]/rawpt[j] < 0.01)continue; // Cut for jets with only very low pT particles
	        if(trackMax[j]/rawpt[j] > 0.98)continue; // Cut for jets where all the pT is taken by one track
			if(jteta[j] < -4.0 || jteta[j] > 4.0) continue; // no accept jets with |eta| > 4

			// Define jet kinematics
		 	float jet_rawpt = rawpt[j];
		 	float jet_eta = jteta[j];
		 	float jet_phi = jtphi[j];

		 	// Apply JEC
		 	JEC.SetJetPT(jet_rawpt); 
		 	JEC.SetJetEta(jet_eta); 
		 	JEC.SetJetPhi(jet_phi);
		 	float jet_pt_corr = JEC.GetCorrectedPT();
	 	
 			//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
			double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr); // Jet weight (specially for MC)
 		 	jet_weight = jet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr, do_jet_smearing, 0.663); // Jet smearing (For systematics)
			
			//leading and subleading
		 	find_leading_subleading(jet_pt_corr,jet_eta,jet_phi,leadrecojet_pt,leadrecojet_eta,leadrecojet_phi,sublrecojet_pt,sublrecojet_eta,sublrecojet_phi); // Find leading and subleading jets

		 	if(jet_eta > jet_eta_min_cut && jet_eta < jet_eta_max_cut){ // Jet eta cut
				
				njets=njets+1;

				jet_eta = jet_eta + boost; // In pPb case, for the center-of-mass correction if needed

		 		if(jet_pt_corr > 40 && jet_pt_corr < jet_pt_max_cut){isjet40included = true;}
		 		if(jet_pt_corr > 50 && jet_pt_corr < jet_pt_max_cut){isjet50included = true;}
		 		if(jet_pt_corr > 60 && jet_pt_corr < jet_pt_max_cut){isjet60included = true;}
		 		if(jet_pt_corr > 80 && jet_pt_corr < jet_pt_max_cut){isjet80included = true;}
		 		if(jet_pt_corr > 100 && jet_pt_corr < jet_pt_max_cut){isjet100included = true;}
		
		 		if(jet_pt_corr > subleading_pT_min) njetssub=njetssub+1;

		 		hist_reco_jet_weighted_nocut->Fill(jet_pt_corr,event_weight*jet_weight); // Fill histogram without any pt cut
		 		
		 		if(jet_pt_corr > jet_pt_min_cut && jet_pt_corr < jet_pt_max_cut){ // Jet pT cut

					njetslead=njetslead+1;
		 			isjetincluded = true;

		 		 	// Fill reco jet QA histograms
		 			double x4D_reco_jet[4]={jet_rawpt,jet_eta,jet_phi,(double) multcentbin}; 
		 			hist_reco_jet->Fill(x4D_reco_jet);
		 			hist_reco_jet_weighted->Fill(x4D_reco_jet,event_weight*jet_weight);
		 			double x4D_reco_jet_corr[4]={jet_pt_corr,jet_eta,jet_phi,(double) multcentbin}; 
		 			hist_reco_jet_corr->Fill(x4D_reco_jet_corr);
		 			hist_reco_jet_corr_weighted->Fill(x4D_reco_jet_corr,event_weight*jet_weight);
		 			// Fill reco jet vectors
		 			TVector3 GoodJets;
		 			GoodJets.SetPtEtaPhi(jet_pt_corr, jet_eta, jet_phi);
		 			jets_reco.push_back(GoodJets);
		 			jet_w_reco.push_back(jet_weight);
		 			
		 			if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
		 				// Psi 2
						double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;	
						double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
		 				Dphi_flat_EP2_inclusive_plus->Fill(deltaphi2PC(jet_phi, Psi2_EP_flat_plus),(double)multcentbin,event_weight*jet_weight);
		 				Dphi_flat_EP2_inclusive_minus->Fill(deltaphi2PC(jet_phi, Psi2_EP_flat_minus),(double)multcentbin,event_weight*jet_weight);

						// Psi 3
						double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;	
						double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
		 				Dphi_flat_EP3_inclusive_plus->Fill(deltaphi2PC(jet_phi, Psi3_EP_flat_plus),(double)multcentbin,event_weight*jet_weight);
		 				Dphi_flat_EP3_inclusive_minus->Fill(deltaphi2PC(jet_phi, Psi3_EP_flat_minus),(double)multcentbin,event_weight*jet_weight);

						// Psi 4
						double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;	
						double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
		 				Dphi_flat_EP4_inclusive_plus->Fill(deltaphi2PC(jet_phi, Psi4_EP_flat_plus),(double)multcentbin,event_weight*jet_weight);
		 				Dphi_flat_EP4_inclusive_minus->Fill(deltaphi2PC(jet_phi, Psi4_EP_flat_minus),(double)multcentbin,event_weight*jet_weight);
		 			}
		 		}
		 	}

			if(is_MC){ 

			 	float ref_pt = refpt[j];
			 	float ref_eta = refeta[j];
			 	float ref_phi = refphi[j];

		 		double refjet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, ref_pt); // Jet weight (specially for MC)

		 		if(ref_pt <= 0) continue;
			 	if(ref_eta < -4.0 || ref_eta > 4.0) continue; // max jet eta
		 		find_leading_subleading(ref_pt,ref_eta,ref_phi,leadrefjet_pt,leadrefjet_eta,leadrefjet_phi,sublrefjet_pt,sublrefjet_eta,sublrefjet_phi); // Find leading and subleading ref jets

		 		if(ref_eta <= jet_eta_min_cut || ref_eta >= jet_eta_max_cut)continue;
				hist_matched_jet_weighted_nocut->Fill(ref_pt,event_weight*refjet_weight); // Fill histogram without any pT cut
				
		 		if(jet_rawpt <= 0) continue;
		 		if(jet_pt_corr <= 0) continue;
				if(ref_pt > jet_pt_min_cut && ref_pt < jet_pt_max_cut){
					double x4D_match_jet[4]={ref_pt,ref_eta,ref_phi,(double) multcentbin}; 
					hist_matched_jet->Fill(x4D_match_jet);	
					hist_matched_jet_weighted->Fill(x4D_match_jet,event_weight*refjet_weight);	
				}

		 		if(jet_eta <= jet_eta_min_cut || jet_eta >= jet_eta_max_cut)continue;
				
				ref_eta = ref_eta + boost;  // In pPb case, for the center-of-mass correction if needed
		 		
		 		double JES_ratio_raw_vs_ref = jet_rawpt/ref_pt;
		 		double JES_ratio_reco_vs_ref = jet_pt_corr/ref_pt;
		 		double JER_ratio_raw_vs_ref = (jet_rawpt - ref_pt)/ref_pt;
		 		double JER_ratio_reco_vs_ref = (jet_pt_corr - ref_pt)/ref_pt;

		 		int refparton;
		 		if(fabs(refparton_flavor[j]) >= 1 && fabs(refparton_flavor[j]) <= 6){
		 			refparton = fabs(refparton_flavor[j]);
		 		}else if(fabs(refparton_flavor[j]) == 21){
		 			refparton = 7;
		 		}else{refparton = 0;}

				//double x5D_JES_ratio_raw_vs_ref[5]={JES_ratio_raw_vs_ref,ref_pt,ref_eta,(double)refparton,(double)multcentbin}; 
				double x5D_JES_ratio_reco_vs_ref[5]={JES_ratio_reco_vs_ref,ref_pt,ref_eta,(double)refparton,(double)multcentbin}; 
				//double x5D_JER_ratio_raw_vs_ref[5]={JER_ratio_raw_vs_ref,ref_pt,ref_eta,(double)refparton,(double)multcentbin}; 
				double x5D_JER_ratio_reco_vs_ref[5]={JER_ratio_reco_vs_ref,ref_pt,ref_eta,(double)refparton,(double)multcentbin}; 

				//hist_jes_raw->Fill(x5D_JES_ratio_raw_vs_ref);
				//hist_jes_raw_weighted->Fill(x5D_JES_ratio_raw_vs_ref,event_weight*refjet_weight*jet_weight);
				//hist_jes_reco->Fill(x5D_JES_ratio_reco_vs_ref);
				hist_jes_reco_weighted->Fill(x5D_JES_ratio_reco_vs_ref,event_weight*refjet_weight*jet_weight);
				//hist_jer_raw->Fill(x5D_JER_ratio_raw_vs_ref);
				//hist_jer_raw_weighted->Fill(x5D_JER_ratio_raw_vs_ref,event_weight*refjet_weight*jet_weight);
				//hist_jer_reco->Fill(x5D_JER_ratio_reco_vs_ref);
				hist_jer_reco_weighted->Fill(x5D_JER_ratio_reco_vs_ref,event_weight*refjet_weight*jet_weight);

		 		int refpartonfromB;
		 		if(fabs(refparton_flavorForB[j]) >= 1 && fabs(refparton_flavorForB[j]) <= 6){
		 			refpartonfromB = fabs(refparton_flavorForB[j]);
		 		}else if(fabs(refparton_flavorForB[j]) == 21){
		 			refpartonfromB = 7;
		 		}else{refpartonfromB = 0;}

				//double x5D_JES_ratio_raw_vs_reffromB[5]={JES_ratio_raw_vs_ref,ref_pt,ref_eta,(double)refpartonfromB,(double)multcentbin}; 
				double x5D_JES_ratio_reco_vs_reffromB[5]={JES_ratio_reco_vs_ref,ref_pt,ref_eta,(double)refpartonfromB,(double)multcentbin}; 
				//double x5D_JER_ratio_raw_vs_reffromB[5]={JER_ratio_raw_vs_ref,ref_pt,ref_eta,(double)refpartonfromB,(double)multcentbin}; 
				double x5D_JER_ratio_reco_vs_reffromB[5]={JER_ratio_reco_vs_ref,ref_pt,ref_eta,(double)refpartonfromB,(double)multcentbin}; 

				//hist_jes_raw_fromB->Fill(x5D_JES_ratio_raw_vs_reffromB);
				//hist_jes_raw_fromB_weighted->Fill(x5D_JES_ratio_raw_vs_reffromB,event_weight*refjet_weight*jet_weight);
				//hist_jes_reco_fromB->Fill(x5D_JES_ratio_reco_vs_reffromB);
				hist_jes_reco_fromB_weighted->Fill(x5D_JES_ratio_reco_vs_reffromB,event_weight*refjet_weight*jet_weight);
				//hist_jer_raw_fromB->Fill(x5D_JER_ratio_raw_vs_reffromB);
				//hist_jer_raw_fromB_weighted->Fill(x5D_JER_ratio_raw_vs_reffromB,event_weight*refjet_weight*jet_weight);
				//hist_jer_reco_fromB->Fill(x5D_JER_ratio_reco_vs_reffromB);
				hist_jer_reco_fromB_weighted->Fill(x5D_JER_ratio_reco_vs_reffromB,event_weight*refjet_weight*jet_weight);

				// Matched for JES studies and parton fraction
				// Fill jet pT for parton fractions --> no cuts applied
				double x4D_parton_ref[4]={ref_pt,ref_eta,(double)refparton,(double)multcentbin}; 
				double x4D_parton_ref_fromB[4]={ref_pt,ref_eta,(double)refpartonfromB,(double)multcentbin}; 
				hist_matched_jet_parton->Fill(x4D_parton_ref,event_weight*refjet_weight*jet_weight);
				hist_matched_jet_parton_fromB->Fill(x4D_parton_ref_fromB,event_weight*refjet_weight*jet_weight);
				
			} // End loop over matched (MC)

		} // End loop over jets

		if(isjet40included) multiplicity_withonejet40->Fill(mult,event_weight);				
		if(isjet50included) multiplicity_withonejet50->Fill(mult,event_weight);
		if(isjet60included) multiplicity_withonejet60->Fill(mult,event_weight);
		if(isjet80included) multiplicity_withonejet80->Fill(mult,event_weight);
		if(isjet100included) multiplicity_withonejet100->Fill(mult,event_weight);

		if(isjetincluded){
			multiplicity_withonejet->Fill(mult);
			multiplicity_withonejet_weighted->Fill(mult,event_weight);
			reco_mult_withonejet->Fill(recomult);
			reco_mult_withonejet_weighted->Fill(recomult, event_weight);
			if(colliding_system=="pPb" && year_of_datataking==2016){
				double x3D_hiHF_onejet[3]={hfplus,hfminus,(double) mult}; /*hfhist_onejet->Fill(x3D_hiHF_onejet);*/ hfhist_onejet_weighted->Fill(x3D_hiHF_onejet,event_weight);
				double x3D_hiHFEta4_onejet[3]={hfplusEta4,hfminusEta4,(double) mult}; /*hfhistEta4_onejet->Fill(x3D_hiHFEta4_onejet);*/ hfhistEta4_onejet_weighted->Fill(x3D_hiHFEta4_onejet,event_weight);
				double x3D_hiZDC_onejet[3]={zdcplus,zdcminus,(double) mult}; /*zdchist_onejet->Fill(x3D_hiZDC_onejet);*/ zdchist_onejet_weighted->Fill(x3D_hiZDC_onejet,event_weight);
			}
		}
		
		if(leadrecojet_eta < jet_eta_min_cut) outsideetarange = true;
		if(leadrecojet_eta > jet_eta_max_cut) outsideetarange = true;
		if(sublrecojet_eta < jet_eta_min_cut) outsideetarange = true;	
		if(sublrecojet_eta > jet_eta_max_cut) outsideetarange = true;	
		
		NJets->Fill(njets);
		NJetsSub->Fill(njetssub);
		NJetsLead->Fill(njetslead);
		if(leadrecojet_pt > leading_pT_min) NJetsLJSLJ->Fill(njetssub);
		
		bool isdijet = false;

		//leading/subleading jets
		if(jetsize > 1 && !outsideetarange){

			Nevents->Fill(6);

			double ljet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrecojet_pt);  // Jet weight (specially for MC)
			//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
	 		ljet_weight = ljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrecojet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

			double sljet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrecojet_pt);  // Jet weight (specially for MC)
			//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
	 		sljet_weight = sljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrecojet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

	 		// Fill leading and subleading pT histograms without cuts
			hist_reco_leadjet_pt_nocut->Fill(leadrecojet_pt);
			hist_reco_leadjet_pt_nocut_weighted->Fill(leadrecojet_pt, event_weight*ljet_weight);
			hist_reco_subljet_pt_nocut->Fill(sublrecojet_pt);
			hist_reco_subljet_pt_nocut_weighted->Fill(sublrecojet_pt, event_weight*sljet_weight);

			//leading/subleading pT cuts
			if(leadrecojet_pt > leading_pT_min && sublrecojet_pt > subleading_pT_min){
				
				Nevents->Fill(7);

				leadrecojet_eta = leadrecojet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				sublrecojet_eta = sublrecojet_eta + boost;  // In pPb case, for the center-of-mass correction if needed

				Nevents->Fill(8);	
				
				// Fill leading/subleading jet quenching quantities
				double delta_phi_reco = fabs(deltaphi(leadrecojet_phi, sublrecojet_phi));
				double delta_phi_reco_2pc = deltaphi2PC(leadrecojet_phi, sublrecojet_phi);
				double Aj_reco = asymmetry(leadrecojet_pt,sublrecojet_pt);
				double Xj_reco = xjvar(leadrecojet_pt,sublrecojet_pt);
				double x4D_reco[4]={Xj_reco,Aj_reco,delta_phi_reco,(double)multcentbin}; hist_reco_lead_reco_subl_quench->Fill(x4D_reco,event_weight*ljet_weight*sljet_weight);
				double x4D_reco2pc[4]={Xj_reco,Aj_reco,delta_phi_reco_2pc,(double)multcentbin}; hist_reco_lead_reco_subl_quench2pc->Fill(x4D_reco2pc,event_weight*ljet_weight*sljet_weight);
				double etadijet = (leadrecojet_eta + sublrecojet_eta)/2.0;
				double etadiff = deltaeta(leadrecojet_eta,sublrecojet_eta);
				double x5D_etaassym_ETEta4[5] = {etadijet,Xj_reco,delta_phi_reco,(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
				hist_etaDijet_ETEta4_reco->Fill(x5D_etaassym_ETEta4,event_weight*ljet_weight*sljet_weight);
				double x5D_etadiff_ETEta4[5] = {etadiff,Xj_reco,delta_phi_reco,(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
				hist_etaDiff_ETEta4_reco->Fill(x5D_etadiff_ETEta4,event_weight*ljet_weight*sljet_weight);

				/*
				double x5D_etaassym_ET[5] = {etadijet,Xj_reco,delta_phi_reco,(double)(hfplus+hfminus),(double)multcentbin};
				hist_etaDijet_ET_reco->Fill(x5D_etaassym_ET,event_weight*ljet_weight*sljet_weight);
				double x5D_etadiff_ET[5] = {etadiff,Xj_reco,delta_phi_reco,(double)(hfplus+hfminus),(double)multcentbin};
				hist_etaDiff_ET_reco->Fill(x5D_etadiff_ET,event_weight*ljet_weight*sljet_weight);
				*/

				// leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
				if(delta_phi_reco > leading_subleading_deltaphi_min){

					double x5D_etaAssym[5] = {0,Xj_reco,(double)(hfplus+hfminus),(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
					if(etadijet > boost){hist_etaAssym_numerator_reco->Fill(x5D_etaAssym,event_weight*ljet_weight*sljet_weight);}else if(etadijet < boost){hist_etaAssym_denominator_reco->Fill(x5D_etaAssym,event_weight*ljet_weight*sljet_weight);}

					// Fill leading and subleading jet vectors
					TVector3 GoodLeadingJets_reco;
					GoodLeadingJets_reco.SetPtEtaPhi(leadrecojet_pt, leadrecojet_eta, leadrecojet_phi);
					lead_jets_reco.push_back(GoodLeadingJets_reco);
					lead_jet_w_reco.push_back(ljet_weight);
					TVector3 GoodSubLeadingJets_reco;
					GoodSubLeadingJets_reco.SetPtEtaPhi(sublrecojet_pt, sublrecojet_eta, sublrecojet_phi);
					subl_jets_reco.push_back(GoodSubLeadingJets_reco);
					subl_jet_w_reco.push_back(sljet_weight);

					if((Xj_reco >= xjmin && Xj_reco <= xjmax) && (Aj_reco >= Ajmin && Aj_reco <= Ajmax)){

						Nevents->Fill(9);
						pass_Aj_or_Xj_reco_cut = true; // if we apply Xj or Aj cuts
						isdijet = true;

						// Fill leading and subleading jet QA histograms
						double x4D_lead[4]={leadrecojet_pt,leadrecojet_eta,leadrecojet_phi,(double) multcentbin}; 
						hist_reco_leadjet->Fill(x4D_lead);	
						hist_reco_leadjet_weighted->Fill(x4D_lead,event_weight*ljet_weight);
						double x4D_sublead[4]={sublrecojet_pt,sublrecojet_eta,sublrecojet_phi,(double) multcentbin}; 
						hist_reco_subljet->Fill(x4D_sublead);	
						hist_reco_subljet_weighted->Fill(x4D_sublead,event_weight*sljet_weight);

		 				if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
		 					// Psi 2
							double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;	
							double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
		 					Dphi_flat_EP2_leading_plus->Fill(deltaphi2PC(leadrecojet_phi, Psi2_EP_flat_plus),(double)multcentbin,event_weight*ljet_weight);
		 					Dphi_flat_EP2_leading_minus->Fill(deltaphi2PC(leadrecojet_phi, Psi2_EP_flat_minus),(double)multcentbin,event_weight*ljet_weight);
		 					Dphi_flat_EP2_subleading_plus->Fill(deltaphi2PC(sublrecojet_phi, Psi2_EP_flat_plus),(double)multcentbin,event_weight*sljet_weight);
		 					Dphi_flat_EP2_subleading_minus->Fill(deltaphi2PC(sublrecojet_phi, Psi2_EP_flat_minus),(double)multcentbin,event_weight*sljet_weight);

							// Psi 3
							double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;	
							double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
		 					Dphi_flat_EP3_leading_plus->Fill(deltaphi2PC(leadrecojet_phi, Psi3_EP_flat_plus),(double)multcentbin,event_weight*ljet_weight);
		 					Dphi_flat_EP3_leading_minus->Fill(deltaphi2PC(leadrecojet_phi, Psi3_EP_flat_minus),(double)multcentbin,event_weight*ljet_weight);
		 					Dphi_flat_EP3_subleading_plus->Fill(deltaphi2PC(sublrecojet_phi, Psi3_EP_flat_plus),(double)multcentbin,event_weight*sljet_weight);
		 					Dphi_flat_EP3_subleading_minus->Fill(deltaphi2PC(sublrecojet_phi, Psi3_EP_flat_minus),(double)multcentbin,event_weight*sljet_weight);

							// Psi 4
							double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;	
							double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
		 					Dphi_flat_EP4_leading_plus->Fill(deltaphi2PC(leadrecojet_phi, Psi4_EP_flat_plus),(double)multcentbin,event_weight*ljet_weight);
		 					Dphi_flat_EP4_leading_minus->Fill(deltaphi2PC(leadrecojet_phi, Psi4_EP_flat_minus),(double)multcentbin,event_weight*ljet_weight);
		 					Dphi_flat_EP4_subleading_plus->Fill(deltaphi2PC(sublrecojet_phi, Psi4_EP_flat_plus),(double)multcentbin,event_weight*sljet_weight);
		 					Dphi_flat_EP4_subleading_minus->Fill(deltaphi2PC(sublrecojet_phi, Psi4_EP_flat_minus),(double)multcentbin,event_weight*sljet_weight);
		 					
		 				}
					}
				}
			}
		}
		
		if(jetsize > 1 && !outsideetarange){
			//leading/subleading pT cuts
			if(is_MC && leadrefjet_pt > leading_pT_min && sublrefjet_pt > subleading_pT_min){
			
				leadrefjet_eta = leadrefjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				sublrefjet_eta = sublrefjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed

				double lrefjet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrefjet_pt);  // Jet weight (specially for MC)
				double slrefjet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrefjet_pt);  // Jet weight (specially for MC)
				// Fill leading/subleading jet quenching quantities
				double delta_phi_ref = fabs(deltaphi(leadrefjet_phi, sublrefjet_phi));
				double delta_phi_ref_2pc = deltaphi2PC(leadrefjet_phi, sublrefjet_phi);
				double Aj_ref = asymmetry(leadrefjet_pt,sublrefjet_pt);
				double Xj_ref = xjvar(leadrefjet_pt,sublrefjet_pt);
				double x4D_ref[4]={Xj_ref,Aj_ref,delta_phi_ref,(double)multcentbin}; hist_ref_lead_ref_subl_quench->Fill(x4D_ref,event_weight*lrefjet_weight*slrefjet_weight);
				double x4D_ref2pc[4]={Xj_ref,Aj_ref,delta_phi_ref_2pc,(double)multcentbin}; hist_ref_lead_ref_subl_quench2pc->Fill(x4D_ref2pc,event_weight*lrefjet_weight*slrefjet_weight);

				double etadijet = (leadrefjet_eta + sublrefjet_eta)/2.0;
				double etadiff = deltaeta(leadrefjet_eta,sublrefjet_eta);
				double x5D_etaassym_ETEta4[5] = {etadijet,Xj_ref,delta_phi_ref,(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
				hist_etaDijet_ETEta4_ref->Fill(x5D_etaassym_ETEta4,event_weight*lrefjet_weight*slrefjet_weight);
				double x5D_etadiff_ETEta4[5] = {etadiff,Xj_ref,delta_phi_ref,(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
				hist_etaDiff_ETEta4_ref->Fill(x5D_etadiff_ETEta4,event_weight*lrefjet_weight*slrefjet_weight);
				/*
				double x5D_etaassym_ET[5] = {etadijet,Xj_ref,delta_phi_ref,(double)(hfplus+hfminus),(double)multcentbin};
				hist_etaDijet_ET_ref->Fill(x5D_etaassym_ET,event_weight*lrefjet_weight*slrefjet_weight);
				double x5D_etadiff_ET[5] = {etadiff,Xj_ref,delta_phi_ref,(double)(hfplus+hfminus),(double)multcentbin};
				hist_etaDiff_ET_ref->Fill(x5D_etadiff_ET,event_weight*lrefjet_weight*slrefjet_weight);
				*/
				// leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
				if(delta_phi_ref > leading_subleading_deltaphi_min){
					double x5D_etaAssym[5] = {0,Xj_ref,(double)(hfplus+hfminus),(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
					if(etadijet > boost){hist_etaAssym_numerator_ref->Fill(x5D_etaAssym,event_weight*lrefjet_weight*slrefjet_weight);}else if(etadijet < boost){hist_etaAssym_denominator_ref->Fill(x5D_etaAssym,event_weight*lrefjet_weight*slrefjet_weight);}
				}
			}
		}
		
		if(isdijet){
			multiplicity_withdijets->Fill(mult);
			multiplicity_withdijets_weighted->Fill(mult,event_weight);
			reco_mult_withdijets->Fill(recomult);
			reco_mult_withdijets_weighted->Fill(recomult, event_weight);
			if(colliding_system=="pPb" && year_of_datataking==2016){
				double x3D_hiHF_dijet[3]={hfplus,hfminus,(double) mult}; /*hfhist_dijet->Fill(x3D_hiHF_dijet);*/ hfhist_dijet_weighted->Fill(x3D_hiHF_dijet,event_weight);
				double x3D_hiHFEta4_dijet[3]={hfplusEta4,hfminusEta4,(double) mult};/*hfhistEta4_dijet->Fill(x3D_hiHFEta4_dijet);*/ hfhistEta4_dijet_weighted->Fill(x3D_hiHFEta4_dijet,event_weight);
				double x3D_hiZDC_dijet[3]={zdcplus,zdcminus,(double) mult}; /*zdchist_dijet->Fill(x3D_hiZDC_dijet);*/ zdchist_dijet_weighted->Fill(x3D_hiZDC_dijet,event_weight);
			}
		}

		// Measure correlations and filling mixing vectors
		// Reco-Reco
		// Inclusive jets
		if(do_inclusejettrack_correlation && pass_Aj_or_Xj_reco_cut){
			correlation(jets_reco, jet_w_reco, tracks_reco, track_w_reco, hist_correlation_signal_jet_reco_track_reco, hist_jet_from_reco_reco_sig, hist_trk_from_reco_reco_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_jet_reco_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_jet_reco_track_reco,JetR,hist_injet_reco_track_reco,do_flow); // calculate correlations
			fillvectors(similar_events, Nev_recoreco, jets_reco, jet_w_reco, tracks_reco, track_w_reco, mult, vertexz, event_weight, ev_jet_vector_reco_reco, jet_weights_reco_reco, ev_track_vector_reco_reco, trk_weights_reco_reco, multvec_reco_reco, vzvec_reco_reco, weights_reco_reco);	// for mixing --> store vectors for mixing
		}

		// Leading/SubLeading jets
		if(do_leading_subleading_jettrack_correlation && pass_Aj_or_Xj_reco_cut){
			correlation(lead_jets_reco, lead_jet_w_reco, tracks_reco, track_w_reco, hist_correlation_signal_lead_jet_reco_track_reco, hist_lead_jet_from_reco_reco_sig, hist_LJ_trk_from_reco_reco_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_reco_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_lead_jet_reco_track_reco,JetR,hist_inLeadjet_reco_track_reco,do_flow); // calculate correlations
			correlation(subl_jets_reco, subl_jet_w_reco, tracks_reco, track_w_reco, hist_correlation_signal_subl_jet_reco_track_reco, hist_subl_jet_from_reco_reco_sig, hist_SLJ_trk_from_reco_reco_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_reco_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_subl_jet_reco_track_reco,JetR,hist_inSubljet_reco_track_reco,do_flow); // calculate correlations
			fillvectors(similar_events, Nev_recoreco_lead, lead_jets_reco, lead_jet_w_reco, tracks_reco, track_w_reco, mult, vertexz, event_weight, ev_jet_vector_leadjet_reco_reco, jet_weights_leadjet_reco_reco, ev_track_vector_leadjet_reco_reco, trk_weights_leadjet_reco_reco, multvec_leadjet_reco_reco, vzvec_leadjet_reco_reco, weights_leadjet_reco_reco);	// for mixing --> store vectors for mixing
			fillvectors(similar_events, Nev_recoreco_subl, subl_jets_reco, subl_jet_w_reco, tracks_reco, track_w_reco, mult, vertexz, event_weight, ev_jet_vector_subleadjet_reco_reco, jet_weights_subleadjet_reco_reco, ev_track_vector_subleadjet_reco_reco, trk_weights_subleadjet_reco_reco, multvec_subleadjet_reco_reco, vzvec_subleadjet_reco_reco, weights_subleadjet_reco_reco); // for mixing --> store vectors for mixing
		}

		if(do_flow){twoparticlecorrelation(tracks_reco, track_w_reco, hist_reco_reco_2pcorrelation_signal, event_weight, mult, sube_tracks_reco, hist_reco_reco_2pcorrelation_signal_subg0, hist_reco_reco_2pcorrelation_signal_subcross);} // calculate 2 particle correlations

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

				// Track/particle QA histogram filling
				double x4D_gen_trk[4]={gtrk_pt,gtrk_eta,gtrk_phi,(double) multcentbin}; 
				hist_gen_trk->Fill(x4D_gen_trk);
				hist_gen_trk_weighted->Fill(x4D_gen_trk,event_weight);

				// Track/particle vector filling
				TVector3 GoodTracks_gen;
				GoodTracks_gen.SetPtEtaPhi(gtrk_pt, gtrk_eta, gtrk_phi);
				tracks_gen.push_back(GoodTracks_gen);
		 		sube_tracks_gen.push_back(gen_trksube->at(j)); // get sube from the tree
		 		double trk_etamix_weight = get_trketamix_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gtrk_eta, false); // weight to deal with Seagull (test)
		 		track_w_gen.push_back(trk_etamix_weight); // save weight to apply in the mixing
		 		
		 		int trackbin = (int) find_my_bin(trk_pt_bins, (float) gtrk_pt);
		 		
		 		if(do_flow){
					//event plane information
					// Psi 2
					double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;	
					double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
					double x3D_gen_plus_dphi2_flat[3]={deltaphi2PC(gtrk_phi, Psi2_EP_flat_plus), (double) trackbin, (double) multcentbin}; 
					Dphi_GEN_EP2_flat_trk_plus->Fill(x3D_gen_plus_dphi2_flat,event_weight);
					double x3D_gen_minus_dphi2_flat[3]={deltaphi2PC(gtrk_phi, Psi2_EP_flat_minus), (double) trackbin, (double) multcentbin}; 
					Dphi_GEN_EP2_flat_trk_minus->Fill(x3D_gen_minus_dphi2_flat,event_weight);
								
					// Psi 3
					double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;	
					double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
					double x3D_gen_plus_dphi3_flat[3]={deltaphi2PC(gtrk_phi, Psi3_EP_flat_plus), (double) trackbin, (double) multcentbin}; 
					Dphi_GEN_EP3_flat_trk_plus->Fill(x3D_gen_plus_dphi3_flat,event_weight);
					double x3D_gen_minus_dphi3_flat[3]={deltaphi2PC(gtrk_phi, Psi3_EP_flat_minus), (double) trackbin, (double) multcentbin}; 
					Dphi_GEN_EP3_flat_trk_minus->Fill(x3D_gen_minus_dphi3_flat,event_weight);
			
					// Psi 4
					double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;	
					double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
					double x3D_gen_plus_dphi4_flat[3]={deltaphi2PC(gtrk_phi, Psi4_EP_flat_plus), (double) trackbin, (double) multcentbin}; 
					Dphi_GEN_EP4_flat_trk_plus->Fill(x3D_gen_plus_dphi4_flat,event_weight);
					double x3D_gen_minus_dphi4_flat[3]={deltaphi2PC(gtrk_phi, Psi4_EP_flat_minus), (double) trackbin, (double) multcentbin}; 
					Dphi_GEN_EP4_flat_trk_minus->Fill(x3D_gen_minus_dphi4_flat,event_weight);
				}
			}

			// Start loop over gen jets
			float leadgenjet_pt=-999, leadgenjet_eta=-999, leadgenjet_phi=-999; // leading jet quantities
			float sublgenjet_pt=-999, sublgenjet_eta=-999, sublgenjet_phi=-999; // subleading jet quantities

			bool isgjetincluded = false;

			for(int j = 0; j < gen_jetsize; j++){
				
				// Define jet kinematics
				float gjet_pt = gen_jtpt[j];
				float gjet_eta = gen_jteta[j];
				float gjet_phi = gen_jtphi[j];
				
		 		if(gjet_eta < -4.0 || gjet_eta > 4.0) continue; // no accept jets with |eta| > 4

				find_leading_subleading(gjet_pt,gjet_eta,gjet_phi,leadgenjet_pt,leadgenjet_eta,leadgenjet_phi,sublgenjet_pt,sublgenjet_eta,sublgenjet_phi); // Find leading and subleading jets

		 		if(gjet_eta < jet_eta_min_cut || gjet_eta > jet_eta_max_cut) continue; // jet eta cut

				gjet_eta = gjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed

				hist_gen_jet_weighted_nocut->Fill(gjet_pt,event_weight); // Fill jet pT without cut

				if(gjet_pt > jet_pt_min_cut && gjet_pt < jet_pt_max_cut){  // Jet pT cut
				
					isgjetincluded=true;
					
					double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gjet_pt); // Jet weight (specially for MC)
					//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
			 		jet_weight = jet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gjet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

			 		// Fill gen jet QA histograms
					double x4D_gen_jet[4]={gjet_pt,gjet_eta,gjet_phi,(double) multcentbin}; 
					hist_gen_jet->Fill(x4D_gen_jet);
					hist_gen_jet_weighted->Fill(x4D_gen_jet,event_weight*jet_weight);

					// Fill gen jet vectors
					TVector3 GoodJets_gen;
					GoodJets_gen.SetPtEtaPhi(gjet_pt, gjet_eta, gjet_phi);
					jets_gen.push_back(GoodJets_gen);
		 			jet_w_gen.push_back(jet_weight);
		 			/*
					for (int k = 0; k < jetsize; k++){ // for delta R different axis calculation (near future, new analysis/studies)

						if(refpt[k] != gen_jtpt[j])continue;
						if(use_WTA){
							if(refeta[k] != gen_jteta_otheraxis[j])continue;
							if(refphi[k] != gen_jtphi_otheraxis[j])continue;
						}else{
							if(refeta[k] != gen_jteta[j])continue;
							if(refphi[k] != gen_jtphi[j])continue;
						}

						double Delta_R = deltaR(gen_jteta[j], gen_jtphi[j], gen_jteta_otheraxis[j], gen_jtphi_otheraxis[j]);

		 				int refpartonfromB;
		 				if(fabs(refparton_flavorForB[k]) >= 1 && fabs(refparton_flavorForB[k]) <= 6){
		 					refpartonfromB = fabs(refparton_flavorForB[k]);
		 				}else if(fabs(refparton_flavorForB[k]) == 21){
		 					refpartonfromB = 7;
		 				}else{refpartonfromB = 0;}

						double x4D_gen_jetaxischeck[4]={Delta_R,refpt[k],(double)refpartonfromB,(double) multcentbin}; 
						genjetaxischeck->Fill(x4D_gen_jetaxischeck,event_weight*jet_weight);

					}
		 			*/
		 			if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
		 				// Psi 2
						double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;	
						double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
		 				Dphi_GEN_flat_EP2_inclusive_plus->Fill(deltaphi2PC(gjet_phi, Psi2_EP_flat_plus),(double)multcentbin,event_weight*jet_weight);
		 				Dphi_GEN_flat_EP2_inclusive_minus->Fill(deltaphi2PC(gjet_phi, Psi2_EP_flat_minus),(double)multcentbin,event_weight*jet_weight);

						// Psi 3
						double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;	
						double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
		 				Dphi_GEN_flat_EP3_inclusive_plus->Fill(deltaphi2PC(gjet_phi, Psi3_EP_flat_plus),(double)multcentbin,event_weight*jet_weight);
		 				Dphi_GEN_flat_EP3_inclusive_minus->Fill(deltaphi2PC(gjet_phi, Psi3_EP_flat_minus),(double)multcentbin,event_weight*jet_weight);

						// Psi 4
						double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;	
						double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
		 				Dphi_GEN_flat_EP4_inclusive_plus->Fill(deltaphi2PC(gjet_phi, Psi4_EP_flat_plus),(double)multcentbin,event_weight*jet_weight);
		 				Dphi_GEN_flat_EP4_inclusive_minus->Fill(deltaphi2PC(gjet_phi, Psi4_EP_flat_minus),(double)multcentbin,event_weight*jet_weight);
		 			}
				}
			}
			
			if(isgjetincluded){gen_mult_withonejet->Fill(genmult); gen_mult_withonejet_weighted->Fill(genmult, event_weight);}

			if(leadgenjet_eta < jet_eta_min_cut) outsidegenetarange = true;
			if(leadgenjet_eta > jet_eta_max_cut) outsidegenetarange = true;
			if(sublgenjet_eta < jet_eta_min_cut) outsidegenetarange = true;	
			if(sublgenjet_eta > jet_eta_max_cut) outsidegenetarange = true;	
		
			bool isgdijet = false;
			
			//leading/subleading jets
			if(gen_jetsize > 1 && !outsidegenetarange){

				double ljet_weight = get_leadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadgenjet_pt); // Jet weight (specially for MC)
				//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
		 		ljet_weight = ljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadgenjet_pt, do_jet_smearing, 0.663); // Jet smearing (For systematics)

				double sljet_weight = get_subleadjetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublgenjet_pt); // Jet weight (specially for MC)
				//resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
		 		sljet_weight = sljet_weight*get_jetpTsmering_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublgenjet_pt, do_jet_smearing, 0.663);  // Jet smearing (For systematics)
				
				// Fill leading and subleading pT histograms without cuts
				hist_gen_leadjet_pt_nocut->Fill(leadgenjet_pt);
				hist_gen_leadjet_pt_nocut_weighted->Fill(leadgenjet_pt, event_weight*ljet_weight);
				hist_gen_subljet_pt_nocut->Fill(sublgenjet_pt);
				hist_gen_subljet_pt_nocut_weighted->Fill(sublgenjet_pt, event_weight*sljet_weight);

				//leading/subleading pT cuts
				if(leadgenjet_pt > leading_pT_min && sublgenjet_pt > subleading_pT_min){ 
				
					leadgenjet_eta = leadgenjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
					sublgenjet_eta = sublgenjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed

					// Fill leading/subleading jet quenching quantities
					double delta_phi_gen = fabs(deltaphi(leadgenjet_phi, sublgenjet_phi));
					double delta_phi_gen_2pc = deltaphi2PC(leadgenjet_phi, sublgenjet_phi);
					double Aj_gen = asymmetry(leadgenjet_pt,sublgenjet_pt);
					double Xj_gen = xjvar(leadgenjet_pt,sublgenjet_pt);
					double x4D_gen[4]={Xj_gen,Aj_gen,delta_phi_gen,(double)multcentbin}; hist_gen_lead_gen_subl_quench->Fill(x4D_gen,event_weight*ljet_weight*sljet_weight);
					double x4D_gen2pc[4]={Xj_gen,Aj_gen,delta_phi_gen_2pc,(double)multcentbin}; hist_gen_lead_gen_subl_quench2pc->Fill(x4D_gen2pc,event_weight*ljet_weight*sljet_weight);

					double etadijet = (leadgenjet_eta + sublgenjet_eta)/2.0;
					double etadiff = deltaeta(leadgenjet_eta,sublgenjet_eta);
					double x5D_etaassym_ETEta4[5] = {etadijet,Xj_gen,delta_phi_gen,(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
					hist_etaDijet_ETEta4_gen->Fill(x5D_etaassym_ETEta4,event_weight*ljet_weight*sljet_weight);
					double x5D_etadiff_ETEta4[5] = {etadiff,Xj_gen,delta_phi_gen,(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
					hist_etaDiff_ETEta4_gen->Fill(x5D_etadiff_ETEta4,event_weight*ljet_weight*sljet_weight);
					/*
					double x5D_etaassym_ET[5] = {etadijet,Xj_gen,delta_phi_gen,(double)(hfplus+hfminus),(double)multcentbin};
					hist_etaDijet_ET_gen->Fill(x5D_etaassym_ET,event_weight*ljet_weight*sljet_weight);
					double x5D_etadiff_ET[5] = {etadiff,Xj_gen,delta_phi_gen,(double)(hfplus+hfminus),(double)multcentbin};
					hist_etaDiff_ET_gen->Fill(x5D_etadiff_ET,event_weight*ljet_weight*sljet_weight);
					*/
					// leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations				
					if(delta_phi_gen > leading_subleading_deltaphi_min){

						double x5D_etaAssym[5] = {0,Xj_gen,(double)(hfplus+hfminus),(double)(hfplusEta4+hfminusEta4),(double)multcentbin};
						if(etadijet > boost){hist_etaAssym_numerator_gen->Fill(x5D_etaAssym);}else if(etadijet < boost){hist_etaAssym_denominator_gen->Fill(x5D_etaAssym);}
			
						// Fill leading and subleading jet vectors
						TVector3 GoodLeadingJets_gen;
						GoodLeadingJets_gen.SetPtEtaPhi(leadgenjet_pt, leadgenjet_eta, leadgenjet_phi);
						lead_jets_gen.push_back(GoodLeadingJets_gen);
						lead_jet_w_gen.push_back(ljet_weight);
						TVector3 GoodSubLeadingJets_gen;
						GoodSubLeadingJets_gen.SetPtEtaPhi(sublgenjet_pt, sublgenjet_eta, sublgenjet_phi);
						subl_jets_gen.push_back(GoodSubLeadingJets_gen);
						subl_jet_w_gen.push_back(sljet_weight);

						if((Xj_gen >= xjmin && Xj_gen < xjmax) && (Aj_gen >= Ajmin && Aj_gen < Ajmax)){

							// Fill leading and subleading jet QA histograms
							double x4D_lead[4]={leadgenjet_pt,leadgenjet_eta,leadgenjet_phi,(double) multcentbin}; 
							hist_gen_leadjet->Fill(x4D_lead);	
							hist_gen_leadjet_weighted->Fill(x4D_lead,event_weight*ljet_weight);
							double x4D_sublead[4]={sublgenjet_pt,sublgenjet_eta,sublgenjet_phi,(double) multcentbin}; 
							hist_gen_subljet->Fill(x4D_sublead);	
							hist_gen_subljet_weighted->Fill(x4D_sublead,event_weight*sljet_weight);

							isgdijet = true;
							pass_Aj_or_Xj_gen_cut = true; // if we apply Xj or Aj cuts
							
		 					if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
		 					
		 						// Psi 2
								double Psi2_EP_flat_plus = (double) EP_Psi2_plus_flat;	
								double Psi2_EP_flat_minus = (double)EP_Psi2_minus_flat;
		 						Dphi_GEN_flat_EP2_leading_plus->Fill(deltaphi2PC(leadrecojet_phi, Psi2_EP_flat_plus),(double)multcentbin,event_weight*ljet_weight);
		 						Dphi_GEN_flat_EP2_leading_minus->Fill(deltaphi2PC(leadrecojet_phi, Psi2_EP_flat_minus),(double)multcentbin,event_weight*ljet_weight);
		 						Dphi_GEN_flat_EP2_subleading_plus->Fill(deltaphi2PC(sublrecojet_phi, Psi2_EP_flat_plus),(double)multcentbin,event_weight*sljet_weight);
		 						Dphi_GEN_flat_EP2_subleading_minus->Fill(deltaphi2PC(sublrecojet_phi, Psi2_EP_flat_minus),(double)multcentbin,event_weight*sljet_weight);

								// Psi 3
								double Psi3_EP_flat_plus = (double) EP_Psi3_plus_flat;	
								double Psi3_EP_flat_minus = (double)EP_Psi3_minus_flat;
		 						Dphi_GEN_flat_EP3_leading_plus->Fill(deltaphi2PC(leadrecojet_phi, Psi3_EP_flat_plus),(double)multcentbin,event_weight*ljet_weight);
		 						Dphi_GEN_flat_EP3_leading_minus->Fill(deltaphi2PC(leadrecojet_phi, Psi3_EP_flat_minus),(double)multcentbin,event_weight*ljet_weight);
		 						Dphi_GEN_flat_EP3_subleading_plus->Fill(deltaphi2PC(sublrecojet_phi, Psi3_EP_flat_plus),(double)multcentbin,event_weight*sljet_weight);
		 						Dphi_GEN_flat_EP3_subleading_minus->Fill(deltaphi2PC(sublrecojet_phi, Psi3_EP_flat_minus),(double)multcentbin,event_weight*sljet_weight);

								// Psi 4
								double Psi4_EP_flat_plus = (double) EP_Psi4_plus_flat;	
								double Psi4_EP_flat_minus = (double)EP_Psi4_minus_flat;
		 						Dphi_GEN_flat_EP4_leading_plus->Fill(deltaphi2PC(leadrecojet_phi, Psi4_EP_flat_plus),(double)multcentbin,event_weight*ljet_weight);
		 						Dphi_GEN_flat_EP4_leading_minus->Fill(deltaphi2PC(leadrecojet_phi, Psi4_EP_flat_minus),(double)multcentbin,event_weight*ljet_weight);
		 						Dphi_GEN_flat_EP4_subleading_plus->Fill(deltaphi2PC(sublrecojet_phi, Psi4_EP_flat_plus),(double)multcentbin,event_weight*sljet_weight);
		 						Dphi_GEN_flat_EP4_subleading_minus->Fill(deltaphi2PC(sublrecojet_phi, Psi4_EP_flat_minus),(double)multcentbin,event_weight*sljet_weight);
		 					
		 					}
						}
					}
				}
			}
			
			if(isgdijet){gen_mult_withdijets->Fill(genmult); gen_mult_withdijets_weighted->Fill(genmult, event_weight);}

			// Measure correlations and fill mixing vectors for inclusive jet+track correlations
			if(do_inclusejettrack_correlation){
				// Reco-Gen
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(jets_reco, jet_w_reco, tracks_gen, track_w_gen, hist_correlation_signal_jet_reco_track_gen, hist_jet_from_reco_gen_sig, hist_trk_from_reco_gen_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_jet_reco_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_jet_reco_track_gen,JetR,hist_injet_reco_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_recogen, jets_reco, jet_w_reco, tracks_gen, track_w_gen, mult, vertexz, event_weight, ev_jet_vector_reco_gen, jet_weights_reco_gen, ev_track_vector_reco_gen, trk_weights_reco_gen, multvec_reco_gen, vzvec_reco_gen, weights_reco_gen);	// for mixing --> store vectors for mixing
				// Gen-Reco			
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(jets_gen, jet_w_gen, tracks_reco, track_w_reco, hist_correlation_signal_jet_gen_track_reco, hist_jet_from_gen_reco_sig, hist_trk_from_gen_reco_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_jet_gen_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_jet_gen_track_reco,JetR,hist_injet_gen_track_reco,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_genreco, jets_gen, jet_w_gen, tracks_reco, track_w_reco, mult, vertexz, event_weight, ev_jet_vector_gen_reco, jet_weights_gen_reco, ev_track_vector_gen_reco, trk_weights_gen_reco, multvec_gen_reco, vzvec_gen_reco, weights_gen_reco);	// for mixing --> store vectors for mixing
				// Gen-Gen
				if(pass_Aj_or_Xj_gen_cut) correlation(jets_gen, jet_w_gen, tracks_gen, track_w_gen, hist_correlation_signal_jet_gen_track_gen, hist_jet_from_gen_gen_sig, hist_trk_from_gen_gen_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_jet_gen_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_jet_gen_track_gen,JetR,hist_injet_gen_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_gengen, jets_gen, jet_w_gen, tracks_gen, track_w_gen, mult, vertexz, event_weight, ev_jet_vector_gen_gen, jet_weights_gen_gen, ev_track_vector_gen_gen, trk_weights_gen_gen, multvec_gen_gen, vzvec_gen_gen, weights_gen_gen);	// for mixing --> store vectors for mixing
			}
			// Measure correlations and fill mixing vectors for (leading/subleading) jet+track correlations
			if(do_leading_subleading_jettrack_correlation){
				// Reco-Gen
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(lead_jets_reco, lead_jet_w_reco, tracks_gen, track_w_gen, hist_correlation_signal_lead_jet_reco_track_gen, hist_lead_jet_from_reco_gen_sig, hist_LJ_trk_from_reco_gen_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_reco_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_lead_jet_reco_track_gen,JetR,hist_inLeadjet_reco_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(subl_jets_reco, subl_jet_w_reco, tracks_gen, track_w_gen, hist_correlation_signal_subl_jet_reco_track_gen, hist_subl_jet_from_reco_gen_sig, hist_SLJ_trk_from_reco_gen_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_reco_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_subl_jet_reco_track_gen,JetR,hist_inSubljet_reco_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_recogen_lead, lead_jets_reco, lead_jet_w_reco, tracks_gen, track_w_gen, mult, vertexz, event_weight, ev_jet_vector_leadjet_reco_gen, jet_weights_leadjet_reco_gen, ev_track_vector_leadjet_reco_gen, trk_weights_leadjet_reco_gen, multvec_leadjet_reco_gen, vzvec_leadjet_reco_gen, weights_leadjet_reco_gen);	// for mixing --> store vectors for mixing
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_recogen_subl, subl_jets_reco, subl_jet_w_reco, tracks_gen, track_w_gen, mult, vertexz, event_weight, ev_jet_vector_subleadjet_reco_gen, jet_weights_subleadjet_reco_gen, ev_track_vector_subleadjet_reco_gen, trk_weights_subleadjet_reco_gen, multvec_subleadjet_reco_gen, vzvec_subleadjet_reco_gen, weights_subleadjet_reco_gen);	// for mixing --> store vectors for mixing
				// Gen-Reco	
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(lead_jets_gen, lead_jet_w_gen, tracks_reco, track_w_reco,hist_correlation_signal_lead_jet_gen_track_reco, hist_lead_jet_from_gen_reco_sig, hist_LJ_trk_from_gen_reco_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_gen_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_lead_jet_gen_track_reco,JetR,hist_inLeadjet_gen_track_reco,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) correlation(subl_jets_gen, subl_jet_w_gen, tracks_reco, track_w_reco,hist_correlation_signal_subl_jet_gen_track_reco, hist_subl_jet_from_gen_reco_sig, hist_SLJ_trk_from_gen_reco_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_gen_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_subl_jet_gen_track_reco,JetR,hist_inSubljet_gen_track_reco,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_genreco_lead, lead_jets_gen, lead_jet_w_gen, tracks_reco, track_w_reco, mult, vertexz, event_weight, ev_jet_vector_leadjet_gen_reco, jet_weights_leadjet_gen_reco, ev_track_vector_leadjet_gen_reco, trk_weights_leadjet_gen_reco, multvec_leadjet_gen_reco, vzvec_leadjet_gen_reco, weights_leadjet_gen_reco);	// for mixing --> store vectors for mixing
				if(pass_Aj_or_Xj_reco_cut && pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_genreco_subl, subl_jets_gen, subl_jet_w_gen, tracks_reco, track_w_reco, mult, vertexz, event_weight, ev_jet_vector_subleadjet_gen_reco, jet_weights_subleadjet_gen_reco, ev_track_vector_subleadjet_gen_reco, trk_weights_subleadjet_gen_reco, multvec_subleadjet_gen_reco, vzvec_subleadjet_gen_reco, weights_subleadjet_gen_reco);	// for mixing --> store vectors for mixing
				// Gen-Gen
				if(pass_Aj_or_Xj_gen_cut) correlation(lead_jets_gen, lead_jet_w_gen, tracks_gen, track_w_gen, hist_correlation_signal_lead_jet_gen_track_gen, hist_lead_jet_from_gen_gen_sig, hist_LJ_trk_from_gen_gen_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_lead_jet_gen_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_lead_jet_gen_track_gen,JetR,hist_inLeadjet_gen_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_gen_cut) correlation(subl_jets_gen, subl_jet_w_gen, tracks_gen, track_w_gen, hist_correlation_signal_subl_jet_gen_track_gen, hist_subl_jet_from_gen_gen_sig, hist_SLJ_trk_from_gen_gen_sig, event_weight, mult, do_rotation, N_of_rot, hist_correlation_rotation_subl_jet_gen_track_gen, sube_tracks_gen, hist_correlation_signal_subg0_subl_jet_gen_track_gen,JetR,hist_inSubljet_gen_track_gen,do_flow); // calculate correlations
				if(pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_gengen_lead, lead_jets_gen, lead_jet_w_gen, tracks_gen, track_w_gen, mult, vertexz, event_weight, ev_jet_vector_leadjet_gen_gen, jet_weights_leadjet_gen_gen, ev_track_vector_leadjet_gen_gen, trk_weights_leadjet_gen_gen, multvec_leadjet_gen_gen, vzvec_leadjet_gen_gen, weights_leadjet_gen_gen);	// for mixing --> store vectors for mixing
				if(pass_Aj_or_Xj_gen_cut) fillvectors(similar_events, Nev_gengen_subl, subl_jets_gen, subl_jet_w_gen, tracks_gen, track_w_gen, mult, vertexz, event_weight, ev_jet_vector_subleadjet_gen_gen, jet_weights_subleadjet_gen_gen, ev_track_vector_subleadjet_gen_gen, trk_weights_subleadjet_gen_gen, multvec_subleadjet_gen_gen, vzvec_subleadjet_gen_gen, weights_subleadjet_gen_gen);	// for mixing --> store vectors for mixing
			}
			
			if(do_flow){twoparticlecorrelation(tracks_gen, track_w_gen, hist_gen_gen_2pcorrelation_signal, event_weight, mult, sube_tracks_gen, hist_gen_gen_2pcorrelation_signal_subg0, hist_gen_gen_2pcorrelation_signal_subcross);} // calculate 2 particle correlations

			
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
		if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_reco_reco, multiplicity_centrality_bins, vzvec_reco_reco, DVz_range, ev_jet_vector_reco_reco, jet_weights_reco_reco, ev_track_vector_reco_reco, trk_weights_reco_reco, hist_correlation_mixing_jet_reco_track_reco, trk_pt_bins, weights_reco_reco, hist_jet_from_reco_reco_mix, hist_trk_from_reco_reco_mix, double_weight_mix, do_flow);
		if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_reco_reco, multiplicity_centrality_bins, vzvec_leadjet_reco_reco, DVz_range, ev_jet_vector_leadjet_reco_reco, jet_weights_leadjet_reco_reco, ev_track_vector_leadjet_reco_reco, trk_weights_leadjet_reco_reco, hist_correlation_mixing_lead_jet_reco_track_reco, trk_pt_bins, weights_leadjet_reco_reco, hist_lead_jet_from_reco_reco_mix, hist_LJ_trk_from_reco_reco_mix, double_weight_mix, do_flow);
		if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_reco_reco, multiplicity_centrality_bins, vzvec_subleadjet_reco_reco, DVz_range, ev_jet_vector_subleadjet_reco_reco, jet_weights_subleadjet_reco_reco, ev_track_vector_subleadjet_reco_reco, trk_weights_subleadjet_reco_reco, hist_correlation_mixing_subl_jet_reco_track_reco, trk_pt_bins, weights_subleadjet_reco_reco, hist_subl_jet_from_reco_reco_mix, hist_SLJ_trk_from_reco_reco_mix, double_weight_mix, do_flow);
		if(do_flow) call_mix_random_2pc(N_ev_mix, Mult_or_Cent_range, multvec_reco_reco, multiplicity_centrality_bins, vzvec_reco_reco, DVz_range, ev_track_vector_reco_reco, trk_weights_reco_reco, hist_reco_reco_2pcorrelation_mixing, trk_pt_bins, weights_reco_reco, double_weight_mix);
		if(is_MC){
			cout << "Running --> RECO-GEN" << endl;
			// RECO-GEN
			if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_reco_gen, multiplicity_centrality_bins, vzvec_reco_gen, DVz_range, ev_jet_vector_reco_gen, jet_weights_reco_gen, ev_track_vector_reco_gen, trk_weights_reco_gen, hist_correlation_mixing_jet_reco_track_gen, trk_pt_bins, weights_reco_gen, hist_jet_from_reco_gen_mix, hist_trk_from_reco_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_reco_gen, multiplicity_centrality_bins, vzvec_leadjet_reco_gen, DVz_range, ev_jet_vector_leadjet_reco_gen, jet_weights_leadjet_reco_gen, ev_track_vector_leadjet_reco_gen, trk_weights_leadjet_reco_gen, hist_correlation_mixing_lead_jet_reco_track_gen, trk_pt_bins, weights_leadjet_reco_gen, hist_lead_jet_from_reco_gen_mix, hist_LJ_trk_from_reco_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_reco_gen, multiplicity_centrality_bins, vzvec_subleadjet_reco_gen, DVz_range, ev_jet_vector_subleadjet_reco_gen, jet_weights_subleadjet_reco_gen, ev_track_vector_subleadjet_reco_gen, trk_weights_subleadjet_reco_gen, hist_correlation_mixing_subl_jet_reco_track_gen, trk_pt_bins, weights_subleadjet_reco_gen, hist_subl_jet_from_reco_gen_mix, hist_SLJ_trk_from_reco_gen_mix, double_weight_mix, do_flow);
			cout << "Running --> GEN-RECO" << endl;
			// GEN-RECO
			if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_gen_reco, multiplicity_centrality_bins, vzvec_gen_reco, DVz_range, ev_jet_vector_gen_reco, jet_weights_gen_reco, ev_track_vector_gen_reco, trk_weights_gen_reco, hist_correlation_mixing_jet_gen_track_reco, trk_pt_bins, weights_gen_reco, hist_jet_from_gen_reco_mix, hist_trk_from_gen_reco_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_gen_reco, multiplicity_centrality_bins, vzvec_leadjet_gen_reco, DVz_range, ev_jet_vector_leadjet_gen_reco, jet_weights_leadjet_gen_reco, ev_track_vector_leadjet_gen_reco, trk_weights_leadjet_gen_reco, hist_correlation_mixing_lead_jet_gen_track_reco, trk_pt_bins, weights_leadjet_gen_reco, hist_lead_jet_from_gen_reco_mix, hist_LJ_trk_from_gen_reco_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_gen_reco, multiplicity_centrality_bins, vzvec_subleadjet_gen_reco, DVz_range, ev_jet_vector_subleadjet_gen_reco, jet_weights_subleadjet_gen_reco, ev_track_vector_subleadjet_gen_reco, trk_weights_subleadjet_gen_reco, hist_correlation_mixing_subl_jet_gen_track_reco, trk_pt_bins, weights_subleadjet_gen_reco, hist_subl_jet_from_gen_reco_mix, hist_SLJ_trk_from_gen_reco_mix, double_weight_mix, do_flow);
			cout << "Running --> GEN-GEN" << endl;
			// GEN-GEN
			if(do_inclusejettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_gen_gen, multiplicity_centrality_bins, vzvec_gen_gen, DVz_range, ev_jet_vector_gen_gen, jet_weights_gen_gen, ev_track_vector_gen_gen, trk_weights_gen_gen, hist_correlation_mixing_jet_gen_track_gen, trk_pt_bins, weights_gen_gen, hist_jet_from_gen_gen_mix, hist_trk_from_gen_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_leadjet_gen_gen, multiplicity_centrality_bins, vzvec_leadjet_gen_gen, DVz_range, ev_jet_vector_leadjet_gen_gen, jet_weights_leadjet_gen_gen, ev_track_vector_leadjet_gen_gen, trk_weights_leadjet_gen_gen, hist_correlation_mixing_lead_jet_gen_track_gen, trk_pt_bins, weights_leadjet_gen_gen, hist_lead_jet_from_gen_gen_mix, hist_LJ_trk_from_gen_gen_mix, double_weight_mix, do_flow);
			if(do_leading_subleading_jettrack_correlation) call_mix_random(N_ev_mix, Mult_or_Cent_range, multvec_subleadjet_gen_gen, multiplicity_centrality_bins, vzvec_subleadjet_gen_gen, DVz_range, ev_jet_vector_subleadjet_gen_gen, jet_weights_subleadjet_gen_gen, ev_track_vector_subleadjet_gen_gen, trk_weights_subleadjet_gen_gen, hist_correlation_mixing_subl_jet_gen_track_gen, trk_pt_bins, weights_subleadjet_gen_gen, hist_subl_jet_from_gen_gen_mix, hist_SLJ_trk_from_gen_gen_mix, double_weight_mix, do_flow);
			if(do_flow) call_mix_random_2pc(N_ev_mix, Mult_or_Cent_range, multvec_gen_gen, multiplicity_centrality_bins, vzvec_gen_gen, DVz_range, ev_track_vector_gen_gen, trk_weights_gen_gen, hist_gen_gen_2pcorrelation_mixing, trk_pt_bins, weights_gen_gen, double_weight_mix);
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
	MyFile->mkdir("QA_histograms"); 
	MyFile->cd("QA_histograms"); 
	w_QA_hist(is_MC,do_leading_subleading_jettrack_correlation); 
	if(do_inclusejettrack_correlation || do_leading_subleading_jettrack_correlation){
		MyFile->mkdir("correlation_reco_reco_histograms"); 
		MyFile->cd("correlation_reco_reco_histograms");
		w_recoreco_hist(do_mixing,do_rotation,do_inclusejettrack_correlation,do_leading_subleading_jettrack_correlation);
	}
	
	if(is_MC){
		if(do_inclusejettrack_correlation || do_leading_subleading_jettrack_correlation){
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
		MyFile->mkdir("JES_JER"); 
		MyFile->cd("JES_JER");  
		w_jes_jer_hist();
	}

	MyFile->mkdir("jetquenching_histograms"); 
	MyFile->cd("jetquenching_histograms");  
	w_jetquenching_hist(is_MC);
	
	if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow){
		MyFile->mkdir("eventplane_histograms"); 
		MyFile->cd("eventplane_histograms");  
		w_ep_hist(is_MC);
	}
	
	if(do_flow){
		MyFile->mkdir("TwoPC_histograms"); 
		MyFile->cd("TwoPC_histograms");  
		w_2pc_hist(is_MC, do_mixing);
	}

	MyFile->mkdir("etaassymetry_histograms"); 
	MyFile->cd("etaassymetry_histograms");  
	w_etassym_hist(is_MC);
	

	MyFile->Close();

	cout << endl;
	cout << "------------------------------------- DONE --------------------------------------" << endl;
	cout << endl;


	sec_end = clock(); // stop time counting
	cout << "========================================" << endl;
	cout << "Total running time: " << (double)(sec_end - sec_start) / CLOCKS_PER_SEC << " [s]" << endl;
	cout << "========================================" << endl;

	print_stop(); // Print time, date and hour when it stops
	
}
