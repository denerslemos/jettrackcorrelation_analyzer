#include "call_libraries.h"  // call libraries from ROOT and C++
#include "vector_definition.h"  // define the vectors for mixing
#include "read_tree.h" // read the TChains
#include "random_mixing.h" // random mixing function
#include "unfold_weight_correction.h" // unfold correction
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
	if(!is_MC) do_pid = false; // MC only
	if(!do_pid) particles = "";
	if(colliding_system!="pPb") do_CM_pPb = false; // Only do center-of-mass for pPb (or future asymmetric systems like pO)

	//print important informations in the output file
	TString data_or_mc;
	if(!is_MC){data_or_mc="Data";}else{data_or_mc="MC";}
	if(colliding_system == "pPb" && is_pgoing && invert_pgoing){data_or_mc+="_invside";}
	TString simev; if(similar_events){simev = "simevs";}else{simev = "";}
	TString ref_sample = "noref"; if(do_mixing && !do_rotation){ref_sample = Form("mix%iMult%iDVz%.1f%s",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data());}else if(!do_mixing && do_rotation){ref_sample = Form("rot%ievs",N_of_rot);}else if(do_mixing && do_rotation){ref_sample = Form("mix%ievsMult%iDVz%.1f%s_rot%ievs",N_ev_mix,Mult_or_Cent_range,DVz_range,simev.Data(),N_of_rot);}
	TString jet_axis; if(use_WTA){jet_axis = "WTA";}else{jet_axis = "ESC";}
	if(doUE_areabased) jet_axis += "_RhoUE_";
	jet_axis += Form("_jid_%s",jet_idselec.Data());
	jet_axis += Form("_meth%i_",jetid_method);
	jet_axis += Form("3rdm%i_c%.1f_",thirdjet_removal_method,thirdjet_removal_cut);
	jet_axis += Form("rem4thjet%i_",do_fourjet_removal);
	jet_axis += Form("splt%i_",splitMCsample);
	if(dodataweightpt1pt2){jet_axis += "DTw_";}else{jet_axis += "MCw_";}
	TString smear;
	if(do_jeu_down && !do_jeu_up){smear += "_jeudown_";}else if(!do_jeu_down && do_jeu_up){smear += "_jeuup_";}
	if(do_jer_down && !do_jer_up){smear += "_jerdown_";}else if(!do_jer_down && do_jer_up){smear += "_jerup_";}
	TString jet_type; if(do_inclusejettrack_correlation) jet_type += "Incl"; if(do_leading_subleading_jettrack_correlation) jet_type += Form("LS_%s",fwdbkw_jettrk_option.Data()); if(do_dijetstudies) jet_type += "Qch";
	TString XjAj; if(do_Xj_or_Ajcut){XjAj = Form("_Ajmin_%.1f_Ajmax_%.1f_Xjmin_%.1f_Xjmax_%.1f",Ajmin,Ajmax,xjmin,xjmax);}else{XjAj = "";}
	TString isflow;	if(do_flow){isflow="flw";}else{isflow="jtsh";}

	// In case of wrong input, printout error message and kill the job
	if(year_of_datataking!=2012 && year_of_datataking!=2016 && year_of_datataking!=2017 && year_of_datataking!=2018){cout << "Data and MC not supported: choose 2012 for pp at 8 TeV, 2016 for pPb at 8.16 TeV, 2017 for pp at 5.02 TeV or XeXe at 5.44 TeV and 2018 for PbPb at 5.02 TeV" << endl; return;}
	if(colliding_system!="pp" && colliding_system!="pPb" && colliding_system!="XeXe" && colliding_system!="PbPb"){cout << "Data and MC not supported: choose pp for proton-proton, pPb for proton-lead, PbPb for lead-lead and XeXe for xenon-xenon" << endl; return;}
	if(sNN_energy_GeV!=5020 && sNN_energy_GeV!=5440 && sNN_energy_GeV!=8000 && sNN_energy_GeV!=8160 && sNN_energy_GeV!=13000){cout << "Data and MC not supported: 5020 for pp 2017 or PbPb 2018, 5440 for XeXe, 8000 for pp 2018, 8160 for pPb 2016" << endl; return;}
	float sqrts = (float)sNN_energy_GeV;

	// Read Jet Energy Correction file
	vector<string> Files;
	Files.push_back(Form("aux_files/%s_%i/JEC/%s",colliding_system.Data(),sNN_energy_GeV,JEC_file.Data()));
	if(!is_MC) Files.push_back(Form("aux_files/%s_%i/JEC/%s",colliding_system.Data(),sNN_energy_GeV,JEC_file_data.Data()));
	JetCorrector JEC(Files);
	JetUncertainty JEU(Form("aux_files/%s_%i/JEC/%s",colliding_system.Data(),sNN_energy_GeV,JEU_file.Data()));
	//JER sys file
	TFile *filejersys = TFile::Open(Form("aux_files/%s_%i/JEC/%s",colliding_system.Data(),sNN_energy_GeV,JER_sys_file.Data()));
	TH1 *resolution_histo = nullptr;
	if(!do_jer_up && !do_jer_down) filejersys->GetObject("JERnominal", resolution_histo);
	if(do_jer_up && !do_jer_down) filejersys->GetObject("JERup", resolution_histo);
	if(!do_jer_up && do_jer_down) filejersys->GetObject("JERdown", resolution_histo);
	TF1* JetSmear = new TF1("JetSmear","sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))",30,800);
	JetSmear->SetParameters(4.25985e-02, 9.51054e-01, 0.0); // fitted from JER
	TF1* JetScaleCorrection = new TF1("JetScaleCorrection","[3] + ([0]-[3]) / ( 1.0 + pow( x/[2],[1] ) )",30,800);
	if(jet_collection == "ak4PFJetAnalyzer" && doUE_areabased) JetScaleCorrection->SetParameters(3.19771e+00, 9.45364e-01, 9.71077e-01, 1.00049e+00);
	if(jet_collection == "akCs4PFJetAnalyzer" && !doUE_areabased) JetScaleCorrection->SetParameters(2.27008e+00, 9.18625e-01, 1.43067e+00, 1.00002e+00);

	// Track or particle efficiency file
	TFile *fileeff = TFile::Open(Form("aux_files/%s_%i/trk_eff_table/%s",colliding_system.Data(),sNN_energy_GeV,trk_eff_file.Data()));
	TFile *fileefftight = TFile::Open(Form("aux_files/%s_%i/trk_eff_table/Hijing_8TeV_MB_eff_v3_tight.root",colliding_system.Data(),sNN_energy_GeV));
	TFile *fileeffloose = TFile::Open(Form("aux_files/%s_%i/trk_eff_table/Hijing_8TeV_MB_eff_v3_loose.root",colliding_system.Data(),sNN_energy_GeV));

	// Weight correction for unfold
	TFile *fileunfoldweight;
	if(dodataweightpt1pt2){fileunfoldweight = TFile::Open(Form("aux_files/%s_%i/Unfolding/weighthistosData.root",colliding_system.Data(),sNN_energy_GeV));}else{fileunfoldweight = TFile::Open(Form("aux_files/%s_%i/Unfolding/weighthistosMC.root",colliding_system.Data(),sNN_energy_GeV));}

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
	TChain *rho_tree = new TChain("hiFJRhoAnalyzer/rhotree"); // jet rho quantities	
	TChain *gen_tree;
	if(is_MC){gen_tree = new TChain("HiGenParticleAna/hi");} // MC gen particles
	TChain *ski_tree = new TChain("skimanalysis/HltTree"); // event filters
	TChain *ep_tree;
	if(colliding_system=="pPb" && year_of_datataking==2016){ep_tree = new TChain("checkflattening/tree");} // event plane for 2016 pPb data
	
	// loop to add all the trees to the chain
	for (std::vector<TString>::iterator listIterator = file_name_vector.begin(); listIterator != file_name_vector.end(); listIterator++){
		TFile *testfile = TFile::Open(*listIterator,"READ");
		if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){ // safety against corrupted files
			cout << "Adding file " << *listIterator << " to the chains" << endl; // adding files to the chains for each step
			hlt_tree->Add(*listIterator);
			trk_tree->Add(*listIterator);
			hea_tree->Add(*listIterator);
			jet_tree->Add(*listIterator);
			ski_tree->Add(*listIterator);
			rho_tree->Add(*listIterator);
			if(is_MC){gen_tree->Add(*listIterator);}
			if(colliding_system=="pPb" && year_of_datataking==2016){ep_tree->Add(*listIterator);}
		}else{cout << "File: " << *listIterator << " failed!" << endl;}
	}
    file_name_vector.clear();
	
	// Connect all chains
	hlt_tree->AddFriend(trk_tree);
	hlt_tree->AddFriend(hea_tree);
	hlt_tree->AddFriend(jet_tree);
	hlt_tree->AddFriend(ski_tree);
	hlt_tree->AddFriend(rho_tree);
	if(is_MC){hlt_tree->AddFriend(gen_tree);}
	if(colliding_system=="pPb" && year_of_datataking==2016){hlt_tree->AddFriend(ep_tree);}

    // Read the desired branchs in the trees
	read_tree(hlt_tree, is_MC, use_WTA, jet_trigger.Data(), colliding_system.Data(), sNN_energy_GeV, year_of_datataking, event_filter_str, jetid_method); // access the tree informations
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

	// Start loop over events
	double nev = (double)nevents;
	for(int i = 0; i < nevents; i++){
		

		hlt_tree->GetEntry(i);
		if(i != 0 && (i % 10000) == 0){double alpha = (double)i; cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;} // % processed
		if(do_quicktest) if(i != 0 && i % 20000 == 0 ) break; // just for quick tests

		bool splitremoval = false;
		if( (i % 2 != 0) && (splitMCsample == 1) ) { splitremoval = true; } else if ( (i % 2 == 0) && (splitMCsample == 2) ) { splitremoval = true; }
		if(is_MC && splitremoval) continue;

		Nevents->Fill(0); // filled after each event cut

		// Booleans to remove events which does not pass the Aj or Xj selection
		bool pass_Aj_or_Xj_reco_cut = true;
		bool pass_Aj_or_Xj_gen_cut = true;
		if(do_Xj_or_Ajcut){pass_Aj_or_Xj_reco_cut = false; pass_Aj_or_Xj_gen_cut = false;}

		// Apply trigger
		if(dojettrigger){if(jet_trigger_bit != 1) continue;}
		Nevents->Fill(1);

		// Apply event filters
		bool event_filter = false;
		for(int ii = 0; ii < event_filter_str.size(); ii++){ if(event_filter_bool[ii] != 1) event_filter = true; }
		if(event_filter) continue;
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
		std::vector<TVector3> alljets_reco;
		std::vector<TVector3> alljets_ref;
		std::vector<TVector3> alljets_gen;

		//apply event selections
		//Vz
		if(colliding_system == "pPb" && is_pgoing && invert_pgoing) { vertexz = -vertexz; } // in pPb we should invert all z quantities
		if(vertexz < vz_cut_min || vertexz > vz_cut_max) continue;
		Nevents->Fill(3);

		//pthat (MC only)
		if(is_MC){if(pthat <= pthatmin || pthat > pthatmax) continue;} //pthat ranges
		Nevents->Fill(4);

		//multiplicity or centrality
		int trksize = (int)ntrk;
		int mult; // centrality or Ntroff
		if(use_centrality){mult = hiBin;}else{mult = get_Ntrkoff(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);}
		int recomult = get_simple_mult_reco(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge); // reco multiplicity (only pt and eta cuts)
		int genmult;
		if(is_MC) genmult = get_simple_mult_gen(colliding_system, sNN_energy_GeV, year_of_datataking, (int) gen_trkpt->size(), gen_trketa, gen_trkpt, gen_trkchg); // gen multiplicity (only pt and eta cuts)
		double mult_corr = get_Ntrkcorr(fileeff, "nominal", use_centrality, mult, colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkphi, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);
		double mult_corr_tight = get_Ntrkcorr(fileefftight, "tight", use_centrality, mult, colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkphi, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);
		double mult_corr_loose = get_Ntrkcorr(fileeffloose, "loose", use_centrality, mult, colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkphi, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);
		multiplicity_nocut->Fill((double) mult);
		multiplicity_corrected->Fill((double) mult_corr);
		multiplicity_tight->Fill((double) mult_corr_tight);
		multiplicity_loose->Fill((double) mult_corr_loose);
		multiplicity2D->Fill((double) mult, (double) mult_corr);

		// invert the HF, ZDC and EP sides for pPb
		if(colliding_system=="pPb" && is_pgoing && invert_pgoing) invert_eventquantities(hfplus, hfminus, hfplusEta4, hfminusEta4, zdcplus, zdcminus, EP_Psi2_plus_flat, EP_Qx2_plus, EP_Qy2_plus, EP_Mult2_plus, EP_Psi3_plus_flat, EP_Qx3_plus, EP_Qy3_plus, EP_Mult3_plus, EP_Psi4_plus_flat, EP_Qx4_plus, EP_Qy4_plus, EP_Mult4_plus, EP_Psi2_minus_flat, EP_Qx2_minus, EP_Qy2_minus, EP_Mult2_minus, EP_Psi3_minus_flat, EP_Qx3_minus, EP_Qy3_minus, EP_Mult3_minus, EP_Psi4_minus_flat, EP_Qx4_minus, EP_Qy4_minus, EP_Mult4_minus);

		if(colliding_system=="pPb" && year_of_datataking==2016){
			hfhistEta4Sum_nocut->Fill(hfplusEta4+hfminusEta4);
			hfhistEta4Plus_nocut->Fill(hfplusEta4);
			hfhistEta4Minus_nocut->Fill(hfminusEta4);
		}
		
		if(mult < multiplicity_centrality_bins[0] || mult > multiplicity_centrality_bins[multiplicity_centrality_bins.size()-1])continue; //centrality of multiplicity range
		double multcentbin = (double) mult;
		Nevents->Fill(5);
	
		//  add extra dependency
		float extra_variable = hfminusEta4;
		double extrabin = (double) extra_variable;

		// event weight(s), this must be applied in all histograms
		double event_weight = get_event_weight(nev,is_MC, use_centrality, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, vertexz, mult, weight, pthat, extra_variable, is_embedded, is_multdep); // get the event weight

		// Fill vertex, pthat (MC) and multiplicity/centrality histograms
		int multat1 = get_Ntrkoff_at1(colliding_system, sNN_energy_GeV, year_of_datataking, trksize, trketa, trkpt, trkcharge, highpur, trkpterr, trkdcaxy, trkdcaxyerr, trkdcaz, trkdcazerr, trkchi2, trkndof, trknlayer, trknhits, trkalgo, trkmva);	 
		multiplicity->Fill(mult);
		multiplicity_weighted->Fill(mult, event_weight);	
		reco_mult->Fill(recomult);
		reco_mult_weighted->Fill(recomult, event_weight);
		gen_mult->Fill(genmult);
		gen_mult_weighted->Fill(genmult, event_weight);
		multiplicity_weighted_at1->Fill(multat1, event_weight);
		double x_vz[3]={(double) vertexz, (double) multcentbin, (double) extrabin}; vzhist->Fill(x_vz); vzhist_weighted->Fill(x_vz,event_weight);
		double x_vxy[4]={(double) vertexx, (double) vertexy, (double) multcentbin, (double) extrabin}; vxyhist->Fill(x_vxy); vxyhist_weighted->Fill(x_vxy,event_weight);
		double x_pthat[3]={(double) pthat, (double) multcentbin, (double) extrabin}; pthathist->Fill(x_pthat); pthathist_weighted->Fill(x_pthat,event_weight);

		if(colliding_system=="pPb" && year_of_datataking==2016){
			// Fill HF and ZDC histograms 
			double x3D_hiHF[3]={hfplus,hfminus,(double) mult}; hfhist->Fill(x3D_hiHF); hfhist_weighted->Fill(x3D_hiHF,event_weight);
			double x3D_hiHFEta4[3]={hfplusEta4,hfminusEta4,(double) mult}; hfhistEta4->Fill(x3D_hiHFEta4); hfhistEta4_weighted->Fill(x3D_hiHFEta4,event_weight);
			double x3D_hiZDC[3]={zdcplus,zdcminus,(double) mult}; zdchist->Fill(x3D_hiZDC); zdchist_weighted->Fill(x3D_hiZDC,event_weight);
			double x2D_hiHFSum[2]={hfplus+hfminus,(double) mult}; hfhistSum_weighted->Fill(x2D_hiHFSum,event_weight);
			double x2D_hiHFEta4Sum[2]={hfplusEta4+hfminusEta4,(double) mult}; hfhistEta4Sum_weighted->Fill(x2D_hiHFEta4Sum,event_weight);	
			// Fill EP histograms		
			if(do_flow) fillEPhistograms((double) multcentbin, (double) extrabin, event_weight, EP_Psi2_plus_flat, EP_Qx2_plus, EP_Qy2_plus, EP_Mult2_plus, EP_Psi3_plus_flat, EP_Qx3_plus, EP_Qy3_plus, EP_Mult3_plus, EP_Psi4_plus_flat, EP_Qx4_plus, EP_Qy4_plus, EP_Mult4_plus, EP_Psi2_minus_flat, EP_Qx2_minus, EP_Qy2_minus, EP_Mult2_minus, EP_Psi3_minus_flat, EP_Qx3_minus, EP_Qy3_minus, EP_Mult3_minus,  EP_Psi4_minus_flat, EP_Qx4_minus, EP_Qy4_minus, EP_Mult4_minus, EP2_plus_flat, EP2_minus_flat, EP3_plus_flat, EP3_minus_flat, EP4_plus_flat, EP4_minus_flat);
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

			// Calculate track efficiency correction
			double trk_weight = 1.0;
			trk_weight = trk_weight*getTrkCorrWeight(fileeff, use_centrality, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, mult, trk_pt, trk_eta, trk_phi);

			trk_eta = trk_eta + boost; // In pPb case, for the center-of-mass correction if needed
			if(colliding_system == "pPb" && is_pgoing && invert_pgoing) trk_eta = -trk_eta; // in case of merging pgoing and Pbgoing use this

			// Fill track QA histograms
			double x_reco_trk[5]={trk_pt,trk_eta,trk_phi,(double) multcentbin,(double) extrabin}; 
			hist_reco_trk->Fill(x_reco_trk);
			hist_reco_trk_corr->Fill(x_reco_trk,trk_weight);
			hist_reco_trk_weighted->Fill(x_reco_trk,trk_weight*event_weight);

			// Fill track vectors
			TVector3 GoodTracks;
			GoodTracks.SetPtEtaPhi(trk_pt, trk_eta, trk_phi);
			tracks_reco.push_back(GoodTracks);
			sube_tracks_reco.push_back(0); // set == 0 because this is only valid for gen
			double trk_etamix_weight = 1.0;// set to one, not useful so far
			track_w_reco.push_back(trk_weight*trk_etamix_weight); // save weight to apply in the mixing

			double trackbin = (double) trk_pt;
			if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filltrkEPhistograms(trk_phi, trackbin, (double) multcentbin, (double) extrabin, trk_weight*event_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_EP2_flat_trk_plus, Dphi_EP2_flat_trk_minus, Dphi_EP3_flat_trk_plus, Dphi_EP3_flat_trk_minus, Dphi_EP4_flat_trk_plus, Dphi_EP4_flat_trk_minus);
			
		} // End loop over tracks

		// Used to find leading, subleading and third jets
		float leadrecojet_pt=-999999999.9, leadrecojet_eta=-999999999.9, leadrecojet_phi=-999999999.9, leadrecojet_mass=-999999999.9, leadrecojet_flavor=-999999999.9;  // leading jet quantities
		float sublrecojet_pt=-999999999.9, sublrecojet_eta=-999999999.9, sublrecojet_phi=-999999999.9, sublrecojet_mass=-999999999.9, sublrecojet_flavor=-999999999.9; // subleading jet quantities
		float thirdrecojet_pt=-999999999.9, thirdrecojet_eta=-999999999.9, thirdrecojet_phi=-999999999.9, thirdrecojet_mass=-999999999.9, thirdrecojet_flavor=-999999999.9; // third jet quantities
		float leadrefjet_pt=-999999999.9, leadrefjet_eta=-999999999.9, leadrefjet_phi=-999999999.9, leadrefjet_mass=-999999999.9, leadrefjet_flavor=-999999999.9; // leading jet ref quantities
		float sublrefjet_pt=-999999999.9, sublrefjet_eta=-999999999.9, sublrefjet_phi=-999999999.9, sublrefjet_mass=-999999999.9, sublrefjet_flavor=-999999999.9; // subleading jet ref quantities
		float thirdrefjet_pt=-999999999.9, thirdrefjet_eta=-999999999.9, thirdrefjet_phi=-999999999.9, thirdrefjet_mass=-999999999.9, thirdrefjet_flavor=-999999999.9; // third jet quantities
		int leadrecojet_index=2000, sublrecojet_index=2000, thirdrecojet_index=2000; // jet reco index
		int leadrefjet_index=2000, sublrefjet_index=2000, thirdrefjet_index=2000; // jet ref index
		float fourthrecojet_pt=-999999999.9, fourthrefjet_pt=-999999999.9; 
				
		bool isjetincluded = false;
		int njets = 0;
	
		int jetsize = (int)nref; // number of jets in an event

		std::vector<int> jetidremovalindex;

		// -> Start loop over jets
		for (int j = 0; j < jetsize; j++){
			
			// control selections and plots
			//if(fabs(jteta[j]) > 5.1) continue;
			double x_trkmaxjet[3]={trackMax[j]/rawpt[j], rawpt[j], (double) multcentbin}; 
			if(trackMax[j]/rawpt[j] >= 0) jettrackmaxptinjethisto->Fill(x_trkmaxjet,event_weight);
			if(trackMax[j]/rawpt[j] >  0) jettrackmaxptinjethisto_no0->Fill(x_trkmaxjet,event_weight);
			if(is_MC && refpt[j] > 0) jettrackmaxptinjethisto_ref->Fill(x_trkmaxjet,event_weight);

			// Define jet kinematics
			float jet_rawpt = rawpt[j];
			float jet_eta = jteta[j];
			float jet_phi = jtphi[j];
			float jet_mass = jtmass[j];
			/*
			if(jetid_method == 1 || jetid_method == 2 ){
	 			double x_NHFbc[4]={(double)jtPfNHF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_NHF_hist_beforeid->Fill(x_NHFbc,event_weight);
	 			double x_NEFbc[4]={(double)jtPfNEF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_NEF_hist_beforeid->Fill(x_NEFbc,event_weight);
	 			double x_CHFbc[4]={(double)jtPfCHF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_CHF_hist_beforeid->Fill(x_CHFbc,event_weight);
	 			double x_MUFbc[4]={(double)jtPfMUF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_MUF_hist_beforeid->Fill(x_MUFbc,event_weight);
	 			double x_CEFbc[4]={(double)jtPfCEF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_CEF_hist_beforeid->Fill(x_CEFbc,event_weight);
	 			double x_CMbc[4]={(double)(jtPfCHM[j] + jtPfCEM[j] + jtPfMUM[j]), (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_chargedmult_hist_beforeid->Fill(x_CMbc,event_weight);
	 			double x_NMbc[4]={(double)(jtPfNHM[j] + jtPfNEM[j]), (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_neutralmult_hist_beforeid->Fill(x_NMbc,event_weight);
	 			double x_Mbc[4]={(double)(jtPfCHM[j] + jtPfCEM[j] + jtPfMUM[j] + jtPfNHM[j] + jtPfNEM[j]), (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_mult_hist_beforeid->Fill(x_Mbc,event_weight);
			}
			*/
			if( jetid_method == 1 || jetid_method == 2 ){ 
				bool isjetid = passJetIDcuts(jet_idselec, year_of_datataking, jet_eta, jtPfNHF[j], jtPfNEF[j], jtPfCHF[j], jtPfMUF[j], jtPfCEF[j], jtPfCHM[j], jtPfCEM[j], jtPfNHM[j], jtPfNEM[j], jtPfMUM[j]);
				if(!isjetid) { 
					jetidremovalindex.push_back(j); 
					if( jetid_method == 1 ) continue; 
				}
			}
			if( jetid_method == 3 || jetid_method == 4 ){ 
				bool isjetid = true;
				if(trackMax[j]/rawpt[j] < 0.01) isjetid = false; 
				if(trackMax[j]/rawpt[j] > 0.98) isjetid = false; 
				if(!isjetid) { 
					jetidremovalindex.push_back(j); 
					if( jetid_method == 3 ) continue; 
				}				
			}
			if(trackMax[j] < trackmaxpt) continue; // Can be use to remove jets from low pT tracks
			/*
			if(jetid_method == 1 || jetid_method == 2 ){
	 			double x_NHFac[4]={(double)jtPfNHF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_NHF_hist_afterid->Fill(x_NHFac,event_weight);
	 			double x_NEFac[4]={(double)jtPfNEF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_NEF_hist_afterid->Fill(x_NEFac,event_weight);
	 			double x_CHFac[4]={(double)jtPfCHF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_CHF_hist_afterid->Fill(x_CHFac,event_weight);
	 			double x_MUFac[4]={(double)jtPfMUF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_MUF_hist_afterid->Fill(x_MUFac,event_weight);
	 			double x_CEFac[4]={(double)jtPfCEF[j], (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_CEF_hist_afterid->Fill(x_CEFac,event_weight);
	 			double x_CMac[4]={(double)(jtPfCHM[j] + jtPfCEM[j] + jtPfMUM[j]), (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_chargedmult_hist_afterid->Fill(x_CMac,event_weight);
	 			double x_NMac[4]={(double)(jtPfNHM[j] + jtPfNEM[j]), (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_neutralmult_hist_afterid->Fill(x_NMac,event_weight);
	 			double x_Mac[4]={(double)(jtPfCHM[j] + jtPfCEM[j] + jtPfMUM[j] + jtPfNHM[j] + jtPfNEM[j]), (double) jet_eta, (double) jet_rawpt,(double) multcentbin}; jet_mult_hist_afterid->Fill(x_Mac,event_weight);
			}
			*/
			// UE subtraction
 			double UE = GetUE(etamin, etamax, rho, jet_eta, JetR);
 			double AverageRho = UE / (JetR * JetR * TMath::Pi());
 			double checkcorrection = jet_rawpt - UE;
 			double x_ue[3]={UE,(double) multcentbin,(double) extrabin}; histo_jetUE->Fill(x_ue,event_weight);
 			double x_rho[3]={AverageRho,(double) multcentbin,(double) extrabin}; histo_jetAverageRho->Fill(x_rho,event_weight);
 			double x_check[3]={checkcorrection,(double) multcentbin,(double) extrabin}; histo_jetcheckcorrection->Fill(x_check,event_weight);

 			// Apply JEC
			JEC.SetJetPT(jet_rawpt); 
			JEC.SetJetEta(jet_eta); 
			JEC.SetJetPhi(jet_phi);
			float jet_pt_corr = JEC.GetCorrectedPT();
						
			if(doUE_areabased) jet_pt_corr = (jet_rawpt - UE)*JEC.GetCorrection();
			
			float JEScorrection = 1.0;
			if(jet_collection.EqualTo("ak4PFJetAnalyzer") && doUE_areabased) JEScorrection = JetScaleCorrection->Eval(jet_pt_corr);
			if(jet_collection.EqualTo("akCs4PFJetAnalyzer") && !doUE_areabased) JEScorrection = JetScaleCorrection->Eval(jet_pt_corr);
			jet_pt_corr = jet_pt_corr * JEScorrection;

			if(is_MC && (!do_jer_up || !do_jer_down)) {
			
				double resolution_factor = resolution_histo->GetBinContent( resolution_histo->GetXaxis()->FindBin(jet_eta) );
				double extraResolution = TMath::Sqrt(TMath::Max(resolution_factor*resolution_factor-1.0,0.0)); // found jet resolution
				double sigma_smear = extraResolution*JetSmear->Eval(jet_pt_corr); // some % worst --> from JetMET
				double mu_smar = 1.0;
				double smear = gRandom->Gaus(mu_smar,sigma_smear);
//				while( smear < 0 ){ smear = gRandom->Gaus(mu_smar,sigma_smear); }
				if(  smear < 0 ) cout << "Negative smear" << endl;
				jet_pt_corr = jet_pt_corr*smear;	
					
			}

			// Apply JEU for systematics in data!
			JEU.SetJetPT(jet_pt_corr);
			JEU.SetJetEta(jet_eta);
			JEU.SetJetPhi(jet_phi);
			if(do_jeu_down && !do_jeu_up){jet_pt_corr = jet_pt_corr * (1 - JEU.GetUncertainty().first);}else if(!do_jeu_down && do_jeu_up){jet_pt_corr = jet_pt_corr * (1 + JEU.GetUncertainty().second);}

			int refpartonfromB = 0; 
			if(is_MC){
				if(fabs(refparton_flavorForB[j]) >= 1 && fabs(refparton_flavorForB[j]) <= 6){
					refpartonfromB = fabs(refparton_flavorForB[j]);
				}else if(fabs(refparton_flavorForB[j]) == 21){
					refpartonfromB = 7;
				}else{
				refpartonfromB = 0;}
			}
			float jet_flavor = (float) refpartonfromB;
			int jet_index = (int) j;

			//leading and subleading
			find_leading_subleading_third(jet_pt_corr,jet_eta,jet_phi,jet_mass,jet_flavor,jet_index,leadrecojet_pt,leadrecojet_eta,leadrecojet_phi,leadrecojet_mass,leadrecojet_flavor,leadrecojet_index,sublrecojet_pt,sublrecojet_eta,sublrecojet_phi,sublrecojet_mass,sublrecojet_flavor,sublrecojet_index,thirdrecojet_pt,thirdrecojet_eta,thirdrecojet_phi,thirdrecojet_mass,thirdrecojet_flavor,thirdrecojet_index,fourthrecojet_pt); // Find leading and subleading jets

			jet_eta = jet_eta + boost; // In pPb case, for the center-of-mass correction if needed
			if(colliding_system == "pPb" && is_pgoing && invert_pgoing)jet_eta = -jet_eta;

			double jet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, jet_pt_corr, jet_eta); // Jet weight (specially for MC)

			// Fill reco jet QA histograms
			double x_reco_jet[5]={jet_rawpt,jet_eta,jet_phi,(double) multcentbin,(double) extrabin}; 
			hist_reco_jet->Fill(x_reco_jet);
			hist_reco_jet_weighted->Fill(x_reco_jet,event_weight*jet_weight);
			double x_reco_jet_corr[5]={jet_pt_corr,jet_eta,jet_phi,(double) multcentbin,(double) extrabin}; 
			hist_reco_jet_corr->Fill(x_reco_jet_corr);
			hist_reco_jet_corr_weighted->Fill(x_reco_jet_corr,event_weight*jet_weight);
		
			TVector3 AllJets;
			AllJets.SetPtEtaPhi(jet_pt_corr, jet_eta, jet_phi);
			alljets_reco.push_back(AllJets);

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
					
				if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filljetEPhistograms(jet_phi, (double) multcentbin, (double) extrabin, event_weight*jet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_flat_EP2_inclusive_plus, Dphi_flat_EP2_inclusive_minus, Dphi_flat_EP3_inclusive_plus, Dphi_flat_EP3_inclusive_minus, Dphi_flat_EP4_inclusive_plus, Dphi_flat_EP4_inclusive_minus);
			
			}

			if(is_MC){ 

				float ref_pt = refpt[j];
				float ref_eta = refeta[j];
				float ref_phi = refphi[j];
				float ref_mass = refmass[j];

    	        if(jet_rawpt < 0.0) continue;
 	            if(jet_pt_corr < 0.0) continue;
 	            //if(fabs(ref_eta) > 5.1) continue;
    	        if(ref_pt < 0.0) continue;

				int jet_index_ref = (int) j;

				find_leading_subleading_third(ref_pt,ref_eta,ref_phi,ref_mass,jet_flavor,jet_index_ref,leadrefjet_pt,leadrefjet_eta,leadrefjet_phi,leadrefjet_mass,leadrefjet_flavor,leadrefjet_index,sublrefjet_pt,sublrefjet_eta,sublrefjet_phi,sublrefjet_mass,sublrefjet_flavor,sublrefjet_index,thirdrefjet_pt,thirdrefjet_eta,thirdrefjet_phi,thirdrefjet_mass,thirdrefjet_flavor,thirdrefjet_index,fourthrefjet_pt); // Find leading and subleading ref jets
				
				float ref_eta_lab = ref_eta;
				ref_eta = ref_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				if(colliding_system == "pPb" && is_pgoing && invert_pgoing)ref_eta = -ref_eta;

				double refjet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, ref_pt, ref_eta); // Jet weight (specially for MC)

				double x_ref_QA[5]={ref_pt, ref_eta, ref_phi, (double)multcentbin, (double) extrabin}; 
				hist_ref_jet_weighted->Fill(x_ref_QA,event_weight*refjet_weight);

				double JES_ratio_reco_vs_ref = jet_pt_corr/ref_pt;
				double x_JES_ratio_reco_vs_reffromB[6]={JES_ratio_reco_vs_ref,ref_pt,ref_eta,(double)refpartonfromB,(double)multcentbin,(double) extrabin}; 
				hist_jes_reco_fromB_weighted->Fill(x_JES_ratio_reco_vs_reffromB,event_weight*refjet_weight*jet_weight);
				double x_JES_ratio_reco_vs_ref[6]={JES_ratio_reco_vs_ref,ref_pt,ref_eta_lab,(double)refpartonfromB,(double)multcentbin,(double) extrabin}; 
				hist_jes_reco_weighted->Fill(x_JES_ratio_reco_vs_ref,event_weight*refjet_weight*jet_weight);
				
				TVector3 AllJetsRef;
				AllJetsRef.SetPtEtaPhi(ref_pt, ref_eta, ref_phi);
				alljets_ref.push_back(AllJetsRef);
				
			} // End if over ref (MC)

		}
		// -> End loop over jets

		// Save quantities for events with inclusive jets		
		if(isjetincluded){
			multiplicity_withonejet_weighted->Fill(mult,event_weight);
			reco_mult_withonejet_weighted->Fill(recomult, event_weight);
			multiplicity_withonejet_weighted_at1->Fill(multat1, event_weight);
			vzhist_jet_weighted->Fill(x_vz, event_weight);
			if(colliding_system=="pPb" && year_of_datataking==2016){
				double x3D_hiHF_onejet[3]={hfplus,hfminus,(double) mult}; hfhist_onejet_weighted->Fill(x3D_hiHF_onejet,event_weight);
				double x3D_hiHFEta4_onejet[3]={hfplusEta4,hfminusEta4,(double) mult}; hfhistEta4_onejet_weighted->Fill(x3D_hiHFEta4_onejet,event_weight);
				double x3D_hiZDC_onejet[3]={zdcplus,zdcminus,(double) mult}; zdchist_onejet_weighted->Fill(x3D_hiZDC_onejet,event_weight);
				double x2D_hiHFSum_onejet[2]={hfplus+hfminus,(double) mult}; hfhistSum_onejet_weighted->Fill(x2D_hiHFSum_onejet,event_weight);
				double x2D_hiHFEta4Sum_onejet[2]={hfplusEta4+hfminusEta4,(double) mult}; hfhistEta4Sum_onejet_weighted->Fill(x2D_hiHFEta4Sum_onejet,event_weight);			
			}
		}
		double xnjets[3]={(double) njets, (double) multcentbin, (double) extrabin}; 
		NJets->Fill(xnjets,event_weight);

		// ----> RECO <------
		// 3rd jet removal
		bool removethirdjet = false;
		if(thirdjet_removal_method == 1){ if(thirdrecojet_pt > thirdjet_removal_cut) removethirdjet = true; }
		if(thirdjet_removal_method == 2){ if(thirdrecojet_pt > thirdjet_removal_cut * sublrecojet_pt) removethirdjet = true; }
		if(thirdjet_removal_method == 3){ if(thirdrecojet_pt > thirdjet_removal_cut * 0.5 * (leadrecojet_pt + sublrecojet_pt)) removethirdjet = true; }

		// 4th jet removal
		bool removefourjet = false;
		if(do_fourjet_removal == 1 && fourthrecojet_pt > 0.0){ if( fourthrecojet_pt+thirdrecojet_pt > sublrecojet_pt ) removefourjet = true; } 
		if(do_fourjet_removal == 2 && fourthrecojet_pt > 0.0){ if( fourthrecojet_pt+thirdrecojet_pt > subleading_pT_min ) removefourjet = true; } 
		if(do_fourjet_removal == 3 && fourthrecojet_pt > 0.0){ if( fabs(fourthrecojet_pt-thirdrecojet_pt) > 10.0 ) removefourjet = true; } 
		if(do_fourjet_removal == 4 && fourthrecojet_pt > 0.0){ if( fabs(fourthrecojet_pt-thirdrecojet_pt) > 25.0 ) removefourjet = true; } 
		if(do_fourjet_removal == 5 && fourthrecojet_pt > 0.0){ if( fabs(fourthrecojet_pt-thirdrecojet_pt) > subleading_pT_min ) removefourjet = true; } 

		// weights for reco and ref
		double ljet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrecojet_pt, leadrecojet_eta);  // Jet weight (specially for MC)
		double sljet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrecojet_pt, sublrecojet_eta);  // Jet weight (specially for MC)
        double lrefjet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadrefjet_pt, leadrefjet_eta);  // Jet weight (specially for MC)
        double slrefjet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublrefjet_pt, sublrefjet_eta);  // Jet weight (specially for MC)

		bool removejetid = true;
		if( ( jetid_method == 2 ||  jetid_method == 4 ) && jetidremovalindex.size() > 0) { 
			bool leadjetid = comparevectorandindexforjetid(jetidremovalindex, leadrecojet_index);
			bool subljetid = comparevectorandindexforjetid(jetidremovalindex, sublrecojet_index);
			if(!leadjetid) removejetid = false;
			if(!subljetid) removejetid = false;
		}

		// reco dijets
		bool tworecojets = twojetfounded(leadrecojet_pt, leading_pT_min, jet_pt_max_cut, sublrecojet_pt, subleading_pT_min);
		bool isrecodijets = isdijet(leadrecojet_pt, leading_pT_min, jet_pt_max_cut, sublrecojet_pt, subleading_pT_min, leadrecojet_phi, sublrecojet_phi, leading_subleading_deltaphi_min, xjmin, xjmax, Ajmin, Ajmax);
		bool goodrecoevent = (leadrecojet_pt > 0.0 && sublrecojet_pt > 0.0 && !removethirdjet && !removefourjet && removejetid);

		// for reco jets
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

		// for ref matched jets
		double leadrefmatchjet_eta = refeta[leadrecojet_index]; // before boost for eta dijet
		double sublrefmatchjet_eta = refeta[sublrecojet_index]; // before boost for eta dijet
		leadrefmatchjet_eta = leadrefmatchjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
		sublrefmatchjet_eta = sublrefmatchjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
		if(colliding_system == "pPb" && is_pgoing && invert_pgoing){
			leadrefmatchjet_eta = -leadrefmatchjet_eta; 
			sublrefmatchjet_eta = -sublrefmatchjet_eta;
		}
		
		// for ref jets
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
		// -> Start of dijet calculations
		if(goodrecoevent){
			Nevents->Fill(6);	

			if(tworecojets){ 	//leading/subleading pT cuts
				Nevents->Fill(7);

				// Fill leading/subleading jet quenching quantities
				double delta_phi_reco = fabs(deltaphi(leadrecojet_phi, sublrecojet_phi));
				double Aj_reco = asymmetry(leadrecojet_pt,sublrecojet_pt);
				double Xj_reco = xjvar(leadrecojet_pt,sublrecojet_pt);
				float ptdijetreco = 0.5*(leadrecojet_pt + sublrecojet_pt);
				double ptdijetbinreco = (double) ptdijetreco;
				float dijetetarecotype = dijetetabin(leadrecojet_eta, sublrecojet_eta, jet_eta_min_cut, jet_eta_max_cut, jet_fwd_eta_min_cut, jet_fwd_eta_max_cut, jet_bkw_eta_min_cut, jet_bkw_eta_max_cut);

				if(do_dijetstudies){
					double x_reco[9]={Xj_reco,Aj_reco,delta_phi_reco,(double)multcentbin,(double)ptdijetbinreco,(double)extrabin,(double)leadrecojet_pt,(double)sublrecojet_pt, (double) dijetetarecotype}; 
					hist_reco_lead_reco_subl_quench->Fill(x_reco, event_weight*ljet_weight*sljet_weight);
					hist_reco_lead_reco_subl_quench_unweighted->Fill(x_reco);
					float correctionmapweightreco = getUnfCorrWeight(fileunfoldweight, leadrecojet_pt, sublrecojet_pt, mult, dijetetarecotype);//do correction here
					if(!is_MC) correctionmapweightreco = 1.0;
					hist_reco_lead_reco_subl_quench_checkmapcorrection->Fill(x_reco, correctionmapweightreco*ljet_weight*sljet_weight);
					if(is_MC && refpt[leadrecojet_index] < 0)hist_fake_lead_reco_subl_quench->Fill(x_reco, event_weight*ljet_weight*sljet_weight);
					if(is_MC && refpt[sublrecojet_index] < 0)hist_reco_lead_fake_subl_quench->Fill(x_reco, event_weight*ljet_weight*sljet_weight);
					if(is_MC && refpt[leadrecojet_index] < 0 && refpt[sublrecojet_index] < 0)hist_fake_lead_fake_subl_quench->Fill(x_reco, event_weight*ljet_weight*sljet_weight);
					filletadijethistograms( (double) sqrts, leadrecojet_eta_lab, sublrecojet_eta_lab, leadrecojet_eta, sublrecojet_eta, leadrecojet_pt, leadrecojet_phi, leadrecojet_mass, sublrecojet_pt, sublrecojet_phi, sublrecojet_mass, (double)multcentbin, (double)extrabin, (double) event_weight*ljet_weight*sljet_weight, hist_etaDijet_reco, hist_etaDijet_CM_reco, hist_yDijet_CM_reco );
					if(colliding_system=="pPb" && year_of_datataking==2016) fillxjEPhistograms_withfake(is_MC, Xj_reco, delta_phi_reco, leadrecojet_phi, refpt[leadrecojet_index], (double) multcentbin, (double)extrabin, (double)ptdijetbinreco, (double) dijetetarecotype, event_weight*ljet_weight*sljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, hist_reco_leadEP_quench_plus, hist_reco_leadEP_quench_minus, hist_fake_leadEP_quench_plus, hist_fake_leadEP_quench_minus);
					if(colliding_system=="pPb" && year_of_datataking==2016) fillxjEPhistograms_withfake(is_MC, Xj_reco, delta_phi_reco, sublrecojet_phi, refpt[sublrecojet_index], (double) multcentbin, (double)extrabin, (double)ptdijetbinreco, (double) dijetetarecotype, event_weight*ljet_weight*sljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, hist_reco_sublEP_quench_plus, hist_reco_sublEP_quench_minus, hist_fake_sublEP_quench_plus, hist_fake_sublEP_quench_minus);
				/*					
					if(is_MC){
						double xjptaveunfreco = TransformToUnfoldingAxis_xjptave(Xj_reco, ptdijetbinreco,xjmax);						
						double x_unf_meas_xjptavereco[4]={xjptaveunfreco, delta_phi_reco, (double)dijetetarecotype,(double)multcentbin}; 
						fhUnfoldingMeasu_xjptave->Fill(x_unf_meas_xjptavereco,event_weight*ljet_weight*sljet_weight);
						double pt1pt2unfreco = TransformToUnfoldingAxis_pt1pt2(leadrecojet_pt, sublrecojet_pt, jet_pt_max_cut);						
						double x_unf_meas_pt1pt2reco[4]={pt1pt2unfreco, delta_phi_reco, (double)dijetetarecotype,(double)multcentbin}; 
						fhUnfoldingMeasu_pt1pt2->Fill(x_unf_meas_pt1pt2reco,event_weight*ljet_weight*sljet_weight);
						double xjpt2unfreco = TransformToUnfoldingAxis_xjptave(Xj_reco, sublrecojet_pt,xjmax);						
						double x_unf_meas_xjpt2reco[5]={xjpt2unfreco, delta_phi_reco, (double)dijetetarecotype,(double)multcentbin, (double)leadrecojet_pt}; 
						fhUnfoldingMeasu_xjpt2->Fill(x_unf_meas_xjpt2reco,event_weight*ljet_weight*sljet_weight);
					}
				*/

				}

				if(isrecodijets){ // leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
					Nevents->Fill(8);
					pass_Aj_or_Xj_reco_cut = true; // if we apply Xj or Aj cuts

					// Sort the vector by Pt
					std::vector<TVector3> sortedVectorReco = alljets_reco;
					std::sort(sortedVectorReco.begin(), sortedVectorReco.end(), sortByPt);
					TVector3 GoodLeadingJets_reco_3vec;
					GoodLeadingJets_reco_3vec.SetPtEtaPhi(leadrecojet_pt,leadrecojet_eta_lab,leadrecojet_phi);
					TVector3 GoodSubLeadingJets_reco_3vec;
					GoodSubLeadingJets_reco_3vec.SetPtEtaPhi(sublrecojet_pt, sublrecojet_eta_lab, sublrecojet_phi);
					TVector3 GoodTrdLeadingJets_reco_3vec;
					if (sortedVectorReco.size() > 2){ GoodTrdLeadingJets_reco_3vec = findClosestVectorInPhi(sortedVectorReco); }
					TVector3 GoodProjLeadingJets_reco_3vec = -(GoodLeadingJets_reco_3vec + GoodSubLeadingJets_reco_3vec);
					TVector3 GoodProjLeadingJets_reco_3vec_diff = (GoodLeadingJets_reco_3vec - GoodSubLeadingJets_reco_3vec);
					bool remove_3rdjetevent = false;
					if(thirdrecojet_pt > 0 && thirdjet_removal_method == 4 && (fabs(deltaphi(GoodTrdLeadingJets_reco_3vec.Phi(),GoodProjLeadingJets_reco_3vec.Phi())) < thirdjet_removal_cut)) remove_3rdjetevent = true;
		
					if(thirdrecojet_pt > 0){
						double delphi_3rdProjReco = fabs(deltaphi(GoodTrdLeadingJets_reco_3vec.Phi(),GoodProjLeadingJets_reco_3vec.Phi()));
						double deleta_3rdProjReco = deltaeta(GoodTrdLeadingJets_reco_3vec.Eta(),GoodProjLeadingJets_reco_3vec.Eta());
						double xjetprojReco[3] = {delphi_3rdProjReco, deleta_3rdProjReco, (double)multcentbin};
						jetjet_Dphi_Deta_reco->Fill(xjetprojReco,event_weight);
						double delphi_3rdProjRecodiff = fabs(deltaphi(GoodTrdLeadingJets_reco_3vec.Phi(),GoodProjLeadingJets_reco_3vec_diff.Phi()));
						double deleta_3rdProjRecodiff = deltaeta(GoodTrdLeadingJets_reco_3vec.Eta(),GoodProjLeadingJets_reco_3vec_diff.Eta());
						double xjetprojRecodiff[3] = {delphi_3rdProjRecodiff, deleta_3rdProjRecodiff, (double)multcentbin};
						jetjet_Dphi_Deta_reco_diff->Fill(xjetprojRecodiff,event_weight);
					}
						
					// Fill leading and subleading jet QA histograms
					double x_lead[5]={leadrecojet_pt,leadrecojet_eta,leadrecojet_phi,(double) multcentbin,(double)extrabin}; 
					hist_reco_leadjet->Fill(x_lead);
					hist_reco_leadjet_weighted->Fill(x_lead,event_weight*ljet_weight);
					double x_sublead[5]={sublrecojet_pt,sublrecojet_eta,sublrecojet_phi,(double) multcentbin,(double)extrabin}; 
					hist_reco_subljet->Fill(x_sublead);
					hist_reco_subljet_weighted->Fill(x_sublead,event_weight*sljet_weight);
					double x_3rdlead[5]={thirdrecojet_pt,thirdrecojet_eta,thirdrecojet_phi,(double) multcentbin,(double)extrabin}; 
					hist_reco_thrdjet_weighted->Fill(x_3rdlead,event_weight);
					double x_3rdproj[5]={GoodProjLeadingJets_reco_3vec.Pt(),GoodProjLeadingJets_reco_3vec.Eta(),GoodProjLeadingJets_reco_3vec.Phi(),(double) multcentbin,(double)extrabin}; 
					hist_reco_thrdproj_weighted->Fill(x_3rdproj,event_weight);
					double x_3rdprojdiff[5]={GoodProjLeadingJets_reco_3vec_diff.Pt(),GoodProjLeadingJets_reco_3vec_diff.Eta(),GoodProjLeadingJets_reco_3vec_diff.Phi(),(double) multcentbin,(double)extrabin}; 
					hist_reco_thrdprojdiff_weighted->Fill(x_3rdprojdiff,event_weight);
						
					bool usethisrecojets = false;
					multiplicity_etadep_weighted->Fill( (double) mult, (double) dijetetarecotype, event_weight);
					if(dijetetarecotype == 0.5 && fwdbkw_jettrk_option == "mid_mid"){usethisrecojets = true;}
					if(dijetetarecotype == 1.5 && fwdbkw_jettrk_option == "mid_fwd"){usethisrecojets = true;}
					if(dijetetarecotype == 2.5 && fwdbkw_jettrk_option == "mid_bkw"){usethisrecojets = true;}
					if(dijetetarecotype == 3.5 && fwdbkw_jettrk_option == "fwd_mid"){usethisrecojets = true;}
					if(dijetetarecotype == 4.5 && fwdbkw_jettrk_option == "bkw_mid"){usethisrecojets = true;}
					if(dijetetarecotype == 5.5 && fwdbkw_jettrk_option == "fwd_fwd"){usethisrecojets = true;}
					if(dijetetarecotype == 6.5 && fwdbkw_jettrk_option == "fwd_bkw"){usethisrecojets = true;}
					if(dijetetarecotype == 7.5 && fwdbkw_jettrk_option == "bkw_fwd"){usethisrecojets = true;}
					if(dijetetarecotype == 8.5 && fwdbkw_jettrk_option == "bkw_bkw"){usethisrecojets = true;}

					// Fill leading and subleading jet vectors for correlations
					TVector3 GoodLeadingJets_reco;
					TVector3 GoodSubLeadingJets_reco;
					if(usethisrecojets){
						GoodLeadingJets_reco.SetPtEtaPhi(leadrecojet_pt, leadrecojet_eta, leadrecojet_phi);
						lead_jets_reco.push_back(GoodLeadingJets_reco);
						lead_jet_w_reco.push_back(ljet_weight);
						GoodSubLeadingJets_reco.SetPtEtaPhi(sublrecojet_pt, sublrecojet_eta, sublrecojet_phi);
						subl_jets_reco.push_back(GoodSubLeadingJets_reco);
						subl_jet_w_reco.push_back(sljet_weight);
					}
					// Fill EP histograms
					if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filljetEPhistograms(leadrecojet_phi, (double) multcentbin, (double) extrabin, event_weight*ljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_flat_EP2_leading_plus, Dphi_flat_EP2_leading_minus, Dphi_flat_EP3_leading_plus, Dphi_flat_EP3_leading_minus, Dphi_flat_EP4_leading_plus, Dphi_flat_EP4_leading_minus);
					if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filljetEPhistograms(sublrecojet_phi, (double) multcentbin, (double) extrabin, event_weight*sljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_flat_EP2_subleading_plus, Dphi_flat_EP2_subleading_minus, Dphi_flat_EP3_subleading_plus, Dphi_flat_EP3_subleading_minus, Dphi_flat_EP4_subleading_plus, Dphi_flat_EP4_subleading_minus);

					multiplicity_withdijets_weighted->Fill(mult,event_weight);
					reco_mult_withdijets_weighted->Fill(recomult, event_weight);
					multiplicity_withdijets_weighted_at1->Fill(multat1, event_weight);
					vzhist_dijet_weighted->Fill(x_vz, event_weight);
					if(colliding_system=="pPb" && year_of_datataking==2016){
						double x3D_hiHF_dijet[3]={hfplus,hfminus,(double) mult}; hfhist_dijet_weighted->Fill(x3D_hiHF_dijet,event_weight);
						double x3D_hiHFEta4_dijet[3]={hfplusEta4,hfminusEta4,(double) mult}; hfhistEta4_dijet_weighted->Fill(x3D_hiHFEta4_dijet,event_weight);
						double x3D_hiZDC_dijet[3]={zdcplus,zdcminus,(double) mult}; zdchist_dijet_weighted->Fill(x3D_hiZDC_dijet,event_weight);
						double x2D_hiHFSum_dijet[2]={hfplus+hfminus,(double) mult}; hfhistSum_dijet_weighted->Fill(x2D_hiHFSum_dijet,event_weight);
						double x2D_hiHFEta4Sum_dijet[2]={hfplusEta4+hfminusEta4,(double) mult}; hfhistEta4Sum_dijet_weighted->Fill(x2D_hiHFEta4Sum_dijet,event_weight);
					}

					if(is_MC && refpt[leadrecojet_index] > 0.0 && refpt[sublrecojet_index] > 0.0){

						// add CM etas here!
						double ref_eta_lead = (double)leadrefmatchjet_eta;  
						double ref_eta_subl = (double)sublrefmatchjet_eta;

						double JES_ratio_reco_vs_ref_leading = leadrecojet_pt/refpt[leadrefjet_index];
						double JES_ratio_reco_vs_ref_subleading = sublrecojet_pt/refpt[sublrefjet_index];
						/*
						if(refpt[sublrefjet_index] > refpt[leadrefjet_index]){
							JES_ratio_reco_vs_ref_leading = leadrecojet_pt/refpt[sublrefjet_index];
							JES_ratio_reco_vs_ref_subleading = sublrecojet_pt/refpt[leadrefjet_index];
						}
						*/
						double x_JES_ratio_reco_vs_ref_leading[6]={JES_ratio_reco_vs_ref_leading,refpt[leadrefjet_index],ref_eta_lead,(double)leadrecojet_flavor,(double)multcentbin,(double) extrabin}; 
						double x_JES_ratio_reco_vs_ref_sleading[6]={JES_ratio_reco_vs_ref_subleading,refpt[sublrefjet_index],ref_eta_subl,(double)sublrecojet_flavor,(double)multcentbin,(double) extrabin}; 

					    hist_leadjes_reco_fromB_weighted->Fill(x_JES_ratio_reco_vs_ref_leading,event_weight*lrefjet_weight*ljet_weight);
						hist_subleadjes_reco_fromB_weighted->Fill(x_JES_ratio_reco_vs_ref_sleading,event_weight*slrefjet_weight*sljet_weight);			 
						if(leadrecojet_index==leadrefjet_index && sublrecojet_index==sublrefjet_index){
							 hist_leadjes_reco_weighted->Fill(x_JES_ratio_reco_vs_ref_leading,event_weight*lrefjet_weight*ljet_weight);
							 hist_subleadjes_reco_weighted->Fill(x_JES_ratio_reco_vs_ref_sleading,event_weight*slrefjet_weight*sljet_weight);
						}
					}	
				}			
			}	
		}
		// -> End of dijet calculations

		// ----> REF <------
		bool removethirdjet_ref = false;
		if(thirdjet_removal_method == 1){ if(thirdrefjet_pt > thirdjet_removal_cut) removethirdjet_ref = true; }
		if(thirdjet_removal_method == 2){ if(thirdrefjet_pt > thirdjet_removal_cut * sublrefjet_pt) removethirdjet_ref = true; }
		if(thirdjet_removal_method == 3){ if(thirdrefjet_pt > thirdjet_removal_cut * 0.5 * (leadrefjet_pt + sublrefjet_pt)) removethirdjet_ref = true; }

		bool removefourjet_ref = false;
		if(do_fourjet_removal == 1 && fourthrefjet_pt > 0.0){ if( fourthrefjet_pt+thirdrefjet_pt > sublrefjet_pt ) removefourjet_ref = true; } 
		if(do_fourjet_removal == 2 && fourthrefjet_pt > 0.0){ if( fourthrefjet_pt+thirdrefjet_pt > subleading_pT_min ) removefourjet_ref = true; } 
		if(do_fourjet_removal == 3 && fourthrefjet_pt > 0.0){ if( fabs(fourthrefjet_pt-thirdrefjet_pt) > 10.0 ) removefourjet_ref = true; } 
		if(do_fourjet_removal == 4 && fourthrefjet_pt > 0.0){ if( fabs(fourthrefjet_pt-thirdrefjet_pt) > 25.0 ) removefourjet_ref = true; } 
		if(do_fourjet_removal == 5 && fourthrefjet_pt > 0.0){ if( fabs(fourthrefjet_pt-thirdrefjet_pt) > subleading_pT_min ) removefourjet_ref = true; } 

		bool tworefjets = twojetfounded(leadrefjet_pt, leading_pT_min, jet_pt_max_cut, sublrefjet_pt, subleading_pT_min);
		bool isrefdijets = isdijet(leadrefjet_pt, leading_pT_min, jet_pt_max_cut, sublrefjet_pt, subleading_pT_min, leadrefjet_phi, sublrefjet_phi, leading_subleading_deltaphi_min, xjmin, xjmax, Ajmin, Ajmax);
		bool goodrefevent = (is_MC && leadrefjet_pt > 0.0 && sublrefjet_pt > 0.0 && !removethirdjet_ref && !removefourjet_ref && removejetid);

		if(goodrefevent){
			//leading/subleading pT cuts
			if(tworefjets){

				// Fill leading/subleading jet quenching quantities
				double delta_phi_ref = fabs(deltaphi(leadrefjet_phi, sublrefjet_phi));
				double Aj_ref = asymmetry(leadrefjet_pt,sublrefjet_pt);
				double Xj_ref = xjvar(leadrefjet_pt,sublrefjet_pt);
				float ptdijetref = 0.5*(leadrefjet_pt + sublrefjet_pt);
				double ptdijetbinref = (double) ptdijetref;
				float dijetetareftype = dijetetabin(leadrefjet_eta, sublrefjet_eta, jet_eta_min_cut, jet_eta_max_cut, jet_fwd_eta_min_cut, jet_fwd_eta_max_cut, jet_bkw_eta_min_cut, jet_bkw_eta_max_cut);

				if(do_dijetstudies){
					double x_ref[9]={Xj_ref,Aj_ref,delta_phi_ref,(double)multcentbin,(double)ptdijetbinref,(double)extrabin,(double)leadrefjet_pt,(double)sublrefjet_pt, (double) dijetetareftype};
					hist_ref_lead_ref_subl_quench_unweighted->Fill(x_ref);	
					hist_ref_lead_ref_subl_quench->Fill(x_ref,event_weight*lrefjet_weight*slrefjet_weight);				
					filletadijethistograms( (double) sqrts, leadrefjet_eta_lab, sublrefjet_eta_lab, leadrefjet_eta, sublrefjet_eta, leadrefjet_pt, leadrefjet_phi, leadrefjet_mass, sublrefjet_pt, sublrefjet_phi, sublrefjet_mass, (double)multcentbin, (double)extrabin, (double) event_weight*lrefjet_weight*slrefjet_weight, hist_etaDijet_ref, hist_etaDijet_CM_ref, hist_yDijet_CM_ref );
					if(colliding_system=="pPb" && year_of_datataking==2016) fillxjEPhistograms_nofake(Xj_ref, delta_phi_ref, leadrefjet_phi, (double) multcentbin, (double)extrabin, (double)ptdijetbinref, (double) dijetetareftype, event_weight*lrefjet_weight*slrefjet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, hist_ref_leadEP_quench_plus, hist_ref_leadEP_quench_minus);
					if(colliding_system=="pPb" && year_of_datataking==2016) fillxjEPhistograms_nofake(Xj_ref, delta_phi_ref, sublrefjet_phi, (double) multcentbin, (double)extrabin, (double)ptdijetbinref, (double) dijetetareftype, event_weight*lrefjet_weight*slrefjet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, hist_ref_sublEP_quench_plus, hist_ref_sublEP_quench_minus);

					/*
					double xjptaveunfref = TransformToUnfoldingAxis_xjptave(Xj_ref, ptdijetbinref,xjmax);						
					double x_unf_meas_xjptaveref[4]={xjptaveunfref, delta_phi_ref, (double)dijetetareftype, (double)multcentbin}; 
					fhUnfoldingTruthRef_xjptave->Fill(x_unf_meas_xjptaveref,event_weight*lrefjet_weight*slrefjet_weight);
					double pt1pt2unfref = TransformToUnfoldingAxis_pt1pt2(leadrefjet_pt, sublrefjet_pt, jet_pt_max_cut);						
					double x_unf_meas_pt1pt2ref[4]={pt1pt2unfref, delta_phi_ref, (double)dijetetareftype, (double)multcentbin}; 
					fhUnfoldingTruthRef_pt1pt2->Fill(x_unf_meas_pt1pt2ref,event_weight*lrefjet_weight*slrefjet_weight);
					double xjpt2unfref = TransformToUnfoldingAxis_xjptave(Xj_ref, sublrefjet_pt,xjmax);						
					double x_unf_meas_xjpt2ref[5]={xjpt2unfref, delta_phi_ref, (double)dijetetareftype, (double)multcentbin, (double)leadrefjet_pt}; 
					fhUnfoldingTruthRef_xjpt2->Fill(x_unf_meas_xjpt2ref,event_weight*lrefjet_weight*slrefjet_weight);
					*/

				}

				if(isrefdijets){

				    // Sort the vector by Pt
				    std::vector<TVector3> sortedVectorRef = alljets_ref;
 	 				std::sort(sortedVectorRef.begin(), sortedVectorRef.end(), sortByPt);
					TVector3 GoodLeadingJets_ref_3vec;
					GoodLeadingJets_ref_3vec.SetPtEtaPhi(leadrefjet_pt,leadrefjet_eta_lab,leadrefjet_phi);
					TVector3 GoodSubLeadingJets_ref_3vec;
					GoodSubLeadingJets_ref_3vec.SetPtEtaPhi(sublrefjet_pt, sublrefjet_eta_lab, sublrefjet_phi);
					TVector3 GoodTrdLeadingJets_ref_3vec;
					if (sortedVectorRef.size() > 2){ GoodTrdLeadingJets_ref_3vec = findClosestVectorInPhi(sortedVectorRef); }else{GoodTrdLeadingJets_ref_3vec.SetPtEtaPhi(thirdrefjet_pt, thirdrefjet_eta, thirdrefjet_phi); }
					TVector3 GoodProjLeadingJets_ref_3vec = -(GoodLeadingJets_ref_3vec + GoodSubLeadingJets_ref_3vec);
					TVector3 GoodProjLeadingJets_ref_3vec_diff = (GoodLeadingJets_ref_3vec - GoodSubLeadingJets_ref_3vec);
					bool remove_3rdjetevent_ref = false;
					if(thirdrefjet_pt > 0 && thirdjet_removal_method == 4 && (fabs(deltaphi(GoodTrdLeadingJets_ref_3vec.Phi(),GoodProjLeadingJets_ref_3vec.Phi())) < thirdjet_removal_cut)) remove_3rdjetevent_ref = true;
			
					if(thirdrefjet_pt > 0){
						double delphi_3rdProjRef = fabs(deltaphi(GoodTrdLeadingJets_ref_3vec.Phi(),GoodProjLeadingJets_ref_3vec.Phi()));
						double deleta_3rdProjRef = deltaeta(GoodTrdLeadingJets_ref_3vec.Eta(),GoodProjLeadingJets_ref_3vec.Eta());
						double xjetprojRef[3] = {delphi_3rdProjRef, deleta_3rdProjRef, (double)multcentbin};
						jetjet_Dphi_Deta_ref->Fill(xjetprojRef,event_weight);
						double delphi_3rdProjRefdiff = fabs(deltaphi(GoodTrdLeadingJets_ref_3vec.Phi(),GoodProjLeadingJets_ref_3vec_diff.Phi()));
						double deleta_3rdProjRefdiff = deltaeta(GoodTrdLeadingJets_ref_3vec.Eta(),GoodProjLeadingJets_ref_3vec_diff.Eta());
						double xjetprojRefdiff[3] = {delphi_3rdProjRefdiff, deleta_3rdProjRefdiff, (double)multcentbin};
						jetjet_Dphi_Deta_ref_diff->Fill(xjetprojRefdiff,event_weight);
					}

						// Fill leading and subleading jet QA histograms
						double x_ref_QA_L[5]={leadrefjet_pt,leadrefjet_eta,leadrefjet_phi,(double)multcentbin,(double) extrabin}; 
						hist_ref_leadjet_weighted->Fill(x_ref_QA_L,event_weight*lrefjet_weight);
						double x_ref_QA_SL[5]={sublrefjet_pt,sublrefjet_eta,sublrefjet_phi,(double)multcentbin,(double) extrabin}; 
						hist_ref_subljet_weighted->Fill(x_ref_QA_SL,event_weight*slrefjet_weight);
						double x_ref_QA_3L[5]={thirdrefjet_pt,thirdrefjet_eta,thirdrefjet_phi,(double)multcentbin,(double) extrabin}; 
						hist_ref_thrdjet_weighted->Fill(x_ref_QA_3L,event_weight);
						double x_ref_QA_3LProj[5]={GoodProjLeadingJets_ref_3vec.Pt(),GoodProjLeadingJets_ref_3vec.Eta(),GoodProjLeadingJets_ref_3vec.Phi(),(double) multcentbin,(double)extrabin}; 
						hist_ref_thrdproj_weighted->Fill(x_ref_QA_3LProj,event_weight);
						double x_ref_QA_3LProjDiff[5]={GoodProjLeadingJets_ref_3vec_diff.Pt(),GoodProjLeadingJets_ref_3vec_diff.Eta(),GoodProjLeadingJets_ref_3vec_diff.Phi(),(double) multcentbin,(double)extrabin}; 
						hist_ref_thrdprojdiff_weighted->Fill(x_ref_QA_3LProjDiff,event_weight);	

				}
			}
		}
	
		// Measure correlations and filling mixing vectors
		// Reco-Reco
		// Inclusive jets
		if(do_inclusejettrack_correlation && pass_Aj_or_Xj_reco_cut && !removethirdjet){
			correlation(jets_reco, jet_w_reco, tracks_reco, track_w_reco, hist_correlation_signal_jet_reco_track_reco, hist_jet_from_reco_reco_sig, hist_trk_from_reco_reco_sig, event_weight, mult, extra_variable, do_rotation, N_of_rot, hist_correlation_rotation_jet_reco_track_reco, sube_tracks_reco, hist_correlation_signal_subg0_jet_reco_track_reco,JetR,hist_injet_reco_track_reco,do_flow); // calculate correlations
			fillvectors(similar_events, Nev_recoreco, jets_reco, jet_w_reco, tracks_reco, track_w_reco, mult, vertexz, event_weight, extra_variable, ev_jet_vector_reco_reco, jet_weights_reco_reco, ev_track_vector_reco_reco, trk_weights_reco_reco, multvec_reco_reco, vzvec_reco_reco, weights_reco_reco,extravec_reco_reco);	// for mixing --> store vectors for mixing
		}

		// Leading/SubLeading jets
		if(do_leading_subleading_jettrack_correlation && pass_Aj_or_Xj_reco_cut && !removethirdjet){
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
				double trk_etamix_weight = 1.0;//get_trketamix_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gtrk_eta, false); // weight to deal with Seagull (test)
				track_w_gen.push_back(trk_etamix_weight); // save weight to apply in the mixing
				
				double trackbingen = (double) gtrk_pt;
				if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filltrkEPhistograms(gtrk_phi, trackbingen, (double) multcentbin, (double) extrabin, event_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_GEN_EP2_flat_trk_plus, Dphi_GEN_EP2_flat_trk_minus, Dphi_GEN_EP3_flat_trk_plus, Dphi_GEN_EP3_flat_trk_minus, Dphi_GEN_EP4_flat_trk_plus, Dphi_GEN_EP4_flat_trk_minus);

			}

			// Start loop over gen jets
			float leadgenjet_pt=-999999999.9, leadgenjet_eta=-999999999.9, leadgenjet_phi=-999999999.9, leadgenjet_mass=-999999999.9, leadgenjet_flavor=-999999999.9; // leading jet quantities
			float sublgenjet_pt=-999999999.9, sublgenjet_eta=-999999999.9, sublgenjet_phi=-999999999.9, sublgenjet_mass=-999999999.9, sublgenjet_flavor=-999999999.9; // subleading jet quantities
			float thirdgenjet_pt=-999999999.9, thirdgenjet_eta=-999999999.9, thirdgenjet_phi=-999999999.9, thirdgenjet_mass=-999999999.9, thirdgenjet_flavor=-999999999.9; // third jet quantities
			int leadgenjet_index=2000, sublgenjet_index=2000, thirdgenjet_index=2000; // jet indexes
			float fourthgenjet_pt=-999999999.9;
			
			bool isgjetincluded = false;

			for(int j = 0; j < gen_jetsize; j++){

				if(fabs(gen_jteta[j]) > 5.1) continue;
				// Define jet kinematics
				float gjet_pt = gen_jtpt[j];
				float gjet_eta = gen_jteta[j];
				float gjet_phi = gen_jtphi[j];
				float gjet_mass = gen_jtmass[j];
				float gjet_flavor = 0.0;
				
				int jet_index_gen = (int) j;

				find_leading_subleading_third(gjet_pt,gjet_eta,gjet_phi,gjet_mass,gjet_flavor,jet_index_gen,leadgenjet_pt,leadgenjet_eta,leadgenjet_phi,leadgenjet_mass,leadgenjet_flavor,leadgenjet_index,sublgenjet_pt,sublgenjet_eta,sublgenjet_phi,sublgenjet_mass,sublgenjet_flavor,sublgenjet_index,thirdgenjet_pt,thirdgenjet_eta,thirdgenjet_phi,thirdgenjet_mass,thirdgenjet_flavor,thirdgenjet_index, fourthgenjet_pt); // Find leading and subleading jets
				gjet_eta = gjet_eta + boost;  // In pPb case, for the center-of-mass correction if needed
				if(colliding_system == "pPb" && is_pgoing && invert_pgoing){gjet_eta = -gjet_eta;}

				double gjet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, gjet_pt, gjet_eta); // Jet weight (specially for MC)

				// Fill gen jet QA histograms
				double x_gen_jet[5]={gjet_pt,gjet_eta,gjet_phi,(double) multcentbin, (double) extrabin}; 
				hist_gen_jet->Fill(x_gen_jet);
				hist_gen_jet_weighted->Fill(x_gen_jet,event_weight*gjet_weight);
				
				TVector3 AllJetsGen;
				AllJetsGen.SetPtEtaPhi(gjet_pt, gjet_eta, gjet_phi);
				alljets_gen.push_back(AllJetsGen);

				if((gjet_pt > jet_pt_min_cut && gjet_pt < jet_pt_max_cut) && (gjet_eta > jet_eta_min_cut && gjet_eta < jet_eta_max_cut)){  // Gen jet pT and eta cut

					isgjetincluded=true;
					// Fill gen jet vectors
					TVector3 GoodJets_gen;
					GoodJets_gen.SetPtEtaPhi(gjet_pt, gjet_eta, gjet_phi);
					jets_gen.push_back(GoodJets_gen);
					jet_w_gen.push_back(gjet_weight);

					if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filljetEPhistograms(gjet_phi, (double) multcentbin, (double) extrabin, event_weight*gjet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_GEN_flat_EP2_inclusive_plus, Dphi_GEN_flat_EP2_inclusive_minus, Dphi_GEN_flat_EP3_inclusive_plus, Dphi_GEN_flat_EP3_inclusive_minus, Dphi_GEN_flat_EP4_inclusive_plus, Dphi_GEN_flat_EP4_inclusive_minus);

				}
			}
			
			if(isgjetincluded){gen_mult_withonejet_weighted->Fill(genmult, event_weight);}

			// ----> GEN <------
			// 3rd jet removal
			bool removethirdjet_gen = false;
			if(thirdjet_removal_method == 1){ if(thirdgenjet_pt > thirdjet_removal_cut) removethirdjet_gen = true; }
			if(thirdjet_removal_method == 2){ if(thirdgenjet_pt > thirdjet_removal_cut * sublgenjet_pt) removethirdjet_gen = true; }
			if(thirdjet_removal_method == 3){ if(thirdgenjet_pt > thirdjet_removal_cut * 0.5 * (leadgenjet_pt + sublgenjet_pt)) removethirdjet_gen = true; }

			// 4th jet removal
			bool removefourjet_gen = false;
			if(do_fourjet_removal == 1 && fourthgenjet_pt > 0.0){ if( fourthgenjet_pt+thirdgenjet_pt > sublgenjet_pt ) removefourjet_gen = true; } 
			if(do_fourjet_removal == 2 && fourthgenjet_pt > 0.0){ if( fourthgenjet_pt+thirdgenjet_pt > subleading_pT_min ) removefourjet_gen = true; } 
			if(do_fourjet_removal == 3 && fourthgenjet_pt > 0.0){ if( fabs(fourthgenjet_pt-thirdgenjet_pt) > 10.0 ) removefourjet_gen = true; } 
			if(do_fourjet_removal == 4 && fourthgenjet_pt > 0.0){ if( fabs(fourthgenjet_pt-thirdgenjet_pt) > 25.0 ) removefourjet_gen = true; } 
			if(do_fourjet_removal == 5 && fourthgenjet_pt > 0.0){ if( fabs(fourthgenjet_pt-thirdgenjet_pt) > subleading_pT_min ) removefourjet_gen = true; } 				

			// weights and dijets
			double gljet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, leadgenjet_pt, leadgenjet_eta); // Jet weight (specially for MC)
			double gsljet_weight = get_jetpT_weight(is_MC, colliding_system.Data(), year_of_datataking, sNN_energy_GeV, sublgenjet_pt, sublgenjet_eta); // Jet weight (specially for MC)
			bool twogenjets = twojetfounded(leadgenjet_pt, leading_pT_min, jet_pt_max_cut, sublgenjet_pt, subleading_pT_min);
			bool isgendijets = isdijet(leadgenjet_pt, leading_pT_min, jet_pt_max_cut, sublgenjet_pt, subleading_pT_min, leadgenjet_phi, sublgenjet_phi, leading_subleading_deltaphi_min, xjmin, xjmax, Ajmin, Ajmax);
			bool goodgenevent = (leadgenjet_pt > 0.0 && sublgenjet_pt > 0.0 && !removethirdjet_gen && !removefourjet_gen && removejetid);

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


			//leading/subleading jets
			if(goodgenevent){
				//leading/subleading pT cuts
				if(twogenjets){ 

					// Fill leading/subleading jet quenching quantities
					double delta_phi_gen = fabs(deltaphi(leadgenjet_phi, sublgenjet_phi));
					double Aj_gen = asymmetry(leadgenjet_pt,sublgenjet_pt);
					double Xj_gen = xjvar(leadgenjet_pt,sublgenjet_pt);
					double ptdijetgen = 0.5*(leadgenjet_pt + sublgenjet_pt);
					double ptdijetbingen = (double) ptdijetgen;
					float dijetetagentype = dijetetabin(leadgenjet_eta, sublgenjet_eta, jet_eta_min_cut, jet_eta_max_cut, jet_fwd_eta_min_cut, jet_fwd_eta_max_cut, jet_bkw_eta_min_cut, jet_bkw_eta_max_cut);
					if(do_dijetstudies){
						double x_gen[9]={Xj_gen,Aj_gen,delta_phi_gen,(double)multcentbin,(double)ptdijetbingen,(double)extrabin,(double)leadgenjet_pt,(double)sublgenjet_pt, (double) dijetetagentype}; 
						hist_gen_lead_gen_subl_quench_unweighted->Fill(x_gen);
						hist_gen_lead_gen_subl_quench->Fill(x_gen,event_weight*gljet_weight*gsljet_weight);
						float correctionmapweightgen = getUnfCorrWeight(fileunfoldweight, leadgenjet_pt, sublgenjet_pt, mult, dijetetagentype);//do correction here
						if(!is_MC) correctionmapweightgen = 1.0;
						hist_gen_lead_gen_subl_quench_checkmapcorrection->Fill(x_gen, correctionmapweightgen*gljet_weight*gsljet_weight);
						filletadijethistograms( (double) sqrts, leadgenjet_eta_lab, sublgenjet_eta_lab, leadgenjet_eta, sublgenjet_eta, leadgenjet_pt, leadgenjet_phi, leadgenjet_mass, sublgenjet_pt, sublgenjet_phi, sublgenjet_mass, (double)multcentbin, (double)extrabin, (double) event_weight*gljet_weight*gsljet_weight, hist_etaDijet_gen, hist_etaDijet_CM_gen, hist_yDijet_CM_gen );
						if(colliding_system=="pPb" && year_of_datataking==2016) fillxjEPhistograms_nofake(Xj_gen, delta_phi_gen, leadgenjet_phi, (double) multcentbin, (double)extrabin, (double)ptdijetbingen, (double) dijetetagentype, event_weight*gljet_weight*gsljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, hist_gen_leadEP_quench_plus, hist_gen_leadEP_quench_minus);
						if(colliding_system=="pPb" && year_of_datataking==2016) fillxjEPhistograms_nofake(Xj_gen, delta_phi_gen, sublgenjet_phi, (double) multcentbin, (double)extrabin, (double)ptdijetbingen, (double) dijetetagentype, event_weight*gljet_weight*gsljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, hist_gen_sublEP_quench_plus, hist_gen_sublEP_quench_minus);
						/*
						double xjptaveunfgen = TransformToUnfoldingAxis_xjptave(Xj_gen, ptdijetbingen,xjmax);						
						double x_unf_meas_xjptavegen[4]={xjptaveunfgen,delta_phi_gen,(double)dijetetagentype,(double)multcentbin}; 
						fhUnfoldingTruthGen_xjptave->Fill(x_unf_meas_xjptavegen,event_weight*gljet_weight*gsljet_weight);
						double pt1pt2unfgen = TransformToUnfoldingAxis_pt1pt2(leadgenjet_pt, sublgenjet_pt, jet_pt_max_cut);						
						double x_unf_meas_pt1pt2gen[4]={pt1pt2unfgen,delta_phi_gen,(double)dijetetagentype,(double)multcentbin}; 
						fhUnfoldingTruthGen_pt1pt2->Fill(x_unf_meas_pt1pt2gen,event_weight*gljet_weight*gsljet_weight);
						double xjpt2unfgen = TransformToUnfoldingAxis_xjptave(Xj_gen, sublgenjet_pt,xjmax);						
						double x_unf_meas_xjpt2gen[5]={xjpt2unfgen,delta_phi_gen,(double)dijetetagentype,(double)multcentbin,(double)leadgenjet_pt}; 
						fhUnfoldingTruthGen_xjpt2->Fill(x_unf_meas_xjpt2gen,event_weight*gljet_weight*gsljet_weight);					
						*/
					}

					// leading/subleading Delta Phi cuts for (leading/subleading)jet+track correlations
					if(isgendijets){
							pass_Aj_or_Xj_gen_cut = true; // if we apply Xj or Aj cuts
							gen_mult_withdijets_weighted->Fill(genmult, event_weight);

						    // Sort the vector by Pt
						    std::vector<TVector3> sortedVectorGen = alljets_gen;
    						std::sort(sortedVectorGen.begin(), sortedVectorGen.end(), sortByPt);
							TVector3 GoodLeadingJets_gen_3vec;
							GoodLeadingJets_gen_3vec.SetPtEtaPhi(leadgenjet_pt,leadgenjet_eta_lab,leadgenjet_phi);
							TVector3 GoodSubLeadingJets_gen_3vec;
							GoodSubLeadingJets_gen_3vec.SetPtEtaPhi(sublgenjet_pt, sublgenjet_eta_lab, sublgenjet_phi);
							TVector3 GoodTrdLeadingJets_gen_3vec;
    						if (sortedVectorGen.size() > 2){ GoodTrdLeadingJets_gen_3vec = findClosestVectorInPhi(sortedVectorGen); }else{ GoodTrdLeadingJets_gen_3vec.SetPtEtaPhi(thirdgenjet_pt, thirdgenjet_eta, thirdgenjet_phi); }
							TVector3 GoodProjLeadingJets_gen_3vec = -(GoodLeadingJets_gen_3vec + GoodSubLeadingJets_gen_3vec);
							TVector3 GoodProjLeadingJets_gen_3vec_diff = (GoodLeadingJets_gen_3vec - GoodSubLeadingJets_gen_3vec);
							bool remove_3rdjetevent_gen = false;
							if(thirdgenjet_pt > 0 && thirdjet_removal_method == 4 && (fabs(deltaphi(GoodTrdLeadingJets_gen_3vec.Phi(),GoodProjLeadingJets_gen_3vec.Phi())) < thirdjet_removal_cut)) remove_3rdjetevent_gen = true;

							if(thirdgenjet_pt > 0){
								double delphi_3rdProjGen = fabs(deltaphi(GoodTrdLeadingJets_gen_3vec.Phi(),GoodProjLeadingJets_gen_3vec.Phi()));
								double deleta_3rdProjGen = deltaeta(GoodTrdLeadingJets_gen_3vec.Eta(),GoodProjLeadingJets_gen_3vec.Eta());
								double xjetprojGen[3] = {delphi_3rdProjGen, deleta_3rdProjGen, (double)multcentbin};
								jetjet_Dphi_Deta_gen->Fill(xjetprojGen,event_weight);
								double delphi_3rdProjGendiff = fabs(deltaphi(GoodTrdLeadingJets_gen_3vec.Phi(),GoodProjLeadingJets_gen_3vec_diff.Phi()));
								double deleta_3rdProjGendiff = deltaeta(GoodTrdLeadingJets_gen_3vec.Eta(),GoodProjLeadingJets_gen_3vec_diff.Eta());
								double xjetprojGendiff[3] = {delphi_3rdProjGendiff, deleta_3rdProjGendiff, (double)multcentbin};
								jetjet_Dphi_Deta_gen_diff->Fill(xjetprojGendiff,event_weight);
							}
							
							
							// Fill leading and subleading jet QA histograms
							double x_lead[5]={leadgenjet_pt,leadgenjet_eta,leadgenjet_phi,(double) multcentbin,(double)extrabin}; 
							hist_gen_leadjet->Fill(x_lead);
							hist_gen_leadjet_weighted->Fill(x_lead,event_weight*gljet_weight);
							double x_sublead[5]={sublgenjet_pt,sublgenjet_eta,sublgenjet_phi,(double) multcentbin,(double)extrabin}; 
							hist_gen_subljet->Fill(x_sublead);
							hist_gen_subljet_weighted->Fill(x_sublead,event_weight*gsljet_weight);
							double x_3rdlead[5]={thirdgenjet_pt,thirdgenjet_eta,thirdgenjet_phi,(double) multcentbin,(double)extrabin}; 
							hist_gen_thrdjet_weighted->Fill(x_3rdlead,event_weight);
							double x_3rdproj[5]={GoodProjLeadingJets_gen_3vec.Pt(),GoodProjLeadingJets_gen_3vec.Eta(),GoodProjLeadingJets_gen_3vec.Phi(),(double) multcentbin,(double)extrabin}; 
							hist_gen_thrdproj_weighted->Fill(x_3rdproj,event_weight);	
							double x_3rdprojdiff[5]={GoodProjLeadingJets_gen_3vec_diff.Pt(),GoodProjLeadingJets_gen_3vec_diff.Eta(),GoodProjLeadingJets_gen_3vec_diff.Phi(),(double) multcentbin,(double)extrabin}; 
							hist_gen_thrdprojdiff_weighted->Fill(x_3rdprojdiff,event_weight);	
													
							bool usethisgenjets = false;
							if(dijetetagentype == 0.5 && fwdbkw_jettrk_option == "mid_mid"){usethisgenjets = true;}
							if(dijetetagentype == 1.5 && fwdbkw_jettrk_option == "mid_fwd"){usethisgenjets = true;}
							if(dijetetagentype == 2.5 && fwdbkw_jettrk_option == "mid_bkw"){usethisgenjets = true;}
							if(dijetetagentype == 3.5 && fwdbkw_jettrk_option == "fwd_mid"){usethisgenjets = true;}
							if(dijetetagentype == 4.5 && fwdbkw_jettrk_option == "bkw_mid"){usethisgenjets = true;}
							if(dijetetagentype == 5.5 && fwdbkw_jettrk_option == "fwd_fwd"){usethisgenjets = true;}
							if(dijetetagentype == 6.5 && fwdbkw_jettrk_option == "fwd_bkw"){usethisgenjets = true;}
							if(dijetetagentype == 7.5 && fwdbkw_jettrk_option == "bkw_fwd"){usethisgenjets = true;}
							if(dijetetagentype == 8.5 && fwdbkw_jettrk_option == "bkw_bkw"){usethisgenjets = true;}

							// Fill leading and subleading jet vectors
							TVector3 GoodLeadingJets_gen;
							TVector3 GoodSubLeadingJets_gen;
							if(usethisgenjets){
								GoodLeadingJets_gen.SetPtEtaPhi(leadgenjet_pt, leadgenjet_eta, leadgenjet_phi);
								lead_jets_gen.push_back(GoodLeadingJets_gen);
								lead_jet_w_gen.push_back(ljet_weight);
								GoodSubLeadingJets_gen.SetPtEtaPhi(sublgenjet_pt, sublgenjet_eta, sublgenjet_phi);
								subl_jets_gen.push_back(GoodSubLeadingJets_gen);
								subl_jet_w_gen.push_back(sljet_weight);
							}
							if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filljetEPhistograms(leadgenjet_phi, (double) multcentbin, (double) extrabin, event_weight*gljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_GEN_flat_EP2_leading_plus, Dphi_GEN_flat_EP2_leading_minus, Dphi_GEN_flat_EP3_leading_plus, Dphi_GEN_flat_EP3_leading_minus, Dphi_GEN_flat_EP4_leading_plus, Dphi_GEN_flat_EP4_leading_minus);
							if(colliding_system=="pPb" && year_of_datataking==2016 && do_flow) filljetEPhistograms(sublgenjet_phi, (double) multcentbin, (double) extrabin, event_weight*gsljet_weight, EP_Psi2_plus_flat, EP_Psi2_minus_flat, EP_Psi3_plus_flat, EP_Psi3_minus_flat, EP_Psi4_plus_flat, EP_Psi4_minus_flat, Dphi_GEN_flat_EP2_subleading_plus, Dphi_GEN_flat_EP2_subleading_minus, Dphi_GEN_flat_EP3_subleading_plus, Dphi_GEN_flat_EP3_subleading_minus, Dphi_GEN_flat_EP4_subleading_plus, Dphi_GEN_flat_EP4_subleading_minus);
					}
				}
			}

			if(do_dijetstudies){
			
				bool missinfo = false;
			
				if(goodrecoevent && tworecojets && removejetid){
					// define variables
					// reco		
					double xjrecoforunfold = xjvar(leadrecojet_pt,sublrecojet_pt);
					double ptaveragerecoforunfold = 0.5*(leadrecojet_pt + sublrecojet_pt);
					double delta_phi_reco_forunfold = fabs(deltaphi(leadrecojet_phi,sublrecojet_phi));
					float dijetetarecotypeforunfold = dijetetabin(leadrecojet_eta, sublrecojet_eta, jet_eta_min_cut, jet_eta_max_cut, jet_fwd_eta_min_cut, jet_fwd_eta_max_cut, jet_bkw_eta_min_cut, jet_bkw_eta_max_cut);
					// ref
					double xjrefforunfold = xjvar(leadrefjet_pt,sublrefjet_pt);
					double ptaveragerefforunfold = 0.5*(leadrefjet_pt + sublrefjet_pt);		
					double delta_phi_ref_forunfold = fabs(deltaphi(leadrefjet_phi,sublrefjet_phi));
					float dijetetareftypeforunfold = dijetetabin(leadrefjet_eta, sublrefjet_eta, jet_eta_min_cut, jet_eta_max_cut, jet_fwd_eta_min_cut, jet_fwd_eta_max_cut, jet_bkw_eta_min_cut, jet_bkw_eta_max_cut);
					// ref matched with reco			
					double xjrefmatchforunfold = xjvar(refpt[leadrecojet_index],refpt[sublrecojet_index]);
					double ptaveragerefmatchforunfold = 0.5*(refpt[leadrecojet_index] + refpt[sublrecojet_index]);		
					double delta_phi_refmatch_forunfold = fabs(deltaphi(refphi[leadrecojet_index],refphi[sublrecojet_index]));
					float dijetetarefmatchtypeforunfold = dijetetabin(leadrefmatchjet_eta, sublrefmatchjet_eta, jet_eta_min_cut, jet_eta_max_cut, jet_fwd_eta_min_cut, jet_fwd_eta_max_cut, jet_bkw_eta_min_cut, jet_bkw_eta_max_cut);
					// matrix reco x ref for leading and subleadings
					double x_unf_meas_ptlead_response[7]={leadrecojet_pt, leadrefjet_pt, delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold, (double)multcentbin}; 
					fhUnfoldingResponse_ptlead->Fill(x_unf_meas_ptlead_response,event_weight*ljet_weight*lrefjet_weight);
					double x_unf_meas_ptsubl_response[7]={sublrecojet_pt, sublrefjet_pt, delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold, (double)multcentbin}; 
					fhUnfoldingResponse_ptsubl->Fill(x_unf_meas_ptsubl_response,event_weight*sljet_weight*slrefjet_weight);
					// matrix reco x ref for Xj
					double x_unf_meas_xj_response[7]={xjrecoforunfold, xjrefforunfold, delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold, (double)multcentbin}; 
					fhUnfoldingResponse_xj->Fill(x_unf_meas_xj_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
					fhUnfoldingResponse_xj_noevtweight->Fill(x_unf_meas_xj_response,ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
					// matrix reco x ref for Xj due corrections			
					float correctionmapweight = getUnfCorrWeight(fileunfoldweight, leadrecojet_pt, sublrecojet_pt, mult, dijetetarecotypeforunfold); //do correction here
					if(!is_MC) correctionmapweight = 1.0;
					
					fhUnfoldingResponse_xj_mapcorrected->Fill(x_unf_meas_xj_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight*correctionmapweight);
					fhUnfoldingResponse_xj_mapcorrected_noevtweight->Fill(x_unf_meas_xj_response,correctionmapweight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);

					// ref matched quantities (matrix + histos)
					if(refpt[leadrecojet_index] > 0.0 && refpt[sublrecojet_index] > 0.0){
						// main observable histogram
						double x_refMatch[9]={xjrefmatchforunfold,xjrefmatchforunfold,delta_phi_refmatch_forunfold,(double)multcentbin,(double) 0.5*(refpt[leadrecojet_index] + refpt[sublrecojet_index]),(double)extrabin,(double)refpt[leadrecojet_index],(double)refpt[sublrecojet_index], (double)dijetetarefmatchtypeforunfold};
						hist_ref_lead_ref_subl_quench_Match_unweighted->Fill(x_refMatch);
						hist_ref_lead_ref_subl_quench_Match->Fill(x_refMatch,event_weight*lrefjet_weight*slrefjet_weight);

						// pt1 && pt2		
						// response
						double x_unf_meas_ptleadMatch_response[7]={leadrecojet_pt, refpt[leadrecojet_index], delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold, (double)multcentbin}; 
						fhUnfoldingResponse_ptleadMatch->Fill(x_unf_meas_ptleadMatch_response,event_weight*ljet_weight*lrefjet_weight);
						double x_unf_meas_ptsublMatch_response[7]={sublrecojet_pt, refpt[sublrecojet_index], delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold, (double)multcentbin}; 
						fhUnfoldingResponse_ptsublMatch->Fill(x_unf_meas_ptsublMatch_response,event_weight*sljet_weight*slrefjet_weight);
				
						//Xj
						// response
						double x_unf_measMatch_xj_response[7]={xjrecoforunfold, xjrefmatchforunfold, delta_phi_reco_forunfold, delta_phi_refmatch_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetarefmatchtypeforunfold, (double)multcentbin}; 
						fhUnfoldingResponse_xjMatch->Fill(x_unf_measMatch_xj_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
						fhUnfoldingResponse_xjMatch_noevtweight->Fill(x_unf_measMatch_xj_response,ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);

						// inverted or subtracted
						if(xjrefmatchforunfold >= 1.0){
				
							double xjrefmatchforunfoldinv = 1.0 / xjrefmatchforunfold;

							double x_refMatchInv[9]={xjrefmatchforunfoldinv, xjrefmatchforunfoldinv, delta_phi_refmatch_forunfold, (double)multcentbin, (double) 0.5*(refpt[leadrecojet_index] + refpt[sublrecojet_index]), (double)extrabin, (double)refpt[sublrecojet_index], (double)refpt[leadrecojet_index], (double)dijetetarefmatchtypeforunfold};
							hist_ref_lead_ref_subl_quench_MatchInv_unweighted->Fill(x_refMatchInv);
							hist_ref_lead_ref_subl_quench_MatchInv->Fill(x_refMatchInv,event_weight*lrefjet_weight*slrefjet_weight);

							double x_unf_measMatchInv_xj_response[7]={xjrecoforunfold, xjrefmatchforunfoldinv, delta_phi_reco_forunfold, delta_phi_refmatch_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetarefmatchtypeforunfold, (double)multcentbin}; 
							fhUnfoldingResponse_xjMatchInv->Fill(x_unf_measMatchInv_xj_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
							fhUnfoldingResponse_xjMatchInv_noevtweight->Fill(x_unf_measMatchInv_xj_response,ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
							fhUnfoldingResponse_xjMatchInv_mapcorrected->Fill(x_unf_measMatchInv_xj_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight*correctionmapweight);
							fhUnfoldingResponse_xjMatchInv_mapcorrected_noevtweight->Fill(x_unf_measMatchInv_xj_response,correctionmapweight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);

						}else{

							double x_refMatchInv[9]={xjrefmatchforunfold, xjrefmatchforunfold, delta_phi_refmatch_forunfold, (double)multcentbin, (double) 0.5*(refpt[leadrecojet_index] + refpt[sublrecojet_index]), (double)extrabin, (double)refpt[leadrecojet_index], (double)refpt[sublrecojet_index], (double)dijetetarefmatchtypeforunfold};
							hist_ref_lead_ref_subl_quench_MatchInv_unweighted->Fill(x_refMatchInv);
							hist_ref_lead_ref_subl_quench_MatchInv->Fill(x_refMatchInv,event_weight*lrefjet_weight*slrefjet_weight);

							double x_unf_measMatchInv_xj_response[7]={xjrecoforunfold, xjrefmatchforunfold, delta_phi_reco_forunfold, delta_phi_refmatch_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetarefmatchtypeforunfold, (double)multcentbin}; 
							fhUnfoldingResponse_xjMatchInv->Fill(x_unf_measMatchInv_xj_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
							fhUnfoldingResponse_xjMatchInv_noevtweight->Fill(x_unf_measMatchInv_xj_response,ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
							fhUnfoldingResponse_xjMatchInv_mapcorrected->Fill(x_unf_measMatchInv_xj_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight*correctionmapweight);
							fhUnfoldingResponse_xjMatchInv_mapcorrected_noevtweight->Fill(x_unf_measMatchInv_xj_response,correctionmapweight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);

						}

					}
/*					
			double xjptaveunfreco_response = TransformToUnfoldingAxis_xjptave(xjrecoforunfold,ptaveragerecoforunfold,xjmax);						
			double xjptaveunfref_response = TransformToUnfoldingAxis_xjptave(xjrefforunfold,ptaveragerefforunfold,xjmax);						
			double x_unf_meas_xjptave_response[7]={xjptaveunfreco_response, xjptaveunfref_response, delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold, (double)multcentbin}; 
			fhUnfoldingResponse_xjptave->Fill(x_unf_meas_xjptave_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);

			double pt1pt2unfreco_response = TransformToUnfoldingAxis_pt1pt2(leadrecojet_pt, sublrecojet_pt, jet_pt_max_cut);						
			double pt1pt2unfref_response = TransformToUnfoldingAxis_pt1pt2(leadrefjet_pt, sublrefjet_pt, jet_pt_max_cut);						
			double x_unf_meas_pt1pt2_response[7]={pt1pt2unfreco_response, pt1pt2unfref_response, delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold,(double)multcentbin}; 
			fhUnfoldingResponse_pt1pt2->Fill(x_unf_meas_pt1pt2_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);

			double xjpt2unfreco_response = TransformToUnfoldingAxis_xjptave(xjrecoforunfold,sublrecojet_pt,xjmax);						
			double xjpt2unfref_response = TransformToUnfoldingAxis_xjptave(xjrefforunfold,sublrefjet_pt,xjmax);						
			double x_unf_meas_xjpt2_response[9]={xjpt2unfreco_response, xjpt2unfref_response, delta_phi_reco_forunfold, delta_phi_ref_forunfold, (double) dijetetarecotypeforunfold, (double)dijetetareftypeforunfold, (double)multcentbin, (double)leadrecojet_pt, (double)leadrefjet_pt}; 
			fhUnfoldingResponse_xjpt2->Fill(x_unf_meas_xjpt2_response,event_weight*ljet_weight*lrefjet_weight*sljet_weight*slrefjet_weight);
*/		
				}
			}

			
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
			// Measure 2 particle correlations and fill mixing vectors for (leading/subleading) track+track flow
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
		if(do_flow){cout << "Running 2PC --> RECO-RECO" << endl; call_mix_random_2pc(N_ev_mix, Mult_or_Cent_range, multvec_reco_reco, multiplicity_centrality_bins, extravec_reco_reco, extra_bins, vzvec_reco_reco, DVz_range, ev_track_vector_reco_reco, trk_weights_reco_reco, hist_reco_reco_2pcorrelation_mixing, trk_pt_bins, weights_reco_reco, double_weight_mix);}
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
	string file_output = Form("%s_%s_%s_%iGeV_%s_%s_%s_Jptmin_%.1f_Jptmax_%.1f_Jetamin_%.1f_Jetamax_%.1f_%s%s%s%s%s_%s_%i",ouputfilename.Data(),colliding_system.Data(),data_or_mc.Data(),sNN_energy_GeV,jet_type.Data(),jet_collection.Data(),jet_trigger.Data(),jet_pt_min_cut,jet_pt_max_cut,jet_eta_min_cut,jet_eta_max_cut,jet_axis.Data(),smear.Data(),XjAj.Data(),ref_sample.Data(),particles.Data(),isflow.Data(),date->GetDate()); // output file
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
	/*
	if(jetid_method == 1 || jetid_method == 2){
		MyFile->mkdir("JetIDVar"); 
		MyFile->cd("JetIDVar"); 	
		w_jetid_hist();
	}
	*/
	if(is_MC && do_dijetstudies){
	
		MyFile->mkdir("Unfolding"); 
		MyFile->cd("Unfolding"); 	
		w_unf_hist();
	
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
