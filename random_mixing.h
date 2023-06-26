#include "call_libraries.h" // call libraries from ROOT and C++
#include "function_definition.h" // call for functions

/*
Function for random mixing
--> Arguments
ntrkoff_int: range of multiplicity or centrality between first and mixed events
nEvt_to_mix: number of events to mix
ev_ntrkoff: vector with multiplicity or centrality for each event
multiplicity_centrality_bins: multiplicity or centrality bins defined at input_variables.h
vtx_z: vector with z vertex values for each event
vzcut: range of Vz between first and mixed events 
Jet_Vector: vector of vectors with TVector3 (pT, eta, phi) of jets (basically vectors of jets for each event)
Jet_W_Vector: vector of vectors with weights from jets (basically vectors of weights for each jet)
Track_Vector: vector of vectors with TVector3 (pT, eta, phi) of tracks (basically vectors of tracks for each event)
Track_W_Vector: vector of vectors with weights from track (basically vectors of weights for each track)
histo: multidimentional histograms with correlations (DeltaPhi, DeltaEta, TrackBin, MultorCentBin)
trk_pt_bin: track bins defined at input_variables.h
event_weight: event weight vector
histo_jet: multidimentional histograms with jet quantities (Pt, Eta, Phi, MultorCentBin)
histo_trk: multidimentional histograms with track quantities (Pt, Eta, Phi, MultorCentBin)
double_weight: boolean to apply or not double weighting
do_flow: correlation not weighted by pT that can be used for flow analysis
*/
void MixEvents_random(int ntrkoff_int, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> ev_extravar, std::vector<double> extravar_bins, std::vector<double> vtx_z, double vzcut, std::vector<std::vector<TVector3>> Jet_Vector, std::vector<std::vector<double>> Jet_W_Vector, std::vector<std::vector<TVector3>> Track_Vector, std::vector<std::vector<double>> Track_W_Vector, THnSparse *histo, std::vector<double> trk_pt_bin, std::vector<double> event_weight, THnSparse *histo_jet, THnSparse *histo_trk, bool double_weight, bool flow){

   int aux_n_evts = (int)ev_ntrkoff.size(); // total number of events
   //TRandom2 *r = new TRandom2(); // random number producer
   std::vector<int> indexes; // vector to make sure we do not have same combinations

   // first loop over all events to define the jet event
   for(int nevt_trg = 0; nevt_trg < aux_n_evts; nevt_trg++){

      if(nevt_trg != 0 && (nevt_trg % 1000) == 0){double alpha = (double)nevt_trg; cout << " Status: " << std::setprecision(3) << ((alpha / aux_n_evts) * 100) << "%" << endl;}
      auto r = new TRandom2();
      int n_associated = 0; // counter used to find the number to mix 
      int n_associated_check = 0; // counter assure we find the required number of events to mix
      std::vector<int> eventcheck; // vector to make sure we do not have same combinations

      std::vector<TVector3> Jet_nevt_trg_vec = Jet_Vector[nevt_trg]; // jet vector for each event
      std::vector<double> Jet_w_nevt_trg_vec = Jet_W_Vector[nevt_trg]; // jet weight vector for each event
      int nMix_nevt_trg = Jet_nevt_trg_vec.size(); // jet vector size
      std::vector<TVector3> Trk_nevt_trg_vec = Track_Vector[nevt_trg]; // jet vector for each event
      std::vector<double> Trk_w_nevt_trg_vec = Track_W_Vector[nevt_trg]; // jet weight vector for each event
      int nMix_nevt_trg_Trk = Trk_nevt_trg_vec.size(); // jet vector size
      double weight_trgev = event_weight[nevt_trg]; // event weight vector

      // while the number of mixed events is less than the required do this loop
      while (n_associated < nEvt_to_mix){ 
         // additional check, if we do not have nEvt_to_mix after check all events, increase the values
         n_associated_check = n_associated_check + 1;
         if (n_associated_check == (aux_n_evts - 1)*10){
               cout << "Number of events for mixing is not enough in range: [0," << n_associated_check << "]"<< endl; 
               cout << "Increasing the ranges: vz in 0.1 cm and multiplicity or centrality in 1" << endl; 
               vzcut = vzcut + 0.1;
               ntrkoff_int = ntrkoff_int+1;
               n_associated_check = 0;
         }
         int nevt_assoc = r->Integer(aux_n_evts-1); // random events between 0 and aux_n_evts-1
         if(nevt_trg == nevt_assoc) continue; // make sure that we do not choose same events
         if(fabs(ev_ntrkoff[nevt_trg] - ev_ntrkoff[nevt_assoc]) > ntrkoff_int) continue; // multiplicity or centrality matching
         if(fabs(vtx_z[nevt_trg] - vtx_z[nevt_assoc]) > vzcut) continue; // vz matching

         // make sure that we do not duplicate the correlation
         bool isduplicated = false; 
         if(eventcheck.size() > 0){for(int i = 0; i < (int)eventcheck.size(); i++){if(nevt_assoc == eventcheck[i]){isduplicated = true; break;}}}
         if(isduplicated == true) continue;

         // remove repeated combinations: (nevt_trg,nevt_assoc) -> (nevt_assoc,nevt_trg)
         int unique_index = (nevt_trg+nevt_assoc)*nevt_trg+(nevt_trg+nevt_assoc)*nevt_assoc+(nevt_trg*nevt_assoc)*(nevt_trg+nevt_assoc);
         bool check_uniqueindex = false; 
         if(indexes.size()>0){ for(int k = 0; k < indexes.size(); k++){ if(unique_index == indexes[k]){ check_uniqueindex = true; } } }
         if(check_uniqueindex) continue;
         indexes.push_back(unique_index);

         // get weights
         double weight_assev = event_weight[nevt_assoc];
         double f_weight = 1.0;
         if(!double_weight){f_weight =  weight_trgev;}else{f_weight = weight_trgev * weight_assev;} //weighting 
         eventcheck.push_back(nevt_assoc); // save all events that pass the requirements
         n_associated = n_associated + 1; // if pass the requirements sum 1
         std::vector<TVector3> Track_nevt_ass_vec = Track_Vector[nevt_assoc]; // get the vector of tracks
         std::vector<double> Track_w_nevt_ass_vec = Track_W_Vector[nevt_assoc]; // get the vector of tracks      
         int nMix_nevt_ass = (int)Track_nevt_ass_vec.size(); // get the vector of tracks size
         std::vector<TVector3> Jet_nevt_ass_vec = Jet_Vector[nevt_assoc]; // get the vector of tracks
         std::vector<double> Jet_w_nevt_ass_vec = Jet_W_Vector[nevt_assoc]; // get the vector of tracks      
         int nMix_nevt_ass_Jet = (int)Jet_nevt_ass_vec.size(); // get the vector of tracks size

         // loop and fill correlation histograms
         for(int imix = 0; imix < nMix_nevt_trg; imix++){
            double jet_weight = Jet_w_nevt_trg_vec[imix];
            for(int iimix = 0; iimix < nMix_nevt_ass; iimix++){
              double trkpt = Track_nevt_ass_vec[iimix].Pt();
              // track efficiency correction for reco
              double trk_weight = Track_w_nevt_ass_vec[iimix];
              // find track and multiplicity bins
              double trkbin = (double) trkpt;
              double multcentbin = (double) ev_ntrkoff[nevt_trg];
              double extrabin = (double) ev_extravar[nevt_trg];
              // Fill correlation histograms
              double x_jet[5]={Jet_nevt_trg_vec[imix].Pt(),Jet_nevt_trg_vec[imix].Eta(),Jet_nevt_trg_vec[imix].Phi(),(double)multcentbin,(double)extrabin}; histo_jet->Fill(x_jet,jet_weight*f_weight);
              double x_trk[5]={Track_nevt_ass_vec[iimix].Pt(),Track_nevt_ass_vec[iimix].Eta(),Track_nevt_ass_vec[iimix].Phi(),(double)multcentbin,(double)extrabin}; histo_trk->Fill(x_trk,trk_weight*f_weight);
              double del_phi = deltaphi2PC(Jet_nevt_trg_vec[imix].Phi(), Track_nevt_ass_vec[iimix].Phi());
              double del_eta = deltaeta(Jet_nevt_trg_vec[imix].Eta(), Track_nevt_ass_vec[iimix].Eta());
              double xcorr[5]={del_phi,del_eta,(double)trkbin,(double)multcentbin,(double)extrabin}; 
              if(flow){histo->Fill(xcorr,jet_weight*trk_weight*f_weight);}else{histo->Fill(xcorr,jet_weight*trk_weight*f_weight*trkpt);}
            }
         } // end of correlation loop
      } // end of while loop (after finding the number of events required)
      eventcheck.clear();
   } // end of all events loop
}

// Function used only to call the mixing in the main code (Arguments the similar as function above)
void call_mix_random(int nEvt_to_mix, int nmix_int, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> ev_extravar, std::vector<double> extravar_bins, std::vector<double> vtx_z, double vzcut, std::vector<std::vector<TVector3>> Jet_Vector, std::vector<std::vector<double>> Jet_W_Vector, std::vector<std::vector<TVector3>> Track_Vector, std::vector<std::vector<double>> Track_W_Vector, THnSparse *histo,  std::vector<double> trk_pt_bin, std::vector<double> event_weight, THnSparse *histo_jet, THnSparse *histo_trk, bool double_weight, bool flow){
   MixEvents_random(nmix_int, nEvt_to_mix, ev_ntrkoff, multiplicity_centrality_bins, ev_extravar, extravar_bins, vtx_z, vzcut, Jet_Vector, Jet_W_Vector, Track_Vector, Track_W_Vector, histo, trk_pt_bin, event_weight, histo_jet, histo_trk, double_weight, flow);  
}


/*
Function for random mixing for 2pc
--> Arguments
ntrkoff_int: range of multiplicity or centrality between first and mixed events
nEvt_to_mix: number of events to mix
ev_ntrkoff: vector with multiplicity or centrality for each event
multiplicity_centrality_bins: multiplicity or centrality bins defined at input_variables.h
vtx_z: vector with z vertex values for each event
vzcut: range of Vz between first and mixed events 
Track_Vector: vector of vectors with TVector3 (pT, eta, phi) of tracks (basically vectors of tracks for each event)
Track_W_Vector: vector of vectors with weights from track (basically vectors of weights for each track)
histo: multidimentional histograms with correlations (DeltaPhi, DeltaEta, TrackBin, MultorCentBin)
trk_pt_bin: track bins defined at input_variables.h
event_weight: event weight vector
double_weight: boolean to apply or not double weighting
*/
void MixEvents_random_2pc(int ntrkoff_int, int nEvt_to_mix, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> ev_extravar, std::vector<double> extravar_bins, std::vector<double> vtx_z, double vzcut, std::vector<std::vector<TVector3>> Track_Vector, std::vector<std::vector<double>> Track_W_Vector, THnSparse *histo, std::vector<double> trk_pt_bin, std::vector<double> event_weight, bool double_weight){

   int aux_n_evts = (int)ev_ntrkoff.size(); // total number of events
   //TRandom2 *r = new TRandom2(); // random number producer
   std::vector<int> indexes; // vector to make sure we do not have same combinations

   // first loop over all events to define the jet event
   for(int nevt_trg = 0; nevt_trg < aux_n_evts; nevt_trg++){
      
      if(nevt_trg != 0 && (nevt_trg % 1000) == 0){double alpha = (double)nevt_trg; cout << " Status: " << std::setprecision(3) << ((alpha / aux_n_evts) * 100) << "%" << endl;}
      auto r = new TRandom2();
      int n_associated = 0; // counter used to find the number to mix 
      int n_associated_check = 0; // counter assure we find the required number of events to mix
      std::vector<int> eventcheck; // vector to make sure we do not have same combinations

      std::vector<TVector3> Trk_nevt_trg_vec = Track_Vector[nevt_trg]; // jet vector for each event
      std::vector<double> Trk_w_nevt_trg_vec = Track_W_Vector[nevt_trg]; // jet weight vector for each event
      int nMix_nevt_trg = Trk_nevt_trg_vec.size(); // jet vector size
      
      double weight_trgev = event_weight[nevt_trg]; // event weight vector

      // while the number of mixed events is less than the required do this loop
      while (n_associated < nEvt_to_mix){ 
      
         // additional check, if we do not have nEvt_to_mix after check all events, increase the values
         n_associated_check = n_associated_check + 1;
         if (n_associated_check == (aux_n_evts - 1)*10){
               cout << "Number of events for mixing is not enough in range: [0," << n_associated_check << "]"<< endl; 
               cout << "Increasing the ranges: vz in 0.1 cm and multiplicity or centrality in 1" << endl; 
               vzcut = vzcut + 0.1;
               ntrkoff_int = ntrkoff_int+1;
               n_associated_check = 0;
         }
         int nevt_assoc = r->Integer(aux_n_evts-1); // random events between 0 and aux_n_evts-1
         if(nevt_trg == nevt_assoc) continue; // make sure that we do not choose same events
         if(fabs(ev_ntrkoff[nevt_trg] - ev_ntrkoff[nevt_assoc]) > ntrkoff_int) continue; // multiplicity or centrality matching
         if(fabs(vtx_z[nevt_trg] - vtx_z[nevt_assoc]) > vzcut) continue; // vz matching

         // make sure that we do not duplicate the correlation
         bool isduplicated = false; 
         if(eventcheck.size() > 0){for(int i = 0; i < (int)eventcheck.size(); i++){if(nevt_assoc == eventcheck[i]){isduplicated = true; break;}}}
         if(isduplicated == true) continue;

         // remove repeated combinations: (nevt_trg,nevt_assoc) -> (nevt_assoc,nevt_trg)
         int unique_index = (nevt_trg+nevt_assoc)*nevt_trg+(nevt_trg+nevt_assoc)*nevt_assoc+(nevt_trg*nevt_assoc)*(nevt_trg+nevt_assoc);
         bool check_uniqueindex = false; 
         if(indexes.size()>0){ for(int k = 0; k < indexes.size(); k++){ if(unique_index == indexes[k]){ check_uniqueindex = true; } } }
         if(check_uniqueindex) continue;
         indexes.push_back(unique_index);
         
         // get weights
         double weight_assev = event_weight[nevt_assoc];
         double f_weight = 1.0;
         if(!double_weight){f_weight =  weight_trgev;}else{f_weight = weight_trgev * weight_assev;} //weighting 

         eventcheck.push_back(nevt_assoc); // save all events that pass the requirements

         n_associated = n_associated + 1; // if pass the requirements sum 1
         
         std::vector<TVector3> Track_nevt_ass_vec = Track_Vector[nevt_assoc]; // get the vector of tracks
         std::vector<double> Track_w_nevt_ass_vec = Track_W_Vector[nevt_assoc]; // get the vector of tracks      
         int nMix_nevt_ass = (int)Track_nevt_ass_vec.size(); // get the vector of tracks size

         // loop and fill correlation histograms
         for(int imix = 0; imix < nMix_nevt_trg; imix++){
            double trkpt1 = Trk_nevt_trg_vec[imix].Pt();
            double trk_weight1 = Trk_w_nevt_trg_vec[imix];
            int trkbin1 = (int) find_my_bin(trk_pt_bins,trkpt1);	
            for(int iimix = 0; iimix < nMix_nevt_ass; iimix++){
              	double trkpt2 = Track_nevt_ass_vec[iimix].Pt();
              	double trk_weight2 = Track_w_nevt_ass_vec[iimix];
              	int trkbin2 = (int) find_my_bin(trk_pt_bins,trkpt2);
              	if( trkbin1 != trkbin2 ) continue; // only same bin to get vn as sqrt of Vn
              	double trkbin = (double) trkpt1;
              	// track efficiency correction for reco
              	double trk_weight = trk_weight1*trk_weight2;
              	// Find track and multiplicity bins
    	        double multcentbin = (double) ev_ntrkoff[nevt_trg];
	        double extrabin = (double) ev_extravar[nevt_trg];              
              	// Fill correlation histograms
              	double del_phi = deltaphi2PC(Trk_nevt_trg_vec[imix].Phi(), Track_nevt_ass_vec[iimix].Phi());
              	double del_eta = deltaeta(Trk_nevt_trg_vec[imix].Eta(), Track_nevt_ass_vec[iimix].Eta());
             
              	if(del_phi == 0 && del_eta == 0 && trkpt1 == trkpt2) continue; // do not fill histograms if particles are identicles
              	double x2pc[5]={del_phi,del_eta,(double)trkbin,(double)multcentbin,(double)extrabin}; histo->Fill(x2pc,trk_weight*f_weight);
              
            }
         } // end of correlation loop
      } // end of while loop (after finding the number of events required)
      eventcheck.clear();
   } // end of all events loop
}

// Function used only to call the mixing in the main code (Arguments the similar as function above)
void call_mix_random_2pc(int nEvt_to_mix, int nmix_int, std::vector<int> ev_ntrkoff, std::vector<double> multiplicity_centrality_bins, std::vector<double> ev_extravar, std::vector<double> extravar_bins, std::vector<double> vtx_z, double vzcut, std::vector<std::vector<TVector3>> Track_Vector, std::vector<std::vector<double>> Track_W_Vector, THnSparse *histo,  std::vector<double> trk_pt_bin, std::vector<double> event_weight, bool double_weight){
   MixEvents_random_2pc(nmix_int, nEvt_to_mix, ev_ntrkoff, multiplicity_centrality_bins, ev_extravar, extravar_bins, vtx_z, vzcut, Track_Vector, Track_W_Vector, histo, trk_pt_bin, event_weight, double_weight);
}
