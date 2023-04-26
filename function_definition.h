#include "call_libraries.h"  // call libraries from ROOT and C++
#include "trk_efficiency_correction.h" // track efficiency correction
#include "weights.h" // weights applied

/*
Find Ntrk offline -> updated for all systems (and easy to update for future systems)
The Ntrk offline is a definition with specific cuts (we should not change it). The track systematics must be applied using the input_variables.h!
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
hp: track high purity workflow
pterr: track pT uncertainty
dcaxy: track DCA in the transverse plane
dcaxyerr: track DCA in the transverse plane uncertainty
dcaz: track DCA in the longitudinal plane
dcazerr: track DCA in the longitudinal plane uncertainty
chi2: track chi2 of reconstruction
ndof: track number of degrees of freedom reconstruction
nlayer: track number of layers with measurements
nhits: track number of hits with measurements
algo: track MVA algorith step
mva: track MVA algorith value [-1,1]
*/
int get_Ntrkoff(TString col_sys, int col_energy, int yearofdatataking, int size, float *eta, float *pt, int *charge, bool *hp, float *pterr, float *dcaxy, float *dcaxyerr,  float *dcaz, float *dcazerr, float* chi2, unsigned char* ndof, unsigned char* nlayer, unsigned char* nhits, int* algo, float* mva){
	int Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){ 
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(hp[ii] == false) continue;
		if(fabs(pterr[ii]/pt[ii]) >= 0.1) continue;
		if(fabs(dcaxy[ii]/dcaxyerr[ii]) >= 3.0) continue;
		if(fabs(dcaz[ii]/dcazerr[ii]) >= 3.0) continue;
		double calomatching = ((pfEcal[ii]+pfHcal[ii])/cosh(eta[ii]))/pt[ii];
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){if(pt[ii] <= 0.4) continue;}
		if(col_sys=="pp" && col_energy==5020 && yearofdatataking==2017){if(pt[ii] <= 0.5) continue;}
		if(col_sys=="pp" && col_energy==13000 && yearofdatataking==2017){if(pt[ii] <= 0.5) continue;}
		if(col_sys=="XeXe" && col_energy==5440 && yearofdatataking==2017){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/ndof[ii])/nlayer[ii] >= 0.15) continue;
		 	if(nhits[ii] < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
		}
		if(col_sys=="PbPb" && col_energy==5020 && yearofdatataking==2018){
			if(pt[ii] <= 0.5) continue; 
			if((chi2[ii]/ndof[ii])/nlayer[ii] >= 0.18) continue;
		 	if(nhits[ii] < 11) continue;
		 	if(pt[ii] > 20.0){if(calomatching<=0.5)continue;} //is this applicable in pp or pPb?
			if(algo[ii]==6 && mva[ii]<0.98) continue;
		}
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}
/*
Find number of tracks for reco -> updated for all systems (and easy to update for future systems)
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
*/
int get_simple_mult_reco(TString col_sys, int col_energy, int yearofdatataking, int size, float *eta, float *pt, int *charge){
	float Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){
		if(fabs(eta[ii]) > 2.4) continue; 
		if(fabs(charge[ii]) == 0)continue;
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){
			if(pt[ii] <= 0.4) continue;
		}else{
			if(pt[ii] <= 0.5) continue;
		}
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}
/*
Find number of tracks for gen -> updated for all systems (and easy to update for future systems)
--> Arguments
col_sys: colliding system
col_energy:  colliding energy
yearofdatataking: year of data-taking
size: track collection size per event
eta: track eta
pt: track pT
charge: track charge
*/
int get_simple_mult_gen(TString col_sys, int col_energy, int yearofdatataking, int size, std::vector<float> *eta, std::vector<float> *pt, std::vector<int> *charge){
	float Ntrk_off = 0;
	for(int ii=0; ii<size; ii++){
		if(fabs(eta->at(ii)) > 2.4) continue; 
		if(fabs(charge->at(ii)) == 0)continue;
		if(col_sys=="pPb" && col_energy==8160 && yearofdatataking==2016){
			if(pt->at(ii) <= 0.4) continue;
		}else{
			if(pt->at(ii) <= 0.5) continue;
		}
		Ntrk_off=Ntrk_off+1;
	}
	return Ntrk_off;
}

/*
Calculate jet Aj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float asymmetry(float pt_leading, float pt_subleading){
	float Avariable = (pt_leading - pt_subleading)/(pt_leading + pt_subleading);
	return Avariable;
}

/*
Calculate jet Xj asymmetry
--> Arguments
pt_leading: leading jet pT
pt_subleading: subleading jet pT
*/
float xjvar(float pt_leading, float pt_subleading){
	float XJvariable = pt_subleading/pt_leading;
	return XJvariable;
}

/*
Calculate Delta Eta
--> Arguments
eta1: eta of first object
eta2: eta of second object
*/
float deltaeta(float eta1, float eta2){
	float deltaEta = ( eta1 - eta2 );
	return deltaEta;
}

/*
Calculate Delta Phi in the range [-pi,pi]
--> Arguments
phi1: phi of first object
phi2: phi of second object
*/
float deltaphi(float phi1, float phi2){
	float deltaPhi = ( phi1 - phi2 );
    while( deltaPhi >  TMath::Pi() ){deltaPhi  += -2*TMath::Pi();}
   	while( deltaPhi < -TMath::Pi() ){deltaPhi  +=  2*TMath::Pi();}
	return deltaPhi;
}

/*
Calculate Delta Phi in the range [-pi/2 , 3/2 pi]
--> Arguments
phi1: phi of first object
phi2: phi of second object
*/
float deltaphi2PC(float phi1, float phi2){     
	float deltaPhi = (phi1 - phi2);
	while( deltaPhi >  1.5*TMath::Pi() ) deltaPhi += -2.*TMath::Pi();
	while( deltaPhi < -0.5*TMath::Pi() ) deltaPhi +=  2.*TMath::Pi();
	return deltaPhi;
}

/*
Calculate delta R (distance)
--> Arguments
eta1: eta of first object
phi1: phi of first object
eta2: eta of second object
phi2: phi of second object
*/
float deltaR(float eta1, float phi1, float eta2, float phi2){
	float deltaR = sqrt(deltaeta(eta1,eta2)*deltaeta(eta1,eta2) + deltaphi(phi1,phi2)*deltaphi(phi1,phi2));
	return deltaR;
}

/*
Find the leading and subleading jets, return jet leading and subleading pt, eta and phi
--> Arguments
pt: jet pT
eta: jet Eta
phi: jet Phi
leadpt: leading jet pT
leadeta: leading jet Eta
leadphi: leading jet Phi
sublpt: subleading jet pT
subleta: subleading jet Eta
sublphi: subleading jet Phi
*/
void find_leading_subleading(float pt, float eta, float phi, float &leadpt, float &leadeta, float &leadphi, float &sublpt, float &subleta, float &sublphi){
	if( pt > leadpt ) {
    	sublpt = leadpt;
        leadpt = pt;
        leadeta = eta;
        leadphi = phi;
    } else if( sublpt < pt) {
    	sublpt = pt;
        subleta = eta;
        sublphi = phi;
    }
}

/*
Find bin dynamically 
--> Arguments
quant_vec: vector with binning
quant: variable
*/
int find_my_bin(std::vector<double> quant_vec, double quant){
	int bin = -999;
	for(int ii = 0; ii < quant_vec.size()-1; ii++) {if(quant >= quant_vec[ii] && quant < quant_vec[ii+1]){bin = ii;} }
    	return bin;
}

/*
Measure the correlation between objects
--> Arguments
jets: vector with jet informations
jet_w: vector with jet weight informations
tracks: vector with track informations
tracks_w: vector with track weight informations
histo_corr: multidimentional histogram for correlations {Delta Phi, Delta Eta, track pT bin, multiplicity or centrality}
histjet: multidimentional histogram for jets in correlations {pT, Eta, Phi}
histtrk: multidimentional histogram for tracks in correlations {pT, Eta, Phi}
event_weight: event weight vector for each event
mult: multiplicity or centrality vector for each event
do_rotation: true means apply/fill rotation method, otherwise use false
N_rot: number of rotation (only use if do_rotation is true)
histo_rot: histogram using rotation method (only use if do_rotation is true)
sube_trk: vector with sube track (MC embedded samples) sube == 0 means PYTHIA embedded tracks while sube > 0 means the other MC (HYDJET, EPOS, ...)
histo_corr_subeg0: if sube > 0 save in this histogram
flow: true for flow measurement false for jet shapes
*/
void correlation(std::vector<TVector3> jets, std::vector<double> jets_w, std::vector<TVector3> tracks, std::vector<double> tracks_w, THnSparse* histo_corr, THnSparse* histjet, THnSparse* histtrk, float event_weight, int mult, bool do_rotation, int N_rot, THnSparse* histo_rot, std::vector<int> sube_trk, THnSparse* histo_corr_subeg0, bool flow){
	// get correlation histograms
	for (int a = 0; a < jets.size(); a++){ // start loop over jets
        	double jet_weight = jets_w[a];
		for (int b = 0; b < tracks.size(); b++){ // start loop over tracks
			double trkpt = tracks[b].Pt();
			double trketa = tracks[b].Eta();
			int subetrk = sube_trk[b];
            // track efficiency correction for reco
            double trk_weight = tracks_w[b];
           	// Find track and multiplicity bins
			int trkbin = (int) find_my_bin(trk_pt_bins,trkpt);
			int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float) mult);
			// Fill jet and track quantities
			double x4D_jet[4]={jets[a].Pt(),jets[a].Eta(),jets[a].Phi(), (double)multcentbin}; histjet->Fill(x4D_jet,jet_weight*event_weight);
			double x4D_trk[4]={tracks[b].Pt(),tracks[b].Eta(),tracks[b].Phi(), (double)multcentbin}; histtrk->Fill(x4D_trk,trk_weight*event_weight);
			// Fill correlation histograms
			double del_phi = deltaphi2PC(jets[a].Phi(), tracks[b].Phi());
			double del_eta = deltaeta(jets[a].Eta(), tracks[b].Eta());
			double x4D[4]={del_phi,del_eta,(double)trkbin, (double)multcentbin}; 
			if(flow){if(subetrk==0){histo_corr->Fill(x4D,jet_weight*trk_weight*event_weight);}else{histo_corr_subeg0->Fill(x4D,jet_weight*trk_weight*event_weight);}
			}else{if(subetrk==0){histo_corr->Fill(x4D,jet_weight*trk_weight*event_weight*trkpt);}else{histo_corr_subeg0->Fill(x4D,jet_weight*trk_weight*event_weight*trkpt);}}
		}
		// get rotation histograms 
		if(do_rotation){
			for(int c = 0; c < N_rot; c++){
				TRandom2 *r = new TRandom2(); // make a random number
				float alpha = r->Uniform(2.0*TMath::Pi()); // make a random number between 0 and 2pi
				float newphi = jets[a].Phi() + alpha; // add to the jet phi and reorder the axis
				if (newphi > 3.0*TMath::Pi()/2.) newphi -= 2.0*TMath::Pi();
				if (newphi < -TMath::Pi()/2.) newphi += 2.0*TMath::Pi();
				float neweta = -1.0*jets[a].Eta(); // invert the jet eta
				for (int d = 0; d < tracks.size(); d++){ // start loop over tracks using the new jet eta and phi (similar as above)
					double trkpt = tracks[d].Pt();
					double trketa = tracks[d].Eta();
           			// track efficiency correction for reco
           			double trk_weight = tracks_w[d];
            		// Find track and multiplicity bins
					int trkbin = (int) find_my_bin(trk_pt_bins,trkpt);
					int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float)mult);
					// Fill correlation histograms     
					double del_phi_rot = deltaphi2PC(newphi, tracks[d].Phi());
					double del_eta_rot = deltaeta(neweta, tracks[d].Eta());
					double x4D_rot[4]={del_phi_rot,del_eta_rot,(double)trkbin,(double)multcentbin}; 
					if(flow){histo_rot->Fill(x4D_rot,jet_weight*trk_weight*event_weight);}else{histo_rot->Fill(x4D_rot,jet_weight*trk_weight*event_weight*trkpt);}
				}
			}
		}//end rotation
	}//end jet track loop
}

/*
Function to fill vectors to be used during the mixing (the ones with &)
--> Arguments
similar_events: true for using tracks only if event has one jet within the jet cut or falt for all tracks
nev: histogram with number of events stored to be used in the mixing
jets: vector with jet informations
jet_weight: vector with jet weight informations
tracks: vector with track informations
trk_weight: vector with track weight informations
mult: multiplicity or centrality
vertexz: Z vertex in centimeters
weight: event weight
ev_jet_vector: vector to be used in the mixing with jet information for each event
ev_jet_weight_vector: vector to be used in the mixing with jet weight information for each event
ev_track_vector: vector to be used in the mixing with track information for each event
ev_trk_weight_vector: vector to be used in the mixing with track weight information for each event
multvec: vector to be used in the mixing with event multiplicity or centrality information
vzvec: vector to be used in the mixing with event Z vertex position information
weightvec: vector to be used in the mixing with event weight information
*/
void fillvectors(bool similar_events, TH1* nev, std::vector<TVector3> jets, std::vector<double> jet_weight, std::vector<TVector3> tracks, std::vector<double> trk_weight, int mult, double vertexz, double weight, std::vector<std::vector<TVector3>> &ev_jet_vector, std::vector<std::vector<double>> &ev_jet_weight_vector, std::vector<std::vector<TVector3>> &ev_track_vector, std::vector<std::vector<double>> &ev_trk_weight_vector, std::vector<int> &multvec, std::vector<double> &vzvec, std::vector<double> &weightvec){
			if(similar_events){
				if(jets.size() > 0 && tracks.size() > 0){
					nev->Fill(0);
					ev_jet_vector.push_back(jets);
					ev_track_vector.push_back(tracks);
					ev_jet_weight_vector.push_back(jet_weight);
					ev_trk_weight_vector.push_back(trk_weight);
					multvec.push_back(mult);
					vzvec.push_back(vertexz);
					weightvec.push_back(weight);
				}
			}else{
				nev->Fill(0);
				ev_jet_vector.push_back(jets);
				ev_track_vector.push_back(tracks);
				ev_jet_weight_vector.push_back(jet_weight);
				ev_trk_weight_vector.push_back(trk_weight);
				multvec.push_back(mult);
				vzvec.push_back(vertexz);
				weightvec.push_back(weight);
			}
}


/*
Measure the 2 particle correlation
--> Arguments
tracks: vector with track informations
tracks_w: vector with track weight informations
histo_corr: multidimentional histogram for correlations {Delta Phi, Delta Eta, track pT bin, multiplicity or centrality}
event_weight: event weight vector for each event
mult: multiplicity or centrality vector for each event
sube_trk: vector with sube track (MC embedded samples) sube == 0 means PYTHIA embedded tracks while sube > 0 means the other MC (HYDJET, EPOS, ...)
histo_corr_subeg0: if sube > 0 save in this histogram
*/
void twoparticlecorrelation(std::vector<TVector3> tracks, std::vector<double> tracks_w, THnSparse* histo_2pcorr, float event_weight, int mult, std::vector<int> sube_trk, THnSparse* histo_2pcorr_subeg0, THnSparse* histo_2pcorr_subeg0_cross){
	// get correlation histograms
	for (int a = 0; a < tracks.size(); a++){ // start loop over tracks
		double trkpt1 = tracks[a].Pt();
        double trk_weight1 = tracks_w[a];
		int subetrk1 = sube_trk[a];
		int trkbin1 = (int) find_my_bin(trk_pt_bins,trkpt1);

		for (int b = 0; b < tracks.size(); b++){ // start loop over tracks+1
			double trkpt2 = tracks[b].Pt();
            		double trk_weight2 = tracks_w[b];
			int subetrk2 = sube_trk[b];
			int trkbin2 = (int) find_my_bin(trk_pt_bins,trkpt2);
			
			if(trkbin1 != trkbin2) continue; // only same bin to get vn as sqrt of Vn
			
			int trkbin = trkbin1;
            
           	// track efficiency correction for reco
            double trk_weight = trk_weight1*trk_weight2;

            // Find track and multiplicity bins
			int multcentbin = (int) find_my_bin(multiplicity_centrality_bins, (float) mult);

			// Fill correlation histograms
			double del_phi = deltaphi2PC(tracks[a].Phi(), tracks[b].Phi());
			double del_eta = deltaeta(tracks[a].Eta(), tracks[b].Eta());
			
			if(del_phi == 0 && del_eta == 0 && trkpt1 == trkpt2) continue; // do not fill histograms if particles are identical
			
			double x4D[4]={del_phi,del_eta,(double)trkbin,(double)multcentbin}; 
			if(subetrk1==0 && subetrk2==0){histo_2pcorr->Fill(x4D,trk_weight*event_weight);
			}else if(subetrk1>0 && subetrk2>0){histo_2pcorr_subeg0->Fill(x4D,trk_weight*event_weight);
			}else{histo_2pcorr_subeg0_cross->Fill(x4D,trk_weight*event_weight);}
			
		} // b loop
	} // a loop
}

