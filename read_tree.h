#include "call_libraries.h"  // call libraries from ROOT and C++

// declare variables

// event quantities
float vertexz; // event z vertex
float hfplus; // event hf + deposit of energy
float hfminus; // event hf - deposit of energy
float zdcplus; // event zdc + deposit of energy
float zdcminus; // event zdc - deposit of energy
int hiBin; // event centrality (used if use_centrality = true in input_variables.h)

//pPb event plane information
float EP_Psi2_plus_flat; // Psi of EP 2 after flattening for HF +
float EP_Psi2_minus_flat; // Psi of EP 2 after flattening for HF -
float EP_Qx2_plus; // Qx of EP 2 no flattening for HF +
float EP_Qx2_minus; // Qx of EP 2 no flattening for HF -
float EP_Qy2_plus; // Qy of EP 2 no flattening for HF +
float EP_Qy2_minus; // Qy of EP 2 no flattening for HF -
float EP_Mult2_plus; // Multiplicity of EP 2 no flattening for HF +
float EP_Mult2_minus; // Multiplicity of EP 2 no flattening for HF -

float EP_Psi3_plus_flat; // Psi of EP 3 after flattening for HF +
float EP_Psi3_minus_flat; // Psi of EP 3 after flattening for HF -
float EP_Qx3_plus; // Qx of EP 3 no flattening for HF +
float EP_Qx3_minus; // Qx of EP 3 no flattening for HF -
float EP_Qy3_plus; // Qy of EP 3 no flattening for HF +
float EP_Qy3_minus; // Qy of EP 3 no flattening for HF -
float EP_Mult3_plus; // Multiplicity of EP 3 no flattening for HF +
float EP_Mult3_minus; // Multiplicity of EP 3 no flattening for HF -

float EP_Psi4_plus_flat; // Psi of EP 4 after flattening for HF +
float EP_Psi4_minus_flat; // Psi of EP 4 after flattening for HF -
float EP_Qx4_plus; // Qx of EP 4 no flattening for HF +
float EP_Qx4_minus; // Qx of EP 4 no flattening for HF -
float EP_Qy4_plus; // Qy of EP 4 no flattening for HF +
float EP_Qy4_minus; // Qy of EP 4 no flattening for HF -
float EP_Mult4_plus; // Multiplicity of EP 4 no flattening for HF +
float EP_Mult4_minus; // Multiplicity of EP 4 no flattening for HF -

// trigger quantities
int jet_trigger_bit; // jet HLT path trigger used for analysis (jet_trigger variable in input_variables.h)

// reco jets
int nref;         // number of jets
float jteta[9999]; // jet eta
float jtphi[9999]; // jet phi
float rawpt[9999]; // jet pT without JEC
float trackMax[9999]; // track maximum pT in a jet
// reco tracks
int ntrk;                 // number of track
float trkpt[9999];       // track pT
float trketa[9999];      // track eta
float trkphi[9999];      // track phi
float trkpterr[9999];    // track pT error (uncertainty)
float trkdcaxy[9999];    // track dxy impact parameter (transverse distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaz[9999];     // track dz impact parameter (longitudinal distance between primary vertex and collision - distance of closest approuch - DCA)
float trkdcaxyerr[9999]; // track dxy error (uncertainty)
float trkdcazerr[9999];  // track dxy error (uncertainty)
float trkchi2[9999];     // track reconstruction chi2 of the fitting
float pfEcal[9999];      // particle flow energy deposit in ECAL
float pfHcal[9999];      // particle flow energy deposit in HCAL
float trkmva[9999];      // track mva for each step
int trkalgo[9999];       // track algorithm/step
unsigned char trkndof[9999];       // track number of degrees of freedom in the fitting 
int trkcharge[9999];     // track charge
unsigned char trknhits[9999];      // number of hits in the tracker
unsigned char trknlayer[9999];     // number of layers with measurement in the tracker
unsigned char trkpixhits[9999];// number of hits in the pixel detector
bool highpur[9999];      // tracker steps MVA selection

// events quantities from gen
float weight; // event weight --> pthat weight
float pthat;  // pthat (initial parton pT)

// gen jets
int ngen;             // number of gen jets
float gen_jtpt[9999];  // gen jet pT
float gen_jteta[9999]; // gen jet eta
float gen_jtphi[9999]; // gen jet phi

// matched jets
float refpt[9999]; // jet pT matched with Gen pT
float refeta[9999]; // jet eta matched with Gen eta
float refphi[9999]; // jet phi matched with Gen phi
int refparton_flavor[9999]; // jet phi matched with Gen phi
int refparton_flavorForB[9999]; // jet phi matched with Gen phi

// gen tracks
std::vector<float> *gen_trkpt = 0;  // gen particle pT
std::vector<float> *gen_trketa = 0; // gen particle eta
std::vector<float> *gen_trkphi = 0; // gen particle phi
std::vector<int> *gen_trkchg = 0;   // gen particle charge
std::vector<int> *gen_trkpid = 0;   // gen particle pid
std::vector<int> *gen_trksube = 0;   // gen particle pid

//All variables listed above are readed in the function bellow
/*
Function to read the Forest/Skim tree
Arguments ->  transfer quantities from trees to our variables
tree: input TChain from jet_analyzer.C file
all the arguments bellow are defined in input_variables.h
is_MC: true -> MC; false -> Data
use_WTA: true -> use WTA (winner-takes-all); false -> use E-Scheme
jet_trigger: string with trigger name
colliding_system: pp, pPb, PbPb, XeXe, ... (future)
colliding_energy: colliding energy in GeV -> 5020, 5440, 8160, 13000, ... 
year_of_datataking: year of data taking
event_filterstr: string of event filters
event_filters: integer (0 or 1) from event filters
*/
void read_tree(TChain *tree, bool is_MC, bool use_WTA, TString jet_trigger, TString colliding_system, int colliding_energy, int year_of_datataking, std::vector<TString> event_filterstr, std::vector<int> event_filters){

    tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files

    // enable branches of interest -> see definition of each variables above

    // event quantities
    tree->SetBranchStatus(Form("%s",jet_trigger.Data()), 1);
    tree->SetBranchStatus("vz", 1);
    
    if(colliding_system=="pPb" && colliding_energy == 8160){

        tree->SetBranchStatus("hiHFplus", 1);
        tree->SetBranchAddress("hiHFplus", &hfplus);
        tree->SetBranchStatus("hiHFminus", 1);
        tree->SetBranchAddress("hiHFminus", &hfminus);
        tree->SetBranchStatus("hiZDCplus", 1);
        tree->SetBranchAddress("hiZDCplus", &zdcplus);
        tree->SetBranchStatus("hiZDCminus", 1);
        tree->SetBranchAddress("hiZDCminus", &zdcminus);

        //EP information

        //angle after the flattening
        tree->SetBranchStatus("epang_HFm2", 1);
        tree->SetBranchAddress("epang_HFm2", &EP_Psi2_minus_flat);
        tree->SetBranchStatus("epang_HFp2", 1);
        tree->SetBranchAddress("epang_HFp2", &EP_Psi2_plus_flat);
        tree->SetBranchStatus("epang_HFm3", 1);
        tree->SetBranchAddress("epang_HFm3", &EP_Psi3_minus_flat);
        tree->SetBranchStatus("epang_HFp3", 1);
        tree->SetBranchAddress("epang_HFp3", &EP_Psi3_plus_flat);
        tree->SetBranchStatus("epang_HFm4", 1);
        tree->SetBranchAddress("epang_HFm4", &EP_Psi4_minus_flat);
        tree->SetBranchStatus("epang_HFp4", 1);
        tree->SetBranchAddress("epang_HFp4", &EP_Psi4_plus_flat);

        //multiplicity
        tree->SetBranchStatus("mult_HFm2", 1);
        tree->SetBranchAddress("mult_HFm2", &EP_Mult2_minus);
        tree->SetBranchStatus("mult_HFp2", 1);
        tree->SetBranchAddress("mult_HFp2", &EP_Mult2_plus);
        tree->SetBranchStatus("mult_HFm3", 1);
        tree->SetBranchAddress("mult_HFm3", &EP_Mult3_minus);
        tree->SetBranchStatus("mult_HFp3", 1);
        tree->SetBranchAddress("mult_HFp3", &EP_Mult3_plus);
        tree->SetBranchStatus("mult_HFm4", 1);
        tree->SetBranchAddress("mult_HFm4", &EP_Mult4_minus);
        tree->SetBranchStatus("mult_HFp4", 1);
        tree->SetBranchAddress("mult_HFp4", &EP_Mult4_plus);

        //Qx
        tree->SetBranchStatus("qx_HFm2", 1);
        tree->SetBranchAddress("qx_HFm2", &EP_Qx2_minus);
        tree->SetBranchStatus("qx_HFp2", 1);
        tree->SetBranchAddress("qx_HFp2", &EP_Qx2_plus);
        tree->SetBranchStatus("qx_HFm3", 1);
        tree->SetBranchAddress("qx_HFm3", &EP_Qx3_minus);
        tree->SetBranchStatus("qx_HFp3", 1);
        tree->SetBranchAddress("qx_HFp3", &EP_Qx3_plus);
        tree->SetBranchStatus("qx_HFm4", 1);
        tree->SetBranchAddress("qx_HFm4", &EP_Qx4_minus);
        tree->SetBranchStatus("qx_HFp4", 1);
        tree->SetBranchAddress("qx_HFp4", &EP_Qx4_plus);

        //Qy
        tree->SetBranchStatus("qx_HFm2", 1);
        tree->SetBranchAddress("qx_HFm2", &EP_Qx2_minus);
        tree->SetBranchStatus("qx_HFp2", 1);
        tree->SetBranchAddress("qx_HFp2", &EP_Qx2_plus);
        tree->SetBranchStatus("qx_HFm3", 1);
        tree->SetBranchAddress("qx_HFm3", &EP_Qx3_minus);
        tree->SetBranchStatus("qx_HFp3", 1);
        tree->SetBranchAddress("qx_HFp3", &EP_Qx3_plus);
        tree->SetBranchStatus("qx_HFm4", 1);
        tree->SetBranchAddress("qx_HFm4", &EP_Qx4_minus);
        tree->SetBranchStatus("qx_HFp4", 1);
        tree->SetBranchAddress("qx_HFp4", &EP_Qx4_plus);

        tree->SetBranchStatus("qy_HFm2", 1);
        tree->SetBranchAddress("qy_HFm2", &EP_Qy2_minus);
        tree->SetBranchStatus("qy_HFp2", 1);
        tree->SetBranchAddress("qy_HFp2", &EP_Qy2_plus);
        tree->SetBranchStatus("qy_HFm3", 1);
        tree->SetBranchAddress("qy_HFm3", &EP_Qy3_minus);
        tree->SetBranchStatus("qy_HFp3", 1);
        tree->SetBranchAddress("qy_HFp3", &EP_Qy3_plus);
        tree->SetBranchStatus("qy_HFm4", 1);
        tree->SetBranchAddress("qy_HFm4", &EP_Qy4_minus);
        tree->SetBranchStatus("qy_HFp4", 1);
        tree->SetBranchAddress("qy_HFp4", &EP_Qy4_plus);
        
    }

    if(colliding_system=="PbPb" || colliding_system=="XeXe") tree->SetBranchStatus("hiBin", 1); //centrality only for PbPb and XeXe
    for(int i = 0; i < event_filterstr.size(); i++) tree->SetBranchStatus(Form("%s",event_filterstr[i].Data()), 1); //event filters

    tree->SetBranchAddress(Form("%s",jet_trigger.Data()), &jet_trigger_bit);
    tree->SetBranchAddress("vz", &vertexz);
    if(colliding_system=="PbPb" || colliding_system=="XeXe") tree->SetBranchAddress("hiBin", &hiBin); //centrality only for PbPb and XeXe
    for(int i = 0; i < event_filters.size(); i++) tree->SetBranchAddress(Form("%s",event_filterstr[i].Data()),&event_filters[i]);

    if(is_MC){
        tree->SetBranchStatus("weight", 1);
        tree->SetBranchStatus("pthat", 1); 
        tree->SetBranchAddress("weight", &weight);
        tree->SetBranchAddress("pthat", &pthat);
    }

    // jet quantities
    tree->SetBranchStatus("nref", 1);
    tree->SetBranchStatus("rawpt", 1);
    tree->SetBranchStatus("trackMax", 1);
    if(use_WTA){
        tree->SetBranchStatus("WTAeta", 1);
        tree->SetBranchStatus("WTAphi", 1);
    }else{
        tree->SetBranchStatus("jteta", 1);
        tree->SetBranchStatus("jtphi", 1);
    }

    tree->SetBranchAddress("nref", &nref);
    tree->SetBranchAddress("rawpt", &rawpt);
    tree->SetBranchAddress("trackMax", &trackMax);
    if(use_WTA){
        tree->SetBranchAddress("WTAeta", &jteta);
        tree->SetBranchAddress("WTAphi", &jtphi);
    }else{
        tree->SetBranchAddress("jteta", &jteta);
        tree->SetBranchAddress("jtphi", &jtphi);
    }

    // gen jet quantities
    if (is_MC){
        tree->SetBranchStatus("ngen", 1);
        tree->SetBranchStatus("genpt", 1);
        if(use_WTA){
            tree->SetBranchStatus("WTAgeneta", 1);
            tree->SetBranchStatus("WTAgenphi", 1);
        }else{
            tree->SetBranchStatus("geneta", 1);
            tree->SetBranchStatus("genphi", 1);
        }

        tree->SetBranchAddress("ngen", &ngen);
        tree->SetBranchAddress("genpt", &gen_jtpt);
        if(use_WTA){
            tree->SetBranchAddress("WTAgeneta", &gen_jteta);
            tree->SetBranchAddress("WTAgenphi", &gen_jtphi);
        }else{
            tree->SetBranchAddress("geneta", &gen_jteta);
            tree->SetBranchAddress("genphi", &gen_jtphi);
        }
    }
    //matching quantities
    if (is_MC){
        tree->SetBranchStatus("refpt", 1);
        tree->SetBranchAddress("refpt", &refpt);
        tree->SetBranchStatus("refeta", 1);
        tree->SetBranchAddress("refeta", &refeta);
        tree->SetBranchStatus("refphi", 1);
        tree->SetBranchAddress("refphi", &refphi);
        tree->SetBranchStatus("refparton_flavor", 1);
        tree->SetBranchAddress("refparton_flavor", &refparton_flavor);
        tree->SetBranchStatus("refparton_flavorForB", 1);
        tree->SetBranchAddress("refparton_flavorForB", &refparton_flavorForB);
    }

    // track quantities
    tree->SetBranchStatus("nTrk", 1);
    tree->SetBranchStatus("trkPt", 1);
    tree->SetBranchStatus("trkEta", 1);
    tree->SetBranchStatus("trkPhi", 1);
    tree->SetBranchStatus("trkPtError", 1);
    tree->SetBranchStatus("trkDxy1", 1);
    tree->SetBranchStatus("trkDxyError1", 1);
    tree->SetBranchStatus("trkDz1", 1);
    tree->SetBranchStatus("trkDzError1", 1);
    if(colliding_system=="PbPb" || colliding_system=="XeXe"){
    	tree->SetBranchStatus("trkChi2", 1);
    	tree->SetBranchStatus("trkNdof", 1);
    	tree->SetBranchStatus("trkNHit", 1);
		tree->SetBranchStatus("trkNlayer", 1);
    }
    if(colliding_system=="pPb" || colliding_system=="XeXe"){
        tree->SetBranchStatus("trkNPixelHit", 1);
    }
    tree->SetBranchStatus("trkCharge", 1);
    tree->SetBranchStatus("highPurity", 1);
    tree->SetBranchStatus("pfEcal", 1);
    tree->SetBranchStatus("pfHcal", 1);

    tree->SetBranchAddress("nTrk", &ntrk);
    tree->SetBranchAddress("trkPt", &trkpt);
    tree->SetBranchAddress("trkEta", &trketa);
    tree->SetBranchAddress("trkPhi", &trkphi);
    tree->SetBranchAddress("trkPtError", &trkpterr);
    tree->SetBranchAddress("trkDxy1", &trkdcaxy);
    tree->SetBranchAddress("trkDxyError1", &trkdcaxyerr);
    tree->SetBranchAddress("trkDz1", &trkdcaz);
    tree->SetBranchAddress("trkDzError1", &trkdcazerr);
    tree->SetBranchAddress("trkCharge", &trkcharge);
    if(colliding_system=="PbPb" || colliding_system=="XeXe"){
    	tree->SetBranchAddress("trkChi2", &trkchi2);
    	tree->SetBranchAddress("trkNdof", &trkndof);
    	tree->SetBranchAddress("trkNHit", &trknhits);
   	tree->SetBranchAddress("trkNlayer", &trknlayer);
    }
    if(colliding_system=="pPb" || colliding_system=="XeXe"){
        tree->SetBranchAddress("trkNPixelHit", &trkpixhits);
    }
    tree->SetBranchAddress("highPurity", &highpur);
    tree->SetBranchAddress("pfEcal", &pfEcal);
    tree->SetBranchAddress("pfHcal", &pfHcal);


    if(colliding_system=="PbPb" && colliding_energy==5020 && year_of_datataking==2018){ //special for 2018 PbPb MC
        tree->SetBranchStatus("trkMVA", 1);
        tree->SetBranchStatus("trkAlgo", 1);
        tree->SetBranchAddress("trkMVA", &trkmva);
        tree->SetBranchAddress("trkAlgo", &trkalgo);
    }

    // gen particle quantities
    if(is_MC){
        tree->SetBranchStatus("pt", 1);
        tree->SetBranchStatus("eta", 1);
        tree->SetBranchStatus("phi", 1);
        tree->SetBranchStatus("chg", 1);
        tree->SetBranchStatus("pdg", 1);
        tree->SetBranchStatus("sube", 1);

        tree->SetBranchAddress("pt", &gen_trkpt);
        tree->SetBranchAddress("eta", &gen_trketa);
        tree->SetBranchAddress("phi", &gen_trkphi);
        tree->SetBranchAddress("chg", &gen_trkchg);
        tree->SetBranchAddress("pdg", &gen_trkpid);
        tree->SetBranchAddress("sube", &gen_trksube);
    }

}
