#include "call_libraries.h"  // important header file

// Call for weights to be applied while running the code

/*
For compatibility between MC RECO and Data
--> Arguments
isMC: true for MC and false for Data
system: colliding system
year: data-taking year
energy: colliding energy
vz: event Vz
weighttree: pt hat weight
leadjetpt: leading jet pT
*/
float get_event_weight(int nevents, bool isMC, bool use_centrality, string system, int year, int energy, float vz, int mult, float weighttree, float leadjetpt, float extraquantity){

	float vzweight = 1.0;
	float multweight = 1.0;
	float evtweight = 1.0;
	float multefficiency = 1.0;
	float jetefficiency = 1.0;		
	float totalweight = 1.0;
	
	// VzWeightFunction is derived from MC vs data event Vz --> MC only --> vzweight
	// MultCentWeightFunction is derived from MC vs data event multiplicity or centrality --> MC only --> multweight
	// MultTriggerWeightFunction is derived from the turn on plots as function of multiplicity --> RECO only
	// JetTriggerWeightFunction is derived from the turn on plots as function of leading jet pT --> RECO only
	// weighttree is the pt hat weight --> MC only
	
	if(isMC && !use_centrality && system == "pp" && energy == 5020 && year == 2017){

		TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol6", -15.0, 15.0);
		VzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10, 3.48615e-09);
		vzweight = VzWeightFunction->Eval(vz);

		TF1 *MultCentWeightFunction = new TF1("MultCentWeightFunction", "pol0", 0.0, 500.0);
		MultCentWeightFunction->SetParameter(0,1.0);
		multweight = MultCentWeightFunction->Eval(mult);

		TF1 *MultTriggerWeightFunction = new TF1("MultTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
		MultTriggerWeightFunction->SetParameter(0,1.0);
		float multtrigweight = 1.0;
		multtrigweight = MultTriggerWeightFunction->Eval(mult);
		multefficiency = 1./multtrigweight;

		TF1 *JetTriggerWeightFunction = new TF1("JetTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
		JetTriggerWeightFunction->SetParameter(0,1.0);
		float jettrigweight = 1.0;
		jettrigweight = JetTriggerWeightFunction->Eval(leadjetpt);
		jetefficiency = 1./jettrigweight;

		evtweight = weighttree;

	}

    if(isMC && !use_centrality && system == "pPb" && energy == 8160 && year == 2016){ // pPb data

		if(leadjetpt > 15.0 && leadjetpt <= 30.){evtweight = 1.0404701e-06 * 961104;}
		else if(leadjetpt > 30. && leadjetpt <= 50.){evtweight = 7.7966624e-08 * 952110;}
		else if(leadjetpt > 50. && leadjetpt <= 80.){evtweight = 1.0016052e-08 * 952554;}
		else if(leadjetpt > 80. && leadjetpt <= 120.){evtweight = 1.3018269e-09 * 996844;}
		else if(leadjetpt > 120. && leadjetpt <= 170.){evtweight = 2.2648493e-10 * 964681;}
		else if(leadjetpt > 170. && leadjetpt <= 220.){evtweight = 4.0879112e-11 * 999260;}
		else if(leadjetpt > 220. && leadjetpt <= 280.){evtweight = 1.1898939e-11 * 964336;}
		else if(leadjetpt > 280. && leadjetpt <= 370.){evtweight = 3.3364433e-12 * 995036;}
		else if(leadjetpt > 370. && leadjetpt <= 460.){evtweight = 7.6612402e-13 * 958160;}
		else if(leadjetpt > 460. && leadjetpt <= 540.){evtweight = 2.1341026e-13 * 981427;}
		else if(leadjetpt > 540.){evtweight = 7.9191586e-14 * 1000000;}
		evtweight = (float) evtweight / nevents;

		// multiplicity weight
		// PYTHIA+EPOS
		/*
    	if(mult < 210.0){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol8", 0, 214);
			EtWeightFunction->SetParameters(-0.0915839,-0.00212374,0.00242633,-8.75813e-05,1.48422e-06,-1.3365e-08,6.57739e-11,-1.67841e-13,1.74335e-16);
			multweight = EtWeightFunction->Eval(mult);
		}else if(mult >= 210 && mult <= 300){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol4",210,320);
			EtWeightFunction->SetParameters(-6.2635,0.135792,-0.000754343,1.4557e-06,-7.25631e-10);
			multweight = EtWeightFunction->Eval(mult);					
		}else{multweight = 1.0;}
		TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1);
		VzWeightFunction->SetParameters(0.843091,-0.0145326,0.00468752,-0.000202477,3.39591e-05,4.51739e-07,-8.6924e-08,-7.52764e-09,1.1983e-09)
		vzweight = VzWeightFunction->Eval(vz);
		*/
		// PYTHIA only
		/*
    	if(mult < 105.0){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol8", 0,109);
			EtWeightFunction->SetParameters(-3.3272,0.641081,-0.0356275,0.00120248,-2.63759e-05,3.62415e-07,-2.97515e-09,1.33199e-11,-2.50344e-14);
			multweight = EtWeightFunction->Eval(mult);
		}else if(mult >= 105.0 && mult <= 250.0){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "expo",105,250);
			EtWeightFunction->SetParameters(9.99045,-0.13026);
			multweight = EtWeightFunction->Eval(mult);					
		}else{multweight = 1.0;}
		TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1);
		VzWeightFunction->SetParameters(0.848574,-0.0182413,0.00402673,8.70187e-06,8.54011e-06,-1.76185e-06,1.39702e-07,2.28844e-09,1.01863e-10);
		vzweight = VzWeightFunction->Eval(vz);
		*/
		// EPb weight
		// PYTHIA+EPOS
		/*
	   	if(extraquantity < 65.0){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol8",0,70);
			EtWeightFunction->SetParameters(0.173544,0.0295155,0.00574286,3.39713e-05,-2.61854e-05,1.12039e-06,-2.08602e-08,1.86687e-10,-6.57817e-13);
			multweight = EtWeightFunction->Eval(extraquantity);
		}else if(extraquantity >= 65.0 && extraquantity < 95.0){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol2",60,120);
			EtWeightFunction->SetParameters(2.56581,-0.0458321,0.000231857);
			multweight = EtWeightFunction->Eval(extraquantity);					
		}else{
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol0",90,200);
			EtWeightFunction->SetParameters(0.29061);	
			multweight = EtWeightFunction->Eval(extraquantity);			
		}
		TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1);
		VzWeightFunction->SetParameters(0.853343,-0.0157977,0.00431199,-0.000162212,3.61976e-05,-1.21137e-07,-9.833e-08,-4.70594e-09,1.05049e-09);
		vzweight = VzWeightFunction->Eval(vz);
		*/
  		// PYTHIA only
  		/*
	   	if(extraquantity < 75.0){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol7",3,80);
			EtWeightFunction->SetParameters(1.45325,0.0357072,-0.011459,0.00062546,-1.64134e-05,2.32377e-07,-1.70596e-09,5.09787e-12);
			multweight = EtWeightFunction->Eval(extraquantity);
		}else{
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol0",75,200);
			EtWeightFunction->SetParameters(0.0727252);
			multweight = EtWeightFunction->Eval(extraquantity);			
		}
		TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1);
		VzWeightFunction->SetParameters(0.852218,-0.0164135,0.0047504,-9.06856e-05,2.02e-05,-7.25535e-07,6.2906e-08,-2.08797e-09,5.13e-10);
		vzweight = VzWeightFunction->Eval(vz);
		*/
		multweight = 1./multweight;
		vzweight = 1./vzweight;
    }
/* --> To recover EPb in MB compared to HM triggers --> not useful, bias still there
    if(!isMC && !use_centrality && system == "pPb" && energy == 8160 && year == 2016){
    	if(mult <= 100.0){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol1", 70.0, 140.0);
			EtWeightFunction->SetParameters(8.15578e-01, -4.86085e-03);
			multweight = EtWeightFunction->Eval(extraquantity);
		}else if(mult > 100. && mult <= 185.){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol4", 70.0, 140.0);
			EtWeightFunction->SetParameters(4.08753e+01, -1.56574e+00, 2.31274e-02, -1.48438e-04, 3.51648e-07);
			multweight = EtWeightFunction->Eval(extraquantity);		
		}else if(mult > 185. && mult <= 250.){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol2", 70.0, 140.0);
			EtWeightFunction->SetParameters(2.19684e+01, -3.16469e-01, 1.69276e-03);
			multweight = EtWeightFunction->Eval(extraquantity);			
		}else if(mult > 250.){
			TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol3", 70.0, 140.0);
			EtWeightFunction->SetParameters(-2.66599e+02, 7.95438e+00, -7.28621e-02, 2.26233e-04);
			multweight = EtWeightFunction->Eval(extraquantity);					
		}
    }
*/
	totalweight = evtweight*multweight*vzweight*multefficiency*jetefficiency;
	return totalweight;

}

/*
For compatibility between MC RECO and Data
--> Arguments
isMC: true for MC and false for Data
system: colliding system
year: data-taking year
energy: colliding energy
jetpt: jet pT weight
*/
float get_jetpT_weight(bool isMC, string system, int year, int energy, float jetpt){

	float jetptweight = 1.0;

	// JetPtWeightFunction is derived from MC vs data jet pT spectra.
/*
	if(isMC && system == "pp" && energy == 5020 && year == 2017){
		TF1 *JetPtWeightFunction = new TF1("JetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from all jets above 120 GeV and JECv6
	    JetPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09);
		jetptweight = JetPtWeightFunction->Eval(jetpt);
	}
*/
	return jetptweight;

}

/*
For compatibility between MC RECO and Data
--> Arguments
isMC: true for MC and false for Data
system: colliding system
year: data-taking year
energy: colliding energy
leadjetpt: leading jet pT weight
*/
float get_leadjetpT_weight(bool isMC, string system, int year, int energy, float leadjetpt){

	float leadjetptweight = 1.0;

	// LeadJetPtWeightFunction is derived from MC vs data leading jet pT spectra.
/*
	if(isMC && system == "pp" && energy == 5020 && year == 2017){
		TF1 *LeadJetPtWeightFunction = new TF1("LeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from leading jets above 120 GeV and JECv6
	    LeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
		leadjetptweight = LeadJetPtWeightFunction->Eval(leadjetpt);
	}
*/
	return leadjetptweight;

}

/*
For compatibility between MC RECO and Data
--> Arguments
isMC: true for MC and false for Data
system: colliding system
year: data-taking year
energy: colliding energy
subleadjetpt: subleading jet pT weight
*/
float get_subleadjetpT_weight(bool isMC, string system, int year, int energy, float subleadjetpt){

	float subleadjetptweight = 1.0;

	// SubLeadJetPtWeightFunction is derived from MC vs data subleading jet pT spectra.
/*
	if(isMC && system == "pp" && energy == 5020 && year == 2017){
		TF1 *SubLeadJetPtWeightFunction = new TF1("SubLeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from leading jets above 120 GeV and JECv6
	    SubLeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
		subleadjetptweight = SubLeadJetPtWeightFunction->Eval(subleadjetpt);
	}
*/
	return subleadjetptweight;

}
