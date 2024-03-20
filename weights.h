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
float get_event_weight(float nevents, bool isMC, bool use_centrality, string system, int year, int energy, float vz, int mult, float weighttree, float leadjetpt, float extraquantity, bool is_embedded, bool is_multdep){

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
		else if(leadjetpt > 120.&& leadjetpt <= 170.){evtweight = 2.2648493e-10 * 964681;}
		else if(leadjetpt > 170. && leadjetpt <= 220.){evtweight = 4.0879112e-11 * 999260;}
		else if(leadjetpt > 220. && leadjetpt <= 280.){evtweight = 1.1898939e-11 * 964336;}
		else if(leadjetpt > 280. && leadjetpt <= 370.){evtweight = 3.3364433e-12 * 995036;}
		else if(leadjetpt > 370. && leadjetpt <= 460.){evtweight = 7.6612402e-13 * 958160;}
		else if(leadjetpt > 460. && leadjetpt <= 540.){evtweight = 2.1341026e-13 * 981427;}
		else if(leadjetpt > 540.){evtweight = 7.9191586e-14 * 1000000;}
		evtweight = (float) evtweight/nevents;
		
		// Vz weighting
        TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol8", -15.1, 15.1);
        VzWeightFunction->SetParameters(0.856516,-0.0159813,0.00436628,-0.00012862,2.61129e-05,-4.16965e-07,1.73711e-08,-3.11953e-09,6.24993e-10);
        vzweight = VzWeightFunction->Eval(vz);
        vzweight = 1./vzweight;
		
		// multiplicity weight
		// PYTHIA+EPOS
        if(is_embedded && is_multdep){
           if(mult < 205){
              TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol8", 10, 210);
              EtWeightFunction->SetParameters(-0.0958084,-0.00115779,0.00233641,-8.42997e-05,1.42306e-06,-1.27251e-08,6.20308e-11,-1.56457e-13,1.60319e-16);
              multweight = EtWeightFunction->Eval(mult);
           }else if(mult >= 205 && mult <= 280){
              TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol1",200,280);
              EtWeightFunction->SetParameters(4.24931,-0.0149767);
              multweight = EtWeightFunction->Eval(mult);
           }else{multweight = 1.0;}
        }
        // PYTHIA only
        if(!is_embedded && is_multdep){
           if(mult < 105){
              TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol8", 10,110);
              EtWeightFunction->SetParameters(-3.15516,0.608522,-0.0334141,0.00112313,-2.46795e-05,3.40076e-07,-2.79767e-09,1.25403e-11,-2.35808e-14);
              multweight = EtWeightFunction->Eval(mult);
           }else if(mult >= 105.0 && mult <= 250.0){
              TF1 *EtWeightFunction = new TF1("EtWeightFunction", "expo",100,250);
              EtWeightFunction->SetParameters(9.97281,-0.13071);
              multweight = EtWeightFunction->Eval(mult);
           }else{multweight = 1.0;}
        }
		/*
		// EPb weight
		// PYTHIA+EPOS
		if(is_embedded && !is_multdep){
		   	if(extraquantity < 65.0){
				TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol8",0,70);
				EtWeightFunction->SetParameters(0.176622,0.0267885,0.00621628,-1.08812e-05,-2.36853e-05,1.04298e-06,-1.95357e-08,1.74939e-10,-6.1565e-13);
				multweight = EtWeightFunction->Eval(extraquantity);
			}else if(extraquantity >= 65.0 && extraquantity < 95.0){
				TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol2",60,120);
				EtWeightFunction->SetParameters(2.71433,-0.0485389,0.000250576);
				multweight = EtWeightFunction->Eval(extraquantity);					
			}else{
				TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol0",90,200);
				EtWeightFunction->SetParameter(0, 0.362791);	
				multweight = EtWeightFunction->Eval(extraquantity);			
			}
		}
  		// PYTHIA only
		if(!is_embedded && !is_multdep){
	   		if(extraquantity < 75.0){
				TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol7",0,80);
				EtWeightFunction->SetParameters(1.43436,0.0374677,-0.0114677,0.000621941,-1.62494e-05,2.29115e-07,-1.67501e-09,4.98366e-12);
				multweight = EtWeightFunction->Eval(extraquantity);
			}else{
				TF1 *EtWeightFunction = new TF1("EtWeightFunction", "pol0",75,200);
				EtWeightFunction->SetParameter(0, 0.0866719);
				multweight = EtWeightFunction->Eval(extraquantity);			
			}
		}
		*/
		multweight = 1./multweight;
    }
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
float get_jetpT_weight(bool isMC, string system, int year, int energy, float jetpt, float jeteta){

	float jetptweight = 1.0;

	// JetPtWeightFunction is derived from MC vs data jet pT spectra for pp at 5.02 TeV
	if(isMC && system == "pp" && energy == 5020 && year == 2017){
		TF1 *JetPtWeightFunction = new TF1("JetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from all jets above 120 GeV and JECv6
	    JetPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09);
		jetptweight = JetPtWeightFunction->Eval(jetpt);
	}

	return jetptweight;

}
