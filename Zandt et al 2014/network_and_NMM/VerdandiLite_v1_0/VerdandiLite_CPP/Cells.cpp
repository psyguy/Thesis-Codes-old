#include "Cells.h"
#include "Synapses.h"
#include "Network.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using namespace std;

#define DEBUG_RUN 0

Neuron::Neuron(int np, int ns): DynSys(np,ns) {
	Id.LocalID = Network::AddNeuron(this);
}

void Neuron::Spike() {	
	#if DEBUG_RUN > 0
	cout << "  Spike from cell " << Id.LocalID << endl;
	#endif
	Network::AddSpike(Id);
}

void Neuron::PreTimeStep() {
	#if DEBUG_RUN > 3
	cout << "      Neuron: determining inputs" << endl;
	#endif

	// Determine input from all synapses:
	PulsesIn = 0;
	I_external = 0;

	int nSynapse = SynapseList.size();


	for(int iSynapse=0; iSynapse<nSynapse; iSynapse++) {
		#if DEBUG_RUN > 4
		cout << "        Checking synapse " << iSynapse << " of " << nSynapse << endl;
		#endif
		PulsesIn += SynapseList.at(iSynapse)->GetPulse(State[0]);
		I_external += SynapseList.at(iSynapse)->GetCurrent(State[0]);
	}
			
	State[0]+=PulsesIn;
}

void Neuron::AddSynapse(Synapse* s) {
	SynapseList.push_back(s);
	//#if DEBUG_INPUT > 3
	//cout << "      Synapse " << s << " added to cell " << Id.LocalID << endl;
	//#endif
}


NeuronTypeLIF::NeuronTypeLIF(): Neuron(4,1) {
/*
	ParamId | Description
	--------+-----------
	   0    | C
	   1    | V_rest
	   2    | V_thres
	   3    | V_reset

	StateId | Description
	--------+----------
	   0    | V
*/
	MirrorVm = 0;
}

void NeuronTypeLIF::VectorField(double* dState, double* StateIn) {
	// dV = -(V-Vrest)/tau
	dState[0] = (I_external-(StateIn[0]-Param[1]))/Param[0];
}

void NeuronTypeLIF::PostTimeStep() {
	// Check for threshold passing:
	if(State[0] > Param[2]) {
		State[0] = Param[3];
		MirrorVm = 30;
		Spike();
	} else {
		MirrorVm = State[0];
	}
}

double NeuronTypeLIF::GetVm() {
	return MirrorVm;
}



NeuronTypeIzhikevich::NeuronTypeIzhikevich(): Neuron(5,2) {
/*
	ParamId | Description
	--------+-----------
	   0    | Injected current (const)
	   1    | a
	   2    | b
	   3    | c
	   4    | d

	StateId | Description
	--------+-----------
	   0    | V
	   1    | u

*/
}

void NeuronTypeIzhikevich::VectorField(double* dState, double* StateIn) {
	// dV = 0.04v^2 + 5v + 140 - u + I
	// du = a(bv-u)
	dState[0] = 0.04*StateIn[0]*StateIn[0] + 5*StateIn[0] + 140 - StateIn[1] + Param[0] + I_external;
	dState[1] = Param[1]*(Param[2]*StateIn[0]-StateIn[1]);
}

void NeuronTypeIzhikevich::PostTimeStep() {
	// if(Vm>=30)
	if(State[0] >= 30) {
		State[0] = Param[3];
		State[1] = State[1]+Param[4];
		Spike();
	}
}


NeuronTypePoisson::NeuronTypePoisson(): Neuron(1,1) {
/*
	ParamId | Description
	--------+-----------
	   0    | Poisson rate (s^-1)

	StateId | Description
	--------+-----------
	   0    | V
*/
}

void NeuronTypePoisson::VectorField(double* dState, double* StateIn) {
	dState[0]=0;
}

void NeuronTypePoisson::PreTimeStep() {
	double r = ((double)rand())/RAND_MAX;
	if(r<DynSys::dt*Param[0]/1000) { // Factor 1000 corrects for rate expressed in s^-1, while time is measured in ms by default
		State[0] = 1;
	} else {
		State[0] = 0;
	}
}

void NeuronTypePoisson::PostTimeStep() {
	if(State[0]>0) Spike();
}

NeuronTypePoissonStep::NeuronTypePoissonStep(): Neuron(3,2) {
/* Same as poisson neuron, but changes rate after time Param[2]
	ParamId | Description
	--------+-----------
	   0    | initial Poisson rate (s^-1)
	   1    | after step Poisson rate
	   2    | Ts step time (s)

	StateId | Description
	--------+-----------
	   0    | V
	   1    | t (time)
*/
}

void NeuronTypePoissonStep::VectorField(double* dState, double* StateIn) {
	dState[0]=0;
	dState[1]=1;
}
void NeuronTypePoissonStep::PreTimeStep() { 
    double r = ((double)rand())/RAND_MAX;
	double rate;		
    if(State[1] < Param[2]*1000) {
         rate = Param[0];
    } else {
         rate = Param[1];
//         cout<<State[1]<<"   "<<Param[1]<<endl;
    }
    if(r<DynSys::dt*rate/1000) { // Factor 1000 corrects for rate expressed in s^-1, while time is measured in ms by default
		State[0] = 1;
	} else {
		State[0] = 0;
	}
}

void NeuronTypePoissonStep::PostTimeStep() {
	if(State[0]>0) Spike();
}



NeuronTypeHH::NeuronTypeHH(): Neuron(13,3) {
/*
	ParamId | Description
	--------+-----------
	   0    | Iinput           % [uF /cm^2]    Injected current (const) 
	   1    | Cm               % [uF / cm^2],  membrane capacitance
	   2    | g_na             % [mS / cm^2],  maximum gate conductances
	   3    | g_naL            % [mS / cm^2],  
	   4    | g_k
	   5    | g_kL
	   6    | g_clL
	   7    | phi              % [1/ms] gate time constant
	   8    | E_k              % [mV]   Nernst potentials 
       9    | E_na
	   10   | E_cl
	   11   | gsyn             % constant synaptic conductance, used for testing
	   12   | Esyn             % reversal potential of constant synapse


	StateId | Description
	--------+-----------
	   0    | V
	   1    | n
	   2    | h

*/
SpikedRecently = 0;
}

#define Iinput Param[0]
#define Cm     Param[1]
#define g_na   Param[2]
#define g_naL  Param[3]
#define g_k    Param[4]
#define g_kL   Param[5]
#define g_clL  Param[6]
#define phi    Param[7]
#define E_k    Param[8]
#define E_na   Param[9]
#define E_cl   Param[10]
//#define SpikeThreshold   Param[11]

void NeuronTypeHH::VectorField(double* dState, double* StateIn) {
    double alpha_n;
    double alpha_m;
    if (fabs(StateIn[0] + 34.0) > 0.001)    {
        alpha_n =0.01 * (StateIn[0]+34.0)/( 1.0 - exp(-0.1 * (StateIn[0]+34.0)) );
    } else {
        alpha_n = 0.1;
    }    
    double beta_n = 0.125 * exp(-(StateIn[0]+44.0)/80.0);
    if (fabs(StateIn[0] + 30.0) > 0.001)    {
        alpha_m = 0.1 * (StateIn[0]+30.0)/( 1.0 - exp(-0.1 * (StateIn[0]+30.0)) );
    } else {
        alpha_m = 1.0;
    } 
	double beta_m = 4.0 * exp(-(StateIn[0]+55.0)/18.0);
	double alpha_h = 0.07 * exp(-(StateIn[0]+44.0)/20.0);
	double beta_h = 1.0/( 1.0 + exp(-0.1 * (StateIn[0]+14.0)) );
	double m_inf = alpha_m/(alpha_m + beta_m);

	double Ina = g_na*(m_inf*m_inf*m_inf)*StateIn[2]*(StateIn[0]-E_na) + g_naL*(StateIn[0]-E_na);              // [mS/cm^2 * mV = uA/cm^2], currents per membrane area
    double Ik =  (g_k*StateIn[1]*StateIn[1]*StateIn[1]*StateIn[1])*(StateIn[0]-E_k) + g_kL*(StateIn[0]-E_k);
    double Icl = g_clL*(StateIn[0]-E_cl);
	
	double Iconstsyn = Param[11]*(Param[12]-State[0]);
	
	dState[0] = (1.0/Cm)*(-Ina -Ik -Icl +Iinput + I_external + Iconstsyn);
	//cout << I_external << "    " << (Ina +Ik +Icl) << endl;
    dState[1] = phi*(alpha_n*(1-StateIn[1])-beta_n*StateIn[1]);
	dState[2] = phi*(alpha_h*(1-StateIn[2])-beta_h*StateIn[2]);
	
}
//void NeuronTypeHH::PreTimeStep() {
//      if(State[0] < 10) {
//		UnderThresholdbeforeStep = true;
//	} else {
//		UnderThresholdbeforeStep = false;
//	}
//		#if DEBUG_RUN > 3
//	cout << "      Neuron: determining inputs" << endl;
//	#endif

	// Determine input from all synapses:
//	PulsesIn = 0;
//	I_external = 0;

//	int nSynapse = SynapseList.size();

	//for(int iSynapse=0; iSynapse<nSynapse; iSynapse++) {
	//	#if DEBUG_RUN > 4
		//cout << "        Checking synapse " << iSynapse << " of " << nSynapse << endl;
		//#endif
		//PulsesIn += SynapseList.at(iSynapse)->GetPulse(State[0]);
		//I_external += SynapseList.at(iSynapse)->GetCurrent(State[0]);
	//}
			
//	State[0]+=PulsesIn;
//}
void NeuronTypeHH::PostTimeStep() {
     if(State[0] < -10){SpikedRecently = 0;}  //oscillation around zero only are denoted as spikes when they fall below -10 mV
     if(State[0] > 0) {
                 if(SpikedRecently==0) {
                                       Spike();
                                       SpikedRecently = 1;
                 }
     }       

}



NeuronTypeTest::NeuronTypeTest(): Neuron(4,1) {
/*
	ParamId | Description
	--------+-----------
	   0    | C
	   1    | V_rest
	   2    | V_thres
	   3    | V_reset

	StateId | Description
	--------+----------
	   0    | V
*/
	MirrorVm = 0;
}

void NeuronTypeTest::VectorField(double* dState, double* StateIn) {
	// dV = -(V-Vrest)/tau
	dState[0] = (I_external-(StateIn[0]-Param[1]))/Param[0];
}

void NeuronTypeTest::PostTimeStep() {
	// Check for threshold passing:
	if(State[0] > Param[2]) {
		State[0] = Param[3];
		MirrorVm = 30;
		Spike();
	} else {
		MirrorVm = State[0];
	}
}

double NeuronTypeTest::GetVm() {
	return MirrorVm;
}
