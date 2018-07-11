#include "Synapses.h"
#include "Network.h"

#include <iostream>

using namespace std;

#define DEBUG 0
#define EXP1 2.71828//182845905


Synapse::Synapse(int np, int ns): DynSys(np,ns) {
	Network::AddSynapse(this);
};

void Synapse::AddSpike(double time, double weight) {
	#if DEBUG > 0
	cout << "Event added at synapse " << this << endl;
	#endif
	Queue.insert(pair<double,double>(time,weight));
	#if DEBUG > 2
	cout << "  Queue contains: " << endl;
	multimap<double,double>::iterator it;
	for(it=Queue.begin(); it!=Queue.end(); it++) {
		cout << "    t = " << (*it).first << endl;
	}
	#endif
}

double SynapseTypeDelta::GetPulse(double Vm) {
	double Weight = 0;

	#if DEBUG > 0
	cout << "        Processing delta synapse " << this << endl;
	#endif
	
	// Processes all spikes that activate the synapse this time step
	double t = DynSys::Time + DynSys::dt;
	if(!Queue.empty()) {
		#if DEBUG > 1
		cout << "          Queue non-empty" << endl;
		#endif
		double EventTime = Queue.begin()->first;
		while(!Queue.empty() && EventTime < t) {
			#if DEBUG > 2
			cout << "            Processing queue item" << endl;
			#endif
			Weight += Queue.begin()->second;
			Queue.erase(Queue.begin());
			//double EventTime = Queue.begin()->first;  %BJZ wordt niet gebruikt
		}
	}
	return Weight;
};



SynapseTypeExponential::SynapseTypeExponential(): Synapse(3,1) {
/*
	ParamId | Description
	--------+-----------
	   0    | Time constant
       1    | E (reversal potential)
	   2    | g (maximal conductance)

	StateId | Description
	--------+-----------
	   0    | s (synaptic gating)
*/	

}

void SynapseTypeExponential::PreTimeStep() {
	// Processes all spikes that activate the synapse this time step
	double t = DynSys::Time + DynSys::dt;
	if(!Queue.empty()) {
		double EventTime = Queue.begin()->first;
		while(!Queue.empty() && EventTime < t) {
			State[0] += Queue.begin()->second;
			Queue.erase(Queue.begin());
		}
	}
}

void SynapseTypeExponential::VectorField(double* dState, double* StateIn) {
	dState[0] = -StateIn[0]/Param[0];
}


double SynapseTypeExponential::GetCurrent(double Vm) {
	// I_ext = s*g*(Vm-E)
	//cout << State[0]*Param[2]*(Param[1]-Vm) << endl;
	return State[0]*Param[2]*(Param[1]-Vm);
}
double SynapseTypeExponential::GetSynState() {
       return State[0];
}

SynapseTypeAlphaLiley::SynapseTypeAlphaLiley(): Synapse(3,2) {
/*
	ParamId | Description
	--------+-----------
	   0    | gamma (Rate constant, 1/s)
       1    | E    (reversal potential)
	   2    | G    (maximal conductance)

	StateId | Description
	--------+-----------
	   0    | s (synaptic gating)
	   1    | first derivative
*/	

}

void SynapseTypeAlphaLiley::PreTimeStep() {
	// Processes all spikes that activate the synapse this time step
	double t = DynSys::Time + DynSys::dt;
	if(!Queue.empty()) {
		double EventTime = Queue.begin()->first;
		while(!Queue.empty() && EventTime < t) {
			State[1] += EXP1*Param[0]*(Queue.begin()->second);
			Queue.erase(Queue.begin());
		}
	}
}

void SynapseTypeAlphaLiley::VectorField(double* dState, double* StateIn) {
     // dI2/d2t = -2 gamma * dI/dt - gamma^2 * I (+ e*gamma*G*input)
     // input are delta pulses, so they are added, instead of integrated over. G is multiplied with at GetCurrent
	dState[0] = State[1];
    dState[1] = -2*Param[0]*State[1] - Param[0]*Param[0]*State[0]; //+EXP1*gamma*G*inp
}


double SynapseTypeAlphaLiley::GetCurrent(double Vm) {
	// I_ext = s*g*(Vm-E)
	//cout << State[0]*Param[2]*(Param[1]-Vm) << endl;
	return State[0]*Param[2]*(Param[1]-Vm);
	
	
	//dh(3) = h(7); 
    //dh(7)  = -2*g*h(7)  - g^2*h(3) + E*G*g*Spikes/s;
	
	
}
double SynapseTypeAlphaLiley::GetSynState() {
	return State[0];	
}


