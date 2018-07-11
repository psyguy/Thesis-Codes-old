#include "Basics.h"
#include "Network.h"

#include <iostream>
#include <stdlib.h>
#include <time.h>


using namespace std;

int main() {
	
	ifstream fInit("Inits.txt");
	int nStep;//, flag;
	fInit >> nStep;
	fInit >> DynSys::dt;
	int OutVm, OutSp, OutSynState;
	fInit >> OutVm;
	fInit >> OutSp;
	fInit >> OutSynState;
	fInit.close();


	ifstream fCell("CellsIn.txt");
	ifstream fSyn("SynIn.txt");
	ifstream fCon("ConIn.txt");
	Network::ReadCells(&fCell);
	Network::ReadSynapses(&fSyn);
	Network::ReadConnections(&fCon);

	fCell.close();
	fSyn.close();
	fCon.close();

	ofstream DumpVm;
	ofstream DumpSp;
	ofstream DumpSynState;
	if(OutVm) {
		DumpVm.open("Vm.txt");
		Network::SetOutVm(&DumpVm);
	}
	if(OutSp) {
		DumpSp.open("Spikes.txt");
		Network::SetOutSpikes(&DumpSp);
	}
	if(OutSynState) {
		DumpSynState.open("SynState.txt");
		Network::SetOutSynState(&DumpSynState);
	}

	// Initialize random seed
	srand(time(NULL));
    int Nupdate = nStep / 100;
	for(int iStep=0; iStep<nStep; iStep++) {
		DynSys::StepAll();

		if(iStep%2==0) Network::DistributeSpike();

		if(iStep%20==0) Network::PrintAll();
		
		if(iStep%Nupdate==0) cout<<iStep<<" of "<<nStep<<endl;
	}	
	DumpVm.close();
	DumpSp.close();
	DumpSynState.close();

	return 0;
}
