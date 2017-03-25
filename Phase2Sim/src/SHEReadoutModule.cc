#include "HGCReadoutModule.hh"

#include "G4UnitsTable.hh"
#include <math.h>
#include <iostream>

using namespace std;
using namespace CLHEP;

SHEReadoutModule::SHEReadoutModule()
{

}

SHEReadoutModule::SHEReadoutModule(TFile* rootFile, RunParams runParams)
	:ReadoutModule(rootFile, runParams)
{
	cout << "### Configuring Shashlik+HE Configuration... Reaodut Module" << endl;

	//
	//	Declare all the readout trees
	//
}

SHEReadoutModule::~SHEReadoutModule()
{

}

//
//	Do in the beginning of the event
//
void SHEReadoutModule::BeginEvent(int)
{

}

//
//	Do at the end of the event
//
void SHEReadoutModule::FinishEvent(int)
{

}

//
//	Book the vectors/Branches we need
//
void SHEReadoutModule::Book(int)
{

}

//
//	Push the Prim Hit
//
void SHEReadoutModule::PushHit(PRIMHit)
{

}

//
//	Push the HE Hit
//
void SHEReadoutModule::PushHit(HEHit)
{

}

//
//	Push the Shashlik Hit
//
void SHEReadoutModule::PushHit(SHHit)
{

}

//
//	Push the Geometry Info
//
void SHEReadoutModule::PushGeomInfo(CMSSHE)
{

}

//
//	Check if this channel has been triggered.
//	for HE Detector
//	
bool SHEReadoutModule::IsChannelTriggered(int, int, HEHit)
{

	return false;
}

//
//	Check if the channel has been triggered for the Shashlik
//	Detector
//
bool SHEReadoutModule::IsChannelTriggered(int, int, SHHit)
{

	return false;
}

//
//	Compute the channel for SH Detector
//
int SHEReadoutModule::ComputeSHChannel(G4ThreeVector, int)
{
	int channel = 0;

	return channel;
}

int SHEReadoutModule::ComputeHEChannel(G4ThreeVector, int)
{
	int channel = 0;

	return channel;
}































