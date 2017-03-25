#include "SHSDCounter.hh"
#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"

#include <iostream>

using namespace std;
using namespace CLHEP;

const double SCINTYIELD_LYSO = 32000.;

SHSDCounter::SHSDCounter(G4String name, RunParams inData, TTree *tree)
	: G4VSensitiveDetector(name),
	_runParams(inData)
{
	
}


//	
//	Constructor
//
SHSDCounter::SHSDCounter(G4String name, RunParams runParams,
		int subDetid, HGCReadoutModule *readout)
	: G4VSensitiveDetector(name), _runParams(runParams),
	_readout(readout)
{
	G4cout << "### Setting up SHSD for ID: " << subDetid << G4endl;

	_id = subDetid;
	_subDName = name;
	_rand.SetSeed(_runParams.seed);
	_readout->Book(_id);
	_readout->Book(_id+1);
	G4cout << "### DONE!" << G4endl;

}

//	Destructor
//
SHSDCounter::~SHSDCounter()
{

}

//	Initialize
//
void SHSDCounter::Initialize(G4HCofThisEvent*)
{
	_readout->BeginEvent(_id);
	_readout->BeginEvent(_id+1);
}

//	Process Hits
//
G4bool SHSDCounter::ProcessHits(G4Step *aStep, G4TouchableHistory*)
{
	if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == 
			"opticalphoton")
	{
		//
		//	For Legacy
		//
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);	
	}

	if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() !=
			"opticalphoton")
	{
//		if (_runParams.verbosity > 0)
//			cout << "### Shashlik Hit" <<endl;

		double thisDep = aStep->GetTotalEnergyDeposit()/MeV;
		eneDep += thisDep;
		int gen = ComputeOP(thisDep, SCINTYIELD_LYSO);
		numPhotons += gen;

		//
		//	Get touchable
		//
		G4TouchableHandle touchable = 
			aStep->GetPreStepPoint()->GetTouchableHandle();
		G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
		
		int layerNum = -1;
		int detNumber = -1;
		G4String physName = 
			aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
		if (physName != "pLastEMAct")
		{
			layerNum = touchable->GetCopyNumber(1);
			detNumber = touchable->GetCopyNumber(2);
		}
		else
		{
			layerNum = touchable->GetCopyNumber();
			detNumber = touchable->GetCopyNumber(1);
		}

		//
		//	Generate a Hit
		//
		SHHit hit;
		hit.pos = pos;
		hit.numPhotons = gen;
		hit.layerNumber = layerNum;
		hit.detID = detNumber;
		_readout->PushHit(hit);
	}

	return true;
}

//	Finalize the event
//
void SHSDCounter::EndOfEvent(G4HCofThisEvent*)
{
	_readout->FinishEvent(_id);
	_readout->FinishEvent(_id+1);
}

int SHSDCounter::ComputeOP(double ene, double opPerMeV)
{
	int genPhotons = 0;
	double mean = ene*opPerMeV;
	double sigma = sqrt(mean);
//	genPhotons = floor(_rand.Gaus(mean, sigma));

	genPhotons = mean;
	return genPhotons;
}




















