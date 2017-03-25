#include "PRIMSD.hh"

#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"

using namespace CLHEP;

//
//	Constructor
//
PRIMSD::PRIMSD(G4String name, RunParams runParams, int id, 
		HGCReadoutModule *readout)
	: G4VSensitiveDetector(name),
	_runParams(runParams), _readout(readout)
{
//	_readout = (HGCReadoutModule*)readout;
	G4cout << "### Setting up PRIMSD for id: " << id << "  " << name << G4endl;
	_id = id;
	_readout->Book(id);
}

//
//	Destructor
//
PRIMSD::~PRIMSD()
{

}

//
//	Pre Event TODO
//
void PRIMSD::Initialize(G4HCofThisEvent*)
{
	_readout->BeginEvent(_id);
}

//
//	End of Event TODO
//
void PRIMSD::EndOfEvent(G4HCofThisEvent*)
{
	_readout->FinishEvent(_id);
}

//
//	Process Hit
//	
G4bool PRIMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	//
	//	Look at Primary Tracks only
	//
	if (aStep->GetTrack()->GetParentID() == 0)
	{
		PRIMHit aHit;
		aHit.pos = aStep->GetPreStepPoint()->GetPosition();
		aHit.mom = aStep->GetPreStepPoint()->GetMomentum();
		aHit.ene = aStep->GetPreStepPoint()->GetKineticEnergy();
		aHit.name = 
			aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
		aHit.id = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
		_readout->PushHit(aHit);
	}

	return true;
}

















