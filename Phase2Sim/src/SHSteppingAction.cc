#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UnitsTable.hh"

#include "SHSteppingAction.hh"

#include <iostream>

using namespace std;
using namespace CLHEP;

SHSteppingAction::SHSteppingAction(RunParams inData)
	: runParams(inData)
{;}

SHSteppingAction::~SHSteppingAction()
{;}

void SHSteppingAction::UserSteppingAction(const G4Step* aStep)
{
	//	To avoid INfinite Internal Reflections
	//
	if (aStep->GetTrack()->GetCurrentStepNumber()>1000)
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
/*	G4String pName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
	if (pName == "opticalphoton")
	{
		//	Just to double check
		//
		if (runParams.verbosity>0)
			cout << pName << endl;

		G4double opEne = aStep->GetTrack()->GetTotalEnergy()/eV;
		G4double hPlanck = 4.136*ev;
		G4double cc = 
	}
	*/
}

