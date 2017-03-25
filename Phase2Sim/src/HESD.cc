




#include "HESD.hh"
#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"

using namespace CLHEP;

const Double_t SCINTYIELD_SCSN81 = 10000;

//
//	COnstructor
//
HESD::HESD(G4String name, RunParams runParams, int subDetid, 
		HGCReadoutModule *readout)
	: G4VSensitiveDetector(name),
	_runParams(runParams), _readout(readout)
{
	_id = subDetid;
	G4cout << "### Setting up the HESD for : " << _id << "  " << name << G4endl;
	_subDName = name;
	_rand.SetSeed(_runParams.seed);
	_readout->Book(_id);

	G4cout << "### Setting up the HESD(Minus Side) for : " << _id+1 
		<< "  " << name << G4endl;
	_readout->Book(_id+1);
	G4cout << "### DONE!" << G4endl;
}

//
//	Destructor
//
HESD::~HESD()
{

}

//
//	Initialize 
//
void HESD::Initialize(G4HCofThisEvent*)
{
	_readout->BeginEvent(_id);
	_readout->BeginEvent(_id+1);
}

//
//	Process Hits
//
G4bool HESD::ProcessHits(G4Step *aStep, G4TouchableHistory*)
{

	if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() == 
			"opticalphoton")
	{
		//
		//	This mustn't happen!!!
		//
		aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	}

	if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName() !=
			"opticalphoton")
	{
		double thisDep = aStep->GetTotalEnergyDeposit()/MeV;
		eneDep += thisDep;
		int gen = ComputeOP(thisDep, SCINTYIELD_SCSN81);
		numPhotons += gen;

		//
		//	Get touchable 
		//
		G4TouchableHandle touchable = 
			aStep->GetPreStepPoint()->GetTouchableHandle();
		G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
		int layerNum = touchable->GetCopyNumber(1);
		int detNumber = touchable->GetCopyNumber(2);

//		cout << layerNum << "  " << detNumber << endl;
//		G4cout << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;

		//
		//	Generate a Hit
		//
		HEHit hit;
		hit.pos = pos;
		hit.numPhotons = gen;
		hit.layerNumber = layerNum;
		hit.detID = detNumber;

		if (_runParams.verbosity>0)
		{
			cout << "### " << detNumber << "  " << layerNum << endl <<
			aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << 
			_id << endl;
		}

		_readout->PushHit(hit);
	}

	return true;
}

//
//	Simulate the #photons for Scintillation
//
int HESD::ComputeOP(double ene, double opPerMeV)
{
	int genPhotons = 0;
	double mean = ene*opPerMeV;
	double sigma = sqrt(mean);
	genPhotons = floor(_rand.Gaus(mean, sigma));

	return genPhotons;
}

//
//	End of Event
//
void HESD::EndOfEvent(G4HCofThisEvent*)
{
	_readout->FinishEvent(_id);
	_readout->FinishEvent(_id+1);
}

