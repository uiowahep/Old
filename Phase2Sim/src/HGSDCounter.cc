#include "HGSDCounter.hh"
#include "G4UnitsTable.hh"
#include "G4TrackStatus.hh"
#include "G4VProcess.hh"

#include "HGSDCounter.hh"

#include <iostream>

using namespace std;
using namespace CLHEP;

//
//	Constructor
//
HGSDCounter::HGSDCounter(G4String  name, RunParams inData, TTree *tree)
	: G4VSensitiveDetector(name), 
	_runParams(inData)
{
	_tree = tree;
	_subDName = name;
//	_tree->Branch("dr", &_dr);
//	_tree->Branch("resp", &_response, "resp/D");
}

//
//	Constructor
//
HGSDCounter::HGSDCounter(G4String name, RunParams runParams, int id, 
		HGCReadoutModule *readout)
	: G4VSensitiveDetector(name), _readout(readout), _runParams(runParams)
{
	G4cout << "### Setting up HGSDCounter for: " << id << "  " << name << G4endl;
//	_tree = _readout->_subDets[id].dataTree;
//	_subDName = _readout->_subDets[id].name;
	_id = id;

	//
	//	Book histos for P and M sides(id differen of 1);
	//
	_readout->Book(_id);
	_readout->Book(id+1);
	G4cout << "### Done" << G4endl;
}

//
//	Desctructor
//
HGSDCounter::~HGSDCounter()
{}

//
//	Initialize
//
void HGSDCounter::Initialize(G4HCofThisEvent*)
{
//	_dx.clear();
//	_dy.clear();
//	_dz.clear();
//	_dr.clear();

	_response=0;
	//
	//	Clean up for both P/M sides
	//
	_readout->BeginEvent(_id);
	_readout->BeginEvent(_id+1);
}

//
//	Process Hits
//
G4bool HGSDCounter::ProcessHits(G4Step *aStep, G4TouchableHistory*)
{
//	cout << "Generic Hit" << endl;
//	G4cout << aStep->GetTrack()->GetParticleDefinition()->GetParticleName()
//		<< G4endl;
//	cout << 
	
	if (aStep->GetTrack()->GetParticleDefinition()->GetPDGCharge()!=0)
	{
		if (_runParams.verbosity>1)
		{
			cout << "### Hit..." << endl;
		}

		//
		//	Get touchable
		//
		G4TouchableHandle touchable = 
			aStep->GetPreStepPoint()->GetTouchableHandle();

		//
		//	Get whatever is neccessary to produce a HIT
		//
		G4ThreeVector dPos = aStep->GetDeltaPosition();
		G4ThreeVector glPrePos = aStep->GetPreStepPoint()->GetPosition();
		G4ThreeVector glPostPos = aStep->GetPostStepPoint()->GetPosition();
		G4int layerNum = touchable->GetCopyNumber(1);
		G4int detNumber = touchable->GetCopyNumber(2);

		//
		//	Generate a Hit:
		//	M detector parts have ID+1
		//	P detector have just ID
		//
		HGCHit hit;
		hit.dPos = dPos;
		hit.glPrePos = glPrePos;
		hit.layerNumber = layerNum;
		hit.detID = detNumber;

		//
		//	Push the hit to the ReadoutModule
		//
		_readout->PushHit(hit);

//		_dx.push_back(dPos.x()/um);
//		_dy.push_back(dPos.y()/um);
//		_dz.push_back(dPos.z()/um);
//		_dr.push_back(dPos.mag()/um);

		Double_t cresponse =  (dPos.mag()/um)*80;
		
		//	
		//	This is just a total deposit...
		//	Check!!!
		//
		_response += cresponse;

		//
		//	Debug
		//
		if (_runParams.verbosity>1)
		{
			G4cout << "### SD#" << _id << ":" << G4endl;
			G4cout << "### glPrePos=" << glPrePos/mm << "  " << glPostPos/mm
				<< "  " << layerNum << "  " << detNumber << G4endl;
		}


/*		if (_runParams.verbosity>1)
		{
			G4String preName = 
				aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
			G4String postName =
				aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
	
			//	Get local/global position for pre Volume
			//
			G4ThreeVector preGlobalPos = 
				aStep->GetPreStepPoint()->GetPosition();
			G4TouchableHandle theTouchable = 
				aStep->GetPreStepPoint()->GetTouchableHandle();
			G4ThreeVector preLocalPos =
				theTouchable->GetHistory()->GetTopTransform().TransformPoint(
				preGlobalPos);


			G4ThreeVector postGlobalPos = 
				aStep->GetPostStepPoint()->GetPosition();
			theTouchable = 
				aStep->GetPostStepPoint()->GetTouchableHandle();
			G4ThreeVector postLocalPos =
				theTouchable->GetHistory()->GetTopTransform().TransformPoint(
				postGlobalPos);

			G4cout << "### A hit in Part: " << _subDName
				<< G4endl
				<< preName << "  " << postName << " " << dPos.mag()/um
				<< G4endl
				<< preGlobalPos/um << "  " << preLocalPos/um
				<< G4endl
				<< postGlobalPos/um << "  " << postLocalPos/um
				<< G4endl;
		}

	*/
	}

	return true;
}

//
//	End of Event
//
void HGSDCounter::EndOfEvent(G4HCofThisEvent*)
{
	//
	//	Finish for both P/M sides
	//
	_readout->FinishEvent(_id);
	_readout->FinishEvent(_id + 1);
//	_tree->Fill();
}
