#include "SHPrimaryGeneratorAction.hh"

#include <iostream>

#include "G4HEPEvtInterface.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "CLHEP/Random/RandFlat.h"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

using namespace CLHEP;
using namespace std;

//	0-7 are old energies(for EM)
//	8-12 are new energies(for HAD)
//
const int numEnergies = 13;
const int numPrims = 3;
G4double primEnergies[numEnergies] = {
	1.*GeV, 2.*GeV, 4.*GeV, 8.*GeV, 16.*GeV, 32.*GeV, 50.*GeV, 60.*GeV,
	20*GeV, 50*GeV, 100*GeV, 200*GeV, 300*GeV
};
G4String primNames[numPrims] = {"e-", "pi-", "mu-"};



//	Constructor
//
SHPrimaryGeneratorAction::SHPrimaryGeneratorAction(RunParams params, 
		HGCReadoutModule *readout)
	: runParams(params)
{
	_readout = readout;
	_rand.SetSeed(runParams.seed);

	//	Define a generator
	//	Random dir(up to dims of Endcap)
	//	Random Energy(up to 200GeV)
	//
	primName = primNames[runParams.iPrim];
	primEnergy = _rand.Uniform(1, 500)*GeV;

	particleGun = new G4ParticleGun(1);
//	particleGun->SetParticleEnergy(primEnergy);
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *particle = particleTable->FindParticle(primName);
//	particleGun->SetParticleDefinition(particle);

	primPos = G4ThreeVector(0*mm, 0*mm, 0.*mm);
//	primDir = G4ThreeVector(1.*m, 0., 0);
	
	G4double dirx = _rand.Uniform(-1600, 1600)*mm;
	G4double diry = _rand.Uniform(-1600, 1600)*mm;
	G4double dirz = 3150*mm;
	primDir = G4ThreeVector(dirx, diry, dirz);
//	primDir = G4ThreeVector(0, 0, -3150*mm);
//	particleGun->SetParticlePosition(primPos);
//	particleGun->SetParticleMomentumDirection(primDir);
//	particleGun->SetParticleEnergy(primEnergy);

	//
	//	Initialize the SIG  and MinBias Files
	//
	_sigFile = new TFile(runParams.sigFileName);
	_sigFile->cd("HtoGGAnalyzer");
	_sigTree = (TTree*)gDirectory->Get("GEN");
	_sigTree->SetBranchAddress("status_SIG", &genSig.status_SIG);
	_sigTree->SetBranchAddress("pdgID_SIG", &genSig.pdgId_SIG);
	_sigTree->SetBranchAddress("pE_SIG", &genSig.pE_SIG);
	_sigTree->SetBranchAddress("pX_SIG", &genSig.pX_SIG);
	_sigTree->SetBranchAddress("pY_SIG", &genSig.pY_SIG);
	_sigTree->SetBranchAddress("pZ_SIG", &genSig.pZ_SIG);
	_sigTree->SetBranchAddress("x_SIG", &genSig.x_SIG);
	_sigTree->SetBranchAddress("y_SIG", &genSig.y_SIG);
	_sigTree->SetBranchAddress("z_SIG", &genSig.z_SIG);

	_sigTree->SetBranchAddress("status", &genSig.status);
	_sigTree->SetBranchAddress("pdgID", &genSig.pdgId);
	_sigTree->SetBranchAddress("pE", &genSig.pE);
	_sigTree->SetBranchAddress("pX", &genSig.pX);
	_sigTree->SetBranchAddress("pY", &genSig.pY);
	_sigTree->SetBranchAddress("pZ", &genSig.pZ);
	_sigTree->SetBranchAddress("x", &genSig.x);
	_sigTree->SetBranchAddress("y", &genSig.y);
	_sigTree->SetBranchAddress("z", &genSig.z);

	_minBiasFile = new TFile(runParams.minBiasFileName);
	_minBiasFile->cd("HtoGGAnalyzer");
	_minBiasTree = (TTree*)gDirectory->Get("GEN");

	_minBiasTree->SetBranchAddress("status", &genMinBias.status);
	_minBiasTree->SetBranchAddress("pdgID", &genMinBias.pdgId);
	_minBiasTree->SetBranchAddress("pE", &genMinBias.pE);
	_minBiasTree->SetBranchAddress("pX", &genMinBias.pX);
	_minBiasTree->SetBranchAddress("pY", &genMinBias.pY);
	_minBiasTree->SetBranchAddress("pZ", &genMinBias.pZ);
	_minBiasTree->SetBranchAddress("x", &genMinBias.x);
	_minBiasTree->SetBranchAddress("y", &genMinBias.y);
	_minBiasTree->SetBranchAddress("z", &genMinBias.z);

}

//	Destructor
//
SHPrimaryGeneratorAction::~SHPrimaryGeneratorAction()
{
	delete particleGun;
}

//	Generate Primary Event
//
void SHPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
	G4double totE = 0;
	G4double totE_PU = 0;

	//
	//	Genearte Primaries from the Sig Tree and from MinBias Tree
	//
	int globalEventToGenerate = anEvent->GetEventID()	+ runParams.startEvent;
	_sigTree->GetEntry(globalEventToGenerate);
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();

	cout << "### Global Start Event Number: " << runParams.startEvent << endl;
	cout << "### Global Event Number: " << globalEventToGenerate << endl;


	//
	//	first Generate the Signal/Signature Verteces
	//
	for (int iPrim=0; iPrim<genSig.pdgId->size(); iPrim++)
	{
		G4ParticleDefinition *particle = 
			particleTable->FindParticle(genSig.pdgId->at(iPrim));
		particleGun->SetParticleDefinition(particle);
		G4double x = genSig.x->at(iPrim)*mm;
		G4double y = genSig.y->at(iPrim)*mm;
		G4double z = genSig.z->at(iPrim)*mm;
		G4double px = genSig.pX->at(iPrim)*GeV;
		G4double py = genSig.pY->at(iPrim)*GeV;
		G4double pz = genSig.pZ->at(iPrim)*GeV;
		G4double ene = genSig.pE->at(iPrim)*GeV;
		particleGun->SetParticlePosition(G4ThreeVector(x,y,z));
		particleGun->SetParticleMomentumDirection(
				G4ThreeVector(px, py, pz).unit());
		particleGun->SetParticleEnergy(ene);
		particleGun->GeneratePrimaryVertex(anEvent);
		totE +=  ene;
		totE_PU += ene;
	}

	int startPUEvent = globalEventToGenerate*140;
	for (int iPUEvent=0; iPUEvent<140; iPUEvent++)
	{
		_minBiasTree->GetEntry(startPUEvent + iPUEvent);
		cout << "### Mixing PU Event " << iPUEvent << "  Global PU Event: "
			<< startPUEvent+iPUEvent << endl;
		for (int iPrim=0; iPrim<genMinBias.pdgId->size(); iPrim++)
		{
			G4ParticleDefinition *particle = 
				particleTable->FindParticle(genMinBias.pdgId->at(iPrim));
			particleGun->SetParticleDefinition(particle);
			G4double x = genMinBias.x->at(iPrim)*mm;
			G4double y = genMinBias.y->at(iPrim)*mm;
			G4double z = genMinBias.z->at(iPrim)*mm;
			G4double px = genMinBias.pX->at(iPrim)*GeV;
			G4double py = genMinBias.pY->at(iPrim)*GeV;
			G4double pz = genMinBias.pZ->at(iPrim)*GeV;
			G4double ene = genMinBias.pE->at(iPrim)*GeV;
			particleGun->SetParticlePosition(G4ThreeVector(x, y, z));
			particleGun->SetParticleMomentumDirection(
					G4ThreeVector(px, py, pz).unit());
			particleGun->SetParticleEnergy(ene);
			particleGun->GeneratePrimaryVertex(anEvent);
			totE_PU += ene;
		}
	}
	
	G4cout << "### Total Energy Per Event: " << totE/GeV << " without PU" << endl
		<< "### Total Energy Per Event: " << totE_PU/GeV << " with PU." << endl;

	/*
	//	All gun's settings have been set up in Constructor
	//
	//particleGun->SetParticleEnergy(10.*keV);
//	G4cout << particleGun->GetParticleEnergy()/GeV << "  " 
//		<< particleGun->GetParticleDefinition()->GetParticleName()
//		<< G4endl; 

//	primEnergy = _rand.Uniform(1, 500)*GeV;
	particleGun->SetParticleEnergy(primEnergy);

	G4double dirx = _rand.Uniform(-1600., 1600.)*mm;
	G4double diry = _rand.Uniform(-1600., 1600.)*mm;
	G4double dirz = 3150*mm;
	primDir = G4ThreeVector(dirx, diry, dirz);
//	primDir = G4ThreeVector(0, 0, -3150*mm);
	particleGun->SetParticleMomentumDirection(primDir);
	particleGun->GeneratePrimaryVertex(anEvent);
	*/
}
