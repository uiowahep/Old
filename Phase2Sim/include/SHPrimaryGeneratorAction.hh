


#ifndef SHPRIMARYGENERATORACTION_H
#define SHPRIMARYGENERATORACTION_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4VPrimaryGenerator.hh"

#include "TRandom.h"

#include "SHDefs.hh"
#include "HGCReadoutModule.hh"

struct genInfo
{
	vector<int> *status;
	vector<int> *pdgId;
	vector<double> *charge;
	vector<double> *pE;
	vector<double> *pX;
	vector<double> *pY;
	vector<double> *pZ;
	vector<double> *x;
	vector<double> *y;
	vector<double> *z;

	vector<int> *status_SIG;
	vector<int> *pdgId_SIG;
	vector<double> *charge_SIG;
	vector<double> *pE_SIG;
	vector<double> *pX_SIG;
	vector<double> *pY_SIG;
	vector<double> *pZ_SIG;
	vector<double> *x_SIG;
	vector<double> *y_SIG;
	vector<double> *z_SIG;
};

class SHPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
		//	Constructor and Destructor
		SHPrimaryGeneratorAction(RunParams params, TTree *tree);
		SHPrimaryGeneratorAction(RunParams, HGCReadoutModule*);
		virtual ~SHPrimaryGeneratorAction();

		//	Standard func	
		void GeneratePrimaries(G4Event *anEvent);

		//	To keep things simplier
		G4ParticleGun *particleGun;
		G4String primName;
		G4double primEnergy;
		G4ThreeVector primPos;
		G4ThreeVector primDir;
		G4int verbosityLevel;

		RunParams runParams;
		TTree *shTree;
		HGCReadoutModule *_readout;
		TRandom _rand;
		TFile *_sigFile;
		TTree *_sigTree;
		genInfo genSig;
		TFile *_minBiasFile;
		TTree *_minBiasTree;
		genInfo genMinBias;
};

#endif
