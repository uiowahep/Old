#ifndef HESD_H
#define HESD_H

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "SHDefs.hh"
#include "HGCReadoutModule.hh"

#include <vector>
#include <string>

#include "TTree.h"
#include "TRandom.h"

using namespace std;

//
//	Define the HE(BH) SD Class
//
class HESD : public G4VSensitiveDetector
{
	public:
		HESD(G4String, RunParams, int subDetid, HGCReadoutModule *readout);
		virtual ~HESD();

		void Initialize(G4HCofThisEvent*);
		G4bool ProcessHits(G4Step*, G4TouchableHistory*);
		void EndOfEvent(G4HCofThisEvent*);
		int ComputeOP(double ene, double opPerMeV);

		int _id;
		G4String _subDName;
		HGCReadoutModule *_readout;
		RunParams _runParams;
		Double_t numPhotons;
		Double_t eneDep;
		TRandom _rand;
};

#endif
