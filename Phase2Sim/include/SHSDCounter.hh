#ifndef CCOUNTERSD_H
#define CCOUNTERSD_H 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "globals.hh"

#include <vector>
#include <string>

#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"

#include "SHDefs.hh"
#include "HGCReadoutModule.hh"

using namespace std;

//	Define the SD Class
//
class SHSDCounter : public G4VSensitiveDetector
{
	public:
		SHSDCounter(G4String, RunParams, TTree*);
		SHSDCounter(G4String, RunParams, int subDetid, HGCReadoutModule *readout);
		virtual ~SHSDCounter();

		void Initialize(G4HCofThisEvent*);
		G4bool ProcessHits(G4Step*, G4TouchableHistory*);
		void EndOfEvent(G4HCofThisEvent*);
		int ComputeOP(double ene, double opPerMeV);

		//	Members
		//
		Double_t numPhotons;
		Double_t eneDep;
		TRandom _rand;
		HGCReadoutModule *_readout;
		int _id;
		G4String _subDName;

		RunParams _runParams;
};

#endif
