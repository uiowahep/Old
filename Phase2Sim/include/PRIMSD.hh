#ifndef PRIMSD_H
#define PRIMSD_H

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "SHDefs.hh"
#include "HGCReadoutModule.hh"

#include <iostream>
#include <vector>

#include "TTree.h"
#include "TH1D.h"

using namespace std;

//
//	SD for primaries
//
class PRIMSD : public G4VSensitiveDetector
{
	public:
		PRIMSD(G4String, RunParams, int , HGCReadoutModule *);
		virtual ~PRIMSD();

		void Initialize(G4HCofThisEvent*);
		G4bool ProcessHits(G4Step*, G4TouchableHistory*);
		void EndOfEvent(G4HCofThisEvent*);

		TTree *tree;
		int _id;
		G4String _subDName;
		RunParams _runParams;
		HGCReadoutModule *_readout;
};

#endif
