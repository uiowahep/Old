#ifndef HGSDCOUNTER_H
#define HGSDCOUNTER_H

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "SHDefs.hh"
#include "HGCReadoutModule.hh"

#include <vector>
#include <string>

#include "TTree.h"
#include "TH1D.h"

using namespace std;

//
//	Define the HGCal SD Class
//
class HGSDCounter : public G4VSensitiveDetector
{
	public:
		HGSDCounter(G4String, RunParams, TTree*);
		HGSDCounter(G4String, RunParams, int subDetid, HGCReadoutModule *readout);
		virtual ~HGSDCounter();

		void Initialize(G4HCofThisEvent*);
		G4bool ProcessHits(G4Step*, G4TouchableHistory*);
		void EndOfEvent(G4HCofThisEvent*);

		TTree *_tree;
		int _id;
		G4String _subDName;
		RunParams _runParams;
		HGCReadoutModule *_readout;
//		vector<double> _dr_EM;
//		vector<double> _dr_FH;
//		vector<double> _dr_BH;
		vector<double> _dr;
		vector<double> _dx;
		vector<double> _dy;
		Double_t _response;
		vector<double> _dz;
};

#endif 
