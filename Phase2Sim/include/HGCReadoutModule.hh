#ifndef HGCREADOUTMODULE_H
#define HGCREADOUTMODULE_H

//
//	ROOT Inclusions
//
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"

//
//	STD inclusions
//
#include <iostream>
#include <vector>
#include <string>

//
//	Phase2Sim inclusions
//
#include "CMSEndcapDefs.hh"
#include "SHDefs.hh"

using namespace std;

const int NUMETAS = 20;

//
//	Declare a struct which will keep info on subDets
//
struct SubDetector
{
	string name;
	int id;
	TTree *dataTree;
	vector<int> *channelIDs;
	vector<Double_t> *responses;
	vector<Double_t> *x;
	vector<Double_t> *y;
	vector<Double_t> *z;
};

//
//	Declare a Struct to record PRIM info
//
struct PrimDetector
{
	string name;
	int id;
	TTree *dataTree;
	vector<Double_t> *x;
	vector<Double_t> *y;
	vector<Double_t> *z;
	vector<Double_t> *px;
	vector<Double_t> *py;
	vector<Double_t> *pz;
	vector<Double_t> *ene;
	vector<string> *pname;
	vector<int> *pid;
};

struct HGMaps
{
	int nEM, n1EM, n2EM, n3EM;
	int nFH;
	int nHE;
	Double_t *pL2R_EM;
	Double_t *pL2Z_EM;
	int *pL2NChs1D_EM;
	Double_t *pL2R_FH;
	Double_t *pL2Z_FH;
	int*pL2NChs1D_FH;
	int numEtas;
	int numPhis;
	Double_t *pEtaMap;
	Double_t *pL2Z_HE;
	TTree *tree;
};

struct SHMaps
{
	int nHE;
	int nSH;
	int numEtas;
	int numPhis;
	Double_t *pL2Z_SH;
	Double_t *pL2R_SH;
	Int_t *pL2NChs1D_SH;

	Double_t *pL2Z_HE;
	Double_t *pEtaMap;
	TTree *tree;
};

class ReadoutModule
{
	public:
		ReadoutModule();
		ReadoutModule(TFile*, RunParams);
		virtual ~ReadoutModule();

		virtual void BeginEvent(int);
		virtual void FinishEvent(int);
		virtual void Book(int);

		TFile *_file;
		vector<SubDetector> _subDets;
		PrimDetector _primDet;
		RunParams _runParams;
};

//
//	SHE Readout Module
//
class SHEReadoutModule : public ReadoutModule
{
	public:
		SHEReadoutModule();
		SHEReadoutModule(TFile*, RunParams);
		virtual ~SHEReadoutModule();

		virtual void BeginEvent(int);
		virtual void FinishEvent(int);
		virtual void Book(int);

		void PushHit(PRIMHit);
		void PushHit(HEHit);
		void PushHit(SHHit);
		void PushGeomInfo(CMSSHE);
		bool IsChannelTriggered(int, int, HEHit);
		bool IsChannelTriggered(int, int, SHHit);

		int ComputeSHChannel(G4ThreeVector, int);
		int ComputeHEChannel(G4ThreeVector, int);

		CMSSHE _SHE;
		SHMaps _SHMaps;
};

//
//	Declare out HGCReadoutModule Class:
//	--	Responsible for ReadOut of all subparts of HGCal Endcap
//		-- EM + FH + HE
//	-- Responsible for Mapping of spacial hit info into the actual channels.
//
class HGCReadoutModule : public ReadoutModule
{
	public:
		HGCReadoutModule();
		HGCReadoutModule(TFile*, RunParams);
		virtual ~HGCReadoutModule();

		virtual void BeginEvent(int);
		virtual void FinishEvent(int);
		virtual void Book(int);

		void PushHit(HGCHit);
		void PushHit(PRIMHit);
		void PushHit(HEHit);
		void PushHit(SHHit);
		void PushGeomInfo(CMSHGCal);
		void PushGeomInfo(CMSSHE);
		bool IsChannelTriggered(int, int, HGCHit);
		bool IsChannelTriggered(int, int, HEHit);
		bool IsChannelTriggered(int, int, SHHit);

		int ComputeEMChannel(G4ThreeVector, int);
		int ComputeFHChannel(G4ThreeVector, int);
		int ComputeHEChannel(G4ThreeVector, int);
		int ComputeSHChannel(G4ThreeVector, int);
		int ComputeSHHEChannel(G4ThreeVector, int);

//		TFile *_file;
//		vector<SubDetector> _subDets;
//		PrimDetector _primDet;
//		RunParams _runParams;
		CMSHGCal _HGCal;
		CMSSHE _SHE;
		HGMaps _GMaps;
		SHMaps _SHMaps;
};

#endif




















