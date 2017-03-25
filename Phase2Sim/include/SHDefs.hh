#ifndef SHDEFS_H
#define SHDEFS_H

#include "globals.hh"
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TApplication.h"
#include "TTree.h"

#include "G4ThreeVector.hh"

using namespace std;

/*
 *	Run Parameters
 *
 */
typedef struct RunParams
{
	int iEnergy;
	int iPrim;
	char *hgCalInputFileName;
	char *lysoInputFileName;
	char *heInputFileName;
	int numLayers;
	int numModules;
	int isInteractive;
	int verbosity;
	int seed;
	int numEvents;
	TFile *rootFile;
	char *endcapConfigFileName;
	int endcapType;
	char *minBiasFileName;
	char *sigFileName;
	int startEvent;
};

//
//	HE Parameters
//
typedef struct HEParameters
{
	int nLayers_Total;
//	int iAbsMaterial;
	G4double fullHEAbsXYZ[3];
	G4double fullHEActXYZ[3];
};


typedef struct HGCal_EMParameters
{
	int nLayers_Total;
	int nLayers_1, nLayers_2, nLayers_3;
	int iAbsMaterial;
//	G4double fullEMModuleXYZ[3];
//	G4double fullEMLayerXYZ[9];
	G4double fullEMAbsXYZ[5];
//	G4double fullEMPadLayerXYZ[3];
	G4double fullEMPadXYZ[5];
	G4double fullEMPadReadoutXYZ[3];
};

typedef struct HGCal_FHParameters
{
	int nLayers_Total;
	int iAbsMaterial;
//	G4double fullFHModuleXYZ[3];
//	G4double fullFHLayerXYZ[3];
	G4double fullFHAbsXYZ[3];
//	G4double fullFHPadLayerXYZ[3];
	G4double fullFHPadXYZ[3];
	G4double fullFHPadReadoutXYZ[3];
};

typedef struct HGCal_BHParameters
{
	int nLayers_Total;
	int iAbsMaterial;
//	G4double fullBHModuleXYZ[3];
//	G4double fullBHLayerXYZ[3];
	G4double fullBHAbsXYZ[3];
//	G4double fullBHPadLayerXYZ[3];
	G4double fullBHPadXYZ[3];
	G4double fullBHPadReadoutXYZ[3];
};

typedef struct HGCalParameters
{
//	G4double fullHGCalXYZ[3];
	HGCal_EMParameters em;
	HGCal_FHParameters fh;
	HGCal_BHParameters bh;
	int emOnOff;
	int fhOnOff;
	int bhOnOff;
};

typedef struct HGCHit
{
	G4ThreeVector dPos;
	G4ThreeVector glPrePos;
	int layerNumber;
	int detID;
};

typedef struct HEHit
{
	G4ThreeVector pos;
	double numPhotons;
	int layerNumber;
	int detID;
};

typedef struct PRIMHit
{
	G4ThreeVector pos;
	G4ThreeVector mom;
	G4double ene;
	G4String name;
	G4int id;
};

typedef struct SHHit
{
	G4ThreeVector pos;
	double numPhotons;
	int layerNumber;
	int detID;
};

#endif
