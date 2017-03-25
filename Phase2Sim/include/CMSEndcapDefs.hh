#ifndef CMSENDCAPDEFS_H
#define CMSENDCAPDEFS_H

//	G4 inclusions
//
#include "globals.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Material.hh"

//	ROOT inclusions
//
#include "TROOT.h"
#include "TFile.h"

//	STD inclusions
//
#include <vector>

//
//	HGCal Detector Description
//

//	Describe the constituents first
//
struct HGC_EM
{
	//	EM module itself
	//	-> 30 Layers(3 kinds of 10 each, different thickness - have to have 
	//	3 solids)
	//		-> Abs (3 kinds as stated above)
	//		-> Si Pads (2 kinds)
	//		-> Electronics Readout
	//
	//	
	//	Important Dimensions, # layers, etc...
	//	Radii are given for the whole EM
	//	Zs are for abs, pad, readout from respective parts;
	//	The rest can be calculated
	//	Part1: 10 layers 0.5X0
	//	Part2: 10 layers 0.8X0
	//	Part3: 10 layers 1.2X0
	//	NOTE: rmin2, rmax2 aren't used - Number of layers + dimensions of 
	//	constituents are enough!!!
	//
	int n1, n2, n3;
	int iMat;
	int onOff;
	G4Material *absMat;
	G4double centerZ;
	G4double totalZ;
	G4double rmin1, rmin2;
	G4double rmax1, rmax2;
	G4double fullAbsZ[3];
	G4double fullPadZ;
	G4double fullReadoutZ;
};

//
//	HGCal::Front H-part
//	rmin2, rmax2 aren't used currently....
//
struct HGC_FH
{
	int n;
	int iMat;
	G4Material *absMat;
	G4double centerZ;
	G4double totalZ;
	G4double rmin1, rmax1;
	G4double rmin2, rmax2;
	G4double fullAbsZ, fullPadZ, fullReadoutZ;
};

//
//	HGCal::Back H-part(HE)
//
//
struct HGC_BH
{
	int n;
	int iMat;
	G4Material *absMat;
	G4double centerZ;
	G4double totalZ;
	G4double rmin1, rmax1;
	G4double rmin2, rmax2;
	G4double fullAbsZ;
	G4double fullActZ;
};

//
//	HGCal::Side Hpart(Side HE)
//
struct HGC_SHE
{
	int n;
	int iMat;
	G4Material *absMat;
	G4double centerZ;
	G4double totalZ;
	double incline;
	G4double rmin1, rmax1;
	G4double rmin2, rmax2;
	G4double fullAbsZ, fullActZ;
};

//	HGC itself
//
struct CMSHGCal
{
	G4double startZ;
	HGC_EM EM;
	HGC_FH FH;
	HGC_BH BH;
	HGC_SHE SHE;
};

//
//	EM(Shashlik)
//
struct SHE_EM
{
	int n;
	int iMat;
	int onOff;
	G4Material *absMat;
	G4double centerZ;
	G4double totalZ;
	G4double rmin1, rmax1;
	G4double rmin2, rmax2;
	G4double fullAbsZ;
	G4double fullActZ;
};

//	
//	FHE - HE part just after Shashlik and goes from 1.47 up to 3 in eta
//	and in z to be determined.
//
struct SHE_HE
{
	int n;
	int iMat;
	G4Material *absMat;
	G4double centerZ;
	G4double totalZ;
	G4double rmin1, rmax1;
	G4double rmin2, rmax2;
	G4double fullAbsZ, fullActZ;
};

//

//
//	CMSSPHE:
//	-- EM is Shashlik 
//	-- HE is the upgraded HE
//
struct CMSSHE
{
	G4double startZ;
	SHE_EM EM;
	SHE_HE FHE;
	SHE_HE MHE1;
	SHE_HE MHE2;
	SHE_HE BHE;
};

//
//	CMS Encap Scenarious as Structs
//
struct CMSEndcap
{
	CMSHGCal cmsHGCal;	
	CMSSHE cmsSHE;
};

#endif
