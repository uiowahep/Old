#ifndef SHDETECTORCONSTRUCTION_H
#define SHDETECTORCONSTRUCTION_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include <vector>

#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "SHDefs.hh"
#include "CMSEndcapDefs.hh"
#include "SHSDCounter.hh"
#include "HGCReadoutModule.hh"

using namespace CLHEP;

class SHDetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		//	Constructor and Destructor
		//
		SHDetectorConstruction(RunParams, const vector<TTree*> &);
		SHDetectorConstruction(RunParams, TTree *tree);
		SHDetectorConstruction(RunParams, TTree *emTree, TTree *fhTree, 
				TTree *bhTree);
		SHDetectorConstruction(RunParams, HGCReadoutModule*);
		~SHDetectorConstruction(); 

		void MyConstructor(RunParams, TTree*);
		G4VPhysicalVolume* Construct();
		G4VPhysicalVolume* BuildGeometry();

		//
		//	BuildEndCapScenario
		//
		int BuildEndcap();
		int BuildHGC();
		int BuildShashlikPlusHE();

		G4double Getr2(G4double r1, G4double z1, G4double z2);
		int ConfigHGC();
		int ConfigShashlikPlusHE();

		//	
		//	Individual Simulations Build Options
		//
		void BuildShashlik();
		void BuildHGCal();
		void BuildHE();
		int ReadHGConfigFile();
		int ReadHEConfigFile();
		void UpdateGeometry();

		RunParams runParams;
		HGCalParameters hgcParameters;
		HEParameters heParameters;
		TTree *_emTree;
		TTree *_fhTree;
		TTree *_bhTree;
		TTree *shTree;
		TTree *_heTree;
		HGCReadoutModule *_readout;

		//	List of Materials/Elements
		//	mAl2O3 Index of 0 and mGaInP index of 8
		G4Element *eLu;
		G4Element *eY;
		G4Element *eSi;
		G4Element *eO;
		G4Element *eH;
		
		//	Materials
		//
		G4Material *mAir;
		G4Material *mFe;
		G4Material *mPb;
		G4Material *mCu;
		G4Material *mSi;
		G4Material *mZn;
		G4Material *mBrass;
		G4Material *mVacuum;
		G4Material *mAl2O3;
		G4Material *mBeO;
		G4Material *mMgO;
		G4Material *mTiO;
		G4Material *mCs3Sb;
		G4Material *mGaP;
		G4Material *mGaAsP;
		G4Material *mGaPCs;
		G4Material *mGaInP;
		G4Material *mW;
		G4Material *mLYSO;
		G4Material *mSiO2;
		G4Material *mG10;
		G4Material *mClm;
		G4Material *mGlass;
		G4Material *mEpoxy;
		
		G4Material *mPOLYETHYLENE;
		G4Material *mEMAbsMat;
		G4Material *mFHAbsMat;
		G4Material *mBHAbsMat;
		G4Material *mHEAbsMat;
		G4Material *mHEScintMain;
		G4Material *mHEScintFront;
		G4Material *mElectronicsMat;

		//
		//	Different Endcap Scenarious
		//
		CMSEndcap cmsEndcap;

		//
		//	Volumes
		//

		//	World
		//	/
		G4Box *solidWorld;
		G4LogicalVolume *logicWorld;
		G4VPhysicalVolume *physWorld;

		//
		//	Shashlik's Part
		//	
		
		//	Shashlik's Container
		//
		G4Box *solidShashlikContainer;
		G4LogicalVolume *logicShashlikContainer;
		G4VPhysicalVolume *physShashlikContainer;

		//	Shashlik Full Length Module
		//
		G4Box *solidShashlik;
		G4LogicalVolume *logicShashlik;
		G4VPhysicalVolume *physShashlik;

		//	1 Layer: Abs + Act Part
		//
		G4Box *solidLayer;
		G4LogicalVolume *logicLayer;
		G4VPhysicalVolume *physLayer;

		//	Absorber
		//
		G4Box *solidAbs;
		G4LogicalVolume *logicAbs;
		G4VPhysicalVolume *physAbs;

		//	Act
		//
		G4Box *solidAct;
		G4LogicalVolume *logicAct;
		G4VPhysicalVolume *physAct;

		//	Fiber: 2 kinds
		//	1) Abs
		//	2) Act
		//
		G4Tubs *solidAbsFiber;
		G4LogicalVolume *logicAbsFiber;
		G4VPhysicalVolume *physAbsFiber;

		G4Tubs *solidActFiber;
		G4LogicalVolume *logicActFiber;
		G4VPhysicalVolume *physActFiber;

		/*
		 *	HGCAL Part
		 *
		 */
		G4Box *solidHGCal;
		G4LogicalVolume *logicHGCal;
		G4VPhysicalVolume *physHGCal;

		//
		//	EM Part Whole Module
		//
		G4Box *solidHG_EMModule;
		G4LogicalVolume *logicHG_EMModule;
		G4VPhysicalVolume *physHG_EMModule;

		//	EM Part: Layer
		//
		G4Box *solidHG_EMLayer[3];
		G4LogicalVolume *logicHG_EMLayer[3];
		G4VPhysicalVolume *physHG_EMLayer;

		//	EM Part: Abs: 3types depending on thickness;
		//
		G4Box *solidHG_EMAbs[3];
		G4LogicalVolume *logicHG_EMAbs[3];
		G4VPhysicalVolume *physHG_EMAbs;

		//	EM Part: Si Pad Layer
		//
		G4Box *solidHG_EMPadLayer;
		G4LogicalVolume *logicHG_EMPadLayer;
		G4VPhysicalVolume *physHG_EMPadLayer;

		//	EM Part: Si Pad: 2 types as in Slide 16
		//	-- 0.9cm2 for first 20layers, and 1.8cm2 for the last 10.
		//	-- the exact dimensions can be varied, layer lay out too....
		//
		G4Box *solidHG_EMPad[2];
		G4LogicalVolume *logicHG_EMPad[2];
		G4VPhysicalVolume *physHG_EMPad;

		G4Box *solidHG_EMPadReadout;
		G4LogicalVolume *logicHG_EMPadReadout;
		G4VPhysicalVolume *physHG_EMPadReadout;

		//
		//	Front Hadronic CAL
		//
		G4Box *solidHG_FHModule;
		G4LogicalVolume *logicHG_FHModule;
		G4VPhysicalVolume *physHG_FHModule;

		//	Fron HCAL: Layer
		//
		G4Box *solidHG_FHLayer;
		G4LogicalVolume *logicHG_FHLayer;
		G4VPhysicalVolume *physHG_FHLayer;

		//	Font HCAL: Abs
		//
		G4Box *solidHG_FHAbs;
		G4LogicalVolume *logicHG_FHAbs;
		G4VPhysicalVolume *physHG_FHAbs;

		//	Front HCAL: Si Pad Layer
		//
		G4Box *solidHG_FHPadLayer;
		G4LogicalVolume *logicHG_FHPadLayer;
		G4VPhysicalVolume *physHG_FHPadLayer;

		//	Front HCAL: Si Pad
		//
		G4Box *solidHG_FHPad;
		G4LogicalVolume *logicHG_FHPad;
		G4VPhysicalVolume *physHG_FHPad;

		//	Front HCAl Readout
		//
		G4Box *solidHG_FHPadReadout;
		G4LogicalVolume *logicHG_FHPadReadout;
		G4VPhysicalVolume *physHG_FHPadReadout;

		//
		//	Back HCAL
		//
		G4Box *solidHG_BHModule;
		G4LogicalVolume *logicHG_BHModule;
		G4VPhysicalVolume *physHG_BHModule;

		//	Back HCAL: Layer
		//
		G4Box *solidHG_BHLayer;
		G4LogicalVolume *logicHG_BHLayer;
		G4VPhysicalVolume *physHG_BHLayer;

		//	Back HCAL: Abs
		//
		G4Box *solidHG_BHAbs;
		G4LogicalVolume *logicHG_BHAbs;
		G4VPhysicalVolume *physHG_BHAbs;

		//	Back HCAL: Si Pad Layer
		//
		G4Box *solidHG_BHPadLayer;
		G4LogicalVolume *logicHG_BHPadLayer;
		G4VPhysicalVolume *physHG_BHPadLayer;

		//	Back HCAL: Si Pad
		//
		G4Box *solidHG_BHPad;
		G4LogicalVolume *logicHG_BHPad;
		G4VPhysicalVolume *physHG_BHPad;

		//	Back HCAL: Pad Readout
		//
		G4Box *solidHG_BHPadReadout;
		G4LogicalVolume *logicHG_BHPadReadout;
		G4VPhysicalVolume *physHG_BHPadReadout;

		
		//
		//	HE Volumes
		//
		G4Box *solidHE_Module;
		G4LogicalVolume *logicHE_Module;
		G4VPhysicalVolume *physHE_Module;

		//	HE Layer
		//
		G4Box *solidHE_Layer;
		G4LogicalVolume *logicHE_Layer;
		G4VPhysicalVolume *physHE_Layer;

		//	Abs
		//
		G4Box *solidHE_Abs;
		G4LogicalVolume *logicHE_Abs;
		G4VPhysicalVolume *physHE_Abs;

		//	Active Part
		//
		G4Box *solidHE_Act;
		G4LogicalVolume *logicHE_Act;
		G4VPhysicalVolume *physHE_Act;	

};
#endif











