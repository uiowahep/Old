/*
 *	SHashlik Simulation
 *	Author: VK
 *	Date: 05/19/204
 *	Description:
 *
 * */


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "Randomize.hh"

#include "QBBC.hh"
//#include "QGSP.hh"

#include "SHDetectorConstruction.hh"
#include "SHPrimaryGeneratorAction.hh"
#include "SHEventAction.hh"
#include "SHRunAction.hh"
#include "SHDefs.hh"
#include "SHPhysicsList.hh"
#include "SHSteppingAction.hh"
#include "HGCReadoutModule.hh"

#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

using namespace std;
using namespace CLHEP;

/*
 *	Main Function:
 *	Input Format:
 *		
 *
 */
int main(int argc, char** argv)
{
	//	Parse all the input
	//
	RunParams runParams;
	runParams.isInteractive = atoi(argv[1]);
	runParams.iEnergy = 0;
	runParams.numModules = 5;
	runParams.numLayers = 28;
	runParams.numEvents = 10;
	runParams.seed = 1;
	runParams.verbosity = 5;
	runParams.hgCalInputFileName = "hgCalInput.txt";
	runParams.iPrim = 0;
	runParams.heInputFileName = "heInput.txt";
	runParams.minBiasFileName = "/Users/vkhristenko/Projects/Simulations/Phase2Sim/Phase2Sim-V3/GEN_PRIM/MinBias.root";
	runParams.sigFileName = "/Users/vkhristenko/Projects/Simulations/Phase2Sim/Phase2Sim-V3/GEN_PRIM/Higgs.root";
	runParams.startEvent = 0;
//	runParams.endcapConfigFileName = "/home/vitya/Projects/Simulations/Phase2Sim/Phase2Sim-V0-build/endcap.hgcal.cfg";
	runParams.endcapConfigFileName = "endcap.she.cfg";
	runParams.endcapType = 2; 

	if (!runParams.isInteractive)
	{
//		runParams.iEnergy = atoi(argv[2]);
//		runParams.numLayers = atoi(argv[3]);
//		runParams.numModules = atoi(argv[4]);
		runParams.numEvents = atoi(argv[2]);
		runParams.seed = atoi(argv[3]);
		runParams.verbosity = atoi(argv[4]);
		runParams.endcapType = atoi(argv[5]);
		runParams.endcapConfigFileName = argv[6];
		runParams.iPrim = atoi(argv[7]);
		runParams.minBiasFileName = argv[8];
		runParams.sigFileName = argv[9];
		runParams.startEvent = atoi(argv[10]);
//		runParams.hgCalInputFileName = argv[8];
//		runParams.iPrim = atoi(argv[9]);
//		runParams.heInputFileName = argv[10];
	}

	//	Initialize Random Number Generator
	//
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
//	CLHEP::HepRandom::setTheSeed(time(NULL));
	CLHEP::HepRandom::setTheSeed(runParams.seed);
	CLHEP::HepRandom::showEngineStatus();

	//
	//	Initialize ROOT Output Format: EM, FH. BH parts are separate
	//
	char rootOutFileName[200];
	sprintf(rootOutFileName, "CMSEndcap_%d_%d_%d_%d_%d.root",
			runParams.endcapType, runParams.verbosity, runParams.seed,
			runParams.numEvents, runParams.startEvent);
	TFile *rootOutFile = new TFile(rootOutFileName, "RECREATE");
	HGCReadoutModule *readout = new HGCReadoutModule(rootOutFile, runParams);

/*
	rootOutFile->mkdir("EM");
	rootOutFile->cd("EM");
	TTree *emTree = new TTree("Events", "Events");
	rootOutFile->cd();
	rootOutFile->mkdir("PRIM");
	rootOutFile->cd("PRIM");
	TTree *primTree = new TTree("Events", "Events");
	rootOutFile->cd();
	runParams.rootFile = rootOutFile;
	vector<TTree*> vTrees;
	vTrees.push_back(emTree);
	*/

	//	Initialize ROOT Output/ EM, FH, BH parts are separate...
	//	As of V5.
	//
/*	char rootOutFile[200];
	sprintf(rootOutFile, "HE_%d_%d_%d_%d.root", runParams.iEnergy,
			runParams.seed, runParams.numEvents, runParams.iPrim);
	TFile *shFile = new TFile(rootOutFile, "RECREATE");
	shFile->mkdir("EM");
	shFile->cd("EM");
	TTree *emTree = new TTree("detData", "detData");
	shFile->cd();
	shFile->mkdir("FH");
	shFile->cd("FH");
	TTree *fhTree = new TTree("detData", "detData");
	shFile->cd();
	shFile->mkdir("BH");
	shFile->cd("BH");
	TTree *bhTree = new TTree("detData", "detData");
	shFile->cd();
	shFile->mkdir("PRIM");
	shFile->cd("PRIM");
	TTree *primTree = new TTree("primData", "primData");
	shFile->cd();
	*/
//	runParams.rootFile = shFile;

	//	Initialize G4 Classes
	//
	G4RunManager *runManager = new G4RunManager;
	SHDetectorConstruction *detector = 
		new SHDetectorConstruction(runParams, readout); 
//		new SHDetectorConstruction(runParams, emTree, fhTree, bhTree);
	runManager->SetUserInitialization(detector);
	
	G4PhysListFactory factory;
	G4VModularPhysicsList *phys = factory.GetReferencePhysList("FTFP_BERT");
//	phys->RegisterPhysics(new SHPhysicsList);
//	SHPhysicsList *phys = new SHPhysicsList();
	runManager->SetUserInitialization(phys);

	SHPrimaryGeneratorAction *primAction = 
		new SHPrimaryGeneratorAction(runParams, readout);
	runManager->SetUserAction(primAction);

	// optional user action Classes
	//
	SHEventAction *eventAction = new SHEventAction();
	SHRunAction *runAction = new SHRunAction();
//	SHSteppingAction *steppingAction = new SHSteppingAction(runParams);
	runManager->SetUserAction(runAction);
	runManager->SetUserAction(eventAction);
//	runManager->SetUserAction(steppingAction);

	// Init G4 Kernel
	//
	runManager->Initialize();
	G4UImanager *UI = G4UImanager::GetUIpointer();

	if (runParams.isInteractive)
	{
		G4VisManager *visManager = 0;
		#ifdef G4VIS_USE
			//	Visualization
			//
			visManager = new G4VisExecutive;
			visManager->Initialize();
		#endif

		#ifdef G4UI_USE
			G4UIExecutive *ui = new G4UIExecutive(argc, argv);
			ui->SessionStart();
			delete ui;
		#endif
	}
	else
	{
		char cmd[200];
		sprintf(cmd, "/run/beamOn %d", runParams.numEvents);
		UI->ApplyCommand(cmd);
	}

	//	Finalize
	//
	delete runManager;
	rootOutFile->Print();
	rootOutFile->Write();
	rootOutFile->Close();

	return 0;
} 
