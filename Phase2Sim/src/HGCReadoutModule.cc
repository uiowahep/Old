//
//	HGCal Readout Module
//

#include "HGCReadoutModule.hh"
#include "G4UnitsTable.hh"

#include <math.h>
#include <iostream>

using namespace std;
using namespace CLHEP;


const Double_t HGCMIPRESPPERUM = 80.;
const Double_t PADZ1 = 9.;
const Double_t PADZ2 = 18;
const Double_t ETAMAX = 3.0;
const Double_t ETAMIN = 1.305;
const Double_t DETA = 0.087;
const Double_t DETA_32T33 = 0.175;
const Double_t DETA_33TEND = 0.128;
const Double_t DPHI = 5;
const Double_t PHIMAX = 360.;
const Double_t PHIMIN = 0;
const Double_t DX_SH = 14;


//
//	Constructor by Default
//
HGCReadoutModule::HGCReadoutModule()
{
	
}

//
//	Specialied Constructor:
//	Detector Output Description:
//	-- HGCal has 3 important sections: 3 SubDetectors
//		-- EM
//		-- FH
//		-- HE
//	
HGCReadoutModule::HGCReadoutModule(TFile *rootFile, RunParams runParams)
	: ReadoutModule(rootFile, runParams)
{
	cout << "### Configuring Readout Module..." << endl;
	//
	//	Declare the names of 3 Readout Parts:
	//	NOTE: SHE and HE go together
	//
	string nameEMP = "EMP"; int idEMP = 0;
	string nameEMM = "EMM"; int idEMM = 1;
	string nameFHP = "FHP"; int idFHP = 2;
	string nameFHM = "FHM"; int idFHM = 3;
	string nameBHP = "BHP"; int idBHP = 4;
	string nameBHM = "BHM"; int idBHM = 5;
	string namePRIM = "PRIM"; int idPRIM = 10;
	string nameGEOM = "GEOM"; int idGEOM = 11;
	string nameSHP = "SHP"; int idSHP = 6;
	string nameSHM = "SHM"; int idSHM = 7;
	string nameHEP = "HEP"; int idHEP = 8;
	string nameHEM = "HEM"; int idHEM = 9;

	_file->mkdir(nameEMP.c_str());
	_file->mkdir(nameEMM.c_str());
	_file->mkdir(nameFHP.c_str());
	_file->mkdir(nameFHM.c_str());
	_file->mkdir(nameBHP.c_str());
	_file->mkdir(nameBHM.c_str());
	_file->mkdir(namePRIM.c_str());
	_file->mkdir(nameGEOM.c_str());
	_file->mkdir(nameSHP.c_str());
	_file->mkdir(nameSHM.c_str());
	_file->mkdir(nameHEP.c_str());
	_file->mkdir(nameHEM.c_str());

	//
	//	Initialize the Trees
	//
	_file->cd(nameEMP.c_str());
//	TTree *tree = new TTree("Events", "Events");
	SubDetector sub;
	sub.name = nameEMP;
	sub.id = idEMP;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);

	_file->cd(nameEMM.c_str());
//	tree = new TTree("Events", "Events");
	sub.name = nameEMM;
	sub.id = idEMM;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);

	_file->cd(nameFHP.c_str());
//	tree = new TTree("Events", "Events");
	sub.name = nameFHP;
	sub.id = idFHP;
	sub.dataTree = new TTree("Events", "Events");
	sub.responses = new vector<Double_t>;
	sub.channelIDs = new vector<int>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);

	_file->cd(nameFHM.c_str());
//	tree = new TTree("Events", "Events");
	sub.name = nameFHM;
	sub.id = idFHM;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);

	_file->cd(nameBHP.c_str());
//	tree = new TTree("Events", "Events");
	sub.name = nameBHP;
	sub.id = idBHP;
	sub.dataTree = new TTree("Events", "Events");
	sub.responses = new vector<Double_t>;
	sub.channelIDs = new vector<int>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);


	_file->cd(nameBHM.c_str());
//	tree = new TTree("Events", "Events");
	sub.name = nameBHM;
	sub.id = idBHM;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);


	//
	//	PRIM Detector
	//
	_file->cd(namePRIM.c_str());
	_primDet.name = namePRIM;
	_primDet.id = idPRIM;
	_primDet.dataTree = new TTree("Events", "Events");
	_primDet.x = new vector<Double_t>;
	_primDet.y = new vector<Double_t>;
	_primDet.z = new vector<Double_t>;
	_primDet.px = new vector<Double_t>;
	_primDet.py = new vector<Double_t>;
	_primDet.pz = new vector<Double_t>;
	_primDet.ene = new vector<Double_t>;
	_primDet.pname = new vector<string>;
	_primDet.pid = new vector<int>;

	//
	//	GEOM "Detector" 
	//	Gotta keep Run Detector info in the ROOT file
	//
	_file->cd(nameGEOM.c_str());
	if (_runParams.endcapType == 1)
		_GMaps.tree = new TTree("GEOM", "GEOM");
	else if (_runParams.endcapType == 2)
		_SHMaps.tree = new TTree("GEOM", "GEOM");
	else
		cout << "### UNKNOWN READOUT CONFIGURATION ID!" << endl;


	//
	//	SHP
	//
	_file->cd(nameSHP.c_str());
	sub.name = nameSHP;
	sub.id = idSHP;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);

	//
	//	SHM
	//
	_file->cd(nameSHM.c_str());
	sub.name = nameSHM;
	sub.id = idSHM;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);

	//
	//	HE for Shashlik + HE
	//
	_file->cd(nameHEP.c_str());
	sub.name = nameHEP;
	sub.id = idHEP;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);

	_file->cd(nameHEM.c_str());
	sub.name = nameHEM;
	sub.id = idHEM;
	sub.dataTree = new TTree("Events", "Events");
	sub.channelIDs = new vector<int>;
	sub.responses = new vector<Double_t>;
	sub.x = new vector<Double_t>;
	sub.y = new vector<Double_t>;
	_subDets.push_back(sub);
	

	_file->cd();

	cout << "### Done Configuring Readout Module..." << endl;
}


//	
//	Destructor
//
HGCReadoutModule::~HGCReadoutModule()
{

}

//
//	Do it in the beginning of each event
//
void HGCReadoutModule::BeginEvent(int id)
{
	if (_runParams.verbosity>1)
		cout << "### ReadoutModule::BeginEvent id=" << id << endl;

	if (_primDet.id == id)
	{
		//
		//	For PRIM Detector clean up
		//
		_primDet.x->clear();
		_primDet.y->clear();
		_primDet.z->clear();
		_primDet.px->clear();
		_primDet.py->clear();
		_primDet.pz->clear();
		_primDet.ene->clear();
		_primDet.pname->clear();
		_primDet.pid->clear();
	}
	else
	{
		//
		//	Endcap SubDetector clean up
		//
		_subDets[id].channelIDs->clear();
		_subDets[id].responses->clear();
		_subDets[id].x->clear();
		_subDets[id].y->clear();
	}

	return;
}

//
//	Do it in the end of each event.
//
void HGCReadoutModule::FinishEvent(int id)
{
	if (_runParams.verbosity>1)
		cout << "### Finishing Event id=" << id << endl;

	if (_primDet.id == id)
		_primDet.dataTree->Fill();
	else
		_subDets[id].dataTree->Fill();

	return;
}

//
//	Book Histograms for that subdetector
//
void HGCReadoutModule::Book(int id)
{
	if (_runParams.verbosity>1)
		cout << "### Booking Histos for id: " << id << endl;

	if (_primDet.id == id)
	{
		//
		//	Do all the branching for PRIM Detector
		//
		TTree *tree = _primDet.dataTree;

		tree->Branch("x", &_primDet.x);
		tree->Branch("y", &_primDet.y);
		tree->Branch("z", &_primDet.z);
		tree->Branch("px", &_primDet.px);
		tree->Branch("py", &_primDet.py);
		tree->Branch("pz", &_primDet.pz);
		tree->Branch("ene", &_primDet.ene);
		tree->Branch("pname", &_primDet.pname);
		tree->Branch("pid", &_primDet.pid);
	}
	else
	{
		//
		//	Do all the branching for Endcap Sub Detector
		//
		TTree *tree = _subDets[id].dataTree;
	
		tree->Branch("channelIDs", &_subDets[id].channelIDs);
		tree->Branch("responses", &_subDets[id].responses);
//		tree->Branch("x", &_subDets[id].x);
//		tree->Branch("y", &_subDets[id].y);
	}
}

//
//	Push a RPIM Hit
//
void HGCReadoutModule::PushHit(PRIMHit aHit)
{
	_primDet.x->push_back(aHit.pos.x()/mm);
	_primDet.y->push_back(aHit.pos.y()/mm);
	_primDet.z->push_back(aHit.pos.z()/mm);
	_primDet.px->push_back(aHit.mom.x()/GeV);
	_primDet.py->push_back(aHit.mom.y()/GeV);
	_primDet.pz->push_back(aHit.mom.z()/GeV);
	_primDet.ene->push_back(aHit.ene/GeV);
	_primDet.pname->push_back(aHit.name.data());
	_primDet.pid->push_back(aHit.id);
}

//
//	Push a Hit for HE
//
void HGCReadoutModule::PushHit(HEHit aHit)
{
	int channel = -1;

	//
	//	Compute the REadout Channel
	//
	if (_runParams.endcapType == 1)
		channel = ComputeHEChannel(aHit.pos, aHit.layerNumber);
	else if	(_runParams.endcapType == 2)
		channel = ComputeSHHEChannel(aHit.pos, aHit.layerNumber);
	else
		cout << "### UNKNOWN READOUT CONFIGURATION ID!" << endl;


	if (channel<0)
		cout << "### NEGATIVE CHANNEL: DETID: " << aHit.detID << endl;

	if (_runParams.verbosity>0)
		cout << "### HEHit Info: " << aHit.detID << "  " 
			<< aHit.pos.x()/mm << "  " << aHit.pos.y()/mm << "  " 
			<< aHit.pos.z()/mm << endl;

	//
	//	Push
	//
	Double_t response = aHit.numPhotons;
	if (!IsChannelTriggered(channel, aHit.detID, aHit))
	{
		_subDets[aHit.detID].channelIDs->push_back(channel);
		_subDets[aHit.detID].responses->push_back(response);
	}

}

void HGCReadoutModule::PushHit(SHHit aHit)
{
	int channel = -1;

	//
	//	Compute the Channel First
	//	
	channel = ComputeSHChannel(aHit.pos, aHit.layerNumber);

	if (channel <= 0)
		cout << "### NEGATIVE OR CHANNEL: DETID: " << aHit.detID << endl;

	if (_runParams.verbosity>1)
		cout << "### Shashlik Hit:Channel: " << channel << endl;

	//
	//	Push
	//
	Double_t response = aHit.numPhotons;
	if (!IsChannelTriggered(channel, aHit.detID, aHit))
	{
		_subDets[aHit.detID].channelIDs->push_back(channel);
		_subDets[aHit.detID].responses->push_back(response);
	}
}


//
//	Push a Hit in the Detector
//
void HGCReadoutModule::PushHit(HGCHit aHit)
{
	int channel = -1;
	
	if (_runParams.verbosity>0)
		cout << "### HGCHit Info: " << aHit.detID << "  " 
			<< aHit.glPrePos.x()/mm << "  " << aHit.glPrePos.y()/mm 
			<< "  " << aHit.glPrePos.z()/mm << endl;

	//
	//	0 for EM part, 2 for FH part, 4 will be for HE(BH) part
	//
	if (aHit.detID == 0 || aHit.detID==1)
		channel = ComputeEMChannel(aHit.glPrePos, aHit.layerNumber);
	else if (aHit.detID == 2 || aHit.detID==3)
		channel = ComputeFHChannel(aHit.glPrePos, aHit.layerNumber);
	else if (aHit.detID == 4)
	{
		return;
	}

	if (channel<0)
		cout << "### NEGATIVE CHANNEL: DETID: " << aHit.detID << endl;


	//
	//	Get the step length - proportional to the response
	//	Fill the appropriate detector information
	//
	Double_t response = aHit.dPos.mag()/um * HGCMIPRESPPERUM;
	if (!IsChannelTriggered(channel, aHit.detID, aHit))
	{
		_subDets[aHit.detID].channelIDs->push_back(channel);
		_subDets[aHit.detID].responses->push_back(response);
		_subDets[aHit.detID].x->push_back(aHit.glPrePos.x()/mm);
		_subDets[aHit.detID].y->push_back(aHit.glPrePos.y()/mm);
	}

}

//
//	Check if the channel has been triggered already...
//	
bool HGCReadoutModule::IsChannelTriggered(int channel, int detID, HGCHit aHit)
{
	if (_subDets[detID].channelIDs->size() == 0)
		return false;

	int ii = 0;
	for (vector<int>::iterator it=_subDets[detID].channelIDs->begin(); 
			it!=_subDets[detID].channelIDs->end(); ++it)
	{
		if (channel == *it)
		{
			Double_t &response = _subDets[detID].responses->at(ii);
			response += aHit.dPos.mag()/um * HGCMIPRESPPERUM;
			return true;
		}
		ii++;
	}
			
	return false;
}

//
//	Check if the channel has been triggered yet...
//
bool HGCReadoutModule::IsChannelTriggered(int channel, int detID,
		HEHit aHit)
{
	if (_subDets[detID].channelIDs->size() == 0)
		return false;

	int ii = 0;
	for (vector<int>::iterator it=_subDets[detID].channelIDs->begin();
			it!=_subDets[detID].channelIDs->end(); ++it)
	{
		if (channel == *it)
		{
			Double_t &response = _subDets[detID].responses->at(ii);
			response += aHit.numPhotons;
			return true;
		}
		ii++;
	}

	return false;
}

//
//	Push Geometry Information
//	and Initialize Geometry Maps
//
void HGCReadoutModule::PushGeomInfo(CMSHGCal hgCal)
{
	_HGCal = hgCal;	

	//
	//	Init GMaps
	//
	double EMstartZ = _HGCal.startZ/mm;
	double FHstartZ = EMstartZ + 2*(_HGCal.EM.centerZ/mm - EMstartZ);
	int EMNLayers = _HGCal.EM.n1 + _HGCal.EM.n2 + _HGCal.EM.n3;
	int FHNLayers = _HGCal.FH.n;
	int HENLayers = _HGCal.SHE.n + _HGCal.BH.n;

	_GMaps.pL2R_EM = new Double_t[EMNLayers];
	_GMaps.pL2Z_EM = new Double_t[EMNLayers];
	_GMaps.pL2NChs1D_EM = new int[EMNLayers];
	_GMaps.pL2R_FH = new Double_t[FHNLayers];
	_GMaps.pL2Z_FH = new Double_t[FHNLayers];
	_GMaps.pL2NChs1D_FH = new int[FHNLayers];
	_GMaps.pEtaMap = new Double_t[NUMETAS];
	_GMaps.pL2Z_HE = new Double_t[HENLayers];
	_GMaps.nEM = EMNLayers;
	_GMaps.n1EM = _HGCal.EM.n1;
	_GMaps.n2EM = _HGCal.EM.n2;
	_GMaps.n3EM = _HGCal.EM.n3;
	_GMaps.nFH = FHNLayers;
	_GMaps.nHE = HENLayers;

	double EMabsDz = 0;
	double dx = 0;
	double dy = 0;
	double FHabsDz = _HGCal.FH.fullAbsZ/mm;
	double padDz = _HGCal.EM.fullPadZ/mm;
	double readoutDz = _HGCal.EM.fullReadoutZ/mm;

	//
	//	Map EM
	//
	for (int iLayer=0; iLayer<EMNLayers; iLayer++)
	{
		if (iLayer<_HGCal.EM.n1)
		{
			dx = PADZ1;
			dy = PADZ1;
			EMabsDz = _HGCal.EM.fullAbsZ[0]/mm;
		}
		else if (iLayer < (_HGCal.EM.n1 + _HGCal.EM.n2)) 
		{
			dx = PADZ1;
			dy = PADZ1;
			EMabsDz = _HGCal.EM.fullAbsZ[1]/mm;
		}
		else
		{
			dx = PADZ2;
			dy = PADZ2;
			EMabsDz = _HGCal.EM.fullAbsZ[2]/mm;
		}

		Double_t dLayer = EMabsDz + padDz + readoutDz;
		Double_t z = EMstartZ + iLayer*dLayer + padDz + EMabsDz;
		Double_t r = _HGCal.EM.rmax1/mm * z/EMstartZ;
		int nChs1D = ceil(2*r/dx);

		//
		//	Assing the radius and the #channels
		//	for that particular layer...
		//
		_GMaps.pL2R_EM[iLayer] = r;
		_GMaps.pL2NChs1D_EM[iLayer] = nChs1D;
		_GMaps.pL2Z_EM[iLayer] = z;
	}

	//
	//	Map FH
	//
	for (int iLayer=0; iLayer<FHNLayers; iLayer++)
	{
		dx = PADZ2;
		dy = PADZ2;

		Double_t dLayer = FHabsDz + padDz + readoutDz;
		Double_t z = FHstartZ + iLayer*dLayer + padDz + FHabsDz;
		Double_t r = _HGCal.FH.rmax1/mm * z/FHstartZ;
		int nChs1D = ceil(2*r/dx);

		//
		//	Assign the radius and the #channels for that 
		//	particular layer to the map...
		//
		_GMaps.pL2R_FH[iLayer] = r;
		_GMaps.pL2Z_FH[iLayer] = z;
		_GMaps.pL2NChs1D_FH[iLayer] = nChs1D;
	}

	//
	//	Map HE:
	//	-- generate Eta Divisions
	//	-- Map layer to Z
	//	The actual NUMBER OF ETAS/PHIS is Hard-Coded for now
	//	Also, this is the number of 
	//
	_GMaps.numEtas = 18;
	_GMaps.numPhis = 72;
	for (int ieta=0; ieta<_GMaps.numEtas-1; ieta++)
	{
		_GMaps.pEtaMap[ieta] = ETAMIN + ieta*DETA;
	}
	_GMaps.pEtaMap[_GMaps.numEtas-1] = 
		_GMaps.pEtaMap[_GMaps.numEtas-2] + DETA_32T33;
	_GMaps.pEtaMap[_GMaps.numEtas] = 
		_GMaps.pEtaMap[_GMaps.numEtas-1] + DETA_33TEND;
	

	//
	//	Branch them now..
	//
	_GMaps.tree->Branch("nEM", &_GMaps.nEM, "nEM/I");
	_GMaps.tree->Branch("n1EM", &_GMaps.n1EM, "n1EM/I");
	_GMaps.tree->Branch("n2EM", &_GMaps.n2EM, "n2EM/I");
	_GMaps.tree->Branch("n3EM", &_GMaps.n3EM, "n3EM/I");
	_GMaps.tree->Branch("L2R_EM", _GMaps.pL2R_EM, "L2R_EM[nEM]/D");
	_GMaps.tree->Branch("L2Z_EM", _GMaps.pL2Z_EM, "L2Z_EM[nEM]/D");
	_GMaps.tree->Branch("L2NChs1D_EM", _GMaps.pL2NChs1D_EM, "L2NChs1D_EM[nEM]/I");
	_GMaps.tree->Branch("nFH", &_GMaps.nFH, "nFH/I");
	_GMaps.tree->Branch("L2R_FH", _GMaps.pL2R_FH, "L2R_FH[nFH]/D");
	_GMaps.tree->Branch("L2Z_FH", _GMaps.pL2Z_FH, "L2Z_FH[nFH]/D");
	_GMaps.tree->Branch("L2NChs1D_FH", _GMaps.pL2NChs1D_FH, "L2NChs1D_FH[nFH]/I");
	_GMaps.tree->Branch("numEtas_HE", &_GMaps.numEtas, "numEtas_HE/I");
	_GMaps.tree->Branch("numPhis_HE", &_GMaps.numPhis, "numPhis_HE/I");
	_GMaps.tree->Branch("pEtaMap",_GMaps.pEtaMap, "pEtaMap[20]/D");
	_GMaps.tree->Branch("nBH", &_GMaps.nHE, "nBH/I");

	_GMaps.tree->Fill();
}


//
//	Compute EM channel
//
int HGCReadoutModule::ComputeEMChannel(G4ThreeVector pos, int ilayer)
{
	int sumChs = 0;
	
	//
	//	Compute the channel#
	//
	for (int i=0; i<ilayer; i++)
		sumChs += _GMaps.pL2NChs1D_EM[i] * _GMaps.pL2NChs1D_EM[i];

	//
	//	Compute channel
	//
	double x = _GMaps.pL2R_EM[ilayer] + pos.x()/mm;
	double y = _GMaps.pL2R_EM[ilayer] + pos.y()/mm;
	double dx = PADZ1;
	double dy = PADZ1;
	if (ilayer>=(_HGCal.EM.n1 + _HGCal.EM.n2))
	{
		dx = PADZ2;
		dy = PADZ2;
	}

	int ix = x/dx;
	int iy = y/dy;

	//
	//	Debug
	//
	if (ix<0 || iy<0)
		cout << "### ERROR: Readout Bins for EM are NEGATIVE" << endl;

	int channel = sumChs + iy*_GMaps.pL2NChs1D_EM[ilayer] + ix;
	return channel;
}

//
//	Compute FH channel
//
int HGCReadoutModule::ComputeFHChannel(G4ThreeVector pos, int ilayer)
{
	int sumChs = 0;

	//
	//	Compute the #channels before the our one
	//
	for (int i=0; i<ilayer; i++)
		sumChs +=_GMaps.pL2NChs1D_FH[i] * _GMaps.pL2NChs1D_FH[i];

	//
	//	Compute the channel
	//
	double x = _GMaps.pL2R_FH[ilayer] + pos.x()/mm;
	double y = _GMaps.pL2R_FH[ilayer] + pos.y()/mm;
	double dx = PADZ2;
	double dy = PADZ2;
	int ix = x/dx;
	int iy = y/dy;

	//
	//	Debug
	//
	if (ix<0 || iy<0)
		cout << "### ERROR: Readout Bins for FH are NEGATIVE!!!" << endl;

	int channel = sumChs  + iy*_GMaps.pL2NChs1D_FH[ilayer] + ix;

	return channel;
}

//
//	Compute HE Channel
//
int HGCReadoutModule::ComputeHEChannel(G4ThreeVector pos, int layer)
{
	int sumChs = _GMaps.numEtas*_GMaps.numPhis*layer;
	double epsilon = 0.0001;
	
	double x = pos.x()/mm;
	double y = pos.y()/mm;
	double z = pos.z()/mm;

	double r = sqrt(x*x + y*y);
//	double phi = pos.phi();
	double phi = 180 + atan2(y,x)*180/3.14159265;
	double eta = abs(pos.eta());

	int iphi = phi/DPHI;
	int ieta = 0;
//	if (eta == _GMaps.pEtaMap[_GMaps.numEtas])
//		ieta = _GMaps.numEtas - 1;
	if (eta > (_GMaps.pEtaMap[_GMaps.numEtas] + epsilon))
		cout << "### ERROR: eta is out of bounds!!!" << ieta << "  "
			<< _GMaps.pEtaMap[_GMaps.numEtas] << endl;
	else if (eta > _GMaps.pEtaMap[_GMaps.numEtas-1])
		ieta = _GMaps.numEtas-1;
	else if (eta > _GMaps.pEtaMap[_GMaps.numEtas-2])
		ieta = _GMaps.numEtas-2;
	else
		ieta = floor((eta - _GMaps.pEtaMap[0])/DETA);

//	cout << "### " << ieta << "  " << eta << "  " << phi << "  " << iphi << endl;
	
	if (iphi<0)
		cout << "### NEGATIVE iphi=" << iphi << "  phi=" << phi << endl;

	int channel = sumChs + iphi*_GMaps.numEtas + ieta;

	return channel;
}

//
//	Compute Shashlik Channel
//
int HGCReadoutModule::ComputeSHChannel(G4ThreeVector pos, int ilayer)
{
	
	int sumChs = _SHMaps.pL2NChs1D_SH[_SHMaps.nSH-1]*
		_SHMaps.pL2NChs1D_SH[_SHMaps.nSH-1]*ilayer;

	double x = _SHMaps.pL2R_SH[_SHMaps.nSH-1] + pos.x()/mm;
	double y = _SHMaps.pL2R_SH[_SHMaps.nSH-1] + pos.y()/mm;
	double dx = DX_SH;
	double dy = DX_SH;
	
	int ix = x/dx;
	int iy = y/dy;

	if (ix<0 || iy<0)
		cout << "### ERROR: ReadoutBins for EM are NEGATIVE: "
			<< ix << "  " << iy << "  " << x << "  " << y << endl;

	int channel = sumChs + iy*_SHMaps.pL2NChs1D_SH[_SHMaps.nSH-1] + ix;

/*	
	int sumChs = _SHMaps.numEtas * _SHMaps.numPhis * ilayer;
	double epsilon = 0.0001;

	double x = pos.x()/mm;
	double y = pos.y()/mm;
	double z = pos.z()/mm;

	double r = sqrt(x*x + y*y);
	double phi = 180. + atan2(y,x)*180/3.14159265;
	double eta = abs(pos.eta());

	int iphi = phi/DPHI;
	int ieta = 0;
	if (eta > (_SHMaps.pEtaMap[_SHMaps.numEtas] + epsilon))
		cout << "### ERROR: eta is out of range: " << eta
			<< "  " << _SHMaps.pEtaMap[_SHMaps.numEtas] << endl;
	else if (eta > _SHMaps.pEtaMap[_SHMaps.numEtas-1])
		ieta = _SHMaps.numEtas - 1;
	else if (eta > _SHMaps.pEtaMap[_SHMaps.numEtas-2])
		ieta = _SHMaps.numEtas - 2;
	else
		ieta = floor((eta - _SHMaps.pEtaMap[0])/DETA);

	if (iphi<0)
		cout << "### NEGATIVE iphi=" << iphi << "  phi=" << phi << endl;

	int channel = sumChs + iphi*_SHMaps.numEtas + ieta;
*/

	return channel;
}

//
//	Compute HE channel but in Shashlik + HE configuration
//	Different Map
//
int HGCReadoutModule::ComputeSHHEChannel(G4ThreeVector pos, int layer)
{
	int sumChs = _SHMaps.numEtas*_SHMaps.numPhis*layer;
	double epsilon = 0.0001;

	double x = pos.x()/mm;
	double y = pos.y()/mm;
	double z = pos.z()/mm;

	double r = sqrt(x*x + y*y);
	double phi = 180. + atan2(y,x)*180/3.14159265;
	double eta = abs(pos.eta());

	int iphi = phi/DPHI;
	int ieta = 0;
//	if (eta == _SHMaps.pEtaMap[_SHMaps.numEtas])
//		ieta = _SHMaps.numEtas-1;
	if (eta > (_SHMaps.pEtaMap[_SHMaps.numEtas] + epsilon))
		cout << "### ERROR: eta is out of range: "	<< eta 
			<< "  " << _SHMaps.pEtaMap[_SHMaps.numEtas] << endl;
	else if (eta > _SHMaps.pEtaMap[_SHMaps.numEtas-1])
		ieta = _SHMaps.numEtas-1;
	else if (eta > _SHMaps.pEtaMap[_SHMaps.numEtas-2])
		ieta = _SHMaps.numEtas-2;
	else 
		ieta = floor((eta - _SHMaps.pEtaMap[0])/DETA);

//	cout << "### " << ieta << "  " << eta << "  " << phi << "  " << iphi << endl;

	if (iphi<0)
		cout << "### NEGATIVE iphi=" << iphi << "  phi=" << phi << endl;

	int channel = sumChs + iphi*_SHMaps.numEtas + ieta;

	return channel;
}

//
//	Check wherther the channel has been triggered
//
bool HGCReadoutModule::IsChannelTriggered(int channel, int detID, 
		SHHit aHit)
{
	if (_subDets[detID].channelIDs->size() == 0)
		return false;

	int ii=0;
	for (vector<int>::iterator it=_subDets[detID].channelIDs->begin();
			it!=_subDets[detID].channelIDs->end(); ++it)
	{
		if (channel == *it)
		{
			Double_t &response = _subDets[detID].responses->at(ii);
			response += aHit.numPhotons;
			return true;
		}
		ii++;
	}

	return false;
}

void HGCReadoutModule::PushGeomInfo(CMSSHE she)
{
	_SHE = she;
	cout << _SHE.startZ/mm << endl;

	//
	//	Init SHMap
	//
	double SHstartZ = _SHE.startZ/mm;
	double HEstartZ = _SHE.startZ/mm + _SHE.EM.totalZ/mm;
	int SHNLayers = _SHE.EM.n;
	int HENLayers = _SHE.FHE.n + _SHE.MHE1.n + _SHE.MHE2.n + _SHE.BHE.n;

	_SHMaps.pL2Z_SH = new Double_t[SHNLayers+1];
	_SHMaps.pL2R_SH = new Double_t[SHNLayers+1];
	_SHMaps.pL2NChs1D_SH = new Int_t[SHNLayers+1];

	_SHMaps.pL2Z_HE = new Double_t[HENLayers];
	_SHMaps.pEtaMap = new Double_t[NUMETAS];
	_SHMaps.nHE = HENLayers;
	_SHMaps.nSH = SHNLayers+1;

	double EMabsZ = _SHE.EM.fullAbsZ/mm;
	double HEabsZ = _SHE.FHE.fullAbsZ/mm;
	double EMactZ = _SHE.EM.fullActZ/mm;

	//
	//	Map Shashlik XY for now
	//
	for (int iLayer=0; iLayer<_SHMaps.nSH; iLayer++)
	{
		Double_t dLayer = EMabsZ + EMactZ;
 		Double_t z = SHstartZ + iLayer*dLayer + EMactZ;
		Double_t r = _SHE.EM.rmax1/mm * z/SHstartZ;
		int nChs1D = ceil(2.*r/DX_SH);
		cout << iLayer << "  " << dLayer << "  " << z << "  " << r  << "  "
			<< nChs1D << endl;

		_SHMaps.pL2R_SH[iLayer] = r;
		_SHMaps.pL2NChs1D_SH[iLayer] = nChs1D;
		_SHMaps.pL2Z_SH[iLayer] = z;
	}

	//
	//	Map HE
	//
	_SHMaps.numEtas = 18;
	_SHMaps.numPhis = 72;
	for (int ieta=0; ieta<_SHMaps.numEtas-1; ieta++)
		_SHMaps.pEtaMap[ieta] = ETAMIN + ieta*DETA;
	_SHMaps.pEtaMap[_SHMaps.numEtas-1] = 
		_SHMaps.pEtaMap[_SHMaps.numEtas-2] + DETA_32T33;
	_SHMaps.pEtaMap[_SHMaps.numEtas] = 
		_SHMaps.pEtaMap[_SHMaps.numEtas-1] + DETA_33TEND;

	for (int ieta=0; ieta<=_SHMaps.numEtas; ieta++)
		cout << "### ieta=" << ieta << "  eta=" << _SHMaps.pEtaMap[ieta]
			<< endl;

	//
	//	Branch now
	//
	_SHMaps.tree->Branch("nSH", &_SHMaps.nSH, "nSH/I");
	_SHMaps.tree->Branch("nHE", &_SHMaps.nHE, "nHE/I");
	_SHMaps.tree->Branch("numEtas_HE", &_SHMaps.numEtas, "numEtas_HE/I");
	_SHMaps.tree->Branch("numPhis_HE", &_SHMaps.numPhis, "numPhis_HE/I");
	_SHMaps.tree->Branch("L2Z_SH", _SHMaps.pL2Z_SH, "L2Z_SH[nSH]/D");
	_SHMaps.tree->Branch("L2R_SH", _SHMaps.pL2R_SH, "L2R_SH[nSH]/D");
	_SHMaps.tree->Branch("L2NChs1D_SH", _SHMaps.pL2NChs1D_SH, 
			"L2NChs1D_SH[nSH]/I");
	_SHMaps.tree->Branch("L2Z_HE", _SHMaps.pL2Z_HE, "L2Z_HE[nHE]");
	_SHMaps.tree->Branch("pEtaMap", _SHMaps.pEtaMap, "pEtaMap[20]/D");

	_SHMaps.tree->Fill();

}








