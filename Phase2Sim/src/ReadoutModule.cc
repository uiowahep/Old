
#include "HGCReadoutModule.hh"

ReadoutModule::ReadoutModule()
{
	_file = NULL;
}

ReadoutModule::ReadoutModule(TFile *rootFile, RunParams runParams)
{
	_file = rootFile;
	_runParams = runParams;
}

ReadoutModule::~ReadoutModule()
{

}

void ReadoutModule::BeginEvent(int)
{
	return;
}

void ReadoutModule::FinishEvent(int)
{
	return;
}

void ReadoutModule::Book(int)
{
	return;
}
