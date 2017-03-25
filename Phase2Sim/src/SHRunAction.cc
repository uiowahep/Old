#include "SHRunAction.hh"

SHRunAction::SHRunAction()
{
	timer = new G4Timer;
}

SHRunAction::~SHRunAction()
{
	delete timer;
}

void SHRunAction::BeginOfRunAction(const G4Run* aRun)
{
	G4cout << "### Run " << aRun->GetRunID() << " start. " << G4endl;
	timer->Start();
}

void SHRunAction::EndOfRunAction(const G4Run *aRun)
{
	timer->Stop();
	G4cout << "Number of Events: " << aRun->GetNumberOfEvent()
		<< " took " << *timer << G4endl;
}
