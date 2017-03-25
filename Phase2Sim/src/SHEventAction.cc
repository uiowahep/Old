#include "SHEventAction.hh"

#include <iostream>

using namespace std;

SHEventAction::SHEventAction()
{
	timer = new G4Timer();
}

SHEventAction::~SHEventAction()
{
	delete timer;
}

void SHEventAction::BeginOfEventAction(const G4Event* anEvent)
{
	G4cout << "### Starting Event: " << anEvent->GetEventID() << endl;
	timer->Start();
}

void SHEventAction::EndOfEventAction(const G4Event*)
{
	timer->Stop();
	G4cout << "### Event Duration: " << *timer << G4endl;
}
