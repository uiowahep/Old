#ifndef SHEVENTACTION_H
#define SHEVENTACTION_H 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4Timer.hh"

class SHEventAction : public G4UserEventAction
{
	public:
		SHEventAction();
		~SHEventAction();

		void BeginOfEventAction(const G4Event*);
		void EndOfEventAction(const G4Event*);

		G4Timer *timer;
};

#endif
