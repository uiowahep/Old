#ifndef SHRUNACTION_H
#define SHRUNACTION_H 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"
#include "G4Timer.hh"

class SHRunAction : public G4UserRunAction
{
	public:
		SHRunAction();
		~SHRunAction();

		void BeginOfRunAction(const G4Run* aRun);
		void EndOfRunAction(const G4Run* aRun);

		G4Timer *timer;
};

#endif
