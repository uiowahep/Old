#ifndef CSTEPPINGACTION_H
#define CSTEPPINGACTION_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

#include "SHDefs.hh"

class SHSteppingAction : public G4UserSteppingAction
{
  public:
    SHSteppingAction(RunParams);
    virtual ~SHSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

	RunParams runParams;
};

#endif

