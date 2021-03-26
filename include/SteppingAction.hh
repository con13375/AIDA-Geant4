// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingAction.hh,v 1.3 1999/12/15 14:49:20 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//////////////////////////////////////////////////////////////////////
//
//  Modified: 31-01-2006   J.L. Tain  --> testSteppingAction.hh
//  Modified: 17-01-2013  A. Algora to run with newer versions
//
//  Get energy deposited in sensitive volume
//
 

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class EventAction;

class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction(EventAction*);
    ~SteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    EventAction* eventAction;
};

#endif
