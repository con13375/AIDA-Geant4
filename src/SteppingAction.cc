// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingAction.cc,v 1.3 1999/12/15 14:49:22 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
/////////////////////////////////////////////////////////////////////
//
//  Modified: 17-01-2013  A. Algora to run with newer versions
// 
//  Collect the energy deposited in the sensitive test volume
///////////////////////////////////////////////////////////////////

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"

SteppingAction::SteppingAction(EventAction* EvAct) 
:G4UserSteppingAction(),eventAction(EvAct)
{ }

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // collect the energy deposited in the sensitive volume
  
  G4TouchableHandle touch = 
    aStep->GetPreStepPoint()->GetTouchableHandle();  
  const G4String currentMaterialName = 
    touch->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();

  const G4String currentPhysicalName 
    = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  //std::cout << " detector1 =" << currentMaterialName << std::endl;
  const G4double particleMass = aStep->GetPreStepPoint()->GetMass();

  if (currentPhysicalName == "DDSD" and 0.510 <= particleMass){// <= 0.511
    G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();

    G4double EdepStep1 = aStep->GetTotalEnergyDeposit();
     std::cout << "A mass of " << particleMass << "MeV at detector1 =" << currentMaterialName << " left Edep =" 
 	       << EdepStep1 << "MeV at (mm): (" << position << ")" << std::endl;
 
    if (EdepStep1 > 0.) eventAction->addEdep(EdepStep1);
   };  

/* G4double EdepStep = aStep->GetTotalEnergyDeposit();
 G4cout << " material =" << currentMaterialName << " Edep =" 
 	       << EdepStep << G4endl;
if (currentMaterialName == "Plastik"){
    G4double EdepStep1 = aStep->GetTotalEnergyDeposit();
    if (EdepStep1 > 0.) eventAction->addEdep(EdepStep1);
   }; 
if (currentMaterialName == "CF4"){
    G4double EdepStep3 = aStep->GetTotalEnergyDeposit();
    if (EdepStep3 > 0.) eventAction->addEdep3(EdepStep3);
   }; 
*/
}


