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

int discretize(G4double x, G4double a, G4double b, G4int N){
  G4double m = N/(b-a);
  if(x==b){
    return N;
  }
  else{
    return std::floor(m*(x-a));
  }
}

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
  G4double particleMass = aStep->GetPreStepPoint()->GetMass();
  G4double particleCharge = aStep->GetPreStepPoint()->GetCharge();

  if (currentPhysicalName == "DDSD" and 0.510 <= particleMass){// <= 0.511
    G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();

    // This part transforms position to strip number (x,y), number of plaque (z)
    G4double detector_XY = 38.15;
    G4double first_pos = 38.7, plaque_sep = 11.6;
    G4double max_Z = first_pos+0.5*plaque_sep;
    G4double min_Z = first_pos-5.5*plaque_sep;
    G4int N_x = discretize(position[0], -detector_XY, detector_XY, 128);
    G4int N_y = discretize(position[1], -detector_XY, detector_XY, 128);
    G4int N_z = discretize(-position[2], -max_Z, -min_Z, 6);

    G4double EdepStep1 = aStep->GetTotalEnergyDeposit();
     std::cout << "##" << "," << particleMass << "," << particleCharge << "," 
 	       << EdepStep1 << "," << position[0] << "," << position[1] << "," << position[2] <<
               "," << N_x+1 << "," << N_y+1 << "," << N_z+1 << std::endl;
 
    eventAction->addEdep(EdepStep1, N_z, N_y, N_x);
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


