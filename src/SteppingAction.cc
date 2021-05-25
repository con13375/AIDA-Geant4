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
  G4double particleEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
  //G4double particleDirection = aStep->GetPreStepPoint()->GetMomentumDirection();

  G4double time_res = 2000; // event window is 2 microseconds; or it could be 20 nanoseconds corresponding to 50MHz
  G4double energy_res = 0; // set at zero for now so no effect

  if (currentPhysicalName == "DDSD"){// and abs(particleCharge) > 0){// <= 0.511
    G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
    G4double time = aStep->GetPostStepPoint()->GetGlobalTime();
    if (eventAction->FIRST == 0){eventAction->previous_time = time; eventAction->first_time = time; eventAction->FIRST += 1;} // only the first time, previous time is set to current time before timestamp check
    if (time-eventAction->previous_time > time_res){eventAction->timestamp += 1;}
    eventAction->previous_time = time;

    // This part transforms position to strip number (x,y), number of plaque (z)
    G4double detector_XY = 38.15;
    G4double first_pos = 38.7, plaque_sep = 11.6;
    G4double max_Z = first_pos+0.5*plaque_sep;
    G4double min_Z = first_pos-5.5*plaque_sep;
    G4int N_x = discretize(position[0], -detector_XY, detector_XY, 128);
    G4int N_y = discretize(position[1], -detector_XY, detector_XY, 128);
    G4int N_z = discretize(-position[2], -max_Z, -min_Z, 6);

    G4double EdepStep1 = aStep->GetTotalEnergyDeposit();
    if (EdepStep1 > energy_res and time-eventAction->first_time < time_res and 0 < N_x and N_x < 129 and 0 < N_y and N_y < 129){
      std::cout <<  std::fixed << std::setprecision(9) 
                << "##" << "," << particleMass << "," << particleCharge << "," << particleEnergy << "," 
 	        << EdepStep1 <<// "," << position[0] << "," << position[1] << "," << position[2] <<
                "," << N_x+1 << "," << N_y+1 << "," << N_z+1 << "," << eventAction->first_time << "," << time-eventAction->first_time << std::endl;
 
      eventAction->addEdep(EdepStep1, N_z, N_y, N_x, time, particleCharge, particleEnergy);//if (timestamp == 0){}
      }
   }; 
  
  G4String ParticleName = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  if(ParticleName == "geantino"){
    G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
    std::cout << "######" << position << "," << currentPhysicalName << "," << currentMaterialName << std::endl;
  }

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


