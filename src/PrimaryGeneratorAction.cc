//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include <math.h>       /* pow */
#include <chrono>
#include <random>       
// std::uniform_int_distribution<int> distrib(a, b)
#include <vector>
//std::vector< int > arr;
//arr.push_back(1);

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
                    = particleTable->FindParticle("chargedgeantino");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(1*eV);    
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double FermiDistribution( G4double a ){
  // particle energy by fermi distribution (pending)
  G4double ep_energy = 3720000; // beta+ decay of 100Sn in eV

  G4double kT = 0.0254; // kT at room temperature in eV
  //kT = 170000; // what if this is the actual temperature
  G4double kT_norm = kT/ep_energy;
  G4double FD;
  FD = 1./(std::exp((a-1)/kT_norm)+1);
  //FD = 1./( pow((a-1, 2.7182) + 1. );
  return FD;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    //fluorine 
    //G4int Z = 50, A = 100;
    //G4double ionCharge   = Z*eplus;
    //G4double excitEnergy = 0.*keV;
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* positron
                    = particleTable->FindParticle("e+");
    //G4ParticleDefinition* ion
    //   = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(positron);//positron);
    //fParticleGun->SetParticleCharge(+Z*eplus);
  }

  // randomized direction  
  G4double theta = 0, cos_phi = 0;
  theta = 2*CLHEP::pi*G4UniformRand();
  cos_phi = 2*G4UniformRand()-1.0;
  G4double px = std::cos(theta)*pow(1-pow(cos_phi,2),0.5);
  G4double py = std::sin(theta)*pow(1-pow(cos_phi,2),0.5);
  G4double pz = cos_phi;

  // randomized plaque
  // random int generator experiment
  // construct a trivial random generator engine from a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_int_distribution<int> distrib(0, 5);
  G4double plaque_nb = distrib(generator);
  //std::cout << plaque_nb << std::endl;

  // randomize position: uniform random depth, circular xy position distribution
  theta = 2*CLHEP::pi*G4UniformRand();
  G4double r_max = 2*mm;
  G4double r_1 = G4UniformRand()*r_max;
  G4double r = r_1*(1-G4UniformRand());
  G4double x = r*std::cos(theta);
  G4double y = r*std::sin(theta);

  G4double first_pos = 38.7*mm, plaque_sep = 11.6*mm, detector_Z = 0.5*mm;
  G4double z = first_pos-plaque_nb*plaque_sep + (0.5-G4UniformRand())*detector_Z;


  // energy distribution part
  G4double end_point_energy = 3.72*MeV; // beta+ decay of 100Sn

  G4double a_x;
  G4int n_x = 1000, n_y = 1000;
  std::vector< G4double > arr;
  for(G4int i=0; i<n_x; i++){
    G4double i_f = i;
    G4double n_x_f = n_x;
    a_x = 1.3*i_f/n_x_f;
    G4double a_E = (a_x+0.511/3.72);
    G4double a_p = pow(pow(a_E,2)+pow(0.511/3.72,2),0.5);
    G4double top = FermiDistribution(a_x)*a_p*a_E*pow(1-a_x,2);
    //std::cout << a_x << "," << a_E << "," << a_p << "," << pow(1-a_x,2) << "," << FermiDistribution(a_x) << "," << top << std::endl;
    //std::cout << top << std::endl;
    for(G4int j=0; j<top*n_y; j++){
      arr.push_back(a_x*end_point_energy);
    }
  }

  std::uniform_int_distribution<int> lucky(0, arr.size());  // which plaque is getting the beta

  G4double rdenergy = arr.at(lucky(generator));
  //for(;;){
  //  rdenergy = 2*MeV*G4UniformRand();
  //  bool result = G4BetaFunction(0,2*MeV,rdenergy,1,49);
  //  if(true == result){
  //    break;
  //  }
  //}

  std::cout << plaque_nb << "," << rdenergy << "," << x << "," << y << "," << z << "," << px << "," << py << "," << pz << std::endl;

  fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));
  fParticleGun->SetParticleEnergy(rdenergy);
  //create vertex
  //
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

