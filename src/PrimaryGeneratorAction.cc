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

#include <math.h>       /* pow */ /* tgamma */
#include <chrono>
#include <string>
#include <random>       
// std::uniform_int_distribution<int> distrib(a, b)
#include <vector>
#include <complex>
#include "complex4.h"
//std::vector< int > arr;
//arr.push_back(1);

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

int discretize2(G4double x, G4double a, G4double b, G4int N){
  G4double slope = N/(b-a);
  if(x==b){
    return N;
  }
  else{
    return std::floor(slope*(x-a));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // default particle kinematic

//  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  fParticleGun->SetParticleDefinition(positron);
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

G4double FermiDistribution(G4int Z, G4double x, G4double EP){
  // particle energy by fermi distribution
  // x runs from 0 to 1 (where 1 represents energy end point)
  G4double x_E = x*EP/0.511+1; // E/mc2, where EP is the endpoint energy
  G4double x_p = pow(pow(x_E,2)-1,0.5); // pc/W

  G4double S;
  G4double A=1./137.0359895;
  G4double y;

//  std::complex<G4double> k;
  S = sqrt(1.- pow((A*Z),2.));
  y = A*Z*x_E/x_p;
  std::complex<G4double> k(S,y);

  G4double FD;
  if(x>0){
    FD = pow(x_p,(2.*(S-1)))*exp(M_PI*y)*pow((abs(cgamma(k))),2.);
  }
  else{
    FD = 0.;
  }

  return FD;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* positron = G4Positron::PositronDefinition();

//  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
//  G4int Z = 50, A = 100;
//  G4double ionCharge   = Z*eplus;
//  G4double excitEnergy = 0.*keV;
//  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
//  G4ParticleDefinition* geantino = particleTable->FindParticle("geantino");
  //G4ParticleDefinition* ion
  //   = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
  fParticleGun->SetParticleDefinition(positron);//ion);//(geantino);//
  //fParticleGun->SetParticleCharge(ionCharge);

  // randomized direction  
  G4double theta = 0, cos_phi = 0;
  theta = 2*CLHEP::pi*G4UniformRand();
  cos_phi = 2*G4UniformRand()-1.0;
  G4double px = std::cos(theta)*pow(1-pow(cos_phi,2),0.5);
  G4double py = std::sin(theta)*pow(1-pow(cos_phi,2),0.5);
  G4double pz = cos_phi;

  // this part chooses a plaque at random
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::uniform_int_distribution<int> distrib(0, 5);
  G4double plaque_nb = distrib(generator);
  //std::cout << plaque_nb << std::endl;

  // randomize position x,y,z
  	// circular xy position distribution with peak at center
  theta = 2*CLHEP::pi*G4UniformRand();
  G4double r_max = 40*mm; // from figure 4.3 in Oscar's thesis, we see the radius covers most of the detector
  G4double r_1 = G4UniformRand()*r_max; // this and the next line are a way to generate a random distribution with a peak in the center
  G4double r = r_1*(1-G4UniformRand());

  // for a gaussian distribution of implantation place (quadrilateral symmetry)
//  std::default_random_engine generator1;
  std::normal_distribution<G4double> distribution(0.0,10.0);
  G4double num1 = distribution(generator), num2 = distribution(generator);
  //std::cout << num1 << "," << num2 << std::endl;
  G4double number_x = std::min(40.0,std::max(-40.0,num1));
  G4double number_y = std::min(40.0,std::max(-40.0,num2));

  G4double x = number_x*mm;//60*mm*(1-2*G4UniformRand());//r*std::cos(theta);//
  G4double y = number_y*mm;//60*mm*(1-2*G4UniformRand());//r*std::sin(theta);//
  // uniform distribution in z within depth of detector from the chosen plaque
  G4double first_pos = 38.7*mm, plaque_sep = 11.6*mm, detector_Z = 0.5*mm;
  G4double z = first_pos-plaque_nb*plaque_sep + (0.5-G4UniformRand())*detector_Z;

 	// This part transforms position to strip number (x,y), number of plaque (z)
  G4double detector_XY = 71.63/2.0*mm;
  	//G4double first_pos = 38.7*mm, plaque_sep = 11.6*mm;
  G4double max_Z = first_pos+0.5*plaque_sep;
  G4double min_Z = first_pos-5.5*plaque_sep;
  G4int nx = discretize2(x, -detector_XY, detector_XY, 128);
  G4int ny = discretize2(y, -detector_XY, detector_XY, 128);
  G4int nz = discretize2(-z, -max_Z, -min_Z, 6);

  // energy distribution part
//  G4double rdenergy = 0; //setting the ion at rest now

  G4double EP = 4.74;
  G4double end_point_energy = EP*MeV; // beta+ decay of 100Sn, but also 3.72 meV??

  G4int Z = 99; // daughter nucleus charge, positive because it's beta+
  G4double a_x;

// I created this 'for' loop as a brute force way to generate numbers from a probability distribution
  const G4int n_x = 1000; // number of bins in the spectrum
  G4int n_box = 100; // the number of boxes per bin in average for the spectrum
  G4int n_y = n_box*n_x; // the total number of boxes for the spectrum
  std::vector< G4double > arr; // this vector will be up to size n_y, will contain all the "energy spectrum boxes" from which a random box will be chosen

  G4double i_f; // not sure if this is necessary but this is a float version
  G4double n_x_f; // float version
  G4double a_E;
  G4double a_p;
  G4double top[n_x];
  G4double accum;
  G4double value;
  accum = 0;
  G4double height;
  for(G4int i=0; i<n_x; i++){
    i_f = i;
    n_x_f = n_x;
    a_x = 1*i_f/n_x_f; // was thinking of putting 1.25 to see what happens if i go overboard but no
    a_E = (a_x+0.511/EP);
    a_p = pow(pow(a_E,2)-pow(0.511/EP,2),0.5);
    value = FermiDistribution(Z, a_x, EP)*a_p*a_E*pow(1-a_x,2);
    top[i] = value;
    accum += value;
    //std::cout << a_x << "," << a_E << "," << a_p << "," << pow(1-a_x,2) << "," << FermiDistribution(a_x) << "," << top << std::endl;
    //std::cout << top << std::endl;
  }
  for(G4int k=0; k<n_x; k++){
    i_f = k;
    n_x_f = n_x;
    a_x = 1*i_f/n_x_f;
    height = top[k]*n_y/accum; // dividing by accum normalizes
    //std::cout << k << "," << height << std::endl;
    for(G4int j=0; abs(j-height) < n_box/100 or j<height; j++){ // because 'height' is not an integer, we could not have known the size of arr from the start, and the condition z allows a window around j == height
      arr.push_back(a_x*end_point_energy); // this makes a vector where the higher the value in array "top", the more frequent an element within it that has the corresponding energy a_x
    }
  }

  std::uniform_int_distribution<int> lucky(0, arr.size());  // which plaque is getting the beta

  G4double rdenergy = arr.at(lucky(generator));

  //if(GunCount == 0.0){x = 0*mm, y = 00*mm;}
  //if(GunCount == 1.0){x = -42.0883*mm, y = -40.4312*mm;}
  //if(GunCount == 2.0){x = -8.3195*mm, y = 39.9027*mm;}
  //if(GunCount == 3.0){x = -34.5304*mm, y = -48.0646*mm;}
  //if(GunCount == 4.0){x = -47.2406*mm, y = 34.5269*mm;}
  //if(GunCount == 5.0){x = (47.5-1.3*1.5)*mm, y = (47.5-1.3*1.5)*mm;}
  //x = 0, y = 0, z =  0*mm;
  
  GunCount += 1;
  std::cout << "Event# " << GunCount << std::endl;
  //std::cout << "#" << "," <<  "Event" << "," << "plaque_nb" << "," << "energy(MeV)" << "," << "x(mm)" << "," << "y(mm)" << "," << "z(mm)" << "," << "px" << "," << "py" << "," << "pz" << "," << "n_x" << "," << "n_y" << "," << "n_z" << std::endl;
  std::cout << "#" << "," << GunCount << "," << rdenergy << "," << nx+1 << "," << ny+1 << "," << nz+1 << std::endl;
  
  // Setting particle gun
  fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));//(x,y,z));//
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px,py,pz));//(1,0,0));//
  fParticleGun->SetParticleEnergy(rdenergy);//(1*GeV);//
  //create vertex
  //
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

