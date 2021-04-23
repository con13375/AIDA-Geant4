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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
 : G4UserEventAction(), 
   fRunAction(runAction),
   fCollID_cryst(-1),
   fCollID_patient(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::addEdep(G4double Edep, G4int N_z, G4int N_y, G4int N_x, G4double t)
{

  TotalEnergyDepositX[128*N_z + N_x] += Edep;
  TotalEnergyDepositY[128*N_z + N_y] += Edep;
  TEDPixel[128*128*N_z + 128*N_x + N_y] += Edep; // set up like this, n_y is the residual i%128, then n_x is the residual i/128%128, and then n_z is i/128/128

  // these next lines save the min and max global time values for energy deposits in the two paradigms (X/Y or pixel)
  minAndMaxTimeX[128*N_z + N_x] = std::min(minAndMaxTimeX[128*N_z + N_x],t);
  minAndMaxTimeX[128*6 + 128*N_z + N_x] = std::max(minAndMaxTimeX[128*N_z + N_x],t);
  minAndMaxTimeY[128*N_z + N_x] = std::min(minAndMaxTimeY[128*N_z + N_x],t);
  minAndMaxTimeY[128*6 + 128*N_z + N_x] = std::max(minAndMaxTimeY[128*N_z + N_x],t);
  minAndMaxTimePixel[128*128*N_z + 128*N_x + N_y] = std::min(minAndMaxTimePixel[128*128*N_z + 128*N_x + N_y],t);
  minAndMaxTimePixel[128*128*6 + 128*128*N_z + 128*N_x + N_y] = std::max(minAndMaxTimePixel[128*128*6 + 128*128*N_z + 128*N_x + N_y],t);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*evt*/)
{
  for(G4int i=0; i<128*6; i++){TotalEnergyDepositX[i] = 0;}
  for(G4int i=0; i<128*6; i++){TotalEnergyDepositY[i] = 0;}
  for(G4int i=0; i<128*128*6; i++){TEDPixel[i] = 0;}
  for(G4int i=0; i<2*128*6; i++){if(i<128*6){minAndMaxTimeX[i] = 1000;} else{minAndMaxTimeX[i] = 0;}} // 1000 is chosen supposing it is bigger than the time values fed into the addEdep function
  for(G4int i=0; i<2*128*6; i++){if(i<128*6){minAndMaxTimeY[i] = 1000;} else{minAndMaxTimeY[i] = 0;}}
  for(G4int i=0; i<2*128*128*6; i++){if(i<128*128*6){minAndMaxTimePixel[i] = 1000;} else{minAndMaxTimePixel[i] = 0;}}
  //std::cout << "##" << "," << "massOfParticle" << "," << "chargeOfParticle" << "," << "energyDeposited(MeV)" << "," << "x(mm)" << "," << "y(mm)" << "," << "z(mm)" << "," << "t(ns)" << "," << "n_x" << "," << "n_y" << "," << "n_z" << std::endl;
  EventCount += 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* /*evt*/ )
{
  //std::cout << "Event" << "," << "x(0)/y(1)" << "," << "energyDep_x" << "," << "N_x" << "," << "n_plaque" << std::endl;// << "," << "t_min" << "," << "t_max" << std::endl;
  for(G4int i=0; i<128*6; i++){
    if(TotalEnergyDepositX[i] > 0){
      std::cout << "##" << "," << EventCount << "," << "0" << "," << TotalEnergyDepositX[i] << "," << i%128+1 << "," << i/128+1 << std::endl;// << "," << minAndMaxTimeX[i] << "," << minAndMaxTimeX[128*6+i] << std::endl;
    }
  }
  //std::cout << "x(0)/y(1)" << "," << "energyDep_Y" << "," << "N_y" << "," << "n_plaque" << "," << "t_min" << "," << "t_max" << std::endl;
  for(G4int i=0; i<128*6; i++){
    if(TotalEnergyDepositY[i] > 0){
      std::cout << "##" << "," << EventCount << "," << "1" << "," << TotalEnergyDepositY[i] << "," << i%128+1 << "," << i/128+1 << std::endl;// << "," << minAndMaxTimeY[i] << "," << minAndMaxTimeY[128*6+i] << std::endl;
    }
  }
  //std::cout << "~" << "," << "energyDep_pixel" << "," << "N_x" << "," << "N_y" << "," << "n_plaque" << "," << "t_min" << "," << "t_max" << std::endl;
  for(G4int i=0; i<128*128*6; i++){
    if(TEDPixel[i] > 0){
      std::cout << "###" << "," << EventCount << "," << TEDPixel[i] << "," << (i/128)%128+1 << "," << i%128+1 << "," << (i/128/128)+1 << std::endl;//<< "," << minAndMaxTimePixel[i] << "," << minAndMaxTimePixel[128*128*6+i] << std::endl;
    }
  }
}  


