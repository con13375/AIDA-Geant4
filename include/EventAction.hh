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
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RunAction;

/// Event action class
///
/// In EndOfEventAction() there is collected information event per event 
/// from Hits Collections, and accumulated statistic for 
/// RunAction::EndOfRunAction().

class EventAction : public G4UserEventAction
{
  public:
    EventAction(RunAction* runAction);
    virtual ~EventAction(); 

    virtual void  BeginOfEventAction(const G4Event*); //virtual
    virtual void  EndOfEventAction(const G4Event*); //virtual
    void addEdep(G4double Edep, G4int N_z, G4int N_y, G4int N_x, G4double t, G4double charge, G4double energy);

    G4int                       timestamp;
    G4int                       FIRST;
    G4double                    previous_time;
    G4double                    first_time;    

  private:
    RunAction*  fRunAction;
    G4double                  TED;
    G4double                  TotalBetaPercent;
    G4double                  BETA_energy;
    G4int                     primary_x;
    G4int                     primary_y;
    G4int                     primary_z;
    G4double                  TotalEnergyDepositX[128*6];
    G4double                  TotalEnergyDepositY[128*6];
    G4double                  TEDPixel[128*128*6];
    G4double                  betaPercent[128*128*6];
//    G4int                     chargePixel[128*128*6];
    G4double                  minAndMaxTimeX[128*6*2];
    G4double                  minAndMaxTimeY[128*6*2];
    G4double                  minAndMaxTimePixel[128*128*6*2];
    G4double                  Beta_total_w;
    G4int                  EventCount;
    G4int fCollID_cryst;
    G4int fCollID_patient;   
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
