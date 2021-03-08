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
/// \file B3DetectorConstruction.cc
/// \brief Implementation of the B3DetectorConstruction class

#include "B3DetectorConstruction.hh"

#include "G4SubtractionSolid.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::B3DetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B3DetectorConstruction::~B3DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();
  
  G4bool isotopes = false;
  
  G4Element*  O = man->FindOrBuildElement("O" , isotopes); 
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);  
  
  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O , 5);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B3DetectorConstruction::Construct()
{  

  // Gamma detector Parameters
  //
  G4double cryst_dX = 0.6*cm, cryst_dY = 0.6*cm, cryst_dZ = 0.3*cm;
  G4int nb_cryst = 32;
  G4int nb_rings = 3;
  //
  G4double dPhi = twopi/nb_cryst, half_dPhi = 0.5*dPhi;
  G4double cosdPhi = std::cos(half_dPhi);
  G4double tandPhi = std::tan(half_dPhi);
  // 
  G4double ring_R1 = 0.5*cryst_dY/tandPhi;
  G4double ring_R2 = (ring_R1+cryst_dZ)/cosdPhi;
  //
  G4double detector_dZ = nb_rings*cryst_dX;
  //
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2SiO5");
        
  //     
  // World
  //
  G4double world_sizeXY = 10*cm;
  G4double world_sizeZ  = 60*cm;
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       1*world_sizeXY, 1*world_sizeXY, 1*world_sizeZ); //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        default_mat,         //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      fCheckOverlaps);       // checking overlaps 
                 
  // Parameters
  G4double plaqueVol_XY = 6.3*cm, plaqueVol_Z = 0.65*cm;
  G4double DDSD_XY = 4*cm, DDSD_Z = 0.5*mm;
  G4double Plastic_XY = 6*cm, Plastic_Z = 1*mm;
  G4double connector_dx = 2.125*cm, connector_dy = 0.6*cm, connector_dz = 0.5*cm;
  G4double kapton_Y = 0.5*mm;

  G4int nb_plaques = 6;
  G4double separation = 10*cm;
  G4double tube_in = 2.3*mm, tube_out = 3*mm, tube_Z = world_sizeZ;

  G4int n_tubes = 4;

  // materials

// define Elements
 
  G4double a = 1.01*g/mole;
  G4Element* elH  = new G4Element("Hydrogen","H" , 1., a); // name, symbol, z, molar density

  a = 12.01*g/mole;
  G4Element* elC = new G4Element("Carbon", "C", 6., a);

  a = 14.01*g/mole;
  G4Element* elN  = new G4Element("Nitrogen","N" , 7., a);

  a = 16.00*g/mole;
  G4Element* elO  = new G4Element("Oxygen"  ,"O" , 8., a);

  //G4NistManager* nist = G4NistManager::Instance();
  //G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* DDSD_mat = nist->FindOrBuildMaterial("G4_Si");
  
  G4double density = 1.06*g/cm3;
  G4int nel = 3;
  G4Material* Noryl = new G4Material("Noryl", density, nel);
  Noryl->AddMaterial(nist->FindOrBuildMaterial("G4_C"), 47.06*perCent);
  Noryl->AddMaterial(nist->FindOrBuildMaterial("G4_H"), 47.06*perCent);
  Noryl->AddMaterial(nist->FindOrBuildMaterial("G4_O"), 5.88*perCent);
  G4Material* mat_tubes = nist->FindOrBuildMaterial("G4_Ti");

  G4double kapton_density = 1.413*g/cm3;
  G4int kapton_nel = 4;
  G4Material* Kapton = new G4Material("Kapton", kapton_density, kapton_nel);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);

  // Titanium cilinders

  G4Tubs* tubes = new G4Tubs("tubes", tube_in, tube_out, tube_Z, 0, twopi);

  G4LogicalVolume* logic_tube = new G4LogicalVolume(tubes, mat_tubes, "tubes");

  G4double tube_pos_XY = Plastic_XY-1.5*tube_out, tube_pos_Z = tube_Z*0;
  for (G4int itubes = 0; itubes < n_tubes; itubes++){
    G4ThreeVector y = G4ThreeVector(std::cos(itubes*pi/2),std::sin(itubes*pi/2),0.);
    G4ThreeVector x = G4ThreeVector(-std::sin(itubes*pi/2),std::cos(itubes*pi/2),0.);
    G4ThreeVector z = G4ThreeVector(0,0,tube_pos_Z);
    G4ThreeVector position = (tube_pos_XY)*y+(tube_pos_XY)*x+z;
    new G4PVPlacement(0,position, logic_tube, "tube", logicWorld, false, itubes, fCheckOverlaps);
    }

  // plaques

  // solids
  G4Box* plaqueVol = new G4Box("plaqueVol", plaqueVol_XY, plaqueVol_XY, plaqueVol_Z);
  G4Box* DDSD = new G4Box("DDSD",DDSD_XY,DDSD_XY,DDSD_Z);
  G4Box* hole = new G4Box("hole",DDSD_XY,DDSD_XY,Plastic_Z);

  G4double mainPlastic_XY = Plastic_XY-3*tube_out, miniPlastic_XY = 1.5*tube_out;
  G4Box* Plastic = new G4Box("Plastic", mainPlastic_XY, mainPlastic_XY, Plastic_Z);
  G4SubtractionSolid* plastic = new G4SubtractionSolid("plastic", Plastic, hole);
  G4Box* miniPlastic = new G4Box("miniPlastic", miniPlastic_XY, miniPlastic_XY,Plastic_Z);
  G4SubtractionSolid* miniplastic = new G4SubtractionSolid("miniplastic", miniPlastic, tubes);
  G4Box* lateralplastic = new G4Box("miniPlastic", mainPlastic_XY, miniPlastic_XY, Plastic_Z);
  G4Box* connector = new G4Box("connector", connector_dx, connector_dy, connector_dz);

  // logics
  G4LogicalVolume* logic_plaqueVol = new G4LogicalVolume(plaqueVol, default_mat, "plaqueVol");
  G4LogicalVolume* logic_DDSD = new G4LogicalVolume(DDSD, DDSD_mat, "DDSD");
  G4LogicalVolume* logic_plastic = new G4LogicalVolume(plastic, Noryl, "plastic");
  G4LogicalVolume* logic_miniplastic = new G4LogicalVolume(miniplastic, Noryl, "plastic");
  G4LogicalVolume* logic_lateralplastic = new G4LogicalVolume(lateralplastic, Noryl, "plastic");
  G4LogicalVolume* logic_connector = new G4LogicalVolume(connector, Noryl, "connector");

  // Assembling plaque
  G4double AIDA_nose_Z = 50*cm;
  G4double Z = Plastic_Z-plaqueVol_Z;
  G4double Z2 = connector_dz+2*Plastic_Z-plaqueVol_Z;
  for (G4int iplaque = 0; iplaque < nb_plaques ; iplaque++) {
    
    G4double dZ = AIDA_nose_Z-(1+iplaque)*separation;

    std::cout << dZ << std::endl;

    G4ThreeVector plaque_center = G4ThreeVector(0,0,dZ);
    new G4PVPlacement(0,G4ThreeVector(0,0,Z),logic_DDSD,"DDSD",
			logic_plaqueVol,false,iplaque,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0,0,Z),logic_plastic,"Plastic",
			logic_plaqueVol,false,iplaque,fCheckOverlaps);

    // Assembling the four connectors and kapton strips
    for (G4int iconn = 0; iconn < 4 ; iconn++) {
      G4double kapton_Z = 0.5*(world_sizeZ+dZ);
      std::cout<<kapton_Z<<std::endl;
      G4Box* strip = new G4Box("kapton", connector_dx, kapton_Y, kapton_Z);
      G4LogicalVolume* logic_strip = new G4LogicalVolume(strip, Kapton, "kapton");

      G4RotationMatrix rotm = G4RotationMatrix();
      rotm.rotateZ(-iconn*pi/2);

      G4ThreeVector y = G4ThreeVector(std::cos(iconn*pi/2),std::sin(-iconn*pi/2),0.);
      G4ThreeVector x = G4ThreeVector(std::sin(iconn*pi/2),std::cos(iconn*pi/2),0.);
      G4ThreeVector z = G4ThreeVector(0,0,Z2);
      G4ThreeVector z_k = G4ThreeVector(0,0,-world_sizeZ+kapton_Z); //kapton Z position wrt world

      G4ThreeVector position = (0.5*DDSD_XY)*y+(Plastic_XY*1-connector_dy*1)*x+z;
      G4ThreeVector position_kapton = (0.5*DDSD_XY)*y+(Plastic_XY*1+iplaque*kapton_Y)*x+z_k;
      std::cout<<position_kapton<<std::endl;
      G4Transform3D transform = G4Transform3D(rotm,position);
      G4Transform3D transform_kapton = G4Transform3D(rotm,position_kapton);
      new G4PVPlacement(transform, logic_connector, "connector", logic_plaqueVol, false,
			4*iplaque+iconn, fCheckOverlaps);
      new G4PVPlacement(transform_kapton, logic_strip, "kapton", logicWorld, false,
			4*iplaque+iconn, fCheckOverlaps);
    }
  // Assembling the plastic with holes
    for (G4int imini = 0; imini < 4 ; imini++) {
      G4RotationMatrix rotm = G4RotationMatrix();
      rotm.rotateZ(-imini*pi/2);

      G4ThreeVector y = G4ThreeVector(std::cos(imini*pi/2),std::sin(-imini*pi/2),0.);
      G4ThreeVector x = G4ThreeVector(std::sin(imini*pi/2),std::cos(imini*pi/2),0.);
      G4ThreeVector z = G4ThreeVector(0,0,Z);
      G4ThreeVector position_mini = (Plastic_XY-1.5*tube_out)*y+(Plastic_XY-1.5*tube_out)*x+z;
      G4ThreeVector position_lateral = (Plastic_XY-1.5*tube_out)*x+z;

      G4Transform3D transform = G4Transform3D(rotm,position_lateral);
      new G4PVPlacement(transform, logic_lateralplastic, "lateralplastic", logic_plaqueVol, false,
			4*iplaque+imini, fCheckOverlaps);
      new G4PVPlacement(0, position_mini, logic_miniplastic, "miniplastic", logic_plaqueVol, false,
			4*iplaque+imini, fCheckOverlaps);
    }  
    new G4PVPlacement(0,G4ThreeVector(0,0,dZ), logic_plaqueVol, "plaque",
			logicWorld, false, iplaque, fCheckOverlaps);
  }

  //
  // ring
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring", ring_R1, ring_R2, 0.5*cryst_dX, 0., twopi);
      
  G4LogicalVolume* logicRing =                         
    new G4LogicalVolume(solidRing,           //its solid
                        default_mat,         //its material
                        "Ring");             //its name
                    
  //     
  // define crystal
  //
  G4double gap = 0.5*mm;        //a gap for wrapping
  G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  G4Box* solidCryst = new G4Box("crystal", dX/2, dY/2, cryst_dZ/2);
                     
  G4LogicalVolume* logicCryst = 
    new G4LogicalVolume(solidCryst,          //its solid
                        cryst_mat,           //its material
                        "CrystalLV");        //its name
               
  // place crystals within a ring 
  //
  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = icrys*dPhi;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg); 
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);     
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);
                                    
    new G4PVPlacement(transform,             //rotation,position
                      logicCryst,            //its logical volume
                      "crystal",             //its name
                      logicRing,             //its mother  volume
                      false,                 //no boolean operation
                      icrys,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }
                                                      
  //
  // full detector
  //
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2, 0.5*detector_dZ, 0., twopi);
      
  G4LogicalVolume* logicDetector =                         
    new G4LogicalVolume(solidDetector,       //its solid
                        default_mat,         //its material
                        "Detector");         //its name
                                 
  // 
  // place rings within detector 
  //
  G4double OG = -0.5*(detector_dZ + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(0,0,OG), //position
                      logicRing,             //its logical volume
                      "ring",                //its name
                      logicDetector,         //its mother  volume
                      false,                 //no boolean operation
                      iring,                 //copy number
                      fCheckOverlaps);       // checking overlaps 
  }
                       
  //
  // place detector in world
  //                    
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-100*cm),         //at (0,0,0)
                    logicDetector,           //its logical volume
                    "Detector",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps 
                 
  //
  // patient
  //
  G4double patient_radius = 0.8*cm;
  G4double patient_dZ = 0.10*cm;  
  G4Material* patient_mat = nist->FindOrBuildMaterial("G4_BRAIN_ICRP");
    
  G4Tubs* solidPatient =
    new G4Tubs("Patient", 0., patient_radius, 0.5*patient_dZ, 0., twopi);
      
  G4LogicalVolume* logicPatient =                         
    new G4LogicalVolume(solidPatient,        //its solid
                        patient_mat,         //its material
                        "PatientLV");        //its name
               
  //
  // place patient in world
  //                    
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-100*cm),         //at (0,0,0)
                    logicPatient,            //its logical volume
                    "Patient",               //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);         // checking overlaps 
                                          
  // Visualization attributes
  //
  logicRing->SetVisAttributes (G4VisAttributes::GetInvisible());
  logicDetector->SetVisAttributes (G4VisAttributes::GetInvisible());    

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl; 

  //always return the physical World
  //
  return physWorld;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B3DetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  
  // declare crystal as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4SDManager::GetSDMpointer()->AddNewDetector(cryst);
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);
  
  // declare patient as a MultiFunctionalDetector scorer
  //  
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  G4SDManager::GetSDMpointer()->AddNewDetector(patient);
  G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(primitiv2);
  SetSensitiveDetector("PatientLV",patient);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
