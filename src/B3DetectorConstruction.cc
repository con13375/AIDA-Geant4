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

#include "G4UnionSolid.hh"
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
  G4double world_sizeZ  = 31*cm;
  
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

  // plaque stuff
  //G4double plaqueVol_XY = 47.5*mm, plaqueVol_Z = 0.8*mm; // empty space mother of plaque, i regret doing this
  G4double DDSD_XY = 4*cm, DDSD_Z = 0.5*mm; //including inactive area
  G4double detector_XY = 3.815*cm; // just active area
  G4double Plastic_XY = 4.75*cm, Plastic_Z = 0.8*mm;

  G4double bit_size = 1.27*mm;
  G4double connector_dx = (bit_size*34+1.55)/2*mm, connector_dy = 2.36*mm, connector_dz = 2.39*mm;
  G4double chole1_dx = (bit_size*34)/2*mm, chole1_dy = connector_dy, chole1_dz = connector_dz-bit_size/2;
  G4double chole2_dx = (bit_size*26)/2*mm, chole2_dy = connector_dy, chole2_dz = connector_dz;
  G4double c_bits_dx = 0.47*bit_size, c_bits_dy = (1.105+1.255)*mm, c_bits_dz = chole1_dz;

  G4double c_copper1_XZ = 0.3*mm;
  G4double c_copper1_Y = connector_dy+3*c_copper1_XZ; //upper
  G4double c_copper2_XZ = c_copper1_XZ;
  G4double c_copper2_Y = connector_dy+1*c_copper1_XZ;  //lower
  G4double c_copper3_XY = c_copper1_XZ;
  G4double c_copper3_Z = 0.5*(connector_dz+bit_size/2+c_copper1_XZ);
  G4double c_copper4_XY = c_copper1_XZ;
  G4double c_copper4_Z = 0.5*(connector_dz-bit_size/2+c_copper1_XZ);
  
  G4int nb_plaques = 6;
  G4double separation = 10*mm;
  G4double AIDA_nose_Z = 5*cm;

  // kapton stuff
  G4double kapton_Y = 0.1*mm;
  G4double cu_X = 0.15*mm, cu_Y = 0.009*mm, cu_sep = 0.635*mm;
  G4int n_cu = 68; // this should be 68, but i put it lower while building to ease loading time

  // tubes, separators and bolts
  G4double tube_in = 0.825*mm, tube_out = 1.5*mm, tube_Z = 0.9*world_sizeZ;
  G4double separator_R1 = tube_out, separator_R2 = 2*mm, separator_R3 = 2.75*mm, separator_Z = 4.15*mm;
  G4double bolt_R1 = tube_out, bolt_R2 = 2*mm, bolt_R3 = 3.5*mm, bolt_Z1 = 2.65*mm, bolt_Z2 = 3.5*mm;
  G4double nut_R = 0.8*mm, nut_Z = 0.5*(bolt_R3-bolt_R1);

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

  a = 28.0855*g/mole;
  G4Element* elSi  = new G4Element("Silicon"  ,"Si" , 14., a);

//  a = 65.38*g/mole;
//  G4Element* elZn  = new G4Element("Zinc"  ,"Zn" , 12., a);
//
//  a = 118.710*g/mole;
//  G4Element* elSn  = new G4Element("Tin"  ,"Sn" , 50., a);
//
//  a = 55.845*g/mole;
//  G4Element* elFe  = new G4Element("Iron"  ,"Fe" , 26., a);
//
//  a = 30.973*g/mole;
//  G4Element* elP  = new G4Element("Phosphorus"  ,"P" , 15., a);
//
//  a = 65.546*g/mole;
//  G4Element* elCu  = new G4Element("Copper"  ,"Cu" , 29., a);

  //G4NistManager* nist = G4NistManager::Instance();
  //G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* DDSD_mat = nist->FindOrBuildMaterial("G4_Si");
  G4Material* Cu = nist->FindOrBuildMaterial("G4_Cu");
  
  // for bolts
  G4double density = 1.06*g/cm3;
  G4int nel = 3;
  G4Material* Noryl = new G4Material("Noryl", density, nel);
  Noryl->AddMaterial(nist->FindOrBuildMaterial("G4_C"), 47.06*perCent);
  Noryl->AddMaterial(nist->FindOrBuildMaterial("G4_H"), 47.06*perCent);
  Noryl->AddMaterial(nist->FindOrBuildMaterial("G4_O"), 5.88*perCent);

  // for connector (LCP), i'm just setting the noryl percentages because I could find them
  density = 1.75*g/cm3;
  nel = 3;
  G4Material* LCP = new G4Material("LCP", density, nel);
  LCP->AddMaterial(nist->FindOrBuildMaterial("G4_C"), 47.06*perCent);
  LCP->AddMaterial(nist->FindOrBuildMaterial("G4_H"), 47.06*perCent);
  LCP->AddMaterial(nist->FindOrBuildMaterial("G4_O"), 5.88*perCent);
  
  //Epoxy (for FR4 )
  density = 1.2*g/cm3;
  nel = 2;
  G4Material* Epoxy = new G4Material("Epoxy" , density, nel);
  Epoxy->AddElement(elH, 2);
  Epoxy->AddElement(elC, 2); 

  //SiO2 (Quarz)
  density= 2.200*g/cm3;
  nel = 2;
  G4Material* SiO2 =  new G4Material("SiO2",density, nel);
  SiO2->AddElement(elSi, 1); 
  SiO2->AddElement(elO , 2); 
 
  //FR4 (Glass + Epoxy)  
  density = 1.86*g/cm3; 
  nel = 2;
  G4Material* FR4 = new G4Material("FR4"  , density, nel); 
  FR4->AddMaterial(Epoxy, 47.2*perCent); 
  FR4->AddMaterial(SiO2, 52.8*perCent);

  //Phosphor Bronze
  density = 8.86*g/cm3; //unverified
  nel = 5;
  G4Material* PBronze = new G4Material("FR4"  , density, nel); 
  PBronze->AddMaterial(nist->FindOrBuildMaterial("G4_Zn"), 9.9*perCent);
  PBronze->AddMaterial(nist->FindOrBuildMaterial("G4_Sn"), 2.2*perCent);
  PBronze->AddMaterial(nist->FindOrBuildMaterial("G4_Fe"), 1.9*perCent);
  PBronze->AddMaterial(nist->FindOrBuildMaterial("G4_P"), 0.03*perCent);
  PBronze->AddMaterial(nist->FindOrBuildMaterial("G4_Cu"), 85.97*perCent);

  G4Material* mat_tubes = nist->FindOrBuildMaterial("G4_Ti");

  G4double kapton_density = 1.413*g/cm3;
  G4int kapton_nel = 4;
  G4Material* Kapton = new G4Material("Kapton", kapton_density, kapton_nel);
  Kapton->AddElement(elO,5);
  Kapton->AddElement(elC,22);
  Kapton->AddElement(elN,2);
  Kapton->AddElement(elH,10);

  // Titanium cilinders, separators, and bolts

  G4Tubs* tubes = new G4Tubs("tubes", tube_in, tube_out, tube_Z, 0, twopi);

  G4Tubs* Sep1 = new G4Tubs("Sep1", separator_R1, separator_R2, separation*0.5, 0, twopi);
  G4Tubs* Sep2 = new G4Tubs("Sep1", separator_R2, separator_R3, separator_Z, 0, twopi);
  G4ThreeVector zTrans(0,0,-0.85);
  G4UnionSolid* separator = new G4UnionSolid("separator", Sep1, Sep2, 0, zTrans);

  G4Tubs* Bolt1 = new G4Tubs("bolt1", bolt_R1, bolt_R2, bolt_Z2, 0, twopi);
  G4Tubs* Bolt2 = new G4Tubs("bolt2", bolt_R2, bolt_R3, bolt_Z1, 0, twopi);
  G4Tubs* nut = new G4Tubs("nut", 0, nut_R, nut_Z, 0, twopi);
  G4UnionSolid* bolt1 = new G4UnionSolid("bolt1", Bolt1, Bolt2, 0, zTrans);
  G4RotationMatrix* yRot = new G4RotationMatrix;
  yRot->rotateY(twopi/4.*rad);
  G4ThreeVector xzTrans(2*nut_Z-bolt_R3,0,-0.85);
  G4SubtractionSolid* bolt = new G4SubtractionSolid("bolt",bolt1,nut,yRot,xzTrans);

  G4LogicalVolume* logic_tube = new G4LogicalVolume(tubes, mat_tubes, "tubes");
  G4LogicalVolume* logic_separator = new G4LogicalVolume(separator, Noryl, "separator");
  G4LogicalVolume* logic_bolt = new G4LogicalVolume(bolt, Noryl, "bolt");
  G4LogicalVolume* logic_nut = new G4LogicalVolume(nut, Cu, "bolt");

  G4double tube_pos_XY = Plastic_XY-2.5*tube_out, tube_pos_Z = tube_Z*0;
  for (G4int itubes = 0; itubes < n_tubes; itubes++){
    G4ThreeVector y = G4ThreeVector(std::cos(itubes*pi/2),std::sin(itubes*pi/2),0.);
    G4ThreeVector x = G4ThreeVector(-std::sin(itubes*pi/2),std::cos(itubes*pi/2),0.);
    G4ThreeVector z = G4ThreeVector(0,0,tube_pos_Z);
    G4ThreeVector position = (tube_pos_XY)*y+(tube_pos_XY)*x+z;
    //std::cout << position << std::endl;
    //std::cout << itubes << std::endl;

    // separators
    for (G4int isep = 0; isep < nb_plaques + 1; isep++){ //nb_plaques has a +1 to put the other endstops
      if(isep == 0 or isep == nb_plaques){ // if end or start, put end stops
        G4int inv = std::cos(pi*isep/nb_plaques); // if at end, this is -1, and rotates and translates accordingly
        G4int oneOrzero = std::cos(pi*isep/nb_plaques/2); // if at end, this is 0, and places at correct plaque
        G4ThreeVector z_i = G4ThreeVector(0,0, AIDA_nose_Z+inv*Plastic_Z+inv*bolt_Z2-(oneOrzero+isep)*(separation+2*Plastic_Z));
        G4RotationMatrix* boltRot = new G4RotationMatrix;
        boltRot -> rotateY(pi*isep/nb_plaques);
        G4ThreeVector position_bolt = position+z_i;
        //std::cout << itubes%2 << std::endl;
        //std::cout << position_bolt << std::endl;
        new G4PVPlacement(yRot,position_bolt+inv*G4ThreeVector(nut_Z-bolt_R3,0,-0.85), logic_nut, "nut",
	  logicWorld, false, 2*itubes+isep, fCheckOverlaps);
        new G4PVPlacement(boltRot,position_bolt, logic_bolt, "bolt", logicWorld, false, 2*itubes+isep, fCheckOverlaps);
      }
      else{ //otherwise, put separators in between
        G4ThreeVector z_i = G4ThreeVector(0,0, AIDA_nose_Z-(0.5+isep)*(separation+2*Plastic_Z));
        G4ThreeVector position_sep = position+z_i;
        //std::cout << itubes%2+10 << std::endl;
        //std::cout << position_sep << std::endl;
        new G4PVPlacement(0,position_sep, logic_separator, "separator", logicWorld, false,2*itubes+isep, fCheckOverlaps);
      }
    }
    new G4PVPlacement(0,position, logic_tube, "tube", logicWorld, false, itubes, fCheckOverlaps);
  }
  // plaques

  // solids
  //G4Box* plaqueVol = new G4Box("plaqueVol", plaqueVol_XY, plaqueVol_XY, plaqueVol_Z);
  G4Box* DDSD = new G4Box("DDSD",DDSD_XY,DDSD_XY,DDSD_Z);
  G4Box* hole1 = new G4Box("hole",DDSD_XY,DDSD_XY,DDSD_Z); // this hole is where the silicon is placed
  G4Box* hole2 = new G4Box("hole",detector_XY,detector_XY,2*Plastic_Z); // this deep hole is where the active silicon is

  G4Box* Plastic = new G4Box("Plastic", Plastic_XY, Plastic_XY, Plastic_Z);

  G4ThreeVector detector_spot(0,0,Plastic_Z-DDSD_Z); // placed just at the surface of the plastic
  G4SubtractionSolid* plastic0 = new G4SubtractionSolid("plastic0", Plastic, hole1, 0, detector_spot); //silicon detector spot
  G4SubtractionSolid* plastic1 = new G4SubtractionSolid("plastic1", plastic0, hole2); //silicon detector hole
  // perforing four holes in the borders of the plastic for the titanium tube to go through
  G4Tubs* solid_tubes = new G4Tubs("tubes", 0, 1.05*tube_out, 1.05*Plastic_Z, 0, twopi);
  G4double th = tube_pos_XY;
  G4ThreeVector tubeHole1(th,th,0);
  G4ThreeVector tubeHole2(th,-th,0);
  G4ThreeVector tubeHole3(-th,th,0);
  G4ThreeVector tubeHole4(-th,-th,0);
  G4SubtractionSolid* plastic2 = new G4SubtractionSolid("plastic2", plastic1, solid_tubes,0,tubeHole1);
  G4SubtractionSolid* plastic3 = new G4SubtractionSolid("plastic3", plastic2, solid_tubes,0,tubeHole2);
  G4SubtractionSolid* plastic4 = new G4SubtractionSolid("plastic4", plastic3, solid_tubes,0,tubeHole3);
  G4SubtractionSolid* plastic = new G4SubtractionSolid("plastic", plastic4, solid_tubes,0,tubeHole4);

  G4Box* connector1 = new G4Box("connector1", connector_dx, connector_dy, connector_dz);
  G4Box* chole_1 = new G4Box("c_hole", chole1_dx, chole1_dy, chole1_dz);
  G4Box* chole_2 = new G4Box("c_hole2", chole2_dx, chole2_dy, chole2_dz);
  G4SubtractionSolid* connector2 = new G4SubtractionSolid("connector", connector1, chole_1);
  G4SubtractionSolid* connector = new G4SubtractionSolid("connector", connector2, chole_2);

// Naming these copper was a mistake because they're actually phosphor bronze
  G4Box* c_copper1 = new G4Box("c_copper1", c_copper1_XZ, c_copper1_Y, c_copper1_XZ);
  G4Box* c_copper2 = new G4Box("c_copper2", c_copper2_XZ, c_copper2_Y, c_copper2_XZ);
  G4Box* c_copper3 = new G4Box("c_copper3", c_copper3_XY, c_copper3_XY, c_copper3_Z);
  G4Box* c_copper4 = new G4Box("c_copper4", c_copper4_XY, c_copper4_XY, c_copper4_Z);
  G4Box* connector_bit1 = new G4Box("connector_bit1", c_bits_dx, c_bits_dy, c_bits_dz);
  G4ThreeVector bit1z(0,0,bit_size/2);
  G4ThreeVector bit2z(0,0,-bit_size/2);
  G4RotationMatrix* xRot = new G4RotationMatrix;
  xRot->rotateX(0*twopi/4.*rad);
  G4SubtractionSolid* connector_bit2 = new G4SubtractionSolid("connector_bit2", connector_bit1, c_copper1,xRot,bit1z);
  G4SubtractionSolid* connector_bit = new G4SubtractionSolid("connector_bit1", connector_bit2, c_copper2,xRot,bit2z);

  // logics
//  G4LogicalVolume* logic_plaqueVol = new G4LogicalVolume(plaqueVol, default_mat, "plaqueVol");
  G4LogicalVolume* logic_DDSD = new G4LogicalVolume(DDSD, DDSD_mat, "DDSD");
  G4LogicalVolume* logic_plastic = new G4LogicalVolume(plastic, FR4, "plastic");
//  G4LogicalVolume* logic_miniplastic = new G4LogicalVolume(miniplastic, FR4, "plastic");
//  G4LogicalVolume* logic_lateralplastic = new G4LogicalVolume(lateralplastic, FR4, "plastic");
  G4LogicalVolume* logic_connector = new G4LogicalVolume(connector, LCP, "connector");
  G4LogicalVolume* logic_bit = new G4LogicalVolume(connector_bit, LCP, "bit");
  G4LogicalVolume* logic_copper1 = new G4LogicalVolume(c_copper1, PBronze, "Cu1");
  G4LogicalVolume* logic_copper2 = new G4LogicalVolume(c_copper2, PBronze, "Cu2");
  G4LogicalVolume* logic_copper3 = new G4LogicalVolume(c_copper3, PBronze, "Cu3");
  G4LogicalVolume* logic_copper4 = new G4LogicalVolume(c_copper4, PBronze, "Cu4");

  // Assembling plaque
  G4double Z = 0;
  G4double Z2 = connector_dz+Plastic_Z;
  for (G4int iplaque = 0; iplaque < nb_plaques ; iplaque++) {
    
    G4double dZ = AIDA_nose_Z-(1+iplaque)*(separation+2*Plastic_Z);
    G4ThreeVector dZ_3V = G4ThreeVector(0,0,dZ); // I'm completely ditching the mother volume of plaque and defining everything wrt world
    //new G4PVPlacement(0,G4ThreeVector(0,0,dZ), logic_plaqueVol, "plaque",
	//		logicWorld, false, iplaque, fCheckOverlaps);
    //std::cout << dZ << std::endl;

    G4ThreeVector plaque_center = dZ_3V;
    std::cout << G4ThreeVector(0,0,Z)+dZ_3V << std::endl;
    new G4PVPlacement(0,G4ThreeVector(0,0,Z)+dZ_3V+detector_spot,logic_DDSD,"DDSD",
			logicWorld,false,iplaque,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0,0,Z)+dZ_3V,logic_plastic,"plastic",
			logicWorld,false,iplaque,fCheckOverlaps);

    // Assembling the four connectors and kapton strips with their copper things
    for (G4int iconn = 0; iconn < 4 ; iconn++) {
      // Kapton strip definitions (here because of changing Z lenght)
      G4double kapton_Z = 0.5*(world_sizeZ+dZ);
      G4Box* strip = new G4Box("kapton", connector_dx, kapton_Y, kapton_Z);
      G4LogicalVolume* logic_strip = new G4LogicalVolume(strip, Kapton, "kapton");
  
      G4Colour kapton_brown(0.6, 0.2, 0.2);
      G4VisAttributes* kapton_color = new G4VisAttributes(kapton_brown);
      logic_strip -> SetVisAttributes(kapton_color);

      // Cu strips in kapton
      G4Box* cu_strip = new G4Box("cu_strip", cu_X, cu_Y, kapton_Z);
      G4LogicalVolume* logic_cu_strip = new G4LogicalVolume(cu_strip, Cu, "cu_strip");

      G4Colour copper_brown(0.4, 0.15, 0.15);
      G4VisAttributes* copper_color = new G4VisAttributes(copper_brown);
      logic_cu_strip -> SetVisAttributes(copper_color);

      G4RotationMatrix rotm = G4RotationMatrix();
      rotm.rotateZ(-iconn*pi/2);

      G4ThreeVector y = G4ThreeVector(std::cos(iconn*pi/2),std::sin(-iconn*pi/2),0.);
      G4ThreeVector x = G4ThreeVector(std::sin(iconn*pi/2),std::cos(iconn*pi/2),0.);
      G4ThreeVector z = G4ThreeVector(0,0,1);
      G4ThreeVector z_k = G4ThreeVector(0,0,-world_sizeZ+kapton_Z+5*mm); //kapton Z position wrt world

      G4ThreeVector position = (0.4*DDSD_XY)*y+(Plastic_XY*1-connector_dy*1)*x+Z2*z+dZ_3V;
      G4ThreeVector position_kapton = (0.4*DDSD_XY)*y+(Plastic_XY*1+kapton_Y+2.1*iplaque*(cu_Y+kapton_Y))*x+z_k;
      //std::cout<<position_kapton<<std::endl;
      G4Transform3D transform = G4Transform3D(rotm,position);
      G4Transform3D transform_kapton = G4Transform3D(rotm,position_kapton);

      // Putting the Cu strips on kapton
      for (G4int icu = 0; icu < n_cu; icu++) {
        G4ThreeVector position_cu = position_kapton+x*(kapton_Y+cu_Y)+y*(-n_cu*0.5*cu_sep+icu*cu_sep);
        G4Transform3D transform_cu = G4Transform3D(rotm,position_cu);
        //std::cout<<position_cu<<std::endl;
        new G4PVPlacement(transform_cu, logic_cu_strip, "cu_strip", logicWorld, false, 4*n_cu*iplaque+n_cu*iconn+icu, fCheckOverlaps);
      }
      
      // Putting 34 connector bits
      G4int n_bits = n_cu/2;
      for (G4int ibit = 0; ibit < n_bits; ibit++) {
        G4double bit_dX = -n_bits*0.5*bit_size+ibit*bit_size;
        G4ThreeVector position_bit = bit_dX*y+position;
        G4Transform3D transform_bit = G4Transform3D(rotm,position_bit);
        new G4PVPlacement(transform_bit, logic_bit, "bit", logicWorld, false,
			  4*n_bits*iplaque+n_bits*iconn+ibit, fCheckOverlaps);

        G4ThreeVector position_cu1 = position_bit+(-3*c_copper1_XZ)*x+bit_size/2*z;
        G4Transform3D transform_cu1 = G4Transform3D(rotm,position_cu1);
        new G4PVPlacement(transform_cu1, logic_copper1, "cu1", logicWorld, false,
			  4*n_bits*iplaque+n_bits*iconn+ibit, fCheckOverlaps);

        G4ThreeVector position_cu2 = position_bit+(-1*c_copper1_XZ)*x-bit_size/2*z;
        G4Transform3D transform_cu2 = G4Transform3D(rotm,position_cu2);
        new G4PVPlacement(transform_cu2, logic_copper2, "cu2", logicWorld, false,
			  4*n_bits*iplaque+n_bits*iconn+ibit, fCheckOverlaps);

        G4ThreeVector position_cu3 = position_bit+(-connector_dy-7*c_copper1_XZ)*x+(-connector_dz+c_copper3_Z)*z;
        G4Transform3D transform_cu3 = G4Transform3D(rotm,position_cu3);
        new G4PVPlacement(transform_cu3, logic_copper3, "cu3", logicWorld, false,
			  4*n_bits*iplaque+n_bits*iconn+ibit, fCheckOverlaps);

        G4ThreeVector position_cu4 = position_bit+(-connector_dy-3*c_copper1_XZ)*x+(-connector_dz+c_copper4_Z)*z;
        G4Transform3D transform_cu4 = G4Transform3D(rotm,position_cu4);
        new G4PVPlacement(transform_cu4, logic_copper4, "cu4", logicWorld, false,
			  4*n_bits*iplaque+n_bits*iconn+ibit, fCheckOverlaps);
      }
      // Putting connector and kapton strip
      new G4PVPlacement(transform, logic_connector, "connector", logicWorld, false,
			4*iplaque+iconn, fCheckOverlaps);
      new G4PVPlacement(transform_kapton, logic_strip, "kapton", logicWorld, false,
			4*iplaque+iconn, fCheckOverlaps);
    }
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
                    G4ThreeVector(0,0,-20*cm),         //at (0,0,0)
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
                    G4ThreeVector(0,0,-20*cm),         //at (0,0,0)
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
  logicWorld->SetVisAttributes (G4VisAttributes::GetInvisible());   
  logicCryst->SetVisAttributes (G4VisAttributes::GetInvisible()); 
  logicPatient->SetVisAttributes (G4VisAttributes::GetInvisible()); 

  // Colors
  G4Colour gray_black(0.2, 0.2, 0.2);
  G4VisAttributes* connector_color = new G4VisAttributes(gray_black);
  logic_connector -> SetVisAttributes(connector_color);
  logic_bit -> SetVisAttributes(connector_color);

//  G4Colour gray_black(0.2, 0.2, 0.2);
  G4VisAttributes* tube_color = new G4VisAttributes(gray_black);
  logic_tube -> SetVisAttributes(tube_color);

  G4Colour light_blue(0.7,0.7,1);
  G4VisAttributes* detector_color = new G4VisAttributes(light_blue);
  logic_DDSD -> SetVisAttributes(detector_color);

  G4Colour bronze(141./256,108./256,41./256);
  G4VisAttributes* bronze_color = new G4VisAttributes(bronze);
  logic_copper1 -> SetVisAttributes(bronze_color);
  logic_copper2 -> SetVisAttributes(bronze_color);
  logic_copper3 -> SetVisAttributes(bronze_color);
  logic_copper4 -> SetVisAttributes(bronze_color);

  G4Colour pcb_Green(0,153./256,77./256);
  G4VisAttributes* plastic_color = new G4VisAttributes(pcb_Green);
  logic_plastic -> SetVisAttributes(plastic_color);

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
