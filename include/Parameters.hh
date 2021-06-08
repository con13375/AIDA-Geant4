// ~~~~~~~~~~~~~~~~ Parameters ~~~~~~~~~~~~~~~~
#ifndef Parameters_h
#define Parameters_h 1
#include <chrono>
#include <random>       
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"



#include "globals.hh"

// ~~~~~~~~~~~~~~~~ World
constexpr G4double world_sizeXY = 100*mm;
constexpr G4double world_sizeZ  = 310*mm;

// ~~~~~~~~~~~~~~~~ plaque stuff
constexpr G4double DDSD_XY_num = 75.60/2.0;
constexpr G4double DDSD_XY = DDSD_XY_num*mm, DDSD_Z = 0.5*mm; //including inactive area
constexpr G4double detector_XY_num = 71.63/2;
constexpr G4double detector_XY = detector_XY_num*mm; // just active area
constexpr G4double Plastic_XY = 47.5*mm, Plastic_Z = 0.8*mm;

// ~~~~~~~~~~~~~~~~ connector stuff
constexpr G4double bit_size = 1.27*mm;
constexpr G4double connector_dx = (bit_size*34+1.55)/2*mm, connector_dy = 2.36*mm, connector_dz = 2.39*mm;
// In the connector I put 34 bits with a plastic connected to pbronce pieces, these parameters are the size of the whole needed in the connector for all 34 bits
constexpr G4double chole1_dx = (bit_size*34)/2*mm, chole1_dy = connector_dy, chole1_dz = connector_dz-1.55/2*mm;
constexpr G4double chole2_dx = (bit_size*26)/2*mm, chole2_dy = connector_dy, chole2_dz = connector_dz; // like shown in the pics, the larger connector piece covers some of the bits (4 on each side), so this deeper hole is smaller in x
constexpr G4double c_bits_dx = 0.47*bit_size, c_bits_dy = (1.105+1.255)*mm, c_bits_dz = 0.97*chole1_dz;

constexpr G4double c_copper1_XZ = 0.3*mm;
constexpr G4double c_copper1_Y = connector_dy+3*c_copper1_XZ; //upper
constexpr G4double c_copper2_XZ = c_copper1_XZ;
constexpr G4double c_copper2_Y = connector_dy+1*c_copper1_XZ;  //lower
constexpr G4double c_copper3_XY = c_copper1_XZ;
constexpr G4double c_copper3_Z = 0.5*(connector_dz+bit_size/2+c_copper1_XZ);
constexpr G4double c_copper4_XY = c_copper1_XZ;
constexpr G4double c_copper4_Z = 0.5*(connector_dz-bit_size/2+c_copper1_XZ);
 
constexpr G4int nb_plaques = 6;
constexpr G4double separation = 10*mm;
constexpr G4double AIDA_nose_Z = 50*mm;

// ~~~~~~~~~~~~~~~~ kapton stuff
constexpr G4double kapton_Y = 0.1*mm;
constexpr G4double cu_X = 0.15*mm, cu_Y = 0.009*mm, cu_sep = 0.635*mm;
constexpr G4int n_cu = 4; // this should be 68, but i put it lower while building to ease loading time

// ~~~~~~~~~~~~~~~~ tubes, separators and bolts
constexpr G4double tube_in = 0.825*mm, tube_out = 1.5*mm, tube_Z = (AIDA_nose_Z + world_sizeZ)/2.0;
constexpr G4double separator_R1 = tube_out, separator_R2 = 2*mm, separator_R3 = 2.75*mm, separator_Z = 4.15*mm;
constexpr G4double bolt_R1 = tube_out, bolt_R2 = 2*mm, bolt_R3 = 3.5*mm, bolt_Z1 = 2.65*mm, bolt_Z2 = 3.5*mm;
constexpr G4double nut_R = 0.8*mm, nut_Z = 0.5*(bolt_R3-bolt_R1);

constexpr G4int n_tubes = 4;

// ~~~~~~~~~~~~~~~~ Al and Mylar
constexpr G4double Al_case_dx = 5*mm;
constexpr G4double Al_case_xy = Plastic_XY+Al_case_dx+3*mm, Al_case_z = (AIDA_nose_Z + world_sizeZ)/2.0;

constexpr G4double Mylar_case_r1 = Al_case_xy*1.4143+1*mm, Mylar_case_dr = 0.03*mm, Mylar_case_z = Al_case_z;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
