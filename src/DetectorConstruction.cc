//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The Geant4 software  is  copyright of the Copyright Holders  of *
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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4GenericTrap.hh"
#include "G4Para.hh"

#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4RotationMatrix.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fTargetMater(0), fLogicTarget(0),
 fDetectorMater(0), fLogicDetector(0), 
 fWorldMater(0), fPhysiWorld(0),
 fDetectorMessenger(0)
{
  fTargetLength      = 1*cm; 
  fTargetRadius      = 0.5*cm;
  fDetectorLength    = 5*cm; 
  fDetectorThickness = 2*cm;
  
  fWorldLength = std::max(fTargetLength,fDetectorLength);
  fWorldRadius = fTargetRadius + fDetectorThickness;
      
  DefineMaterials();
    
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // build materials
  //
  fDetectorMater = 
  new G4Material("Germanium", 32, 72.61*g/mole, 5.323*g/cm3);
  

  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
  //
  G4int ncomponents; G4double fractionmass;      
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2,
                      kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
  //
  fWorldMater = Air20;
  
  // or use G4 materials data base
  //
  G4NistManager* man = G4NistManager::Instance();  
  fTargetMater = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");
                   
  ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

/*  
  // World
  //
  // (re) compute World dimensions if necessary
  fWorldLength = std::max(fTargetLength,fDetectorLength);
  fWorldRadius = fTargetRadius + fDetectorThickness;
    
  G4Tubs*
  sWorld = new G4Tubs("World",                                 //name
                 0.,fWorldRadius, 0.5*fWorldLength, 0.,twopi); //dimensions  
                   
  G4LogicalVolume*
  lWorld = new G4LogicalVolume(sWorld,                  //shape
                             fWorldMater,               //material
                             "World");                  //name

  fPhysiWorld = new G4PVPlacement(0,                    //no rotation
                            G4ThreeVector(),            //at (0,0,0)
                            lWorld,                     //logical volume
                            "World",                    //name
                            0,                          //mother volume
                            false,                      //no boolean operation
                            0);                         //copy number

  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    lWorld,                  //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

*/
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 10*m, env_sizeZ = 10*m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* sWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* lWorld =                         
    new G4LogicalVolume(sWorld,              //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  fPhysiWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      lWorld,                //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    lWorld,                  //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //
  // Shrould and shells
  //

  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  G4ThreeVector pos_shrould = G4ThreeVector();

  //Main Block
  G4VSolid* main_block =
    new G4Tubs("Main Block",                            // Name
                0, 100.33*mm, 88.9*mm, 0, twopi);       // Size

  //G4RotationMatrix* rs = new G4RotationMatrix();
  //rs->rotateY(90.*deg);

  G4VSolid* shave_block =
    new G4Box("Shave Block",
              34.415*mm, 100.33*mm, 90.*mm);

  G4ThreeVector pos_shave_block1(65.915*mm, 0, 0);
  G4ThreeVector pos_shave_block2(-65.915*mm, 0, 0);
  
  G4VSolid* inside_cyl =
    new G4Tubs("Inside Cylinder",
               0, 26.046*mm, 80.01*mm, 0, twopi);

  G4ThreeVector pos_inside_cyl1(0, -64.146*mm, -3.81*mm);
  G4ThreeVector pos_inside_cyl2(0, 64.146*mm, -3.81*mm);

  G4VSolid * inside_block =
    new G4Box("Inside Block",
              31.115*mm, 60.146*mm, 80.01*mm); // The walls are actually thinner.
                                               // Set x to 31.115*mm, but 26.046 matches the shape
                                               // of the inside cylinders. y dimension 64.146 appears
                                               // to have trouble ->60
  G4ThreeVector pos_inside_block(0, 0, -3.81*mm);

  // Side Rectangle
  G4VSolid* side_rect =
    new G4Box("Side Rectangle",
              3.175*mm, 10.5*mm, 31.75*mm);
  G4ThreeVector pos_side_rect(0, 100.3*mm, 0);

  // Top Stuff geometry
  G4VSolid* top_cyl =
    new G4Tubs("Top Cylinder",
               0, 6.0325*mm, 16*mm, 0, twopi);

  G4ThreeVector pos_top_cyl1(9.017*mm, -70.358*mm, 80.9*mm);
  G4ThreeVector pos_top_cyl2(0, -36.5*mm, 80.9*mm);
  G4ThreeVector pos_top_cyl3(0, 36.5*mm, 80.9*mm);
  G4ThreeVector pos_top_cyl4(-9.017*mm, 70.369*mm, 80.9*mm);

  // Center hole
  G4VSolid* cent_hole =
    new G4Tubs("Center Hole",
               0, 10.1205*mm, 32.*mm, 0, twopi);
  G4RotationMatrix* rm = new G4RotationMatrix();

  //rm->rotateZ(90.*deg);
  rm->rotateY(90.*deg);
  G4ThreeVector Tr(0, 0, 0);
  //MIC_____________________________
  G4VSolid* MIC_body =
    new G4Box("MIC body",
      98.044*mm, 185.039*mm, 12.7*mm);
  G4ThreeVector pos_MIC_body(0, 0, 0);//height of 36.703 + 12.7 = 49.403

  G4VSolid* MIC_cyl = 
    new G4Tubs("Subt Cyl",
      0, 38.227*mm, 12.8*mm, 0, twopi);
  G4ThreeVector pos_MIC_cyl1(0, 63.627*mm, 0);
  G4ThreeVector pos_MIC_cyl2(0, -63.627*mm, 0);

  G4VSolid* MIC_rect =
    new G4Box("Subt Rect",
      38.227*mm, 63.627*mm, 12.8*mm);
  G4ThreeVector pos_MIC_rect(0, 0, 0);

  //MIC Subtraction
  G4SubtractionSolid* m1 =
    new G4SubtractionSolid("M1", MIC_body, MIC_cyl, 0, pos_MIC_cyl1);
  G4SubtractionSolid* m2 =
    new G4SubtractionSolid("M2", m1, MIC_cyl, 0, pos_MIC_cyl2);
  G4SubtractionSolid* m3 =
    new G4SubtractionSolid("M3", m2, MIC_rect, 0, pos_MIC_rect);

  // Subtraction solids
  G4SubtractionSolid* s1 =
    new G4SubtractionSolid("S1", main_block, shave_block,  0, pos_shave_block1);
  
  G4SubtractionSolid* s2 =
    new G4SubtractionSolid("S2", s1, shave_block, 0, pos_shave_block2);
 
  G4SubtractionSolid* s3 =
    new G4SubtractionSolid("S3", s2, inside_cyl, 0, pos_inside_cyl1);
  
  G4SubtractionSolid* s4 =
    new G4SubtractionSolid("S4", s3, inside_cyl, 0, pos_inside_cyl2);
  
  G4SubtractionSolid* s5 =
    new G4SubtractionSolid("S5", s4, inside_block, 0, pos_inside_block);
  
  G4SubtractionSolid* s6 =
    new G4SubtractionSolid("S6", s5, side_rect, 0, pos_side_rect);
  
  G4SubtractionSolid* s7 =
    new G4SubtractionSolid("S7", s6, top_cyl, 0, pos_top_cyl1);
  
  G4SubtractionSolid* s8 =
    new G4SubtractionSolid("S8", s7, top_cyl, 0, pos_top_cyl2);
  
  G4SubtractionSolid* s9 =
    new G4SubtractionSolid("S9", s8, top_cyl, 0, pos_top_cyl3);
  
  G4SubtractionSolid* s10 =
    new G4SubtractionSolid("S10", s9, top_cyl, 0, pos_top_cyl4);
  
  G4SubtractionSolid* s11 =
    new G4SubtractionSolid("S11", s10, cent_hole, rm, Tr);

  //Logic Volume

  G4LogicalVolume* logical_shrould =
    new G4LogicalVolume(s11,                  // Its solid volume
                        aluminum,             // Its material
                        "Logical Shrould");   // Its name
  
  G4LogicalVolume* logical_MIC =
    new G4LogicalVolume(m3,
      aluminum,
      "Logical MIC");
  new G4PVPlacement(0,                        // No rotation
                    G4ThreeVector(0, 0, 0),   // At position
                    logical_shrould,          // Its logical volume
                    "Physical Shrould",       // Its name
                    logicEnv,                 // Its mother volume
                    true,                     // Boolean operation
                    0,                        // Copy number
                    checkOverlaps); // Overlaps checking

  new G4PVPlacement(0,
    G4ThreeVector(0, 0, 138.303*mm), // 49.403 + 88.9 = 138.303
    logical_MIC,
    "Main Insulator Cap",
    logicEnv,
    true,
    0,
    checkOverlaps);
  //
  // Neutron moderator/reflector
  //

  //
  // Outer Box
  G4Material* box_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4ThreeVector pos1 = G4ThreeVector();
  G4double boxsize = .75*m;

  G4Box* outer_boxshape = 
    new G4Box("Outer_Box", boxsize, boxsize, boxsize);
  
  // G4LogicalVolume* outer_boxlogic = new G4LogicalVolume(outer_boxshape,box_mat,"Outer_Box");               

  //
  // Inner Box
  G4Material* innerbox_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4double innerboxsize = .70*m;

  G4Box* inner_boxshape =
    new G4Box("Inner_Box", innerboxsize, innerboxsize, innerboxsize);
  // G4LogicalVolume* inner_boxlogic = new G4LogicalVolume(inner_boxshape,innerbox_mat,"Inner_Box");

  //
  // Beam Holes
  G4double hole_rad = .05*m;
  G4double hole_len = .8*m;

  G4Tubs* beamhole_shape =
    new G4Tubs("Beam_hole",
               0.*m,
               hole_rad,
               hole_len,
               0.*deg,
               360.*deg);

  // G4LogicalVolume* beamhole_logic = new G4LogicalVolume(beamhole_shape,innerbox_mat,"Beam_hole");
  // G4Tubs beamhole_shape("Beam_hole",0.*m,hole_rad,hole_len,0.*deg,360.*deg);
  // G4Box inner_boxshape("Inner_Box",innerboxsize,innerboxsize,innerboxsize);
  // G4Box outer_boxshape("Outer_Box",boxsize,boxsize,boxsize);  

  G4SubtractionSolid* shell_noholes =
    new G4SubtractionSolid("Shell_noholes",
                           outer_boxshape,
                           inner_boxshape);

  G4SubtractionSolid* shell_withholes =
    new G4SubtractionSolid("ShellSol",shell_noholes,beamhole_shape);

  G4LogicalVolume* shell_withholesL =
    new G4LogicalVolume(shell_noholes,box_mat,"ShellLog");

  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    shell_withholesL,        //its logical volume
                    "Shell",                 //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  // Rotation
  G4RotationMatrix* rotate_object =          // Define the rotation matrix
    new G4RotationMatrix();

  rotate_object->rotateZ(pi/8);         // Perform the rotation operations
  rotate_object->rotateY(pi/2);         // around the specified axes using
                                        // the right-hand rule for each axis.
                                        // Rotations are performed in backward
                                        // order: Y and then Z in this case.

  //
  // Target
  //
  G4Tubs* 
  sTarget = new G4Tubs("Target",                                   //name
                  0., fTargetRadius, 0.5*fTargetLength, 0.,twopi); //dimensions


  fLogicTarget = new G4LogicalVolume(sTarget,           //shape
                             fTargetMater,              //material
                             "Target");                 //name
                               
           new G4PVPlacement(rotate_object,             //rotation
                            G4ThreeVector(0, 1*m, 1*m), //location
                            fLogicTarget,               //logical volume
                            "Target",                   //name
                            lWorld,                     //mother  volume
                            false,                      //no boolean operation
                            0);                         //copy number

  //
  // Detector
  //
  G4Tubs* 
  sDetector = new G4Tubs("Detector",  
                fTargetRadius, fWorldRadius, 0.5*fDetectorLength, 0.,twopi);


  fLogicDetector = new G4LogicalVolume(sDetector,       //shape
                             fDetectorMater,            //material
                             "Detector");               //name
                               
           new G4PVPlacement(rotate_object,             //no rotation
                            G4ThreeVector(0, 2*m, 2*m), //location
                            fLogicDetector,             //logical volume
                            "Detector",                 //name
                            lWorld,                     //mother  volume
                            false,                      //no boolean operation
                            0);                         //copy number


  PrintParameters();
  
  //always return the root volume
  //
  return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
         << " Radius = " << G4BestUnit(fTargetRadius,"Length")  
         << " Material = " << fTargetMater->GetName();
  G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
         << " Tickness = " << G4BestUnit(fDetectorThickness,"Length")  
         << " Material = " << fDetectorMater->GetName() << G4endl;          
  G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fTargetMater = pttoMaterial;
    if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    fDetectorMater = pttoMaterial;
    if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
  fTargetRadius = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
  fTargetLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
  fDetectorThickness = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
  fDetectorLength = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
  return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
  return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
  return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
  return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
  return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
  return fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
  return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
  return fLogicDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
