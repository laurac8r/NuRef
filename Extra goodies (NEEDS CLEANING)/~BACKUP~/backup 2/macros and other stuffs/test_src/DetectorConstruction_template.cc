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
// $Id: HFNG_model_DetectorConstruction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file HFNG_model_DetectorConstruction.cc
/// \brief Implementation of the HFNG_model_DetectorConstruction class

#include "HFNG_model_DetectorConstruction.hh"

#include "G4RunManager.hh"
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

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4RotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HFNG_model_DetectorConstruction::HFNG_model_DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HFNG_model_DetectorConstruction::~HFNG_model_DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* HFNG_model_DetectorConstruction::Construct()
{  
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
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.8*world_sizeXY, 0.8*world_sizeXY, 0.8*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no Boolean operation
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
                    logicWorld,              //its mother  volume
                    false,                   //no Boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //
  // Geometry of interest
  //

  // Defining the material for the item
  //

  // Defining elements and compounds using the internal Geant4 database
  G4Material* aluminum = nist->FindOrBuildMaterial("G4_Al");
  G4Material* titanium = nist->FindOrBuildMaterial("G4_Ti");
  G4Material* polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* calcium_oxide = nist->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
  
  // Defining an element using its isotopic composition
  G4Isotope* boron_10 =
    new G4Isotope("B-10",                 // Name
                  5,                      // Atomic (proton) number
                  10,                     // Nucleon (mass) number
                  10.012936992*g/mole);   // Molar mass (atomic mass in g/mole)

  G4Isotope* boron_11 =
    new G4Isotope("B-11",                 // Name
                  5,                      // Atomic (proton) number
                  11,                     // Nucleon (mass) number
                  11.009305406*g/mole);   // Molar mass (atomic mass in g/mole)

  G4Element* boron =
    new G4Element("Natural Boron",        // Name
                  "B",                    // Element symbol
                  2);                     // Number of isotopes

  boron->AddIsotope(boron_10,             // First isotope
                        19.9*perCent);    // Relative abundance (mole percent composition)
                                          // Source: Wikipedia

  boron->AddIsotope(boron_11,             // Second isotope
                        80.1*perCent);    // Relative abundance (mole percent composition)
                                          // Source: Wikipedia

  // Defining elements using their atomic numbers and molar masses
  G4Element* hydrogen =
    new G4Element("Hydrogen",             // Name
                  "H",                    // Element symbol
                  1,                      // Atomic number
                  1.00794*g/mole);        // Molar mass

  G4Element* oxygen =
    new G4Element("Oxygen",               // Name
                  "O",                    // Element symbol
                  8,                      // Atomic number
                  15.9994*g/mole);        // Molar mass

  // Defining a compound using its molecular stoichiometry
  G4Material* boric_acid =
    new G4Material("Boric Acid",          // Name
                   2.46*g/cm3,            // Density
                   3);                    // Number of components

  boric_acid->AddElement(boron,           // Name of first element in the compound
                         1);              // Number of atoms of the element in a molecule

  boric_acid->AddElement(hydrogen,        // Name of second element in the compound
                         3);              // Number of atoms of the element in a molecule

  boric_acid->AddElement(oxygen,          // Name of third element in the compound
                         3);              // Number of atoms of the element in a molecule

  // Defining a mixture using its components
  // The composite sample with an arbitrary composition is chosen here.
  G4double dens_comp = 10*g/cm3;

  G4double fract_mass_boric_acid = 10*perCent;
  G4double fract_mass_steel = 20*perCent;
  G4double fract_mass_epoxy = 70*perCent;

  G4Material* composite =
    new G4Material("Composite",           // Name
                   dens_comp,             // Density
                   3);                    // Number of components

  composite->AddMaterial(boric_acid,                // Name
                         fract_mass_boric_acid);    // Mass fraction

  // NOTE: The steel material NEEDS to be redefined. Iron is chosen as a
  // placeholder below just to make the overall code compile correctly.
  // Steel's elemental composition must be defined from scratch for complete
  // accuracy.
  G4Material* steel = nist->FindOrBuildMaterial("G4_Fe");

  composite->AddMaterial(steel,                     // Name
                         fract_mass_steel);         // Mass fraction

  // NOTE: The epoxy material NEEDS to be redefined. Polyethylene is
  // chosen as a placeholder below just to make the overall code compile
  // correctly. Epoxy's elemental composition must be defined from scratch
  // for complete accuracy.
  G4Material* epoxy = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

  composite->AddMaterial(epoxy,                     // Name
                         fract_mass_epoxy);         // Mass fraction

  // Place SOLID VOLUME code for item of interest here.
  //

  // Various shapes for the solid volume are placed here as example references.
  //
  // NOTE: Each shape requires pointers to the header file for that shape. A
  // call to the header file MUST be made at the beginning of this source file.
  // The calls to the shapes in the examples below have been included.

  // Box shape
  G4double x_half_dimension = 0.5*cm;
  G4double y_half_dimension = 0.5*cm;
  G4double z_half_dimension = 0.5*cm;

  G4Box* box =
    new G4Box("Box",                  // Its name
              x_half_dimension,       // Half the intended x-dimension
              y_half_dimension,       // Half the intended y-dimension
              z_half_dimension);      // Half the intended z-dimension

  // Cylinderical (tubular) shape
  G4double r_min = 0;
  G4double r_max = 1*cm;

  G4double half_length = 0.5*cm;

  G4double start_angle = 0;
  G4double end_angle = twopi;

  G4Tubs* cylinder =
    new G4Tubs("Cylinder",            // Its name
               r_min,                 // Inner radius
               r_max,                 // Outer radius
               half_length,           // Half the intended length
               start_angle,           // Start angle of the cylindrical segment
               end_angle);            // End angle of the cylindrical segment

  // Conical shape
  G4double r_min_bottom = 0;
  G4double r_max_bottom = 0.4*cm;

  G4double r_min_top = 0;
  G4double r_max_top = 0.1*cm;

  G4double z_dim = 1*cm;

  start_angle = 0;
  end_angle = twopi;

  G4Cons* cone =
    new G4Cons("Cone",                // Its name
               r_min_bottom,          // First inner radius
               r_max_bottom,          // First outer radius
               r_min_top,             // Second inner radius
               r_max_top,             // Second outer radius
               z_dim,                 // Half the intended length
               start_angle,           // Start angle of the conical segment
               end_angle);            // End angle of the conical segment


  // Spherical shape
  r_min = 0;
  r_max = 0.6*cm;
  G4double phi_min = 0;
  G4double phi_max = twopi;
  G4double theta_min = 0;
  G4double theta_max = twopi;

  G4Sphere* sphere =
    new G4Sphere("Sphere",            // Its name
                 r_min,               // Inner radius
                 r_max,               // Outer radius
                 phi_min,             // Start phi angle
                 phi_max,             // End phi angle
                 theta_min,           // Start theta angle
                 theta_max);          // End theta angle

  // Triangular prism shape (uses general trapezoid constructor)

  // This example shape is not very useful for this project because
  // Boolean operations do not support it, but it is included as an
  // example of a more complex solid if the need arises in future
  // projects.

  // Note: Two vertices are collapsed into one vertex for each face
  // to create the triangular top and bottom faces.
  //
  // Vertices' coordinates are indicated for the xy-plane of the
  // respective face they lie in.
  G4ThreeVector pos_prism = G4ThreeVector();

  G4double half_height = 2.5*cm;

  std::vector<G4TwoVector> vertices;

  // Vertex positions are chosen to make an equilateral triangle for
  // each face.
  vertices.push_back(G4TwoVector(-3*mm, -3*mm));    // Bottom (-z) face
  vertices.push_back(G4TwoVector(-3*mm,  3*mm));
  vertices.push_back(G4TwoVector(2.2*mm, 0));
  vertices.push_back(G4TwoVector(2.2*mm, 0));
  vertices.push_back(G4TwoVector(-3*mm, -3*mm));     // Top (+z) face
  vertices.push_back(G4TwoVector(-3*mm,  3*mm));
  vertices.push_back(G4TwoVector(2.2*mm, 0));
  vertices.push_back(G4TwoVector(2.2*mm, 0));

  G4GenericTrap* tri_prism =
    new G4GenericTrap("Triangular Prism",     // Its name
                      half_height,            // Half the height
                      vertices);              // Vertices

  // Parallelepiped shape (this implementation mimics two triangular prisms
  // placed back-to-back)
  G4double half_x_dim = 0.15*cm;
  G4double half_y_dim = 0.10*cm;
  G4double half_z_dim = 1*cm;
  G4double alpha = pi/4;
  G4double theta = 0;
  G4double phi = 0;

  G4Para* parallelepiped =
    new G4Para("Parallelepiped",              // Name
               half_x_dim,                    // Half of intended length
               half_y_dim,                    // Half of intended width
               half_z_dim,                    // Half of intended height
               alpha,                         // Clockwise twist of yx-planes
               theta,                         // Polar angle between yx-planes
               phi);                          // Azimuthal angle between yx-planes

  // Perform an example of a Boolean subtraction operation of the parallelepiped from
  // the box with multiple rotations and a translation operation. The rotation and
  // translation order is irrelevant.  
  G4RotationMatrix* rotate_object = new G4RotationMatrix();   // Define the rotation matrix
  rotate_object->rotateY(pi/8);         // Perform the rotation operations
  rotate_object->rotateX(pi/2);         // around the specified axes using
                                        // the right-hand rule for each axis.
                                        // Rotations are performed in backward
                                        // order: X and then Y in this case.
  
  G4ThreeVector translate_object(0.5*cm, 0, 0);      // Define the translation vector

  G4SubtractionSolid* cut_box =
    new G4SubtractionSolid("Cut Box",           // Name
                           box,                 // Box
                           parallelepiped,      // Parallelepiped
                           rotate_object,       // Rotation
                           translate_object);   // Translation

  // Place LOGICAL VOLUME code for item of interest here.
  //
  G4LogicalVolume* logical_object =
    new G4LogicalVolume(cut_box,                // Its solid volume
                        titanium,               // Its material
                        "Logical Object");      // Its name

  // Place PHYSICAL VOLUME code for item of interest here.
  //
  new G4PVPlacement(0,
                    G4ThreeVector(),
                    logical_object,
                    "Physical Object",
                    logicEnv,
                    true,
                    0,
                    checkOverlaps);

  // Set the logical_object as the scoring volume.
  //
  fScoringVolume = logical_object;

  //
  // Always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
