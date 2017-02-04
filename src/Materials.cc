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

#include "DetectorConstruction.hh"

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
  G4double env_sizeXY = 10 * m, env_sizeZ = 10 * m;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ = 1.2*env_sizeZ;
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
  G4Material* calcium_oxide = nist->FindOrBuildMaterial("G4_CALCIUM_OXIDE");
  G4Material* iron = nist->FindOrBuildMaterial("G4_Fe");
  G4Material* chromium = nist->FindOrBuildMaterial("G4_Cr");

  // Defining an element using its isotopic composition
  G4Isotope* boron_10 =
    new G4Isotope("B-10",                 // Name
      5,                      // Atomic (proton) number
      10,                     // Nucleon (mass) number
      10.012936992*g / mole);   // Molar mass (atomic mass in g/mole)

  G4Isotope* boron_11 =
    new G4Isotope("B-11",                 // Name
      5,                      // Atomic (proton) number
      11,                     // Nucleon (mass) number
      11.009305406*g / mole);   // Molar mass (atomic mass in g/mole)

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
      1.00794*g / mole);        // Molar mass

  G4Element* oxygen =
    new G4Element("Oxygen",               // Name
      "O",                    // Element symbol
      8,                      // Atomic number
      15.9994*g / mole);        // Molar mass

  G4Element* nitrogen =
    new G4Element("Nitrogen",
      "N",
      7,
      14.007*g / mole);
  G4Element* carbon =
    new G4Element("Carbon",
      "C",
      6,
      12.0107*g / mole);
  

                    // Defining a compound using its molecular stoichiometry
  G4Material* boric_acid =
    new G4Material("Boric Acid",          // Name
      2.46*g / cm3,            // Density
      3);                    // Number of components

  boric_acid->AddElement(boron,           // Name of first element in the compound
    1);              // Number of atoms of the element in a molecule

  boric_acid->AddElement(hydrogen,        // Name of second element in the compound
    3);              // Number of atoms of the element in a molecule

  boric_acid->AddElement(oxygen,          // Name of third element in the compound
    3);              // Number of atoms of the element in a molecule

             // Defining a mixture using its components
             // The composite sample with an arbitrary composition is chosen here.
  
  G4Material* epoxy_resin =
    new G4Material("Epoxy Resin",          // Name
      1.1628*g / cm3,            // Density
      3);
  epoxy_resin->AddElement(carbon,
    21);
  epoxy_resin->AddElement(hydrogen,
    24);
  epoxy_resin->AddElement(oxygen,
    4);

  G4Material* epoxy_hardener =
    new G4Material("Epoxy Hardener",          // Name
      0.922*g / cm3,            // Density
      3);
  epoxy_hardener->AddElement(carbon,
    10);
  epoxy_hardener->AddElement(hydrogen,
    22);
  epoxy_hardener->AddElement(nitrogen,
    2);
  
  G4double dens_comp = 10 * g / cm3;
  G4double dens_steel = 7.9 * g / cm3;

  G4double fract_mass_boric_acid = 50 * perCent;
  G4double fract_mass_iron = 89.5 * perCent;
  G4double fract_mass_epoxy_resin = 30.6 * perCent;
  G4double fract_mass_epoxy_hardener = 19.4 * perCent;
  G4double fract_mass_chromium = 10.5 * perCent;

  G4Material* steel =
    new G4Material("Absorber",           // Name
      dens_steel,             // Density
      2);
  steel->AddMaterial(iron,                // Name
    fract_mass_iron);
  steel->AddMaterial(chromium,                // Name
    fract_mass_chromium);


  G4Material* absorber =
    new G4Material("Absorber",           // Name
      dens_comp,             // Density
      3);                    // Number of components

  absorber->AddMaterial(boric_acid,                // Name
    fract_mass_boric_acid);    // Mass fraction

                   // NOTE: The steel material NEEDS to be redefined. Iron is chosen as a
                   // placeholder below just to make the overall code compile correctly.
                   // Steel's elemental composition must be defined from scratch for complete
                   // accuracy.

  absorber->AddMaterial(epoxy_resin,                     // Name
    fract_mass_epoxy_resin);         // Mass fraction

                   // NOTE: The epoxy material NEEDS to be redefined. Polyethylene is
                   // chosen as a placeholder below just to make the overall code compile
                   // correctly. Epoxy's elemental composition must be defined from scratch
                   // for complete accuracy.

  absorber->AddMaterial(epoxy_hardener,                     // Name
    fract_mass_epoxy_hardener);         // Mass fraction

  G4Material* test1 =
    new G4Material("test1",           // Name
      dens_comp,             // Density
      4);                    // Number of components
  test1->AddMaterial(steel,                // Name
    33.3 * perCent);

  test1->AddMaterial(boric_acid,                // Name
    33.3 * perCent);    // Mass fraction

  test1->AddMaterial(epoxy_resin,                     // Name
    20.6 * perCent);         // Mass fraction

  test1->AddMaterial(epoxy_hardener,                     // Name
    12.8 * perCent);         // Mass fraction

  G4Material* test2 =
    new G4Material("test2",           // Name
      dens_comp,             // Density
      3);                    // Number of components
  
  test2->AddMaterial(boric_acid,                // Name
    25.0 * perCent);    // Mass fraction

  test2->AddMaterial(epoxy_resin,                     // Name
    45.8 * perCent);         // Mass fraction

  test2->AddMaterial(epoxy_hardener,                     // Name
    29.2 * perCent);         // Mass fraction

  G4Material* test3 =
    new G4Material("test3",           // Name
      dens_comp,             // Density
      4);                    // Number of components
  test3->AddMaterial(steel,                // Name
    33.0 * perCent);

  test3->AddMaterial(boric_acid,                // Name
    17.0 * perCent);    // Mass fraction

  test3->AddMaterial(epoxy_resin,                     // Name
    30.6 * perCent);         // Mass fraction

  test3->AddMaterial(epoxy_hardener,                     // Name
    19.4 * perCent);         // Mass fraction

                   // Place SOLID VOLUME code for item of interest here.
                   //

                   // Various shapes for the solid volume are placed here as example references.
                   //
                   // NOTE: Each shape requires pointers to the header file for that shape. A
                   // call to the header file MUST be made at the beginning of this source file.
                   // The calls to the shapes in the examples below have been included.

                   // Box shape
  G4double x_half_dimension = 3.2*cm;
  G4double y_half_dimension = 3.2*cm;
  G4double z_half_dimension = 0.65*cm;

  G4Box* box =
    new G4Box("Box",                  // Its name
      x_half_dimension,       // Half the intended x-dimension
      y_half_dimension,       // Half the intended y-dimension
      z_half_dimension);      // Half the intended z-dimension

  
                 // Place LOGICAL VOLUME code for item of interest here.
                 //
  G4LogicalVolume* logical_object =
    new G4LogicalVolume(box,                // Its solid volume
      absorber,               // Its material
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