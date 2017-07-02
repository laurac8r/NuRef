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
#include "globals.hh"

#include "G4GDMLParser.hh"
#include "DetectorMessenger.hh"

#include "G4SDManager.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSCellFlux.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fWorldMater(0), fPhysiWorld(0),
 fDetectorMessenger(0),fScoringVolume(0),
 fScoringVolumeVec(0)
{
  fReadFile ="test3.gdml";
  fWritingChoice=1;
  DefineMaterials();

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4VPhysicalVolume* fWorldPhysVol;
    fParser.Read(fReadFile,false);
    fWorldPhysVol = fParser.GetWorldVolume();
    fScoringVolume = fParser.GetVolume("Target_logical");
    fScoringVolumeVec.push_back(fParser.GetVolume("TargetV4_Assem-1_Target_Faceplate2-2_logical"));
    fScoringVolumeVec.push_back(fParser.GetVolume("Target_logical"));

    ConstructSDandField();

  return fWorldPhysVol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials(){
  // build materials

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

  // G4NistManager* nist = G4NistManager::Instance();

  // G4Material* fWorldMater = nist->FindOrBuildMaterial("G4_Al");

  // or use G4 materials data base

  // G4NistManager* mat_manager = G4NistManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  return fPhysiWorld;
}


// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// SetReadFile
//
void DetectorConstruction::SetReadFile( const G4String& File )
{
  fReadFile=File;
  fWritingChoice=0;
}

void DetectorConstruction::ConstructSDandField()
{
  // Create a new sensitive detector (SD) manager pointer.
  G4SDManager* SD_manager = G4SDManager::GetSDMpointer();

  // Set the verbose level of the SD manager to 1.
  SD_manager->SetVerboseLevel(1);

  // Create a reference to the SD cache pointer.
  G4MultiFunctionalDetector* &sensitiveDetector =
    fSensitiveDetectorCache.Get();

  if (!sensitiveDetector)
  {
    sensitiveDetector = new G4MultiFunctionalDetector("MyDetector");

    G4VPrimitiveScorer* flux_scorer;

    G4String filterName, particleName;
    G4SDParticleFilter* neutronFilter =
      new G4SDParticleFilter(filterName = "neutronFilter",
                             particleName = "neutron");

    flux_scorer = new G4PSCellFlux("neutron cell flux");
    flux_scorer->SetFilter(neutronFilter);
    sensitiveDetector->RegisterPrimitive(flux_scorer);
  }

  SD_manager->AddNewDetector(sensitiveDetector);
  fScoringVolume->SetSensitiveDetector(sensitiveDetector);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
