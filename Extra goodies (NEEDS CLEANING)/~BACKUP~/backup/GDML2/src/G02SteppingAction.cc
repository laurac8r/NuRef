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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "G02SteppingAction.hh"
#include "G02DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G02EventAction.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4GDMLParser.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G02SteppingAction::G02SteppingAction(G02EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),fScoringVolume1(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G02SteppingAction::~G02SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G02SteppingAction::UserSteppingAction(const G4Step* step)
{
   const G02DetectorConstruction* detectorConstruction
      = static_cast<const G02DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
   if (!fScoringVolume) { 
    fScoringVolume = detectorConstruction->ReturnVolume();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  if (volume == fScoringVolume){ 
  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
//G4cout<<"Edeposit in this step is:"<<G4BestUnit(edepStep,"Energy")<<G4endl;
 fEventAction->AddEdep(edepStep);  
}
// for vector alternative method
const G02DetectorConstruction* detectorConstruction1
      = static_cast<const G02DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
   if (fScoringVolume1.empty()) { 
    fScoringVolume1 = detectorConstruction1->GetScoringVolume1();   
  }
 // get volume of the current step
  G4LogicalVolume* volume1 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
for(int i =0;i<2;i++){
  if (volume1 == fScoringVolume1[i]) 
{
//G4cout<<"volume is "<<volume1->GetName()<<G4endl;
//getchar();
 G4double edepStep1 = step->GetTotalEnergyDeposit();
 fEventAction->AddEdep1(edepStep1,i); 
//G4cout<<"Edeposit in this step is:"<<G4BestUnit(edepStep1,"Energy")<<G4endl;
}
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

