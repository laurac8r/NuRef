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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// $Id: TrackingAction.cc 69099 2013-04-18 12:25:19Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"
#include "G4LogicalVolume.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4ParticleTypes.hh"
#include "G4HadronicProcessType.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* event, HistoManager* histo)
:G4UserTrackingAction(), fDetector(det), fEventAction(event), fHisto(histo),
  fScoringVolume(0),fScoringVolume1(0.)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{  
  //count secondary particles
  if (track->GetTrackID() == 1) return;  
  G4String name   = track->GetDefinition()->GetParticleName();
  G4double energy = track->GetKineticEnergy();
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  run->ParticleCount(name,energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
 // Check the status of exiting particles only. Return only if the particle is
 // still within the world.
 G4StepStatus status = track->GetStep()->GetPostStepPoint()->GetStepStatus();
 if (status != fWorldBoundary) return; 

 // Get the particle associated with the track. Grab the name and energy as well.
 const G4ParticleDefinition* particle = track->GetParticleDefinition();
 G4String name   = particle->GetParticleName();
 G4double energy = track->GetKineticEnergy();
 
 // // Add the energy flow to the event action.
 // fEventAction->AddEflow(energy);  
 
 // Get the current run of the simulation..
 Run* run = static_cast<Run*>(
              G4RunManager::GetRunManager()->GetNonConstCurrentRun());

 // Get the geometry of the world
 const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

 // Add a flux modifier to the neutron energy.
 run->ParticleFlux(name,energy);               
 
 // Create a flux histogram for every energy flow.
 //
 G4AnalysisManager* analysis = G4AnalysisManager::Instance();

 G4LogicalVolume* track_volume = track->GetVolume()->GetLogicalVolume();

 // If the scoring volume vector is empty, use the base logical volume.
 if (fScoringVolume1.empty()){fScoringVolume1 = detectorConstruction->GetScoringVolume1();}

 // Loop through the scoring volumes.
 // for(int ih=0; ih < fScoringVolume1.size(); ih++)
 for(int ih=0; ih < 2; ih++)
  {
    
    // Check if we are in the scoring volume and particle is a neutron.
    if ( track_volume == fScoringVolume1[ih] && particle == G4Neutron::Neutron() )
      {
        analysis->FillH1(ih+6,energy);
      }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// TrackingAction::TrackingAction(EventAction* event, DetectorConstruction* det)
// :G4UserTrackingAction(), fEventAction(event), fDetector(det)
// { }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void TrackingAction::PreUserTrackingAction(const G4Track* track)
// {
//   Run* run = static_cast<Run*>(
//         G4RunManager::GetRunManager()->GetNonConstCurrentRun());

//   const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>
//       (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  
//   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

//   //which volume ?
//   //
//   G4LogicalVolume* lVolume = track->GetVolume()->GetLogicalVolume();
//   G4int iVol = 0;
//   // if (lVolume == fDetector->GetLogicTarget())   iVol = 1;
//   // if (lVolume == fDetector->GetLogicDetector()) iVol = 2;
    
//   //secondary particles only
//   if (track->GetTrackID() == 1) return;
  
//   const G4ParticleDefinition* particle = track->GetParticleDefinition();  
//   G4String name   = particle->GetParticleName();
//   G4int pid       = particle->GetPDGEncoding();
//   G4int Z         = particle->GetAtomicNumber();
//   G4int A         = particle->GetAtomicMass();
//   G4double charge = particle->GetPDGCharge();    
//   G4double energy = track->GetKineticEnergy();
//   G4double time   = track->GetGlobalTime();
//   G4double weight = track->GetWeight();
  
//   run->ParticleCount(name,energy,iVol);
  
//   //Radioactive decay products
//   G4int processType = track->GetCreatorProcess()->GetProcessSubType();
//   if (processType == fRadioactiveDecay) {
//     //fill ntuple id = 3
//     G4int id = 3;
//     analysisManager->FillNtupleDColumn(id,0, double(pid));
//     analysisManager->FillNtupleDColumn(id,1, double(Z));
//     analysisManager->FillNtupleDColumn(id,2, double(A));
//     analysisManager->FillNtupleDColumn(id,3, energy);
//     analysisManager->FillNtupleDColumn(id,4, time/s);
//     analysisManager->FillNtupleDColumn(id,5, weight);
//     analysisManager->AddNtupleRow(id);
    
//     if (charge < 3.) {   
//       //fill ntuple id = 0
//       id = 0;
//       analysisManager->FillNtupleDColumn(id,0, double(pid));
//       analysisManager->FillNtupleDColumn(id,1, energy);
//       analysisManager->FillNtupleDColumn(id,2, time/s);
//       analysisManager->FillNtupleDColumn(id,3, weight);
//       analysisManager->AddNtupleRow(id);
    
//       analysisManager->FillH1(6, energy, weight);
//     }                        
//   }
  
//   //all unstable ions produced in target
//   G4bool unstableIon = ((charge > 2.) && !(particle->GetPDGStable()));
//   if ((unstableIon) && (iVol == 1)) {
//     //fill ntuple id = 1
//     G4int id = 1;
//     analysisManager->FillNtupleDColumn(id,0, double(pid));
//     analysisManager->FillNtupleDColumn(id,1, time/s);
//     analysisManager->FillNtupleDColumn(id,2, weight);
//     analysisManager->AddNtupleRow(id);  
//   }
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void TrackingAction::PostUserTrackingAction(const G4Track* )
// { }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

