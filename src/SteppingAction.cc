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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 71404 2013-06-14 16:56:38Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
//#include "GeneralParticleSource.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "G4Event.hh"
#include "HistoManager.hh"
#include "G4Step.hh"
#include "G4LogicalVolume.hh"
#include "G4GDMLParser.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "G4TouchableHandle.hh"

#include "G4Neutron.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event, HistoManager* histo)
: G4UserSteppingAction(), fDetector(det), fEventAction(event), fHisto(histo),
  fScoringVolume(0),fScoringVolume1(0.)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  //get scoring volume
  //
  const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  if (fScoringVolume1.empty()){fScoringVolume1 = detectorConstruction->GetScoringVolume1();}
  
  // Get the post and pre step points for a track.
  G4StepPoint* post = aStep->GetPostStepPoint();
  G4StepPoint* pre = aStep->GetPreStepPoint();

  // Obtain the name of the particle.
  G4String particle = aStep->GetTrack()->GetDefinition()->GetParticleName();
 
  G4Track * theTrack = aStep->GetTrack();

  // Get the process object.
  const G4VProcess* process = post->GetProcessDefinedStep();

  // Get the process name.
  G4String processName = process->GetProcessName();

  // check if it is alive
  if( theTrack->GetTrackStatus()!=fAlive && particle =="neutron")
    {
      G4cout<<"Neutron is dead"<<G4endl;
      
      G4cout<<"The process in the step is: "<<processName<<G4endl;}

      // Touch the volume for binning.
      // G4TouchableHandle touch_handle = post->GetTouchableHandle();

      // // Obtain the name of the volume for binning
      G4String name = pre->GetPhysicalVolume()->GetLogicalVolume()->GetName();
      // G4LogicalVolume* step_vol = pre->GetPhysicalVolume()->GetLogicalVolume();

      // Obtain the initial kinetic energy of the current event.
      const G4Event* event = static_cast<const G4Event*>(G4RunManager::GetRunManager()->GetCurrentEvent());
      G4PrimaryVertex* primaryVertex = event->GetPrimaryVertex();
      G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
      G4double initial_ke = primaryParticle->GetKineticEnergy();

      // Obtain the kinetic energy of the particle in the current step
      G4double current_ke = aStep->GetTrack()->GetKineticEnergy();

      G4double deposited_ke = initial_ke - current_ke;

      if (processName == "Transportation") { analysisManager->FillH1(8, current_ke); }

      // Check if the pre-step point in the track is the volume of interest for binning.
      // If so, bin the particle if it is a neutron.
      // for(int ih=0; ih < fScoringVolume1.size(); ih++)
      for(int ih=0; ih < 2; ih++)
        {
          //Check if we are in the scorin volume and particle is neutron
          // if ( step_vol == fScoringVolume1[i] && particle == "neutron" )
          if ( name == fScoringVolume1[ih]->GetName() && particle == "neutron" )
            {
              if (processName == "nCapture") { analysisManager->FillH1(ih+4, deposited_ke); }

              // Bin any neutron that have lost kinetic energy and thus have undergone scattering.
              if (deposited_ke > 0) { analysisManager->FillH1(ih+2, current_ke); }

              //Check if the particle just entered the volume, ie it its first step in the volume
              if (pre->GetStepStatus() == fGeomBoundary)
                {

                  // Fill the histogram for the kinetic energy.
                  analysisManager->FillH1(ih, current_ke);
                  //G4cout<<"vol is: "<<fScoringVolume1[i]->GetName()<<" Energy is: "<<kinetic_energy<<G4endl;


                  //G4cout<<"E is: "<<kinetic_energy<<G4endl;
                  //G4double energy_deposited = aStep->GetTotalEnergyDeposit();
                  //analysisManager->FillH1(3,energy_deposited);
                }
          }
  }
}
// void SteppingAction::UserSteppingAction(const G4Step* aStep)
// {
//   Run* run = static_cast<Run*>(
//         G4RunManager::GetRunManager()->GetNonConstCurrentRun());    
  
//   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

// //get scoring volume
// //
// const DetectorConstruction* detectorConstruction
//       = static_cast<const DetectorConstruction*>
//         (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//    if (!fScoringVolume) { 
//     fScoringVolume = detectorConstruction->ReturnVolume();   
//   }

//   // get volume of the current step
//   //
//   G4LogicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()
//                              ->GetVolume()->GetLogicalVolume();
  
//  // check if we are in scoring volume
//   if (volume == fScoringVolume){ 

// // collect energy deposited in this step
//   G4double edepStep = aStep->GetTotalEnergyDeposit();
// //G4cout<<"Edeposit in this step is:"<<G4BestUnit(edepStep,"Energy")<<G4endl;
//  fEventAction->AddEdep(edepStep);  
// }


// //vector alternative method for more than one volume 
// //
// const DetectorConstruction* detectorConstruction1
//       = static_cast<const DetectorConstruction*>
//         (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//    if (fScoringVolume1.empty()) { 
//     fScoringVolume1 = detectorConstruction1->GetScoringVolume1();   
//   }
//  // get volume of the current step
//   G4LogicalVolume* volume1 
//     = aStep->GetPreStepPoint()->GetTouchableHandle()
//       ->GetVolume()->GetLogicalVolume();
      
//   // check if we are in scoring volume
// for(int i =0;i<2;i++){
//   if (volume1 == fScoringVolume1[i]) 
// {
// //G4cout<<"volume is "<<volume1->GetName()<<G4endl;
// //getchar();
//  G4double edepStep1 = aStep->GetTotalEnergyDeposit();
//  fEventAction->AddEdep1(edepStep1,i); 
// //G4cout<<"Edeposit in this step is:"<<G4BestUnit(edepStep1,"Energy")<<G4endl;
// }
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......