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

SteppingAction::SteppingAction(RunAction* runAction, EventAction* eventAction)
: G4UserSteppingAction(),
  fRunAction(runAction), fEventAction(eventAction),
  fScoringVolume(0), fScoringVolumeVec(0.)
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

  if (fScoringVolumeVec.empty()){fScoringVolumeVec = detectorConstruction->GetScoringVolumeVec();}

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
  // G4String processName = process->GetProcessDefinedStep()->GetName;

  // Touch the volume for binning.
  // G4TouchableHandle touch_handle = post->GetTouchableHandle();

  // // Obtain the name of the volume for binning
  G4String name = pre->GetPhysicalVolume()->GetLogicalVolume()->GetName();
  // G4LogicalVolume* step_vol = pre->GetPhysicalVolume()->GetLogicalVolume();

  // G4cout << "The name of the logical volume of particle death is: " << name << G4endl << G4endl;

  // Obtain the initial kinetic energy of the current event.
  const G4Event* current_event = static_cast<const G4Event*>(G4RunManager::GetRunManager()->GetCurrentEvent());
  G4PrimaryVertex* primaryVertex = current_event->GetPrimaryVertex();
  G4PrimaryParticle* primaryParticle = primaryVertex->GetPrimary();
  G4double initial_ke = primaryParticle->GetKineticEnergy();

  // Obtain the kinetic energy of the particle in the current step
  G4double current_ke = aStep->GetTrack()->GetKineticEnergy();

  // Obtain the deposited kinetic energy from the particle in the current step.
  G4double deposited_ke = aStep->GetTotalEnergyDeposit();

  // Add the current-step energy deposition to the event action energy
  //  deposition counter if and only if there is some energy deposition.
  if (deposited_ke > 0) {fEventAction->AddEngDep(deposited_ke);}

  // G4cout << "Initial kinetic energy of particle: " << initial_ke << G4endl << G4endl;
  // G4cout << "Current kinetic energy of particle: " << current_ke << G4endl << G4endl;
  // G4cout << "Deposited energy from particle kinetic energy: " << deposited_ke << G4endl << G4endl;

  // Check if the pre-step point in the track is the volume of interest for binning.
  // If so, bin the particle if it is a neutron.
  // for(int ih=0; ih < fScoringVolumeVec.size(); ih++)
  for(int ih=0; ih < 2; ih++)
    {
      // Check if we are in the scoring volume and the particle is a neutron.
      // if ( step_vol == fScoringVolumeVec[i] && particle == "neutron" )
      if ( name == fScoringVolumeVec[ih]->GetName() && particle == "neutron" )
        {
          // Check if the neutron just entered the volume. That is, if its
          //  first step is in the volume.
          if (pre->GetStepStatus() == fGeomBoundary)
            {
              // Bin the kinetic energy of the neutron at the scoring volume
              //  boundary.
              if (current_ke > 0) {analysisManager->FillH1(ih+1, current_ke);}

              //G4cout<<"vol is: "<<fScoringVolumeVec[i]->GetName()<<" Energy is: "<<kinetic_energy<<G4endl;

              //G4cout<<"Edeposit in this step is:"<<G4BestUnit(edepStep,"Energy")<<G4endl;

              // Add the current-step energy deposition in the logical volume
              //  of interest to the event action energy deposition counter
              //  array if and only if there is some energy deposition.
              if (deposited_ke > 0) {fEventAction->AddEngDepArr(deposited_ke, ih);}
            }

          // Bin the current kinetic energy of any neutron that have lost kinetic
          //  energy and thus have undergone scattering.
          if (deposited_ke > 0) {analysisManager->FillH1(ih+3, current_ke);}

          // Bin the kinetic energy the neutron has deposited in the
          //  scoring volume by the if it is non-zero.
          if (deposited_ke > 0) {analysisManager->FillH1(ih+5, deposited_ke);}
        }
    }

  // If the particle is alive and is a neutron, check if the particle has been
  //  transported out of this world or if a neutron capture has occured and
  //  either bin the current kinetic energy (if it is non-zero) or increment
  //  a neutron capture counter, respectively.
  if( theTrack->GetTrackStatus() != fAlive && particle == "neutron")
    {
      // G4cout<<"Neutron is dead"<<G4endl;

      // G4cout<<"The process in the step is: "<<processName<<G4endl;

      // If the neutron has been transported out of this world, bin the current
      //  kinetic energy.
      if (processName == "Transportation" && pre->GetStepStatus() == fWorldBoundary )
        {
          analysisManager->FillH1(0, current_ke);
        }

      // Increment the neutron capture counter if the neutron has been captured
      if (processName == "nCapture") {fRunAction->IncrementNeutCap();}
    }
}

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


// //vector alternative method for more than one volume
// //
// const DetectorConstruction* detectorConstruction1
//       = static_cast<const DetectorConstruction*>
//         (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
//    if (fScoringVolumeVec.empty()) {
//     fScoringVolumeVec = detectorConstruction1->GetScoringVolume1();
//   }
//  // get volume of the current step
//   G4LogicalVolume* volume1
//     = aStep->GetPreStepPoint()->GetTouchableHandle()
//       ->GetVolume()->GetLogicalVolume();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
