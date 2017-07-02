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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 70756 2013-06-05 12:20:06Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4AccumulableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"

#include "Randomize.hh"
#include <iomanip>
#include <assert.h>

#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(const G4String& outputFile)
: G4UserRunAction(), fOutputFileSpec(outputFile)
// fEdep(0.),
// fEdep2(0.)
{

// // add new units for dose
//   //
//   const G4double milligray = 1.e-3*gray;
//   const G4double microgray = 1.e-6*gray;
//   const G4double nanogray  = 1.e-9*gray;
//   const G4double picogray  = 1.e-12*gray;
//
//   new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
//   new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
//   new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
//   new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);
//
//   // Register accumulable to the accumulable manager
//   G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
//   accumulableManager->RegisterAccumulable(fEdep);
//   accumulableManager->RegisterAccumulable(fEdep2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void RunAction::BeginOfRunAction(const G4Run* aRun)
G4Run* RunAction::GenerateRun()
{
  // Histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open a ROOT file if the analysis manager is active.
  if (analysisManager->IsActive()) analysisManager->OpenFile();

  // Save the random number store status.
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  if (isMaster) G4Random::showEngineStatus();

  // // reset accumulables to their initial values
  // G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  // accumulableManager->Reset();

  // Generate a new run.
  return new Run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{

  // If a visul manager exists.
  //
  if (G4VVisManager::GetConcreteInstance())
  {
    // Refreshes the visual viewer
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  // Print the number of events that occurred in the run.
  G4cout << "Number of Events Processed: " << aRun->GetNumberOfEvent() << G4endl;

  // Cast the G4Run pointer onto a Run pointer. The cast pointer is expected to
  //   exist. An error is thrown and the program exits if the cast fails.
  const Run* theRun = dynamic_cast<const Run*>(aRun);
  assert (0 != theRun);

  // G4int nofEvents = aRun->GetNumberOfEvent();
  // if (nofEvents == 0) return;

  // Get an instance of the analysis manager to save the histograms.
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // If the analysis manager is active, write the histograms to the ROOT file
  //   and close it.
  if (analysisManager->IsActive())
  {
    analysisManager->Write();
    analysisManager->CloseFile();
  }

  // Dump the fluence data into the output file.
  theRun->DumpData(fOutputFileSpec);
}

  // // Merge accumulables
  // G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  // accumulableManager->Merge();

  // // Compute dose = total energy deposit in a run and its variance
  // //
  // G4double edep  = fEdep.GetValue();
  // G4double edep2 = fEdep2.GetValue();
  //
  // G4double rms = edep2 - edep*edep/nofEvents;
  // if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;
  //
  // const DetectorConstruction* detectorConstruction
  //  = static_cast<const DetectorConstruction*>
  //    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());//you have to et the mass by callin the detector construction
  //
  // G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  // //G4cout<<"mass is "<<G4BestUnit(mass,"Mass")<<G4endl;
  // G4double dose = edep/mass;
  // G4double rmsDose = rms/mass;
  //
  // // Run conditions
  //
  // G4cout
  //    << G4endl
  //    << " The run consists of " << nofEvents<<" particle(s)"
  //    << G4endl
  //    << " Cumulated dose per run, in scoring volume : "
  //    << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
  //    << G4endl
  //    << "------------------------------------------------------------"
  //    << G4endl
  //    << G4endl;

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
void RunAction::AddEngDep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}
