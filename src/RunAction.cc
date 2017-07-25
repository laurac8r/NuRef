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
// #include "DetectorConstruction.hh"
// #include "PrimaryGeneratorAction.hh"
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
: G4UserRunAction(),
  fOutputFileSpec(outputFile), fEngDep(0.), fEngDepSqr(0.), fNeutCap(0),
{
  // Create an instance of an accumulable manager.
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();

  // Register the energy deposit accumulable and the square of that quantity for
  //   the case of only a single scoring volume.
  accumulableManager->RegisterAccumulable(fEngDep);
  accumulableManager->RegisterAccumulable(fEngDepSqr);

  // Iterate through and register all the elements in the vectors containing the
  //   energy deposition accumulables and the squares of those quantities in the
  //   case of multiple scoring volumes.
  for (int i = 0; i < 2; i++)
  {
    accumulableManager->RegisterAccumulable(fEngDepArr[i]);
    accumulableManager->RegisterAccumulable(fEngDepSqrArr[i]);

    // accumulableManager->RegisterAccumulable(fEngDepArr.at(i));
    // accumulableManager->RegisterAccumulable(fEngDepSqrArr.at(i));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  // Histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open a ROOT file if the analysis manager is active.
  if (analysisManager->IsActive()) analysisManager->OpenFile();

  // Save the random number store status.
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  if (isMaster) G4Random::showEngineStatus();

  // Reset accumulables held in the accumulable manager to their initial values.
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{

  // If a visual manager exists.
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

  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0) return;

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

  // Merge the energy deposition accumulables.
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double engDep  = fEngDep.GetValue();
  G4double engDepSqr = fEngDepSqr.GetValue();

  G4double rms = engDepSqr - engDep*engDep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  // Run conditions
  //
  G4cout << G4endl << G4endl << " The run consists of " << nofEvents
         << " particle(s)" << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEngDep(G4double& engDep)
{
  // Accumulate the energy deposition from the event calling this method for
  //   the case of only one scoring volume.
  fEngDep  += engDep;
  fEngDepSqr += engDep * engDep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::AddEngDepArr(G4double (&engDep)[2])
{
  // Iterate through all the scoring volumes and accumulate the energy
  //   deposition from the event action calling this method.
  for (int i = 0; i < 2; i++)
  {
    // Grab the energy deposition for the scoring volume of interest. Throw an
    //   error if we are out of range of the energy deposition vector number
    //   of elements. (Use of ".at(i)" in place of "[i]" ensures this error
    //   checking.)
    G4double engDepForI = engDep[i];
    // engDepForI = engDep.at(i);

    // Accumulate the energy deposition from the event action calling this
    //   method for the case of multiple scoring volumes.
    fEngDepArr[i]  += engDepForI;
    fEngDepSqrArr[i] += engDepForI * engDepForI;
  }
}
