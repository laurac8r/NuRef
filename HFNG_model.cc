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
/// \file HFNG_model.cc
/// \brief Main program of the radioactivedecay/HFNG_model example
//
//
// $Id: HFNG_model.cc 98257 2016-07-04 17:39:46Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
  #include "G4MTRunManager.hh"
#else
  #include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4ScoringManager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
//#include "PhysicsList.hh"
#include "QGSP_BERT_HP.hh"
// #include "Shielding.hh"
#include "ActionInitialization.hh"
// #include "SteppingVerbose.hh"

#ifdef G4VIS_USE
  #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
  #include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Define default arguments for the executable.
  G4String macroFile = "None";
  G4String startingSeed = "1";
  G4String outputFile = "fluence.csv";

  // Parse the arguments for the executable depending on the number of
  // arguments passed to the executable.
  if (argc > 1) macroFile = argv[1];
  if (argc > 2) startingSeed = argv[2];
  if (argc > 3) outputFile = argv[3];

  // Print the arguments used for the executable.
  G4cout << "Starting run with... " << G4endl << G4endl;
  G4cout << "Macro file    : " << macroFile << G4endl << G4endl;
  G4cout << "Starting seed : " << startingSeed << G4endl << G4endl;
  G4cout << "Output file   : " << outputFile << G4endl << G4endl;

  // Initiate the multithreaded run manager if Geant4 is run in multithreaded
  //   mode. Otherwise, initiate the regular run manager.
  #ifdef G4MULTITHREADED

    G4MTRunManager* runManager = new G4MTRunManager;

  #else

    G4RunManager* runManager = new G4RunManager;

  #endif

  // Initiate the random engine. POSSIBLY UPDATE ???
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // // Convert the starting seed to an integer and feed it to the random engine.
  // unsigned startingSeedInt;
  // std::istringstream is(startingSeed);
  // is >> startingSeedInt;
  // G4Random::setTheSeed(startingSeedInt);

  // Instantiate the scoring manager. NECESSARY ???
  G4ScoringManager::GetScoringManager();

  // Instantiate the geometry.
  DetectorConstruction* det = new DetectorConstruction;
  runManager->SetUserInitialization(det);

  // Instantiate the physics list. POSSIBLY UPDATE ???
  runManager->SetUserInitialization(new QGSP_BERT_HP);

  // Set the action initialization.
  runManager->SetUserInitialization(new ActionInitialization(det, outputFile));

  // Initialize the visual and UI managers with or without the macro file.
  //  If initialized with the macro file, the simulation is run in interactive
  //  mode.
  //
  // Note: This gives possibility to import the geometry from different GDML
  //       files.
  #ifdef G4VIS_USE

    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

  #endif

  // Open an UI session: It will stay open until the user types "exit" or
  // simply uses a keyboard interrupt.
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  #ifdef G4UI_USE

  // Run the default macro file if no macro has been specified as an argument
  //  to the executable. The executable is set to use an UI and visualization
  //  is turned on.
  //
  //  aka., run in interactive mode
  if (argc == 1)
  {

      // Create a new UI executive object to start the simulation session if the
      //  executable has been run in interactive mode.
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    #ifdef G4VIS_USE

      // Run the visualization macro.
      UImanager->ApplyCommand("/control/execute vis.mac");

    #endif

      // Start the simulation session.
      ui->SessionStart();

      // Delete the UI executive object after the session is over.
      delete ui;
  }

  // Run the executable with the specified macro provided as an argument to the
  //   executable if one or more arguments have been provided by the user.
  //
  //  aka., run in batch mode
  else
  {
      // Create a new UI executive object to start the simulation session if
      //  the executable has been run in interactive mode.
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);

      // Define command prefix to use to execute the macro file.
      G4String command_prefix = "/control/execute ";

      // Apply the full command, including the macro file name.
      UImanager->ApplyCommand(command_prefix + macroFile);

      // Start the simulation session.
      ui->SessionStart();

      // Delete the UI executive object after the session is over.
      delete ui;
  }

  #endif

  // Terminate the visualization if it was used.
  #ifdef G4VIS_USE

    delete visManager;

  #endif

  // Terminate the simulation.
  delete runManager;

  // Exit the main function with no errors.
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
