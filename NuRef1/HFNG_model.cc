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
#include "Randomize.hh"

#include "DetectorConstruction.hh"
//#include "PhysicsList.hh"
#include "QGSP_BERT.hh"
#include "Shielding.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

#ifdef G4VIS_USE
 #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {

  G4RunManager* runManager = new G4RunManager;

  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // set mandatory initialization classes
  DetectorConstruction* det= new DetectorConstruction;
  runManager->SetUserInitialization(det);

  G4VUserPhysicsList* phyList = new QGSP_BERT;
  runManager->SetUserInitialization(phyList);

  runManager->SetUserInitialization(new ActionInitialization(det));

  // Initialisation of runManager via macro for the interactive mode
  // This gives possibility to give different names for GDML file to READ
 
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

 // Open a UI session: will stay there until the user types "exit"
  //
  G4UImanager* UImanager = G4UImanager::GetUIpointer(); 

  if ( argc==1 )   // Automatically run default macro for writing... 
{
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
    ui->SessionStart();
    delete ui;
#endif
  } else {            // Provides macro in input

#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4String command = "/control/execute "; 
    G4String fileName = argv[1]; 
    UImanager->ApplyCommand(command+fileName); 
    ui->SessionStart();
      ui->SessionStart();
      delete ui;
#endif
}

#ifdef G4VIS_USE
     delete runManager;
#endif

  // job termination
  //
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
