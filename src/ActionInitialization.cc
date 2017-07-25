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
// $Id: ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "globals.hh"

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "HistoManager.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"

// #include "TrackingAction.hh"
// #include "SteppingVerbose.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::ActionInitialization(DetectorConstruction* detector,
                                           const G4String& outputFile)
 : G4VUserActionInitialization(),
   fDetector(detector), fOutputFileSpec(outputFile)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitialization::~ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction(fOutputFileSpec));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  // Define and set a new run action.
  RunAction* runAction = new RunAction(fOutputFileSpec);
  SetUserAction(runAction);

  // Define and set a new primary generator action.
  SetUserAction(new PrimaryGeneratorAction());

  // Define and set a new event action.
  EventAction* eventAction = new EventAction(runAction);
  SetUserAction(eventAction);

  // Define a new histogram manager. The setup of new histograms is done with this
  //  HistoManager class instantiation.
  // HistoManager* histoManager = new HistoManager();

  // SetUserAction(new TrackingAction(fDetector, event, histo));

  // Define and set a new stepping action.
  SetUserAction(new SteppingAction(runAction, eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4VSteppingVerbose* ActionInitialization::InitializeSteppingVerbose() const
// {
//   return new SteppingVerbose();
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
