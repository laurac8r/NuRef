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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 76293 2013-11-08 13:11:23Z gcosmo $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "RunAction.hh"

#include "G4ios.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction(const G4String& outputFile))
: G4UserEventAction(),
  fRunAction(runAction), fEngDep(0.), fEngDepArr{0.}
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction()
{
  // Initialize the accumulated energy deposition variabl.
  // fEngDep = fEngDepArr[0] = fEngDepArr[1] = 0.;
  fEngDep 0.;
  fEngDepArr = {0.};

  G4cout << G4endl << G4endl << "The accumulated energy deposition vector is "
         << fEngDepArr << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{

  // Obtain the event ID so we can add it's energy deposition to the
  //   accumulated energy deposition variables.
  G4int eventID = event->GetEventID();

  // Print energy deposited in each event.
  fRunAction->AddEngDep(fEngDep);

  // If the accumulated energy deposition is nonzero for the single logical
  //   scoring volume at the end of the event action, print it.
  if (fEngDep > 0)
  {
    G4cout << G4endl << G4endl << "Total energy deposition for this event "
           << eventID << " in the single logical scoring volume is:"
           << G4BestUnit(fEngDep, "Energy") << G4endl;
  }

  // Using vector container
  fRunAction->AddEngDepArr(fEngDepArr);

  // Print energy deposited in each event
  for(int i = 0; i<2; i++)
   {
    if(fEngDepArr[i]>0){
    G4cout<<"Total Edep for this event "<< eventID<<" for logical volume having the ID "<<i<<" is:"<<G4BestUnit(fEngDepArr[i],"Energy")<<G4endl;}
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
