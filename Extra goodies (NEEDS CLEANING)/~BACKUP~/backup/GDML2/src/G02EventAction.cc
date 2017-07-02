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
// $Id: G02EventAction.cc 75214 2013-10-29 16:04:42Z gcosmo $
//
/// \file G02EventAction.cc
/// \brief Implementation of the G02EventAction class

#include "G02EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G02RunAction.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G02EventAction::G02EventAction(G02RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction),
  fEdep(0.),fEdep1{0.}
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G02EventAction::~G02EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G02EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
fEdep1[0]=fEdep1[1] = {0.};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G02EventAction::EndOfEventAction(const G4Event* event)
{
G4int eventID = event->GetEventID();
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);
  //print energy deposited in each event
/*
if(fEdep>0){ 
G4cout<<"Total Edep for this event "<< eventID<<" is:"<<G4BestUnit(fEdep,"Energy")<<G4endl;}
*/
//using vector container

//fRunAction->AddEdep1(fEdep1);
  //print energy deposited in each event
for(int i = 0;i<2;i++)
{
if(fEdep1[i]>0){ 
G4cout<<"Total Edep for this event "<< eventID<<" for logical volume having the ID "<<i<<" is:"<<G4BestUnit(fEdep1[i],"Energy")<<G4endl;}
}  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
