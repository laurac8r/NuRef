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
/// \file RunAction.hh
/// \brief Definition of the RunAction class
//
// $Id: RunAction.hh 66241 2012-12-13 18:34:42Z gunter $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Accumulable.hh"

#include <vector>

#include "Run.hh"

// class Run;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RunAction : public G4UserRunAction
{
  public:
    // Define the constructor for the run action.
    RunAction(const G4String&);

    // Define the destructor for the run action.
    virtual ~RunAction();

    // Define a method to generate a new run.
    virtual G4Run* GenerateRun() {return new Run;};

    // Define the beginning of the run action, in which initialization of the
    //   energy deposition accumulables and the histogram analysis manager
    //   occurs.
    virtual void BeginOfRunAction(const G4Run*);

    // Define the end of the run action, in which saving of the ROOT file,
    //   dumping of the fluence data, merging of the energy deposition
    //   accumulables, and calculation of the root-mean-squared (RMS) energy
    //   deposition occurs.
    virtual void EndOfRunAction(const G4Run*);

    // Define a method to increment the neutron capture counter.
    void IncrementNeutCap() {fNeutCap++;}

    void AddEngDep (G4double&);
    void AddEngDepArr (G4double (&engDep)[2]);

  private:
    // Create ordinary accumulables for energy deposition for the case of only
    //   one scoring logical volume.
    G4Accumulable<G4double> fEngDep;
    G4Accumulable<G4double> fEngDepSqr;

    // Define vectors to store the accumulables for the case of multiple
    //   scoring logical volumes. CHANGE THE DEFAULT NUMBER OF SCORING
    //   VOLUMES TO THAT DEFINED BY SOME CUSTOM METHOD CALLING A PRIVATE
    //   VARIABLE IN THE DETECTOR CONSTRUCTION CLASS.
    std::vector< G4Accumulable<G4double> > fEngDepArr;
    std::vector< G4Accumulable<G4double> > fEngDepSqrArr;

    // Define a double-precision floating-point number to store the
    //  accumulated number of neutron captures.
    G4int fNeutCap;

    // // Reserve the number of scoring logical volumes for memory for the energy
    // //   deposit arrays.
    // fEngDepArr.reserve(2);
    // fEngDepSqrArr.reserve(2);

    G4String fOutputFileSpec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
