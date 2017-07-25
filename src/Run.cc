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
/// \file Run.cc
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $

#include "Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"

#include "G4SystemOfUnits.hh"
#include <assert.h>

// #include "DetectorConstruction.hh"
// #include "PrimaryGeneratorAction.hh"
// #include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run()
: G4Run(), fluenceMap()
{
  // Create a hits map for each primitive scorer.
  fluenceMap[0] = new G4THitsMap<G4double>("Scoring Volume",
                                           "Neutron Cell Fluence [cm^-2]");
  // // The following is just an example of how to add another scorer.
  // fluenceMap[1] = new G4THitsMap<G4double>("Scoring Volume",
  //                                          "Neutron Surface Fluence [cm^-2]")
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run()
{
  // Create an iterator to clean up the fluence map.
  std::map< G4int, G4THitsMap<G4double>* >::iterator iter = fluenceMap.begin();

  // Clean up the fluence map.
  while (iter != fluenceMap.end())
  {
    delete iter->second;
    iter++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* anEvent)
{
  // Get the hits collection.
  G4HCofThisEvent* eventHitCollection = anEvent->GetHCofThisEvent();

  // Exit this RecordEvent instance if there is no event hit collection.
  if (!eventHitCollection) return;

  // Create an iterator for updating the fluence map.
  std::map< G4int, G4THitsMap<G4double>* >::iterator iter = fluenceMap.begin();

  // Iterate through and update the fluence map.
  while (iter != fluenceMap.end())
  {
    // Get the ID of the fluence entry in the map picked by the iterator.
    G4int id = iter->first;

    // Get the hit collection corresponding to the fluence entry ID.
    G4THitsMap<G4double>* eventHitsMap =
      dynamic_cast< G4THitsMap<G4double>* >(eventHitCollection->GetHC(id));

    // Expect the event hits map to exist. Non-existant if dynamic cast fails.
    //   Used as a fail-safe.
    assert (0 != eventHitsMap);

    // Accumulate (add) event data into our fluence map.
    *(iter->second) += *eventHitsMap;

    // Increment the iterator.
    iter++;
  }

  // Record the event.
  G4Run::RecordEvent(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::DumpData(G4String &outputFileSpec) const
{
  // Create a vector for the titles of the hit maps.
  std::vector<G4String> titles;
  titles.push_back("Length");

  // Output fluence map - energy binning on x-axis, theta on y.
  std::map< G4int, std::vector<G4double> > output;

  // Define the number of energy bins based on the number of scorers.
  G4int numEnergyBins = fluenceMap.size();

  // Initialize the theta and energy bin iterators.
  G4int i(0), j(0);

  // Initialize fluence to 0 in all bins of the output map. Uses theta and energy
  //   iteration.
  for (i=0; i < numThetaBins; i++)
  {
    for (j=0; j < numEnergyBins; j++)
    {
      output[i].push_back(0);
    }
  }

  // Reinitialize the energy bin iterator to zero.
  j = 0;

  // Create an iterator to fill the output fluence map.
  std::map< G4int, G4THitsMap<G4double>* >::const_iterator iter = fluenceMap.begin();

  // Iterate through and fill the output fluence map.
  while (iter != fluenceMap.end())
  {
    // Call the hit map value from the output fluence map.
    G4THitsMap<G4double>* hitMap = iter->second;

    // Add the hit map name to the titles vector.
    titles.push_back(hitMap->GetName());

    // Create a pointer to the hit map value pointer for dereferencing.
    std::map< G4int, G4double* >* myMap = hitMap->GetMap();

    // Iterate through the theta and energy bins and fill the output fluence map.
    for (i=0; i < numThetaBins; i++)
    {
      // Create a fluence pointer to the theta- and energy-binned hit map pointer
      //   selected by the iterator.
      G4double* fluence = (*myMap)[i];

      // Fill the selected element of the output fluence map only if the fluence
      //   pointer is not null. Use the dereference of the fluence pointer.
      if (0 != fluence) output[i][j] = (*fluence);
    }

    // Increment the energy bin and output fluence map iterators.
    j++;
    iter++;
  }

  Print(titles, output, outputFileSpec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Print(const std::vector<G4String>& titles,
                const std::map< G4int, std::vector<G4double> >& myMap,
                G4String &outputFileSpec) const
{
  // Create file to output the fluence map to.
  std::ofstream outFile(outputFileSpec);

  // Create a constant iterator for the titles vector.
  std::vector<G4String>::const_iterator titlesIter = titles.begin();

  // Iterate through the titles vector and print each title to the screen.
  while (titlesIter != titles.end())
  {
    // Print to the screen with proper spacing in between each title using the
    //   dereference of the titles iterator.
    G4cout << std::setw(8) << *titlesIter << " ";

    // Increment the titles iterator.
    titlesIter++;
  }

  // Create a blank line to seperate the titles from the fluence data.
  G4cout << G4endl << G4endl;

  // Create a constant iterator for the fluence map data, as well as for each
  //   scorer.
  std::map< G4int, std::vector<G4double> >::const_iterator iter = myMap.begin();

  // Iterate through fluence map and print out the fluence data.
  while (iter != myMap.end())
  {
    // Print out the integer ID number of the scorer of interest to the screen.
    G4cout << std::setw(8) << std::setprecision(3) << iter->first << " ";

    // Create a constant iterator for the energy bins in the fluence data for
    //   the scorer of interest.
    std::vector<G4double>::const_iterator energyBinIter = iter->second.begin();

    // Iterator through the energy bins in the fluence map and output the
    //   fluence data.
    while (energyBinIter != iter->second.end())
    {
      // Extract the fluence value using the dereference of the energy bin
      //   iterator.
      G4double value = *energyBinIter;

      // Print the fluence values to the screen.
      G4cout << std::setw(10) << std::setprecision(5) << value * cm * cm << " ";

      // Print the fluence values to the output file.
      outFile << value * cm * cm;

      // Icrement the energy bin iterator.
      energyBinIter++;
    }

    // Append a blank line to the output file and the screen.
    outFile << G4endl;
    G4cout << G4endl;

    // Increment the fluence map iterator.
    iter++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* run)
{
  // Cast the G4Run object to the Run object.
  const Run* localRun = static_cast<const Run*>(run);

  // Create a local map reference object using the pointer to the fluence map.
  const std::map< G4int, G4THitsMap<G4double>* >& localMap = localRun->fluenceMap;

  // Create a constant iterator for the local fluence map.
  std::map< G4int, G4THitsMap<G4double>* >::const_iterator iter = localMap.begin();

  // Iterate through and merge all the event fluence data into the run fluence
  //   data using the dereference of the element of the run fluence map
  //   (corresponding to the iterator's first element at the iterator position)
  //   to equal to the dereference of the iterator's second element.
  for( /* Iterator already initialized */ ; iter != localMap.end(); ++iter )

      ( *(fluenceMap[iter->first]) ) += ( *(iter->second) );


  // Eneable Geant4 to output the overall summary information for the run.
  G4Run::Merge(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
