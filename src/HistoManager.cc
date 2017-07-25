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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 98265 2016-07-04 17:47:54Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "Run.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// HistoManager::HistoManager(DetectorConstruction* det)
  // , fDetector(det), fScoringVolume1(0.)
HistoManager::HistoManager()
: fFileName("HFNG_histograms")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // // Obtain the run manager.
  // Run* run = static_cast<Run*>(
  //       G4RunManager::GetRunManager()->GetNonConstCurrentRun());

  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     // Enable activation of histograms

  //get scoring volume
  //
  // const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>
        // (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // G4LogicalVolume* fScoringVolume1 = detectorConstruction->GetScoringVolume();

 // if (fScoringVolume1.empty()){fScoringVolume1 = detectorConstruction->GetScoringVolume1();}
//if (fScoringVolume1.empty()){fScoringVolume1 = fDetector->GetScoringVolume1();}
//G4cout<<"First scoring volume is: "<<fScoringVolume1[0]<<G4endl;
//fScoringVolume1 = fDetector->GetScoringVolume1();
//G4cout<<"Scoring volume size is: "<<fScoringVolume1.size()<<G4endl;
//getchar();

  // fScoringVolume1 = detectorConstruction->GetScoringVolume1();

  //for(int j = 0; j < 3; j++)
   // {
      // Retrieve the name of the current scoring volume.
  //  if (!fScoringVolume) {
   //  fScoringVolume = detectorConstruction->ReturnVolume();
  // }
   //   G4String volume_name = fScoringVolume1[j]->GetName();


      // Define an array to store the title of each histogram.
      const G4String title[] =
              {
                "Neutrons Leaving World"                          // ID = 0
                "Neutrons Entering Volume 1",                     // ID = 1
                "Neutrons Entering Volume 2",                     // ID = 2
                "Neutrons Scattering in Volume 1"                 // ID = 3
                "Neutrons Scattering in Volume 2"                 // ID = 4
                "Energy Deposition in Volume 1",                  // ID = 5
                "Energy Deposition in Volume 2",                  // ID = 6
                "Neutron Fluence Into Volume 1",                  // ID = 7
                "Neutron Fluence Into Volume 2",                  // ID = 8
              };

      // Define the maximum number of histograms as the number of elements in
      //  the title array. The numerator gives the total number of bytes used to
      //  hold the array in memory. The denominator gives the number of bytes
      //  used to hold the one element of the array in memory (the first).
      const G4int kMaxHisto = sizeof(title) / sizeof(*title);

      // Initialize an array to hold the histogram IDs.
      G4String id[kMaxHisto];

      // Fill the histogram ID array.
      for (int i=0; i<kMaxHisto; i++)
      {
        // // Create a temporary, dynamic pointer to hold the string version of
        // //  the histogram ID.
        // std::string* temp_ptr = new std::string;
        //
        // // Set the temporary pointer's dereference value equal to the string
        // //  version of the histogram ID currently chosen by the for loop
        // //  iterator.
        // *temp_ptr = std::to_string(i);
        //
        // // Dynamically cast the temporary pointer to a G4String reference
        // //  object and set the casted reference equal to the reference of the
        // //  ID array element currently chosen by the for loop iterator.
        // &id[i] = dynamic_cast<G4String*> (temp_ptr);

        id[i] = std::to_string(i);

        // Release the memory held by the temporary pointer.
        // delete temp_ptr;
      }
      // const G4String id[] = {"0", "1", "2","3", "4", "5", "6", "7", "8"};

      // Bin parameters
      G4int nbins = 500;    // For a set energy window, affects energy width
                            //   of each bin
      G4double x_min = 0.;  // Minimum value on the x-axis (MeV)
      G4double y_min = 5.;  // Minimum value on the y-axis (MeV)

      // Iterate through and create all histograms.
      for (G4int k=0; k<kMaxHisto; k++)
      {
        G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, x_min, y_min);
        analysisManager->SetH1Activation(ih, true);
      }
    //}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
