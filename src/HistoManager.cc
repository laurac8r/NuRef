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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("HFNG")
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
  // Obtain the run manager.
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 

  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     // Enable activation of histograms

  //get scoring volume
  //
  const DetectorConstruction* detectorConstruction = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  // G4LogicalVolume* fScoringVolume1 = detectorConstruction->GetScoringVolume();

  if (fScoringVolume1.empty()){fScoringVolume1 = detectorConstruction->GetScoringVolume1();}
  
  // fScoringVolume1 = detectorConstruction->GetScoringVolume1();

  for(int j = 0; j < fScoringVolume1.size(); j++)
    {
      // Retrieve the name of the current scoring volume.
      G4String volume_name = fScoringVolume1[j]->GetName();

      // Define histograms start values.
      const G4int kMaxHisto = 3;
      // const G4String id[] = {"0", "1", "2"};
      const G4String id[] = {std::to_string(j*10), std::to_string(j*10 + 1), 
                             std::to_string(j*10 + 2), std::to_string(j*10 + 3),
                             std::to_string(j*10 + 4)};
      const G4String title[] = 
                    {
                      "All Neutrons Entering the Volume " + volume_name,            // ID = 10*j + 0
                      "Neutrons Scattered in the Volume " + volume_name,            // ID = 10*j + 1
                      "Energy Absorbed in the Volume " + volume_name,               // ID = 10*j + 2
                      "Fission Reactions Occurring in the Volume " + volume_name,   // ID = 10*j + 3
                      "Capture Reactions Occurring in the Volume " + volume_name    // ID = 10*j + 4
                    };

      // Bin parameters               
      G4int nbins = 500;
      G4double x_min = 0.;  // Minimum value on the x-axis
      G4double x_max = 5.;  // Minimum value on the y-axis

      // Create all histograms
      for (G4int k=0; k<kMaxHisto; k++)
        {
          G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, x_min, x_max);
          analysisManager->SetH1Activation(ih, true);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
