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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  Detector Construction Class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4GDMLParser.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4MultiFunctionalDetector;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    // G4LogicalVolume* ReturnVolume() const {return fScoringVolume;}
    G4LogicalVolume* GetScoringVolume() const {return fScoringVolume;}

    std::vector<G4LogicalVolume*> GetScoringVolumeVec() const
    {
      return fScoringVolumeVec;
    }

    // Return the number of elements in the scoring volume vector.
    G4int GetScoreVolVecSize() const {return fScoringVolumeVec.size();}

    // Reading GDML
    //
    // void SetReadFile( const G4String& File);
    void SetReadFile( const G4String&);

    // Construct the sensitive detectors and electromagnetic fields through them,
    //  if desired.
    void ConstructSDandField();

  private:

    G4GDMLParser       fParser;
    G4Material*        fWorldMater;
    G4VPhysicalVolume* fPhysiWorld;
    G4int fWritingChoice;
    G4String fReadFile;
    DetectorMessenger* fDetectorMessenger;
    G4LogicalVolume* fScoringVolume;
    std::vector<G4LogicalVolume*> fScoringVolumeVec;

    void DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();

    // Use of G4Cache for multi-threaded operation. Allows for storage of
    // a seperate detector pointer per thread.
    const G4Cache<G4MultiFunctionalDetector*> fSensitiveDetectorCache;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif
