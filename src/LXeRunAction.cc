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
// $Id: LXeRunAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/LXe/src/LXeRunAction.cc
/// \brief Implementation of the LXeRunAction class
//
//
#include "LXeRunAction.hh"
#include "LXeRecorderBase.hh"
#include "G4SystemOfUnits.hh"

#include "LXeAnalysis.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeRunAction::LXeRunAction(LXeRecorderBase* r) : fRecorder(r) {}

/*LXeRunAction::LXeRunAction(LXeRecorderBase* r, G4String outputfile) : fRecorder(r) {
outputfilename = outputfile;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeRunAction::~LXeRunAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeRunAction::BeginOfRunAction(const G4Run* aRun){
  if(fRecorder)fRecorder->RecordBeginOfRun(aRun);
 // KL COMMENTED TO REMOVE WAVELENGTH HISTORGRAM
  G4AnalysisManager *AnalysisManE = G4AnalysisManager::Instance();

  AnalysisManE->OpenFile("newLAr");
  AnalysisManE->CreateNtuple("NewLAr", "dEdx");
  AnalysisManE->CreateNtupleDColumn("energy");
  AnalysisManE->CreateNtupleDColumn("dEdx");
  AnalysisManE->FinishNtuple();

  //AnalysisManE->CreateH1("0","-dE/dx of E", 100, 0., 1.);
//

/*  G4AnalysisManager *AnalysisMan = G4AnalysisManager::Instance();

  AnalysisMan->OpenFile("OutputData");

  AnalysisMan->CreateH1("0", "Hit X Position", 300, -0.15*m, 0.15*m);
  AnalysisMan->CreateH1("1", "Hit Y Position", 300, -0.15*m, 0.15*m);
  AnalysisMan->CreateH1("2", "Hit Z Position", 300, -0.15*m, 0.15*m);

  AnalysisMan->CreateNtuple("HitPosition", "Hits Position");
  AnalysisMan->CreateNtupleDColumn("Xpos");
  AnalysisMan->CreateNtupleDColumn("Ypos");
  AnalysisMan->CreateNtupleDColumn("Zpos");
  AnalysisMan->CreateNtupleDColumn("TrackID");
  AnalysisMan->FinishNtuple(); */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeRunAction::EndOfRunAction(const G4Run* aRun){
  if(fRecorder)fRecorder->RecordEndOfRun(aRun);

  G4AnalysisManager *AnalysisMan = G4AnalysisManager::Instance();
  AnalysisMan->Write();
  AnalysisMan->CloseFile();
}
