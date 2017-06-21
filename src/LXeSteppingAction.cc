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
// $Id: LXeSteppingAction.cc 73915 2013-09-17 07:32:26Z gcosmo $
//
/// \file optical/LXe/src/LXeSteppingAction.cc
/// \brief Implementation of the LXeSteppingAction class
//
//
#include "LXeSteppingAction.hh"
#include "LXeEventAction.hh"
#include "LXeTrackingAction.hh"
#include "LXeTrajectory.hh"
#include "LXePMTSD.hh"
#include "LXeUserTrackInformation.hh"
#include "LXeUserEventInformation.hh"
#include "LXeSteppingMessenger.hh"
#include "LXeRecorderBase.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"

#include "LXeAnalysis.hh"

//static double RtotalStepLength = 0;
static double totalStepLength = 0;
static double TOTtotalEnergyDeposit = 0;
static double energyparticle = 0;
static int stepCounter = 1;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingAction::LXeSteppingAction(LXeRecorderBase* r)
  : fRecorder(r),fOneStepPrimaries(false)
{
  fSteppingMessenger = new LXeSteppingMessenger(this);

  fExpectedNextStatus = Undefined;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeSteppingAction::~LXeSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeSteppingAction::UserSteppingAction(const G4Step *theStep){
       G4AnalysisManager *AnalysisManE = G4AnalysisManager::Instance();
       G4Track* theTrack = theStep->GetTrack();
       //G4Track* theCurrentTrack = theStep->GetTrack();

  if (theStep->GetTrack()->GetParentID() == 0) {

    G4StepPoint *endPoint = theStep->GetPostStepPoint();
    G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();

    if(theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "scintillator"){
    //if(procName == "muBrem" || procName == "muIoni" || procName == "muPairProd"){

        //  fRunAction->CountProcesses(procName);
      G4StepPoint *startPoint = theStep->GetPreStepPoint();
      G4double E0 = startPoint->GetKineticEnergy();
      G4double E2 = endPoint->GetKineticEnergy();
      G4double totalEnergyDeposit = E0 - E2;
      energyparticle = (E0 + energyparticle*(stepCounter-1))/stepCounter;
      G4double stepLength = theStep->GetStepLength();
      totalStepLength = stepLength + totalStepLength;
      TOTtotalEnergyDeposit = totalEnergyDeposit + TOTtotalEnergyDeposit;
      G4double deDx = 0.;
      //deDx = (totalEnergyDeposit / MeV) * (cm / totalStepLength);
      deDx = (TOTtotalEnergyDeposit / MeV) * (cm / totalStepLength);
      //totalStepLength = stepLength + totalStepLength;
      //TOTtotalEnergyDeposit = totalEnergyDeposit + TOTtotalEnergyDeposit;
      stepCounter = stepCounter + 1;

      //if (totalEnergyDeposit > 0.)
      /*  deDx = (totalEnergyDeposit / MeV) * (cm / ItotalStepLength);
        meandeDx = (deDx + meandeDx*(stepCounter-1)) / stepCounter;
        stepCounter +=1 ;

             G4cout << " Radiative processes StepLength " << RtotalStepLength/cm <<G4endl;
             G4cout << " Ionisation StepLength " << ItotalStepLength/cm <<G4endl;
             G4cout << " stepLength " << stepLength/cm <<G4endl;
             G4cout << " deDx " << G4endl;
             G4cout << "    " << deDx << G4endl;
             G4cout << "    " << meandeDx  << G4endl;
             G4cout << "    " << meandeDx / 1.3954 << G4endl;
             G4cout << " meandeDx  " <<  G4endl;
             G4cout << " stepCounter  " << stepCounter << G4endl;
             G4cout << " **************" << G4endl; */
    //G4int id = 0;
    //BEA if ((procName == "muIoni") || (procName == "muPairProd") || (procName == "muBrem") || (procName == "muonNuclear"))
    //BEA  id = 1;
  //  fHistoManager->FillHisto(id, deDx);
      //if(totalStepLength > 8*cm){

    /*  G4cout << "Muon Energy: "<< energyparticle / MeV << G4endl;
      G4cout << "process name: "<< procName << G4endl;
      G4cout << "Current step length " << stepLength / cm << G4endl;
      G4cout << "Total step length " << totalStepLength / cm << G4endl;
      G4cout << "Current energy deposit " << totalEnergyDeposit / MeV << G4endl;
      G4cout << "Total energy deposit " << TOTtotalEnergyDeposit / MeV << G4endl;
      G4cout << "deDx " << deDx << G4endl;
      G4cout << "step count " << stepCounter << G4endl;
      G4cout << "-----" << G4endl;*/

      if(totalStepLength > 1*m){
      if(TOTtotalEnergyDeposit > 0.0 && totalStepLength > 0.0) {

        AnalysisManE->FillNtupleDColumn(0, energyparticle / MeV);
        AnalysisManE->FillNtupleDColumn(1, deDx);
        AnalysisManE->AddNtupleRow();

        totalStepLength = 0;
        TOTtotalEnergyDeposit = 0;
        energyparticle = 0;
        stepCounter = 1;
        theTrack->SetTrackStatus(fStopAndKill);
        G4cout << "Muon dead" << G4endl;
        G4cout << "-----" << G4endl;


      }}
        //  RtotalStepLength = 0;
}}
/*

    if(theTrack->GetParentID()==0){
       if(theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "scintillator"){
         if(procName == "muIoni"){
            G4cout << "-----" << G4endl;
            G4cout << "Particle's definition: " << theCurrentTrack->GetParticleDefinition()->GetParticleName() << G4endl;
            //if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Scintillation"){
            //  G4cout << "Process: " << theCurrentTrack-> GetCreatorProcess()->GetProcessName() << G4endl;
            //}
            G4cout << "Muon Energy: "<< (E0 + E2) / 2.0 << G4endl;
            G4cout << "E0-E2: " << E0 - E2<< G4endl;
            G4cout << "Delta Energy: " << theStep->GetDeltaEnergy() << G4endl;
            G4cout << "Step Length: " << stepLength << G4endl;
            G4cout << "Step Length*: " << theCurrentTrack->GetStepLength() << G4endl;
            G4cout << "Energy Deposited to Medium: " << theStep->GetTotalEnergyDeposit() << G4endl;
            G4cout << "-----" << G4endl;
          }
}
}*/




//////////////////////////////////////////////////////////////////////////////


// -------------------- ADD MUON ENERGY TO OUTPUTDATA.ROOT --------------------
/*
   if(theTrack->GetParentID()==0){
	  if(theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName() == "scintillator"){
         G4cout << "-----" << G4endl;
         G4cout << "Particle's definition: " << theCurrentTrack->GetParticleDefinition()->GetParticleName() << G4endl;
         G4cout << "Muon Energy: "<< (E0 + E2) / 2.0 << G4endl;
         G4cout << "E0-E2: " << E0 - E2<< G4endl;
         G4cout << "Delta Energy: " << theStep->GetDeltaEnergy() << G4endl;
         G4cout << "Step Length: " << stepLength << G4endl;
         G4cout << "Step Length*: " << theCurrentTrack->GetStepLength() << G4endl;
         G4cout << "Energy Deposited to Medium: " << theStep->GetTotalEnergyDeposit() << G4endl;
         G4cout << "-----" << G4endl;
         //AnalysisManE->FillH1(0,theStep->GetTotalEnergyDeposit()/theCurrentTrack->GetStepLength());
        // AnalysisManE->FillNtupleDColumn(0, theStep->GetPreStepPoint()->GetKineticEnergy()/MeV);
         //AnalysisManE->FillNtupleDColumn(1, theStep->GetDeltaEnergy()/MeV);
         //AnalysisManE->FillNtupleDColumn(2, theStep->GetStepLength()/cm);
         //AnalysisManE->AddNtupleRow();
         //if(theStep->GetDeltaEnergy()!=0){
           //if(theStep->GetPreStepPoint()->GetKineticEnergy()!=0){
            // G4cout << "E = " << theStep->GetPreStepPoint()->GetKineticEnergy() << G4endl;
             //G4cout << "EDep/stepLen = " << theStep->GetDeltaEnergy()/theStep->GetStepLength() << G4endl;
	   //  G4cout << "stepLen = " << theStep->GetStepLength() << G4endl;
      //       G4cout << "----------" << G4endl;
           //}
         //}
        }
      }*/


// -------------------- PRINT TO SCREEN  --------------------
/*       if (theTrack->GetCurrentStepNumber()==1){
         if(theTrack->GetParentID()==1){
          if(theCurrentTrack->GetCreatorProcess()->GetProcessName()!="Scintillation"){
             if(theCurrentTrack->GetCreatorProcess()->GetProcessName()!="Cerenkov"){
//               if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="muPairProd"){
                 G4cout << "-----" << G4endl;
                 G4cout << "Parent Particle's ID: " << theTrack->GetParentID() << G4endl;
                 G4cout << "Particle's definition: " << theCurrentTrack->GetParticleDefinition()->GetParticleName() << G4endl;
                 G4cout << "Process name: " <<  theCurrentTrack->GetCreatorProcess()->GetProcessName() << G4endl;
                 G4StepPoint* thePreStepPoint = theStep->GetPreStepPoint();
                 G4cout << "Prestep Point: " << thePreStepPoint->GetPosition().getX()<<","<<thePreStepPoint->GetPosition().getY()<<","<<thePreStepPoint->GetPosition().getZ() << G4endl;
                 G4cout << "-----" << G4endl;
               }
             }
           }


// ----------------------- COUNT SCINT AND CEREN PHOTONS -----------------------
           if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Scintillation"){
              scintphoton += 1;
              if (scintphoton % 10000 == 0){
                 G4cout << " ==> Number of Scintillated  Photons: " << scintphoton << G4endl;
              }
           }
           if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov"){
              cerenphoton += 1;
              if (cerenphoton % 10000 == 0){
                G4cout << " ==> Number of Cerenkov  Photons: " << cerenphoton << G4endl;
              }
            }
         }
	}
*/
// ----------  WHEN SECONDARY IS CREATED NOT SCINT OR CEREN ----------



/*       if (theTrack->GetCurrentStepNumber()==1){
         if(theTrack->GetParentID()==1){
           G4cout << "Process name: " <<  theCurrentTrack->GetCreatorProcess()->GetProcessName() << G4endl;
           G4cout << "Particles definition: " <<  theCurrentTrack->GetParticleDefinition()->GetParticleName() << G4endl;

           if (theCurrentTrack->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
             createdphoton += 1;
             G4cout << "Number of created photons: " <<  createdphoton << G4endl;
             if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Scintillation"){
                   scintphoton += 1;
                   //G4cout << "Process: " << theCurrentTrack-> GetCreatorProcess()->GetProcessName() << G4endl;
                   G4cout << "Number of scintillated photons: " <<  scintphoton << G4endl;
               }
             if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Cerenkov"){
                cerenphoton += 1;
                G4cout << "Number of cerenkov photons: " <<  cerenphoton << G4endl;
              }
            }
          }
        }
*/

  if ( theTrack->GetCurrentStepNumber() == 1 ) fExpectedNextStatus = Undefined;

/*KL PRINT ENERGY TO SCREEN
  if (theTrack->GetCurrentStepNumber() == 1){
    if(theStep->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
      // G4cout << theStep->GetTrack()->GetKineticEnergy() << G4endl;
    }
  }
*/ //KL

 //KL SAVE WAVELENGTHS TO HISTORGRAM
 //G4AnalysisManager *AnalysisManE = G4AnalysisManager::Instance();
/*  double theWavelength;
    //G4Track* theCurrentTrack = theStep->GetTrack();
    G4StepPoint* thePreStepPoint = theStep->GetPreStepPoint();
      if (theCurrentTrack->GetCurrentStepNumber() == 1) {
        if (theCurrentTrack->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
          if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Scintillation"){
            if(theTrack->GetParentID()==1){
              theWavelength = 1.2398/(theCurrentTrack->GetTotalEnergy()/keV);
              //G4cout << "-----" << G4endl;
              //G4cout << "Prestep Point: " << thePreStepPoint->GetPosition().getX()<<","<<thePreStepPoint->GetPosition().getY()<<","<<thePreStepPoint->GetPosition().getZ() << G4endl;
              //G4cout << "Position: " << theCurrentTrack->GetPosition().getX()<<","<<theCurrentTrack->GetPosition().getY()<<","<<theCurrentTrack->GetPosition().getZ() << G4endl;
              //G4cout << "Wavelength: " << theWavelength << G4endl;
              //G4cout << "Creator Process: " << theCurrentTrack->GetCreatorProcess()->GetProcessName() << G4endl;
              //AnalysisManE->FillH1(0,theWavelength);
            }
            //else{
              //G4cout << "Parent ID: " <<  theTrack->GetParentID() << G4endl;
            //}
          }
        }
      } */
    //G4cout << "Length: " << step->GetStepLength()/mm << " (mm). ";
    return;
 //KL

  LXeUserTrackInformation* trackInformation
    =(LXeUserTrackInformation*)theTrack->GetUserInformation();
  LXeUserEventInformation* eventInformation
    =(LXeUserEventInformation*)G4EventManager::GetEventManager()
    ->GetConstCurrentEvent()->GetUserInformation();

  G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm
      = theStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

  if(theTrack->GetParentID()==0){
    //This is a primary track

    G4TrackVector* fSecondary=fpSteppingManager->GetfSecondary();
    G4int tN2ndariesTot = fpSteppingManager->GetfN2ndariesAtRestDoIt()
      + fpSteppingManager->GetfN2ndariesAlongStepDoIt()
      + fpSteppingManager->GetfN2ndariesPostStepDoIt();

    //If we havent already found the conversion position and there were
    //secondaries generated, then search for it
    if(!eventInformation->IsConvPosSet() && tN2ndariesTot>0 ){
      for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
          lp1<(*fSecondary).size(); lp1++){
        const G4VProcess* creator=(*fSecondary)[lp1]->GetCreatorProcess();
        if(creator){
          G4String creatorName=creator->GetProcessName();
          if(creatorName=="phot"||creatorName=="compt"||creatorName=="conv"){
            //since this is happening before the secondary is being tracked
            //the Vertex position has not been set yet(set in initial step)
            eventInformation->SetConvPos((*fSecondary)[lp1]->GetPosition());
          }
        }
      }
    }

    if(fOneStepPrimaries&&thePrePV->GetName()=="scintillator")
      theTrack->SetTrackStatus(fStopAndKill);
  }

  if(!thePostPV){//out of world
    fExpectedNextStatus=Undefined;
    return;
  }
//tracking photon
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
    //Optical photon only

    //if(thePrePV->GetName()=="Slab")
      //force drawing of photons in WLS slab
      //trackInformation->SetForceDrawTrajectory(true);
    /*else if(thePostPV->GetName()=="expHall")
      //Kill photons entering expHall from something other than Slab
      theTrack->SetTrackStatus(fStopAndKill);
      */

    //Was the photon absorbed by the absorption process
    if(thePostPoint->GetProcessDefinedStep()->GetProcessName()
       =="OpAbsorption"){
      eventInformation->IncAbsorption();
      trackInformation->AddTrackStatusFlag(absorbed);
    }

    boundaryStatus=boundary->GetStatus();

    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check
    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      if(fExpectedNextStatus==StepTooSmall){
        if(boundaryStatus!=StepTooSmall){
          G4ExceptionDescription ed;
          ed << "LXeSteppingAction::UserSteppingAction(): "
                << "No reallocation step after reflection!"
                << G4endl;
          G4Exception("LXeSteppingAction::UserSteppingAction()", "LXeExpl01",
          FatalException,ed,
          "Something is wrong with the surface normal or geometry");
        }
      }
      fExpectedNextStatus=Undefined;
      switch(boundaryStatus){
      case Absorption:
      {
        //FillHistogram(theStep);
        FillNtuple(theStep);
        trackInformation->AddTrackStatusFlag(boundaryAbsorbed);
        eventInformation->IncBoundaryAbsorption();
        break;
      }
      case Detection: //Note, this assumes that the volume causing detection
                      //is the photocathode because it is the only one with
                      //non-zero efficiency
        {
          //FillHistogram(theStep);
          FillNtuple(theStep);
        //Triger sensitive detector manually since photon is
        //absorbed but status was Detection
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdName="/LXeDet/pmtSD";
        LXePMTSD* pmtSD = (LXePMTSD*)SDman->FindSensitiveDetector(sdName);
        if(pmtSD)pmtSD->ProcessHits_constStep(theStep,NULL);
        trackInformation->AddTrackStatusFlag(hitPMT);
        break;
        }
      case FresnelReflection: FillNtuple(theStep);
      case TotalInternalReflection: FillNtuple(theStep);
      case LambertianReflection: FillNtuple(theStep);
      case LobeReflection: FillNtuple(theStep);
      case SpikeReflection: FillNtuple(theStep);
      case BackScattering:
      {
        //FillHistogram(theStep);
        FillNtuple(theStep);
        trackInformation->IncReflections();
        fExpectedNextStatus=StepTooSmall;
        break;
      }
      default:
        break;
      }
      //if(thePostPV->GetName()=="sphere")
        //trackInformation->AddTrackStatusFlag(hitSphere);
    }
  }

  if(fRecorder)fRecorder->RecordStep(theStep);
}
 // KL
/*void LXeSteppingAction::FillHistogram(const G4Step *theStep) {
  G4AnalysisManager *AnalysisMan = G4AnalysisManager::Instance();
// KL ONLY COLLECT SCINTILLATION
  G4Track* theCurrentTrack = theStep->GetTrack();
  G4Track* theTrack = theStep->GetTrack();
  if (theCurrentTrack->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Scintillation"){
      if(theTrack->GetParentID()==1){
        AnalysisMan->FillH1(0, theTrack->GetPosition().getX());
        AnalysisMan->FillH1(1, theTrack->GetPosition().getY());
        AnalysisMan->FillH1(2, theTrack->GetPosition().getZ());
      }
    }
  }
}

void LXeSteppingAction::FillNtuple(const G4Step *theStep) {
  G4AnalysisManager *AnalysisMan = G4AnalysisManager::Instance();
  G4Track* theTrack = theStep->GetTrack();
  G4Track* theCurrentTrack = theStep->GetTrack();
  if (theCurrentTrack->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    if(theCurrentTrack->GetCreatorProcess()->GetProcessName()=="Scintillation"){
      if(theTrack->GetParentID()==1){
        AnalysisMan->FillNtupleDColumn(0, theTrack->GetPosition().getX());
        AnalysisMan->FillNtupleDColumn(1, theTrack->GetPosition().getY());
        AnalysisMan->FillNtupleDColumn(2, theTrack->GetPosition().getZ());
        AnalysisMan->FillNtupleDColumn(3, theTrack->GetTrackID());
        AnalysisMan->AddNtupleRow();
      }
    }
  }
}*/
//
