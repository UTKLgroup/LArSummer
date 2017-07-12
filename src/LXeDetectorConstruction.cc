#include "LXeDetectorConstruction.hh"
#include "LXePMTSD.hh"
#include "LXeScintSD.hh"
#include "LXeDetectorMessenger.hh"
#include "LXeMainVolume.hh"
#include "LXeWLSSlab.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
//G4bool LXeDetectorConstruction::fSphereOn = true;

LXeDetectorConstruction::LXeDetectorConstruction() : fLAr_mt(NULL), fMPTPStyrene(NULL) {
  fExperimentalHall_box = NULL;
  fExperimentalHall_log = NULL;
  fExperimentalHall_phys = NULL;

  fLXe = fAl = fAir = fVacuum = fGlass = NULL;
  fPstyrene = fPMMA = fPethylene1 = fPethylene2 = NULL;

  fN = fO = fC = fH = NULL;

  SetDefaults();

  fDetectorMessenger = new LXeDetectorMessenger(this);
  fUserLimits = new G4UserLimits();
}

LXeDetectorConstruction::~LXeDetectorConstruction() {}

void LXeDetectorConstruction::DefineMaterials() {
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4int polyPMMA = 1;
  G4int nC_PMMA = 3+2*polyPMMA;
  G4int nH_PMMA = 6+2*polyPMMA;

  G4int polyeth = 1;
  G4int nC_eth = 2*polyeth;
  G4int nH_eth = 4*polyeth;

  //***Elements
  fH = new G4Element("H", "H", z=1., a=1.01*g/mole);
  fC = new G4Element("C", "C", z=6., a=12.01*g/mole);
  fN = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  fO = new G4Element("O"  , "O", z=8., a= 16.00*g/mole);

  //***Materials
  // Liquid Xenon
  fLXe = new G4Material("LXe",z=54.,a=131.29*g/mole,density=3.020*g/cm3);
  // Liquid Argon
  fLAr = new G4Material("LAr",z=18.,a=39.948*g/mole,density=1.3954*g/cm3);
  //Hydrogen
  //fH = new G4Material("Hi",z=1.,a=1.01*g/mole,density=0.000083748*g/cm3);

  //Aluminum
  fAl = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.6989*g/cm3);
  //Vacuum
  fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
                          density=universe_mean_density,kStateGas,0.1*kelvin,
                          1.e-19*pascal);
  //Air//changed to Hydrogen
  fAir = new G4Material("Air", density= 1.29*mg/cm3, 2);
  fAir->AddElement(fN, 70*perCent);
  fAir->AddElement(fO, 30*perCent);
  //Glass
  fGlass = new G4Material("Glass", density=1.032*g/cm3,2);
  fGlass->AddElement(fC,91.533*perCent);
  fGlass->AddElement(fH,8.467*perCent);
  //Polystyrene
  fPstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  fPstyrene->AddElement(fC, 8);
  fPstyrene->AddElement(fH, 8);
  //Fiber(PMMA)
  fPMMA = new G4Material("PMMA", density=1190*kg/m3,3);
  fPMMA->AddElement(fH,nH_PMMA);
  fPMMA->AddElement(fC,nC_PMMA);
  fPMMA->AddElement(fO,2);
  //Cladding(polyethylene)
  fPethylene1 = new G4Material("Pethylene1", density=1200*kg/m3,2);
  fPethylene1->AddElement(fH,nH_eth);
  fPethylene1->AddElement(fC,nC_eth);
  //Double cladding(flourinated polyethylene)
  fPethylene2 = new G4Material("Pethylene2", density=1400*kg/m3,2);
  fPethylene2->AddElement(fH,nH_eth);
  fPethylene2->AddElement(fC,nC_eth);

  //***Material properties tables
  // LXe
  G4double lxe_Energy[]    = { 7.0*eV , 7.07*eV, 7.14*eV };
  const G4int lxenum = sizeof(lxe_Energy)/sizeof(G4double);
  G4double lxe_SCINT[] = { 0.1, 1.0, 0.1 };
  assert(sizeof(lxe_SCINT) == sizeof(lxe_Energy));
  G4double lxe_RIND[]  = { 1.59 , 1.57, 1.54 };
  assert(sizeof(lxe_RIND) == sizeof(lxe_Energy));
  G4double lxe_ABSL[]  = { 35.*cm, 35.*cm, 35.*cm};
  assert(sizeof(lxe_ABSL) == sizeof(lxe_Energy));
  fLXe_mt = new G4MaterialPropertiesTable();
  fLXe_mt->AddProperty("FASTCOMPONENT", lxe_Energy, lxe_SCINT, lxenum);
  fLXe_mt->AddProperty("SLOWCOMPONENT", lxe_Energy, lxe_SCINT, lxenum);
  fLXe_mt->AddProperty("RINDEX",        lxe_Energy, lxe_RIND,  lxenum);
  fLXe_mt->AddProperty("ABSLENGTH",     lxe_Energy, lxe_ABSL,  lxenum);
  fLXe_mt->AddConstProperty("SCINTILLATIONYIELD",12000./MeV);
  fLXe_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  fLXe_mt->AddConstProperty("FASTTIMECONSTANT",20.*ns);
  fLXe_mt->AddConstProperty("SLOWTIMECONSTANT",45.*ns);
  fLXe_mt->AddConstProperty("YIELDRATIO",1.0);
  fLXe->SetMaterialPropertiesTable(fLXe_mt);
  // Set the Birks Constant for the LXe scintillator
  fLXe->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // LAr
  // Data taken from LArSOFT
  G4double LAr_ScintEnergy_fast[] =
//KL   {6.0*eV, 6.7*eV, 7.1*eV, 7.4*eV, 7.7*eV, 7.9*eV, 8.1*eV, 8.4*eV,
//KL    8.5*eV, 8.6*eV, 8.8*eV, 9.0*eV, 9.1*eV, 9.4*eV, 9.8*eV, 10.4*eV,
//KL    10.7*eV};
  {8.95*eV,9.14*eV,9.27*eV,9.36*eV,9.43*eV,9.48*eV,9.52*eV,9.57*eV,9.61*eV,
  9.64*eV,9.67*eV,9.71*eV,9.75*eV,9.81*eV,9.93*eV,9.98*eV,10.03*eV,10.05*eV,
  10.07*eV,10.11*eV,10.13*eV,10.16*eV,10.20*eV,10.23*eV,10.28*eV,10.33*eV,
  10.41*eV,10.52*eV,10.64*eV,10.84*eV};
  const G4int LArScintNum = sizeof(LAr_ScintEnergy_fast)/sizeof(G4double);
  G4double LAr_Scint_fast[] =
  {0.02,0.06,0.12,0.20,0.28,0.36,0.43,0.52,0.60,0.67,0.73,0.82,0.89,0.96,0.97,
  0.91,0.83,0.75,0.68,0.62,0.55,0.47,0.38,0.29,0.21,0.14,0.08,0.02,0.01,0.0};
//KL  {0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0};
  assert(sizeof(LAr_Scint_fast) == sizeof(LAr_ScintEnergy_fast));

  G4double LAr_ScintEnergy_slow[] =
  {8.97*eV,9.11*eV,9.22*eV,9.32*eV,9.37*eV,9.42*eV,9.48*eV,9.52*eV,9.56*eV,
  9.58*eV,9.62*eV,9.68*eV,9.73*eV,9.84*eV,9.90*eV,9.93*eV,9.95*eV,10.00*eV,
  10.03*eV,10.06*eV,10.10*eV,10.15*eV,10.21*eV,10.26*eV,10.33*eV,10.45*eV};

  G4double LAr_Scint_slow[] =
  {0.06,0.10,0.18,0.24,0.30,0.39,0.48,0.58,0.65,0.72,0.80,0.88,
  0.94,0.94,0.89,0.82,0.74,0.67,0.57,0.48,0.40,0.32,0.25,0.15,0.07,0.01};

  G4double LAr_RIndEnergy[] =
  {1.900*eV,  2.934*eV,  3.592*eV,  5.566*eV,  6.694*eV,  7.540*eV,  8.574*eV,  9.044*eV,
  9.232*eV,  9.420*eV,  9.514*eV,  9.608*eV,  9.702*eV,  9.796*eV,  9.890*eV,  9.984*eV,
  10.08*eV,  10.27*eV,  10.45*eV,  10.74*eV,  10.92*eV};
  const G4int LArRIndNum = sizeof(LAr_RIndEnergy)/sizeof(G4double);
  G4double LAr_RIND[]  =
  {1.232,  1.236,  1.240,  1.261,  1.282,  1.306,  1.353,  1.387,  1.404,  1.423,  1.434,
  1.446,  1.459,  1.473,  1.488,  1.505,  1.524,  1.569,  1.627,  1.751,  1.879};
  assert(sizeof(LAr_RIND) == sizeof(LAr_RIndEnergy));

  G4double LAr_AbsLEnergy[] = {4*eV, 5*eV, 6*eV, 7*eV, 8*eV, 9*eV, 10*eV, 11*eV};
  const G4int LArAbsLNum = sizeof(LAr_AbsLEnergy)/sizeof(G4double);
  G4double LAr_ABSL[]  = {2000.*cm, 2000.*cm, 2000.*cm, 2000.*cm, 2000.*cm, 2000.*cm, 2000.*cm, 2000.*cm};
  assert(sizeof(LAr_ABSL) == sizeof(LAr_AbsLEnergy));

  fLAr_mt = new G4MaterialPropertiesTable();
  fLAr_mt->AddProperty("FASTCOMPONENT", LAr_ScintEnergy_fast, LAr_Scint_fast, LArScintNum);
  fLAr_mt->AddProperty("SLOWCOMPONENT", LAr_ScintEnergy_fast, LAr_Scint_fast, LArScintNum);
  fLAr_mt->AddProperty("RINDEX",        LAr_RIndEnergy, LAr_RIND,  LArRIndNum);
  fLAr_mt->AddProperty("ABSLENGTH",     LAr_AbsLEnergy, LAr_ABSL,  LArAbsLNum);
  //BEA fLAr_mt->AddConstProperty("SCINTILLATIONYIELD",24000./MeV);
  fLAr_mt->AddConstProperty("SCINTILLATIONYIELD",2./MeV);
  fLAr_mt->AddConstProperty("RESOLUTIONSCALE",1.0);
  fLAr_mt->AddConstProperty("FASTTIMECONSTANT",6.*ns);
  fLAr_mt->AddConstProperty("SLOWTIMECONSTANT",1590.*ns);
  fLAr_mt->AddConstProperty("YIELDRATIO",0.3);
  fLAr->SetMaterialPropertiesTable(fLAr_mt);
  // Set the Birks Constant for the LAr scintillator
  fLAr->GetIonisation()->SetBirksConstant(0.69*mm/MeV);

  // Glass
  G4double glass_RIND[]={1.49,1.49,1.49};
  assert(sizeof(glass_RIND) == sizeof(lxe_Energy));
  G4double glass_AbsLength[]={420.*cm,420.*cm,420.*cm};
  assert(sizeof(glass_AbsLength) == sizeof(lxe_Energy));
  G4MaterialPropertiesTable *glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH",lxe_Energy,glass_AbsLength,lxenum); // Fix lxe_Energy
  glass_mt->AddProperty("RINDEX",lxe_Energy,glass_RIND,lxenum);
  fGlass->SetMaterialPropertiesTable(glass_mt);

  // Vacuum
  G4double vacuum_Energy[]={2.0*eV,7.0*eV,7.14*eV};
  const G4int vacnum = sizeof(vacuum_Energy)/sizeof(G4double);
  G4double vacuum_RIND[]={1.,1.,1.};
  assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));
  G4MaterialPropertiesTable *vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND,vacnum);
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  fAir->SetMaterialPropertiesTable(vacuum_mt);//Give air the same rindex

  G4double wls_Energy[] = {2.00*eV,2.87*eV,2.90*eV,3.47*eV};
  const G4int wlsnum = sizeof(wls_Energy)/sizeof(G4double);

  G4double rIndexPstyrene[]={ 1.5, 1.5, 1.5, 1.5};
  assert(sizeof(rIndexPstyrene) == sizeof(wls_Energy));
  G4double absorption1[]={2.*cm, 2.*cm, 2.*cm, 2.*cm};
  assert(sizeof(absorption1) == sizeof(wls_Energy));
  G4double scintilFast[]={0.00, 0.00, 1.00, 1.00};
  assert(sizeof(scintilFast) == sizeof(wls_Energy));
  fMPTPStyrene = new G4MaterialPropertiesTable();
  fMPTPStyrene->AddProperty("RINDEX",wls_Energy,rIndexPstyrene,wlsnum);
  fMPTPStyrene->AddProperty("ABSLENGTH",wls_Energy,absorption1,wlsnum);
  fMPTPStyrene->AddProperty("FASTCOMPONENT",wls_Energy, scintilFast,wlsnum);
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  fMPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  fPstyrene->SetMaterialPropertiesTable(fMPTPStyrene);
  // Set the Birks Constant for the Polystyrene scintillator
  fPstyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double RefractiveIndexFiber[]={ 1.60, 1.60, 1.60, 1.60};
  assert(sizeof(RefractiveIndexFiber) == sizeof(wls_Energy));
  G4double AbsFiber[]={9.00*m,9.00*m,0.1*mm,0.1*mm};
  assert(sizeof(AbsFiber) == sizeof(wls_Energy));
  G4double EmissionFib[]={1.0, 1.0, 0.0, 0.0};
  assert(sizeof(EmissionFib) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* fiberProperty = new G4MaterialPropertiesTable();
  fiberProperty->AddProperty("RINDEX",wls_Energy,RefractiveIndexFiber,wlsnum);
  fiberProperty->AddProperty("WLSABSLENGTH",wls_Energy,AbsFiber,wlsnum);
  fiberProperty->AddProperty("WLSCOMPONENT",wls_Energy,EmissionFib,wlsnum);
  fiberProperty->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
  fPMMA->SetMaterialPropertiesTable(fiberProperty);

  G4double RefractiveIndexClad1[]={ 1.49, 1.49, 1.49, 1.49};
  assert(sizeof(RefractiveIndexClad1) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* clad1Property = new G4MaterialPropertiesTable();
  clad1Property->AddProperty("RINDEX",wls_Energy,RefractiveIndexClad1,wlsnum);
  clad1Property->AddProperty("ABSLENGTH",wls_Energy,AbsFiber,wlsnum);
  fPethylene1->SetMaterialPropertiesTable(clad1Property);

  G4double RefractiveIndexClad2[]={ 1.42, 1.42, 1.42, 1.42};
  assert(sizeof(RefractiveIndexClad2) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* clad2Property = new G4MaterialPropertiesTable();
  clad2Property->AddProperty("RINDEX",wls_Energy,RefractiveIndexClad2,wlsnum);
  clad2Property->AddProperty("ABSLENGTH",wls_Energy,AbsFiber,wlsnum);
  fPethylene2->SetMaterialPropertiesTable(clad2Property);
}

G4VPhysicalVolume* LXeDetectorConstruction::Construct(){

  if (fExperimentalHall_phys) {
     G4GeometryManager::GetInstance()->OpenGeometry();
     G4PhysicalVolumeStore::GetInstance()->Clean();
     G4LogicalVolumeStore::GetInstance()->Clean();
     G4SolidStore::GetInstance()->Clean();
     G4LogicalSkinSurface::CleanSurfaceTable();
     G4LogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  return ConstructDetector();
}

G4VPhysicalVolume* LXeDetectorConstruction::ConstructDetector()
{
  //The experimental hall walls are all 1m away from housing walls
  //Bea changed to 10m
  G4double expHall_x = fScint_x+fD_mtl+1.*m;
  G4double expHall_y = fScint_y+fD_mtl+1.*m;
  G4double expHall_z = fScint_z+fD_mtl+1.*m;

  //Create experimental hall

  fExperimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box,
                                             fVacuum,"expHall_log",0,0,0);
  fExperimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                              fExperimentalHall_log,"expHall",0,false,0);

  fExperimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  fUserLimits->SetMaxAllowedStep(1.0*cm);
  fExperimentalHall_log->SetUserLimits(fUserLimits);
  //Place the main volume
  if(fMainVolumeOn){
    fMainVolume
      = new LXeMainVolume(0,G4ThreeVector(),fExperimentalHall_log,false,0,this, fUserLimits);
  }

/*
  //Place the WLS slab
  if(fWLSslab){
    G4VPhysicalVolume* slab = new LXeWLSSlab(0,G4ThreeVector(0.,0.,
                                             -fScint_z/2.-fSlab_z-1.*cm),
                                             fExperimentalHall_log,false,0,
                                             this);

    //Surface properties for the WLS slab
    G4OpticalSurface* scintWrap = new G4OpticalSurface("ScintWrap");

    new G4LogicalBorderSurface("ScintWrap", slab,
                               fExperimentalHall_phys,
                               scintWrap);

    scintWrap->SetType(dielectric_metal);
    scintWrap->SetFinish(polished);
    scintWrap->SetModel(glisur);

    G4double pp[] = {2.0*eV, 3.5*eV};
    const G4int num = sizeof(pp)/sizeof(G4double);
    G4double reflectivity[] = {1., 1.};
    assert(sizeof(reflectivity) == sizeof(pp));
    G4double efficiency[] = {0.0, 0.0};
    assert(sizeof(efficiency) == sizeof(pp));

    G4MaterialPropertiesTable* scintWrapProperty
      = new G4MaterialPropertiesTable();

    scintWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
    scintWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,num);
    scintWrap->SetMaterialPropertiesTable(scintWrapProperty);
  }
*/
  return fExperimentalHall_phys;
}

void LXeDetectorConstruction::ConstructSDandField() {

  if (!fMainVolume) return;

  // PMT SD

  if (!fPmt_SD.Get()) {
    //Created here so it exists as pmts are being placed
    G4cout << "Construction /LXeDet/pmtSD" << G4endl;
    LXePMTSD* pmt_SD = new LXePMTSD("/LXeDet/pmtSD");
    fPmt_SD.Put(pmt_SD);

    pmt_SD->InitPMTs((fNx*fNy+fNx*fNz+fNy*fNz)*2); //let pmtSD know # of pmts
    pmt_SD->SetPmtPositions(fMainVolume->GetPmtPositions());
  }

  //sensitive detector is not actually on the photocathode.
  //processHits gets done manually by the stepping action.
  //It is used to detect when photons hit and get absorbed&detected at the
  //boundary to the photocathode (which doesnt get done by attaching it to a
  //logical volume.
  //It does however need to be attached to something or else it doesnt get
  //reset at the begining of events

  SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fPmt_SD.Get());

  // Scint SD

  if (!fScint_SD.Get()) {
    G4cout << "Construction /LXeDet/scintSD" << G4endl;
    LXeScintSD* scint_SD = new LXeScintSD("/LXeDet/scintSD");
    fScint_SD.Put(scint_SD);
  }
  SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
}

void LXeDetectorConstruction::SetDimensions(G4ThreeVector dims) {
  this->fScint_x=dims[0];
  this->fScint_y=dims[1];
  this->fScint_z=dims[2];
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void LXeDetectorConstruction::SetHousingThickness(G4double d_mtl) {
  this->fD_mtl=d_mtl;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void LXeDetectorConstruction::SetNX(G4int nx) {
  this->fNx=nx;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void LXeDetectorConstruction::SetNY(G4int ny) {
  this->fNy=ny;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void LXeDetectorConstruction::SetNZ(G4int nz) {
  this->fNz=nz;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void LXeDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt) {
  this->fOuterRadius_pmt=outerRadius_pmt;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


void LXeDetectorConstruction::SetDefaults() {

  //Resets to default values
  fD_mtl=0.0635*cm;

  fScint_x = 10.0*cm;
  fScint_y = 10.0*cm;
//BEA fScint_y = 17.8*cm;
//KL  fScint_z = 22.6*cm;
  fScint_z = 1.0*m;

  fNx = 2;
  fNy = 2;
  fNz = 3;

//KL  fOuterRadius_pmt = 2.3*cm; // comment out to remove PMT

  //fSphereOn = true;
  fRefl=0.0;

  //fNfibers=15;
  //fWLSslab=false;
  fMainVolumeOn=true;
  fMainVolume=NULL;
  fSlab_z=2.5*mm;

  G4UImanager::GetUIpointer()
    ->ApplyCommand("/LXe/detector/scintYieldFactor 1.");

  //KL if(fLAr_mt)fLAr_mt->AddConstProperty("SCINTILLATIONYIELD",12000./MeV);
  if(fMPTPStyrene)fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);

  G4RunManager::GetRunManager()->ReinitializeGeometry();
}


/*
void LXeDetectorConstruction::SetSphereOn(G4bool b) {
  fSphereOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
*/


void LXeDetectorConstruction::SetHousingReflectivity(G4double r) {
  fRefl=r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void LXeDetectorConstruction::SetWLSSlabOn(G4bool b) {
  fWLSslab=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainVolumeOn(G4bool b) {
  fMainVolumeOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void LXeDetectorConstruction::SetNFibers(G4int n) {
  fNfibers=n;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeDetectorConstruction::SetMainScintYield(G4double y) {
  fLAr_mt->AddConstProperty("SCINTILLATIONYIELD",y/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*
void LXeDetectorConstruction::SetWLSScintYield(G4double y) {
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",y/MeV);
}
*/
