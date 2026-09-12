#include "Actions.hh"

#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4GenericMessenger.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <numeric>

namespace scint {
namespace {

/// Photon energy (Geant4 units) to wavelength in nm.
///
/// The constant is CODATA hc in eV nm, the same value MaterialLibrary.cc uses to
/// convert the material tables the other way. They must agree, or a wavelength
/// reported here would not line up with the table that produced the photon.
constexpr G4double kHcEvNm = 1239.841984;

G4double WavelengthNm(G4double energy) { return kHcEvNm / (energy / eV); }

/// Locate the boundary process once; Geant4 offers no direct accessor.
G4OpBoundaryProcess* BoundaryProcess() {
  static G4OpBoundaryProcess* cached = nullptr;
  if (cached != nullptr) return cached;
  G4ProcessManager* manager = G4OpticalPhoton::Definition()->GetProcessManager();
  if (manager == nullptr) return nullptr;
  G4ProcessVector* processes = manager->GetProcessList();
  for (G4int i = 0; i < static_cast<G4int>(processes->size()); ++i) {
    if (auto* boundary = dynamic_cast<G4OpBoundaryProcess*>((*processes)[i])) {
      cached = boundary;
      return cached;
    }
  }
  return nullptr;
}

EventAction* CurrentEventAction() {
  return const_cast<EventAction*>(static_cast<const EventAction*>(
      G4RunManager::GetRunManager()->GetUserEventAction()));
}

}  // namespace

// --------------------------------------------------------------------------- //
// Primary generator
// --------------------------------------------------------------------------- //

PrimaryGeneratorAction::PrimaryGeneratorAction() {
  fGun = new G4ParticleGun(1);
  fGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
  fGun->SetParticleMomentumDirection({0.0, 0.0, 1.0});
  DefineCommands();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fGun;
  delete fMessenger;
}

void PrimaryGeneratorAction::DefineCommands() {
  fMessenger = new G4GenericMessenger(this, "/scint/source/", "Primary particle source");
  fMessenger->DeclarePropertyWithUnit("energy", "keV", fEnergy, "gamma energy");
  fMessenger->DeclarePropertyWithUnit("offsetZ", "mm", fOffsetZ,
                                      "start position along z, relative to the origin");
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  fGun->SetParticleEnergy(fEnergy);
  fGun->SetParticlePosition({0.0, 0.0, fOffsetZ});
  fGun->GeneratePrimaryVertex(event);
}

// --------------------------------------------------------------------------- //
// Event
// --------------------------------------------------------------------------- //

void EventAction::BeginOfEventAction(const G4Event*) { fData.Reset(); }

void EventAction::EndOfEventAction(const G4Event* event) {
  auto* analysis = G4AnalysisManager::Instance();

  G4double firstTime = -1.0;
  G4double meanTime = -1.0;
  if (!fData.detectionTimes.empty()) {
    firstTime = *std::min_element(fData.detectionTimes.begin(), fData.detectionTimes.end());
    meanTime = std::accumulate(fData.detectionTimes.begin(), fData.detectionTimes.end(), 0.0) /
               static_cast<G4double>(fData.detectionTimes.size());
  }

  G4int column = 0;
  analysis->FillNtupleIColumn(column++, event->GetEventID());
  analysis->FillNtupleDColumn(column++, fData.edep / keV);
  analysis->FillNtupleIColumn(column++, fData.generated);
  analysis->FillNtupleIColumn(column++, fData.scintillation);
  analysis->FillNtupleIColumn(column++, fData.cherenkov);
  analysis->FillNtupleIColumn(column++, fData.detected);
  analysis->FillNtupleIColumn(column++, fData.killedAtStepLimit);
  analysis->FillNtupleDColumn(column++, firstTime < 0 ? -1.0 : firstTime / ns);
  analysis->FillNtupleDColumn(column++, meanTime < 0 ? -1.0 : meanTime / ns);
  // Mean wavelength emitted and mean wavelength detected. The gap between them
  // is bulk self-absorption: the blue side of the emission band runs into the
  // crystal's own absorption edge, so the light that survives to the readout is
  // redder than the light that was emitted. -1 marks "no photons".
  analysis->FillNtupleDColumn(
      column++, fData.generated > 0
                    ? fData.generatedWavelengthSum / static_cast<G4double>(fData.generated)
                    : -1.0);
  analysis->FillNtupleDColumn(
      column++, fData.detected > 0
                    ? fData.detectedWavelengthSum / static_cast<G4double>(fData.detected)
                    : -1.0);
  analysis->AddNtupleRow();
}

// --------------------------------------------------------------------------- //
// Stacking: count optical photons at creation
// --------------------------------------------------------------------------- //

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track) {
  if (track->GetDefinition() == G4OpticalPhoton::Definition() && track->GetParentID() > 0) {
    EventData& data = fEventAction->Data();
    data.generated += 1;
    data.generatedWavelengthSum += WavelengthNm(track->GetTotalEnergy());
    // Optical photons are not all scintillation photons. G4OpticalPhysics also
    // enables Cerenkov, and in NaI (n ~ 1.85) the Cerenkov threshold for
    // electrons is about 95 keV, which Compton electrons from a 662 keV gamma
    // pass easily. Counting "optical photons" as "scintillation photons" would
    // therefore overstate the light yield -- by a little here, by much more in a
    // dense high-index crystal or at higher energy. They are separated so that
    // the photon-budget gate tests what it claims to test.
    const G4VProcess* creator = track->GetCreatorProcess();
    if (creator != nullptr) {
      const G4String& name = creator->GetProcessName();
      if (name == "Scintillation") {
        data.scintillation += 1;
      } else if (name == "Cerenkov") {
        data.cherenkov += 1;
      }
    }
  }
  return fUrgent;
}

// --------------------------------------------------------------------------- //
// Stepping: energy deposition and photon detection
// --------------------------------------------------------------------------- //

SteppingAction::SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {
  DefineCommands();
}

SteppingAction::~SteppingAction() { delete fMessenger; }

void SteppingAction::DefineCommands() {
  fMessenger = new G4GenericMessenger(this, "/scint/optical/", "Optical photon transport");
  fMessenger->DeclareProperty(
      "maxSteps", fMaxOpticalSteps,
      "kill an optical photon after this many steps (0 = no limit); the number "
      "killed is recorded per event so the bias can be measured");
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  EventData& data = fEventAction->Data();

  // The step cap, applied before anything else so that a capped photon cannot
  // also be counted as detected in the same step.
  if (fMaxOpticalSteps > 0 &&
      step->GetTrack()->GetDefinition() == G4OpticalPhoton::Definition() &&
      step->GetTrack()->GetCurrentStepNumber() >= fMaxOpticalSteps) {
    data.killedAtStepLimit += 1;
    const_cast<G4Track*>(step->GetTrack())->SetTrackStatus(fStopAndKill);
    return;
  }

  const G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if (volume != nullptr && volume->GetName() == "Crystal" &&
      step->GetTrack()->GetDefinition() != G4OpticalPhoton::Definition()) {
    data.edep += step->GetTotalEnergyDeposit();
  }

  if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::Definition()) return;
  if (step->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) return;

  G4OpBoundaryProcess* boundary = BoundaryProcess();
  if (boundary == nullptr) return;
  if (boundary->GetStatus() != Detection) return;

  data.detected += 1;
  data.detectionTimes.push_back(step->GetPostStepPoint()->GetGlobalTime());
  data.detectionEnergies.push_back(step->GetTrack()->GetTotalEnergy());
  data.detectedWavelengthSum += WavelengthNm(step->GetTrack()->GetTotalEnergy());
}

// --------------------------------------------------------------------------- //
// Run
// --------------------------------------------------------------------------- //

RunAction::RunAction() {
  DefineCommands();
  auto* analysis = G4AnalysisManager::Instance();
  analysis->SetDefaultFileType("root");
  analysis->SetVerboseLevel(0);
  analysis->CreateNtuple("events", "Per-event scintillation tallies");
  analysis->CreateNtupleIColumn("event");
  analysis->CreateNtupleDColumn("edep_keV");
  analysis->CreateNtupleIColumn("photons_generated");
  analysis->CreateNtupleIColumn("photons_scintillation");
  analysis->CreateNtupleIColumn("photons_cherenkov");
  analysis->CreateNtupleIColumn("photons_detected");
  analysis->CreateNtupleIColumn("photons_killed_steplimit");
  analysis->CreateNtupleDColumn("first_detection_ns");
  analysis->CreateNtupleDColumn("mean_detection_ns");
  analysis->CreateNtupleDColumn("mean_generated_nm");
  analysis->CreateNtupleDColumn("mean_detected_nm");
  analysis->FinishNtuple();
}

RunAction::~RunAction() { delete fMessenger; }

void RunAction::DefineCommands() {
  fMessenger = new G4GenericMessenger(this, "/scint/output/", "Run output");
  fMessenger->DeclareProperty("file", fOutputFile, "output file stem (without extension)");
  fMessenger->DeclareProperty("format", fOutputFormat, "root or csv");
}

void RunAction::BeginOfRunAction(const G4Run*) {
  auto* analysis = G4AnalysisManager::Instance();
  analysis->SetDefaultFileType(fOutputFormat);
  analysis->OpenFile(fOutputFile);
}

void RunAction::EndOfRunAction(const G4Run* run) {
  auto* analysis = G4AnalysisManager::Instance();
  analysis->Write();
  analysis->CloseFile();

  const auto* detector = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4cout << "\n=== run summary ===\n"
         << "events: " << run->GetNumberOfEvent() << "\n"
         << "geometry: " << detector->Describe() << "\n"
         << "output: " << fOutputFile << "." << fOutputFormat << "\n"
         << "===================" << G4endl;
}

// --------------------------------------------------------------------------- //

void ActionInitialization::Build() const {
  SetUserAction(new PrimaryGeneratorAction());
  auto* eventAction = new EventAction();
  SetUserAction(eventAction);
  SetUserAction(new StackingAction(eventAction));
  SetUserAction(new SteppingAction(eventAction));
  SetUserAction(new RunAction());
}

}  // namespace scint
