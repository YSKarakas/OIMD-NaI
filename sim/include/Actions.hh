// User actions: primary generation, photon tallying and output.
//
// The division of labour is the conventional one. StackingAction counts optical
// photons as they are created; SteppingAction records energy deposition in the
// crystal and photons that reach detection status at the readout surface;
// EventAction owns the tallies and fills one ntuple row per event; RunAction
// opens and closes the output file.

#ifndef SCINT_ACTIONS_HH
#define SCINT_ACTIONS_HH 1

#include "G4String.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserStackingAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "EventData.hh"

class G4GenericMessenger;
class G4ParticleGun;

namespace scint {

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event* event) override;

 private:
  void DefineCommands();
  G4ParticleGun* fGun = nullptr;
  G4GenericMessenger* fMessenger = nullptr;
  G4double fEnergy = 662.0;      ///< keV
  G4double fOffsetZ = -100.0;    ///< mm, relative to the crystal entrance face
};

class EventAction : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;
  EventData& Data() { return fData; }

 private:
  EventData fData;
};

class StackingAction : public G4UserStackingAction {
 public:
  explicit StackingAction(EventAction* eventAction) : fEventAction(eventAction) {}
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;

 private:
  EventAction* fEventAction;
};

class SteppingAction : public G4UserSteppingAction {
 public:
  explicit SteppingAction(EventAction* eventAction) : fEventAction(eventAction) {}
  void UserSteppingAction(const G4Step* step) override;

 private:
  EventAction* fEventAction;
};

class RunAction : public G4UserRunAction {
 public:
  RunAction();
  ~RunAction() override;
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  void DefineCommands();
  G4GenericMessenger* fMessenger = nullptr;
  G4String fOutputFile = "output";
  G4String fOutputFormat = "csv";
};

class ActionInitialization : public G4VUserActionInitialization {
 public:
  void Build() const override;
};

}  // namespace scint

#endif
