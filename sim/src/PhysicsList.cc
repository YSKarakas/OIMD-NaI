#include "PhysicsList.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4OpticalParameters.hh"
#include "G4OpticalPhysics.hh"
#include "G4SystemOfUnits.hh"

namespace scint {

PhysicsList::PhysicsList() {
  RegisterPhysics(new G4EmStandardPhysics_option4());
  RegisterPhysics(new G4OpticalPhysics());

  auto* optical = G4OpticalParameters::Instance();
  // Track optical photons before the rest of the event. With tens of thousands
  // of photons per gamma this keeps the stack shallow; it does not change physics.
  optical->SetScintTrackSecondariesFirst(true);

  SetDefaultCutValue(0.1 * mm);
}

}  // namespace scint
