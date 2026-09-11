// Physics for a gamma-excited scintillator with optical photon transport.
//
// G4EmStandardPhysics_option4 is chosen because it uses the most accurate
// low-energy electromagnetic models available, which matters at the tens-to-
// hundreds of keV where photoelectric absorption dominates in a high-Z crystal.
// Hadronic physics is deliberately absent: nothing here produces hadrons.

#ifndef SCINT_PHYSICS_LIST_HH
#define SCINT_PHYSICS_LIST_HH 1

#include "G4VModularPhysicsList.hh"

namespace scint {

class PhysicsList : public G4VModularPhysicsList {
 public:
  PhysicsList();
};

}  // namespace scint

#endif
