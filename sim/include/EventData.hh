// Per-event tallies shared between the stacking, stepping and event actions.
#ifndef SCINT_EVENT_DATA_HH
#define SCINT_EVENT_DATA_HH 1

#include "G4Types.hh"

#include <vector>

namespace scint {

struct EventData {
  G4double edep = 0.0;                    ///< energy deposited in the crystal
  G4int generated = 0;                    ///< optical photons created
  G4int detected = 0;                     ///< photons registering at the readout
  std::vector<G4double> detectionTimes;   ///< arrival times of detected photons
  std::vector<G4double> detectionEnergies;///< photon energies at detection

  void Reset() {
    edep = 0.0;
    generated = 0;
    detected = 0;
    detectionTimes.clear();
    detectionEnergies.clear();
  }
};

}  // namespace scint

#endif
