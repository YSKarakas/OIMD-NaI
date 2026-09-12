// Per-event tallies shared between the stacking, stepping and event actions.
#ifndef SCINT_EVENT_DATA_HH
#define SCINT_EVENT_DATA_HH 1

#include "G4Types.hh"

#include <vector>

namespace scint {

struct EventData {
  G4double edep = 0.0;                    ///< energy deposited in the crystal
  G4int generated = 0;                    ///< optical photons created (all processes)
  G4int scintillation = 0;                ///< of those, from Scintillation
  G4int cherenkov = 0;                    ///< of those, from Cerenkov
  G4int detected = 0;                     ///< photons registering at the readout
  std::vector<G4double> detectionTimes;   ///< arrival times of detected photons
  std::vector<G4double> detectionEnergies;///< photon energies at detection

  /// Running sums of wavelength, in nm, so that the mean emitted and mean
  /// detected wavelength can be compared without storing one entry per photon
  /// (there are ~27,000 per event at 662 keV).  Their difference is the
  /// self-absorption of the crystal made visible: light that leaves the crystal
  /// is redder than the light it emitted.
  G4double generatedWavelengthSum = 0.0;
  G4double detectedWavelengthSum = 0.0;

  void Reset() {
    edep = 0.0;
    generated = 0;
    scintillation = 0;
    cherenkov = 0;
    detected = 0;
    detectionTimes.clear();
    detectionEnergies.clear();
    generatedWavelengthSum = 0.0;
    detectedWavelengthSum = 0.0;
  }
};

}  // namespace scint

#endif
