// Building Geant4 materials, with optical properties, from a plain-text spec.
//
// Every optical property a scintillator simulation needs must be supplied by the
// user -- Geant4 derives none of them from chemistry. The spec file is therefore
// the physical input to the study, and it is kept human-readable and diffable so
// that a reviewer can see exactly which numbers produced a result, and where
// each came from.
//
// Format (see materials/ for examples):
//
//     name          = NaI_Tl
//     density_g_cm3 = 3.67
//     composition   = Na:0.153373922,I:0.846626078   # mass fractions
//     scintillation_yield_per_MeV = 38000
//     resolution_scale            = 1.0
//     time_constant_1_ns          = 250
//
//     [RINDEX]                 # wavelength_nm  value
//     300  1.85
//     ...
//     [ABSLENGTH_mm]
//     ...
//     [SCINTILLATIONCOMPONENT1]
//     ...
//
// Wavelengths are given in nm because that is how optical data is published;
// they are converted to photon energies and re-sorted internally, since Geant4
// requires monotonically increasing energy.

#ifndef SCINT_MATERIAL_LIBRARY_HH
#define SCINT_MATERIAL_LIBRARY_HH 1

#include "G4Material.hh"
#include "G4String.hh"

#include <map>
#include <vector>

namespace scint {

/// Photon energy / value pairs, sorted by increasing energy as Geant4 requires.
struct Spectrum {
  std::vector<G4double> energies;  ///< eV
  std::vector<G4double> values;
  bool empty() const { return energies.empty(); }
};

struct MaterialSpec {
  G4String name;
  G4double density_g_cm3 = 0.0;
  std::map<G4String, G4double> mass_fractions;
  std::map<G4String, G4double> scalars;   ///< e.g. scintillation_yield_per_MeV
  std::map<G4String, Spectrum> spectra;   ///< section name -> spectrum
  G4String source_path;
};

/// Parse a material spec file. Throws std::runtime_error with the offending line.
MaterialSpec ReadMaterialSpec(const G4String& path);

/// Construct the Geant4 material and attach its optical property table.
///
/// Returns the existing material if one of that name was already built, so that
/// repeated calls within a run are harmless.
G4Material* BuildMaterial(const MaterialSpec& spec);

/// Convert a wavelength in nm to a photon energy in Geant4 units.
G4double WavelengthToEnergy(G4double wavelength_nm);

}  // namespace scint

#endif
