// Orthogonal parameterisation of scintillator surface treatment and wrapping.
//
// Geant4 encodes surface treatment, wrapping material and optical coupling in a
// single flat enum (polishedteflonair, groundvm2000glue, ...). That is compact
// but it hides the structure of the parameter space, and a light-collection
// study needs to vary the three axes independently. This header exposes them as
// three enumerations and maps them onto the Geant4 finish, refusing the
// combinations Geant4 has no look-up table for rather than silently substituting
// a neighbour.
//
// The supported grid is 3 treatments x 6 air-coupled wrappings, plus the 3 x 2
// glue-coupled cases Geant4 provides, i.e. 24 combinations in total.

#ifndef SCINT_SURFACE_CATALOGUE_HH
#define SCINT_SURFACE_CATALOGUE_HH 1

#include "G4OpticalSurface.hh"
#include "G4String.hh"

#include <optional>
#include <vector>

namespace scint {

/// Mechanical treatment of the crystal face.
enum class Treatment {
  Polished,  ///< mechanically polished
  Etched,    ///< chemically etched
  Ground     ///< rough-cut / ground
};

/// Reflector applied outside the crystal face.
enum class Wrapping {
  None,      ///< bare surface facing air
  Lumirror,  ///< Lumirror specular reflector
  Teflon,    ///< PTFE tape, diffuse
  TiO,       ///< titanium dioxide paint, diffuse
  Tyvek,     ///< Tyvek, diffuse
  ESR        ///< 3M ESR / VM2000 specular film
};

/// Optical coupling between crystal and reflector.
enum class Coupling {
  Air,  ///< air gap
  Glue  ///< index-matched bond (Geant4: "meltmount")
};

struct SurfaceSpec {
  Treatment treatment = Treatment::Polished;
  Wrapping wrapping = Wrapping::None;
  Coupling coupling = Coupling::Air;
};

/// Result of resolving a specification against Geant4's look-up tables.
struct FinishResolution {
  bool supported = false;
  G4OpticalSurfaceFinish finish = polished;
  G4String name;    ///< the Geant4 enumerator name, for the provenance record
  G4String reason;  ///< why an unsupported combination was refused
};

/// Map a (treatment, wrapping, coupling) triple onto a Geant4 surface finish.
///
/// Never falls back to an approximate neighbour: an unsupported combination
/// returns supported == false with an explanation. Silently substituting a
/// different wrapping would corrupt exactly the comparison this study makes.
FinishResolution ResolveFinish(const SurfaceSpec& spec);

/// Every combination Geant4 supports, in a stable order.
std::vector<SurfaceSpec> SupportedSpecs();

// ---- string conversion (for macros, configuration files and reports) ----

std::optional<Treatment> ParseTreatment(const G4String& text);
std::optional<Wrapping> ParseWrapping(const G4String& text);
std::optional<Coupling> ParseCoupling(const G4String& text);
std::optional<G4OpticalSurfaceModel> ParseSurfaceModel(const G4String& text);

G4String ToString(Treatment value);
G4String ToString(Wrapping value);
G4String ToString(Coupling value);
G4String ToString(G4OpticalSurfaceModel value);

}  // namespace scint

#endif
