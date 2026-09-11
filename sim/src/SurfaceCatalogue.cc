#include "SurfaceCatalogue.hh"

#include <algorithm>
#include <array>
#include <map>

namespace scint {
namespace {

struct Key {
  Treatment treatment;
  Wrapping wrapping;
  Coupling coupling;

  bool operator<(const Key& other) const {
    return std::tie(treatment, wrapping, coupling) <
           std::tie(other.treatment, other.wrapping, other.coupling);
  }
};

// The full set Geant4 ships look-up tables for. Teflon, TiO and Tyvek exist only
// in air-coupled form; Geant4 provides no glue-coupled tables for them, and this
// table therefore omits them rather than inventing a substitute.
const std::map<Key, std::pair<G4OpticalSurfaceFinish, const char*>>& Table() {
  static const std::map<Key, std::pair<G4OpticalSurfaceFinish, const char*>> table = {
      // --- polished -------------------------------------------------------
      {{Treatment::Polished, Wrapping::None, Coupling::Air}, {polishedair, "polishedair"}},
      {{Treatment::Polished, Wrapping::Lumirror, Coupling::Air}, {polishedlumirrorair, "polishedlumirrorair"}},
      {{Treatment::Polished, Wrapping::Teflon, Coupling::Air}, {polishedteflonair, "polishedteflonair"}},
      {{Treatment::Polished, Wrapping::TiO, Coupling::Air}, {polishedtioair, "polishedtioair"}},
      {{Treatment::Polished, Wrapping::Tyvek, Coupling::Air}, {polishedtyvekair, "polishedtyvekair"}},
      {{Treatment::Polished, Wrapping::ESR, Coupling::Air}, {polishedvm2000air, "polishedvm2000air"}},
      {{Treatment::Polished, Wrapping::Lumirror, Coupling::Glue}, {polishedlumirrorglue, "polishedlumirrorglue"}},
      {{Treatment::Polished, Wrapping::ESR, Coupling::Glue}, {polishedvm2000glue, "polishedvm2000glue"}},
      // --- etched ---------------------------------------------------------
      {{Treatment::Etched, Wrapping::None, Coupling::Air}, {etchedair, "etchedair"}},
      {{Treatment::Etched, Wrapping::Lumirror, Coupling::Air}, {etchedlumirrorair, "etchedlumirrorair"}},
      {{Treatment::Etched, Wrapping::Teflon, Coupling::Air}, {etchedteflonair, "etchedteflonair"}},
      {{Treatment::Etched, Wrapping::TiO, Coupling::Air}, {etchedtioair, "etchedtioair"}},
      {{Treatment::Etched, Wrapping::Tyvek, Coupling::Air}, {etchedtyvekair, "etchedtyvekair"}},
      {{Treatment::Etched, Wrapping::ESR, Coupling::Air}, {etchedvm2000air, "etchedvm2000air"}},
      {{Treatment::Etched, Wrapping::Lumirror, Coupling::Glue}, {etchedlumirrorglue, "etchedlumirrorglue"}},
      {{Treatment::Etched, Wrapping::ESR, Coupling::Glue}, {etchedvm2000glue, "etchedvm2000glue"}},
      // --- ground ---------------------------------------------------------
      {{Treatment::Ground, Wrapping::None, Coupling::Air}, {groundair, "groundair"}},
      {{Treatment::Ground, Wrapping::Lumirror, Coupling::Air}, {groundlumirrorair, "groundlumirrorair"}},
      {{Treatment::Ground, Wrapping::Teflon, Coupling::Air}, {groundteflonair, "groundteflonair"}},
      {{Treatment::Ground, Wrapping::TiO, Coupling::Air}, {groundtioair, "groundtioair"}},
      {{Treatment::Ground, Wrapping::Tyvek, Coupling::Air}, {groundtyvekair, "groundtyvekair"}},
      {{Treatment::Ground, Wrapping::ESR, Coupling::Air}, {groundvm2000air, "groundvm2000air"}},
      {{Treatment::Ground, Wrapping::Lumirror, Coupling::Glue}, {groundlumirrorglue, "groundlumirrorglue"}},
      {{Treatment::Ground, Wrapping::ESR, Coupling::Glue}, {groundvm2000glue, "groundvm2000glue"}},
  };
  return table;
}

}  // namespace

FinishResolution ResolveFinish(const SurfaceSpec& spec) {
  FinishResolution out;
  const auto it = Table().find({spec.treatment, spec.wrapping, spec.coupling});
  if (it == Table().end()) {
    out.supported = false;
    out.reason = "Geant4 provides no look-up table for " + ToString(spec.treatment) + " + " +
                 ToString(spec.wrapping) + " + " + ToString(spec.coupling) +
                 ". Glue coupling exists only for lumirror and ESR.";
    return out;
  }
  out.supported = true;
  out.finish = it->second.first;
  out.name = it->second.second;
  return out;
}

std::vector<SurfaceSpec> SupportedSpecs() {
  std::vector<SurfaceSpec> specs;
  specs.reserve(Table().size());
  for (const auto& [key, value] : Table()) {
    specs.push_back({key.treatment, key.wrapping, key.coupling});
  }
  return specs;
}

// ---- string conversion ------------------------------------------------------

namespace {

G4String Lower(const G4String& text) {
  G4String out = text;
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return out;
}

}  // namespace

std::optional<Treatment> ParseTreatment(const G4String& text) {
  const G4String key = Lower(text);
  if (key == "polished") return Treatment::Polished;
  if (key == "etched") return Treatment::Etched;
  if (key == "ground") return Treatment::Ground;
  return std::nullopt;
}

std::optional<Wrapping> ParseWrapping(const G4String& text) {
  const G4String key = Lower(text);
  if (key == "none" || key == "bare" || key == "air") return Wrapping::None;
  if (key == "lumirror") return Wrapping::Lumirror;
  if (key == "teflon" || key == "ptfe") return Wrapping::Teflon;
  if (key == "tio" || key == "tio2") return Wrapping::TiO;
  if (key == "tyvek") return Wrapping::Tyvek;
  // "vm2000" is Geant4's name for the 3M film everyone else calls ESR.
  if (key == "esr" || key == "vm2000") return Wrapping::ESR;
  return std::nullopt;
}

std::optional<Coupling> ParseCoupling(const G4String& text) {
  const G4String key = Lower(text);
  if (key == "air") return Coupling::Air;
  if (key == "glue" || key == "meltmount" || key == "glued") return Coupling::Glue;
  return std::nullopt;
}

std::optional<G4OpticalSurfaceModel> ParseSurfaceModel(const G4String& text) {
  const G4String key = Lower(text);
  if (key == "unified") return unified;
  if (key == "lut" || key == "lbnl" || key == "lbnl_lut") return LUT;
  if (key == "davis" || key == "davis_lut") return DAVIS;
  if (key == "glisur") return glisur;
  return std::nullopt;
}

G4String ToString(Treatment value) {
  switch (value) {
    case Treatment::Polished: return "polished";
    case Treatment::Etched: return "etched";
    case Treatment::Ground: return "ground";
  }
  return "unknown";
}

G4String ToString(Wrapping value) {
  switch (value) {
    case Wrapping::None: return "none";
    case Wrapping::Lumirror: return "lumirror";
    case Wrapping::Teflon: return "teflon";
    case Wrapping::TiO: return "tio";
    case Wrapping::Tyvek: return "tyvek";
    case Wrapping::ESR: return "esr";
  }
  return "unknown";
}

G4String ToString(Coupling value) {
  switch (value) {
    case Coupling::Air: return "air";
    case Coupling::Glue: return "glue";
  }
  return "unknown";
}

G4String ToString(G4OpticalSurfaceModel value) {
  switch (value) {
    case glisur: return "glisur";
    case unified: return "unified";
    case LUT: return "lut";
    case DAVIS: return "davis";
    case dichroic: return "dichroic";
  }
  return "unknown";
}

}  // namespace scint
