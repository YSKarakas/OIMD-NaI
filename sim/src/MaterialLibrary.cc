#include "MaterialLibrary.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <algorithm>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace scint {
namespace {

// CODATA: hc expressed for the nm <-> eV conversion used throughout optics.
constexpr G4double kHcEvNm = 1239.841984;

G4String Trim(const G4String& text) {
  const std::size_t first = text.find_first_not_of(" \t\r\n");
  if (first == G4String::npos) return "";
  const std::size_t last = text.find_last_not_of(" \t\r\n");
  return text.substr(first, last - first + 1);
}

G4String StripComment(const G4String& line) {
  const std::size_t hash = line.find('#');
  return hash == G4String::npos ? line : G4String(line.substr(0, hash));
}

std::vector<G4String> Split(const G4String& text, char sep) {
  std::vector<G4String> parts;
  std::stringstream ss(text);
  G4String item;
  while (std::getline(ss, item, sep)) {
    const G4String trimmed = Trim(item);
    if (!trimmed.empty()) parts.push_back(trimmed);
  }
  return parts;
}

/// Sort a spectrum by increasing photon energy, as G4MaterialPropertiesTable requires.
void SortByEnergy(Spectrum& spectrum) {
  std::vector<std::size_t> order(spectrum.energies.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return spectrum.energies[a] < spectrum.energies[b];
  });
  std::vector<G4double> energies, values;
  energies.reserve(order.size());
  values.reserve(order.size());
  for (const std::size_t index : order) {
    energies.push_back(spectrum.energies[index]);
    values.push_back(spectrum.values[index]);
  }
  spectrum.energies = std::move(energies);
  spectrum.values = std::move(values);
}

}  // namespace

G4double WavelengthToEnergy(G4double wavelength_nm) {
  if (wavelength_nm <= 0.0) {
    throw std::runtime_error("wavelength must be positive, got " + std::to_string(wavelength_nm));
  }
  return (kHcEvNm / wavelength_nm) * eV;
}

MaterialSpec ReadMaterialSpec(const G4String& path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error("cannot open material spec: " + path);

  MaterialSpec spec;
  spec.source_path = path;
  G4String section;  // empty => scalar section
  G4String raw;
  std::size_t lineNumber = 0;

  while (std::getline(file, raw)) {
    ++lineNumber;
    const G4String line = Trim(StripComment(raw));
    if (line.empty()) continue;

    if (line.front() == '[' && line.back() == ']') {
      section = line.substr(1, line.size() - 2);
      spec.spectra[section];  // ensure the section exists even if empty
      continue;
    }

    if (section.empty()) {
      const std::size_t equals = line.find('=');
      if (equals == G4String::npos) {
        throw std::runtime_error(path + ":" + std::to_string(lineNumber) +
                                 ": expected 'key = value', got: " + line);
      }
      const G4String key = Trim(line.substr(0, equals));
      const G4String value = Trim(line.substr(equals + 1));
      if (key == "name") {
        spec.name = value;
      } else if (key == "density_g_cm3") {
        spec.density_g_cm3 = std::stod(value);
      } else if (key == "composition") {
        for (const auto& token : Split(value, ',')) {
          const auto pieces = Split(token, ':');
          if (pieces.size() != 2) {
            throw std::runtime_error(path + ":" + std::to_string(lineNumber) +
                                     ": expected 'Symbol:fraction', got: " + token);
          }
          spec.mass_fractions[pieces[0]] = std::stod(pieces[1]);
        }
      } else {
        spec.scalars[key] = std::stod(value);
      }
      continue;
    }

    // Spectrum row: wavelength_nm value
    std::istringstream row(line);
    G4double wavelength = 0.0, value = 0.0;
    if (!(row >> wavelength >> value)) {
      throw std::runtime_error(path + ":" + std::to_string(lineNumber) +
                               ": expected 'wavelength_nm value', got: " + line);
    }
    spec.spectra[section].energies.push_back(WavelengthToEnergy(wavelength));
    spec.spectra[section].values.push_back(value);
  }

  if (spec.name.empty()) throw std::runtime_error(path + ": missing 'name'");
  if (spec.density_g_cm3 <= 0.0) throw std::runtime_error(path + ": missing or invalid 'density_g_cm3'");
  if (spec.mass_fractions.empty()) throw std::runtime_error(path + ": missing 'composition'");

  const G4double total = std::accumulate(
      spec.mass_fractions.begin(), spec.mass_fractions.end(), 0.0,
      [](G4double sum, const auto& entry) { return sum + entry.second; });
  if (std::abs(total - 1.0) > 1e-6) {
    throw std::runtime_error(path + ": mass fractions sum to " + std::to_string(total) +
                             ", expected 1 within 1e-6");
  }

  for (auto& [name, spectrum] : spec.spectra) SortByEnergy(spectrum);
  return spec;
}

G4Material* BuildMaterial(const MaterialSpec& spec) {
  if (G4Material* existing = G4Material::GetMaterial(spec.name, false)) return existing;

  auto* nist = G4NistManager::Instance();
  auto* material = new G4Material(spec.name, spec.density_g_cm3 * g / cm3,
                                  static_cast<G4int>(spec.mass_fractions.size()));
  for (const auto& [symbol, fraction] : spec.mass_fractions) {
    G4Element* element = nist->FindOrBuildElement(symbol);
    if (element == nullptr) {
      throw std::runtime_error(spec.source_path + ": unknown element symbol '" + symbol + "'");
    }
    material->AddElement(element, fraction);
  }

  auto* table = new G4MaterialPropertiesTable();

  // Spectra. Section names carry their unit suffix where one is needed, because
  // a silently wrong length unit is the easiest way to get a plausible but
  // wrong light-collection efficiency.
  for (const auto& [section, spectrum] : spec.spectra) {
    if (spectrum.empty()) continue;
    if (section == "ABSLENGTH_mm") {
      std::vector<G4double> scaled;
      scaled.reserve(spectrum.values.size());
      for (const G4double value : spectrum.values) scaled.push_back(value * mm);
      table->AddProperty("ABSLENGTH", spectrum.energies, scaled);
    } else {
      table->AddProperty(section, spectrum.energies, spectrum.values);
    }
  }

  // Scalars.
  for (const auto& [key, value] : spec.scalars) {
    if (key == "scintillation_yield_per_MeV") {
      table->AddConstProperty("SCINTILLATIONYIELD", value / MeV);
    } else if (key == "resolution_scale") {
      table->AddConstProperty("RESOLUTIONSCALE", value);
    } else if (key == "time_constant_1_ns") {
      table->AddConstProperty("SCINTILLATIONTIMECONSTANT1", value * ns);
    } else if (key == "time_constant_2_ns") {
      table->AddConstProperty("SCINTILLATIONTIMECONSTANT2", value * ns);
    } else if (key == "rise_time_1_ns") {
      table->AddConstProperty("SCINTILLATIONRISETIME1", value * ns);
    } else if (key == "yield_ratio_1") {
      table->AddConstProperty("SCINTILLATIONYIELD1", value);
    } else if (key == "yield_ratio_2") {
      table->AddConstProperty("SCINTILLATIONYIELD2", value);
    } else if (key == "birks_mm_per_MeV") {
      material->GetIonisation()->SetBirksConstant(value * mm / MeV);
    }
    // Unrecognised scalars are left in the spec: they document provenance
    // (measurement temperature, source of a number) without affecting physics.
  }

  material->SetMaterialPropertiesTable(table);
  return material;
}

}  // namespace scint
