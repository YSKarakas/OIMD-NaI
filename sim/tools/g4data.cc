// g4data — dump authoritative atomic and photon-attenuation data from Geant4 itself.
//
// Rationale: the EPICS2017 cross-section files shipped with Geant4 are licensed
// "NOT FOR COMMERCIAL USE AND MUST BE USED WITHIN GEANT4", so they must not be
// vendored into a redistributable Python package. Querying them through Geant4's
// own G4EmCalculator keeps us inside that licence, guarantees the screening layer
// and the transport simulation share one source of truth, and removes any need to
// download external cross-section tables.
//
// Usage:
//   g4data version
//   g4data elements --out <file.csv>
//   g4data attenuation --name <label> --density <g/cm3> \
//                      --massfrac "Sym:frac,Sym:frac,..." \
//                      --energies "e1,e2,..."   (MeV) \
//                      --out <file.csv>
//
// Mass fractions are supplied by the caller (scint.materials computes them from
// stoichiometry in exact arithmetic); this tool only builds the material and
// evaluates cross sections.

#include "G4EmCalculator.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4Version.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Gamma.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

// ---------- tiny argument parser ----------
std::map<G4String, G4String> ParseArgs(int argc, char** argv, int from) {
  std::map<G4String, G4String> opts;
  for (int i = from; i + 1 < argc; i += 2) {
    G4String key = argv[i];
    if (key.rfind("--", 0) != 0) {
      std::cerr << "g4data: expected an option starting with '--', got '" << key << "'\n";
      std::exit(2);
    }
    opts[key.substr(2)] = argv[i + 1];
  }
  return opts;
}

G4String Require(const std::map<G4String, G4String>& o, const G4String& k) {
  auto it = o.find(k);
  if (it == o.end()) {
    std::cerr << "g4data: missing required option --" << k << "\n";
    std::exit(2);
  }
  return it->second;
}

std::vector<G4String> Split(const G4String& s, char sep) {
  std::vector<G4String> out;
  std::stringstream ss(s);
  G4String item;
  while (std::getline(ss, item, sep)) {
    if (!item.empty()) out.push_back(item);
  }
  return out;
}

// ---------- minimal Geant4 scaffolding ----------
class World : public G4VUserDetectorConstruction {
 public:
  G4VPhysicalVolume* Construct() override {
    auto* air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
    auto* box = new G4Box("World", 1 * m, 1 * m, 1 * m);
    auto* log = new G4LogicalVolume(box, air, "World");
    return new G4PVPlacement(nullptr, {}, log, "World", nullptr, false, 0);
  }
};

class Gun : public G4VUserPrimaryGeneratorAction {
 public:
  Gun() : fGun(new G4ParticleGun(1)) { fGun->SetParticleDefinition(G4Gamma::Definition()); }
  ~Gun() override { delete fGun; }
  void GeneratePrimaries(G4Event* e) override { fGun->GeneratePrimaryVertex(e); }

 private:
  G4ParticleGun* fGun;
};

class Actions : public G4VUserActionInitialization {
 public:
  void Build() const override { SetUserAction(new Gun()); }
};

class Physics : public G4VModularPhysicsList {
 public:
  Physics() { RegisterPhysics(new G4EmStandardPhysics_option4()); }
};

// ---------- commands ----------

// Strip the CVS-style "$Name: ... $" wrapper Geant4 keeps in its version tag.
G4String CleanTag(const G4String& raw) {
  const std::size_t start = raw.find(':');
  const std::size_t end = raw.rfind('$');
  if (start == G4String::npos || end == G4String::npos || end <= start + 1) return raw;
  G4String inner = raw.substr(start + 1, end - start - 1);
  const std::size_t first = inner.find_first_not_of(" \t");
  const std::size_t last = inner.find_last_not_of(" \t");
  if (first == G4String::npos) return raw;
  return inner.substr(first, last - first + 1);
}

// geant4-config --version cannot distinguish a beta from a release: the beta of
// a series already carries the target G4VERSION_NUMBER. The tag string is the
// only reliable discriminator, so it is exported for the provenance record.
int PrintVersion() {
  std::cout << "{\n"
            << "  \"version_number\": " << G4VERSION_NUMBER << ",\n"
            << "  \"tag\": \"" << CleanTag(G4VERSION_TAG) << "\",\n"
            << "  \"build_tag\": \"" << CleanTag(G4Version) << "\",\n"
            << "  \"date\": \"" << G4Date << "\",\n"
            << "  \"reference_tag\": " << G4VERSION_REFERENCE_TAG << "\n"
            << "}\n";
  return 0;
}

int DumpElements(const std::map<G4String, G4String>& opts) {
  const G4String out = Require(opts, "out");
  auto* nist = G4NistManager::Instance();
  const std::vector<G4String>& names = nist->GetNistElementNames();

  std::ofstream f(out);
  if (!f) {
    std::cerr << "g4data: cannot write " << out << "\n";
    return 1;
  }
  f << "# Atomic data exported from Geant4 G4NistManager\n";
  f << "# Source: Geant4 " << G4Version << " NIST element database\n";
  f << "z,symbol,atomic_mass_amu\n";
  f << std::setprecision(10);
  for (std::size_t z = 1; z < names.size(); ++z) {
    const G4String& sym = names[z];
    if (sym.empty()) continue;
    f << z << "," << sym << "," << nist->GetAtomicMassAmu(static_cast<G4int>(z)) << "\n";
  }
  std::cout << "g4data: wrote " << out << " (" << names.size() - 1 << " elements)\n";
  return 0;
}

int DumpAttenuation(const std::map<G4String, G4String>& opts) {
  const G4String label = Require(opts, "name");
  const G4double density = std::stod(Require(opts, "density"));
  const G4String massfrac = Require(opts, "massfrac");
  const G4String energies = Require(opts, "energies");
  const G4String out = Require(opts, "out");

  auto* nist = G4NistManager::Instance();

  std::vector<std::pair<G4String, G4double>> comp;
  G4double fracSum = 0.0;
  for (const auto& tok : Split(massfrac, ',')) {
    auto kv = Split(tok, ':');
    if (kv.size() != 2) {
      std::cerr << "g4data: malformed mass-fraction entry '" << tok << "'\n";
      return 2;
    }
    const G4double v = std::stod(kv[1]);
    comp.emplace_back(kv[0], v);
    fracSum += v;
  }
  if (std::abs(fracSum - 1.0) > 1e-6) {
    std::cerr << "g4data: mass fractions sum to " << fracSum << ", expected 1 (tolerance 1e-6)\n";
    return 2;
  }

  auto* mat = new G4Material(label, density * g / cm3, static_cast<G4int>(comp.size()));
  for (const auto& [sym, frac] : comp) {
    G4Element* el = nist->FindOrBuildElement(sym);
    if (!el) {
      std::cerr << "g4data: unknown element symbol '" << sym << "'\n";
      return 2;
    }
    mat->AddElement(el, frac);
  }

  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
  rm->SetUserInitialization(new World());
  rm->SetUserInitialization(new Physics());
  rm->SetUserInitialization(new Actions());
  rm->Initialize();
  rm->BeamOn(0);  // builds physics tables; required before G4EmCalculator queries

  G4EmCalculator calc;
  std::ofstream f(out);
  if (!f) {
    std::cerr << "g4data: cannot write " << out << "\n";
    return 1;
  }
  f << "# Photon attenuation computed by Geant4 " << G4Version << "\n";
  f << "# Physics list: G4EmStandardPhysics_option4\n";
  f << "# Material: " << label << ", density " << density << " g/cm3\n";
  f << "# attenuation_length is the total mean free path (1/mu).\n";
  f << "energy_MeV,attenuation_length_cm,mu_per_cm,mu_over_rho_cm2_per_g\n";
  f << std::setprecision(10);
  for (const auto& tok : Split(energies, ',')) {
    const G4double e = std::stod(tok) * MeV;
    const G4double lambda = calc.ComputeGammaAttenuationLength(e, mat);
    const G4double lambda_cm = lambda / cm;
    const G4double mu = (lambda_cm > 0.0) ? 1.0 / lambda_cm : 0.0;
    f << e / MeV << "," << lambda_cm << "," << mu << "," << mu / density << "\n";
  }
  std::cout << "g4data: wrote " << out << "\n";
  delete rm;
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: g4data <version|elements|attenuation> [--opt value ...]\n";
    return 2;
  }
  const G4String cmd = argv[1];
  if (cmd == "version") return PrintVersion();
  const auto opts = ParseArgs(argc, argv, 2);
  if (cmd == "elements") return DumpElements(opts);
  if (cmd == "attenuation") return DumpAttenuation(opts);
  std::cerr << "g4data: unknown command '" << cmd << "'\n";
  return 2;
}
