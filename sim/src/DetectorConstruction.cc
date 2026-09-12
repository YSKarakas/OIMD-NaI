#include "DetectorConstruction.hh"

#include "MaterialLibrary.hh"

#include "G4Box.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace scint {
namespace {

/// Photon energy range over which flat optical properties are defined.
/// 1.2 eV (about 1030 nm) to 6.2 eV (about 200 nm) brackets every scintillator
/// emission band of interest here.
constexpr G4double kEnergyMin = 1.2 * eV;
constexpr G4double kEnergyMax = 6.2 * eV;

void AddFlatProperty(G4MaterialPropertiesTable* table, const G4String& key, G4double value) {
  const std::vector<G4double> energies = {kEnergyMin, kEnergyMax};
  const std::vector<G4double> values = {value, value};
  table->AddProperty(key, energies, values);
}

/// Read a "wavelength_nm value" table, e.g. a photodetector efficiency curve.
Spectrum ReadSpectrumFile(const G4String& path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error("cannot open spectrum file: " + path);
  Spectrum spectrum;
  G4String line;
  while (std::getline(file, line)) {
    const std::size_t hash = line.find('#');
    if (hash != G4String::npos) line = line.substr(0, hash);
    std::istringstream row(line);
    G4double wavelength = 0.0, value = 0.0;
    if (row >> wavelength >> value) {
      spectrum.energies.push_back(WavelengthToEnergy(wavelength));
      spectrum.values.push_back(value);
    }
  }
  if (spectrum.empty()) throw std::runtime_error("no usable rows in spectrum file: " + path);
  // Wavelength-ordered input becomes energy-reversed; Geant4 needs increasing energy.
  if (spectrum.energies.front() > spectrum.energies.back()) {
    std::reverse(spectrum.energies.begin(), spectrum.energies.end());
    std::reverse(spectrum.values.begin(), spectrum.values.end());
  }
  return spectrum;
}

}  // namespace

DetectorConstruction::DetectorConstruction() { DefineCommands(); }

DetectorConstruction::~DetectorConstruction() {
  delete fCrystalMessenger;
  delete fSurfaceMessenger;
  delete fCouplingMessenger;
  delete fReadoutMessenger;
}

void DetectorConstruction::DefineCommands() {
  fCrystalMessenger = new G4GenericMessenger(this, "/scint/crystal/", "Crystal geometry");
  fCrystalMessenger->DeclareProperty("shape", fShape, "box or cylinder");
  fCrystalMessenger->DeclarePropertyWithUnit("diameter", "mm", fDiameter, "cylinder diameter");
  fCrystalMessenger->DeclarePropertyWithUnit("width", "mm", fWidth, "box width (x)");
  fCrystalMessenger->DeclarePropertyWithUnit("height", "mm", fHeight, "box height (y)");
  fCrystalMessenger->DeclarePropertyWithUnit("length", "mm", fLength, "length along z");
  fCrystalMessenger->DeclareProperty("materialSpec", fMaterialSpec,
                                     "path to a material spec file");

  fSurfaceMessenger = new G4GenericMessenger(this, "/scint/surface/", "Wrapped surface");
  fSurfaceMessenger->DeclareProperty("treatment", fTreatment, "polished, etched or ground");
  fSurfaceMessenger->DeclareProperty("wrapping", fWrapping,
                                     "none, lumirror, teflon, tio, tyvek or esr");
  fSurfaceMessenger->DeclareProperty("coupling", fWrapCoupling, "air or glue");
  fSurfaceMessenger->DeclareProperty("model", fSurfaceModel, "unified, lut or davis");
  fSurfaceMessenger->DeclareProperty("sigmaAlpha", fSigmaAlpha,
                                     "micro-facet spread in degrees (UNIFIED model only)");

  fCouplingMessenger = new G4GenericMessenger(this, "/scint/coupling/", "Readout coupling");
  fCouplingMessenger->DeclareProperty("type", fCouplingType, "air or grease");
  fCouplingMessenger->DeclareProperty("rindex", fCouplingIndex, "coupling refractive index");
  fCouplingMessenger->DeclarePropertyWithUnit("thickness", "mm", fCouplingThickness,
                                              "coupling layer thickness");

  fReadoutMessenger = new G4GenericMessenger(this, "/scint/readout/", "Photodetector");
  fReadoutMessenger->DeclarePropertyWithUnit("diameter", "mm", fReadoutDiameter,
                                             "active diameter; 0 matches the crystal face");
  fReadoutMessenger->DeclarePropertyWithUnit("thickness", "mm", fReadoutThickness,
                                             "photodetector thickness");
  fReadoutMessenger->DeclareProperty("efficiency", fReadoutEfficiency,
                                     "flat detection efficiency used when no PDE file is given");
  fReadoutMessenger->DeclareProperty("pdeFile", fPdeFile, "wavelength_nm / efficiency table");
}

G4Material* DetectorConstruction::BuildCouplingMaterial() {
  auto* nist = G4NistManager::Instance();
  if (fCouplingType == "air") {
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    if (air->GetMaterialPropertiesTable() == nullptr) {
      auto* table = new G4MaterialPropertiesTable();
      AddFlatProperty(table, "RINDEX", 1.0);
      air->SetMaterialPropertiesTable(table);
    }
    return air;
  }
  if (fCouplingType != "grease") {
    throw std::runtime_error("unknown coupling type '" + fCouplingType + "' (expected air or grease)");
  }
  if (G4Material* existing = G4Material::GetMaterial("OpticalGrease", false)) return existing;

  // Silicone grease, approximated by its elemental composition; only the
  // refractive index matters optically and that is supplied explicitly.
  auto* grease = new G4Material("OpticalGrease", 1.06 * g / cm3, 4);
  grease->AddElement(nist->FindOrBuildElement("C"), 6);
  grease->AddElement(nist->FindOrBuildElement("H"), 18);
  grease->AddElement(nist->FindOrBuildElement("O"), 1);
  grease->AddElement(nist->FindOrBuildElement("Si"), 2);
  auto* table = new G4MaterialPropertiesTable();
  AddFlatProperty(table, "RINDEX", fCouplingIndex);
  AddFlatProperty(table, "ABSLENGTH", 10.0 * m);
  grease->SetMaterialPropertiesTable(table);
  return grease;
}

G4Material* DetectorConstruction::BuildDetectorMaterial() {
  if (G4Material* existing = G4Material::GetMaterial("DetectorWindow", false)) return existing;
  auto* nist = G4NistManager::Instance();
  G4Material* silica = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  auto* window = new G4Material("DetectorWindow", silica->GetDensity(), 1);
  window->AddMaterial(silica, 1.0);
  auto* table = new G4MaterialPropertiesTable();
  AddFlatProperty(table, "RINDEX", 1.52);  // borosilicate window
  AddFlatProperty(table, "ABSLENGTH", 1.0 * m);
  window->SetMaterialPropertiesTable(table);
  return window;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
  if (fMaterialSpec.empty()) {
    throw std::runtime_error(
        "no crystal material specified: set /scint/crystal/materialSpec before /run/initialize");
  }

  const MaterialSpec spec = ReadMaterialSpec(fMaterialSpec);
  G4Material* crystalMaterial = BuildMaterial(spec);
  G4Material* couplingMaterial = BuildCouplingMaterial();
  G4Material* detectorMaterial = BuildDetectorMaterial();

  auto* nist = G4NistManager::Instance();
  G4Material* worldMaterial = nist->FindOrBuildMaterial("G4_AIR");
  if (worldMaterial->GetMaterialPropertiesTable() == nullptr) {
    auto* table = new G4MaterialPropertiesTable();
    AddFlatProperty(table, "RINDEX", 1.0);
    worldMaterial->SetMaterialPropertiesTable(table);
  }

  const G4bool isCylinder = (fShape == "cylinder");
  if (!isCylinder && fShape != "box") {
    throw std::runtime_error("unknown crystal shape '" + fShape + "' (expected box or cylinder)");
  }

  const G4double halfLength = 0.5 * fLength;
  const G4double transverseHalf =
      isCylinder ? 0.5 * fDiameter : 0.5 * std::max(fWidth, fHeight);
  const G4double stackLength = fLength + fCouplingThickness + fReadoutThickness;
  const G4double worldHalf = std::max(transverseHalf, 0.5 * stackLength) * 3.0 + 10.0 * cm;

  auto* worldSolid = new G4Box("World", worldHalf, worldHalf, worldHalf);
  auto* worldLogical = new G4LogicalVolume(worldSolid, worldMaterial, "World");
  worldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto* worldPhysical =
      new G4PVPlacement(nullptr, {}, worldLogical, "World", nullptr, false, 0, true);

  // --- crystal -------------------------------------------------------------
  G4VSolid* crystalSolid = nullptr;
  if (isCylinder) {
    crystalSolid = new G4Tubs("Crystal", 0.0, 0.5 * fDiameter, halfLength, 0.0, CLHEP::twopi);
  } else {
    crystalSolid = new G4Box("Crystal", 0.5 * fWidth, 0.5 * fHeight, halfLength);
  }
  auto* crystalLogical = new G4LogicalVolume(crystalSolid, crystalMaterial, "Crystal");
  auto* crystalPhysical =
      new G4PVPlacement(nullptr, {}, crystalLogical, "Crystal", worldLogical, false, 0, true);

  // --- coupling layer ------------------------------------------------------
  const G4double couplingZ = halfLength + 0.5 * fCouplingThickness;
  G4VSolid* couplingSolid = nullptr;
  if (isCylinder) {
    couplingSolid = new G4Tubs("Coupling", 0.0, 0.5 * fDiameter, 0.5 * fCouplingThickness, 0.0,
                               CLHEP::twopi);
  } else {
    couplingSolid = new G4Box("Coupling", 0.5 * fWidth, 0.5 * fHeight, 0.5 * fCouplingThickness);
  }
  auto* couplingLogical = new G4LogicalVolume(couplingSolid, couplingMaterial, "Coupling");
  auto* couplingPhysical = new G4PVPlacement(nullptr, {0, 0, couplingZ}, couplingLogical,
                                             "Coupling", worldLogical, false, 0, true);

  // --- photodetector -------------------------------------------------------
  const G4double detectorZ = halfLength + fCouplingThickness + 0.5 * fReadoutThickness;
  const G4double activeRadius =
      (fReadoutDiameter > 0.0) ? 0.5 * fReadoutDiameter : 0.5 * fDiameter;
  G4VSolid* detectorSolid = nullptr;
  if (isCylinder) {
    detectorSolid = new G4Tubs("Photodetector", 0.0, activeRadius, 0.5 * fReadoutThickness, 0.0,
                               CLHEP::twopi);
  } else {
    const G4double halfX = (fReadoutDiameter > 0.0) ? activeRadius : 0.5 * fWidth;
    const G4double halfY = (fReadoutDiameter > 0.0) ? activeRadius : 0.5 * fHeight;
    detectorSolid = new G4Box("Photodetector", halfX, halfY, 0.5 * fReadoutThickness);
  }
  auto* detectorLogical =
      new G4LogicalVolume(detectorSolid, detectorMaterial, DetectorVolumeName());
  auto* detectorPhysical = new G4PVPlacement(nullptr, {0, 0, detectorZ}, detectorLogical,
                                             DetectorVolumeName(), worldLogical, false, 0, true);

  AttachSurfaces(crystalPhysical, couplingPhysical, detectorPhysical, worldPhysical);

  // Echo the configuration Geant4 actually built. This is not decoration: a
  // macro command that fails to bind executes silently and changes nothing, so
  // the only safe assumption is that the geometry must state what it became.
  G4cout << "[geometry] " << Describe() << G4endl;
  G4cout << "[material] " << fMaterialSpec << G4endl;
  return worldPhysical;
}

void DetectorConstruction::AttachSurfaces(G4VPhysicalVolume* crystal, G4VPhysicalVolume* coupling,
                                          G4VPhysicalVolume* detector, G4VPhysicalVolume* world) {
  // --- wrapped faces: crystal against the world ----------------------------
  const auto treatment = ParseTreatment(fTreatment);
  const auto wrapping = ParseWrapping(fWrapping);
  const auto wrapCoupling = ParseCoupling(fWrapCoupling);
  const auto model = ParseSurfaceModel(fSurfaceModel);
  if (!treatment) throw std::runtime_error("unknown surface treatment '" + fTreatment + "'");
  if (!wrapping) throw std::runtime_error("unknown wrapping '" + fWrapping + "'");
  if (!wrapCoupling) throw std::runtime_error("unknown surface coupling '" + fWrapCoupling + "'");
  if (!model) throw std::runtime_error("unknown surface model '" + fSurfaceModel + "'");

  // An unwrapped surface is handled analytically whatever model was asked for.
  // Geant4 ships no look-up table for a bare surface and does not say so: it
  // hangs inside G4OpBoundaryProcess::DielectricLUT() instead, at full CPU and
  // without an error. Fresnel and Snell describe a bare dielectric interface
  // exactly, so the UNIFIED model is not an approximation here -- it is the
  // right calculation. The substitution is echoed in the geometry line so that
  // a reader of the output can see that this configuration used a different
  // surface model from the wrapped ones.
  const bool bare = (*wrapping == Wrapping::None);
  G4OpticalSurfaceModel effectiveModel = *model;
  FinishResolution resolution;
  if (bare) {
    resolution = ResolveBareFinish({*treatment, *wrapping, *wrapCoupling});
    effectiveModel = unified;
  } else {
    resolution = ResolveFinish({*treatment, *wrapping, *wrapCoupling});
  }
  if (!resolution.supported) throw std::runtime_error(resolution.reason);
  fResolvedFinish = resolution.name;
  fEffectiveModel = (effectiveModel == unified)  ? "unified"
                    : (effectiveModel == DAVIS)  ? "davis"
                    : (effectiveModel == LUT)    ? "lut"
                                                 : "glisur";
  if (effectiveModel != *model) {
    fEffectiveModel += " (requested '" + fSurfaceModel + "', substituted)";
  }

  auto* wrapped = new G4OpticalSurface("WrappedSurface");
  wrapped->SetType(dielectric_LUT);
  wrapped->SetModel(effectiveModel);
  wrapped->SetFinish(resolution.finish);
  if (effectiveModel == unified) {
    // The look-up-table finishes carry their own roughness; only UNIFIED takes
    // an explicit micro-facet spread, and it must be given in radians.
    wrapped->SetType(dielectric_dielectric);
    wrapped->SetSigmaAlpha(fSigmaAlpha * deg);
  } else if (effectiveModel == DAVIS) {
    wrapped->SetType(dielectric_LUTDAVIS);
  }
  new G4LogicalBorderSurface("CrystalToWorld", crystal, world, wrapped);

  // --- readout face: crystal to coupling -----------------------------------
  auto* readoutFace = new G4OpticalSurface("CrystalToCoupling");
  readoutFace->SetType(dielectric_dielectric);
  readoutFace->SetModel(unified);
  readoutFace->SetFinish(polished);
  new G4LogicalBorderSurface("CrystalToCoupling", crystal, coupling, readoutFace);

  // --- detection surface ---------------------------------------------------
  // A photon absorbed here is recorded as detected, so the efficiency curve is
  // the photodetector's PDE. Reflectivity is zero: anything not detected is lost.
  auto* detection = new G4OpticalSurface("Detection");
  detection->SetType(dielectric_metal);
  detection->SetModel(unified);
  detection->SetFinish(polished);
  auto* table = new G4MaterialPropertiesTable();
  AddFlatProperty(table, "REFLECTIVITY", 0.0);
  if (!fPdeFile.empty()) {
    const Spectrum pde = ReadSpectrumFile(fPdeFile);
    table->AddProperty("EFFICIENCY", pde.energies, pde.values);
  } else {
    AddFlatProperty(table, "EFFICIENCY", fReadoutEfficiency);
  }
  detection->SetMaterialPropertiesTable(table);
  new G4LogicalBorderSurface("CouplingToDetector", coupling, detector, detection);
}

G4String DetectorConstruction::Describe() const {
  std::ostringstream out;
  out << "crystal=" << fShape;
  if (fShape == "cylinder") {
    out << " d=" << fDiameter / mm << "mm";
  } else {
    out << " " << fWidth / mm << "x" << fHeight / mm << "mm";
  }
  out << " L=" << fLength / mm << "mm"
      << " surface=" << fTreatment << "/" << fWrapping << "/" << fWrapCoupling
      << " (" << fResolvedFinish << ", model=" << fEffectiveModel << ")"
      << " coupling=" << fCouplingType << " n=" << fCouplingIndex
      << " readout=" << (fReadoutDiameter > 0 ? fReadoutDiameter / mm : fDiameter / mm) << "mm";
  return out.str();
}

}  // namespace scint
