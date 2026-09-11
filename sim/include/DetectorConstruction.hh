// Parametric scintillator + coupling + photodetector geometry.
//
// Layout along +z, all volumes sharing a transverse profile (box or cylinder):
//
//     world (air)
//       |-- crystal           centred on the origin
//       |-- coupling layer    against the crystal's +z face
//       `-- photodetector     against the coupling layer
//
// Optical surfaces:
//   crystal -> world        the wrapped surface (treatment + wrapping + coupling),
//                           resolved through SurfaceCatalogue. Because the
//                           coupling layer physically occupies the +z face, this
//                           border covers every face except the readout face --
//                           which is what a wrapped, one-end-read crystal is.
//   crystal -> coupling     dielectric-dielectric, polished
//   coupling -> detector    dielectric-metal carrying the photon detection
//                           efficiency; a photon absorbed here counts as detected.
//
// Geometry is set from a macro before /run/initialize via G4GenericMessenger.

#ifndef SCINT_DETECTOR_CONSTRUCTION_HH
#define SCINT_DETECTOR_CONSTRUCTION_HH 1

#include "G4VUserDetectorConstruction.hh"
#include "G4String.hh"
#include "G4Types.hh"

#include "SurfaceCatalogue.hh"

class G4GenericMessenger;
class G4LogicalVolume;
class G4Material;
class G4VPhysicalVolume;

namespace scint {

class DetectorConstruction : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  /// Name of the logical volume that registers photon detections.
  static G4String DetectorVolumeName() { return "Photodetector"; }

  /// One-line description of the resolved configuration, for the run record.
  G4String Describe() const;

 private:
  void DefineCommands();
  G4Material* BuildCouplingMaterial();
  G4Material* BuildDetectorMaterial();
  void AttachSurfaces(G4VPhysicalVolume* crystal, G4VPhysicalVolume* coupling,
                      G4VPhysicalVolume* detector, G4VPhysicalVolume* world);

  // One messenger per command directory. G4GenericMessenger keys its property
  // map by the declared name but looks values up by the command's last path
  // token, so a nested name such as "crystal/shape" is registered under a key
  // that can never be found again -- the command then executes silently and
  // changes nothing. Flat names under a per-directory messenger avoid this.
  G4GenericMessenger* fCrystalMessenger = nullptr;
  G4GenericMessenger* fSurfaceMessenger = nullptr;
  G4GenericMessenger* fCouplingMessenger = nullptr;
  G4GenericMessenger* fReadoutMessenger = nullptr;

  // --- crystal ---
  G4String fShape = "cylinder";      ///< box | cylinder
  G4double fDiameter = 76.2;         ///< mm, cylinder
  G4double fWidth = 50.0;            ///< mm, box (x)
  G4double fHeight = 50.0;           ///< mm, box (y)
  G4double fLength = 76.2;           ///< mm, along z
  G4String fMaterialSpec;            ///< path to a material spec file

  // --- surface ---
  G4String fTreatment = "polished";
  G4String fWrapping = "teflon";
  G4String fWrapCoupling = "air";
  G4String fSurfaceModel = "lut";
  G4double fSigmaAlpha = 0.0;        ///< degrees, UNIFIED model only

  // --- optical coupling to the readout ---
  G4String fCouplingType = "grease"; ///< air | grease
  G4double fCouplingIndex = 1.465;
  G4double fCouplingThickness = 0.1; ///< mm

  // --- readout ---
  G4double fReadoutDiameter = 0.0;   ///< mm; 0 means "match the crystal face"
  G4double fReadoutThickness = 1.0;  ///< mm
  G4double fReadoutEfficiency = 1.0; ///< flat PDE when no spectrum is supplied
  G4String fPdeFile;                 ///< optional wavelength/efficiency table

  G4String fResolvedFinish;          ///< filled in during Construct()
};

}  // namespace scint

#endif
