#!/usr/bin/env python
"""Print the Layer 0 validation report.

Two independent checks:
  1. Shower quantities against the PDG Atomic and Nuclear Properties tables.
  2. Photon attenuation against Geant4, compared with a published compilation.

Run from the repository root:  ./.venv/bin/python scripts/validate_layer0.py
"""

from __future__ import annotations

from scint.attenuation import attenuation_lengths_cm
from scint.geant4 import version_info
from scint.materials import Material, attenuation_depth_cm

PDG = {
    "BGO": dict(formula="Bi4Ge3O12", rho=7.130, X0_cm=1.118, RM_cm=2.259, lamI_cm=22.32),
    "NaI": dict(formula="NaI", rho=3.667, X0_cm=2.588, RM_cm=4.105, lamI_cm=42.16),
}

BENCHMARK = [
    ("Pr:LuAG", "Lu3Al5O12", 6.73, 0.12, 0.65),
    ("Ce:GAGG", "Gd3Al2Ga3O12", 6.63, 0.18, 0.91),
    ("LaBr3:Ce", "LaBr3", 5.08, 0.33, 1.54),
    ("CeBr3", "CeBr3", 5.18, 0.31, 1.48),
    ("LYSO", "Lu1.8Y0.2SiO5", 7.20, 0.12, 0.67),
    ("CsI:Tl", "CsI", 4.51, 0.23, 1.25),
    ("NaI:Tl", "NaI", 3.67, 0.35, 1.76),
    ("BGO", "Bi4Ge3O12", 7.13, 0.08, 0.43),
]


def main() -> None:
    info = version_info()
    tag = info.get("tag") or info.get("config_version") or "unknown"
    print(f"Geant4: {tag}   (reported by geant4-config as {info.get('config_version')})")
    if info.get("is_prerelease"):
        print()
        print("!! This is a PRE-RELEASE Geant4 build. Results from it are for")
        print("!! development only and must not be published. Rerun in the pinned")
        print("!! container before quoting any number from this report.")
    print()
    print("=" * 78)
    print("1. Shower quantities vs PDG Atomic and Nuclear Properties (2024)")
    print("=" * 78)
    print(f"{'material':<8} {'quantity':<20} {'computed':>10} {'PDG':>10} {'deviation':>10}")
    for name, ref in PDG.items():
        material = Material(name, ref["formula"], ref["rho"])
        for label, got, expected in [
            ("radiation length", material.radiation_length_cm, ref["X0_cm"]),
            ("Moliere radius", material.moliere_radius_cm, ref["RM_cm"]),
            ("interaction length", material.nuclear_interaction_length_cm, ref["lamI_cm"]),
        ]:
            dev = (got - expected) / expected * 100
            print(f"{name:<8} {label:<20} {got:10.3f} {expected:10.3f} {dev:+9.1f}%")

    print()
    print("=" * 78)
    print("2. 88% attenuation depth: Geant4 vs arXiv:2505.06929 Table 1")
    print("=" * 78)
    print(f"{'material':<10} {'Zeff':>6} {'X0/cm':>7} | "
          f"{'100 keV':>8} {'ref':>5} {'dev':>7} | {'200 keV':>8} {'ref':>5} {'dev':>7}")
    for name, formula, rho, ref100, ref200 in BENCHMARK:
        material = Material(name, formula, rho)
        lengths = attenuation_lengths_cm(material, [0.1, 0.2])
        d100 = attenuation_depth_cm(lengths[0.1], 0.88)
        d200 = attenuation_depth_cm(lengths[0.2], 0.88)
        dev100 = (d100 - ref100) / ref100 * 100
        dev200 = (d200 - ref200) / ref200 * 100
        flag = "  <-- see note in benchmark CSV" if max(abs(dev100), abs(dev200)) > 10 else ""
        print(f"{name:<10} {material.z_eff:6.1f} {material.radiation_length_cm:7.3f} | "
              f"{d100:8.3f} {ref100:5.2f} {dev100:+6.1f}% | "
              f"{d200:8.3f} {ref200:5.2f} {dev200:+6.1f}%{flag}")
    print()
    print("Z_eff uses the power-law convention with exponent 3.5.")


if __name__ == "__main__":
    main()
