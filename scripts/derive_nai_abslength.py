"""Derive the NaI(Tl) bulk absorption length at 365 nm from a published datum.

THE DATUM.  Mao, Zhang and Zhu, IEEE Trans. Nucl. Sci. 55, 2425-2431 (2008),
Section III, state for their sample set:

    "the values of the cut-off wavelength, at which the transmittance data show
     50% of that at 800 nm, are 140 nm, 280 nm, 293 nm, 315 nm, 318 nm, 342 nm,
     358 nm, 365 nm and 390 nm for [LiF], CsI, [BaF2], BGO, CsI(Na), PWO,
     CsI(Tl), NaI(Tl) and LSO/LYSO respectively."

so for NaI(Tl):  T(365 nm) = 0.50 * T(800 nm).

THE SAMPLE.  Section II: "The NaI:Tl sample is a cylinder of 1.5 X0 long".
X0(NaI) = 9.49 g/cm^2 (PDG), rho = 3.667 g/cm^3  ->  X0 = 2.588 cm, L = 3.881 cm.

THE MODEL.  Their Eq. (1) is the transmittance of a slab with two parallel faces,
allowing multiple internal bounces and no internal absorption:

    T_theo = (1 - R)^2 / (1 - R^2) = (1 - R) / (1 + R),   R = ((n-1)/(n+1))^2

Adding bulk absorption with a = exp(-L / L_abs) between the faces:

    T = (1 - R)^2 a / (1 - R^2 a^2)

We invert this for a, hence for L_abs, at 365 nm.  At 800 nm the paper states
the measured transmittance "approaches the theoretical limit", so we take
T(800) = T_theo(800); the sensitivity to that assumption is reported below.

WHAT THIS IS NOT.  One anchor point, not a curve, and not a measurement we made.
It fixes the absorption length where the crystal is already strongly absorbing.
The wavelength dependence away from 365 nm is a model (see scint/optical.py).
"""

import math

# --- dispersion models (see scint/optical.py for provenance) ------------------

def n_li1976(lam_nm: float) -> float:
    l2 = (lam_nm / 1000.0) ** 2
    return math.sqrt(1 + 0.478 + 1.532 * l2 / (l2 - 0.170**2) + 4.27 * l2 / (l2 - 86.21**2))

def n_jellison_endpoints(lam_nm: float) -> float:
    # single-term Sellmeier through the two values quoted in the abstract of
    # Jellison et al., J. Appl. Phys. 111, 043521 (2012)
    A, C2 = 2.041276, 0.027186
    l2 = (lam_nm / 1000.0) ** 2
    return math.sqrt(1 + A * l2 / (l2 - C2))


def reflectance(n: float) -> float:
    return ((n - 1.0) / (n + 1.0)) ** 2

def t_theoretical(n: float) -> float:
    r = reflectance(n)
    return (1.0 - r) / (1.0 + r)

def abslength_from_transmittance(t_meas: float, n: float, length_mm: float) -> float:
    """Invert T = (1-R)^2 a / (1 - R^2 a^2) for L_abs, with a = exp(-L/L_abs)."""
    r = reflectance(n)
    # (T R^2) a^2 + (1-R)^2 a - T = 0
    A = t_meas * r * r
    B = (1.0 - r) ** 2
    C = -t_meas
    a = (-B + math.sqrt(B * B - 4 * A * C)) / (2 * A)
    if not 0.0 < a < 1.0:
        raise ValueError(f"unphysical transmission factor a={a}")
    return -length_mm / math.log(a)


X0_G_CM2_PDG = 9.49      # PDG, NaI
DENSITY = 3.667          # g/cm^3
LENGTH_MM = 10.0 * 1.5 * X0_G_CM2_PDG / DENSITY

def main() -> None:
    print(f"sample length L = 1.5 X0 = {LENGTH_MM:.3f} mm\n")
    for label, n_of in (("Li 1976", n_li1976), ("Jellison 2012 (endpoints)", n_jellison_endpoints)):
        n800, n365 = n_of(800.0), n_of(365.0)
        t800 = t_theoretical(n800)
        print(f"--- refractive index: {label}")
        print(f"    n(800) = {n800:.5f}   T_theo(800) = {t800:.5f}")
        print(f"    n(365) = {n365:.5f}   T_theo(365) = {t_theoretical(n365):.5f}")
        for assumed_t800, note in ((t800, "T(800) = theoretical limit"),
                                   (0.98 * t800, "T(800) 2 % below theory"),
                                   (0.95 * t800, "T(800) 5 % below theory")):
            t365 = 0.5 * assumed_t800
            l_abs = abslength_from_transmittance(t365, n365, LENGTH_MM)
            print(f"    T(365) = {t365:.5f}  ->  L_abs(365 nm) = {l_abs:7.2f} mm   [{note}]")
        print()

    # sensitivity to the sample length, i.e. to which X0 the authors used
    print("--- sensitivity to the assumed sample length")
    for x0 in (9.49, 9.67):
        L = 10.0 * 1.5 * x0 / DENSITY
        n365 = n_li1976(365.0)
        t365 = 0.5 * t_theoretical(n_li1976(800.0))
        print(f"    X0 = {x0:.2f} g/cm^2 -> L = {L:.3f} mm -> "
              f"L_abs(365 nm) = {abslength_from_transmittance(t365, n365, L):.2f} mm")

if __name__ == "__main__":
    main()
