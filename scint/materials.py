"""Layer 0: deterministic material properties derived from stoichiometry alone.

Nothing in this module is a simulation. Given a chemical formula and a density,
every quantity here follows from exact arithmetic plus published closed-form
expressions, so results are reproducible to machine precision.

Atomic masses are read from ``data/derived/atomic_data.csv``, which is exported
from Geant4's own NIST element database by ``sim/tools/g4data.cc``. That keeps
the screening layer and the transport simulation on a single source of truth and
avoids vendoring Geant4's licence-restricted cross-section files.

Formulae used (all from the Particle Data Group, *Review of Particle Physics*,
"Passage of particles through matter"):

* Radiation length, Dahl approximation::

      X0 = 716.4 g/cm^2 * A / (Z (Z+1) ln(287 / sqrt(Z)))

  combined for mixtures as ``1/X0 = sum_j w_j / X0_j`` over mass fractions w_j.

* Critical energy for solids and liquids::  Ec = 610 MeV / (Z + 1.24)

* Moliere radius::  RM = X0 * Es / Ec  with  Es = 21.2052 MeV,
  combined for mixtures as ``1/RM = (1/Es) sum_j w_j Ec_j / X0_j``.

* Nuclear interaction length, approximated as ``lambda_I = 35.0 * A^(1/3) g/cm^2``,
  combined as ``1/lambda = sum_j w_j / lambda_j``. Measured against PDG for BGO
  and NaI this holds to within 2 %, but it remains an approximation and can be
  worse for other compositions; take Geant4 values where the number matters.

* Effective atomic number, power-law definition::

      Z_eff = (sum_i alpha_i Z_i^m)^(1/m)

  with alpha_i the fraction of total electrons contributed by element i. There is
  no single accepted exponent m; values between 2.94 and 3.5 appear in the
  literature. The exponent is therefore an explicit parameter and is reported
  alongside every value it produces.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

# PDG constants
RADIATION_LENGTH_CONSTANT_G_PER_CM2 = 716.4
ES_MEV = 21.2052
CRITICAL_ENERGY_NUMERATOR_MEV = 610.0  # solids and liquids
CRITICAL_ENERGY_OFFSET = 1.24
NUCLEAR_INTERACTION_CONSTANT = 35.0

DEFAULT_ZEFF_EXPONENT = 3.5

_REPO_ROOT = Path(__file__).resolve().parent.parent
_ATOMIC_DATA = _REPO_ROOT / "data" / "derived" / "atomic_data.csv"

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)|(\()|(\))(\d*\.?\d*)")


class FormulaError(ValueError):
    """Raised when a chemical formula cannot be parsed unambiguously."""


@lru_cache(maxsize=1)
def atomic_masses() -> dict[str, float]:
    """Return {symbol: standard atomic mass in amu} from the Geant4 NIST export.

    Raises FileNotFoundError with a pointer to the generating command if the
    export has not been produced yet -- the data is deliberately not committed,
    because it is derived and must stay in step with the installed Geant4.
    """
    if not _ATOMIC_DATA.exists():
        raise FileNotFoundError(
            f"{_ATOMIC_DATA} is missing. Generate it with:\n"
            "    ./build/sim/g4data elements --out data/derived/atomic_data.csv"
        )
    masses: dict[str, float] = {}
    with _ATOMIC_DATA.open() as fh:
        rows = csv.DictReader(line for line in fh if not line.startswith("#"))
        for row in rows:
            masses[row["symbol"]] = float(row["atomic_mass_amu"])
    return masses


@lru_cache(maxsize=1)
def atomic_numbers() -> dict[str, int]:
    """Return {symbol: Z} from the Geant4 NIST export."""
    if not _ATOMIC_DATA.exists():
        raise FileNotFoundError(
            f"{_ATOMIC_DATA} is missing. Generate it with:\n"
            "    ./build/sim/g4data elements --out data/derived/atomic_data.csv"
        )
    numbers: dict[str, int] = {}
    with _ATOMIC_DATA.open() as fh:
        rows = csv.DictReader(line for line in fh if not line.startswith("#"))
        for row in rows:
            numbers[row["symbol"]] = int(row["z"])
    return numbers


def parse_formula(formula: str) -> dict[str, float]:
    """Parse a chemical formula into {element symbol: atom count}.

    Supports nested parentheses and fractional subscripts, so multicomponent
    garnets write naturally::

        >>> sorted(parse_formula("Gd3Al2Ga3O12").items())
        [('Al', 2.0), ('Ga', 3.0), ('Gd', 3.0), ('O', 12.0)]
        >>> sorted(parse_formula("(Lu0.5Gd2.5)(Al2Ga3)O12").items())
        [('Al', 2.0), ('Ga', 3.0), ('Gd', 2.5), ('Lu', 0.5), ('O', 12.0)]

    An omitted subscript means 1, matching chemical convention.
    """
    text = formula.replace(" ", "")
    if not text:
        raise FormulaError("empty formula")

    stack: list[dict[str, float]] = [{}]
    pos = 0
    while pos < len(text):
        match = _TOKEN.match(text, pos)
        if match is None or match.end() == pos:
            raise FormulaError(f"cannot parse {formula!r} at position {pos}: {text[pos:]!r}")
        symbol, count, open_paren, close_paren, group_count = match.groups()
        if open_paren:
            stack.append({})
        elif close_paren:
            if len(stack) == 1:
                raise FormulaError(f"unbalanced ')' in {formula!r}")
            group = stack.pop()
            multiplier = float(group_count) if group_count else 1.0
            for element, n in group.items():
                stack[-1][element] = stack[-1].get(element, 0.0) + n * multiplier
        else:
            n = float(count) if count else 1.0
            stack[-1][symbol] = stack[-1].get(symbol, 0.0) + n
        pos = match.end()

    if len(stack) != 1:
        raise FormulaError(f"unbalanced '(' in {formula!r}")

    known = atomic_numbers()
    unknown = sorted(set(stack[0]) - set(known))
    if unknown:
        raise FormulaError(f"unknown element symbol(s) in {formula!r}: {', '.join(unknown)}")
    if not stack[0]:
        raise FormulaError(f"no elements found in {formula!r}")
    return stack[0]


@dataclass(frozen=True)
class Material:
    """A compound with a known stoichiometry and measured density.

    ``density_g_cm3`` is an input, not a prediction: it is either a measured
    value or one computed from a crystallographic unit cell elsewhere.
    ``zeff_exponent`` is carried on the object so that any Z_eff it reports can
    always be traced to the convention that produced it.
    """

    name: str
    formula: str
    density_g_cm3: float
    zeff_exponent: float = DEFAULT_ZEFF_EXPONENT
    composition: dict[str, float] = field(init=False)

    def __post_init__(self) -> None:
        if self.density_g_cm3 <= 0:
            raise ValueError(f"{self.name}: density must be positive, got {self.density_g_cm3}")
        object.__setattr__(self, "composition", parse_formula(self.formula))

    # ---- basic stoichiometry -------------------------------------------------

    @property
    def molar_mass(self) -> float:
        """Formula mass in amu (equivalently g/mol)."""
        masses = atomic_masses()
        return sum(n * masses[el] for el, n in self.composition.items())

    @property
    def mass_fractions(self) -> dict[str, float]:
        """Mass fraction per element, summing to 1."""
        masses = atomic_masses()
        total = self.molar_mass
        return {el: n * masses[el] / total for el, n in self.composition.items()}

    @property
    def electron_fractions(self) -> dict[str, float]:
        """Fraction of the compound's electrons contributed by each element."""
        numbers = atomic_numbers()
        per_element = {el: n * numbers[el] for el, n in self.composition.items()}
        total = sum(per_element.values())
        return {el: v / total for el, v in per_element.items()}

    @property
    def electron_density_per_cm3(self) -> float:
        """Electron number density, in electrons per cm^3."""
        avogadro = 6.02214076e23  # exact, SI definition
        numbers = atomic_numbers()
        electrons_per_formula = sum(n * numbers[el] for el, n in self.composition.items())
        return avogadro * self.density_g_cm3 * electrons_per_formula / self.molar_mass

    # ---- derived shower / interaction quantities -----------------------------

    @property
    def z_eff(self) -> float:
        """Effective atomic number under the power-law convention.

        The exponent is ``self.zeff_exponent``; report it with the value.
        """
        m = self.zeff_exponent
        numbers = atomic_numbers()
        return sum(a * numbers[el] ** m for el, a in self.electron_fractions.items()) ** (1.0 / m)

    @property
    def radiation_length_g_cm2(self) -> float:
        """Radiation length in g/cm^2 (PDG Dahl approximation, mixture rule)."""
        masses, numbers = atomic_masses(), atomic_numbers()
        inverse = 0.0
        for el, w in self.mass_fractions.items():
            z, a = numbers[el], masses[el]
            x0 = RADIATION_LENGTH_CONSTANT_G_PER_CM2 * a / (z * (z + 1) * _ln(287.0 / z**0.5))
            inverse += w / x0
        return 1.0 / inverse

    @property
    def radiation_length_cm(self) -> float:
        """Radiation length in cm."""
        return self.radiation_length_g_cm2 / self.density_g_cm3

    @property
    def moliere_radius_g_cm2(self) -> float:
        """Moliere radius in g/cm^2 (PDG mixture rule)."""
        masses, numbers = atomic_masses(), atomic_numbers()
        inverse = 0.0
        for el, w in self.mass_fractions.items():
            z, a = numbers[el], masses[el]
            x0 = RADIATION_LENGTH_CONSTANT_G_PER_CM2 * a / (z * (z + 1) * _ln(287.0 / z**0.5))
            ec = CRITICAL_ENERGY_NUMERATOR_MEV / (z + CRITICAL_ENERGY_OFFSET)
            inverse += w * ec / x0
        return ES_MEV / inverse

    @property
    def moliere_radius_cm(self) -> float:
        """Moliere radius in cm."""
        return self.moliere_radius_g_cm2 / self.density_g_cm3

    @property
    def nuclear_interaction_length_cm(self) -> float:
        """Approximate nuclear interaction length in cm.

        Uses ``lambda = 35 A^(1/3) g/cm^2``. Agrees with PDG to within 2 % for
        BGO and NaI (see tests/test_materials.py), but it is an approximation.
        """
        masses = atomic_masses()
        inverse = sum(
            w / (NUCLEAR_INTERACTION_CONSTANT * masses[el] ** (1.0 / 3.0))
            for el, w in self.mass_fractions.items()
        )
        return 1.0 / inverse / self.density_g_cm3

    # ---- interop -------------------------------------------------------------

    def massfrac_arg(self, precision: int = 9) -> str:
        """Format mass fractions for the ``g4data attenuation --massfrac`` option.

        Renormalises after rounding so the string always sums to 1 within the
        tolerance g4data enforces.
        """
        fractions = self.mass_fractions
        items = sorted(fractions.items())
        rounded = {el: round(v, precision) for el, v in items}
        drift = round(1.0 - sum(rounded.values()), precision)
        if drift:
            heaviest = max(rounded, key=lambda el: rounded[el])
            rounded[heaviest] = round(rounded[heaviest] + drift, precision)
        return ",".join(f"{el}:{rounded[el]:.{precision}f}" for el, _ in items)

    def summary(self) -> dict[str, float | str]:
        """A flat, report-friendly dictionary of every Layer 0 quantity."""
        return {
            "name": self.name,
            "formula": self.formula,
            "density_g_cm3": self.density_g_cm3,
            "molar_mass_amu": self.molar_mass,
            "z_eff": self.z_eff,
            "z_eff_exponent": self.zeff_exponent,
            "electron_density_per_cm3": self.electron_density_per_cm3,
            "radiation_length_cm": self.radiation_length_cm,
            "moliere_radius_cm": self.moliere_radius_cm,
            "nuclear_interaction_length_cm": self.nuclear_interaction_length_cm,
        }


def attenuation_depth_cm(attenuation_length_cm: float, fraction: float) -> float:
    """Depth at which ``fraction`` of an incident beam has interacted.

    Inverts ``1 - exp(-z/lambda) = fraction``. With fraction = 0.88 this
    reproduces the "88 % attenuation depth" convention used in the benchmark
    compilation in ``data/benchmark_scintillators.csv``.
    """
    if not 0.0 < fraction < 1.0:
        raise ValueError(f"fraction must lie strictly between 0 and 1, got {fraction}")
    return -attenuation_length_cm * _ln(1.0 - fraction)


def _ln(x: float) -> float:
    from math import log

    return log(x)
