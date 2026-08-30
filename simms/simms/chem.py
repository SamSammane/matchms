"""Molecular formulas, monoisotopic masses and coarse isotope patterns.

Self-contained (numpy only) so that spectrum simulation does not depend on
heavyweight chemistry toolkits. Masses are monoisotopic, abundances are
natural terrestrial values.
"""

from __future__ import annotations

import re
from typing import Dict, List, Tuple

PROTON = 1.00727646688
ELECTRON = 0.00054857990

# element -> list of (isotope mass, natural abundance), sorted by mass.
ISOTOPES: Dict[str, List[Tuple[float, float]]] = {
    "H": [(1.0078250319, 0.999885), (2.0141017780, 0.000115)],
    "C": [(12.0, 0.9893), (13.0033548378, 0.0107)],
    "N": [(14.0030740052, 0.99632), (15.0001088984, 0.00368)],
    "O": [(15.9949146221, 0.99757), (16.9991315000, 0.00038), (17.9991604000, 0.00205)],
    "S": [(31.9720706900, 0.9493), (32.9714585000, 0.0076), (33.9678668300, 0.0429), (35.9670808800, 0.0002)],
    "P": [(30.9737615100, 1.0)],
    "F": [(18.9984032200, 1.0)],
    "Cl": [(34.9688527100, 0.7578), (36.9659026000, 0.2422)],
    "Br": [(78.9183376000, 0.5069), (80.9162910000, 0.4931)],
    "I": [(126.9044680000, 1.0)],
    "Si": [(27.9769265327, 0.92230), (28.9764947200, 0.04683), (29.9737702200, 0.03087)],
    "Na": [(22.9897692800, 1.0)],
    "K": [(38.9637069000, 0.932581), (39.9639986700, 0.000117), (40.9618259700, 0.067302)],
    "Li": [(6.0151223000, 0.0759), (7.0160040000, 0.9241)],
    "B": [(10.0129370000, 0.199), (11.0093055000, 0.801)],
    "Se": [(73.9224766000, 0.0089), (75.9192141000, 0.0937), (76.9199146000, 0.0763),
           (77.9173095000, 0.2377), (79.9165218000, 0.4961), (81.9167000000, 0.0873)],
    "Fe": [(53.9396148000, 0.05845), (55.9349421000, 0.91754), (56.9353987000, 0.02119), (57.9332805000, 0.00282)],
}

MONOISOTOPIC: Dict[str, float] = {
    element: max(isos, key=lambda iso: iso[1])[0] for element, isos in ISOTOPES.items()
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")

# adduct -> (mass delta added to neutral M, charge). Sign of charge sets polarity.
ADDUCTS: Dict[str, Tuple[float, int]] = {
    "[M+H]+": (PROTON, 1),
    "[M+2H]2+": (2 * PROTON, 2),
    "[M+3H]3+": (3 * PROTON, 3),
    "[M+Na]+": (22.98976928 - ELECTRON, 1),
    "[M+K]+": (38.9637069 - ELECTRON, 1),
    "[M+NH4]+": (18.03437413, 1),
    "[M-H2O+H]+": (PROTON - 18.0105646863, 1),
    "[M]+": (-ELECTRON, 1),
    "[M-H]-": (-PROTON, -1),
    "[M-2H]2-": (-2 * PROTON, -2),
    "[M+Cl]-": (34.96885271 + ELECTRON, -1),
    "[M+HCOO]-": (44.99820285 + PROTON, -1),
    "[M]-": (ELECTRON, -1),
}


def parse_formula(formula: str) -> Dict[str, int]:
    """Parse a plain molecular formula (e.g. ``C16H30N2O3``) into element counts."""
    formula = formula.strip()
    if not formula:
        raise ValueError("empty molecular formula")
    counts: Dict[str, int] = {}
    pos = 0
    for match in _FORMULA_TOKEN.finditer(formula):
        if match.start() != pos:
            raise ValueError(f"cannot parse formula {formula!r} at position {pos}")
        pos = match.end()
        element, digits = match.groups()
        if element not in ISOTOPES:
            raise ValueError(f"unsupported element {element!r} in formula {formula!r}")
        counts[element] = counts.get(element, 0) + (int(digits) if digits else 1)
    if pos != len(formula):
        raise ValueError(f"cannot parse formula {formula!r} at position {pos}")
    return counts


def monoisotopic_mass(formula: str) -> float:
    """Monoisotopic (most abundant isotope) mass of a neutral formula."""
    return sum(MONOISOTOPIC[el] * n for el, n in parse_formula(formula).items())


def adduct_mz(formula: str, adduct: str = "[M+H]+") -> float:
    """m/z of a formula for a named adduct."""
    if adduct not in ADDUCTS:
        raise ValueError(f"unsupported adduct {adduct!r}; known: {sorted(ADDUCTS)}")
    delta, charge = ADDUCTS[adduct]
    return (monoisotopic_mass(formula) + delta) / abs(charge)


# A coarse pattern maps integer nominal-mass shift -> (abundance, average mass).
_Pattern = Dict[int, Tuple[float, float]]


def _convolve(a: _Pattern, b: _Pattern, prune: float) -> _Pattern:
    out: Dict[int, Tuple[float, float]] = {}
    for shift_a, (ab_a, mass_a) in a.items():
        for shift_b, (ab_b, mass_b) in b.items():
            abundance = ab_a * ab_b
            if abundance < prune:
                continue
            shift = shift_a + shift_b
            mass = mass_a + mass_b
            if shift in out:
                prev_ab, prev_mass = out[shift]
                total = prev_ab + abundance
                out[shift] = (total, (prev_mass * prev_ab + mass * abundance) / total)
            else:
                out[shift] = (abundance, mass)
    return out


def _element_pattern(element: str) -> _Pattern:
    isos = ISOTOPES[element]
    base_nominal = round(isos[0][0])
    return {round(mass) - base_nominal: (abundance, mass) for mass, abundance in isos}


def isotope_pattern(formula: str, prune: float = 1e-5, max_peaks: int = 12) -> List[Tuple[float, float]]:
    """Coarse isotope pattern of a neutral formula.

    Returns a list of (mass, relative abundance) with the most abundant
    peak normalized to 1.0, sorted by mass.
    """
    counts = parse_formula(formula)
    pattern: _Pattern = {0: (1.0, 0.0)}
    for element, n in counts.items():
        elem = _element_pattern(element)
        # exponentiation by squaring on the convolution
        power = elem
        remaining = n
        while remaining:
            if remaining & 1:
                pattern = _convolve(pattern, power, prune)
            remaining >>= 1
            if remaining:
                power = _convolve(power, power, prune)
    peaks = [(mass, abundance) for abundance, mass in pattern.values()]
    peaks.sort(key=lambda p: p[0])
    top = max(abundance for _, abundance in peaks)
    peaks = [(mass, abundance / top) for mass, abundance in peaks if abundance / top > prune]
    if len(peaks) > max_peaks:
        peaks = sorted(peaks, key=lambda p: p[1], reverse=True)[:max_peaks]
        peaks.sort(key=lambda p: p[0])
    return peaks


def isotope_envelope_mz(formula: str, adduct: str = "[M+H]+", prune: float = 1e-4,
                        max_peaks: int = 8) -> List[Tuple[float, float]]:
    """Isotope pattern shifted onto the m/z axis for a given adduct."""
    if adduct not in ADDUCTS:
        raise ValueError(f"unsupported adduct {adduct!r}")
    delta, charge = ADDUCTS[adduct]
    z = abs(charge)
    return [((mass + delta) / z, abundance)
            for mass, abundance in isotope_pattern(formula, prune=prune, max_peaks=max_peaks)]


# Senko's averagine: average elemental composition of peptides per
# 111.1254 Da of monoisotopic mass, used to model isotope envelopes when
# only a mass (not a formula) is known.
_AVERAGINE = {"C": 4.9384, "H": 7.7583, "N": 1.3577, "O": 1.4773, "S": 0.0417}
_AVERAGINE_UNIT = 111.1254


def averagine_formula(mass: float) -> Dict[str, int]:
    """Estimate an elemental composition for a peptide-like neutral mass."""
    if mass <= 0:
        raise ValueError("mass must be positive")
    scale = mass / _AVERAGINE_UNIT
    counts = {el: max(0, round(n * scale)) for el, n in _AVERAGINE.items()}
    counts = {el: n for el, n in counts.items() if n > 0}
    # top up hydrogens so the model mass tracks the requested mass
    model_mass = sum(MONOISOTOPIC[el] * n for el, n in counts.items())
    counts["H"] = counts.get("H", 0) + max(0, round((mass - model_mass) / MONOISOTOPIC["H"]))
    return counts


def averagine_envelope_mz(neutral_mass: float, charge: int = 1, prune: float = 1e-3,
                          max_peaks: int = 8) -> List[Tuple[float, float]]:
    """Approximate isotope envelope (m/z, rel abundance) for a neutral mass
    observed at the given positive charge, via the averagine model."""
    if charge < 1:
        raise ValueError("charge must be >= 1")
    counts = averagine_formula(neutral_mass)
    formula = "".join(f"{el}{n}" for el, n in sorted(counts.items()))
    pattern = isotope_pattern(formula, prune=prune, max_peaks=max_peaks)
    if not pattern:
        return [((neutral_mass + charge * PROTON) / charge, 1.0)]
    # anchor the monoisotopic peak on the requested mass; keep isotope spacings
    mono = pattern[0][0]
    return [((neutral_mass + (mass - mono) + charge * PROTON) / charge, abundance)
            for mass, abundance in pattern]
