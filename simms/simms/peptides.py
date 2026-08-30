"""In-silico proteolytic digestion and theoretical peptide fragment spectra.

Replaces the CLI-accessible part of what OpenMS's removed MSSimulator /
TheoreticalSpectrumGenerator provided, in pure Python: tryptic digestion of
FASTA entries and b/y (optionally a) fragment ion spectra at configurable
charge states.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import numpy as np

from .chem import PROTON

WATER = 18.0105646863
AMMONIA = 17.0265491015
CO = 27.9949146221

# monoisotopic residue masses (peptide-bond residues, i.e. minus water)
RESIDUE_MASS: Dict[str, float] = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406,
    "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
    "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841,
    "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}

ENZYMES = {
    # regex matching the cleavage position (after the match start)
    "trypsin": re.compile(r"(?<=[KR])(?!P)"),
    "trypsin/p": re.compile(r"(?<=[KR])"),
    "lys-c": re.compile(r"(?<=K)"),
    "arg-c": re.compile(r"(?<=R)"),
    "chymotrypsin": re.compile(r"(?<=[FWYL])(?!P)"),
    "no cleavage": None,
}


def peptide_mass(sequence: str) -> float:
    """Monoisotopic neutral mass of a peptide."""
    try:
        return sum(RESIDUE_MASS[aa] for aa in sequence) + WATER
    except KeyError as err:
        raise ValueError(f"unknown amino acid {err.args[0]!r} in {sequence!r}") from None


def precursor_mz(sequence: str, charge: int) -> float:
    if charge < 1:
        raise ValueError("charge must be >= 1")
    return (peptide_mass(sequence) + charge * PROTON) / charge


def digest(sequence: str, enzyme: str = "trypsin", missed_cleavages: int = 0,
           min_length: int = 6, max_length: int = 40) -> List[str]:
    """Digest a protein sequence into peptides."""
    enzyme = enzyme.lower()
    if enzyme not in ENZYMES:
        raise ValueError(f"unsupported enzyme {enzyme!r}; known: {sorted(ENZYMES)}")
    sequence = sequence.strip().upper()
    pattern = ENZYMES[enzyme]
    if pattern is None:
        fragments = [sequence]
    else:
        fragments = [f for f in pattern.split(sequence) if f]
    peptides = []
    for start in range(len(fragments)):
        for extra in range(missed_cleavages + 1):
            end = start + extra + 1
            if end > len(fragments):
                break
            peptide = "".join(fragments[start:end])
            if min_length <= len(peptide) <= max_length and all(aa in RESIDUE_MASS for aa in peptide):
                peptides.append(peptide)
    # stable de-duplication
    seen = set()
    unique = []
    for peptide in peptides:
        if peptide not in seen:
            seen.add(peptide)
            unique.append(peptide)
    return unique


@dataclass
class FragmentIon:
    mz: float
    intensity: float
    ion_type: str  # e.g. "b3", "y5++"


def fragment_ions(sequence: str, ion_types: str = "by", max_fragment_charge: int = 1) -> List[FragmentIon]:
    """Theoretical fragment ions for a peptide.

    Intensities follow a simple deterministic heuristic: y ions are stronger
    than b ions, mid-sequence fragments stronger than termini — enough for
    realistic-looking, reproducible test spectra.
    """
    if not sequence:
        raise ValueError("empty peptide sequence")
    masses = [RESIDUE_MASS[aa] for aa in sequence]
    n = len(masses)
    prefix = np.cumsum(masses)
    total = prefix[-1]
    ions: List[FragmentIon] = []
    for i in range(1, n):  # fragment lengths 1..n-1
        # bell-shaped positional weight peaking mid-sequence
        weight = 0.25 + 0.75 * np.sin(np.pi * i / n)
        b_neutral = prefix[i - 1]
        y_neutral = total - prefix[i - 1] + WATER
        a_neutral = b_neutral - CO
        for z in range(1, max_fragment_charge + 1):
            suffix = "+" * z
            if "b" in ion_types:
                ions.append(FragmentIon((b_neutral + z * PROTON) / z, 55.0 * weight / z, f"b{i}{suffix}"))
            if "y" in ion_types:
                ions.append(FragmentIon((y_neutral + z * PROTON) / z, 100.0 * weight / z, f"y{n - i}{suffix}"))
            if "a" in ion_types:
                ions.append(FragmentIon((a_neutral + z * PROTON) / z, 20.0 * weight / z, f"a{i}{suffix}"))
    ions.sort(key=lambda ion: ion.mz)
    return ions


# immonium ion m/z for residues that commonly produce them under CID/HCD
IMMONIUM = {
    "L": 86.0964, "I": 86.0964, "P": 70.0651, "V": 72.0808, "F": 120.0808,
    "Y": 136.0757, "W": 159.0917, "H": 110.0713, "R": 129.1135, "K": 101.1073,
    "Q": 101.0709, "E": 102.0550, "M": 104.0528, "N": 87.0553, "D": 88.0393,
}

_BASIC_WEIGHT = {"R": 1.0, "K": 0.8, "H": 0.5}


def proton_mobility(sequence: str, precursor_charge: int) -> str:
    """Mobile-proton classification: 'mobile', 'partial', or 'nonmobile'.

    Charges beyond the count of basic residues are mobile and drive
    charge-directed backbone fragmentation; sequestered protons favor
    charge-remote pathways (enhanced cleavage C-terminal to D/E).
    """
    basic = sum(_BASIC_WEIGHT.get(aa, 0.0) for aa in sequence)
    if precursor_charge > basic + 0.5:
        return "mobile"
    if precursor_charge > basic - 0.5:
        return "partial"
    return "nonmobile"


def realistic_fragment_ions(sequence: str, precursor_charge: int = 2,
                            collision_energy: float = 25.0,
                            max_fragment_charge: Optional[int] = None) -> List[FragmentIon]:
    """Fragment ions with intensities from a mobile-proton heuristic model.

    Captures the qualitative rules of peptide CID/HCD spectra:

    - y ions dominate b ions; b1 is essentially absent, b2 is prominent;
    - cleavage N-terminal to proline strongly enhances the corresponding y ion;
    - under non-mobile conditions, cleavage C-terminal to Asp/Glu is enhanced;
    - higher collision energy shifts intensity toward low-mass fragments;
    - abundant fragments carry M+1 isotope peaks (averagine carbon count);
    - fragments containing S/T/E/D show -H2O and R/K/N/Q show -NH3 losses;
    - immonium ions appear at high collision energy.
    """
    if not sequence:
        raise ValueError("empty peptide sequence")
    if max_fragment_charge is None:
        max_fragment_charge = min(2, max(1, precursor_charge - 1))
    masses = [RESIDUE_MASS[aa] for aa in sequence]
    n = len(masses)
    prefix = np.cumsum(masses)
    total = prefix[-1]
    precursor_mass = total + WATER
    mobility = proton_mobility(sequence, precursor_charge)
    ce_ratio = collision_energy / 25.0

    ions: List[FragmentIon] = []

    def add(mz: float, intensity: float, label: str, fragment_mass: float,
            residues: str) -> None:
        if intensity < 0.05:
            return
        # collision energy: high CE erodes large fragments, low CE keeps them
        size_fraction = fragment_mass / precursor_mass
        intensity *= float(np.exp(-(ce_ratio - 1.0) * size_fraction * 1.5))
        if intensity < 0.05:
            return
        ions.append(FragmentIon(mz, intensity, label))
        # M+1 isotope for abundant, heavier fragments (averagine carbons)
        if intensity > 3.0 and fragment_mass > 250:
            n_carbon = fragment_mass / 111.1254 * 4.9384
            iso_frac = min(0.85, n_carbon * 0.0107 + fragment_mass / 111.1254 * 1.3577 * 0.0037)
            ions.append(FragmentIon(mz + 1.00336 / _label_charge(label),
                                    intensity * iso_frac, label + "+i"))
        # neutral losses
        if intensity > 2.0:
            if any(aa in residues for aa in "STED"):
                ions.append(FragmentIon(mz - WATER / _label_charge(label),
                                        intensity * 0.15, label + "-H2O"))
            if any(aa in residues for aa in "RKNQ"):
                ions.append(FragmentIon(mz - AMMONIA / _label_charge(label),
                                        intensity * 0.12, label + "-NH3"))

    def _label_charge(label: str) -> int:
        return max(1, label.count("+")) if "+" in label else 1

    for i in range(1, n):
        weight = 0.25 + 0.75 * float(np.sin(np.pi * i / n))
        b_res, y_res = sequence[:i], sequence[i:]
        b_neutral = float(prefix[i - 1])
        y_neutral = float(total - prefix[i - 1] + WATER)

        b_base, y_base = 40.0 * weight, 100.0 * weight
        if sequence[i] == "P":       # X|P: proline effect, strong y / weak b
            y_base *= 3.0
            b_base *= 0.3
        if sequence[i - 1] == "P":   # P|X cleavage is suppressed
            y_base *= 0.3
            b_base *= 0.3
        if mobility == "nonmobile" and sequence[i - 1] in "DE":
            b_base *= 2.5            # charge-remote D/E enhancement
            y_base *= 2.0
        elif mobility == "partial" and sequence[i - 1] in "DE":
            b_base *= 1.5
        if i == 1:
            b_base *= 0.02           # b1 practically never observed
        elif i == 2:
            b_base *= 1.5            # b2 prominent
        if mobility == "nonmobile":
            b_base *= 0.6            # overall weaker backbone fragmentation
            y_base *= 0.6

        for z in range(1, max_fragment_charge + 1):
            # multiply-charged fragments need enough residues to hold charge
            charge_capacity = 1 + sum(1 for aa in y_res if aa in "RKH")
            z_scale = 1.0 if z == 1 else (0.35 if z <= charge_capacity else 0.02)
            suffix = "+" * z
            add((y_neutral + z * PROTON) / z, y_base * z_scale / z,
                f"y{n - i}{suffix}", y_neutral, y_res)
            charge_capacity_b = 1 + sum(1 for aa in b_res if aa in "RKH")
            z_scale_b = 1.0 if z == 1 else (0.3 if z <= charge_capacity_b else 0.02)
            add((b_neutral + z * PROTON) / z, b_base * z_scale_b / z,
                f"b{i}{suffix}", b_neutral, b_res)

    # precursor-related ions (weak survivor + water loss)
    prec_mz = (precursor_mass + precursor_charge * PROTON) / precursor_charge
    survivor = max(0.0, 12.0 * (1.6 - ce_ratio))
    if survivor >= 0.05:
        ions.append(FragmentIon(prec_mz, survivor, "M"))
        if any(aa in sequence for aa in "STED"):
            ions.append(FragmentIon(prec_mz - WATER / precursor_charge,
                                    survivor * 0.4, "M-H2O"))

    # immonium ions show up under energetic activation
    if collision_energy >= 25:
        for aa in sorted(set(sequence) & set(IMMONIUM)):
            ions.append(FragmentIon(IMMONIUM[aa], 3.0 * ce_ratio, f"imm({aa})"))

    top = max(ion.intensity for ion in ions)
    ions = [FragmentIon(ion.mz, ion.intensity * 100.0 / top, ion.ion_type)
            for ion in ions]
    ions.sort(key=lambda ion: ion.mz)
    return ions


def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)
