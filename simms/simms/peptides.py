"""In-silico proteolytic digestion and theoretical peptide fragment spectra.

Replaces the CLI-accessible part of what OpenMS's removed MSSimulator /
TheoreticalSpectrumGenerator provided, in pure Python: tryptic digestion of
FASTA entries and b/y (optionally a) fragment ion spectra at configurable
charge states.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

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
