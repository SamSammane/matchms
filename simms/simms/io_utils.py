"""Loading and saving spectra across formats, built on matchms.

Read: mgf, msp, mzML, mzXML, json, pickle (matchms dispatch) plus MassBank
record .txt files/directories (this package). Write: mgf, msp, json, pickle
(matchms) plus MassBank record directories.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import List, Optional, Sequence

from matchms import Spectrum
from matchms.importing.load_spectra import load_spectra

from . import massbank

READ_EXTENSIONS = {".mgf", ".msp", ".mzml", ".mzxml", ".json", ".pickle", ".txt"}
WRITE_EXTENSIONS = {".mgf", ".msp", ".json", ".pickle"}

_INT_KEYS = ("ms_level", "variant")
_FLOAT_KEYS = ("retention_time", "parent_mass", "precursor_mz")


def _coerce_numeric_metadata(spectrum: Optional[Spectrum]) -> Optional[Spectrum]:
    """Text formats (mgf/msp) stringify metadata; restore numeric types."""
    if spectrum is None:
        return None
    for key, caster in [(k, int) for k in _INT_KEYS] + [(k, float) for k in _FLOAT_KEYS]:
        value = spectrum.get(key)
        if isinstance(value, str):
            try:
                spectrum.set(key, caster(float(value)) if caster is int else caster(value))
            except ValueError:
                pass
    return spectrum


def load_any(path: str) -> List[Spectrum]:
    """Load all spectra from one file or MassBank record directory."""
    p = Path(path)
    if p.is_dir():
        spectra = []
        for record_path in sorted(p.glob("**/MSBNK-*.txt")):
            spectrum = massbank.record_to_spectrum(massbank.load_record(str(record_path)))
            if spectrum is not None:
                spectra.append(spectrum)
        return spectra
    suffix = p.suffix.lower()
    if suffix == ".txt":
        spectrum = massbank.record_to_spectrum(massbank.load_record(str(p)))
        return [spectrum] if spectrum is not None else []
    spectra = [_coerce_numeric_metadata(s) for s in load_spectra(str(p))]
    return [s for s in spectra if s is not None]


def load_many(paths: Sequence[str]) -> List[Spectrum]:
    spectra: List[Spectrum] = []
    for path in paths:
        spectra.extend(load_any(path))
    return spectra


def save_any(spectra: Sequence[Spectrum], path: str,
             export_style: str = "matchms") -> int:
    """Save spectra to mgf/msp/json/pickle, or a directory of MassBank records."""
    spectra = [s for s in spectra if s is not None]
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == "" or suffix == ".massbank":
        directory = p if suffix == "" else p.with_suffix("")
        directory.mkdir(parents=True, exist_ok=True)
        for i, spectrum in enumerate(spectra):
            accession = f"MSBNK-SIMMS-SIM{i + 1:06d}"
            text = massbank.spectrum_to_record_text(spectrum, accession)
            (directory / f"{accession}.txt").write_text(text, encoding="utf-8")
        return len(spectra)
    if suffix not in WRITE_EXTENSIONS:
        raise ValueError(
            f"unsupported output format {suffix!r}; supported: "
            f"{sorted(WRITE_EXTENSIONS)} or a directory for MassBank records")
    if p.exists():
        p.unlink()
    from matchms.exporting.save_spectra import save_spectra as _save
    if suffix in {".mgf", ".msp"}:
        _save(spectra, str(p), export_style=export_style)
    else:
        _save(spectra, str(p))
    return len(spectra)
