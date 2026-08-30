"""Combining spectral libraries across files and formats."""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from matchms import Spectrum

from .io_utils import load_many, save_any


def _dedupe(spectra: Sequence[Spectrum], key: str) -> List[Spectrum]:
    """Keep, per metadata key value, the spectrum with the most peaks."""
    best: Dict[str, Spectrum] = {}
    keyless: List[Spectrum] = []
    for spectrum in spectra:
        value = spectrum.metadata.get(key)
        if not value:
            keyless.append(spectrum)
            continue
        current = best.get(value)
        if current is None or len(spectrum.peaks) > len(current.peaks):
            best[value] = spectrum
    return list(best.values()) + keyless


def merge_files(inputs: Sequence[str], output: str,
                dedupe_key: Optional[str] = None,
                min_peaks: int = 0,
                ms_level: Optional[int] = None,
                ionmode: Optional[str] = None,
                export_style: str = "matchms") -> Dict[str, int]:
    """Merge spectra from any mix of supported input files into one output."""
    spectra = load_many(inputs)
    loaded = len(spectra)
    if min_peaks:
        spectra = [s for s in spectra if len(s.peaks) >= min_peaks]
    if ms_level is not None:
        spectra = [s for s in spectra if s.metadata.get("ms_level") == ms_level]
    if ionmode:
        spectra = [s for s in spectra if str(s.metadata.get("ionmode", "")).lower() == ionmode.lower()]
    if dedupe_key:
        spectra = _dedupe(spectra, dedupe_key)
    written = save_any(spectra, output, export_style=export_style)
    return {"inputs": len(inputs), "loaded": loaded, "written": written}
