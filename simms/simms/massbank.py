"""Reading, sampling, writing and lightly validating MassBank record files.

The MassBank-data repository (https://github.com/MassBank/MassBank-data) is a
tree of contributor directories holding one ``MSBNK-*.txt`` record per
spectrum. This module turns those records into matchms Spectrum objects (for
use as templates in simulation), samples them efficiently without loading the
whole 139k-record collection, and can write simulated spectra back out in
MassBank record format.
"""

from __future__ import annotations

import datetime as _dt
import os
import random
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np
from matchms import Spectrum

RECORD_GLOB = "MSBNK-*.txt"

# record field -> matchms-style metadata key
_SIMPLE_FIELDS = {
    "ACCESSION": "accession",
    "RECORD_TITLE": "record_title",
    "CH$FORMULA": "formula",
    "CH$EXACT_MASS": "parent_mass",
    "CH$SMILES": "smiles",
    "CH$IUPAC": "inchi",
    "AC$INSTRUMENT": "instrument",
    "AC$INSTRUMENT_TYPE": "instrument_type",
}


@dataclass
class MassBankRecord:
    path: Optional[str] = None
    metadata: Dict[str, str] = field(default_factory=dict)
    peaks: List[Tuple[float, float, int]] = field(default_factory=list)  # (mz, intensity, rel_int)


def parse_record(text: str, path: Optional[str] = None) -> MassBankRecord:
    record = MassBankRecord(path=path)
    meta = record.metadata
    in_peaks = False
    for raw_line in text.splitlines():
        if raw_line.startswith("//"):
            break
        if in_peaks:
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    mz = float(parts[0])
                    intensity = float(parts[1])
                    rel = int(float(parts[2])) if len(parts) > 2 else 0
                except ValueError:
                    continue
                record.peaks.append((mz, intensity, rel))
            continue
        if ":" not in raw_line:
            continue
        key, _, value = raw_line.partition(":")
        key = key.strip()
        value = value.strip()
        if key == "PK$PEAK":
            in_peaks = True
            continue
        if key in _SIMPLE_FIELDS:
            meta[_SIMPLE_FIELDS[key]] = value
        elif key == "CH$NAME":
            meta.setdefault("compound_name", value)
        elif key == "CH$LINK":
            link_type, _, link_value = value.partition(" ")
            if link_type == "INCHIKEY":
                meta["inchikey"] = link_value.strip()
        elif key == "AC$MASS_SPECTROMETRY":
            sub_key, _, sub_value = value.partition(" ")
            sub_value = sub_value.strip()
            if sub_key == "MS_TYPE":
                meta["ms_type"] = sub_value
            elif sub_key == "ION_MODE":
                meta["ionmode"] = sub_value.lower()
            elif sub_key == "COLLISION_ENERGY":
                meta["collision_energy"] = sub_value
            elif sub_key == "FRAGMENTATION_MODE":
                meta["fragmentation_mode"] = sub_value
        elif key == "AC$CHROMATOGRAPHY":
            sub_key, _, sub_value = value.partition(" ")
            if sub_key == "RETENTION_TIME":
                meta["retention_time"] = sub_value.strip()
        elif key == "MS$FOCUSED_ION":
            sub_key, _, sub_value = value.partition(" ")
            sub_value = sub_value.strip()
            if sub_key == "PRECURSOR_M/Z":
                meta["precursor_mz"] = sub_value
            elif sub_key == "PRECURSOR_TYPE":
                meta["adduct"] = sub_value
        elif key == "PK$SPLASH":
            meta["splash"] = value
    return record


def load_record(path: str) -> MassBankRecord:
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        return parse_record(handle.read(), path=str(path))


def retention_time_seconds(value: Optional[str]) -> Optional[float]:
    """Parse a MassBank RETENTION_TIME value ('11.9 min', '714 sec', '11.9') to seconds."""
    if not value:
        return None
    match = re.match(r"([0-9.]+)\s*(min|sec|s)?", value.strip(), re.IGNORECASE)
    if not match:
        return None
    number = float(match.group(1))
    unit = (match.group(2) or "min").lower()
    return number * 60.0 if unit == "min" else number


def record_to_spectrum(record: MassBankRecord) -> Optional[Spectrum]:
    """Convert a parsed record into a matchms Spectrum (None if peak-less)."""
    if not record.peaks:
        return None
    mz = np.array([p[0] for p in record.peaks], dtype=float)
    intensities = np.array([p[1] for p in record.peaks], dtype=float)
    order = np.argsort(mz)
    metadata = dict(record.metadata)
    for numeric_key in ("precursor_mz", "parent_mass"):
        if numeric_key in metadata:
            try:
                metadata[numeric_key] = float(metadata[numeric_key])
            except ValueError:
                del metadata[numeric_key]
    rt = retention_time_seconds(metadata.get("retention_time"))
    if rt is not None:
        metadata["retention_time"] = rt
    else:
        metadata.pop("retention_time", None)
    # translate MS_TYPE into ms_level and drop the raw key: matchms
    # harmonization would otherwise map ms_type onto ms_level and clash.
    ms_type = metadata.pop("ms_type", "")
    if ms_type.upper().startswith("MS") and ms_type[2:].isdigit():
        metadata["ms_level"] = int(ms_type[2:])
    elif ms_type.upper() == "MS":
        metadata["ms_level"] = 1
    return Spectrum(mz=mz[order], intensities=intensities[order],
                    metadata=metadata, metadata_harmonization=False)


def iter_record_paths(repo: str, contributors: Optional[List[str]] = None) -> Iterator[Path]:
    root = Path(repo)
    if not root.is_dir():
        raise FileNotFoundError(f"MassBank data directory not found: {repo}")
    dirs = sorted(d for d in root.iterdir()
                  if d.is_dir() and not d.name.startswith("."))
    if contributors:
        wanted = {c.lower() for c in contributors}
        dirs = [d for d in dirs if d.name.lower() in wanted]
    for directory in dirs:
        yield from sorted(directory.glob(RECORD_GLOB))


def _matches(record: MassBankRecord, filters: Dict[str, str]) -> bool:
    for key, expected in filters.items():
        actual = record.metadata.get(key, "")
        if expected.lower() not in str(actual).lower():
            return False
    return True


def sample_records(repo: str, n: int, seed: int = 0,
                   filters: Optional[Dict[str, str]] = None,
                   contributors: Optional[List[str]] = None,
                   scan_limit: Optional[int] = None) -> List[MassBankRecord]:
    """Reservoir-sample n records matching the filters.

    Filters are substring matches on parsed metadata keys (e.g.
    ``ms_type=MS2``, ``ionmode=positive``, ``instrument_type=QTOF``).
    Filtering requires parsing each candidate, so with heavy filters over the
    full repository expect a full scan; ``scan_limit`` caps the number of
    files inspected for faster (biased) sampling.
    """
    rng = random.Random(seed)
    paths = list(iter_record_paths(repo, contributors))
    rng.shuffle(paths)
    if scan_limit:
        paths = paths[:scan_limit]
    chosen: List[MassBankRecord] = []
    if not filters:
        for path in paths[: n * 3]:  # small headroom for unparseable/peak-less records
            record = load_record(str(path))
            if record.peaks:
                chosen.append(record)
            if len(chosen) >= n:
                break
        return chosen[:n]
    for path in paths:
        record = load_record(str(path))
        if record.peaks and _matches(record, filters):
            chosen.append(record)
            if len(chosen) >= n:
                break
    return chosen


def repo_stats(repo: str, per_contributor_limit: Optional[int] = None) -> Dict[str, object]:
    """Counts of records per contributor (fast, filename-based)."""
    root = Path(repo)
    contributors = {}
    total = 0
    for directory in sorted(d for d in root.iterdir() if d.is_dir() and not d.name.startswith(".")):
        count = sum(1 for _ in directory.glob(RECORD_GLOB))
        if count:
            contributors[directory.name] = count
            total += count
        if per_contributor_limit and len(contributors) >= per_contributor_limit:
            break
    return {"path": str(root), "total_records": total, "contributors": contributors}


REQUIRED_FIELDS = ["ACCESSION", "RECORD_TITLE", "DATE", "AUTHORS", "LICENSE",
                   "CH$NAME", "AC$INSTRUMENT", "AC$INSTRUMENT_TYPE",
                   "AC$MASS_SPECTROMETRY: MS_TYPE", "AC$MASS_SPECTROMETRY: ION_MODE",
                   "PK$NUM_PEAK", "PK$PEAK"]


def validate_record_text(text: str) -> List[str]:
    """Light structural validation of a MassBank record. Returns issues found.

    This is not a replacement for the official Java Validator from
    MassBank-cli-tools; it catches the structural mistakes that break parsers.
    """
    issues = []
    for required in REQUIRED_FIELDS:
        if not re.search(rf"^{re.escape(required)}", text, re.MULTILINE):
            issues.append(f"missing required field {required}")
    if not text.rstrip().endswith("//"):
        issues.append("record must end with //")
    record = parse_record(text)
    declared = re.search(r"PK\$NUM_PEAK:\s*(\d+)", text)
    if declared and int(declared.group(1)) != len(record.peaks):
        issues.append(f"PK$NUM_PEAK is {declared.group(1)} but {len(record.peaks)} peak rows found")
    mzs = [p[0] for p in record.peaks]
    if mzs != sorted(mzs):
        issues.append("peaks are not sorted by m/z")
    return issues


def spectrum_to_record_text(spectrum: Spectrum, accession: str,
                            contributor: str = "SIMMS",
                            license_name: str = "CC BY") -> str:
    """Serialize a (simulated) spectrum as a MassBank-format record."""
    meta = spectrum.metadata
    mz = spectrum.peaks.mz
    intensities = spectrum.peaks.intensities
    rel = np.round(intensities * 999.0 / intensities.max()).astype(int) if len(intensities) else []
    name = meta.get("compound_name", "Simulated compound")
    instrument_type = meta.get("instrument_type", "LC-ESI-QTOF")
    level = meta.get("ms_level", 2)
    try:
        level = int(level)
    except (TypeError, ValueError):
        level = 2
    ms_type = meta.get("ms_type") or ("MS" if level <= 1 else f"MS{level}")
    ion_mode = str(meta.get("ionmode", "positive")).upper()
    lines = [
        f"ACCESSION: {accession}",
        f"RECORD_TITLE: {name}; {instrument_type}; {ms_type}; simulated",
        f"DATE: {_dt.date.today().strftime('%Y.%m.%d')}",
        f"AUTHORS: simms simulated data generator",
        f"LICENSE: {license_name}",
        f"COMMENT: SIMULATED record generated by simms"
        + (f" from {meta['accession']}" if meta.get("accession") else ""),
        f"CH$NAME: {name}",
        "CH$COMPOUND_CLASS: N/A; Simulated",
    ]
    if meta.get("formula"):
        lines.append(f"CH$FORMULA: {meta['formula']}")
    if meta.get("parent_mass"):
        lines.append(f"CH$EXACT_MASS: {meta['parent_mass']}")
    if meta.get("smiles"):
        lines.append(f"CH$SMILES: {meta['smiles']}")
    if meta.get("inchi"):
        lines.append(f"CH$IUPAC: {meta['inchi']}")
    if meta.get("inchikey"):
        lines.append(f"CH$LINK: INCHIKEY {meta['inchikey']}")
    lines += [
        f"AC$INSTRUMENT: simms in-silico instrument",
        f"AC$INSTRUMENT_TYPE: {instrument_type}",
        f"AC$MASS_SPECTROMETRY: MS_TYPE {ms_type}",
        f"AC$MASS_SPECTROMETRY: ION_MODE {ion_mode}",
    ]
    if meta.get("collision_energy"):
        lines.append(f"AC$MASS_SPECTROMETRY: COLLISION_ENERGY {meta['collision_energy']}")
    if meta.get("precursor_mz"):
        lines.append(f"MS$FOCUSED_ION: PRECURSOR_M/Z {meta['precursor_mz']}")
    if meta.get("adduct"):
        lines.append(f"MS$FOCUSED_ION: PRECURSOR_TYPE {meta['adduct']}")
    lines.append(f"PK$NUM_PEAK: {len(mz)}")
    lines.append("PK$PEAK: m/z int. rel.int.")
    for i in range(len(mz)):
        lines.append(f"  {mz[i]:.4f} {intensities[i]:.1f} {int(rel[i])}")
    lines.append("//")
    return "\n".join(lines) + "\n"
