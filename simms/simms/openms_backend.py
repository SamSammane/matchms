"""Optional OpenMS TOPP tool backend.

When an OpenMS installation is on PATH (or pointed to by OPENMS_BIN), simms
exposes its ~150 TOPP tools for native-format work that the Python stack
does not cover: mzML/featureXML/consensusXML merging (FileMerger), FASTA ->
theoretical mzML (IDFileConverter), format conversion (FileConverter),
feature linking, decoy generation, etc. Nothing here requires OpenMS at
import time; discovery is lazy and absence is reported, not fatal.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from typing import Dict, List, Optional, Sequence

KNOWN_TOOLS = [
    # generation / synthetic data
    "IDFileConverter", "DecoyDatabase", "Digestor", "DigestorMotif", "RNADigestor",
    "OpenSwathDecoyGenerator", "AssayGeneratorMetabo", "SeedListGenerator",
    "MassCalculator", "Resampler",
    # combination / linking
    "FileMerger", "IDMerger", "SpectraMerger", "QCMerger",
    "FeatureLinkerUnlabeled", "FeatureLinkerUnlabeledQT", "FeatureLinkerUnlabeledKD",
    "ConsensusID",
    # conversion / processing / validation
    "FileConverter", "FileInfo", "FileFilter", "PeakPickerHiRes",
    "FeatureFinderCentroided", "FuzzyDiff", "MzMLSplitter",
]


def _search_dirs() -> List[str]:
    dirs = []
    env_bin = os.environ.get("OPENMS_BIN")
    if env_bin:
        dirs.append(env_bin)
    dirs.extend((os.environ.get("PATH") or "").split(os.pathsep))
    return dirs


def find_tool(name: str) -> Optional[str]:
    env_bin = os.environ.get("OPENMS_BIN")
    if env_bin:
        candidate = os.path.join(env_bin, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return shutil.which(name)


def discover() -> Dict[str, str]:
    """Map of available TOPP tool name -> executable path."""
    found = {}
    for tool in KNOWN_TOOLS:
        path = find_tool(tool)
        if path:
            found[tool] = path
    return found


def run_tool(name: str, args: Sequence[str]) -> int:
    """Run a TOPP tool, passing args through verbatim. Returns its exit code."""
    path = find_tool(name)
    if path is None:
        raise FileNotFoundError(
            f"OpenMS tool {name!r} not found on PATH or in $OPENMS_BIN. "
            "Install OpenMS or set OPENMS_BIN to its bin directory.")
    completed = subprocess.run([path, *args])
    return completed.returncode


def merge_mzml(inputs: Sequence[str], output: str) -> int:
    """Merge mzML (or featureXML/consensusXML/TraML) files via FileMerger."""
    return run_tool("FileMerger", ["-in", *inputs, "-out", output])


def fasta_to_mzml(fasta: str, output: str, extra_args: Sequence[str] = ()) -> int:
    """Generate theoretical spectra from a FASTA database via IDFileConverter."""
    return run_tool("IDFileConverter", ["-in", fasta, "-out", output, *extra_args])
