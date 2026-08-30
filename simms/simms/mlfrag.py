"""ML-predicted peptide fragment intensities via Koina (Prosit models).

Koina (https://koina.wilhelmlab.org) serves the Prosit deep-learning models
behind a Triton inference HTTP API. This module implements the wire format
for the intensity models (``Prosit_2020_intensity_HCD`` by default), with a
local response cache so repeated simulations do not re-query the service.

The network is optional infrastructure: when Koina is unreachable (offline
machines, restricted network policies) prediction raises
:class:`KoinaUnavailable` with actionable advice, and the caller can fall
back to the built-in mobile-proton model (``--fragment-model realistic``).
"""

from __future__ import annotations

import hashlib
import json
import os
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

DEFAULT_MODEL = "Prosit_2020_intensity_HCD"
DEFAULT_URL = "https://koina.wilhelmlab.org"
CACHE_DIR = Path(os.environ.get("SIMMS_CACHE",
                                Path.home() / ".cache" / "simms")) / "koina"

# Prosit supports plain sequences up to 30 residues, charge 1-6
MAX_SEQUENCE_LENGTH = 30


class KoinaUnavailable(RuntimeError):
    """Raised when the Koina inference service cannot be reached."""


@dataclass
class FragmentPrediction:
    peptide: str
    charge: int
    collision_energy: float
    mz: List[float]
    intensities: List[float]
    annotations: List[str]


def _cache_key(model: str, peptide: str, charge: int, ce: float) -> Path:
    digest = hashlib.sha256(f"{model}|{peptide}|{charge}|{ce:.2f}".encode()).hexdigest()
    return CACHE_DIR / model / f"{digest}.json"


def _build_request(peptides: Sequence[str], charges: Sequence[int],
                   collision_energies: Sequence[float]) -> Dict:
    n = len(peptides)
    return {
        "id": "simms",
        "inputs": [
            {"name": "peptide_sequences", "shape": [n, 1], "datatype": "BYTES",
             "data": list(peptides)},
            {"name": "precursor_charges", "shape": [n, 1], "datatype": "INT32",
             "data": [int(c) for c in charges]},
            {"name": "collision_energies", "shape": [n, 1], "datatype": "FP32",
             "data": [float(ce) for ce in collision_energies]},
        ],
    }


def _parse_response(payload: Dict, peptides: Sequence[str], charges: Sequence[int],
                    collision_energies: Sequence[float]) -> List[FragmentPrediction]:
    outputs = {out["name"]: out for out in payload.get("outputs", [])}
    for required in ("intensities", "mz", "annotation"):
        if required not in outputs:
            raise ValueError(f"Koina response missing output {required!r}")
    n = len(peptides)
    width = outputs["intensities"]["shape"][1]
    intensities = outputs["intensities"]["data"]
    mz = outputs["mz"]["data"]
    annotations = outputs["annotation"]["data"]
    predictions = []
    for i in range(n):
        row = slice(i * width, (i + 1) * width)
        row_int = intensities[row]
        row_mz = mz[row]
        row_ann = annotations[row]
        keep_mz, keep_int, keep_ann = [], [], []
        for m, inten, ann in zip(row_mz, row_int, row_ann):
            # negative intensity marks an impossible ion for this peptide
            if inten is None or inten <= 0 or m is None or m <= 0:
                continue
            keep_mz.append(float(m))
            keep_int.append(float(inten))
            if isinstance(ann, bytes):
                ann = ann.decode("utf-8", "replace")
            keep_ann.append(str(ann))
        predictions.append(FragmentPrediction(
            peptide=peptides[i], charge=int(charges[i]),
            collision_energy=float(collision_energies[i]),
            mz=keep_mz, intensities=keep_int, annotations=keep_ann))
    return predictions


def _post(url: str, body: Dict, timeout: float) -> Dict:
    request = urllib.request.Request(
        url, data=json.dumps(body).encode("utf-8"),
        headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def predict(peptides: Sequence[str], charges: Sequence[int],
            collision_energies: Sequence[float],
            model: str = DEFAULT_MODEL,
            url: Optional[str] = None,
            timeout: float = 30.0,
            use_cache: bool = True,
            _post_fn=None) -> List[FragmentPrediction]:
    """Predict fragment intensity spectra for (peptide, charge, CE) triples.

    Results are cached under ``~/.cache/simms/koina`` (override with the
    ``SIMMS_CACHE`` environment variable). ``_post_fn`` exists for tests.
    """
    if not (len(peptides) == len(charges) == len(collision_energies)):
        raise ValueError("peptides, charges and collision_energies must align")
    for peptide in peptides:
        if len(peptide) > MAX_SEQUENCE_LENGTH:
            raise ValueError(
                f"{peptide!r} exceeds Prosit's {MAX_SEQUENCE_LENGTH}-residue limit")

    url = (url or os.environ.get("SIMMS_KOINA_URL") or DEFAULT_URL).rstrip("/")
    poster = _post_fn or _post

    results: List[Optional[FragmentPrediction]] = [None] * len(peptides)
    missing: List[int] = []
    if use_cache:
        for i, (peptide, charge, ce) in enumerate(zip(peptides, charges, collision_energies)):
            cache_path = _cache_key(model, peptide, int(charge), float(ce))
            if cache_path.is_file():
                data = json.loads(cache_path.read_text())
                results[i] = FragmentPrediction(**data)
            else:
                missing.append(i)
    else:
        missing = list(range(len(peptides)))

    if missing:
        batch_peptides = [peptides[i] for i in missing]
        batch_charges = [charges[i] for i in missing]
        batch_ces = [collision_energies[i] for i in missing]
        endpoint = f"{url}/v2/models/{model}/infer"
        try:
            payload = poster(endpoint, _build_request(batch_peptides, batch_charges,
                                                      batch_ces), timeout)
        except (urllib.error.URLError, OSError, TimeoutError) as err:
            raise KoinaUnavailable(
                f"Koina inference service unreachable at {endpoint} ({err}). "
                "Check network access (SIMMS_KOINA_URL overrides the host) or "
                "use --fragment-model realistic for the built-in offline model."
            ) from err
        predicted = _parse_response(payload, batch_peptides, batch_charges, batch_ces)
        for index, prediction in zip(missing, predicted):
            results[index] = prediction
            if use_cache:
                cache_path = _cache_key(model, prediction.peptide,
                                        prediction.charge, prediction.collision_energy)
                cache_path.parent.mkdir(parents=True, exist_ok=True)
                cache_path.write_text(json.dumps(prediction.__dict__))
    return [r for r in results if r is not None]
