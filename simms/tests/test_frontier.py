"""Tests for profile mode, Prosit/Koina ML intensities, and DIA acquisition."""

import json
import math

import numpy as np
import pytest

from simms import mlfrag
from simms.mlfrag import (FragmentPrediction, KoinaUnavailable, _build_request,
                          _parse_response, predict)
from simms.realism import (REALISM_PRESETS, RealismConfig, centroids_to_profile,
                           peak_sigma_mz)
from simms.simulate import simulate_lcms_run, theoretical_peptide_spectrum


# ------------------------------------------------------------ profile mode

def test_peak_sigma_scales_with_mz():
    sigma_200 = peak_sigma_mz(np.array([200.0]), 30000)[0]
    sigma_800 = peak_sigma_mz(np.array([800.0]), 30000)[0]
    assert sigma_200 == pytest.approx(200.0 / 30000 / 2.3548, rel=1e-4)
    # Orbitrap-like: resolution halves at 4x m/z, so sigma grows 8x
    assert sigma_800 / sigma_200 == pytest.approx(8.0, rel=1e-3)


def test_centroids_to_profile_shape_and_area():
    mz = np.array([500.0])
    intensities = np.array([1000.0])
    grid, profile = centroids_to_profile(mz, intensities, 30000)
    assert profile.max() == pytest.approx(1000.0, rel=1e-6)
    assert grid[np.argmax(profile)] == pytest.approx(500.0, abs=1e-4)
    # measured FWHM equals mz / R(mz)
    half = grid[profile >= 500.0]
    fwhm = half.max() - half.min()
    expected = 500.0 / (30000 * math.sqrt(200.0 / 500.0))
    assert fwhm == pytest.approx(expected, rel=0.15)


def test_centroids_to_profile_merges_unresolved_doublet():
    # two peaks closer than the peak width fuse into one profile maximum
    mz = np.array([500.000, 500.002])
    intensities = np.array([800.0, 900.0])
    grid, profile = centroids_to_profile(mz, intensities, 20000)
    maxima = np.sum((profile[1:-1] > profile[:-2]) & (profile[1:-1] > profile[2:])
                    & (profile[1:-1] > 100))
    assert maxima == 1
    assert profile.max() > 900.0  # overlap sums


def test_centroids_to_profile_empty():
    grid, profile = centroids_to_profile(np.array([]), np.array([]), 30000)
    assert grid.size == 0 and profile.size == 0


def test_lcms_profile_mode_roundtrip(tmp_path):
    from pyteomics import mzml as pymzml
    spectra = [theoretical_peptide_spectrum("ELVISLIVESK", 2, model="simple")]
    config = RealismConfig(profile_mode=True, resolving_power=60000,
                           chemical_noise_peaks=0, contaminants=False,
                           charge_envelope=False, drift_ppm=0,
                           spray_instability_cv=0, tailing_tau_factor=0)
    out = tmp_path / "profile.mzML"
    result = simulate_lcms_run(spectra, str(out), gradient_seconds=30,
                               ms1_interval_seconds=5, realism=config, seed=1)
    assert result["profile_mode"] is True
    scans = list(pymzml.read(str(out)))
    busy = max(scans, key=lambda s: len(s["m/z array"]))
    assert "profile spectrum" in busy
    assert "centroid spectrum" not in busy
    # profile data has many more points than the centroid count
    assert len(busy["m/z array"]) > 50


# ------------------------------------------------------------- Koina/Prosit

def _fake_koina_payload(n, width=6):
    # two valid ions per peptide; the rest masked with -1 (impossible ions)
    intensities, mz, annotation = [], [], []
    for i in range(n):
        intensities += [0.2, 1.0] + [-1.0] * (width - 2)
        mz += [175.119, 401.256] + [-1.0] * (width - 2)
        annotation += ["y1+1", "y3+1"] + [""] * (width - 2)
    return {"outputs": [
        {"name": "intensities", "shape": [n, width], "data": intensities},
        {"name": "mz", "shape": [n, width], "data": mz},
        {"name": "annotation", "shape": [n, width], "data": annotation},
    ]}


def test_build_request_wire_format():
    body = _build_request(["PEPTIDEK"], [2], [27.0])
    names = {i["name"]: i for i in body["inputs"]}
    assert names["peptide_sequences"]["datatype"] == "BYTES"
    assert names["peptide_sequences"]["shape"] == [1, 1]
    assert names["precursor_charges"]["data"] == [2]
    assert names["collision_energies"]["data"] == [27.0]


def test_parse_response_masks_invalid_ions():
    predictions = _parse_response(_fake_koina_payload(2), ["AAK", "AAR"], [2, 2],
                                  [25.0, 25.0])
    assert len(predictions) == 2
    assert predictions[0].mz == [175.119, 401.256]
    assert predictions[0].intensities == [0.2, 1.0]
    assert predictions[0].annotations == ["y1+1", "y3+1"]


def test_predict_uses_cache(tmp_path, monkeypatch):
    monkeypatch.setattr(mlfrag, "CACHE_DIR", tmp_path / "koina")
    calls = []

    def fake_post(url, body, timeout):
        calls.append(url)
        n = body["inputs"][0]["shape"][0]
        return _fake_koina_payload(n)

    first = predict(["PEPTIDEK"], [2], [25.0], _post_fn=fake_post)
    assert len(calls) == 1
    assert first[0].peptide == "PEPTIDEK"
    # second call is served from cache: no new network hit
    second = predict(["PEPTIDEK"], [2], [25.0],
                     _post_fn=lambda *a: (_ for _ in ()).throw(AssertionError("no")))
    assert second[0].mz == first[0].mz
    assert len(calls) == 1


def test_predict_raises_koina_unavailable(tmp_path, monkeypatch):
    monkeypatch.setattr(mlfrag, "CACHE_DIR", tmp_path / "koina")

    def down(url, body, timeout):
        raise OSError("connection refused")

    with pytest.raises(KoinaUnavailable, match="realistic"):
        predict(["PEPTIDEK"], [2], [25.0], _post_fn=down)


def test_predict_rejects_overlong_sequence():
    with pytest.raises(ValueError, match="30"):
        predict(["A" * 31], [2], [25.0], _post_fn=lambda *a: None)


def test_prosit_spectrum_via_mocked_transport(tmp_path, monkeypatch):
    monkeypatch.setattr(mlfrag, "CACHE_DIR", tmp_path / "koina")
    monkeypatch.setattr(mlfrag, "_post",
                        lambda url, body, timeout: _fake_koina_payload(
                            body["inputs"][0]["shape"][0]))
    spectrum = theoretical_peptide_spectrum("PEPTIDEK", 2, model="prosit",
                                            collision_energy=27)
    assert spectrum.metadata["fragment_model"] == "prosit"
    assert list(spectrum.peaks.mz) == [175.119, 401.256]
    assert "y3+1" in spectrum.metadata["ion_annotations"]


# --------------------------------------------------------------------- DIA

def _dia_run(tmp_path, **kwargs):
    spectra = [theoretical_peptide_spectrum("ELVISLIVESK", 2, model="simple"),
               theoretical_peptide_spectrum("SETPEPTIDER", 2, model="simple")]
    out = tmp_path / "dia.mzML"
    defaults = dict(gradient_seconds=30, ms1_interval_seconds=5,
                    acquisition="dia", dia_range=(400.0, 700.0),
                    dia_window=25.0, dia_overlap=1.0,
                    realism=REALISM_PRESETS["none"], seed=2)
    defaults.update(kwargs)
    return out, simulate_lcms_run(spectra, str(out), **defaults)


def test_dia_full_cycle_every_ms1(tmp_path):
    _, result = _dia_run(tmp_path)
    assert result["acquisition"] == "dia"
    assert result["dia_windows_per_cycle"] == 12  # (700-400)/25
    assert result["ms2_scans"] == result["ms1_scans"] * 12


def test_dia_isolation_windows_in_mzml(tmp_path):
    from pyteomics import mzml as pymzml
    path, _ = _dia_run(tmp_path)
    ms2 = [s for s in pymzml.read(str(path)) if s["ms level"] == 2]
    window = ms2[0]["precursorList"]["precursor"][0]["isolationWindow"]
    assert window["isolation window target m/z"] == pytest.approx(412.5)
    assert window["isolation window lower offset"] == pytest.approx(13.0)
    assert window["isolation window upper offset"] == pytest.approx(13.0)


def test_dia_multiplexing_counts(tmp_path):
    # a wide single window over both precursors multiplexes them
    _, result = _dia_run(tmp_path, dia_range=(400.0, 1000.0), dia_window=600.0)
    assert result["dia_windows_per_cycle"] == 1
    assert result["multiplexed_ms2_scans"] > 0


def test_dia_rejects_unknown_scheme(tmp_path):
    with pytest.raises(ValueError, match="acquisition"):
        _dia_run(tmp_path, acquisition="swath")


def test_dia_cli_end_to_end(tmp_path, capsys):
    from simms.cli import main
    mgf = tmp_path / "in.mgf"
    code = main(["generate", "peptides", "--sequences", "ELVISLIVESK",
                 "--fragment-model", "simple", "--out", str(mgf)])
    assert code == 0
    capsys.readouterr()
    out = tmp_path / "run.mzML"
    code = main(["generate", "lcms-run", "-i", str(mgf), "--gradient", "30",
                 "--ms1-interval", "5", "--acquisition", "dia",
                 "--dia-range", "400", "700", "--realism", "none",
                 "--out", str(out), "--json"])
    assert code == 0
    result = json.loads(capsys.readouterr().out)
    assert result["acquisition"] == "dia"
    assert out.exists()
