"""Tests for the instrument/chromatography realism layer."""

import numpy as np
import pytest

from simms.chem import averagine_envelope_mz, averagine_formula
from simms.peptides import proton_mobility, realistic_fragment_ions
from simms.realism import (REALISM_PRESETS, RealismConfig, SprayStability,
                           apply_saturation, calibration_drift_ppm,
                           emg_profile)
from simms.simulate import simulate_lcms_run, theoretical_peptide_spectrum


# ------------------------------------------------------------- averagine

def test_averagine_formula_scales_with_mass():
    small = averagine_formula(500.0)
    large = averagine_formula(5000.0)
    assert large["C"] > small["C"] * 8
    assert set(small) <= {"C", "H", "N", "O", "S"}


def test_averagine_envelope_spacing_and_growth():
    env1 = averagine_envelope_mz(1000.0, 1)
    assert env1[1][0] - env1[0][0] == pytest.approx(1.00336, abs=0.001)
    env2 = averagine_envelope_mz(1000.0, 2)
    assert env2[1][0] - env2[0][0] == pytest.approx(0.50168, abs=0.001)
    # M+1/M ratio grows with mass (more carbons)
    ratio_1k = averagine_envelope_mz(1000.0, 1)[1][1]
    ratio_3k = averagine_envelope_mz(3000.0, 1)[1][1]
    assert ratio_3k > ratio_1k > 0.3


# -------------------------------------------------- realistic fragments

def test_proton_mobility_classes():
    assert proton_mobility("ELVISLIVESK", 3) == "mobile"
    assert proton_mobility("DGDGRDRAAK", 1) == "nonmobile"


def test_realistic_fragments_b1_suppressed_b2_present():
    ions = {i.ion_type: i.intensity for i in realistic_fragment_ions("ELVISLIVESK", 2)}
    assert ions.get("b1+", 0) < 1.0
    assert ions["b2+"] > 20.0
    assert ions["y6+"] > ions["b6+"]  # y dominates b


def test_realistic_fragments_proline_effect():
    # cleavage N-terminal to P boosts y and suppresses b at that position;
    # spectra are peak-normalized, so compare the y/b ratio at the site
    with_p = {i.ion_type: i.intensity for i in realistic_fragment_ions("AAAPAAAK", 2)}
    without = {i.ion_type: i.intensity for i in realistic_fragment_ions("AAALAAAK", 2)}
    ratio_p = with_p["y5+"] / with_p["b3+"]
    ratio_l = without["y5+"] / without["b3+"]
    assert ratio_p > 5.0 * ratio_l


def test_realistic_fragments_neutral_losses_and_immonium():
    labels = [i.ion_type for i in realistic_fragment_ions("SETPEPTIDER", 2, 30.0)]
    assert any(l.endswith("-H2O") for l in labels)
    assert any(l.endswith("-NH3") for l in labels)
    assert any(l.startswith("imm(") for l in labels)
    assert any(l.endswith("+i") for l in labels)  # fragment isotopes


def test_collision_energy_shifts_to_low_mass():
    def weighted_mean_mz(ce):
        ions = realistic_fragment_ions("ELVISLIVESK", 2, ce)
        mz = np.array([i.mz for i in ions])
        w = np.array([i.intensity for i in ions])
        return float((mz * w).sum() / w.sum())
    assert weighted_mean_mz(40) < weighted_mean_mz(15) - 30


def test_theoretical_spectrum_models():
    simple = theoretical_peptide_spectrum("PEPTIDE", model="simple")
    realistic = theoretical_peptide_spectrum("PEPTIDE", model="realistic")
    assert len(simple.peaks) == 12
    assert len(realistic.peaks) > len(simple.peaks)
    with pytest.raises(ValueError):
        theoretical_peptide_spectrum("PEPTIDE", model="nope")


# ------------------------------------------------------- realism pieces

def test_emg_degenerates_to_gaussian():
    t = np.linspace(0, 100, 500)
    gauss = emg_profile(t, 50, 5, 0.0)
    assert t[np.argmax(gauss)] == pytest.approx(50, abs=0.5)
    # symmetric: equal mass on both sides of apex
    left, right = gauss[t < 50].sum(), gauss[t > 50].sum()
    assert left == pytest.approx(right, rel=0.02)


def test_emg_tails_to_the_right():
    t = np.linspace(0, 200, 2000)
    profile = emg_profile(t, 50, 5, 10.0)
    apex_t = t[np.argmax(profile)]
    mean_t = float((t * profile).sum() / profile.sum())
    assert mean_t > apex_t + 2.0  # right-skewed
    assert profile.max() == pytest.approx(1.0)


def test_spray_stability_reproducible_and_correlated():
    a = SprayStability(0.1, np.random.default_rng(1))
    b = SprayStability(0.1, np.random.default_rng(1))
    seq_a = [a.next_factor() for _ in range(50)]
    seq_b = [b.next_factor() for _ in range(50)]
    assert seq_a == seq_b
    log = np.log(seq_a)
    lag1 = np.corrcoef(log[:-1], log[1:])[0, 1]
    assert lag1 > 0.5  # AR(1) memory


def test_calibration_drift_bounds():
    for t in np.linspace(0, 600, 50):
        assert abs(calibration_drift_ppm(t, 600, 6.0)) <= 3.0 + 1e-9
    assert calibration_drift_ppm(100, 600, 0.0) == 0.0


def test_saturation_clips_and_compresses():
    intensities = np.array([1e3, 4e8, 9.9e8, 5e9, 1e12])
    out = apply_saturation(intensities, 1e9)
    assert out[0] == 1e3                      # untouched below the knee
    assert out.max() == 1e9                   # hard clip at full scale
    assert out[2] < intensities[2]            # compressed above the knee
    np.testing.assert_array_equal(apply_saturation(intensities, None), intensities)


# --------------------------------------------------------- full LC-MS run

def _run(tmp_path, realism, seed=3, **kwargs):
    spectra = [theoretical_peptide_spectrum("ELVISLIVESK", 2, model="realistic"),
               theoretical_peptide_spectrum("SETPEPTIDER", 2, model="realistic")]
    out = tmp_path / "run.mzML"
    result = simulate_lcms_run(spectra, str(out), gradient_seconds=60,
                               peak_fwhm_seconds=8, ms1_interval_seconds=2,
                               realism=realism, seed=seed, **kwargs)
    return out, result


def test_lcms_charge_envelope_expands_species(tmp_path):
    _, with_env = _run(tmp_path, REALISM_PRESETS["default"])
    _, without = _run(tmp_path, REALISM_PRESETS["none"])
    assert without["precursor_species"] == 2
    assert with_env["precursor_species"] > 4  # 2-3 charge states per peptide


def test_lcms_contaminants_in_every_ms1(tmp_path):
    from pyteomics import mzml as pymzml
    path, _ = _run(tmp_path, REALISM_PRESETS["default"])
    ms1 = [s for s in pymzml.read(str(path)) if s["ms level"] == 1]
    hits = sum(1 for s in ms1 if np.any(np.abs(s["m/z array"] - 371.10124) < 0.05))
    assert hits == len(ms1)


def test_lcms_none_preset_is_sparse(tmp_path):
    from pyteomics import mzml as pymzml
    path, result = _run(tmp_path, REALISM_PRESETS["none"])
    assert result["chimeric_ms2_scans"] == 0
    sparse_min = min(len(s["m/z array"]) for s in pymzml.read(str(path))
                     if s["ms level"] == 1)
    path2, _ = _run(tmp_path, REALISM_PRESETS["default"])
    busy_min = min(len(s["m/z array"]) for s in pymzml.read(str(path2))
                   if s["ms level"] == 1)
    # default preset carries a persistent background (contaminants + chemical
    # noise) in every scan; the none preset has only real signal
    assert busy_min >= sparse_min + 20


def test_lcms_realism_run_is_deterministic(tmp_path):
    _, r1 = _run(tmp_path, REALISM_PRESETS["high"], seed=11)
    _, r2 = _run(tmp_path, REALISM_PRESETS["high"], seed=11)
    assert r1 == {**r2, "output": r1["output"]}


def test_lcms_dynamic_exclusion_limits_retriggers(tmp_path):
    tight = RealismConfig(dynamic_exclusion_seconds=1.0, charge_envelope=False,
                          chimeras=False, contaminants=False,
                          chemical_noise_peaks=0, spray_instability_cv=0,
                          drift_ppm=0, tailing_tau_factor=0)
    loose = RealismConfig(dynamic_exclusion_seconds=30.0, charge_envelope=False,
                          chimeras=False, contaminants=False,
                          chemical_noise_peaks=0, spray_instability_cv=0,
                          drift_ppm=0, tailing_tau_factor=0)
    _, many = _run(tmp_path, tight)
    _, few = _run(tmp_path, loose)
    assert many["ms2_scans"] > few["ms2_scans"]
