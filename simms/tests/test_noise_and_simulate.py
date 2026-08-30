import numpy as np
import pytest

from simms.noise import NOISE_PRESETS, NoiseModel
from simms.simulate import (isotope_spectrum, simulate_variants,
                            theoretical_peptide_spectrum)


def test_noise_none_is_identity():
    rng = np.random.default_rng(0)
    mz = np.array([100.0, 200.0, 300.0])
    intensities = np.array([1.0, 2.0, 3.0])
    out_mz, out_int = NOISE_PRESETS["none"].apply(mz, intensities, rng)
    np.testing.assert_array_equal(out_mz, mz)
    np.testing.assert_array_equal(out_int, intensities)


def test_noise_is_seed_reproducible():
    model = NoiseModel(mz_ppm=10, intensity_cv=0.2, dropout=0.1, noise_peaks=5)
    mz = np.linspace(100, 1000, 50)
    intensities = np.ones(50)
    a = model.apply(mz, intensities, np.random.default_rng(123))
    b = model.apply(mz, intensities, np.random.default_rng(123))
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])


def test_noise_never_drops_all_peaks():
    model = NoiseModel(dropout=1.0, noise_peaks=0, mz_ppm=0, intensity_cv=0)
    mz = np.array([100.0, 200.0])
    intensities = np.array([1.0, 5.0])
    out_mz, out_int = model.apply(mz, intensities, np.random.default_rng(0))
    assert len(out_mz) == 1
    assert out_mz[0] == 200.0  # the base peak survives


def test_noise_adds_noise_peaks():
    model = NoiseModel(mz_ppm=0, intensity_cv=0, dropout=0, noise_peaks=10)
    mz = np.array([100.0, 500.0])
    intensities = np.array([1.0, 2.0])
    out_mz, _ = model.apply(mz, intensities, np.random.default_rng(0))
    assert len(out_mz) == 12
    assert np.all(np.diff(out_mz) >= 0)  # sorted


def test_theoretical_peptide_spectrum():
    spectrum = theoretical_peptide_spectrum("PEPTIDE", charge=2)
    assert spectrum.metadata["precursor_mz"] == pytest.approx(400.6873, abs=1e-3)
    assert spectrum.metadata["ms_level"] == 2
    assert len(spectrum.peaks) == 12  # b1-b6 + y1-y6


def test_isotope_spectrum():
    spectrum = isotope_spectrum("C9H8O4", compound_name="aspirin")
    assert spectrum.metadata["compound_name"] == "aspirin"
    assert spectrum.peaks.mz[0] == pytest.approx(181.0495, abs=1e-3)
    assert spectrum.peaks.intensities.max() == pytest.approx(100.0)


def test_simulate_variants_reproducible_and_annotated():
    template = theoretical_peptide_spectrum("ELVISLIVESK")
    noise = NOISE_PRESETS["default"]
    run1 = simulate_variants([template], noise, n_variants=3, seed=99)
    run2 = simulate_variants([template], noise, n_variants=3, seed=99)
    assert len(run1) == 3
    for a, b in zip(run1, run2):
        np.testing.assert_array_equal(a.peaks.mz, b.peaks.mz)
    assert run1[0].metadata["variant"] == 0
    assert run1[2].metadata["variant"] == 2
    different = simulate_variants([template], noise, n_variants=1, seed=100)
    assert not np.array_equal(run1[0].peaks.mz, different[0].peaks.mz)
