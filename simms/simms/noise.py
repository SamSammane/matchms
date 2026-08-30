"""Configurable noise models for turning ideal spectra into realistic ones.

All randomness flows through a caller-supplied ``numpy.random.Generator`` so
every simulation is reproducible from a seed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np


@dataclass
class NoiseModel:
    """Parameters of the peak-level noise model.

    mz_ppm: gaussian sigma of mass error, in ppm.
    intensity_cv: coefficient of variation of multiplicative intensity noise.
    dropout: probability that any real peak is missing from the spectrum.
    noise_peaks: number of spurious peaks to add.
    noise_intensity: mean intensity of spurious peaks, as a fraction of the
        base peak (drawn from an exponential distribution).
    mz_range: m/z window in which spurious peaks appear; defaults to the
        span of the real peaks.
    """

    mz_ppm: float = 5.0
    intensity_cv: float = 0.15
    dropout: float = 0.05
    noise_peaks: int = 5
    noise_intensity: float = 0.01
    mz_range: Optional[Tuple[float, float]] = None

    def apply(self, mz: np.ndarray, intensities: np.ndarray,
              rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
        mz = np.asarray(mz, dtype=float)
        intensities = np.asarray(intensities, dtype=float)
        if mz.size == 0:
            return mz, intensities

        keep = rng.random(mz.size) >= self.dropout
        if not keep.any():  # never drop everything
            keep[int(np.argmax(intensities))] = True
        mz = mz[keep]
        intensities = intensities[keep]

        if self.mz_ppm > 0:
            mz = mz * (1.0 + rng.normal(0.0, self.mz_ppm * 1e-6, size=mz.size))
        if self.intensity_cv > 0:
            factors = rng.normal(1.0, self.intensity_cv, size=intensities.size)
            intensities = intensities * np.clip(factors, 0.05, None)

        if self.noise_peaks > 0:
            low, high = self.mz_range if self.mz_range else (float(mz.min()), float(mz.max()))
            if high <= low:
                high = low + 1.0
            base_peak = float(intensities.max())
            noise_mz = rng.uniform(low, high, size=self.noise_peaks)
            noise_int = rng.exponential(self.noise_intensity * base_peak, size=self.noise_peaks)
            mz = np.concatenate([mz, noise_mz])
            intensities = np.concatenate([intensities, noise_int])

        order = np.argsort(mz)
        return mz[order], intensities[order]


NOISE_PRESETS = {
    "none": NoiseModel(mz_ppm=0, intensity_cv=0, dropout=0, noise_peaks=0),
    "clean-orbitrap": NoiseModel(mz_ppm=2, intensity_cv=0.08, dropout=0.02, noise_peaks=3, noise_intensity=0.005),
    "default": NoiseModel(),
    "noisy-qtof": NoiseModel(mz_ppm=10, intensity_cv=0.25, dropout=0.10, noise_peaks=15, noise_intensity=0.02),
    "harsh": NoiseModel(mz_ppm=25, intensity_cv=0.5, dropout=0.25, noise_peaks=40, noise_intensity=0.05),
}
