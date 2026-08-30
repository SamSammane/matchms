"""Instrument- and chromatography-realism models for LC-MS run simulation.

Everything here is deterministic given a numpy Generator, and configured
through :class:`RealismConfig` (presets: none / default / high).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.special import erfc


@dataclass
class RealismConfig:
    """Knobs for physically-motivated LC-MS realism.

    tailing_tau_factor: EMG exponential tail constant as a multiple of the
        gaussian sigma (0 = symmetric gaussian; ~1 = typical tailing).
    rt_broadening: fractional peak-width growth across the gradient
        (late-eluting peaks are wider).
    charge_envelope: distribute peptide precursors over neighboring charge
        states (z-1, z, z+1) as real electrospray does.
    isotope_envelopes: model MS1 isotope envelopes from the molecular
        formula, or the averagine approximation when only a mass is known.
    contaminants: include ever-present background ions (polysiloxanes,
        phthalates) in every MS1 scan.
    chemical_noise_peaks: count of low-intensity chemical-noise peaks per
        MS1 scan (0 disables).
    chemical_noise_level: their mean intensity, as a fraction of the most
        abundant compound apex.
    spray_instability_cv: per-scan correlated intensity fluctuation
        (AR(1) electrospray flicker) applied to all species together.
    drift_ppm: peak-to-peak amplitude of slow mass-calibration drift over
        the run (sinusoidal), on top of per-peak noise.
    saturation: detector full-scale; intensities are clipped there
        (None disables).
    isolation_window: MS2 isolation window full width in m/z; co-eluting
        precursors inside it produce chimeric MS2 spectra.
    chimeras: whether to mix co-isolated fragments into MS2 scans.
    dynamic_exclusion_seconds: re-trigger lockout per precursor.
    """

    tailing_tau_factor: float = 1.0
    rt_broadening: float = 0.5
    charge_envelope: bool = True
    isotope_envelopes: bool = True
    contaminants: bool = True
    chemical_noise_peaks: int = 30
    chemical_noise_level: float = 2e-4
    spray_instability_cv: float = 0.08
    drift_ppm: float = 3.0
    saturation: Optional[float] = 5e9
    isolation_window: float = 1.6
    chimeras: bool = True
    dynamic_exclusion_seconds: float = 20.0


REALISM_PRESETS: Dict[str, RealismConfig] = {
    "none": RealismConfig(tailing_tau_factor=0.0, rt_broadening=0.0,
                          charge_envelope=False, isotope_envelopes=False,
                          contaminants=False, chemical_noise_peaks=0,
                          spray_instability_cv=0.0, drift_ppm=0.0,
                          saturation=None, chimeras=False,
                          dynamic_exclusion_seconds=20.0),
    "default": RealismConfig(),
    "high": RealismConfig(tailing_tau_factor=1.8, rt_broadening=0.8,
                          chemical_noise_peaks=120, chemical_noise_level=5e-4,
                          spray_instability_cv=0.15, drift_ppm=8.0,
                          saturation=1e9, isolation_window=2.0,
                          dynamic_exclusion_seconds=12.0),
}


# ubiquitous ESI background ions: polysiloxanes and phthalates
CONTAMINANT_IONS: List[Tuple[float, float]] = [
    (371.10124, 1.0),    # polysiloxane [Si(CH3)2O]5 + H
    (445.12003, 0.6),    # polysiloxane [Si(CH3)2O]6 + H
    (536.16537, 0.25),   # polysiloxane cluster
    (149.02332, 0.8),    # phthalate fragment
    (279.15909, 0.5),    # dibutyl phthalate + H
    (391.28429, 0.4),    # diisooctyl phthalate + H
    (413.26623, 0.3),    # diisooctyl phthalate + Na
    (122.09643, 0.35),   # triethylamine + H
    (102.12773, 0.3),    # triethylamine-related
]


def emg_profile(t: np.ndarray, rt: float, sigma: float, tau: float) -> np.ndarray:
    """Exponentially modified gaussian elution profile, peak-normalized.

    tau -> 0 degenerates to a pure gaussian.
    """
    t = np.asarray(t, dtype=float)
    if tau < 1e-6:
        profile = np.exp(-0.5 * ((t - rt) / sigma) ** 2)
    else:
        arg = (sigma / tau - (t - rt) / sigma) / np.sqrt(2.0)
        # log-space to avoid overflow of exp() before erfc() cancels it
        log_amp = (sigma * sigma) / (2 * tau * tau) - (t - rt) / tau
        profile = np.where(arg < 25,
                           np.exp(np.clip(log_amp, -700, 700)) * erfc(arg),
                           0.0)
        # for very negative arg (far into the tail) erfc ~ 2, exp decays: fine
    peak = profile.max()
    if peak <= 0:
        return np.zeros_like(t)
    return profile / peak


def calibration_drift_ppm(t: float, gradient_seconds: float, drift_ppm: float,
                          phase: float = 0.0) -> float:
    """Slow sinusoidal mass-calibration drift (ppm) at time t."""
    if drift_ppm <= 0:
        return 0.0
    return 0.5 * drift_ppm * float(np.sin(2 * np.pi * t / gradient_seconds + phase))


class SprayStability:
    """AR(1) log-intensity flicker shared by all species in a scan."""

    def __init__(self, cv: float, rng: np.random.Generator, rho: float = 0.9):
        self.cv = cv
        self.rho = rho
        self.rng = rng
        self._state = 0.0

    def next_factor(self) -> float:
        if self.cv <= 0:
            return 1.0
        self._state = self.rho * self._state + self.rng.normal(0.0, self.cv)
        return float(np.exp(self._state))


def chemical_noise(rng: np.random.Generator, n_peaks: int, level: float,
                   mz_range: Tuple[float, float]) -> Tuple[np.ndarray, np.ndarray]:
    """Low-intensity chemical noise peaks for one scan."""
    if n_peaks <= 0:
        return np.array([]), np.array([])
    mz = rng.uniform(mz_range[0], mz_range[1], size=n_peaks)
    intensities = rng.exponential(level, size=n_peaks)
    return mz, intensities


def apply_saturation(intensities: np.ndarray, saturation: Optional[float]) -> np.ndarray:
    """Detector saturation: soft-compress the top decade, hard-clip at full scale."""
    if saturation is None:
        return intensities
    knee = 0.5 * saturation
    out = np.asarray(intensities, dtype=float).copy()
    over = out > knee
    # exponential approach to full scale above the knee: always below the
    # input, asymptotic to saturation
    span = saturation - knee
    out[over] = knee + span * (1.0 - np.exp(-(out[over] - knee) / span))
    return np.minimum(out, saturation)
