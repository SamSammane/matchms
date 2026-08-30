"""Simulation of spectra and whole LC-MS runs.

Generation modes:

- theoretical peptide fragment spectra (simple deterministic model, or a
  mobile-proton realistic model with neutral losses / immonium ions),
- isotope-pattern MS1 spectra from molecular formulas,
- template-based simulation: real MassBank spectra perturbed by a noise
  model into any number of simulated variants.

The LC-MS run simulator places compounds on a gradient with exponentially
modified gaussian elution, models MS1 isotope envelopes (formula-based or
averagine), electrospray charge-state envelopes, background contaminant
ions, chemical noise, spray instability, calibration drift, detector
saturation, and data-dependent MS2 with dynamic exclusion and co-isolation
chimeras — then writes a standards-compliant indexed mzML via psims.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from matchms import Spectrum

from . import chem, peptides
from .noise import NoiseModel
from .realism import (REALISM_PRESETS, RealismConfig, SprayStability,
                      apply_saturation, calibration_drift_ppm, chemical_noise,
                      emg_profile, CONTAMINANT_IONS)


def theoretical_peptide_spectrum(sequence: str, charge: int = 2,
                                 ion_types: str = "by",
                                 max_fragment_charge: int = 1,
                                 model: str = "simple",
                                 collision_energy: float = 25.0) -> Spectrum:
    """Theoretical MS2 spectrum of a peptide.

    model="simple": clean b/y ladder with a positional intensity heuristic.
    model="realistic": mobile-proton intensity model with isotope peaks,
    neutral losses, immonium ions and collision-energy dependence.
    """
    if model == "realistic":
        ions = peptides.realistic_fragment_ions(
            sequence, precursor_charge=charge, collision_energy=collision_energy)
    elif model == "simple":
        ions = peptides.fragment_ions(sequence, ion_types=ion_types,
                                      max_fragment_charge=max_fragment_charge)
    else:
        raise ValueError(f"unknown fragment model {model!r}; use simple|realistic")
    mz = np.array([ion.mz for ion in ions])
    intensities = np.array([ion.intensity for ion in ions])
    order = np.argsort(mz)
    metadata = {
        "compound_name": sequence,
        "peptide_sequence": sequence,
        "charge": charge,
        "precursor_mz": peptides.precursor_mz(sequence, charge),
        "parent_mass": peptides.peptide_mass(sequence),
        "ms_level": 2,
        "ionmode": "positive",
        "simulated": f"theoretical-peptide-{model}",
        "fragment_model": model,
        "collision_energy": collision_energy,
        "ion_annotations": ",".join(ions[i].ion_type for i in order),
    }
    return Spectrum(mz=mz[order], intensities=intensities[order],
                    metadata=metadata, metadata_harmonization=False)


def isotope_spectrum(formula: str, adduct: str = "[M+H]+",
                     compound_name: Optional[str] = None) -> Spectrum:
    envelope = chem.isotope_envelope_mz(formula, adduct=adduct)
    mz = np.array([peak[0] for peak in envelope])
    intensities = np.array([peak[1] * 100.0 for peak in envelope])
    charge = chem.ADDUCTS[adduct][1]
    metadata = {
        "compound_name": compound_name or formula,
        "formula": formula,
        "adduct": adduct,
        "charge": charge,
        "ionmode": "positive" if charge > 0 else "negative",
        "precursor_mz": float(mz[np.argmax(intensities)]),
        "parent_mass": chem.monoisotopic_mass(formula),
        "ms_level": 1,
        "simulated": "isotope-pattern",
    }
    return Spectrum(mz=mz, intensities=intensities, metadata=metadata,
                    metadata_harmonization=False)


def simulate_variants(templates: Sequence[Spectrum], noise: NoiseModel,
                      n_variants: int = 1, seed: int = 0) -> List[Spectrum]:
    """Perturb template spectra into simulated variants (reproducible by seed)."""
    rng = np.random.default_rng(seed)
    out: List[Spectrum] = []
    for template in templates:
        if template is None:
            continue
        for variant in range(n_variants):
            mz, intensities = noise.apply(template.peaks.mz,
                                          template.peaks.intensities, rng)
            metadata = dict(template.metadata)
            # matchms harmonization stores MassBank accessions as spectrum_id
            template_id = metadata.pop("accession", None) or metadata.pop("spectrum_id", None)
            metadata["simulated"] = "massbank-template" if template_id else "variant"
            metadata["variant"] = variant
            if template_id:
                metadata["template_accession"] = template_id
            out.append(Spectrum(mz=mz, intensities=intensities, metadata=metadata,
                                metadata_harmonization=False))
    return out


# ------------------------------------------------------------------ LC-MS


@dataclass
class PrecursorSpecies:
    """One observable precursor: a compound at one charge state."""
    compound_index: int
    charge: int
    mz: float
    weight: float                      # fraction of the compound's signal
    envelope: List[Tuple[float, float]]


@dataclass
class SimulatedCompound:
    name: str
    rt_seconds: float
    sigma_seconds: float
    tau_seconds: float
    abundance: float
    ms2_mz: np.ndarray
    ms2_intensities: np.ndarray
    species: List[PrecursorSpecies] = field(default_factory=list)


def _neutral_mass(metadata: Dict, precursor_mz: float, charge: int) -> Optional[float]:
    mass = metadata.get("parent_mass")
    if isinstance(mass, (int, float)) and mass > 0:
        return float(mass)
    if precursor_mz and charge:
        return float(precursor_mz * charge - charge * chem.PROTON)
    return None


def _charge_states(base_charge: int, is_peptide: bool,
                   realism: RealismConfig) -> List[Tuple[int, float]]:
    if not realism.charge_envelope or not is_peptide or base_charge < 2:
        return [(base_charge, 1.0)]
    states = [(base_charge, 1.0), (base_charge + 1, 0.35)]
    if base_charge - 1 >= 1:
        states.append((base_charge - 1, 0.25))
    return states


def _compound_from_spectrum(index: int, spectrum: Spectrum,
                            rng: np.random.Generator,
                            gradient_seconds: float,
                            peak_fwhm_seconds: float,
                            realism: RealismConfig) -> Optional[SimulatedCompound]:
    metadata = spectrum.metadata
    mz_arr = spectrum.peaks.mz
    if len(mz_arr) == 0:
        return None
    precursor = metadata.get("precursor_mz")
    if not isinstance(precursor, (int, float)):
        precursor = float(mz_arr.max())
    charge = metadata.get("charge", 1)
    try:
        charge = abs(int(charge)) or 1
    except (TypeError, ValueError):
        charge = 1

    rt = metadata.get("retention_time")
    if not isinstance(rt, (int, float)) or not 0 < rt < gradient_seconds:
        rt = float(rng.uniform(0.05, 0.95) * gradient_seconds)
    sigma0 = peak_fwhm_seconds / 2.3548
    sigma = sigma0 * (1.0 + realism.rt_broadening * rt / gradient_seconds)
    tau = realism.tailing_tau_factor * sigma

    neutral = _neutral_mass(metadata, float(precursor), charge)
    formula = metadata.get("formula")
    is_peptide = bool(metadata.get("peptide_sequence"))

    species: List[PrecursorSpecies] = []
    for z, weight in _charge_states(charge, is_peptide, realism):
        if neutral is not None:
            species_mz = (neutral + z * chem.PROTON) / z
        elif z == charge:
            species_mz = float(precursor)
        else:
            continue  # cannot place other charge states without a mass
        envelope: List[Tuple[float, float]] = []
        if realism.isotope_envelopes:
            if formula:
                try:
                    adduct = metadata.get("adduct", "")
                    adduct = adduct if adduct in chem.ADDUCTS and abs(chem.ADDUCTS[adduct][1]) == z else None
                    if adduct:
                        envelope = chem.isotope_envelope_mz(formula, adduct=adduct)
                    else:
                        pattern = chem.isotope_pattern(formula)
                        mono = pattern[0][0]
                        envelope = [(species_mz + (m - mono) / z, a) for m, a in pattern]
                except ValueError:
                    envelope = []
            elif neutral is not None:
                envelope = chem.averagine_envelope_mz(neutral, z)
        if not envelope:
            envelope = [(species_mz, 1.0)]
        species.append(PrecursorSpecies(index, z, species_mz, weight, envelope))

    total_weight = sum(s.weight for s in species)
    for s in species:
        s.weight /= total_weight

    return SimulatedCompound(
        name=str(metadata.get("compound_name", f"compound{index}")),
        rt_seconds=float(rt), sigma_seconds=float(sigma), tau_seconds=float(tau),
        abundance=float(rng.uniform(1e5, 1e8)),
        ms2_mz=np.asarray(mz_arr, dtype=float),
        ms2_intensities=np.asarray(spectrum.peaks.intensities, dtype=float),
        species=species,
    )


def _apply_drift(mz: np.ndarray, ppm: float) -> np.ndarray:
    return mz * (1.0 + ppm * 1e-6) if ppm else mz


def simulate_lcms_run(spectra: Sequence[Spectrum], output_path: str,
                      gradient_seconds: float = 600.0,
                      peak_fwhm_seconds: float = 10.0,
                      ms1_interval_seconds: float = 1.0,
                      dda_top_n: int = 3,
                      ms2_trigger_fraction: float = 0.05,
                      noise: Optional[NoiseModel] = None,
                      realism: Optional[RealismConfig] = None,
                      seed: int = 0) -> Dict[str, object]:
    """Simulate a DDA LC-MS/MS run from a list of (fragment) spectra.

    Each input spectrum becomes one eluting compound (EMG elution profile,
    one or more charge states, isotope envelopes). MS1 scans are emitted
    every ``ms1_interval_seconds``; the ``dda_top_n`` most intense precursor
    species above threshold trigger MS2 scans (with co-isolation chimeras
    and dynamic exclusion). The run is written to ``output_path`` as
    indexed mzML.
    """
    from psims.mzml import MzMLWriter  # deferred: import cost & optional at runtime

    if realism is None:
        realism = REALISM_PRESETS["default"]

    rng = np.random.default_rng(seed)
    compounds = [c for c in (_compound_from_spectrum(i, s, rng, gradient_seconds,
                                                     peak_fwhm_seconds, realism)
                             for i, s in enumerate(spectra) if s is not None)
                 if c is not None]
    if not compounds:
        raise ValueError("no usable input spectra for LC-MS simulation")

    times = np.arange(0.0, gradient_seconds, ms1_interval_seconds)
    # precompute peak-normalized elution profiles (compounds x scans)
    profiles = np.vstack([
        c.abundance * emg_profile(times, c.rt_seconds, c.sigma_seconds, c.tau_seconds)
        for c in compounds
    ])

    all_species = [s for c in compounds for s in c.species]
    max_abundance = max(c.abundance for c in compounds)
    contaminant_base = 0.005 * max_abundance
    spray = SprayStability(realism.spray_instability_cv, rng)
    drift_phase = float(rng.uniform(0, 2 * np.pi))
    mz_lo = max(50.0, min(s.mz for s in all_species) - 100.0)
    mz_hi = max(s.mz for s in all_species) + 200.0

    scans: List[Dict[str, object]] = []
    scan_id = 0
    ms1_count = ms2_count = chimera_count = 0
    last_triggered: Dict[Tuple[int, int], float] = {}

    for scan_index, t in enumerate(times):
        levels = profiles[:, scan_index]
        spray_factor = spray.next_factor()
        drift = calibration_drift_ppm(t, gradient_seconds, realism.drift_ppm, drift_phase)

        mz_list: List[float] = []
        int_list: List[float] = []
        species_level: List[float] = []
        for s in all_species:
            level = levels[s.compound_index] * s.weight * spray_factor
            species_level.append(level)
            if level < 1.0:
                continue
            for iso_mz, iso_ab in s.envelope:
                mz_list.append(iso_mz)
                int_list.append(level * iso_ab)
        if realism.contaminants:
            jitter = rng.normal(1.0, 0.05, size=len(CONTAMINANT_IONS))
            for (c_mz, c_rel), j in zip(CONTAMINANT_IONS, jitter):
                mz_list.append(c_mz)
                int_list.append(contaminant_base * c_rel * max(0.1, j) * spray_factor)
        noise_mz, noise_int = chemical_noise(
            rng, realism.chemical_noise_peaks,
            realism.chemical_noise_level * max_abundance, (mz_lo, mz_hi))
        mz_arr = np.concatenate([np.asarray(mz_list), noise_mz])
        int_arr = np.concatenate([np.asarray(int_list), noise_int])
        mz_arr = _apply_drift(mz_arr, drift)
        if noise is not None and mz_arr.size:
            mz_arr, int_arr = noise.apply(mz_arr, int_arr, rng)
        int_arr = apply_saturation(int_arr, realism.saturation)
        order = np.argsort(mz_arr)
        scan_id += 1
        ms1_count += 1
        ms1_id = f"scan={scan_id}"
        scans.append({"id": ms1_id, "ms_level": 1, "time": t,
                      "mz": mz_arr[order], "intensities": int_arr[order]})

        # --- data-dependent MS2 on the top-N species of this MS1 scan
        eligible = []
        for s, level in zip(all_species, species_level):
            compound = compounds[s.compound_index]
            if level < ms2_trigger_fraction * compound.abundance * s.weight:
                continue
            key = (s.compound_index, s.charge)
            if t - last_triggered.get(key, -1e12) < realism.dynamic_exclusion_seconds:
                continue
            eligible.append((s, level))
        eligible.sort(key=lambda pair: pair[1], reverse=True)

        for s, level in eligible[:dda_top_n]:
            compound = compounds[s.compound_index]
            scale = level / (compound.abundance * s.weight)
            frag_mz = compound.ms2_mz.copy()
            frag_int = compound.ms2_intensities * level
            is_chimera = False
            if realism.chimeras:
                half_window = realism.isolation_window / 2.0
                for other, other_level in zip(all_species, species_level):
                    # co-isolation matters only when comparable to the target
                    if other is s or other_level < max(1.0, 0.02 * level):
                        continue
                    if abs(other.mz - s.mz) <= half_window:
                        other_compound = compounds[other.compound_index]
                        frag_mz = np.concatenate([frag_mz, other_compound.ms2_mz])
                        frag_int = np.concatenate([
                            frag_int, other_compound.ms2_intensities * other_level])
                        is_chimera = True
            if is_chimera:
                chimera_count += 1
            frag_mz = _apply_drift(frag_mz, drift)
            if noise is not None and frag_mz.size:
                frag_mz, frag_int = noise.apply(frag_mz, frag_int, rng)
            frag_int = apply_saturation(frag_int, realism.saturation)
            if frag_mz.size == 0:
                continue
            order = np.argsort(frag_mz)
            scan_id += 1
            ms2_count += 1
            scans.append({
                "id": f"scan={scan_id}", "ms_level": 2,
                "time": t + ms1_interval_seconds * 0.3,
                "mz": frag_mz[order], "intensities": frag_int[order],
                "precursor": {"mz": s.mz, "charge": s.charge,
                              "intensity": level, "parent_id": ms1_id},
            })
            last_triggered[(s.compound_index, s.charge)] = t

    with MzMLWriter(open(output_path, "wb"), close=True) as writer:
        writer.controlled_vocabularies()
        writer.file_description(["MS1 spectrum", "MSn spectrum", "centroid spectrum"])
        writer.software_list([{"id": "simms", "version": "0.1.0",
                               "params": ["custom unreleased software tool"]}])
        instrument = writer.InstrumentConfiguration(
            id="IC1", component_list=writer.ComponentList([
                writer.Source(order=1, params=["electrospray ionization"]),
                writer.Analyzer(order=2, params=["quadrupole"]),
                writer.Detector(order=3, params=["electron multiplier"]),
            ]))
        writer.instrument_configuration_list([instrument])
        writer.data_processing_list([
            writer.DataProcessing([
                writer.ProcessingMethod(order=1, software_reference="simms",
                                        params=["Conversion to mzML"])], id="DP1"),
        ])
        with writer.run(id="simulated_run", instrument_configuration="IC1"):
            with writer.spectrum_list(count=len(scans)):
                for scan in scans:
                    params = [
                        {"ms level": scan["ms_level"]},
                        "centroid spectrum",
                        "MS1 spectrum" if scan["ms_level"] == 1 else "MSn spectrum",
                    ]
                    precursor_information = None
                    if scan["ms_level"] == 2:
                        prec = scan["precursor"]
                        precursor_information = {
                            "mz": prec["mz"], "intensity": prec["intensity"],
                            "charge": prec["charge"], "scan_id": prec["parent_id"],
                            "activation": ["HCD", {"collision energy": 25.0}],
                        }
                    writer.write_spectrum(
                        scan["mz"], scan["intensities"],
                        id=scan["id"], params=params,
                        scan_start_time=scan["time"] / 60.0,
                        precursor_information=precursor_information,
                    )

    return {
        "output": output_path,
        "compounds": len(compounds),
        "precursor_species": len(all_species),
        "ms1_scans": ms1_count,
        "ms2_scans": ms2_count,
        "chimeric_ms2_scans": chimera_count,
        "total_scans": len(scans),
        "gradient_seconds": gradient_seconds,
        "seed": seed,
    }
