"""Simulation of spectra and whole LC-MS runs.

Three generation modes:

- theoretical peptide fragment spectra (from sequences or FASTA digestion),
- isotope-pattern MS1 spectra from molecular formulas,
- template-based simulation: real MassBank spectra perturbed by a noise
  model into any number of simulated variants.

Plus a full LC-MS run simulator that places compounds on a chromatographic
gradient (gaussian elution profiles), synthesizes MS1 scans with isotope
envelopes and data-dependent MS2 scans, and writes a standards-compliant
mzML file via psims.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from matchms import Spectrum

from . import chem, peptides
from .noise import NoiseModel


def theoretical_peptide_spectrum(sequence: str, charge: int = 2,
                                 ion_types: str = "by",
                                 max_fragment_charge: int = 1) -> Spectrum:
    ions = peptides.fragment_ions(sequence, ion_types=ion_types,
                                  max_fragment_charge=max_fragment_charge)
    mz = np.array([ion.mz for ion in ions])
    intensities = np.array([ion.intensity for ion in ions])
    metadata = {
        "compound_name": sequence,
        "peptide_sequence": sequence,
        "charge": charge,
        "precursor_mz": peptides.precursor_mz(sequence, charge),
        "parent_mass": peptides.peptide_mass(sequence),
        "ms_level": 2,
        "ionmode": "positive",
        "simulated": "theoretical-peptide",
        "ion_annotations": ",".join(ion.ion_type for ion in ions),
    }
    return Spectrum(mz=mz, intensities=intensities, metadata=metadata,
                    metadata_harmonization=False)


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


@dataclass
class SimulatedCompound:
    """One compound placed on the simulated gradient."""
    name: str
    precursor_mz: float
    rt_seconds: float
    peak_sigma_seconds: float
    abundance: float
    envelope: List[Tuple[float, float]]  # MS1 isotope envelope (mz, rel abundance)
    ms2_mz: np.ndarray
    ms2_intensities: np.ndarray
    charge: int = 1


def _compound_from_spectrum(spectrum: Spectrum, rng: np.random.Generator,
                            gradient_seconds: float,
                            peak_fwhm_seconds: float) -> Optional[SimulatedCompound]:
    metadata = spectrum.metadata
    precursor = metadata.get("precursor_mz")
    if precursor is None:
        mz_arr = spectrum.peaks.mz
        if len(mz_arr) == 0:
            return None
        precursor = float(mz_arr.max())
    rt = metadata.get("retention_time")
    if not isinstance(rt, (int, float)) or not 0 < rt < gradient_seconds:
        rt = float(rng.uniform(0.05, 0.95) * gradient_seconds)
    envelope: List[Tuple[float, float]] = []
    formula = metadata.get("formula")
    adduct = metadata.get("adduct", "[M+H]+")
    if formula:
        try:
            if adduct not in chem.ADDUCTS:
                adduct = "[M+H]+"
            envelope = chem.isotope_envelope_mz(formula, adduct=adduct)
        except ValueError:
            envelope = []
    if not envelope:
        envelope = [(float(precursor), 1.0)]
    sigma = peak_fwhm_seconds / 2.3548
    charge = metadata.get("charge", 1)
    try:
        charge = abs(int(charge)) or 1
    except (TypeError, ValueError):
        charge = 1
    return SimulatedCompound(
        name=str(metadata.get("compound_name", "unknown")),
        precursor_mz=float(precursor),
        rt_seconds=float(rt),
        peak_sigma_seconds=float(sigma),
        abundance=float(rng.uniform(1e5, 1e8)),
        envelope=envelope,
        ms2_mz=np.asarray(spectrum.peaks.mz, dtype=float),
        ms2_intensities=np.asarray(spectrum.peaks.intensities, dtype=float),
        charge=charge,
    )


def simulate_lcms_run(spectra: Sequence[Spectrum], output_path: str,
                      gradient_seconds: float = 600.0,
                      peak_fwhm_seconds: float = 10.0,
                      ms1_interval_seconds: float = 1.0,
                      dda_top_n: int = 3,
                      ms2_trigger_fraction: float = 0.05,
                      noise: Optional[NoiseModel] = None,
                      seed: int = 0) -> Dict[str, object]:
    """Simulate a DDA LC-MS/MS run from a list of (fragment) spectra.

    Each input spectrum becomes one eluting compound. MS1 scans are emitted
    every ``ms1_interval_seconds`` containing the isotope envelopes of all
    compounds scaled by their gaussian elution profile; the ``dda_top_n``
    most intense precursors above threshold trigger MS2 scans using the
    compound's fragment spectrum (noise applied per scan when a NoiseModel
    is given). The run is written to ``output_path`` as indexed mzML.
    """
    from psims.mzml import MzMLWriter  # deferred: import cost & optional at runtime

    rng = np.random.default_rng(seed)
    compounds = [c for c in (_compound_from_spectrum(s, rng, gradient_seconds, peak_fwhm_seconds)
                             for s in spectra if s is not None) if c is not None]
    if not compounds:
        raise ValueError("no usable input spectra for LC-MS simulation")

    times = np.arange(0.0, gradient_seconds, ms1_interval_seconds)
    scans: List[Dict[str, object]] = []
    scan_id = 0
    ms1_count = 0
    ms2_count = 0
    active_exclusion: Dict[str, float] = {}

    for t in times:
        profile = np.array([
            c.abundance * np.exp(-0.5 * ((t - c.rt_seconds) / c.peak_sigma_seconds) ** 2)
            for c in compounds
        ])
        mz_list: List[float] = []
        int_list: List[float] = []
        for compound, level in zip(compounds, profile):
            if level < 1.0:
                continue
            for iso_mz, iso_ab in compound.envelope:
                mz_list.append(iso_mz)
                int_list.append(level * iso_ab)
        mz_arr = np.array(mz_list)
        int_arr = np.array(int_list)
        if noise is not None and mz_arr.size:
            mz_arr, int_arr = noise.apply(mz_arr, int_arr, rng)
        order = np.argsort(mz_arr)
        scan_id += 1
        ms1_count += 1
        ms1_id = f"scan={scan_id}"
        scans.append({
            "id": ms1_id, "ms_level": 1, "time": t,
            "mz": mz_arr[order], "intensities": int_arr[order],
        })

        # data-dependent MS2 on the top-N precursors of this MS1 scan
        eligible = [(compound, level) for compound, level in zip(compounds, profile)
                    if level >= ms2_trigger_fraction * compound.abundance
                    and active_exclusion.get(compound.name, -1e9) + 3 * compound.peak_sigma_seconds < t]
        eligible.sort(key=lambda pair: pair[1], reverse=True)
        for compound, level in eligible[:dda_top_n]:
            frag_mz = compound.ms2_mz
            frag_int = compound.ms2_intensities * (level / compound.abundance)
            if noise is not None and frag_mz.size:
                frag_mz, frag_int = noise.apply(frag_mz, frag_int, rng)
            if frag_mz.size == 0:
                continue
            scan_id += 1
            ms2_count += 1
            scans.append({
                "id": f"scan={scan_id}", "ms_level": 2,
                "time": t + ms1_interval_seconds * 0.3,
                "mz": frag_mz, "intensities": frag_int,
                "precursor": {"mz": compound.precursor_mz, "charge": compound.charge,
                              "intensity": level, "parent_id": ms1_id},
            })
            active_exclusion[compound.name] = t

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
        "ms1_scans": ms1_count,
        "ms2_scans": ms2_count,
        "total_scans": len(scans),
        "gradient_seconds": gradient_seconds,
        "seed": seed,
    }
