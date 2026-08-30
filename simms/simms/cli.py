"""simms command-line interface.

Designed to be driven both by humans and by coding agents (Claude Code,
Cursor, Codex): every command takes explicit options, honors --seed for
reproducibility, and with --json prints a single machine-readable result
object to stdout.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any, Dict, List, Optional, Sequence

from . import __version__
from .noise import NOISE_PRESETS, NoiseModel


def _add_noise_options(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("noise model")
    group.add_argument("--noise-preset", choices=sorted(NOISE_PRESETS),
                       default="none", help="base noise preset (default: none)")
    group.add_argument("--mz-ppm", type=float, help="override: gaussian m/z error sigma in ppm")
    group.add_argument("--intensity-cv", type=float, help="override: intensity noise CV")
    group.add_argument("--dropout", type=float, help="override: peak dropout probability")
    group.add_argument("--noise-peaks", type=int, help="override: number of spurious peaks")
    group.add_argument("--noise-peak-intensity", type=float,
                       help="override: spurious peak mean intensity as fraction of base peak")


def _noise_from_args(args: argparse.Namespace) -> NoiseModel:
    preset = NOISE_PRESETS[args.noise_preset]
    return NoiseModel(
        mz_ppm=preset.mz_ppm if args.mz_ppm is None else args.mz_ppm,
        intensity_cv=preset.intensity_cv if args.intensity_cv is None else args.intensity_cv,
        dropout=preset.dropout if args.dropout is None else args.dropout,
        noise_peaks=preset.noise_peaks if args.noise_peaks is None else args.noise_peaks,
        noise_intensity=(preset.noise_intensity if args.noise_peak_intensity is None
                         else args.noise_peak_intensity),
    )


def _add_realism_options(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("instrument & chromatography realism")
    group.add_argument("--realism", choices=["none", "default", "high"],
                       default="default", help="realism preset (default: default)")
    group.add_argument("--tailing-tau", type=float,
                       help="override: EMG tail constant as multiple of peak sigma")
    group.add_argument("--rt-broadening", type=float,
                       help="override: fractional peak-width growth across the gradient")
    group.add_argument("--drift-ppm", type=float,
                       help="override: peak-to-peak mass calibration drift over the run")
    group.add_argument("--chemical-noise", type=int,
                       help="override: chemical noise peaks per MS1 scan")
    group.add_argument("--spray-cv", type=float,
                       help="override: spray instability CV (correlated per-scan flicker)")
    group.add_argument("--saturation", type=float,
                       help="override: detector full scale (intensity clip level)")
    group.add_argument("--isolation-window", type=float,
                       help="override: MS2 isolation window width in m/z")
    group.add_argument("--exclusion", type=float,
                       help="override: dynamic exclusion time in seconds")
    group.add_argument("--no-charge-envelope", action="store_true",
                       help="disable electrospray charge-state envelopes")
    group.add_argument("--no-isotope-envelopes", action="store_true",
                       help="disable MS1 isotope envelopes (single precursor peak)")
    group.add_argument("--no-contaminants", action="store_true",
                       help="disable background contaminant ions")
    group.add_argument("--no-chimeras", action="store_true",
                       help="disable co-isolation chimeric MS2 spectra")


def _realism_from_args(args: argparse.Namespace):
    from dataclasses import replace
    from .realism import REALISM_PRESETS
    config = REALISM_PRESETS[args.realism]
    overrides = {}
    for arg_name, field_name in [
            ("tailing_tau", "tailing_tau_factor"),
            ("rt_broadening", "rt_broadening"),
            ("drift_ppm", "drift_ppm"),
            ("chemical_noise", "chemical_noise_peaks"),
            ("spray_cv", "spray_instability_cv"),
            ("saturation", "saturation"),
            ("isolation_window", "isolation_window"),
            ("exclusion", "dynamic_exclusion_seconds")]:
        value = getattr(args, arg_name, None)
        if value is not None:
            overrides[field_name] = value
    if args.no_charge_envelope:
        overrides["charge_envelope"] = False
    if args.no_isotope_envelopes:
        overrides["isotope_envelopes"] = False
    if args.no_contaminants:
        overrides["contaminants"] = False
    if args.no_chimeras:
        overrides["chimeras"] = False
    return replace(config, **overrides)


def _parse_filters(pairs: Optional[List[str]]) -> Dict[str, str]:
    filters = {}
    for pair in pairs or []:
        if "=" not in pair:
            raise SystemExit(f"--filter must be key=value, got {pair!r}")
        key, _, value = pair.partition("=")
        filters[key.strip()] = value.strip()
    return filters


def _emit(args: argparse.Namespace, result: Dict[str, Any]) -> None:
    if getattr(args, "json", False):
        print(json.dumps(result, indent=2, default=str))
    else:
        for key, value in result.items():
            print(f"{key}: {value}")


# ---------------------------------------------------------------- commands

def cmd_describe(args: argparse.Namespace) -> int:
    import shutil
    from . import openms_backend
    manifest = {
        "tool": "simms",
        "version": __version__,
        "purpose": "generate, mutate, combine and validate simulated mass spectrometry data",
        "reproducibility": "all randomness is controlled by --seed",
        "commands": {
            "generate peptides": "theoretical peptide fragment spectra from sequences or FASTA digestion",
            "generate isotopes": "MS1 isotope-pattern spectrum from a molecular formula",
            "generate from-massbank": "sample real MassBank records and emit noisy simulated variants",
            "generate lcms-run": "simulate a full DDA LC-MS/MS run and write indexed mzML",
            "merge": "combine spectra from any mix of supported files into one library",
            "convert": "convert between spectral file formats",
            "clean": "apply matchms default cleaning filters",
            "massbank stats": "record counts per contributor in a MassBank-data checkout",
            "massbank validate": "structural validation of MassBank record files",
            "openms list": "discover installed OpenMS TOPP tools",
            "openms run": "pass a command through to an OpenMS TOPP tool",
            "pipeline run": "execute a YAML pipeline chaining any of the above",
        },
        "read_formats": ["mgf", "msp", "mzML", "mzXML", "json", "pickle",
                         "MassBank record .txt / directory"],
        "write_formats": ["mgf", "msp", "json", "pickle",
                          "MassBank record directory", "mzML (lcms-run)"],
        "noise_presets": sorted(NOISE_PRESETS),
        "realism_presets": ["none", "default", "high"],
        "realism_features": [
            "EMG (tailing) elution profiles with gradient peak broadening",
            "MS1 isotope envelopes from formula or averagine model",
            "electrospray charge-state envelopes",
            "background contaminant ions (polysiloxanes, phthalates)",
            "per-scan chemical noise and AR(1) spray instability",
            "sinusoidal mass calibration drift",
            "detector saturation (soft knee + hard clip)",
            "DDA dynamic exclusion and co-isolation chimeric MS2",
            "mobile-proton MS2 intensity model with neutral losses, "
            "immonium ions, fragment isotopes and collision-energy dependence",
        ],
        "openms_tools_available": sorted(openms_backend.discover()),
    }
    print(json.dumps(manifest, indent=2))
    return 0


def cmd_generate_peptides(args: argparse.Namespace) -> int:
    from . import io_utils, peptides, simulate
    sequences: List[str] = []
    if args.sequences:
        sequences.extend(s.strip().upper() for s in args.sequences.split(",") if s.strip())
    if args.fasta:
        for _, protein in peptides.read_fasta(args.fasta):
            sequences.extend(peptides.digest(protein, enzyme=args.enzyme,
                                             missed_cleavages=args.missed_cleavages,
                                             min_length=args.min_length,
                                             max_length=args.max_length))
    if not sequences:
        raise SystemExit("no peptides: provide --sequences and/or --fasta")
    charges = [int(c) for c in args.charges.split(",")]
    spectra = [simulate.theoretical_peptide_spectrum(
                   seq, charge=z, ion_types=args.ion_types,
                   max_fragment_charge=args.max_fragment_charge,
                   model=args.fragment_model,
                   collision_energy=args.collision_energy)
               for seq in sequences for z in charges]
    noise = _noise_from_args(args)
    if args.variants > 0:
        spectra = simulate.simulate_variants(spectra, noise, n_variants=args.variants,
                                             seed=args.seed)
    written = io_utils.save_any(spectra, args.out)
    _emit(args, {"command": "generate peptides", "peptides": len(sequences),
                 "charges": charges, "variants": args.variants,
                 "fragment_model": args.fragment_model,
                 "spectra_written": written, "out": args.out, "seed": args.seed})
    return 0


def cmd_generate_isotopes(args: argparse.Namespace) -> int:
    from . import io_utils, simulate
    spectra = [simulate.isotope_spectrum(formula, adduct=args.adduct)
               for formula in args.formula]
    noise = _noise_from_args(args)
    if args.variants > 0:
        spectra = simulate.simulate_variants(spectra, noise, n_variants=args.variants,
                                             seed=args.seed)
    written = io_utils.save_any(spectra, args.out)
    _emit(args, {"command": "generate isotopes", "formulas": args.formula,
                 "adduct": args.adduct, "spectra_written": written,
                 "out": args.out, "seed": args.seed})
    return 0


def cmd_generate_from_massbank(args: argparse.Namespace) -> int:
    from . import io_utils, massbank, simulate
    records = massbank.sample_records(
        args.massbank, args.n, seed=args.seed,
        filters=_parse_filters(args.filter),
        contributors=args.contributors.split(",") if args.contributors else None,
        scan_limit=args.scan_limit)
    templates = [massbank.record_to_spectrum(r) for r in records]
    templates = [t for t in templates if t is not None]
    if not templates:
        raise SystemExit("no MassBank records matched the given filters")
    if args.variants > 0:
        noise = _noise_from_args(args)
        spectra = simulate.simulate_variants(templates, noise,
                                             n_variants=args.variants, seed=args.seed)
    else:
        spectra = templates
    written = io_utils.save_any(spectra, args.out)
    _emit(args, {"command": "generate from-massbank",
                 "templates_sampled": len(templates), "variants_per_template": args.variants,
                 "spectra_written": written, "out": args.out, "seed": args.seed})
    return 0


def cmd_generate_lcms_run(args: argparse.Namespace) -> int:
    from . import io_utils, massbank, simulate
    spectra = []
    if args.inputs:
        spectra.extend(io_utils.load_many(args.inputs))
    if args.massbank:
        records = massbank.sample_records(args.massbank, args.n, seed=args.seed,
                                          filters=_parse_filters(args.filter))
        spectra.extend(s for s in (massbank.record_to_spectrum(r) for r in records)
                       if s is not None)
    if not spectra:
        raise SystemExit("no input spectra: provide -i files and/or --massbank")
    noise = _noise_from_args(args) if args.noise_preset != "none" or args.mz_ppm else None
    result = simulate.simulate_lcms_run(
        spectra, args.out, gradient_seconds=args.gradient,
        peak_fwhm_seconds=args.peak_fwhm, ms1_interval_seconds=args.ms1_interval,
        dda_top_n=args.top_n, noise=noise, realism=_realism_from_args(args),
        seed=args.seed)
    result["command"] = "generate lcms-run"
    result["realism"] = args.realism
    _emit(args, result)
    return 0


def cmd_merge(args: argparse.Namespace) -> int:
    from . import merging
    result = merging.merge_files(args.inputs, args.out, dedupe_key=args.dedupe_key,
                                 min_peaks=args.min_peaks, ms_level=args.ms_level,
                                 ionmode=args.ionmode, export_style=args.export_style)
    result["command"] = "merge"
    result["out"] = args.out
    _emit(args, result)
    return 0


def cmd_convert(args: argparse.Namespace) -> int:
    from . import io_utils
    spectra = io_utils.load_any(args.input)
    written = io_utils.save_any(spectra, args.out, export_style=args.export_style)
    _emit(args, {"command": "convert", "input": args.input, "out": args.out,
                 "spectra_written": written})
    return 0


def cmd_clean(args: argparse.Namespace) -> int:
    from matchms.filtering import default_filters, normalize_intensities
    from . import io_utils
    spectra = io_utils.load_any(args.input)
    cleaned = []
    for spectrum in spectra:
        spectrum = default_filters(spectrum)
        if spectrum is not None and args.normalize:
            spectrum = normalize_intensities(spectrum)
        if spectrum is not None:
            cleaned.append(spectrum)
    written = io_utils.save_any(cleaned, args.out)
    _emit(args, {"command": "clean", "input": args.input, "loaded": len(spectra),
                 "written": written, "out": args.out})
    return 0


def cmd_massbank_stats(args: argparse.Namespace) -> int:
    from . import massbank
    result = massbank.repo_stats(args.massbank)
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"path: {result['path']}")
        print(f"total_records: {result['total_records']}")
        for name, count in sorted(result["contributors"].items(),
                                  key=lambda kv: kv[1], reverse=True):
            print(f"  {name}: {count}")
    return 0


def cmd_massbank_validate(args: argparse.Namespace) -> int:
    from . import massbank
    failures = 0
    report = []
    for path in args.files:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            issues = massbank.validate_record_text(handle.read())
        report.append({"file": path, "ok": not issues, "issues": issues})
        failures += bool(issues)
    _emit(args, {"command": "massbank validate", "files": len(args.files),
                 "failed": failures, "report": report})
    return 1 if failures else 0


def cmd_openms_list(args: argparse.Namespace) -> int:
    from . import openms_backend
    found = openms_backend.discover()
    _emit(args, {"command": "openms list", "available": found,
                 "known_tools": openms_backend.KNOWN_TOOLS,
                 "hint": "install OpenMS or set OPENMS_BIN to enable native TOPP backends"})
    return 0


def cmd_openms_run(args: argparse.Namespace) -> int:
    from . import openms_backend
    return openms_backend.run_tool(args.tool, args.tool_args)


def cmd_pipeline_run(args: argparse.Namespace) -> int:
    from . import pipeline
    result = pipeline.run_pipeline(args.file, dry_run=args.dry_run,
                                   quiet_steps=args.json)
    _emit(args, result)
    return 0 if result["ok"] else 1


# ---------------------------------------------------------------- parser

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="simms",
        description="Simulated mass spectrometry test data: generate, mutate, "
                    "combine and validate spectra and LC-MS runs from the CLI.")
    parser.add_argument("--version", action="version", version=f"simms {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    def common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--seed", type=int, default=0, help="random seed (default 0)")
        p.add_argument("--json", action="store_true", help="print machine-readable JSON result")

    p = sub.add_parser("describe", help="machine-readable capability manifest")
    p.set_defaults(func=cmd_describe)

    generate = sub.add_parser("generate", help="generate simulated spectra / runs")
    gen_sub = generate.add_subparsers(dest="mode", required=True)

    p = gen_sub.add_parser("peptides", help="theoretical peptide fragment spectra")
    p.add_argument("--sequences", help="comma-separated peptide sequences")
    p.add_argument("--fasta", help="FASTA file to digest in silico")
    p.add_argument("--enzyme", default="trypsin", help="digestion enzyme (default trypsin)")
    p.add_argument("--missed-cleavages", type=int, default=0)
    p.add_argument("--min-length", type=int, default=6)
    p.add_argument("--max-length", type=int, default=40)
    p.add_argument("--charges", default="2", help="comma-separated precursor charges")
    p.add_argument("--ion-types", default="by",
                   help="fragment ion series for the simple model, e.g. by, aby")
    p.add_argument("--max-fragment-charge", type=int, default=1)
    p.add_argument("--fragment-model", choices=["simple", "realistic"],
                   default="realistic",
                   help="intensity model: simple b/y ladder, or mobile-proton "
                        "realistic model with losses/immonium/isotopes (default)")
    p.add_argument("--collision-energy", type=float, default=25.0,
                   help="normalized collision energy for the realistic model")
    p.add_argument("--variants", type=int, default=0,
                   help="noisy variants per spectrum (0 = ideal spectra)")
    p.add_argument("--out", required=True)
    common(p)
    _add_noise_options(p)
    p.set_defaults(func=cmd_generate_peptides)

    p = gen_sub.add_parser("isotopes", help="isotope-pattern MS1 spectra from formulas")
    p.add_argument("--formula", nargs="+", required=True, help="molecular formula(s), e.g. C9H8O4")
    p.add_argument("--adduct", default="[M+H]+")
    p.add_argument("--variants", type=int, default=0)
    p.add_argument("--out", required=True)
    common(p)
    _add_noise_options(p)
    p.set_defaults(func=cmd_generate_isotopes)

    p = gen_sub.add_parser("from-massbank",
                           help="simulate variants of real MassBank spectra")
    p.add_argument("--massbank", required=True, help="path to a MassBank-data checkout")
    p.add_argument("--n", type=int, default=10, help="number of template records to sample")
    p.add_argument("--filter", action="append",
                   help="metadata substring filter key=value (repeatable), "
                        "e.g. ms_type=MS2 ionmode=positive instrument_type=QTOF")
    p.add_argument("--contributors", help="comma-separated contributor directories")
    p.add_argument("--scan-limit", type=int, help="cap files scanned when filtering")
    p.add_argument("--variants", type=int, default=1,
                   help="simulated variants per template (0 = raw copies)")
    p.add_argument("--out", required=True)
    common(p)
    _add_noise_options(p)
    p.set_defaults(func=cmd_generate_from_massbank)

    p = gen_sub.add_parser("lcms-run", help="simulate a DDA LC-MS/MS run as mzML")
    p.add_argument("-i", "--inputs", nargs="*", default=[],
                   help="spectral files whose entries become eluting compounds")
    p.add_argument("--massbank", help="also sample compounds from a MassBank-data checkout")
    p.add_argument("--n", type=int, default=10, help="records to sample from --massbank")
    p.add_argument("--filter", action="append", help="MassBank metadata filter key=value")
    p.add_argument("--gradient", type=float, default=600.0, help="gradient length in seconds")
    p.add_argument("--peak-fwhm", type=float, default=10.0, help="chromatographic FWHM in seconds")
    p.add_argument("--ms1-interval", type=float, default=1.0, help="MS1 scan interval in seconds")
    p.add_argument("--top-n", type=int, default=3, help="DDA precursors per MS1 scan")
    p.add_argument("--out", required=True, help="output mzML path")
    common(p)
    _add_noise_options(p)
    _add_realism_options(p)
    p.set_defaults(func=cmd_generate_lcms_run)

    p = sub.add_parser("merge", help="combine spectral files into one library")
    p.add_argument("-i", "--inputs", nargs="+", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--dedupe-key", help="metadata key to deduplicate on, e.g. inchikey")
    p.add_argument("--min-peaks", type=int, default=0)
    p.add_argument("--ms-level", type=int)
    p.add_argument("--ionmode", choices=["positive", "negative"])
    p.add_argument("--export-style", default="matchms",
                   choices=["matchms", "massbank", "nist", "riken", "gnps"])
    common(p)
    p.set_defaults(func=cmd_merge)

    p = sub.add_parser("convert", help="convert between spectral formats")
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--export-style", default="matchms",
                   choices=["matchms", "massbank", "nist", "riken", "gnps"])
    common(p)
    p.set_defaults(func=cmd_convert)

    p = sub.add_parser("clean", help="apply matchms default cleaning filters")
    p.add_argument("-i", "--input", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--normalize", action="store_true", help="also normalize intensities")
    common(p)
    p.set_defaults(func=cmd_clean)

    mb = sub.add_parser("massbank", help="work with a MassBank-data checkout")
    mb_sub = mb.add_subparsers(dest="mb_command", required=True)
    p = mb_sub.add_parser("stats", help="record counts per contributor")
    p.add_argument("--massbank", required=True)
    common(p)
    p.set_defaults(func=cmd_massbank_stats)
    p = mb_sub.add_parser("validate", help="structurally validate record files")
    p.add_argument("files", nargs="+")
    common(p)
    p.set_defaults(func=cmd_massbank_validate)

    om = sub.add_parser("openms", help="OpenMS TOPP tool backend")
    om_sub = om.add_subparsers(dest="om_command", required=True)
    p = om_sub.add_parser("list", help="discover installed TOPP tools")
    common(p)
    p.set_defaults(func=cmd_openms_list)
    p = om_sub.add_parser("run", help="run a TOPP tool, passing args through")
    p.add_argument("tool")
    p.add_argument("tool_args", nargs=argparse.REMAINDER)
    p.set_defaults(func=cmd_openms_run)

    pl = sub.add_parser("pipeline", help="run a YAML pipeline of simms commands")
    pl_sub = pl.add_subparsers(dest="pl_command", required=True)
    p = pl_sub.add_parser("run")
    p.add_argument("file")
    p.add_argument("--dry-run", action="store_true")
    common(p)
    p.set_defaults(func=cmd_pipeline_run)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except BrokenPipeError:
        return 0  # stdout closed early (e.g. piped through head)
    except (FileNotFoundError, ValueError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
