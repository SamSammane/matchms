# simms — guidance for coding agents

This file is for AI coding agents (Claude Code, Cursor, Codex, …) driving the
`simms` CLI or editing this package.

## Driving the CLI

- Start with `simms describe` — it prints a JSON capability manifest:
  commands, supported read/write formats, noise presets, and which OpenMS
  TOPP tools are installed.
- Add `--json` to any command for a single machine-readable result object on
  stdout. For `pipeline run --json`, step logs go to stderr; stdout carries
  only the final JSON.
- Always pass `--seed <int>` when generating data that must be reproducible.
  Identical seed + arguments ⇒ identical output, byte-for-byte spectra.
- Exit codes: 0 success, 1 validation/pipeline failure, 2 usage or missing
  file errors (message on stderr).

## Recipes

Simulated MS2 library from real MassBank spectra:

```bash
simms generate from-massbank --massbank <MassBank-data checkout> \
  --n 50 --filter ms_type=MS2 --variants 4 --noise-preset default \
  --seed 1 --out sim_library.mgf --json
```

Full simulated LC-MS/MS run for pipeline testing (valid indexed mzML,
readable by pyteomics, OpenMS, ProteoWizard):

```bash
simms generate lcms-run -i sim_library.mgf --gradient 600 --out run.mzML --seed 1 --json
```

Realism defaults to the `default` preset (EMG tailing, isotope/charge
envelopes, contaminants, chemical noise, spray flicker, calibration drift,
saturation, chimeric MS2, dynamic exclusion). Use `--realism none` for
idealized textbook output, `--realism high` for a stress-test instrument,
and the per-knob overrides listed in `simms generate lcms-run --help` to
isolate a single effect. Peptide MS2 spectra use the mobile-proton
intensity model by default; `--fragment-model simple` gives the plain b/y
ladder. The JSON result reports `precursor_species` (charge-envelope
expansion) and `chimeric_ms2_scans` so agents can assert on realism
behavior without parsing the mzML.

Combine heterogeneous libraries: `simms merge -i a.mgf b.msp c.json -o all.mgf --json`.

Multi-step dataset builds belong in a YAML pipeline (`simms pipeline run f.yaml`)
rather than ad-hoc shell chains — the pipeline file documents the dataset and
re-runs deterministically.

## Filter keys for --filter (from-massbank)

Substring matches on parsed record metadata: `ms_type` (MS, MS2, MS3…),
`ionmode` (positive/negative), `instrument_type` (e.g. LC-ESI-QTOF, EI-B),
`formula`, `inchikey`, `compound_name`, `collision_energy`. Repeatable;
all must match. Heavy filters scan many files — bound work with
`--scan-limit` when speed matters more than unbiased sampling.

## Editing this package

- Layout: `simms/chem.py` (formulas/isotopes), `peptides.py` (digestion +
  fragment ions), `noise.py` (noise models), `massbank.py` (record
  parse/sample/write/validate), `simulate.py` (variants + LC-MS runs),
  `io_utils.py` (format dispatch), `merging.py`, `openms_backend.py`,
  `pipeline.py`, `cli.py`.
- matchms metadata harmonization renames keys at Spectrum creation
  (`accession` → `spectrum_id`, `ms_type` → `ms_level`). Never store a key
  that harmonization maps onto another key you also set; check
  `matchms/data/known_key_conversions.csv` before adding metadata keys.
- Text formats stringify metadata; `io_utils._coerce_numeric_metadata`
  restores int/float types on load. Extend it when adding numeric keys.
- Run `python -m pytest tests/ -q` before committing; tests must stay green
  without a MassBank checkout (those tests skip) and with one.
- Keep every new capability exposed through `cli.py`, covered by an
  end-to-end test in `tests/test_cli_end_to_end.py`, and listed in the
  `describe` manifest.
