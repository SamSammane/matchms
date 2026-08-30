# Development guide

Architecture, design decisions, and how to extend simms.

## Module map

```
simms/
├── chem.py            isotope tables, formula parsing, monoisotopic masses,
│                      isotope patterns (convolution), adduct m/z, averagine
├── peptides.py        residue masses, digestion, b/y/a fragment m/z,
│                      simple + mobile-proton intensity models, immonium table
├── noise.py           NoiseModel (ppm error, intensity CV, dropout, spurious
│                      peaks) + presets; all randomness via a passed Generator
├── massbank.py        MassBank record parse/sample/stats/validate/write;
│                      record -> matchms Spectrum
├── mlfrag.py          Koina (Prosit) client: Triton wire format, response
│                      masking, local cache, KoinaUnavailable
├── realism.py         RealismConfig + presets, EMG profiles, contaminant
│                      table, chemical noise, spray AR(1), calibration drift,
│                      saturation, profile-mode rendering
├── simulate.py        spectrum-level generation (peptide/isotope/variants)
│                      and the LC-MS run engine (DDA + DIA) -> mzML via psims
├── io_utils.py        format dispatch on top of matchms + MassBank dirs;
│                      numeric metadata coercion on load
├── merging.py         multi-file merge with filters and dedup
├── openms_backend.py  TOPP tool discovery and passthrough (optional)
├── pipeline.py        YAML pipeline runner (re-enters the CLI per step)
└── cli.py             argparse surface; every capability is exposed here
```

Dependency policy: `chem`, `peptides`, `noise`, `realism`, `mlfrag` are
numpy/scipy/stdlib-only. matchms appears only in `simulate`, `io_utils`,
`merging`, `massbank`; psims only inside `simulate_lcms_run` (deferred
import). This keeps the physics testable without the I/O stack.

## Design decisions worth knowing

**Seed discipline.** No module calls `np.random` globally. Every stochastic
function takes a `numpy.random.Generator` (or a seed at the CLI boundary).
A whole LC-MS run is reproducible byte-for-byte from `--seed`; tests assert
this.

**matchms metadata harmonization.** matchms renames keys at `Spectrum`
creation (`accession` → `spectrum_id`, `ms_type` → `ms_level`, see
`matchms/data/known_key_conversions.csv`). Two rules follow: never store a
key that harmonization maps onto another key you also set, and never rely
on a key surviving a save/load round-trip without checking the conversion
table. This bit us once (`template_accession` logic); the fix and the rule
live in `simulate.simulate_variants`.

**Text formats stringify.** mgf/msp round-trips turn numbers into strings.
`io_utils._coerce_numeric_metadata` restores `ms_level`, `variant` (int)
and `retention_time`, `parent_mass`, `precursor_mz` (float) on load.
Extend it when adding numeric metadata keys.

**psims specifics.** `scan_start_time` is a float in minutes.
`precursor_information` accepts `isolation_window_args={"lower": off,
"target": mz, "upper": off}` (offsets, not absolute bounds). The
`centroided=` flag controls the centroid/profile CV param — do not also
put "centroid spectrum" in `params` or it duplicates.

**Koina client.** `mlfrag.predict` takes `_post_fn` for tests — unit tests
mock the transport and assert on the exact request/response wire format
rather than hitting the network. The cache key is
sha256(model|peptide|charge|ce) under `~/.cache/simms/koina`
(`SIMMS_CACHE` overrides). Never fall back silently from prosit to another
model: fidelity substitution must be the caller's explicit choice.

**Namespace collision.** Running pytest from a parent directory can
resolve `import simms` to the project directory (a namespace package)
instead of the package. `tests/conftest.py` pins the package directory
first on `sys.path`; keep it.

## The LC-MS run engine in one pass

`simulate_lcms_run` does, in order per scan tick:

1. Look up each compound's EMG elution level (precomputed matrix,
   compounds × scans).
2. Assemble MS1: per precursor species (compound × charge state), scale
   its isotope envelope by level × spray factor; append contaminants and
   chemical noise; apply calibration drift, the per-peak noise model,
   and saturation.
3. DDA: rank eligible species (threshold + dynamic exclusion), take
   top-N, build each MS2 from the compound's fragment array plus any
   co-isolated species (relative-intensity floor), same noise chain.
   DIA: iterate the fixed window plan, multiplex every species inside
   each window, always emitting the scan.
4. At write time only: optionally render centroids to profile mode, then
   hand psims the arrays with the right CV params and precursor/isolation
   metadata.

Adding a new physical effect usually means: a knob on `RealismConfig`
(+ presets), a pure function in `realism.py`, one call site in
`simulate_lcms_run`, a CLI flag in `_add_realism_options` /
`_realism_from_args`, a line in the `describe` manifest, and a test.

## Testing

```
python -m pytest tests/ -q        # 74 tests, ~17 s
```

- `test_chem.py`, `test_peptides.py` — literature-value chemistry.
- `test_noise_and_simulate.py` — noise invariants, seed reproducibility.
- `test_massbank.py` — record parsing, validation, round-trip.
- `test_realism.py` — each realism effect in isolation plus written-mzML
  assertions (XIC skew, contaminant presence, preset contrast).
- `test_frontier.py` — profile mode (FWHM vs resolution, doublet
  merging), Koina wire format/caching/failure modes, DIA cycles,
  isolation windows, multiplexing.
- `test_cli_end_to_end.py` — every command through `cli.main`, including
  runs against a real MassBank-data checkout when present (skipped
  otherwise).

Rules: tests must stay green both with and without the MassBank checkout
and without network. Anything touching randomness asserts determinism.
New CLI capability ⇒ end-to-end test + `describe` manifest entry.

## Release checklist

1. `python -m pytest tests/ -q` from the package dir *and* the repo root.
2. Bump `__version__` in `simms/__init__.py` and `pyproject.toml`.
3. Update README / docs / AGENTS.md for new flags.
4. Verify `simms describe` reflects reality (it is the agent contract).
