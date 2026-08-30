# CLI reference

Every command supports `--json` (machine-readable result on stdout) and,
where randomness is involved, `--seed` (default 0; identical seed +
arguments ⇒ identical output). Exit codes: `0` success, `1` validation or
pipeline failure, `2` usage/missing-file/network errors (message on
stderr).

```
simms describe
simms generate peptides | isotopes | from-massbank | lcms-run
simms merge | convert | clean
simms massbank stats | validate
simms openms list | run
simms pipeline run
```

---

## simms describe

Prints the JSON capability manifest: commands, read/write formats, noise
and realism presets, acquisition schemes, fragment models, and which
OpenMS TOPP tools are installed. Agents should call this first.

---

## simms generate peptides

Theoretical peptide MS2 spectra from explicit sequences and/or in-silico
digestion of a FASTA file.

| Option | Default | Meaning |
|---|---|---|
| `--sequences` | — | comma-separated peptide sequences |
| `--fasta` | — | FASTA file to digest |
| `--enzyme` | trypsin | trypsin, trypsin/p, lys-c, arg-c, chymotrypsin, "no cleavage" |
| `--missed-cleavages` | 0 | allowed missed cleavages |
| `--min-length` / `--max-length` | 6 / 40 | peptide length bounds |
| `--charges` | 2 | comma-separated precursor charges |
| `--fragment-model` | realistic | `simple` (b/y ladder), `realistic` (mobile-proton model), `prosit` (Koina ML prediction) |
| `--collision-energy` | 25 | NCE for realistic/prosit models |
| `--ion-types` | by | series for the simple model (`by`, `aby`) |
| `--max-fragment-charge` | 1 | fragment charges for the simple model |
| `--variants` | 0 | noisy variants per spectrum (0 = ideal) |
| `--out` | required | `.mgf` `.msp` `.json` `.pickle` or a directory (MassBank records) |
| noise options | | see **Noise model** below |

```bash
simms generate peptides --fasta proteome.fasta --missed-cleavages 1 \
    --charges 2,3 --fragment-model prosit --collision-energy 28 \
    --out library.mgf --json
```

Prosit notes: requires network to https://koina.wilhelmlab.org on first
use per (peptide, charge, NCE); responses are cached in
`~/.cache/simms/koina` (override root with `SIMMS_CACHE`; override host
with `SIMMS_KOINA_URL`). Sequences ≤ 30 residues. On an unreachable
service the command exits 2 with guidance — the offline fallback is
`--fragment-model realistic`, chosen by you, never silently.

## simms generate isotopes

MS1 isotope-pattern spectra from molecular formulas.

| Option | Default | Meaning |
|---|---|---|
| `--formula` | required | one or more formulas (e.g. `C9H8O4`) |
| `--adduct` | `[M+H]+` | any of the supported adducts (see `chem.ADDUCTS`) |
| `--variants` | 0 | noisy variants per spectrum |
| `--out` | required | output file |

## simms generate from-massbank

Sample real records from a MassBank-data checkout and emit simulated
variants.

| Option | Default | Meaning |
|---|---|---|
| `--massbank` | required | path to a MassBank-data checkout |
| `--n` | 10 | template records to sample |
| `--filter` | — | repeatable `key=value` substring filters: `ms_type`, `ionmode`, `instrument_type`, `formula`, `inchikey`, `compound_name`, `collision_energy` |
| `--contributors` | — | comma-separated contributor directories |
| `--scan-limit` | — | cap files scanned when filtering (speed over unbiased sampling) |
| `--variants` | 1 | simulated variants per template (0 = raw copies) |
| `--out` | required | output file |

Each variant keeps `template_accession` — the ground-truth link back to
the real record.

## simms generate lcms-run

Simulate a full LC-MS/MS run and write indexed mzML.

| Option | Default | Meaning |
|---|---|---|
| `-i/--inputs` | — | spectral files; each spectrum becomes an eluting compound |
| `--massbank` + `--n` + `--filter` | — | additionally sample compounds from MassBank-data |
| `--gradient` | 600 | gradient length, seconds |
| `--peak-fwhm` | 10 | chromatographic peak FWHM, seconds |
| `--ms1-interval` | 1 | MS1 scan interval, seconds |
| `--acquisition` | dda | `dda` (top-N) or `dia` (SWATH-style windows) |
| `--top-n` | 3 | DDA precursors per MS1 scan |
| `--dia-range` | 400 1000 | DIA precursor m/z range |
| `--dia-window` | 25 | DIA isolation window width, m/z |
| `--dia-overlap` | 1 | overlap between adjacent windows, m/z |
| `--out` | required | output mzML path |

**Realism** (`--realism none|default|high` preset, then overrides):

| Override | Effect |
|---|---|
| `--tailing-tau` | EMG tail constant as multiple of peak sigma |
| `--rt-broadening` | fractional peak-width growth along the gradient |
| `--drift-ppm` | peak-to-peak calibration drift over the run |
| `--chemical-noise` | chemical-noise peaks per MS1 scan |
| `--spray-cv` | AR(1) spray instability CV |
| `--saturation` | detector full scale (soft knee + hard clip) |
| `--isolation-window` | DDA isolation window width, m/z |
| `--exclusion` | dynamic exclusion, seconds |
| `--no-charge-envelope` | single charge state per compound |
| `--no-isotope-envelopes` | single precursor peak per species |
| `--no-contaminants` | drop background ions |
| `--no-chimeras` | disable co-isolation mixing (DDA) |
| `--profile` | write profile-mode peak shapes |
| `--resolving-power` | resolving power at m/z 200 (implies `--profile`) |

JSON result fields: `compounds`, `precursor_species`, `ms1_scans`,
`ms2_scans`, `total_scans`, `profile_mode`, and per scheme
`chimeric_ms2_scans` (DDA) or `dia_windows_per_cycle` +
`multiplexed_ms2_scans` (DIA).

```bash
# DIA run with profile peaks at 60k resolution
simms generate lcms-run -i library.mgf --gradient 900 \
    --acquisition dia --dia-range 400 1200 --dia-window 20 \
    --profile --resolving-power 60000 --realism high \
    --out dia_run.mzML --seed 3 --json
```

---

## simms merge

Combine spectra from any mix of readable files into one library.

| Option | Default | Meaning |
|---|---|---|
| `-i/--inputs` | required | input files (mgf/msp/mzML/mzXML/json/pickle/MassBank .txt or dir) |
| `-o/--out` | required | output file |
| `--dedupe-key` | — | metadata key to dedup on (keeps most peaks), e.g. `inchikey` |
| `--min-peaks` | 0 | drop spectra with fewer peaks |
| `--ms-level` | — | keep only this MS level |
| `--ionmode` | — | `positive` or `negative` |
| `--export-style` | matchms | metadata style: matchms, massbank, nist, riken, gnps |

## simms convert

`-i input -o output [--export-style ...]` — same formats as merge. An
extensionless output (or `.massbank`) writes a directory of MassBank
records named `MSBNK-SIMMS-SIM*.txt`.

## simms clean

`-i input -o output [--normalize]` — matchms `default_filters`
(+ intensity normalization) over a library.

---

## simms massbank stats / validate

`stats --massbank PATH` — record counts per contributor.
`validate FILE...` — structural validation of MassBank records (required
fields, peak-count consistency, m/z ordering, terminator); exit 1 if any
file fails.

## simms openms list / run

`list` — discover installed OpenMS TOPP tools (PATH or `$OPENMS_BIN`).
`run TOOL -- args...` — verbatim passthrough, e.g.

```bash
simms openms run FileMerger -- -in a.mzML b.mzML -out merged.mzML
```

## simms pipeline run

`simms pipeline run pipeline.yaml [--dry-run] [--json]` — execute steps
declaratively; with `--json`, step logs go to stderr and stdout carries
only the final result object. Step forms:

```yaml
name: my-dataset
steps:
  - generate peptides --sequences ELVISLIVESK --out a.mgf     # string
  - ["merge", "-i", "a.mgf", "-o", "b.mgf"]                   # token list
  - run: generate lcms-run                                    # mapping
    args: {inputs: [b.mgf], gradient: 300, out: run.mzML, acquisition: dia}
```

---

## Environment variables

| Variable | Meaning |
|---|---|
| `OPENMS_BIN` | directory containing OpenMS TOPP executables |
| `SIMMS_KOINA_URL` | Koina host override for `--fragment-model prosit` |
| `SIMMS_CACHE` | cache root (default `~/.cache/simms`) |
