# simms — simulated mass spectrometry test data from the CLI

> Full documentation: [docs/APPROACH.md](docs/APPROACH.md) (scientific
> models), [docs/CLI.md](docs/CLI.md) (complete command reference),
> [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) (architecture & contributing),
> [AGENTS.md](AGENTS.md) (driving simms from coding agents).

`simms` is a single command-line tool that generates, mutates, combines and
validates simulated mass spectrometry data. It unifies the capabilities of
four ecosystems:

| Ecosystem | What simms uses it for |
|---|---|
| **MassBank-data** | source library of 139k+ real reference spectra to sample as simulation templates; MassBank record read/write/validate |
| **matchms** | spectral I/O (mgf/msp/mzML/mzXML/json/pickle), metadata harmonization, cleaning filters, merging |
| **psims / pyteomics** | writing standards-compliant indexed mzML for simulated LC-MS runs |
| **OpenMS** (optional) | native TOPP tool backend (FileMerger, IDFileConverter, FeatureLinker…) when an OpenMS install is on PATH |

Everything is reproducible: all randomness flows from `--seed`. Every command
supports `--json` for machine-readable output, which makes `simms` easy to
drive from coding agents (Claude Code, Cursor, Codex — see `AGENTS.md`).

## Install

```bash
pip install -e .          # from this directory
simms describe            # machine-readable capability manifest
```

## Generating simulated data

```bash
# Theoretical peptide fragment spectra (b/y ions), with noisy variants
simms generate peptides --sequences ELVISLIVESK,PEPTIDEK --charges 2,3 \
    --variants 3 --noise-preset noisy-qtof --seed 42 --out peptides.mgf

# Digest a FASTA in silico first (trypsin, missed cleavages, length limits)
simms generate peptides --fasta proteins.fasta --missed-cleavages 1 --out digest.msp

# Isotope-pattern MS1 spectra from molecular formulas
simms generate isotopes --formula C9H8O4 C16H30N2O3 --adduct "[M+H]+" --out iso.msp

# Sample real MassBank records and emit noisy simulated variants
simms generate from-massbank --massbank /path/to/MassBank-data \
    --n 20 --filter ms_type=MS2 --filter ionmode=positive \
    --variants 5 --noise-preset default --seed 7 --out simulated.mgf

# Simulate a full DDA LC-MS/MS run (MS1 isotope envelopes + top-N MS2) as mzML
simms generate lcms-run -i simulated.mgf --gradient 600 --peak-fwhm 10 \
    --ms1-interval 1 --top-n 3 --noise-preset clean-orbitrap --out run.mzML

# SWATH-style DIA instead of DDA, with profile-mode peaks at 60k resolution
simms generate lcms-run -i simulated.mgf --acquisition dia \
    --dia-range 400 1200 --dia-window 20 --profile --resolving-power 60000 \
    --out dia_run.mzML

# ML-predicted fragment intensities (Prosit via the Koina service, cached)
simms generate peptides --fasta proteome.fasta --fragment-model prosit \
    --collision-energy 28 --out prosit_library.mgf
```

Noise presets: `none`, `clean-orbitrap`, `default`, `noisy-qtof`, `harsh` —
each overridable per parameter (`--mz-ppm`, `--intensity-cv`, `--dropout`,
`--noise-peaks`, `--noise-peak-intensity`).

### Realism model

Peptide MS2 intensities use a **mobile-proton model** by default
(`--fragment-model realistic`): y-over-b dominance, suppressed b1/prominent
b2, the proline effect, enhanced D/E cleavage for non-mobile precursors,
-H2O/-NH3 neutral losses, immonium ions, fragment M+1 isotopes, and
collision-energy dependence (`--collision-energy`, higher CE shifts
intensity to low-mass fragments). `--fragment-model simple` restores the
clean b/y ladder.

`generate lcms-run` layers instrument and chromatography physics,
controlled by `--realism none|default|high` plus per-knob overrides:

| Feature | Knob |
|---|---|
| EMG (tailing) elution, width growth along gradient | `--tailing-tau`, `--rt-broadening` |
| MS1 isotope envelopes (formula or averagine from mass) | `--no-isotope-envelopes` |
| Electrospray charge-state envelopes (z−1, z, z+1) | `--no-charge-envelope` |
| Background contaminant ions (polysiloxanes, phthalates) | `--no-contaminants` |
| Chemical noise per scan / spray instability (AR(1)) | `--chemical-noise`, `--spray-cv` |
| Mass calibration drift over the run | `--drift-ppm` |
| Detector saturation (soft knee, hard clip) | `--saturation` |
| DDA dynamic exclusion | `--exclusion` |
| Co-isolation chimeric MS2 spectra | `--isolation-window`, `--no-chimeras` |
| Profile-mode peak shapes (resolution-scaled FWHM) | `--profile`, `--resolving-power` |
| DIA acquisition (SWATH windows, multiplexed MS2) | `--acquisition dia`, `--dia-range`, `--dia-window`, `--dia-overlap` |
| Prosit ML fragment intensities (Koina, cached) | `--fragment-model prosit` |

## Combining and converting

```bash
# Merge any mix of formats into one library, with filtering and dedup
simms merge -i a.mgf b.msp c.json run.mzML -o combined.mgf \
    --min-peaks 5 --ionmode positive --dedupe-key inchikey

# Convert between formats (mgf/msp/json/pickle, or a MassBank record dir)
simms convert -i combined.mgf -o combined.msp --export-style nist
simms convert -i simulated.mgf -o records_dir     # writes MSBNK-SIMMS-*.txt

# Clean with matchms default filters
simms clean -i combined.mgf -o cleaned.msp --normalize
```

## MassBank utilities

```bash
simms massbank stats --massbank /path/to/MassBank-data
simms massbank validate records_dir/MSBNK-SIMMS-SIM000001.txt
```

## OpenMS backend (optional)

```bash
simms openms list                       # discover installed TOPP tools
simms openms run FileMerger -- -in a.mzML b.mzML -out merged.mzML
simms openms run IDFileConverter -- -in proteins.fasta -out theoretical.mzML
```

Set `OPENMS_BIN=/path/to/OpenMS/bin` if the tools are not on PATH.

## Pipelines

Chain any commands declaratively in YAML and run them as one unit:

```yaml
name: simulated-test-dataset
steps:
  - generate peptides --sequences ELVISLIVESK --variants 2 --noise-preset default --out p1.mgf
  - run: generate from-massbank
    args: {massbank: /data/MassBank-data, n: 10, variants: 3, out: p2.mgf, filter: ms_type=MS2}
  - merge -i p1.mgf p2.mgf -o dataset.mgf
  - clean -i dataset.mgf -o dataset.msp --normalize
  - generate lcms-run -i dataset.msp --gradient 300 --out dataset_run.mzML
```

```bash
simms pipeline run pipeline.yaml --json     # step logs on stderr, JSON result on stdout
simms pipeline run pipeline.yaml --dry-run  # show resolved commands without running
```

## Tests

```bash
python -m pytest tests/ -v
```

The suite validates chemistry against literature values (aspirin masses, the
Br₂ isotope triplet, canonical PEPTIDE b/y ions), MassBank parsing and
round-tripping, seed reproducibility, and full mzML round-trips via
pyteomics. When a MassBank-data checkout exists at `/home/user/massbank-data`
the end-to-end sampling tests run against the real collection.
