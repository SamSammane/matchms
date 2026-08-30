# Approach: how simms simulates mass spectrometry data

This document explains the scientific models behind simms — what is
simulated, how, and where the boundaries of realism lie.

## Why simulate

Testing LC-MS/MS software (search engines, feature finders, library-matching
tools, QC pipelines) requires data whose ground truth is known. Real data
has no ground truth; hand-crafted fixtures are trivial to pass. OpenMS
removed its MSSimulator in version 3.1, leaving no maintained CLI
simulator in the open-source ecosystem. simms fills that gap with a
generative model that is **seeded** (every output is exactly reproducible),
**layered** (each physical effect can be enabled in isolation), and
**grounded** (real MassBank reference spectra can serve as templates).

## The generation stack

### 1. Molecules and masses (`chem.py`)

Element isotope tables (masses and natural abundances) drive everything.
A molecular formula is parsed into element counts; the **isotope pattern**
is computed by convolving per-element distributions with
exponentiation-by-squaring, aggregating into coarse (nominal-mass) bins
with abundance-weighted centroid masses. Validated against literature:
aspirin monoisotopic 180.04226 Da, [M+H]+ 181.04954, and the 1:2:1
M/M+2/M+4 triplet of Br₂ compounds.

When only a mass is known (peptides without an explicit formula), the
**averagine model** (Senko et al.) estimates an elemental composition from
the average composition of peptides per 111.1254 Da, anchoring the
monoisotopic peak on the true mass. The M+1/M ratio it yields (~0.55 at
1000 Da, ~0.8 at 1500 Da) matches published envelopes.

### 2. Peptides (`peptides.py`)

In-silico digestion (trypsin with the proline rule, Lys-C, Arg-C,
chymotrypsin; configurable missed cleavages and length bounds) produces
peptides from FASTA. Fragment m/z values are exact monoisotopic b/y (and
a) ion masses.

Fragment **intensities** come from one of three models:

- **simple** — a positional bell curve with y > b. Clean, minimal,
  useful when tests must enumerate peaks.
- **realistic** — a heuristic implementation of the **mobile proton
  model**: precursor charge is weighed against the basic residue content
  (R > K > H) to classify proton mobility. Mobile protons give
  charge-directed backbone cleavage (strong y series); sequestered protons
  shift intensity to charge-remote pathways (enhanced cleavage C-terminal
  to Asp/Glu). Overlaid rules: b1 suppression, b2 prominence, the proline
  effect (enhanced y for cleavage N-terminal to Pro, suppressed cleavage
  C-terminal to Pro), −H₂O losses from S/T/E/D-containing fragments,
  −NH₃ from R/K/N/Q, immonium ions above ~25 NCE, fragment M+1 isotopes
  scaled by averagine carbon count, and collision-energy dependence
  (higher NCE erodes large fragments; verified: NCE 15→40 lowers the
  intensity-weighted mean fragment m/z by ~95).
- **prosit** — deep-learning predictions from the **Prosit** model
  (Gessulat et al., Nat. Methods 2019) served by the
  [Koina](https://koina.wilhelmlab.org) inference API
  (`Prosit_2020_intensity_HCD`). simms speaks Koina's Triton wire format,
  masks the impossible-ion sentinel values, and caches every
  (peptide, charge, NCE) prediction locally so repeated simulations are
  offline after the first call. When the service is unreachable the CLI
  fails with a clear message and the `realistic` model as the suggested
  fallback — no silent substitution.

### 3. Template-based simulation (`massbank.py` + `noise.py`)

Instead of inventing spectra, simms can sample **real reference spectra**
from a MassBank-data checkout (139k+ records), filtered by MS level, ion
mode, instrument type, etc., and perturb them through a parametric noise
model: gaussian ppm mass error, multiplicative intensity noise, peak
dropout (the base peak always survives), and spurious peaks with
exponentially distributed intensities. Each template yields any number of
variants; the record's accession is preserved as `template_accession` —
the ground-truth link.

### 4. Chromatography and the instrument (`realism.py` + `simulate.py`)

A full LC-MS/MS run places each input spectrum on a gradient as an
eluting compound:

- **Elution**: exponentially modified gaussian (EMG) profiles — the
  standard model for tailing chromatographic peaks — with the tail
  constant τ expressed as a multiple of the gaussian σ, and peak width
  growing linearly along the gradient. Computed in log-space to avoid
  overflow.
- **MS1 content**: per compound, one or more **precursor species**
  (electrospray charge envelope z−1/z/z+1 for peptides), each with an
  isotope envelope from the formula or averagine. On top: ever-present
  **background contaminants** (polysiloxane 371.101/445.120, phthalates
  149.023/279.159/391.284, …), per-scan **chemical noise**, and an AR(1)
  log-normal **spray instability** factor shared by all ions in a scan
  (correlated flicker, as electrospray actually behaves).
- **Mass axis**: slow sinusoidal **calibration drift** across the run,
  plus the per-peak noise model.
- **Detector**: **saturation** with a soft knee (exponential approach to
  full scale) and a hard clip.
- **Peak shapes**: optionally **profile mode** — each centroid is
  rendered as a gaussian whose FWHM follows Orbitrap-like resolution
  scaling, R(m/z) = R₂₀₀·√(200/m/z); overlapping peaks sum on a shared
  grid, so unresolved doublets merge exactly as a detector records them.
  Verified: measured FWHM in the written mzML matches m/z ÷ R(m/z) to
  5 decimals.

### 5. Acquisition schemes

- **DDA**: top-N precursor selection per MS1 scan with an intensity
  threshold, per-precursor **dynamic exclusion** in seconds, and
  **co-isolation chimeras** — co-eluting species inside the isolation
  window (and above a relative-intensity floor) mix their fragments into
  the triggered MS2.
- **DIA** (SWATH-style): every MS1 is followed by a complete cycle of
  fixed-width isolation windows over a configurable precursor range, with
  configurable overlap. Each window's MS2 **multiplexes** the fragments
  of every species currently eluting inside it, plus its own chemical
  background; windows fire even when empty, as a real method does.
  Isolation window target and offsets are written into the mzML
  precursor element.

### 6. Output

Runs are written as **indexed mzML** via psims — controlled vocabulary,
instrument/software/data-processing metadata, selected-ion and
isolation-window precursor records — and are read back cleanly by
pyteomics (verified in the test suite), OpenMS, and ProteoWizard-based
tools. Spectral libraries go to mgf/msp/json/pickle via matchms, or to
MassBank record format (which round-trips through simms's own validator).

## Validation strategy

Every physical claim above is asserted by a test, not by inspection:
literature masses, isotope-ratio values, XIC tail skewness for EMG,
FWHM-vs-resolution for profile mode, CE-dependent fragment mass shift,
seed determinism of entire runs, and mzML round-trips through an
independent reader (pyteomics).

## Known limits

- Fragment intensity heuristics are qualitative, not learned; use the
  Prosit backend when fidelity to real HCD spectra matters.
- Elution is EMG-only (no fronting), ion mobility is not modeled, and
  MS2 is modeled for protonated species (no ETD/EThcD).
- Profile mode renders gaussian peaks; real Orbitrap peak shapes have
  mild asymmetry that is not reproduced.
- The mzML writes centroid-derived profile points, not raw transients.
