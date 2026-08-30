"""simms: simulated mass spectrometry test data toolkit.

A command-line orchestrator that combines:

- pure-Python simulation (theoretical peptide fragment spectra, isotope
  patterns, configurable noise models),
- the MassBank-data repository as a source of real reference spectra,
- matchms for spectral I/O, cleaning and merging,
- psims/pyteomics for writing and reading mzML,
- OpenMS TOPP tools as an optional native backend when installed.
"""

__version__ = "0.2.0"
