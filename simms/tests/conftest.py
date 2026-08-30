"""Ensure the real simms package wins over the outer project directory.

When pytest runs from the matchms repository root, ``import simms`` would
otherwise resolve to this project directory itself (a namespace package)
instead of the ``simms/simms`` package.
"""

import pathlib
import sys

_pkg_root = str(pathlib.Path(__file__).resolve().parents[1])
if _pkg_root not in sys.path:
    sys.path.insert(0, _pkg_root)
else:
    sys.path.remove(_pkg_root)
    sys.path.insert(0, _pkg_root)

# drop a half-imported namespace 'simms' so the regular package is found
_mod = sys.modules.get("simms")
if _mod is not None and getattr(_mod, "__file__", None) is None:
    del sys.modules["simms"]
