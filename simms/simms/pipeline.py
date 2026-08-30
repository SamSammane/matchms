"""YAML pipeline runner: chain any simms commands declaratively.

Pipeline file format::

    name: my-simulated-dataset
    steps:
      - generate from-massbank --massbank /data/massbank --n 20 --out step1.mgf
      - run: generate peptides
        args: {sequences: "PEPTIDE,ELVISLIVESK", out: step2.mgf}
      - ["merge", "-i", "step1.mgf", "step2.mgf", "-o", "combined.mgf"]

Each step is one simms command line, given as a string (shlex-split), a
token list, or a mapping with ``run`` (subcommand words) and ``args``
(mapping rendered to --key value options; true renders as a bare flag).
"""

from __future__ import annotations

import shlex
from typing import Any, Dict, List

import yaml


def _step_to_argv(step: Any) -> List[str]:
    if isinstance(step, str):
        return shlex.split(step)
    if isinstance(step, list):
        return [str(token) for token in step]
    if isinstance(step, dict):
        run = step.get("run")
        if not run:
            raise ValueError(f"pipeline step mapping needs a 'run' key: {step!r}")
        argv = shlex.split(run) if isinstance(run, str) else [str(t) for t in run]
        for key, value in (step.get("args") or {}).items():
            option = f"--{str(key).replace('_', '-')}"
            if value is True:
                argv.append(option)
            elif value is False or value is None:
                continue
            elif isinstance(value, list):
                argv.append(option)
                argv.extend(str(item) for item in value)
            else:
                argv.extend([option, str(value)])
        return argv
    raise ValueError(f"unsupported pipeline step type: {type(step).__name__}")


def run_pipeline(path: str, dry_run: bool = False,
                 quiet_steps: bool = False) -> Dict[str, Any]:
    import contextlib
    import sys

    from . import cli  # deferred to avoid a circular import

    with open(path, "r", encoding="utf-8") as handle:
        spec = yaml.safe_load(handle)
    if not isinstance(spec, dict) or "steps" not in spec:
        raise ValueError("pipeline file must be a mapping with a 'steps' list")
    results = []
    for index, step in enumerate(spec["steps"], start=1):
        argv = _step_to_argv(step)
        if argv and argv[0] == "simms":
            argv = argv[1:]
        if dry_run:
            results.append({"step": index, "argv": argv, "status": "dry-run"})
            continue
        if quiet_steps:
            # keep stdout clean for the final JSON result; step logs go to stderr
            with contextlib.redirect_stdout(sys.stderr):
                code = cli.main(argv)
        else:
            code = cli.main(argv)
        results.append({"step": index, "argv": argv,
                        "status": "ok" if code == 0 else f"exit {code}"})
        if code != 0:
            break
    return {"pipeline": spec.get("name", path), "steps_run": len(results),
            "results": results,
            "ok": all(r["status"] in ("ok", "dry-run") for r in results)}
