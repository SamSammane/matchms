"""End-to-end CLI tests exercising generation, merging, conversion,
pipelines and mzML round-tripping through the public entry point."""

import json
import os

import pytest

from simms.cli import main
from simms.io_utils import load_any

MASSBANK_REPO = "/home/user/massbank-data"
needs_massbank = pytest.mark.skipif(not os.path.isdir(MASSBANK_REPO),
                                    reason="massbank-data checkout not present")


def run(capsys, *argv):
    code = main(list(argv))
    out = capsys.readouterr().out
    return code, out


def test_generate_peptides_to_mgf(tmp_path, capsys):
    out = tmp_path / "peps.mgf"
    code, stdout = run(capsys, "generate", "peptides", "--sequences", "PEPTIDE,ELVISLIVESK",
                       "--charges", "2,3", "--out", str(out), "--json")
    assert code == 0
    result = json.loads(stdout)
    assert result["spectra_written"] == 4
    spectra = load_any(str(out))
    assert len(spectra) == 4
    assert spectra[0].metadata["ms_level"] == 2


def test_generate_peptides_from_fasta(tmp_path, capsys):
    fasta = tmp_path / "prot.fasta"
    fasta.write_text(">p1\nMKWVTFISLLFLFSSAYSRGVFRR\n")
    out = tmp_path / "digest.msp"
    code, stdout = run(capsys, "generate", "peptides", "--fasta", str(fasta),
                       "--missed-cleavages", "1", "--out", str(out), "--json")
    assert code == 0
    assert json.loads(stdout)["spectra_written"] > 0


def test_generate_isotopes_with_variants(tmp_path, capsys):
    out = tmp_path / "iso.mgf"
    code, stdout = run(capsys, "generate", "isotopes", "--formula", "C9H8O4",
                       "--variants", "3", "--noise-preset", "default",
                       "--seed", "5", "--out", str(out), "--json")
    assert code == 0
    assert json.loads(stdout)["spectra_written"] == 3


def test_merge_and_convert_roundtrip(tmp_path, capsys):
    a = tmp_path / "a.mgf"
    b = tmp_path / "b.msp"
    run(capsys, "generate", "peptides", "--sequences", "PEPTIDE", "--out", str(a))
    run(capsys, "generate", "isotopes", "--formula", "C6H12O6", "--out", str(b))
    merged = tmp_path / "merged.mgf"
    code, stdout = run(capsys, "merge", "-i", str(a), str(b), "-o", str(merged), "--json")
    assert code == 0
    assert json.loads(stdout)["written"] == 2
    converted = tmp_path / "merged.json"
    code, _ = run(capsys, "convert", "-i", str(merged), "-o", str(converted))
    assert code == 0
    assert len(load_any(str(converted))) == 2


def test_massbank_record_export_and_validate(tmp_path, capsys):
    mgf = tmp_path / "x.mgf"
    run(capsys, "generate", "isotopes", "--formula", "C9H8O4", "--out", str(mgf))
    records_dir = tmp_path / "records"
    code, _ = run(capsys, "convert", "-i", str(mgf), "-o", str(records_dir))
    assert code == 0
    record_files = sorted(records_dir.glob("MSBNK-SIMMS-*.txt"))
    assert len(record_files) == 1
    code, stdout = run(capsys, "massbank", "validate", str(record_files[0]), "--json")
    assert code == 0
    assert json.loads(stdout)["failed"] == 0


def test_lcms_run_mzml_roundtrip(tmp_path, capsys):
    mgf = tmp_path / "compounds.mgf"
    run(capsys, "generate", "peptides", "--sequences", "PEPTIDE,ELVISLIVESK",
        "--out", str(mgf))
    mzml_path = tmp_path / "run.mzML"
    code, stdout = run(capsys, "generate", "lcms-run", "-i", str(mgf),
                       "--gradient", "60", "--ms1-interval", "5",
                       "--out", str(mzml_path), "--seed", "1", "--json")
    assert code == 0
    result = json.loads(stdout)
    assert result["ms1_scans"] == 12
    assert result["ms2_scans"] > 0

    from pyteomics import mzml as pymzml
    scans = list(pymzml.read(str(mzml_path)))
    assert len(scans) == result["total_scans"]
    ms2 = [s for s in scans if s["ms level"] == 2]
    precursor = ms2[0]["precursorList"]["precursor"][0]
    assert "selectedIonList" in precursor


def test_lcms_run_deterministic(tmp_path, capsys):
    mgf = tmp_path / "c.mgf"
    run(capsys, "generate", "peptides", "--sequences", "PEPTIDE", "--out", str(mgf))
    out1, out2 = tmp_path / "r1.mzML", tmp_path / "r2.mzML"
    _, s1 = run(capsys, "generate", "lcms-run", "-i", str(mgf), "--gradient", "30",
                "--out", str(out1), "--seed", "9", "--json")
    _, s2 = run(capsys, "generate", "lcms-run", "-i", str(mgf), "--gradient", "30",
                "--out", str(out2), "--seed", "9", "--json")
    r1, r2 = json.loads(s1), json.loads(s2)
    assert r1["total_scans"] == r2["total_scans"]


def test_pipeline_run(tmp_path, capsys):
    mgf1 = tmp_path / "s1.mgf"
    merged = tmp_path / "merged.msp"
    pipeline_file = tmp_path / "p.yaml"
    pipeline_file.write_text(f"""
name: test-pipeline
steps:
  - generate peptides --sequences PEPTIDE --out {mgf1}
  - run: merge
    args: {{inputs: ["{mgf1}"], out: "{merged}"}}
""")
    code, stdout = run(capsys, "pipeline", "run", str(pipeline_file), "--json")
    assert code == 0
    result = json.loads(stdout)
    assert result["ok"] is True
    assert result["steps_run"] == 2
    assert merged.exists()


def test_pipeline_dry_run(tmp_path, capsys):
    pipeline_file = tmp_path / "p.yaml"
    pipeline_file.write_text("steps:\n  - generate peptides --sequences X --out nope.mgf\n")
    code, stdout = run(capsys, "pipeline", "run", str(pipeline_file), "--dry-run", "--json")
    assert code == 0
    assert json.loads(stdout)["results"][0]["status"] == "dry-run"


def test_openms_list_reports_gracefully(capsys):
    code, stdout = run(capsys, "openms", "list", "--json")
    assert code == 0
    result = json.loads(stdout)
    assert "available" in result and "known_tools" in result


def test_describe_manifest(capsys):
    code, stdout = run(capsys, "describe")
    assert code == 0
    manifest = json.loads(stdout)
    assert "generate lcms-run" in manifest["commands"]


@needs_massbank
def test_generate_from_real_massbank(tmp_path, capsys):
    out = tmp_path / "sim.mgf"
    code, stdout = run(capsys, "generate", "from-massbank", "--massbank", MASSBANK_REPO,
                       "--n", "3", "--variants", "2", "--noise-preset", "default",
                       "--seed", "11", "--out", str(out), "--json")
    assert code == 0
    result = json.loads(stdout)
    assert result["spectra_written"] == 6
    spectra = load_any(str(out))
    assert all(s.metadata.get("simulated") for s in spectra)
    assert all(s.metadata.get("template_accession", "").startswith("MSBNK-")
               for s in spectra)


@needs_massbank
def test_massbank_stats(capsys):
    code, stdout = run(capsys, "massbank", "stats", "--massbank", MASSBANK_REPO, "--json")
    assert code == 0
    result = json.loads(stdout)
    assert result["total_records"] > 100000
