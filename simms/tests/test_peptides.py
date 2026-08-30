import pytest

from simms.peptides import (digest, fragment_ions, peptide_mass, precursor_mz,
                            read_fasta)


def test_peptide_mass_known_value():
    # PEPTIDE monoisotopic neutral mass, standard reference value
    assert peptide_mass("PEPTIDE") == pytest.approx(799.35997, abs=1e-3)


def test_precursor_mz_elvislivesk():
    assert precursor_mz("ELVISLIVESK", 2) == pytest.approx(615.3712, abs=1e-3)


def test_fragment_ions_by_series():
    ions = fragment_ions("PEPTIDE", ion_types="by")
    by_type = {ion.ion_type: ion.mz for ion in ions}
    # canonical values for PEPTIDE
    assert by_type["b2+"] == pytest.approx(227.1026, abs=1e-3)
    assert by_type["y1+"] == pytest.approx(148.0604, abs=1e-3)
    assert by_type["y6+"] == pytest.approx(703.3145, abs=1e-3)
    assert len([i for i in ions if i.ion_type.startswith("b")]) == 6
    assert len([i for i in ions if i.ion_type.startswith("y")]) == 6


def test_fragment_ions_sorted_by_mz():
    ions = fragment_ions("ELVISLIVESK")
    mzs = [ion.mz for ion in ions]
    assert mzs == sorted(mzs)


def test_trypsin_digest():
    peptides = digest("MKWVTFISLLRPEPTIDEK", enzyme="trypsin",
                      missed_cleavages=0, min_length=1, max_length=50)
    assert "MK" in peptides
    # no cleavage between R and P (proline rule)
    assert "WVTFISLLRPEPTIDEK" in peptides


def test_trypsin_missed_cleavages():
    no_missed = digest("AAAKBBBRCCCK".replace("B", "G"), missed_cleavages=0,
                       min_length=1, max_length=50)
    one_missed = digest("AAAKBBBRCCCK".replace("B", "G"), missed_cleavages=1,
                        min_length=1, max_length=50)
    assert len(one_missed) > len(no_missed)


def test_read_fasta(tmp_path):
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">prot1 description\nPEPTIDEK\nAAAR\n>prot2\nGGGK\n")
    entries = list(read_fasta(str(fasta)))
    assert entries == [("prot1 description", "PEPTIDEKAAAR"), ("prot2", "GGGK")]
