import pytest

from simms.chem import (adduct_mz, isotope_envelope_mz, isotope_pattern,
                        monoisotopic_mass, parse_formula)


def test_parse_formula():
    assert parse_formula("C9H8O4") == {"C": 9, "H": 8, "O": 4}
    assert parse_formula("CH4") == {"C": 1, "H": 4}
    assert parse_formula("C2H4Br2") == {"C": 2, "H": 4, "Br": 2}


def test_parse_formula_rejects_garbage():
    with pytest.raises(ValueError):
        parse_formula("C9H8O4X")
    with pytest.raises(ValueError):
        parse_formula("")


def test_monoisotopic_mass_aspirin():
    assert monoisotopic_mass("C9H8O4") == pytest.approx(180.04226, abs=1e-4)


def test_adduct_mz_protonated_aspirin():
    assert adduct_mz("C9H8O4", "[M+H]+") == pytest.approx(181.04954, abs=1e-4)


def test_adduct_mz_deprotonated():
    assert adduct_mz("C9H8O4", "[M-H]-") == pytest.approx(179.03498, abs=1e-4)


def test_isotope_pattern_bromine_triplet():
    # Br2 gives the characteristic ~1:2:1 M/M+2/M+4 triplet
    peaks = isotope_pattern("C2H4Br2")
    intense = sorted(peaks, key=lambda p: p[1], reverse=True)[:3]
    intense.sort(key=lambda p: p[0])
    assert intense[1][1] == pytest.approx(1.0)
    assert intense[0][1] == pytest.approx(0.51, abs=0.03)
    assert intense[2][1] == pytest.approx(0.49, abs=0.03)
    assert intense[1][0] - intense[0][0] == pytest.approx(2.0, abs=0.01)


def test_isotope_envelope_mz_spacing():
    envelope = isotope_envelope_mz("C50H100", "[M+H]+")
    assert envelope[1][0] - envelope[0][0] == pytest.approx(1.00336, abs=0.001)
    envelope_2plus = isotope_envelope_mz("C50H100", "[M+2H]2+")
    assert envelope_2plus[1][0] - envelope_2plus[0][0] == pytest.approx(0.5017, abs=0.001)
