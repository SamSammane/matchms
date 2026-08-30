import textwrap

import pytest

from simms.massbank import (parse_record, record_to_spectrum,
                            retention_time_seconds, spectrum_to_record_text,
                            validate_record_text)

SAMPLE_RECORD = textwrap.dedent("""\
    ACCESSION: MSBNK-TEST-XX000001
    RECORD_TITLE: Testol; LC-ESI-QTOF; MS2
    DATE: 2024.01.01
    AUTHORS: Test Author
    LICENSE: CC BY
    CH$NAME: Testol
    CH$FORMULA: C9H8O4
    CH$EXACT_MASS: 180.04226
    CH$SMILES: CC(=O)OC1=CC=CC=C1C(=O)O
    CH$LINK: INCHIKEY BSYNRYMUTXBXSQ-UHFFFAOYSA-N
    AC$INSTRUMENT: Test QTOF
    AC$INSTRUMENT_TYPE: LC-ESI-QTOF
    AC$MASS_SPECTROMETRY: MS_TYPE MS2
    AC$MASS_SPECTROMETRY: ION_MODE POSITIVE
    AC$CHROMATOGRAPHY: RETENTION_TIME 5.5 min
    MS$FOCUSED_ION: PRECURSOR_M/Z 181.0495
    MS$FOCUSED_ION: PRECURSOR_TYPE [M+H]+
    PK$NUM_PEAK: 3
    PK$PEAK: m/z int. rel.int.
      95.0491 10000 250
      135.0441 40000 999
      163.0390 20000 500
    //
    """)


def test_parse_record():
    record = parse_record(SAMPLE_RECORD)
    assert record.metadata["accession"] == "MSBNK-TEST-XX000001"
    assert record.metadata["formula"] == "C9H8O4"
    assert record.metadata["ms_type"] == "MS2"
    assert record.metadata["ionmode"] == "positive"
    assert record.metadata["inchikey"] == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert record.metadata["adduct"] == "[M+H]+"
    assert len(record.peaks) == 3
    assert record.peaks[1] == (135.0441, 40000.0, 999)


def test_retention_time_parsing():
    assert retention_time_seconds("5.5 min") == pytest.approx(330.0)
    assert retention_time_seconds("30 sec") == pytest.approx(30.0)
    assert retention_time_seconds("11.9") == pytest.approx(714.0)
    assert retention_time_seconds(None) is None
    assert retention_time_seconds("N/A") is None


def test_record_to_spectrum():
    spectrum = record_to_spectrum(parse_record(SAMPLE_RECORD))
    assert spectrum is not None
    assert spectrum.metadata["precursor_mz"] == pytest.approx(181.0495)
    assert spectrum.metadata["ms_level"] == 2
    assert spectrum.metadata["retention_time"] == pytest.approx(330.0)
    assert list(spectrum.peaks.mz) == sorted(spectrum.peaks.mz)


def test_validate_ok_record():
    assert validate_record_text(SAMPLE_RECORD) == []


def test_validate_catches_problems():
    broken = SAMPLE_RECORD.replace("PK$NUM_PEAK: 3", "PK$NUM_PEAK: 5")
    issues = validate_record_text(broken)
    assert any("PK$NUM_PEAK" in issue for issue in issues)
    truncated = SAMPLE_RECORD.replace("//", "").replace("ACCESSION: ", "XACCESSION: ")
    issues = validate_record_text(truncated)
    assert any("ACCESSION" in issue for issue in issues)
    assert any("//" in issue for issue in issues)


def test_roundtrip_spectrum_to_record():
    spectrum = record_to_spectrum(parse_record(SAMPLE_RECORD))
    text = spectrum_to_record_text(spectrum, "MSBNK-SIMMS-SIM000001")
    assert validate_record_text(text) == []
    reparsed = record_to_spectrum(parse_record(text))
    assert len(reparsed.peaks) == 3
    assert reparsed.metadata["formula"] == "C9H8O4"
