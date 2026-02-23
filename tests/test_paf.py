"""Tests for PAF format output."""

from rusty_dot._rusty_dot import py_coords_to_paf


def test_paf_basic():
    """Test basic PAF output format."""
    matches = [(0, 10, 5, 15)]
    lines = py_coords_to_paf(matches, 'queryA', 100, 'targetB', 200)
    assert len(lines) == 1
    fields = lines[0].split('\t')
    assert fields[0] == 'queryA'
    assert fields[1] == '100'
    assert fields[2] == '0'
    assert fields[3] == '10'
    assert fields[4] == '+'
    assert fields[5] == 'targetB'
    assert fields[6] == '200'
    assert fields[7] == '5'
    assert fields[8] == '15'


def test_paf_empty_matches():
    """Test that empty match list returns empty PAF."""
    lines = py_coords_to_paf([], 'q', 100, 't', 200)
    assert lines == []


def test_paf_multiple_matches():
    """Test multiple matches produce multiple PAF lines."""
    matches = [(0, 10, 0, 10), (20, 30, 25, 35)]
    lines = py_coords_to_paf(matches, 'q', 50, 't', 60)
    assert len(lines) == 2


def test_paf_match_length():
    """Test that match length field is correct."""
    matches = [(0, 15, 0, 15)]
    lines = py_coords_to_paf(matches, 'q', 50, 't', 50)
    fields = lines[0].split('\t')
    # Field 9 = number of residue matches = 15
    assert fields[9] == '15'
    # Field 10 = alignment block length = 15
    assert fields[10] == '15'
