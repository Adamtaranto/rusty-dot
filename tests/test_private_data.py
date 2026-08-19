"""Opt-in benchmarks against private local assemblies (never committed).

These tests exercise the real-data hot paths with the git-ignored fungal
assemblies under ``.plans/data/asm``.  They are skipped unless the files
exist AND ``RUSTY_DOT_PRIVATE_DATA=1`` is set, so ordinary test runs (local
and CI) never touch them:

    RUSTY_DOT_PRIVATE_DATA=1 pytest tests/test_private_data.py -m private_data -s

Timings are printed with ``-s`` for before/after comparisons in PR
descriptions; nothing derived from the private data is written to the repo.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
import time

import pytest

_REPO = Path(__file__).resolve().parents[1]
_ASM_DIR = _REPO / '.plans' / 'data' / 'asm'
_ASSEMBLIES = sorted(_ASM_DIR.glob('*.fasta')) if _ASM_DIR.is_dir() else []

sys.path.insert(0, str(_REPO / 'app'))

pytestmark = [
    pytest.mark.private_data,
    pytest.mark.skipif(
        os.environ.get('RUSTY_DOT_PRIVATE_DATA') != '1' or len(_ASSEMBLIES) < 2,
        reason='set RUSTY_DOT_PRIVATE_DATA=1 with assemblies in .plans/data/asm',
    ),
]


def _timed(label: str, fn, *args, **kwargs):
    """Run *fn*, print its wall-clock time, and return its result."""
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    print(f'{label}: {time.perf_counter() - t0:.2f} s')
    return result


def test_parse_real_assemblies():
    """Byte-parse both real assemblies and report throughput."""
    from core.fasta import parse_fasta_bytes

    for path in _ASSEMBLIES[:2]:
        data = path.read_bytes()
        fi = _timed(
            f'parse {path.name} ({len(data) / 1e6:.0f} MB)', parse_fasta_bytes, data
        )
        assert fi.records
        print(f'  {len(fi.records)} contigs, {fi.total_length / 1e6:.1f} Mb')


def test_kmer_pipeline_real_assemblies():
    """Full k-mer pipeline on real data: index build, matches, reorder."""
    from core.fasta import parse_fasta_bytes

    from rusty_dot.paf_io import CrossIndex

    query = parse_fasta_bytes(_ASSEMBLIES[0].read_bytes())
    target = parse_fasta_bytes(_ASSEMBLIES[1].read_bytes())

    cross = CrossIndex(k=21)

    def _build():
        for name, seq in query.records:
            cross.add_sequence(name, seq, group='query')
        for name, seq in target.records:
            cross.add_sequence(name, seq, group='target')

    _timed('index build', _build)
    _timed('compute_matches', cross.compute_matches, 'query', 'target', True)
    paf_lines = cross.get_paf([('query', 'target')])
    assert paf_lines
    print(f'  {len(paf_lines)} PAF records')
    _timed(
        'reorder_for_colinearity',
        cross.reorder_for_colinearity,
        'query',
        'target',
    )
