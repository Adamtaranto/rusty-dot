"""Golden parity tests for the rolling-hash k-mer matching engine.

These tests validate the Rust ``SequenceIndex`` matching output against an
independent, brute-force Python reference that enumerates shared k-mers by exact
string comparison.  The reference is deliberately naive (O(n·m)) so it shares no
logic with the ntHash-based engine — if the two agree across forward, reverse
complement, palindrome, inversion, repeat, ``N``-containing, and randomised
inputs, the engine reproduces exact k-mer search.
"""

import random

import pytest

from rusty_dot._rusty_dot import SequenceIndex

_COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def _revcomp(seq: str) -> str:
    """Reverse complement of an uppercase ACGT string."""
    return seq.translate(_COMPLEMENT)[::-1]


def _kmer_positions(seq: str, k: int) -> dict[str, list[int]]:
    """Map each ACGT-only k-mer to its sorted 0-based start positions."""
    positions: dict[str, list[int]] = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if all(b in 'ACGT' for b in kmer):
            positions.setdefault(kmer, []).append(i)
    return positions


def _naive_forward(query: str, target: str, k: int) -> set[tuple[int, int, int, int]]:
    """Brute-force set of forward (query_start, query_end, target_start, target_end)."""
    q_pos = _kmer_positions(query, k)
    t_pos = _kmer_positions(target, k)
    hits = set()
    for kmer, qs in q_pos.items():
        ts = t_pos.get(kmer)
        if ts is None:
            continue
        for qp in qs:
            for tp in ts:
                hits.add((qp, qp + k, tp, tp + k))
    return hits


def _naive_stranded(
    query: str, target: str, k: int
) -> set[tuple[int, int, int, int, str]]:
    """Brute-force set of both-strand hits with strand labels."""
    q_pos = _kmer_positions(query, k)
    t_pos = _kmer_positions(target, k)
    hits: set[tuple[int, int, int, int, str]] = set()
    # Forward
    for kmer, qs in q_pos.items():
        ts = t_pos.get(kmer)
        if ts is not None:
            for qp in qs:
                for tp in ts:
                    hits.add((qp, qp + k, tp, tp + k, '+'))
    # Reverse complement: query k-mer whose RC occurs in the target
    for kmer, qs in q_pos.items():
        ts = t_pos.get(_revcomp(kmer))
        if ts is not None:
            for qp in qs:
                for tp in ts:
                    hits.add((qp, qp + k, tp, tp + k, '-'))
    return hits


def _engine_forward(query: str, target: str, k: int):
    idx = SequenceIndex(k=k)
    idx.add_sequence('q', query)
    idx.add_sequence('t', target)
    return set(tuple(m) for m in idx.compare_sequences('q', 't', False))


def _engine_stranded(query: str, target: str, k: int):
    idx = SequenceIndex(k=k)
    idx.add_sequence('q', query)
    idx.add_sequence('t', target)
    return set(tuple(m) for m in idx.compare_sequences_stranded('q', 't', False))


# --- Targeted fixtures -----------------------------------------------------

_FIXTURES = [
    # (name, query, target, k)
    ('identical', 'ACGTACGTACGT', 'ACGTACGTACGT', 4),
    ('offset', 'ACGTACGTACGT', 'TACGTACGTACG', 4),
    ('pure_rc', 'AAACAAACAAAC', 'GTTTGTTTGTTT', 4),
    ('palindrome', 'ACGTACGT', 'ACGTACGT', 4),
    ('inversion', 'AAAACCCCGGGG', 'CCCCGGGGTTTT', 5),
    ('tandem_repeat', 'ATATATATATAT', 'ATATATATATAT', 3),
    ('with_n', 'ACGTNACGTACG', 'ACGTACGTNACG', 4),
    ('no_shared', 'AAAAAAAA', 'CCCCCCCC', 4),
    ('k_equals_len', 'ACGTACGT', 'ACGTACGT', 8),
]


@pytest.mark.parametrize('name,query,target,k', _FIXTURES)
def test_forward_matches_naive_reference(name, query, target, k):
    """Forward k-mer hits match the brute-force reference on targeted fixtures."""
    assert _engine_forward(query, target, k) == _naive_forward(query, target, k), name


@pytest.mark.parametrize('name,query,target,k', _FIXTURES)
def test_stranded_matches_naive_reference(name, query, target, k):
    """Both-strand k-mer hits match the brute-force reference on targeted fixtures."""
    assert _engine_stranded(query, target, k) == _naive_stranded(query, target, k), name


# --- Randomised fuzzing -----------------------------------------------------


@pytest.mark.parametrize('seed', range(25))
def test_random_sequences_match_naive_reference(seed):
    """Randomised query/target pairs agree with the reference for both strands."""
    rng = random.Random(seed)
    k = rng.randint(3, 12)
    query = ''.join(rng.choice('ACGT') for _ in range(rng.randint(k, 80)))
    # Build a target that shares real homology: embed part of the query and its
    # reverse complement, plus random padding, so both strands are exercised.
    piece = query[: rng.randint(1, max(1, len(query) // 2))]
    parts = [
        ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, 10))),
        piece,
        ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, 10))),
        _revcomp(piece),
        ''.join(rng.choice('ACGT') for _ in range(rng.randint(0, 10))),
    ]
    target = ''.join(parts)

    assert _engine_forward(query, target, k) == _naive_forward(query, target, k)
    assert _engine_stranded(query, target, k) == _naive_stranded(query, target, k)
