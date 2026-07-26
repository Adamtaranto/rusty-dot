"""Deterministic synthetic DNA generators for benchmarks.

The real assemblies used during development are private and git-ignored, so the
committed benchmarks feed on reproducible synthetic sequences instead.  All
generators are seeded, so CodSpeed's instruction-count simulation sees identical
inputs on every run.
"""

from __future__ import annotations

import random

_BASES = 'ACGT'
_COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def revcomp(seq: str) -> str:
    """Reverse complement of an uppercase ACGT string."""
    return seq.translate(_COMPLEMENT)[::-1]


def random_dna(length: int, seed: int) -> str:
    """Generate a deterministic random DNA string.

    Parameters
    ----------
    length : int
        Number of bases to generate.
    seed : int
        Seed for the pseudo-random generator (fixed output per seed).

    Returns
    -------
    str
        An uppercase ACGT string of the requested length.
    """
    rng = random.Random(seed)
    return ''.join(rng.choice(_BASES) for _ in range(length))


def homologous_pair(
    length: int,
    seed: int = 0,
    divergence: float = 0.02,
    n_inversions: int = 2,
) -> tuple[str, str]:
    """Build two related sequences sharing homology, with planted inversions.

    The target is a mutated copy of the query (point substitutions at rate
    *divergence*) in which *n_inversions* evenly spaced segments are replaced by
    their reverse complement.  This exercises both forward and reverse-strand
    matching on realistic, collinear-with-inversions input.

    Parameters
    ----------
    length : int
        Length of each sequence in bases.
    seed : int, optional
        Seed controlling the base sequence and mutations.  Default ``0``.
    divergence : float, optional
        Per-base substitution probability applied to the target.  Default
        ``0.02``.
    n_inversions : int, optional
        Number of reverse-complement segments to plant in the target.  Default
        ``2``.

    Returns
    -------
    tuple of str
        ``(query, target)`` — two uppercase ACGT strings of equal length.
    """
    rng = random.Random(seed)
    query = ''.join(rng.choice(_BASES) for _ in range(length))

    # Mutate a copy to produce the target.
    target = list(query)
    for i in range(length):
        if rng.random() < divergence:
            target[i] = rng.choice(_BASES)

    # Plant evenly spaced reverse-complement segments.
    if n_inversions > 0 and length > 0:
        seg = max(1, length // (n_inversions * 4))
        for j in range(n_inversions):
            start = (j * 2 + 1) * (length // (n_inversions * 2 + 1))
            end = min(length, start + seg)
            if end > start:
                target[start:end] = list(revcomp(''.join(target[start:end])))

    return query, ''.join(target)


def multi_contig_group(
    n_contigs: int,
    contig_length: int,
    seed: int,
) -> list[tuple[str, str]]:
    """Generate a list of ``(name, sequence)`` contigs for one group.

    Parameters
    ----------
    n_contigs : int
        Number of contigs to generate.
    contig_length : int
        Length of each contig in bases.
    seed : int
        Base seed; contig ``i`` uses seed ``seed + i``.

    Returns
    -------
    list of tuple of str
        ``[(name, sequence), ...]`` with names ``contig_0``, ``contig_1``, ...
    """
    return [
        (f'contig_{i}', random_dna(contig_length, seed + i)) for i in range(n_contigs)
    ]
