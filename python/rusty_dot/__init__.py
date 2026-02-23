"""
rusty-dot: Fast dot plot comparisons of DNA sequences using an FM-Index.

This package provides:
- Rust-backed FM-index for fast k-mer lookup (via PyO3 bindings)
- Python classes for building and querying sequence indexes
- Dotplot visualization tools
- PAF file I/O with CIGAR string parsing
- Gravity-based contig reordering for maximising dotplot collinearity

Examples
--------
Basic usage:

>>> from rusty_dot import SequenceIndex
>>> idx = SequenceIndex(k=10)
>>> idx.add_sequence("seq1", "ACGTACGTACGT")
>>> idx.add_sequence("seq2", "TACGTACGTACG")
>>> matches = idx.compare_sequences_stranded("seq1", "seq2")
"""

from rusty_dot._rusty_dot import (  # noqa: F401
    SequenceIndex,
    py_build_kmer_set,
    py_coords_to_paf,
    py_find_kmer_coords,
    py_load_index,
    py_merge_kmer_runs,
    py_read_fasta,
    py_save_index,
)
from rusty_dot.dotplot import DotPlotter  # noqa: F401
from rusty_dot.paf_io import (  # noqa: F401
    PafAlignment,
    PafRecord,
    compute_gravity_contigs,
    parse_paf_file,
)

__version__ = '0.1.0'
__all__ = [
    'SequenceIndex',
    'DotPlotter',
    'PafRecord',
    'PafAlignment',
    'parse_paf_file',
    'compute_gravity_contigs',
    'py_read_fasta',
    'py_build_kmer_set',
    'py_find_kmer_coords',
    'py_merge_kmer_runs',
    'py_coords_to_paf',
    'py_save_index',
    'py_load_index',
]
