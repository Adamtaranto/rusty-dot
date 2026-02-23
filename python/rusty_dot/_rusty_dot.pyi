"""
Type stubs for the rusty_dot Rust extension module.

All functions in this module are implemented in Rust via PyO3.
"""

from __future__ import annotations

class SequenceIndex:
    """FM-index backed sequence comparison engine.

    Builds FM-indexes and k-mer sets for DNA sequences read from FASTA
    files or provided directly, then supports fast pairwise k-mer
    coordinate lookup with optional sequential run merging.

    Parameters
    ----------
    k : int
        The k-mer length to use for indexing and comparison.

    Examples
    --------
    >>> idx = SequenceIndex(k=10)
    >>> idx.add_sequence("seq1", "ACGTACGTACGT")
    >>> idx.add_sequence("seq2", "TACGTACGTACG")
    >>> matches = idx.compare_sequences("seq1", "seq2")
    """

    def __init__(self, k: int) -> None: ...

    def add_sequence(self, name: str, seq: str) -> None:
        """Add a single sequence to the index."""
        ...

    def load_fasta(self, path: str) -> list[str]:
        """Load all sequences from a FASTA or gzipped FASTA file."""
        ...

    def sequence_names(self) -> list[str]:
        """Get the list of sequence names in the index."""
        ...

    def get_kmer_set(self, name: str) -> set[str]:
        """Get the k-mer set for a named sequence."""
        ...

    def get_sequence_length(self, name: str) -> int:
        """Get the length of a named sequence."""
        ...

    def compare_sequences(
        self,
        query_name: str,
        target_name: str,
        merge: bool = True,
    ) -> list[tuple[int, int, int, int]]:
        """Find shared k-mer matches between two sequences."""
        ...

    def get_paf(
        self,
        query_name: str,
        target_name: str,
        merge: bool = True,
    ) -> list[str]:
        """Get PAF-formatted alignment strings for a sequence pair."""
        ...

    def precompute_all_pairs(self, merge: bool = True) -> list[tuple[str, str]]:
        """Pre-calculate pairwise comparisons for all sequence pairs."""
        ...

    def save(self, path: str) -> None:
        """Save the index to a binary file."""
        ...

    def load(self, path: str) -> None:
        """Load sequences from a saved index file."""
        ...

    @property
    def k(self) -> int:
        """The k-mer length used for indexing."""
        ...

    def __len__(self) -> int: ...
    def __repr__(self) -> str: ...


def py_read_fasta(path: str) -> dict[str, str]: ...
def py_build_kmer_set(seq: str, k: int) -> set[str]: ...
def py_find_kmer_coords(seq: str, kmers: list[str]) -> dict[str, list[int]]: ...
def py_merge_kmer_runs(
    kmer_coords: dict[str, list[int]],
    query_kmer_positions: dict[str, list[int]],
    k: int,
) -> list[tuple[int, int, int, int]]: ...
def py_coords_to_paf(
    matches: list[tuple[int, int, int, int]],
    query_name: str,
    query_len: int,
    target_name: str,
    target_len: int,
) -> list[str]: ...
def py_save_index(path: str, sequences: dict[str, str], k: int) -> None: ...
def py_load_index(path: str) -> tuple[dict[str, list[str]], int]: ...
