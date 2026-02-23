"""Type stubs for the rusty_dot Rust extension module.

All functions and classes in this module are implemented in Rust via PyO3.
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

    def __init__(self, k: int) -> None:
        """Initialise an empty SequenceIndex.

        Parameters
        ----------
        k : int
            The k-mer length to use for indexing and comparison.
            Must be greater than zero.

        Raises
        ------
        ValueError
            If ``k`` is zero.
        """
        ...

    def add_sequence(self, name: str, seq: str) -> None:
        """Add a single sequence to the index.

        Builds an FM-index and collects the k-mer set for the given
        sequence, storing both for subsequent comparisons.

        Parameters
        ----------
        name : str
            Unique identifier for the sequence.
        seq : str
            DNA sequence string. Uppercase is recommended; lowercase
            input is accepted and treated as uppercase.

        Raises
        ------
        ValueError
            If the FM-index cannot be built (e.g., invalid characters).
        """
        ...

    def load_fasta(self, path: str) -> list[str]:
        """Load all sequences from a FASTA or gzipped FASTA file.

        Parses the file with needletail (automatic gzip detection) and
        calls ``add_sequence`` for every record found.

        Parameters
        ----------
        path : str
            Path to a FASTA (``.fa``, ``.fasta``) or gzipped FASTA
            (``.fa.gz``, ``.fasta.gz``) file.

        Returns
        -------
        list[str]
            Sequence names (record identifiers) that were added, in the
            order they appear in the file.

        Raises
        ------
        ValueError
            If the file cannot be opened or parsed.
        """
        ...

    def sequence_names(self) -> list[str]:
        """Return all sequence names currently held in the index.

        Returns
        -------
        list[str]
            Unordered list of sequence identifiers.
        """
        ...

    def get_kmer_set(self, name: str) -> set[str]:
        """Return the set of unique k-mers for a named sequence.

        Only k-mers composed entirely of the characters A, C, G, T are
        included (k-mers containing N or other IUPAC codes are excluded).

        Parameters
        ----------
        name : str
            The sequence identifier.

        Returns
        -------
        set[str]
            The set of unique k-mer strings (length ``k``) found in the
            sequence.

        Raises
        ------
        KeyError
            If ``name`` is not present in the index.
        """
        ...

    def get_sequence_length(self, name: str) -> int:
        """Return the length of a named sequence in base pairs.

        Parameters
        ----------
        name : str
            The sequence identifier.

        Returns
        -------
        int
            The sequence length in base pairs (not including any
            internal sentinel character).

        Raises
        ------
        KeyError
            If ``name`` is not present in the index.
        """
        ...

    def compare_sequences(
        self,
        query_name: str,
        target_name: str,
        merge: bool = True,
    ) -> list[tuple[int, int, int, int]]:
        """Find shared k-mer matches between two sequences.

        Intersects the k-mer sets of the two sequences and looks up the
        coordinates of each shared k-mer in both FM-indexes.  Uses the
        smaller k-mer set as the probe for efficiency.

        When ``merge=True`` (default) consecutive k-mer hits on the same
        co-linear diagonal are merged into a single coordinate block.
        Results for each ``(query_name, target_name, merge)`` combination
        are cached; repeated calls with the same arguments are free.

        Parameters
        ----------
        query_name : str
            Name of the query sequence (defines the y-axis in a dotplot).
        target_name : str
            Name of the target sequence (defines the x-axis in a dotplot).
        merge : bool, optional
            When ``True`` (default), merge consecutive co-linear k-mer
            hits into contiguous coordinate blocks.  When ``False``,
            every individual k-mer hit is returned as its own block.

        Returns
        -------
        list[tuple[int, int, int, int]]
            List of ``(query_start, query_end, target_start, target_end)``
            coordinate tuples.  All coordinates are 0-based; end
            positions are exclusive.

        Raises
        ------
        KeyError
            If either ``query_name`` or ``target_name`` is not present in
            the index.
        """
        ...

    def get_paf(
        self,
        query_name: str,
        target_name: str,
        merge: bool = True,
    ) -> list[str]:
        """Return PAF-formatted alignment strings for a sequence pair.

        Calls ``compare_sequences`` internally (using the cache if
        available) and formats the result as 12-column PAF lines.

        Parameters
        ----------
        query_name : str
            Name of the query sequence.
        target_name : str
            Name of the target sequence.
        merge : bool, optional
            Whether to merge consecutive co-linear k-mer runs before
            generating PAF lines.  Default is ``True``.

        Returns
        -------
        list[str]
            PAF lines as tab-separated strings.  Each line has 12 fields:
            query name, query length, query start, query end, strand,
            target name, target length, target start, target end,
            residue matches, alignment block length, mapping quality.

        Raises
        ------
        KeyError
            If either sequence name is not present in the index.
        """
        ...

    def precompute_all_pairs(self, merge: bool = True) -> list[tuple[str, str]]:
        """Pre-calculate comparisons for every ordered sequence pair.

        Iterates over all ``(i, j)`` pairs where ``i != j`` and calls
        ``compare_sequences`` for each, populating the cache.  Subsequent
        individual calls to ``compare_sequences`` for any pair will then
        be served from the cache.

        Parameters
        ----------
        merge : bool, optional
            Whether to merge consecutive co-linear k-mer runs.  Default
            is ``True``.

        Returns
        -------
        list[tuple[str, str]]
            Ordered list of ``(query_name, target_name)`` pairs that were
            computed.
        """
        ...

    def save(self, path: str) -> None:
        """Serialise the index to a binary file.

        Stores the original sequence bytes and k-mer sets using bincode.
        The FM-index is rebuilt from the sequence bytes when the file is
        loaded, so the on-disk format is compact and version-independent.

        Parameters
        ----------
        path : str
            Destination file path.  The file is created or overwritten.

        Raises
        ------
        ValueError
            If the file cannot be created or serialisation fails.
        """
        ...

    def load(self, path: str) -> None:
        """Load sequences from a previously serialised index file.

        Deserialises sequence bytes and k-mer sets from a file written by
        ``save``, then rebuilds the FM-index for each sequence in memory.

        Parameters
        ----------
        path : str
            Path to a binary index file produced by ``save``.

        Raises
        ------
        ValueError
            If the file cannot be read, deserialisation fails, or the
            k-mer length stored in the file does not match ``self.k``.
        """
        ...

    @property
    def k(self) -> int:
        """The k-mer length used for this index.

        Returns
        -------
        int
            The k-mer length supplied at construction time.
        """
        ...

    def __len__(self) -> int:
        """Return the number of sequences currently held in the index.

        Returns
        -------
        int
            Count of indexed sequences.
        """
        ...

    def __repr__(self) -> str:
        """Return a concise string representation of the index.

        Returns
        -------
        str
            A string of the form ``SequenceIndex(k=<k>, sequences=<n>)``.
        """
        ...


def py_read_fasta(path: str) -> dict[str, str]:
    """Read sequences from a FASTA or gzipped FASTA file.

    Uses needletail for fast parsing.  Gzip decompression is applied
    automatically when the file extension ends in ``.gz``.  All sequences
    are returned in uppercase.

    Parameters
    ----------
    path : str
        Path to a FASTA (``.fa``, ``.fasta``) or gzipped FASTA
        (``.fa.gz``, ``.fasta.gz``) file.

    Returns
    -------
    dict[str, str]
        Mapping of sequence name (record identifier, whitespace-split at
        the first space) to the full uppercase sequence string.

    Raises
    ------
    ValueError
        If the file cannot be opened or its contents cannot be parsed.
    """
    ...


def py_build_kmer_set(seq: str, k: int) -> set[str]:
    """Build the complete set of unique k-mers in a DNA sequence.

    Slides a window of length ``k`` along ``seq`` and collects every
    k-mer that consists exclusively of the characters A, C, G, T.
    K-mers that contain N or other ambiguous bases are silently skipped.

    Parameters
    ----------
    seq : str
        The DNA sequence string.  Uppercase is recommended.
    k : int
        The k-mer length.  Must be greater than zero.

    Returns
    -------
    set[str]
        Set of unique k-mer strings of length ``k`` found in ``seq``.
        Returns an empty set when ``k > len(seq)``.

    Raises
    ------
    ValueError
        If ``k`` is zero.
    """
    ...


def py_find_kmer_coords(seq: str, kmers: list[str]) -> dict[str, list[int]]:
    """Find all start positions of each k-mer in a sequence.

    Builds an FM-index for ``seq`` and performs a backward search for
    every k-mer in ``kmers``.  Only k-mers that actually occur in the
    sequence are included in the result.

    Parameters
    ----------
    seq : str
        The DNA sequence to search in.
    kmers : list[str]
        K-mer strings to look up.

    Returns
    -------
    dict[str, list[int]]
        Mapping of k-mer string to the list of 0-based start positions
        where it occurs in ``seq``.  K-mers with no occurrences are
        omitted.

    Raises
    ------
    ValueError
        If the FM-index cannot be constructed for ``seq``.
    """
    ...


def py_merge_kmer_runs(
    kmer_coords: dict[str, list[int]],
    query_kmer_positions: dict[str, list[int]],
    k: int,
) -> list[tuple[int, int, int, int]]:
    """Merge co-linear k-mer hits into contiguous coordinate blocks.

    For every ``(query_pos, target_pos)`` pair derived from the two
    mappings, this function groups hits by diagonal
    (``target_pos - query_pos``) and, within each diagonal, merges
    consecutive hits where the query position increments by exactly one.

    Parameters
    ----------
    kmer_coords : dict[str, list[int]]
        Mapping of k-mer string to list of 0-based start positions in the
        target sequence.
    query_kmer_positions : dict[str, list[int]]
        Mapping of k-mer string to list of 0-based start positions in the
        query sequence.
    k : int
        The k-mer length, used to compute the (exclusive) end coordinates
        of each merged block.

    Returns
    -------
    list[tuple[int, int, int, int]]
        List of ``(query_start, query_end, target_start, target_end)``
        tuples.  Coordinates are 0-based; end positions are exclusive.
        Results are ordered by diagonal and then by query start position.
    """
    ...


def py_coords_to_paf(
    matches: list[tuple[int, int, int, int]],
    query_name: str,
    query_len: int,
    target_name: str,
    target_len: int,
) -> list[str]:
    """Convert coordinate tuples to PAF (Pairwise mApping Format) strings.

    Each match is formatted as a 12-column, tab-separated PAF line.
    The strand is always reported as ``+`` and mapping quality as ``255``
    (missing).

    Parameters
    ----------
    matches : list[tuple[int, int, int, int]]
        List of ``(query_start, query_end, target_start, target_end)``
        coordinate tuples (0-based, end exclusive).
    query_name : str
        Name of the query sequence (column 1).
    query_len : int
        Total length of the query sequence in base pairs (column 2).
    target_name : str
        Name of the target sequence (column 6).
    target_len : int
        Total length of the target sequence in base pairs (column 7).

    Returns
    -------
    list[str]
        PAF-formatted lines, one per match.  Returns an empty list when
        ``matches`` is empty.
    """
    ...


def py_save_index(path: str, sequences: dict[str, str], k: int) -> None:
    """Build and serialise an index collection to a binary file.

    Constructs FM-indexes and k-mer sets for all provided sequences and
    writes them to a compact bincode file.

    Parameters
    ----------
    path : str
        Destination file path.  The file is created or overwritten.
    sequences : dict[str, str]
        Mapping of sequence name to DNA sequence string.
    k : int
        The k-mer length to use when building k-mer sets.

    Raises
    ------
    ValueError
        If the file cannot be written or serialisation fails.
    """
    ...


def py_load_index(path: str) -> tuple[dict[str, list[str]], int]:
    """Load a previously serialised index collection from a binary file.

    Reads a bincode file written by ``py_save_index`` and returns the
    k-mer sets and k-mer length stored in it.  The FM-indexes themselves
    are not returned; use :class:`SequenceIndex` with ``load`` for a
    fully functional in-memory index.

    Parameters
    ----------
    path : str
        Path to a binary index file produced by ``py_save_index`` or
        :meth:`SequenceIndex.save`.

    Returns
    -------
    tuple[dict[str, list[str]], int]
        A 2-tuple ``(kmer_sets, k)`` where *kmer_sets* maps each sequence
        name to its list of k-mer strings and *k* is the k-mer length
        used when the index was built.

    Raises
    ------
    ValueError
        If the file cannot be read or deserialisation fails.
    """
    ...
