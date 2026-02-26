"""Type stubs for the rusty_dot Rust extension module.

All functions and classes in this module are implemented in Rust via PyO3.
"""

from __future__ import annotations

class SequenceIndex:
    """FM-index backed sequence comparison engine.

    Each sequence added to the index receives its **own independent FM-index**
    built by rust-bio.  The FM-index is constructed once per sequence and
    cannot be extended after construction, so adding more sequences never
    modifies an existing FM-index — it only creates a new one for the
    newly added sequence.

    The index behaves as a **dictionary of per-sequence FM-indexes**:

    * :meth:`add_sequence` and :meth:`load_fasta` **add** new entries to
      the collection; calling either method multiple times accumulates
      sequences rather than replacing the collection.
    * If a sequence name already exists in the index, a ``UserWarning`` is
      emitted and the existing entry is **overwritten** with a new FM-index
      for the new sequence.
    * :meth:`load_fasta` raises ``ValueError`` if the FASTA file itself
      contains duplicate sequence names.
    * Pairwise comparisons (:meth:`compare_sequences`,
      :meth:`compare_sequences_stranded`) always operate on exactly two
      independent FM-indexes.

    The k-mer length ``k`` is fixed at construction time and applies to all
    sequences held in the index.

    Parameters
    ----------
    k : int
        The k-mer length to use for indexing and comparison.

    Examples
    --------
    Build an index from individual sequences:

    >>> idx = SequenceIndex(k=10)
    >>> idx.add_sequence("seq1", "ACGTACGTACGT")
    >>> idx.add_sequence("seq2", "TACGTACGTACG")
    >>> idx.sequence_names()
    ['seq1', 'seq2']

    Accumulate sequences from multiple FASTA files:

    >>> idx = SequenceIndex(k=15)
    >>> idx.load_fasta("assembly_a.fasta")   # adds seqs from file A
    >>> idx.load_fasta("assembly_b.fasta")   # adds seqs from file B, keeps file A seqs
    >>> matches = idx.compare_sequences("seq1", "seq2")

    Overwrite a sequence by re-using its name:

    >>> idx = SequenceIndex(k=10)
    >>> idx.add_sequence("seq1", "ACGTACGT")
    >>> idx.add_sequence("seq1", "GGGGGGGG")  # silently replaces the previous seq1
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

        Builds a **new independent FM-index** for ``seq`` using rust-bio and
        stores it alongside the k-mer set and raw sequence bytes.  Each call
        creates a separate FM-index for that sequence only — the rust-bio
        FM-index cannot be extended after construction, so adding sequences
        never modifies existing FM-indexes.

        Calling ``add_sequence`` does **not** affect any other sequence already
        in the index.  Sequences accumulate: calling this method *N* times with
        *N* distinct names results in an index holding *N* independent
        FM-indexes.

        If a sequence named ``name`` already exists in the index, a
        ``UserWarning`` is emitted and the existing entry is **overwritten**
        with a new FM-index for the new ``seq``.

        Parameters
        ----------
        name : str
            Unique identifier for the sequence.  Re-using an existing name
            emits a :class:`UserWarning` and replaces that sequence (and its
            FM-index).
        seq : str
            DNA sequence string. Uppercase is recommended; lowercase
            input is accepted and treated as uppercase.

        Raises
        ------
        ValueError
            If the FM-index cannot be built (e.g., invalid characters).

        Warns
        -----
        UserWarning
            If ``name`` already exists in the index (the existing entry is
            overwritten).
        """
        ...

    def load_fasta(self, path: str) -> list[str]:
        """Load all sequences from a FASTA or gzipped FASTA file.

        Parses the file with needletail (automatic gzip detection) and
        builds a fresh **independent FM-index** for each record.

        Sequences already in the index are **preserved** — ``load_fasta``
        only adds new entries (or overwrites entries whose name already exists
        in the index).  Calling ``load_fasta`` on two separate files
        accumulates all sequences from both files in the same index.

        If the FASTA file contains **duplicate sequence names** (two records
        with the same identifier), a ``ValueError`` is raised before any
        sequences are added to the index.

        If a record's name **already exists in the index**, a ``UserWarning``
        is emitted and the existing entry is overwritten with the new sequence.

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
            If the file cannot be opened or parsed, or if the file contains
            duplicate sequence names.

        Warns
        -----
        UserWarning
            For each record whose name already exists in the index (those
            entries are overwritten).
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

    def get_paf_all(self, merge: bool = True) -> list[str]:
        """Return PAF-formatted strings for every ordered sequence pair.

        Calls :meth:`get_paf` for every ``(i, j)`` pair where ``i != j``,
        populating the comparison cache as a side-effect.

        Parameters
        ----------
        merge : bool, optional
            Whether to merge consecutive co-linear k-mer runs before
            generating PAF lines.  Default is ``True``.

        Returns
        -------
        list[str]
            All PAF lines for every pairwise comparison, one line per match.
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

    def compare_sequences_stranded(
        self,
        query_name: str,
        target_name: str,
        merge: bool = True,
    ) -> list[tuple[int, int, int, int, str]]:
        """Find shared k-mer matches between two sequences on both strands.

        In addition to the forward (``+``) strand matches returned by
        :meth:`compare_sequences`, this method also searches for k-mers in
        the query whose reverse complement appears in the target, reporting
        those as ``"-"`` strand matches.

        Parameters
        ----------
        query_name : str
            Name of the query sequence.
        target_name : str
            Name of the target sequence.
        merge : bool, optional
            Whether to merge co-linear k-mer runs.  Forward runs are merged
            by diagonal; reverse runs are merged by anti-diagonal.
            Default is ``True``.

        Returns
        -------
        list[tuple[int, int, int, int, str]]
            List of ``(query_start, query_end, target_start, target_end, strand)``
            tuples.  Coordinates are 0-based; end positions are exclusive.
            ``strand`` is ``"+"`` for forward matches and ``"-"`` for
            reverse-complement matches.

        Raises
        ------
        KeyError
            If either sequence name is not present in the index.
        """
        ...

    def optimal_contig_order(
        self,
        query_names: list[str],
        target_names: list[str],
    ) -> tuple[list[str], list[str]]:
        """Return query and target contig names sorted for maximum collinearity.

        Uses the gravity-centre algorithm: for each query contig the gravity
        is the weighted mean of target mid-point positions (normalised by
        total target span) across all matches.  Query contigs with no matches
        are placed at the end.  The same algorithm is applied symmetrically to
        reorder the target contigs.

        Parameters
        ----------
        query_names : list[str]
            Names of the query sequences to reorder.
        target_names : list[str]
            Names of the target sequences to use as the reference axis.

        Returns
        -------
        tuple[list[str], list[str]]
            ``(sorted_query_names, sorted_target_names)`` ordered by ascending
            gravity centre.

        Raises
        ------
        KeyError
            If any sequence name is not present in the index.
        """
        ...

    def save(self, path: str) -> None:
        """Serialise the index to a binary file.

        Stores the original sequence bytes and k-mer sets using postcard.
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

def py_merge_rev_runs(
    target_rev_coords: dict[str, list[int]],
    query_kmer_positions: dict[str, list[int]],
    k: int,
) -> list[tuple[int, int, int, int]]:
    """Merge reverse-complement k-mer hits into contiguous anti-diagonal blocks.

    For every ``(query_pos, target_rev_pos)`` pair derived from the two
    mappings, this function groups hits by anti-diagonal
    (``query_pos + target_rev_pos``) and, within each anti-diagonal, merges
    consecutive hits where ``query_pos`` increments by 1 and
    ``target_rev_pos`` decrements by 1.

    Parameters
    ----------
    target_rev_coords : dict[str, list[int]]
        Mapping of k-mer string to the 0-based start positions of its
        **reverse complement** in the target sequence (as returned by
        ``find_rev_coords_in_index``).
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
        tuples representing merged ``-``-strand match regions.  Coordinates
        are 0-based; end positions are exclusive.  ``target_start`` and
        ``target_end`` are the forward-strand boundaries of the RC region
        on the target.
    """
    ...

def py_merge_rev_fwd_runs(
    target_rev_coords: dict[str, list[int]],
    query_kmer_positions: dict[str, list[int]],
    k: int,
) -> list[tuple[int, int, int, int]]:
    """Merge RC k-mer hits co-linear on a forward diagonal into contiguous blocks.

    Handles the case where sequential query k-mers have their reverse complements
    at **increasing** positions in the target — i.e., as the query position
    advances by 1, the forward-strand target position of the RC match also
    advances by 1 (``t_rc - q = constant``).  This is complementary to
    :func:`py_merge_rev_runs`, which handles the anti-diagonal case
    (``t_rc`` decreases as ``q`` increases).

    This pattern arises in inverted-repeat contexts where both repeat arms
    run in the same left-to-right direction along the forward strand:

    * :func:`py_merge_rev_runs` — anti-diagonal (``q + t_rc = const``):
      standard inverted repeat / reverse-complement alignment.
    * :func:`py_merge_rev_fwd_runs` — forward diagonal (``t_rc - q = const``):
      inverted-repeat case where both arms advance in the same direction.

    Parameters
    ----------
    target_rev_coords : dict[str, list[int]]
        Mapping of k-mer to the 0-based start positions of its **reverse
        complement** in the target sequence (as returned by
        ``find_rev_coords_in_index``).
    query_kmer_positions : dict[str, list[int]]
        Mapping of k-mer to its 0-based start positions in the query sequence.
    k : int
        The k-mer length, used to compute (exclusive) end coordinates.

    Returns
    -------
    list[tuple[int, int, int, int]]
        List of ``(query_start, query_end, target_start, target_end)`` tuples
        representing merged ``-``-strand match regions where RC target positions
        advance together with query positions.  Coordinates are 0-based with
        end positions exclusive.
    """
    ...

def py_merge_runs(
    kmer_coords: dict[str, list[int]],
    query_kmer_positions: dict[str, list[int]],
    k: int,
    strand: str,
) -> list[tuple[int, int, int, int, str]]:
    """Merge k-mer hits into contiguous blocks for either strand orientation.

    A unified entry-point that dispatches to the forward or reverse-complement
    merge algorithm based on ``strand``.

    For ``strand="+"`` the forward (co-linear diagonal) algorithm is used:
    ``kmer_coords`` holds positions of each k-mer in the target and
    ``query_kmer_positions`` holds positions of the same k-mer in the query.

    For ``strand="-"`` **both** the anti-diagonal algorithm
    (:func:`py_merge_rev_runs`) and the co-diagonal algorithm
    (:func:`py_merge_rev_fwd_runs`) are applied to ``kmer_coords`` (positions
    of the reverse complement of each k-mer in the target).  Results from
    both algorithms are combined and deduplicated, so every valid RC alignment
    — whether the RC target positions increase or decrease as the query
    advances — is reported exactly once.

    Parameters
    ----------
    kmer_coords : dict[str, list[int]]
        Mapping of k-mer to 0-based target positions.  For ``strand="-"``,
        these are the positions of the RC of each k-mer in the target.
    query_kmer_positions : dict[str, list[int]]
        Mapping of k-mer to 0-based query positions.
    k : int
        The k-mer length, used to compute (exclusive) end coordinates.
    strand : str
        Orientation of the match: ``"+"`` for forward (co-linear diagonal)
        or ``"-"`` for reverse-complement (both anti-diagonal and co-diagonal
        patterns are merged and returned together, deduplicated).

    Returns
    -------
    list[tuple[int, int, int, int, str]]
        List of ``(query_start, query_end, target_start, target_end, strand)``
        5-tuples.  Coordinates are 0-based; end positions are exclusive.
        ``strand`` echoes the input argument.

    Raises
    ------
    ValueError
        If ``strand`` is neither ``"+"`` nor ``"-"``.
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
    writes them to a compact postcard binary file.

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

    Reads a postcard binary file written by ``py_save_index`` and returns the
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
