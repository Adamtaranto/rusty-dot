"""PAF (Pairwise mApping Format) file I/O with CIGAR string support.

This module provides:
- :class:`PafRecord` — a dataclass representing one PAF alignment line.
- :func:`parse_paf_file` — a generator that yields :class:`PafRecord` objects.
- :class:`PafAlignment` — a container that loads PAF records from a file or
  a list and provides a :meth:`~PafAlignment.reorder_contigs` method to
  reorganise sequences for maximum collinearity, and a
  :meth:`~PafAlignment.filter_by_min_length` method to discard short
  alignments.
- :class:`CrossIndex` — multi-group sequence index for cross-group pairwise
  comparisons (e.g. two or more genome assemblies).  Sequences are organised
  into named groups; alignments are computed between non-self group pairs.
  Compatible with :class:`~rusty_dot.dotplot.DotPlotter`.

CIGAR support
-------------
If a PAF record carries the optional ``cg:Z:<cigar>`` SAM-like tag, the
following alignment statistics are derived automatically:

``alignment_length``
    Total span on the target (M + D + N operations).
``n_matches``
    Number of exact sequence matches (``=`` operations).  If no ``=``
    operations are present, the ``residue_matches`` field from column 10
    is used as a fallback.
``n_mismatches``
    Number of sequence mismatches (``X`` operations).
``n_gaps``
    Number of gap-opening events (count of distinct ``I`` and ``D``
    operations).
``n_gap_bases``
    Total number of inserted or deleted bases.

Examples
--------
>>> from rusty_dot.paf_io import PafAlignment
>>> aln = PafAlignment.from_file("matches.paf")
>>> q_order, t_order = aln.reorder_contigs(aln.query_names, aln.target_names)
"""

from __future__ import annotations

from dataclasses import dataclass, field
import logging
from pathlib import Path
import re
from typing import Any, Generator, Iterable

from rusty_dot._rusty_dot import SequenceIndex

_log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# CIGAR parsing
# ---------------------------------------------------------------------------

_CIGAR_RE = re.compile(r'(\d+)([MIDNSHP=X])')


def _parse_cigar(cigar: str) -> dict[str, int]:
    """Parse a CIGAR string and return a dict of operation counts.

    Parameters
    ----------
    cigar : str
        CIGAR string (e.g. ``"10M2I5D"``).

    Returns
    -------
    dict[str, int]
        Mapping of CIGAR operation character to total count of *bases*
        covered by that operation.

    Raises
    ------
    ValueError
        If the CIGAR string contains invalid characters.
    """
    counts: dict[str, int] = {}
    for length_str, op in _CIGAR_RE.findall(cigar):
        counts[op] = counts.get(op, 0) + int(length_str)
    return counts


def _cigar_stats(cigar: str, residue_matches: int) -> dict[str, int]:
    """Derive alignment statistics from a CIGAR string.

    Parameters
    ----------
    cigar : str
        CIGAR string (e.g. ``"10=2X3I4D"``).
    residue_matches : int
        Fallback value for ``n_matches`` when no ``=`` operations exist.

    Returns
    -------
    dict[str, int]
        Keys: ``alignment_length``, ``n_matches``, ``n_mismatches``,
        ``n_gaps``, ``n_gap_bases``.
    """
    ops = _parse_cigar(cigar)
    # Alignment length = bases consumed on the target
    alignment_length = (
        ops.get('M', 0)
        + ops.get('D', 0)
        + ops.get('N', 0)
        + ops.get('=', 0)
        + ops.get('X', 0)
    )
    n_matches = ops.get('=', residue_matches)
    n_mismatches = ops.get('X', 0)
    n_gap_bases = ops.get('I', 0) + ops.get('D', 0)
    # Count distinct gap-opening events (each I or D run counts as one event).
    n_gaps = sum(1 for _, op in _CIGAR_RE.findall(cigar) if op in ('I', 'D'))
    return {
        'alignment_length': alignment_length,
        'n_matches': n_matches,
        'n_mismatches': n_mismatches,
        'n_gaps': n_gaps,
        'n_gap_bases': n_gap_bases,
    }


# ---------------------------------------------------------------------------
# PafRecord dataclass
# ---------------------------------------------------------------------------


@dataclass
class PafRecord:
    """A single PAF alignment record.

    The twelve required PAF columns are represented as typed attributes.
    Optional SAM-like tags (e.g. ``tp:A:P``, ``cg:Z:10M``) are stored in
    :attr:`tags`.  If a ``cg:Z:`` tag is present, CIGAR-derived alignment
    statistics are populated automatically.

    Parameters
    ----------
    query_name : str
        Query sequence name (column 1).
    query_len : int
        Query sequence length (column 2).
    query_start : int
        Query start position, 0-based (column 3).
    query_end : int
        Query end position, exclusive (column 4).
    strand : str
        Relative strand: ``"+"`` or ``"-"`` (column 5).
    target_name : str
        Target sequence name (column 6).
    target_len : int
        Target sequence length (column 7).
    target_start : int
        Target start position, 0-based (column 8).
    target_end : int
        Target end position, exclusive (column 9).
    residue_matches : int
        Number of residue matches (column 10).
    alignment_block_len : int
        Number of bases in the alignment block (column 11).
    mapping_quality : int
        Mapping quality (0–255; 255 = missing) (column 12).
    tags : dict[str, Any]
        Optional SAM-like tags decoded as ``{tag_name: value}``.
    cigar : str or None
        CIGAR string from ``cg:Z:`` tag, or ``None`` if absent.
    alignment_length : int or None
        Target-span alignment length derived from CIGAR, or ``None``.
    n_matches : int or None
        Count of exact-match bases (``=`` ops) from CIGAR; falls back to
        ``residue_matches`` when only ``M`` ops are present.
    n_mismatches : int or None
        Count of mismatch bases (``X`` ops) from CIGAR, or ``None``.
    n_gaps : int or None
        Total number of gap bases (``I`` + ``D`` bases) from CIGAR.
    n_gap_bases : int or None
        Same as ``n_gaps`` (alias kept for clarity).
    """

    query_name: str
    query_len: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_len: int
    target_start: int
    target_end: int
    residue_matches: int
    alignment_block_len: int
    mapping_quality: int
    tags: dict[str, Any] = field(default_factory=dict)
    cigar: str | None = None
    alignment_length: int | None = None
    n_matches: int | None = None
    n_mismatches: int | None = None
    n_gaps: int | None = None
    n_gap_bases: int | None = None

    @property
    def query_aligned_len(self) -> int:
        """Return the aligned length on the query sequence.

        Returns
        -------
        int
            ``query_end - query_start``.
        """
        return self.query_end - self.query_start

    @property
    def target_aligned_len(self) -> int:
        """Return the aligned length on the target sequence.

        Returns
        -------
        int
            ``target_end - target_start``.
        """
        return self.target_end - self.target_start

    @classmethod
    def from_line(cls, line: str) -> 'PafRecord':
        """Parse a single PAF text line into a :class:`PafRecord`.

        Parameters
        ----------
        line : str
            A single PAF record line (tab-separated, trailing newline optional).

        Returns
        -------
        PafRecord
            The parsed record.

        Raises
        ------
        ValueError
            If the line has fewer than 12 tab-separated fields.
        """
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 12:
            raise ValueError(
                f'PAF line has {len(fields)} fields; expected at least 12: {line!r}'
            )
        tags: dict[str, Any] = {}
        cigar: str | None = None
        for tag_field in fields[12:]:
            parts = tag_field.split(':', 2)
            if len(parts) == 3:
                tag_name, tag_type, tag_value = parts
                if tag_type == 'i':
                    tags[tag_name] = int(tag_value)
                elif tag_type == 'f':
                    tags[tag_name] = float(tag_value)
                else:
                    tags[tag_name] = tag_value
                if tag_name == 'cg' and tag_type == 'Z':
                    cigar = tag_value

        residue_matches = int(fields[9])
        stats: dict[str, int] = {}
        if cigar is not None:
            stats = _cigar_stats(cigar, residue_matches)

        return cls(
            query_name=fields[0],
            query_len=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            strand=fields[4],
            target_name=fields[5],
            target_len=int(fields[6]),
            target_start=int(fields[7]),
            target_end=int(fields[8]),
            residue_matches=residue_matches,
            alignment_block_len=int(fields[10]),
            mapping_quality=int(fields[11]),
            tags=tags,
            cigar=cigar,
            alignment_length=stats.get('alignment_length'),
            n_matches=stats.get('n_matches'),
            n_mismatches=stats.get('n_mismatches'),
            n_gaps=stats.get('n_gaps'),
            n_gap_bases=stats.get('n_gap_bases'),
        )

    def to_line(self) -> str:
        """Serialise this record back to a PAF-format string (no trailing newline).

        Returns
        -------
        str
            Tab-separated PAF line with the 12 required columns.  Optional
            tags are not included.
        """
        return '\t'.join(
            str(v)
            for v in [
                self.query_name,
                self.query_len,
                self.query_start,
                self.query_end,
                self.strand,
                self.target_name,
                self.target_len,
                self.target_start,
                self.target_end,
                self.residue_matches,
                self.alignment_block_len,
                self.mapping_quality,
            ]
        )


# ---------------------------------------------------------------------------
# Generator helper
# ---------------------------------------------------------------------------


def parse_paf_file(path: str | Path) -> Generator[PafRecord, None, None]:
    """Yield :class:`PafRecord` objects from a PAF file.

    Lines beginning with ``#`` are treated as comments and skipped.  Empty
    lines are also skipped.

    Parameters
    ----------
    path : str or Path
        Path to the PAF file.

    Yields
    ------
    PafRecord
        One record per non-comment, non-empty line.

    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    ValueError
        If a line cannot be parsed as a PAF record.
    """
    path = Path(path)
    with path.open('r', encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            yield PafRecord.from_line(line)


# ---------------------------------------------------------------------------
# Gravity-based contig reordering (pure Python)
# ---------------------------------------------------------------------------


def compute_gravity_contigs(
    records: Iterable[PafRecord],
    query_names: list[str],
    target_names: list[str],
) -> tuple[list[str], list[str]]:
    """Return query and target contig names sorted by gravity centre.

    For each query contig the gravity centre is the weighted mean of target
    mid-point positions (normalised by the total target span) across all
    alignment records that involve that contig.  Target contigs are sorted
    symmetrically against the query axis.

    Contigs with no alignment records receive a gravity of ``float("inf")``
    and are placed at the end of the sorted list.

    Parameters
    ----------
    records : iterable of PafRecord
        Alignment records to use for computing gravity centres.
    query_names : list[str]
        The query contig names to reorder.
    target_names : list[str]
        The target contig names to reorder.

    Returns
    -------
    tuple[list[str], list[str]]
        ``(sorted_query_names, sorted_target_names)`` ordered by ascending
        gravity centre.
    """
    query_set = set(query_names)
    target_set = set(target_names)

    # Collect all records into a list and build sequence-length maps from them.
    q_len_map: dict[str, int] = {}
    t_len_map: dict[str, int] = {}
    all_records: list[PafRecord] = []
    for rec in records:
        all_records.append(rec)
        q_len_map[rec.query_name] = rec.query_len
        t_len_map[rec.target_name] = rec.target_len

    # Build cumulative target offsets using actual sequence lengths.
    t_offsets: dict[str, int] = {}
    t_off = 0
    for t in target_names:
        t_offsets[t] = t_off
        t_off += t_len_map.get(t, 1)
    total_target_len = max(t_off, 1)

    # Build cumulative query offsets using actual sequence lengths.
    q_offsets_real: dict[str, int] = {}
    q_off = 0
    for q in query_names:
        q_offsets_real[q] = q_off
        q_off += q_len_map.get(q, 1)
    total_query_len = max(q_off, 1)

    # Accumulate weighted positions.
    q_weight: dict[str, float] = dict.fromkeys(query_names, 0.0)
    q_wpos: dict[str, float] = dict.fromkeys(query_names, 0.0)
    t_weight: dict[str, float] = dict.fromkeys(target_names, 0.0)
    t_wpos: dict[str, float] = dict.fromkeys(target_names, 0.0)

    for rec in all_records:
        if rec.query_name not in query_set or rec.target_name not in target_set:
            continue
        size = float(rec.alignment_block_len or (rec.query_end - rec.query_start))
        if size <= 0:
            continue

        # Target gravity from query's perspective.
        t_mid = (
            t_offsets.get(rec.target_name, 0)
            + (rec.target_start + rec.target_end) / 2.0
        )
        q_weight[rec.query_name] += size
        q_wpos[rec.query_name] += size * t_mid

        # Query gravity from target's perspective.
        q_mid = (
            q_offsets_real.get(rec.query_name, 0)
            + (rec.query_start + rec.query_end) / 2.0
        )
        t_weight[rec.target_name] += size
        t_wpos[rec.target_name] += size * q_mid

    def _gravity(name: str, wt: dict, wp: dict, total: float) -> float:
        w = wt.get(name, 0.0)
        return (wp.get(name, 0.0) / w / total) if w > 0 else float('inf')

    def _sort_key_with_len(
        name: str,
        wt: dict,
        wp: dict,
        total: float,
        len_map: dict,
    ) -> tuple:
        g = _gravity(name, wt, wp, total)
        if g == float('inf'):
            # Unmatched: sort after matched (1 > 0), then by descending length
            return (1, -len_map.get(name, 0))
        return (0, g)

    sorted_q = sorted(
        query_names,
        key=lambda n: _sort_key_with_len(
            n, q_weight, q_wpos, total_target_len, q_len_map
        ),
    )
    sorted_t = sorted(
        target_names,
        key=lambda n: _sort_key_with_len(
            n, t_weight, t_wpos, total_query_len, t_len_map
        ),
    )
    return sorted_q, sorted_t


# ---------------------------------------------------------------------------
# PafAlignment class
# ---------------------------------------------------------------------------


class PafAlignment:
    """A collection of PAF alignment records with contig-ordering utilities.

    Can be constructed from a file path or an iterable of :class:`PafRecord`
    objects.  Provides :meth:`reorder_contigs` to sort query and target
    sequence names so that a subsequent dotplot shows maximum collinearity.

    Parameters
    ----------
    records : list of PafRecord
        The alignment records.

    Examples
    --------
    Load from a file and reorder contigs:

    >>> aln = PafAlignment.from_file("alignments.paf")
    >>> q_order, t_order = aln.reorder_contigs(aln.query_names, aln.target_names)
    """

    def __init__(self, records: list[PafRecord]) -> None:
        self.records: list[PafRecord] = records
        # Custom group assignments.  None means use the default (query_names
        # → 'a', target_names → 'b') which is computed lazily from records.
        self._groups: dict[str, list[str]] | None = None

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_file(cls, path: str | Path) -> 'PafAlignment':
        """Load records from a PAF file.

        Parameters
        ----------
        path : str or Path
            Path to the PAF file.

        Returns
        -------
        PafAlignment
            New instance with all records loaded.
        """
        return cls(list(parse_paf_file(path)))

    @classmethod
    def from_records(cls, records: Iterable[PafRecord]) -> 'PafAlignment':
        """Construct from an iterable of :class:`PafRecord` objects.

        Parameters
        ----------
        records : iterable of PafRecord
            Source records.

        Returns
        -------
        PafAlignment
            New instance.
        """
        return cls(list(records))

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def query_names(self) -> list[str]:
        """Return a deduplicated list of query sequence names (insertion order).

        Returns
        -------
        list[str]
            Unique query names in the order first seen.
        """
        seen: dict[str, None] = {}
        for rec in self.records:
            seen[rec.query_name] = None
        return list(seen)

    @property
    def target_names(self) -> list[str]:
        """Return a deduplicated list of target sequence names (insertion order).

        Returns
        -------
        list[str]
            Unique target names in the order first seen.
        """
        seen: dict[str, None] = {}
        for rec in self.records:
            seen[rec.target_name] = None
        return list(seen)

    def sequence_names(self) -> list[str]:
        """Return a deduplicated list of all query and target sequence names.

        The list contains each name at most once, in the order it was first
        encountered (queries before targets within each record).  This method
        makes :class:`PafAlignment` compatible with :class:`~rusty_dot.dotplot.DotPlotter`.

        Returns
        -------
        list[str]
            All unique sequence names across query and target fields.
        """
        seen: dict[str, None] = {}
        for rec in self.records:
            seen[rec.query_name] = None
            seen[rec.target_name] = None
        return list(seen)

    def get_sequence_length(self, name: str) -> int:
        """Return the length of a sequence by name as stored in PAF records.

        Looks up *name* in the ``query_name`` and ``target_name`` fields of
        every record and returns the corresponding ``query_len`` or
        ``target_len``.  This method makes :class:`PafAlignment` compatible
        with :class:`~rusty_dot.dotplot.DotPlotter`.

        Parameters
        ----------
        name : str
            Sequence name to look up.

        Returns
        -------
        int
            Length of the sequence.

        Raises
        ------
        KeyError
            If *name* is not found in any record.
        """
        for rec in self.records:
            if rec.query_name == name:
                return rec.query_len
            if rec.target_name == name:
                return rec.target_len
        raise KeyError(f'Sequence {name!r} not found in PAF records.')

    def __len__(self) -> int:
        """Return the number of records.

        Returns
        -------
        int
            Record count.
        """
        return len(self.records)

    def __repr__(self) -> str:
        """Return a concise string representation.

        Returns
        -------
        str
            ``PafAlignment(records=<n>, queries=<q>, targets=<t>)``.
        """
        return (
            f'PafAlignment(records={len(self.records)}, '
            f'queries={len(self.query_names)}, '
            f'targets={len(self.target_names)})'
        )

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def filter_by_query(self, names: Iterable[str]) -> 'PafAlignment':
        """Return a new :class:`PafAlignment` containing only the given query names.

        Parameters
        ----------
        names : iterable of str
            Query names to keep.

        Returns
        -------
        PafAlignment
            Filtered alignment.
        """
        keep = set(names)
        return PafAlignment([r for r in self.records if r.query_name in keep])

    def filter_by_target(self, names: Iterable[str]) -> 'PafAlignment':
        """Return a new :class:`PafAlignment` containing only the given target names.

        Parameters
        ----------
        names : iterable of str
            Target names to keep.

        Returns
        -------
        PafAlignment
            Filtered alignment.
        """
        keep = set(names)
        return PafAlignment([r for r in self.records if r.target_name in keep])

    def filter_by_min_length(self, min_length: int) -> 'PafAlignment':
        """Return a new :class:`PafAlignment` keeping only records of sufficient length.

        Filters on the query aligned length (``query_end - query_start``), which
        equals the alignment block span for both merged k-mer runs and PAF
        alignments imported from a file.

        Parameters
        ----------
        min_length : int
            Minimum alignment length (inclusive).  Records with a query aligned
            length strictly less than ``min_length`` are discarded.

        Returns
        -------
        PafAlignment
            Filtered alignment containing only records with
            ``query_aligned_len >= min_length``.
        """
        return PafAlignment(
            [r for r in self.records if r.query_aligned_len >= min_length]
        )

    # ------------------------------------------------------------------
    # Group management
    # ------------------------------------------------------------------

    @property
    def groups(self) -> dict[str, list[str]]:
        """Return the current group assignments.

        If groups have not been set explicitly via :meth:`set_groups` or
        :meth:`rename_group`, returns the default: all query sequence names
        in group ``'a'`` and all target sequence names in group ``'b'``.

        Returns
        -------
        dict[str, list[str]]
            Mapping of group label → list of sequence names.
        """
        if self._groups is not None:
            return dict(self._groups)
        return {'a': self.query_names, 'b': self.target_names}

    def set_groups(self, groups: dict[str, list[str]]) -> None:
        """Set custom group assignments for sequence names.

        Parameters
        ----------
        groups : dict[str, list[str]]
            Mapping of group label → list of sequence names belonging to
            that group.

        Warns
        -----
        Logs a warning for every sequence name that appears in more than one
        group.
        """
        seen: dict[str, str] = {}
        for group, names in groups.items():
            for name in names:
                if name in seen:
                    _log.warning(
                        'PafAlignment.set_groups: sequence %r is assigned to '
                        'both group %r and group %r',
                        name,
                        seen[name],
                        group,
                    )
                else:
                    seen[name] = group
        self._groups = {g: list(ns) for g, ns in groups.items()}

    def rename_group(self, old_name: str, new_name: str) -> None:
        """Rename a group label.

        If custom groups have not been set yet, the default assignment
        (``'a'`` → query names, ``'b'`` → target names) is materialised
        first.

        Parameters
        ----------
        old_name : str
            Current group label.
        new_name : str
            New group label.

        Raises
        ------
        KeyError
            If *old_name* is not a known group.
        ValueError
            If *new_name* already exists as a different group label.
        """
        current = self._groups if self._groups is not None else self.groups
        if old_name not in current:
            raise KeyError(f"Group {old_name!r} not found.")
        if new_name in current and new_name != old_name:
            raise ValueError(f"Group {new_name!r} already exists.")
        self._groups = {
            (new_name if k == old_name else k): v for k, v in current.items()
        }
        _log.info('PafAlignment: renamed group %r → %r', old_name, new_name)

    # ------------------------------------------------------------------
    # Contig reordering
    # ------------------------------------------------------------------

    def reorder_contigs(
        self,
        query_names: list[str] | None = None,
        target_names: list[str] | None = None,
        query_group: str | None = None,
        target_group: str | None = None,
    ) -> tuple[list[str], list[str]]:
        """Sort query and target contigs to maximise collinearity in the dotplot.

        Uses the gravity-centre algorithm: each contig is assigned a gravity
        equal to the weighted mean position of its alignment blocks on the
        opposing axis.  Contigs are then sorted by ascending gravity.

        Parameters
        ----------
        query_names : list[str] or None, optional
            Query contigs to reorder.  Ignored when *query_group* is given.
            Defaults to :attr:`query_names`.
        target_names : list[str] or None, optional
            Target contigs to reorder.  Ignored when *target_group* is given.
            Defaults to :attr:`target_names`.
        query_group : str or None, optional
            Group label whose members are used as query contigs.  When
            provided, the corresponding entry in :attr:`groups` is used and
            *query_names* is ignored.
        target_group : str or None, optional
            Group label whose members are used as target contigs.  When
            provided, the corresponding entry in :attr:`groups` is used and
            *target_names* is ignored.

        Returns
        -------
        tuple[list[str], list[str]]
            ``(sorted_query_names, sorted_target_names)``.

        Raises
        ------
        KeyError
            If a supplied group label is not present in :attr:`groups`.
        """
        current_groups = self.groups
        if query_group is not None:
            if query_group not in current_groups:
                raise KeyError(f"Group {query_group!r} not found.")
            q = current_groups[query_group]
        else:
            q = query_names if query_names is not None else self.query_names

        if target_group is not None:
            if target_group not in current_groups:
                raise KeyError(f"Group {target_group!r} not found.")
            t = current_groups[target_group]
        else:
            t = target_names if target_names is not None else self.target_names

        return compute_gravity_contigs(self.records, q, t)


# ---------------------------------------------------------------------------
# CrossIndex: multi-group pairwise sequence comparisons
# ---------------------------------------------------------------------------


class CrossIndex:
    """Multi-group sequence index for cross-group pairwise comparisons.

    Sequences are organised into named groups (e.g. ``'assembly_a'``,
    ``'assembly_b'``).  Each sequence is stored in a shared
    :class:`~rusty_dot.SequenceIndex` under a ``group:name`` internal key,
    which keeps names unique even when the same sequence identifier appears
    in multiple groups.

    **Workflow**

    Loading sequences and computing matches are separate, explicit steps:

    1. Load sequences via :meth:`add_sequence` or :meth:`load_fasta`.
       Each sequence addition is logged at ``DEBUG`` level.  A ``WARNING``
       is emitted if the name already exists in the same group (FM-index
       overwritten) or in a different group.
    2. Call :meth:`compute_matches` to find k-mer matches between groups.
       This must be done before calling :meth:`reorder_contigs` or
       :meth:`reorder_for_colinearity`.
    3. Inspect :attr:`computed_group_pairs` to verify which pairs have been
       computed.

    **Alignment scope by number of groups**

    * **2 groups** — :meth:`compute_matches` compares the two groups by
      default.
    * **3+ groups** — all non-self ordered pairs by default.  Use the
      *query_group* / *target_group* arguments to restrict to a specific pair.

    **DotPlotter compatibility**

    ``CrossIndex`` exposes :meth:`get_sequence_length`,
    :meth:`compare_sequences_stranded`, and :meth:`sequence_names` so that it
    can be passed directly to :class:`~rusty_dot.dotplot.DotPlotter`::

        cross = CrossIndex(k=15)
        cross.load_fasta("assembly_a.fasta", group="a")
        cross.load_fasta("assembly_b.fasta", group="b")
        cross.compute_matches()

        from rusty_dot.dotplot import DotPlotter
        plotter = DotPlotter(cross)
        plotter.plot(
            query_names=cross.sequence_names(group="a"),
            target_names=cross.sequence_names(group="b"),
            output_path="cross_plot.png",
        )

    Parameters
    ----------
    k : int
        K-mer length to use for indexing and comparison.

    Examples
    --------
    >>> from rusty_dot.paf_io import CrossIndex
    >>> cross = CrossIndex(k=10)
    >>> cross.load_fasta("genome_a.fasta", group="a")
    >>> cross.load_fasta("genome_b.fasta", group="b")
    >>> cross.compute_matches()
    >>> paf_lines = cross.get_paf()
    """

    def __init__(self, k: int) -> None:
        """Initialise an empty CrossIndex.

        Parameters
        ----------
        k : int
            K-mer length to use when building the sequence index.
        """
        self._k: int = k
        self._index: SequenceIndex = SequenceIndex(k=k)
        # group_label -> ordered list of original (un-prefixed) sequence names
        self._groups: dict[str, list[str]] = {}
        # group_label -> internal prefix used in _index (supports rename_group)
        self._internal_group: dict[str, str] = {}
        # (query_group, target_group) -> list[PafRecord] from compute_matches
        self._records_by_pair: dict[tuple[str, str], list[PafRecord]] = {}

    @property
    def _paf_records(self) -> list[PafRecord]:
        """All cached PAF records across every computed group pair (flat list).

        Note: this property rebuilds the list on every access.  Avoid calling
        it repeatedly in a tight loop; assign to a local variable instead.
        """
        result: list[PafRecord] = []
        for recs in self._records_by_pair.values():
            result.extend(recs)
        return result

    @property
    def computed_group_pairs(self) -> list[tuple[str, str]]:
        """Group pairs for which k-mer matches have been computed.

        Returns
        -------
        list[tuple[str, str]]
            List of ``(query_group, target_group)`` pairs that have had
            :meth:`compute_matches` run on them.  Use this to confirm that
            the required pair is ready before calling
            :meth:`reorder_contigs`.
        """
        return list(self._records_by_pair.keys())

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _make_internal(self, group: str, name: str) -> str:
        """Format an internal (prefixed) name for use in SequenceIndex.

        Uses :attr:`_internal_group` to resolve the actual internal prefix,
        which may differ from *group* after a :meth:`rename_group` call.
        The ``get(group, group)`` fallback is safe because every call to
        :meth:`add_sequence` or :meth:`load_fasta` immediately registers the
        new group in ``_internal_group`` before any internal name is formed.
        """
        prefix = self._internal_group.get(group, group)
        return f'{prefix}:{name}'

    @staticmethod
    def _split_internal(internal: str) -> tuple[str, str]:
        """Split ``'group:name'`` into ``(group, name)``."""
        group, _, name = internal.partition(':')
        return group, name

    # ------------------------------------------------------------------
    # Adding sequences
    # ------------------------------------------------------------------

    def _check_name_collision(self, name: str, group: str) -> None:
        """Emit warnings if *name* already exists in the same or another group."""
        if name in self._groups.get(group, []):
            _log.warning(
                "CrossIndex: sequence %r already exists in group %r; "
                "its FM-index will be overwritten",
                name,
                group,
            )
        else:
            for other_g, other_names in self._groups.items():
                if other_g != group and name in other_names:
                    _log.warning(
                        "CrossIndex: sequence %r already exists in group %r; "
                        "adding the same name to group %r may cause confusion",
                        name,
                        other_g,
                        group,
                    )

    def add_sequence(self, name: str, seq: str, group: str = 'a') -> None:
        """Add a single sequence to the specified group.

        Logs a ``DEBUG``-level message for every sequence loaded, and a
        ``WARNING`` if *name* already exists in the same group (its FM-index
        will be overwritten) or in a different group (potential confusion).

        Parameters
        ----------
        name : str
            Sequence identifier (must be unique within the group).
        seq : str
            DNA sequence string.
        group : str, optional
            Group label.  Any non-empty string is accepted; ``':'`` is
            forbidden because it is used as an internal separator.
            Default is ``'a'``.

        Raises
        ------
        ValueError
            If *group* contains ``':'``.
        """
        if ':' in group:
            raise ValueError(f"Group name must not contain ':', got {group!r}")
        self._check_name_collision(name, group)
        _log.debug(
            "CrossIndex: adding sequence %r (len=%d) to group %r",
            name,
            len(seq),
            group,
        )
        if group not in self._groups:
            self._groups[group] = []
            self._internal_group[group] = group
        internal = self._make_internal(group, name)
        self._index.add_sequence(internal, seq)
        if name not in self._groups[group]:
            self._groups[group].append(name)

    def load_fasta(self, path: str, group: str = 'a') -> list[str]:
        """Load sequences from a FASTA file into the specified group.

        Logs an ``INFO``-level message when opening the file, a ``DEBUG``
        message for each sequence loaded (including sequence name and length),
        and a ``WARNING`` if any sequence name already exists in the same or
        another group.

        Parameters
        ----------
        path : str
            Path to a FASTA (``.fa`` / ``.fasta``) or gzipped FASTA file.
        group : str, optional
            Group label.  Default is ``'a'``.

        Returns
        -------
        list[str]
            The original (un-prefixed) sequence names that were loaded, in
            file order.

        Raises
        ------
        ValueError
            If *group* contains ``':'``, or if the file cannot be parsed, or
            if the FASTA file contains duplicate sequence names.
        """
        if ':' in group:
            raise ValueError(f"Group name must not contain ':', got {group!r}")
        from rusty_dot._rusty_dot import py_read_fasta

        _log.info("CrossIndex: loading sequences from %r into group %r", path, group)
        seqs = py_read_fasta(path)
        if group not in self._groups:
            self._groups[group] = []
            self._internal_group[group] = group
        names: list[str] = []
        for name, seq in seqs.items():
            self._check_name_collision(name, group)
            _log.debug(
                "CrossIndex: adding sequence %r (len=%d) to group %r",
                name,
                len(seq),
                group,
            )
            internal = self._make_internal(group, name)
            self._index.add_sequence(internal, seq)
            if name not in self._groups[group]:
                self._groups[group].append(name)
            names.append(name)
        _log.info(
            "CrossIndex: loaded %d sequence(s) from %r into group %r",
            len(names),
            path,
            group,
        )
        return names

    # ------------------------------------------------------------------
    # Properties and name helpers
    # ------------------------------------------------------------------

    @property
    def group_names(self) -> list[str]:
        """Return the list of group labels in insertion order.

        Returns
        -------
        list[str]
            Group labels.
        """
        return list(self._groups.keys())

    def sequence_names(self, group: str | None = None) -> list[str]:
        """Return internal (``group:name``) identifiers suitable for DotPlotter.

        Parameters
        ----------
        group : str or None, optional
            If given, return only names from that group.  If ``None``
            (default), return names from all groups.

        Returns
        -------
        list[str]
            Internal ``group:name`` strings in current :attr:`contig_order`.
        """
        if group is not None:
            return [self._make_internal(group, n) for n in self._groups.get(group, [])]
        result: list[str] = []
        for g, names in self._groups.items():
            result.extend(self._make_internal(g, n) for n in names)
        return result

    @property
    def contig_order(self) -> dict[str, list[str]]:
        """Current contig order per group as original (un-prefixed) names.

        Returns
        -------
        dict[str, list[str]]
            Mapping of group label → ordered list of sequence names.
            Updated in-place by :meth:`reorder_by_length` and
            :meth:`reorder_for_colinearity`.
        """
        return {g: list(names) for g, names in self._groups.items()}

    # ------------------------------------------------------------------
    # Backward-compatible properties (two-group a/b model)
    # ------------------------------------------------------------------

    @property
    def query_names(self) -> list[str]:
        """Un-prefixed names for group ``'a'`` (backward compatible).

        Returns
        -------
        list[str]
        """
        return list(self._groups.get('a', []))

    @property
    def target_names(self) -> list[str]:
        """Un-prefixed names for group ``'b'`` (backward compatible).

        Returns
        -------
        list[str]
        """
        return list(self._groups.get('b', []))

    # ------------------------------------------------------------------
    # DotPlotter-compatible interface
    # ------------------------------------------------------------------

    def get_sequence_length(self, name: str) -> int:
        """Return the length of the sequence identified by its internal name.

        Parameters
        ----------
        name : str
            Internal (``group:name``) identifier.

        Returns
        -------
        int
            Sequence length in bases.
        """
        return self._index.get_sequence_length(name)

    def compare_sequences_stranded(
        self, name1: str, name2: str, merge: bool = True
    ) -> list:
        """Compare two sequences by their internal names, returning stranded matches.

        Parameters
        ----------
        name1 : str
            Internal name of the query sequence.
        name2 : str
            Internal name of the target sequence.
        merge : bool, optional
            Whether to merge consecutive co-linear k-mer runs.
            Default is ``True``.

        Returns
        -------
        list of (int, int, int, int, str)
            List of ``(query_start, query_end, target_start, target_end, strand)``
            tuples.
        """
        return self._index.compare_sequences_stranded(name1, name2, merge)

    # ------------------------------------------------------------------
    # Contig reordering
    # ------------------------------------------------------------------

    def rename_group(self, old_name: str, new_name: str) -> None:
        """Rename a group label without re-indexing sequences.

        The internal prefix used in the underlying :class:`SequenceIndex` is
        preserved; only the public-facing group label is changed.

        Parameters
        ----------
        old_name : str
            Current group label to rename.
        new_name : str
            New group label.  Must not contain ``':'``.

        Raises
        ------
        KeyError
            If *old_name* is not a known group.
        ValueError
            If *new_name* contains ``':'`` or already exists as a group label.
        """
        if old_name not in self._groups:
            raise KeyError(f"Group {old_name!r} not found.")
        if ':' in new_name:
            raise ValueError(f"Group name must not contain ':', got {new_name!r}")
        if new_name in self._groups and new_name != old_name:
            raise ValueError(f"Group {new_name!r} already exists.")
        # Rebuild _groups preserving insertion order
        self._groups = {
            (new_name if k == old_name else k): v for k, v in self._groups.items()
        }
        # Update the internal-prefix mapping
        self._internal_group[new_name] = self._internal_group.pop(old_name)
        _log.info('CrossIndex: renamed group %r → %r', old_name, new_name)

    def set_group_members(self, group: str, names: list[str]) -> None:
        """Assign a custom list of sequence names to an existing group.

        Only the logical membership list is updated; sequences already indexed
        are not moved or removed from the underlying
        :class:`~rusty_dot.SequenceIndex`.

        Parameters
        ----------
        group : str
            Group label to update.  The group must already exist.
        names : list[str]
            New ordered list of un-prefixed sequence names for the group.

        Raises
        ------
        KeyError
            If *group* is not a known group.

        Warns
        -----
        Logs a warning for every name that is also present in another group.
        """
        if group not in self._groups:
            raise KeyError(f"Group {group!r} not found.")
        for n in names:
            for other_g, other_ns in self._groups.items():
                if other_g != group and n in other_ns:
                    _log.warning(
                        "CrossIndex: sequence %r is assigned to both group %r and group %r",
                        n,
                        other_g,
                        group,
                    )
        self._groups[group] = list(names)

    def reorder_by_length(self, group: str | None = None) -> None:
        """Reorder contigs within one or all groups by descending sequence length.

        Updates :attr:`contig_order` in-place.

        Parameters
        ----------
        group : str or None, optional
            Group to reorder.  If ``None`` (default), all groups are reordered.
        """
        groups_to_sort = [group] if group is not None else list(self._groups.keys())
        for g in groups_to_sort:
            self._groups[g].sort(
                key=lambda n: self._index.get_sequence_length(
                    self._make_internal(g, n)
                ),
                reverse=True,
            )

    def reorder_for_colinearity(self, query_group: str, target_group: str) -> None:
        """Reorder sequences in two groups to maximise dotplot collinearity.

        Uses the gravity-centre algorithm via
        :meth:`~rusty_dot.SequenceIndex.optimal_contig_order`.  Updates
        :attr:`contig_order` in-place for both groups.

        .. note::
            :meth:`compute_matches` must be called for ``(query_group,
            target_group)`` before calling this method.

        Parameters
        ----------
        query_group : str
            Group label for the query (y-axis / rows).
        target_group : str
            Group label for the target (x-axis / columns).

        Raises
        ------
        KeyError
            If either group label is not present in the index.
        ValueError
            If :meth:`compute_matches` has not been called for this group
            pair.
        """
        pair = (query_group, target_group)
        if pair not in self._records_by_pair:
            raise ValueError(
                f"No matches computed for group pair {pair!r}. "
                "Call compute_matches() for this pair first."
            )
        q_internal = [
            self._make_internal(query_group, n) for n in self._groups[query_group]
        ]
        t_internal = [
            self._make_internal(target_group, n) for n in self._groups[target_group]
        ]
        sorted_q_int, sorted_t_int = self._index.optimal_contig_order(
            q_internal, t_internal
        )
        self._groups[query_group] = [self._split_internal(n)[1] for n in sorted_q_int]
        self._groups[target_group] = [self._split_internal(n)[1] for n in sorted_t_int]

    # ------------------------------------------------------------------
    # PAF output and match computation
    # ------------------------------------------------------------------

    def _get_default_group_pairs(self) -> list[tuple[str, str]]:
        """Return default group pairs for alignment.

        * 2 groups → one pair between the two groups.
        * 3+ groups → all non-self ordered pairs.
        """
        groups = list(self._groups.keys())
        if len(groups) == 2:
            return [(groups[0], groups[1])]
        return [
            (a, b) for i, a in enumerate(groups) for j, b in enumerate(groups) if i != j
        ]

    def compute_matches(
        self,
        query_group: str | None = None,
        target_group: str | None = None,
        merge: bool = True,
    ) -> None:
        """Compute k-mer matches between groups and cache the results.

        This is the primary computation step and must be called **before**
        :meth:`reorder_contigs` or :meth:`reorder_for_colinearity`.  Matches
        are computed only between groups — not within a single group.

        When *query_group* and *target_group* are both ``None``:

        * **2 groups** — the single cross-group pair is used.
        * **3+ groups** — all non-self ordered pairs are computed.

        The computed records are stored internally, keyed by
        ``(query_group, target_group)``, and the pair is added to
        :attr:`computed_group_pairs`.

        Parameters
        ----------
        query_group : str or None, optional
            Group label for query sequences.  When ``None`` (default) the
            groups are auto-detected (see above).
        target_group : str or None, optional
            Group label for target sequences.  When ``None`` (default) the
            groups are auto-detected.
        merge : bool, optional
            Whether to merge consecutive co-linear k-mer runs into single
            alignment blocks.  Default is ``True``.

        Raises
        ------
        ValueError
            If group auto-detection fails (≠2 groups, no explicit params), or
            if only one of *query_group* / *target_group* is supplied.
        KeyError
            If an explicit group label is not present in the index.
        """
        if query_group is None and target_group is None:
            pairs = self._get_default_group_pairs()
        elif (query_group is None) ^ (target_group is None):
            raise ValueError(
                'Provide both query_group and target_group, or neither.'
            )
        else:
            if query_group not in self._groups:
                raise KeyError(f"Group {query_group!r} not found.")
            if target_group not in self._groups:
                raise KeyError(f"Group {target_group!r} not found.")
            pairs = [(query_group, target_group)]

        for qg, tg in pairs:
            q_seqs = self._groups.get(qg, [])
            t_seqs = self._groups.get(tg, [])
            _log.info(
                "CrossIndex.compute_matches: computing matches between "
                "group %r (%d sequence(s)) and group %r (%d sequence(s))",
                qg,
                len(q_seqs),
                tg,
                len(t_seqs),
            )
            pair_records: list[PafRecord] = []
            for q_orig in q_seqs:
                for t_orig in t_seqs:
                    q_int = self._make_internal(qg, q_orig)
                    t_int = self._make_internal(tg, t_orig)
                    _log.debug(
                        "CrossIndex.compute_matches: comparing %r (group %r) "
                        "vs %r (group %r)",
                        q_orig,
                        qg,
                        t_orig,
                        tg,
                    )
                    lines = self._index.get_paf(q_int, t_int, merge)
                    for line in lines:
                        fields = line.split('\t')
                        fields[0] = q_orig
                        fields[5] = t_orig
                        pair_records.append(PafRecord.from_line('\t'.join(fields)))
            self._records_by_pair[(qg, tg)] = pair_records
            _log.info(
                "CrossIndex.compute_matches: stored %d record(s) for pair (%r, %r)",
                len(pair_records),
                qg,
                tg,
            )

    def get_paf(
        self,
        group_pairs: list[tuple[str, str]] | None = None,
        merge: bool = True,
    ) -> list[str]:
        """Return PAF lines for cross-group sequence comparisons.

        Parameters
        ----------
        group_pairs : list of (str, str) or None, optional
            Explicit list of ``(query_group, target_group)`` pairs to compare.
            If ``None`` (default):

            * 2 groups → the single cross-group pair.
            * 3+ groups → all non-self ordered pairs.
        merge : bool, optional
            Whether to merge consecutive co-linear k-mer runs before
            generating PAF lines.  Default is ``True``.

        Returns
        -------
        list[str]
            PAF-formatted lines (12 tab-separated columns each).
        """
        if group_pairs is None:
            group_pairs = self._get_default_group_pairs()

        paf_lines: list[str] = []
        for query_group, target_group in group_pairs:
            q_seqs = self._groups.get(query_group, [])
            t_seqs = self._groups.get(target_group, [])
            _log.info(
                'CrossIndex.get_paf: on-demand computation of %d x %d alignments '
                'between group %r and group %r '
                '(tip: call compute_matches() first to pre-cache results)',
                len(q_seqs),
                len(t_seqs),
                query_group,
                target_group,
            )
            for q_orig in q_seqs:
                for t_orig in t_seqs:
                    q_int = self._make_internal(query_group, q_orig)
                    t_int = self._make_internal(target_group, t_orig)
                    _log.debug(
                        "CrossIndex.get_paf: comparing %r (group %r) vs %r (group %r)",
                        q_orig,
                        query_group,
                        t_orig,
                        target_group,
                    )
                    lines = self._index.get_paf(q_int, t_int, merge)
                    for line in lines:
                        fields = line.split('\t')
                        fields[0] = q_orig
                        fields[5] = t_orig
                        paf_lines.append('\t'.join(fields))
        return paf_lines

    def run_merge(
        self,
        group_pairs: list[tuple[str, str]] | None = None,
    ) -> None:
        """Compute merged alignments and store the result as :attr:`_paf_records`.

        .. deprecated::
            Use :meth:`compute_matches` instead.  ``run_merge`` now delegates
            to ``compute_matches`` and is retained only for backward
            compatibility.

        Parameters
        ----------
        group_pairs : list of (str, str) or None, optional
            Group pairs to compare (same semantics as :meth:`compute_matches`).
            Defaults to all cross-group pairs.
        """
        if group_pairs is None:
            self.compute_matches(merge=True)
        else:
            for qg, tg in group_pairs:
                self.compute_matches(query_group=qg, target_group=tg, merge=True)

    # ------------------------------------------------------------------
    # Backward-compatible API (two-group a/b model)
    # ------------------------------------------------------------------

    def get_paf_all(self, merge: bool = True) -> list[str]:
        """Return PAF lines for all cross-group comparisons.

        Backward-compatible wrapper around :meth:`get_paf`.  When a group
        ``'b'`` is present, computes ``a`` vs ``b`` alignments; otherwise
        performs all-vs-all within group ``'a'``.

        Parameters
        ----------
        merge : bool, optional
            Whether to merge consecutive co-linear k-mer runs.
            Default is ``True``.

        Returns
        -------
        list[str]
            PAF-formatted lines.
        """
        if 'b' in self._groups and self._groups['b']:
            return self.get_paf(group_pairs=[('a', 'b')], merge=merge)
        # Single group or no group 'b': all-vs-all within group 'a'
        names_a = self._groups.get('a', [])
        _log.info(
            'CrossIndex: computing all-vs-all pairwise alignments '
            'within group a (%d sequences)',
            len(names_a),
        )
        paf_lines: list[str] = []
        for i, q_orig in enumerate(names_a):
            for j, t_orig in enumerate(names_a):
                if i == j:
                    continue
                q_int = self._make_internal('a', q_orig)
                t_int = self._make_internal('a', t_orig)
                lines = self._index.get_paf(q_int, t_int, merge)
                for line in lines:
                    fields = line.split('\t')
                    fields[0] = q_orig
                    fields[5] = t_orig
                    paf_lines.append('\t'.join(fields))
        return paf_lines

    def reorder_contigs(
        self,
        query_names: list[str] | None = None,
        target_names: list[str] | None = None,
        query_group: str | None = None,
        target_group: str | None = None,
    ) -> tuple[list[str], list[str]]:
        """Sort contigs for maximum collinearity.

        .. note::
            :meth:`compute_matches` must be called for the relevant group pair
            before calling this method.

        When *query_group* and *target_group* are not provided the method
        auto-detects the two groups to compare:

        * If there are **exactly two groups**, those two groups are used
          (regardless of their labels) and an info-level log message records
          which groups were selected and the order of the returned tuple.
        * Otherwise a :exc:`ValueError` is raised and the caller must supply
          explicit group labels via *query_group* / *target_group*.

        Parameters
        ----------
        query_names : list[str] or None, optional
            Explicit un-prefixed names within *query_group* to reorder.
            Defaults to all sequences in *query_group*.
        target_names : list[str] or None, optional
            Explicit un-prefixed names within *target_group* to reorder.
            Defaults to all sequences in *target_group*.
        query_group : str or None, optional
            Group label for the query (first element of the returned tuple).
            When ``None`` the group is auto-detected (requires exactly two
            groups).
        target_group : str or None, optional
            Group label for the target (second element of the returned tuple).
            When ``None`` the group is auto-detected (requires exactly two
            groups).

        Returns
        -------
        tuple[list[str], list[str]]
            ``(sorted_query_names, sorted_target_names)`` — both using
            original un-prefixed names.  The log output names the groups
            in the same order as the tuple elements.

        Raises
        ------
        ValueError
            If groups cannot be auto-detected (i.e. there are not exactly two
            groups and no explicit group labels were supplied), if only one of
            *query_group* / *target_group* is given, or if
            :meth:`compute_matches` has not been called for the resolved group
            pair.
        KeyError
            If an explicitly supplied group label is not present.
        """
        groups = list(self._groups.keys())

        if query_group is None and target_group is None:
            if len(groups) == 2:
                query_group, target_group = groups[0], groups[1]
                _log.info(
                    'CrossIndex.reorder_contigs: auto-selected groups '
                    '%r (query / first) and %r (target / second)',
                    query_group,
                    target_group,
                )
            else:
                raise ValueError(
                    'reorder_contigs requires exactly two groups when query_group '
                    'and target_group are not specified; '
                    f'found {len(groups)} group(s): {groups!r}. '
                    'Provide query_group and target_group explicitly, or use '
                    'reorder_for_colinearity for full control.'
                )
        elif (query_group is None) ^ (target_group is None):
            raise ValueError(
                'Provide both query_group and target_group, or neither.'
            )
        else:
            _log.info(
                'CrossIndex.reorder_contigs: using groups '
                '%r (query / first) and %r (target / second)',
                query_group,
                target_group,
            )

        pair = (query_group, target_group)
        if pair not in self._records_by_pair:
            raise ValueError(
                f"No matches computed for group pair {pair!r}. "
                "Call compute_matches() for this pair first."
            )

        q_names = (
            query_names
            if query_names is not None
            else list(self._groups[query_group])
        )
        t_names = (
            target_names
            if target_names is not None
            else list(self._groups[target_group])
        )
        q_internal = [self._make_internal(query_group, n) for n in q_names]
        t_internal = [self._make_internal(target_group, n) for n in t_names]
        sorted_q_int, sorted_t_int = self._index.optimal_contig_order(
            q_internal, t_internal
        )
        sorted_q = [self._split_internal(n)[1] for n in sorted_q_int]
        sorted_t = [self._split_internal(n)[1] for n in sorted_t_int]
        return sorted_q, sorted_t

    # ------------------------------------------------------------------
    # Dunder methods
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        """Return a concise machine-readable representation.

        Returns
        -------
        str
            ``CrossIndex(k=<k>, groups={<label>=<n>, ...})``.
        """
        group_info = ', '.join(f'{g}={len(names)}' for g, names in self._groups.items())
        return f'CrossIndex(k={self._k}, groups={{{group_info}}})'

    def __str__(self) -> str:
        """Return a human-readable stats summary.

        Returns
        -------
        str
            Multi-line summary of groups, sequence counts, computed pairs,
            and cached PAF record count.
        """
        n_total = sum(len(v) for v in self._groups.values())
        lines = [f'CrossIndex (k={self._k})']
        lines.append(f'  Total sequences : {n_total}')
        for g, names in self._groups.items():
            lines.append(f'  Group {g!r:12s}: {len(names):>6d} sequences')
        if self._records_by_pair:
            for (qg, tg), recs in self._records_by_pair.items():
                lines.append(
                    f'  Computed pair   : ({qg!r}, {tg!r}) → {len(recs)} record(s)'
                )
        else:
            lines.append('  Computed pairs  : none (call compute_matches() first)')
        lines.append(f'  PAF records     : {len(self._paf_records)}')
        return '\n'.join(lines)
