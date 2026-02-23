"""PAF (Pairwise mApping Format) file I/O with CIGAR string support.

This module provides:
- :class:`PafRecord` — a dataclass representing one PAF alignment line.
- :func:`parse_paf_file` — a generator that yields :class:`PafRecord` objects.
- :class:`PafAlignment` — a container that loads PAF records from a file or
  a list and provides a :meth:`~PafAlignment.reorder_contigs` method to
  reorganise sequences for maximum collinearity.

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

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Generator, Iterable


# ---------------------------------------------------------------------------
# CIGAR parsing
# ---------------------------------------------------------------------------

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


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
        ops.get("M", 0)
        + ops.get("D", 0)
        + ops.get("N", 0)
        + ops.get("=", 0)
        + ops.get("X", 0)
    )
    n_matches = ops.get("=", residue_matches)
    n_mismatches = ops.get("X", 0)
    n_gap_bases = ops.get("I", 0) + ops.get("D", 0)
    # Count distinct gap-opening events (each I or D run counts as one event).
    n_gaps = sum(
        1
        for _, op in _CIGAR_RE.findall(cigar)
        if op in ("I", "D")
    )
    return {
        "alignment_length": alignment_length,
        "n_matches": n_matches,
        "n_mismatches": n_mismatches,
        "n_gaps": n_gaps,
        "n_gap_bases": n_gap_bases,
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
    def from_line(cls, line: str) -> "PafRecord":
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
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12:
            raise ValueError(
                f"PAF line has {len(fields)} fields; expected at least 12: {line!r}"
            )
        tags: dict[str, Any] = {}
        cigar: str | None = None
        for tag_field in fields[12:]:
            parts = tag_field.split(":", 2)
            if len(parts) == 3:
                tag_name, tag_type, tag_value = parts
                if tag_type == "i":
                    tags[tag_name] = int(tag_value)
                elif tag_type == "f":
                    tags[tag_name] = float(tag_value)
                else:
                    tags[tag_name] = tag_value
                if tag_name == "cg" and tag_type == "Z":
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
            alignment_length=stats.get("alignment_length"),
            n_matches=stats.get("n_matches"),
            n_mismatches=stats.get("n_mismatches"),
            n_gaps=stats.get("n_gaps"),
            n_gap_bases=stats.get("n_gap_bases"),
        )

    def to_line(self) -> str:
        """Serialise this record back to a PAF-format string (no trailing newline).

        Returns
        -------
        str
            Tab-separated PAF line with the 12 required columns.  Optional
            tags are not included.
        """
        return "\t".join(
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
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
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
    q_weight: dict[str, float] = {n: 0.0 for n in query_names}
    q_wpos: dict[str, float] = {n: 0.0 for n in query_names}
    t_weight: dict[str, float] = {n: 0.0 for n in target_names}
    t_wpos: dict[str, float] = {n: 0.0 for n in target_names}

    for rec in all_records:
        if rec.query_name not in query_set or rec.target_name not in target_set:
            continue
        size = float(rec.alignment_block_len or (rec.query_end - rec.query_start))
        if size <= 0:
            continue

        # Target gravity from query's perspective.
        t_mid = t_offsets.get(rec.target_name, 0) + (rec.target_start + rec.target_end) / 2.0
        q_weight[rec.query_name] += size
        q_wpos[rec.query_name] += size * t_mid

        # Query gravity from target's perspective.
        q_mid = q_offsets_real.get(rec.query_name, 0) + (rec.query_start + rec.query_end) / 2.0
        t_weight[rec.target_name] += size
        t_wpos[rec.target_name] += size * q_mid

    def _gravity(name: str, wt: dict, wp: dict, total: float) -> float:
        w = wt.get(name, 0.0)
        return (wp.get(name, 0.0) / w / total) if w > 0 else float("inf")

    sorted_q = sorted(
        query_names,
        key=lambda n: _gravity(n, q_weight, q_wpos, total_target_len),
    )
    sorted_t = sorted(
        target_names,
        key=lambda n: _gravity(n, t_weight, t_wpos, total_query_len),
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

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_file(cls, path: str | Path) -> "PafAlignment":
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
    def from_records(cls, records: Iterable[PafRecord]) -> "PafAlignment":
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
            f"PafAlignment(records={len(self.records)}, "
            f"queries={len(self.query_names)}, "
            f"targets={len(self.target_names)})"
        )

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def filter_by_query(self, names: Iterable[str]) -> "PafAlignment":
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

    def filter_by_target(self, names: Iterable[str]) -> "PafAlignment":
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

    # ------------------------------------------------------------------
    # Contig reordering
    # ------------------------------------------------------------------

    def reorder_contigs(
        self,
        query_names: list[str] | None = None,
        target_names: list[str] | None = None,
    ) -> tuple[list[str], list[str]]:
        """Sort query and target contigs to maximise collinearity in the dotplot.

        Uses the gravity-centre algorithm: each contig is assigned a gravity
        equal to the weighted mean position of its alignment blocks on the
        opposing axis.  Contigs are then sorted by ascending gravity.

        Parameters
        ----------
        query_names : list[str] or None, optional
            Query contigs to reorder.  Defaults to :attr:`query_names`.
        target_names : list[str] or None, optional
            Target contigs to reorder.  Defaults to :attr:`target_names`.

        Returns
        -------
        tuple[list[str], list[str]]
            ``(sorted_query_names, sorted_target_names)``.
        """
        q = query_names if query_names is not None else self.query_names
        t = target_names if target_names is not None else self.target_names
        return compute_gravity_contigs(self.records, q, t)
