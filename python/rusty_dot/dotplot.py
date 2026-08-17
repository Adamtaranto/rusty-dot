"""
Dotplot visualization module for rusty-dot.

Provides the DotPlotter class for generating all-vs-all dotplots from
DNA sequence comparison data.

Reference: https://github.com/rrwick/Autocycler/blob/b0523350898faac71686251ec58f7d83bc2b1c28/src/dotplot.rs
"""

from __future__ import annotations

import dataclasses
import logging
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from matplotlib.collections import LineCollection, PatchCollection
import matplotlib.colors as mcolors
import matplotlib.figure
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.paf_io import CrossIndex, PafAlignment

if TYPE_CHECKING:
    from rusty_dot.annotation import GffAnnotation
    from rusty_dot.paf_io import CrossIndex

_log = logging.getLogger(__name__)


def _resolve_rasterized(
    n_segments: int,
    rasterized: Union[bool, str],
    threshold: int,
) -> bool:
    """Decide whether a match layer should be rasterised.

    Parameters
    ----------
    n_segments : int
        Number of line segments about to be drawn in this layer.
    rasterized : bool or str
        ``True`` / ``False`` force the choice; ``'auto'`` rasterises only when
        *n_segments* exceeds *threshold* (keeping small plots true vector while
        bounding the size of very dense ones).
    threshold : int
        Segment count above which ``'auto'`` switches to rasterised.

    Returns
    -------
    bool
        Whether to set ``rasterized=True`` on the layer.

    Raises
    ------
    ValueError
        If *rasterized* is a string other than ``'auto'``.
    """
    if isinstance(rasterized, bool):
        return rasterized
    if rasterized == 'auto':
        return n_segments > threshold
    raise ValueError(f"rasterized must be True, False, or 'auto', got {rasterized!r}")


def _chain_blocks(
    blocks: list[tuple[int, int, int, int, str]],
    chain_gap: int,
) -> list[tuple[int, int, int, int, str]]:
    """Chain co-linear match blocks separated by small gaps into single blocks.

    Blocks that lie on the same diagonal (forward strand) or anti-diagonal
    (reverse-complement strand) and whose query coordinates are separated by no
    more than *chain_gap* base pairs are merged into one block.  This collapses
    diagonals that exact k-mer matching has broken into many pieces (e.g. by
    SNPs) into a small number of long lines, greatly reducing the number of
    segments that must be drawn.

    Indels shift the diagonal, so blocks either side of an indel stay separate
    (a deliberate v1 choice — only strictly co-linear runs are joined).

    Parameters
    ----------
    blocks : list of tuple
        ``(query_start, query_end, target_start, target_end, strand)`` tuples,
        with ``target_start < target_end`` and ``strand`` one of ``'+'`` / ``'-'``.
    chain_gap : int
        Maximum gap (bp, measured on the query axis) between consecutive blocks
        on the same diagonal that will still be joined.  ``0`` disables chaining
        and returns *blocks* unchanged.

    Returns
    -------
    list of tuple
        The chained blocks, in the same 5-tuple format.
    """
    if chain_gap <= 0 or len(blocks) <= 1:
        return blocks

    arr = np.asarray(
        [(qs, qe, ts, te) for (qs, qe, ts, te, _s) in blocks], dtype=np.int64
    )
    strands = np.array([s for (_qs, _qe, _ts, _te, s) in blocks])
    qs, qe, ts, te = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]

    # Diagonal key: forward strand is constant t_start - q_start; reverse strand
    # is constant q_start + t_end (as query advances, target retreats).
    is_rev = strands == '-'
    diag = np.where(is_rev, qs + te, ts - qs)
    strand_code = is_rev.astype(np.int64)

    # Sort by (strand, diagonal, query_start) so each co-linear run is contiguous.
    order = np.lexsort((qs, diag, strand_code))
    qs_o, qe_o, ts_o, te_o = qs[order], qe[order], ts[order], te[order]
    diag_o, strand_o = diag[order], strand_code[order]

    # Index the (strand, diagonal) groups, then shift every coordinate by
    # group_index * BASE so distinct groups occupy non-overlapping value ranges
    # separated by more than chain_gap.  A single global sorted-interval merge is
    # then correct both within a group and across group boundaries (no running
    # state leaks between groups).
    grp_change = np.ones(len(qs_o), dtype=bool)
    grp_change[1:] = (strand_o[1:] != strand_o[:-1]) | (diag_o[1:] != diag_o[:-1])
    base_group = np.cumsum(grp_change) - 1
    span = int(max(qe_o.max(), qs_o.max())) + chain_gap + 1
    qs_shift = qs_o + base_group * span
    qe_shift = qe_o + base_group * span

    # Sorted-interval merge: a new chain starts at index 0 or wherever the next
    # start lies more than chain_gap beyond the running max end so far.
    running_end = np.maximum.accumulate(qe_shift)
    is_new = np.ones(len(qs_o), dtype=bool)
    is_new[1:] = qs_shift[1:] > running_end[:-1] + chain_gap
    starts = np.flatnonzero(is_new)

    # Reduce each chain: query span = [min qs, max qe]; target span = [min ts, max te].
    chained_qs = qs_o[starts]
    chained_qe = np.maximum.reduceat(qe_o, starts)
    chained_ts = np.minimum.reduceat(ts_o, starts)
    chained_te = np.maximum.reduceat(te_o, starts)
    chained_strand = strand_o[starts]

    result: list[tuple[int, int, int, int, str]] = []
    for i in range(len(starts)):
        strand = '-' if chained_strand[i] == 1 else '+'
        result.append(
            (
                int(chained_qs[i]),
                int(chained_qe[i]),
                int(chained_ts[i]),
                int(chained_te[i]),
                strand,
            )
        )
    return result


class DotPlotter:
    """Generate all-vs-all dotplots for sets of DNA sequences.

    Accepts a :class:`~rusty_dot.SequenceIndex` (single sequence collection),
    a :class:`~rusty_dot.paf_io.CrossIndex` (multi-group collection), or a
    :class:`~rusty_dot.paf_io.PafAlignment` loaded from an external aligner
    such as minimap2.

    When a ``PafAlignment`` is passed as *index*, sequence lengths are read
    from the PAF records and alignments are rendered directly — no k-mer index
    is required::

        from rusty_dot.paf_io import PafAlignment
        from rusty_dot.dotplot import DotPlotter

        aln = PafAlignment.from_file("alignments.paf")
        q_order, t_order = aln.reorder_contigs()

        plotter = DotPlotter(aln)
        plotter.plot(
            query_names=q_order,
            target_names=t_order,
            output_path="dotplot.png",
        )

    When using a ``CrossIndex``, use *query_group* and *target_group* in
    :meth:`plot` and :meth:`plot_single` so that sequence names are resolved
    automatically and pre-computed merged alignments are used for rendering::

        cross = CrossIndex(k=15)
        cross.load_fasta("assembly_a.fasta", group="a")
        cross.load_fasta("assembly_b.fasta", group="b")
        cross.compute_matches()  # pre-compute merged alignments

        plotter = DotPlotter(cross)
        plotter.plot(
            query_group="a",   # sequence names looked up from group 'a'
            target_group="b",  # sequence names looked up from group 'b'
            output_path="cross_plot.png",
        )

    To colour alignments by sequence identity, supply a
    :class:`~rusty_dot.paf_io.PafAlignment` and set
    ``color_by_identity=True``::

        from rusty_dot.paf_io import PafAlignment
        aln = PafAlignment.from_file("alignments.paf")
        plotter = DotPlotter(aln)
        fig = plotter.plot(color_by_identity=True, identity_palette="viridis")
        cbar = plotter.plot_identity_colorbar(palette="viridis")

    Parameters
    ----------
    index : SequenceIndex, CrossIndex, or PafAlignment
        A populated index or alignment collection.  When a
        :class:`~rusty_dot.paf_io.PafAlignment` is supplied, it is used both
        to resolve sequence lengths and as the source of alignment segments.
    paf_alignment : PafAlignment, optional
        Pre-loaded PAF alignments used as the data source when
        ``color_by_identity=True`` and *index* is a ``SequenceIndex`` or
        ``CrossIndex``.  When *index* is already a ``PafAlignment`` this
        argument is ignored.  When ``None`` (default) and *index* is not a
        ``PafAlignment``, k-mer matches from *index* are used for plotting.

    Examples
    --------
    >>> from rusty_dot import SequenceIndex
    >>> from rusty_dot.dotplot import DotPlotter
    >>> idx = SequenceIndex(k=10)
    >>> idx.add_sequence("seq1", "ACGTACGTACGT" * 10)
    >>> idx.add_sequence("seq2", "TACGTACGTACG" * 10)
    >>> plotter = DotPlotter(idx)
    >>> fig = plotter.plot(output_path="dotplot.png")  # save to file
    >>> fig = plotter.plot()  # display inline in Jupyter, no file saved
    """

    def __init__(
        self,
        index: Union[SequenceIndex, 'CrossIndex', 'PafAlignment'],
        paf_alignment: Optional['PafAlignment'] = None,
    ) -> None:
        """Initialise the DotPlotter.

        Parameters
        ----------
        index : SequenceIndex, CrossIndex, or PafAlignment
            A populated index or alignment collection.  When a
            :class:`~rusty_dot.paf_io.PafAlignment` is supplied, it is used
            both to resolve sequence lengths and as the source of alignment
            segments.
        paf_alignment : PafAlignment, optional
            Pre-loaded PAF alignments.  Used for identity-based colouring
            when *index* is a ``SequenceIndex`` or ``CrossIndex``.  When
            *index* is already a ``PafAlignment`` this argument is ignored.
            When ``None`` (default), k-mer matches from *index* are used.
        """
        self.index = index
        # When a PafAlignment is passed as the primary index, use it for
        # rendering alignment segments (the explicit paf_alignment kwarg is
        # then redundant and is ignored to avoid confusion).
        if isinstance(index, PafAlignment):
            self.paf_alignment: Optional[PafAlignment] = index
        else:
            self.paf_alignment = paf_alignment

    def _index_is_paf(self) -> bool:
        """Return ``True`` when *index* is a :class:`~rusty_dot.paf_io.PafAlignment`.

        Helper used by :meth:`_plot_panel` to decide whether to draw from PAF
        records or from the k-mer engine.
        """
        return isinstance(self.index, PafAlignment)

    def _index_is_cross(self) -> bool:
        """Return ``True`` when *index* is a :class:`~rusty_dot.paf_io.CrossIndex`."""
        from rusty_dot.paf_io import CrossIndex

        return isinstance(self.index, CrossIndex)

    @staticmethod
    def _strip_group_prefix(name: str) -> str:
        """Strip the ``'group:'`` prefix from a CrossIndex internal name.

        For a ``SequenceIndex`` or ``PafAlignment`` name that contains no
        ``':'``, the name is returned unchanged.

        Parameters
        ----------
        name : str
            Sequence name, possibly in ``'group:name'`` format.

        Returns
        -------
        str
            The un-prefixed name.
        """
        _, sep, suffix = name.partition(':')
        return suffix if sep else name

    def _get_paf_override(
        self,
        query_group: Optional[str],
        target_group: Optional[str],
    ) -> Optional['PafAlignment']:
        """Return a :class:`PafAlignment` built from pre-computed CrossIndex records.

        Parameters
        ----------
        query_group : str or None
        target_group : str or None

        Returns
        -------
        PafAlignment or None
            ``None`` when no pre-computed records exist for the pair or when
            either group is ``None``.
        """
        if query_group is None or target_group is None or not self._index_is_cross():
            return None
        cross = self.index  # type: ignore[assignment]
        pair = (query_group, target_group)
        if pair in cross.computed_group_pairs:
            records = cross.get_records_for_pair(query_group, target_group)
            paf = PafAlignment(records)
            _log.debug(
                'DotPlotter: using %d pre-computed record(s) from pair %r',
                len(records),
                pair,
            )
            return paf
        _log.debug(
            'DotPlotter: no pre-computed records for pair %r; '
            'will compute k-mer matches on demand',
            pair,
        )
        return None

    def _resolve_group_names(
        self,
        query_group: Optional[str],
        target_group: Optional[str],
        query_names: Optional[list[str]],
        target_names: Optional[list[str]],
    ) -> tuple[Optional[list[str]], Optional[list[str]], Optional['PafAlignment']]:
        """Resolve sequence name lists and an optional cached PAF alignment.

        When *query_group* / *target_group* are provided and *index* is a
        :class:`~rusty_dot.paf_io.CrossIndex`:

        * Query and target name lists are populated from the group's internal
          (``'group:name'``) identifiers, overriding any explicitly provided
          *query_names* / *target_names*.
        * If the group pair has pre-computed matches (via
          :meth:`~rusty_dot.paf_io.CrossIndex.compute_matches`), a
          :class:`~rusty_dot.paf_io.PafAlignment` is built from those records
          and returned so that :meth:`_plot_panel` can use them directly.

        Parameters
        ----------
        query_group : str or None
            Group label for query sequences.
        target_group : str or None
            Group label for target sequences.
        query_names : list[str] or None
            Caller-supplied query names (used when groups are not provided).
        target_names : list[str] or None
            Caller-supplied target names.

        Returns
        -------
        tuple of (list[str] | None, list[str] | None, PafAlignment | None)
            Resolved ``(query_names, target_names, paf_override)``.
            *paf_override* is ``None`` when no pre-computed records are found.

        Raises
        ------
        ValueError
            If *query_group* / *target_group* are provided but *index* is not
            a ``CrossIndex``.
        """
        if query_group is None and target_group is None:
            return query_names, target_names, None

        if not self._index_is_cross():
            raise ValueError(
                'query_group and target_group can only be used when index is '
                'a CrossIndex.'
            )

        cross = self.index  # type: ignore[assignment]

        if query_group is not None:
            query_names = cross.sequence_names(group=query_group)
        if target_group is not None:
            target_names = cross.sequence_names(group=target_group)

        paf_override = self._get_paf_override(query_group, target_group)
        return query_names, target_names, paf_override

    def plot(
        self,
        query_names: Optional[list[str]] = None,
        target_names: Optional[list[str]] = None,
        query_group: Optional[str] = None,
        target_group: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize_per_panel: float = 4.0,
        dot_size: float = 0.5,
        dot_color: str = 'blue',
        rc_color: str = 'red',
        merge: bool = True,
        title: Optional[str] = None,
        dpi: int = 150,
        scale_sequences: bool = True,
        format: Optional[str] = None,
        min_length: int = 0,
        color_by_identity: bool = False,
        identity_palette: str = 'viridis',
        annotation: Optional['GffAnnotation'] = None,
        chain_gap: int = 0,
        rasterized: Union[bool, str] = 'auto',
        rasterization_threshold: int = 50_000,
        reverse_contigs: Optional[set[str]] = None,
        hide_internal_axes: bool = False,
    ) -> matplotlib.figure.Figure:
        """Plot an all-vs-all dotplot grid.

        If both ``query_names`` and ``target_names`` are provided, the plot
        will show each query sequence (rows) against each target sequence
        (columns). If only one set is provided, or neither, all pairwise
        combinations within the available sequences are plotted.

        When *index* is a :class:`~rusty_dot.paf_io.CrossIndex`, use
        *query_group* and *target_group* to specify which groups supply the
        query and target sequences.  The corresponding internal
        (``'group:name'``) identifiers are looked up automatically and used
        for sequence-length queries and k-mer comparisons.  If
        :meth:`~rusty_dot.paf_io.CrossIndex.compute_matches` has already been
        called for that pair, the pre-computed merged alignments are used for
        rendering rather than recomputing on the fly.

        The figure is always returned so it can be displayed inline in a
        Jupyter notebook.  When ``output_path`` is provided the figure is
        also saved to disk.

        Parameters
        ----------
        query_names : list[str], optional
            Sequence names for the y-axis (rows). If ``None``, uses all
            sequences in the index.  Ignored when *query_group* is provided
            and *index* is a ``CrossIndex``.
        target_names : list[str], optional
            Sequence names for the x-axis (columns). If ``None``, uses all
            sequences in the index.  Ignored when *target_group* is provided
            and *index* is a ``CrossIndex``.
        query_group : str or None, optional
            Group label whose sequences are used as query (rows).  When
            provided, the group's sequences are looked up from *index*
            (which must be a ``CrossIndex``) and *query_names* is ignored.
        target_group : str or None, optional
            Group label whose sequences are used as target (columns).  When
            provided, the group's sequences are looked up from *index*
            (which must be a ``CrossIndex``) and *target_names* is ignored.
        output_path : str or Path, optional
            Output image file path.  When ``None`` (default) the figure is
            not saved to disk.  Use a ``.svg`` extension (or set
            ``format='svg'``) to produce an SVG vector image.
        figsize_per_panel : float, optional
            Base size in inches for each subplot panel when
            ``scale_sequences=False``.  When ``scale_sequences=True`` this
            value sets the size of the *longest* sequence axis and all
            other axes are scaled proportionally.  Default is ``4.0``.
        dot_size : float, optional
            Size of each dot in the scatter plot. Default is ``0.5``.
        dot_color : str, optional
            Colour for forward-strand (``+``) match lines. Default is ``"blue"``.
        rc_color : str, optional
            Colour for reverse-complement (``-``) strand match lines.
            Default is ``"red"``.
        merge : bool, optional
            Whether to merge sequential k-mer runs before plotting.
            Default is ``True``.
        title : str, optional
            Overall figure title. If ``None``, no title is added.
        dpi : int, optional
            Resolution of the output image. Default is ``150``.  For vector
            formats this only affects any rasterised match layer (see
            *rasterized*); axes and labels remain resolution-independent.  Raise
            it (e.g. ``dpi=300``) for a higher-resolution raster (PNG) figure.
        scale_sequences : bool, optional
            When ``True`` (default), subplot widths and heights are
            proportional to the lengths of the corresponding sequences so that
            relative sequence sizes are preserved.  When ``False``, every
            panel has the same fixed size.
        format : str, optional
            Output image format (e.g. ``'png'``, ``'svg'``, ``'pdf'``).
            When ``None`` (default), the format is inferred from the
            ``output_path`` file extension.
        min_length : int, optional
            Minimum alignment length to display.  Matches shorter than this
            value are not drawn.  Applies to merged k-mer runs and pre-computed
            PAF alignments.  Default is ``0`` (no filtering).
        color_by_identity : bool, optional
            When ``True``, alignments are coloured by sequence identity using
            the *identity_palette* colormap.  Requires a
            :class:`~rusty_dot.paf_io.PafAlignment` to be supplied as
            ``paf_alignment`` to :meth:`__init__`; if no PAF alignment is
            available a warning is logged and the default strand colours are
            used instead.  Default is ``False``.
        identity_palette : str, optional
            Matplotlib colormap name used to map identity values (0–1) to
            colours when ``color_by_identity=True``.  Default is
            ``'viridis'``.
        annotation : GffAnnotation, optional
            Feature annotations to overlay on self-vs-self diagonal panels.
            Each feature is drawn as a coloured square at its genomic
            position.  Sequence names in *annotation* that are absent from
            the index emit a warning.  Default is ``None``.
        chain_gap : int, optional
            When greater than ``0``, co-linear match blocks on the same diagonal
            separated by up to *chain_gap* bp are chained into single lines
            before drawing, greatly reducing the number of segments (and thus
            render time and file size) for dense plots.  Default is ``0`` (off).
        rasterized : bool or str, optional
            Whether to rasterise the match layer.  ``'auto'`` (default) keeps it
            true vector — infinitely zoomable in SVG/PDF — until a panel's
            segment count exceeds *rasterization_threshold*, above which that
            layer is rasterised at *dpi* to bound file size.  ``True`` / ``False``
            force the behaviour.  Axes, ticks and labels always stay vector.
        rasterization_threshold : int, optional
            Segment count per strand/panel above which ``rasterized='auto'``
            rasterises the layer.  Default is ``50_000``.
        reverse_contigs : set[str] or None, optional
            Un-prefixed query (row) contig names to render reverse-complemented
            so reverse-oriented contigs read along the main diagonal.  When
            ``None`` (default) the set is pulled automatically from the index:
            :meth:`~rusty_dot.paf_io.CrossIndex.reversed_contigs` for the
            *query_group* of a ``CrossIndex``, or
            :attr:`~rusty_dot.paf_io.PafAlignment.reversed_contigs` for a
            ``PafAlignment`` (both populated by a prior ``reorder`` call).  Pass
            an explicit set (including ``set()`` to disable) to override.
        hide_internal_axes : bool, optional
            When ``True``, internal panel boundaries are removed so the grid
            reads as one continuous plot: inter-panel gaps collapse to zero,
            ticks and spines shared between adjacent panels are hidden, and
            only the outer frame with its tick labels remains.  Default is
            ``False``.

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.  In a Jupyter notebook the figure is
            displayed inline automatically; call ``matplotlib.pyplot.close``
            on the returned object when it is no longer needed.

        Raises
        ------
        ValueError
            If *query_group* / *target_group* are provided but *index* is
            not a ``CrossIndex``.
        """
        # Resolve group names and optional pre-computed PAF records.
        query_names, target_names, paf_override = self._resolve_group_names(
            query_group, target_group, query_names, target_names
        )

        all_names = self.index.sequence_names()
        if not all_names:
            raise ValueError('No sequences in the index.')

        if query_names is None:
            query_names = sorted(all_names)
        if target_names is None:
            target_names = sorted(all_names)

        # Use the per-call override if available, otherwise fall back to the
        # paf_alignment set at construction time.
        effective_paf = paf_override if paf_override is not None else self.paf_alignment

        # Resolve the set of reverse-oriented query contigs (un-prefixed names).
        # An explicit argument wins; otherwise auto-pull from the index.
        if reverse_contigs is not None:
            reverse_set = set(reverse_contigs)
        elif isinstance(self.index, CrossIndex) and query_group is not None:
            reverse_set = self.index.reversed_contigs(query_group)
        elif isinstance(self.index, PafAlignment):
            reverse_set = set(self.index.reversed_contigs)
        else:
            reverse_set = set()

        # Warn about annotation sequences missing from the index.
        if annotation is not None:
            index_seqs = set(all_names)
            for ann_seq in annotation.sequence_names():
                if ann_seq not in index_seqs:
                    _log.warning(
                        'Annotation contains features for sequence %r which is '
                        'not present in the index. These features will not be '
                        'plotted.',
                        ann_seq,
                    )

        nrows = len(query_names)
        ncols = len(target_names)

        flush_kw: dict[str, float] = (
            {'wspace': 0.0, 'hspace': 0.0} if hide_internal_axes else {}
        )
        if scale_sequences:
            q_lens = [self.index.get_sequence_length(n) for n in query_names]
            t_lens = [self.index.get_sequence_length(n) for n in target_names]
            max_len = max(max(q_lens), max(t_lens), 1)
            col_widths = [figsize_per_panel * (seq_len / max_len) for seq_len in t_lens]
            row_heights = [
                figsize_per_panel * (seq_len / max_len) for seq_len in q_lens
            ]
            fig_w = sum(col_widths)
            fig_h = sum(row_heights)
            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(fig_w, fig_h),
                squeeze=False,
                gridspec_kw={
                    'width_ratios': col_widths,
                    'height_ratios': row_heights,
                    **flush_kw,
                },
            )
        else:
            fig_w = figsize_per_panel * ncols
            fig_h = figsize_per_panel * nrows
            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(fig_w, fig_h),
                squeeze=False,
                gridspec_kw=flush_kw if hide_internal_axes else None,
            )

        for row_idx, q_name in enumerate(query_names):
            for col_idx, t_name in enumerate(target_names):
                ax = axes[row_idx][col_idx]
                self._plot_panel(
                    ax,
                    q_name,
                    t_name,
                    dot_size=dot_size,
                    dot_color=dot_color,
                    rc_color=rc_color,
                    merge=merge,
                    min_length=min_length,
                    # Sequence name labels: y-label on leftmost column only;
                    # column (x) labels are shown as titles on the top row.
                    show_xlabel=False,
                    show_ylabel=(col_idx == 0),
                    color_by_identity=color_by_identity,
                    identity_palette=identity_palette,
                    paf_alignment_override=effective_paf,
                    chain_gap=chain_gap,
                    rasterized=rasterized,
                    rasterization_threshold=rasterization_threshold,
                    reverse_query=self._strip_group_prefix(q_name) in reverse_set,
                )

                # Column label at top of each column (top row only), rotated.
                # Use the display name (strip group prefix for CrossIndex).
                if row_idx == 0:
                    ax.set_title(
                        self._strip_group_prefix(t_name),
                        fontsize=8,
                        rotation=45,
                        ha='left',
                        va='bottom',
                    )

                # Suppress redundant tick labels on internal panels.
                if row_idx < nrows - 1:
                    ax.tick_params(axis='x', labelbottom=False)
                if col_idx > 0:
                    ax.tick_params(axis='y', labelleft=False)

                # Remove internal ticks and spines so the grid reads as one
                # continuous plot, keeping the outer frame intact.
                if hide_internal_axes:
                    if row_idx < nrows - 1:
                        ax.tick_params(axis='x', bottom=False)
                        ax.spines['bottom'].set_visible(False)
                    if row_idx > 0:
                        ax.spines['top'].set_visible(False)
                    if col_idx > 0:
                        ax.tick_params(axis='y', left=False)
                        ax.spines['left'].set_visible(False)
                    if col_idx < ncols - 1:
                        ax.spines['right'].set_visible(False)

                # Annotation squares on self-vs-self (diagonal) panels.
                if annotation is not None and q_name == t_name:
                    self._draw_annotation_squares(ax, q_name, annotation)

        if title:
            fig.suptitle(title, fontsize=14, y=1.01)

        if hide_internal_axes:
            # tight_layout() would reinsert inter-panel gaps; keep the panels
            # flush and just leave margins for the outer labels and titles.
            fig.subplots_adjust(
                left=0.1, right=0.98, bottom=0.06, top=0.92, wspace=0.0, hspace=0.0
            )
        else:
            plt.tight_layout()
        if output_path is not None:
            plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight', format=format)
        return fig

    def _plot_panel(
        self,
        ax: plt.Axes,
        query_name: str,
        target_name: str,
        dot_size: float = 0.5,
        dot_color: str = 'blue',
        rc_color: str = 'red',
        merge: bool = True,
        min_length: int = 0,
        show_xlabel: bool = True,
        show_ylabel: bool = True,
        color_by_identity: bool = False,
        identity_palette: str = 'viridis',
        paf_alignment_override: Optional['PafAlignment'] = None,
        chain_gap: int = 0,
        rasterized: Union[bool, str] = 'auto',
        rasterization_threshold: int = 50_000,
        reverse_query: bool = False,
    ) -> None:
        """Render a single comparison panel onto the given Axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        query_name : str
            Name of the query sequence (y-axis).  For a
            :class:`~rusty_dot.paf_io.CrossIndex` this is the internal
            (``'group:name'``) identifier; the group prefix is stripped for
            axis labels and PAF record lookup.
        target_name : str
            Name of the target sequence (x-axis).  Same note as *query_name*.
        dot_size : float, optional
            Marker size. Default is ``0.5``.
        dot_color : str, optional
            Marker colour for forward-strand (``+``) matches. Default is ``"blue"``.
        rc_color : str, optional
            Marker colour for reverse-complement (``-``) matches. Default is ``"red"``.
        merge : bool, optional
            Whether to merge sequential runs. Default is ``True``.
        min_length : int, optional
            Minimum alignment length to display.  Matches shorter than this
            value are skipped.  Default is ``0`` (no filtering).
        show_xlabel : bool, optional
            Whether to render the target sequence name as an x-axis label.
            Default is ``True``.
        show_ylabel : bool, optional
            Whether to render the query sequence name as a y-axis label.
            Default is ``True``.
        color_by_identity : bool, optional
            When ``True``, colour alignments by sequence identity using
            *identity_palette*.  Requires a PAF alignment (either
            ``paf_alignment_override`` or ``self.paf_alignment``) to be set;
            if not, a warning is logged and strand colours are used instead.
            Default is ``False``.
        identity_palette : str, optional
            Matplotlib colormap name for identity-based colouring.
            Default is ``'viridis'``.
        paf_alignment_override : PafAlignment or None, optional
            Pre-computed PAF alignments to use for this panel.  When
            provided, this takes precedence over ``self.paf_alignment`` for
            record lookup.  Typically supplied from pre-computed
            :class:`~rusty_dot.paf_io.CrossIndex` records.
            Default is ``None``.
        chain_gap : int, optional
            When greater than ``0``, co-linear match blocks on the same diagonal
            separated by up to *chain_gap* bp (on the query axis) are chained
            into a single line before drawing, reducing the segment count.
            Ignored for identity-coloured rendering.  Default is ``0`` (off).
        rasterized : bool or str, optional
            Controls whether the match layer is rasterised.  ``'auto'`` (default)
            keeps it true vector when the segment count is at or below
            *rasterization_threshold* and rasterises it otherwise; ``True`` /
            ``False`` force the choice.  Axes, ticks and labels always stay vector.
        rasterization_threshold : int, optional
            Segment count above which ``rasterized='auto'`` switches a layer to
            rasterised.  Default is ``50_000``.
        reverse_query : bool, optional
            When ``True`` the query (row) contig is rendered reverse-complemented:
            every match's query coordinates are mirrored (``q → q_len - q``) and
            its strand colour flipped, so a reverse-oriented contig reads along
            the main diagonal.  The underlying records are not modified.
            Default is ``False``.
        """
        q_len = self.index.get_sequence_length(query_name)
        t_len = self.index.get_sequence_length(target_name)

        # Display names: strip 'group:' prefix for CrossIndex internal names.
        display_q = self._strip_group_prefix(query_name)
        display_t = self._strip_group_prefix(target_name)

        # Effective PAF alignment: per-call override takes precedence.
        effective_paf = (
            paf_alignment_override
            if paf_alignment_override is not None
            else self.paf_alignment
        )

        # Determine rendering mode:
        # • When index is a PafAlignment (no k-mer index), always draw from
        #   PAF records using strand colours unless color_by_identity is set.
        # • When a paf_alignment_override is provided (CrossIndex pre-computed
        #   records), use it for rendering.
        # • When index has a k-mer engine but color_by_identity is requested
        #   without a PafAlignment, fall back to k-mer matches with a warning.
        use_paf = (
            color_by_identity
            or self._index_is_paf()
            or (paf_alignment_override is not None)
        )

        if use_paf and effective_paf is None:
            _log.warning(
                'color_by_identity=True requires a PafAlignment; k-mer matches '
                'are always 100% identity. Pass paf_alignment= to DotPlotter '
                'to enable identity colouring.'
            )
            use_paf = False

        if use_paf and color_by_identity:
            # Identity-coloured PAF records: one segment per record with a
            # per-segment colour.  Chaining is not applied here because each
            # record carries its own identity value.
            records = self._records_for_pair(effective_paf, display_q, display_t)
            if reverse_query:
                # Mirror the query coordinates and flip the strand so the
                # contig renders reverse-complemented (originals untouched).
                records = [
                    dataclasses.replace(
                        rec,
                        query_start=q_len - rec.query_end,
                        query_end=q_len - rec.query_start,
                        strand='-' if rec.strand == '+' else '+',
                    )
                    for rec in records
                ]
            self._draw_identity_records(
                ax,
                records,
                identity_palette=identity_palette,
                dot_size=dot_size,
                min_length=min_length,
                rasterized=rasterized,
                rasterization_threshold=rasterization_threshold,
            )
        else:
            # Strand-coloured rendering, from either pre-computed PAF records or
            # the k-mer engine.  Both produce (qs, qe, ts, te, strand) blocks.
            if use_paf:
                records = self._records_for_pair(effective_paf, display_q, display_t)
                blocks = [
                    (
                        rec.query_start,
                        rec.query_end,
                        rec.target_start,
                        rec.target_end,
                        rec.strand,
                    )
                    for rec in records
                ]
            else:
                blocks = [
                    (qs, qe, ts, te, strand)
                    for qs, qe, ts, te, strand in self.index.compare_sequences_stranded(
                        query_name, target_name, merge
                    )
                ]
            if reverse_query:
                # Mirror the query coordinates and flip the strand so the
                # contig renders reverse-complemented (originals untouched).
                blocks = [
                    (
                        q_len - qe,
                        q_len - qs,
                        ts,
                        te,
                        '-' if strand == '+' else '+',
                    )
                    for qs, qe, ts, te, strand in blocks
                ]
            if chain_gap > 0:
                blocks = _chain_blocks(blocks, chain_gap)
            self._draw_stranded_blocks(
                ax,
                blocks,
                dot_color=dot_color,
                rc_color=rc_color,
                dot_size=dot_size,
                min_length=min_length,
                rasterized=rasterized,
                rasterization_threshold=rasterization_threshold,
            )

        ax.set_xlim(0, t_len)
        ax.set_ylim(0, q_len)
        ax.invert_yaxis()
        if show_xlabel:
            ax.set_xlabel(display_t, fontsize=8)
        if show_ylabel:
            ax.set_ylabel(display_q, fontsize=8)
        ax.tick_params(axis='both', labelsize=6)
        ax.set_aspect('auto')

    @staticmethod
    def _draw_stranded_blocks(
        ax: plt.Axes,
        blocks: list[tuple[int, int, int, int, str]],
        dot_color: str,
        rc_color: str,
        dot_size: float,
        min_length: int,
        rasterized: Union[bool, str],
        rasterization_threshold: int,
    ) -> None:
        """Draw strand-coloured match blocks as vectorised LineCollections.

        Forward (``+``) blocks are drawn as diagonal segments and reverse (``-``)
        blocks as anti-diagonal segments, each strand batched into a single
        ``LineCollection`` built from a NumPy ``(N, 2, 2)`` array.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        blocks : list of tuple
            ``(query_start, query_end, target_start, target_end, strand)`` blocks.
        dot_color, rc_color : str
            Colours for forward and reverse-complement matches.
        dot_size : float
            Line width in points.
        min_length : int
            Skip blocks whose query length is below this (``0`` = keep all).
        rasterized : bool or str
            Passed to :func:`_resolve_rasterized` per strand collection.
        rasterization_threshold : int
            Segment-count threshold for ``rasterized='auto'``.
        """
        if not blocks:
            return
        arr = np.asarray(
            [(qs, qe, ts, te) for (qs, qe, ts, te, _s) in blocks], dtype=float
        )
        strand = np.array([s for (_qs, _qe, _ts, _te, s) in blocks])
        lengths = arr[:, 1] - arr[:, 0]
        keep = lengths >= min_length if min_length > 0 else np.ones(len(arr), bool)

        fwd = arr[keep & (strand != '-')]
        rev = arr[keep & (strand == '-')]

        if len(fwd):
            # Forward: (t_start, q_start) -> (t_end, q_end).
            fwd_seg = np.stack(
                [
                    np.column_stack([fwd[:, 2], fwd[:, 0]]),
                    np.column_stack([fwd[:, 3], fwd[:, 1]]),
                ],
                axis=1,
            )
            ax.add_collection(
                LineCollection(
                    fwd_seg,
                    colors=dot_color,
                    linewidths=dot_size,
                    alpha=0.7,
                    rasterized=_resolve_rasterized(
                        len(fwd_seg), rasterized, rasterization_threshold
                    ),
                )
            )
        if len(rev):
            # Reverse complement: (t_end, q_start) -> (t_start, q_end).
            rev_seg = np.stack(
                [
                    np.column_stack([rev[:, 3], rev[:, 0]]),
                    np.column_stack([rev[:, 2], rev[:, 1]]),
                ],
                axis=1,
            )
            ax.add_collection(
                LineCollection(
                    rev_seg,
                    colors=rc_color,
                    linewidths=dot_size,
                    alpha=0.7,
                    rasterized=_resolve_rasterized(
                        len(rev_seg), rasterized, rasterization_threshold
                    ),
                )
            )

    @staticmethod
    def _draw_identity_records(
        ax: plt.Axes,
        records: list,
        identity_palette: str,
        dot_size: float,
        min_length: int,
        rasterized: Union[bool, str],
        rasterization_threshold: int,
    ) -> None:
        """Draw PAF records coloured by per-record identity.

        Each record becomes one segment; colours come from the record identity
        (``residue_matches / alignment_block_len``) mapped through
        *identity_palette*.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to draw on.
        records : list of PafRecord
            Records for this panel's sequence pair.
        identity_palette : str
            Matplotlib colormap name for identity colouring.
        dot_size : float
            Line width in points.
        min_length : int
            Skip records whose query aligned length is below this.
        rasterized : bool or str
            Passed to :func:`_resolve_rasterized`.
        rasterization_threshold : int
            Segment-count threshold for ``rasterized='auto'``.
        """
        cmap = plt.get_cmap(identity_palette)
        norm = mcolors.Normalize(vmin=0, vmax=1)
        segments: list[list[tuple[float, float]]] = []
        colors: list = []
        for rec in records:
            if min_length > 0 and rec.query_aligned_len < min_length:
                continue
            identity = (
                rec.residue_matches / rec.alignment_block_len
                if rec.alignment_block_len > 0
                else 1.0
            )
            color = cmap(norm(identity))
            if rec.strand == '-':
                xs = (rec.target_end, rec.target_start)
            else:
                xs = (rec.target_start, rec.target_end)
            segments.append([(xs[0], rec.query_start), (xs[1], rec.query_end)])
            colors.append(color)
        if segments:
            ax.add_collection(
                LineCollection(
                    segments,
                    colors=colors,
                    linewidths=dot_size,
                    alpha=0.7,
                    rasterized=_resolve_rasterized(
                        len(segments), rasterized, rasterization_threshold
                    ),
                )
            )

    def _records_for_pair(
        self,
        effective_paf: 'PafAlignment',
        display_q: str,
        display_t: str,
    ) -> list:
        """Return the PAF records for one ``(query, target)`` pair.

        A ``(query_name, target_name) -> list[PafRecord]`` index is built once
        per unique ``PafAlignment`` object and cached on the plotter, so panel
        rendering does not rescan every record for each of the ``N × N`` panels.

        Parameters
        ----------
        effective_paf : PafAlignment
            The alignment whose records should be indexed.
        display_q : str
            Query sequence display name (group prefix already stripped).
        display_t : str
            Target sequence display name (group prefix already stripped).

        Returns
        -------
        list of PafRecord
            Records matching the requested pair (empty list if none).
        """
        cache = getattr(self, '_paf_record_index_cache', None)
        # Rebuild the index only when the PafAlignment object changes identity
        # (accessing ``.records`` can itself be an O(records) rebuild).
        if cache is None or cache[0] is not effective_paf:
            index: dict[tuple[str, str], list] = {}
            for r in effective_paf.records:
                index.setdefault((r.query_name, r.target_name), []).append(r)
            cache = (effective_paf, index)
            self._paf_record_index_cache = cache
        return cache[1].get((display_q, display_t), [])

    def _draw_annotation_squares(
        self,
        ax: plt.Axes,
        seq_name: str,
        annotation: 'GffAnnotation',
    ) -> None:
        """Overlay annotation feature squares on a self-vs-self panel.

        Each feature ``[start, end)`` is drawn as a filled square at position
        ``(start, start)`` to ``(end, end)`` in the dotplot coordinate system.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes of the self-vs-self panel.
        seq_name : str
            Sequence name whose features should be drawn.
        annotation : GffAnnotation
            The annotation object providing features and colours.
        """
        features = annotation.get_features_for_sequence(seq_name)
        # Batch all feature squares into a single PatchCollection rather than
        # adding one artist per feature.
        rects = []
        facecolors = []
        for feat in features:
            width = feat.end - feat.start
            rects.append(mpatches.Rectangle((feat.start, feat.start), width, width))
            facecolors.append(annotation.get_color(feat.feature_type))
        if rects:
            ax.add_collection(
                PatchCollection(
                    rects,
                    facecolors=facecolors,
                    edgecolors='none',
                    alpha=0.35,
                    match_original=False,
                )
            )

    def plot_annotation_legend(
        self,
        annotation: 'GffAnnotation',
        output_path: Optional[Union[str, Path]] = None,
        figsize: tuple[float, float] = (3.0, 4.0),
        dpi: int = 150,
        format: Optional[str] = None,
    ) -> matplotlib.figure.Figure:
        """Render the annotation feature-type legend as a standalone figure.

        Produces a figure containing only a colour legend that maps each
        feature type to its assigned colour.  This is intended to be
        displayed alongside dotplots produced with an *annotation* argument.

        Parameters
        ----------
        annotation : GffAnnotation
            The annotation object whose feature-type colours are displayed.
        output_path : str or Path, optional
            Output image file path.  When ``None`` (default) the figure is
            not saved to disk.
        figsize : tuple[float, float], optional
            Figure size as ``(width, height)`` in inches.
            Default is ``(3.0, 4.0)``.
        dpi : int, optional
            Output image resolution. Default is ``150``.
        format : str, optional
            Output image format (e.g. ``'png'``, ``'svg'``, ``'pdf'``).
            When ``None`` (default), the format is inferred from the
            ``output_path`` file extension.

        Returns
        -------
        matplotlib.figure.Figure
            A figure containing only the legend.
        """
        handles = [
            mpatches.Patch(
                facecolor=annotation.get_color(ft),
                edgecolor='none',
                label=ft,
            )
            for ft in annotation.feature_types()
        ]
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_visible(False)
        fig.legend(handles=handles, loc='center', fontsize=10, frameon=True)
        plt.tight_layout()
        if output_path is not None:
            plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight', format=format)
        return fig

    def plot_single(
        self,
        query_name: str,
        target_name: str,
        query_group: Optional[str] = None,
        target_group: Optional[str] = None,
        output_path: Optional[Union[str, Path]] = None,
        figsize: tuple[float, float] = (6.0, 6.0),
        dot_size: float = 0.5,
        dot_color: str = 'blue',
        rc_color: str = 'red',
        merge: bool = True,
        title: Optional[str] = None,
        dpi: int = 150,
        format: Optional[str] = None,
        min_length: int = 0,
        color_by_identity: bool = False,
        identity_palette: str = 'viridis',
        annotation: Optional['GffAnnotation'] = None,
        annotation_track_size: float = 0.4,
        chain_gap: int = 0,
        rasterized: Union[bool, str] = 'auto',
        rasterization_threshold: int = 50_000,
    ) -> matplotlib.figure.Figure:
        """Plot a single pairwise dotplot.

        When *annotation* is provided, a linear annotation track is drawn
        below the x-axis (target sequence features) and to the left of the
        y-axis (query sequence features).

        When *index* is a :class:`~rusty_dot.paf_io.CrossIndex`, supply
        *query_group* and *target_group* to have the sequence names resolved
        to internal (``'group:name'``) identifiers automatically, and to
        render from pre-computed records when available.

        Parameters
        ----------
        query_name : str
            Name of the query sequence (y-axis).  When *query_group* is
            provided and *index* is a ``CrossIndex``, this is treated as an
            un-prefixed name and the internal identifier is looked up.
        target_name : str
            Name of the target sequence (x-axis).  Same note as *query_name*.
        query_group : str or None, optional
            Group label for the query sequence.  When provided and *index* is
            a ``CrossIndex``, the internal name is resolved as
            ``'{query_group}:{query_name}'``.
        target_group : str or None, optional
            Group label for the target sequence.  When provided and *index*
            is a ``CrossIndex``, the internal name is resolved as
            ``'{target_group}:{target_name}'``.
        output_path : str or Path, optional
            Output image file path.  When ``None`` (default) the figure is
            not saved to disk.  Use a ``.svg`` extension (or set
            ``format='svg'``) to produce an SVG vector image.
        figsize : tuple[float, float], optional
            Figure size as (width, height) in inches for the main dotplot
            panel.  When annotation tracks are added the overall figure will
            be slightly larger.  Default is ``(6, 6)``.
        dot_size : float, optional
            Marker/line size for each match. Default is ``0.5``.
        dot_color : str, optional
            Colour for forward-strand (``+``) matches. Default is ``"blue"``.
        rc_color : str, optional
            Colour for reverse-complement (``-``) matches. Default is ``"red"``.
        merge : bool, optional
            Whether to merge sequential k-mer runs. Default is ``True``.
        title : str, optional
            Plot title. If ``None``, a default title is used.
        dpi : int, optional
            Output image resolution. Default is ``150``.
        format : str, optional
            Output image format (e.g. ``'png'``, ``'svg'``, ``'pdf'``).
            When ``None`` (default), the format is inferred from the
            ``output_path`` file extension.
        min_length : int, optional
            Minimum alignment length to display.  Matches shorter than this
            value are not drawn.  Applies to merged k-mer runs and pre-computed
            PAF alignments.  Default is ``0`` (no filtering).
        color_by_identity : bool, optional
            When ``True``, alignments are coloured by sequence identity using
            the *identity_palette* colormap.  Requires a
            :class:`~rusty_dot.paf_io.PafAlignment` to be supplied as
            ``paf_alignment`` to :meth:`__init__`; if no PAF alignment is
            available a warning is logged and the default strand colours are
            used instead.  Default is ``False``.
        identity_palette : str, optional
            Matplotlib colormap name used to map identity values (0–1) to
            colours when ``color_by_identity=True``.  Default is
            ``'viridis'``.
        annotation : GffAnnotation, optional
            Feature annotations to display as linear tracks flanking the
            dotplot.  Target features are drawn below the x-axis; query
            features are drawn to the left of the y-axis.  Sequence names
            in *annotation* absent from the index emit a warning.
            Default is ``None``.
        annotation_track_size : float, optional
            Height/width in inches of each annotation track.
            Default is ``0.4``.
        chain_gap : int, optional
            When greater than ``0``, chain co-linear match blocks on the same
            diagonal separated by up to *chain_gap* bp into single lines before
            drawing.  Default is ``0`` (off).  See :meth:`plot`.
        rasterized : bool or str, optional
            Whether to rasterise the match layer; ``'auto'`` (default) keeps it
            true vector until the segment count exceeds
            *rasterization_threshold*.  See :meth:`plot`.
        rasterization_threshold : int, optional
            Segment count above which ``rasterized='auto'`` rasterises the
            layer.  Default is ``50_000``.

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.  In a Jupyter notebook the figure is
            displayed inline automatically; call ``matplotlib.pyplot.close``
            on the returned object when it is no longer needed.

        Raises
        ------
        ValueError
            If *query_group* / *target_group* are provided but *index* is
            not a ``CrossIndex``.
        """
        # Resolve group-prefixed names for CrossIndex.
        if query_group is not None or target_group is not None:
            if not self._index_is_cross():
                raise ValueError(
                    'query_group and target_group can only be used when index '
                    'is a CrossIndex.'
                )
            cross = self.index  # type: ignore[assignment]
            if query_group is not None:
                query_name = cross.make_internal_name(query_group, query_name)
            if target_group is not None:
                target_name = cross.make_internal_name(target_group, target_name)

        # Use pre-computed records when available (via shared helper).
        paf_override = self._get_paf_override(query_group, target_group)
        effective_paf = paf_override if paf_override is not None else self.paf_alignment

        import matplotlib.gridspec as gridspec

        if annotation is not None:
            # Warn about annotation sequences not in the index.
            index_seqs = set(self.index.sequence_names())
            for ann_seq in annotation.sequence_names():
                if ann_seq not in index_seqs:
                    _log.warning(
                        'Annotation contains features for sequence %r which is '
                        'not present in the index. These features will not be '
                        'plotted.',
                        ann_seq,
                    )
            x_feats = annotation.get_features_for_sequence(target_name)
            y_feats = annotation.get_features_for_sequence(query_name)
            has_tracks = True
        else:
            x_feats = []
            y_feats = []
            has_tracks = False

        if has_tracks:
            fw, fh = figsize
            ts = annotation_track_size
            # GridSpec layout:
            #   rows: [main (fh), x-track (ts)]
            #   cols: [y-track (ts), main (fw)]
            total_w = fw + ts
            total_h = fh + ts
            fig = plt.figure(figsize=(total_w, total_h))
            gs = gridspec.GridSpec(
                2,
                2,
                width_ratios=[ts, fw],
                height_ratios=[fh, ts],
                hspace=0.02,
                wspace=0.02,
            )
            main_ax = fig.add_subplot(gs[0, 1])
            y_track_ax = fig.add_subplot(gs[0, 0], sharey=main_ax)
            x_track_ax = fig.add_subplot(gs[1, 1], sharex=main_ax)
            corner_ax = fig.add_subplot(gs[1, 0])
            corner_ax.set_visible(False)
        else:
            fig, main_ax = plt.subplots(figsize=figsize)

        self._plot_panel(
            main_ax,
            query_name,
            target_name,
            dot_size=dot_size,
            dot_color=dot_color,
            rc_color=rc_color,
            merge=merge,
            min_length=min_length,
            color_by_identity=color_by_identity,
            identity_palette=identity_palette,
            paf_alignment_override=effective_paf,
            chain_gap=chain_gap,
            rasterized=rasterized,
            rasterization_threshold=rasterization_threshold,
        )

        if has_tracks:
            t_len = self.index.get_sequence_length(target_name)

            # Hide main-axis tick labels that duplicate the track labels.
            plt.setp(main_ax.get_xticklabels(), visible=False)
            plt.setp(main_ax.get_yticklabels(), visible=False)

            # Display names for track axis labels (strip group prefix).
            display_t = self._strip_group_prefix(target_name)
            display_q = self._strip_group_prefix(query_name)

            # ---- x-annotation track (below x-axis: target features) ----
            x_track_ax.set_xlim(0, t_len)
            x_track_ax.set_ylim(0, 1)
            x_track_ax.set_yticks([])
            x_track_ax.tick_params(axis='x', labelsize=6)
            x_track_ax.set_xlabel(display_t, fontsize=8)
            for feat in x_feats:
                rect = mpatches.Rectangle(
                    (feat.start, 0.1),
                    feat.end - feat.start,
                    0.8,
                    facecolor=annotation.get_color(feat.feature_type),  # type: ignore[union-attr]
                    edgecolor='none',
                )
                x_track_ax.add_patch(rect)

            # ---- y-annotation track (left of y-axis: query features) ----
            # The main axes y-axis is inverted, so sharey keeps inversion.
            y_track_ax.set_xlim(1, 0)  # reversed so features face main plot
            y_track_ax.set_xticks([])
            y_track_ax.tick_params(axis='y', labelsize=6)
            y_track_ax.set_ylabel(display_q, fontsize=8)
            for feat in y_feats:
                rect = mpatches.Rectangle(
                    (0.1, feat.start),
                    0.8,
                    feat.end - feat.start,
                    facecolor=annotation.get_color(feat.feature_type),  # type: ignore[union-attr]
                    edgecolor='none',
                )
                y_track_ax.add_patch(rect)

        # Title: use display names (strip group prefix for CrossIndex).
        if title is None:
            dq = self._strip_group_prefix(query_name)
            dt = self._strip_group_prefix(target_name)
            title = f'{dq} vs {dt}'
        main_ax.set_title(title, fontsize=10)

        if has_tracks:
            fig.subplots_adjust(hspace=0.02, wspace=0.02)
        else:
            plt.tight_layout()
        if output_path is not None:
            plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight', format=format)
        return fig

    def plot_identity_colorbar(
        self,
        palette: str = 'viridis',
        figsize: tuple[float, float] = (1.5, 4.0),
        output_path: Optional[Union[str, Path]] = None,
        dpi: int = 150,
        format: Optional[str] = None,
    ) -> matplotlib.figure.Figure:
        """Render the identity colour scale as a standalone figure.

        Produces a figure containing only a vertical colorbar that maps
        identity values (0–100 %) to colours from *palette*.  This is
        intended to be displayed alongside a dotplot produced with
        ``color_by_identity=True``.

        Parameters
        ----------
        palette : str, optional
            Matplotlib colormap name.  Should match the *identity_palette*
            used when calling :meth:`plot` or :meth:`plot_single`.
            Default is ``'viridis'``.
        figsize : tuple[float, float], optional
            Figure size as ``(width, height)`` in inches.
            Default is ``(1.5, 4.0)``.
        output_path : str or Path, optional
            Output image file path.  When ``None`` (default) the figure is
            not saved to disk.
        dpi : int, optional
            Output image resolution. Default is ``150``.
        format : str, optional
            Output image format (e.g. ``'png'``, ``'svg'``, ``'pdf'``).
            When ``None`` (default), the format is inferred from the
            ``output_path`` file extension.

        Returns
        -------
        matplotlib.figure.Figure
            A figure containing only the colorbar.
        """
        norm = mcolors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap=plt.get_cmap(palette), norm=norm)
        sm.set_array([])
        fig, ax = plt.subplots(figsize=figsize)
        cb = fig.colorbar(sm, ax=ax, orientation='vertical')
        cb.set_label('Identity', fontsize=10)
        cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
        cb.set_ticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax.set_visible(False)
        plt.tight_layout()
        if output_path is not None:
            plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight', format=format)
        return fig
