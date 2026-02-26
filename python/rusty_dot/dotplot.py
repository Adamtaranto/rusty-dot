"""
Dotplot visualization module for rusty-dot.

Provides the DotPlotter class for generating all-vs-all dotplots from
DNA sequence comparison data.

Reference: https://github.com/rrwick/Autocycler/blob/b0523350898faac71686251ec58f7d83bc2b1c28/src/dotplot.rs
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

import matplotlib
import matplotlib.cm

matplotlib.use('Agg')
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

from rusty_dot._rusty_dot import SequenceIndex

if TYPE_CHECKING:
    from rusty_dot.paf_io import CrossIndex, PafAlignment

_log = logging.getLogger(__name__)



class DotPlotter:
    """Generate all-vs-all dotplots for sets of DNA sequences.

    Accepts a :class:`~rusty_dot.SequenceIndex` (single sequence collection),
    a :class:`~rusty_dot.paf_io.CrossIndex` (multi-group collection), or a
    :class:`~rusty_dot.paf_io.PafAlignment` (alignments loaded from a PAF
    file).

    When using a ``CrossIndex``, pass group-specific names via ``query_names``
    and ``target_names``::

        cross = CrossIndex(k=15)
        cross.load_fasta("assembly_a.fasta", group="a")
        cross.load_fasta("assembly_b.fasta", group="b")

        plotter = DotPlotter(cross)
        plotter.plot(
            query_names=cross.sequence_names(group="a"),
            target_names=cross.sequence_names(group="b"),
            output_path="cross_plot.png",
        )

    When using a ``PafAlignment``, alignments are drawn directly from the
    PAF records and identity data is available for colour-mapping::

        from rusty_dot.paf_io import PafAlignment
        aln = PafAlignment.from_file("alignments.paf")
        plotter = DotPlotter(aln)
        plotter.plot(
            query_names=aln.query_names,
            target_names=aln.target_names,
            color_by_identity=True,
            output_path="dotplot.png",
        )
        plotter.plot_identity_scale("identity_scale.png")

    Parameters
    ----------
    index : SequenceIndex, CrossIndex, or PafAlignment
        A populated index or PAF alignment instance.

    Examples
    --------
    >>> from rusty_dot import SequenceIndex
    >>> from rusty_dot.dotplot import DotPlotter
    >>> idx = SequenceIndex(k=10)
    >>> idx.add_sequence("seq1", "ACGTACGTACGT" * 10)
    >>> idx.add_sequence("seq2", "TACGTACGTACG" * 10)
    >>> plotter = DotPlotter(idx)
    >>> plotter.plot(output_path="dotplot.png")
    """

    def __init__(self, index: Union[SequenceIndex, 'CrossIndex', 'PafAlignment']) -> None:
        """Initialise the DotPlotter.

        Parameters
        ----------
        index : SequenceIndex, CrossIndex, or PafAlignment
            A populated index or PAF alignment instance.
        """
        self.index = index
        # Pre-build a length map when the source is a PafAlignment so that
        # _get_sequence_length() does not need to re-scan records each call.
        self._paf_len_map: dict[str, int] = {}
        if self._is_paf_source:
            for rec in index.records:  # type: ignore[union-attr]
                self._paf_len_map[rec.query_name] = rec.query_len
                self._paf_len_map[rec.target_name] = rec.target_len

    # ------------------------------------------------------------------
    # Internal helpers for dual-source (k-mer / PAF) support
    # ------------------------------------------------------------------

    @property
    def _is_paf_source(self) -> bool:
        """Return ``True`` when the index is a :class:`~rusty_dot.paf_io.PafAlignment`."""
        from rusty_dot.paf_io import PafAlignment

        return isinstance(self.index, PafAlignment)

    def _sequence_names(self) -> list[str]:
        """Return all sequence names, handling both index and PAF sources."""
        if self._is_paf_source:
            return list(self._paf_len_map.keys())
        return self.index.sequence_names()  # type: ignore[union-attr]

    def _get_sequence_length(self, name: str) -> int:
        """Return sequence length, handling both index and PAF sources."""
        if self._is_paf_source:
            if name not in self._paf_len_map:
                raise KeyError(
                    f'Sequence {name!r} not found in the PAF alignment records.'
                )
            return self._paf_len_map[name]
        return self.index.get_sequence_length(name)  # type: ignore[union-attr]

    def _get_paf_matches(
        self,
        query_name: str,
        target_name: str,
        min_length: int,
    ) -> list[tuple[int, int, int, int, str, Optional[float]]]:
        """Return alignment tuples from PAF records for the given sequence pair.

        Parameters
        ----------
        query_name : str
            Query sequence name.
        target_name : str
            Target sequence name.
        min_length : int
            Minimum query-aligned length.

        Returns
        -------
        list of (q_start, q_end, t_start, t_end, strand, identity_or_None)
        """
        result: list[tuple[int, int, int, int, str, Optional[float]]] = []
        for rec in self.index.records:  # type: ignore[union-attr]
            if rec.query_name != query_name or rec.target_name != target_name:
                continue
            if min_length > 0 and rec.query_aligned_len < min_length:
                continue
            identity: Optional[float] = (
                rec.residue_matches / rec.alignment_block_len
                if rec.alignment_block_len > 0
                else None
            )
            result.append(
                (
                    rec.query_start,
                    rec.query_end,
                    rec.target_start,
                    rec.target_end,
                    rec.strand,
                    identity,
                )
            )
        return result

    def plot(
        self,
        query_names: Optional[list[str]] = None,
        target_names: Optional[list[str]] = None,
        output_path: Union[str, Path] = 'dotplot.png',
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
        palette: str = 'plasma',
    ) -> None:
        """Plot an all-vs-all dotplot grid.

        If both ``query_names`` and ``target_names`` are provided, the plot
        will show each query sequence (rows) against each target sequence
        (columns). If only one set is provided, or neither, all pairwise
        combinations within the available sequences are plotted.

        Parameters
        ----------
        query_names : list[str], optional
            Sequence names for the y-axis (rows). If ``None``, uses all
            sequences in the index.
        target_names : list[str], optional
            Sequence names for the x-axis (columns). If ``None``, uses all
            sequences in the index.
        output_path : str or Path, optional
            Output image file path. Default is ``"dotplot.png"``. Use a
            ``.svg`` extension (or set ``format='svg'``) to produce an SVG
            vector image.
        figsize_per_panel : float, optional
            Base size in inches for each subplot panel when
            ``scale_sequences=False``.  When ``scale_sequences=True`` this
            value sets the size of the *longest* sequence axis and all
            other axes are scaled proportionally.  Default is ``4.0``.
        dot_size : float, optional
            Size of each dot in the scatter plot. Default is ``0.5``.
        dot_color : str, optional
            Colour for forward-strand (``+``) match lines. Default is ``"blue"``.
            Ignored when ``color_by_identity=True``.
        rc_color : str, optional
            Colour for reverse-complement (``-``) strand match lines.
            Default is ``"red"``.  Ignored when ``color_by_identity=True``.
        merge : bool, optional
            Whether to merge sequential k-mer runs before plotting.
            Default is ``True``.  Ignored when the source is a
            :class:`~rusty_dot.paf_io.PafAlignment`.
        title : str, optional
            Overall figure title. If ``None``, no title is added.
        dpi : int, optional
            Resolution of the output image. Default is ``150``.
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
            When ``True``, colour each alignment by its sequence identity using
            the colourmap specified by *palette*.  Identity data is only
            available when the source is a
            :class:`~rusty_dot.paf_io.PafAlignment`; a warning is emitted and
            the parameter has no effect for k-mer match sources.
            Default is ``False``.
        palette : str, optional
            Matplotlib colourmap name used for identity colouring (e.g.
            ``'plasma'``, ``'viridis'``, ``'RdYlGn'``).  Only used when
            ``color_by_identity=True``.  Default is ``'plasma'``.
        """
        all_names = self._sequence_names()
        if not all_names:
            raise ValueError('No sequences in the index.')

        if color_by_identity and not self._is_paf_source:
            _log.warning(
                'color_by_identity=True has no effect when the source is k-mer '
                'matches; individual k-mer matches are necessarily 100%% identity. '
                'Load alignments from a PAF file to enable identity colouring.'
            )

        if query_names is None:
            query_names = sorted(all_names)
        if target_names is None:
            target_names = sorted(all_names)

        nrows = len(query_names)
        ncols = len(target_names)

        if scale_sequences:
            q_lens = [self._get_sequence_length(n) for n in query_names]
            t_lens = [self._get_sequence_length(n) for n in target_names]
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
                    color_by_identity=color_by_identity,
                    palette=palette,
                    # Only label the leftmost column (y) and bottom row (x)
                    show_xlabel=(row_idx == nrows - 1),
                    show_ylabel=(col_idx == 0),
                )

        if title:
            fig.suptitle(title, fontsize=14, y=1.01)

        plt.tight_layout()
        plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight', format=format)
        plt.close(fig)

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
        palette: str = 'plasma',
    ) -> None:
        """Render a single comparison panel onto the given Axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes to draw on.
        query_name : str
            Name of the query sequence (y-axis).
        target_name : str
            Name of the target sequence (x-axis).
        dot_size : float, optional
            Marker size. Default is ``0.5``.
        dot_color : str, optional
            Marker colour for forward-strand (``+``) matches. Default is ``"blue"``.
            Ignored when ``color_by_identity=True``.
        rc_color : str, optional
            Marker colour for reverse-complement (``-``) matches. Default is ``"red"``.
            Ignored when ``color_by_identity=True``.
        merge : bool, optional
            Whether to merge sequential runs. Default is ``True``.  Ignored
            when the source is a :class:`~rusty_dot.paf_io.PafAlignment`.
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
            When ``True``, colour each alignment by its sequence identity.
            Only meaningful for :class:`~rusty_dot.paf_io.PafAlignment` sources.
            Default is ``False``.
        palette : str, optional
            Matplotlib colourmap name for identity colouring.
            Default is ``'plasma'``.
        """
        q_len = self._get_sequence_length(query_name)
        t_len = self._get_sequence_length(target_name)

        # Build colourmap / norm once per panel if needed.
        cmap = None
        norm = None
        if color_by_identity and self._is_paf_source:
            cmap = matplotlib.colormaps[palette]
            norm = mcolors.Normalize(vmin=0.0, vmax=1.0)

        if self._is_paf_source:
            # Draw from PAF records; identity data available.
            matches_with_id = self._get_paf_matches(query_name, target_name, min_length)
            for q_start, q_end, t_start, t_end, strand, identity in matches_with_id:
                if strand == '-':
                    xs = [t_end, t_start]
                else:
                    xs = [t_start, t_end]
                if cmap is not None and norm is not None and identity is not None:
                    color: Union[str, tuple] = cmap(norm(identity))
                elif strand == '-':
                    color = rc_color
                else:
                    color = dot_color
                ax.plot(
                    xs,
                    [q_start, q_end],
                    color=color,
                    linewidth=dot_size,
                    alpha=0.7,
                )
        else:
            matches = self.index.compare_sequences_stranded(query_name, target_name, merge)  # type: ignore[union-attr]

            # Draw match lines/dots; RC matches are drawn as anti-diagonal lines.
            for q_start, q_end, t_start, t_end, strand in matches:
                if min_length > 0 and (q_end - q_start) < min_length:
                    continue
                if strand == '-':
                    # Reverse complement: as query advances (q_start→q_end) the
                    # target position retreats (t_end→t_start).
                    xs = [t_end, t_start]
                    color = rc_color
                else:
                    xs = [t_start, t_end]
                    color = dot_color
                ax.plot(
                    xs,
                    [q_start, q_end],
                    color=color,
                    linewidth=dot_size,
                    alpha=0.7,
                )

        ax.set_xlim(0, t_len)
        ax.set_ylim(0, q_len)
        ax.invert_yaxis()
        if show_xlabel:
            ax.set_xlabel(target_name, fontsize=8)
        if show_ylabel:
            ax.set_ylabel(query_name, fontsize=8)
        ax.tick_params(axis='both', labelsize=6)
        ax.set_aspect('auto')

    def plot_single(
        self,
        query_name: str,
        target_name: str,
        output_path: Union[str, Path] = 'dotplot.png',
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
        palette: str = 'plasma',
    ) -> None:
        """Plot a single pairwise dotplot.

        Parameters
        ----------
        query_name : str
            Name of the query sequence (y-axis).
        target_name : str
            Name of the target sequence (x-axis).
        output_path : str or Path, optional
            Output image file path. Default is ``"dotplot.png"``. Use a
            ``.svg`` extension (or set ``format='svg'``) to produce an SVG
            vector image.
        figsize : tuple[float, float], optional
            Figure size as (width, height) in inches. Default is ``(6, 6)``.
        dot_size : float, optional
            Marker/line size for each match. Default is ``0.5``.
        dot_color : str, optional
            Colour for forward-strand (``+``) matches. Default is ``"blue"``.
            Ignored when ``color_by_identity=True``.
        rc_color : str, optional
            Colour for reverse-complement (``-``) matches. Default is ``"red"``.
            Ignored when ``color_by_identity=True``.
        merge : bool, optional
            Whether to merge sequential k-mer runs. Default is ``True``.
            Ignored when the source is a
            :class:`~rusty_dot.paf_io.PafAlignment`.
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
            When ``True``, colour each alignment by its sequence identity using
            the colourmap specified by *palette*.  Identity data is only
            available when the source is a
            :class:`~rusty_dot.paf_io.PafAlignment`; a warning is emitted and
            the parameter has no effect for k-mer match sources.
            Default is ``False``.
        palette : str, optional
            Matplotlib colourmap name used for identity colouring.
            Default is ``'plasma'``.
        """
        if color_by_identity and not self._is_paf_source:
            _log.warning(
                'color_by_identity=True has no effect when the source is k-mer '
                'matches; individual k-mer matches are necessarily 100%% identity. '
                'Load alignments from a PAF file to enable identity colouring.'
            )
        fig, ax = plt.subplots(figsize=figsize)
        self._plot_panel(
            ax,
            query_name,
            target_name,
            dot_size=dot_size,
            dot_color=dot_color,
            rc_color=rc_color,
            merge=merge,
            min_length=min_length,
            color_by_identity=color_by_identity,
            palette=palette,
        )
        if title is None:
            title = f'{query_name} vs {target_name}'
        ax.set_title(title, fontsize=10)
        plt.tight_layout()
        plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight', format=format)
        plt.close(fig)

    def plot_identity_scale(
        self,
        output_path: Union[str, Path] = 'identity_scale.png',
        palette: str = 'plasma',
        figsize: tuple[float, float] = (1.5, 4.0),
        dpi: int = 150,
        format: Optional[str] = None,
    ) -> None:
        """Render a standalone identity colour scale as a separate figure.

        Produces a vertical colourbar ranging from 0 % (low identity) to
        100 % (identical), using the requested colourmap.  Useful as a
        companion figure to a dotplot generated with ``color_by_identity=True``.

        Parameters
        ----------
        output_path : str or Path, optional
            Output image file path.  Default is ``"identity_scale.png"``.
        palette : str, optional
            Matplotlib colourmap name.  Should match the *palette* used in the
            corresponding :meth:`plot` or :meth:`plot_single` call.
            Default is ``'plasma'``.
        figsize : tuple[float, float], optional
            Figure size as (width, height) in inches.  Default is ``(1.5, 4.0)``.
        dpi : int, optional
            Output image resolution.  Default is ``150``.
        format : str, optional
            Output image format (e.g. ``'png'``, ``'svg'``, ``'pdf'``).
            When ``None`` (default), the format is inferred from
            *output_path*.
        """
        fig, ax = plt.subplots(figsize=figsize)
        cmap = matplotlib.colormaps[palette]
        norm = mcolors.Normalize(vmin=0.0, vmax=1.0)
        cb = fig.colorbar(
            matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=ax,
            orientation='vertical',
        )
        cb.set_label('Identity', fontsize=10)
        cb.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
        cb.set_ticklabels(['0%', '25%', '50%', '75%', '100%'])
        plt.tight_layout()
        plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight', format=format)
        plt.close(fig)
