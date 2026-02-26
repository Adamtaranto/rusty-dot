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

import matplotlib.colors as mcolors
import matplotlib.figure
import matplotlib.pyplot as plt

from rusty_dot._rusty_dot import SequenceIndex

if TYPE_CHECKING:
    from rusty_dot.paf_io import CrossIndex, PafAlignment

_log = logging.getLogger(__name__)


class DotPlotter:
    """Generate all-vs-all dotplots for sets of DNA sequences.

    Accepts either a :class:`~rusty_dot.SequenceIndex` (single sequence
    collection) or a :class:`~rusty_dot.paf_io.CrossIndex` (multi-group
    collection).  When using a ``CrossIndex``, pass group-specific names via
    ``query_names`` and ``target_names``::

        cross = CrossIndex(k=15)
        cross.load_fasta("assembly_a.fasta", group="a")
        cross.load_fasta("assembly_b.fasta", group="b")

        plotter = DotPlotter(cross)
        plotter.plot(
            query_names=cross.sequence_names(group="a"),
            target_names=cross.sequence_names(group="b"),
            output_path="cross_plot.png",
        )

    To colour alignments by sequence identity, supply a
    :class:`~rusty_dot.paf_io.PafAlignment` and set
    ``color_by_identity=True``::

        from rusty_dot.paf_io import PafAlignment
        aln = PafAlignment.from_file("alignments.paf")
        plotter = DotPlotter(idx, paf_alignment=aln)
        fig = plotter.plot(color_by_identity=True, identity_palette="viridis")
        cbar = plotter.plot_identity_colorbar(palette="viridis")

    Parameters
    ----------
    index : SequenceIndex or CrossIndex
        A populated index instance.
    paf_alignment : PafAlignment, optional
        Pre-loaded PAF alignments used as the data source when
        ``color_by_identity=True``.  When ``None`` (default) k-mer matches
        from *index* are used for plotting.

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
        index: Union[SequenceIndex, 'CrossIndex'],
        paf_alignment: Optional['PafAlignment'] = None,
    ) -> None:
        """Initialise the DotPlotter.

        Parameters
        ----------
        index : SequenceIndex or CrossIndex
            A populated index instance.
        paf_alignment : PafAlignment, optional
            Pre-loaded PAF alignments.  Required for identity-based colouring.
            When ``None`` (default), k-mer matches from *index* are used.
        """
        self.index = index
        self.paf_alignment = paf_alignment

    def plot(
        self,
        query_names: Optional[list[str]] = None,
        target_names: Optional[list[str]] = None,
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
    ) -> matplotlib.figure.Figure:
        """Plot an all-vs-all dotplot grid.

        If both ``query_names`` and ``target_names`` are provided, the plot
        will show each query sequence (rows) against each target sequence
        (columns). If only one set is provided, or neither, all pairwise
        combinations within the available sequences are plotted.

        The figure is always returned so it can be displayed inline in a
        Jupyter notebook.  When ``output_path`` is provided the figure is
        also saved to disk.

        Parameters
        ----------
        query_names : list[str], optional
            Sequence names for the y-axis (rows). If ``None``, uses all
            sequences in the index.
        target_names : list[str], optional
            Sequence names for the x-axis (columns). If ``None``, uses all
            sequences in the index.
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

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.  In a Jupyter notebook the figure is
            displayed inline automatically; call ``matplotlib.pyplot.close``
            on the returned object when it is no longer needed.
        """
        all_names = self.index.sequence_names()
        if not all_names:
            raise ValueError('No sequences in the index.')

        if query_names is None:
            query_names = sorted(all_names)
        if target_names is None:
            target_names = sorted(all_names)

        nrows = len(query_names)
        ncols = len(target_names)

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
                    # Only label the leftmost column (y) and bottom row (x)
                    show_xlabel=(row_idx == nrows - 1),
                    show_ylabel=(col_idx == 0),
                    color_by_identity=color_by_identity,
                    identity_palette=identity_palette,
                )

        if title:
            fig.suptitle(title, fontsize=14, y=1.01)

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
            *identity_palette*.  Requires ``self.paf_alignment`` to be set;
            if not, a warning is logged and strand colours are used instead.
            Default is ``False``.
        identity_palette : str, optional
            Matplotlib colormap name for identity-based colouring.
            Default is ``'viridis'``.
        """
        q_len = self.index.get_sequence_length(query_name)
        t_len = self.index.get_sequence_length(target_name)

        if color_by_identity and self.paf_alignment is None:
            _log.warning(
                'color_by_identity=True requires a PafAlignment; k-mer matches '
                'are always 100% identity. Pass paf_alignment= to DotPlotter '
                'to enable identity colouring.'
            )
            color_by_identity = False

        if color_by_identity:
            # Use PAF records for this sequence pair, coloured by identity.
            cmap = plt.get_cmap(identity_palette)
            norm = mcolors.Normalize(vmin=0, vmax=1)
            records = [
                r
                for r in self.paf_alignment.records  # type: ignore[union-attr]
                if r.query_name == query_name and r.target_name == target_name
            ]
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
                    xs = [rec.target_end, rec.target_start]
                else:
                    xs = [rec.target_start, rec.target_end]
                ax.plot(
                    xs,
                    [rec.query_start, rec.query_end],
                    color=color,
                    linewidth=dot_size,
                    alpha=0.7,
                )
        else:
            # Draw match lines/dots from k-mer index; RC matches are drawn as
            # anti-diagonal lines.
            matches = self.index.compare_sequences_stranded(query_name, target_name, merge)
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
    ) -> matplotlib.figure.Figure:
        """Plot a single pairwise dotplot.

        Parameters
        ----------
        query_name : str
            Name of the query sequence (y-axis).
        target_name : str
            Name of the target sequence (x-axis).
        output_path : str or Path, optional
            Output image file path.  When ``None`` (default) the figure is
            not saved to disk.  Use a ``.svg`` extension (or set
            ``format='svg'``) to produce an SVG vector image.
        figsize : tuple[float, float], optional
            Figure size as (width, height) in inches. Default is ``(6, 6)``.
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

        Returns
        -------
        matplotlib.figure.Figure
            The generated figure.  In a Jupyter notebook the figure is
            displayed inline automatically; call ``matplotlib.pyplot.close``
            on the returned object when it is no longer needed.
        """
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
            identity_palette=identity_palette,
        )
        if title is None:
            title = f'{query_name} vs {target_name}'
        ax.set_title(title, fontsize=10)
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
