"""
Dotplot visualization module for rusty-dot.

Provides the DotPlotter class for generating all-vs-all dotplots from
DNA sequence comparison data.

Reference: https://github.com/rrwick/Autocycler/blob/b0523350898faac71686251ec58f7d83bc2b1c28/src/dotplot.rs
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from rusty_dot._rusty_dot import SequenceIndex


class DotPlotter:
    """Generate all-vs-all dotplots for sets of DNA sequences.

    Parameters
    ----------
    index : SequenceIndex
        A SequenceIndex with the sequences to plot.

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

    def __init__(self, index: SequenceIndex) -> None:
        """Initialise the DotPlotter.

        Parameters
        ----------
        index : SequenceIndex
            A populated SequenceIndex instance.
        """
        self.index = index

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
        scale_sequences: bool = False,
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
            Output image file path. Default is ``"dotplot.png"``.
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
            When ``True``, subplot widths and heights are proportional to
            the lengths of the corresponding sequences so that relative
            sequence sizes are preserved.  When ``False`` (default) every
            panel has the same fixed size.
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
                )

        if title:
            fig.suptitle(title, fontsize=14, y=1.01)

        plt.tight_layout()
        plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight')
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
        """
        q_len = self.index.get_sequence_length(query_name)
        t_len = self.index.get_sequence_length(target_name)

        matches = self.index.compare_sequences_stranded(query_name, target_name, merge)

        # Draw match lines/dots; RC matches are drawn as anti-diagonal lines.
        for q_start, q_end, t_start, t_end, strand in matches:
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
        ax.set_xlabel(target_name, fontsize=8)
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
    ) -> None:
        """Plot a single pairwise dotplot.

        Parameters
        ----------
        query_name : str
            Name of the query sequence (y-axis).
        target_name : str
            Name of the target sequence (x-axis).
        output_path : str or Path, optional
            Output image file path. Default is ``"dotplot.png"``.
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
        )
        if title is None:
            title = f'{query_name} vs {target_name}'
        ax.set_title(title, fontsize=10)
        plt.tight_layout()
        plt.savefig(str(output_path), dpi=dpi, bbox_inches='tight')
        plt.close(fig)
