"""Plot configuration state for the browser app.

A single dataclass captures every display option that can change without
recomputing an alignment.  The app re-renders from cached alignment results
whenever this config changes, and remembers it across sub-plot drill-down.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

#: Contig-ordering choices exposed in the UI, mapped to display labels.
ORDER_CHOICES: dict[str, str] = {
    'input': 'Input order',
    'length': 'Sort by size',
    'colinearity': 'Maximise colinearity',
}


@dataclass
class PlotConfig:
    """Display options for the dotplot, independent of the alignment.

    Attributes
    ----------
    contig_order : str
        One of ``'input'``, ``'length'`` or ``'colinearity'`` (keys of
        :data:`ORDER_CHOICES`).
    auto_reverse : bool
        Flip contigs detected as reverse-oriented so colinear alignments
        render as forward diagonals.
    hide_internal_axes : bool
        Remove internal ticks/spines so the grid reads as one continuous
        plot.
    nature : bool
        Apply the Nature-journal rcParams style while rendering.
    dot_size : float
        Line width for match segments.
    min_length : int
        Minimum match length (bp) to draw.
    title : str or None
        Optional figure title.
    """

    contig_order: str = 'input'
    auto_reverse: bool = False
    hide_internal_axes: bool = False
    nature: bool = False
    dot_size: float = 0.5
    min_length: int = 0
    title: str | None = field(default=None)

    def plot_kwargs(self) -> dict[str, Any]:
        """Translate this config into ``DotPlotter.plot`` keyword arguments.

        Returns
        -------
        dict[str, Any]
            Keyword arguments accepted by
            :meth:`rusty_dot.dotplot.DotPlotter.plot`.
        """
        order = None if self.contig_order == 'input' else self.contig_order
        return {
            'contig_order': order,
            'auto_reverse': self.auto_reverse,
            'hide_internal_axes': self.hide_internal_axes,
            'dot_size': self.dot_size,
            'min_length': self.min_length,
            'title': self.title,
        }
