"""
Nature-journal-style plot formatting for rusty-dot.

Provides an opt-in matplotlib style tuned to Nature journal figure
guidelines: Helvetica/Arial sans-serif fonts at 5-7 pt, 0.5 pt line and
axes widths, no top/right spines, outward ticks, and 300 dpi tight-bbox
figure export defaults.

The style is never applied automatically — wrap plotting calls in
:func:`nature_style` or call :func:`use_nature_style` to opt in.

Examples
--------
Scoped application via the context manager:

>>> from rusty_dot import DotPlotter, nature_style
>>> with nature_style():  # doctest: +SKIP
...     fig = plotter.plot(output_path='dotplot.png')

Global application:

>>> from rusty_dot import use_nature_style
>>> use_nature_style()  # doctest: +SKIP
"""

from __future__ import annotations

from contextlib import AbstractContextManager
from typing import Any, cast

import matplotlib

NATURE_RC: dict[str, Any] = {
    # Fonts: Helvetica/Arial with sensible sans-serif fallbacks.
    'font.family': 'sans-serif',
    'font.sans-serif': [
        'Helvetica',
        'Arial',
        'Helvetica Neue',
        'Liberation Sans',
        'Nimbus Sans',
        'DejaVu Sans',
    ],
    # Font sizes: Nature figures use 5-7 pt text.
    'font.size': 7.0,
    'axes.titlesize': 7.0,
    'axes.labelsize': 7.0,
    'xtick.labelsize': 6.0,
    'ytick.labelsize': 6.0,
    'legend.fontsize': 6.0,
    'legend.title_fontsize': 7.0,
    'figure.titlesize': 7.0,
    # Line and axes widths: 0.5 pt throughout.
    'lines.linewidth': 0.5,
    'axes.linewidth': 0.5,
    'grid.linewidth': 0.5,
    'patch.linewidth': 0.5,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.5,
    'ytick.minor.width': 0.5,
    # Spines: hide top and right.
    'axes.spines.top': False,
    'axes.spines.right': False,
    # Ticks point outward.
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    # Figure export: 300 dpi with a tight bounding box.  Note that explicit
    # keyword arguments always win over rcParams, so methods that pass
    # ``dpi=`` / ``bbox_inches=`` to ``savefig`` explicitly (e.g.
    # ``DotPlotter.plot(dpi=...)``) are unaffected by these two defaults.
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
}
"""rcParams implementing Nature journal figure-formatting guidelines.

Keys mirror :data:`matplotlib.rcParams`: Helvetica/Arial sans-serif fonts
at 5-7 pt, 0.5 pt line and axes widths, no top/right spines, outward
ticks, and 300 dpi tight-bbox savefig defaults.
"""


def nature_style() -> AbstractContextManager[None]:
    """Return a context manager that applies the Nature style temporarily.

    Inside the ``with`` block :data:`matplotlib.rcParams` reflects
    :data:`NATURE_RC`; on exit the previous rcParams are restored,
    even if an exception was raised.

    Returns
    -------
    contextlib.AbstractContextManager
        A :func:`matplotlib.rc_context` context manager configured with
        :data:`NATURE_RC`.

    Examples
    --------
    >>> import matplotlib
    >>> from rusty_dot.style import nature_style
    >>> with nature_style():
    ...     assert matplotlib.rcParams['savefig.dpi'] == 300
    """
    # Cast: matplotlib stubs key rc dicts by a Literal of every rcParam name,
    # which a dict[str, Any] constant cannot satisfy directly.
    return matplotlib.rc_context(rc=cast('dict[Any, Any]', NATURE_RC))


def use_nature_style() -> None:
    """Apply the Nature style globally to :data:`matplotlib.rcParams`.

    Unlike :func:`nature_style`, the change persists for the rest of the
    session (or until rcParams are reset, e.g. via
    :func:`matplotlib.rcdefaults`).

    Returns
    -------
    None
        This function mutates :data:`matplotlib.rcParams` in place.

    Examples
    --------
    >>> import matplotlib
    >>> from rusty_dot.style import use_nature_style
    >>> use_nature_style()
    >>> matplotlib.rcParams['axes.spines.top']
    False
    >>> matplotlib.rcdefaults()  # reset to matplotlib library defaults
    """
    matplotlib.rcParams.update(cast('dict[Any, Any]', NATURE_RC))
