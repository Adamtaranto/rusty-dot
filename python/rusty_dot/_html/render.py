"""Assemble the single-file interactive HTML dotplot report.

:func:`render_html_report` renders a matplotlib figure to SVG in memory,
loads the packaged template/CSS/JS assets via :mod:`importlib.resources`,
and splices everything into one self-contained HTML document — no external
requests, no vendored libraries, no template engine (plain ``str.replace``
on unique slot tokens).
"""

from __future__ import annotations

from importlib.resources import files
import io
import json
import logging
from pathlib import Path
from typing import Any, Union

import matplotlib.figure

_log = logging.getLogger(__name__)

# Slot tokens replaced in template.html.  Chosen to be valid nowhere in
# HTML/CSS/JS syntax so a stray token is easy to spot in output.
_SLOT_TITLE = '__RD_TITLE__'
_SLOT_CSS = '__RD_CSS__'
_SLOT_SVG = '__RD_SVG__'
_SLOT_DATA = '__RD_DATA__'
_SLOT_JS = '__RD_JS__'


def _load_asset(name: str) -> str:
    """Read a packaged asset file from ``rusty_dot._html``.

    Parameters
    ----------
    name : str
        Asset filename (e.g. ``'template.html'``).

    Returns
    -------
    str
        The asset's text content.
    """
    return files('rusty_dot._html').joinpath(name).read_text(encoding='utf-8')


def _figure_to_svg(fig: matplotlib.figure.Figure) -> str:
    """Render a figure to an SVG fragment string.

    The XML prolog and DOCTYPE emitted by matplotlib are stripped so the
    ``<svg>`` element can be embedded inline in the HTML document.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to render.  Artists tagged with ``set_gid`` keep their
        ids in the SVG output.

    Returns
    -------
    str
        The ``<svg>...</svg>`` markup.
    """
    buf = io.StringIO()
    # bbox_inches='tight' matches the raster save path so HTML output frames
    # the grid the same way as PNG/SVG output.
    fig.savefig(buf, format='svg', bbox_inches='tight')
    svg = buf.getvalue()
    start = svg.find('<svg')
    if start == -1:  # pragma: no cover - matplotlib always emits <svg>
        raise ValueError('matplotlib SVG output contained no <svg> element')
    return svg[start:]


def render_html_report(
    fig: matplotlib.figure.Figure,
    payload: dict[str, Any],
    output_path: Union[str, Path],
    title: str = 'rusty-dot report',
) -> Path:
    """Write a self-contained interactive HTML report for a dotplot figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The dotplot grid figure.  Panel axes and match collections must be
        gid-tagged (``rd-panel-*`` / ``rd-matches-*``) so the embedded JS can
        find them; :meth:`~rusty_dot.dotplot.DotPlotter.plot` does this when
        HTML output is requested.
    payload : dict
        JSON-ready match payload from
        :func:`rusty_dot._html.serialize.build_panel_payload`.
    output_path : str or Path
        Destination ``.html`` file path.
    title : str, optional
        Document and header title.  Default is ``'rusty-dot report'``.

    Returns
    -------
    pathlib.Path
        The path the report was written to.
    """
    output_path = Path(output_path)
    svg = _figure_to_svg(fig)
    data = json.dumps(payload, separators=(',', ':'))
    # A '</script>' inside JSON text would terminate the data block early;
    # escaping '</' is invisible to JSON.parse but safe in HTML.
    data = data.replace('</', '<\\/')

    html = _load_asset('template.html')
    html = html.replace(_SLOT_TITLE, title)
    html = html.replace(_SLOT_CSS, _load_asset('report.css'))
    html = html.replace(_SLOT_JS, _load_asset('report.js'))
    html = html.replace(_SLOT_DATA, data)
    html = html.replace(_SLOT_SVG, svg)

    output_path.write_text(html, encoding='utf-8')
    _log.info(
        'Wrote interactive HTML report to %s (%d panel(s), %.1f kB)',
        output_path,
        len(payload.get('panels', {})),
        output_path.stat().st_size / 1024,
    )
    return output_path
