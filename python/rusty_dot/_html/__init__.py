"""Interactive HTML report generation for rusty-dot dotplots.

This private package renders a :class:`~rusty_dot.dotplot.DotPlotter` figure
as a single self-contained HTML file: the matplotlib figure is embedded as
inline SVG, match coordinates (and optionally sequences) are embedded as a
JSON payload, and vanilla CSS/JS assets bundled with the package provide
panel selection, scroll zoom and click-to-inspect behaviour.
"""

from rusty_dot._html.render import render_html_report
from rusty_dot._html.serialize import build_panel_payload

__all__ = ['build_panel_payload', 'render_html_report']
