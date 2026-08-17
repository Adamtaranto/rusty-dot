"""rusty-dot assembly comparison — browser app (Shiny for Python / Shinylive).

Runs entirely client-side: uploads never leave the browser.  Under Pyodide
the rusty-dot wasm wheel bundled in ``wheels/`` is installed at startup;
run natively (``shiny run app/app.py``) it uses the installed rusty-dot.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import logging
from pathlib import Path
import sys
import tempfile

from core.align import (
    AVAILABLE_METHODS,
    METHOD_LABELS,
    paf_alignment_from_text,
    paf_text_from_alignment,
)
from core.cache import QUERY_GROUP, TARGET_GROUP, SessionCache
from core.export import reordered_fasta_text
from core.fasta import FastaInput, parse_fasta_bytes
from core.panels import panel_pair, resolve_orders
from core.state import ORDER_CHOICES, PlotConfig
import matplotlib  # noqa: F401  (ensures shinylive bundles the pyodide package)
import numpy  # noqa: F401
from shiny import App, reactive, render, req, ui

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('rusty_dot_app')

APP_DIR = Path(__file__).parent

_boot_done = False


async def ensure_rusty_dot() -> None:
    """Make :mod:`rusty_dot` importable, installing the wasm wheel if needed.

    Under Pyodide (Shinylive) the wheel bundled in the app's ``wheels/``
    directory is installed with micropip on first call.  Natively this is a
    no-op when rusty-dot is already installed.

    Raises
    ------
    RuntimeError
        If rusty-dot is unavailable and cannot be installed (no bundled
        wheel under Pyodide, or not pip-installed natively).
    """
    global _boot_done
    if _boot_done or importlib.util.find_spec('rusty_dot') is not None:
        _boot_done = True
        return
    if sys.platform == 'emscripten':
        import micropip  # noqa: PLC0415 - pyodide-only module

        wheels = sorted(APP_DIR.glob('wheels/*.whl'))
        if not wheels:
            raise RuntimeError(
                'No rusty-dot wasm wheel bundled with the app '
                '(expected app/wheels/*.whl at export time).'
            )
        logger.info('Installing bundled wheel %s', wheels[-1].name)
        await micropip.install(f'emfs:{wheels[-1]}')
        _boot_done = True
        return
    raise RuntimeError(
        'rusty-dot is not installed — run `pip install rusty-dot` '
        '(or `maturin develop` from the repo).'
    )


def _method_choices() -> dict[str, str]:
    """Return UI choices for the method selector, marking unimplemented ones.

    Returns
    -------
    dict[str, str]
        Method key -> label, with unavailable methods suffixed.
    """
    return {
        key: (label if key in AVAILABLE_METHODS else f'{label} — coming soon')
        for key, label in METHOD_LABELS.items()
    }


# --- W2: interactive plot ---------------------------------------------------
# Script injected into the generated HTML report before embedding: posts a
# message to the parent app window when a dotplot panel is double-clicked.
# app/www/bridge.js listens for it and forwards to the 'panel_dblclick' input.
_PANEL_DBLCLICK_JS = """
<script>
(function () {
  'use strict';
  document.addEventListener('dblclick', function (ev) {
    var el = ev.target;
    var g = el && el.closest ? el.closest('g[id^="rd-panel-"]') : null;
    if (!g) { return; }
    var m = /^rd-panel-(\\d+)-(\\d+)$/.exec(g.id);
    if (!m) { return; }
    window.parent.postMessage(
      {type: 'rd-panel-dblclick',
       row: parseInt(m[1], 10),
       col: parseInt(m[2], 10)},
      '*'
    );
  });
})();
</script>
"""


def inject_panel_bridge(html: str) -> str:
    """Insert the panel double-click bridge script into a report document.

    Parameters
    ----------
    html : str
        Full HTML report text produced by ``DotPlotter.to_html``.

    Returns
    -------
    str
        The report with the bridge script spliced in before ``</body>``
        (appended at the end if no closing body tag is found).
    """
    idx = html.rfind('</body>')
    if idx == -1:
        return html + _PANEL_DBLCLICK_JS
    return html[:idx] + _PANEL_DBLCLICK_JS + html[idx:]


# --- end W2 ------------------------------------------------------------------


app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h5('Input'),
        ui.input_file(
            'query_fasta',
            'Query assembly (FASTA / .gz)',
            accept=['.fa', '.fasta', '.fna', '.gz'],
        ),
        ui.input_file(
            'target_fasta',
            'Target / reference assembly (FASTA / .gz)',
            accept=['.fa', '.fasta', '.fna', '.gz'],
        ),
        ui.input_select('method', 'Alignment method', choices=_method_choices()),
        ui.panel_conditional(
            "input.method === 'kmer'",
            ui.input_slider('k', 'k-mer size', min=8, max=64, value=21, step=1),
            ui.input_checkbox('merge', 'Merge adjacent matches', True),
        ),
        ui.panel_conditional(
            "input.method === 'paf'",
            ui.input_file('paf_file', 'PAF file', accept=['.paf', '.txt']),
        ),
        ui.input_action_button('run', 'Run comparison', class_='btn-primary'),
        ui.hr(),
        ui.h5('Plot options'),
        # --- W2: interactive plot ---
        ui.input_checkbox('interactive', 'Interactive plot (zoom & drill-down)', True),
        ui.input_select('contig_order', 'Contig order', choices=ORDER_CHOICES),
        ui.input_checkbox('auto_reverse', 'Auto-flip reversed contigs', False),
        ui.input_checkbox('hide_internal_axes', 'Hide internal axes', False),
        ui.input_checkbox('nature', 'Nature journal style', False),
        ui.input_numeric('min_length', 'Min match length (bp)', 0, min=0),
        ui.input_numeric('dot_size', 'Line width', 0.5, min=0.1, max=5, step=0.1),
        ui.hr(),
        ui.h5('Downloads'),
        ui.download_button('dl_svg', 'Plot (SVG)'),
        ui.download_button('dl_pdf', 'Plot (PDF)'),
        ui.download_button('dl_paf', 'Alignment (PAF)'),
        ui.download_button('dl_fasta', 'Reordered query (FASTA)'),
        width=320,
    ),
    ui.output_ui('status'),
    # --- W2: interactive plot ---
    ui.output_ui('plot_area'),
    ui.head_content(
        ui.include_css(APP_DIR / 'www' / 'app.css'),
        ui.include_js(APP_DIR / 'www' / 'bridge.js'),
    ),
    title='rusty-dot · assembly comparison',
    fillable=True,
)


def server(input, output, session) -> None:  # noqa: A002, D103
    cache = SessionCache()
    ready = reactive.value(False)
    boot_error = reactive.value('')
    # (kind, alignment-object, {'query': FastaInput|None, 'target': ...})
    result = reactive.value(None)

    @reactive.effect
    async def _boot():
        try:
            await ensure_rusty_dot()
            ready.set(True)
            logger.info('rusty-dot ready (platform=%s)', sys.platform)
        except RuntimeError as exc:
            boot_error.set(str(exc))

    def _parse_upload(file_input, label: str) -> FastaInput:
        files = file_input()
        if not files:
            raise ValueError(f'Please upload a {label} assembly.')
        raw = Path(files[0]['datapath']).read_bytes()
        return parse_fasta_bytes(raw)

    @reactive.effect
    @reactive.event(input.run)
    async def _run():
        req(ready())
        method = input.method()
        if method not in AVAILABLE_METHODS:
            ui.notification_show(
                f'{METHOD_LABELS[method]} is not available yet.',
                type='warning',
            )
            return
        try:
            with ui.Progress(min=0, max=3) as progress:
                if method == 'kmer':
                    progress.set(0, message='Parsing assemblies…')
                    query = _parse_upload(input.query_fasta, 'query')
                    target = _parse_upload(input.target_fasta, 'target')
                    progress.set(1, message=f'Building k-mer index (k={input.k()})…')
                    index = cache.kmer_index(
                        input.k(), query, target, merge=input.merge()
                    )
                    progress.set(3, message='Done')
                    result.set(('kmer', index, {'query': query, 'target': target}))
                elif method == 'paf':
                    progress.set(0, message='Parsing PAF…')
                    files = input.paf_file()
                    if not files:
                        raise ValueError('Please upload a PAF file.')
                    text = Path(files[0]['datapath']).read_text()
                    alignment = paf_alignment_from_text(text)
                    progress.set(3, message='Done')
                    result.set(('paf', alignment, {'query': None, 'target': None}))
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=8)

    @reactive.calc
    def config() -> PlotConfig:
        return PlotConfig(
            contig_order=input.contig_order(),
            auto_reverse=input.auto_reverse(),
            hide_internal_axes=input.hide_internal_axes(),
            nature=input.nature(),
            dot_size=input.dot_size() or 0.5,
            min_length=int(input.min_length() or 0),
        )

    # --- W2: interactive plot ------------------------------------------------
    # Focused (query, target) contig pair for the drill-down view, or None
    # for the full overview grid.
    focus = reactive.value(None)

    @reactive.effect
    def _reset_focus_on_new_result():
        result()
        focus.set(None)

    @reactive.calc
    def layout():
        """Explicit plotted axis orders for the current result + config.

        Single source of truth for row/column order: the plot call, the
        panel double-click mapping and the FASTA download all use it.
        """
        res = result()
        req(res)
        cfg = config()
        kind, obj, meta = res
        if kind == 'kmer':
            records = obj.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
            q_in = list(meta['query'].names)
            t_in = list(meta['target'].names)
            lengths = {n: len(s) for n, s in meta['target'].records}
            lengths.update({n: len(s) for n, s in meta['query'].records})
        else:
            records = obj.records
            q_in = list(obj.query_names)
            t_in = list(obj.target_names)
            lengths = {n: obj.get_sequence_length(n) for n in (*q_in, *t_in)}
        q_order, t_order, reversed_q = resolve_orders(
            cfg.contig_order, records, q_in, t_in, lengths
        )
        reverse = reversed_q if cfg.auto_reverse else set()
        return {'query_names': q_order, 'target_names': t_order, 'reverse': reverse}

    def make_figure(res, cfg: PlotConfig, lay, pair=None, output_path=None):
        from rusty_dot import DotPlotter
        from rusty_dot.paf_io import PafAlignment
        from rusty_dot.style import nature_style

        kind, obj, _meta = res
        q_names = [pair[0]] if pair is not None else lay['query_names']
        t_names = [pair[1]] if pair is not None else lay['target_names']
        kwargs = cfg.plot_kwargs()
        # Ordering is fully explicit (see layout()); never reorder in-plot.
        kwargs.update(contig_order=None, auto_reverse=False)
        kwargs['reverse_contigs'] = set(lay['reverse'])
        if pair is not None and cfg.title is None:
            kwargs['title'] = f'{pair[0]} vs {pair[1]}'
        if kind == 'kmer':
            # Internal 'group:name' identifiers keep name resolution
            # unambiguous; the cached cross-group records are handed to the
            # plotter directly so nothing is recomputed.
            paf = PafAlignment(obj.get_records_for_pair(QUERY_GROUP, TARGET_GROUP))
            plotter = DotPlotter(obj, paf_alignment=paf)
            kwargs['query_names'] = [
                obj.make_internal_name(QUERY_GROUP, n) for n in q_names
            ]
            kwargs['target_names'] = [
                obj.make_internal_name(TARGET_GROUP, n) for n in t_names
            ]
        else:
            plotter = DotPlotter(obj)
            kwargs['query_names'] = list(q_names)
            kwargs['target_names'] = list(t_names)
        style = nature_style() if cfg.nature else contextlib.nullcontext()
        with style:
            if output_path is not None:
                return plotter.to_html(output_path, **kwargs)
            return plotter.plot(**kwargs)

    def _report_html(pair=None) -> str:
        """Render the interactive HTML report and return it as a string."""
        import matplotlib.pyplot as plt

        res = result()
        req(res)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'report.html'
            fig = make_figure(res, config(), layout(), pair=pair, output_path=path)
            html = path.read_text(encoding='utf-8')
        plt.close(fig)
        return inject_panel_bridge(html)

    @reactive.calc
    def overview_html() -> str:
        # Cached by reactive.calc: entering/leaving the drill-down view does
        # not invalidate it, so 'Back to overview' re-embeds the same string.
        return _report_html(None)

    @reactive.calc
    def focus_html() -> str:
        pair = focus()
        req(pair)
        return _report_html(pair)

    @render.ui
    def report_frame():
        html = focus_html() if focus() is not None else overview_html()
        return ui.tags.iframe(
            srcdoc=html,
            class_='rd-report-frame',
            sandbox='allow-scripts',
            title='Interactive dotplot report',
        )

    @render.ui
    def plot_area():
        if result() is None:
            return None
        pair = focus()
        toolbar = []
        if pair is not None:
            toolbar.append(
                ui.input_action_button(
                    'back_overview',
                    'Back to overview',
                    class_='btn-primary btn-sm rd-back-btn',
                )
            )
            toolbar.append(ui.span(f'{pair[0]} vs {pair[1]}', class_='rd-focus-label'))
        if input.interactive():
            toolbar.append(
                ui.div(
                    ui.tags.b('Navigate: '),
                    'scroll = pan up/down · Shift+scroll = pan left/right · '
                    'Cmd/Ctrl+scroll = zoom · drag = zoom to region · '
                    'click panel = focus · click match = details · '
                    'Esc = reset · double-click panel = standalone view',
                    class_='rd-nav-hint',
                )
            )
            body = ui.output_ui('report_frame')
        else:
            body = ui.output_plot('dotplot', height='72vh')
        return ui.div(
            ui.div(*toolbar, class_='rd-plot-toolbar') if toolbar else None,
            body,
            class_='rd-plot-area',
        )

    @reactive.effect
    @reactive.event(input.panel_dblclick)
    def _drill_down():
        if result() is None or focus() is not None:
            # No result yet, or already in the standalone single-pair view
            # (whose only panel is (0, 0) and must not remap).
            return
        info = input.panel_dblclick()
        lay = layout()
        try:
            q, t = panel_pair(
                lay['query_names'],
                lay['target_names'],
                int(info['row']),
                int(info['col']),
            )
        except (KeyError, TypeError, ValueError, IndexError):
            logger.warning('Ignoring unmappable panel_dblclick payload: %r', info)
            return
        focus.set((q, t))

    @reactive.effect
    @reactive.event(input.back_overview)
    def _back_to_overview():
        focus.set(None)

    # --- end W2 --------------------------------------------------------------

    @render.ui
    def status():
        if boot_error():
            return ui.div(boot_error(), class_='rd-status rd-status-error')
        if not ready():
            return ui.div(
                'Loading rusty-dot (first visit compiles the WASM runtime — '
                'this can take a few seconds)…',
                class_='rd-status',
            )
        if result() is None:
            return ui.div(
                'Upload two assemblies (or a PAF file) and press '
                '"Run comparison". Everything runs in your browser — '
                'files are never uploaded to a server. Large genomes are '
                'limited by browser memory (≈2 GB); bacterial/fungal-scale '
                'assemblies work best.',
                class_='rd-status',
            )
        return None

    @render.plot
    def dotplot():
        # --- W2: interactive plot --- (static fallback; honours drill-down)
        res = result()
        req(res)
        return make_figure(res, config(), layout(), pair=focus())

    def _figure_bytes(fmt: str) -> bytes:
        # --- W2: interactive plot --- (exports the currently shown view)
        res = result()
        req(res)
        import matplotlib.pyplot as plt

        fig = make_figure(res, config(), layout(), pair=focus())
        buf = io.BytesIO()
        fig.savefig(buf, format=fmt, bbox_inches='tight')
        plt.close(fig)
        return buf.getvalue()

    @render.download_button(filename='dotplot.svg')
    def dl_svg():
        yield _figure_bytes('svg')

    @render.download_button(filename='dotplot.pdf')
    def dl_pdf():
        yield _figure_bytes('pdf')

    @render.download_button(filename='alignment.paf')
    def dl_paf():
        res = result()
        req(res)
        kind, obj, _meta = res
        if kind == 'kmer':
            lines = obj.get_paf(group_pairs=[(QUERY_GROUP, TARGET_GROUP)], merge=True)
            yield '\n'.join(lines) + '\n'
        else:
            yield paf_text_from_alignment(obj)

    @render.download_button(filename='query_reordered.fasta')
    def dl_fasta():
        # --- W2: interactive plot --- (uses the same explicit layout as
        # the plot, so the exported FASTA matches what is displayed)
        res = result()
        req(res)
        kind, obj, meta = res
        lay = layout()
        order = lay['query_names']
        reverse = lay['reverse']
        if kind == 'kmer':
            with tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp) / 'query_reordered.fasta'
                obj.write_fasta(out, group=QUERY_GROUP, order=order, reverse=reverse)
                yield out.read_bytes()
        else:
            # Coordinate-only alignment (PAF import / external tool): export
            # from the query assembly sequences when they are available —
            # either attached to the result or from the sidebar upload.
            query = meta.get('query')
            if not isinstance(query, FastaInput):
                try:
                    query = _parse_upload(input.query_fasta, 'query')
                except ValueError:
                    query = None
            if query is None:
                ui.notification_show(
                    'Reordered FASTA export needs sequences — upload the '
                    'query assembly in the sidebar (PAF files carry '
                    'coordinates only).',
                    type='warning',
                    duration=8,
                )
                req(False)
            yield reordered_fasta_text(query.records, order, reverse)


app = App(app_ui, server, static_assets=APP_DIR / 'www')
