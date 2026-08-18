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
import uuid

from core.align import (
    AVAILABLE_METHODS,
    BIOWASM_TOOLS,
    METHOD_LABELS,
    MINIMAP2_PRESETS,
    alignment_from_tool_output,
    build_tool_args,
    fasta_text,
    paf_alignment_from_text,
    paf_text_from_alignment,
)
from core.cache import QUERY_GROUP, TARGET_GROUP, SessionCache
from core.export import reordered_fasta_text
from core.fasta import FastaInput, parse_fasta_bytes
from core.panels import panel_pair, resolve_orders
from core.state import ORDER_CHOICES, PlotConfig
from core.wheels import pick_wasm_wheel, runtime_platform_tag
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

        wheel = pick_wasm_wheel(
            list(APP_DIR.glob('wheels/*.whl')), runtime_platform_tag()
        )
        logger.info('Installing bundled wheel %s', wheel.name)
        await micropip.install(f'emfs:{wheel}')
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
        # PAF import is an input *mode* (the input_mode radio), not an
        # alignment method — with a PAF there is nothing left to align.
        if key != 'paf'
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


_HIDE_REPORT_HEADER_CSS = '<style>#rd-header{display:none}</style>'


def strip_report_header(html: str) -> str:
    """Hide the report's built-in title/navigation header for embedding.

    The standalone HTML report carries its own header with navigation hints
    (``#rd-header`` in ``_html/template.html``); inside the app that
    duplicates the app-level hint bar, so the embed hides it with CSS.
    Standalone ``to_html`` exports are unaffected.

    Parameters
    ----------
    html : str
        Full HTML report text produced by ``DotPlotter.to_html``.

    Returns
    -------
    str
        The report with a header-hiding style spliced in before ``</head>``
        (prepended if no closing head tag is found).
    """
    idx = html.find('</head>')
    if idx == -1:
        return _HIDE_REPORT_HEADER_CSS + html
    return html[:idx] + _HIDE_REPORT_HEADER_CSS + html[idx:]


# --- end W2 ------------------------------------------------------------------


app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h5('Input'),
        ui.input_radio_buttons(
            'input_mode',
            None,
            {'fasta': 'Assemblies (FASTA)', 'paf': 'Alignment (PAF)'},
            selected='fasta',
            inline=True,
        ),
        ui.panel_conditional(
            "input.input_mode === 'fasta'",
            ui.input_file(
                'query_fasta',
                'Query assembly (FASTA / .gz)',
                accept=['.fa', '.fasta', '.fna', '.gz'],
            ),
            ui.input_checkbox('self_align', 'Align assembly to itself', False),
            ui.panel_conditional(
                '!input.self_align',
                ui.input_file(
                    'target_fasta',
                    'Target / reference assembly (FASTA / .gz)',
                    accept=['.fa', '.fasta', '.fna', '.gz'],
                ),
            ),
            ui.input_select('method', 'Alignment method', choices=_method_choices()),
            ui.panel_conditional(
                "input.method === 'kmer'",
                ui.input_slider('k', 'k-mer size', min=8, max=64, value=21, step=1),
                ui.input_checkbox('merge', 'Merge adjacent matches', True),
            ),
            # --- W1: biowasm aligner options ---
            ui.panel_conditional(
                "input.method === 'minimap2'",
                ui.input_select(
                    'mm2_preset',
                    'Preset (-x)',
                    choices={p: p for p in MINIMAP2_PRESETS},
                    selected='asm20',
                ),
                ui.input_numeric(
                    'mm2_k', 'k-mer size (-k, 0 = preset default)', 0, min=0, max=28
                ),
                ui.input_numeric(
                    'mm2_w', 'Minimizer window (-w, 0 = preset default)', 0, min=0
                ),
            ),
            ui.panel_conditional(
                "input.method === 'lastz'",
                ui.input_numeric('lastz_step', 'Seed step (--step)', 1, min=1),
                ui.input_checkbox('lastz_gapped', 'Gapped extension', True),
                ui.input_checkbox(
                    'lastz_notransition', 'Exact seeds only (--notransition)', False
                ),
            ),
            ui.panel_conditional(
                "input.method === 'nucmer'",
                ui.input_numeric('nucmer_l', 'Min match length (-l)', 20, min=1),
                ui.input_numeric('nucmer_c', 'Min cluster length (-c)', 65, min=1),
                ui.input_checkbox(
                    'nucmer_maxmatch', 'Use all matches (--maxmatch)', False
                ),
            ),
            # --- end W1 ---
        ),
        ui.panel_conditional(
            "input.input_mode === 'paf'",
            ui.input_file('paf_file', 'PAF file', accept=['.paf', '.txt']),
            ui.input_file(
                'paf_query_fasta',
                'Query assembly (optional — enables the reordered-FASTA download)',
                accept=['.fa', '.fasta', '.fna', '.gz'],
            ),
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
        # Identity colouring only makes sense for tool/PAF alignments (k-mer
        # matches are always 100% identity); output.result_kind is a hidden
        # text output that tracks the current result.
        ui.panel_conditional(
            "output.result_kind === 'paf'",
            ui.input_checkbox('color_by_identity', 'Colour by % identity', False),
            ui.panel_conditional(
                'input.color_by_identity',
                ui.input_select(
                    'identity_palette',
                    'Identity palette',
                    choices=['viridis', 'plasma', 'cividis', 'coolwarm'],
                ),
            ),
        ),
        ui.input_numeric('min_length', 'Min match length (bp)', 0, min=0),
        ui.input_numeric('dot_size', 'Line width', 0.5, min=0.1, max=5, step=0.1),
        ui.hr(),
        ui.h5('Downloads'),
        ui.output_ui('downloads'),
        width=320,
    ),
    ui.div(ui.output_text('result_kind'), class_='rd-hidden'),
    ui.output_ui('status'),
    # --- W2: interactive plot ---
    ui.output_ui('plot_area'),
    ui.head_content(
        ui.include_css(APP_DIR / 'www' / 'app.css'),
        ui.include_js(APP_DIR / 'www' / 'bridge.js'),
    ),
    # W1: biowasm aligner bridge (Aioli loaded lazily from the CDN on use).
    ui.head_content(ui.include_js(APP_DIR / 'www' / 'aligners.js', method='inline')),
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

    def _parse_inputs(progress=None) -> tuple[FastaInput, FastaInput]:
        """Parse the query (and target, or reuse query when self-aligning)."""
        if progress is not None:
            progress.set(0, message='Parsing query assembly…')
        query = _parse_upload(input.query_fasta, 'query')
        if input.self_align():
            return query, query
        if progress is not None:
            progress.set(1, message='Parsing target assembly…')
        return query, _parse_upload(input.target_fasta, 'target')

    @reactive.effect
    @reactive.event(input.run)
    async def _run():
        req(ready())
        mode = input.input_mode()
        try:
            if mode == 'paf':
                with ui.Progress(min=0, max=3) as progress:
                    progress.set(1, message='Parsing PAF…')
                    files = input.paf_file()
                    if not files:
                        raise ValueError('Please upload a PAF file.')
                    text = Path(files[0]['datapath']).read_text()
                    alignment = paf_alignment_from_text(text)
                    progress.set(3, message=f'{len(alignment)} alignment(s) loaded')
                    result.set(('paf', alignment, {'query': None, 'target': None}))
                return
            method = input.method()
            if method not in AVAILABLE_METHODS:
                ui.notification_show(
                    f'{METHOD_LABELS[method]} is not available yet.',
                    type='warning',
                )
                return
            if method != 'kmer':
                return  # biowasm tools are handled by _run_biowasm
            with ui.Progress(min=0, max=4) as progress:
                query, target = _parse_inputs(progress)
                progress.set(
                    2,
                    message=(
                        f'Building k-mer index (k={input.k()}, '
                        f'{len(query.records)}×{len(target.records)} contigs)…'
                    ),
                )
                index = cache.kmer_index(input.k(), query, target, merge=input.merge())
                progress.set(3, message='Rendering dotplot…')
                result.set(('kmer', index, {'query': query, 'target': target}))
                progress.set(4, message='Done')
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=8)

    # --- W1: biowasm aligners ---
    # Runs minimap2 / LASTZ / nucmer in a browser WebWorker via biowasm
    # (Aioli).  Python builds the CLI args and plain-FASTA payload, sends a
    # 'rd_run_aligner' custom message to www/aligners.js, and receives the
    # tool's text output back through the 'aligner_result' input.

    # All biowasm tools share one 2 GB wasm heap cap; warn above ~200 MB.
    _BIOWASM_SIZE_WARN = 200 * 1024 * 1024
    _NOTIF_ID = 'rd_aligner_progress'
    # In-flight request: {'request_id', 'method', 'params', 'query', 'target'}
    aligner_pending = reactive.value(None)

    def _tool_params(method: str) -> dict:
        if method == 'minimap2':
            return {
                'preset': input.mm2_preset(),
                'k': int(input.mm2_k() or 0),
                'w': int(input.mm2_w() or 0),
            }
        if method == 'lastz':
            return {
                'step': int(input.lastz_step() or 1),
                'gapped': bool(input.lastz_gapped()),
                'notransition': bool(input.lastz_notransition()),
            }
        return {
            'l': int(input.nucmer_l() or 20),
            'c': int(input.nucmer_c() or 65),
            'maxmatch': bool(input.nucmer_maxmatch()),
        }

    @reactive.effect
    @reactive.event(input.run)
    async def _run_biowasm():
        req(ready())
        # Every Run click supersedes any in-flight aligner request: a late
        # 'aligner_result' for an invalidated request_id must not clobber
        # whatever the user asked for afterwards (including k-mer/PAF runs,
        # which the sibling _run effect handles on this same event).
        if aligner_pending() is not None:
            aligner_pending.set(None)
            ui.notification_remove(_NOTIF_ID)
        if input.input_mode() != 'fasta':
            return
        method = input.method()
        if method not in BIOWASM_TOOLS:
            return
        try:
            query, target = _parse_inputs()
            params = _tool_params(method)
            args = build_tool_args(method, params)
            query_text = fasta_text(query.records)
            target_text = fasta_text(target.records)
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=8)
            return
        cached = cache.get_paf(method, params, query.digest, target.digest)
        if cached is not None:
            logger.info('%s alignment cache hit', method)
            result.set(('paf', cached, {'query': query, 'target': target}))
            return
        if len(query_text) + len(target_text) > _BIOWASM_SIZE_WARN:
            ui.notification_show(
                'Combined input exceeds ~200 MB — browser aligners share a '
                '2 GB memory cap and may fail on inputs this large.',
                type='warning',
                duration=12,
            )
        request_id = uuid.uuid4().hex
        aligner_pending.set(
            {
                'request_id': request_id,
                'method': method,
                'params': params,
                'query': query,
                'target': target,
            }
        )
        ui.notification_show(
            f'Running {METHOD_LABELS[method]} in your browser (the first '
            'run downloads the tool from the biowasm CDN)…',
            id=_NOTIF_ID,
            duration=None,
        )
        await session.send_custom_message(
            'rd_run_aligner',
            {
                'tool': method,
                'args': args,
                'query_fasta': query_text,
                'target_fasta': target_text,
                'request_id': request_id,
            },
        )

    @reactive.effect
    @reactive.event(input.aligner_result)
    def _on_aligner_result():
        res = input.aligner_result()
        info = aligner_pending()
        if not res or not info or res.get('request_id') != info['request_id']:
            return  # stale or unsolicited result
        aligner_pending.set(None)
        ui.notification_remove(_NOTIF_ID)
        method = info['method']
        if res.get('error'):
            ui.notification_show(
                f'{METHOD_LABELS[method]} failed: {res["error"]}',
                type='error',
                duration=12,
            )
            return
        try:
            alignment = alignment_from_tool_output(method, res.get('output') or '')
        except ValueError as exc:
            ui.notification_show(
                f'Could not parse {METHOD_LABELS[method]} output: {exc}',
                type='error',
                duration=12,
            )
            return
        cache.put_paf(
            method,
            info['params'],
            alignment,
            info['query'].digest,
            info['target'].digest,
        )
        ui.notification_show(
            f'{METHOD_LABELS[method]}: {len(alignment)} alignment(s).',
            type='message',
            duration=5,
        )
        result.set(
            ('paf', alignment, {'query': info['query'], 'target': info['target']})
        )

    @reactive.effect
    @reactive.event(input.aligner_progress)
    def _on_aligner_progress():
        prog = input.aligner_progress()
        info = aligner_pending()
        if not prog or not info or prog.get('request_id') != info['request_id']:
            return  # stale or pre-warm activity, not the run in flight
        stage_msgs = {
            'loading-aioli': 'Loading the biowasm runtime…',
            'initialising-tool': (
                f'Downloading {METHOD_LABELS[info["method"]]} from the '
                'biowasm CDN (first use only)…'
            ),
            'aligning': f'Running {METHOD_LABELS[info["method"]]}…',
            'reading-output': 'Reading alignment output…',
        }
        msg = stage_msgs.get(prog.get('stage'))
        if msg:
            ui.notification_show(msg, id=_NOTIF_ID, duration=None)

    # --- end W1: biowasm aligners ---

    @reactive.calc
    def config() -> PlotConfig:
        return PlotConfig(
            contig_order=input.contig_order(),
            auto_reverse=input.auto_reverse(),
            hide_internal_axes=input.hide_internal_axes(),
            nature=input.nature(),
            dot_size=input.dot_size() or 0.5,
            min_length=int(input.min_length() or 0),
            color_by_identity=bool(input.color_by_identity()),
            identity_palette=input.identity_palette() or 'viridis',
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
        if kind == 'kmer':
            # k-mer matches are exact: identity is uniformly 100%, so
            # identity colouring would be a misleading single-colour plot.
            kwargs['color_by_identity'] = False
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
        return strip_report_header(inject_panel_bridge(html))

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

    @output(suspend_when_hidden=False)
    @render.text
    def result_kind():
        # Hidden output driving panel_conditional visibility (e.g. the
        # identity-colouring controls, shown only for tool/PAF results).
        # suspend_when_hidden=False: the output lives in a display:none div,
        # which Shiny would otherwise never update.
        res = result()
        return res[0] if res else ''

    def _has_sequences(res) -> bool:
        """Whether a reordered-FASTA export is possible for this result."""
        kind, _obj, meta = res
        return (
            kind == 'kmer'
            or isinstance(meta.get('query'), FastaInput)
            or bool(input.paf_query_fasta())
        )

    @render.ui
    def downloads():
        res = result()
        if res is None:
            return ui.div(
                ui.tags.button(
                    'Plot (SVG)', class_='btn rd-dl-disabled', disabled=True
                ),
                ui.tags.button(
                    'Plot (PDF)', class_='btn rd-dl-disabled', disabled=True
                ),
                ui.tags.button(
                    'Alignment (PAF)', class_='btn rd-dl-disabled', disabled=True
                ),
                ui.tags.button(
                    'Reordered query (FASTA)',
                    class_='btn rd-dl-disabled',
                    disabled=True,
                ),
                ui.div('Run a comparison first.', class_='rd-dl-note'),
            )
        parts = [
            ui.download_button('dl_svg', 'Plot (SVG)'),
            ui.download_button('dl_pdf', 'Plot (PDF)'),
            ui.download_button('dl_paf', 'Alignment (PAF)'),
        ]
        if _has_sequences(res):
            parts.append(ui.download_button('dl_fasta', 'Reordered query (FASTA)'))
        else:
            parts += [
                ui.tags.button(
                    'Reordered query (FASTA)',
                    class_='btn rd-dl-disabled',
                    disabled=True,
                ),
                ui.div(
                    'FASTA export needs sequences — upload the query '
                    'assembly in the sidebar.',
                    class_='rd-dl-note',
                ),
            ]
        return ui.div(*parts)

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
                    query = _parse_upload(input.paf_query_fasta, 'query')
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
