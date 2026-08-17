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
from core.fasta import FastaInput, parse_fasta_bytes
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
            ui.input_checkbox('nucmer_maxmatch', 'Use all matches (--maxmatch)', False),
        ),
        # --- end W1 ---
        ui.input_action_button('run', 'Run comparison', class_='btn-primary'),
        ui.hr(),
        ui.h5('Plot options'),
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
    ui.output_plot('dotplot', height='75vh'),
    ui.head_content(ui.include_css(APP_DIR / 'www' / 'app.css')),
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
        method = input.method()
        if method not in BIOWASM_TOOLS:
            return
        try:
            query = _parse_upload(input.query_fasta, 'query')
            target = _parse_upload(input.target_fasta, 'target')
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
        )

    def make_figure(res, cfg: PlotConfig):
        from rusty_dot import DotPlotter
        from rusty_dot.style import nature_style

        kind, obj, _meta = res
        kwargs = cfg.plot_kwargs()
        if kind == 'kmer':
            kwargs.update(query_group=QUERY_GROUP, target_group=TARGET_GROUP)
        else:
            # Keep PAF queries on rows and targets on columns (the default
            # would plot the union of all names on both axes).
            kwargs.update(query_names=obj.query_names, target_names=obj.target_names)
        style = nature_style() if cfg.nature else contextlib.nullcontext()
        with style:
            return DotPlotter(obj).plot(**kwargs)

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
        res = result()
        req(res)
        return make_figure(res, config())

    def _figure_bytes(fmt: str) -> bytes:
        res = result()
        req(res)
        fig = make_figure(res, config())
        buf = io.BytesIO()
        fig.savefig(buf, format=fmt, bbox_inches='tight')
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
        res = result()
        req(res)
        kind, obj, meta = res
        if kind != 'kmer':
            ui.notification_show(
                'Reordered FASTA export needs sequences — run the k-mer '
                'method (PAF files carry coordinates only).',
                type='warning',
                duration=8,
            )
            req(False)
        cfg = config()
        query: FastaInput = meta['query']
        if cfg.contig_order == 'length':
            order = [
                name
                for name, _ in sorted(
                    query.records, key=lambda r: len(r[1]), reverse=True
                )
            ]
            reverse: set[str] = set()
        elif cfg.contig_order == 'colinearity':
            order, _t_order = obj.reorder_contigs(
                query_group=QUERY_GROUP, target_group=TARGET_GROUP
            )
            reverse = obj.reversed_contigs(QUERY_GROUP) if cfg.auto_reverse else set()
        else:
            order = query.names
            reverse = set()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / 'query_reordered.fasta'
            obj.write_fasta(out, group=QUERY_GROUP, order=order, reverse=reverse)
            yield out.read_bytes()


app = App(app_ui, server, static_assets=APP_DIR / 'www')
