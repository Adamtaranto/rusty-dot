"""rusty-dot assembly comparison — browser app (Shiny for Python / Shinylive).

Runs entirely client-side: uploads never leave the browser.  Under Pyodide
the rusty-dot wasm wheel bundled in ``wheels/`` is installed at startup;
run natively (``shiny run app/app.py``) it uses the installed rusty-dot.
"""

from __future__ import annotations

from html import escape as _esc
import importlib.util
import io
import logging
from pathlib import Path
import sys
import tempfile
import time
import uuid

from core.align import (
    AVAILABLE_METHODS,
    BIOWASM_TOOLS,
    METHOD_LABELS,
    MINIMAP2_PRESET_DEFAULTS,
    MINIMAP2_PRESETS,
    alignment_from_tool_output,
    build_tool_args,
    fasta_text,
    paf_alignment_from_text,
    paf_text_from_alignment,
)
from core.annotation_colors import (
    assign_shared_colors,
    color_map_for,
    display_name,
    normalise_type,
)
from core.annotation_state import (
    ANNOTATION_ROLES,
    apply_annotation_config,
    apply_feature_overrides,
    build_feature_rows,
    merge_annotations,
    replace_source,
    type_slug_map,
)
from core.cache import QUERY_GROUP, TARGET_GROUP, SessionCache
from core.export import reordered_fasta_text
from core.fasta import FastaInput, content_digest, parse_fasta_bytes
from core.genbank import parse_genbank_bytes
from core.panels import (
    filter_by_min_length,
    has_self_pair,
    nav_tips,
    panel_pair,
    resolve_orders,
)
from core.state import ORDER_CHOICES, PlotConfig
from core.validate import validate_annotation_names, validate_query_names
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


def _lbl(text: str, tip: str):
    """Input label with a hoverable info icon explaining the setting.

    Parameters
    ----------
    text : str
        The visible label text.
    tip : str
        One-line explanation shown in a tooltip on the icon.

    Returns
    -------
    htmltools.Tag
        Label span suitable as any input's ``label`` argument.
    """
    return ui.span(
        text,
        ' ',
        ui.tooltip(ui.tags.span('ⓘ', class_='rd-info'), tip, placement='right'),
    )


# --- W2: interactive plot ---------------------------------------------------
# Script injected into the generated HTML report before embedding: posts a
# message to the parent app window when a dotplot panel is double-clicked.
# app/www/bridge.js listens for it and forwards to the 'panel_dblclick' input.
_PANEL_DBLCLICK_JS = """
<script>
(function () {
  'use strict';
  // Tell report.js that double-click drills down here, so it defers the
  // single-click focus zoom briefly and a double-click cancels it —
  // without this, double-clicking a panel zooms in and then swaps views.
  window.RD_DBLCLICK_DRILLDOWN = true;
  document.addEventListener('dblclick', function (ev) {
    // Walk ancestors testing the exact panel id.  closest() with a prefix
    // selector would stop at the axes background group, whose id merely
    // starts the same way, and the drill-down would silently never fire.
    var node = ev.target;
    var m = null;
    while (node && node.nodeType === 1) {
      m = /^rd-panel-(\\d+)-(\\d+)$/.exec(node.id || '');
      if (m) { break; }
      node = node.parentNode;
    }
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


def debounce(delay_secs: float):
    """Debounce a reactive calculation (standard Shiny recipe).

    Wraps a ``reactive.calc``-decorated function so downstream consumers
    only see a new value once the source has been stable for *delay_secs*.
    Used for numeric inputs whose spinner arrows would otherwise trigger a
    full plot re-render for every 0.1 step.

    Parameters
    ----------
    delay_secs : float
        Quiet period before the new value propagates.

    Returns
    -------
    Callable
        Decorator for a ``reactive.calc`` function.
    """

    def wrapper(f):
        when = reactive.value(None)
        trigger = reactive.value(0)

        @reactive.effect(priority=102)
        def _primer():
            try:
                f()
            except Exception:  # noqa: BLE001 - value read only to register deps
                pass
            with reactive.isolate():
                when.set(time.monotonic() + delay_secs)

        @reactive.effect(priority=101)
        def _timer():
            deadline = when()
            if deadline is None:
                return
            now = time.monotonic()
            if now >= deadline:
                when.set(None)
                with reactive.isolate():
                    trigger.set(trigger() + 1)
            else:
                reactive.invalidate_later(deadline - now)

        @reactive.calc
        @reactive.event(trigger, ignore_none=False)
        def debounced():
            with reactive.isolate():
                return f()

        return debounced

    return wrapper


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
            {
                'fasta': 'Assemblies (FASTA)',
                'genbank': 'Assemblies (GenBank)',
                'paf': 'Alignment (PAF)',
            },
            selected='fasta',
            inline=True,
        ),
        # Everything except the PAF branch shares the self-align checkbox,
        # the method selector and the per-method parameter panels; only the
        # upload widgets differ between FASTA and GenBank.
        ui.panel_conditional(
            "input.input_mode !== 'paf'",
            ui.panel_conditional(
                "input.input_mode === 'fasta'",
                ui.input_file(
                    'query_fasta',
                    'Query assembly (FASTA / .gz)',
                    accept=['.fa', '.fasta', '.fna', '.gz'],
                ),
            ),
            ui.panel_conditional(
                "input.input_mode === 'genbank'",
                ui.input_file(
                    'query_gbk',
                    _lbl(
                        'Query assembly (GenBank / .gz)',
                        'Sequences and annotations are read from the same '
                        'file; any GFF uploaded below is merged with them.',
                    ),
                    accept=['.gb', '.gbk', '.gbff', '.genbank', '.gz'],
                ),
            ),
            ui.input_checkbox(
                'self_align',
                _lbl(
                    'Align assembly to itself',
                    'Compare one assembly against itself to reveal repeats '
                    'and segmental duplications — no second upload needed. '
                    'Only the query assembly is used.',
                ),
                False,
            ),
            ui.panel_conditional(
                "!input.self_align && input.input_mode === 'fasta'",
                ui.input_file(
                    'target_fasta',
                    'Target / reference assembly (FASTA / .gz)',
                    accept=['.fa', '.fasta', '.fna', '.gz'],
                ),
            ),
            ui.panel_conditional(
                "!input.self_align && input.input_mode === 'genbank'",
                ui.input_file(
                    'target_gbk',
                    'Target / reference assembly (GenBank / .gz)',
                    accept=['.gb', '.gbk', '.gbff', '.genbank', '.gz'],
                ),
            ),
            ui.input_select(
                'method',
                _lbl(
                    'Alignment method',
                    'k-mer: fast exact matching, fully offline; minimap2 / '
                    'nucmer: full aligners compiled to WebAssembly (fetched '
                    'from the biowasm CDN at runtime).',
                ),
                choices=_method_choices(),
            ),
            ui.panel_conditional(
                "input.method === 'kmer'",
                ui.input_slider(
                    'k',
                    _lbl(
                        'k-mer size',
                        'Exact-match seed length; larger k gives fewer, more '
                        'specific matches.',
                    ),
                    min=8,
                    max=64,
                    value=21,
                    step=1,
                ),
                ui.input_checkbox(
                    'merge',
                    _lbl(
                        'Merge adjacent matches',
                        'Join runs of adjacent k-mer hits into longer '
                        'diagonal segments.',
                    ),
                    True,
                ),
                # Compute-time filter: repeat-rich assembly pairs can produce
                # millions of short match blocks that exhaust browser memory;
                # dropping them natively keeps real genomes workable.
                ui.input_numeric(
                    'kmer_min_block',
                    _lbl(
                        'Min match block (bp, 0 = keep all)',
                        'Drop merged match blocks shorter than this before '
                        'they reach the plot; keeps repeat-rich genome pairs '
                        'from producing millions of records.',
                    ),
                    50,
                    min=0,
                ),
            ),
            # --- W1: biowasm aligner options ---
            ui.panel_conditional(
                "input.method === 'minimap2'",
                ui.input_select(
                    'mm2_preset',
                    _lbl(
                        'Preset (-x)',
                        'Divergence preset: asm5 ≲1% sequence divergence, '
                        'asm10 ≲5%, asm20 ≲10% (safe default across '
                        'strains/isolates).',
                    ),
                    choices={p: p for p in MINIMAP2_PRESETS},
                    selected='asm20',
                ),
                # k/w/m are pre-filled with the selected preset's actual
                # values (asm20 defaults) and refreshed by _sync_mm2_defaults
                # whenever the preset changes.
                ui.input_numeric(
                    'mm2_k',
                    _lbl(
                        'k-mer size (-k)',
                        'Minimizer k-mer length; smaller values increase '
                        'sensitivity for diverged repeats.',
                    ),
                    MINIMAP2_PRESET_DEFAULTS['asm20']['k'],
                    min=1,
                    max=28,
                ),
                ui.input_numeric(
                    'mm2_w',
                    _lbl(
                        'Minimizer window (-w)',
                        'Minimizer window size; smaller windows increase seed density.',
                    ),
                    MINIMAP2_PRESET_DEFAULTS['asm20']['w'],
                    min=1,
                ),
                ui.input_numeric(
                    'mm2_m',
                    _lbl(
                        'Min chaining score (-m)',
                        'Discard chains scoring below this; raise (e.g. 200) '
                        'to filter weak or noisy repeat matches.',
                    ),
                    MINIMAP2_PRESET_DEFAULTS['asm20']['m'],
                    min=1,
                ),
                ui.input_checkbox(
                    'mm2_c',
                    _lbl(
                        'Base-level alignment (-c)',
                        'Compute exact CIGARs and identity tags (cg/NM/de): '
                        'enables gap-compressed identity colouring and '
                        'aligned-sequence display, but is slower on large '
                        'assemblies.',
                    ),
                    False,
                ),
                ui.input_checkbox(
                    'mm2_p',
                    _lbl(
                        'Retain all chains (-P)',
                        'Keep every chain instead of dropping secondary '
                        'matches, so each repeat copy maps — larger output.',
                    ),
                    False,
                ),
                ui.panel_conditional(
                    'input.self_align',
                    ui.input_checkbox(
                        'mm2_d',
                        _lbl(
                            'Skip self-diagonal matches (-D)',
                            'Drop the trivial full-length match of each '
                            'contig against itself on the main diagonal.',
                        ),
                        False,
                    ),
                ),
            ),
            ui.panel_conditional(
                "input.method === 'nucmer'",
                # mummer's -l 20 -c 65 defaults are tuned for small regions;
                # for assembly-scale dotplots they mostly add noise clusters
                # and minutes of runtime.  -l 100 -c 200 is the conventional
                # whole-genome comparison setting.
                ui.input_numeric(
                    'nucmer_l',
                    _lbl(
                        'Min match length (-l)',
                        'Minimum exact-match anchor length; larger values '
                        'run faster with fewer spurious matches.',
                    ),
                    100,
                    min=1,
                ),
                ui.input_numeric(
                    'nucmer_c',
                    _lbl(
                        'Min cluster length (-c)',
                        'Minimum combined anchor length for a cluster to be '
                        'reported as an alignment.',
                    ),
                    200,
                    min=1,
                ),
                ui.input_checkbox(
                    'nucmer_maxmatch',
                    _lbl(
                        'Use all matches (--maxmatch)',
                        'Use every anchor match regardless of uniqueness — '
                        'required to see all copies of repetitive regions.',
                    ),
                    False,
                ),
                ui.input_checkbox(
                    'nucmer_nosimplify',
                    _lbl(
                        'Keep repeat-induced alignments (--nosimplify)',
                        'Keep shadowed clusters instead of simplifying them '
                        'away — needed to find inexact repeats in '
                        'self-alignments.',
                    ),
                    False,
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
        ui.input_checkbox(
            'interactive',
            _lbl(
                'Interactive plot (zoom & drill-down)',
                'Render a zoomable HTML report with clickable matches '
                'instead of a static image.',
            ),
            True,
        ),
        ui.input_select(
            'contig_order',
            _lbl(
                'Contig order',
                'How contigs are arranged along the axes: upload order, '
                'longest first, or reordered to maximise colinearity with '
                'the other assembly.',
            ),
            choices=ORDER_CHOICES,
        ),
        ui.input_numeric(
            'min_contig_len',
            _lbl(
                'Min contig length (bp, 0 = keep all)',
                'Leave contigs shorter than this out of the plot. One panel '
                'per contig means a few chromosomes can be buried under '
                'hundreds of short scaffolds. Excluded contigs are still '
                'written to the reordered-FASTA download.',
            ),
            value=0,
            min=0,
            step=1000,
        ),
        ui.input_checkbox(
            'auto_reverse',
            _lbl(
                'Auto-flip reversed contigs',
                'Display contigs that align mostly in reverse on their '
                'reverse strand so synteny reads as a forward diagonal.',
            ),
            False,
        ),
        ui.input_checkbox(
            'hide_internal_axes',
            _lbl(
                'Hide internal axes',
                'Remove the internal panel borders and ticks so the '
                'all-vs-all grid reads as one continuous plot.',
            ),
            False,
        ),
        # Identity colouring only makes sense for tool/PAF alignments (k-mer
        # matches are always 100% identity); output.result_kind is a hidden
        # text output that tracks the current result.
        ui.panel_conditional(
            "output.result_kind === 'paf'",
            ui.input_checkbox(
                'color_by_identity',
                _lbl(
                    'Colour by % identity',
                    'Colour each match by its percent identity instead of '
                    'forward/reverse strand (aligner and PAF results only).',
                ),
                False,
            ),
            ui.panel_conditional(
                'input.color_by_identity',
                ui.input_select(
                    'identity_palette',
                    _lbl(
                        'Identity palette',
                        'Colour map used for the identity scale.',
                    ),
                    choices=['viridis', 'plasma', 'cividis', 'coolwarm'],
                ),
            ),
        ),
        ui.input_numeric(
            'min_length',
            _lbl(
                'Min match length (bp)',
                'Hide matches shorter than this many base pairs — applied '
                'instantly, without re-rendering.',
            ),
            0,
            min=0,
        ),
        ui.input_numeric(
            'dot_size',
            _lbl(
                'Line width',
                'Stroke width of the match segments — applied instantly, '
                'without re-rendering.',
            ),
            0.5,
            min=0.1,
            max=5,
            step=0.1,
        ),
        ui.hr(),
        # --- GFF annotations -------------------------------------------------
        ui.h5('Annotations (GFF3)'),
        # Same query-then-target order as the assembly uploads above.
        ui.input_file(
            'query_gff',
            'Query annotations (.gff / .gff3 / .gz)',
            accept=['.gff', '.gff3', '.gz'],
        ),
        # Hidden alongside the target assembly when self-aligning: both axes
        # are then the query assembly, so there is no target to annotate.
        # PAF input has no self-align notion but does have both roles, so it
        # keeps the upload regardless of a stale checkbox value.
        ui.panel_conditional(
            "input.input_mode === 'paf' || !input.self_align",
            ui.input_file(
                'target_gff',
                'Target annotations (.gff / .gff3 / .gz)',
                accept=['.gff', '.gff3', '.gz'],
            ),
        ),
        ui.panel_conditional(
            "output.gff_mode === 'plain' || output.gff_mode === 'self'",
            ui.input_action_button(
                'clear_gff',
                'Clear annotations',
                class_='btn-outline-secondary btn-sm',
            ),
        ),
        # Static, so pressing Run never destroys them.  Only their visibility
        # is reactive (output.gff_mode); panel_conditional toggles CSS display
        # and leaves the inputs bound, so they can be read unconditionally.
        ui.panel_conditional(
            "output.gff_mode === 'plain' || output.gff_mode === 'self'",
            ui.input_checkbox('gff_tracks', 'Side tracks in focused pair view', True),
        ),
        ui.panel_conditional(
            "output.gff_mode === 'self'",
            ui.tooltip(
                ui.input_checkbox('gff_diagonal', 'Shade features on diagonal', True),
                'Draws each feature as a square on the diagonal of '
                'self-comparison panels.',
            ),
        ),
        ui.output_ui('gff_controls'),
        ui.hr(),
        ui.h5('Downloads'),
        ui.output_ui('downloads'),
        width=320,
    ),
    # Pure-CSS busy pill: html.shiny-busy (set by Shiny while any
    # computation runs) makes it visible, with spinner + pulse animations.
    ui.div('Processing…', class_='rd-busy-pill'),
    # Fixed memory note (bottom-right; hidden while the readout is empty,
    # e.g. on native runs where the wasm heap does not exist).
    ui.div(ui.output_text('app_memory'), class_='rd-mem-fixed'),
    ui.div(ui.output_text('result_kind'), class_='rd-hidden'),
    ui.div(ui.output_text('gff_mode'), class_='rd-hidden'),
    ui.output_ui('status'),
    # --- W2: interactive plot ---
    ui.output_ui('plot_area'),
    ui.output_ui('aligner_log_ui'),
    ui.head_content(
        ui.include_css(APP_DIR / 'www' / 'app.css'),
        ui.include_js(APP_DIR / 'www' / 'bridge.js'),
        # Custom Shiny binding for native <input type="color"> pickers.
        ui.include_js(APP_DIR / 'www' / 'color-input.js', method='inline'),
        # Delegated events for the drill-down Annotations table (one
        # listener instead of ~1200 Shiny-bound inputs).
        ui.include_js(APP_DIR / 'www' / 'feature-table.js', method='inline'),
        # Hold the sidebar's scroll position across dynamic-UI re-renders.
        ui.include_js(APP_DIR / 'www' / 'sidebar-scroll.js', method='inline'),
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

    def _parse_seq_upload(role: str) -> FastaInput:
        """Parse one assembly upload for the current input mode.

        GenBank carries its annotations in the same file, so parsing also
        registers them as an annotation source for *role* — merged with
        any GFF the user uploads separately.
        """
        if input.input_mode() != 'genbank':
            fasta_input = input.query_fasta if role == 'query' else input.target_fasta
            return _parse_upload(fasta_input, role)

        gbk_input = input.query_gbk if role == 'query' else input.target_gbk
        files = gbk_input()
        if not files:
            raise ValueError(f'Please upload a {role} assembly.')
        parsed = parse_genbank_bytes(Path(files[0]['datapath']).read_bytes())
        # digest is over the raw upload, so an unchanged file re-run keeps the
        # existing annotation source (and the user's annotation choices).
        _set_ann_source(
            role,
            'genbank',
            files[0]['name'],
            parsed.gff_text,
            key=(parsed.fasta.digest, files[0]['name']),
        )
        return parsed.fasta

    def _parse_inputs(progress=None) -> tuple[FastaInput, FastaInput]:
        """Parse the query (and target, or reuse query when self-aligning)."""
        if progress is not None:
            progress.set(0, message='Parsing query assembly…')
        query = _parse_seq_upload('query')
        if input.self_align():
            return query, query
        if progress is not None:
            progress.set(1, message='Parsing target assembly…')
        return query, _parse_seq_upload('target')

    # The canonical CSR k-mer index costs ~13-16 bytes/bp (both strands);
    # a real 74 Mb pair completes in the browser (verified), so the guard
    # sits at ~80 Mb combined under Pyodide's 2 GB heap.  Exceeding the
    # heap does not fail gracefully — the wasm allocator aborts and takes
    # the whole Python runtime down — so refuse up front with a way
    # forward, and warn earlier that big k-mer runs take minutes.
    _KMER_HARD_LIMIT = 80 * 1024 * 1024
    _KMER_WARN_LIMIT = 40 * 1024 * 1024

    def _check_kmer_memory(query: FastaInput, target: FastaInput) -> None:
        if sys.platform != 'emscripten':
            return  # native runs are bounded by system RAM, not the wasm heap
        total = query.total_length + (0 if target is query else target.total_length)
        if total > _KMER_HARD_LIMIT:
            raise ValueError(
                f'Combined assemblies are {total / 1e6:.0f} Mb — the k-mer '
                'index needs more memory than the browser allows above '
                '~80 Mb and would crash the app. Use minimap2 (or nucmer) '
                'for assemblies this size.'
            )
        if total > _KMER_WARN_LIMIT:
            ui.notification_show(
                'Large input for the k-mer method — expect a few minutes of '
                'processing; minimap2 is much faster at this scale.',
                type='warning',
                duration=10,
            )

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
                    query = None
                    if input.paf_query_fasta():
                        progress.set(2, message='Parsing query assembly…')
                        query = _parse_upload(input.paf_query_fasta, 'query')
                        for warning in validate_query_names(
                            query.names,
                            alignment.query_names,
                            alignment.target_names,
                        ):
                            ui.notification_show(warning, type='warning', duration=12)
                    progress.set(3, message=f'{len(alignment)} alignment(s) loaded')
                    result.set(('paf', alignment, {'query': query, 'target': None}))
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
                _check_kmer_memory(query, target)
                progress.set(
                    2,
                    message=(
                        f'Building k-mer index (k={input.k()}, '
                        f'{len(query.records)}×{len(target.records)} contigs)…'
                    ),
                )
                index = cache.kmer_index(
                    input.k(),
                    query,
                    target,
                    merge=input.merge(),
                    min_block_len=int(input.kmer_min_block() or 0),
                )
                progress.set(3, message='Rendering dotplot…')
                result.set(('kmer', index, {'query': query, 'target': target}))
                progress.set(4, message='Done')
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=8)

    # --- W1: biowasm aligners ---
    # Runs minimap2 / nucmer in a browser WebWorker via biowasm
    # (Aioli).  Python builds the CLI args and plain-FASTA payload, sends a
    # 'rd_run_aligner' custom message to www/aligners.js, and receives the
    # tool's text output back through the 'aligner_result' input.

    # All biowasm tools share one 2 GB wasm heap cap; warn above ~200 MB.
    _BIOWASM_SIZE_WARN = 200 * 1024 * 1024
    _NOTIF_ID = 'rd_aligner_progress'
    # Extra budget granted each time the user chooses to keep waiting.
    _TIMEOUT_EXTEND_MS = 300_000
    # In-flight request: {'request_id', 'method', 'params', 'query', 'target'}
    aligner_pending = reactive.value(None)
    # Rolling log of completed tool runs: [{'tool', 'cmd', 'stderr', 'error'}]
    aligner_log = reactive.value([])
    _ALIGNER_LOG_MAX = 10
    # Dataset digests already shipped to aligners.js this session.  FASTA
    # payloads are sent once per upload and referenced by digest afterwards,
    # so re-runs and tool switches never re-copy whole genomes across the
    # Pyodide/JS boundary.
    sent_datasets: set[str] = set()

    def _log_run(tool: str, cmd: str, stderr: str, error: str | None) -> None:
        entries = list(aligner_log())
        entries.append({'tool': tool, 'cmd': cmd, 'stderr': stderr, 'error': error})
        aligner_log.set(entries[-_ALIGNER_LOG_MAX:])

    async def _cancel_pending(reason: str) -> None:
        """Cancel the in-flight aligner run (JS drops/terminates it)."""
        info = aligner_pending()
        if info is None:
            return
        aligner_pending.set(None)
        ui.modal_remove()
        ui.notification_remove(_NOTIF_ID)
        await session.send_custom_message(
            'rd_cancel_aligner', {'request_id': info['request_id']}
        )
        ui.notification_show(
            f'{METHOD_LABELS[info["method"]]} run cancelled ({reason}).',
            type='message',
            duration=5,
        )

    @reactive.effect
    @reactive.event(input.method)
    async def _cancel_on_method_change():
        info = aligner_pending()
        if info is not None and info['method'] != input.method():
            await _cancel_pending('method changed')

    @reactive.effect
    @reactive.event(input.mm2_preset, ignore_init=True)
    def _sync_mm2_defaults():
        # Selecting a preset means adopting its parameters: refresh the
        # k/w/m inputs to the preset's actual values so what is shown is
        # always what runs.
        defaults = MINIMAP2_PRESET_DEFAULTS.get(input.mm2_preset())
        if defaults is None:
            return
        ui.update_numeric('mm2_k', value=defaults['k'])
        ui.update_numeric('mm2_w', value=defaults['w'])
        ui.update_numeric('mm2_m', value=defaults['m'])

    def _tool_params(method: str) -> dict:
        if method == 'minimap2':
            return {
                'preset': input.mm2_preset(),
                'k': int(input.mm2_k() or 0),
                'w': int(input.mm2_w() or 0),
                'm': int(input.mm2_m() or 0),
                'c': bool(input.mm2_c()),
                'P': bool(input.mm2_p()),
                # -D only applies when a sequence can meet itself; storing
                # the effective value keeps the cache key honest when the
                # checkbox stays ticked but self-align is turned off.
                'D': bool(input.mm2_d()) and bool(input.self_align()),
            }
        return {
            'l': int(input.nucmer_l() or 100),
            'c': int(input.nucmer_c() or 200),
            'maxmatch': bool(input.nucmer_maxmatch()),
            'nosimplify': bool(input.nucmer_nosimplify()),
        }

    async def _send_dataset(data: FastaInput) -> None:
        """Ship a parsed assembly to aligners.js once, keyed by digest."""
        if data.digest in sent_datasets:
            return
        await session.send_custom_message(
            'rd_mount_fasta',
            {'dataset_id': data.digest, 'text': fasta_text(data.records)},
        )
        sent_datasets.add(data.digest)

    @reactive.effect
    @reactive.event(input.run)
    async def _run_biowasm():
        req(ready())
        # Every Run click supersedes any in-flight aligner request: cancel
        # it JS-side too so the superseded wasm computation actually stops
        # instead of burning CPU ahead of the run the user wants.
        if aligner_pending() is not None:
            await _cancel_pending('superseded by a new run')
        if input.input_mode() not in ('fasta', 'genbank'):
            return
        method = input.method()
        if method not in BIOWASM_TOOLS:
            return
        try:
            query, target = _parse_inputs()
            params = _tool_params(method)
            args = build_tool_args(method, params)
        except ValueError as exc:
            ui.notification_show(str(exc), type='error', duration=8)
            return
        cached = cache.get_paf(method, params, query.digest, target.digest)
        if cached is not None:
            logger.info('%s alignment cache hit', method)
            result.set(('paf', cached, {'query': query, 'target': target}))
            return
        if query.total_length + target.total_length > _BIOWASM_SIZE_WARN:
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
        try:
            await _send_dataset(query)
            await _send_dataset(target)
        except ValueError as exc:
            aligner_pending.set(None)
            ui.notification_remove(_NOTIF_ID)
            ui.notification_show(str(exc), type='error', duration=8)
            return
        await session.send_custom_message(
            'rd_run_aligner',
            {
                'tool': method,
                'args': args,
                'query_id': query.digest,
                'target_id': target.digest,
                'request_id': request_id,
            },
        )

    @reactive.effect
    @reactive.event(input.aligner_result)
    def _on_aligner_result():
        res = input.aligner_result()
        info = aligner_pending()
        if not res or res.get('cancelled'):
            return  # cancelled runs were already reported when cancelled
        if not info or res.get('request_id') != info['request_id']:
            return  # stale or unsolicited result
        aligner_pending.set(None)
        # A result can land while the "still running?" modal is open — the
        # run finished on its own.  Take the modal down with it.
        ui.modal_remove()
        ui.notification_remove(_NOTIF_ID)
        method = info['method']
        _log_run(
            method,
            res.get('cmd') or '',
            res.get('stderr') or '',
            res.get('error') or None,
        )
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
            'mounting-data': 'Mounting assemblies into the tool sandbox…',
            'aligning': f'Running {METHOD_LABELS[info["method"]]}…',
            'reading-output': 'Reading alignment output…',
        }
        msg = stage_msgs.get(prog.get('stage'))
        if msg:
            ui.notification_show(msg, id=_NOTIF_ID, duration=None)

    @reactive.effect
    @reactive.event(input.aligner_timeout)
    def _on_aligner_timeout():
        """Offer to keep waiting when a run outlives its watchdog budget.

        The tool has not failed and has not been stopped — aligners.js only
        notifies — so the choice really is "wait longer" vs "give up", and
        it can be offered again each time the extended budget expires.
        """
        ev = input.aligner_timeout()
        info = aligner_pending()
        if not ev or not info or ev.get('request_id') != info['request_id']:
            return  # stale: this run already finished or was superseded
        minutes = max(1, round((ev.get('total_ms') or 0) / 60000))
        label = METHOD_LABELS[info['method']]
        ui.modal_show(
            ui.modal(
                ui.p(
                    f'{label} has been running for about {minutes} minute'
                    f'{"s" if minutes != 1 else ""} and has not finished yet.'
                ),
                ui.p(
                    'It is still working in the background — nothing has been '
                    'lost. Large genomes (especially minimap2 with base-level '
                    'alignment) can take a while in the browser.',
                    class_='text-muted',
                ),
                title='Still aligning',
                easy_close=False,
                footer=ui.TagList(
                    ui.input_action_button(
                        'aligner_give_up', 'Cancel run', class_='btn-outline-secondary'
                    ),
                    ui.input_action_button(
                        'aligner_wait',
                        f'Wait another {_TIMEOUT_EXTEND_MS // 60000} minutes',
                        class_='btn-primary',
                    ),
                ),
            )
        )

    @reactive.effect
    @reactive.event(input.aligner_wait, ignore_init=True)
    async def _on_aligner_wait():
        info = aligner_pending()
        ui.modal_remove()
        if info is None:
            return  # finished while the modal was up
        await session.send_custom_message(
            'rd_extend_aligner',
            {'request_id': info['request_id'], 'ms': _TIMEOUT_EXTEND_MS},
        )
        ui.notification_show(
            f'Still running {METHOD_LABELS[info["method"]]} — waiting another '
            f'{_TIMEOUT_EXTEND_MS // 60000} minutes…',
            id=_NOTIF_ID,
            duration=None,
        )

    @reactive.effect
    @reactive.event(input.aligner_give_up, ignore_init=True)
    async def _on_aligner_give_up():
        await _cancel_pending('timed out')

    # --- end W1: biowasm aligners ---

    # Numeric spinner arrows fire one input event per 0.1 step; without a
    # quiet period every step triggers a full plot re-render.  Only the
    # settled value should reach config().
    @debounce(0.7)
    @reactive.calc
    def dot_size_settled() -> float:
        return input.dot_size() or 0.5

    @debounce(0.7)
    @reactive.calc
    def min_length_settled() -> int:
        return int(input.min_length() or 0)

    # --- GFF annotations -----------------------------------------------------
    # Annotation sources per role (upload-derived state, deliberately outside
    # PlotConfig).  A role can hold both a GenBank-derived set and a GFF
    # upload; they are merged, tagged by the file they came from.  Each entry
    # is {'kind': 'gff'|'genbank', 'filename': str, 'annotation': GffAnnotation}.
    ann_sources = {role: reactive.value(()) for role in ANNOTATION_ROLES}

    def _set_ann_source(role: str, kind: str, filename: str, ann, key=None) -> None:
        """Add, replace or clear one annotation source for *role*.

        *ann* may be a parsed GffAnnotation or GFF3 text.  Every record is
        tagged with *filename* so the drill-down can say where a feature
        came from once several sources are merged.

        *key* identifies the upload — normally ``(content_digest, filename)``.
        When it matches what is already stored this is a no-op: re-running an
        alignment re-parses the same GenBank file every time, and blindly
        re-setting the reactive value would rebuild the feature-type controls
        and wipe the user's toggles, colours and per-feature overrides.
        """
        from rusty_dot.annotation import GffAnnotation  # noqa: PLC0415

        entries = ann_sources[role]()
        # Cheap identity check first, so an unchanged GenBank upload does not
        # even pay for re-parsing its GFF text.
        if key is not None:
            current = next((e for e in entries if e['kind'] == kind), None)
            if current is not None and current.get('key') == key:
                return
        if isinstance(ann, str):
            ann = GffAnnotation.from_text(ann) if ann.strip() else None
        updated = replace_source(entries, kind, filename, ann, key)
        if updated is None:
            return  # nothing changed; leave the controls and overrides alone
        if ann is not None and len(ann):
            for rec in ann.records:
                rec.source_file = filename
        ann_sources[role].set(updated)
        _reset_feature_overrides()

    def _parse_gff_upload(file_input, role: str) -> None:
        from rusty_dot.annotation import GffAnnotation  # noqa: PLC0415

        files = file_input()
        if not files:
            _set_ann_source(role, 'gff', '', None)
            return
        raw = Path(files[0]['datapath']).read_bytes()
        try:
            ann = GffAnnotation.from_bytes(raw)
        except ValueError as exc:
            _set_ann_source(role, 'gff', '', None)
            ui.notification_show(
                f'Could not parse {role} GFF: {exc}', type='error', duration=10
            )
            return
        if len(ann) == 0:
            _set_ann_source(role, 'gff', '', None)
            ui.notification_show(
                f'No features found in the {role} GFF file.',
                type='warning',
                duration=8,
            )
            return
        _set_ann_source(
            role,
            'gff',
            files[0]['name'],
            ann,
            key=(content_digest(raw), files[0]['name']),
        )
        kinds = {e['kind'] for e in ann_sources[role]()}
        if 'genbank' in kinds:
            # Overlapping duplicates are likely, but silently deduping would
            # discard features the user deliberately supplied.  Say so and
            # let them toggle sources off instead.
            ui.notification_show(
                f'{role.capitalize()} now has annotations from both a GenBank '
                'file and a GFF; features present in both will be drawn twice.',
                type='warning',
                duration=10,
            )
        ui.notification_show(
            f'{role.capitalize()} annotations: {len(ann)} feature(s), '
            f'{len(ann.feature_types())} type(s).',
            type='message',
            duration=5,
        )

    def gff_raw_for(role: str):
        """Return the merged annotation for *role* across all its sources."""
        return merge_annotations([e['annotation'] for e in ann_sources[role]()])

    # Per-feature state from the drill-down Annotations tab, keyed by the
    # positional uids build_feature_rows hands out.  Cleared whenever a
    # role's sources change, since re-uploading renumbers those uids.
    feature_hidden = reactive.value(frozenset())
    feature_colors = reactive.value({})

    def _reset_feature_overrides() -> None:
        feature_hidden.set(frozenset())
        feature_colors.set({})

    @reactive.effect
    @reactive.event(input.query_gff)
    def _on_query_gff():
        req(ready())
        _parse_gff_upload(input.query_gff, 'query')

    @reactive.effect
    @reactive.event(input.target_gff)
    def _on_target_gff():
        req(ready())
        _parse_gff_upload(input.target_gff, 'target')

    _MIN_LEN_NOTIF_ID = 'rd_min_contig_len'

    @reactive.effect
    def _report_excluded_contigs():
        """Say how many contigs the length filter is hiding.

        Silently dropping panels looks like missing data, so keep a
        standing note while the filter is active.  Uses a fixed id so
        adjusting the threshold replaces the note instead of stacking up.
        """
        if result() is None:
            return
        lay = layout()
        dropped = len(lay['excluded_query']) + len(lay['excluded_target'])
        min_len = max(0, int(input.min_contig_len() or 0))
        if not min_len:
            ui.notification_remove(_MIN_LEN_NOTIF_ID)
            return
        if not dropped:
            ui.notification_show(
                f'No contigs are shorter than {min_len:,} bp; nothing excluded.',
                id=_MIN_LEN_NOTIF_ID,
                type='message',
                duration=6,
            )
            return
        ui.notification_show(
            f'Hiding {dropped} contig(s) shorter than {min_len:,} bp. '
            'They are still included in the reordered-FASTA download.',
            id=_MIN_LEN_NOTIF_ID,
            type='message',
            duration=8,
        )

    @reactive.effect
    def _warn_on_annotation_name_mismatch():
        """Warn when a GFF annotates contigs the assembly does not have.

        Keyed on both the result and the annotation sources, so it fires
        whichever order the user supplies them in — a GFF can be uploaded
        long before the first run, and a new assembly can be run against an
        already-loaded GFF.  The plotter only logs about this, which the
        browser user never sees.
        """
        res = result()
        if res is None:
            return
        _kind, _obj, meta = res
        for role in ANNOTATION_ROLES:
            ann = gff_raw_for(role)
            fasta = meta.get(role)
            if ann is None or not isinstance(fasta, FastaInput):
                continue
            for warning in validate_annotation_names(
                fasta.names, ann.sequence_names(), role
            ):
                ui.notification_show(warning, type='warning', duration=12)

    @reactive.effect
    @reactive.event(input.clear_gff, ignore_init=True)
    async def _on_clear_gff():
        """Drop every uploaded GFF and reset the file inputs.

        Shiny's file-input binding has a no-op ``setValue``, so the widget
        cannot be cleared from the server: without the custom message the
        filename would linger, and re-selecting the same file would fire no
        change event and so never reload.
        """
        for role in ANNOTATION_ROLES:
            _set_ann_source(role, 'gff', '', None)
        await session.send_custom_message(
            'rd_clear_file_inputs', {'ids': ['query_gff', 'target_gff']}
        )
        ui.notification_show(
            'Cleared uploaded annotations.', type='message', duration=4
        )

    def _read_dynamic(input_id: str, default):
        """Read a dynamically rendered input, tolerating its absence."""
        try:
            return input[input_id]()
        except Exception:  # noqa: BLE001 - silent until the control renders
            # Expected before gff_controls() has rendered; a genuine id typo
            # looks identical, so leave a breadcrumb rather than a silent
            # fallback to the default for every type.
            logger.debug('dynamic input %r not available; using default', input_id)
            return default

    @reactive.calc
    def gff_type_index():
        """Union of feature types across roles, with shared colours.

        Returns ``(rows, slugs, shared)`` where *rows* is one entry per
        normalised type -- ``{'key', 'label', 'roles', 'color'}`` -- and
        *slugs* maps normalised key to input-id fragment.  Colours are
        assigned over the union so the same type never gets two colours on
        the two axes (see core/annotation_colors.py).
        """
        types_by_role = {
            role: list(ann.feature_types())
            for role in ANNOTATION_ROLES
            if (ann := gff_raw_for(role)) is not None
        }
        if not types_by_role:
            return [], {}, {}
        shared = assign_shared_colors(types_by_role)
        spellings: dict[str, list[str]] = {}
        roles_with: dict[str, list[str]] = {}
        for role, fts in types_by_role.items():
            for ft in fts:
                key = normalise_type(ft)
                spellings.setdefault(key, []).append(ft)
                if role not in roles_with.setdefault(key, []):
                    roles_with[key].append(role)
        slugs = type_slug_map(sorted(spellings))
        rows = [
            {
                'key': key,
                'label': display_name(spellings[key]),
                'roles': roles_with[key],
                'color': shared[key],
            }
            for key in sorted(spellings)
        ]
        return rows, slugs, shared

    @render.ui
    def gff_controls():
        rows, slugs, _shared = gff_type_index()
        if not rows:
            return None
        controls = []
        for row in rows:
            badges = ''.join(
                f'<span class="rd-gff-role" title="{r} annotations">'
                f'{r[0].upper()}</span>'
                for r in row['roles']
            )
            controls.append(
                ui.div(
                    ui.input_checkbox(f'gtyp_{slugs[row["key"]]}', row['label'], True),
                    ui.HTML(badges),
                    ui.HTML(
                        f'<input type="color" class="rd-color-input" '
                        f'id="gcol_{slugs[row["key"]]}" value="{row["color"]}">'
                    ),
                    class_='rd-gff-type-row',
                )
            )
        # Deliberately depends on gff_type_index() alone.  The gff_diagonal /
        # gff_tracks toggles live in the static sidebar because anything that
        # made this output re-render — self_panels() used to, via result() —
        # rebuilds every checkbox and colour picker below at its default,
        # silently discarding the user's choices on each run.
        return ui.div(
            ui.div(
                ui.h6('Feature types'),
                *controls,
                class_='rd-gff-section',
            ),
        )

    @reactive.calc
    def annotations():
        """Return (query_ann, target_ann) with the user's type/colour choices."""
        rows, slugs, shared = gff_type_index()
        # One control per *normalised* type, shared by both roles.
        chosen = {
            row['key']: (
                bool(_read_dynamic(f'gtyp_{slugs[row["key"]]}', True)),
                str(_read_dynamic(f'gcol_{slugs[row["key"]]}', '') or ''),
            )
            for row in rows
        }
        hidden = feature_hidden()
        feat_colors = feature_colors()
        result = {}
        for role in ANNOTATION_ROLES:
            ann = gff_raw_for(role)
            if ann is None:
                result[role] = None
                continue
            # Per-feature choices first: their uids are positions in the
            # *unfiltered* record list, so filtering by type first would
            # renumber them.  The type filter runs second and therefore
            # wins — a type switched off hides its features regardless.
            ann = apply_feature_overrides(ann, hidden, feat_colors, role)
            if ann is None:
                result[role] = None
                continue
            fts = list(ann.feature_types())
            enabled = {ft: chosen.get(normalise_type(ft), (True, ''))[0] for ft in fts}
            picked = {key: color for key, (_on, color) in chosen.items() if color}
            colors = color_map_for(fts, {**shared, **picked})
            result[role] = apply_annotation_config(ann, enabled, colors)
        return result['query'], result['target']

    # --- end GFF annotations -------------------------------------------------

    @reactive.calc
    def config() -> PlotConfig:
        return PlotConfig(
            contig_order=input.contig_order(),
            auto_reverse=input.auto_reverse(),
            hide_internal_axes=input.hide_internal_axes(),
            dot_size=dot_size_settled(),
            min_length=min_length_settled(),
            color_by_identity=bool(input.color_by_identity()),
            identity_palette=input.identity_palette() or 'viridis',
        )

    @reactive.calc
    def render_config() -> PlotConfig:
        """Structural plot options — the ones that genuinely need a re-render.

        Display-only options (line width, min match length) are applied
        client-side inside the embedded report via ``rd_display_opts``
        messages, so they are deliberately absent here: the report HTML is
        rendered with ``min_length=0`` and the default line width, and the
        client owns both from then on.  The static plot and the SVG/PDF
        downloads keep using the full :func:`config` (server-side
        semantics unchanged there).
        """
        return PlotConfig(
            contig_order=input.contig_order(),
            auto_reverse=input.auto_reverse(),
            hide_internal_axes=input.hide_internal_axes(),
            color_by_identity=bool(input.color_by_identity()),
            identity_palette=input.identity_palette() or 'viridis',
        )

    @reactive.effect
    async def _send_display_opts():
        # Forward the debounced display options to the embedded report
        # (bridge.js relays them into the iframe), where they are applied
        # instantly — no matplotlib re-render in Pyodide.
        await session.send_custom_message(
            'rd_display_opts',
            {
                'dot_size': float(dot_size_settled()),
                'min_length': int(min_length_settled()),
            },
        )

    # --- W2: interactive plot ------------------------------------------------
    # Focused (query, target) contig pair for the drill-down view, or None
    # for the full overview grid.
    focus = reactive.value(None)
    # Per-result memos, keyed by the alignment object's identity: computed
    # contig orders per ordering mode (gravity is seconds-long on real
    # assemblies) and the DotPlotter/PafAlignment pair (whose construction
    # copies the full record list).  Both are invalidated with the result.
    order_cache: dict = {}
    figure_ctx_cache: dict = {}
    # (query_name, target_name) -> [PafRecord], filled lazily from the
    # current result's records for aligned-sequence lookups.  The grouping
    # mirrors DotPlotter._records_for_pair, so a payload segment index maps
    # straight onto the per-pair list.
    paf_pair_index: dict = {}

    @reactive.effect
    def _reset_focus_on_new_result():
        result()
        focus.set(None)
        order_cache.clear()
        figure_ctx_cache.clear()
        paf_pair_index.clear()

    @reactive.calc
    def ordering_config() -> tuple[str, bool]:
        """Return only the config fields that affect contig ordering.

        ``layout()`` depends on this instead of ``config()`` so display-only
        edits (line width, min length, palettes, styling) never re-run the
        gravity ordering — on real assemblies that reorder is the expensive
        step.  Reads the raw inputs directly: going through ``config()``
        would re-invalidate ``layout()`` on every config recompute.
        """
        return (
            input.contig_order(),
            input.auto_reverse(),
            max(0, int(input.min_contig_len() or 0)),
        )

    def _axis_inputs(res) -> tuple[list[str], list[str], dict[str, int]]:
        """Return ``(query_names, target_names, lengths)`` before ordering.

        Shared by :func:`layout` and :func:`grid_panel_count` so the two
        cannot disagree about which contigs the grid contains.
        """
        kind, obj, meta = res
        if kind == 'kmer':
            q_in = list(meta['query'].names)
            t_in = list(meta['target'].names)
            lengths = {n: len(s) for n, s in meta['target'].records}
            lengths.update({n: len(s) for n, s in meta['query'].records})
        else:
            q_in = list(obj.query_names)
            t_in = list(obj.target_names)
            lengths = {n: obj.get_sequence_length(n) for n in (*q_in, *t_in)}
        return q_in, t_in, lengths

    @reactive.calc
    def grid_panel_count() -> int:
        """Count the panels the grid will draw, from the result alone.

        Deliberately does NOT go through ``layout()``: that would pull
        ``input.contig_order`` / ``input.auto_reverse`` into ``plot_area``,
        and re-rendering ``plot_area`` rebuilds the report ``<iframe>``,
        throwing away its client-side state (zoom, highlight bands, pushed
        display options, the selected match).  Reading ``result()`` — which
        ``plot_area`` already depends on — costs no extra invalidation.

        This is an upper bound: the ``colinearity`` ordering modes delegate
        to ``compute_gravity_contigs``, which may return fewer contigs.  The
        worst case is one over-advertised navigation tip.
        """
        res = result()
        if res is None:
            return 0
        q_in, t_in, lengths = _axis_inputs(res)
        min_len = max(0, int(input.min_contig_len() or 0))
        q_keep, _ = filter_by_min_length(q_in, lengths, min_len)
        t_keep, _ = filter_by_min_length(t_in, lengths, min_len)
        # Same empty-axis fallback as layout(); the two must agree or the
        # hint advertises a click-to-focus the report has disabled.
        return (len(q_keep) or len(q_in)) * (len(t_keep) or len(t_in))

    @reactive.calc
    def self_panels() -> bool:
        """Whether the panel grid contains a self-comparison panel.

        Gates the diagonal-shading control, which can only ever draw on
        such a panel.  Judged on the whole grid rather than the focused
        pair, so the control does not appear and vanish as the user drills
        in and out.  Returns ``False`` rather than blocking before the
        first run, so the feature-type list stays visible while a plot is
        being set up.
        """
        if input.input_mode() != 'paf' and input.self_align():
            return True
        if result() is None:
            return False
        lay = layout()
        return has_self_pair(lay['query_names'], lay['target_names'])

    @reactive.calc
    def layout():
        """Explicit plotted axis orders for the current result + config.

        Single source of truth for row/column order: the plot call, the
        panel double-click mapping and the FASTA download all use it.
        """
        res = result()
        req(res)
        contig_order, auto_reverse, min_len = ordering_config()
        kind, obj, _meta = res
        # Memoise per ordering mode: re-selecting e.g. 'maximise colinearity'
        # after trying another mode must not recompute the gravity sort.
        # The length filter is part of the key because it changes which
        # contigs the ordering runs over.  auto_reverse stays outside it —
        # it only selects whether the cached reversed set is applied.
        cache_key = (id(obj), contig_order, min_len)
        cached = order_cache.get(cache_key)
        if cached is None:
            q_in, t_in, lengths = _axis_inputs(res)
            records = (
                obj.get_records_for_pair(QUERY_GROUP, TARGET_GROUP)
                if kind == 'kmer'
                else obj.records
            )
            q_keep, q_drop = filter_by_min_length(q_in, lengths, min_len)
            t_keep, t_drop = filter_by_min_length(t_in, lengths, min_len)
            # A threshold above every contig would leave nothing to draw;
            # keep that axis whole rather than rendering an empty grid.
            if not q_keep:
                q_keep, q_drop = q_in, []
            if not t_keep:
                t_keep, t_drop = t_in, []
            cached = (
                *resolve_orders(contig_order, records, q_keep, t_keep, lengths),
                q_drop,
                t_drop,
            )
            order_cache[cache_key] = cached
        else:
            logger.info('contig-order cache hit for mode %r', contig_order)
        q_order, t_order, reversed_q, q_drop, t_drop = cached
        reverse = reversed_q if auto_reverse else set()
        # Return copies so downstream mutation cannot poison the memo.
        return {
            'query_names': list(q_order),
            'target_names': list(t_order),
            'reverse': set(reverse),
            # Excluded by the length filter.  Kept so the FASTA export can
            # stay complete: CrossIndex.write_fasta writes exactly the names
            # it is given, so anything omitted here is silently lost.
            'excluded_query': list(q_drop),
            'excluded_target': list(t_drop),
        }

    def make_figure(res, cfg: PlotConfig, lay, pair=None, output_path=None):
        from rusty_dot import DotPlotter
        from rusty_dot.paf_io import PafAlignment

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
        # Identity colouring is unreadable without a key.
        kwargs['identity_colorbar'] = bool(kwargs.get('color_by_identity'))
        kwargs['reverse_contigs'] = set(lay['reverse'])
        if pair is not None and cfg.title is None:
            kwargs['title'] = f'{pair[0]} vs {pair[1]}'
        # GFF annotations: diagonal shading on self panels plus side tracks
        # in the focused (1x1) drill-down view.  Reading annotations() here
        # keeps every figure consumer reactive to toggle/colour changes
        # without touching layout()'s dependencies.
        ann_q, ann_t = annotations()
        if ann_q is not None or ann_t is not None:
            # Static inputs now, so read them directly.  The self_panels()
            # conjunct still matters: a hidden checkbox keeps its last value,
            # so a box ticked during a self-comparison must not shade the
            # diagonal of a later cross-assembly run.
            if self_panels() and bool(input.gff_diagonal()):
                # On a self-comparison both roles describe the same
                # sequences, so merge them: a user who uploaded the
                # annotation under either role — or split it across two
                # files — gets everything shaded.
                kwargs['annotation'] = merge_annotations([ann_q, ann_t])
            kwargs['annotation_query'] = ann_q
            kwargs['annotation_target'] = ann_t
            kwargs['annotation_tracks'] = bool(input.gff_tracks())
        # Memoise the plotter per result: for the k-mer path, building the
        # PafAlignment copies the full record list — pointless to repeat on
        # every re-render of the same result.
        plotter = figure_ctx_cache.get(id(obj))
        if plotter is None:
            if kind == 'kmer':
                # Internal 'group:name' identifiers keep name resolution
                # unambiguous; the cached cross-group records are handed to
                # the plotter directly so nothing is recomputed.
                paf = PafAlignment(obj.get_records_for_pair(QUERY_GROUP, TARGET_GROUP))
                plotter = DotPlotter(obj, paf_alignment=paf)
            else:
                plotter = DotPlotter(obj)
            figure_ctx_cache[id(obj)] = plotter
        if kind == 'kmer':
            kwargs['query_names'] = [
                obj.make_internal_name(QUERY_GROUP, n) for n in q_names
            ]
            kwargs['target_names'] = [
                obj.make_internal_name(TARGET_GROUP, n) for n in t_names
            ]
        else:
            kwargs['query_names'] = list(q_names)
            kwargs['target_names'] = list(t_names)
        if output_path is not None:
            return plotter.to_html(output_path, **kwargs)
        return plotter.plot(**kwargs)

    def _report_html(pair=None) -> str:
        """Render the interactive HTML report and return it as a string.

        Uses :func:`render_config` (structural options only): line width
        and min match length are applied client-side in the embedded
        report, so changing them never re-runs this seconds-long render.
        """
        import matplotlib.pyplot as plt

        res = result()
        req(res)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'report.html'
            fig = make_figure(
                res, render_config(), layout(), pair=pair, output_path=path
            )
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
    def aligner_log_ui():
        entries = aligner_log()
        if not entries:
            return None
        blocks = []
        for e in entries:
            text = f'$ {e["cmd"]}\n{e["stderr"] or "(no tool output)"}'
            if e.get('error'):
                text += f'\nERROR: {e["error"]}'
            blocks.append(ui.tags.pre(text, class_='rd-log'))
        return ui.accordion(
            ui.accordion_panel(
                f'Aligner log ({len(entries)} run(s))', *blocks, value='log'
            ),
            open=False,
            class_='rd-log-accordion',
        )

    @output(suspend_when_hidden=False)
    @render.ui
    def report_frame():
        # suspend_when_hidden=False: switching to the drill-down's
        # Annotations tab hides this output, and Shiny would re-render it
        # on the way back — a seconds-long matplotlib pass under Pyodide
        # for a view that has not changed.
        #
        # Rendering the report re-runs layout (contig ordering) and the
        # matplotlib figure — seconds-long on real assemblies, so show what
        # is happening rather than freezing silently.
        with ui.Progress(min=0, max=2) as progress:
            progress.set(1, message='Computing layout & rendering report…')
            html = focus_html() if focus() is not None else overview_html()
            progress.set(2, message='Done')
        return ui.tags.iframe(
            srcdoc=html,
            class_='rd-report-frame',
            sandbox='allow-scripts',
            title='Interactive dotplot report',
        )

    def _nav_hint(focused: bool, multi_panel: bool):
        """Build the navigation-tips box shown under the interactive report.

        Parameters
        ----------
        focused : bool
            Whether the focused single-pair view is active; panel-click tips
            do not apply there (click-to-focus is disabled on single-panel
            reports).
        multi_panel : bool
            Whether the grid has more than one panel.  See
            :func:`core.panels.nav_tips`.

        Returns
        -------
        htmltools.Tag
            The hint ``div`` with each action term in bold.
        """
        tips = nav_tips(focused, multi_panel)
        parts: list = [ui.tags.b('Navigate: ')]
        for i, (action, effect) in enumerate(tips):
            if i:
                parts.append(' · ')
            parts += [ui.tags.b(action), f' = {effect}']
        return ui.div(*parts, class_='rd-nav-hint')

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
            body = ui.output_ui('report_frame')
        else:
            body = ui.output_plot('dotplot', height='72vh')
        hint = (
            _nav_hint(pair is not None, grid_panel_count() > 1)
            if input.interactive()
            else None
        )
        if pair is None or not feature_rows():
            return ui.div(
                ui.div(*toolbar, class_='rd-plot-toolbar') if toolbar else None,
                body,
                hint,
                class_='rd-plot-area',
            )
        # Tabs only in the drill-down, so the overview path is untouched.
        return ui.div(
            ui.div(*toolbar, class_='rd-plot-toolbar') if toolbar else None,
            ui.navset_tab(
                ui.nav_panel('Plot', body, hint),
                ui.nav_panel('Annotations', ui.output_ui('annotation_table')),
                id='drill_tabs',
            ),
            class_='rd-plot-area',
        )

    @reactive.calc
    def feature_rows() -> list[dict]:
        """Every feature on the focused pair's sequences, as table rows."""
        pair = focus()
        if pair is None:
            return []
        # Both roles are always listed: on a self panel the two axes carry
        # the same sequence but may still hold different uploads.
        return build_feature_rows(
            gff_raw_for('query'), pair[0], 'query'
        ) + build_feature_rows(gff_raw_for('target'), pair[1], 'target')

    @reactive.effect
    @reactive.event(input.feature_table_change)
    def _on_feature_table_change():
        """Fold one delta from the annotations table into the override state.

        The table ships raw HTML inputs and posts deltas through
        www/feature-table.js rather than binding one Shiny input per
        feature: the real GFFs here run to ~600 features per sequence, and
        ~1200 bound inputs visibly freezes Pyodide.
        """
        ev = input.feature_table_change()
        if not ev:
            return
        kind = ev.get('kind')
        uids = ev.get('uids') or ([ev['uid']] if ev.get('uid') else [])
        if not uids:
            return
        if kind in ('vis', 'bulk'):
            visible = bool(ev.get('value'))
            hidden = set(feature_hidden())
            hidden.difference_update(uids) if visible else hidden.update(uids)
            feature_hidden.set(frozenset(hidden))
        elif kind == 'color':
            colors = dict(feature_colors())
            value = ev.get('value') or ''
            for uid in uids:
                if value:
                    colors[uid] = value
                else:
                    colors.pop(uid, None)  # reset to the type colour
            feature_colors.set(colors)

    @render.ui
    def annotation_table():
        rows = feature_rows()
        pair = focus()
        if not rows or pair is None:
            return None
        hidden = feature_hidden()
        overrides = feature_colors()
        _type_rows, _slugs, shared = gff_type_index()

        def type_color(ft: str) -> str:
            return shared.get(normalise_type(ft), '#888888')

        body = []
        for r in rows:
            checked = '' if r['uid'] in hidden else ' checked'
            color = overrides.get(r['uid']) or type_color(r['type'])
            meta = ' · '.join(f'{k}={v}' for k, v in list(r['attributes'].items())[:6])
            body.append(
                '<tr data-uid="{uid}" data-type="{type}" data-source="{src}">'
                '<td><input type="checkbox" data-uid="{uid}" data-kind="vis"{ck}></td>'
                '<td><input type="color" data-uid="{uid}" data-kind="color" '
                'value="{color}" class="rd-color-input"></td>'
                '<td class="rd-ft-role">{role}</td>'
                '<td class="rd-ft-type">{type}</td>'
                '<td class="rd-ft-name">{name}</td>'
                '<td class="rd-ft-num">{start:,}</td>'
                '<td class="rd-ft-num">{end:,}</td>'
                '<td class="rd-ft-num">{length:,}</td>'
                '<td>{strand}</td>'
                '<td class="rd-ft-src">{src}</td>'
                '<td class="rd-ft-attrs" title="{meta}">{meta}</td>'
                '</tr>'.format(
                    uid=_esc(r['uid']),
                    ck=checked,
                    color=_esc(color),
                    role=_esc(r['role']),
                    type=_esc(r['type']),
                    name=_esc(r['name'] or ''),
                    start=r['start'],
                    end=r['end'],
                    length=r['length'],
                    strand=_esc(r['strand']),
                    src=_esc(r['source_file'] or r['source'] or ''),
                    meta=_esc(meta),
                )
            )
        header = (
            '<thead><tr><th>Show</th><th>Colour</th><th>Axis</th><th>Type</th>'
            '<th>Name</th><th>Start</th><th>End</th><th>Length</th>'
            '<th>Strand</th><th>Source</th><th>Attributes</th></tr></thead>'
        )
        return ui.div(
            ui.div(
                ui.HTML(
                    f'<b>{len(rows)}</b> feature(s) on '
                    f'<b>{_esc(pair[0])}</b> and <b>{_esc(pair[1])}</b>. '
                    'Coordinates are 1-based inclusive, as in the source file.'
                ),
                class_='rd-ft-caption',
            ),
            ui.div(
                ui.HTML(
                    '<input type="search" id="rd-ft-filter" '
                    'placeholder="Filter by type, name or source…">'
                    '<button type="button" data-bulk="show">Show all</button>'
                    '<button type="button" data-bulk="hide">Hide all</button>'
                    '<button type="button" data-bulk="reset">Reset colours</button>'
                ),
                class_='rd-ft-tools',
            ),
            ui.HTML(
                f'<div class="rd-ft-scroll"><table id="rd-feature-table">'
                f'{header}<tbody>{"".join(body)}</tbody></table></div>'
            ),
            class_='rd-ft-panel',
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

    def _pair_records(obj, q: str, t: str) -> list:
        """Return the PAF records for one pair, building the index lazily."""
        if not paf_pair_index:
            for rec in obj.records:
                paf_pair_index.setdefault((rec.query_name, rec.target_name), []).append(
                    rec
                )
        return paf_pair_index.get((q, t), [])

    def _record_for_segment(obj, q: str, t: str, idx: int):
        """Return the PAF record behind identity-layer payload row *idx*."""
        recs = _pair_records(obj, q, t)
        return recs[idx] if 0 <= idx < len(recs) else None

    def _genomic_coords(info: dict, q: str, lay, qlen: int):
        """Map a clicked segment's payload coords to genomic values.

        The payload's query side is mirrored on reverse-oriented contigs;
        the target side is always genomic.

        Returns
        -------
        tuple or None
            ``(q_start, q_end, t_start, t_end, strand)`` or ``None`` when
            the message lacks usable coordinates.
        """
        try:
            gqs, gqe = int(info['qs']), int(info['qe'])
            ts_, te_ = int(info['ts']), int(info['te'])
        except (KeyError, TypeError, ValueError):
            return None
        strand = info.get('strand', '+')
        if strand not in ('+', '-'):
            strand = '+'
        if q in lay['reverse']:
            gqs, gqe = qlen - gqe, qlen - gqs
            strand = '-' if strand == '+' else '+'
        return gqs, gqe, ts_, te_, strand

    def _preview_slice(seq: str, start: int, end: int, minus: bool) -> str:
        """Return an alignment-oriented preview of ``seq[start:end]``.

        Clipped to 20,000 bases *before* any copy or revcomp so megabase
        matches never materialise a full slice for the preview.  On the
        minus strand the alignment-oriented sequence begins at the genomic
        end, so the window is taken from there.
        """
        from rusty_dot.alignment_view import revcomp

        cap = 20_000
        n = end - start
        if n <= cap:
            s = seq[start:end]
            return revcomp(s) if minus else s
        if minus:
            s = revcomp(seq[end - cap : end])
        else:
            s = seq[start : start + cap]
        return s + f'… [truncated at {cap:,} bases]'

    def _sequence_for(meta: dict, name: str, side: str):
        """Look up a full sequence by name, preferring *side*'s FASTA.

        Falls back to the other assembly so self-align mode (one file) and
        PAF uploads with a single FASTA still resolve.
        """
        order = ('query', 'target') if side == 'query' else ('target', 'query')
        for key in order:
            fasta = meta.get(key)
            if fasta is not None:
                seq = next((s for n, s in fasta.records if n == name), None)
                if seq is not None:
                    return seq
        return None

    def _match_context(info):
        """Resolve a clicked segment to sequences and genomic coordinates.

        Shared by the match-detail and copy-request handlers.

        Returns
        -------
        dict or None
            ``{'q', 't', 'rec', 'qseq', 'tseq', 'gqs', 'gqe', 'ts', 'te',
            'strand'}`` with genomic (unmirrored) coordinates, or ``None``
            when the panel or sequences cannot be resolved.  ``rec`` is
            the exact PAF record when one matches (its CIGAR, if any,
            enables the gapped view).
        """
        res = result()
        if not info or res is None:
            return None
        kind, obj, meta = res
        pair = focus()
        if pair is not None:
            q, t = pair
            lay = layout()
        else:
            try:
                lay = layout()
                q, t = panel_pair(
                    lay['query_names'],
                    lay['target_names'],
                    int(info['row']),
                    int(info['col']),
                )
            except (KeyError, TypeError, ValueError, IndexError):
                return None
        qseq = _sequence_for(meta, q, 'query')
        tseq = _sequence_for(meta, t, 'target')
        if qseq is None or tseq is None:
            return None
        rec = None
        if kind == 'paf' and info.get('layer') == 'identity':
            try:
                rec = _record_for_segment(obj, q, t, int(info['idx']))
            except (TypeError, ValueError):
                rec = None
        # Strand-coloured (fwd/rev) layers have no stable index->record
        # mapping (blocks may be chained), but an unchained block's coords
        # match its record exactly — recover the CIGAR that way so runs
        # with -c get the gapped view without identity colouring on.
        if rec is None and kind == 'paf':
            g = _genomic_coords(info, q, lay, len(qseq))
            if g is not None:
                gqs0, gqe0, ts0, te0, gstrand0 = g
                rec = next(
                    (
                        r
                        for r in _pair_records(obj, q, t)
                        if r.query_start == gqs0
                        and r.query_end == gqe0
                        and r.target_start == ts0
                        and r.target_end == te0
                        and r.strand == gstrand0
                    ),
                    None,
                )
        if rec is not None:
            gqs, gqe = rec.query_start, rec.query_end
            ts_, te_ = rec.target_start, rec.target_end
            strand = rec.strand
        else:
            g = _genomic_coords(info, q, lay, len(qseq))
            if g is None:
                return None
            gqs, gqe, ts_, te_, strand = g
        return {
            'q': q,
            't': t,
            'rec': rec,
            'qseq': qseq,
            'tseq': tseq,
            'gqs': gqs,
            'gqe': gqe,
            'ts': ts_,
            'te': te_,
            'strand': strand,
        }

    @reactive.effect
    @reactive.event(input.match_select)
    async def _on_match_select():
        """Serve the sequence preview for a clicked match.

        The report posts 'rd-match-select'; the reply travels back through
        bridge.js as an 'rd-seq-response' echoing row/col/layer/idx so the
        report can discard stale responses.  When the record has a CIGAR
        (minimap2 ``-c``) the reply is a gapped alignment view; otherwise
        the raw query and target slices are shown unaligned.  Only the
        clipped preview is sent — full sequences are fetched separately on
        copy (`_on_copy_request`), so click latency stays flat regardless
        of match size.
        """
        from rusty_dot.alignment_view import aligned_text

        info = input.match_select()
        ctx = _match_context(info)
        reply = {k: (info or {}).get(k) for k in ('row', 'col', 'layer', 'idx')}
        if ctx is None:
            reply['error'] = 'no-sequences'
            await session.send_custom_message('rd_match_seq', reply)
            return
        rec = ctx['rec']
        # Full sequences are available for on-demand copy in either branch.
        reply['copy'] = True
        if rec is not None and rec.cigar is not None:
            view = aligned_text(rec, ctx['qseq'], ctx['tseq'])
            reply['aligned'] = True
            reply['text'] = (
                f'{rec.query_name}:{rec.query_start:,}-{rec.query_end:,} '
                f'({rec.strand}) vs '
                f'{rec.target_name}:{rec.target_start:,}-{rec.target_end:,}'
                f'\n\n{view["text"]}'
            )
        else:
            reply['aligned'] = False
            reply['text'] = (
                f'{ctx["q"]}:{ctx["gqs"]:,}-{ctx["gqe"]:,} ({ctx["strand"]}) '
                f'vs {ctx["t"]}:{ctx["ts"]:,}-{ctx["te"]:,}\n'
                'No CIGAR for this match — sequences shown unaligned; run '
                'minimap2 with base-level alignment (-c) for a gapped '
                'alignment view.\n\n'
                f'>query {ctx["q"]}:{ctx["gqs"]:,}-{ctx["gqe"]:,} '
                f'({ctx["strand"]})\n'
                f'{_preview_slice(ctx["qseq"], ctx["gqs"], ctx["gqe"], ctx["strand"] == "-")}\n'
                f'>target {ctx["t"]}:{ctx["ts"]:,}-{ctx["te"]:,}\n'
                f'{_preview_slice(ctx["tseq"], ctx["ts"], ctx["te"], False)}'
            )
        await session.send_custom_message('rd_match_seq', reply)

    @reactive.effect
    @reactive.event(input.copy_request)
    async def _on_copy_request():
        """Serve one full match sequence for a copy-button press.

        Shipping megabases of sequence with every click made match details
        slow; instead the detail reply carries only the preview and the
        full sequence crosses the boundary once, here, when actually
        requested.  The reply echoes the segment key plus ``side`` so the
        report caches it client-side (repeat copies are then instant).
        """
        from rusty_dot.alignment_view import revcomp

        info = input.copy_request()
        side = (info or {}).get('side')
        reply = {k: (info or {}).get(k) for k in ('row', 'col', 'layer', 'idx')}
        reply['side'] = side
        ctx = _match_context(info)
        if ctx is None or side not in ('query', 'target'):
            reply['error'] = 'no-sequences'
        elif side == 'query':
            seq = ctx['qseq'][ctx['gqs'] : ctx['gqe']]
            reply['seq'] = revcomp(seq) if ctx['strand'] == '-' else seq
        else:
            reply['seq'] = ctx['tseq'][ctx['ts'] : ctx['te']]
        await session.send_custom_message('rd_copy_seq', reply)

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
    def app_memory():
        # suspend_when_hidden=False: the fixed note is CSS-hidden while its
        # text is empty, and Shiny would otherwise never render into a
        # hidden output — leaving it permanently empty (and hidden).
        # Refresh every 5 s.  The wasm linear memory only ever grows, so
        # this reports the high-water mark of the Python runtime's heap —
        # the resource the ~2 GB browser cap actually constrains.  Native
        # runs (and non-Pyodide environments) show nothing.
        reactive.invalidate_later(5)
        if sys.platform != 'emscripten':
            return ''
        try:
            import pyodide_js  # noqa: PLC0415 - pyodide-only module

            used = int(pyodide_js._module.HEAPU8.length)
        except Exception:  # noqa: BLE001 - internals unavailable -> hide
            return ''
        return f'App memory: {used / 1048576:.0f} MB of ~2048 MB wasm heap'

    @output(suspend_when_hidden=False)
    @render.text
    def result_kind():
        # Hidden output driving panel_conditional visibility (e.g. the
        # identity-colouring controls, shown only for tool/PAF results).
        # suspend_when_hidden=False: the output lives in a display:none div,
        # which Shiny would otherwise never update.
        res = result()
        return res[0] if res else ''

    @output(suspend_when_hidden=False)
    @render.text
    def gff_mode():
        """Visibility state for the static annotation toggles.

        ``''`` no annotations · ``'plain'`` cross-comparison ·
        ``'self'`` at least one self-comparison panel.  Conditions on this
        output test *positively* (``=== 'self'``), because before the first
        report it is ``undefined`` client-side and a negative test would
        flash the controls on a fresh page.

        Checks ``ann_sources`` rather than ``gff_type_index()`` — same truth
        value, without re-running the shared-colour assignment on each run.
        """
        if not any(gff_raw_for(role) is not None for role in ANNOTATION_ROLES):
            return ''
        return 'self' if self_panels() else 'plain'

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
        with ui.Progress(min=0, max=2) as progress:
            progress.set(1, message='Computing layout & rendering plot…')
            fig = make_figure(res, config(), layout(), pair=focus())
            progress.set(2, message='Done')
        return fig

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
        # Contigs the length filter left out of the plot still belong in the
        # export — they are simply unordered, appended after the plotted
        # ones.  CrossIndex.write_fasta writes exactly the names it is
        # given, so omitting them here would quietly drop sequence.
        order = lay['query_names'] + lay['excluded_query']
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
