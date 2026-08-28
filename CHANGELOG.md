# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — library

- Interactive single-file HTML dotplot reports (`DotPlotter.to_html()`, or an
  `.html` output path): click a panel to focus it, scroll/drag to pan and
  zoom, drag a box to zoom to a region, click a match for coordinates,
  identity and sequence details.
- Match selection in reports: a selected match stays highlighted after its
  detail bar is closed, and its query/target ranges are drawn as translucent
  bands through every panel in its grid row and column. `Shift`+drag
  box-selects every match intersecting the box (testing the drawn line, not
  its bounding box) without opening details; `Esc` clears selections first
  and only then resets/exits.
- Annotation support: GFF3 parsing upgrades (percent-decoded attributes,
  `feature_id`/`parent`/`name`, `from_text`/`from_bytes` with gzip
  detection, embedded `##FASTA`, multi-part grouping); diagonal feature
  squares drawn behind the matches with an automatic type legend; lane-packed
  side annotation tracks on focused plots (strand arrows, multi-part CDS
  connectors); clickable features in reports, with track features banding
  their row/column — bands survive into saved SVG/PNG/PDF figures.
  `GffFeature` gained `color` and `source_file`; the report payload gained a
  `tracks` key and per-segment `uid`s.
- Identity and aligned sequences: `PafRecord.identity` (gap-compressed `de`
  tag, else CIGAR-derived, else the BLAST-style estimate);
  `alignment_view.aligned_text()` for gapped query/match/target views;
  `to_line()` preserves optional SAM-style tags; identity colouring uses the
  best available metric and `plot(identity_colorbar=True)` appends a colour
  key.
- Plotting: `contig_order='length'|'colinearity'` with `auto_reverse=True`,
  `hide_internal_axes=True`, opt-in `nature_style()` styling, and a
  reordering/reorientation tutorial notebook.
- `SequenceIndex.compare_pairs_stranded`: batched both-strand comparison in
  one native call (GIL released) with a `min_block_len` filter;
  `CrossIndex.compute_matches` gains the same option and computes the whole
  grid in one batched call.

### Added — browser app

- Fully client-side assembly comparison (`app/`) built with Shiny for Python
  and deployed as a static Shinylive/Pyodide site alongside the docs. Upload
  FASTA/FASTA.gz assemblies, align with rusty-dot's k-mer engine or in-browser
  biowasm aligners (minimap2 2.22, nucmer/MUMmer4), reconfigure the dotplot
  without recomputing, and download SVG/PDF plots, PAF, and a reordered query
  FASTA. Includes a loading splash with staged progress, background aligner
  pre-download, and a wasm-heap memory readout.
- Input modes: PAF as a top-level mode (with contig-name validation against
  an optional query assembly), self-alignment ("Align assembly to itself",
  which clears/hides target annotations), and GenBank assemblies parsed by a
  pure-Python flat-file parser (gzip, multi-record, full location grammar) so
  sequences and annotations come from one upload.
- Annotations: GFF3/GenBank uploads per role with merged multi-source
  handling, per-feature-type visibility toggles and colour pickers gated
  behind an **Apply changes** button (pre-first-plot edits apply on render),
  self-panel shading, side tracks in the drill-down, a **Clear annotations**
  reset, and warnings when uploaded sequence names match nothing.
- Drill-down **Annotations** tab: lists every feature on both focused
  sequences (contig column included, self-panels deduplicated) with
  per-feature show/hide and colour overrides held until applied, header
  sorting, multiple per-column text filters with an Apply button, and bulk
  actions. Rows are built client-side, so switching tabs is instant. Click a
  row to select it (⌘/Ctrl-click for several); selected features are banded
  on the Plot tab, and the bands survive re-renders and exports.
- Interactive drill-down: double-click a panel for the focused single-pair
  view with contig-name axis labels and bp/Kbp/Mbp ticks; clicking an
  identity-coloured match shows the gapped CIGAR alignment (minimap2 gained
  a **Base-level alignment (`-c`)** option), other matches show their raw
  sequences; copy buttons fetch and cache full sequences on demand.
- Fullscreen plot mode: an expand button at the top-left of the plot panel
  requests true browser fullscreen (falling back to a window-filling
  overlay), scales the figure to the viewport, hides the hints, memory
  readout and tab strip, and keeps **Back to overview** working without
  leaving fullscreen. With a selection active, `Esc` deselects first; a
  second `Esc` exits.
- Manual dark-mode toggle in the header (defaults to the system preference).
- Aligner UX: per-stage progress messages surfaced at the right end of the
  header bar, a collapsible aligner log with exact command lines, in-flight
  runs cancelled on method switch or re-run, per-tool timeouts with a
  **Still aligning** wait/cancel prompt, repeat-genome options (minimap2
  `-m`/`-P`/`-D`, nucmer `--nosimplify`), presets pre-filling their real
  `-k`/`-w`/`-m` values, hoverable ⓘ help on every option, and a size guard
  that steers >80 Mb k-mer runs to minimap2 instead of exhausting the wasm
  heap.
- A **Min contig length** setting hides short contigs from the plot (the
  reordered-FASTA download keeps them).
- Cmd/Ctrl+click adds matches to (or removes them from) the selection one
  by one, in the overview and the drill-down. Selections — matches and the
  selected annotation feature — persist across overview <-> drill-down
  swaps: the app holds them in view-independent terms (names + data
  coordinates) and re-applies them to each rebuilt report, keeping
  selections on pairs a narrower view does not show.
- Upload-size guidance, from empirically measured ceilings (synthetic
  pairs in Chrome): beyond ~80 Mb combined the k-mer method is removed
  from the selector (90 Mb completes at a 2.9 GB heap peak; 100 Mb
  silently kills the Python runtime), with a pointer to the tutorial
  notebooks for running the library locally; above ~200 MB the app
  suggests uploading a precomputed PAF instead (minimap2 completed at
  200 MB, nucmer at 250 MB; 300 MB crashed the browser tab for both).

### Changed

- Releases publish real binary wheels: maturin builds for Linux
  (x86_64/aarch64), macOS (arm64/x86_64) and Windows (x64) across CPython
  3.12–3.14 and 3.14t, plus an sdist. The free-threaded wheel declares
  `gil_used = false` and CI fails if importing it re-enables the GIL. A
  PEP 783 `pyemscripten_2026_0_wasm32` wheel is published to PyPI for
  browser runtimes; the Shinylive app keeps its separate legacy
  `emscripten_3_1_58_wasm32` build (Pyodide 0.27.7).
- App identity: restyled around a CSS custom-property theme layer (teal
  accent, light/dark parity) with a bundled display font, a
  "rusty·dot — live assembly comparison" wordmark, a dot-grid page
  background, a translucent sidebar, and a matching loading splash. The
  progress popup, title-bar pulse and corner "Processing…" pill were
  replaced by the single header task status plus spinner.
- Feature-type colours are assigned once across both uploads
  (case-insensitive with aliases, conventional colours for common types),
  shown in a single **Feature types** list with `Q`/`T` badges;
  **Shade features on diagonal** only appears when a self-comparison panel
  exists.
- Navigation tips moved below the plot with bold action terms, tailored to
  the view; report copy buttons are labelled **Fetch …**/**Copy …** by what
  the next press does.
- Assembly-scale aligner defaults (nucmer `-l 100 -c 200`); docs build moved
  from MkDocs + Material to Zensical with pre-rendered notebooks and a CI
  docs-build job; the docs Web App page opens the app in its own window.

### Fixed

- Drag-zoom in the embedded overview lands exactly where the box was
  drawn. A figure taller than the iframe kept its intrinsic size and
  scrolled, so the zoomed viewBox filled a mostly off-screen element and
  appeared shifted up; embedded reports now fit the SVG to the iframe
  viewport (as fullscreen already did).
- The memory readout and docs no longer claim a 2 GB wasm cap: the
  Pyodide heap measurably grows to 4096 MB (~4.0 GB allocatable before a
  clean MemoryError), and the first-load download is 26 MB compressed
  (37 MB of assets), replacing the stale 30/45 MB claims.
- Sequence previews no longer go stale after re-running the aligner with
  different settings: the per-pair record cache is cleared on every new
  result, so a run without `-c` stops showing the previous run's base-level
  alignments.
- The focused-view title is centred over the dot-plot panel rather than the
  whole figure (the identity colorbar and legend shift the panel
  off-centre).
- The aligner log no longer covers the report's bottom detail bar in the
  drill-down (tab panes now form the same flex column as the overview).
- Annotations survive app startup in the browser: the wasm-wheel install and
  the annotation priming raced, caching an import error for the whole
  session.
- Panel-id collisions no longer break click-to-focus, dimming or
  double-click drill-down (axes backgrounds use their own `rd-plotbg-`
  prefix and panel ids are matched exactly).
- Feature-type toggles and colours survive "Run comparison"; annotation
  edits no longer rebuild the figure once per click; toggling six types
  costs one redraw.
- Unstranded track features no longer render as full-height blocks (the
  corner radius is now computed per axis).
- Extreme-ratio focused views keep the title, y label and ticks legible
  without distorting bp-per-inch parity; short contig row labels overhang
  instead of collapsing to an ellipsis.
- The sidebar keeps its scroll position through file-picker interactions;
  file inputs reset properly after **Clear annotations**.
- The deploy pipeline no longer bundles stale wasm wheels (build dirs are
  cleared and the app selects the wheel by the running Pyodide's platform
  tag).

### Performance

- The k-mer index was restructured into a canonical-hash CSR table
  (~13–16 bytes/bp vs ~100) with a sorted two-pointer match walk: on a real
  37 Mb × 37 Mb fungal pair, index build 6.4 s → 3.4 s, matching
  20 s → 3.5 s, peak memory roughly quartered — and the k-mer method now
  completes in the browser instead of exhausting the 2 GB wasm heap.
  `reorder_for_colinearity` reuses cached match records (72 s → 11 s).
- Browser app: `parse_fasta_bytes` rewritten single-pass (~7.5× faster);
  line width and min match length apply instantly inside the report over
  postMessage instead of re-rendering matplotlib; contig orders and plotter
  construction are memoised per result; assemblies ship to the biowasm
  aligners once and are referenced by digest; session caches are bounded
  LRUs.
- Match details ship only the clipped preview (full sequences are fetched on
  demand and cached), so megabase matches never materialise full slices for
  display.

### Removed

- LASTZ as an in-browser alignment method (its wasm build was impractically
  slow at assembly scale); precomputed LASTZ alignments still import as PAF.

## [0.1.0] - 2026-02-25

### Added

- Initial release: rolling-hash ntHash k-mer indexing, both-strand k-mer
  matching, strand-aware run merging, PAF output, `CrossIndex` cross-assembly
  comparisons, collinearity contig ordering with reordered/reoriented FASTA
  export, and matplotlib dotplot visualisation via `DotPlotter`.

[Unreleased]: https://github.com/Adamtaranto/rusty-dot/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Adamtaranto/rusty-dot/releases/tag/v0.1.0
