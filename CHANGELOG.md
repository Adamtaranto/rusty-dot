# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added (GFF annotations)

- GFF3 parsing upgrades: `parse_attributes()` with percent-decoding,
  `GffFeature.feature_id`/`.parent`/`.name` properties,
  `GffAnnotation.from_text()`/`.from_bytes()` constructors (gzip detected
  by magic bytes; `from_file` also reads `.gz`), embedded `##FASTA`
  section handling, and multi-part feature grouping via `iter_groups()`.
- Dotplot annotation overlays reworked: diagonal feature squares now draw
  **behind** the alignment segments and mirror correctly on
  reverse-displayed contigs; a feature-type colour legend is added
  automatically (`annotation_legend=False` to disable).
- Side annotation tracks on focused single-pair plots
  (`plot(annotation_tracks=True)` with `annotation_query=` /
  `annotation_target=`): lane-packed feature shapes with strand arrows
  for gene/mRNA/exon/CDS/ORF, rounded rectangles for unstranded types,
  and connector lines joining multi-part CDS groups.  `plot_single` now
  renders its tracks through the same engine.
- Interactive HTML reports: diagonal annotation features are clickable —
  the detail bar shows the feature's name, type, coordinates, strand,
  parent and source.
- README section and dotplot-tutorial examples for annotation overlays.
- Browser app: GFF3 annotation uploads for the query and/or target
  assembly (`.gff`/`.gff3`/`.gz`) with per-feature-type visibility
  toggles and colour pickers (native colour inputs via a small Shiny
  binding). Features shade self-vs-self panels behind the alignments
  (including cross-group self-alignments), the focused drill-down view
  adds side annotation tracks, and features are clickable in the
  interactive report.

### Performance

- The per-sequence k-mer index was restructured from two hash maps
  (~100 bytes per base pair) into a single canonical-hash CSR table
  (~13-16 bytes/bp, both strands included), and match computation now
  intersects indexes with a sorted two-pointer walk emitting raw hit pairs
  instead of building per-k-mer String maps.  On a real 37 Mb × 37 Mb
  fungal assembly pair: index build 6.4 s → 3.4 s, match computation
  20 s → 3.5 s, peak memory roughly quartered — and the k-mer method now
  completes in the browser (~5 min upload → interactive plot) where it
  previously aborted the 2 GB wasm heap fatally. Match output is unchanged
  (byte-verified, parity-tested against exact search).

- `SequenceIndex.compare_pairs_stranded`: batched both-strand comparison of
  many sequence pairs in a single native call (GIL released, parallel on
  native builds), with an optional `min_block_len` filter that drops short
  match blocks before they are materialised as Python objects — repeat-rich
  assembly pairs can otherwise produce millions of records.
  `CrossIndex.compute_matches` gains the same `min_block_len` option and now
  computes the whole Q×T grid through one batched call.
- `CrossIndex.reorder_for_colinearity` reuses the match records cached by
  `compute_matches` instead of recomputing the full match grid (72 s → 11 s
  on a real 37 Mb × 37 Mb fungal assembly pair, and 0.2 s when combined with
  `min_block_len=50`).
- Browser app: `parse_fasta_bytes` rewritten as a single-pass byte-level
  parser — ~7.5× faster on a 37 Mb assembly (754 ms → 101 ms) with peak
  allocations reduced from ~4× to ~1.5× the file size.
- Browser app: display-only plot edits (line width, min length, palettes)
  no longer re-run the contig-ordering computation; assemblies are shipped
  to the biowasm aligners once per upload and referenced by digest instead
  of re-copied across the Pyodide/JS boundary on every run; session caches
  are bounded LRUs; and the k-mer method gains a compute-time "min match
  block" filter (default 50 bp) that keeps repeat-rich genome pairs from
  producing millions of match records.
- Browser app: assembly-scale aligner defaults — LASTZ now defaults to
  `--step=20 --notransition` and nucmer to `-l 100 -c 200` (the previous
  small-region defaults were impractically slow on whole assemblies);
  minimap2 preset help documents the asm5/asm10/asm20 divergence bands.

### Added (app UX)

- Browser app: switching alignment method (or re-running) while a biowasm
  aligner is working now cancels the in-flight run — the tool's WebWorker
  is terminated so the computation actually stops — and run timeouts are
  per-tool (5 min minimap2, 10 min LASTZ/nucmer).
- Browser app: collapsible "Aligner log" panel showing each tool run's
  exact command line and stderr output.
- Browser app: progress messages while rusty-dot itself is working
  (parsing, index build, contig ordering, report rendering), plus a
  "mounting assemblies" stage for the biowasm tools.
- Browser app: the k-mer method now refuses inputs beyond ~80 Mb combined
  with guidance to use minimap2 (and warns above 40 Mb that big runs take
  minutes) — exceeding the 2 GB browser heap previously crashed the whole
  app irrecoverably mid-run.

### Fixed

- Interactive reports: clicking the plot in a single-panel report (the
  app's drill-down view) no longer recentres it — click-to-focus only
  applies when there are multiple panels to choose from.
- Focused single-pair views now keep the sequences' true aspect ratio:
  the annotation-track layout sizes its figure proportionally to the two
  sequence lengths, and trackless single-pair plots enforce exact
  bp-per-inch parity on both axes.
- Browser app: double-clicking a sub-plot to open the standalone
  single-pair view no longer zooms into the panel first and then swaps
  views — the single-click focus zoom now waits briefly and is cancelled
  by a double-click. Standalone HTML reports (no drill-down) keep
  instant clicks.
- Browser app: the embedded interactive report no longer shows its own
  navigation-hint header (the app-level hint bar already covers it);
  standalone HTML report exports keep their header.
- Browser app: download buttons are disabled until a comparison has produced
  content, and the reordered-FASTA download is offered only when sequences
  are available.

- Browser app: the deploy pipeline could bundle stale wasm wheels restored
  from the CI cargo cache, and the app installed whichever wheel sorted
  last — a mismatched Emscripten tag that micropip rejects, breaking app
  startup. The workflows now clear `target/wheels` before building and only
  publish the matching tag, and the app selects the bundled wheel by the
  running Pyodide's platform tag.

### Added

- Browser app: PAF is now a top-level input mode (assemblies vs alignment)
  instead of an alignment-method entry — selecting PAF hides the method and
  FASTA controls, with an optional query-assembly upload to enable the
  reordered-FASTA download.
- Browser app: self-alignment option ("Align assembly to itself") — upload
  one assembly instead of the same file twice.
- Browser app: per-stage progress messages for every alignment method,
  including live biowasm stages (runtime load, tool download, aligning).
- `DotPlotter.plot(identity_colorbar=True)`: append an identity colour key
  (0-100 %) at the right of the figure when colouring by identity; the
  browser app shows the key automatically.
- Browser app: uploading a query assembly alongside a PAF now validates the
  contig names against the PAF query column, warning when nothing matches
  (with a swapped-query/target hint), when contigs or sequences are missing
  on either side, and when names appear in both PAF columns.
- Browser app: numeric plot options (line width, min match length) are
  debounced — holding the spinner arrows re-renders once with the settled
  value instead of once per 0.1 step.
- Browser app: colour alignments by percent identity (palette selectable)
  for methods that report identity (minimap2, LASTZ, nucmer, PAF import);
  k-mer matches are exact so the option is hidden there.
- Browser app: loading splash screen with staged progress messages while the
  Pyodide runtime downloads and boots (previously a blank page for many
  seconds), and biowasm aligner binaries now pre-download in the background
  as soon as the tool is selected instead of on the first Run.
- Browser app (`app/`): fully client-side assembly comparison built with
  Shiny for Python and deployed as a static Shinylive/Pyodide site alongside
  the docs at <https://adamtaranto.github.io/rusty-dot/app/>. Upload
  FASTA/FASTA.gz assemblies, align with rusty-dot's k-mer engine or import
  PAF, reconfigure the dotplot without recomputing, and download SVG/PDF
  plots, PAF, and a reordered query FASTA.
- Browser app: additional in-browser aligners via biowasm — minimap2 2.22,
  LASTZ 1.04.52, and nucmer (MUMmer4); these fetch tool binaries from the
  biowasm CDN at runtime. NCBI BLAST is not offered — no production
  WebAssembly build of BLAST+ exists; nucmer is the substitute.
- Browser app: interactive in-app dotplot with sub-plot drill-down — click a
  sub-panel of the all-vs-all grid to focus a single query/target pair.
- Docs page for the browser app (`docs/webapp.md`) with an embedded live app,
  capabilities table, and memory/network limits; linked from the README,
  docs landing page, and site navigation.
- Pyodide wheel pipeline: the `wasm-build` CI job now targets Pyodide 0.27
  (CPython 3.12 + Emscripten 3.1.58) with proper side-module link flags and a
  Node smoke test, producing a wheel loadable in the browser; the docs deploy
  workflow builds the wheel, exports the app, and publishes both with the site.
- Interactive single-file HTML dotplot reports (`DotPlotter.to_html()`, or an
  `.html` output path): click a sub-panel to focus it, scroll to zoom, and
  click a match to show its details.
- Opt-in Nature-journal plot styling via the `nature_style()` context manager
  in `rusty_dot.style`.
- `plot(hide_internal_axes=True)` to suppress internal axis spines and ticks
  for continuous grid plots.
- Plot-time contig ordering with `plot(contig_order='length'|'colinearity')`
  and automatic flipping of reverse-oriented contigs with `auto_reverse=True`.
- Tutorial notebook on reordering and reorienting one assembly to match
  another.

## [0.1.0] - 2026-02-25

### Added

- Initial release: rolling-hash ntHash k-mer indexing, both-strand k-mer
  matching, strand-aware run merging, PAF output, `CrossIndex` cross-assembly
  comparisons, collinearity contig ordering with reordered/reoriented FASTA
  export, and matplotlib dotplot visualisation via `DotPlotter`.

[Unreleased]: https://github.com/Adamtaranto/rusty-dot/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Adamtaranto/rusty-dot/releases/tag/v0.1.0
