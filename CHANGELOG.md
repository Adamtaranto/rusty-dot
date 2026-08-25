# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added (contig length filter)

- A **Min contig length** setting leaves contigs shorter than the given
  size out of the plot. One panel per contig means a couple of
  chromosomes can be buried under hundreds of short scaffolds, each drawn
  as a sliver. Excluded contigs are still written to the
  reordered-FASTA download — the plot is filtered, the data is not — and
  a threshold that would empty an axis leaves that axis whole instead of
  rendering an empty grid.

### Changed (sidebar layout and plot navigation)

- Ticking **Align assembly to itself** now hides the target annotation
  upload as well as the target assembly: both axes are the query assembly,
  so there is no target to annotate. PAF input keeps the upload, having
  both roles but no self-align notion.
- Multi-panel plots: clicking a panel centres it as before, but **any**
  subsequent click in the figure returns to the full view. Previously the
  second click had to land on the same panel — every listener was attached
  per panel and stopped propagation, so no click elsewhere could reset.
  Clicks that a match, annotation or side-track claims still show their
  details rather than throwing the zoom away.

### Added (annotation handling)

- A **Clear annotations** button in the GFF section drops every uploaded
  annotation and resets the file inputs. Shiny cannot clear a file input
  from the server, so this also resets the widget itself — without which
  the filename would linger and re-picking the same file would never
  reload it.
- Uploading a GFF whose sequence names are absent from the corresponding
  assembly now warns, naming the offending contigs (or, when nothing
  matches at all, suggesting the likely cause: wrong file, or assigned to
  the wrong role). Previously the features simply never appeared, which
  reads as a broken plot — the library logged a warning the browser user
  never saw.

### Fixed (sidebar and multi-panel labels)

- Choosing a GFF file no longer scrolls the sidebar back to the top. The
  jump happens during the file-picker interaction, before any re-render,
  so the position is captured when the picker is opened and re-asserted
  once the panel settles; a real scroll gesture cancels it.
- Multi-panel plots: per-row contig labels are angled rather than
  vertical, and shrink to fit their row. A vertical label is as tall as
  the contig name is long, so an assembly mixing one large chromosome with
  small contigs had the short rows' labels overrun into their neighbours'
  and smear together.

### Fixed (single-panel plots and annotation controls)

- Clicking a dot plot no longer dims the panel it just selected. The axes
  background was gid-tagged `rd-panel-<r>-<c>-bg`, and matplotlib nests
  gid'd artists, so it matched every `g[id^="rd-panel-"]` selector: the
  report counted it as an extra panel (defeating the guard that disables
  click-to-focus on single-panel reports) and dimmed the clicked panel's
  parent. On a multi-panel grid, clicking blank area dimmed *every* panel.
  The background now uses its own `rd-plotbg-` prefix, and the report
  matches panel ids on their exact shape rather than the prefix alone.
- Double-clicking a plot opens the drill-down again. The same id collision
  made the app's bridge resolve to the background group, whose id failed
  the strict panel regex, so no message was ever posted and the drill-down
  silently never opened — on any grid size, though it was most visible on
  a single-panel plot where nearly all of the panel is background.
- Feature-type toggles and colours survive pressing "Run comparison". The
  controls are rendered dynamically, and the block also drew the
  "Shade features on diagonal" checkbox, whose visibility depends on the
  result — so every run rebuilt the whole list at its defaults. The two
  global toggles are now static (shown or hidden through a conditional
  panel), leaving the type list dependent only on the uploaded
  annotations. Re-running in GenBank mode no longer re-registers an
  unchanged annotation source either, which also preserves the
  Annotations tab's per-feature show/hide and colour overrides.
- The navigation hint no longer advertises "click panel = focus" on a
  single-panel plot, where click-to-focus is deliberately inert.

### Added (GenBank input)

- Browser app: **Assemblies (GenBank)** input mode reads sequences *and*
  annotations from the same upload, so annotated GenBank assemblies no
  longer have to be split into FASTA + GFF by hand. `app/core/genbank.py`
  is a pure-Python flat-file parser (Biopython is far too large for the
  Pyodide bundle) handling gzip, multi-record files, CRLF, multi-line
  qualifiers, and the location grammar — `join`/`order`/`complement`
  including nesting, `<a..>b` fuzzy bounds, `a^b` sites, bare positions,
  and remote `ACC:1..10` references (skipped with a warning rather than
  failing the upload). `/translation` is dropped and the whole-record
  `source` feature is skipped by default.
- A role's annotations are now a list of sources tagged by filename and
  kind, so a GenBank file and a GFF upload coexist and are merged.
  Overlaps are not silently deduped — that would discard features the user
  deliberately supplied — so the app warns instead.

### Added (drill-down annotations)

- Browser app: an **Annotations** tab in the focused-pair view lists every
  feature on the two sequences with its type, name, coordinates (1-based
  inclusive, as in the source file), strand, source file and remaining
  attributes, plus a per-feature show/hide toggle and colour override, a
  client-side filter, and bulk show/hide/reset actions.
- Interactive report: clicking a side-track feature bands the matching
  column (x track) or row (y track) behind the matches, and shows the
  feature in the detail bar. Bands toggle, several can be active at once,
  `Shift`+click bands every part of a multi-part feature (a spliced CDS),
  and `Esc` clears them. The detail bar now also reports the source file.
- Aligner runs that outlive their time budget now prompt instead of
  failing: a **Still aligning** dialog offers "Wait another 5 minutes"
  (re-arming the watchdog, and offerable again) or "Cancel run". The
  budget became a watchdog that only notifies, so `minimap2 -c` on large
  genomes is no longer killed mid-run.
- `GffFeature` gained `color` (per-feature colour override, honoured by the
  side tracks and the diagonal squares) and `source_file` (provenance for
  merged uploads). Both default to unset, so existing callers are
  unaffected.
- HTML report payload gained an optional top-level `tracks` key describing
  every drawn side-track part, with its SVG gid and a `group` shared by the
  parts of one feature.

### Changed (annotation colours)

- Browser app: feature-type colours are now assigned once across *both*
  uploads instead of per file. `GffAnnotation` colours types by their index
  into that file's own sorted type list, so the same type landed on
  different palette entries in the query and target — `gene` green on one
  track and blue on the other — and `CDS` vs `cds` were unrelated keys.
- Types are matched case-insensitively with a small alias table, so
  `repeat_region` / `TE` / `LTR` share one colour. `gene`, `CDS`, `mRNA`,
  `exon` and `repeat` get reserved conventional colours (green / yellow /
  maroon / grey / red); the rest walk a 24-entry palette and then
  deterministic hash-derived colours, so nothing silently repeats.
- The two per-role colour pickers collapse into one **Feature types** list
  with `Q`/`T` badges — once the colour is shared, two controls for one
  value would need a feedback-breaking two-way sync.
- **Shade features on diagonal** is only offered when the panel grid
  actually contains a self-comparison panel; it could never draw anything
  on a plain cross-assembly comparison. On a self-alignment the query and
  target annotations are now merged for shading, instead of the target's
  being ignored whenever a query GFF existed.
- Interactive report: the sequence copy buttons say what the next press
  will do — **Fetch query seq** / **Fetch target seq** until the sequence
  has been fetched from the app, **Copy …** afterwards — replacing the
  mystery "Press again to copy" state.

### Fixed (annotation tracks)

- Unstranded features (`repeat_region`, `TE`, `TIR`, `TRF_SSR`, `flank`,
  `Pseudogene`, …) rendered as full-height blocks covering the entire
  y-axis annotation band. `FancyBboxPatch` applies `rounding_size` on x and
  `rounding_size × mutation_aspect` on y, but the y-orientation branch
  passed the along-sequence radius (in bases) as `rounding_size` even
  though that axis is measured in lane units, then inverted the aspect on
  top — on a 236 kb contig, a ~1200-lane corner radius in x and ~10⁷ bp in
  y. The radius is now expressed per axis. The x track is unchanged.

### Changed (documentation build)

- The docs site is now built with [Zensical](https://zensical.org) instead of
  MkDocs + Material (the same authors' successor project). `mkdocs.yml` is
  replaced by `zensical.toml`; the site URLs, the `site/` output directory and
  the `/app/` Shinylive export path are unchanged.
- Tutorial notebooks are pre-rendered to Markdown by the new
  `scripts/notebooks_to_md.py` (Zensical has no notebook support). The
  `.ipynb` files remain the source of truth and are still not executed at
  build time; the generated `.md` and image directories are gitignored.
- `docs` extra: `zensical`, `mkdocstrings-python` and `nbconvert` replace
  `mkdocs`, `mkdocs-material`, `mkdocstrings[python]` and `mkdocs-jupyter`.
- CI gained a `Docs Build` job, so a broken docs build now fails on the pull
  request rather than after merge.
- `docs/api/style.md` was orphaned (not in the nav) and is now published as
  **Plot Style**; a few leftover Sphinx `:class:` roles in the API pages were
  turned into working links.

### Added (browser-app build docs)

- `environment-wasm.yml`: conda environment for building and serving the
  browser app, pinning Emscripten 3.1.58 from conda-forge (no `emsdk` clone
  needed) alongside `shiny`, `shinylive` and `maturin`. Kept separate from
  `environment.yml`, whose `python-freethreading` pin the wasm build cannot
  use.

### Fixed (docs)

- `app/README.md` cited the wrong pinned Rust nightly (`nightly-2025-02-17`);
  the toolchain required by `.cargo/config.toml` and CI is
  `nightly-2025-05-01`. The local build recipe is now complete and verified
  end to end: nightly install, Emscripten setup, the version pins that must
  agree (Pyodide 0.27.x / CPython 3.12 / Emscripten 3.1.58), the wheel-tag
  check, `scripts/add_loading_splash.py`, and a troubleshooting table.

### Added (identity & aligned sequences)

- Browser app: minimap2 gained a **Base-level alignment (`-c`)** option that
  adds `cg:Z` CIGAR, `NM:i` and `de:f` tags to the PAF output. (`-L` was
  evaluated and is unnecessary: it only affects SAM/BAM CIGAR storage, not
  PAF stdout.)
- Browser app: clicking an identity-coloured match now fetches and displays
  the gapped query/match/target alignment in the detail bar, built
  server-side from the record's CIGAR and the uploaded assemblies (sequences
  are sliced on demand — nothing is pre-extracted or embedded in the
  report). Requires a minimap2 run with `-c` (or an uploaded PAF with `cg`
  tags) plus the source FASTAs.
- `PafRecord.identity`: fraction identity choosing the most accurate
  available metric — `1 - de` (gap-compressed, from `-c`/`--cs` output),
  else gap-compressed identity derived from the CIGAR, else the BLAST-style
  `residue_matches / alignment_block_len` estimate from the required PAF
  columns.
- `rusty_dot.alignment_view.aligned_text()`: render a gapped pairwise
  alignment view (query / match line / target) from a PAF record's CIGAR.
- Browser app: matches without a CIGAR (minimap2 without `-c`, nucmer,
  k-mer) show the raw query and target sequences as two unaligned lines in
  the detail bar instead of nothing, each truncated at 20,000 bases.
- Browser app: sequence copy buttons in the match detail bar copy the full
  (untruncated) sequence of the selected match to the clipboard. They are
  labelled **Fetch …** until the sequence has been retrieved and **Copy …**
  afterwards (see *Changed (annotation colours)*).

### Performance (match detail)

- Browser app: clicking a match no longer ships the full query/target
  sequences with the preview (up to ~0.7 s on 8 Mbp matches, ~60–80 ms per
  Mbp); the detail reply now carries only the clipped preview, and the copy
  buttons fetch the full sequence on demand, caching it client-side so
  repeat copies are instant. Preview slices and the gapped-alignment
  builder also clip to the 20,000-column window before any copy/revcomp,
  so megabase matches never materialise full slices for display.

### Changed (identity & PAF round-trip)

- Identity colouring (`color_by_identity=True`) now uses
  `PafRecord.identity`, preferring the gap-compressed `de` tag over the PAF
  column 10/11 estimate when base-level alignment output is available.
- `PafRecord.to_line()` now preserves optional SAM-style tags (with their
  original type characters), so exported and cached PAF retains
  `cg`/`NM`/`de` instead of being truncated to the 12 required columns.

### Fixed (focused-view canvas)

- Focused single-pair (drill-down) views with extreme sequence-length
  ratios (tiny query vs long target) no longer let the title overlap the
  plot area or truncate the y-axis contig label: the canvas reserves
  absolute inch-based margins around the untouched panel (bp-per-inch
  parity preserved), the title keeps a fixed physical offset, an
  over-tall rotated y label falls back to horizontal, and ultra-thin
  axes show a single end tick instead of overlapping position labels.

### Changed (aligner options)

- Browser app: the minimap2 `-k`/`-w`/`-m` inputs are pre-filled with the
  selected preset's actual values (asm5/asm10: k=19 w=19; asm20: k=19
  w=10; m=40) and refresh when the preset changes, replacing the
  "0 = preset default" sentinel.

### Removed

- Browser app: LASTZ was removed as an alignment method (its wasm build was
  impractically slow on assembly-scale inputs compared with minimap2 and
  nucmer). Precomputed LASTZ alignments can still be imported as PAF.

### Added (app UI batch)

- Browser app: every aligner and plot setting gained a hoverable ⓘ info
  icon explaining what the option does.
- Browser app: repeat-genome aligner options — minimap2 `-m` (min chaining
  score), `-P` (retain all chains) and `-D` (skip self-diagonal matches;
  offered in self-alignment mode), and nucmer `--nosimplify` alongside the
  existing `--maxmatch`.
- Focused single-pair (drill-down) views now label the axes with the
  contig names — left of the y axis, below the x axis, offset past the
  GFF side annotation tracks when present — and axis ticks read in
  bp/Kbp/Mbp units instead of scientific notation.
- Multi-panel grids use one bp/Kbp/Mbp unit across every contig (chosen
  from the longest, announced once per axis as "Position (Mbp)") and
  angle the x tick labels 45° so long position labels no longer overlap.

### Changed (app UI batch)

- Browser app: the navigation tips moved from above the plot to below it,
  with each action term (scroll, drag, Esc, …) in bold; the focused view
  omits the panel-click tips that don't apply there.

### Performance (app rendering)

- Browser app: line width and min match length now update **instantly**
  inside the embedded interactive report (an injected CSS rule and
  per-segment hiding, driven over postMessage) instead of re-rendering
  the whole matplotlib figure in Pyodide — previously a multi-second
  wait per tweak on real assemblies. The static (non-interactive) plot
  and SVG/PDF downloads still honour both options server-side.
- Browser app: computed contig orders are cached per ordering mode for
  the current result, so re-selecting "maximise colinearity" after
  trying another mode no longer recomputes the gravity sort; the
  plotter construction (which copies the full record list) is also
  memoised per result.

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

### Changed (app restyle)

- Browser app: restyled around a CSS custom-property theme layer (teal
  accent matching the docs site) with light/dark parity and smooth
  transitions on buttons and inputs, plus a wordmark treatment for the
  header.
- Browser app: an animated "Processing…" pill appears whenever the app is
  computing (pure CSS on Shiny's busy state — covers aligner runs, index
  builds and report rendering).
- Browser app: live wasm-heap memory readout in the sidebar footer
  (refreshes every 5 s; Pyodide only) so users can see how close they are
  to the ~2 GB browser cap.
- Browser app: the Nature-journal styling option was removed from the app
  (the library's `nature_style()` context manager is unchanged).
- Docs: the Web App page now opens the app in its own window via a
  launch button instead of embedding it in an 800 px iframe.

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
