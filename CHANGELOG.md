# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Browser app: the deploy pipeline could bundle stale wasm wheels restored
  from the CI cargo cache, and the app installed whichever wheel sorted
  last — a mismatched Emscripten tag that micropip rejects, breaking app
  startup. The workflows now clear `target/wheels` before building and only
  publish the matching tag, and the app selects the bundled wheel by the
  running Pyodide's platform tag.

### Added

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
