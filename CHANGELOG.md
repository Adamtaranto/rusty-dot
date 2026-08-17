# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

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
