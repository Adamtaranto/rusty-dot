# Web App

**rusty-dot** ships with a browser app for comparing two genome assemblies —
no installation required. The app is fully client-side: your FASTA files are
parsed, aligned, and plotted **entirely in your browser** (via WebAssembly and
[Pyodide](https://pyodide.org)) and **never leave your machine**.

<p style="margin: 1.2rem 0;">
  <a href="../app/" target="_blank" rel="noopener" class="md-button md-button--primary">
    Launch the app in a new window &rarr;
  </a>
</p>

The app opens in its own window — it needs the full viewport for the
interactive dotplot, and running it standalone keeps its memory separate
from the docs site.

## Capabilities

| Category | Options |
|---|---|
| Input | Upload query and target assemblies as FASTA or gzipped FASTA (`.fa`, `.fasta`, `.fna`, or any of these gzipped), align an assembly to itself ("Align assembly to itself"), or import a pre-computed PAF file |
| Alignment — offline | **k-mer matching** with rusty-dot's rolling-hash engine (the Rust crate compiled to a wasm wheel, running in Pyodide); **PAF import** (pure Python). Both work fully offline once the app has loaded. The k-mer method also offers a compute-time "min match block" filter (default 50 bp) that drops short blocks before they are materialised — changing it re-runs the comparison. |
| Alignment — via biowasm | **minimap2 2.22** (preset `asm20` by default), **LASTZ 1.04.52** (`--step=20 --notransition` by default), and **nucmer (MUMmer4)** (`-l 100 -c 200` by default) run in-browser via the [biowasm](https://biowasm.com) CDN — assembly-scale defaults; the tool binaries are fetched from biowasm.com at runtime, so these need a network connection. Runs can be cancelled mid-flight (switch method or re-run), and each run's exact command line and stderr appear in the collapsible **Aligner log**. |
| Plot reconfiguration (no recompute) | Sort contigs by size or maximise colinearity, auto-flip reverse-oriented contigs, hide internal axes, minimum alignment length filter, identity colouring (tool/PAF results only — k-mer matches are exact, so the option is hidden there) |
| Annotations | Upload GFF3 tracks for either assembly; per-type toggles and colours, diagonal shading, side tracks in the focused pair view, clickable features |
| Downloads | Plot as SVG or PDF, alignments as PAF, and the reordered/reoriented query assembly as FASTA |
| Diagnostics | Live wasm-heap memory readout in the sidebar; animated processing indicator while anything is computing |

## Limits and notes

- **Memory**: browser WebAssembly is capped at roughly 2 GB, so the app is
  best suited to viral, bacterial, and fungal-scale assemblies. For large
  plant or animal genomes, use rusty-dot [locally](installation.md) instead.
  The sidebar shows the live wasm-heap usage against this cap.
- **k-mer input limit**: the k-mer method refuses inputs beyond **~80 Mb
  combined** (and warns above 40 Mb that big runs take minutes) — its index
  would exhaust the browser heap beyond that. minimap2 / LASTZ / nucmer
  handle larger inputs; they warn above ~200 MB combined.
- **First load**: the app downloads roughly 45 MB of Python/wasm runtime
  assets on first visit; they are cached by the browser, so subsequent loads
  are fast.
- **Why no BLAST?** No production WebAssembly build of NCBI BLAST+ exists, so
  BLAST cannot run in the browser. MUMmer4's `nucmer` is offered as the
  closest substitute for sensitive genome-vs-genome alignment.

## Run it locally

The app is a standard Shiny for Python app and can be run natively (with the
installed rusty-dot instead of the wasm wheel) or exported as a static
Shinylive site. See
[`app/README.md`](https://github.com/Adamtaranto/rusty-dot/blob/main/app/README.md)
in the repository for the development and export recipes.
