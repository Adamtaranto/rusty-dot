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
| Input | Upload query and target assemblies as FASTA or gzipped FASTA (`.fa`, `.fasta`, `.fna`, or any of these gzipped), or import a pre-computed PAF file |
| Alignment — offline | **k-mer matching** with rusty-dot's rolling-hash engine (the Rust crate compiled to a wasm wheel, running in Pyodide); **PAF import** (pure Python). Both work fully offline once the app has loaded. |
| Alignment — via biowasm | **minimap2 2.22**, **LASTZ 1.04.52**, and **nucmer (MUMmer4)** run in-browser via the [biowasm](https://biowasm.com) CDN. These fetch the tool binaries from biowasm.com at runtime, so they need a network connection. |
| Plot reconfiguration (no recompute) | Sort contigs by size or maximise colinearity, auto-flip reverse-oriented contigs, hide internal axes, minimum alignment length filter, identity colouring |
| Annotations | Upload GFF3 tracks for either assembly; per-type toggles and colours, diagonal shading, side tracks in the focused pair view, clickable features |
| Downloads | Plot as SVG or PDF, alignments as PAF, and the reordered/reoriented query assembly as FASTA |

## Limits and notes

- **Memory**: browser WebAssembly is capped at roughly 2 GB, so the app is
  best suited to viral, bacterial, and fungal-scale assemblies. For large
  plant or animal genomes, use rusty-dot [locally](installation.md) instead.
- **First load**: the app downloads a ~30 MB Python/wasm runtime on first
  visit; it is cached by the browser, so subsequent loads are fast.
- **Why no BLAST?** No production WebAssembly build of NCBI BLAST+ exists, so
  BLAST cannot run in the browser. MUMmer4's `nucmer` is offered as the
  closest substitute for sensitive genome-vs-genome alignment.

## Run it locally

The app is a standard Shiny for Python app and can be run natively (with the
installed rusty-dot instead of the wasm wheel) or exported as a static
Shinylive site. See
[`app/README.md`](https://github.com/Adamtaranto/rusty-dot/blob/main/app/README.md)
in the repository for the development and export recipes.
