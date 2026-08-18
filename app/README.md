# rusty-dot browser app

A fully client-side assembly-comparison app (Shiny for Python, deployed as a
static Shinylive/Pyodide site). Uploads never leave the browser.

## Layout

- `app.py` — UI + server wiring (thin; logic lives in `core/`)
- `core/` — pure-Python logic, unit-tested natively (`tests/test_app_*.py`)
- `www/` — custom CSS/JS
- `wheels/` — **not in git**; the rusty-dot wasm wheel is dropped here at
  export time (CI artifact from the `wasm-build` job) and installed from the
  Pyodide virtual filesystem at app startup
- `requirements.txt` — packages micropip-installs at startup (Pyodide builds)

## Run natively (development)

```bash
pip install shiny rusty-dot   # or maturin develop from the repo root
shiny run --launch-browser app/app.py
```

## Build the static site

```bash
pip install shinylive
# 1. build or download the wasm wheel (CI job "Wasm (Pyodide) Wheel Build"):
mkdir -p app/wheels && cp rusty_dot-*-cp312-cp312-emscripten_3_1_58_wasm32.whl app/wheels/
# 2. export and serve (service worker => must be http, not file://):
shinylive export app site/app
python3 -m http.server --directory site 8741
# open http://127.0.0.1:8741/app/
```

The wasm wheel must target the Pyodide release bundled by shinylive
(currently Pyodide 0.27.x = CPython 3.12 + Emscripten 3.1.58, built with the
pinned `nightly-2025-02-17` Rust toolchain); see the `wasm-build` job in
`.github/workflows/ci.yml` for the exact recipe.

## Usage notes

- **Input modes**: upload two assemblies (FASTA/FASTA.gz) *or* a
  precomputed PAF alignment. In PAF mode the aligner options are hidden; an
  optional query-assembly upload enables the reordered-FASTA download.
- **Self-alignment**: tick "Align assembly to itself" to compare one
  assembly against itself without uploading it twice.
- **Identity colouring**: for methods that report identity (minimap2,
  LASTZ, nucmer, PAF import) alignments can be coloured by % identity
  instead of forward/reverse strand.

## Alignment methods

| Method | Runs | Status |
|---|---|---|
| k-mer matching | rusty-dot wasm wheel (Pyodide) | ✅ |
| PAF import | pure Python | ✅ |
| minimap2 / LASTZ / nucmer | biowasm (Aioli WebWorker; fetched from the biowasm CDN at runtime) | ✅ |
| BLAST | — | not possible: no production WASM build of BLAST+ exists |
