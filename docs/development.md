# Development Guide

This page explains how to set up a local development environment for rusty-dot so you can edit both the Rust extension and the Python package and run the full test suite.

## Prerequisites

### Install Rust

rusty-dot requires a working Rust toolchain (stable channel).
Install it with [rustup](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
rustup update stable
```

Verify the installation:

```bash
rustc --version
cargo --version
```

### Create a conda environment

A minimal conda environment file is provided at `environment.yml`.
It pins Python to 3.13 and installs the key Python dependencies.

```bash
conda env create -f environment.yml
conda activate rustydot
```

Alternatively, use any Python ≥ 3.9 virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
```

## Fork and clone the repository

[Fork rusty-dot](https://github.com/Adamtaranto/rusty-dot/fork) on GitHub, then clone your fork:

```bash
git clone https://github.com/<your-username>/rusty-dot.git
cd rusty-dot
git remote add upstream https://github.com/Adamtaranto/rusty-dot.git
```

## Install in development mode

[maturin](https://www.maturin.rs/) is the build backend for the Rust extension.
Install it and then build the package in editable mode:

```bash
pip install maturin
# Build Rust extension + install package in editable mode
maturin develop --extras dev,docs
```

The `--extras dev,docs` flag installs all optional development and documentation
dependencies declared in `pyproject.toml`.

Verify the install:

```python
import rusty_dot
print(rusty_dot.__version__)
```

!!! note
    Always re-run `maturin develop` after modifying any Rust source files
    (`src/*.rs`).  Pure Python changes in `python/rusty_dot/` take effect
    immediately without a rebuild.

## Installing pre-commit hooks

rusty-dot uses [pre-commit](https://pre-commit.com/) to enforce code quality
on every commit.  Install the hooks once after cloning:

```bash
pip install pre-commit   # already included in the dev extras
pre-commit install                 # install commit-stage hooks
pre-commit install --hook-type pre-push   # install push-stage hooks
```

The hooks include:

| Hook | Stage | What it checks |
|------|-------|---------------|
| `check-ast` / `check-yaml` / `check-toml` | commit | File syntax |
| `trailing-whitespace` / `end-of-file-fixer` | commit | Whitespace |
| `cargo fmt --check` | commit | Rust formatting |
| `cargo check` | commit | Rust compilation |
| `cargo clippy -- -D warnings` | commit | Rust linting |
| `ruff format` / `ruff check` | commit | Python formatting & linting |
| `pydocstyle` (`.pyi` stubs only) | commit | Docstring style |
| `cargo test --lib` | push | Rust unit tests |
| `pytest tests/ -x -q` | push | Python test suite |

Run all hooks manually at any time:

```bash
pre-commit run --all-files
```

## Running Python tests

```bash
pytest tests/ -v
```

To run a specific test file or test function:

```bash
pytest tests/test_index.py -v
pytest tests/test_index.py::test_get_paf_all_returns_paf_lines -v
```

## Python code style — ruff

rusty-dot uses [ruff](https://docs.astral.sh/ruff/) for Python linting and
formatting (configured in `pyproject.toml`).

**Check for issues:**

```bash
ruff check python/ tests/
```

**Auto-fix fixable issues:**

```bash
ruff check --fix python/ tests/
```

**Format code:**

```bash
ruff format python/ tests/
```

## Rust code quality

### Format

```bash
cargo fmt
```

Check formatting without modifying files (also used by the pre-commit hook):

```bash
cargo fmt --all -- --check
```

### Lint (Clippy)

```bash
cargo clippy -- -D warnings
```

### Compile check

Quickly verify the crate compiles without producing a binary:

```bash
cargo check
```

### Rust unit tests

Run the Rust-side library tests:

```bash
cargo test --lib
```

Run all Rust tests (including integration tests, if any):

```bash
cargo test
```

## Building the documentation locally

Install the docs dependencies (included in `pip install maturin --extras dev,docs`):

```bash
mkdocs serve
```

Open <http://127.0.0.1:8000> in your browser.  The site rebuilds automatically
when you save a documentation file.

To build a static site:

```bash
mkdocs build
```
