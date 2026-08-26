#!/usr/bin/env bash
# Set up the cross-build toolchain for the PEP 783 wasm wheel published to PyPI
# (publish.yml).  This is NOT the toolchain for the browser app's wheel -- that
# one targets the Pyodide release shinylive bundles (0.27.7 = CPython 3.12 +
# Emscripten 3.1.58) and is built with a plain `maturin build`; see
# app/README.md.
#
# Everything about this build is dictated by the target Pyodide release, so the
# versions are queried from pyodide-build's cross-build environment rather than
# pinned here.  Run it from a CPython environment whose version matches the
# target Pyodide's CPython (3.14 for Pyodide 314.x) -- pyodide-build requires
# that.
#
# Two things this handles that a plain `emsdk install` + `rustup target add`
# does not:
#
#   1. Emscripten must be Pyodide's *patched* build.  A stock emscripten
#      rejects Rust's mangled symbols when linking a side module
#      ("emcc: error: invalid export name: _ZN..."); Pyodide carries
#      emsdk/patches/0002-Don-t-check-exports-for-being-valid-C-C-identifiers.
#      `pyodide build` installs and patches its own emsdk, so do not install
#      emscripten from conda-forge or setup-emsdk for this build, and do not
#      leave a stock emcc on PATH.
#   2. Pyodide ABIs built with WebAssembly exception handling (0.28/0.29, i.e.
#      Emscripten 4.0.9) need a Rust std built the same way; `rustup target
#      add wasm32-unknown-emscripten` gives the wrong one, so a prebuilt
#      sysroot replaces it.  Newer ABIs (Emscripten 5.x with Rust >= 1.93)
#      need no such sysroot and this is skipped.
#
# Usage:
#   scripts/setup_wasm_toolchain.sh [pyodide-version] [rust-toolchain|auto]
#
# Prints the toolchain to use as RUSTUP_TOOLCHAIN on the last line.
set -euo pipefail

# Defaults target the PyPI wheel: Pyodide 314.x = CPython 3.14 + Emscripten
# 5.0.3, whose wheels carry the PEP 783 pyemscripten_2026_0_wasm32 tag.
PYODIDE_VERSION="${1:-314.0.6}"
# "auto" uses whatever `pyodide config get rust_toolchain` reports.  Override
# it only when that toolchain cannot build this crate -- e.g. Pyodide 0.29
# suggests nightly-2025-02-01 (rustc 1.86), but bio 4.x calls is_multiple_of,
# stabilised in 1.87, so it would need nightly-2025-06-27 (the next toolchain
# with a matching emcc 4.0.9 wasm-eh sysroot).
RUST_TOOLCHAIN="${2:-auto}"

pip install pyodide-build
pyodide xbuildenv install "$PYODIDE_VERSION"

EMSCRIPTEN_VERSION="$(pyodide config get emscripten_version)"
SYSROOT_URL="$(pyodide config get rust_emscripten_target_url)"
if [ "$RUST_TOOLCHAIN" = "auto" ]; then
  RUST_TOOLCHAIN="$(pyodide config get rust_toolchain)"
else
  # Overriding the toolchain invalidates the URL pyodide-build reports (the
  # sysroot is built for one exact emcc/toolchain pair), so point at the
  # release for the pair actually being used.
  SYSROOT_TAG="emcc-${EMSCRIPTEN_VERSION}_${RUST_TOOLCHAIN}"
  SYSROOT_URL="https://github.com/pyodide/rust-emscripten-wasm-eh-sysroot/releases/download/${SYSROOT_TAG}/${SYSROOT_TAG}.tar.bz2"
fi

echo "Pyodide $PYODIDE_VERSION -> Emscripten $EMSCRIPTEN_VERSION, CPython $(pyodide config get python_version), Rust $RUST_TOOLCHAIN"

rustup toolchain install "$RUST_TOOLCHAIN" --profile minimal

if [ -z "$SYSROOT_URL" ]; then
  # Newer Pyodide ABIs (Emscripten 5.x with Rust >= 1.93) need no custom
  # sysroot: upstream rust-std is built with wasm exception handling.
  echo "No custom sysroot needed; using rustup's wasm32-unknown-emscripten"
  rustup target add wasm32-unknown-emscripten --toolchain "$RUST_TOOLCHAIN"
  echo "$RUST_TOOLCHAIN"
  exit 0
fi

# Replace the stock wasm32-unknown-emscripten std with the wasm-eh one built
# for this exact (emcc, toolchain) pair.  Mirrors what pyodide-build does for
# its own recipe builds.
RUSTLIB="$(dirname "$(dirname "$(rustup which --toolchain "$RUST_TOOLCHAIN" rustc)")")/lib/rustlib"
TOKEN="$RUSTLIB/wasm32-unknown-emscripten_install-url.txt"

if [ -f "$TOKEN" ] && [ "$(cat "$TOKEN")" = "$SYSROOT_URL" ]; then
  echo "wasm-eh sysroot already installed"
else
  echo "Installing wasm-eh sysroot from $SYSROOT_URL"
  rm -rf "$RUSTLIB/wasm32-unknown-emscripten"
  curl -fsSL "$SYSROOT_URL" | tar -xj -C "$RUSTLIB"
  echo "$SYSROOT_URL" > "$TOKEN"
fi

echo "$RUST_TOOLCHAIN"
