# Installation

## Requirements

- Python 3.9 or later
- A working Rust toolchain (for building from source)
- [maturin](https://www.maturin.rs/) ≥ 1.0

## Install from source

Clone the repository and build the Rust extension with [maturin](https://www.maturin.rs/):

```bash
git clone https://github.com/Adamtaranto/rusty-dot.git
cd rusty-dot
pip install maturin
maturin develop --release
```

The `--release` flag enables full Rust compiler optimisations, which is strongly recommended for any non-trivial dataset.

## Install Python dependencies

rusty-dot depends on:

| Package | Purpose |
|---------|---------|
| `matplotlib ≥ 3.5` | Dotplot visualisation |
| `numpy ≥ 1.21` | Array operations used by matplotlib |

These are declared as package dependencies and will be installed automatically by pip.

## Optional: documentation dependencies

To build the documentation locally:

```bash
pip install rusty-dot[docs]
mkdocs serve
```

## Verify the installation

```python
import rusty_dot
print(rusty_dot.__version__)  # 0.1.0

from rusty_dot import SequenceIndex
idx = SequenceIndex(k=10)
idx.add_sequence("test", "ACGTACGTACGT")
print(idx)  # SequenceIndex(k=10, sequences=1)
```
