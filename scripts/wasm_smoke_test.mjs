// Smoke-test the wasm (Pyodide) wheel under Node.
//
// Usage: node scripts/wasm_smoke_test.mjs target/wheels/rusty_dot-*.whl
//
// Requires the `pyodide` npm package at the version matching the wheel's
// platform tag (Pyodide 0.27.x for emscripten_3_1_58 wheels):
//   npm install pyodide@0.27.7
//
// The test loads Pyodide, installs the wheel from the virtual filesystem via
// micropip (which validates the platform tag exactly as the browser would),
// imports rusty_dot, and runs a tiny k-mer comparison through the Rust core.
import { readFileSync } from 'node:fs';
import { basename } from 'node:path';
import { loadPyodide } from 'pyodide';

const wheelPath = process.argv[2];
if (!wheelPath) {
  console.error('usage: node scripts/wasm_smoke_test.mjs <wheel.whl>');
  process.exit(2);
}

const py = await loadPyodide();
// matplotlib/numpy are imported by rusty_dot's __init__; load the prebuilt
// Pyodide packages rather than letting micropip try to resolve them.
await py.loadPackage(['micropip', 'numpy', 'matplotlib']);

const wheelName = basename(wheelPath);
py.FS.writeFile(`/tmp/${wheelName}`, readFileSync(wheelPath));
const micropip = py.pyimport('micropip');
await micropip.install(`emfs:/tmp/${wheelName}`);

await py.runPythonAsync(`
import rusty_dot

idx = rusty_dot.SequenceIndex(8)
idx.add_sequence('a', 'ACGTTGCAAGGCTTAA' * 32)
idx.add_sequence('b', 'ACGTTGCAAGGCTTAA' * 32)
matches = idx.compare_sequences_stranded('a', 'b')
assert matches, 'expected self-similarity matches from the Rust core'

# The pure-Python layer on top of the Rust index must work too.
ci = rusty_dot.CrossIndex(8)
ci.add_sequence('q1', 'ACGTTGCAAGGCTTAA' * 32, group='q')
ci.add_sequence('t1', 'ACGTTGCAAGGCTTAA' * 32, group='t')
ci.compute_matches(query_group='q', target_group='t')
print('rusty_dot wasm OK:', len(matches), 'match blocks')
`);
console.log(`smoke test passed for ${wheelName}`);
