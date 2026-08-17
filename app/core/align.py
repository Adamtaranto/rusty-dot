"""Alignment-method registry and result construction.

Two method families exist:

* **In-Pyodide** — k-mer matching via the rusty-dot Rust core, and PAF
  import (parsed directly from uploaded text).  Implemented here.
* **biowasm tools** — minimap2 / LASTZ / nucmer run in an Aioli WebWorker
  on the JS side; their text output is handed back to Python and converted
  to ``PafRecord`` objects by parsers in this module.  (Tool execution and
  output parsers land with the biowasm work unit; the registry already
  declares them so the UI can present a stable method list.)

NCBI BLAST is deliberately absent: no production WASM build of BLAST+
exists, so it cannot run client-side.  ``nucmer`` (MUMmer4) is offered as
the closest substitute.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from rusty_dot.paf_io import PafAlignment

#: Method key -> UI label for every planned alignment method.
METHOD_LABELS: dict[str, str] = {
    'kmer': 'k-mer matching (rusty-dot)',
    'paf': 'Import PAF file',
    'minimap2': 'minimap2 (biowasm)',
    'lastz': 'LASTZ (biowasm)',
    'nucmer': 'nucmer / MUMmer4 (biowasm)',
}

#: Methods implemented so far; the UI greys out the rest.
AVAILABLE_METHODS: frozenset[str] = frozenset({'kmer', 'paf'})


def paf_alignment_from_text(text: str) -> PafAlignment:
    """Build a ``PafAlignment`` from PAF text.

    Parameters
    ----------
    text : str
        PAF-format content (e.g. an uploaded ``.paf`` file or minimap2
        stdout).  Blank lines and ``#`` comment lines are skipped.

    Returns
    -------
    rusty_dot.paf_io.PafAlignment
        Alignment holding one record per PAF line.

    Raises
    ------
    ValueError
        If no PAF records are found, or a line is not valid PAF (propagated
        from ``PafRecord.from_line`` with the line number added).
    """
    from rusty_dot.paf_io import PafAlignment, PafRecord

    records = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        try:
            records.append(PafRecord.from_line(line))
        except (ValueError, IndexError) as exc:
            raise ValueError(f'Invalid PAF at line {lineno}: {exc}') from exc
    if not records:
        raise ValueError('No PAF records found in input')
    return PafAlignment.from_records(records)


def paf_text_from_alignment(alignment: PafAlignment) -> str:
    """Serialise a ``PafAlignment`` back to PAF text.

    Parameters
    ----------
    alignment : rusty_dot.paf_io.PafAlignment
        Alignment to serialise.

    Returns
    -------
    str
        One PAF line per record, newline-terminated.
    """
    return '\n'.join(rec.to_line() for rec in alignment.records) + '\n'
