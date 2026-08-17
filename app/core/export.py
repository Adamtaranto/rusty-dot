"""Reordered / reoriented FASTA export from in-memory sequences.

The k-mer method stores sequences inside the ``CrossIndex`` (which has its
own ``write_fasta``), but PAF-import and external-tool methods only carry
coordinates — for those the app keeps the parsed upload
(:class:`core.fasta.FastaInput`) and exports directly from its records
using this module.
"""

from __future__ import annotations

from typing import Iterable, Sequence


def reordered_fasta_text(
    records: Sequence[tuple[str, str]],
    order: Iterable[str],
    reverse: set[str],
    line_width: int = 60,
) -> str:
    """Serialise sequence records to FASTA in a given order and orientation.

    Contigs named in *order* are written first, in that order; any remaining
    contigs from *records* follow in their original order (so an ordering
    derived from alignments — which may not cover every contig — still
    exports the complete assembly).  Contigs named in *reverse* are written
    reverse-complemented with a ``reverse_complement`` note in the header,
    matching :meth:`rusty_dot.paf_io.CrossIndex.write_fasta` output.

    Parameters
    ----------
    records : sequence of (str, str)
        ``(name, sequence)`` pairs, e.g. :attr:`core.fasta.FastaInput.records`.
    order : iterable of str
        Contig names in the desired output order.  Names not present in
        *records* are ignored.
    reverse : set[str]
        Contig names to reverse-complement.
    line_width : int, optional
        Wrap sequence lines at this many bases; ``0`` or negative writes
        each sequence on one line.  Default is ``60``.

    Returns
    -------
    str
        The FASTA text.
    """
    from rusty_dot.paf_io import reverse_complement  # noqa: PLC0415 - lazy

    seq_map = dict(records)
    names = [n for n in order if n in seq_map]
    seen = set(names)
    names.extend(n for n, _ in records if n not in seen)

    chunks: list[str] = []
    for name in names:
        seq = seq_map[name]
        if name in reverse:
            seq = reverse_complement(seq)
            chunks.append(f'>{name} reverse_complement\n')
        else:
            chunks.append(f'>{name}\n')
        if line_width and line_width > 0:
            for i in range(0, len(seq), line_width):
                chunks.append(seq[i : i + line_width] + '\n')
        else:
            chunks.append(seq + '\n')
    return ''.join(chunks)
