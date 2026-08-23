"""Render a gapped pairwise alignment view from a PAF record's CIGAR.

Used by the web app to display the aligned sequences of a clicked match:
the record supplies coordinates, strand and the ``cg:Z`` CIGAR, and the
caller supplies the full forward-orientation query and target sequences.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .paf_io import _CIGAR_RE

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from .paf_io import PafRecord

_COMPLEMENT = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def revcomp(seq: str) -> str:
    """Return the reverse complement of a nucleotide sequence.

    Parameters
    ----------
    seq : str
        Nucleotide sequence (ACGTN, either case).

    Returns
    -------
    str
        The reverse complement, preserving case.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def clip_sequence(seq: str, max_len: int = 20_000) -> str:
    """Return *seq*, truncated with a notice when longer than *max_len*.

    Parameters
    ----------
    seq : str
        Sequence to clip.
    max_len : int
        Maximum number of bases to keep.

    Returns
    -------
    str
        The sequence, with ``… [truncated at N bases]`` appended when it
        was longer than *max_len*.
    """
    if len(seq) <= max_len:
        return seq
    return seq[:max_len] + f'… [truncated at {max_len:,} bases]'


def aligned_text(
    rec: PafRecord,
    query_seq: str,
    target_seq: str,
    *,
    width: int = 80,
    max_cols: int = 20_000,
) -> dict[str, Any]:
    """Build a text rendering of one alignment from its CIGAR.

    The query slice is reverse-complemented for ``-`` strand records, per
    PAF ``cg`` semantics (the CIGAR describes the query in alignment
    orientation against the forward target).

    Parameters
    ----------
    rec : PafRecord
        Record with a non-``None`` :attr:`~PafRecord.cigar`.  Coordinates
        must be the original genomic (forward-strand) values.
    query_seq : str
        Full query sequence in forward orientation.
    target_seq : str
        Full target sequence in forward orientation.
    width : int
        Wrap width in alignment columns per block.
    max_cols : int
        Truncate the rendering after this many alignment columns.

    Returns
    -------
    dict
        ``{'text': str, 'truncated': bool, 'cols': int}`` where ``text``
        holds blocks of three lines (query / match / target) separated by
        blank lines, and ``cols`` is the number of columns rendered.

    Raises
    ------
    ValueError
        If the record has no CIGAR.
    """
    if rec.cigar is None:
        raise ValueError('PafRecord has no CIGAR (cg:Z tag)')
    # Each alignment column consumes at most one base per side, so only
    # the first max_cols bases of each slice can ever be rendered — clip
    # before any copy/revcomp so megabase matches stay cheap.  In
    # alignment orientation a '-' strand query starts at its genomic end.
    q_need = min(rec.query_end - rec.query_start, max_cols)
    if rec.strand == '-':
        q = revcomp(query_seq[rec.query_end - q_need : rec.query_end])
    else:
        q = query_seq[rec.query_start : rec.query_start + q_need]
    t_need = min(rec.target_end - rec.target_start, max_cols)
    t = target_seq[rec.target_start : rec.target_start + t_need]
    q_buf: list[str] = []
    m_buf: list[str] = []
    t_buf: list[str] = []
    qi = ti = cols = 0
    truncated = False
    for length_str, op in _CIGAR_RE.findall(rec.cigar):
        n = int(length_str)
        if cols + n > max_cols:
            n = max_cols - cols
            truncated = True
        if op in 'M=X':
            qs, ts = q[qi : qi + n], t[ti : ti + n]
            if op == '=':
                mid = '|' * n
            elif op == 'X':
                mid = ' ' * n
            else:
                mid = ''.join(
                    '|' if a.upper() == b.upper() else ' ' for a, b in zip(qs, ts)
                )
            qi += n
            ti += n
        elif op == 'I':
            qs, mid, ts = q[qi : qi + n], ' ' * n, '-' * n
            qi += n
        elif op in 'DN':
            qs, mid, ts = '-' * n, ' ' * n, t[ti : ti + n]
            ti += n
        else:
            # S/H/P should not appear in PAF cg tags; skip defensively.
            continue
        q_buf.append(qs)
        m_buf.append(mid)
        t_buf.append(ts)
        cols += n
        if truncated:
            break
    qa, ma, ta = ''.join(q_buf), ''.join(m_buf), ''.join(t_buf)
    blocks = [
        '\n'.join((qa[i : i + width], ma[i : i + width], ta[i : i + width]))
        for i in range(0, len(qa), width)
    ]
    text = '\n\n'.join(blocks)
    if truncated:
        text += f'\n\n[truncated at {max_cols:,} alignment columns]'
    return {'text': text, 'truncated': truncated, 'cols': cols}
