"""Parity tests for the Rust and Python gravity-order implementations.

The collinearity sort is implemented twice: in Rust
(:meth:`rusty_dot.SequenceIndex.optimal_contig_order`, used by ``CrossIndex``)
and in pure Python (:func:`rusty_dot.paf_io.compute_gravity_contigs`, used by
``PafAlignment``).  Both must produce identical orderings for the same matches.

To make the two engines genuinely comparable, the PAF records fed to the Python
path are derived directly from the Rust index's own stranded matches, so both
sides operate on exactly the same alignment blocks.
"""

from rusty_dot._rusty_dot import SequenceIndex
from rusty_dot.paf_io import PafRecord, compute_gravity_contigs


def _records_from_index(
    idx: SequenceIndex,
    query_names: list[str],
    target_names: list[str],
) -> list[PafRecord]:
    """Build PAF records from the index's own merged stranded matches.

    Parameters
    ----------
    idx : SequenceIndex
        Populated index.
    query_names, target_names : list[str]
        Sequence names to compare pairwise.

    Returns
    -------
    list[PafRecord]
        One record per merged match block across every (query, target) pair.
    """
    records: list[PafRecord] = []
    for q in query_names:
        q_len = idx.get_sequence_length(q)
        for t in target_names:
            t_len = idx.get_sequence_length(t)
            for qs, qe, ts, te, strand in idx.compare_sequences_stranded(q, t, True):
                block = max(qe - qs, te - ts)
                line = (
                    f'{q}\t{q_len}\t{qs}\t{qe}\t{strand}\t'
                    f'{t}\t{t_len}\t{ts}\t{te}\t{block}\t{block}\t255'
                )
                records.append(PafRecord.from_line(line))
    return records


def _make_index() -> SequenceIndex:
    """Build a multi-sequence index with matches spread across targets."""
    idx = SequenceIndex(k=11)
    a = (
        'ACGTTGCAAGGCCTTAGCTAGGATCCGATCGATTACGGCATGCATTGCACGTAGCTAGCATCG'
        'TTAGGCATCCGATTGACCGATACGGATTCAGCTAGGCATTACGGATCCGATTAGCACGTATGC'
    )
    b = (
        'GGATCCTTACGAGCATTGCACCGATTAGGCATCGATCGATTAGCACGTATGCATTGCAAGGCC'
        'TTAGCTAGGATCCGATCGATTACGGCATGCATTGCACGTAGCTAGCATCGTTAGGCATCCGAT'
    )
    idx.add_sequence('tA', a)
    idx.add_sequence('tB', b)
    idx.add_sequence('qA', a)  # maps to tA
    idx.add_sequence('qB', b)  # maps to tB
    idx.add_sequence('qAB', a[:60] + b[:60])  # maps to both, tA dominant
    return idx


def _assert_parity(idx, query_names, target_names):
    rust_q, rust_t = idx.optimal_contig_order(list(query_names), list(target_names))
    records = _records_from_index(idx, query_names, target_names)
    py_q, py_t, _reversed = compute_gravity_contigs(
        records, list(query_names), list(target_names)
    )
    assert rust_q == py_q, f'query order differs: rust={rust_q} python={py_q}'
    assert rust_t == py_t, f'target order differs: rust={rust_t} python={py_t}'


def test_parity_two_targets_multi_contig():
    """Rust and Python agree on a multi-target, multi-contig layout."""
    idx = _make_index()
    _assert_parity(idx, ['qB', 'qAB', 'qA'], ['tB', 'tA'])


def test_parity_with_unmatched_contig():
    """Parity holds when an unmatched contig must sort last by length."""
    idx = _make_index()
    idx.add_sequence('qNone', 'A' * 150)  # no shared k-mers → unmatched
    _assert_parity(idx, ['qNone', 'qA', 'qB'], ['tA', 'tB'])


def test_parity_single_target():
    """Parity holds in the degenerate single-target case."""
    idx = _make_index()
    _assert_parity(idx, ['qB', 'qA', 'qAB'], ['tA'])
