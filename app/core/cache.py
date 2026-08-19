"""Session-scoped caches for alignment results.

Alignments are expensive (k-mer indexing in wasm, or an external wasm tool
run); plot reconfiguration is cheap.  The cache keys every result by the
content digests of the uploaded files plus the method parameters, so
changing a display option never recomputes an alignment, and re-running
with identical inputs is instant.
"""

from __future__ import annotations

from collections import OrderedDict
import logging
from typing import TYPE_CHECKING, Any

from .fasta import FastaInput

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from rusty_dot import CrossIndex
    from rusty_dot.paf_io import PafAlignment

logger = logging.getLogger(__name__)

#: Group labels used for the two assemblies inside every ``CrossIndex``.
QUERY_GROUP = 'query'
TARGET_GROUP = 'target'


def _params_key(params: dict[str, Any]) -> tuple[tuple[str, Any], ...]:
    """Return a hashable, order-independent form of a parameter dict.

    Parameters
    ----------
    params : dict[str, Any]
        Method parameters (must have hashable values).

    Returns
    -------
    tuple[tuple[str, Any], ...]
        Sorted ``(name, value)`` pairs.
    """
    return tuple(sorted(params.items()))


class SessionCache:
    """Per-session cache of k-mer indices and alignment results.

    One instance lives for the lifetime of a browser session.  K-mer
    ``CrossIndex`` objects are cached per ``(k, query digest, target
    digest)`` so per-contig k-mer indices are reused across runs; imported
    or tool-generated ``PafAlignment`` objects are cached per ``(method,
    params, digests)``.

    Both caches are bounded LRUs: a ``CrossIndex`` holds full sequence
    copies plus the match cache, so on real assemblies a single entry can
    reach hundreds of megabytes — unbounded growth would exhaust the wasm
    heap long before the user noticed.
    """

    #: Maximum retained k-mer indices (each can hold two full assemblies).
    KMER_MAX = 3
    #: Maximum retained alignment results (coordinate records only).
    PAF_MAX = 8

    def __init__(self) -> None:
        """Initialise empty caches."""
        self._kmer: OrderedDict[tuple[Any, ...], CrossIndex] = OrderedDict()
        self._paf: OrderedDict[tuple[Any, ...], PafAlignment] = OrderedDict()

    @staticmethod
    def _touch(store: OrderedDict, key: tuple[Any, ...], maxsize: int) -> None:
        """Mark *key* most-recently-used and evict the oldest beyond *maxsize*."""
        store.move_to_end(key)
        while len(store) > maxsize:
            evicted, _ = store.popitem(last=False)
            logger.info('Evicted cache entry %s', evicted)

    def kmer_index(
        self,
        k: int,
        query: FastaInput,
        target: FastaInput,
        merge: bool = True,
        min_block_len: int = 0,
    ) -> CrossIndex:
        """Return a ``CrossIndex`` with matches computed, building on miss.

        Parameters
        ----------
        k : int
            K-mer size.
        query : FastaInput
            Parsed query assembly (rows / y-axis).
        target : FastaInput
            Parsed target assembly (columns / x-axis).
        merge : bool, optional
            Merge adjacent k-mer hits into runs.  Default is ``True``.
        min_block_len : int, optional
            Drop match blocks shorter than this many bases at compute time
            (filtered natively, before records are materialised).  Repeat-rich
            assembly pairs can otherwise produce millions of records that
            exhaust the wasm heap.  Default is ``0`` (keep all).

        Returns
        -------
        rusty_dot.CrossIndex
            Index containing both assemblies (groups ``'query'`` and
            ``'target'``) with the stranded match cache populated.
        """
        from rusty_dot import CrossIndex

        key = (k, merge, min_block_len, query.digest, target.digest)
        cached = self._kmer.get(key)
        if cached is not None:
            logger.info('k-mer index cache hit for %s', key)
            self._touch(self._kmer, key, self.KMER_MAX)
            return cached
        logger.info(
            'Building k-mer index (k=%d) for %d query / %d target contigs',
            k,
            len(query.records),
            len(target.records),
        )
        index = CrossIndex(k)
        for name, seq in query.records:
            index.add_sequence(name, seq, group=QUERY_GROUP)
        for name, seq in target.records:
            index.add_sequence(name, seq, group=TARGET_GROUP)
        index.compute_matches(
            query_group=QUERY_GROUP,
            target_group=TARGET_GROUP,
            merge=merge,
            min_block_len=min_block_len,
        )
        self._kmer[key] = index
        self._touch(self._kmer, key, self.KMER_MAX)
        return index

    def get_paf(
        self, method: str, params: dict[str, Any], *digests: str
    ) -> PafAlignment | None:
        """Look up a cached ``PafAlignment`` for a method/params/inputs combo.

        Parameters
        ----------
        method : str
            Alignment method name (e.g. ``'paf'``, ``'minimap2'``).
        params : dict[str, Any]
            Method parameters used to produce the alignment.
        *digests : str
            Content digests of every input file involved.

        Returns
        -------
        rusty_dot.paf_io.PafAlignment or None
            The cached alignment, or ``None`` on miss.
        """
        key = (method, _params_key(params), *digests)
        cached = self._paf.get(key)
        if cached is not None:
            self._touch(self._paf, key, self.PAF_MAX)
        return cached

    def put_paf(
        self,
        method: str,
        params: dict[str, Any],
        alignment: PafAlignment,
        *digests: str,
    ) -> None:
        """Store a ``PafAlignment`` under a method/params/inputs key.

        Parameters
        ----------
        method : str
            Alignment method name.
        params : dict[str, Any]
            Method parameters used to produce the alignment.
        alignment : rusty_dot.paf_io.PafAlignment
            The result to cache.
        *digests : str
            Content digests of every input file involved.
        """
        key = (method, _params_key(params), *digests)
        self._paf[key] = alignment
        self._touch(self._paf, key, self.PAF_MAX)
        logger.info('Cached %s alignment (%d record(s))', method, len(alignment))
