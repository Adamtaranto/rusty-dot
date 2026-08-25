"""GenBank flat-file parsing from in-memory bytes.

The browser app accepts annotated GenBank assemblies directly, so users do
not have to split them into FASTA + GFF by hand.  Parsing is pure Python
for the same reason as :mod:`core.fasta`: the wasm build has no native
readers, and Biopython is far too large to add to the Pyodide bundle for
what amounts to a few hundred lines of flat-file scanning.

The result is a :class:`GenbankInput` carrying

* a :class:`~core.fasta.FastaInput` built from the ORIGIN blocks, so
  everything downstream (aligners, the k-mer index, the dataset cache)
  consumes exactly what it already consumes for a FASTA upload; and
* GFF3 **text** for the FEATURES tables, handed to
  :meth:`rusty_dot.annotation.GffAnnotation.from_text` rather than
  constructing ``GffFeature`` objects here.  That keeps this module free
  of any ``rusty_dot`` import (matching the rest of ``app/core``), makes
  the tests plain string comparisons, and inherits attribute
  percent-decoding and 1-based → 0-based coordinate conversion from the
  already-tested parser.

Only what a dotplot needs is kept.  ``/translation`` qualifiers are
dropped (megabytes on a real genome, useless here) and the whole-record
``source`` feature is skipped by default, since it spans the entire
sequence and would shade the entire diagonal.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import logging

from core.fasta import FastaInput, content_digest

logger = logging.getLogger(__name__)

_GZIP_MAGIC = b'\x1f\x8b'

#: Qualifiers never worth carrying into a dotplot annotation.
_SKIP_QUALIFIERS = frozenset({'translation'})

#: Qualifiers copied into the GFF attributes column, in this order.
_KEPT_QUALIFIERS = ('locus_tag', 'product', 'note', 'db_xref', 'gene')

#: Qualifiers whose continuation lines join without a space (sequence-like).
_UNSPACED_QUALIFIERS = frozenset({'translation'})

#: Column where a FEATURES key starts, and where its location begins.
_FEATURE_KEY_COL = 5
_FEATURE_LOC_COL = 21

#: Characters that must be percent-encoded in a GFF3 attribute value.
_GFF_ESCAPES = {
    '%': '%25',
    ';': '%3B',
    '=': '%3D',
    '&': '%26',
    ',': '%2C',
    '\t': '%09',
    '\n': '%0A',
    '\r': '%0D',
}


@dataclass(frozen=True)
class Location:
    """A parsed GenBank feature location.

    Attributes
    ----------
    parts : tuple[tuple[int, int], ...]
        0-based half-open ``(start, end)`` intervals, in the order they
        appear in the location string.
    strand : str
        ``'+'`` or ``'-'``.
    fuzzy : bool
        Whether any bound was marked ``<`` or ``>`` (partial feature).
    ordered : bool
        ``True`` for ``order(...)`` (parts whose join is not asserted) as
        opposed to ``join(...)``.  Drawn identically; recorded so it can
        be surfaced as an attribute.
    """

    parts: tuple[tuple[int, int], ...]
    strand: str
    fuzzy: bool
    ordered: bool


@dataclass(frozen=True)
class GenbankInput:
    """A parsed GenBank upload.

    Attributes
    ----------
    fasta : FastaInput
        Sequences from the ORIGIN blocks, keyed by the same names used in
        *gff_text* so the annotation lines resolve against the contigs.
    gff_text : str
        GFF3 text for every feature, ready for
        :meth:`~rusty_dot.annotation.GffAnnotation.from_text`.  Empty
        (header only) when the file carries no usable features.
    n_features : int
        Number of GenBank features emitted (a multi-part feature counts
        once, though it produces one GFF line per part).
    """

    fasta: FastaInput
    gff_text: str
    n_features: int


def _decompress(data: bytes) -> bytes:
    """Transparently gunzip *data* when it carries the gzip magic bytes."""
    if data[:2] != _GZIP_MAGIC:
        return data
    try:
        return gzip.decompress(data)
    except (OSError, EOFError) as exc:
        raise ValueError(f'Could not decompress gzip data: {exc}') from exc


def _split_records(lines: list[str]) -> list[tuple[int, list[str]]]:
    """Split a flat file into per-record line blocks.

    Returns ``(first_line_number, lines)`` per record.  A record runs from
    its ``LOCUS`` line to the terminating ``//``; a final record with no
    ``//`` is accepted rather than silently dropped.
    """
    records: list[tuple[int, list[str]]] = []
    current: list[str] = []
    start_line = 0
    for i, line in enumerate(lines, 1):
        if line.startswith('LOCUS') and not current:
            start_line = i
        if line.strip() == '//':
            if current:
                records.append((start_line, current))
                current = []
            continue
        if start_line:
            current.append(line)
    if current:
        records.append((start_line, current))
    return records


def _record_name(lines: list[str], line_no: int) -> str:
    """Pick the sequence name for a record.

    VERSION is preferred over ACCESSION over the LOCUS name, because that
    is what tends to appear in companion files.  This is the single source
    of truth for the name: it labels both the FASTA record and the
    ``seqname`` column of every GFF line emitted for the record, and a
    mismatch between the two would make the features silently vanish from
    the plot (the plotter only logs a warning for unknown sequences).
    """
    accession = locus = None
    for line in lines:
        if line.startswith('VERSION'):
            parts = line.split()
            if len(parts) > 1:
                return parts[1]
        elif line.startswith('ACCESSION') and accession is None:
            parts = line.split()
            if len(parts) > 1:
                accession = parts[1]
        elif line.startswith('LOCUS') and locus is None:
            parts = line.split()
            if len(parts) > 1:
                locus = parts[1]
    name = accession or locus
    if not name:
        raise ValueError(f'GenBank record at line {line_no} has no usable name')
    return name


def _declared_length(lines: list[str]) -> int | None:
    """Read the sequence length declared on the LOCUS line, if present."""
    for line in lines:
        if not line.startswith('LOCUS'):
            continue
        tokens = line.split()
        for i, tok in enumerate(tokens):
            if tok in ('bp', 'aa') and i and tokens[i - 1].isdigit():
                return int(tokens[i - 1])
        return None
    return None


def parse_location(text: str) -> Location:
    """Parse a GenBank location string.

    Handles ``a..b``, a bare position, ``a^b`` (a point between two
    bases), ``<a`` / ``b>`` fuzzy bounds, and nested ``complement(...)``,
    ``join(...)`` and ``order(...)`` operators.  Remote references
    (``ACC:1..10``) are skipped: their bases are not in this file.

    Parameters
    ----------
    text : str
        The location field, with whitespace already collapsed.

    Returns
    -------
    Location
        Parts as 0-based half-open intervals, in file order.

    Raises
    ------
    ValueError
        If the location cannot be parsed.
    """
    state = {'fuzzy': False, 'ordered': False, 'remote': 0}
    parts = _parse_location_expr(text.strip(), state)
    if state['remote']:
        logger.warning(
            'Skipped %d remote location reference(s) in %r', state['remote'], text
        )
    elif not parts:
        # No parts and nothing skipped means the string itself was empty
        # or degenerate (e.g. 'join()') -- that is a real parse failure.
        raise ValueError(f'Location {text!r} selects no bases in this record')
    return Location(
        parts=tuple(parts),
        strand='-' if state.get('complement') else '+',
        fuzzy=bool(state['fuzzy']),
        ordered=bool(state['ordered']),
    )


def _split_top_level(text: str) -> list[str]:
    """Split a comma list, ignoring commas inside nested parentheses."""
    out: list[str] = []
    depth = 0
    start = 0
    for i, ch in enumerate(text):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        elif ch == ',' and depth == 0:
            out.append(text[start:i])
            start = i + 1
    out.append(text[start:])
    return [p for p in (s.strip() for s in out) if p]


def _parse_location_expr(text: str, state: dict) -> list[tuple[int, int]]:
    """Recursive-descent worker for :func:`parse_location`."""
    text = text.strip()
    if not text:
        return []
    for op in ('complement', 'join', 'order'):
        prefix = op + '('
        if text.startswith(prefix) and text.endswith(')'):
            inner = text[len(prefix) : -1]
            if op == 'complement':
                # Nested complement(join(...)) flips the whole group; a
                # doubly-complemented location is back to the plus strand.
                state['complement'] = not state.get('complement', False)
            else:
                state['ordered'] = state['ordered'] or op == 'order'
            parts: list[tuple[int, int]] = []
            for piece in _split_top_level(inner):
                parts.extend(_parse_location_expr(piece, state))
            return parts
    if ':' in text:  # ACC:1..10 -- bases live in another record
        state['remote'] += 1
        return []
    return [_parse_span(text, state)]


def _parse_bound(token: str, state: dict) -> int:
    """Parse one 1-based bound, recording ``<``/``>`` fuzziness."""
    token = token.strip()
    if token[:1] in ('<', '>'):
        state['fuzzy'] = True
        token = token[1:]
    if not token.isdigit():
        raise ValueError(f'Unparsable location bound {token!r}')
    return int(token)


def _parse_span(text: str, state: dict) -> tuple[int, int]:
    """Parse a single ``a..b`` / ``a^b`` / ``n`` span to a half-open pair."""
    if '..' in text:
        lo_s, hi_s = text.split('..', 1)
        lo, hi = _parse_bound(lo_s, state), _parse_bound(hi_s, state)
    elif '^' in text:
        # A site *between* two bases has no width; represent it as the
        # single base before the join so it remains drawable.
        lo_s, _hi_s = text.split('^', 1)
        lo = _parse_bound(lo_s, state)
        hi = lo
    else:
        lo = hi = _parse_bound(text, state)
    if hi < lo:
        raise ValueError(f'Location {text!r} ends before it starts')
    return lo - 1, hi  # 1-based inclusive -> 0-based half-open


def _feature_blocks(lines: list[str]) -> list[tuple[str, str, list[str], int]]:
    """Extract ``(key, location, qualifier_lines, line_no)`` per feature.

    The FEATURES table runs until ORIGIN / CONTIG / BASE COUNT or the end
    of the record.  A new feature starts wherever the key column is
    non-blank; a location continues across lines that are neither a new
    key nor a ``/qualifier``.
    """
    blocks: list[tuple[str, str, list[str], int]] = []
    in_features = False
    key = loc = None
    quals: list[str] = []
    start_line = 0

    def flush() -> None:
        if key:
            blocks.append((key, loc or '', quals[:], start_line))

    for i, raw in enumerate(lines, 1):
        line = raw.rstrip('\n').rstrip('\r')
        if not in_features:
            in_features = line.startswith('FEATURES')
            continue
        if line[:1] not in ('', ' ') or line.startswith(
            ('ORIGIN', 'CONTIG', 'BASE COUNT')
        ):
            break
        if not line.strip():
            continue
        new_key = line[_FEATURE_KEY_COL:_FEATURE_LOC_COL].strip()
        rest = line[_FEATURE_LOC_COL:].strip()
        if new_key:
            flush()
            key, loc, quals, start_line = new_key, rest, [], i
        elif rest.startswith('/'):
            quals.append(rest)
        elif quals:
            quals.append(rest)  # continuation of the previous qualifier
        elif key:
            loc = (loc or '') + rest  # continuation of the location
    flush()
    return blocks


def _parse_qualifiers(qual_lines: list[str]) -> list[tuple[str, str]]:
    """Parse ``/key="value"`` lines into ordered ``(key, value)`` pairs.

    A list rather than a dict so repeated qualifiers (``/db_xref``) all
    survive.  Values may span lines; sequence-like qualifiers join without
    a separator, everything else with a space.
    """
    pairs: list[tuple[str, str]] = []
    key: str | None = None
    chunks: list[str] = []
    quoted = False

    def flush() -> None:
        if key is None:
            return
        joiner = '' if key in _UNSPACED_QUALIFIERS else ' '
        value = joiner.join(c for c in chunks if c)
        pairs.append((key, value.strip().strip('"')))

    for line in qual_lines:
        if line.startswith('/'):
            if quoted:
                # An unterminated quote swallowing the next /qualifier is
                # far worse than ending it here.
                flush()
            elif key is not None:
                flush()
            body = line[1:]
            if '=' in body:
                key, first = body.split('=', 1)
                quoted = first.startswith('"') and not (
                    first.endswith('"') and len(first) > 1
                )
                chunks = [first]
            else:
                key, chunks, quoted = body, [], False
        elif key is not None:
            chunks.append(line)
            if quoted and line.rstrip().endswith('"'):
                quoted = False
    flush()
    return pairs


def _origin_sequence(lines: list[str]) -> str:
    """Extract the ORIGIN sequence block, stripped of offsets and spaces."""
    out: list[str] = []
    in_origin = False
    for raw in lines:
        line = raw.strip()
        if not in_origin:
            in_origin = line.startswith('ORIGIN')
            continue
        if line.startswith('//'):
            break
        out.append(''.join(c for c in line if c.isalpha()))
    return ''.join(out).upper()


def _gff_escape(value: str) -> str:
    """Percent-encode the characters GFF3 reserves in attribute values."""
    return ''.join(_GFF_ESCAPES.get(ch, ch) for ch in value)


def _attributes_for(
    key: str, index: int, quals: list[tuple[str, str]], loc: Location
) -> str:
    """Build the GFF3 column-9 string for one GenBank feature.

    Both ``ID`` and ``Parent`` are always emitted: ``GffAnnotation.iter_groups``
    only joins the parts of a multi-part feature when a record carries
    *both*, so a ``join(...)`` CDS would otherwise draw as unrelated
    fragments with no connector.  When the record has no natural parent
    name a synthetic one is used -- it is never rendered, and it is far
    less invasive than relaxing the library's grouping rule.
    """
    qmap: dict[str, list[str]] = {}
    for qk, qv in quals:
        qmap.setdefault(qk, []).append(qv)

    feature_id = f'{key}_{index}'
    parent = (qmap.get('locus_tag') or qmap.get('gene') or [f'{feature_id}_parent'])[0]
    name = (qmap.get('gene') or qmap.get('product') or qmap.get('locus_tag') or [None])[
        0
    ]

    attrs = [f'ID={_gff_escape(feature_id)}', f'Parent={_gff_escape(parent)}']
    if name:
        attrs.append(f'Name={_gff_escape(name)}')
    for qk in _KEPT_QUALIFIERS:
        for qv in qmap.get(qk, []):
            if qv:
                attrs.append(f'{qk}={_gff_escape(qv)}')
    if loc.fuzzy:
        attrs.append('partial=true')
    if loc.ordered:
        attrs.append('order=true')
    return ';'.join(attrs)


def parse_genbank_bytes(data: bytes, include_source: bool = False) -> GenbankInput:
    """Parse a GenBank flat file (optionally gzipped) from raw bytes.

    Parameters
    ----------
    data : bytes
        Raw uploaded file content.
    include_source : bool, optional
        Keep the whole-record ``source`` feature.  Off by default: it
        spans the entire sequence, so it shades the whole diagonal and
        fills a track lane end to end.

    Returns
    -------
    GenbankInput
        Sequences and GFF3 text for the file.

    Raises
    ------
    ValueError
        On unreadable input: bad gzip, non-UTF-8 bytes, no records, a
        record with no sequence, duplicate record names, or an unparsable
        feature location.
    """
    raw = _decompress(data)
    try:
        text = raw.decode('utf-8')
    except UnicodeDecodeError as exc:
        raise ValueError(f'GenBank file is not UTF-8 text: {exc}') from exc

    blocks = _split_records(text.splitlines())
    if not blocks:
        raise ValueError('No GenBank records found (no LOCUS line).')

    records: list[tuple[str, str]] = []
    gff_lines = ['##gff-version 3']
    seen: set[str] = set()
    n_features = 0

    for line_no, lines in blocks:
        name = _record_name(lines, line_no)
        if name in seen:
            raise ValueError(
                f'Duplicate record name {name!r} (line {line_no}); '
                'record names must be unique.'
            )
        seen.add(name)

        seq = _origin_sequence(lines)
        if not seq:
            raise ValueError(
                f'GenBank record {name!r} (line {line_no}) has no ORIGIN sequence.'
            )
        declared = _declared_length(lines)
        if declared is not None and declared != len(seq):
            # Warn, do not raise: a truncated or hand-edited file is still
            # worth plotting, and the real sequence is what we drew from.
            logger.warning(
                'GenBank record %r declares %d bp but ORIGIN holds %d; '
                'using the ORIGIN sequence.',
                name,
                declared,
                len(seq),
            )
        records.append((name, seq))

        for index, (key, loc_text, quals, feat_line) in enumerate(
            _feature_blocks(lines), 1
        ):
            if key == 'source' and not include_source:
                continue
            try:
                loc = parse_location(loc_text)
            except ValueError as exc:
                raise ValueError(
                    f'Unparsable location {loc_text!r} in feature {key!r} at '
                    f'line {line_no + feat_line - 1}: {exc}'
                ) from exc
            if not loc.parts:
                continue  # entirely remote: its bases are in another record
            pairs = [
                (k, v) for k, v in _parse_qualifiers(quals) if k not in _SKIP_QUALIFIERS
            ]
            attributes = _attributes_for(key, index, pairs, loc)
            n_features += 1
            for start, end in loc.parts:
                gff_lines.append(
                    '\t'.join(
                        (
                            name,
                            'GenBank',
                            key,
                            str(start + 1),  # GFF is 1-based inclusive
                            str(max(end, start + 1)),
                            '.',
                            loc.strand,
                            '.',
                            attributes,
                        )
                    )
                )

    logger.info(
        'Parsed %d GenBank record(s), %d residues, %d feature(s)',
        len(records),
        sum(len(s) for _, s in records),
        n_features,
    )
    return GenbankInput(
        fasta=FastaInput(digest=content_digest(data), records=tuple(records)),
        gff_text='\n'.join(gff_lines) + '\n' if n_features else '',
        n_features=n_features,
    )
