"""FASTA parsing from in-memory bytes.

The browser app receives uploads as raw bytes and the wasm build of
rusty-dot excludes its native FASTA reader (needletail), so parsing happens
here in pure Python.  Plain and gzip-compressed FASTA are supported.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
import logging

logger = logging.getLogger(__name__)

_GZIP_MAGIC = b'\x1f\x8b'


@dataclass(frozen=True)
class FastaInput:
    """A parsed FASTA upload.

    Attributes
    ----------
    digest : str
        Stable content digest of the raw uploaded bytes, used as a cache key.
    records : tuple[tuple[str, str], ...]
        Parsed ``(name, sequence)`` pairs in file order.
    """

    digest: str
    records: tuple[tuple[str, str], ...]

    @property
    def names(self) -> list[str]:
        """Return contig names in file order.

        Returns
        -------
        list[str]
            Contig names in the order they appear in the file.
        """
        return [name for name, _ in self.records]

    @property
    def total_length(self) -> int:
        """Return the total number of residues across all contigs.

        Returns
        -------
        int
            Sum of sequence lengths.
        """
        return sum(len(seq) for _, seq in self.records)


def content_digest(data: bytes) -> str:
    """Compute a short stable digest of raw file bytes.

    Parameters
    ----------
    data : bytes
        Raw file content.

    Returns
    -------
    str
        First 16 hex characters of the SHA-256 digest.
    """
    return hashlib.sha256(data).hexdigest()[:16]


#: Whitespace bytes deleted from sequence blocks in one C-level pass.
_SEQ_WHITESPACE = b'\r\n\t\x0b\x0c '


def _lineno(data: bytes, offset: int) -> int:
    """Return the 1-based line number of a byte offset (error paths only).

    Parameters
    ----------
    data : bytes
        The full decompressed input.
    offset : int
        Byte offset into *data*.

    Returns
    -------
    int
        1-based line number containing *offset*.
    """
    return data.count(b'\n', 0, offset) + 1


def parse_fasta_bytes(data: bytes) -> FastaInput:
    """Parse FASTA (optionally gzip-compressed) from raw bytes.

    Handles multi-line sequences, CRLF line endings, blank lines, and
    old-style ``;`` comment lines.  All whitespace inside sequence data is
    removed.

    The parser works directly on the raw bytes in per-record blocks — record
    boundaries are located with :meth:`bytes.find` and whitespace is stripped
    with a single :meth:`bytes.translate` pass per record — instead of
    decoding and splitting the whole file into per-line strings.  On
    multi-megabyte uploads under Pyodide this cuts peak memory from roughly
    four times the file size to about two, and avoids millions of short-lived
    line objects.

    Parameters
    ----------
    data : bytes
        Raw content of a ``.fasta``/``.fa``/``.fna`` file, or the same
        gzip-compressed (``.gz``).

    Returns
    -------
    FastaInput
        Parsed records plus a content digest of the *raw* (possibly
        compressed) bytes.

    Raises
    ------
    ValueError
        If the file is empty, is not valid FASTA (no ``>`` header, sequence
        data before the first header, duplicate or empty contig names, or a
        record with no sequence), or cannot be decoded as UTF-8.
    """
    digest = content_digest(data)
    if data[:2] == _GZIP_MAGIC:
        try:
            data = gzip.decompress(data)
        except (OSError, EOFError) as exc:
            raise ValueError(f'Could not decompress gzip input: {exc}') from exc

    # Find the first meaningful line: skip blank lines and ';' comments.
    pos = 0
    size = len(data)
    while pos < size:
        nl = data.find(b'\n', pos)
        end = size if nl == -1 else nl
        stripped = data[pos:end].strip()
        if stripped and not stripped.startswith(b';'):
            break
        pos = end + 1
    if pos >= size:
        raise ValueError('No FASTA records found in input')
    if data[pos : pos + 1] != b'>':
        nl = data.find(b'\n', pos)
        first_line = data[pos : size if nl == -1 else nl]
        try:
            first_line.decode('utf-8')
        except UnicodeDecodeError as exc:
            raise ValueError('Input is not UTF-8 text — is this a FASTA file?') from exc
        raise ValueError(
            f'Sequence data before first ">" header at line {_lineno(data, pos)}'
        )

    # Record boundaries: every '>' at the start of a line.  A '\n>' scan
    # covers CRLF ('\r\n>') and blank lines before headers too.
    offsets: list[int] = [pos]
    search = pos
    while True:
        nxt = data.find(b'\n>', search)
        if nxt == -1:
            break
        offsets.append(nxt + 1)
        search = nxt + 1
    offsets.append(size)

    records: list[tuple[str, str]] = []
    seen: set[str] = set()
    for start, stop in zip(offsets, offsets[1:]):
        block = data[start:stop]
        nl = block.find(b'\n')
        header = (block if nl == -1 else block[:nl]).strip()
        seq_block = b'' if nl == -1 else block[nl + 1 :]
        try:
            header_text = header.decode('utf-8')
        except UnicodeDecodeError as exc:
            raise ValueError('Input is not UTF-8 text — is this a FASTA file?') from exc
        name = header_text[1:].split(maxsplit=1)[0] if len(header_text) > 1 else ''
        if not name:
            raise ValueError(f'Empty contig name at line {_lineno(data, start)}')
        if name in seen:
            raise ValueError(
                f'Duplicate contig name {name!r} at line {_lineno(data, start)}'
            )
        seen.add(name)
        if b';' in seq_block:
            # Rare path: drop old-style ';' comment lines inside the record.
            seq_block = b'\n'.join(
                ln for ln in seq_block.split(b'\n') if not ln.lstrip().startswith(b';')
            )
        seq_bytes = seq_block.translate(None, _SEQ_WHITESPACE)
        if not seq_bytes:
            raise ValueError(f'Contig {name!r} has no sequence')
        try:
            seq = seq_bytes.decode('utf-8')
        except UnicodeDecodeError as exc:
            raise ValueError('Input is not UTF-8 text — is this a FASTA file?') from exc
        records.append((name, seq))

    logger.info(
        'Parsed %d FASTA record(s), %d residues total (digest %s)',
        len(records),
        sum(len(s) for _, s in records),
        digest,
    )
    return FastaInput(digest=digest, records=tuple(records))
