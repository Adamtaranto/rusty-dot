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


def parse_fasta_bytes(data: bytes) -> FastaInput:
    """Parse FASTA (optionally gzip-compressed) from raw bytes.

    Handles multi-line sequences, CRLF line endings, and blank lines.
    Sequence characters are kept as-is apart from surrounding whitespace.

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
    try:
        text = data.decode('utf-8')
    except UnicodeDecodeError as exc:
        raise ValueError('Input is not UTF-8 text — is this a FASTA file?') from exc

    names: list[str] = []
    seqs: list[list[str]] = []
    seen: set[str] = set()
    for lineno, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith('>'):
            name = line[1:].split(maxsplit=1)[0] if len(line) > 1 else ''
            if not name:
                raise ValueError(f'Empty contig name at line {lineno}')
            if name in seen:
                raise ValueError(f'Duplicate contig name {name!r} at line {lineno}')
            seen.add(name)
            names.append(name)
            seqs.append([])
        elif line.startswith(';'):
            # Old-style FASTA comment line.
            continue
        else:
            if not names:
                raise ValueError(
                    f'Sequence data before first ">" header at line {lineno}'
                )
            seqs[-1].append(line)

    if not names:
        raise ValueError('No FASTA records found in input')
    records: list[tuple[str, str]] = []
    for name, parts in zip(names, seqs):
        seq = ''.join(parts)
        if not seq:
            raise ValueError(f'Contig {name!r} has no sequence')
        records.append((name, seq))
    logger.info(
        'Parsed %d FASTA record(s), %d residues total (digest %s)',
        len(records),
        sum(len(s) for _, s in records),
        digest,
    )
    return FastaInput(digest=digest, records=tuple(records))
