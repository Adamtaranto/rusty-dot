"""Alignment-method registry and result construction.

Two method families exist:

* **In-Pyodide** — k-mer matching via the rusty-dot Rust core, and PAF
  import (parsed directly from uploaded text).  Implemented here.
* **biowasm tools** — minimap2 / LASTZ / nucmer run in an Aioli WebWorker
  on the JS side (``app/www/aligners.js``); their text output is handed
  back to Python and converted to ``PafRecord`` objects by parsers in this
  module.  Python builds the CLI argument list (:func:`build_tool_args`)
  and the plain-FASTA payload (:func:`fasta_text`); JS only mounts files
  and runs the tool.  The mounted file names are fixed:
  :data:`TARGET_FILENAME` and :data:`QUERY_FILENAME`.

NCBI BLAST is deliberately absent: no production WASM build of BLAST+
exists, so it cannot run client-side.  ``nucmer`` (MUMmer4) is offered as
the closest substitute.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Iterable

if TYPE_CHECKING:  # pragma: no cover - only for type checkers
    from rusty_dot.paf_io import PafAlignment, PafRecord

#: Method key -> UI label for every planned alignment method.
METHOD_LABELS: dict[str, str] = {
    'kmer': 'k-mer matching (rusty-dot)',
    'paf': 'Import PAF file',
    'minimap2': 'minimap2 (biowasm)',
    'lastz': 'LASTZ (biowasm)',
    'nucmer': 'nucmer / MUMmer4 (biowasm)',
}

#: Methods implemented so far; the UI greys out the rest.
AVAILABLE_METHODS: frozenset[str] = frozenset(
    {'kmer', 'paf', 'minimap2', 'lastz', 'nucmer'}
)

#: Methods that run in a biowasm (Aioli) WebWorker on the JS side.
BIOWASM_TOOLS: frozenset[str] = frozenset({'minimap2', 'lastz', 'nucmer'})

#: File names under which aligners.js mounts the two assemblies.
TARGET_FILENAME = 'target.fa'
QUERY_FILENAME = 'query.fa'

#: minimap2 presets exposed in the UI (assembly-to-assembly modes).
MINIMAP2_PRESETS: tuple[str, ...] = ('asm5', 'asm10', 'asm20')

#: Column order requested from LASTZ ``--format=general:`` output.  The
#: parser (:func:`lastz_general_to_records`) indexes columns by this order,
#: so the two must stay in sync.
LASTZ_GENERAL_FIELDS: tuple[str, ...] = (
    'name1',  # target name
    'zstart1',  # target start (0-based)
    'end1',  # target end (exclusive)
    'length1',  # target length
    'name2',  # query name
    'strand2',  # query strand (+/-)
    'zstart2+',  # query start on the forward strand (0-based)
    'end2+',  # query end on the forward strand (exclusive)
    'length2',  # query length
    'nmatch',  # matched bases
)


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


def fasta_text(records: Iterable[tuple[str, str]]) -> str:
    """Reconstruct plain FASTA text from parsed ``(name, sequence)`` records.

    Used to hand assemblies to the biowasm tools: Python owns gzip handling
    and FASTA parsing (``core.fasta``), so JS only ever sees plain text.

    Parameters
    ----------
    records : iterable of tuple[str, str]
        Parsed ``(name, sequence)`` pairs (e.g. ``FastaInput.records``).

    Returns
    -------
    str
        FASTA text with one header line and one (unwrapped) sequence line
        per record.

    Raises
    ------
    ValueError
        If *records* is empty, or a record has an empty name or sequence.
    """
    parts: list[str] = []
    for name, seq in records:
        if not name or not seq:
            raise ValueError(f'FASTA record with empty name or sequence: {name!r}')
        parts.append(f'>{name}\n{seq}\n')
    if not parts:
        raise ValueError('No FASTA records to serialise')
    return ''.join(parts)


def build_tool_args(tool: str, params: dict[str, Any]) -> list[str]:
    """Build the CLI argument list for a biowasm tool run.

    The returned list does *not* include the program name (aligners.js
    prepends it) but does include the mounted input file names
    (:data:`TARGET_FILENAME`, :data:`QUERY_FILENAME`).

    Parameters
    ----------
    tool : str
        One of :data:`BIOWASM_TOOLS`.
    params : dict[str, Any]
        Tool parameters from the UI:

        * ``minimap2`` — ``preset`` (one of :data:`MINIMAP2_PRESETS`),
          optional ``k`` and ``w`` (0 or ``None`` = use the preset default).
        * ``lastz`` — ``step`` (seed step, >= 1), ``notransition`` (bool),
          ``gapped`` (bool).
        * ``nucmer`` — ``l`` (min match length), ``c`` (min cluster
          length), ``maxmatch`` (bool).

    Returns
    -------
    list[str]
        CLI arguments in tool order (options first, then positional
        inputs; minimap2/nucmer take ``target query``, LASTZ takes
        ``target[multiple] query``).

    Raises
    ------
    ValueError
        If *tool* is unknown or a parameter value is invalid.
    """
    if tool == 'minimap2':
        preset = params.get('preset', 'asm20')
        if preset not in MINIMAP2_PRESETS:
            raise ValueError(f'Unknown minimap2 preset: {preset!r}')
        args = ['-x', preset]
        for flag, key in (('-k', 'k'), ('-w', 'w')):
            value = params.get(key) or 0
            if value < 0:
                raise ValueError(f'minimap2 {flag} must be >= 0, got {value}')
            if value:
                args += [flag, str(int(value))]
        # minimap2 writes PAF to stdout by default: minimap2 [opts] target query
        return [*args, TARGET_FILENAME, QUERY_FILENAME]
    if tool == 'lastz':
        raw_step = params.get('step')
        step = 1 if raw_step is None else int(raw_step)
        if step < 1:
            raise ValueError(f'LASTZ step must be >= 1, got {step}')
        args = [
            # [multiple] treats the target file as a multi-sequence database.
            f'{TARGET_FILENAME}[multiple]',
            QUERY_FILENAME,
            '--format=general:' + ','.join(LASTZ_GENERAL_FIELDS),
        ]
        if step != 1:
            args.append(f'--step={step}')
        if params.get('notransition'):
            args.append('--notransition')
        if not params.get('gapped', True):
            args.append('--nogapped')
        return args
    if tool == 'nucmer':
        raw_l = params.get('l')
        raw_c = params.get('c')
        min_match = 20 if raw_l is None else int(raw_l)
        min_cluster = 65 if raw_c is None else int(raw_c)
        if min_match < 1 or min_cluster < 1:
            raise ValueError(
                f'nucmer -l/-c must be >= 1, got l={min_match}, c={min_cluster}'
            )
        args = []
        if params.get('maxmatch'):
            args.append('--maxmatch')
        args += ['-l', str(min_match), '-c', str(min_cluster)]
        # nucmer writes <prefix>.delta (default prefix "out"); aligners.js
        # reads out.delta back after the run.
        return [*args, TARGET_FILENAME, QUERY_FILENAME]
    raise ValueError(f'Unknown biowasm tool: {tool!r}')


def lastz_general_to_records(text: str) -> list[PafRecord]:
    """Convert LASTZ ``--format=general:`` output to PAF records.

    Expects the column order in :data:`LASTZ_GENERAL_FIELDS` (the format
    requested by :func:`build_tool_args`).  LASTZ reports query
    coordinates on the forward strand for the ``zstart2+``/``end2+``
    columns, matching the PAF convention for ``-`` strand records.

    Parameters
    ----------
    text : str
        LASTZ general-format output (header lines start with ``#``).

    Returns
    -------
    list[rusty_dot.paf_io.PafRecord]
        One record per alignment line.

    Raises
    ------
    ValueError
        If a line has the wrong number of columns, non-numeric
        coordinates, an invalid strand, or no alignments are found.
    """
    from rusty_dot.paf_io import PafRecord

    n_cols = len(LASTZ_GENERAL_FIELDS)
    records: list[PafRecord] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) != n_cols:
            raise ValueError(
                f'LASTZ output line {lineno} has {len(fields)} columns; '
                f'expected {n_cols}: {line!r}'
            )
        try:
            t_start, t_end, t_len = (int(fields[i]) for i in (1, 2, 3))
            q_start, q_end, q_len = (int(fields[i]) for i in (6, 7, 8))
            n_match = int(fields[9])
        except ValueError as exc:
            raise ValueError(f'Invalid LASTZ output at line {lineno}: {exc}') from exc
        strand = fields[5]
        if strand not in ('+', '-'):
            raise ValueError(f'Invalid strand {strand!r} at LASTZ output line {lineno}')
        records.append(
            PafRecord(
                query_name=fields[4],
                query_len=q_len,
                query_start=q_start,
                query_end=q_end,
                strand=strand,
                target_name=fields[0],
                target_len=t_len,
                target_start=t_start,
                target_end=t_end,
                residue_matches=n_match,
                alignment_block_len=max(t_end - t_start, q_end - q_start),
                mapping_quality=255,
            )
        )
    if not records:
        raise ValueError('LASTZ produced no alignments')
    return records


def nucmer_delta_to_records(text: str) -> list[PafRecord]:
    """Convert a nucmer ``.delta`` file to PAF records.

    The delta format holds a file-path header line and a program line,
    then per-sequence-pair sections starting ``>ref query rlen qlen``.
    Each alignment in a section is a coordinate line ``rstart rend qstart
    qend err sim stp`` (1-based inclusive; reverse strand = start > end)
    followed by indel offsets terminated by ``0``.  Only the coordinate
    line is used; indel offsets are skipped.

    Parameters
    ----------
    text : str
        Content of a nucmer ``.delta`` output file.

    Returns
    -------
    list[rusty_dot.paf_io.PafRecord]
        One record per alignment, with 0-based half-open forward-strand
        coordinates as required by the PAF spec.

    Raises
    ------
    ValueError
        If a header or coordinate line is malformed, or no alignments are
        found.
    """
    from rusty_dot.paf_io import PafRecord

    records: list[PafRecord] = []
    ref = qry = ''
    rlen = qlen = 0
    in_section = False
    skipping_indels = False
    for lineno, line in enumerate(text.splitlines(), start=1):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            fields = line[1:].split()
            if len(fields) != 4:
                raise ValueError(
                    f'Invalid delta sequence header at line {lineno}: {line!r}'
                )
            ref, qry = fields[0], fields[1]
            try:
                rlen, qlen = int(fields[2]), int(fields[3])
            except ValueError as exc:
                raise ValueError(
                    f'Invalid delta sequence header at line {lineno}: {exc}'
                ) from exc
            in_section = True
            skipping_indels = False
            continue
        if not in_section:
            # File-path header line and the NUCMER/PROMER program line.
            continue
        if skipping_indels:
            if line == '0':
                skipping_indels = False
            continue
        fields = line.split()
        if len(fields) != 7:
            raise ValueError(
                f'Invalid delta coordinate line at line {lineno}: {line!r}'
            )
        try:
            r1, r2, q1, q2, errs = (int(fields[i]) for i in range(5))
        except ValueError as exc:
            raise ValueError(
                f'Invalid delta coordinate line at line {lineno}: {exc}'
            ) from exc
        # 1-based inclusive -> 0-based half-open; start > end marks the
        # reverse strand on that axis.  PAF wants both axes forward with a
        # strand flag, so normalise each axis and XOR the orientations.
        strand = '+'
        if r1 > r2:
            r1, r2 = r2, r1
            strand = '-' if strand == '+' else '+'
        if q1 > q2:
            q1, q2 = q2, q1
            strand = '-' if strand == '+' else '+'
        t_start, t_end = r1 - 1, r2
        q_start, q_end = q1 - 1, q2
        block = max(t_end - t_start, q_end - q_start)
        records.append(
            PafRecord(
                query_name=qry,
                query_len=qlen,
                query_start=q_start,
                query_end=q_end,
                strand=strand,
                target_name=ref,
                target_len=rlen,
                target_start=t_start,
                target_end=t_end,
                residue_matches=max(block - errs, 0),
                alignment_block_len=block,
                mapping_quality=255,
            )
        )
        skipping_indels = True
    if not records:
        raise ValueError('nucmer produced no alignments')
    return records


def alignment_from_tool_output(tool: str, text: str) -> PafAlignment:
    """Parse a biowasm tool's text output into a ``PafAlignment``.

    Parameters
    ----------
    tool : str
        One of :data:`BIOWASM_TOOLS`.
    text : str
        The tool's output: PAF text (minimap2), general-format text
        (LASTZ), or delta-file content (nucmer).

    Returns
    -------
    rusty_dot.paf_io.PafAlignment
        Parsed alignment.

    Raises
    ------
    ValueError
        If *tool* is unknown or the output cannot be parsed / holds no
        alignments.
    """
    from rusty_dot.paf_io import PafAlignment

    if tool == 'minimap2':
        return paf_alignment_from_text(text)
    if tool == 'lastz':
        return PafAlignment.from_records(lastz_general_to_records(text))
    if tool == 'nucmer':
        return PafAlignment.from_records(nucmer_delta_to_records(text))
    raise ValueError(f'Unknown biowasm tool: {tool!r}')
