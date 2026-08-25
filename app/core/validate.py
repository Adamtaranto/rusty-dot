"""Cross-checks between an uploaded assembly and a pre-computed PAF.

When a user pairs a PAF alignment with an optional query assembly (for the
reordered-FASTA export), the contig names must line up with the PAF's
*query* column — a mismatch silently produces empty or wrong exports.
These checks turn the common mistakes (wrong file, swapped query/target,
partial assemblies) into explicit warnings.
"""

from __future__ import annotations

from collections.abc import Iterable


def validate_query_names(
    assembly_names: Iterable[str],
    paf_query_names: Iterable[str],
    paf_target_names: Iterable[str],
) -> list[str]:
    """Check uploaded query-assembly contig names against a PAF's name sets.

    Parameters
    ----------
    assembly_names : Iterable[str]
        Contig names from the uploaded query assembly FASTA.
    paf_query_names : Iterable[str]
        Names from the PAF query-name column (column 1).
    paf_target_names : Iterable[str]
        Names from the PAF target-name column (column 6).

    Returns
    -------
    list[str]
        Human-readable warnings, empty when everything lines up.  Covers:
        no assembly name in the PAF query column (with a hint when they
        match the *target* column instead — swapped inputs), assembly
        contigs absent from the PAF, PAF query names absent from the
        assembly, and names that appear in both the PAF query and target
        columns (ambiguous role).
    """
    assembly = set(assembly_names)
    paf_q = set(paf_query_names)
    paf_t = set(paf_target_names)
    warnings: list[str] = []

    def _preview(names: set[str], limit: int = 5) -> str:
        shown = ', '.join(sorted(names)[:limit])
        return shown + (', …' if len(names) > limit else '')

    matched = assembly & paf_q
    if assembly and not matched:
        msg = (
            'None of the uploaded assembly contig names appear in the PAF '
            'query column — the reordered-FASTA export will be empty.'
        )
        if assembly & paf_t:
            msg += (
                ' They DO match the PAF target column: the assembly looks '
                'like the target/reference, or the PAF was generated with '
                'query and target swapped.'
            )
        warnings.append(msg)
        # The finer-grained checks below would only repeat this.
        return warnings

    missing_from_paf = assembly - paf_q
    if missing_from_paf:
        warnings.append(
            f'{len(missing_from_paf)} assembly contig(s) have no alignments '
            f'in the PAF ({_preview(missing_from_paf)}); they will appear '
            'unplaced in the reordered FASTA.'
        )

    missing_from_assembly = paf_q - assembly
    if missing_from_assembly:
        warnings.append(
            f'{len(missing_from_assembly)} PAF query name(s) are not in the '
            f'uploaded assembly ({_preview(missing_from_assembly)}); their '
            'sequences cannot be exported.'
        )

    ambiguous = paf_q & paf_t
    if ambiguous:
        warnings.append(
            f'{len(ambiguous)} name(s) appear in BOTH the PAF query and '
            f'target columns ({_preview(ambiguous)}); their role is '
            'ambiguous and reordering/orientation may be wrong.'
        )
    return warnings


def validate_annotation_names(
    assembly_names: Iterable[str],
    annotation_names: Iterable[str],
    role: str,
) -> list[str]:
    """Check a GFF's sequence names against the assembly it annotates.

    A GFF whose ``seqname`` column does not match the FASTA headers is the
    single most common annotation mistake — a different assembly version,
    or names rewritten by an intermediate tool.  Nothing errors: the
    features simply never appear, which is easy to misread as "the plot is
    broken".

    Parameters
    ----------
    assembly_names : Iterable[str]
        Contig names from the assembly FASTA for this role.
    annotation_names : Iterable[str]
        Distinct ``seqname`` values from the parsed GFF.
    role : str
        ``'query'`` or ``'target'``, used in the message.

    Returns
    -------
    list[str]
        Human-readable warnings, empty when every GFF sequence is present.
    """
    assembly = set(assembly_names)
    annotated = set(annotation_names)
    if not assembly or not annotated:
        return []

    def _preview(names: set[str], limit: int = 5) -> str:
        shown = ', '.join(sorted(names)[:limit])
        return shown + (', …' if len(names) > limit else '')

    missing = annotated - assembly
    if not missing:
        return []
    if missing == annotated:
        # Nothing lines up at all: almost always the wrong file or the
        # wrong role, so say so rather than listing every name.
        return [
            f'None of the {role} GFF sequence names match the {role} '
            f'assembly ({_preview(annotated)}); no features will be drawn. '
            'Check the GFF belongs to this assembly and is assigned to the '
            'right role.'
        ]
    return [
        f'{len(missing)} sequence name(s) in the {role} GFF are not in the '
        f'{role} assembly ({_preview(missing)}); features on them will not '
        'be drawn.'
    ]
