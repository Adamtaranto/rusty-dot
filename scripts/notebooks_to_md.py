"""Render the tutorial notebooks to Markdown for the documentation build.

Zensical has no notebook support (the MkDocs build used ``mkdocs-jupyter``),
so ``docs/tutorials/*.ipynb`` are converted to sibling ``.md`` files before
``zensical build`` runs.  The notebooks stay the source of truth; the
generated Markdown and extracted images are build products and are
gitignored.

Notebooks are **not** executed — as under ``mkdocs-jupyter``'s
``execute: false``, only outputs already stored in the ``.ipynb`` are
rendered.

Run from the repository root::

    python scripts/notebooks_to_md.py
"""

from __future__ import annotations

from pathlib import Path
import sys

TUTORIALS = Path('docs/tutorials')


def convert(notebook: Path) -> bool:
    """Convert one notebook to Markdown beside itself.

    Skips the conversion when the existing ``.md`` is newer than the
    notebook, so repeated local builds are cheap.

    Parameters
    ----------
    notebook : pathlib.Path
        Path to a ``.ipynb`` file.

    Returns
    -------
    bool
        ``True`` if the Markdown was (re)written, ``False`` if the existing
        output was already up to date.
    """
    from nbconvert import MarkdownExporter
    from nbconvert.writers import FilesWriter
    import nbformat

    target = notebook.with_suffix('.md')
    if target.exists() and target.stat().st_mtime >= notebook.stat().st_mtime:
        return False

    nb = nbformat.read(notebook, as_version=4)
    # Extract images into a per-notebook directory; the default puts them
    # loose in docs/tutorials/, where names like output_6_1.png collide
    # between notebooks.
    body, resources = MarkdownExporter().from_notebook_node(
        nb, resources={'output_files_dir': f'{notebook.stem}_files'}
    )
    if not body.lstrip().startswith('#'):
        # Nav titles and the table of contents need a top-level heading.
        body = f'# {notebook.stem.replace("_", " ").title()}\n\n{body}'

    # Images are written to <name>_files/ next to the Markdown, which is
    # where the exporter's relative links already point.
    resources['output_extension'] = '.md'
    writer = FilesWriter(build_directory=str(notebook.parent))
    writer.write(body, resources, notebook_name=notebook.stem)
    return True


def main() -> int:
    """Convert every tutorial notebook, reporting what changed.

    Returns
    -------
    int
        Process exit status: ``0`` on success, ``1`` if the tutorials
        directory is missing.
    """
    if not TUTORIALS.is_dir():
        print(f'not found: {TUTORIALS} (run from the repository root)', file=sys.stderr)
        return 1
    notebooks = sorted(TUTORIALS.glob('*.ipynb'))
    for notebook in notebooks:
        action = 'converted' if convert(notebook) else 'up to date'
        print(f'{action}: {notebook}')
    print(f'{len(notebooks)} notebook(s) processed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
