"""Selection of the bundled rusty-dot wasm wheel at app startup.

The app ships one (ideally) wasm wheel in ``app/wheels/``.  CI caching has
previously leaked wheels built for other Emscripten versions into the
bundle, and micropip hard-rejects a wheel whose platform tag does not match
the running Pyodide — so the wheel is chosen by matching the interpreter's
platform instead of taking whatever sorts last.
"""

from __future__ import annotations

from pathlib import Path
import sysconfig


def runtime_platform_tag() -> str:
    """Return the wheel platform tag for the running interpreter.

    Returns
    -------
    str
        Platform tag in wheel-filename form (e.g.
        ``'emscripten_3_1_58_wasm32'``), derived from
        ``sysconfig.get_platform()`` (``'emscripten-3.1.58-wasm32'`` under
        Pyodide 0.27).
    """
    return sysconfig.get_platform().replace('-', '_').replace('.', '_')


def pick_wasm_wheel(wheels: list[Path], platform_tag: str) -> Path:
    """Pick the bundled wheel matching the running Pyodide's platform.

    Parameters
    ----------
    wheels : list[pathlib.Path]
        Candidate ``.whl`` paths (contents of ``app/wheels/``).
    platform_tag : str
        Wheel platform tag of the running interpreter, as returned by
        :func:`runtime_platform_tag`.

    Returns
    -------
    pathlib.Path
        The single wheel whose filename carries *platform_tag*.  If several
        match, the lexically last (typically highest version) wins.

    Raises
    ------
    RuntimeError
        If *wheels* is empty or none of them matches *platform_tag* (a
        mis-built bundle; installing a mismatched wheel would make micropip
        fail with a much less helpful error).
    """
    if not wheels:
        raise RuntimeError(
            'No rusty-dot wasm wheel bundled with the app '
            '(expected app/wheels/*.whl at export time).'
        )
    matching = sorted(w for w in wheels if platform_tag in w.name)
    if not matching:
        names = ', '.join(w.name for w in wheels)
        raise RuntimeError(
            f'No bundled wheel matches this Pyodide ({platform_tag}); '
            f'found: {names}. The site was built against a different '
            'Pyodide/Emscripten version.'
        )
    return matching[-1]
