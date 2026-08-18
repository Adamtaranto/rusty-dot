"""Patch a shinylive-exported index.html with a loading splash screen.

The bare shinylive export shows a blank page for the many seconds it takes
to download and boot the Pyodide runtime, which reads as a broken site.
This script injects a self-contained overlay (spinner + staged status
messages + app title) that removes itself once the Shiny app inside the
viewer iframe has actually rendered.

Usage
-----
    python scripts/add_loading_splash.py site/app/index.html

Run after ``shinylive export``; it is idempotent (a marker comment guards
against double insertion) and exits non-zero if the file does not look like
a shinylive export.
"""

from __future__ import annotations

from pathlib import Path
import sys

MARKER = '<!-- rd-loading-splash -->'

SPLASH = f"""{MARKER}
<style>
  #rd-splash {{
    position: fixed; inset: 0; z-index: 9999;
    display: flex; flex-direction: column; align-items: center;
    justify-content: center; gap: 1rem;
    background: #fff; color: #1f3a5f;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto,
      Helvetica, Arial, sans-serif;
    transition: opacity 0.4s ease;
  }}
  #rd-splash.rd-hide {{ opacity: 0; pointer-events: none; }}
  #rd-splash h1 {{ font-size: 1.4rem; font-weight: 600; margin: 0; }}
  #rd-splash .rd-sub {{ font-size: 0.95rem; color: #48627f; }}
  #rd-splash .rd-spinner {{
    width: 42px; height: 42px; border-radius: 50%;
    border: 4px solid #d6e2f0; border-top-color: #2c6bb3;
    animation: rd-spin 0.9s linear infinite;
  }}
  @keyframes rd-spin {{ to {{ transform: rotate(360deg); }} }}
  @media (prefers-color-scheme: dark) {{
    #rd-splash {{ background: #14181d; color: #cfe1f5; }}
    #rd-splash .rd-sub {{ color: #8fa8c4; }}
    #rd-splash .rd-spinner {{ border-color: #2a3947; border-top-color: #5b9bd5; }}
  }}
</style>
<div id="rd-splash">
  <div class="rd-spinner"></div>
  <h1>rusty-dot &middot; assembly comparison</h1>
  <div class="rd-sub" id="rd-splash-msg">Starting&hellip;</div>
  <div class="rd-sub" style="font-size: 0.8rem">
    Everything runs in your browser &mdash; your files never leave your machine.
  </div>
</div>
<script>
  (function () {{
    'use strict';
    var msgs = [
      [0, 'Downloading the Python runtime (\\u224830 MB on first visit, cached afterwards)\\u2026'],
      [8, 'Starting Python in your browser\\u2026'],
      [16, 'Installing rusty-dot and matplotlib\\u2026'],
      [30, 'Almost there \\u2014 first visits take the longest\\u2026'],
      [60, 'Still working. Slow connection? The runtime is \\u224830 MB\\u2026'],
    ];
    var start = Date.now();
    var msgEl = document.getElementById('rd-splash-msg');
    var timer = setInterval(function () {{
      var t = (Date.now() - start) / 1000;
      for (var i = msgs.length - 1; i >= 0; i--) {{
        if (t >= msgs[i][0]) {{ msgEl.textContent = msgs[i][1]; break; }}
      }}
    }}, 1000);
    function appReady() {{
      // The shinylive viewer hosts the app in a same-origin iframe; the app
      // is usable once its own widgets exist in that frame.
      var frames = document.querySelectorAll('iframe');
      for (var i = 0; i < frames.length; i++) {{
        try {{
          var doc = frames[i].contentDocument;
          if (doc && doc.getElementById('run')) return true;
        }} catch (e) {{ /* cross-origin frame: not ours */ }}
      }}
      return false;
    }}
    var poll = setInterval(function () {{
      // Give up (and get out of the way) after 5 minutes regardless.
      if (appReady() || Date.now() - start > 300000) {{
        clearInterval(poll);
        clearInterval(timer);
        var el = document.getElementById('rd-splash');
        el.classList.add('rd-hide');
        setTimeout(function () {{ el.remove(); }}, 500);
      }}
    }}, 400);
  }})();
</script>
"""


def patch(index_html: Path) -> bool:
    """Insert the splash into *index_html* (idempotent).

    Parameters
    ----------
    index_html : pathlib.Path
        Path to a shinylive-exported ``index.html``.

    Returns
    -------
    bool
        ``True`` if the file was modified, ``False`` if the splash was
        already present.

    Raises
    ------
    ValueError
        If the file does not look like a shinylive export (no ``</body>``
        or no shinylive script reference).
    """
    text = index_html.read_text()
    if MARKER in text:
        return False
    if '</body>' not in text or 'shinylive' not in text:
        raise ValueError(f'{index_html} does not look like a shinylive export')
    text = text.replace('<title>Shiny App</title>', '<title>rusty-dot</title>')
    index_html.write_text(text.replace('</body>', SPLASH + '</body>'))
    return True


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('usage: add_loading_splash.py <exported-index.html>')
    target = Path(sys.argv[1])
    changed = patch(target)
    print(f'{"patched" if changed else "already patched"}: {target}')
