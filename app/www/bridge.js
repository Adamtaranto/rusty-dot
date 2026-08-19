/* rusty-dot app <-> embedded report bridge (W2).
 *
 * The interactive dotplot report is embedded in a sandboxed srcdoc iframe.
 * Two message flows cross that boundary:
 *
 * child -> parent
 *   'rd-panel-dblclick' {row, col}  — forwarded to Shiny as the
 *       'panel_dblclick' input so the app swaps to the single-pair view.
 *   'rd-report-ready'               — the report finished wiring; reply with
 *       the current display options (the iframe element is replaced on every
 *       re-render, so state must be re-sent each time).
 *
 * parent -> child
 *   'rd-display-opts' {dot_size, min_length} — display-only options applied
 *       client-side inside the report (CSS stroke width + per-segment
 *       hiding), so changing them never re-renders matplotlib in Pyodide.
 *       Values arrive from the server via the 'rd_display_opts' custom
 *       message. The child has an opaque origin (sandbox without
 *       allow-same-origin), so targetOrigin must be '*'.
 */
(function () {
  'use strict';

  var lastOpts = { dot_size: 0.5, min_length: 0 };

  function reportFrame() {
    return document.querySelector('iframe.rd-report-frame');
  }

  function pushOptsToReport() {
    var frame = reportFrame();
    if (frame && frame.contentWindow) {
      frame.contentWindow.postMessage(
        {
          type: 'rd-display-opts',
          dot_size: lastOpts.dot_size,
          min_length: lastOpts.min_length,
        },
        '*'
      );
    }
  }

  window.addEventListener('message', function (ev) {
    var msg = ev && ev.data;
    if (!msg) {
      return;
    }
    if (msg.type === 'rd-report-ready') {
      // Only trust the ping if it came from the current report iframe.
      var frame = reportFrame();
      if (frame && ev.source === frame.contentWindow) {
        pushOptsToReport();
      }
      return;
    }
    if (msg.type !== 'rd-panel-dblclick') {
      return;
    }
    var row = Number(msg.row);
    var col = Number(msg.col);
    if (!Number.isInteger(row) || !Number.isInteger(col) || row < 0 || col < 0) {
      return;
    }
    if (window.Shiny && typeof window.Shiny.setInputValue === 'function') {
      window.Shiny.setInputValue(
        'panel_dblclick',
        { row: row, col: col },
        { priority: 'event' }
      );
    }
  });

  function register() {
    if (window.Shiny && window.Shiny.addCustomMessageHandler) {
      window.Shiny.addCustomMessageHandler('rd_display_opts', function (opts) {
        if (!opts) {
          return;
        }
        if (typeof opts.dot_size === 'number' && opts.dot_size > 0) {
          lastOpts.dot_size = opts.dot_size;
        }
        if (typeof opts.min_length === 'number' && opts.min_length >= 0) {
          lastOpts.min_length = opts.min_length;
        }
        pushOptsToReport();
      });
    } else {
      setTimeout(register, 100);
    }
  }

  register();
})();
