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
 *   'rd-match-select' {row, col, layer, idx, qs, qe, ts, te, strand} — a
 *       clicked match wants its sequences (aligned when a CIGAR exists,
 *       plain slices otherwise); forwarded to Shiny as the
 *       'match_select' input.
 *
 * parent -> child
 *   'rd-display-opts' {dot_size, min_length} — display-only options applied
 *       client-side inside the report (CSS stroke width + per-segment
 *       hiding), so changing them never re-renders matplotlib in Pyodide.
 *       Values arrive from the server via the 'rd_display_opts' custom
 *       message. The child has an opaque origin (sandbox without
 *       allow-same-origin), so targetOrigin must be '*'.
 *   'rd-seq-response' {row, col, layer, idx, text?, error?} — the aligned
 *       sequences (or an error code) for a prior 'rd-match-select'; arrives
 *       from the server via the 'rd_match_seq' custom message.
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
    if (msg.type === 'rd-match-select') {
      var frame2 = reportFrame();
      if (!frame2 || ev.source !== frame2.contentWindow) {
        return;
      }
      var mrow = Number(msg.row);
      var mcol = Number(msg.col);
      var midx = Number(msg.idx);
      if (
        !Number.isInteger(mrow) || mrow < 0 ||
        !Number.isInteger(mcol) || mcol < 0 ||
        !Number.isInteger(midx) || midx < 0 ||
        typeof msg.layer !== 'string'
      ) {
        return;
      }
      var payload = { row: mrow, col: mcol, layer: msg.layer, idx: midx };
      ['qs', 'qe', 'ts', 'te'].forEach(function (k) {
        var v = Number(msg[k]);
        if (Number.isFinite(v) && v >= 0) {
          payload[k] = v;
        }
      });
      if (msg.strand === '+' || msg.strand === '-') {
        payload.strand = msg.strand;
      }
      if (window.Shiny && typeof window.Shiny.setInputValue === 'function') {
        window.Shiny.setInputValue('match_select', payload, {
          priority: 'event',
        });
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
      window.Shiny.addCustomMessageHandler('rd_match_seq', function (msg) {
        var frame = reportFrame();
        if (frame && frame.contentWindow && msg) {
          frame.contentWindow.postMessage(
            Object.assign({ type: 'rd-seq-response' }, msg),
            '*'
          );
        }
      });
    } else {
      setTimeout(register, 100);
    }
  }

  register();
})();
