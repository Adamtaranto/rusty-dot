/* rusty-dot app <-> embedded report bridge (W2).
 *
 * The interactive dotplot report is embedded in a sandboxed srcdoc iframe.
 * A small script injected into the report posts a message to the parent
 * window when a panel is double-clicked; this listener forwards it to the
 * Shiny server as the 'panel_dblclick' input so the app can swap to a
 * standalone single-pair view.
 */
(function () {
  'use strict';
  window.addEventListener('message', function (ev) {
    var msg = ev && ev.data;
    if (!msg || msg.type !== 'rd-panel-dblclick') {
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
})();
