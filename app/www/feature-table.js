// Delegated event handling for the drill-down Annotations table.
//
// The table is plain HTML with raw <input> controls rather than one Shiny
// input per feature: the real GFFs this app targets run to ~600 features
// per sequence, and ~1200 bound inputs (each with its own binding and
// message round-trip) visibly freezes Pyodide.  Instead a single listener
// on the container posts *deltas*:
//
//   Shiny.setInputValue('feature_table_change',
//       {kind: 'vis'|'color', uid, value}, {priority: 'event'})
//   Shiny.setInputValue('feature_table_change',
//       {kind: 'bulk'|'color', uids: [...], value}, {priority: 'event'})
//
// Filtering is purely client-side and never reaches the server.
//
// This mirrors the approach already used for the display-options bridge:
// heavy state stays in the browser, only what changed crosses over.

(function () {
  'use strict';

  var TABLE_ID = 'rd-feature-table';

  function send(payload) {
    if (!window.Shiny || !window.Shiny.setInputValue) return;
    window.Shiny.setInputValue('feature_table_change', payload, {
      priority: 'event',
    });
  }

  function table() {
    return document.getElementById(TABLE_ID);
  }

  function visibleRows() {
    var t = table();
    if (!t) return [];
    return Array.prototype.filter.call(t.tBodies[0].rows, function (row) {
      return row.style.display !== 'none';
    });
  }

  function uidsOf(rows) {
    return rows.map(function (r) {
      return r.getAttribute('data-uid');
    });
  }

  // --- per-row toggles and colour pickers ---------------------------------
  document.addEventListener('change', function (ev) {
    var el = ev.target;
    if (!el || !el.getAttribute) return;
    var uid = el.getAttribute('data-uid');
    var kind = el.getAttribute('data-kind');
    if (!uid || !kind) return;
    if (!el.closest || !el.closest('#' + TABLE_ID)) return;
    if (kind === 'vis') {
      send({ kind: 'vis', uid: uid, value: el.checked });
    } else if (kind === 'color') {
      send({ kind: 'color', uid: uid, value: el.value });
    }
  });

  // --- bulk buttons -------------------------------------------------------
  // Bulk actions apply to the rows currently *visible* under the filter, so
  // "Hide all" after filtering to one type does what it looks like it does.
  document.addEventListener('click', function (ev) {
    var btn = ev.target;
    if (!btn || !btn.getAttribute) return;
    var action = btn.getAttribute('data-bulk');
    if (!action) return;
    var rows = visibleRows();
    if (!rows.length) return;
    var uids = uidsOf(rows);
    if (action === 'show' || action === 'hide') {
      var show = action === 'show';
      rows.forEach(function (r) {
        var box = r.querySelector('input[data-kind="vis"]');
        if (box) box.checked = show;
      });
      send({ kind: 'bulk', uids: uids, value: show });
    } else if (action === 'reset') {
      // Empty value means "drop the override and fall back to the type
      // colour" -- otherwise there is no way back from a per-feature pick.
      // Repaint the swatches here too: edits are held until Apply, so the
      // table does not re-render on this and the server cannot put them
      // back.  Each input carries its type colour for exactly this.
      rows.forEach(function (r) {
        var picker = r.querySelector('input[data-kind="color"]');
        if (picker && picker.dataset.typeColor) {
          picker.value = picker.dataset.typeColor;
        }
      });
      send({ kind: 'color', uids: uids, value: '' });
    }
  });

  // --- client-side filter -------------------------------------------------
  document.addEventListener('input', function (ev) {
    var box = ev.target;
    if (!box || box.id !== 'rd-ft-filter') return;
    var needle = box.value.trim().toLowerCase();
    var t = table();
    if (!t) return;
    Array.prototype.forEach.call(t.tBodies[0].rows, function (row) {
      var hay = (
        row.getAttribute('data-type') +
        ' ' +
        row.getAttribute('data-source') +
        ' ' +
        row.textContent
      ).toLowerCase();
      row.style.display = !needle || hay.indexOf(needle) !== -1 ? '' : 'none';
    });
  });
})();
