// The drill-down Annotations table, built and driven client-side.
//
// The server ships a static shell (header, tools, empty tbody) and posts
// the row data through the 'rd_feature_rows' custom message; this file
// builds the body.  Rendering ~600 rows per sequence as raw HTML through
// a @render.ui made the Plot <-> Annotations tab switch stall, since
// Shiny re-inserted and re-scanned the whole blob on every re-show.
//
// The rows carry raw <input> controls rather than one Shiny input per
// feature: ~1200 bound inputs (each with its own binding and message
// round-trip) visibly freezes Pyodide.  Instead a single listener on the
// container posts *deltas*:
//
//   Shiny.setInputValue('feature_table_change',
//       {kind: 'vis'|'color', uid, value}, {priority: 'event'})
//   Shiny.setInputValue('feature_table_change',
//       {kind: 'bulk'|'color', uids: [...], value}, {priority: 'event'})
//
// Filtering, sorting and row selection are purely client-side and never
// reach the server; selected rows are announced to bridge.js through a
// 'rd-ft-selection' CustomEvent so the report can band them.

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

  function esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }

  function fmt(n) {
    // Same thousands separators the server used to format with; the
    // numeric sort strips them back out.
    return Number(n).toLocaleString('en-US');
  }

  // --- building the table body from the pushed payload --------------------

  var lastPayload = null;

  function rowHtml(r, idx) {
    return (
      '<tr data-uid="' + esc(r.uid) + '" data-idx="' + idx + '"' +
      ' data-type="' + esc(r.type) + '" data-source="' + esc(r.source) + '"' +
      ' data-role="' + esc(r.role) + '" data-contig="' + esc(r.contig) + '"' +
      // 0-based start for the report payload; the cell shows 1-based.
      ' data-start="' + (Number(r.start) - 1) + '"' +
      ' data-end="' + Number(r.end) + '"' +
      ' data-strand="' + esc(r.strand) + '">' +
      '<td><input type="checkbox" data-uid="' + esc(r.uid) +
      '" data-kind="vis"' + (r.visible ? ' checked' : '') + '></td>' +
      '<td><input type="color" data-uid="' + esc(r.uid) +
      '" data-kind="color" value="' + esc(r.color) +
      '" data-type-color="' + esc(r.type_color) +
      '" class="rd-color-input"></td>' +
      '<td class="rd-ft-role">' + esc(r.role) + '</td>' +
      '<td class="rd-ft-contig">' + esc(r.contig) + '</td>' +
      '<td class="rd-ft-type">' + esc(r.type) + '</td>' +
      '<td class="rd-ft-name">' + esc(r.name) + '</td>' +
      '<td class="rd-ft-num">' + fmt(r.start) + '</td>' +
      '<td class="rd-ft-num">' + fmt(r.end) + '</td>' +
      '<td class="rd-ft-num">' + fmt(r.length) + '</td>' +
      '<td>' + esc(r.strand) + '</td>' +
      '<td class="rd-ft-src">' + esc(r.source) + '</td>' +
      '<td class="rd-ft-attrs" title="' + esc(r.attrs) + '">' +
      esc(r.attrs) + '</td>' +
      '</tr>'
    );
  }

  function fillColumnSelect(t) {
    var sel = document.getElementById('rd-ft-filter-col');
    if (!sel) return;
    var opts = '<option value="-1">Any column</option>';
    Array.prototype.forEach.call(t.tHead.rows[0].cells, function (th, i) {
      var kind = th.getAttribute('data-sort');
      // Checkbox and swatch columns hold no text to filter on.
      if (kind === 'check' || kind === 'color') return;
      opts += '<option value="' + i + '">' + esc(th.textContent) + '</option>';
    });
    sel.innerHTML = opts;
  }

  function resetSortIndicators(t) {
    sortCol = null;
    sortDir = 0;
    Array.prototype.forEach.call(t.tHead.rows[0].cells, function (cell) {
      cell.setAttribute('aria-sort', 'none');
      cell.classList.remove('rd-sort-asc', 'rd-sort-desc');
    });
  }

  function populate() {
    var t = table();
    if (!t || !lastPayload) return;
    var rows = lastPayload.rows || [];
    t.tBodies[0].innerHTML = rows.map(rowHtml).join('');
    resetSortIndicators(t);
    fillColumnSelect(t);
    // Announce the emptied selection so bridge.js drops any held
    // highlights -- new rows mean the old uids no longer apply.
    clearSelection(true);
    applyFilters();
    var caption = document.querySelector('.rd-ft-caption');
    if (caption && lastPayload.pair) {
      caption.innerHTML =
        '<b>' + rows.length + '</b> feature(s) on <b>' +
        esc(lastPayload.pair[0]) + '</b> and <b>' +
        esc(lastPayload.pair[1]) + '</b>. ' +
        'Coordinates are 1-based inclusive, as in the source file. ' +
        'A feature loaded from two files (e.g. GenBank + GFF) lists once ' +
        'per source.';
    }
  }

  // The shell and the row payload race: the message can arrive before the
  // shell's DOM exists, and the shell can be re-inserted (plot_area
  // re-renders) without a fresh message.  The observer catches both: any
  // insertion that leaves an empty tbody next to a held payload gets
  // filled.
  new MutationObserver(function () {
    var t = table();
    if (t && lastPayload && t.tBodies[0].rows.length === 0) populate();
  }).observe(document.documentElement, { childList: true, subtree: true });

  function register() {
    if (window.Shiny && window.Shiny.addCustomMessageHandler) {
      window.Shiny.addCustomMessageHandler('rd_feature_rows', function (msg) {
        if (!msg || !Array.isArray(msg.rows)) return;
        lastPayload = msg;
        if (table()) populate();
      });
    } else {
      setTimeout(register, 100);
    }
  }

  register();

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

  // --- column sorting -----------------------------------------------------
  // Sorting is done here rather than on the server for the same reason the
  // table posts deltas: re-rendering hundreds of rows is the cost this
  // design exists to avoid, and it would throw away edits the user has not
  // applied yet.  Moving the existing <tr> elements keeps every checkbox,
  // swatch and data attribute attached to its row.

  var sortCol = null; // th currently sorted by
  var sortDir = 0; // 1 asc, -1 desc, 0 file order

  function cellKey(row, index, kind) {
    var cell = row.cells[index];
    if (!cell) return '';
    if (kind === 'check') {
      var box = cell.querySelector('input[type="checkbox"]');
      return box && box.checked ? 1 : 0;
    }
    if (kind === 'color') {
      var picker = cell.querySelector('input[type="color"]');
      return picker ? picker.value.toLowerCase() : '';
    }
    if (kind === 'num') {
      // Cells are formatted with thousands separators.
      var n = parseFloat(cell.textContent.replace(/,/g, ''));
      return isNaN(n) ? -Infinity : n;
    }
    return cell.textContent.trim().toLowerCase();
  }

  function applySort(th) {
    var t = table();
    if (!t) return;
    var index = th.cellIndex;
    var kind = th.getAttribute('data-sort') || 'text';

    if (sortCol === th) {
      // asc -> desc -> back to the order the file gave us.
      sortDir = sortDir === 1 ? -1 : sortDir === -1 ? 0 : 1;
    } else {
      sortCol = th;
      sortDir = 1;
    }

    var rows = Array.prototype.slice.call(t.tBodies[0].rows);
    if (sortDir === 0) {
      sortCol = null;
      rows.sort(function (a, b) {
        return (+a.getAttribute('data-idx')) - (+b.getAttribute('data-idx'));
      });
    } else {
      rows.sort(function (a, b) {
        var ka = cellKey(a, index, kind);
        var kb = cellKey(b, index, kind);
        if (ka < kb) return -sortDir;
        if (ka > kb) return sortDir;
        // Ties keep file order, so repeated sorts are stable and the
        // result never depends on the previous sort.
        return (+a.getAttribute('data-idx')) - (+b.getAttribute('data-idx'));
      });
    }

    var frag = document.createDocumentFragment();
    rows.forEach(function (r) {
      frag.appendChild(r);
    });
    t.tBodies[0].appendChild(frag);

    Array.prototype.forEach.call(t.tHead.rows[0].cells, function (cell) {
      var active = cell === th && sortDir !== 0;
      cell.setAttribute(
        'aria-sort',
        active ? (sortDir === 1 ? 'ascending' : 'descending') : 'none'
      );
      cell.classList.toggle('rd-sort-asc', active && sortDir === 1);
      cell.classList.toggle('rd-sort-desc', active && sortDir === -1);
    });
  }

  document.addEventListener('click', function (ev) {
    var th = ev.target.closest ? ev.target.closest('th[data-sort]') : null;
    if (!th || !th.closest('#' + TABLE_ID)) return;
    applySort(th);
  });

  document.addEventListener('keydown', function (ev) {
    if (ev.key !== 'Enter' && ev.key !== ' ') return;
    var th = ev.target.closest ? ev.target.closest('th[data-sort]') : null;
    if (!th || !th.closest('#' + TABLE_ID)) return;
    ev.preventDefault();
    applySort(th);
  });

  // --- per-column filters -------------------------------------------------
  // Several filters can be active at once, ANDed together; each targets one
  // column (or the whole row).  They take effect on "Apply filters" -- a
  // long table repainting on every keystroke is the lag the button avoids.

  var filters = []; // {colIndex, label, needle}

  function renderChips() {
    var box = document.getElementById('rd-ft-chips');
    if (!box) return;
    box.innerHTML = filters
      .map(function (f, i) {
        return (
          '<span class="rd-ft-chip">' + esc(f.label) + ': “' +
          esc(f.needle) + '”' +
          '<button type="button" data-chip="' + i +
          '" title="Remove this filter" aria-label="Remove filter">' +
          '&times;</button></span>'
        );
      })
      .join('');
  }

  function applyFilters() {
    var t = table();
    if (!t) return;
    Array.prototype.forEach.call(t.tBodies[0].rows, function (row) {
      var shown = filters.every(function (f) {
        var hay;
        if (f.colIndex < 0) {
          hay =
            row.getAttribute('data-type') + ' ' +
            row.getAttribute('data-source') + ' ' +
            row.textContent;
        } else {
          var cell = row.cells[f.colIndex];
          hay = cell ? cell.textContent : '';
        }
        return hay.toLowerCase().indexOf(f.needle) !== -1;
      });
      row.style.display = shown ? '' : 'none';
    });
  }

  function addFilterFromInput() {
    var input = document.getElementById('rd-ft-filter');
    var sel = document.getElementById('rd-ft-filter-col');
    if (!input || !sel) return;
    var needle = input.value.trim().toLowerCase();
    if (!needle) return;
    filters.push({
      colIndex: Number(sel.value),
      label: sel.options[sel.selectedIndex].textContent,
      needle: needle,
    });
    input.value = '';
    renderChips();
  }

  document.addEventListener('click', function (ev) {
    var el = ev.target;
    if (!el || !el.getAttribute) return;
    if (el.id === 'rd-ft-filter-add') {
      addFilterFromInput();
      return;
    }
    if (el.id === 'rd-ft-filter-apply') {
      // Text still sitting in the box counts too, so typing one filter
      // and pressing Apply works without an explicit Add.
      addFilterFromInput();
      applyFilters();
      return;
    }
    var chip = el.getAttribute('data-chip');
    if (chip !== null && el.closest('#rd-ft-chips')) {
      filters.splice(Number(chip), 1);
      renderChips();
      applyFilters(); // removing a filter is itself the action
    }
  });

  document.addEventListener('keydown', function (ev) {
    if (ev.key !== 'Enter') return;
    var el = ev.target;
    if (!el || el.id !== 'rd-ft-filter') return;
    ev.preventDefault();
    addFilterFromInput();
    applyFilters();
  });

  // --- row selection ------------------------------------------------------
  // Click selects a row, ⌘/Ctrl-click toggles it into a multi-selection.
  // bridge.js relays the selection into the report iframe, which bands the
  // matching features on the plot.

  var selectedUids = new Set();

  function announceSelection() {
    var t = table();
    var features = [];
    if (t) {
      Array.prototype.forEach.call(t.tBodies[0].rows, function (row) {
        if (!selectedUids.has(row.getAttribute('data-uid'))) return;
        features.push({
          uid: row.getAttribute('data-uid'),
          role: row.getAttribute('data-role'),
          seqname: row.getAttribute('data-contig'),
          start: Number(row.getAttribute('data-start')),
          end: Number(row.getAttribute('data-end')),
          type: row.getAttribute('data-type'),
          strand: row.getAttribute('data-strand'),
        });
      });
    }
    document.dispatchEvent(
      new CustomEvent('rd-ft-selection', { detail: { features: features } })
    );
  }

  function clearSelection(announce) {
    selectedUids.clear();
    var t = table();
    if (t) {
      Array.prototype.forEach.call(
        t.querySelectorAll('tr.rd-row-selected'),
        function (row) {
          row.classList.remove('rd-row-selected');
        }
      );
    }
    if (announce !== false) announceSelection();
  }

  document.addEventListener('click', function (ev) {
    var el = ev.target;
    if (!el || !el.closest) return;
    // The row's own controls keep their meaning; only bare-cell clicks
    // select.
    if (el.closest('input, button, select, a, label')) return;
    var row = el.closest('#' + TABLE_ID + ' tbody tr');
    if (!row) return;
    var uid = row.getAttribute('data-uid');
    if (!uid) return;
    if (ev.metaKey || ev.ctrlKey) {
      if (selectedUids.has(uid)) {
        selectedUids.delete(uid);
        row.classList.remove('rd-row-selected');
      } else {
        selectedUids.add(uid);
        row.classList.add('rd-row-selected');
      }
    } else if (selectedUids.has(uid)) {
      // A plain click on an already-selected row clears the whole
      // selection (Escape does the same).
      clearSelection(false);
    } else {
      clearSelection(false);
      selectedUids.add(uid);
      row.classList.add('rd-row-selected');
    }
    announceSelection();
  });

  document.addEventListener('keydown', function (ev) {
    if (ev.key !== 'Escape' || !selectedUids.size) return;
    clearSelection();
  });
})();
