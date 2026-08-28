/* rusty-dot interactive report behaviours (vanilla JS, no dependencies).
 *
 * Wiring contract with the Python side:
 *  - The figure is inline SVG inside #rd-figure.
 *  - Each dotplot panel's axes group is <g id="rd-panel-<row>-<col>">.
 *  - Each match layer is <g id="rd-matches-<row>-<col>-<layer>"> where
 *    <layer> is 'fwd', 'rev' or 'identity'.  The nth drawable child
 *    (<path> or <use>) of that group, in document order, corresponds to the
 *    nth entry of payload.panels['rd-panel-<row>-<col>'].segments[<layer>].
 *  - Match metadata lives in <script type="application/json" id="rd-data">.
 *
 * Behaviours:
 *  1. Click a panel        -> zoom the SVG viewBox to that panel and dim the
 *                             others; any subsequent click in the figure
 *                             (or Esc) returns to the full view.
 *  2. Wheel over the SVG   -> vertical pan; Shift+wheel -> horizontal pan;
 *                             Cmd/Ctrl+wheel (incl. trackpad pinch, which
 *                             browsers report as ctrl+wheel) -> uniform
 *                             aspect-preserving zoom centred on the cursor.
 *                             Pans are clamped to the figure bounds, and at
 *                             full view the page's natural scroll is kept.
 *  3. Drag inside a panel  -> rubber-band rectangle; on release, zoom to
 *                             that region (constrained to the panel; drags
 *                             under 5px count as plain clicks; Esc cancels).
 *  4. Click a match line   -> show coords/strand/identity (and sequence when
 *                             embedded) in the fixed bottom detail bar.
 *  5. Click a side-track   -> band the matching column (x track) or row (y
 *     feature                 track) behind the matches, and show the
 *                             feature in the detail bar.  Shift+click bands
 *                             every part of a multi-part feature.  Bands
 *                             toggle, stack, and clear on Esc.
 */
(function () {
  'use strict';

  var dataEl = document.getElementById('rd-data');
  var svg = document.querySelector('#rd-figure svg');
  if (!dataEl || !svg) return;

  var payload;
  try {
    payload = JSON.parse(dataEl.textContent);
  } catch (err) {
    console.error('rusty-dot: failed to parse embedded payload', err);
    return;
  }

  // ---------------------------------------------------------------------
  // viewBox helpers
  // ---------------------------------------------------------------------

  function getViewBox() {
    var vb = svg.getAttribute('viewBox');
    if (vb) {
      var p = vb.split(/[\s,]+/).map(Number);
      return { x: p[0], y: p[1], w: p[2], h: p[3] };
    }
    // Fall back to the declared width/height (strip pt/px units).
    var w = parseFloat(svg.getAttribute('width')) || svg.clientWidth;
    var h = parseFloat(svg.getAttribute('height')) || svg.clientHeight;
    return { x: 0, y: 0, w: w, h: h };
  }

  function setViewBox(vb) {
    svg.setAttribute(
      'viewBox',
      vb.x + ' ' + vb.y + ' ' + vb.w + ' ' + vb.h
    );
  }

  var homeViewBox = getViewBox();
  setViewBox(homeViewBox); // ensure the attribute exists for later edits

  /* Expand a viewBox to the figure's own aspect ratio (centred) so target
   * regions of arbitrary shape render without letterboxing surprises: every
   * viewBox we set matches the element's intrinsic aspect exactly. */
  function fitAspect(vb) {
    var ratio = homeViewBox.w / homeViewBox.h;
    var out = { x: vb.x, y: vb.y, w: vb.w, h: vb.h };
    if (out.w / out.h > ratio) {
      var newH = out.w / ratio;
      out.y -= (newH - out.h) / 2;
      out.h = newH;
    } else {
      var newW = out.h * ratio;
      out.x -= (newW - out.w) / 2;
      out.w = newW;
    }
    return out;
  }

  /* Escape text destined for innerHTML (sequence names come from user
   * FASTA headers and must not inject markup). */
  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, function (c) {
      return (
        { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]
      );
    });
  }

  /* Bounding box of an element expressed in the root SVG's viewBox space
   * (getBBox alone ignores ancestor transforms). */
  function bboxInRootSpace(el) {
    var box = el.getBBox();
    var toRoot = svg.getScreenCTM().inverse().multiply(el.getScreenCTM());
    var corners = [
      [box.x, box.y],
      [box.x + box.width, box.y],
      [box.x, box.y + box.height],
      [box.x + box.width, box.y + box.height],
    ];
    var xs = [];
    var ys = [];
    corners.forEach(function (c) {
      var pt = new DOMPoint(c[0], c[1]).matrixTransform(toRoot);
      xs.push(pt.x);
      ys.push(pt.y);
    });
    var x0 = Math.min.apply(null, xs);
    var y0 = Math.min.apply(null, ys);
    return {
      x: x0,
      y: y0,
      w: Math.max.apply(null, xs) - x0,
      h: Math.max.apply(null, ys) - y0,
    };
  }

  /* Pointer event position in viewBox coordinates. */
  function eventPoint(evt) {
    var pt = new DOMPoint(evt.clientX, evt.clientY);
    return pt.matrixTransform(svg.getScreenCTM().inverse());
  }

  // ---------------------------------------------------------------------
  // 1. Panel selection / zoom
  // ---------------------------------------------------------------------

  /* Panel groups are matched on the EXACT id shape, not the prefix alone.
   * matplotlib nests every gid'd artist in its own <g>, so a panel group
   * can contain descendants whose ids merely start the same way; matching
   * loosely once made the axes background masquerade as a second panel,
   * which dimmed the plot it was supposed to select and broke drill-down. */
  var PANEL_ID_RE = /^rd-panel-\d+-\d+$/;

  /* Nearest ancestor (inclusive) that is a real panel group, or null. */
  function closestPanel(el) {
    var node = el;
    while (node && node.nodeType === 1) {
      if (node === svg) return null; // never escape the figure
      if (PANEL_ID_RE.test(node.id || '')) return node;
      node = node.parentNode;
    }
    return null;
  }

  var panelGroups = Array.prototype.filter.call(
    svg.querySelectorAll('g[id^="rd-panel-"]'),
    function (g) {
      return PANEL_ID_RE.test(g.id);
    }
  );
  var selectedPanel = null;

  function resetView() {
    selectedPanel = null;
    setViewBox(homeViewBox);
    panelGroups.forEach(function (g) {
      g.classList.remove('rd-dim');
    });
  }

  function selectPanel(panel) {
    // No same-panel toggle here: while a panel is selected the delegated
    // handler below resets on any click, so this is only ever reached from
    // the unfocused state.
    selectedPanel = panel;
    var box = bboxInRootSpace(panel);
    var pad = Math.max(box.w, box.h) * 0.03;
    setViewBox(
      fitAspect({
        x: box.x - pad,
        y: box.y - pad,
        w: box.w + 2 * pad,
        h: box.h + 2 * pad,
      })
    );
    panelGroups.forEach(function (g) {
      g.classList.toggle('rd-dim', g !== panel);
    });
  }

  /* When an embedding app drills down on double-click (it sets
   * window.RD_DBLCLICK_DRILLDOWN from its injected bridge script), defer
   * the single-click focus zoom briefly so the first click of a
   * double-click never fires it — otherwise a double-click zooms into the
   * panel and then swaps views, a disorienting in-and-out sequence.
   * Standalone reports have no drill-down, so clicks stay instant. */
  var DBLCLICK_GRACE_MS = 300;
  var pendingPanelTimer = null;

  function cancelPendingPanelClick() {
    if (pendingPanelTimer !== null) {
      clearTimeout(pendingPanelTimer);
      pendingPanelTimer = null;
    }
  }

  document.addEventListener('dblclick', cancelPendingPanelClick, true);

  /* Click-to-focus only makes sense with several panels to choose from:
   * in a single-panel report (the drill-down view) it would just recentre
   * the one visible plot, so leave clicks alone there.
   *
   * One delegated listener on the SVG rather than one per panel, so that
   * while a panel is focused a click *anywhere* in the figure returns to
   * the full view — having to find and re-click the panel you zoomed into
   * is a poor way out of a zoom.  Clicks the match / annotation / track
   * handlers claim stop propagation and never reach here, so inspecting a
   * match inside a focused panel still shows its details instead of
   * throwing the zoom away. */
  if (panelGroups.length > 1) {
    svg.addEventListener('click', function (evt) {
      // A completed drag-zoom releases a click too; swallow that one.
      if (consumeDragClick()) return;
      if (selectedPanel !== null) {
        cancelPendingPanelClick();
        resetView();
        return;
      }
      var panel = closestPanel(evt.target);
      if (!panel) return; // click landed on figure margin, not a panel
      if (!window.RD_DBLCLICK_DRILLDOWN) {
        selectPanel(panel);
        return;
      }
      cancelPendingPanelClick();
      pendingPanelTimer = setTimeout(function () {
        pendingPanelTimer = null;
        selectPanel(panel);
      }, DBLCLICK_GRACE_MS);
    });
  }

  // ---------------------------------------------------------------------
  // 2. Wheel: pan (Shift = horizontal), Cmd/Ctrl+wheel = zoom
  // ---------------------------------------------------------------------

  /* Clamp a candidate viewBox to the figure bounds; zooming out to (or
   * past) full size snaps back to the home view. */
  function clampView(vb) {
    if (vb.w >= homeViewBox.w || vb.h >= homeViewBox.h) {
      return { x: homeViewBox.x, y: homeViewBox.y, w: homeViewBox.w, h: homeViewBox.h };
    }
    vb.x = clampVal(vb.x, homeViewBox.x, homeViewBox.x + homeViewBox.w - vb.w);
    vb.y = clampVal(vb.y, homeViewBox.y, homeViewBox.y + homeViewBox.h - vb.h);
    return vb;
  }

  svg.addEventListener(
    'wheel',
    function (evt) {
      var vb = getViewBox();
      if (evt.metaKey || evt.ctrlKey) {
        // Zoom (trackpad pinch arrives as ctrl+wheel).  Same factor on
        // both axes keeps the aspect ratio intact; anchoring on the
        // pointer keeps the hovered feature under the cursor.
        evt.preventDefault();
        var factor = Math.exp(evt.deltaY * 0.002);
        var pt = eventPoint(evt);
        vb.x = pt.x - (pt.x - vb.x) * factor;
        vb.y = pt.y - (pt.y - vb.y) * factor;
        vb.w *= factor;
        vb.h *= factor;
        setViewBox(clampView(vb));
        return;
      }
      // Pan: vertical by default, horizontal with Shift, clamped to the
      // figure.  When the pan cannot move (full view, or already at the
      // relevant edge), fall through WITHOUT preventDefault so the page's
      // natural scrolling still works.
      // Shift+wheel is reported via deltaX by some browsers/trackpads.
      var delta = evt.deltaY !== 0 ? evt.deltaY : evt.deltaX;
      if (evt.shiftKey) {
        var stepX = delta * (vb.w / (svg.clientWidth || homeViewBox.w));
        var newX = clampVal(
          vb.x + stepX,
          homeViewBox.x,
          Math.max(homeViewBox.x, homeViewBox.x + homeViewBox.w - vb.w)
        );
        if (newX === vb.x) return;
        evt.preventDefault();
        vb.x = newX;
      } else {
        var stepY = delta * (vb.h / (svg.clientHeight || homeViewBox.h));
        var newY = clampVal(
          vb.y + stepY,
          homeViewBox.y,
          Math.max(homeViewBox.y, homeViewBox.y + homeViewBox.h - vb.h)
        );
        if (newY === vb.y) return;
        evt.preventDefault();
        vb.y = newY;
      }
      setViewBox(vb);
    },
    { passive: false }
  );

  // ---------------------------------------------------------------------
  // 3. Drag-rectangle region zoom within a panel
  // ---------------------------------------------------------------------

  var DRAG_THRESHOLD_PX = 5; // below this a mousedown+up is a plain click
  var drag = null; // active drag state, or null
  var suppressNextClick = false;

  /* Consume the click that the browser fires right after a completed
   * drag-zoom, so it does not also trigger panel/match click actions. */
  function consumeDragClick() {
    if (suppressNextClick) {
      suppressNextClick = false;
      return true;
    }
    return false;
  }

  function clampVal(v, lo, hi) {
    return Math.min(hi, Math.max(lo, v));
  }

  /* Current selection rectangle in viewBox coordinates, clamped to the
   * panel the drag started in. */
  function dragRegion(evt) {
    var p1 = eventPoint(evt);
    var b = drag.panelBox;
    var x0 = clampVal(Math.min(drag.p0.x, p1.x), b.x, b.x + b.w);
    var x1 = clampVal(Math.max(drag.p0.x, p1.x), b.x, b.x + b.w);
    var y0 = clampVal(Math.min(drag.p0.y, p1.y), b.y, b.y + b.h);
    var y1 = clampVal(Math.max(drag.p0.y, p1.y), b.y, b.y + b.h);
    return { x: x0, y: y0, w: x1 - x0, h: y1 - y0 };
  }

  function cancelDrag() {
    if (drag && drag.rect && drag.rect.parentNode) {
      drag.rect.parentNode.removeChild(drag.rect);
    }
    drag = null;
  }

  svg.addEventListener('mousedown', function (evt) {
    if (evt.button !== 0) return;
    suppressNextClick = false;
    // Only drags starting inside a panel arm region select.
    var panel = closestPanel(evt.target);
    if (!panel) return;
    var pm = panel.id.match(/^rd-panel-(\d+)-(\d+)$/);
    drag = {
      cx0: evt.clientX,
      cy0: evt.clientY,
      p0: eventPoint(evt),
      panelBox: bboxInRootSpace(panel),
      active: false,
      rect: null,
      // Shift turns the gesture into a box select of match segments in
      // the origin panel; without it the drag zooms as before.
      shift: evt.shiftKey,
      row: pm ? parseInt(pm[1], 10) : -1,
      col: pm ? parseInt(pm[2], 10) : -1,
    };
  });

  document.addEventListener('mousemove', function (evt) {
    if (!drag) return;
    var dx = evt.clientX - drag.cx0;
    var dy = evt.clientY - drag.cy0;
    if (!drag.active) {
      if (Math.sqrt(dx * dx + dy * dy) < DRAG_THRESHOLD_PX) return;
      // Passed the threshold: this gesture is a region select, not a click.
      drag.active = true;
      drag.rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      drag.rect.setAttribute(
        'class',
        drag.shift ? 'rd-select-rect' : 'rd-drag-rect'
      );
      svg.appendChild(drag.rect);
    }
    evt.preventDefault();
    var r = dragRegion(evt);
    drag.rect.setAttribute('x', r.x);
    drag.rect.setAttribute('y', r.y);
    drag.rect.setAttribute('width', r.w);
    drag.rect.setAttribute('height', r.h);
  });

  function rectsIntersect(a, b) {
    return a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h;
  }

  function segBox(entry) {
    if (!entry.bbox) entry.bbox = bboxInRootSpace(entry.el);
    return entry.bbox;
  }

  /* Liang-Barsky line/rect clip: does the segment (x1,y1)-(x2,y2) touch
   * rect r?  Needed because a long diagonal's bounding box covers most of
   * its panel -- bbox overlap alone would select it from a small box
   * dragged in the empty corner. */
  function segIntersectsRect(x1, y1, x2, y2, r) {
    var dx = x2 - x1;
    var dy = y2 - y1;
    var p = [-dx, dx, -dy, dy];
    var q = [x1 - r.x, r.x + r.w - x1, y1 - r.y, r.y + r.h - y1];
    var t0 = 0;
    var t1 = 1;
    for (var i = 0; i < 4; i++) {
      if (p[i] === 0) {
        if (q[i] < 0) return false;
        continue;
      }
      var t = q[i] / p[i];
      if (p[i] < 0) {
        if (t > t1) return false;
        if (t > t0) t0 = t;
      } else {
        if (t < t0) return false;
        if (t < t1) t1 = t;
      }
    }
    return t0 <= t1;
  }

  /* A match is drawn along one diagonal of its bbox: forward matches run
   * top-left to bottom-right (query and target advance together, y grows
   * downward), reverse matches along the other diagonal. */
  function segTouchesRegion(entry, region) {
    var b = segBox(entry);
    if (!rectsIntersect(region, b)) return false;
    if (entry.strand === '-') {
      return segIntersectsRect(b.x, b.y + b.h, b.x + b.w, b.y, region);
    }
    return segIntersectsRect(b.x, b.y, b.x + b.w, b.y + b.h, region);
  }

  document.addEventListener('mouseup', function (evt) {
    if (!drag) return;
    var wasActive = drag.active;
    var region = wasActive ? dragRegion(evt) : null;
    var isSelect = drag.shift;
    var row = drag.row;
    var col = drag.col;
    cancelDrag();
    if (!wasActive) return; // plain click: let the click handlers run
    suppressNextClick = true;
    if (!(region.w > 0 && region.h > 0)) return;
    if (isSelect) {
      // Box select: highlight every match segment intersecting the box
      // (bbox approximation) without opening the detail bar or fetching
      // sequences.  An empty box clears the selection.
      var hits = segmentRegistry.filter(function (e) {
        return (
          e.row === row &&
          e.col === col &&
          !e.el.classList.contains('rd-len-hidden') &&
          segTouchesRegion(e, region)
        );
      });
      closeDetailBar();
      setMatchSelection(hits);
    } else {
      setViewBox(fitAspect(region));
    }
  });

  document.addEventListener('keydown', function (evt) {
    if (evt.key === 'Escape') {
      cancelDrag();
      // First Escape clears selections/bands (and resets the view); only
      // an Escape with nothing left to clear asks the embedding app to
      // leave fullscreen.
      var hadSelection = hasAnySelection();
      resetView();
      clearAllSelection();
      clearBands();
      if (!hadSelection && window.parent !== window) {
        window.parent.postMessage({ type: 'rd-esc' }, '*');
      }
    }
  });

  // ---------------------------------------------------------------------
  // 4. Match click -> detail bar
  // ---------------------------------------------------------------------

  var detail = document.getElementById('rd-detail');
  var detailCoords = document.getElementById('rd-detail-coords');
  var detailSeq = document.getElementById('rd-detail-seq');
  var detailActions = document.getElementById('rd-detail-actions');
  var copyQueryBtn = document.getElementById('rd-copy-query');
  var copyTargetBtn = document.getElementById('rd-copy-target');
  // Match selection persists after the detail bar closes and can hold
  // several segments (Shift+drag box select).  Keyed "row:col:layer:idx"
  // onto segmentRegistry entries.  Annotation/track features keep their
  // own single slot so a feature click never silently drops a match
  // selection's row/column highlights (and vice versa).
  var selectedSegs = {};
  var selectedFeatureEl = null;
  // Key ("row:col:layer:idx") of the aligned-sequence request in flight;
  // guards against stale 'rd-seq-response' messages after another click.
  var pendingSeqKey = null;
  // Copy-button state for the selected match.  Full sequences are NOT
  // shipped with the preview (that made clicks slow on megabase
  // matches); they are fetched from the embedding app on the first copy
  // press and cached here, so repeat copies are instant.
  var copyState = {
    query: { avail: false, seq: null, pending: false },
    target: { avail: false, seq: null, pending: false },
  };
  // Segment identity of the current selection, echoed in copy requests.
  var currentSeg = null;

  var SIDE_LABELS = { query: 'query seq', target: 'target seq' };

  function copyBtn(side) {
    return side === 'target' ? copyTargetBtn : copyQueryBtn;
  }

  // The button says what the *next* press will do: before the sequence has
  // been fetched from the embedding app a press only fetches it.
  function labelFor(side) {
    return (copyState[side].seq ? 'Copy ' : 'Fetch ') + SIDE_LABELS[side];
  }

  function setCopyState(qAvail, tAvail, qSeq, tSeq) {
    copyState.query = { avail: !!qAvail, seq: qSeq || null, pending: false };
    copyState.target = { avail: !!tAvail, seq: tSeq || null, pending: false };
    ['query', 'target'].forEach(function (side) {
      var btn = copyBtn(side);
      btn.hidden = !copyState[side].avail;
      btn.disabled = false;
      btn.textContent = labelFor(side);
    });
    detailActions.hidden = !(copyState.query.avail || copyState.target.avail);
  }

  function copyText(text, btn, side) {
    var done = function (ok) {
      btn.textContent = ok ? 'Copied!' : 'Copy failed';
      setTimeout(function () {
        // Recomputed at revert time: the state may have gained a cached
        // sequence since the press, flipping 'Fetch' to 'Copy'.
        btn.textContent = labelFor(side);
      }, 1200);
    };
    // execCommand works on user activation even in a sandboxed
    // (opaque-origin) frame; the async clipboard API is the fallback.
    var ta = document.createElement('textarea');
    ta.value = text;
    ta.style.position = 'fixed';
    ta.style.opacity = '0';
    document.body.appendChild(ta);
    ta.select();
    var ok = false;
    try {
      ok = document.execCommand('copy');
    } catch (e) {
      ok = false;
    }
    document.body.removeChild(ta);
    if (!ok && navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(
        function () {
          done(true);
        },
        function () {
          done(false);
        }
      );
    } else {
      done(ok);
    }
  }

  function wireCopy(btn, side) {
    btn.addEventListener('click', function () {
      var st = copyState[side];
      if (!st.avail || st.pending) return;
      if (st.seq) {
        copyText(st.seq, btn, side);
        return;
      }
      if (!currentSeg || window.parent === window) return;
      st.pending = true;
      btn.disabled = true;
      btn.textContent = 'Fetching…';
      window.parent.postMessage(
        Object.assign({ type: 'rd-copy-request', side: side }, currentSeg),
        '*'
      );
    });
  }

  wireCopy(copyQueryBtn, 'query');
  wireCopy(copyTargetBtn, 'target');

  function segKey(seg) {
    return seg ? seg.row + ':' + seg.col + ':' + seg.layer + ':' + seg.idx : null;
  }

  function handleCopyResponse(msg) {
    var side = msg.side === 'target' ? 'target' : 'query';
    var btn = copyBtn(side);
    var key = msg.row + ':' + msg.col + ':' + msg.layer + ':' + msg.idx;
    if (key !== segKey(currentSeg)) return; // stale: selection changed
    copyState[side].pending = false;
    btn.disabled = false;
    if (typeof msg.seq !== 'string' || !msg.seq) {
      btn.textContent = 'Copy failed';
      setTimeout(function () {
        btn.textContent = labelFor(side); // still 'Fetch …': nothing cached
      }, 1200);
      return;
    }
    // Cache before touching the label, so every labelFor() below reads
    // 'Copy …'.
    copyState[side].seq = msg.seq;
    // The fetch usually completes inside the click's transient-activation
    // window, so the async clipboard write succeeds directly; if it
    // doesn't, the sequence is now cached and the next press copies
    // synchronously -- which is exactly what the 'Copy …' label promises.
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(msg.seq).then(
        function () {
          btn.textContent = 'Copied!';
          setTimeout(function () {
            btn.textContent = labelFor(side);
          }, 1200);
        },
        function () {
          btn.textContent = labelFor(side);
        }
      );
    } else {
      btn.textContent = labelFor(side);
    }
  }

  /* Close the bar without touching the selection: the × keeps the match
   * highlighted so its row/column stay marked across the grid. */
  function closeDetailBar() {
    detail.hidden = true;
    pendingSeqKey = null;
    currentSeg = null;
    setCopyState(false, false, null, null);
  }

  function clearMatchSelection() {
    Object.keys(selectedSegs).forEach(function (k) {
      selectedSegs[k].el.classList.remove('rd-selected-match');
    });
    selectedSegs = {};
    updateRowColHighlights();
  }

  /* Replace the match selection with the given registry entries. */
  function setMatchSelection(entries) {
    clearMatchSelection();
    entries.forEach(function (e) {
      e.el.classList.add('rd-selected-match');
      selectedSegs[e.key] = e;
    });
    updateRowColHighlights();
  }

  function clearAllSelection() {
    closeDetailBar();
    clearMatchSelection();
    if (selectedFeatureEl) {
      selectedFeatureEl.classList.remove('rd-selected-match');
      selectedFeatureEl = null;
    }
    announceFeatureSelection(null);
    notifySelectionState();
  }

  /* Anything Escape would clear before it should mean "exit fullscreen". */
  function hasAnySelection() {
    return (
      Object.keys(selectedSegs).length > 0 ||
      !!selectedFeatureEl ||
      Object.keys(activeBands).length > 0
    );
  }

  /* Keep the embedding app abreast of whether anything is selected: its
   * own Escape handling clears selections first and exits fullscreen only
   * when there is nothing left to clear. */
  function notifySelectionState() {
    if (window.parent === window) return;
    window.parent.postMessage(
      { type: 'rd-selection-state', any: hasAnySelection() },
      '*'
    );
  }

  /* Every (query, target) pair this report displays.  Sent with selection
   * updates so the app can keep selections on pairs a narrower view (the
   * drill-down) does not show, and merge in the changes for those it does. */
  var reportPairs = Object.keys(payload.panels || {}).map(function (gid) {
    var p = payload.panels[gid];
    return [p.query, p.target];
  });

  /* Mirror the match selection to the app in view-independent terms so it
   * survives the iframe being rebuilt (overview <-> drill-down swaps). */
  function announceMatchSelection() {
    if (window.parent === window) return;
    var matches = [];
    Object.keys(selectedSegs).forEach(function (k) {
      var e = selectedSegs[k];
      if (e.q === undefined) return;
      matches.push({
        q: e.q, t: e.t, layer: e.layer,
        qs: e.qs, qe: e.qe, ts: e.ts, te: e.te,
      });
    });
    window.parent.postMessage(
      { type: 'rd-match-selection', pairs: reportPairs, matches: matches },
      '*'
    );
  }

  /* Same for the (single) selected annotation/track feature. */
  function announceFeatureSelection(feat) {
    if (window.parent === window) return;
    window.parent.postMessage(
      {
        type: 'rd-feature-selection',
        feature: feat
          ? {
              seqname: feat.seqname,
              start: feat.start,
              end: feat.end,
              type: feat.type,
              strand: feat.strand,
            }
          : null,
      },
      '*'
    );
  }

  document
    .getElementById('rd-detail-close')
    .addEventListener('click', closeDetailBar);

  function showMatchDetail(panelGid, layer, idx, el, additive) {
    var panel = payload.panels[panelGid];
    if (!panel) return;
    var seg = panel.segments[layer][idx];
    if (!seg) return;

    var pm = panelGid.match(/^rd-panel-(\d+)-(\d+)$/);
    if (!pm) return;
    var row = parseInt(pm[1], 10);
    var col = parseInt(pm[2], 10);
    var key = row + ':' + col + ':' + layer + ':' + idx;
    if (additive && selectedSegs[key]) {
      // Cmd/Ctrl-click on a selected match drops just that match.
      selectedSegs[key].el.classList.remove('rd-selected-match');
      delete selectedSegs[key];
      updateRowColHighlights();
      if (!detail.hidden && currentSeg && segKey(currentSeg) === key) {
        closeDetailBar();
      }
      return;
    }
    // Clicking the sole-selected match while its bar is open deselects;
    // with the bar closed (after ×) the same click just reopens the bar.
    if (!additive && !detail.hidden && selectedSegs[key] &&
        Object.keys(selectedSegs).length === 1) {
      clearAllSelection();
      return;
    }

    var strand = layer === 'rev' ? '-' : '+';
    var identity = null;
    if (layer === 'identity') {
      identity = seg[4];
      if (seg.length > 5) strand = seg[5];
    }

    var html =
      '<b>' + escapeHtml(panel.query) + '</b>:' + seg[0].toLocaleString() +
      '–' + seg[1].toLocaleString() +
      ' &times; <b>' + escapeHtml(panel.target) + '</b>:' + seg[2].toLocaleString() +
      '–' + seg[3].toLocaleString() +
      ' &middot; strand ' + escapeHtml(strand) +
      ' &middot; ' + (seg[1] - seg[0]).toLocaleString() + ' bp';
    if (identity !== null) {
      html += ' &middot; identity ' + (identity * 100).toFixed(1) + '%';
    }
    detailCoords.innerHTML = html;

    var seq =
      payload.has_sequences && panel.seqs && panel.seqs[layer]
        ? panel.seqs[layer][idx]
        : null;
    detailSeq.classList.remove('rd-aligned');
    currentSeg = null;
    setCopyState(false, false, null, null);
    if (seq) {
      detailSeq.textContent = seq;
      detailSeq.hidden = false;
      setCopyState(true, false, seq, null);
    } else if (window.parent && window.parent !== window) {
      // Embedded report: ask the embedding app for this match's
      // sequences (a CIGAR-based alignment when available, otherwise the
      // raw query/target slices).  The segment coordinates travel along
      // so non-identity layers need no server-side record lookup.
      pendingSeqKey = key;
      currentSeg = {
        row: row,
        col: col,
        layer: layer,
        idx: idx,
        qs: seg[0],
        qe: seg[1],
        ts: seg[2],
        te: seg[3],
        strand: strand,
      };
      detailSeq.textContent = 'Fetching sequences…';
      detailSeq.hidden = false;
      window.parent.postMessage(
        Object.assign({ type: 'rd-match-select' }, currentSeg),
        '*'
      );
    } else {
      detailSeq.textContent = '';
      detailSeq.hidden = true;
    }
    detail.hidden = false;

    if (selectedFeatureEl) {
      selectedFeatureEl.classList.remove('rd-selected-match');
      selectedFeatureEl = null;
      announceFeatureSelection(null);
    }
    var entry = registryByKey[key] || {
      el: el, key: key, row: row, col: col, layer: layer,
      q: panel.query, t: panel.target,
      qs: seg[0], qe: seg[1], ts: seg[2], te: seg[3],
    };
    if (additive) {
      // Cmd/Ctrl-click adds to the selection instead of replacing it.
      entry.el.classList.add('rd-selected-match');
      selectedSegs[entry.key] = entry;
      updateRowColHighlights();
    } else {
      setMatchSelection([entry]);
    }
  }

  /* Wire up every match group: index its drawable children in document
   * order (the serialisation contract) and overlay each with an invisible
   * wide-stroke clone so hairline matches are easy to click.  Each segment
   * is also registered with its query-side length so the embedding app can
   * filter by minimum match length client-side (behaviour 6). */
  var segmentRegistry = []; // {el, hit, qlen, row, col, layer, idx, key, bbox}
  var registryByKey = {};
  var matchGroups = Array.prototype.slice.call(
    svg.querySelectorAll('g[id^="rd-matches-"]')
  );
  matchGroups.forEach(function (group) {
    // gid: rd-matches-<row>-<col>-<layer>
    var m = group.id.match(/^rd-matches-(\d+)-(\d+)-(\w+)$/);
    if (!m) return;
    var panelGid = 'rd-panel-' + m[1] + '-' + m[2];
    var layer = m[3];
    var row = parseInt(m[1], 10);
    var col = parseInt(m[2], 10);
    var panel = payload.panels[panelGid];
    var segs = panel && panel.segments ? panel.segments[layer] : null;

    var children = Array.prototype.slice.call(
      group.querySelectorAll('path, use')
    );
    children.forEach(function (el, idx) {
      var onClick = function (evt) {
        // A completed drag-zoom releases a click too; swallow it fully
        // (stopPropagation so the panel handler does not fire either).
        if (consumeDragClick()) {
          evt.stopPropagation();
          return;
        }
        evt.stopPropagation();
        showMatchDetail(panelGid, layer, idx, el, evt.metaKey || evt.ctrlKey);
      };
      el.addEventListener('click', onClick);
      // Invisible widened hit target stacked on top of the original.
      var hit = el.cloneNode(false);
      hit.removeAttribute('style');
      hit.removeAttribute('id');
      hit.setAttribute('class', 'rd-hit');
      hit.addEventListener('click', onClick);
      group.appendChild(hit);
      if (segs && segs[idx]) {
        var entry = {
          el: el,
          hit: hit,
          // Query-side span, matching the server-side min_length semantics.
          qlen: segs[idx][1] - segs[idx][0],
          row: row,
          col: col,
          layer: layer,
          idx: idx,
          // Which bbox diagonal the segment is drawn along (box select
          // tests the actual line, not the bbox): identity segments carry
          // their strand at [5], the fwd/rev layers imply it.
          strand:
            layer === 'rev'
              ? '-'
              : layer === 'identity' && segs[idx].length > 5
                ? segs[idx][5]
                : '+',
          key: row + ':' + col + ':' + layer + ':' + idx,
          // View-independent identity: display names plus data coords, so
          // a selection made in the overview can be re-applied in the
          // drill-down (whose grid position differs) and back.
          q: panel.query,
          t: panel.target,
          qs: segs[idx][0],
          qe: segs[idx][1],
          ts: segs[idx][2],
          te: segs[idx][3],
          // Root-space bbox, filled lazily on the first box select.  The
          // screen-CTM composition cancels the current viewBox, so one
          // computation stays valid at any zoom.
          bbox: null,
        };
        segmentRegistry.push(entry);
        registryByKey[entry.key] = entry;
      }
    });
  });

  /* Axes backgrounds carry their own gid prefix (deliberately outside
   * rd-panel- so they never masquerade as panels); selection bands are
   * inserted right after them so they paint behind the match strokes. */
  var plotBgs = [];
  Array.prototype.forEach.call(
    svg.querySelectorAll('[id^="rd-plotbg-"]'),
    function (el) {
      var m = el.id.match(/^rd-plotbg-(\d+)-(\d+)$/);
      if (m) {
        plotBgs.push({
          el: el,
          row: parseInt(m[1], 10),
          col: parseInt(m[2], 10),
          box: null, // root-space bbox, lazy (viewBox-independent)
        });
      }
    }
  );

  function bgBox(bg) {
    if (!bg.box) bg.box = bboxInRootSpace(bg.el);
    return bg.box;
  }

  /* Merge overlapping/touching 1-D intervals ({a, b} with a <= b). */
  function mergeRanges(ranges) {
    ranges.sort(function (u, v) {
      return u.a - v.a;
    });
    var out = [];
    ranges.forEach(function (r) {
      var last = out[out.length - 1];
      if (last && r.a <= last.b) {
        if (r.b > last.b) last.b = r.b;
      } else {
        out.push({ a: r.a, b: r.b });
      }
    });
    return out;
  }

  /* Selection bands: each selected match paints its query range as a
   * horizontal band through every panel in its grid row and its target
   * range as a vertical band through every panel in its column (the same
   * look as the drill-down track bands).  Rows share the query axis and
   * columns the target axis, so the segment's own root-space extent is
   * valid across its whole row/column.  Overlapping ranges are merged so
   * a dense box-selection stays a handful of rects.  On a 1x1 drilldown
   * this degenerates to a crosshair through the single panel. */
  var selectionBandEls = [];

  function updateRowColHighlights() {
    selectionBandEls.forEach(function (el) {
      if (el.parentNode) el.parentNode.removeChild(el);
    });
    selectionBandEls = [];
    var rowRanges = {}; // row -> [{a, b}] y extents (query)
    var colRanges = {}; // col -> [{a, b}] x extents (target)
    Object.keys(selectedSegs).forEach(function (k) {
      var e = selectedSegs[k];
      if (!e.el || !e.el.getBBox) return;
      var b = segBox(e);
      (rowRanges[e.row] = rowRanges[e.row] || []).push({ a: b.y, b: b.y + b.h });
      (colRanges[e.col] = colRanges[e.col] || []).push({ a: b.x, b: b.x + b.w });
    });
    Object.keys(rowRanges).forEach(function (r) {
      rowRanges[r] = mergeRanges(rowRanges[r]);
    });
    Object.keys(colRanges).forEach(function (c) {
      colRanges[c] = mergeRanges(colRanges[c]);
    });
    plotBgs.forEach(function (bg) {
      var pb = bgBox(bg);
      (rowRanges[bg.row] || []).forEach(function (r) {
        addSelectionBand(bg, pb.x, r.a, pb.w, r.b - r.a);
      });
      (colRanges[bg.col] || []).forEach(function (r) {
        addSelectionBand(bg, r.a, pb.y, r.b - r.a, pb.h);
      });
    });
    notifySelectionState();
    announceMatchSelection();
  }

  function addSelectionBand(bg, x, y, w, h) {
    var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rect.setAttribute('class', 'rd-selband');
    rect.setAttribute('x', x);
    rect.setAttribute('y', y);
    rect.setAttribute('width', w);
    rect.setAttribute('height', h);
    // Right after the background: behind the match strokes by document
    // order, same trick as the track bands (.rd-band).
    bg.el.parentNode.insertBefore(rect, bg.el.nextSibling);
    selectionBandEls.push(rect);
  }

  // ---------------------------------------------------------------------
  // 5. Annotation feature click -> detail bar
  // ---------------------------------------------------------------------

  function showAnnotationDetail(panelGid, idx, el) {
    var panel = payload.panels[panelGid];
    if (!panel || !panel.annotations) return;
    showFeatureDetail(panel.annotations[idx], el);
  }

  /* Render one GFF feature into the detail bar.  Shared by the diagonal
   * squares and the side tracks, which carry the same field set. */
  function showFeatureDetail(feat, el) {
    if (!feat) return;

    var label = feat.name || feat.id || feat.type;
    var html =
      '<b>' + escapeHtml(String(label)) + '</b>' +
      ' &middot; ' + escapeHtml(feat.type) +
      ' &middot; ' + escapeHtml(feat.seqname) + ':' +
      feat.start.toLocaleString() + '–' + feat.end.toLocaleString() +
      ' &middot; strand ' + escapeHtml(feat.strand) +
      ' &middot; ' + (feat.end - feat.start).toLocaleString() + ' bp';
    if (feat.parent) {
      html += ' &middot; parent ' + escapeHtml(feat.parent);
    }
    if (feat.source && feat.source !== '.') {
      html += ' &middot; source ' + escapeHtml(feat.source);
    }
    if (feat.source_file) {
      html += ' &middot; file ' + escapeHtml(feat.source_file);
    }
    detailCoords.innerHTML = html;
    detailSeq.textContent = '';
    detailSeq.hidden = true;
    detail.hidden = false;

    // A feature click replaces any match selection (and its row/column
    // marks) but lives in its own slot; Escape clears it.
    clearMatchSelection();
    if (selectedFeatureEl) selectedFeatureEl.classList.remove('rd-selected-match');
    selectedFeatureEl = el;
    el.classList.add('rd-selected-match');
    announceFeatureSelection(feat);
    notifySelectionState();
  }

  /* Wire up every diagonal-annotation group: the n-th drawable child of
   * 'rd-annot-<row>-<col>' corresponds to panels[gid].annotations[n]
   * (serialisation contract — patch draw order equals SVG child order). */
  var annotRegistry = []; // {el, feat} — for re-applying a restored selection
  var annotGroups = Array.prototype.slice.call(
    svg.querySelectorAll('g[id^="rd-annot-"]')
  );
  annotGroups.forEach(function (group) {
    var m = group.id.match(/^rd-annot-(\d+)-(\d+)$/);
    if (!m) return;
    var panelGid = 'rd-panel-' + m[1] + '-' + m[2];
    var feats = (payload.panels[panelGid] || {}).annotations || [];

    var children = Array.prototype.slice.call(
      group.querySelectorAll('path, use')
    );
    children.forEach(function (el, idx) {
      el.classList.add('rd-annot');
      if (feats[idx]) annotRegistry.push({ el: el, feat: feats[idx] });
      el.addEventListener('click', function (evt) {
        if (consumeDragClick()) {
          evt.stopPropagation();
          return;
        }
        evt.stopPropagation();
        showAnnotationDetail(panelGid, idx, el);
      });
    });
  });

  // ---------------------------------------------------------------------
  // 5b. Side-track feature click -> highlight band across the plot
  // ---------------------------------------------------------------------
  //
  // Clicking a feature in the x track paints a translucent column behind
  // the matches; the y track paints a row.  A multi-part feature (a
  // spliced CDS) is drawn as one patch per part, so Shift+click on any
  // part bands every segment of the group.  Bands toggle, several can be
  // active at once, and Escape clears them.
  //
  // No bp->pixel arithmetic is involved: the track axes are created with
  // sharex/sharey against the main panel, so a clicked feature's own
  // bounding box already gives the band's extent along the sequence axis,
  // and the panel background rect gives the other one.  That also keeps
  // this correct under bbox_inches='tight', which shifts everything a
  // Python-side computation would have had to predict.

  var trackData = payload.tracks || null;
  var activeBands = {}; // gid -> <rect>
  var bandEntries = {}; // gid -> {axis, payload entry}
  var bandLayer = null;

  /* Tell an embedding app which features are banded, so a figure it saves
   * can carry the same highlights.  The bands here are DOM rects measured
   * off the rendered SVG; a saved figure is drawn from coordinates, so the
   * app gets the feature's own numbers rather than pixel geometry. */
  function publishBands() {
    if (window.parent === window) return; // standalone report: no app
    var bands = Object.keys(activeBands).map(function (gid) {
      var held = bandEntries[gid] || {};
      var entry = held.entry || {};
      return {
        gid: gid,
        axis: held.axis,
        seqname: entry.seqname,
        start: entry.start,
        end: entry.end,
        color: entry.color || patchFill(document.getElementById(gid)),
      };
    });
    window.parent.postMessage({ type: 'rd-bands', bands: bands }, '*');
    notifySelectionState();
  }

  function panelBackground() {
    return svg.querySelector('[id="rd-plotbg-0-0"]');
  }

  /* Lazily create the band group, inserted directly after the panel
   * background so document order puts bands behind the match strokes. */
  function ensureBandLayer() {
    if (bandLayer && bandLayer.parentNode) return bandLayer;
    var bg = panelBackground();
    if (!bg) return null;
    var host = bg.parentNode;
    bandLayer = document.createElementNS('http://www.w3.org/2000/svg', 'g');
    bandLayer.setAttribute('id', 'rd-bands');
    host.insertBefore(bandLayer, bg.nextSibling);
    return bandLayer;
  }

  function clearBands() {
    Object.keys(activeBands).forEach(function (gid) {
      var rect = activeBands[gid];
      if (rect && rect.parentNode) rect.parentNode.removeChild(rect);
      var el = document.getElementById(gid);
      if (el) el.classList.remove('rd-track-active');
    });
    activeBands = {};
    bandEntries = {};
    publishBands();
  }

  /* Colour a track patch is actually painted with.  A gid-tagged artist is
   * wrapped in a <g> by matplotlib's SVG backend and the fill sits on the
   * child <path>'s style, so the group's own computed fill is the useless
   * default black. */
  function patchFill(el) {
    var kid = el.querySelector('path, use, rect, polygon') || el;
    var fill =
      kid.getAttribute('fill') || window.getComputedStyle(kid).fill || '';
    if (!fill || fill === 'none' || fill === 'rgb(0, 0, 0)') return '#888888';
    // getComputedStyle reports 'rgb(r, g, b)'.  SVG takes that happily, but
    // the app forwards this colour to matplotlib for saved figures, which
    // does not, so normalise to hex here rather than parsing CSS there.
    var m = /^rgba?\((\d+),\s*(\d+),\s*(\d+)/.exec(fill);
    if (m) {
      return (
        '#' +
        [m[1], m[2], m[3]]
          .map(function (v) {
            return ('0' + Number(v).toString(16)).slice(-2);
          })
          .join('')
      );
    }
    return fill;
  }

  function addBand(gid, entry, axis) {
    var el = document.getElementById(gid);
    var bg = panelBackground();
    var layer = ensureBandLayer();
    if (!el || !bg || !layer) return;
    var feat = bboxInRootSpace(el);
    var panel = bboxInRootSpace(bg);
    var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    if (axis === 'x') {
      rect.setAttribute('x', feat.x);
      rect.setAttribute('width', feat.w);
      rect.setAttribute('y', panel.y);
      rect.setAttribute('height', panel.h);
    } else {
      rect.setAttribute('x', panel.x);
      rect.setAttribute('width', panel.w);
      rect.setAttribute('y', feat.y);
      rect.setAttribute('height', feat.h);
    }
    rect.setAttribute('fill', (entry && entry.color) || patchFill(el));
    rect.setAttribute('class', 'rd-band');
    layer.appendChild(rect);
    activeBands[gid] = rect;
    bandEntries[gid] = { axis: axis, entry: entry };
    el.classList.add('rd-track-active');
  }

  function removeBand(gid) {
    var rect = activeBands[gid];
    if (rect && rect.parentNode) rect.parentNode.removeChild(rect);
    delete activeBands[gid];
    delete bandEntries[gid];
    var el = document.getElementById(gid);
    if (el) el.classList.remove('rd-track-active');
  }

  function toggleBands(gids, entries, axis) {
    // All-on -> all-off, otherwise fill in the missing ones, so a
    // Shift+click on a partly-banded group completes it rather than
    // flip-flopping each segment independently.
    var allOn = gids.every(function (g) {
      return activeBands[g];
    });
    gids.forEach(function (gid, i) {
      if (allOn) removeBand(gid);
      else if (!activeBands[gid]) addBand(gid, entries[i], axis);
    });
    publishBands();
  }

  // Bands driven by the embedding app's annotations-table selection
  // ('rd-highlight-features' messages).  Kept in their own set so
  // replacing the selection never removes a band the user click-toggled
  // on a track glyph themselves.
  var selectionGids = {};
  var lastSelection = [];

  function entryMatches(entry, f) {
    // uid is decisive when it matches; on a mismatch the coordinates
    // still decide, because a self panel draws the query's features on
    // the x axis under target-role uids.
    if (f.uid && entry.uid && entry.uid === f.uid) return true;
    return (
      entry.seqname === f.seqname &&
      Number(entry.start) === Number(f.start) &&
      Number(entry.end) === Number(f.end) &&
      entry.type === f.type &&
      entry.strand === f.strand
    );
  }

  /* Whether the report currently has real layout.  Inside a hidden tab
   * pane (display:none) every box measures zero, so bands drawn then
   * would be degenerate -- and the hidden iframe's window keeps its old
   * size, so size-based signals never fire. */
  function panelVisible() {
    var bg = panelBackground();
    if (!bg) return false;
    var r = bg.getBoundingClientRect();
    return r.width > 0 && r.height > 0;
  }

  function applySelectionHighlights(features) {
    lastSelection = Array.isArray(features) ? features : [];
    Object.keys(selectionGids).forEach(function (gid) {
      removeBand(gid);
    });
    selectionGids = {};
    if (!panelVisible()) {
      // Hidden: hold the selection; the visibility observer below
      // redraws it as soon as the pane is shown again.
      publishBands();
      return;
    }
    if (trackData && Array.isArray(features)) {
      features.forEach(function (f) {
        if (!f) return;
        // Both axes: a self panel draws the same sequence twice, and the
        // role-prefixed uids plus seqnames keep cross pairs from
        // over-matching.
        ['x', 'y'].forEach(function (axis) {
          (trackData[axis] || []).forEach(function (entry) {
            if (!entryMatches(entry, f)) return;
            if (activeBands[entry.gid]) return; // user's own band wins
            addBand(entry.gid, entry, axis);
            selectionGids[entry.gid] = true;
          });
        });
      });
    }
    publishBands();
  }

  // A selection arriving while the report's tab pane is hidden cannot be
  // drawn (see panelVisible); redraw it when the pane is shown again.
  // IntersectionObserver is the one signal that reliably reports the
  // display:none -> shown transition here: the hidden iframe's window
  // keeps its stale size, so neither the resize event nor a
  // ResizeObserver ever fires inside it.
  var selectionRedrawTimer = null;

  function redrawSelectionSoon() {
    if (!lastSelection.length) return;
    clearTimeout(selectionRedrawTimer);
    selectionRedrawTimer = setTimeout(function () {
      if (panelVisible() && !Object.keys(selectionGids).length) {
        applySelectionHighlights(lastSelection);
      }
    }, 100);
  }

  window.addEventListener('resize', redrawSelectionSoon);
  if (window.IntersectionObserver) {
    new IntersectionObserver(function (entries) {
      if (entries.some(function (e) { return e.isIntersecting; })) {
        redrawSelectionSoon();
      }
    }).observe(svg);
  }

  if (trackData) {
    ['x', 'y'].forEach(function (axis) {
      var entries = trackData[axis] || [];
      var byGroup = {};
      entries.forEach(function (e) {
        (byGroup[e.group] = byGroup[e.group] || []).push(e);
      });
      entries.forEach(function (entry) {
        var el = document.getElementById(entry.gid);
        if (!el) return;
        el.classList.add('rd-track-feature');
        el.addEventListener('click', function (evt) {
          if (consumeDragClick()) {
            evt.stopPropagation();
            return;
          }
          evt.stopPropagation();
          var group = evt.shiftKey ? byGroup[entry.group] : [entry];
          toggleBands(
            group.map(function (e) {
              return e.gid;
            }),
            group,
            axis
          );
          showFeatureDetail(entry, el);
        });
      });
    });
  }

  // ---------------------------------------------------------------------
  // 6. Embedded display options (line width, min match length)
  // ---------------------------------------------------------------------
  //
  // When the report is embedded by the rusty-dot app it is rendered with
  // min_length=0 and the default line width; the app then drives both
  // options client-side via 'rd-display-opts' messages, so a cosmetic
  // change never re-renders matplotlib.  Line width is one injected CSS
  // rule (!important beats the per-path inline styles); min length toggles
  // a hiding class on each segment path and its .rd-hit clone using the
  // query-side spans registered above.  Standalone reports simply never
  // receive these messages.

  var displayStyle = document.createElement('style');
  document.head.appendChild(displayStyle);
  var currentMinLength = 0;

  function applyDisplayOpts(opts) {
    if (typeof opts.dot_size === 'number' && opts.dot_size > 0) {
      displayStyle.textContent =
        'g[id^="rd-matches-"] > path:not(.rd-hit) { stroke-width: ' +
        opts.dot_size +
        'px !important; }';
    }
    if (
      typeof opts.min_length === 'number' &&
      opts.min_length >= 0 &&
      opts.min_length !== currentMinLength
    ) {
      currentMinLength = opts.min_length;
      segmentRegistry.forEach(function (seg) {
        var hide = seg.qlen < currentMinLength;
        seg.el.classList.toggle('rd-len-hidden', hide);
        seg.hit.classList.toggle('rd-len-hidden', hide);
      });
    }
  }

  /* Re-apply a selection the app held across an iframe rebuild (an
   * overview <-> drill-down swap): resolve the view-independent match
   * identities against this report's registry and the stored annotation
   * feature against the diagonal squares and track entries.  Anything the
   * current view does not display simply stays stored in the app. */
  function applyRestoredSelection(msg) {
    var matches = Array.isArray(msg.matches) ? msg.matches : [];
    if (matches.length) {
      var byIdentity = {};
      segmentRegistry.forEach(function (e) {
        var k = [e.q, e.t, e.layer, e.qs, e.qe, e.ts, e.te].join(' ');
        (byIdentity[k] = byIdentity[k] || []).push(e);
      });
      var entries = [];
      matches.forEach(function (m) {
        var found = byIdentity[[m.q, m.t, m.layer, m.qs, m.qe, m.ts, m.te].join(' ')];
        if (found) entries.push.apply(entries, found);
      });
      if (entries.length) setMatchSelection(entries);
    }
    var f = msg.feature;
    if (f && typeof f === 'object') {
      var el = findFeatureElement(f);
      if (el) {
        if (selectedFeatureEl) {
          selectedFeatureEl.classList.remove('rd-selected-match');
        }
        selectedFeatureEl = el;
        el.classList.add('rd-selected-match');
        notifySelectionState();
      }
    }
  }

  function sameFeature(a, b) {
    return (
      a && b &&
      a.seqname === b.seqname &&
      a.type === b.type &&
      a.start === b.start &&
      a.end === b.end &&
      a.strand === b.strand
    );
  }

  function findFeatureElement(f) {
    for (var i = 0; i < annotRegistry.length; i++) {
      if (sameFeature(annotRegistry[i].feat, f)) return annotRegistry[i].el;
    }
    if (trackData) {
      var axes = ['x', 'y'];
      for (var a = 0; a < axes.length; a++) {
        var entries = trackData[axes[a]] || [];
        for (var j = 0; j < entries.length; j++) {
          if (sameFeature(entries[j], f)) {
            var el = document.getElementById(entries[j].gid);
            if (el) return el;
          }
        }
      }
    }
    return null;
  }

  window.addEventListener('message', function (ev) {
    var msg = ev && ev.data;
    if (!msg) return;
    if (msg.type === 'rd-display-opts') {
      applyDisplayOpts(msg);
      return;
    }
    if (msg.type === 'rd-highlight-features') {
      applySelectionHighlights(msg.features);
      return;
    }
    if (msg.type === 'rd-fullscreen') {
      // Embedding app entered/left fullscreen: scale the figure to fill
      // the iframe viewport (CSS keys off this class).
      document.body.classList.toggle('rd-fs', !!msg.on);
      return;
    }
    if (msg.type === 'rd-restore-selection') {
      applyRestoredSelection(msg);
      return;
    }
    if (msg.type === 'rd-clear-selection') {
      // The app caught Escape while something was selected: clear here,
      // exactly as an in-report Escape would.
      cancelDrag();
      resetView();
      clearAllSelection();
      clearBands();
      return;
    }
    if (msg.type === 'rd-seq-response') {
      // Sequences for the last clicked segment; ignore responses for
      // anything but the request currently in flight.
      var key = msg.row + ':' + msg.col + ':' + msg.layer + ':' + msg.idx;
      if (key !== pendingSeqKey) return;
      pendingSeqKey = null;
      if (typeof msg.text === 'string' && msg.text) {
        detailSeq.textContent = msg.text;
        // Gapped alignments need column-preserving layout; plain
        // query/target previews soft-wrap instead.
        detailSeq.classList.toggle('rd-aligned', !!msg.aligned);
        detailSeq.hidden = false;
        // Sequences exist server-side; they are fetched (and cached) on
        // the first copy press rather than shipped with the preview.
        setCopyState(!!msg.copy, !!msg.copy, null, null);
      } else {
        detailSeq.textContent = '';
        detailSeq.hidden = true;
      }
      return;
    }
    if (msg.type === 'rd-copy-response') {
      handleCopyResponse(msg);
    }
  });

  // Tell the embedding app the report is wired (the iframe element is
  // replaced on every re-render, so the app re-sends the current options
  // in response to this ping).
  if (window.parent && window.parent !== window) {
    window.parent.postMessage({ type: 'rd-report-ready' }, '*');
  }
})();
