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
 *                             others; click again (or Esc) resets.
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

  var panelGroups = Array.prototype.slice.call(
    svg.querySelectorAll('g[id^="rd-panel-"]')
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
    if (selectedPanel === panel) {
      resetView();
      return;
    }
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
   * the one visible plot, so leave clicks alone there. */
  if (panelGroups.length > 1) {
    panelGroups.forEach(function (panel) {
      panel.addEventListener('click', function (evt) {
        // A completed drag-zoom releases a click too; swallow that one.
        if (consumeDragClick()) return;
        // Match clicks are handled (and stopped) by the match handler below.
        evt.stopPropagation();
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
    var panel =
      evt.target && evt.target.closest
        ? evt.target.closest('g[id^="rd-panel-"]')
        : null;
    if (!panel) return;
    drag = {
      cx0: evt.clientX,
      cy0: evt.clientY,
      p0: eventPoint(evt),
      panelBox: bboxInRootSpace(panel),
      active: false,
      rect: null,
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
      drag.rect.setAttribute('class', 'rd-drag-rect');
      svg.appendChild(drag.rect);
    }
    evt.preventDefault();
    var r = dragRegion(evt);
    drag.rect.setAttribute('x', r.x);
    drag.rect.setAttribute('y', r.y);
    drag.rect.setAttribute('width', r.w);
    drag.rect.setAttribute('height', r.h);
  });

  document.addEventListener('mouseup', function (evt) {
    if (!drag) return;
    var wasActive = drag.active;
    var region = wasActive ? dragRegion(evt) : null;
    cancelDrag();
    if (!wasActive) return; // plain click: let the click handlers run
    suppressNextClick = true;
    if (region.w > 0 && region.h > 0) {
      setViewBox(fitAspect(region));
    }
  });

  document.addEventListener('keydown', function (evt) {
    if (evt.key === 'Escape') {
      cancelDrag();
      resetView();
      hideDetail();
      clearBands();
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
  var selectedMatch = null;
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

  function hideDetail() {
    detail.hidden = true;
    pendingSeqKey = null;
    currentSeg = null;
    setCopyState(false, false, null, null);
    if (selectedMatch) {
      selectedMatch.classList.remove('rd-selected-match');
      selectedMatch = null;
    }
  }

  document
    .getElementById('rd-detail-close')
    .addEventListener('click', hideDetail);

  function showMatchDetail(panelGid, layer, idx, el) {
    var panel = payload.panels[panelGid];
    if (!panel) return;
    var seg = panel.segments[layer][idx];
    if (!seg) return;

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
      var pm = panelGid.match(/^rd-panel-(\d+)-(\d+)$/);
      if (pm) {
        var row = parseInt(pm[1], 10);
        var col = parseInt(pm[2], 10);
        pendingSeqKey = row + ':' + col + ':' + layer + ':' + idx;
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
    } else {
      detailSeq.textContent = '';
      detailSeq.hidden = true;
    }
    detail.hidden = false;

    if (selectedMatch) selectedMatch.classList.remove('rd-selected-match');
    selectedMatch = el;
    el.classList.add('rd-selected-match');
  }

  /* Wire up every match group: index its drawable children in document
   * order (the serialisation contract) and overlay each with an invisible
   * wide-stroke clone so hairline matches are easy to click.  Each segment
   * is also registered with its query-side length so the embedding app can
   * filter by minimum match length client-side (behaviour 6). */
  var segmentRegistry = []; // {el, hit, qlen}
  var matchGroups = Array.prototype.slice.call(
    svg.querySelectorAll('g[id^="rd-matches-"]')
  );
  matchGroups.forEach(function (group) {
    // gid: rd-matches-<row>-<col>-<layer>
    var m = group.id.match(/^rd-matches-(\d+)-(\d+)-(\w+)$/);
    if (!m) return;
    var panelGid = 'rd-panel-' + m[1] + '-' + m[2];
    var layer = m[3];
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
        showMatchDetail(panelGid, layer, idx, el);
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
        // Query-side span, matching the server-side min_length semantics.
        segmentRegistry.push({
          el: el,
          hit: hit,
          qlen: segs[idx][1] - segs[idx][0],
        });
      }
    });
  });

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

    if (selectedMatch) selectedMatch.classList.remove('rd-selected-match');
    selectedMatch = el;
    el.classList.add('rd-selected-match');
  }

  /* Wire up every diagonal-annotation group: the n-th drawable child of
   * 'rd-annot-<row>-<col>' corresponds to panels[gid].annotations[n]
   * (serialisation contract — patch draw order equals SVG child order). */
  var annotGroups = Array.prototype.slice.call(
    svg.querySelectorAll('g[id^="rd-annot-"]')
  );
  annotGroups.forEach(function (group) {
    var m = group.id.match(/^rd-annot-(\d+)-(\d+)$/);
    if (!m) return;
    var panelGid = 'rd-panel-' + m[1] + '-' + m[2];

    var children = Array.prototype.slice.call(
      group.querySelectorAll('path, use')
    );
    children.forEach(function (el, idx) {
      el.classList.add('rd-annot');
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
  var bandLayer = null;

  function panelBackground() {
    return svg.querySelector('[id="rd-panel-0-0-bg"]');
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
    el.classList.add('rd-track-active');
  }

  function removeBand(gid) {
    var rect = activeBands[gid];
    if (rect && rect.parentNode) rect.parentNode.removeChild(rect);
    delete activeBands[gid];
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

  window.addEventListener('message', function (ev) {
    var msg = ev && ev.data;
    if (!msg) return;
    if (msg.type === 'rd-display-opts') {
      applyDisplayOpts(msg);
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
