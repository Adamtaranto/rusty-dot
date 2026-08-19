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

  panelGroups.forEach(function (panel) {
    panel.addEventListener('click', function (evt) {
      // A completed drag-zoom releases a click too; swallow that one.
      if (consumeDragClick()) return;
      // Match clicks are handled (and stopped) by the match handler below.
      evt.stopPropagation();
      selectPanel(panel);
    });
  });

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
    }
  });

  // ---------------------------------------------------------------------
  // 4. Match click -> detail bar
  // ---------------------------------------------------------------------

  var detail = document.getElementById('rd-detail');
  var detailCoords = document.getElementById('rd-detail-coords');
  var detailSeq = document.getElementById('rd-detail-seq');
  var selectedMatch = null;

  function hideDetail() {
    detail.hidden = true;
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
    if (seq) {
      detailSeq.textContent = seq;
      detailSeq.hidden = false;
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
   * wide-stroke clone so hairline matches are easy to click. */
  var matchGroups = Array.prototype.slice.call(
    svg.querySelectorAll('g[id^="rd-matches-"]')
  );
  matchGroups.forEach(function (group) {
    // gid: rd-matches-<row>-<col>-<layer>
    var m = group.id.match(/^rd-matches-(\d+)-(\d+)-(\w+)$/);
    if (!m) return;
    var panelGid = 'rd-panel-' + m[1] + '-' + m[2];
    var layer = m[3];

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
    });
  });

  // ---------------------------------------------------------------------
  // 5. Annotation feature click -> detail bar
  // ---------------------------------------------------------------------

  function showAnnotationDetail(panelGid, idx, el) {
    var panel = payload.panels[panelGid];
    if (!panel || !panel.annotations) return;
    var feat = panel.annotations[idx];
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
})();
