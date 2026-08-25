// Keep the sidebar's scroll position when a file is chosen in it.
//
// Picking a GFF sends the sidebar back to the top, losing the user's place
// right after an action whose result (the feature-type list) is further
// down the panel.
//
// Instrumenting the scroll container shows the jump has already happened by
// the time the input's `change` event fires — it comes out of the file
// picker interaction itself, not out of the dependent outputs re-rendering
// afterwards.  So the position has to be captured when the picker is
// opened, and re-asserted once the file comes back.
//
// The hold spans the upload and the re-renders it triggers, because the
// panel keeps growing as those outputs fill in and a single restore would
// be clamped against a height that is still settling.  Any real scroll
// gesture cancels it, so this can never fight the user for the scrollbar.

(function () {
  'use strict';

  // Long enough to cover the upload plus the dependent re-renders; short
  // enough that a stuck hold would be imperceptible.
  var HOLD_MS = 2000;

  /* Nearest scrollable ancestor, or null.  The scroller is bslib's
   * <aside class="sidebar"> on desktop but the inner .sidebar-content once
   * the layout collapses, so find it rather than hard-coding either. */
  function scrollerFor(el) {
    var node = el;
    while (node && node !== document.documentElement) {
      if (node.scrollHeight > node.clientHeight + 2) {
        var overflowY = window.getComputedStyle(node).overflowY;
        if (overflowY === 'auto' || overflowY === 'scroll') return node;
      }
      node = node.parentElement;
    }
    return null;
  }

  function sidebarFileInput(el) {
    if (!el || el.type !== 'file' || !el.closest) return null;
    return el.closest('.sidebar') ? el : null;
  }

  var mark = null; // {el, top} captured when the picker was opened
  var hold = null; // {el, top, until} active restore

  function release() {
    hold = null;
  }

  // Any genuine scroll gesture wins; the hold exists only to undo the
  // programmatic jump, never to pin the panel against the user.
  ['wheel', 'touchmove', 'keydown'].forEach(function (name) {
    document.addEventListener(name, release, { passive: true, capture: true });
  });

  function tick() {
    if (!hold) return;
    if (Date.now() > hold.until) {
      hold = null;
      return;
    }
    var el = hold.el;
    var want = Math.min(hold.top, el.scrollHeight - el.clientHeight);
    if (el.scrollTop !== want) el.scrollTop = want;
    requestAnimationFrame(tick);
  }

  // Opening the picker is the last moment the user's position is still on
  // screen, so record it here rather than on `change`.
  document.addEventListener(
    'pointerdown',
    function (ev) {
      var input = sidebarFileInput(ev.target);
      if (!input) return;
      var scroller = scrollerFor(input);
      if (scroller && scroller.scrollTop > 0) {
        mark = { el: scroller, top: scroller.scrollTop };
      }
    },
    true
  );

  document.addEventListener(
    'click',
    function (ev) {
      var input = sidebarFileInput(ev.target);
      if (!input) return;
      var scroller = scrollerFor(input);
      if (scroller && scroller.scrollTop > 0) {
        mark = { el: scroller, top: scroller.scrollTop };
      }
    },
    true
  );

  document.addEventListener(
    'change',
    function (ev) {
      var input = sidebarFileInput(ev.target);
      if (!input) return;
      var scroller = scrollerFor(input);
      if (!scroller) return;
      // Prefer the pre-picker position; fall back to the current one when
      // the picker never moved it (nothing to undo, but the re-renders
      // below can still clamp it).
      var top = mark && mark.el === scroller ? mark.top : scroller.scrollTop;
      mark = null;
      if (top <= 0) return;
      hold = { el: scroller, top: top, until: Date.now() + HOLD_MS };
      requestAnimationFrame(tick);
    },
    true
  );
})();
