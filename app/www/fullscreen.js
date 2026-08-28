// Fullscreen toggle for the plot area.  The button lives inside the
// server-rendered plot toolbar (re-created on every plot_area render), so
// the wiring is delegated and the on/off state is a class on <html> --
// both survive re-renders, and app.css derives the overlay layout, the
// hidden chrome (hints, tabs, memory note) and the button's icon from
// that class alone.
//
// Entering also requests native browser fullscreen (hiding the URL bar
// and tabs).  When that is unavailable -- e.g. the shinylive viewer
// iframe lacks allowfullscreen -- the request fails silently and the
// CSS overlay still fills the browser window.
(function () {
  'use strict';

  function fsActive() {
    return document.documentElement.classList.contains('rd-fullscreen');
  }

  /* Toggle the state class and tell the bridge, which relays it into the
   * report iframe so the figure can scale to the viewport. */
  function applyState(on) {
    document.documentElement.classList.toggle('rd-fullscreen', on);
    document.dispatchEvent(
      new CustomEvent('rd-fs-change', { detail: { on: on } })
    );
  }

  /* The drill-down tab strip is hidden in fullscreen, so make sure the
   * Plot pane (always the first tab) is the visible one before entering. */
  function activatePlotTab() {
    var link = document.querySelector('.rd-plot-area .nav-tabs .nav-link');
    if (link && !link.classList.contains('active')) link.click();
  }

  function setFullscreen(on) {
    applyState(on);
    if (on) {
      activatePlotTab();
      if (!document.fullscreenElement && document.documentElement.requestFullscreen) {
        document.documentElement.requestFullscreen().catch(function () {});
      }
    } else if (document.fullscreenElement && document.exitFullscreen) {
      document.exitFullscreen().catch(function () {});
    }
  }

  document.addEventListener('click', function (ev) {
    var btn = ev.target.closest && ev.target.closest('.rd-fs-btn');
    if (!btn) return;
    setFullscreen(!fsActive());
  });

  // Escape with something selected clears the selection instead of
  // exiting; only an Escape with nothing selected leaves fullscreen.
  // The report announces its selection state ('rd-selection-state') and
  // the annotations table announces its rows ('rd-ft-selection').
  var reportHasSelection = false;
  var tableSelectionCount = 0;

  function anySelection() {
    return reportHasSelection || tableSelectionCount > 0;
  }

  function clearReportSelection() {
    var frame = document.querySelector('iframe.rd-report-frame');
    if (frame && frame.contentWindow) {
      frame.contentWindow.postMessage({ type: 'rd-clear-selection' }, '*');
    }
  }

  window.addEventListener('message', function (ev) {
    var msg = ev && ev.data;
    if (!msg) return;
    if (msg.type === 'rd-selection-state') {
      reportHasSelection = !!msg.any;
    } else if (msg.type === 'rd-report-ready') {
      // A fresh iframe starts with nothing selected.
      reportHasSelection = false;
    } else if (msg.type === 'rd-esc' && fsActive()) {
      // Escape inside the report with nothing left to clear.
      setFullscreen(false);
    }
  });

  document.addEventListener('rd-ft-selection', function (ev) {
    var features = ev && ev.detail && ev.detail.features;
    tableSelectionCount = Array.isArray(features) ? features.length : 0;
  });

  // Native fullscreen exits on its own Escape (or browser UI).  If that
  // Escape had a selection to clear, keep the overlay in place (the
  // browser's exit itself cannot be prevented) so the next Escape is the
  // one that actually leaves fullscreen; otherwise drop the overlay so
  // the two states never diverge.
  document.addEventListener('fullscreenchange', function () {
    if (document.fullscreenElement || !fsActive()) return;
    if (anySelection()) {
      clearReportSelection();
    } else {
      applyState(false);
    }
  });

  // Capture phase so this sees the selection state before the table's own
  // Escape handler clears it.  Only reaches this document when focus is
  // outside the sandboxed report iframe; Escape inside the iframe is
  // handled by the report itself (clearing first, then 'rd-esc').
  document.addEventListener(
    'keydown',
    function (ev) {
      if (ev.key !== 'Escape' || !fsActive()) return;
      if (anySelection()) {
        clearReportSelection();
        // The annotations table clears its own rows on this same keydown.
      } else {
        setFullscreen(false);
      }
    },
    true
  );
})();
