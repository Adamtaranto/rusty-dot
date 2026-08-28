// Mirror ui.Progress messages into the header's task-status slot.
//
// The server keeps using ui.Progress unchanged -- its messages flush
// synchronously mid-computation, which a reactive output cannot do -- but
// the popup card it draws (a .shiny-progress-notification inside the
// notification panel) duplicates the busy pill, so app.css hides it and
// this observer copies its text into #rd-task-status at the right end of
// the header bar instead.  Ordinary ui.notification_show toasts have no
// .shiny-progress-notification node and are untouched.

(function () {
  'use strict';

  var scheduled = false;

  function sync() {
    scheduled = false;
    var slot = document.getElementById('rd-task-status');
    if (!slot) return;
    var cards = document.querySelectorAll('.shiny-progress-notification');
    if (!cards.length) {
      if (slot.textContent) slot.textContent = '';
      return;
    }
    // Several progress bars can stack (nested renders); the newest one is
    // the task actually running now.
    var card = cards[cards.length - 1];
    var message = card.querySelector('.progress-message');
    var detail = card.querySelector('.progress-detail');
    var text = (
      (message ? message.textContent : '') +
      ' ' +
      (detail ? detail.textContent : '')
    ).trim();
    // The placeholder text before the first real .set() call.
    if (text === 'message') text = '';
    if (slot.textContent !== text) slot.textContent = text;
  }

  function schedule() {
    if (scheduled) return;
    scheduled = true;
    window.requestAnimationFrame(sync);
  }

  new MutationObserver(schedule).observe(document.documentElement, {
    childList: true,
    subtree: true,
    characterData: true,
  });
})();
