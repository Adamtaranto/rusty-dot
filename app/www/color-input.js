// Minimal Shiny input binding for native <input type="color"> pickers.
//
// Shiny for Python ships no colour input; this binding turns any
// <input type="color" class="rd-color-input" id="..."> into a regular
// Shiny input whose value is the '#rrggbb' hex string.  Change events are
// debounced by the 250 ms rate policy so dragging inside the OS colour
// dialog doesn't re-render the plot per intermediate colour.

(function () {
  'use strict';

  function register() {
    if (!window.Shiny || !window.Shiny.InputBinding) {
      setTimeout(register, 100);
      return;
    }
    var binding = new window.Shiny.InputBinding();
    Object.assign(binding, {
      find: function (scope) {
        return scope.querySelectorAll ?
          scope.querySelectorAll('input.rd-color-input') :
          [];
      },
      getId: function (el) {
        return el.id;
      },
      getValue: function (el) {
        return el.value;
      },
      setValue: function (el, value) {
        el.value = value;
      },
      subscribe: function (el, callback) {
        el.addEventListener('input', el._rdColorHandler = function () {
          callback(true); // true -> apply the debounce rate policy
        });
        el.addEventListener('change', el._rdColorChange = function () {
          callback(false);
        });
      },
      unsubscribe: function (el) {
        el.removeEventListener('input', el._rdColorHandler);
        el.removeEventListener('change', el._rdColorChange);
      },
      getRatePolicy: function () {
        return { policy: 'debounce', delay: 250 };
      },
    });
    window.Shiny.inputBindings.register(binding, 'rusty-dot.colorInput');
  }

  register();
})();
