// biowasm aligner bridge for the rusty-dot browser app.
//
// Message protocol (Shiny custom messages from Python):
//   'rd_mount_fasta'  {dataset_id, text}
//       Store a plain-FASTA payload under a stable id (the upload's content
//       digest).  Python sends each dataset once per session, so re-runs and
//       tool switches never re-copy multi-megabyte genomes across the
//       Pyodide/JS boundary again.
//   'rd_run_aligner'  {tool, args, query_id, target_id, request_id}
//       Run a tool against two previously mounted datasets.
//   'rd_cancel_aligner' {request_id}
//       Cancel the in-flight/queued run with that id: queued runs are
//       dropped before executing, and a best-effort attempt is made to
//       terminate the tool's WebWorker so the wasm computation actually
//       stops (falling back to ignore-the-result if the worker handle is
//       not reachable through Aioli's public surface).
//   'rd_extend_aligner' {request_id, ms}
//       Re-arm the run watchdog for that id, giving the tool another
//       `ms` milliseconds before the app is prompted again.
//
// Results return via
//   Shiny.setInputValue('aligner_result',
//       {request_id, tool, output, error, cancelled, stderr, cmd},
//       {priority: 'event'})
// On failure `output` is null and `error` carries a human-readable message.
// `stderr` (the tool's log) and `cmd` (the exact command line) are always
// included so the app can show a run log.
//
// A long run notifies rather than fails:
//   Shiny.setInputValue('aligner_timeout',
//       {request_id, tool, waited_ms, total_ms}, {priority: 'event'})
// The run keeps going; the app decides whether to extend the watchdog
// ('rd_extend_aligner') or cancel it ('rd_cancel_aligner').
//
// Aioli runs in its own WebWorker, so it never contends with the Pyodide
// runtime that executes the Shiny app itself.

(function () {
  'use strict';

  var AIOLI_URL = 'https://biowasm.com/cdn/v3/aioli/3.2.3/aioli.js';

  // tool key -> {module: biowasm CDN module spec, program: CLI program
  // name, outputFile: file to read back instead of stdout (or null)}
  var TOOLS = {
    minimap2: { module: 'minimap2/2.22', program: 'minimap2', outputFile: null },
    nucmer: {
      // Object spec: Aioli's 3-part string form is tool/PROGRAM/version,
      // so the object form is used to keep the fields unambiguous.
      module: { tool: 'mummer4', version: '4.0.0rc1', program: 'nucmer' },
      program: 'nucmer',
      outputFile: 'out.delta',
    },
  };

  // Python owns these names (core/align.py TARGET_FILENAME / QUERY_FILENAME);
  // build_tool_args embeds them in the CLI args, and we mount under the same
  // names so the args resolve inside the Aioli filesystem.
  var QUERY_FILENAME = 'query.fa';
  var TARGET_FILENAME = 'target.fa';

  var aioliScriptPromise = null; // loads aioli.js once
  var cliPromises = {}; // tool key -> Promise<Aioli instance>
  var runQueue = Promise.resolve(); // serialises runs (shared out.delta etc.)
  var datasets = {}; // dataset_id -> FASTA text (session-lived)
  // Which datasets are mounted on each tool's CLI instance; torn down with it.
  var mountedOn = {}; // tool key -> 'queryId|targetId'
  var cancelledIds = {}; // request_id -> true (queued runs drop themselves)
  var cancelResolvers = {}; // request_id -> resolve fn unblocking the queue
  var currentRun = null; // {requestId, tool} while a tool is executing

  // A hung init (bad CDN asset, unreachable biowasm) must reject rather
  // than block the run queue forever -- there is no useful "wait longer"
  // story for a fetch that is never going to arrive.
  var INIT_TIMEOUT_MS = 180000;

  // Runs are different: real assemblies can legitimately keep a wasm
  // aligner busy for many minutes (minimap2 -c especially), and the
  // WebWorker keeps computing regardless of what this side does with the
  // promise.  So the per-tool budget is a *watchdog*: when it expires the
  // app is told, the run continues, and the app either extends the budget
  // or cancels for real.
  var RUN_TIMEOUT_MS = {
    minimap2: 300000,
    nucmer: 600000,
  };

  // request_id -> {timer, tool, startedAt, waited}
  var watchdogs = {};

  function armWatchdog(requestId, tool, ms) {
    if (!requestId) return; // pre-warm runs are not watched
    var w = watchdogs[requestId];
    if (!w) {
      w = watchdogs[requestId] = { timer: null, tool: tool, startedAt: Date.now() };
    }
    clearTimeout(w.timer);
    w.timer = setTimeout(function () {
      // Notify only.  The entry is kept so 'rd_extend_aligner' can re-arm
      // it; the run itself is untouched and still racing to finish.
      w.timer = null;
      if (!window.Shiny || !window.Shiny.setInputValue) return;
      window.Shiny.setInputValue(
        'aligner_timeout',
        {
          request_id: requestId,
          tool: w.tool,
          waited_ms: ms,
          total_ms: Date.now() - w.startedAt,
        },
        { priority: 'event' }
      );
    }, ms);
  }

  function extendWatchdog(msg) {
    if (!msg || !msg.request_id) return;
    var w = watchdogs[msg.request_id];
    if (!w) return; // already finished or cancelled
    armWatchdog(msg.request_id, w.tool, msg.ms || RUN_TIMEOUT_MS[w.tool] || 300000);
  }

  function clearWatchdog(requestId) {
    var w = watchdogs[requestId];
    if (!w) return;
    clearTimeout(w.timer);
    delete watchdogs[requestId];
  }

  function withTimeout(promise, what, ms) {
    return Promise.race([
      promise,
      new Promise(function (resolve, reject) {
        setTimeout(function () {
          reject(new Error(what + ' timed out after ' + ms / 1000 + 's'));
        }, ms);
      }),
    ]);
  }

  function loadAioliScript() {
    if (window.Aioli) return Promise.resolve();
    if (!aioliScriptPromise) {
      aioliScriptPromise = new Promise(function (resolve, reject) {
        var s = document.createElement('script');
        s.src = AIOLI_URL;
        s.onload = function () {
          resolve();
        };
        s.onerror = function () {
          aioliScriptPromise = null; // allow retry on the next run
          reject(
            new Error(
              'Could not load Aioli from the biowasm CDN (' +
                AIOLI_URL +
                ') — check your network connection.'
            )
          );
        };
        document.head.appendChild(s);
      });
    }
    return aioliScriptPromise;
  }

  function getCli(tool) {
    if (!cliPromises[tool]) {
      mountedOn[tool] = null; // fresh instance, nothing mounted yet
      cliPromises[tool] = withTimeout(
        loadAioliScript().then(function () {
          // printInterleaved:false keeps stdout separate from stderr —
          // minimap2 writes its PAF output to stdout but also logs
          // progress to stderr, and the default interleaved string is
          // unparseable.
          return new window.Aioli([TOOLS[tool].module], {
            printInterleaved: false,
          });
        }),
        'Initialising ' + tool + ' from biowasm',
        INIT_TIMEOUT_MS
      ).catch(function (err) {
          cliPromises[tool] = null; // allow retry on the next run
          throw new Error(
            'Could not initialise ' + tool + ' from biowasm: ' + errText(err)
          );
        });
    }
    return cliPromises[tool];
  }

  // Best-effort worker termination.  Aioli does not officially expose its
  // WebWorker; probe the common handles and fall back to just abandoning
  // the instance (its result is ignored and the next run re-initialises).
  function teardownTool(tool) {
    var pending = cliPromises[tool];
    cliPromises[tool] = null;
    mountedOn[tool] = null;
    if (!pending) return;
    pending.then(
      function (cli) {
        try {
          var worker =
            (cli && cli.worker) ||
            (cli && cli._worker) ||
            (cli && cli.config && cli.config.worker);
          if (worker && typeof worker.terminate === 'function') {
            worker.terminate();
            return;
          }
        } catch (err) {
          // fall through to the console note below
        }
        // eslint-disable-next-line no-console
        console.warn(
          'rusty-dot: could not reach the ' +
            tool +
            ' worker to terminate it; its result will be ignored but the ' +
            'computation may continue until it finishes or times out.'
        );
      },
      function () {
        // init already failed; nothing to terminate
      }
    );
  }

  function errText(err) {
    if (err === undefined || err === null) return 'unknown error';
    if (err instanceof Error) return err.message;
    return String(err);
  }

  async function runAligner(msg) {
    var spec = TOOLS[msg.tool];
    if (!spec) throw new Error('Unknown aligner tool: ' + msg.tool);
    var queryText = datasets[msg.query_id];
    var targetText = datasets[msg.target_id];
    if (typeof queryText !== 'string' || typeof targetText !== 'string') {
      throw new Error(
        'Assembly data not found in the browser session — please press ' +
          'Run again.'
      );
    }
    if (!window.Aioli) sendProgress(msg.request_id, msg.tool, 'loading-aioli');
    if (!cliPromises[msg.tool]) {
      sendProgress(msg.request_id, msg.tool, 'initialising-tool');
    }
    var CLI = await getCli(msg.tool);
    if (spec.outputFile) {
      // Remove any output left by a previous run on this cached CLI
      // instance: if the tool fails, a stale file must not be read back
      // and reported as this run's result.
      try {
        await CLI.fs.unlink(spec.outputFile);
      } catch (err) {
        // Not present (first run) — fine.
      }
    }
    var mountKey = msg.query_id + '|' + msg.target_id;
    if (mountedOn[msg.tool] !== mountKey) {
      sendProgress(msg.request_id, msg.tool, 'mounting-data');
      await CLI.mount([
        { name: QUERY_FILENAME, data: queryText },
        { name: TARGET_FILENAME, data: targetText },
      ]);
      mountedOn[msg.tool] = mountKey;
    }
    var cmd = [spec.program].concat(msg.args || []).join(' ');
    sendProgress(msg.request_id, msg.tool, 'aligning');
    var res;
    currentRun = { requestId: msg.request_id, tool: msg.tool };
    armWatchdog(msg.request_id, msg.tool, RUN_TIMEOUT_MS[msg.tool] || 300000);
    try {
      // Deliberately not raced against a timer: the run ends only by
      // completing, erroring, or being cancelled.  Rejecting here would
      // report a failure while the worker kept burning CPU on a run that
      // was about to succeed.
      res = await CLI.exec(cmd);
    } finally {
      clearWatchdog(msg.request_id);
      currentRun = null;
    }
    sendProgress(msg.request_id, msg.tool, 'reading-output');
    var stdout = res && typeof res === 'object' ? res.stdout : res;
    var stderr = res && typeof res === 'object' ? res.stderr : '';
    if (spec.outputFile) {
      var fileText;
      try {
        fileText = await CLI.fs.readFile(spec.outputFile, { encoding: 'utf8' });
      } catch (err) {
        throw new Error(
          msg.tool +
            ' did not produce ' +
            spec.outputFile +
            ' (tool output: ' +
            String(stderr || stdout).slice(0, 500) +
            ')'
        );
      }
      return { output: fileText, stderr: String(stderr || ''), cmd: cmd };
    }
    return {
      output: typeof stdout === 'string' ? stdout : String(stdout),
      stderr: String(stderr || ''),
      cmd: cmd,
    };
  }

  function sendResult(payload) {
    window.Shiny.setInputValue('aligner_result', payload, { priority: 'event' });
  }

  // Stage updates for the in-flight run ('loading-aioli',
  // 'initialising-tool', 'mounting-data', 'aligning', 'reading-output').
  // Pre-warm runs have no request id and stay silent.
  function sendProgress(requestId, tool, stage) {
    if (!requestId || !window.Shiny || !window.Shiny.setInputValue) return;
    window.Shiny.setInputValue(
      'aligner_progress',
      { request_id: requestId, tool: tool, stage: stage },
      { priority: 'event' }
    );
  }

  function handleMount(msg) {
    if (msg && msg.dataset_id && typeof msg.text === 'string') {
      datasets[msg.dataset_id] = msg.text;
    }
  }

  function handleCancel(msg) {
    if (!msg || !msg.request_id) return;
    cancelledIds[msg.request_id] = true;
    clearWatchdog(msg.request_id);
    if (currentRun && currentRun.requestId === msg.request_id) {
      // The tool is executing right now: kill its worker so the CPU work
      // actually stops (best-effort; see teardownTool).
      teardownTool(currentRun.tool);
      currentRun = null;
    }
    // Unblock the queue immediately — a terminated worker's exec promise
    // never settles, and the next run must not wait behind it.
    var release = cancelResolvers[msg.request_id];
    if (release) release({ cancelled: true });
  }

  function handleMessage(msg) {
    // Serialise runs: concurrent nucmer runs would race on out.delta, and
    // biowasm downloads are friendlier one at a time.
    runQueue = runQueue
      .then(function () {
        if (cancelledIds[msg.request_id]) {
          // Cancelled while queued: drop before doing any work, so a
          // superseded run cannot burn CPU behind the one the user wants.
          delete cancelledIds[msg.request_id];
          return { cancelled: true };
        }
        var cancelPromise = new Promise(function (resolve) {
          cancelResolvers[msg.request_id] = resolve;
        });
        var runPromise = runAligner(msg);
        // The losing promise of the race may reject later (e.g. the
        // timeout of a terminated run); swallow it so it never surfaces
        // as an unhandled rejection.
        runPromise.catch(function () {});
        return Promise.race([runPromise, cancelPromise]).finally(function () {
          delete cancelResolvers[msg.request_id];
        });
      })
      .then(
        function (res) {
          if (res && res.cancelled) {
            sendResult({
              request_id: msg.request_id,
              tool: msg.tool,
              output: null,
              error: null,
              cancelled: true,
              stderr: '',
              cmd: '',
            });
            return;
          }
          sendResult({
            request_id: msg.request_id,
            tool: msg.tool,
            output: res.output,
            error: null,
            cancelled: false,
            stderr: res.stderr,
            cmd: res.cmd,
          });
        },
        function (err) {
          var wasCancelled = !!cancelledIds[msg.request_id];
          delete cancelledIds[msg.request_id];
          sendResult({
            request_id: msg.request_id,
            tool: msg.tool,
            output: null,
            error: wasCancelled ? null : errText(err),
            cancelled: wasCancelled,
            stderr: '',
            cmd: '',
          });
        }
      );
  }

  function register() {
    // Inside the shinylive iframe `Shiny` is the in-frame global; it may
    // not exist yet when this script runs, so poll briefly.
    if (window.Shiny && window.Shiny.addCustomMessageHandler) {
      window.Shiny.addCustomMessageHandler('rd_mount_fasta', handleMount);
      window.Shiny.addCustomMessageHandler('rd_run_aligner', handleMessage);
      window.Shiny.addCustomMessageHandler('rd_cancel_aligner', handleCancel);
      window.Shiny.addCustomMessageHandler('rd_extend_aligner', extendWatchdog);
    } else {
      setTimeout(register, 100);
    }
  }

  // Pre-warm: start downloading a biowasm tool as soon as the user selects
  // it in the method menu, so the (one-time, ~seconds) CDN fetch overlaps
  // with choosing files and options instead of delaying the first Run.
  // getCli() caches the instance, so Run reuses whatever this started.
  document.addEventListener('change', function (ev) {
    var el = ev.target;
    if (!el || el.id !== 'method') return;
    if (TOOLS[el.value]) {
      getCli(el.value).catch(function () {
        // Ignore here: cliPromises resets itself on failure, and the real
        // run reports CDN errors to the user properly.
      });
    }
  });

  register();
})();
