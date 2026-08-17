// biowasm aligner bridge for the rusty-dot browser app.
//
// Listens for the Shiny custom message 'rd_run_aligner' with payload
//   {tool, args, query_fasta, target_fasta, request_id}
// loads Aioli lazily from the biowasm CDN on first use, mounts the two
// FASTA strings, runs the tool in its WebWorker, and returns the text
// output via
//   Shiny.setInputValue('aligner_result',
//                       {request_id, tool, output, error},
//                       {priority: 'event'})
// On any failure (CDN unreachable, tool error, empty output) `output` is
// null and `error` carries a human-readable message.
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
    lastz: { module: 'lastz/1.04.52', program: 'lastz', outputFile: null },
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

  // A hung init/exec (bad CDN asset, wedged WebWorker) must reject rather
  // than block the run queue forever.
  var RUN_TIMEOUT_MS = 180000;

  function withTimeout(promise, what) {
    return Promise.race([
      promise,
      new Promise(function (resolve, reject) {
        setTimeout(function () {
          reject(
            new Error(what + ' timed out after ' + RUN_TIMEOUT_MS / 1000 + 's')
          );
        }, RUN_TIMEOUT_MS);
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
      cliPromises[tool] = withTimeout(
        loadAioliScript().then(function () {
          // printInterleaved:false keeps stdout separate from stderr —
          // minimap2/lastz write their PAF/general output to stdout but
          // also log progress to stderr, and the default interleaved
          // string is unparseable.
          return new window.Aioli([TOOLS[tool].module], {
            printInterleaved: false,
          });
        }),
        'Initialising ' + tool + ' from biowasm'
      ).catch(function (err) {
          cliPromises[tool] = null; // allow retry on the next run
          throw new Error(
            'Could not initialise ' + tool + ' from biowasm: ' + errText(err)
          );
        });
    }
    return cliPromises[tool];
  }

  function errText(err) {
    if (err === undefined || err === null) return 'unknown error';
    if (err instanceof Error) return err.message;
    return String(err);
  }

  async function runAligner(msg) {
    var spec = TOOLS[msg.tool];
    if (!spec) throw new Error('Unknown aligner tool: ' + msg.tool);
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
    await CLI.mount([
      { name: QUERY_FILENAME, data: msg.query_fasta },
      { name: TARGET_FILENAME, data: msg.target_fasta },
    ]);
    var cmd = [spec.program].concat(msg.args || []).join(' ');
    var res = await withTimeout(CLI.exec(cmd), msg.tool + ' run');
    var stdout = res && typeof res === 'object' ? res.stdout : res;
    var stderr = res && typeof res === 'object' ? res.stderr : '';
    if (spec.outputFile) {
      try {
        return await CLI.fs.readFile(spec.outputFile, { encoding: 'utf8' });
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
    }
    return typeof stdout === 'string' ? stdout : String(stdout);
  }

  function sendResult(payload) {
    window.Shiny.setInputValue('aligner_result', payload, { priority: 'event' });
  }

  function handleMessage(msg) {
    // Serialise runs: concurrent nucmer runs would race on out.delta, and
    // biowasm downloads are friendlier one at a time.
    runQueue = runQueue
      .then(function () {
        return runAligner(msg);
      })
      .then(
        function (output) {
          sendResult({
            request_id: msg.request_id,
            tool: msg.tool,
            output: output,
            error: null,
          });
        },
        function (err) {
          sendResult({
            request_id: msg.request_id,
            tool: msg.tool,
            output: null,
            error: errText(err),
          });
        }
      );
  }

  function register() {
    // Inside the shinylive iframe `Shiny` is the in-frame global; it may
    // not exist yet when this script runs, so poll briefly.
    if (window.Shiny && window.Shiny.addCustomMessageHandler) {
      window.Shiny.addCustomMessageHandler('rd_run_aligner', handleMessage);
    } else {
      setTimeout(register, 100);
    }
  }

  register();
})();
