#!/usr/bin/env node

import { createServer } from 'node:http';
import { spawn, spawnSync } from 'node:child_process';
import {
  existsSync,
  mkdirSync,
  readFileSync,
  rmSync,
  statSync,
  writeFileSync,
} from 'node:fs';
import path from 'node:path';
import process from 'node:process';
import { fileURLToPath } from 'node:url';

const rootDir = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const defaultManifest = 'tests/expected/molstar-reference/reference-fixtures.txt';
const defaultMolstarDir = 'artifacts/molstar';
const defaultOutDir = '/private/tmp/molfig-molstar-browser-reference';
const defaultFormats = ['obj', 'stl'];
const defaultChromePath = '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';

function printHelp() {
  console.log(`Usage: node scripts/molstar-browser-reference-convert.mjs [options]

Run the pinned Mol* exporter inside a real Chrome/WebGL page instead of the
headless-gl CommonJS runtime. This is intended for Web App parity probes,
especially cylinder-impostor geometry that requires GL_EXT_frag_depth.

Options:
  --dry-run, --check       Validate inputs, bundle readiness, and output plan only.
  --build-only             Build the browser harness bundle, then exit.
  --manifest <path>        Reference fixture manifest. Default: ${defaultManifest}
  --molstar-dir <path>     Pinned Mol* checkout/build directory. Default: ${defaultMolstarDir}
  --out-dir <path>         Output directory. Default: ${defaultOutDir}
  --formats <list>         Comma-separated targets: obj,stl,report. Default: ${defaultFormats.join(',')}
  --chrome <path>          Chrome/Chromium executable. Default: ${defaultChromePath}
  --headed                 Show the browser window. Default is headless Chrome.
  --keep-profile           Keep the temporary Chrome user-data-dir.
  --render-object-report   Print per-render-object draw-count diagnostics.
  --compare-references     Diff generated OBJ/MTL/STL bytes against contract references.
  --compare-molfig          Diff generated browser OBJ/STL against molfig exports.
  --molfig-diff <path|cargo>
                          molfig-diff command for --compare-molfig. Default: built binary or cargo.
  --debug-stl-facet <n[,n]>
                          Dump browser render-object raw/centered vertices for STL sparse facet(s).
  --timeout-ms <ms>        Browser export timeout per fixture. Default: 120000.
  --help, -h               Show this help.
`);
}

function parseArgs(argv) {
  const args = {
    dryRun: false,
    buildOnly: false,
    manifest: defaultManifest,
    molstarDir: defaultMolstarDir,
    outDir: defaultOutDir,
    formats: [...defaultFormats],
    formatsFromCli: false,
    chromePath: defaultChromePath,
    headed: false,
    keepProfile: false,
    renderObjectReport: false,
    compareReferences: false,
    compareMolfig: false,
    molfigDiff: undefined,
    debugStlFacets: [],
    timeoutMs: 120000,
    help: false,
  };

  for (let i = 0; i < argv.length; i++) {
    const arg = argv[i];
    if (arg === '--dry-run' || arg === '--check') {
      args.dryRun = true;
    } else if (arg === '--build-only') {
      args.buildOnly = true;
    } else if (arg === '--help' || arg === '-h') {
      args.help = true;
    } else if (arg === '--manifest') {
      args.manifest = requireValue(argv, ++i, arg);
    } else if (arg === '--molstar-dir') {
      args.molstarDir = requireValue(argv, ++i, arg);
    } else if (arg === '--out-dir') {
      args.outDir = requireValue(argv, ++i, arg);
    } else if (arg === '--formats') {
      args.formats = parseFormats(requireValue(argv, ++i, arg), arg);
      args.formatsFromCli = true;
    } else if (arg === '--chrome') {
      args.chromePath = requireValue(argv, ++i, arg);
    } else if (arg === '--headed') {
      args.headed = true;
    } else if (arg === '--keep-profile') {
      args.keepProfile = true;
    } else if (arg === '--render-object-report') {
      args.renderObjectReport = true;
    } else if (arg === '--compare-references') {
      args.compareReferences = true;
    } else if (arg === '--compare-molfig') {
      args.compareMolfig = true;
    } else if (arg === '--molfig-diff') {
      args.molfigDiff = requireValue(argv, ++i, arg);
    } else if (arg === '--debug-stl-facet') {
      args.debugStlFacets.push(...parseNonNegativeIntegers(requireValue(argv, ++i, arg), arg));
    } else if (arg === '--timeout-ms') {
      args.timeoutMs = requirePositiveInteger(requireValue(argv, ++i, arg), arg);
    } else {
      throw new Error(`Unknown option: ${arg}`);
    }
  }

  validateFormats(args.formats, '--formats');
  args.debugStlFacets = [...new Set(args.debugStlFacets)].sort((a, b) => a - b);
  return args;
}

function requireValue(argv, index, flag) {
  const value = argv[index];
  if (!value || value.startsWith('--')) throw new Error(`${flag} requires a value.`);
  return value;
}

function requirePositiveInteger(value, label) {
  const parsed = Number(value);
  if (!Number.isInteger(parsed) || parsed <= 0) throw new Error(`${label}: expected a positive integer`);
  return parsed;
}

function parseNonNegativeIntegers(value, label) {
  const values = String(value).split(',').map(part => part.trim()).filter(Boolean);
  if (values.length === 0) throw new Error(`${label}: expected at least one integer`);
  return values.map(part => {
    const parsed = Number(part);
    if (!Number.isInteger(parsed) || parsed < 0) throw new Error(`${label}: expected a non-negative integer, got ${part}`);
    return parsed;
  });
}

function parseFormats(value, label) {
  const formats = String(value).split(',').map(s => s.trim()).filter(Boolean);
  validateFormats(formats, label);
  return formats;
}

function validateFormats(formats, label) {
  if (formats.length === 0) throw new Error(`${label}: expected at least one format`);
  const seen = new Set();
  for (const format of formats) {
    if (format !== 'obj' && format !== 'stl' && format !== 'report') {
      throw new Error(`${label}: unsupported format '${format}'. Expected obj, stl, or report.`);
    }
    if (seen.has(format)) throw new Error(`${label}: duplicate format '${format}'`);
    seen.add(format);
  }
}

function resolveRepoPath(value) {
  return path.resolve(rootDir, value);
}

function resolveInputPath(value) {
  return path.isAbsolute(value) ? path.resolve(value) : resolveRepoPath(value);
}

function displayPath(value) {
  const absolute = path.resolve(value);
  const relative = path.relative(rootDir, absolute);
  return relative && !relative.startsWith('..') && !path.isAbsolute(relative) ? relative : absolute;
}

function materializeGeneratedFixture(contract, contractPath, absFixturePath) {
  const generator = String(contract.fixture_generator ?? '');
  if (!generator) return;
  const absGeneratorPath = resolveRepoPath(generator);
  if (!existsSync(absGeneratorPath) || !statSync(absGeneratorPath).isFile()) {
    throw new Error(`${contractPath}: fixture generator not found: ${generator}`);
  }
  mkdirSync(path.dirname(absFixturePath), { recursive: true });
  const result = spawnSync(process.execPath, [absGeneratorPath, absFixturePath], {
    cwd: rootDir,
    encoding: 'utf8',
    maxBuffer: 20 * 1024 * 1024,
  });
  if (result.status !== 0) {
    throw new Error(`${contractPath}: fixture generator failed.\n${formatProcessOutput(result)}`);
  }
}

function loadPlan(args) {
  const manifestPath = resolveRepoPath(args.manifest);
  if (!existsSync(manifestPath)) throw new Error(`Manifest not found: ${args.manifest}`);

  const entries = readFileSync(manifestPath, 'utf8')
    .split(/\r?\n/)
    .map((line, index) => parseManifestEntry(line, index + 1, args.manifest))
    .filter(Boolean);
  if (entries.length === 0) throw new Error(`Manifest has no contract entries: ${args.manifest}`);

  return entries.map((entry, index) => {
    const absContractPath = resolveRepoPath(entry.contractPath);
    if (!existsSync(absContractPath)) throw new Error(`Contract not found: ${entry.contractPath}`);
    const contract = parseContract(readFileSync(absContractPath, 'utf8'), entry.contractPath);
    const fixture = String(contract.fixture ?? '');
    if (!fixture) throw new Error(`${entry.contractPath}: missing fixture=...`);
    const absFixturePath = resolveRepoPath(fixture);
    materializeGeneratedFixture(contract, entry.contractPath, absFixturePath);
    if (!existsSync(absFixturePath)) throw new Error(`${entry.contractPath}: fixture not found: ${fixture}`);
    if (!statSync(absFixturePath).isFile()) throw new Error(`${entry.contractPath}: fixture is not a file: ${fixture}`);

    const options = contract.options ? parseJson(contract.options, `${entry.contractPath}: options`) : {};
    const objReference = String(contract.obj_reference ?? '');
    const stlReference = String(contract.stl_reference ?? '');
    const browserReportReference = String(contract.browser_report_reference ?? '');
    const browserExpectation = contract.browser_expectation
      ? parseJson(contract.browser_expectation, `${entry.contractPath}: browser_expectation`)
      : undefined;
    const stem = contractOutputStem(entry.contractPath);
    const objExportBasename = validateBasename(
      String(contract.obj_export_basename ?? (objReference ? path.basename(objReference, path.extname(objReference)) : stem)),
      `${entry.contractPath}: obj_export_basename`,
    );
    const formats = args.formatsFromCli
      ? args.formats
      : (entry.formats ? entry.formats.filter(format => format === 'obj' || format === 'stl' || format === 'report') : args.formats);
    if (formats.length === 0) throw new Error(`${entry.contractPath}: no browser targets selected`);

    return {
      id: String(index),
      contractPath: entry.contractPath,
      contract,
      fixture,
      absFixturePath,
      inputFormat: normalizeInputFormat(options.format ?? contract.input_format, fixture),
      options,
      formats,
      stem,
      objReference,
      stlReference,
      browserReportReference,
      browserExpectation,
      objExportBasename,
      objMtllib: String(contract.obj_mtllib ?? `${objExportBasename}.mtl`),
    };
  });
}

function parseManifestEntry(rawLine, lineNumber, label) {
  const line = stripManifestComment(rawLine).trim();
  if (!line) return undefined;

  let contractPath;
  let formats;
  for (const field of line.split(/\s+/)) {
    const eq = field.indexOf('=');
    if (eq < 0) {
      if (contractPath) throw new Error(`${label}:${lineNumber}: unexpected manifest token '${field}'`);
      contractPath = field;
      continue;
    }
    const key = field.slice(0, eq);
    const value = field.slice(eq + 1);
    if (!value) throw new Error(`${label}:${lineNumber}: manifest field '${key}' requires a value`);
    if (key === 'contract' || key === 'path') {
      if (contractPath) throw new Error(`${label}:${lineNumber}: duplicate contract path`);
      contractPath = value;
    } else if (key === 'formats') {
      formats = parseFormatsForManifest(value, `${label}:${lineNumber}: formats`);
    } else if (key === 'tag' || key === 'tags' || key === 'note') {
      // Reserved manifest metadata.
    } else {
      throw new Error(`${label}:${lineNumber}: unknown manifest field '${key}'`);
    }
  }
  if (!contractPath) throw new Error(`${label}:${lineNumber}: missing contract path`);
  return { contractPath, formats };
}

function parseFormatsForManifest(value, label) {
  const formats = String(value).split(',').map(s => s.trim()).filter(Boolean);
  if (formats.length === 0) throw new Error(`${label}: expected at least one format`);
  for (const format of formats) {
    if (format !== 'json' && format !== 'obj' && format !== 'stl' && format !== 'report') {
      throw new Error(`${label}: unsupported format '${format}'. Expected json, obj, stl, or report.`);
    }
  }
  return formats;
}

function stripManifestComment(line) {
  const hash = line.indexOf('#');
  return hash < 0 ? line : line.slice(0, hash);
}

function parseContract(text, label) {
  const values = {};
  for (const [index, rawLine] of text.split(/\r?\n/).entries()) {
    const line = rawLine.trim();
    if (!line || line.startsWith('#')) continue;
    const eq = line.indexOf('=');
    if (eq < 0) throw new Error(`${label}:${index + 1}: expected key=value`);
    values[line.slice(0, eq)] = line.slice(eq + 1);
  }
  return values;
}

function parseJson(value, label) {
  try {
    return JSON.parse(value);
  } catch (error) {
    throw new Error(`${label}: invalid JSON: ${error.message}`);
  }
}

function normalizeInputFormat(format, fixture) {
  const lower = String(format || path.extname(fixture).slice(1)).toLowerCase();
  if (lower === 'mmcif') return 'cif';
  if (lower === 'binarycif') return 'bcif';
  if (lower === 'pdb' || lower === 'cif' || lower === 'bcif') return lower;
  throw new Error(`Unsupported input format '${lower}' for ${fixture}`);
}

function contractOutputStem(contractPath) {
  let base = path.basename(contractPath);
  base = base.replace(/\.contract$/, '');
  base = base.replace(/\.reference$/, '');
  base = base.replace(/\.(ply|obj|stl|json)$/i, '');
  return base;
}

function validateBasename(value, label) {
  if (!value || value !== path.basename(value) || value.includes('/') || value.includes('\\')) {
    throw new Error(`${label}: expected a filename stem without path separators, got ${value}`);
  }
  return value;
}

function plannedOutputPaths(args, item) {
  const outDir = resolveInputPath(args.outDir);
  return {
    obj: path.join(outDir, `${item.stem}.obj`),
    mtl: path.join(outDir, item.objMtllib),
    stl: path.join(outDir, `${item.stem}.stl`),
    report: path.join(outDir, `${item.stem}.browser-report.json`),
  };
}

function printPlan(args, plan, bundlePath) {
  console.log(`Mol* browser reference conversion ${args.dryRun ? 'dry run' : args.buildOnly ? 'bundle build' : 'plan'}`);
  console.log(`manifest: ${args.manifest}`);
  console.log(`molstar: ${args.molstarDir}`);
  console.log(`bundle: ${displayPath(bundlePath)}`);
  console.log(`chrome: ${args.chromePath}`);
  console.log(`fixtures: ${plan.length}`);
  for (const item of plan) {
    const output = plannedOutputPaths(args, item);
    console.log(`- ${item.fixture} input=${item.inputFormat} targets=${item.formats.join(',')} stem=${item.stem} obj-basename=${item.objExportBasename}`);
    console.log(`  exports obj=${displayPath(output.obj)} mtl=${displayPath(output.mtl)} stl=${displayPath(output.stl)}`);
  }
}

function inspectBrowserPrerequisites(args) {
  const molstarDir = resolveRepoPath(args.molstarDir);
  const packageJsonPath = path.join(molstarDir, 'package.json');
  const esbuildBin = path.join(molstarDir, 'node_modules/esbuild/bin/esbuild');
  const sourceFiles = [
    'mol-plugin/context.ts',
    'mol-plugin/spec.ts',
    'mol-plugin-state/actions/file.ts',
    'extensions/geo-export/obj-exporter.ts',
    'extensions/geo-export/stl-exporter.ts',
  ];
  const missing = [];
  if (!existsSync(molstarDir)) missing.push(`${args.molstarDir}/`);
  if (!existsSync(packageJsonPath)) missing.push(`${args.molstarDir}/package.json`);
  if (!existsSync(esbuildBin)) missing.push(`${args.molstarDir}/node_modules/esbuild/bin/esbuild`);
  for (const file of sourceFiles) {
    if (!existsSync(path.join(molstarDir, 'src', file))) missing.push(`${args.molstarDir}/src/${file}`);
  }
  if (!existsSync(args.chromePath)) missing.push(args.chromePath);
  return { molstarDir, packageJsonPath, esbuildBin, missing };
}

function ensureBrowserPrerequisites(args, inspection, { requireChrome }) {
  const missing = requireChrome ? inspection.missing : inspection.missing.filter(item => item !== args.chromePath);
  if (missing.length === 0) return;
  throw new Error([
    'Mol* browser reference exporter cannot run from the current local checkout state.',
    ...missing.map(item => `- missing ${item}`),
  ].join('\n'));
}

function buildHarnessBundle(args, inspection) {
  const outDir = path.join(resolveInputPath(args.outDir), '.browser-harness');
  mkdirSync(outDir, { recursive: true });
  const bundlePath = path.join(outDir, 'molfig-molstar-browser-reference.js');
  const result = spawnSync(inspection.esbuildBin, [
    'scripts/molstar-browser-reference-harness.ts',
    '--bundle',
    '--format=iife',
    '--global-name=MolfigMolstarBrowserReference',
    `--outfile=${bundlePath}`,
    '--tsconfig=artifacts/molstar/tsconfig.json',
    '--define:process.env.NODE_ENV="production"',
    '--define:process.env.DEBUG=false',
    '--define:__MOLSTAR_PLUGIN_VERSION__="5.9.0"',
    '--define:__MOLSTAR_BUILD_TIMESTAMP__=0',
    '--log-level=warning',
  ], {
    cwd: rootDir,
    encoding: 'utf8',
    maxBuffer: 20 * 1024 * 1024,
  });
  if (result.status !== 0) {
    throw new Error(`Mol* browser harness bundle failed.\n${formatProcessOutput(result)}`);
  }
  return bundlePath;
}

function formatProcessOutput(result) {
  const lines = [];
  if (result.error) lines.push(`error: ${result.error.message}`);
  if (result.stdout?.trim()) lines.push(`stdout:\n${result.stdout.trim()}`);
  if (result.stderr?.trim()) lines.push(`stderr:\n${result.stderr.trim()}`);
  if (lines.length === 0) lines.push(`exit status ${result.status}`);
  return lines.join('\n');
}

async function runBrowserConversion(args, plan, bundlePath) {
  const server = await startServer(args, plan, bundlePath);
  const profileDir = path.join(resolveInputPath(args.outDir), '.chrome-profile');
  if (existsSync(profileDir) && !args.keepProfile) rmSync(profileDir, { recursive: true, force: true });
  mkdirSync(profileDir, { recursive: true });

  let chrome;
  let cdp;
  const comparisonFailures = [];
  try {
    const debugPort = await getFreePort();
    chrome = launchChrome(args, debugPort, profileDir, server.url);
    const target = await waitForChromeTarget(debugPort, args.timeoutMs);
    cdp = await CdpClient.connect(target.webSocketDebuggerUrl);
    await cdp.call('Runtime.enable');
    await cdp.call('Page.enable');
    await waitForHarness(cdp, args.timeoutMs);

    for (const item of plan) {
      console.log(`Converting ${item.fixture}`);
      const result = await runFixtureInBrowser(cdp, args, item, server.url);
      writeBrowserReport(args, item, result);
      const contractChecks = validateBrowserContract(item, result);
      for (const check of contractChecks) console.log(check.message);
      comparisonFailures.push(...contractChecks.filter(check => !check.ok));
      if (args.renderObjectReport) printBrowserRenderObjectReport(item, result);
      if (args.debugStlFacets.length > 0) printBrowserStlFacetDebug(item, result);
      printBrowserPassSummary(item, result);
      const molfigContext = args.compareMolfig
        ? {
            optionsPath: writeMolfigDiffOptions(args, item),
            command: resolveMolfigDiffCommand(args),
          }
        : undefined;
      if (molfigContext && args.debugStlFacets.length > 0) {
        for (const check of compareBrowserStlFacetDebugToMolfig(molfigContext.command, item, molfigContext.optionsPath, result)) {
          console.log(check.message);
          if (!check.ok) comparisonFailures.push(check);
        }
      }
      const comparisonChecks = [];
      if (args.compareReferences) comparisonChecks.push(...compareBrowserOutputsToReferences(args, item));
      if (args.compareMolfig) comparisonChecks.push(...compareBrowserOutputsToMolfig(args, item, molfigContext));
      for (const check of comparisonChecks) console.log(check.message);
      comparisonFailures.push(...comparisonChecks.filter(check => !check.ok));
    }

    if (comparisonFailures.length > 0) {
      throw new Error(`Browser comparison failed: ${comparisonFailures.map(check => check.label).join(', ')}`);
    }
  } finally {
    await cdp?.close();
    if (chrome) {
      chrome.kill('SIGTERM');
      await waitForExit(chrome, 2000).catch(() => chrome.kill('SIGKILL'));
    }
    await server.close();
    if (!args.keepProfile) rmSync(profileDir, { recursive: true, force: true });
  }
}

function launchChrome(args, debugPort, profileDir, url) {
  const chromeArgs = [
    `--remote-debugging-port=${debugPort}`,
    `--user-data-dir=${profileDir}`,
    '--no-first-run',
    '--no-default-browser-check',
    '--disable-background-networking',
    '--disable-sync',
    '--disable-extensions',
    '--disable-popup-blocking',
    '--autoplay-policy=no-user-gesture-required',
    '--window-size=1024,1024',
  ];
  if (!args.headed) chromeArgs.push('--headless=new');
  chromeArgs.push(url);

  const chrome = spawn(args.chromePath, chromeArgs, { stdio: ['ignore', 'pipe', 'pipe'] });
  chrome.stdout.on('data', data => process.stdout.write(String(data)));
  chrome.stderr.on('data', data => {
    const text = String(data);
    if (!text.includes('DevTools listening')) process.stderr.write(text);
  });
  return chrome;
}

async function runFixtureInBrowser(cdp, args, item, baseUrl) {
  const expression = `window.molfigBrowserReferenceExport(${JSON.stringify({
    id: item.id,
    fixtureUrl: `${baseUrl}/fixture/${encodeURIComponent(item.id)}`,
    fileName: path.basename(item.fixture),
    formats: item.formats,
    objExportBasename: item.objExportBasename,
    dataFormat: dataFormatProviderId(item),
    structurePreset: item.options['structure-preset'] ?? item.options.structure_preset,
    representation: String(item.options.representation ?? 'default'),
    theme: browserTheme(item.options),
    sizeThresholds: browserObjectOption(item.options, 'viewer-size-thresholds', 'viewerSizeThresholds'),
    gaussianSurfaceParams: browserObjectOption(item.options, 'gaussian-surface-params', 'gaussianSurfaceParams'),
    molecularSurfaceParams: browserObjectOption(item.options, 'molecular-surface-params', 'molecularSurfaceParams'),
    renderObjectReport: args.renderObjectReport,
    exporterOptions: {
      primitivesQuality: String(item.options['export-primitives-quality'] ?? item.options.export_primitives_quality ?? 'auto'),
    },
    debugStlFacets: args.debugStlFacets,
  })})`;
  const response = await cdp.call('Runtime.evaluate', {
    expression,
    awaitPromise: true,
    returnByValue: true,
  }, args.timeoutMs);
  if (response.exceptionDetails) {
    throw new Error(formatRuntimeException(response.exceptionDetails));
  }
  if (response.result?.subtype === 'error') {
    throw new Error(response.result.description || 'browser export failed');
  }
  return response.result.value;
}

function browserTheme(options) {
  const theme = options.theme ?? options['viewer-theme'] ?? options.viewer_theme;
  return theme && typeof theme === 'object' && !Array.isArray(theme) ? theme : undefined;
}

function browserObjectOption(options, ...keys) {
  for (const key of keys) {
    const value = options[key];
    if (value !== undefined) {
      if (!value || typeof value !== 'object' || Array.isArray(value)) {
        throw new Error(`browser option '${key}' must be an object`);
      }
      return value;
    }
  }
  return undefined;
}

function dataFormatProviderId(item) {
  const provider = item.options['data-format'] ?? item.options.data_format;
  if (provider !== undefined) return String(provider);
  return 'auto';
}

function formatRuntimeException(exception) {
  const text = exception.exception?.description || exception.text || 'browser runtime exception';
  const stack = exception.stackTrace?.callFrames
    ?.slice(0, 8)
    ?.map(frame => `    at ${frame.functionName || '<anonymous>'} (${frame.url}:${frame.lineNumber + 1}:${frame.columnNumber + 1})`)
    ?.join('\n');
  return stack ? `${text}\n${stack}` : text;
}

function printBrowserRenderObjectReport(item, result) {
  console.log(`Mol* browser structures for ${item.stem}: ${JSON.stringify(result.structures)}`);
  console.log(`Mol* browser render objects for ${item.stem}: ${result.renderObjects.length}`);
  for (const object of result.renderObjects) {
    const sphere = object.boundingSphere
      ? ` sphereCenter=${formatNumberArray(object.boundingSphere.center)} sphereRadius=${formatNumber(object.boundingSphere.radius)} sphereExtrema=${object.boundingSphere.extremaCount}`
      : '';
    const preset = object.tag
      ? ` component=${object.component ?? '<unknown>'} tag=${object.tag} representation=${object.representation ?? '<unknown>'} colorTheme=${object.colorTheme ?? '<unknown>'} carbonColorTheme=${object.carbonColorTheme ?? '<none>'} representationOrder=${object.representationOrder ?? '<unknown>'}`
      : '';
    const exportCount = object.stlTriangleCount === undefined ? '' : ` stlTriangleCount=${object.stlTriangleCount}`;
    const primitiveCount = object.primitiveCount === undefined ? '' : ` primitiveCount=${object.primitiveCount}`;
    const cylinderCaps = object.cylinderCapHistogram === undefined ? '' : ` cylinderCaps=${JSON.stringify(object.cylinderCapHistogram)}`;
    const cylinderSamples = object.cylinderSamples === undefined ? '' : ` cylinderSamples=${JSON.stringify(object.cylinderSamples)}`;
    const meshVertexSamples = object.meshVertexSamples === undefined ? '' : ` meshVertexSamples=${JSON.stringify(object.meshVertexSamples)}`;
    console.log(`  [${object.index}] type=${object.type} geometry=${object.geometry} visible=${object.visible} drawCount=${object.drawCount} vertexCount=${object.vertexCount} groupCount=${object.groupCount} instanceCount=${object.instanceCount}${preset}${exportCount}${primitiveCount}${cylinderCaps}${cylinderSamples}${meshVertexSamples}${sphere}`);
  }
  console.log(`  totalDrawCount=${result.totalDrawCount} visibleDrawCount=${result.visibleDrawCount} hiddenDrawCount=${result.hiddenDrawCount} exportableDrawCount=${result.exportableDrawCount}`);
  console.log(`  webgl webgl2=${result.webgl.webgl2} fragDepth=${result.webgl.fragDepth} textureFloat=${result.webgl.textureFloat}`);
}

function writeBrowserReport(args, item, result) {
  const output = plannedOutputPaths(args, item);
  const structures = item.browserExpectation
    ? result.structures
    : result.structures.map(({ polymerResidueCount, sizeClass, sizeThresholds, ...legacy }) => legacy);
  const report = {
    ...(item.options['structure-preset'] || item.options.structure_preset || item.browserExpectation
      ? { structures }
      : {}),
    renderObjects: result.renderObjects.map(object => ({
      index: object.index,
      type: object.type,
      geometry: object.geometry,
      visible: object.visible,
      drawCount: object.drawCount,
      vertexCount: object.vertexCount,
      groupCount: object.groupCount,
      instanceCount: object.instanceCount,
      component: object.component,
      tag: object.tag,
      representation: object.representation,
      colorTheme: object.colorTheme,
      carbonColorTheme: object.carbonColorTheme,
      representationOrder: object.representationOrder,
      ...(item.browserExpectation ? {
        visuals: object.visuals,
        surfaceParams: object.surfaceParams,
      } : {}),
      stlTriangleCount: object.stlTriangleCount,
      primitiveCount: object.primitiveCount,
      cylinderCapHistogram: object.cylinderCapHistogram,
      cylinderSamples: object.cylinderSamples,
      meshVertexSamples: object.meshVertexSamples,
      boundingSphere: object.boundingSphere && {
        center: object.boundingSphere.center,
        radius: object.boundingSphere.radius,
        extremaCount: object.boundingSphere.extremaCount,
        extrema: object.boundingSphere.extrema,
      },
    })),
    sceneBoundingSphere: result.sceneBoundingSphere,
    totalDrawCount: result.totalDrawCount,
    visibleDrawCount: result.visibleDrawCount,
    hiddenDrawCount: result.hiddenDrawCount,
    exportableDrawCount: result.exportableDrawCount,
    uploads: result.uploads,
    webgl: result.webgl,
  };
  mkdirSync(path.dirname(output.report), { recursive: true });
  writeFileSync(output.report, `${JSON.stringify(report, null, 2)}\n`);
}

function validateBrowserContract(item, result) {
  const expected = item.browserExpectation;
  if (!expected) return [];
  const failures = [];
  const checks = [];
  const record = (ok, label, actual, wanted) => {
    const message = `- ${ok ? 'PASS' : 'FAIL'} ${item.stem}: ${label}${ok ? '' : `; expected=${JSON.stringify(wanted)} actual=${JSON.stringify(actual)}`}`;
    const check = { ok, label: `${item.stem}: ${label}`, message };
    checks.push(check);
    if (!ok) failures.push(check);
  };

  if (expected.renderObjectCount !== undefined) {
    record(result.renderObjects.length === expected.renderObjectCount, 'renderObjectCount', result.renderObjects.length, expected.renderObjectCount);
  }
  for (const field of ['totalDrawCount', 'visibleDrawCount', 'hiddenDrawCount', 'exportableDrawCount']) {
    if (expected[field] !== undefined) record(result[field] === expected[field], field, result[field], expected[field]);
  }
  if (expected.uploads) record(isObjectSubset(result.uploads, expected.uploads), 'uploads', result.uploads, expected.uploads);
  if (expected.structure) {
    const actual = result.structures[expected.structure.index ?? 0];
    record(Boolean(actual), 'structure exists', actual, expected.structure);
    if (actual) validateExpectedFields(record, 'structure', actual, expected.structure, new Set(['index']));
  }
  for (const [index, objectExpectation] of (expected.renderObjects ?? []).entries()) {
    const actual = findExpectedRenderObject(result.renderObjects, objectExpectation);
    record(Boolean(actual), `renderObjects[${index}] exists`, actual, objectExpectation);
    if (!actual) continue;
    validateExpectedFields(record, `renderObjects[${index}]`, actual, objectExpectation, new Set([
      'minDrawCount', 'minVertexCount', 'minGroupCount',
    ]));
    for (const [field, actualField] of [
      ['minDrawCount', 'drawCount'],
      ['minVertexCount', 'vertexCount'],
      ['minGroupCount', 'groupCount'],
    ]) {
      if (objectExpectation[field] !== undefined) {
        record(actual[actualField] >= objectExpectation[field], `renderObjects[${index}].${field}`, actual[actualField], `>= ${objectExpectation[field]}`);
      }
    }
  }
  if (failures.length === 0) checks.unshift({
    ok: true,
    label: `${item.stem}: browser expectation`,
    message: `- PASS ${item.stem}: browser expectation (${checks.length} assertions)`,
  });
  return checks;
}

function findExpectedRenderObject(objects, expected) {
  return objects.find(object => {
    if (expected.index !== undefined && object.index !== expected.index) return false;
    if (expected.tag !== undefined && object.tag !== expected.tag) return false;
    if (expected.representation !== undefined && object.representation !== expected.representation) return false;
    if (expected.geometry !== undefined && object.geometry !== expected.geometry) return false;
    return true;
  });
}

function validateExpectedFields(record, prefix, actual, expected, ignored = new Set()) {
  for (const [key, wanted] of Object.entries(expected)) {
    if (ignored.has(key)) continue;
    const got = actual?.[key];
    if (wanted && typeof wanted === 'object' && !Array.isArray(wanted)) {
      record(isObjectSubset(got, wanted), `${prefix}.${key}`, got, wanted);
    } else {
      record(JSON.stringify(got) === JSON.stringify(wanted), `${prefix}.${key}`, got, wanted);
    }
  }
}

function isObjectSubset(actual, expected) {
  if (!actual || typeof actual !== 'object') return false;
  return Object.entries(expected).every(([key, value]) => {
    const got = actual[key];
    if (value && typeof value === 'object' && !Array.isArray(value)) return isObjectSubset(got, value);
    return JSON.stringify(got) === JSON.stringify(value);
  });
}

function formatNumber(value) {
  return Number.isFinite(value) ? Number(value).toExponential(17) : String(value);
}

function formatNumberArray(values) {
  return `[${values.map(formatNumber).join(',')}]`;
}

function printBrowserPassSummary(item, result) {
  const parts = item.formats.map(format => {
    if (format === 'report') return `REPORT ${result.renderObjects.length} render objects`;
    const upload = result.uploads[format];
    return `${format.toUpperCase()} ${upload?.byteLength ?? 0} bytes`;
  });
  console.log(`- PASS ${item.stem}: browser ${parts.join('; ')}; exportableDrawCount=${result.exportableDrawCount}; fragDepth=${result.webgl.fragDepth}`);
}

function printBrowserStlFacetDebug(item, result) {
  const facets = result.debug?.stlFacets ?? [];
  for (const facet of facets) {
    console.log(`Mol* browser STL facet debug for ${item.stem}: ${JSON.stringify(facet)}`);
  }
}

function compareBrowserOutputsToReferences(args, item) {
  const output = plannedOutputPaths(args, item);
  const checks = [];

  if (item.formats.includes('obj')) {
    if (!item.objReference) throw new Error(`${item.contractPath}: missing obj_reference for --compare-references`);
    checks.push(compareTextFiles(resolveRepoPath(item.objReference), output.obj, `${item.stem}: obj`));
    checks.push(compareTextFiles(resolveObjMtlReferencePath(item), output.mtl, `${item.stem}: mtl`));
  }
  if (item.formats.includes('stl')) {
    if (!item.stlReference) throw new Error(`${item.contractPath}: missing stl_reference for --compare-references`);
    checks.push(compareBinaryFiles(resolveRepoPath(item.stlReference), output.stl, `${item.stem}: stl`));
  }
  if (item.browserReportReference) {
    checks.push(compareTextFiles(
      resolveRepoPath(item.browserReportReference),
      output.report,
      `${item.stem}: browser report`,
    ));
  }

  return checks;
}

function compareBrowserOutputsToMolfig(args, item, molfigContext) {
  const output = plannedOutputPaths(args, item);
  const optionsPath = molfigContext?.optionsPath ?? writeMolfigDiffOptions(args, item);
  const command = molfigContext?.command ?? resolveMolfigDiffCommand(args);
  const checks = [];

  for (const format of item.formats) {
    if (format !== 'obj' && format !== 'stl') continue;
    checks.push(runMolfigDiff(command, format, item.absFixturePath, optionsPath, output[format], `${item.stem}: molfig-vs-browser ${format}`));
  }
  return checks;
}

function compareBrowserStlFacetDebugToMolfig(command, item, optionsPath, browserResult) {
  const facets = browserResult.debug?.stlFacets ?? [];
  const checks = [];
  for (const facet of facets) {
    if (!facet?.found) continue;
    const result = spawnSync(command.command, [
      ...command.args,
      '--stl-export-facet-context',
      String(facet.stlFacet),
      'stl',
      item.absFixturePath,
      optionsPath,
    ], {
      cwd: rootDir,
      encoding: 'utf8',
      maxBuffer: 20 * 1024 * 1024,
    });
    const label = `${item.stem}: browser-vs-molfig STL facet ${facet.stlFacet} center`;
    if (result.status !== 0) {
      checks.push({
        ok: true,
        label,
        message: `- DIAG ${label}: molfig facet context unavailable: ${formatMolfigDiffOutput(result)}`,
      });
      continue;
    }
    const parsed = parseJsonOutput(result.stdout);
    if (!parsed) {
      checks.push({
        ok: true,
        label,
        message: `- DIAG ${label}: molfig facet context was not JSON: ${formatMolfigDiffOutput(result)}`,
      });
      continue;
    }
    checks.push({
      ok: true,
      label,
      message: formatStlFacetCenterDiagnostic(label, facet, parsed),
    });
  }
  return checks;
}

function parseJsonOutput(output) {
  const text = String(output ?? '').trim();
  if (!text) return undefined;
  try {
    return JSON.parse(text.split(/\r?\n/).at(-1));
  } catch {
    return undefined;
  }
}

function formatStlFacetCenterDiagnostic(label, browserFacet, molfigFacet) {
  const centerDelta = vectorDelta(browserFacet.centerOffset, molfigFacet.vertex_offset);
  const boxMinDelta = vectorDelta(browserFacet.boxMin, molfigFacet.export_box_min);
  const boxMaxDelta = vectorDelta(browserFacet.boxMax, molfigFacet.export_box_max);
  const boxMinPointDeltas = pointDeltas(browserFacet.boxMinPoints, molfigFacet.export_box_min_points);
  const boxMaxPointDeltas = pointDeltas(browserFacet.boxMaxPoints, molfigFacet.export_box_max_points);
  const boxMinIndicesMatch = JSON.stringify(browserFacet.boxMinIndices) === JSON.stringify(molfigFacet.export_box_min_indices);
  const boxMaxIndicesMatch = JSON.stringify(browserFacet.boxMaxIndices) === JSON.stringify(molfigFacet.export_box_max_indices);
  const browserBits = browserFacet.centeredVertexBits;
  const molfigBits = molfigFacet.target_face?.stl_vertex_bits;
  const bitsMatch = JSON.stringify(browserBits) === JSON.stringify(molfigBits);
  return `- DIAG ${label}: center_offset_delta=${formatVector(centerDelta)}; box_min_delta=${formatVector(boxMinDelta)}; box_max_delta=${formatVector(boxMaxDelta)}; box_min_indices_match=${boxMinIndicesMatch}; box_max_indices_match=${boxMaxIndicesMatch}; box_min_point_deltas=${formatPointDeltas(boxMinPointDeltas)}; box_max_point_deltas=${formatPointDeltas(boxMaxPointDeltas)}; stl_vertex_bits_match=${bitsMatch}; browser_center_offset=${formatVector(browserFacet.centerOffset)}; molfig_vertex_offset=${formatVector(molfigFacet.vertex_offset)}`;
}

function vectorDelta(reference, generated) {
  if (!Array.isArray(reference) || !Array.isArray(generated) || reference.length < 3 || generated.length < 3) {
    return undefined;
  }
  return [
    Number(generated[0]) - Number(reference[0]),
    Number(generated[1]) - Number(reference[1]),
    Number(generated[2]) - Number(reference[2]),
  ];
}

function formatVector(value) {
  if (!Array.isArray(value)) return 'n/a';
  return `[${value.map(component => Number(component).toExponential(9)).join(',')}]`;
}

function pointDeltas(referencePoints, generatedPoints) {
  if (!Array.isArray(referencePoints) || !Array.isArray(generatedPoints)) return undefined;
  return referencePoints.map((reference, index) => vectorDelta(reference, generatedPoints[index]));
}

function formatPointDeltas(value) {
  if (!Array.isArray(value)) return 'n/a';
  return `[${value.map(formatVector).join(',')}]`;
}

function writeMolfigDiffOptions(args, item) {
  const dir = path.join(resolveInputPath(args.outDir), '.molfig-options');
  mkdirSync(dir, { recursive: true });
  const optionsPath = path.join(dir, `${item.stem}.json`);
  writeFileSync(optionsPath, `${JSON.stringify(item.options, null, 2)}\n`);
  return optionsPath;
}

function resolveMolfigDiffCommand(args) {
  if (args.molfigDiff) {
    if (args.molfigDiff === 'cargo') return { command: 'cargo', args: ['run', '--quiet', '--bin', 'molfig-diff', '--'] };
    return { command: resolveInputPath(args.molfigDiff), args: [] };
  }

  for (const candidate of [
    'target/debug/molfig-diff',
    'target/release/molfig-diff',
  ]) {
    const abs = resolveRepoPath(candidate);
    if (existsSync(abs)) return { command: abs, args: [] };
  }
  return { command: 'cargo', args: ['run', '--quiet', '--bin', 'molfig-diff', '--'] };
}

function runMolfigDiff(command, format, fixturePath, optionsPath, referencePath, label) {
  const result = spawnSync(command.command, [
    ...command.args,
    format,
    fixturePath,
    optionsPath,
    referencePath,
  ], {
    cwd: rootDir,
    encoding: 'utf8',
    maxBuffer: 20 * 1024 * 1024,
  });
  const ok = result.status === 0;
  const output = formatMolfigDiffOutput(result);
  return {
    ok,
    label,
    message: `- ${ok ? 'PASS' : 'FAIL'} ${label}: ${output}`,
  };
}

function formatMolfigDiffOutput(result) {
  const output = [result.stdout, result.stderr]
    .filter(Boolean)
    .map(part => part.trim())
    .filter(Boolean)
    .join('\n');
  if (output) return output.replace(/\n/g, '\n  ');
  if (result.error) return result.error.message;
  return `exit status ${result.status}`;
}

function resolveObjMtlReferencePath(item) {
  const mtllib = item.objMtllib;
  if (path.isAbsolute(mtllib) || mtllib.includes('/') || mtllib.includes('\\')) return resolveInputPath(mtllib);
  if (!item.objReference) return resolveRepoPath(mtllib);
  return path.join(path.dirname(resolveRepoPath(item.objReference)), mtllib);
}

function compareTextFiles(referencePath, generatedPath, label) {
  const reference = readFileSync(referencePath, 'utf8');
  const generated = readFileSync(generatedPath, 'utf8');
  if (reference === generated) {
    return { ok: true, label, message: `- PASS ${label}: text match (${Buffer.byteLength(reference)} bytes)` };
  }
  const referenceLines = reference.split(/\r?\n/);
  const generatedLines = generated.split(/\r?\n/);
  const lineCount = Math.max(referenceLines.length, generatedLines.length);
  let firstLine = 0;
  while (firstLine < lineCount && referenceLines[firstLine] === generatedLines[firstLine]) firstLine += 1;
  return {
    ok: false,
    label,
    message: `- FAIL ${label}: first difference at line ${firstLine + 1}; reference=${quoteLine(referenceLines[firstLine])}, generated=${quoteLine(generatedLines[firstLine])}; reference_lines=${referenceLines.length}, generated_lines=${generatedLines.length}`,
  };
}

function compareBinaryFiles(referencePath, generatedPath, label) {
  const reference = readFileSync(referencePath);
  const generated = readFileSync(generatedPath);
  if (reference.equals(generated)) {
    return { ok: true, label, message: `- PASS ${label}: byte match (${reference.length} bytes)` };
  }
  const first = firstDifferentByte(reference, generated);
  const referenceByte = first < reference.length ? reference[first] : '<eof>';
  const generatedByte = first < generated.length ? generated[first] : '<eof>';
  return {
    ok: false,
    label,
    message: `- FAIL ${label}: first difference at byte ${first}; reference=${referenceByte}, generated=${generatedByte}; reference_len=${reference.length}, generated_len=${generated.length}${stlFacetContext(reference, generated, first)}`,
  };
}

function firstDifferentByte(reference, generated) {
  const len = Math.min(reference.length, generated.length);
  for (let i = 0; i < len; i++) {
    if (reference[i] !== generated[i]) return i;
  }
  return len;
}

function quoteLine(line) {
  if (line === undefined) return '"<eof>"';
  return JSON.stringify(line.length > 160 ? `${line.slice(0, 157)}...` : line);
}

function stlFacetContext(reference, generated, firstDiff) {
  if (!looksLikeBinaryStl(reference) || !looksLikeBinaryStl(generated) || firstDiff < 84) return '';
  const facet = Math.floor((firstDiff - 84) / 50);
  const offset = (firstDiff - 84) % 50;
  const field = offset < 48 ? stlFacetField(offset) : 'attribute byte count';
  return `; stl_context=facet ${facet} ${field} byte ${offset % 4}`;
}

function looksLikeBinaryStl(bytes) {
  if (bytes.length < 84) return false;
  const facets = bytes.readUInt32LE(80);
  return 84 + facets * 50 === bytes.length;
}

function stlFacetField(offset) {
  const names = [
    'normal.x', 'normal.y', 'normal.z',
    'v1.x', 'v1.y', 'v1.z',
    'v2.x', 'v2.y', 'v2.z',
    'v3.x', 'v3.y', 'v3.z',
  ];
  return names[Math.floor(offset / 4)] ?? 'unknown';
}

async function startServer(args, plan, bundlePath) {
  const uploads = new Map();
  const server = createServer(async (req, res) => {
    try {
      const url = new URL(req.url || '/', 'http://127.0.0.1');
      if (req.method === 'GET' && url.pathname === '/') {
        sendText(res, 200, browserHtml(), 'text/html; charset=utf-8');
        return;
      }
      if (req.method === 'GET' && url.pathname === '/bundle.js') {
        sendBytes(res, 200, readFileSync(bundlePath), 'application/javascript');
        return;
      }
      if (req.method === 'GET' && url.pathname.startsWith('/fixture/')) {
        const id = decodeURIComponent(url.pathname.slice('/fixture/'.length));
        const item = plan.find(item => item.id === id);
        if (!item) return sendText(res, 404, `unknown fixture id ${id}`);
        sendBytes(res, 200, readFileSync(item.absFixturePath), 'application/octet-stream');
        return;
      }
      if (req.method === 'POST' && url.pathname.startsWith('/upload/')) {
        const [, , id, format] = url.pathname.split('/').map(decodeURIComponent);
        const item = plan.find(item => item.id === id);
        if (!item) return sendText(res, 404, `unknown upload id ${id}`);
        if (!['obj', 'mtl', 'stl'].includes(format)) return sendText(res, 400, `unsupported upload format ${format}`);
        const chunks = [];
        for await (const chunk of req) chunks.push(chunk);
        const body = Buffer.concat(chunks);
        const output = plannedOutputPaths(args, item);
        const outputPath = output[format];
        mkdirSync(path.dirname(outputPath), { recursive: true });
        writeFileSync(outputPath, body);
        uploads.set(`${id}:${format}`, body.length);
        sendText(res, 200, 'ok');
        return;
      }
      sendText(res, 404, 'not found');
    } catch (error) {
      sendText(res, 500, error.stack || error.message || String(error));
    }
  });

  await new Promise((resolve, reject) => {
    server.once('error', reject);
    server.listen(0, '127.0.0.1', resolve);
  });
  const address = server.address();
  const url = `http://127.0.0.1:${address.port}`;
  return {
    url,
    uploads,
    close: () => new Promise(resolve => server.close(resolve)),
  };
}

function browserHtml() {
  return `<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>molfig Mol* browser reference exporter</title>
  <style>
    html, body, #app { width: 100%; height: 100%; margin: 0; overflow: hidden; }
    #app { position: relative; }
  </style>
</head>
<body>
  <div id="app"></div>
  <script src="/bundle.js"></script>
</body>
</html>`;
}

function sendText(res, status, text, contentType = 'text/plain; charset=utf-8') {
  sendBytes(res, status, Buffer.from(text), contentType);
}

function sendBytes(res, status, bytes, contentType) {
  res.writeHead(status, {
    'content-type': contentType,
    'content-length': bytes.length,
    'cache-control': 'no-store',
  });
  res.end(bytes);
}

async function getFreePort() {
  const server = createServer();
  await new Promise((resolve, reject) => {
    server.once('error', reject);
    server.listen(0, '127.0.0.1', resolve);
  });
  const port = server.address().port;
  await new Promise(resolve => server.close(resolve));
  return port;
}

async function waitForChromeTarget(debugPort, timeoutMs) {
  const deadline = Date.now() + timeoutMs;
  let lastError;
  while (Date.now() < deadline) {
    try {
      const targets = await fetchJson(`http://127.0.0.1:${debugPort}/json/list`);
      const page = targets.find(target => target.type === 'page' && target.webSocketDebuggerUrl);
      if (page) return page;
    } catch (error) {
      lastError = error;
    }
    await sleep(250);
  }
  throw new Error(`Timed out waiting for Chrome DevTools target${lastError ? `: ${lastError.message}` : ''}`);
}

async function waitForHarness(cdp, timeoutMs) {
  const deadline = Date.now() + timeoutMs;
  while (Date.now() < deadline) {
    const response = await cdp.call('Runtime.evaluate', {
      expression: 'typeof window.molfigBrowserReferenceExport === "function"',
      returnByValue: true,
    }, 5000);
    if (response.result?.value === true) return;
    await sleep(250);
  }
  throw new Error('Timed out waiting for browser reference harness');
}

async function fetchJson(url) {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`${url}: ${response.status} ${await response.text()}`);
  return response.json();
}

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

function waitForExit(child, timeoutMs) {
  return new Promise((resolve, reject) => {
    const timer = setTimeout(() => reject(new Error('process exit timeout')), timeoutMs);
    child.once('exit', code => {
      clearTimeout(timer);
      resolve(code);
    });
  });
}

class CdpClient {
  constructor(ws) {
    this.ws = ws;
    this.nextId = 1;
    this.pending = new Map();
    ws.addEventListener('message', event => this.onMessage(event));
    ws.addEventListener('close', () => {
      for (const { reject } of this.pending.values()) reject(new Error('CDP connection closed'));
      this.pending.clear();
    });
  }

  static async connect(url) {
    const ws = new WebSocket(url);
    await new Promise((resolve, reject) => {
      ws.addEventListener('open', resolve, { once: true });
      ws.addEventListener('error', reject, { once: true });
    });
    return new CdpClient(ws);
  }

  call(method, params = {}, timeoutMs = 30000) {
    const id = this.nextId++;
    const message = JSON.stringify({ id, method, params });
    const promise = new Promise((resolve, reject) => {
      const timer = setTimeout(() => {
        this.pending.delete(id);
        reject(new Error(`CDP ${method} timed out after ${timeoutMs}ms`));
      }, timeoutMs);
      this.pending.set(id, { resolve, reject, timer, method });
    });
    this.ws.send(message);
    return promise;
  }

  onMessage(event) {
    const message = JSON.parse(event.data);
    if (message.id === undefined) return;
    const pending = this.pending.get(message.id);
    if (!pending) return;
    this.pending.delete(message.id);
    clearTimeout(pending.timer);
    if (message.error) {
      pending.reject(new Error(`CDP ${pending.method} failed: ${message.error.message}`));
    } else {
      pending.resolve(message.result);
    }
  }

  close() {
    this.ws.close();
    return Promise.resolve();
  }
}

async function main() {
  try {
    const args = parseArgs(process.argv.slice(2));
    if (args.help) {
      printHelp();
      return;
    }
    const plan = loadPlan(args);
    const inspection = inspectBrowserPrerequisites(args);
    ensureBrowserPrerequisites(args, inspection, { requireChrome: !args.dryRun && !args.buildOnly });
    const bundlePath = buildHarnessBundle(args, inspection);
    printPlan(args, plan, bundlePath);
    if (args.dryRun || args.buildOnly) return;
    await runBrowserConversion(args, plan, bundlePath);
  } catch (error) {
    console.error(`error: ${error.stack || error.message || error}`);
    process.exitCode = 1;
  }
}

main();
