#!/usr/bin/env node

import { createRequire } from 'node:module';
import { spawnSync } from 'node:child_process';
import {
  existsSync,
  mkdirSync,
  readFileSync,
  statSync,
  writeFileSync,
} from 'node:fs';
import path from 'node:path';
import process from 'node:process';
import { fileURLToPath } from 'node:url';

const rootDir = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const molstarReferenceCommit = '1b8117d3f10f7c978aabb5a0d3d47370635aefe4';
const defaultManifest = 'tests/expected/molstar-reference/reference-fixtures.txt';
const defaultArtifactManifest = 'tests/expected/reference-artifacts.json';
const defaultMolstarDir = 'artifacts/molstar';
const defaultOutDir = 'tests/expected/molstar-reference';
const defaultFormats = ['json', 'obj', 'stl'];

function printHelp() {
  console.log(`Usage: node scripts/molstar-reference-convert.mjs [options]

Replay molfig reference fixture contracts through the pinned local Mol* build.

Options:
  --dry-run, --check       Validate the manifest, contracts, fixtures, existing artifacts,
                           local Mol* build/source readiness, and output plan only.
  --manifest <path>        Reference fixture manifest. Default: ${defaultManifest}
  --artifact-manifest <path>
                           Reference artifact/export metadata. Default: ${defaultArtifactManifest}
  --molstar-dir <path>     Pinned Mol* checkout/build directory. Default: ${defaultMolstarDir}
  --out-dir <path>         Output directory for generated Mol* references. Default: ${defaultOutDir}
  --formats <list>         Comma-separated export formats: json,obj,stl. Default: ${defaultFormats.join(',')}
  --scene-source <source>  Override contract options.scene-source for all fixtures:
                           manual, data-format, or open-files.
  --force-cylinder-impostors
                           Override Mol* cylinder-impostor support checks for
                           headless export diagnostics. This matches WebGL
                           browsers with EXT_frag_depth more closely when
                           headless-gl lacks that extension.
  --runtime-module-dir <path>
                           Extra package prefix or node_modules directory used to resolve
                           runtime-only dependencies such as gl and pngjs. Can be repeated.
  --render-object-report   Print Mol* render-object draw-count diagnostics during conversion.
  --no-build-from-source   Require a prebuilt Mol* CommonJS tree; do not try a local
                           source-tree build fallback before conversion.
  --help, -h               Show this help.

Contract options may set scene-source to manual, data-format, or open-files.
The open-files source replays Mol*'s OpenFiles StateAction path used by the
default Web App drag-and-drop/import handler.

Real conversion uses an already-built local CommonJS Mol* tree at
artifacts/molstar/lib/commonjs. If that tree is missing, the script can build it
from artifacts/molstar/src when local build dependencies are present. Runtime
dependencies are resolved from artifacts/molstar/node_modules plus any
--runtime-module-dir entries. This script never installs or downloads
dependencies.`);
}

function parseArgs(argv) {
  const args = {
    dryRun: false,
    manifest: defaultManifest,
    artifactManifest: defaultArtifactManifest,
    molstarDir: defaultMolstarDir,
    outDir: defaultOutDir,
    formats: [...defaultFormats],
    formatsFromCli: false,
    sceneSource: undefined,
    forceCylinderImpostors: false,
    runtimeModuleDirs: [],
    buildFromSource: true,
    renderObjectReport: false,
    help: false,
  };

  for (let i = 0; i < argv.length; i++) {
    const arg = argv[i];
    if (arg === '--dry-run' || arg === '--check') {
      args.dryRun = true;
    } else if (arg === '--help' || arg === '-h') {
      args.help = true;
    } else if (arg === '--manifest') {
      args.manifest = requireValue(argv, ++i, arg);
    } else if (arg === '--artifact-manifest') {
      args.artifactManifest = requireValue(argv, ++i, arg);
    } else if (arg === '--molstar-dir') {
      args.molstarDir = requireValue(argv, ++i, arg);
    } else if (arg === '--out-dir') {
      args.outDir = requireValue(argv, ++i, arg);
    } else if (arg === '--formats') {
      args.formats = parseFormats(requireValue(argv, ++i, arg), arg);
      args.formatsFromCli = true;
    } else if (arg === '--scene-source') {
      args.sceneSource = validateSceneSource(requireValue(argv, ++i, arg), arg);
    } else if (arg === '--force-cylinder-impostors') {
      args.forceCylinderImpostors = true;
    } else if (arg === '--runtime-module-dir') {
      args.runtimeModuleDirs.push(requireValue(argv, ++i, arg));
    } else if (arg === '--render-object-report') {
      args.renderObjectReport = true;
    } else if (arg === '--no-build-from-source') {
      args.buildFromSource = false;
    } else {
      throw new Error(`Unknown option: ${arg}`);
    }
  }

  validateFormats(args.formats, '--formats');
  return args;
}

function loadArtifactManifest(args) {
  const manifestPath = resolveRepoPath(args.artifactManifest);
  if (!existsSync(manifestPath)) {
    return { path: args.artifactManifest, entriesByContract: new Map(), entriesByName: new Map() };
  }

  const raw = JSON.parse(readFileSync(manifestPath, 'utf8'));
  if (raw.schema !== 1) {
    throw new Error(`${args.artifactManifest}: unsupported schema ${raw.schema}`);
  }
  if (raw.molstar_reference_commit !== molstarReferenceCommit) {
    throw new Error(`${args.artifactManifest}: molstar_reference_commit expected ${molstarReferenceCommit}, got ${raw.molstar_reference_commit}`);
  }
  if (!Array.isArray(raw.artifacts)) {
    throw new Error(`${args.artifactManifest}: expected artifacts array`);
  }

  const entriesByContract = new Map();
  const entriesByName = new Map();
  for (const [index, entry] of raw.artifacts.entries()) {
    const label = `${args.artifactManifest}:artifacts[${index}]`;
    const name = requireString(entry.name, `${label}.name`);
    const contract = requireString(entry.contract, `${label}.contract`);
    if (entriesByName.has(name)) throw new Error(`${label}: duplicate artifact name ${name}`);
    if (entriesByContract.has(contract)) throw new Error(`${label}: duplicate contract ${contract}`);
    entriesByName.set(name, entry);
    entriesByContract.set(contract, entry);
  }

  return { path: args.artifactManifest, entriesByContract, entriesByName };
}

function requireString(value, label) {
  if (typeof value !== 'string' || value.length === 0) {
    throw new Error(`${label}: expected a non-empty string`);
  }
  return value;
}

function requireNonNegativeInteger(value, label) {
  if (!Number.isInteger(value) || value < 0) {
    throw new Error(`${label}: expected a non-negative integer`);
  }
  return value;
}

function parseFormats(value, label) {
  const formats = String(value).split(',').map(s => s.trim()).filter(Boolean);
  validateFormats(formats, label);
  return formats;
}

function validateFormats(formats, label) {
  if (formats.length === 0) {
    throw new Error(`${label}: expected at least one format.`);
  }
  const seen = new Set();
  for (const format of formats) {
    if (format !== 'json' && format !== 'obj' && format !== 'stl') {
      throw new Error(`${label}: unsupported format '${format}'. Expected json, obj, or stl.`);
    }
    if (seen.has(format)) throw new Error(`${label}: duplicate format '${format}'.`);
    seen.add(format);
  }
}

function validateSceneSource(value, label) {
  const sceneSource = normalizeSceneSource(value);
  if (sceneSource !== 'manual' && sceneSource !== 'data-format' && sceneSource !== 'open-files') {
    throw new Error(`${label}: unsupported scene source '${value}'. Expected manual, data-format, or open-files.`);
  }
  return sceneSource;
}

function requireValue(argv, index, flag) {
  const value = argv[index];
  if (!value || value.startsWith('--')) throw new Error(`${flag} requires a value.`);
  return value;
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

function loadPlan(args) {
  const manifestPath = resolveRepoPath(args.manifest);
  if (!existsSync(manifestPath)) throw new Error(`Manifest not found: ${args.manifest}`);

  const manifestEntries = readFileSync(manifestPath, 'utf8')
    .split(/\r?\n/)
    .map((line, index) => parseManifestEntry(line, index + 1, args.manifest))
    .filter(Boolean);

  if (manifestEntries.length === 0) {
    throw new Error(`Manifest has no contract entries: ${args.manifest}`);
  }

  return manifestEntries.map(entry => {
    const contractPath = entry.contractPath;
    const absContractPath = resolveRepoPath(contractPath);
    if (!existsSync(absContractPath)) throw new Error(`Contract not found: ${contractPath}`);
    const contract = parseContract(readFileSync(absContractPath, 'utf8'), contractPath);
    assertEqual(String(contract.molstar_reference_commit ?? ''), molstarReferenceCommit, `${contractPath}: molstar_reference_commit`);
    const fixture = String(contract.fixture ?? '');
    if (!fixture) throw new Error(`${contractPath}: missing fixture=...`);
    const absFixturePath = resolveRepoPath(fixture);
    if (!existsSync(absFixturePath)) throw new Error(`${contractPath}: fixture not found: ${fixture}`);
    if (!statSync(absFixturePath).isFile()) throw new Error(`${contractPath}: fixture is not a file: ${fixture}`);

    let options = {};
    if (contract.options) {
      try {
        options = JSON.parse(contract.options);
      } catch (error) {
        throw new Error(`${contractPath}: invalid options JSON: ${error.message}`);
      }
    }
    if (args.sceneSource !== undefined) {
      options['scene-source'] = args.sceneSource;
    }
    if (args.forceCylinderImpostors) {
      options['force-cylinder-impostors'] = true;
    }

    const inputFormat = normalizeInputFormat(options.format ?? contract.input_format, fixture);
    if (contract.input_format) {
      const declaredInputFormat = normalizeInputFormat(contract.input_format, fixture);
      if (declaredInputFormat !== inputFormat) {
        throw new Error(`${contractPath}: input_format=${contract.input_format} disagrees with options.format=${options.format}`);
      }
    }
    const stem = contractOutputStem(contractPath);
    const objReference = String(contract.obj_reference ?? '');
    const objExportBasename = contractObjExportBasename(contract, objReference, stem, contractPath);
    const objMtllib = contractObjMtllib(contract, objExportBasename, contractPath);
    return {
      contractPath,
      contract,
      fixture,
      absFixturePath,
      options,
      inputFormat,
      formats: args.formatsFromCli ? args.formats : entry.formats ?? args.formats,
      stem,
      objReference,
      objExportBasename,
      objMtllib,
      stlReference: String(contract.stl_reference ?? ''),
      jsonReference: String(contract.json_reference ?? ''),
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
      if (contractPath) {
        throw new Error(`${label}:${lineNumber}: unexpected manifest token '${field}'`);
      }
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
      formats = parseFormats(value, `${label}:${lineNumber}: formats`);
    } else if (key === 'tag' || key === 'tags' || key === 'note') {
      // Reserved manifest metadata; exporter planning does not need it yet.
    } else {
      throw new Error(`${label}:${lineNumber}: unknown manifest field '${key}'`);
    }
  }

  if (!contractPath) throw new Error(`${label}:${lineNumber}: missing contract path`);
  return { contractPath, formats };
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

function contractObjExportBasename(contract, objReference, stem, contractPath) {
  const basename = String(contract.obj_export_basename ?? '').trim();
  if (basename) return validateBasename(basename, `${contractPath}: obj_export_basename`);
  if (objReference) {
    const referenceBase = path.basename(objReference, path.extname(objReference));
    if (referenceBase) return validateBasename(referenceBase, `${contractPath}: obj_reference basename`);
  }
  return validateBasename(stem, `${contractPath}: derived obj export basename`);
}

function validateBasename(value, label) {
  if (value !== path.basename(value) || value.includes('/') || value.includes('\\')) {
    throw new Error(`${label}: expected a filename stem without path separators, got ${value}`);
  }
  return value;
}

function contractObjMtllib(contract, objExportBasename, contractPath) {
  const mtllib = String(contract.obj_mtllib ?? `${objExportBasename}.mtl`).trim();
  if (!mtllib) throw new Error(`${contractPath}: obj_mtllib must not be empty`);
  if (mtllib !== path.basename(mtllib) || mtllib.includes('/') || mtllib.includes('\\')) {
    throw new Error(`${contractPath}: obj_mtllib must be a filename without path separators, got ${mtllib}`);
  }
  return mtllib;
}

function inspectMolstar(args) {
  const molstarDir = resolveRepoPath(args.molstarDir);
  const commonjsDir = path.join(molstarDir, 'lib/commonjs');
  const packageJsonPath = path.join(molstarDir, 'package.json');
  const packageLockPath = path.join(molstarDir, 'package-lock.json');
  const nodeModulesPath = path.join(molstarDir, 'node_modules');
  const missing = [];

  if (!existsSync(molstarDir)) missing.push(`${args.molstarDir}/`);
  if (!existsSync(commonjsDir)) missing.push(`${args.molstarDir}/lib/commonjs/`);

  const requiredFiles = [
    'mol-plugin/headless-plugin-context.js',
    'mol-plugin/spec.js',
    'mol-plugin-state/actions/file.js',
    'mol-plugin-state/transforms/data.js',
    'mol-plugin-state/transforms/model.js',
    'mol-plugin-state/transforms/representation.js',
    'extensions/geo-export/obj-exporter.js',
    'extensions/geo-export/stl-exporter.js',
    'mol-repr/structure/visual/util/common.js',
    'mol-math/geometry.js',
    'mol-task/index.js',
    'mol-util/data-source.js',
    'mol-util/assets.js',
    'mol-util/nodejs-shims.js',
    'mol-util/zip/zip.js',
  ];

  for (const file of requiredFiles) {
    const candidate = path.join(commonjsDir, file);
    if (!existsSync(candidate)) missing.push(`${args.molstarDir}/lib/commonjs/${file}`);
  }

  const sourceFiles = [
    'mol-plugin/headless-plugin-context.ts',
    'mol-plugin/spec.ts',
    'mol-plugin-state/actions/file.ts',
    'mol-plugin-state/transforms/data.ts',
    'mol-plugin-state/transforms/model.ts',
    'mol-plugin-state/transforms/representation.ts',
    'extensions/geo-export/obj-exporter.ts',
    'extensions/geo-export/stl-exporter.ts',
    'mol-repr/structure/visual/util/common.ts',
    'mol-math/geometry.ts',
    'mol-task/index.ts',
    'mol-util/data-source.ts',
    'mol-util/assets.ts',
    'mol-util/nodejs-shims.ts',
    'mol-util/zip/zip.ts',
  ].map(file => ({
    file,
    present: existsSync(path.join(molstarDir, 'src', file)),
  }));

  let packageJson = {};
  let packageLock = {};
  if (existsSync(packageJsonPath)) {
    packageJson = JSON.parse(readFileSync(packageJsonPath, 'utf8'));
  }
  if (existsSync(packageLockPath)) {
    packageLock = JSON.parse(readFileSync(packageLockPath, 'utf8'));
  }
  const dependencyResolvers = runtimeDependencyResolvers(args, packageJsonPath);
  const resolveDependency = name => {
    for (const resolver of dependencyResolvers) {
      try {
        return {
          present: true,
          resolvedBy: resolver.label,
          path: resolver.require.resolve(name),
        };
      } catch {
        // Try the next configured module root.
      }
    }
    return { present: false, resolvedBy: '', path: '' };
  };
  const npmCachePath = inspectNpmCachePath();
  const dependencyState = name => {
    const resolved = resolveDependency(name);
    return {
      name,
      present: resolved.present,
      resolvedBy: resolved.resolvedBy,
      path: resolved.path,
      lockPresent: Boolean(packageLock.packages?.[`node_modules/${name}`]),
      cache: inspectNpmCachePackage(name),
    };
  };
  const runtimeDependencies = ['gl', 'pngjs', 'jpeg-js'].map(dependencyState);
  const buildDependencies = ['concurrently', 'typescript', 'tsc-alias', 'cpx2'].map(dependencyState);
  const missingSources = sourceFiles.filter(file => !file.present);
  const missingBuildDependencies = buildDependencies.filter(dep => !dep.present);
  const missingRuntimeDependencies = runtimeDependencies.filter(dep => !dep.present);

  return {
    molstarDir,
    commonjsDir,
    missing,
    commonjsReady: missing.length === 0,
    packageVersion: packageJson.version,
    packageLockPresent: existsSync(packageLockPath),
    nodeModulesPresent: existsSync(nodeModulesPath),
    runtimeModuleDirs: args.runtimeModuleDirs.map(resolveInputPath),
    toolchain: inspectToolchain(),
    npmCachePath,
    tsconfigCommonjsPresent: existsSync(path.join(molstarDir, 'tsconfig.commonjs.json')),
    buildScript: packageJson.scripts?.['build:lib'],
    sourceFiles,
    missingSources,
    runtimeDependencies,
    buildDependencies,
    missingBuildDependencies,
    missingRuntimeDependencies,
    sourceTreeBuildPossible: existsSync(packageJsonPath)
      && existsSync(path.join(molstarDir, 'tsconfig.commonjs.json'))
      && missingSources.length === 0
      && missingBuildDependencies.length === 0,
  };
}

function runtimeDependencyResolvers(args, packageJsonPath) {
  const resolvers = [];
  if (existsSync(packageJsonPath)) {
    resolvers.push({
      label: displayPath(path.dirname(packageJsonPath)),
      require: createRequire(packageJsonPath),
    });
  }
  for (const dir of args.runtimeModuleDirs) {
    const absDir = resolveInputPath(dir);
    const base = path.basename(absDir) === 'node_modules'
      ? path.join(path.dirname(absDir), 'package.json')
      : path.join(absDir, 'package.json');
    resolvers.push({
      label: displayPath(absDir),
      require: createRequire(base),
    });
  }
  return resolvers;
}

function inspectNpmCachePath() {
  const result = spawnSync('npm', ['config', 'get', 'cache'], {
    encoding: 'utf8',
    timeout: 2500,
  });
  if (result.status !== 0) return { checked: false, error: formatProcessOutput(result) };
  return { checked: true, path: result.stdout.trim() };
}

function inspectToolchain() {
  const pythonCommands = [];
  if (process.env.PYTHON) pythonCommands.push(process.env.PYTHON);
  pythonCommands.push('python', 'python3');
  return {
    node: {
      executable: process.execPath,
      version: process.version,
    },
    npm: commandProbe('npm', ['--version']),
    python: pythonCommands.map(command => ({
      command,
      version: commandProbe(command, ['--version']),
      distutils: commandProbe(command, ['-c', 'import distutils.version; print("distutils-ok")']),
    })),
  };
}

function commandProbe(command, args) {
  const result = spawnSync(command, args, {
    encoding: 'utf8',
    timeout: 2500,
  });
  if (result.error) return { ok: false, detail: result.error.message };
  if (result.status !== 0) return { ok: false, detail: formatProcessOutput(result) };
  const output = `${result.stdout || ''}${result.stderr || ''}`.trim();
  return { ok: true, detail: output };
}

function inspectNpmCachePackage(name) {
  const result = spawnSync('npm', ['cache', 'ls', '--parseable', name], {
    encoding: 'utf8',
    timeout: 2500,
  });
  if (result.error) return { checked: false, present: false, error: result.error.message };
  if (result.status !== 0) return { checked: false, present: false, error: formatProcessOutput(result) };
  return { checked: true, present: result.stdout.trim().length > 0 };
}

function printPlan(args, plan, inspection, artifactManifest) {
  console.log(`Mol* reference conversion ${args.dryRun ? 'dry run' : 'plan'}`);
  console.log(`manifest: ${args.manifest}`);
  console.log(`artifact manifest: ${artifactManifest.entriesByContract.size === 0 ? `${artifactManifest.path} (not present)` : artifactManifest.path}`);
  console.log(`fixtures: ${plan.length}`);
  for (const item of plan) {
    const paths = plannedOutputPaths(args, item);
    const sceneSource = normalizeSceneSource(item.options['scene-source'] ?? item.options.scene_source ?? 'manual');
    console.log(`- ${item.fixture} input=${item.inputFormat} scene=${sceneSource} targets=${item.formats.join(',')} stem=${item.stem} obj-basename=${item.objExportBasename}`);
    console.log(`  exports obj=${path.relative(rootDir, paths.objPath)} mtl=${path.relative(rootDir, paths.mtlPath)} stl=${path.relative(rootDir, paths.stlPath)} json=${path.relative(rootDir, paths.jsonPath)}`);
  }
  if (inspection.missing.length) {
    console.log('');
    console.log('Mol* CommonJS prebuilt tree is not ready:');
    for (const missing of inspection.missing) console.log(`- missing ${missing}`);
    console.log('');
    console.log('Observed local Mol* checkout state:');
    console.log(`- package version: ${inspection.packageVersion ?? '<missing package.json>'}`);
    console.log(`- package-lock.json: ${inspection.packageLockPresent ? 'present' : 'missing'}`);
    console.log(`- node_modules: ${inspection.nodeModulesPresent ? 'present' : 'missing'}`);
    if (inspection.runtimeModuleDirs.length) {
      console.log(`- runtime module dirs: ${inspection.runtimeModuleDirs.map(displayPath).join(', ')}`);
    }
    printToolchainStatus(inspection.toolchain);
    console.log(`- npm cache: ${formatNpmCachePath(inspection.npmCachePath)}`);
    console.log(`- tsconfig.commonjs.json: ${inspection.tsconfigCommonjsPresent ? 'present' : 'missing'}`);
    console.log(`- build:lib script: ${inspection.buildScript ?? '<missing>'}`);
    if (inspection.missingSources.length === 0) {
      console.log('- TypeScript source entrypoints: present');
    } else {
      console.log('- TypeScript source entrypoints missing:');
      for (const source of inspection.missingSources) console.log(`  - src/${source.file}`);
    }
    printDependencyStatus('build dependencies', inspection.buildDependencies);
    printDependencyStatus('runtime dependencies', inspection.runtimeDependencies);
    console.log('');
    if (inspection.sourceTreeBuildPossible && inspection.missingRuntimeDependencies.length === 0 && args.buildFromSource) {
      console.log('Source-tree fallback: available. A real conversion will run');
      console.log(`npm --offline run build:lib in ${args.molstarDir} before loading Mol*.`);
    } else {
      console.log('Source-tree fallback: blocked.');
      for (const reason of sourceTreeBuildBlockers(args, inspection)) {
        console.log(`- ${reason}`);
      }
    }
  } else {
    console.log(`Mol* CommonJS prebuilt tree: ${displayPath(inspection.commonjsDir)}`);
    if (inspection.runtimeModuleDirs.length) {
      console.log(`runtime module dirs: ${inspection.runtimeModuleDirs.map(displayPath).join(', ')}`);
    }
    printToolchainStatus(inspection.toolchain);
    printDependencyStatus('runtime dependencies', inspection.runtimeDependencies);
  }
}

function printToolchainStatus(toolchain) {
  console.log(`- node: ${toolchain.node.version} (${toolchain.node.executable})`);
  console.log(`- npm: ${formatProbe(toolchain.npm)}`);
  const seen = new Set();
  for (const python of toolchain.python) {
    if (seen.has(python.command)) continue;
    seen.add(python.command);
    console.log(`- ${python.command}: ${formatProbe(python.version)}; distutils=${formatProbe(python.distutils)}`);
  }
}

function formatProbe(probe) {
  return probe.ok ? probe.detail : `unavailable (${probe.detail})`;
}

function printDependencyStatus(label, dependencies) {
  const missing = dependencies.filter(dep => !dep.present);
  if (missing.length === 0) {
    console.log(`- ${label}: installed`);
    for (const dep of dependencies) {
      if (dep.resolvedBy) console.log(`  - ${dep.name}: ${dep.resolvedBy}`);
    }
    return;
  }

  console.log(`- ${label}: missing ${missing.map(dep => dep.name).join(', ')}`);
  for (const dep of missing) {
    console.log(`  - ${dep.name}: package-lock=${dep.lockPresent ? 'present' : 'missing'}, npm-cache=${formatNpmCachePackage(dep.cache)}`);
  }
  for (const dep of dependencies.filter(dep => dep.present)) {
    if (dep.resolvedBy) console.log(`  - ${dep.name}: ${dep.resolvedBy}`);
  }
}

function formatNpmCachePath(cache) {
  if (!cache.checked) return `unavailable (${cache.error || 'npm cache inspection failed'})`;
  return cache.path || '<empty>';
}

function formatNpmCachePackage(cache) {
  if (!cache.checked) return `unknown (${cache.error || 'inspection failed'})`;
  return cache.present ? 'present' : 'missing';
}

function sourceTreeBuildBlockers(args, inspection) {
  const blockers = [];
  if (!inspection.commonjsReady) {
    if (!args.buildFromSource) blockers.push('--no-build-from-source was set');
    if (!inspection.packageVersion) blockers.push(`${args.molstarDir}/package.json is missing or unreadable`);
    if (!inspection.tsconfigCommonjsPresent) blockers.push(`${args.molstarDir}/tsconfig.commonjs.json is missing`);
    for (const source of inspection.missingSources) blockers.push(`${args.molstarDir}/src/${source.file} is missing`);
    for (const dep of inspection.missingBuildDependencies) {
      blockers.push(`build dependency '${dep.name}' is not installed in ${args.molstarDir}/node_modules`);
    }
  }
  const runtimeSearchRoots = [
    `${args.molstarDir}/node_modules`,
    ...args.runtimeModuleDirs.map(dir => {
      const absDir = resolveInputPath(dir);
      return path.basename(absDir) === 'node_modules'
        ? displayPath(absDir)
        : `${displayPath(absDir)}/node_modules`;
    }),
  ].join(' or ');
  for (const dep of inspection.missingRuntimeDependencies) {
    blockers.push(`runtime dependency '${dep.name}' is not resolved from ${runtimeSearchRoots}`);
  }
  if (blockers.length === 0 && inspection.commonjsReady) {
    blockers.push('prebuilt CommonJS tree is already available');
  }
  return blockers;
}

function validateExistingReferences(args, plan, artifactManifest) {
  console.log('');
  console.log('Existing Mol* reference artifact validation:');
  let checked = 0;
  const matchedArtifactContracts = new Set();

  for (const item of plan) {
    const artifactMetadata = validateArtifactManifestEntry(args, item, artifactManifest);
    if (artifactMetadata) matchedArtifactContracts.add(artifactMetadata.contract);
    const stats = {};

    if (item.formats.includes('obj') || item.formats.includes('json')) {
      stats.obj = validateObjReference(item, artifactMetadata);
    }

    if (item.formats.includes('stl') || item.formats.includes('json')) {
      stats.stl = validateStlReference(item);
    }

    if (stats.obj && stats.stl) {
      validateObjStlSparseSlotDriftReference(item, stats.obj, stats.stl);
    }

    if (item.formats.includes('json')) {
      validateJsonReference(item, stats);
    }

    console.log(`- PASS ${item.stem}: ${referencePassSummary(item, stats)}`);
    checked++;
  }

  if (checked === 0) throw new Error('Reference manifest did not contain any fixtures.');
  for (const contractPath of artifactManifest.entriesByContract.keys()) {
    if (!matchedArtifactContracts.has(contractPath)) {
      throw new Error(`${artifactManifest.path}: artifact contract is not present in ${args.manifest}: ${contractPath}`);
    }
  }
}

function referencePassSummary(item, stats) {
  const parts = [item.formats.join('/')];
  if (stats.obj) {
    parts.push(`OBJ ${stats.obj.vertex_count}v/${stats.obj.normal_count}vn/${stats.obj.face_count}f ${stats.obj.fnv1a64}`);
  }
  if (stats.stl) {
    parts.push(`STL ${stats.stl.facet_count} facets ${stats.stl.fnv1a64}`);
  }
  if (item.formats.includes('json')) {
    parts.push(`JSON ${contractValue(item, 'json_fnv1a64')}`);
  }
  return parts.join('; ');
}

function validateArtifactManifestEntry(args, item, artifactManifest) {
  const entry = artifactManifest.entriesByContract.get(item.contractPath)
    ?? artifactManifest.entriesByName.get(item.stem);
  if (!entry) return undefined;

  assertEqual(entry.name, item.stem, `${artifactManifest.path}: ${item.contractPath}: name`);
  assertEqual(entry.contract, item.contractPath, `${artifactManifest.path}: ${item.contractPath}: contract`);
  assertEqual(entry.fixture, item.fixture, `${artifactManifest.path}: ${item.contractPath}: fixture`);

  const exports = entry.exports ?? {};
  const paths = plannedOutputPaths(args, item);
  validateArtifactExportPath(args, item, artifactManifest, exports, 'obj', item.objReference, paths.objPath);
  validateArtifactExportPath(args, item, artifactManifest, exports, 'stl', item.stlReference, paths.stlPath);
  validateArtifactExportPath(args, item, artifactManifest, exports, 'json', item.jsonReference, paths.jsonPath);
  validateArtifactExportPath(args, item, artifactManifest, exports, 'mtl', path.join(path.dirname(item.objReference), item.objMtllib), paths.mtlPath, item.formats.includes('obj'));

  const obj = entry.obj ?? {};
  if (item.formats.includes('obj') || obj.export_basename !== undefined) {
    assertEqual(requireString(obj.export_basename, `${artifactManifest.path}: ${item.stem}.obj.export_basename`), item.objExportBasename, `${artifactManifest.path}: ${item.contractPath}: obj.export_basename`);
  }
  if (item.formats.includes('obj') || obj.mtllib !== undefined) {
    assertEqual(requireString(obj.mtllib, `${artifactManifest.path}: ${item.stem}.obj.mtllib`), item.objMtllib, `${artifactManifest.path}: ${item.contractPath}: obj.mtllib`);
  }

  const mtl = entry.mtl ?? {};
  if (mtl.byte_len !== undefined) {
    assertEqual(requireNonNegativeInteger(mtl.byte_len, `${artifactManifest.path}: ${item.stem}.mtl.byte_len`), contractUsize(item, 'mtl_byte_len'), `${artifactManifest.path}: ${item.contractPath}: mtl.byte_len`);
  }
  if (mtl.fnv1a64 !== undefined) {
    assertEqual(requireString(mtl.fnv1a64, `${artifactManifest.path}: ${item.stem}.mtl.fnv1a64`), contractValue(item, 'mtl_fnv1a64'), `${artifactManifest.path}: ${item.contractPath}: mtl.fnv1a64`);
  }
  if (mtl.material_count !== undefined) {
    assertEqual(requireNonNegativeInteger(mtl.material_count, `${artifactManifest.path}: ${item.stem}.mtl.material_count`), contractUsize(item, 'mtl_material_count'), `${artifactManifest.path}: ${item.contractPath}: mtl.material_count`);
  }
  if (mtl.materials !== undefined) {
    assertStringArray(mtl.materials, `${artifactManifest.path}: ${item.stem}.mtl.materials`);
    assertStringArray(obj.materials, `${artifactManifest.path}: ${item.stem}.obj.materials`);
    assertEqual(JSON.stringify(mtl.materials), JSON.stringify(obj.materials), `${artifactManifest.path}: ${item.contractPath}: mtl.materials`);
  }

  const stl = entry.stl ?? {};
  if (stl.header !== undefined) {
    assertEqual(requireString(stl.header, `${artifactManifest.path}: ${item.stem}.stl.header`), contractValue(item, 'stl_header'), `${artifactManifest.path}: ${item.contractPath}: stl.header`);
  }
  if (stl.nonzero_attribute_count !== undefined) {
    assertEqual(requireNonNegativeInteger(stl.nonzero_attribute_count, `${artifactManifest.path}: ${item.stem}.stl.nonzero_attribute_count`), contractUsize(item, 'stl_nonzero_attribute_count'), `${artifactManifest.path}: ${item.contractPath}: stl.nonzero_attribute_count`);
  }

  return entry;
}

function validateArtifactExportPath(args, item, artifactManifest, exports, format, expectedPath, plannedPath, required = item.formats.includes(format)) {
  const actualPath = exports[format];
  if (!required && actualPath === undefined) return;
  assertEqual(requireString(actualPath, `${artifactManifest.path}: ${item.stem}.exports.${format}`), expectedPath, `${artifactManifest.path}: ${item.contractPath}: exports.${format}`);
  assertEqual(path.relative(rootDir, plannedPath), actualPath, `${artifactManifest.path}: ${item.contractPath}: planned ${format} export`);
}

function validateObjReference(item, artifactMetadata) {
  const objPath = contractValue(item, 'obj_reference');
  const obj = readFileSync(resolveRepoPath(objPath), 'utf8');
  const objBytes = Buffer.from(obj);
  const objStats = objExportStats(obj);
  failIfObjStatsInvalid(item, objStats);
  if (objStats.mtllibs.length !== 1) {
    throw new Error(`${item.contractPath}: expected exactly one OBJ mtllib line, got ${objStats.mtllibs.length}`);
  }
  assertEqual(objStats.firstContentLine, `mtllib ${item.objMtllib}`, `${item.contractPath}: first OBJ content line`);
  assertEqual(objStats.mtllib, `mtllib ${item.objMtllib}`, `${item.contractPath}: obj_mtllib`);
  validateObjMaterials(item, objPath, objStats, artifactMetadata);
  assertEqual(objBytes.length, contractUsize(item, 'obj_byte_len'), `${item.contractPath}: obj_byte_len`);
  assertEqual(objStats.fnv1a64, contractValue(item, 'obj_fnv1a64'), `${item.contractPath}: obj_fnv1a64`);
  assertEqual(objStats.vertex_count, contractUsize(item, 'obj_vertex_count'), `${item.contractPath}: obj_vertex_count`);
  assertEqual(objStats.normal_count, contractUsize(item, 'obj_normal_count'), `${item.contractPath}: obj_normal_count`);
  assertEqual(objStats.face_count, contractUsize(item, 'obj_face_count'), `${item.contractPath}: obj_face_count`);
  assertVec3Close(objStats.min, contractVec3(item, 'obj_bounds_min'), 0.0001, `${item.contractPath}: obj_bounds_min`);
  assertVec3Close(objStats.max, contractVec3(item, 'obj_bounds_max'), 0.0001, `${item.contractPath}: obj_bounds_max`);
  return objStats;
}

function failIfObjStatsInvalid(item, objStats) {
  if (objStats.unsupportedLines.length) {
    throw new Error(`${item.contractPath}: unsupported OBJ lines: ${formatLineExamples(objStats.unsupportedLines)}`);
  }
  if (objStats.invalidVectors.length) {
    throw new Error(`${item.contractPath}: invalid OBJ vector lines: ${formatLineExamples(objStats.invalidVectors)}`);
  }
  if (objStats.invalidFaces.length) {
    throw new Error(`${item.contractPath}: invalid OBJ face lines: ${formatLineExamples(objStats.invalidFaces)}`);
  }
  if (objStats.nonTriFaceCount !== 0) {
    throw new Error(`${item.contractPath}: expected only triangle OBJ faces, got ${objStats.nonTriFaceCount} non-triangle faces`);
  }
  if (objStats.faceVertexNormalMismatchCount !== 0) {
    throw new Error(`${item.contractPath}: expected OBJ face vertex and normal indices to match, got ${objStats.faceVertexNormalMismatchCount} mismatches`);
  }
  if (objStats.faceIndexErrors.length) {
    throw new Error(`${item.contractPath}: OBJ face index out of range: ${formatLineExamples(objStats.faceIndexErrors)}`);
  }
  if (objStats.vertex_count === 0 || objStats.normal_count === 0 || objStats.face_count === 0) {
    throw new Error(`${item.contractPath}: expected non-empty OBJ vertices, normals, and faces`);
  }
}

function validateObjMaterials(item, objPath, objStats, artifactMetadata) {
  const metadataObj = artifactMetadata?.obj;
  if (metadataObj) {
    assertStringArray(metadataObj.materials, `${item.contractPath}: obj.materials`);
    assertEqual(JSON.stringify(objStats.materials), JSON.stringify(metadataObj.materials), `${item.contractPath}: obj.materials`);
    if (metadataObj.material_switch_count !== undefined) {
      assertEqual(objStats.usemtl.length, metadataObj.material_switch_count, `${item.contractPath}: obj.material_switch_count`);
    }
  }

  const mtlReference = path.join(path.dirname(objPath), item.objMtllib);
  const absMtlReference = resolveRepoPath(mtlReference);
  if (!existsSync(absMtlReference)) {
    throw new Error(`${item.contractPath}: OBJ mtllib references missing material file: ${mtlReference}`);
  }

  const mtl = readFileSync(absMtlReference, 'utf8');
  const mtlStats = mtlExportStats(mtl);
  assertEqual(mtlStats.byte_len, contractUsize(item, 'mtl_byte_len'), `${item.contractPath}: mtl_byte_len`);
  assertEqual(mtlStats.fnv1a64, contractValue(item, 'mtl_fnv1a64'), `${item.contractPath}: mtl_fnv1a64`);
  assertEqual(mtlStats.materials.length, contractUsize(item, 'mtl_material_count'), `${item.contractPath}: mtl_material_count`);
  assertEqual(JSON.stringify(mtlStats.materials), JSON.stringify(objStats.materials), `${item.contractPath}: MTL newmtl materials`);
  if (artifactMetadata?.mtl) {
    if (artifactMetadata.mtl.materials !== undefined) {
      assertEqual(JSON.stringify(mtlStats.materials), JSON.stringify(artifactMetadata.mtl.materials), `${item.contractPath}: artifact MTL materials`);
    }
  }
}

function assertStringArray(value, label) {
  if (!Array.isArray(value) || value.some(item => typeof item !== 'string')) {
    throw new Error(`${label}: expected an array of strings`);
  }
}

function validateStlReference(item) {
  const stlPath = contractValue(item, 'stl_reference');
  const stl = readFileSync(resolveRepoPath(stlPath));
  const stlStats = stlExportStats(stl);
  if (stlStats.invalidFloatCount !== 0) {
    throw new Error(`${item.contractPath}: STL contains ${stlStats.invalidFloatCount} non-finite float values`);
  }
  assertEqual(stlStats.byte_len, contractUsize(item, 'stl_byte_len'), `${item.contractPath}: stl_byte_len`);
  assertEqual(stlStats.fnv1a64, contractValue(item, 'stl_fnv1a64'), `${item.contractPath}: stl_fnv1a64`);
  assertEqual(stlStats.facet_count, contractUsize(item, 'stl_facet_count'), `${item.contractPath}: stl_facet_count`);
  assertEqual(stlStats.header, contractValue(item, 'stl_header'), `${item.contractPath}: stl_header`);
  assertEqual(stlStats.nonzero_attribute_count, contractUsize(item, 'stl_nonzero_attribute_count'), `${item.contractPath}: stl_nonzero_attribute_count`);
  assertVec3Close(stlStats.min, contractVec3(item, 'stl_bounds_min'), 0.0001, `${item.contractPath}: stl_bounds_min`);
  assertVec3Close(stlStats.max, contractVec3(item, 'stl_bounds_max'), 0.0001, `${item.contractPath}: stl_bounds_max`);
  return stlStats;
}

function validateObjStlSparseSlotDriftReference(item, objStats, stlStats) {
  if (item.contract.obj_stl_sparse_slot_total_components === undefined) return;

  assertEqual(stlStats.facet_count, objStats.face_count * 3, `${item.contractPath}: OBJ face count to STL sparse slot count`);
  const drift = objStlSparseSlotDrift(objStats.vertices, objStats.faces, stlStats.bytes);
  assertEqual(drift.total_components, contractUsize(item, 'obj_stl_sparse_slot_total_components'), `${item.contractPath}: obj_stl_sparse_slot_total_components`);
  assertEqual(drift.rounded_mismatch_count, contractUsize(item, 'obj_stl_sparse_slot_rounding_mismatch_count'), `${item.contractPath}: obj_stl_sparse_slot_rounding_mismatch_count`);
  assertNumberClose(drift.max_delta, contractNumber(item, 'obj_stl_sparse_slot_max_delta'), 0.0000001, `${item.contractPath}: obj_stl_sparse_slot_max_delta`);

  const first = drift.first_rounded_mismatch;
  if (!first) {
    throw new Error(`${item.contractPath}: obj_stl_sparse_slot_first_mismatch is pinned but no OBJ/STL sparse-slot mismatch was found`);
  }
  const expectedFirst = contractUsizeTuple(item, 'obj_stl_sparse_slot_first_mismatch', 4);
  assertEqual(JSON.stringify([first.face_index, first.stl_facet_index, first.vertex_slot, first.axis]), JSON.stringify(expectedFirst), `${item.contractPath}: obj_stl_sparse_slot_first_mismatch`);
}

function validateJsonReference(item, stats) {
  const jsonPath = contractValue(item, 'json_reference');
  const json = readFileSync(resolveRepoPath(jsonPath), 'utf8');
  const jsonBytes = Buffer.from(json);
  const parsed = parseJsonReference(json, item);
  assertEqual(parsed.schema, 1, `${item.contractPath}: JSON schema`);
  assertEqual(parsed.name, item.stem, `${item.contractPath}: JSON name`);
  assertEqual(parsed.molstar_reference_commit, molstarReferenceCommit, `${item.contractPath}: JSON molstar_reference_commit`);
  assertEqual(parsed.fixture, item.fixture, `${item.contractPath}: JSON fixture`);
  assertEqual(parsed.options?.format, item.inputFormat, `${item.contractPath}: JSON options.format`);
  assertEqual(parsed.options?.representation, String(item.options.representation ?? 'molstar'), `${item.contractPath}: JSON options.representation`);
  assertEqual(parsed.options?.assembly, String(item.options.assembly ?? '1'), `${item.contractPath}: JSON options.assembly`);
  assertEqual(parsed.options?.['sphere-detail'], Number(item.options['sphere-detail'] ?? item.options.sphere_detail ?? 1), `${item.contractPath}: JSON options.sphere-detail`);
  assertEqual(parsed.source?.obj, item.objReference, `${item.contractPath}: JSON source.obj`);
  assertEqual(parsed.source?.stl, item.stlReference, `${item.contractPath}: JSON source.stl`);
  assertEqual(jsonBytes.length, contractUsize(item, 'json_byte_len'), `${item.contractPath}: json_byte_len`);
  assertEqual(fnv1a64(jsonBytes), contractValue(item, 'json_fnv1a64'), `${item.contractPath}: json_fnv1a64`);
  const generatedSummary = `${referenceSummaryJson(item, stats)}\n`;
  const jsonDiff = diffText(json, generatedSummary, `${item.contractPath}: json summary`);
  if (jsonDiff) throw new Error(jsonDiff);
}

function parseJsonReference(json, item) {
  try {
    return JSON.parse(json);
  } catch (error) {
    throw new Error(`${item.contractPath}: invalid JSON summary: ${error.message}`);
  }
}

function contractValue(item, key) {
  const value = item.contract[key];
  if (value === undefined) throw new Error(`${item.contractPath}: missing ${key}=...`);
  return String(value);
}

function contractUsize(item, key) {
  const value = Number(contractValue(item, key));
  if (!Number.isInteger(value) || value < 0) {
    throw new Error(`${item.contractPath}: invalid integer ${key}=${contractValue(item, key)}`);
  }
  return value;
}

function contractNumber(item, key) {
  const value = Number(contractValue(item, key));
  if (!Number.isFinite(value)) {
    throw new Error(`${item.contractPath}: invalid number ${key}=${contractValue(item, key)}`);
  }
  return value;
}

function contractUsizeTuple(item, key, count) {
  const values = contractValue(item, key).split(',').map(value => Number(value));
  if (values.length !== count || values.some(value => !Number.isInteger(value) || value < 0)) {
    throw new Error(`${item.contractPath}: invalid integer tuple ${key}=${contractValue(item, key)}`);
  }
  return values;
}

function contractVec3(item, key) {
  const values = contractValue(item, key).split(',').map(value => Number(value));
  if (values.length !== 3 || values.some(value => !Number.isFinite(value))) {
    throw new Error(`${item.contractPath}: invalid vec3 ${key}=${contractValue(item, key)}`);
  }
  return values;
}

function assertEqual(actual, expected, label) {
  if (actual !== expected) throw new Error(`${label}: expected ${expected}, got ${actual}`);
}

function assertVec3Close(actual, expected, tolerance, label) {
  for (let axis = 0; axis < 3; axis++) {
    if (Math.abs(actual[axis] - expected[axis]) > tolerance) {
      throw new Error(`${label}[${axis}]: expected ${expected[axis]}, got ${actual[axis]}, tolerance ${tolerance}`);
    }
  }
}

function assertNumberClose(actual, expected, tolerance, label) {
  if (Math.abs(actual - expected) > tolerance) {
    throw new Error(`${label}: expected ${expected}, got ${actual}, tolerance ${tolerance}`);
  }
}

function diffText(reference, generated, label) {
  if (reference === generated) return undefined;
  const referenceLines = reference.split(/\r?\n/);
  const generatedLines = generated.split(/\r?\n/);
  const limit = Math.min(referenceLines.length, generatedLines.length);
  let firstDiff = 0;
  while (firstDiff < limit && referenceLines[firstDiff] === generatedLines[firstDiff]) firstDiff++;
  const referenceLine = firstDiff < referenceLines.length ? JSON.stringify(referenceLines[firstDiff]) : '<eof>';
  const generatedLine = firstDiff < generatedLines.length ? JSON.stringify(generatedLines[firstDiff]) : '<eof>';
  return `FAIL ${label}: first difference at line ${firstDiff + 1}; reference=${referenceLine}, generated=${generatedLine}; reference_lines=${referenceLines.length}, generated_lines=${generatedLines.length}`;
}

function loadMolstar(args, inspection) {
  inspection = ensureMolstarCommonjs(args, inspection);

  if (inspection.missingRuntimeDependencies.length) {
    throw new Error(conversionBlockedMessage(args, inspection));
  }

  installHeadlessDomShim();

  const packageJsonPath = path.join(inspection.molstarDir, 'package.json');
  installMolstarVersionShim(packageJsonPath);
  const molstarRequire = createRequire(packageJsonPath);
  const dependencyResolvers = runtimeDependencyResolvers(args, packageJsonPath);
  const fromCommonjs = rel => molstarRequire(path.join(inspection.commonjsDir, rel));
  const loadDependency = name => {
    const errors = [];
    for (const resolver of dependencyResolvers) {
      try {
        return resolver.require(name);
      } catch (error) {
        errors.push(`${resolver.label}: ${error.message}`);
      }
    }
    throw new Error(`Missing Mol* runtime dependency '${name}'. Install it in ${args.molstarDir} or pass --runtime-module-dir. (${errors.join(' | ')})`);
  };

  const fsModule = molstarRequire('node:fs');
  const gl = loadDependency('gl');
  const pngjs = loadDependency('pngjs');
  const jpegjs = loadDependency('jpeg-js');

  const { setFSModule } = fromCommonjs('mol-util/data-source.js');
  setFSModule(fsModule);
  patchHeadlessCanvas3D(fromCommonjs('mol-canvas3d/canvas3d.js'));

  return {
    externalModules: { gl, pngjs, 'jpeg-js': jpegjs },
    ...fromCommonjs('mol-plugin/headless-plugin-context.js'),
    ...fromCommonjs('mol-plugin/spec.js'),
    ...fromCommonjs('mol-plugin-state/actions/file.js'),
    ...fromCommonjs('mol-plugin-state/transforms/data.js'),
    ...fromCommonjs('mol-plugin-state/transforms/model.js'),
    ...fromCommonjs('mol-plugin-state/transforms/representation.js'),
    ...fromCommonjs('extensions/geo-export/obj-exporter.js'),
    ...fromCommonjs('extensions/geo-export/stl-exporter.js'),
    structureVisualCommon: fromCommonjs('mol-repr/structure/visual/util/common.js'),
    ...fromCommonjs('mol-math/geometry.js'),
    ...fromCommonjs('mol-task/index.js'),
    ...fromCommonjs('mol-util/assets.js'),
    ...fromCommonjs('mol-util/nodejs-shims.js'),
  };
}

function patchHeadlessCanvas3D(canvas3dModule) {
  const canvas3d = canvas3dModule?.Canvas3D;
  if (!canvas3d || canvas3d.__molfigHeadlessCanvasPatch) return;
  const create = canvas3d.create;
  canvas3d.create = (ctx, props, attribs) => {
    if (!ctx.canvas) {
      const canvas = document.createElement('canvas');
      canvas.width = 800;
      canvas.height = 800;
      ctx = { ...ctx, canvas };
    }
    return create(ctx, props, attribs);
  };
  canvas3d.__molfigHeadlessCanvasPatch = true;
}

function installHeadlessDomShim() {
  if (globalThis.document && globalThis.window) return;

  class HeadlessElement {
    constructor(tagName = 'div') {
      this.tagName = String(tagName).toUpperCase();
      this.children = [];
      this.style = {};
      this.parentElement = null;
      this.name = '';
      this.width = 0;
      this.height = 0;
    }

    appendChild(child) {
      if (child && typeof child === 'object') {
        child.parentElement = this;
        this.children.push(child);
      }
      return child;
    }

    removeChild(child) {
      const index = this.children.indexOf(child);
      if (index >= 0) this.children.splice(index, 1);
      if (child && typeof child === 'object') child.parentElement = null;
      return child;
    }

    setAttribute(name, value) {
      this[name] = String(value);
    }

    getAttribute(name) {
      return this[name] ?? null;
    }

    addEventListener() {}
    removeEventListener() {}
    requestFullscreen() { return Promise.resolve(); }
    getBoundingClientRect() {
      return { left: 0, top: 0, right: this.width, bottom: this.height, width: this.width, height: this.height };
    }
  }

  class HeadlessCanvasElement extends HeadlessElement {
    constructor() {
      super('canvas');
      this.width = 800;
      this.height = 800;
    }

    getContext() {
      return null;
    }

    toDataURL() {
      return 'data:image/png;base64,';
    }
  }

  class HeadlessDocument {
    constructor() {
      this.body = new HeadlessElement('body');
      this.head = new HeadlessElement('head');
      this.documentElement = new HeadlessElement('html');
      this.scrollingElement = this.documentElement;
      this.fullscreenElement = null;
    }

    createElement(tagName) {
      return String(tagName).toLowerCase() === 'canvas'
        ? new HeadlessCanvasElement()
        : new HeadlessElement(tagName);
    }

    createElementNS(_namespace, tagName) {
      return this.createElement(tagName);
    }

    createEvent() {
      return { initMouseEvent() {} };
    }

    getElementById() {
      return null;
    }

    getElementsByTagName(tagName) {
      switch (String(tagName).toLowerCase()) {
        case 'body': return [this.body];
        case 'head': return [this.head];
        case 'html': return [this.documentElement];
        default: return [];
      }
    }

    getElementsByClassName() {
      return [];
    }

    addEventListener() {}
    removeEventListener() {}
    exitFullscreen() { return Promise.resolve(); }
  }

  const document = globalThis.document ?? new HeadlessDocument();
  const window = globalThis.window ?? globalThis;
  window.document = document;
  window.devicePixelRatio ??= 1;
  window.innerWidth ??= 800;
  window.innerHeight ??= 800;
  window.performance ??= { now: () => Date.now() };
  window.requestAnimationFrame ??= (callback) => setTimeout(() => callback(Date.now()), 16);
  window.cancelAnimationFrame ??= (id) => clearTimeout(id);
  window.addEventListener ??= () => {};
  window.removeEventListener ??= () => {};
  window.open ??= () => null;
  window.location ??= { search: '', href: '', toString: () => '' };
  globalThis.window = window;
  globalThis.document = document;
  globalThis.navigator ??= { userAgent: 'node', vendor: '', maxTouchPoints: 0 };
  globalThis.Window ??= class Window {};
  globalThis.Document ??= HeadlessDocument;
  globalThis.HTMLElement ??= HeadlessElement;
  globalThis.HTMLCanvasElement ??= HeadlessCanvasElement;
  globalThis.requestAnimationFrame ??= window.requestAnimationFrame;
  globalThis.cancelAnimationFrame ??= window.cancelAnimationFrame;
}

function installMolstarVersionShim(packageJsonPath) {
  if (typeof globalThis.__MOLSTAR_PLUGIN_VERSION__ !== 'undefined') return;
  try {
    const packageJson = JSON.parse(readFileSync(packageJsonPath, 'utf8'));
    if (packageJson.version) {
      globalThis.__MOLSTAR_PLUGIN_VERSION__ = String(packageJson.version);
    }
  } catch {
    // Keep Mol*'s own '(development)' fallback if package metadata is unavailable.
  }
}

function ensureMolstarCommonjs(args, inspection) {
  if (inspection.commonjsReady) return inspection;

  const blockers = sourceTreeBuildBlockers(args, inspection);
  if (blockers.length) {
    throw new Error(conversionBlockedMessage(args, inspection));
  }

  console.log('');
  console.log(`Mol* CommonJS prebuilt tree is missing; building from local source tree in ${args.molstarDir}`);
  const result = spawnSync('npm', ['--offline', 'run', 'build:lib'], {
    cwd: inspection.molstarDir,
    encoding: 'utf8',
    maxBuffer: 20 * 1024 * 1024,
  });

  if (result.status !== 0) {
    throw new Error(`Mol* source-tree build failed while running 'npm --offline run build:lib' in ${args.molstarDir}.\n${formatProcessOutput(result)}`);
  }

  const nextInspection = inspectMolstar(args);
  if (!nextInspection.commonjsReady) {
    throw new Error(`Mol* source-tree build completed but required CommonJS files are still missing.\n${conversionBlockedMessage(args, nextInspection)}`);
  }
  return nextInspection;
}

function conversionBlockedMessage(args, inspection) {
  const lines = ['Mol* reference exporter cannot run from the current local checkout state.'];
  if (inspection.commonjsReady) {
    lines.push(`Prebuilt CommonJS tree is present at ${displayPath(inspection.commonjsDir)}.`);
  } else {
    lines.push(`Prebuilt CommonJS tree is missing at ${args.molstarDir}/lib/commonjs.`);
  }
  lines.push('Runtime toolchain:');
  lines.push(`- node: ${inspection.toolchain.node.version} (${inspection.toolchain.node.executable})`);
  lines.push(`- npm: ${formatProbe(inspection.toolchain.npm)}`);
  for (const python of inspection.toolchain.python) {
    lines.push(`- ${python.command}: ${formatProbe(python.version)}; distutils=${formatProbe(python.distutils)}`);
  }
  lines.push('Blockers:');
  const blockers = sourceTreeBuildBlockers(args, inspection);
  if (blockers.length === 0) {
    blockers.push('required CommonJS files are missing even though source-tree prerequisites are present');
  }
  for (const blocker of blockers) {
    lines.push(`- ${blocker}`);
  }
  lines.push('No dependency install was attempted. Restore local/cache dependencies first, then rerun conversion.');
  return lines.join('\n');
}

function formatProcessOutput(result) {
  const stdout = result.stdout ? result.stdout.trim() : '';
  const stderr = result.stderr ? result.stderr.trim() : '';
  const lines = [];
  if (result.error) lines.push(`error: ${result.error.message}`);
  if (stdout) lines.push(`stdout:\n${tailLines(stdout, 40)}`);
  if (stderr) lines.push(`stderr:\n${tailLines(stderr, 40)}`);
  if (lines.length === 0) lines.push(`exit status ${result.status}`);
  return lines.join('\n');
}

function tailLines(text, maxLines) {
  const lines = text.split(/\r?\n/);
  return lines.slice(Math.max(0, lines.length - maxLines)).join('\n');
}

async function runConversion(args, plan, molstar) {
  const outDir = resolveRepoPath(args.outDir);
  mkdirSync(outDir, { recursive: true });

  for (const item of plan) {
    const outputPaths = plannedOutputPaths(args, item);
    console.log(`Converting ${item.fixture}`);
    const plugin = new molstar.HeadlessPluginContext(
      molstar.externalModules,
      molstar.DefaultPluginSpec(),
      { width: 800, height: 800 },
    );

    try {
      await plugin.init();
      configureSceneRuntime(plugin, molstar, item.options);
      await buildScene(plugin, molstar, item);
      const renderObjects = plugin.canvas3d?.getRenderObjects() ?? [];
      if (renderObjects.length === 0) throw new Error(`${item.fixture}: Mol* produced no render objects`);
      if (args.renderObjectReport) printRenderObjectReport(item, renderObjects);

      const sphere = plugin.canvas3d.boundingSphereVisible;
      const box = molstar.Box3D.fromSphere3D(molstar.Box3D(), sphere);
      const stats = {};

      if (item.formats.includes('obj')) {
        const exporter = new molstar.ObjExporter(item.objExportBasename, box);
        configureExporter(exporter, item.options);
        const data = await exportWith(plugin, molstar, renderObjects, exporter);
        mkdirSync(path.dirname(outputPaths.objPath), { recursive: true });
        mkdirSync(path.dirname(outputPaths.mtlPath), { recursive: true });
        writeFileSync(outputPaths.objPath, data.obj);
        writeFileSync(outputPaths.mtlPath, data.mtl);
        stats.obj = objExportStats(data.obj);
      }

      if (item.formats.includes('stl')) {
        const exporter = new molstar.StlExporter(box);
        configureExporter(exporter, item.options);
        const data = await exportWith(plugin, molstar, renderObjects, exporter);
        mkdirSync(path.dirname(outputPaths.stlPath), { recursive: true });
        writeFileSync(outputPaths.stlPath, data.stl);
        stats.stl = stlExportStats(data.stl);
      }

      if (item.formats.includes('json')) {
        if (!stats.obj) {
          if (!item.objReference) throw new Error(`${item.contractPath}: json summary requires obj export or obj_reference=...`);
          stats.obj = objExportStats(readFileSync(resolveRepoPath(item.objReference), 'utf8'));
        }
        if (!stats.stl) {
          if (!item.stlReference) throw new Error(`${item.contractPath}: json summary requires stl export or stl_reference=...`);
          stats.stl = stlExportStats(readFileSync(resolveRepoPath(item.stlReference)));
        }
        const summary = referenceSummaryJson(item, stats);
        mkdirSync(path.dirname(outputPaths.jsonPath), { recursive: true });
        writeFileSync(outputPaths.jsonPath, `${summary}\n`);
      }
    } finally {
      await plugin.clear?.();
      plugin.dispose();
    }
  }

  console.log(`Mol* reference outputs written to ${path.relative(rootDir, outDir)}`);
}

function configureExporter(exporter, options) {
  const primitivesQuality = options['export-primitives-quality'] ?? options.export_primitives_quality;
  if (primitivesQuality === undefined) return;
  if (!['auto', 'high', 'medium', 'low'].includes(String(primitivesQuality))) {
    throw new Error(`export-primitives-quality: unsupported value '${primitivesQuality}'. Expected auto, high, medium, or low.`);
  }
  exporter.options.primitivesQuality = String(primitivesQuality);
}

function configureSceneRuntime(plugin, molstar, options) {
  if (optionBoolean(options, 'force-cylinder-impostors', 'force_cylinder_impostors')) {
    assertCylinderImpostorRuntimeSupport(plugin);
    forceCylinderImpostorSupport(molstar);
  } else {
    restoreCylinderImpostorSupport(molstar);
  }
}

function assertCylinderImpostorRuntimeSupport(plugin) {
  if (plugin.canvas3d?.webgl?.extensions?.fragDepth) return;
  throw new Error('force-cylinder-impostors requires a headless WebGL context with GL_EXT_frag_depth; the current runtime would fail while creating Mol* cylinder shaders.');
}

function forceCylinderImpostorSupport(molstar) {
  const common = molstar.structureVisualCommon;
  if (!common || typeof common.checkCylinderImpostorSupport !== 'function') {
    throw new Error('Mol* structure visual common module is not available for cylinder-impostor patching');
  }
  if (common.__molfigForceCylinderImpostors) return;
  common.__molfigOriginalCheckCylinderImpostorSupport = common.checkCylinderImpostorSupport;
  common.checkCylinderImpostorSupport = () => true;
  common.__molfigForceCylinderImpostors = true;
}

function restoreCylinderImpostorSupport(molstar) {
  const common = molstar.structureVisualCommon;
  if (!common?.__molfigForceCylinderImpostors) return;
  common.checkCylinderImpostorSupport = common.__molfigOriginalCheckCylinderImpostorSupport;
  delete common.__molfigOriginalCheckCylinderImpostorSupport;
  delete common.__molfigForceCylinderImpostors;
}

function plannedOutputPaths(args, item) {
  const outDir = resolveRepoPath(args.outDir);
  const useContractReferences = path.resolve(outDir) === path.resolve(resolveRepoPath(defaultOutDir));
  const objPath = useContractReferences && item.objReference
    ? resolveRepoPath(item.objReference)
    : path.join(outDir, `${item.stem}.obj`);
  const mtlPath = useContractReferences && item.objReference
    ? path.join(path.dirname(objPath), item.objMtllib)
    : path.join(outDir, item.objMtllib);
  const stlPath = useContractReferences && item.stlReference
    ? resolveRepoPath(item.stlReference)
    : path.join(outDir, `${item.stem}.stl`);
  const jsonPath = useContractReferences && item.jsonReference
    ? resolveRepoPath(item.jsonReference)
    : path.join(outDir, `${item.stem}.summary.json`);

  return { objPath, mtlPath, stlPath, jsonPath };
}

async function buildScene(plugin, molstar, item) {
  const sceneSource = normalizeSceneSource(item.options['scene-source'] ?? item.options.scene_source ?? 'manual');
  if (sceneSource === 'data-format') {
    await buildSceneFromDataFormat(plugin, molstar, item);
    return;
  }
  if (sceneSource === 'open-files') {
    await buildSceneFromOpenFilesAction(plugin, molstar, item);
    return;
  }
  if (sceneSource !== 'manual') {
    throw new Error(`${item.contractPath}: unsupported scene-source '${sceneSource}'. Expected manual, data-format, or open-files.`);
  }

  const update = plugin.build();
  let state = update.toRoot().apply(molstar.RawData, {
    data: readFixtureData(item),
    label: path.basename(item.fixture),
  });

  if (item.inputFormat === 'pdb') {
    state = state.apply(molstar.TrajectoryFromPDB);
  } else {
    state = state.apply(molstar.ParseCif).apply(molstar.TrajectoryFromMmCif, cifParams(item.options));
  }

  const structure = state
    .apply(molstar.ModelFromTrajectory, { modelIndex: Number(item.options['model-index'] ?? 0) })
    .apply(molstar.StructureFromModel, { type: structureType(item.options.assembly) });

  const preset = item.options['representation-preset'] ?? item.options.representation_preset;
  if (preset !== undefined) {
    const structureRef = structure.ref;
    await update.commit();
    await plugin.builders.structure.representation.applyPreset(
      structureRef,
      String(preset),
      representationPresetParams(item.options),
    );
    plugin.canvas3d?.commit(true);
    return;
  }

  applyRepresentations(structure, molstar, item.options);
  await update.commit();
  plugin.canvas3d?.commit(true);
}

async function buildSceneFromDataFormat(plugin, molstar, item) {
  ensureNoCifBlockSelection(item, 'data-format');

  const file = new molstar.File_(
    [readFixtureData(item)],
    path.basename(item.fixture),
  );
  const asset = molstar.Asset.File(file);
  const { data, fileInfo } = await plugin.builders.data.readFile({
    file: asset,
    isBinary: item.inputFormat === 'bcif',
  });

  const providerId = dataFormatProviderId(item);
  const provider = providerId === 'auto'
    ? plugin.dataFormats.auto(fileInfo, data.cell?.obj)
    : plugin.dataFormats.get(providerId);
  if (!provider) {
    throw new Error(`${item.contractPath}: Mol* data-format provider '${providerId}' was not found for ${item.fixture}`);
  }

  const parsed = await provider.parse(plugin, data);
  if (!parsed?.trajectory) {
    throw new Error(`${item.contractPath}: Mol* data-format provider '${providerId}' did not return a trajectory`);
  }

  await plugin.builders.structure.hierarchy.applyPreset(
    parsed.trajectory,
    String(item.options['hierarchy-preset'] ?? item.options.hierarchy_preset ?? 'default'),
    hierarchyPresetParams(item.options),
  );
  plugin.canvas3d?.commit(true);
}

async function buildSceneFromOpenFilesAction(plugin, molstar, item) {
  ensureNoCifBlockSelection(item, 'open-files');
  if (!molstar.OpenFiles) {
    throw new Error(`${item.contractPath}: Mol* OpenFiles StateAction is not available in the loaded CommonJS tree`);
  }

  const file = new molstar.File_(
    [readFixtureData(item)],
    path.basename(item.fixture),
  );
  const asset = molstar.Asset.File(file);
  await plugin.runTask(plugin.state.data.applyAction(molstar.OpenFiles, {
    files: [asset],
    format: openFilesFormatParam(item),
    visuals: true,
  }));
  plugin.canvas3d?.commit(true);
}

function normalizeSceneSource(value) {
  const normalized = String(value).trim().toLowerCase().replace(/_/g, '-');
  if (normalized === 'open-files-action' || normalized === 'web-app-open-files') return 'open-files';
  return normalized;
}

function ensureNoCifBlockSelection(item, sceneSource) {
  if ((item.options['block-header'] ?? item.options.block_header ?? item.options['block-index'] ?? item.options.block_index) !== undefined) {
    throw new Error(`${item.contractPath}: scene-source=${sceneSource} does not support CIF block selection yet; use scene-source=manual.`);
  }
}

function readFixtureData(item) {
  if (item.inputFormat === 'bcif') return readFileSync(item.absFixturePath);
  return readFileSync(item.absFixturePath, 'utf8');
}

function openFilesFormatParam(item) {
  const provider = item.options['open-files-format']
    ?? item.options.open_files_format
    ?? item.options['data-format']
    ?? item.options.data_format;
  if (provider === undefined || String(provider) === 'auto') {
    return { name: 'auto', params: {} };
  }
  return { name: 'specific', params: String(provider) };
}

function dataFormatProviderId(item) {
  const provider = item.options['data-format'] ?? item.options.data_format;
  if (provider !== undefined) return String(provider);
  switch (item.inputFormat) {
    case 'pdb': return 'pdb';
    case 'cif':
    case 'bcif': return 'mmcif';
    default: return 'auto';
  }
}

function cifParams(options) {
  const params = {};
  if (options['block-header']) params.blockHeader = String(options['block-header']);
  if (options['block-index'] !== undefined) params.blockIndex = Number(options['block-index']);
  return params;
}

function structureType(assembly) {
  if (assembly === 'asymmetric-unit' || assembly === 'none') return { name: 'model', params: {} };
  const id = assembly === undefined || assembly === null || assembly === '' ? '1' : String(assembly);
  return { name: 'assembly', params: { id, dynamicBonds: false } };
}

function hierarchyPresetParams(options) {
  const params = {
    showUnitcell: false,
    structure: structureType(options.assembly),
  };
  const modelIndex = options['model-index'] ?? options.model_index;
  if (modelIndex !== undefined) params.model = { modelIndex: Number(modelIndex) };
  const representationPreset = options['representation-preset'] ?? options.representation_preset;
  if (representationPreset !== undefined) params.representationPreset = String(representationPreset);
  const representationParams = representationPresetParams(options);
  if (Object.keys(representationParams).length > 0) params.representationPresetParams = representationParams;
  return params;
}

function applyRepresentations(structure, molstar, options) {
  const representation = String(options.representation ?? 'molstar');
  if (representation === 'spacefill') {
    addRepresentation(structure, molstar, 'all', 'spacefill', 'element-symbol', options);
  } else if (representation === 'ball-and-stick') {
    addRepresentation(structure, molstar, 'all', 'ball-and-stick', 'element-symbol', options);
  } else {
    addRepresentation(structure, molstar, 'polymer', 'cartoon', 'sequence-id', options);
    if (representation === 'molstar') {
      addRepresentation(structure, molstar, 'ligand', 'ball-and-stick', 'element-symbol', options);
      addRepresentation(structure, molstar, 'branched', 'ball-and-stick', 'element-symbol', options);
      addRepresentation(structure, molstar, 'ion', 'spacefill', 'element-symbol', options);
    }
  }
}

function addRepresentation(structure, molstar, component, typeName, colorName, options) {
  structure
    .apply(molstar.StructureComponent, {
      type: { name: 'static', params: component },
      nullIfEmpty: true,
      label: component,
    })
    .apply(molstar.StructureRepresentation3D, {
      type: { name: typeName, params: representationParams(typeName, options) },
      colorTheme: { name: colorName, params: {} },
      sizeTheme: { name: 'physical', params: {} },
    });
}

function representationParams(typeName, options) {
  const params = { alpha: 1 };
  copyStringParam(params, options, 'quality', 'quality');
  copyNumberParam(params, options, 'sphere-detail', 'detail');
  copyNumberParam(params, options, 'sphere_detail', 'detail');
  copyNumberParam(params, options, 'linear-segments', 'linearSegments');
  copyNumberParam(params, options, 'linear_segments', 'linearSegments');
  copyNumberParam(params, options, 'radial-segments', 'radialSegments');
  copyNumberParam(params, options, 'radial_segments', 'radialSegments');

  if (typeName === 'cartoon') {
    copyNumberParam(params, options, 'sheet-arrow-factor', 'arrowFactor');
    copyNumberParam(params, options, 'sheet_arrow_factor', 'arrowFactor');
    copyBooleanParam(params, options, 'tubular-helices', 'tubularHelices');
    copyBooleanParam(params, options, 'tubular_helices', 'tubularHelices');
    copyBooleanParam(params, options, 'round-cap', 'roundCap');
    copyBooleanParam(params, options, 'round_cap', 'roundCap');
    copyStringParam(params, options, 'helix-profile', 'helixProfile');
    copyStringParam(params, options, 'helix_profile', 'helixProfile');
  }

  return params;
}

function representationPresetParams(options) {
  const params = {};
  copyStringParam(params, options, 'quality', 'quality');
  return params;
}

function copyNumberParam(params, options, from, to) {
  if (options[from] === undefined) return;
  const value = Number(options[from]);
  if (Number.isFinite(value)) params[to] = value;
}

function copyBooleanParam(params, options, from, to) {
  if (options[from] === undefined) return;
  params[to] = optionBoolean(options, from);
}

function copyStringParam(params, options, from, to) {
  if (options[from] === undefined) return;
  const value = String(options[from]);
  if (value) params[to] = value;
}

function optionBoolean(options, ...keys) {
  for (const key of keys) {
    if (options[key] === undefined) continue;
    const value = options[key];
    return value === true || value === 'true' || value === 1 || value === '1';
  }
  return false;
}

async function exportWith(plugin, molstar, renderObjects, exporter) {
  return plugin.runTask(molstar.Task.create('Export Mol* reference geometry', async ctx => {
    for (let i = 0; i < renderObjects.length; i++) {
      await ctx.update({ message: `Exporting object ${i + 1}/${renderObjects.length}` });
      await exporter.add(renderObjects[i], plugin.canvas3d.webgl, ctx);
    }
    return exporter.getData();
  }));
}

function printRenderObjectReport(item, renderObjects) {
  console.log(`Mol* render objects for ${item.stem}: ${renderObjects.length}`);
  let totalDrawCount = 0;
  let visibleDrawCount = 0;
  let hiddenDrawCount = 0;
  let exportableDrawCount = 0;
  for (let index = 0; index < renderObjects.length; index++) {
    const renderObject = renderObjects[index];
    const values = renderObject.values ?? {};
    const drawCount = valueCellValue(values.drawCount) ?? 0;
    const vertexCount = valueCellValue(values.uVertexCount) ?? valueCellValue(values.vertexCount) ?? 0;
    const groupCount = valueCellValue(values.uGroupCount) ?? valueCellValue(values.groupCount) ?? 0;
    const instanceCount = valueCellValue(values.instanceCount) ?? valueCellValue(values.uInstanceCount) ?? 0;
    const geometryType = valueCellValue(values.dGeometryType) ?? renderObject.type ?? '<unknown>';
    const visible = renderObject.state?.visible ?? true;
    const numericDrawCount = Number(drawCount) || 0;
    totalDrawCount += numericDrawCount;
    if (visible) {
      visibleDrawCount += numericDrawCount;
      if (numericDrawCount > 0 && Number(instanceCount) !== 0) {
        exportableDrawCount += numericDrawCount;
      }
    } else {
      hiddenDrawCount += numericDrawCount;
    }
    console.log(`  [${index}] type=${renderObject.type} geometry=${geometryType} visible=${visible} drawCount=${drawCount} vertexCount=${vertexCount} groupCount=${groupCount} instanceCount=${instanceCount}`);
  }
  console.log(`  totalDrawCount=${totalDrawCount} visibleDrawCount=${visibleDrawCount} hiddenDrawCount=${hiddenDrawCount} exportableDrawCount=${exportableDrawCount}`);
}

function valueCellValue(cell) {
  return cell?.ref?.value;
}

function referenceSummaryJson(item, stats) {
  return `{
  "schema": 1,
  "name": ${JSON.stringify(item.stem)},
  "molstar_reference_commit": ${JSON.stringify(molstarReferenceCommit)},
  "fixture": ${JSON.stringify(item.fixture)},
  "options": {
    "format": ${JSON.stringify(item.inputFormat)},
    "representation": ${JSON.stringify(String(item.options.representation ?? 'molstar'))},
    "assembly": ${JSON.stringify(String(item.options.assembly ?? '1'))},
    "sphere-detail": ${Number(item.options['sphere-detail'] ?? item.options.sphere_detail ?? 1)}
  },
  "source": {
    "obj": ${JSON.stringify(item.objReference || `${item.stem}.obj`)},
    "stl": ${JSON.stringify(item.stlReference || `${item.stem}.stl`)}
  },
  "obj": {
    "byte_len": ${stats.obj.byte_len},
    "fnv1a64": ${JSON.stringify(stats.obj.fnv1a64)},
    "vertex_count": ${stats.obj.vertex_count},
    "normal_count": ${stats.obj.normal_count},
    "face_count": ${stats.obj.face_count},
    "bounds": {
      "min": ${formatVec3(stats.obj.min)},
      "max": ${formatVec3(stats.obj.max)}
    }
  },
  "stl": {
    "byte_len": ${stats.stl.byte_len},
    "fnv1a64": ${JSON.stringify(stats.stl.fnv1a64)},
    "facet_count": ${stats.stl.facet_count},
    "header": ${JSON.stringify(stats.stl.header)},
    "bounds": {
      "min": ${formatVec3(stats.stl.min)},
      "max": ${formatVec3(stats.stl.max)}
    }
  }
}`;
}

function objExportStats(obj) {
  const bytes = Buffer.from(obj);
  const stats = {
    byte_len: bytes.length,
    fnv1a64: fnv1a64(bytes),
    mtllib: undefined,
    mtllibs: [],
    firstContentLine: undefined,
    usemtl: [],
    materials: [],
    vertex_count: 0,
    normal_count: 0,
    face_count: 0,
    vertices: [],
    faces: [],
    nonTriFaceCount: 0,
    faceVertexNormalMismatchCount: 0,
    unsupportedLines: [],
    invalidVectors: [],
    invalidFaces: [],
    faceIndexErrors: [],
    faceRefs: [],
    min: [Infinity, Infinity, Infinity],
    max: [-Infinity, -Infinity, -Infinity],
  };
  let lineNumber = 0;
  for (const line of obj.split(/\r?\n/)) {
    lineNumber++;
    if (!line.trim()) continue;
    if (stats.firstContentLine === undefined) stats.firstContentLine = line;
    if (line.startsWith('mtllib ')) {
      stats.mtllibs.push(line);
      if (stats.mtllib === undefined) stats.mtllib = line;
    } else if (line.startsWith('usemtl ')) {
      const material = line.slice('usemtl '.length).trim();
      stats.usemtl.push(material);
      if (!stats.materials.includes(material)) stats.materials.push(material);
    } else if (line.startsWith('v ')) {
      stats.vertex_count++;
      const values = line.slice(2).trim().split(/\s+/).slice(0, 3).map(Number);
      if (values.length !== 3 || values.some(value => !Number.isFinite(value))) {
        pushLineExample(stats.invalidVectors, lineNumber, line);
        continue;
      }
      stats.vertices.push(values);
      for (let axis = 0; axis < 3; axis++) {
        stats.min[axis] = Math.min(stats.min[axis], values[axis]);
        stats.max[axis] = Math.max(stats.max[axis], values[axis]);
      }
    } else if (line.startsWith('vn ')) {
      stats.normal_count++;
      const values = line.slice(3).trim().split(/\s+/).slice(0, 3).map(Number);
      if (values.length !== 3 || values.some(value => !Number.isFinite(value))) {
        pushLineExample(stats.invalidVectors, lineNumber, line);
      }
    } else if (line.startsWith('f ')) {
      stats.face_count++;
      const refs = line.slice(2).trim().split(/\s+/);
      if (refs.length !== 3) stats.nonTriFaceCount++;
      const face = [];
      for (const ref of refs) {
        const match = /^([1-9]\d*)\/\/([1-9]\d*)$/.exec(ref);
        if (!match) {
          pushLineExample(stats.invalidFaces, lineNumber, line);
          continue;
        }
        const vertexIndex = Number(match[1]);
        const normalIndex = Number(match[2]);
        if (vertexIndex !== normalIndex) stats.faceVertexNormalMismatchCount++;
        stats.faceRefs.push({ lineNumber, line, vertexIndex, normalIndex });
        face.push(vertexIndex - 1);
      }
      if (face.length === 3) stats.faces.push(face);
    } else {
      pushLineExample(stats.unsupportedLines, lineNumber, line);
    }
  }
  for (const ref of stats.faceRefs) {
    if (ref.vertexIndex > stats.vertex_count || ref.normalIndex > stats.normal_count) {
      pushLineExample(stats.faceIndexErrors, ref.lineNumber, ref.line);
    }
  }
  delete stats.faceRefs;
  return stats;
}

function mtlExportStats(mtl) {
  const bytes = Buffer.from(mtl);
  const stats = {
    byte_len: bytes.length,
    fnv1a64: fnv1a64(bytes),
    materials: [],
  };
  for (const line of mtl.split(/\r?\n/)) {
    if (!line.startsWith('newmtl ')) continue;
    const material = line.slice('newmtl '.length).trim();
    if (stats.materials.includes(material)) {
      throw new Error(`MTL newmtl material appears more than once: ${material}`);
    }
    stats.materials.push(material);
  }
  return stats;
}

function stlExportStats(stl) {
  const bytes = Buffer.from(stl);
  if (bytes.length < 84) {
    throw new Error(`STL file is too short: expected at least 84 bytes, got ${bytes.length}`);
  }
  const facetCount = bytes.readUInt32LE(80);
  const expectedByteLength = 84 + facetCount * 50;
  if (bytes.length !== expectedByteLength) {
    throw new Error(`STL byte length does not match facet count: got ${bytes.length}, expected ${expectedByteLength}`);
  }
  const stats = {
    bytes,
    byte_len: bytes.length,
    fnv1a64: fnv1a64(bytes),
    facet_count: facetCount,
    header: bytes.subarray(0, 80).toString('utf8').replace(/\0+$/g, ''),
    nonzero_attribute_count: 0,
    invalidFloatCount: 0,
    min: [Infinity, Infinity, Infinity],
    max: [-Infinity, -Infinity, -Infinity],
  };
  for (let facet = 0; facet < facetCount; facet++) {
    const base = 84 + facet * 50 + 12;
    const attributeByteCount = bytes.readUInt16LE(84 + facet * 50 + 48);
    if (attributeByteCount !== 0) stats.nonzero_attribute_count++;
    for (let vertex = 0; vertex < 3; vertex++) {
      for (let axis = 0; axis < 3; axis++) {
        const value = bytes.readFloatLE(base + (vertex * 3 + axis) * 4);
        if (!Number.isFinite(value)) {
          stats.invalidFloatCount++;
          continue;
        }
        stats.min[axis] = Math.min(stats.min[axis], value);
        stats.max[axis] = Math.max(stats.max[axis], value);
      }
    }
  }
  return stats;
}

function objStlSparseSlotDrift(objVertices, objFaces, stl) {
  const drift = {
    total_components: 0,
    rounded_mismatch_count: 0,
    first_rounded_mismatch: undefined,
    max_delta: 0,
    max_delta_component: undefined,
  };
  for (let faceIndex = 0; faceIndex < objFaces.length; faceIndex++) {
    const face = objFaces[faceIndex];
    const stlFacetIndex = faceIndex * 3;
    for (let vertexSlot = 0; vertexSlot < 3; vertexSlot++) {
      const objVertex = objVertices[face[vertexSlot]];
      const stlVertex = stlFacetVertex(stl, stlFacetIndex, vertexSlot);
      for (let axis = 0; axis < 3; axis++) {
        drift.total_components++;
        const objValue = objVertex[axis];
        const stlValue = stlVertex[axis];
        const roundedStlValue = molstarObjRoundedCoordinate(stlValue);
        const component = {
          face_index: faceIndex,
          stl_facet_index: stlFacetIndex,
          vertex_slot: vertexSlot,
          axis,
          obj_value: objValue,
          stl_value: stlValue,
          rounded_stl_value: roundedStlValue,
        };
        const delta = Math.abs(objValue - stlValue);
        if (delta > drift.max_delta) {
          drift.max_delta = delta;
          drift.max_delta_component = component;
        }
        if (Math.abs(roundedStlValue - objValue) > 0.00001) {
          drift.rounded_mismatch_count++;
          if (!drift.first_rounded_mismatch) drift.first_rounded_mismatch = component;
        }
      }
    }
  }
  return drift;
}

function stlFacetVertex(stl, facetIndex, vertexIndex) {
  const base = 84 + facetIndex * 50 + 12 + vertexIndex * 12;
  return [
    stl.readFloatLE(base),
    stl.readFloatLE(base + 4),
    stl.readFloatLE(base + 8),
  ];
}

function molstarObjRoundedCoordinate(value) {
  const rounded = Math.floor(value * 1000 + 0.5) / 1000;
  return rounded === 0 ? 0 : rounded;
}

function pushLineExample(examples, lineNumber, line) {
  if (examples.length < 5) examples.push({ lineNumber, line });
}

function formatLineExamples(examples) {
  return examples.map(example => `${example.lineNumber}:${JSON.stringify(example.line)}`).join(', ');
}

function fnv1a64(bytes) {
  let hash = 0xcbf29ce484222325n;
  for (const byte of bytes) {
    hash ^= BigInt(byte);
    hash = (hash * 0x100000001b3n) & 0xffffffffffffffffn;
  }
  return hash.toString(16).padStart(16, '0');
}

function formatVec3(values) {
  return `[${values.map(value => formatFloat(value)).join(', ')}]`;
}

function formatFloat(value) {
  const rounded = Math.round(Number(value) * 10000) / 10000;
  if (Math.abs(rounded) < 0.00005) return '0';
  if (Math.abs(Math.round(rounded) - rounded) < 0.00005) return String(Math.round(rounded));
  return rounded.toFixed(4).replace(/0+$/g, '').replace(/\.$/, '');
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  if (args.help) {
    printHelp();
    return;
  }

  const plan = loadPlan(args);
  const artifactManifest = loadArtifactManifest(args);
  const inspection = inspectMolstar(args);
  if (args.dryRun) {
    printPlan(args, plan, inspection, artifactManifest);
    validateExistingReferences(args, plan, artifactManifest);
    return;
  }

  printPlan(args, plan, inspection, artifactManifest);
  const molstar = loadMolstar(args, inspection);
  await runConversion(args, plan, molstar);
  if (path.resolve(resolveRepoPath(args.outDir)) === path.resolve(resolveRepoPath(defaultOutDir))) {
    validateExistingReferences(args, plan, artifactManifest);
  }
}

main().catch(error => {
  console.error(`error: ${error.stack || error.message}`);
  process.exitCode = 1;
});
