#!/usr/bin/env node

import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { inflateSync } from 'node:zlib';

const scriptDir = dirname(fileURLToPath(import.meta.url));
const pluginRoot = resolve(scriptDir, '..');
const repoRoot = resolve(pluginRoot, '..');

function fail(message) {
  throw new Error(message);
}

function run(command, args, cwd) {
  const result = spawnSync(command, args, {
    cwd,
    encoding: 'utf8',
    maxBuffer: 64 * 1024 * 1024,
  });
  if (result.status !== 0) {
    process.stderr.write(result.stdout ?? '');
    process.stderr.write(result.stderr ?? '');
    fail(`${command} exited with status ${result.status}`);
  }
}

function generateReferences(outputDir) {
  // The first structure loaded into a fresh Mol* canvas can be captured while
  // WebGL postprocessing is still warming up. Render the XYZ fixture twice and
  // use the second image; both entries share a stem, so the settled pass
  // replaces the warm-up capture deterministically.
  mkdirSync(outputDir, { recursive: true });
  const manifest = join(outputDir, 'illustrative-validation-fixtures.txt');
  writeFileSync(manifest, [
    'contract=tests/expected/molstar-reference/molstar-illustrative-benzene.reference.contract formats=report',
    'contract=tests/expected/molstar-reference/molstar-illustrative-benzene.reference.contract formats=report',
    'contract=tests/expected/molstar-reference/molstar-illustrative-pdb.reference.contract formats=report',
    '',
  ].join('\n'));
  run(process.execPath, [
    join(scriptDir, 'molstar-browser-reference-convert.mjs'),
    '--manifest',
    manifest,
    '--out-dir',
    outputDir,
    '--formats',
    'report',
    '--capture-images',
  ], pluginRoot);

  run('typst', [
    'compile',
    '--root',
    repoRoot,
    '--ppi',
    '72',
    join(pluginRoot, 'tests/api/illustrative-rendering-reference.typ'),
    join(outputDir, 'molfig-{p}.png'),
  ], repoRoot);

  return {
    xyz: {
      molstar: join(outputDir, 'molstar-illustrative-benzene.png'),
      molfig: join(outputDir, 'molfig-1.png'),
      report: join(outputDir, 'molstar-illustrative-benzene.browser-report.json'),
    },
    pdb: {
      molstar: join(outputDir, 'molstar-illustrative-pdb.png'),
      molfig: join(outputDir, 'molfig-2.png'),
      report: join(outputDir, 'molstar-illustrative-pdb.browser-report.json'),
    },
  };
}

function assertExact(label, actual, expected) {
  if (actual !== expected) fail(`${label}: expected ${JSON.stringify(expected)}, got ${JSON.stringify(actual)}`);
}

function validateMolstarRenderingReport(path) {
  const report = JSON.parse(readFileSync(path, 'utf8'));
  const rendering = report.rendering;
  if (!rendering) fail(`${path}: browser report is missing realized rendering state`);
  assertExact(`${path}: ignoreLight`, rendering.ignoreLight, true);
  assertExact(`${path}: camera.mode`, rendering.camera?.mode, 'perspective');
  assertExact(`${path}: camera.fov`, rendering.camera?.fov, Math.PI / 4);
  assertExact(`${path}: viewport.width`, rendering.viewport?.width, 1024);
  assertExact(`${path}: viewport.height`, rendering.viewport?.height, 937);

  const outline = rendering.postprocessing?.outline;
  assertExact(`${path}: outline.name`, outline?.name, 'on');
  assertExact(`${path}: outline.scale`, outline?.params?.scale, 1);
  assertExact(`${path}: outline.color`, outline?.params?.color, 0);
  assertExact(`${path}: outline.threshold`, outline?.params?.threshold, 0.33);
  assertExact(`${path}: outline.includeTransparent`, outline?.params?.includeTransparent, true);

  const occlusion = rendering.postprocessing?.occlusion;
  assertExact(`${path}: occlusion.name`, occlusion?.name, 'on');
  assertExact(`${path}: occlusion.multiScale`, occlusion?.params?.multiScale?.name, 'off');
  assertExact(`${path}: occlusion.radius`, occlusion?.params?.radius, 5);
  assertExact(`${path}: occlusion.bias`, occlusion?.params?.bias, 0.8);
  assertExact(`${path}: occlusion.blurKernelSize`, occlusion?.params?.blurKernelSize, 15);
  assertExact(`${path}: occlusion.blurDepthBias`, occlusion?.params?.blurDepthBias, 0.5);
  assertExact(`${path}: occlusion.samples`, occlusion?.params?.samples, 32);
  assertExact(`${path}: occlusion.resolutionScale`, occlusion?.params?.resolutionScale, 1);
  assertExact(`${path}: occlusion.color`, occlusion?.params?.color, 0);
  assertExact(`${path}: occlusion.transparentThreshold`, occlusion?.params?.transparentThreshold, 0.4);
  assertExact(`${path}: shadow.name`, rendering.postprocessing?.shadow?.name, 'off');
}

function paeth(a, b, c) {
  const p = a + b - c;
  const pa = Math.abs(p - a);
  const pb = Math.abs(p - b);
  const pc = Math.abs(p - c);
  return pa <= pb && pa <= pc ? a : pb <= pc ? b : c;
}

function decodePng(path) {
  const data = readFileSync(path);
  const signature = Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]);
  if (!data.subarray(0, 8).equals(signature)) fail(`${path}: invalid PNG signature`);

  let width;
  let height;
  let bitDepth;
  let colorType;
  let interlace;
  const idat = [];
  for (let offset = 8; offset < data.length;) {
    const length = data.readUInt32BE(offset);
    const type = data.toString('ascii', offset + 4, offset + 8);
    const body = data.subarray(offset + 8, offset + 8 + length);
    if (type === 'IHDR') {
      width = body.readUInt32BE(0);
      height = body.readUInt32BE(4);
      bitDepth = body[8];
      colorType = body[9];
      interlace = body[12];
    } else if (type === 'IDAT') {
      idat.push(body);
    } else if (type === 'IEND') {
      break;
    }
    offset += 12 + length;
  }
  if (bitDepth !== 8 || ![2, 6].includes(colorType) || interlace !== 0) {
    fail(`${path}: expected non-interlaced 8-bit RGB/RGBA PNG, got depth=${bitDepth}, type=${colorType}, interlace=${interlace}`);
  }

  const channels = colorType === 2 ? 3 : 4;
  const stride = width * channels;
  const raw = inflateSync(Buffer.concat(idat));
  if (raw.length !== height * (stride + 1)) fail(`${path}: unexpected inflated byte length`);
  const pixels = Buffer.alloc(width * height * channels);
  for (let y = 0, source = 0, target = 0; y < height; y += 1) {
    const filter = raw[source++];
    for (let x = 0; x < stride; x += 1, source += 1, target += 1) {
      const left = x >= channels ? pixels[target - channels] : 0;
      const up = y > 0 ? pixels[target - stride] : 0;
      const upperLeft = y > 0 && x >= channels ? pixels[target - stride - channels] : 0;
      const value = raw[source];
      if (filter === 0) pixels[target] = value;
      else if (filter === 1) pixels[target] = (value + left) & 0xff;
      else if (filter === 2) pixels[target] = (value + up) & 0xff;
      else if (filter === 3) pixels[target] = (value + Math.floor((left + up) / 2)) & 0xff;
      else if (filter === 4) pixels[target] = (value + paeth(left, up, upperLeft)) & 0xff;
      else fail(`${path}: unsupported PNG filter ${filter}`);
    }
  }
  return { width, height, channels, pixels };
}

function percentile(sorted, fraction) {
  return sorted[Math.min(sorted.length - 1, Math.floor(sorted.length * fraction))];
}

function isBackgroundPixel(r, g, b) {
  return r >= 250 && g >= 249 && b >= 248;
}

function pixelRgb(image, x, y) {
  if (x < 0 || x >= image.width || y < 0 || y >= image.height) return [252, 251, 250];
  const offset = (y * image.width + x) * image.channels;
  return [image.pixels[offset], image.pixels[offset + 1], image.pixels[offset + 2]];
}

function spatialComparison(referencePath, candidatePath) {
  const reference = decodePng(referencePath);
  const candidate = decodePng(candidatePath);
  if (reference.width !== candidate.width || reference.height !== candidate.height) {
    fail(`spatial comparison dimensions differ (${reference.width}x${reference.height} vs ${candidate.width}x${candidate.height})`);
  }

  let best;
  for (let shiftY = -2; shiftY <= 2; shiftY += 1) {
    for (let shiftX = -2; shiftX <= 2; shiftX += 1) {
      let union = 0;
      let intersection = 0;
      let disagreement = 0;
      let rgbDelta = 0;
      let luminanceDelta = 0;
      for (let y = 0; y < reference.height * 0.85; y += 1) {
        for (let x = Math.floor(reference.width * 0.2); x < reference.width * 0.8; x += 1) {
          const referenceRgb = pixelRgb(reference, x, y);
          const candidateRgb = pixelRgb(candidate, x + shiftX, y + shiftY);
          const referenceForeground = !isBackgroundPixel(...referenceRgb);
          const candidateForeground = !isBackgroundPixel(...candidateRgb);
          if (!referenceForeground && !candidateForeground) continue;
          union += 1;
          if (referenceForeground && candidateForeground) intersection += 1;
          else disagreement += 1;
          rgbDelta += Math.abs(referenceRgb[0] - candidateRgb[0]);
          rgbDelta += Math.abs(referenceRgb[1] - candidateRgb[1]);
          rgbDelta += Math.abs(referenceRgb[2] - candidateRgb[2]);
          const referenceLuminance = referenceRgb[0] * 0.2126 + referenceRgb[1] * 0.7152 + referenceRgb[2] * 0.0722;
          const candidateLuminance = candidateRgb[0] * 0.2126 + candidateRgb[1] * 0.7152 + candidateRgb[2] * 0.0722;
          luminanceDelta += Math.abs(referenceLuminance - candidateLuminance);
        }
      }
      const comparison = {
        shiftX,
        shiftY,
        foregroundJaccard: intersection / union,
        foregroundDisagreement: disagreement / union,
        rgbMae: rgbDelta / (union * 3 * 255),
        luminanceMae: luminanceDelta / (union * 255),
      };
      const score = comparison.foregroundDisagreement + comparison.rgbMae + comparison.luminanceMae;
      if (!best || score < best.score) best = { ...comparison, score };
    }
  }
  return best;
}

function metrics(path) {
  const image = decodePng(path);
  const colors = new Map();
  const darkColors = new Map();
  const greenChannels = [];
  const moleculeLuminance = [];
  let greenPixels = 0;
  let darkPixels = 0;
  let blackPixels = 0;
  let backgroundPixels = 0;
  let moleculePixels = 0;
  let moleculeMinX = image.width;
  let moleculeMinY = image.height;
  let moleculeMaxX = -1;
  let moleculeMaxY = -1;
  for (let i = 0, pixel = 0; i < image.pixels.length; i += image.channels, pixel += 1) {
    const r = image.pixels[i];
    const g = image.pixels[i + 1];
    const b = image.pixels[i + 2];
    const a = image.channels === 4 ? image.pixels[i + 3] : 255;
    if (a === 0) continue;
    const background = isBackgroundPixel(r, g, b);
    if (background) backgroundPixels += 1;
    const x = pixel % image.width;
    const y = Math.floor(pixel / image.width);
    // Mol* draws its orientation axes at bottom-left. Restrict the molecule
    // mask to the centered viewport region so that UI decoration cannot make
    // a structurally wrong rendering pass the geometry checks.
    const molecule = !background && x >= image.width * 0.2 && x < image.width * 0.8 && y < image.height * 0.85;
    if (molecule) {
      moleculePixels += 1;
      moleculeMinX = Math.min(moleculeMinX, x);
      moleculeMinY = Math.min(moleculeMinY, y);
      moleculeMaxX = Math.max(moleculeMaxX, x);
      moleculeMaxY = Math.max(moleculeMaxY, y);
      moleculeLuminance.push(r * 0.2126 + g * 0.7152 + b * 0.0722);
      const green = g >= 48 && g > r * 2.2 && g > b * 1.18;
      if (green) {
        greenPixels += 1;
        greenChannels.push(g);
        const key = `${r},${g},${b}`;
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
      if (r <= 35 && g <= 100 && b <= 80) {
        darkPixels += 1;
        const key = `${r},${g},${b}`;
        darkColors.set(key, (darkColors.get(key) ?? 0) + 1);
      }
      if (r <= 4 && g <= 4 && b <= 4) blackPixels += 1;
    }
  }
  greenChannels.sort((a, b) => a - b);
  moleculeLuminance.sort((a, b) => a - b);
  if (greenPixels === 0) fail(`${path}: no molecular green pixels detected`);
  const [dominantGreen, dominantGreenCount] = [...colors.entries()].sort((a, b) => b[1] - a[1])[0];
  const [dominantDark, dominantDarkCount] = [...darkColors.entries()].sort((a, b) => b[1] - a[1])[0];
  return {
    path,
    width: image.width,
    height: image.height,
    dominantGreen,
    dominantGreenCount,
    greenPixels,
    greenP10: percentile(greenChannels, 0.10),
    greenP25: percentile(greenChannels, 0.25),
    greenP50: percentile(greenChannels, 0.50),
    greenP75: percentile(greenChannels, 0.75),
    greenP90: percentile(greenChannels, 0.90),
    greenMean: greenChannels.reduce((sum, value) => sum + value, 0) / greenChannels.length,
    luminanceP10: percentile(moleculeLuminance, 0.10),
    luminanceP50: percentile(moleculeLuminance, 0.50),
    luminanceP90: percentile(moleculeLuminance, 0.90),
    darkPixels,
    dominantDark,
    dominantDarkCount,
    blackPixels,
    darkToGreen: darkPixels / greenPixels,
    darkFraction: darkPixels / (image.width * image.height),
    darkMoleculeFraction: darkPixels / moleculePixels,
    blackMoleculeFraction: blackPixels / moleculePixels,
    backgroundFraction: backgroundPixels / (image.width * image.height),
    moleculePixels,
    moleculeBounds: {
      x: moleculeMinX,
      y: moleculeMinY,
      width: moleculeMaxX - moleculeMinX + 1,
      height: moleculeMaxY - moleculeMinY + 1,
    },
  };
}

function assertPair(label, reference, candidate) {
  if (reference.width !== candidate.width || reference.height !== candidate.height) {
    fail(`${label}: dimensions differ (${reference.width}x${reference.height} vs ${candidate.width}x${candidate.height})`);
  }
  const referenceDominant = reference.dominantGreen.split(',').map(Number);
  const candidateDominant = candidate.dominantGreen.split(',').map(Number);
  const dominantDelta = Math.max(...referenceDominant.map((value, index) => Math.abs(value - candidateDominant[index])));
  if (dominantDelta !== 0) {
    fail(`${label}: dominant material color differs (${reference.dominantGreen} vs ${candidate.dominantGreen})`);
  }
  for (const percentile of ['greenP10', 'greenP25', 'greenP50', 'greenP75', 'greenP90']) {
    if (Math.abs(reference[percentile] - candidate[percentile]) > 2) {
      fail(`${label}: ${percentile} differs by more than 2 levels (${reference[percentile]} vs ${candidate[percentile]})`);
    }
  }
  if (Math.abs(reference.greenMean - candidate.greenMean) > 1.6) {
    fail(`${label}: mean green level differs by more than 1.6 (${reference.greenMean} vs ${candidate.greenMean})`);
  }
  if (Math.abs(reference.darkToGreen - candidate.darkToGreen) > 0.006) {
    fail(`${label}: dark-to-green ratio differs by more than 0.006 (${reference.darkToGreen} vs ${candidate.darkToGreen})`);
  }
  if (Math.abs(reference.darkMoleculeFraction - candidate.darkMoleculeFraction) > 0.004) {
    fail(`${label}: molecular dark-pixel coverage differs by more than 0.004 (${reference.darkMoleculeFraction} vs ${candidate.darkMoleculeFraction})`);
  }
  if (candidate.blackMoleculeFraction > reference.blackMoleculeFraction + 0.002) {
    fail(`${label}: pure-black coverage exceeds Mol* by more than 0.002 (${reference.blackMoleculeFraction} vs ${candidate.blackMoleculeFraction})`);
  }
  if (candidate.darkPixels === 0) {
    fail(`${label}: illustrative outline is absent`);
  }
  if (Math.abs(reference.backgroundFraction - candidate.backgroundFraction) > 0.012) {
    fail(`${label}: background fraction differs by more than 0.012 (${reference.backgroundFraction} vs ${candidate.backgroundFraction})`);
  }
  for (const field of ['x', 'y', 'width', 'height']) {
    if (Math.abs(reference.moleculeBounds[field] - candidate.moleculeBounds[field]) > 1) {
      fail(`${label}: molecule bound ${field} differs by more than 1 px (${reference.moleculeBounds[field]} vs ${candidate.moleculeBounds[field]})`);
    }
  }
  const spatial = spatialComparison(reference.path, candidate.path);
  if (spatial.foregroundJaccard < 0.72) {
    fail(`${label}: foreground Jaccard similarity is below 0.72 (${spatial.foregroundJaccard})`);
  }
  if (spatial.rgbMae > 0.16 || spatial.luminanceMae > 0.16) {
    fail(`${label}: spatial color error is too high (${JSON.stringify(spatial)})`);
  }
  return spatial;
}

if (process.argv[2] === '--metrics') {
  for (const path of process.argv.slice(3)) console.log(JSON.stringify(metrics(resolve(path)), null, 2));
  process.exit(0);
}

if (process.argv[2] === '--compare') {
  const reference = metrics(resolve(process.argv[3]));
  const candidate = metrics(resolve(process.argv[4]));
  console.log(JSON.stringify({ reference, candidate, spatial: spatialComparison(reference.path, candidate.path) }, null, 2));
  process.exit(0);
}

if (process.argv[2] === '--spatial') {
  const referencePath = resolve(process.argv[3]);
  for (const candidatePath of process.argv.slice(4)) {
    console.log(JSON.stringify({
      candidate: resolve(candidatePath),
      spatial: spatialComparison(referencePath, resolve(candidatePath)),
    }));
  }
  process.exit(0);
}

const outputDir = process.argv[2]
  ? resolve(process.argv[2])
  : mkdtempSync(join(tmpdir(), 'molfig-illustrative-qa.'));
const pairs = generateReferences(outputDir);
for (const [label, paths] of Object.entries(pairs)) {
  validateMolstarRenderingReport(paths.report);
  const reference = metrics(paths.molstar);
  const candidate = metrics(paths.molfig);
  const spatial = assertPair(label, reference, candidate);
  console.log(JSON.stringify({ label, reference, candidate, spatial }, null, 2));
}
console.log(`Illustrative rendering validation passed; artifacts: ${outputDir}`);
