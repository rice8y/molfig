import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const crateRoot = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const dataRoot = resolve(crateRoot, '../package/examples/data');
const manifestPath = join(dataRoot, 'XYZ_VALIDATION.json');
const manifest = JSON.parse(readFileSync(manifestPath, 'utf8'));

assertEqual(manifest.schema, 1, 'manifest schema');
assertEqual(manifest.coordinate_unit, 'angstrom', 'coordinate unit');
assert(Array.isArray(manifest.records) && manifest.records.length >= 4, 'expected at least four real PubChem XYZ records');

const seenFiles = new Set();
const seenCids = new Set();
for (const record of manifest.records) {
  assert(!seenFiles.has(record.file), `duplicate file ${record.file}`);
  assert(!seenCids.has(record.cid), `duplicate CID ${record.cid}`);
  seenFiles.add(record.file);
  seenCids.add(record.cid);

  const bytes = readFileSync(join(dataRoot, record.file));
  const text = bytes.toString('utf8');
  assert(text.endsWith('\n'), `${record.file}: final newline`);
  assert(!text.includes('\r'), `${record.file}: LF line endings`);
  assertEqual(createHash('sha256').update(bytes).digest('hex'), record.sha256, `${record.file}: sha256`);

  const lines = text.slice(0, -1).split('\n');
  assertEqual(lines.length, record.atom_count + 2, `${record.file}: exact line count`);
  assertEqual(lines[0], String(record.atom_count), `${record.file}: canonical atom-count line`);
  assert(lines[1].includes(`PubChem CID ${record.cid}`), `${record.file}: CID in comment`);
  assert(lines[1].includes(record.conformer_id), `${record.file}: conformer id in comment`);

  const elements = {};
  const positions = [];
  for (const [index, line] of lines.slice(2).entries()) {
    const fields = line.trim().split(/\s+/);
    assertEqual(fields.length, 4, `${record.file}: atom ${index + 1} field count`);
    assert(/^[A-Z][a-z]?$/.test(fields[0]), `${record.file}: atom ${index + 1} element symbol`);
    for (const coordinate of fields.slice(1)) {
      assert(/^-?\d+\.\d{4}$/.test(coordinate), `${record.file}: atom ${index + 1} canonical coordinate ${coordinate}`);
    }
    const xyz = fields.slice(1).map(Number);
    assert(xyz.every(Number.isFinite), `${record.file}: atom ${index + 1} finite coordinates`);
    elements[fields[0]] = (elements[fields[0]] ?? 0) + 1;
    positions.push(xyz);
  }

  const sortedElements = Object.fromEntries(Object.entries(elements).sort(([a], [b]) => a.localeCompare(b)));
  const sortedExpectedElements = Object.fromEntries(Object.entries(record.elements).sort(([a], [b]) => a.localeCompare(b)));
  assertEqual(sortedElements, sortedExpectedElements, `${record.file}: element composition`);
  assertEqual(parseFormula(record.formula), sortedExpectedElements, `${record.file}: molecular formula`);
  const bounds = {
    min: [0, 1, 2].map(axis => Math.min(...positions.map(position => position[axis]))),
    max: [0, 1, 2].map(axis => Math.max(...positions.map(position => position[axis]))),
  };
  assertEqual(bounds, record.bounds, `${record.file}: coordinate bounds`);

  for (let a = 0; a < positions.length; a++) {
    for (let b = a + 1; b < positions.length; b++) {
      const distance = Math.hypot(...positions[a].map((value, axis) => value - positions[b][axis]));
      assert(distance > 0.2, `${record.file}: atoms ${a + 1} and ${b + 1} are implausibly coincident`);
    }
  }

  assertEqual(record.source_bonds.length, record.source_bond_count, `${record.file}: source SDF bond count`);
  const sourceEndpoints = [];
  const seenSourceEndpoints = new Set();
  for (const [bondIndex, bond] of record.source_bonds.entries()) {
    assert(Array.isArray(bond) && bond.length === 3, `${record.file}: source bond ${bondIndex + 1} shape`);
    const [a, b, order] = bond;
    assert(Number.isInteger(a) && Number.isInteger(b) && a >= 1 && a < b && b <= record.atom_count,
      `${record.file}: source bond ${bondIndex + 1} canonical endpoints`);
    assert(Number.isInteger(order) && order >= 1 && order <= 3,
      `${record.file}: source bond ${bondIndex + 1} order`);
    const key = `${a}-${b}`;
    assert(!seenSourceEndpoints.has(key), `${record.file}: duplicate source bond ${key}`);
    seenSourceEndpoints.add(key);
    sourceEndpoints.push([a, b]);
  }
  const inferredEndpoints = inferMolstarBondEndpoints(lines.slice(2));
  assertEqual(inferredEndpoints, sourceEndpoints, `${record.file}: Mol* XYZ bond inference vs source SDF connectivity`);
}

console.log(`PubChem XYZ contract OK (${manifest.records.length} records)`);

function assert(condition, label) {
  if (!condition) throw new Error(label);
}

function assertEqual(actual, expected, label) {
  const actualJson = JSON.stringify(actual);
  const expectedJson = JSON.stringify(expected);
  if (actualJson !== expectedJson) {
    throw new Error(`${label}\nactual:   ${actualJson}\nexpected: ${expectedJson}`);
  }
}

function parseFormula(formula) {
  const elements = {};
  let consumed = '';
  for (const match of formula.matchAll(/([A-Z][a-z]?)(\d*)/g)) {
    consumed += match[0];
    elements[match[1]] = (elements[match[1]] ?? 0) + Number(match[2] || 1);
  }
  assertEqual(consumed, formula, `invalid molecular formula ${formula}`);
  return Object.fromEntries(Object.entries(elements).sort(([a], [b]) => a.localeCompare(b)));
}

function inferMolstarBondEndpoints(atomLines) {
  const atoms = atomLines.map(line => {
    const [element, x, y, z] = line.trim().split(/\s+/);
    return { element, position: [Number(x), Number(y), Number(z)] };
  });
  const endpoints = [];
  for (let a = 0; a < atoms.length; a++) {
    for (let b = a + 1; b < atoms.length; b++) {
      if (atoms[a].element === 'H' && atoms[b].element === 'H') continue;
      const distance = Math.hypot(...atoms[a].position.map((value, axis) => value - atoms[b].position[axis]));
      if (distance > 0 && distance <= molstarPairingThreshold(atoms[a].element, atoms[b].element)) {
        endpoints.push([a + 1, b + 1]);
      }
    }
  }
  return endpoints;
}

function molstarPairingThreshold(elementA, elementB) {
  const elementIndex = { H: 0, C: 6, N: 7, O: 8 };
  const elementThreshold = { H: 1.42, C: 1.75, N: 1.6, O: 1.52 };
  const a = elementIndex[elementA];
  const b = elementIndex[elementB];
  assert(a !== undefined && b !== undefined, `missing Mol* threshold for ${elementA}-${elementB}`);
  const key = (a + b) * (a + b + 1) / 2 + Math.max(a, b);
  const pairThreshold = new Map([
    [0, 0.8], [27, 1.2], [35, 1.15], [44, 1.1], [84, 1.75],
    [98, 1.6], [112, 1.6], [113, 1.59], [129, 1.45], [144, 1.6],
  ]).get(key);
  return pairThreshold ?? (elementThreshold[elementA] + elementThreshold[elementB]) / 1.95;
}
