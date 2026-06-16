import { readFileSync } from 'node:fs';
import { createHash } from 'node:crypto';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = join(dirname(fileURLToPath(import.meta.url)), '..');
const contract = JSON.parse(readFileSync(join(root, 'tests/expected/mesh-contract.json'), 'utf8'));

function parsePdb(text) {
  const atoms = [];
  const bonds = new Set();
  let model = 1;
  for (const line of text.split(/\r?\n/)) {
    const record = line.slice(0, 6).trim();
    if (record === 'MODEL') {
      model = Number(line.slice(10, 14).trim()) || model;
    } else if (record === 'ATOM' || record === 'HETATM') {
      atoms.push({
        serial: Number(line.slice(6, 11)),
        atom: line.slice(12, 16).trim(),
        altLoc: normalizeOptional(line.slice(16, 17).trim()),
        residue: line.slice(17, 20).trim(),
        chain: line.slice(21, 22).trim(),
        model,
        x: Number(line.slice(30, 38)),
        y: Number(line.slice(38, 46)),
        z: Number(line.slice(46, 54)),
        element: line.slice(76, 78).trim()
      });
    } else if (record === 'CONECT') {
      const from = Number(line.slice(6, 11));
      for (let i = 11; i < line.length; i += 5) {
        const to = Number(line.slice(i, i + 5));
        if (Number.isFinite(to) && to > 0) {
          bonds.add([from, to].sort((a, b) => a - b).join('-'));
        }
      }
    }
  }
  return summarize(atoms, bonds.size);
}

function parseCif(text) {
  const lines = text.split(/\r?\n/);
  const atoms = [];
  const bonds = new Set();
  for (let i = 0; i < lines.length; i++) {
    if (lines[i].trim() !== 'loop_') continue;
    const fields = [];
    let j = i + 1;
    while (j < lines.length && lines[j].trim().startsWith('_')) {
      fields.push(lines[j].trim());
      j++;
    }
    const rows = [];
    while (j < lines.length && lines[j].trim() && !lines[j].trim().startsWith('#')) {
      rows.push(lines[j].trim().split(/\s+/));
      j++;
    }
    if (fields.includes('_atom_site.id')) {
      for (const row of rows) {
        const get = name => row[fields.indexOf(name)];
        const getOptional = (name, fallback = '') => {
          const index = fields.indexOf(name);
          return index >= 0 ? row[index] : fallback;
        };
        atoms.push({
          serial: Number(get('_atom_site.id')),
          atom: get('_atom_site.label_atom_id'),
          altLoc: normalizeOptional(getOptional('_atom_site.label_alt_id')),
          residue: get('_atom_site.label_comp_id'),
          chain: get('_atom_site.label_asym_id'),
          model: Number(getOptional('_atom_site.pdbx_PDB_model_num', '1')) || 1,
          x: Number(get('_atom_site.Cartn_x')),
          y: Number(get('_atom_site.Cartn_y')),
          z: Number(get('_atom_site.Cartn_z')),
          element: get('_atom_site.type_symbol')
        });
      }
    }
    if (fields.includes('_struct_conn.ptnr1_label_atom_id')) {
      const key = (...names) => names.map(name => fields.indexOf(name));
      const [aChain, aSeq, aAtom, bChain, bSeq, bAtom] = key(
        '_struct_conn.ptnr1_label_asym_id',
        '_struct_conn.ptnr1_label_seq_id',
        '_struct_conn.ptnr1_label_atom_id',
        '_struct_conn.ptnr2_label_asym_id',
        '_struct_conn.ptnr2_label_seq_id',
        '_struct_conn.ptnr2_label_atom_id'
      );
      for (const row of rows) {
        bonds.add([
          `${row[aChain]}:${row[aSeq]}:${row[aAtom]}`,
          `${row[bChain]}:${row[bSeq]}:${row[bAtom]}`
        ].sort().join('-'));
      }
    }
  }
  return summarize(atoms, bonds.size);
}

function normalizeOptional(value) {
  return value === '.' || value === '?' ? '' : value;
}

function summarize(atoms, bondCount) {
  const byElementRaw = {};
  for (const atom of atoms) byElementRaw[atom.element] = (byElementRaw[atom.element] ?? 0) + 1;
  const byElement = Object.fromEntries(Object.entries(byElementRaw).sort(([a], [b]) => a.localeCompare(b)));
  const altLocs = [];
  for (const atom of atoms) {
    if (atom.altLoc && !altLocs.includes(atom.altLoc)) altLocs.push(atom.altLoc);
  }
  return {
    atomCount: atoms.length,
    bondCount,
    models: new Set(atoms.map(atom => atom.model)).size,
    residues: [...new Set(atoms.map(atom => atom.residue))],
    chains: [...new Set(atoms.map(atom => atom.chain))],
    elements: byElement,
    bboxAngstrom: {
      min: [
        Math.min(...atoms.map(atom => atom.x)),
        Math.min(...atoms.map(atom => atom.y)),
        Math.min(...atoms.map(atom => atom.z))
      ],
      max: [
        Math.max(...atoms.map(atom => atom.x)),
        Math.max(...atoms.map(atom => atom.y)),
        Math.max(...atoms.map(atom => atom.z))
      ]
    },
    altLocs
  };
}

function assertDeepEqual(actual, expected, label) {
  const a = JSON.stringify(actual);
  const e = JSON.stringify(expected);
  if (a !== e) throw new Error(`${label}\nactual:   ${a}\nexpected: ${e}`);
}

for (const [name, fixture] of Object.entries(contract.fixtures)) {
  for (const input of fixture.inputs) {
    const text = readFileSync(join(root, input), 'utf8');
    const summary = input.endsWith('.pdb') ? parsePdb(text) : parseCif(text);
    const expected = fixture.parse;
    assertDeepEqual(summary.atomCount, expected.atomCount, `${name} atom count in ${input}`);
    assertDeepEqual(summary.bondCount, expected.bondCount, `${name} bond count in ${input}`);
    if ('models' in expected) assertDeepEqual(summary.models, expected.models, `${name} model count in ${input}`);
    assertDeepEqual(summary.residues, expected.residues, `${name} residues in ${input}`);
    assertDeepEqual(summary.chains, expected.chains, `${name} chains in ${input}`);
    assertDeepEqual(summary.elements, expected.elements, `${name} elements in ${input}`);
    assertDeepEqual(summary.bboxAngstrom, expected.bboxAngstrom, `${name} bbox in ${input}`);
    if ('altLocs' in expected) assertDeepEqual(summary.altLocs, expected.altLocs, `${name} altLocs in ${input}`);
  }

  for (const input of fixture.binaryInputs ?? []) {
    const bytes = readFileSync(join(root, input));
    const hash = createHash('sha256').update(bytes).digest('hex');
    const checksumPath = `${input}.sha256`;
    const expectedHash = readFileSync(join(root, checksumPath), 'utf8').trim().split(/\s+/)[0];
    assertDeepEqual(bytes.length > 0, true, `${name} non-empty binary fixture in ${input}`);
    assertDeepEqual(hash, expectedHash, `${name} checksum in ${checksumPath}`);
  }
}

console.log('Fixture contract OK');
