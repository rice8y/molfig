// Convert the official QMEAN example, without inventing coordinates or scores.
// Run: node package/examples/data/convert-qmean.mjs [--check]
import assert from 'node:assert/strict';
import { readFileSync, writeFileSync } from 'node:fs';

const file = name => new URL(name, import.meta.url);
const project = JSON.parse(readFileSync(file('QMEAN.json'), 'utf8'));
assert.equal(project.method, 'QMEAN');
const model = project.models.model_001;
const amino = Object.fromEntries('A ALA,R ARG,N ASN,D ASP,C CYS,Q GLN,E GLU,G GLY,H HIS,I ILE,L LEU,K LYS,M MET,F PHE,P PRO,S SER,T THR,W TRP,Y TYR,V VAL'.split(',').map(pair => pair.split(' ')));
assert.deepEqual(Object.keys(model.chains), ['A']);
const sequence = model.chains.A.seqres;
const scores = model.scores.local_scores.A;
assert.equal(sequence.length, scores.length);
const lines = readFileSync(file('qmean-model_001_processed.pdb'), 'utf8').split(/\r?\n/);
assert(!lines.some(line => line.startsWith('MODEL ')), 'Expected one unnumbered PDB model');
const atoms = lines.filter(line => line.startsWith('ATOM  ') || line.startsWith('HETATM'));
assert(atoms.length > 0);
const residues = new Map();
const atomIds = new Set();
const atomRows = atoms.map(line => {
  const field = (start, end) => line.slice(start, end).trim();
  const id = field(6, 11);
  const seq = Number(field(22, 26));
  const comp = field(17, 20);
  assert.equal(field(0, 6), 'ATOM');
  assert.equal(field(21, 22), 'A');
  assert.equal(field(26, 27), '', 'Insertion codes require explicit remapping');
  assert.equal(field(16, 17), '', 'Alternate locations require explicit remapping');
  assert(!atomIds.has(id));
  atomIds.add(id);
  assert.equal(comp, amino[sequence[seq - 1]]);
  const score = scores[seq - 1];
  assert(Number.isFinite(score) && score >= 0 && score <= 1);
  // The official processed PDB stores the same scores rounded to 2 decimals.
  assert(Math.abs(Number(field(60, 66)) - score) <= 0.00501);
  for (const [start, end] of [[30, 38], [38, 46], [46, 54]]) {
    assert(Number.isFinite(Number(field(start, end))));
  }
  assert(/^[A-Z][A-Z]?$/.test(field(76, 78)));
  residues.set(seq, comp);
  return ['ATOM', id, field(76, 78), field(12, 16), '.', comp, 'A', '1', seq, '?',
    field(30, 38), field(38, 46), field(46, 54), field(54, 60), field(60, 66),
    seq, comp, 'A', field(12, 16), '1'].join(' ');
});
assert.equal(residues.size, scores.filter(score => score !== null).length);
let output = '# Derived from the official SWISS-MODEL QMEAN example model_001.\n'
  + '# Coordinates preserved from qmean-model_001_processed.pdb.\n'
  + '# Full-precision local scores preserved from QMEAN.json.\n'
  + '# Conversion only; no synthetic annotations. See README.md. CC BY-SA 4.0.\n'
  + 'data_qmean_model_001\n#\n';
function loop(category, columns, rows) {
  output += 'loop_\n' + columns.map(column => `_${category}.${column}`).join('\n') + '\n'
    + rows.map(row => Array.isArray(row) ? row.join(' ') : row).join('\n') + '\n#\n';
}
loop('entity', ['id', 'type', 'pdbx_description'], [['1', 'polymer', "'Official QMEAN example model_001'"]]);
loop('entity_poly', ['entity_id', 'type', 'pdbx_seq_one_letter_code', 'pdbx_seq_one_letter_code_can', 'pdbx_strand_id'], [['1', "'polypeptide(L)'", sequence, sequence, 'A']]);
loop('entity_poly_seq', ['entity_id', 'num', 'mon_id', 'hetero'], [...sequence].map((letter, i) => ['1', i + 1, amino[letter], 'n']));
loop('struct_asym', ['id', 'entity_id'], [['A', '1']]);
loop('atom_site', ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id', 'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv', 'auth_seq_id', 'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num'], atomRows);
loop('ma_qa_metric', ['id', 'mode', 'name', 'type'], [['1', 'local', 'QMEAN', 'other']]);
loop('ma_qa_metric_local', ['ordinal_id', 'model_id', 'label_asym_id', 'label_seq_id', 'label_comp_id', 'metric_id', 'metric_value'], [...residues].map(([seq, comp], i) => [i + 1, '1', 'A', seq, comp, '1', scores[seq - 1]]));
const target = file('qmean-model_001.cif');
if (process.argv.includes('--check')) {
  assert.equal(readFileSync(target, 'utf8'), output, 'Regenerate qmean-model_001.cif');
} else {
  writeFileSync(target, output);
}
console.log(`QMEAN: ${atoms.length} atoms, ${residues.size} scored residues; coordinates, sequence and scores verified.`);
