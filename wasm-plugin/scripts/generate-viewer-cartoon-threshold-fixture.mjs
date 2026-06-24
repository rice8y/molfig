#!/usr/bin/env node

import { mkdirSync, writeFileSync } from 'node:fs';
import path from 'node:path';
import process from 'node:process';

const output = process.argv[2];
if (!output) {
  throw new Error('usage: generate-viewer-cartoon-threshold-fixture.mjs <output.cif>');
}

const waterCount = 50_001;
const lipidResidueCount = 10_001;
const lines = [
  'data_viewer_cartoon_thresholds',
  'loop_',
  '_entity.id',
  '_entity.type',
  '1 water',
  '2 non-polymer',
  '#',
  'loop_',
  '_struct_asym.id',
  '_struct_asym.entity_id',
  'W 1',
  'L 2',
  '#',
  'loop_',
  '_chem_comp.id',
  '_chem_comp.name',
  '_chem_comp.type',
  '_chem_comp.mon_nstd_flag',
  "HOH WATER non-polymer n",
  "LPX 'LIPID THRESHOLD COMPONENT' lipid n",
  '#',
  'loop_',
  '_chem_comp_bond.comp_id',
  '_chem_comp_bond.atom_id_1',
  '_chem_comp_bond.atom_id_2',
  '_chem_comp_bond.value_order',
  'LPX C1 C2 sing',
  '#',
  'loop_',
  '_atom_site.group_PDB',
  '_atom_site.id',
  '_atom_site.type_symbol',
  '_atom_site.label_atom_id',
  '_atom_site.label_comp_id',
  '_atom_site.label_asym_id',
  '_atom_site.label_entity_id',
  '_atom_site.label_seq_id',
  '_atom_site.auth_seq_id',
  '_atom_site.Cartn_x',
  '_atom_site.Cartn_y',
  '_atom_site.Cartn_z',
  '_atom_site.occupancy',
  '_atom_site.pdbx_PDB_model_num',
];

let atomId = 1;
for (let i = 0; i < waterCount; i++, atomId++) {
  const x = (i % 250) * 3;
  const y = Math.floor(i / 250) * 3;
  lines.push(`HETATM ${atomId} O O HOH W 1 . ${i + 1} ${x}.000 ${y}.000 0.000 1.00 1`);
}

for (let i = 0; i < lipidResidueCount; i++) {
  const x = (i % 200) * 5;
  const y = 700 + Math.floor(i / 200) * 5;
  const authSeq = i + 1;
  lines.push(`HETATM ${atomId++} C C1 LPX L 2 . ${authSeq} ${x}.000 ${y}.000 0.000 1.00 1`);
  lines.push(`HETATM ${atomId++} C C2 LPX L 2 . ${authSeq} ${x + 1}.400 ${y}.000 0.000 1.00 1`);
}
lines.push('#', '');

mkdirSync(path.dirname(path.resolve(output)), { recursive: true });
writeFileSync(output, lines.join('\n'));
console.log(`generated ${output}: ${waterCount} water atoms, ${lipidResidueCount * 2} lipid atoms`);
