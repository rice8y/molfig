import { createHash } from 'node:crypto';
import { readFileSync, readdirSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const here = dirname(fileURLToPath(import.meta.url));

const atoms = [
  ['ATOM', '1', 'N', 'N', '.', 'ALA', 'A', '1', '0.000', '0.000', '0.000', '1.00', '11.00', '0', '1', '1'],
  ['ATOM', '2', 'C', 'CA', '.', 'ALA', 'A', '1', '1.450', '0.050', '0.100', '1.00', '11.00', '1', '1', '1'],
  ['ATOM', '3', 'C', 'C', '.', 'ALA', 'A', '1', '2.050', '1.390', '0.000', '1.00', '11.00', '0', '1', '1'],
  ['ATOM', '4', 'O', 'O', '.', 'ALA', 'A', '1', '1.500', '2.460', '0.100', '1.00', '11.00', '0', '1', '1'],
  ['ATOM', '5', 'C', 'CB', 'A', 'ALA', 'A', '1', '1.950', '-0.820', '1.240', '0.60', '14.00', '0', '1', '1'],
  ['ATOM', '6', 'C', 'CB', 'B', 'ALA', 'A', '1', '1.900', '-0.900', '-1.120', '0.40', '18.00', '0', '1', '1'],
  ['ATOM', '7', 'N', 'N', '.', 'SER', 'A', '2', '3.220', '1.330', '-0.160', '1.00', '12.00', '0', '1', '1'],
  ['ATOM', '8', 'C', 'CA', '.', 'SER', 'A', '2', '3.900', '2.570', '-0.280', '1.00', '12.00', '0', '1', '1'],
  ['ATOM', '9', 'C', 'C', '.', 'SER', 'A', '2', '5.370', '2.470', '-0.020', '1.00', '12.00', '0', '1', '1'],
  ['ATOM', '10', 'O', 'O', '.', 'SER', 'A', '2', '6.030', '3.480', '0.040', '1.00', '12.00', '0', '1', '1'],
  ['ATOM', '11', 'N', 'N', '.', 'GLY', 'A', '3', '5.900', '1.240', '0.120', '1.00', '13.00', '0', '1', '1'],
  ['ATOM', '12', 'C', 'CA', '.', 'GLY', 'A', '3', '7.300', '1.020', '0.420', '1.00', '13.00', '0', '1', '1'],
  ['ATOM', '13', 'C', 'C', '.', 'GLY', 'A', '3', '8.100', '2.230', '0.020', '1.00', '13.00', '0', '1', '1'],
  ['ATOM', '14', 'O', 'O', '.', 'GLY', 'A', '3', '7.700', '3.350', '0.240', '1.00', '13.00', '0', '1', '1'],
].map(row => {
  const [group, id, typeSymbol, labelAtomId, labelAltId, labelCompId, labelAsymId, labelSeqId, x, y, z, occupancy, bIso, charge, modelNum, ihmModelId] = row;
  const insertionCodeBySeq = { '1': 'A', '2': '.', '3': 'B' };
  return [
    group, id, typeSymbol, labelAtomId, `${labelAtomId}X`, labelAltId,
    labelCompId, `${labelCompId}X`, labelAsymId, 'X', '1', labelSeqId,
    `${Number(labelSeqId) + 100}`, insertionCodeBySeq[labelSeqId] ?? '.',
    x, y, z, occupancy, bIso, charge, modelNum, ihmModelId,
  ];
});

const atomFields = [
  'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'auth_atom_id',
  'label_alt_id', 'label_comp_id', 'auth_comp_id', 'label_asym_id',
  'auth_asym_id', 'label_entity_id', 'label_seq_id', 'auth_seq_id',
  'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
  'B_iso_or_equiv', 'pdbx_formal_charge', 'pdbx_PDB_model_num', 'ihm_model_id',
];

const structConns = [
  ['covale1', 'covale', 'A', '1', 'N', 'A', '1', 'CA'],
  ['covale2', 'covale', 'A', '1', 'CA', 'A', '1', 'C'],
  ['covale3', 'covale', 'A', '1', 'CA', 'A', '1', 'CB'],
  ['covale4', 'covale', 'A', '1', 'C', 'A', '1', 'O'],
  ['covale5', 'covale', 'A', '1', 'C', 'A', '2', 'N'],
  ['covale6', 'covale', 'A', '2', 'N', 'A', '2', 'CA'],
  ['covale7', 'covale', 'A', '2', 'CA', 'A', '2', 'C'],
  ['covale8', 'covale', 'A', '2', 'C', 'A', '2', 'O'],
  ['covale9', 'covale', 'A', '2', 'C', 'A', '3', 'N'],
  ['covale10', 'covale', 'A', '3', 'N', 'A', '3', 'CA'],
  ['covale11', 'covale', 'A', '3', 'CA', 'A', '3', 'C'],
  ['covale12', 'covale', 'A', '3', 'C', 'A', '3', 'O'],
].map(row => {
  const insertionCodeBySeq = { '1': 'A', '2': '.', '3': 'B' };
  const [id, connTypeId, p1Chain, p1Seq, p1Atom, p2Chain, p2Seq, p2Atom] = row;
  return [
    id, connTypeId, p1Chain, p1Seq, insertionCodeBySeq[p1Seq] ?? '.', p1Atom,
    p2Chain, p2Seq, insertionCodeBySeq[p2Seq] ?? '.', p2Atom,
  ];
});

const connFields = [
  'id', 'conn_type_id', 'ptnr1_label_asym_id', 'ptnr1_label_seq_id',
  'pdbx_ptnr1_PDB_ins_code', 'ptnr1_label_atom_id', 'ptnr2_label_asym_id',
  'ptnr2_label_seq_id', 'pdbx_ptnr2_PDB_ins_code', 'ptnr2_label_atom_id',
];

const chemCompAtoms = [
  ['ALA', 'N', 'N', 'N', '0', 'N', 'N', 'N', '0.000', '0.000', '0.000', '0.000', '0.000', '0.000'],
  ['ALA', 'CA', 'CA', 'C', '1', 'N', 'N', 'N', '1.450', '0.050', '0.100', '1.450', '0.050', '0.100'],
  ['ALA', 'C', 'C', 'C', '0', 'N', 'N', 'N', '2.050', '1.390', '0.000', '2.050', '1.390', '0.000'],
  ['ALA', 'O', 'O', 'O', '0', 'N', 'N', 'N', '1.500', '2.460', '0.100', '1.500', '2.460', '0.100'],
  ['ALA', 'CB', 'CB', 'C', '0', 'N', 'N', 'N', '1.950', '-0.820', '1.240', '1.950', '-0.820', '1.240'],
];

const chemCompAtomFields = [
  'comp_id', 'atom_id', 'alt_atom_id', 'type_symbol', 'charge',
  'pdbx_aromatic_flag', 'pdbx_leaving_atom_flag', 'pdbx_stereo_config',
  'model_Cartn_x', 'model_Cartn_y', 'model_Cartn_z',
  'pdbx_model_Cartn_x_ideal', 'pdbx_model_Cartn_y_ideal', 'pdbx_model_Cartn_z_ideal',
];

const assemblyOps = [
  ['1', 'identity', '1.000000', '0.000000', '0.000000', '0.00000', '0.000000', '1.000000', '0.000000', '0.00000', '0.000000', '0.000000', '1.000000', '0.00000'],
  ['2', 'rotation', '-1.000000', '0.000000', '0.000000', '12.00000', '0.000000', '1.000000', '0.000000', '0.00000', '0.000000', '0.000000', '1.000000', '0.00000'],
];

const opFields = [
  'id', 'type', 'matrix[1][1]', 'matrix[1][2]', 'matrix[1][3]', 'vector[1]',
  'matrix[2][1]', 'matrix[2][2]', 'matrix[2][3]', 'vector[2]',
  'matrix[3][1]', 'matrix[3][2]', 'matrix[3][3]', 'vector[3]',
];

function stringData(values) {
  const unique = [];
  const map = new Map();
  const indices = [];
  for (const value of values) {
    if (!map.has(value)) {
      map.set(value, unique.length);
      unique.push(value);
    }
    indices.push(map.get(value));
  }

  let stringData = '';
  const offsets = [0];
  for (const value of unique) {
    stringData += value;
    offsets.push(stringData.length);
  }

  return {
    encoding: [{
      kind: 'StringArray',
      dataEncoding: [{ kind: 'ByteArray', type: 4 }],
      stringData,
      offsetEncoding: [{ kind: 'ByteArray', type: 4 }],
      offsets: Uint8Array.from(offsets),
    }],
    data: Uint8Array.from(indices),
  };
}

function intData(values) {
  const bytes = Buffer.alloc(values.length * 4);
  values.forEach((value, index) => bytes.writeInt32LE(Number(value), index * 4));
  return {
    encoding: [{ kind: 'ByteArray', type: 3 }],
    data: bytes,
  };
}

function floatData(values) {
  const bytes = Buffer.alloc(values.length * 4);
  values.forEach((value, index) => bytes.writeFloatLE(Number(value), index * 4));
  return {
    encoding: [{ kind: 'ByteArray', type: 32 }],
    data: bytes,
  };
}

function category(name, fields, rows, typed = {}) {
  return {
    name,
    rowCount: rows.length,
    columns: fields.map((field, index) => ({
      name: field,
      data: (typed[field] || stringData)(rows.map(row => row[index])),
    })),
  };
}

const atomSiteTyped = {
  id: intData,
  label_seq_id: intData,
  Cartn_x: floatData,
  Cartn_y: floatData,
  Cartn_z: floatData,
  occupancy: floatData,
  B_iso_or_equiv: floatData,
  pdbx_formal_charge: intData,
  pdbx_PDB_model_num: intData,
  ihm_model_id: intData,
};

const chemCompTyped = {
  formula_weight: floatData,
  pdbx_formal_charge: intData,
};

const chemCompAtomTyped = {
  charge: intData,
  model_Cartn_x: floatData,
  model_Cartn_y: floatData,
  model_Cartn_z: floatData,
  pdbx_model_Cartn_x_ideal: floatData,
  pdbx_model_Cartn_y_ideal: floatData,
  pdbx_model_Cartn_z_ideal: floatData,
};

const chemCompBondTyped = {
  pdbx_ordinal: intData,
};

const chemCompAngleTyped = {
  value_angle: floatData,
  value_angle_esd: floatData,
};

const anisotropTyped = {
  id: intData,
  'U[1][1]': floatData,
  'U[1][2]': floatData,
  'U[1][3]': floatData,
  'U[2][1]': floatData,
  'U[2][2]': floatData,
  'U[2][3]': floatData,
  'U[3][1]': floatData,
  'U[3][2]': floatData,
  'U[3][3]': floatData,
};

const ihmSphereFields = [
  'id', 'entity_id', 'asym_id', 'seq_id_begin', 'seq_id_end',
  'Cartn_x', 'Cartn_y', 'Cartn_z', 'object_radius',
];

const ihmSphereTyped = {
  id: intData,
  seq_id_begin: intData,
  seq_id_end: intData,
  Cartn_x: floatData,
  Cartn_y: floatData,
  Cartn_z: floatData,
  object_radius: floatData,
};

const ihmGaussianFields = [
  'id', 'entity_id', 'asym_id', 'seq_id_begin', 'seq_id_end',
  'mean_Cartn_x', 'mean_Cartn_y', 'mean_Cartn_z', 'weight',
  'covariance_matrix[1][1]', 'covariance_matrix[1][2]', 'covariance_matrix[1][3]',
  'covariance_matrix[2][1]', 'covariance_matrix[2][2]', 'covariance_matrix[2][3]',
  'covariance_matrix[3][1]', 'covariance_matrix[3][2]', 'covariance_matrix[3][3]',
];

const ihmGaussianTyped = {
  id: intData,
  seq_id_begin: intData,
  seq_id_end: intData,
  mean_Cartn_x: floatData,
  mean_Cartn_y: floatData,
  mean_Cartn_z: floatData,
  weight: floatData,
  'covariance_matrix[1][1]': floatData,
  'covariance_matrix[1][2]': floatData,
  'covariance_matrix[1][3]': floatData,
  'covariance_matrix[2][1]': floatData,
  'covariance_matrix[2][2]': floatData,
  'covariance_matrix[2][3]': floatData,
  'covariance_matrix[3][1]': floatData,
  'covariance_matrix[3][2]': floatData,
  'covariance_matrix[3][3]': floatData,
};

const ihmModelListTyped = {
  model_id: intData,
  assembly_id: intData,
  protocol_id: intData,
  representation_id: intData,
};

const ihmModelGroupTyped = {
  id: intData,
};

const ihmModelGroupLinkTyped = {
  model_id: intData,
  group_id: intData,
};

const ihmCrossLinkRestraintTyped = {
  id: intData,
  group_id: intData,
  seq_id_1: intData,
  seq_id_2: intData,
  distance_threshold: floatData,
  psi: floatData,
  sigma_1: floatData,
  sigma_2: floatData,
};

function str(value) {
  const bytes = Buffer.from(value);
  if (bytes.length < 32) return Buffer.concat([Buffer.from([0xa0 | bytes.length]), bytes]);
  if (bytes.length < 256) return Buffer.concat([Buffer.from([0xd9, bytes.length]), bytes]);
  throw new Error(`string too long: ${value}`);
}

function bin(value) {
  const bytes = Buffer.from(value);
  if (bytes.length < 256) return Buffer.concat([Buffer.from([0xc4, bytes.length]), bytes]);
  throw new Error(`binary too long: ${bytes.length}`);
}

function uint(value) {
  if (value < 128) return Buffer.from([value]);
  if (value < 256) return Buffer.from([0xcc, value]);
  const out = Buffer.alloc(3);
  out[0] = 0xcd;
  out.writeUInt16BE(value, 1);
  return out;
}

function array(values) {
  const head = values.length < 16 ? Buffer.from([0x90 | values.length]) : Buffer.from([0xdc, values.length >> 8, values.length & 0xff]);
  return Buffer.concat([head, ...values.map(pack)]);
}

function object(value) {
  const entries = Object.entries(value);
  const head = entries.length < 16 ? Buffer.from([0x80 | entries.length]) : Buffer.from([0xde, entries.length >> 8, entries.length & 0xff]);
  return Buffer.concat([head, ...entries.flatMap(([key, val]) => [str(key), pack(val)])]);
}

function pack(value) {
  if (typeof value === 'string') return str(value);
  if (typeof value === 'number') return uint(value);
  if (value instanceof Uint8Array) return bin(value);
  if (Array.isArray(value)) return array(value);
  if (value && typeof value === 'object') return object(value);
  throw new Error(`unsupported value: ${value}`);
}

const file = {
  version: '0.3.0',
  encoder: 'molfig fixture generator',
  dataBlocks: [{
    header: 'assembly_altloc_helix',
    categories: [
      category('_entry', ['id'], [['assembly-altloc-helix']]),
      category('_exptl', ['method'], [['X-RAY DIFFRACTION']]),
      category('_entity', ['id', 'type', 'pdbx_description'], [['1', 'polymer', 'test peptide']]),
      category('_entity_poly', [
        'entity_id', 'type', 'pdbx_seq_one_letter_code', 'nstd_linkage', 'nstd_monomer',
      ], [['1', 'polypeptide(L)', 'ASG', 'no', 'no']]),
      category('_entity_poly_seq', ['entity_id', 'num', 'mon_id', 'hetero'], [
        ['1', '1', 'ALA', 'n'],
        ['1', '2', 'SER', 'n'],
        ['1', '3', 'GLY', 'n'],
      ]),
      category('_struct_asym', ['id', 'entity_id', 'details'], [['A', '1', 'test asym']]),
      category('_pdbx_entity_branch', ['entity_id', 'type'], [['2', 'oligosaccharide']]),
      category('_pdbx_entity_branch_link', [
        'link_id', 'details', 'entity_id', 'entity_branch_list_num_1',
        'entity_branch_list_num_2', 'comp_id_1', 'comp_id_2', 'atom_id_1',
        'leaving_atom_id_1', 'atom_stereo_config_1', 'atom_id_2',
        'leaving_atom_id_2', 'atom_stereo_config_2', 'value_order',
      ], [[
        '1', 'test glycosidic link', '2', '1', '2', 'NAG', 'MAN', 'C1',
        'O1', 'n', 'O4', 'HO4', 'n', 'sing',
      ]], {
        link_id: intData,
        entity_branch_list_num_1: intData,
        entity_branch_list_num_2: intData,
      }),
      category('_pdbx_branch_scheme', [
        'entity_id', 'hetero', 'asym_id', 'mon_id', 'num', 'pdb_asym_id',
        'pdb_seq_num', 'pdb_mon_id', 'auth_asym_id', 'auth_seq_num', 'auth_mon_id',
      ], [['2', 'n', 'B', 'NAG', '1', 'B', '101', 'NAG', 'BA', '501', 'NAG']], {
        num: intData,
      }),
      category('_pdbx_nonpoly_scheme', [
        'asym_id', 'entity_id', 'mon_id', 'pdb_strand_id', 'ndb_seq_num',
        'pdb_seq_num', 'auth_seq_num', 'pdb_mon_id', 'auth_mon_id', 'pdb_ins_code',
      ], [['L', '3', 'HEM', 'L', '1', '701', '9001', 'HEM', 'HEM', '.']]),
      category('_pdbx_poly_seq_scheme', [
        'asym_id', 'entity_id', 'seq_id', 'mon_id', 'ndb_seq_num', 'pdb_seq_num',
        'auth_seq_num', 'pdb_mon_id', 'auth_mon_id', 'pdb_strand_id', 'pdb_ins_code',
        'hetero',
      ], [['A', '1', '1', 'ALA', '1', '1', '10', 'ALA', 'ALA', 'A', '.', 'n']], {
        seq_id: intData,
      }),
      category('_chem_comp', [
        'id', 'name', 'type', 'formula', 'formula_weight', 'one_letter_code',
        'three_letter_code', 'mon_nstd_flag', 'pdbx_synonyms', 'pdbx_formal_charge',
        'pdbx_initial_date', 'pdbx_modified_date', 'pdbx_ambiguous_flag',
        'pdbx_release_status',
      ], [
        ['ALA', 'ALANINE', 'l-peptide linking', 'C3 H7 N O2', '89.09', 'A', 'ALA', 'n', 'ALANINE', '0', '1999-07-08', '2024-01-01', 'N', 'REL'],
        ['SER', 'SERINE', 'l-peptide linking', 'C3 H7 N O3', '105.09', 'S', 'SER', 'n', 'SERINE', '0', '1999-07-08', '2024-01-01', 'N', 'REL'],
        ['GLY', 'GLYCINE', 'l-peptide linking', 'C2 H5 N O2', '75.07', 'G', 'GLY', 'n', 'GLYCINE', '0', '1999-07-08', '2024-01-01', 'N', 'REL'],
      ], chemCompTyped),
      category('_chem_comp_atom', chemCompAtomFields, chemCompAtoms, chemCompAtomTyped),
      category('_chem_comp_bond', [
        'comp_id', 'atom_id_1', 'atom_id_2', 'value_order', 'pdbx_aromatic_flag',
        'pdbx_stereo_config', 'pdbx_ordinal',
      ], [
        ['ALA', 'N', 'CA', 'SING', 'N', 'N', '101'],
        ['ALA', 'CA', 'C', 'sing', 'N', 'N', '102'],
        ['ALA', 'C', 'O', 'DOUB', 'N', 'N', '103'],
        ['ALA', 'CA', 'CB', 'delo', 'N', 'N', '104'],
      ], chemCompBondTyped),
      category('_chem_comp_angle', [
        'comp_id', 'atom_id_1', 'atom_id_2', 'atom_id_3', 'value_angle', 'value_angle_esd',
      ], [
        ['ALA', 'N', 'CA', 'C', '111.0', '1.5'],
        ['ALA', 'CA', 'C', 'O', '120.5', '2.0'],
      ], chemCompAngleTyped),
      category('_atom_site', atomFields, atoms, atomSiteTyped),
      category('_atom_site_anisotrop', [
        'id', 'U[1][1]', 'U[1][2]', 'U[1][3]', 'U[2][1]',
        'U[2][2]', 'U[2][3]', 'U[3][1]', 'U[3][2]', 'U[3][3]',
      ], [['2', '0.10', '0.01', '0.02', '0.01', '0.11', '0.03', '0.02', '0.03', '0.12']], anisotropTyped),
      category('_struct_conn', connFields, structConns),
      category('_struct_conf', [
        'conf_type_id', 'id', 'beg_label_comp_id', 'beg_label_asym_id',
        'beg_label_seq_id', 'pdbx_beg_PDB_ins_code', 'end_label_comp_id',
        'end_label_asym_id', 'end_label_seq_id', 'pdbx_end_PDB_ins_code',
        'pdbx_PDB_helix_class', 'details',
      ], [['HELX_P', 'H1', 'ALA', 'A', '1', 'A', 'GLY', 'A', '3', 'B', '1', 'molfig-test-helix']]),
      category('_struct_sheet_range', [
        'sheet_id', 'id', 'beg_label_asym_id', 'beg_label_seq_id',
        'pdbx_beg_PDB_ins_code', 'end_label_asym_id', 'end_label_seq_id',
        'pdbx_end_PDB_ins_code',
      ], [['S1', '1', 'A', '2', '.', 'A', '3', 'B']]),
      category('_pdbx_struct_assembly', ['id', 'details'], [['1', 'author_defined_test_dimer']]),
      category('_pdbx_struct_assembly_gen', ['assembly_id', 'oper_expression', 'asym_id_list'], [['1', '1,2', 'A']]),
      category('_pdbx_struct_oper_list', opFields, assemblyOps),
    ],
  }],
};

const ihmFile = {
  version: '0.3.0',
  encoder: 'molfig fixture generator',
  dataBlocks: [{
    header: 'ihm_only',
    categories: [
      category('_ihm_model_list', [
        'model_id', 'model_name', 'assembly_id', 'protocol_id', 'representation_id',
      ], [['1', 'model one', '1', '1', '1']], ihmModelListTyped),
      category('_ihm_model_group', ['id', 'name', 'details'], [['1', 'ensemble', 'test group']], ihmModelGroupTyped),
      category('_ihm_model_group_link', ['model_id', 'group_id'], [['1', '1']], ihmModelGroupLinkTyped),
      category('_ihm_sphere_obj_site', ihmSphereFields, [
        ['1', '1', 'A', '1', '10', '0.0', '0.0', '0.0', '2.0'],
      ], ihmSphereTyped),
      category('_ihm_gaussian_obj_site', ihmGaussianFields, [
        ['1', '1', 'A', '11', '20', '4.0', '0.0', '0.0', '1.0', '4.0', '0.0', '0.0', '0.0', '1.0', '0.0', '0.0', '0.0', '1.0'],
      ], ihmGaussianTyped),
      category('_ihm_cross_link_restraint', [
        'id', 'group_id', 'entity_id_1', 'entity_id_2', 'asym_id_1', 'asym_id_2',
        'comp_id_1', 'comp_id_2', 'seq_id_1', 'seq_id_2', 'atom_id_1', 'atom_id_2',
        'restraint_type', 'conditional_crosslink_flag', 'model_granularity',
        'distance_threshold', 'psi', 'sigma_1', 'sigma_2',
      ], [[
        '1', '1', '1', '1', 'A', 'A', 'ALA', 'GLY', '1', '10', 'CA', 'CA',
        'upper bound', 'ALL', 'by-residue', '25.0', '0.1', '1.5', '2.5',
      ]], ihmCrossLinkRestraintTyped),
    ],
  }],
};

function writeBcif(filename, value) {
  const bytes = pack(value);
  const out = join(here, filename);
  writeFileSync(out, bytes);
  writeFileSync(`${out}.sha256`, `${createHash('sha256').update(bytes).digest('hex')}  ${filename}\n`);
}

function tokenizeCif(text) {
  const tokens = [];
  let i = 0;
  while (i < text.length) {
    const ch = text[i];
    if (/\s/.test(ch)) {
      i++;
      continue;
    }
    if (ch === '#') {
      while (i < text.length && text[i] !== '\n') i++;
      continue;
    }
    if ((ch === '"' || ch === "'")) {
      const quote = ch;
      i++;
      let value = '';
      while (i < text.length && text[i] !== quote) value += text[i++];
      if (text[i] === quote) i++;
      tokens.push(value);
      continue;
    }
    let value = '';
    while (i < text.length && !/\s/.test(text[i]) && text[i] !== '#') value += text[i++];
    tokens.push(value);
  }
  return tokens;
}

function splitCifField(field) {
  const dot = field.indexOf('.');
  if (!field.startsWith('_') || dot < 0) throw new Error(`unsupported CIF field ${field}`);
  return [field.slice(0, dot), field.slice(dot + 1)];
}

function inferColumnData(values) {
  if (values.length === 0) return stringData;
  if (values.some(value => value === '.' || value === '?')) return stringData;
  if (!values.every(value => /^[-+]?(?:\d+|\d*\.\d+)(?:[eE][-+]?\d+)?$/.test(value))) {
    return stringData;
  }
  return values.some(value => /[.eE]/.test(value)) ? floatData : intData;
}

function cifCategoriesFromText(text) {
  const tokens = tokenizeCif(text);
  const categories = [];
  let i = 0;
  while (i < tokens.length) {
    const token = tokens[i];
    if (token.startsWith('data_')) {
      i++;
      continue;
    }
    if (token === 'loop_') {
      i++;
      const fields = [];
      while (i < tokens.length && tokens[i].startsWith('_')) fields.push(tokens[i++]);
      if (fields.length === 0) throw new Error('loop_ without fields');
      const rows = [];
      while (
        i < tokens.length &&
        tokens[i] !== 'loop_' &&
        !tokens[i].startsWith('data_') &&
        !tokens[i].startsWith('_')
      ) {
        const row = tokens.slice(i, i + fields.length);
        if (row.length !== fields.length) throw new Error(`truncated row for ${fields[0]}`);
        rows.push(row);
        i += fields.length;
      }
      const [name] = splitCifField(fields[0]);
      const columnNames = fields.map(field => {
        const [categoryName, columnName] = splitCifField(field);
        if (categoryName !== name) {
          throw new Error(`mixed category loop ${fields[0]} and ${field}`);
        }
        return columnName;
      });
      const typed = {};
      columnNames.forEach((columnName, index) => {
        typed[columnName] = inferColumnData(rows.map(row => row[index]));
      });
      categories.push(category(name, columnNames, rows, typed));
      continue;
    }
    if (token.startsWith('_')) {
      const [name, column] = splitCifField(token);
      const value = tokens[i + 1];
      if (value === undefined) throw new Error(`missing value for ${token}`);
      const existing = categories.find(category => category.name === name && category.rowCount === 1);
      if (existing) {
        existing.columns.push({
          name: column,
          data: inferColumnData([value])([value]),
        });
      } else {
        categories.push(category(name, [column], [[value]], { [column]: inferColumnData([value]) }));
      }
      i += 2;
      continue;
    }
    throw new Error(`unexpected CIF token ${token}`);
  }
  return categories;
}

function writeBcifFromCif(filename) {
  const cifPath = join(here, '..', 'cif', filename);
  const text = readFileSync(cifPath, 'utf8');
  const header = filename.replace(/\.cif$/, '').replaceAll('-', '_');
  writeBcif(filename.replace(/\.cif$/, '.bcif'), {
    version: '0.3.0',
    encoder: 'molfig generic CIF fixture generator',
    dataBlocks: [{
      header,
      categories: cifCategoriesFromText(text),
    }],
  });
}

writeBcif('assembly-altloc-helix.bcif', file);
writeBcif('ihm-only.bcif', ihmFile);

for (const filename of readdirSync(join(here, '..', 'cif')).filter(name => name.endsWith('.cif')).sort()) {
  if (filename === 'assembly-altloc-helix.cif') continue;
  writeBcifFromCif(filename);
}
