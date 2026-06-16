// Compile-time smoke test for chem component dictionary metadata.

#import "../../../package/lib.typ" as molfig

#let ligand-cif = "data_demo\nloop_\n_chem_comp_bond.comp_id\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n_chem_comp_bond.pdbx_stereo_config\n_chem_comp_bond.pdbx_ordinal\nLIG C1 C2 DOUB N N 7\nLIG C2 C3 delo N E 8\n#\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nHETATM 1 C C1 LIG A 1 0.0 0.0 0.0\nHETATM 2 C C2 LIG A 1 20.0 0.0 0.0\nHETATM 3 C C3 LIG A 1 40.0 0.0 0.0\n#\n"

#let info = molfig.info(
  ligand-cif,
  format: "cif",
  infer-bonds: false,
  assembly: "asymmetric-unit",
)

#assert.eq(info.atom_count, 3)
#assert.eq(info.bond_count, 2)
#assert.eq(info.chem_comp_bond_count, 2)
#assert.eq(info.bond_metadata.chem_comp, 2)
#assert.eq(info.bond_metadata.aromatic, 1)
#assert.eq(info.bond_metadata.resonance, 1)
