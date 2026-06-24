use crate::chemistry::{infer_bonds, infer_element_from_name, normalize_element};
use crate::model::{
    entity_type_from_component, Assembly, AssemblyGenerator, Atom, AtomSiteAnisotrop,
    AtomSiteColumnPresence, Bond, BondFlags, BondMetadata, BondSource, ChemicalComponent,
    ChemicalComponentAngle, ChemicalComponentAtom, ChemicalComponentBond, CoarseGaussian,
    CoarseSphere, Entity, EntityIndexMap, EntityPoly, EntityPolySeq, Entry, Experiment,
    GlobalModelTransform, IhmCrossLinkRestraint, IhmModelGroup, IhmModelGroupLink, IhmModelList,
    IndexPairBonds, Molecule, PartialChargeData, PdbxBranchScheme, PdbxEntityBranch,
    PdbxEntityBranchLink, PdbxMolecule, PdbxNonpolyScheme, PdbxPolySeqScheme,
    QualityAssessmentData, SecondaryRange, SourceCategory, SourceData, StructAsym,
    StructConnMetadata, Transform, Vec3,
};
use crate::options::{InputFormat, MeshOptions};

#[cfg(test)]
thread_local! {
    static SOURCE_PARSE_COUNT: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

mod pdb;

pub(crate) use pdb::parse_pdb;

mod cif_table;

use cif_table::{
    cif_tables, cif_tokens, clean_nonempty_at, is_present_cif_value, normalized_cif_column_pair,
    CifTable,
};

pub fn parse_molecule(data: &[u8], format: InputFormat) -> Result<Molecule, String> {
    parse_molecule_with_options(
        data,
        &MeshOptions {
            format,
            ..MeshOptions::default()
        },
    )
}

pub(crate) struct ParsedMolecule {
    pub(crate) molecule: Molecule,
    pub(crate) available_alt_locs: Vec<String>,
}

pub(crate) fn parse_molecule_with_options(
    data: &[u8],
    options: &MeshOptions,
) -> Result<Molecule, String> {
    let molecule = parse_molecule_source(data, options)?;
    Ok(apply_molecule_options(molecule, options))
}

pub(crate) fn parse_molecule_with_options_and_metadata(
    data: &[u8],
    options: &MeshOptions,
) -> Result<ParsedMolecule, String> {
    let molecule = parse_molecule_source(data, options)?;
    let available_alt_locs = unique_alt_locs(&molecule.atoms);
    Ok(ParsedMolecule {
        molecule: apply_molecule_options(molecule, options),
        available_alt_locs,
    })
}

fn parse_molecule_source(data: &[u8], options: &MeshOptions) -> Result<Molecule, String> {
    #[cfg(test)]
    SOURCE_PARSE_COUNT.with(|count| count.set(count.get() + 1));

    let format = options.format;
    let detected = match format {
        InputFormat::Auto if looks_like_binary_cif(data) => InputFormat::BinaryCif,
        InputFormat::Auto => {
            let text = std::str::from_utf8(data)
                .map_err(|_| "input must be UTF-8 PDB/mmCIF data or BinaryCIF".to_string())?;
            if looks_like_pdb(text) {
                InputFormat::Pdb
            } else {
                InputFormat::Cif
            }
        }
        f => f,
    };
    match detected {
        InputFormat::Auto => unreachable!(),
        InputFormat::Pdb => {
            let text = std::str::from_utf8(data)
                .map_err(|_| "PDB input must be UTF-8 text".to_string())?;
            parse_pdb(text)
        }
        InputFormat::Cif => {
            let text = std::str::from_utf8(data)
                .map_err(|_| "mmCIF input must be UTF-8 text".to_string())?;
            parse_cif(text)
        }
        InputFormat::BinaryCif => binary::parse_binary_cif(data, options),
    }
}

#[cfg(test)]
pub(crate) fn reset_source_parse_count_for_test() {
    SOURCE_PARSE_COUNT.with(|count| count.set(0));
}

#[cfg(test)]
pub(crate) fn source_parse_count_for_test() -> usize {
    SOURCE_PARSE_COUNT.with(std::cell::Cell::get)
}

fn apply_molecule_options(mut molecule: Molecule, options: &MeshOptions) -> Molecule {
    molecule = select_alt_loc(molecule, &options.alt_loc);
    if let Some(id) = &options.assembly {
        molecule = apply_assembly(molecule, id);
    }
    if options.infer_bonds && molecule.bonds.is_empty() {
        let bonds = infer_bonds(&molecule.atoms);
        let bond_metadata = bonds
            .iter()
            .map(|bond| {
                molecule
                    .atoms
                    .get(bond.a)
                    .zip(molecule.atoms.get(bond.b))
                    .map(|(a, b)| BondMetadata::computed_for_atoms(a, b))
                    .unwrap_or_else(BondMetadata::computed)
            })
            .collect();
        molecule.bonds = bonds;
        molecule.bond_metadata = bond_metadata;
    }
    molecule.refresh_topology_metadata();
    molecule
}

fn looks_like_binary_cif(data: &[u8]) -> bool {
    data.first()
        .is_some_and(|b| matches!(*b, 0x83 | 0xde | 0xdf))
        && std::str::from_utf8(data).is_err()
}

fn looks_like_pdb(text: &str) -> bool {
    text.lines()
        .take(64)
        .any(|line| line.starts_with("ATOM  ") || line.starts_with("HETATM"))
}

fn parse_cif(text: &str) -> Result<Molecule, String> {
    let tokens = cif_tokens(text);
    let tables = cif_tables(&tokens)?;
    if tables.is_empty() && tokens.iter().any(|t| t == "atom_site") {
        return Err("mmCIF atom_site data could not be read".to_string());
    }
    let name = cif_data_block_name(&tokens);
    parse_cif_tables(&tables, SourceData::mmcif(name, source_categories(&tables)))
}

fn parse_cif_tables(tables: &[CifTable], source_data: SourceData) -> Result<Molecule, String> {
    let mut atoms = Vec::new();
    let mut atom_site_columns = AtomSiteColumnPresence {
        occupancy_defined: true,
        b_iso_defined: false,
        xyz_defined: false,
    };
    let mut assemblies = Vec::new();
    let entries = parse_cif_entries(tables);
    let experiments = parse_cif_experiments(tables);
    let mut entities = parse_cif_entities(tables);
    let entity_polymers = parse_cif_entity_polymers(tables);
    let entity_poly_seq = parse_cif_entity_poly_seq(tables);
    let pdbx_entity_branch = parse_cif_pdbx_entity_branch(tables);
    let pdbx_entity_branch_links = parse_cif_pdbx_entity_branch_links(tables);
    let pdbx_branch_scheme = parse_cif_pdbx_branch_scheme(tables);
    let pdbx_nonpoly_scheme = parse_cif_pdbx_nonpoly_scheme(tables);
    let pdbx_poly_seq_scheme = parse_cif_pdbx_poly_seq_scheme(tables);
    let ihm_model_list = parse_cif_ihm_model_list(tables);
    let ihm_model_groups = parse_cif_ihm_model_groups(tables);
    let ihm_model_group_links = parse_cif_ihm_model_group_links(tables);
    let ihm_cross_link_restraints = parse_cif_ihm_cross_link_restraints(tables);
    let struct_asym = parse_cif_struct_asym(tables);
    let pdbx_molecule = parse_cif_pdbx_molecule(tables);

    for table in tables {
        if table.name != "atom_site" {
            continue;
        }
        let idx = |name: &str| table.header_index(name);
        let ix = idx("_atom_site.Cartn_x").ok_or("mmCIF atom_site loop is missing Cartn_x")?;
        let iy = idx("_atom_site.Cartn_y").ok_or("mmCIF atom_site loop is missing Cartn_y")?;
        let iz = idx("_atom_site.Cartn_z").ok_or("mmCIF atom_site loop is missing Cartn_z")?;
        let iel = idx("_atom_site.type_symbol");
        let ilabel_name = idx("_atom_site.label_atom_id");
        let iauth_name = idx("_atom_site.auth_atom_id");
        let (ilabel_name, iauth_name) = normalized_cif_column_pair(table, ilabel_name, iauth_name);
        let ilabel_res = idx("_atom_site.label_comp_id");
        let iauth_res = idx("_atom_site.auth_comp_id");
        let (ilabel_res, iauth_res) = normalized_cif_column_pair(table, ilabel_res, iauth_res);
        let ilabel_chain = idx("_atom_site.label_asym_id");
        let iauth_chain = idx("_atom_site.auth_asym_id");
        let (ilabel_chain, iauth_chain) =
            normalized_cif_column_pair(table, ilabel_chain, iauth_chain);
        let ilabel_seq = idx("_atom_site.label_seq_id");
        let iauth_seq = idx("_atom_site.auth_seq_id");
        let (ilabel_seq, iauth_seq) = normalized_cif_column_pair(table, ilabel_seq, iauth_seq);
        let iicode = idx("_atom_site.pdbx_PDB_ins_code");
        let ialt = idx("_atom_site.label_alt_id");
        let iocc = idx("_atom_site.occupancy");
        let ib = idx("_atom_site.B_iso_or_equiv");
        let icharge = idx("_atom_site.pdbx_formal_charge");
        let igroup = idx("_atom_site.group_PDB");
        let iid = idx("_atom_site.id");
        let imodel =
            idx("_atom_site.ihm_model_id").or_else(|| idx("_atom_site.pdbx_PDB_model_num"));
        let ientity = idx("_atom_site.label_entity_id");
        atom_site_columns = AtomSiteColumnPresence {
            occupancy_defined: true,
            b_iso_defined: ib.is_some(),
            xyz_defined: true,
        };

        for row_index in table.row_indices() {
            let x = table
                .float_at(row_index, ix)
                .ok_or_else(|| format!("invalid Cartn_x: {}", table.raw_at(row_index, ix)))?;
            let y = table
                .float_at(row_index, iy)
                .ok_or_else(|| format!("invalid Cartn_y: {}", table.raw_at(row_index, iy)))?;
            let z = table
                .float_at(row_index, iz)
                .ok_or_else(|| format!("invalid Cartn_z: {}", table.raw_at(row_index, iz)))?;
            let name = ilabel_name
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_else(|| "?".to_string());
            let auth_name = iauth_name
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_else(|| name.clone());
            let raw_type_symbol = iel.and_then(|n| clean_nonempty_at(table, row_index, n));
            let type_symbol = raw_type_symbol
                .as_deref()
                .map(normalize_type_symbol_molstar)
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| normalize_type_symbol_molstar(&infer_element_from_name(&name)));
            let element = raw_type_symbol
                .map(normalize_element)
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| normalize_element(infer_element_from_name(&name)));
            let id = iid
                .and_then(|n| table.usize_at(row_index, n))
                .unwrap_or(atoms.len() + 1);
            let group = igroup
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_else(|| "ATOM".to_string());
            let occupancy = iocc
                .and_then(|n| table.float_at(row_index, n))
                .unwrap_or(1.0);
            let b_iso = ib.and_then(|n| table.float_at(row_index, n)).unwrap_or(0.0);
            let formal_charge = icharge
                .and_then(|n| table.i32_at(row_index, n))
                .unwrap_or(0);
            let chain = ilabel_chain
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_default();
            let auth_chain = iauth_chain
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_else(|| chain.clone());
            let residue = ilabel_res
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_default();
            let auth_residue = iauth_res
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_else(|| residue.clone());
            let residue_seq = ilabel_seq
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_default();
            let auth_residue_seq = iauth_seq
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_else(|| residue_seq.clone());
            let entity_id = ientity
                .and_then(|n| clean_nonempty_at(table, row_index, n))
                .or_else(|| struct_asym_entity_id(&struct_asym, &chain))
                .or_else(|| (entities.len() == 1).then(|| entities[0].id.clone()))
                .or_else(|| (!chain.is_empty()).then(|| chain.clone()))
                .unwrap_or_default();
            atoms.push(Atom {
                id,
                source_index: atoms.len(),
                model_num: imodel.and_then(|n| table.i32_at(row_index, n)).unwrap_or(1),
                name,
                type_symbol,
                auth_name,
                element,
                chain,
                auth_chain,
                entity_id,
                residue,
                auth_residue,
                group_pdb: group.clone(),
                residue_seq,
                auth_residue_seq,
                insertion_code: iicode
                    .map(|n| table.clean_at(row_index, n))
                    .unwrap_or_default(),
                alt_id: ialt
                    .map(|n| table.clean_at(row_index, n))
                    .unwrap_or_default(),
                occupancy,
                b_iso,
                formal_charge,
                position: Vec3 { x, y, z },
                het: group.eq_ignore_ascii_case("HETATM"),
                operator_name: String::new(),
            });
        }
    }

    assemblies.extend(parse_cif_assemblies(tables));
    let chemical_components = parse_cif_chemical_components(tables);
    let chemical_component_atoms = parse_cif_chemical_component_atoms(tables);
    let chemical_component_bonds = parse_cif_chemical_component_bonds(tables);
    let chemical_component_angles = parse_cif_chemical_component_angles(tables);
    let atom_site_anisotrop = parse_cif_atom_site_anisotrop(tables);
    let coarse_spheres = parse_ihm_sphere_obj_site(tables);
    let coarse_gaussians = parse_ihm_gaussian_obj_site(tables);
    entities = complete_cif_entities(entities, &atoms, &coarse_spheres, &coarse_gaussians)?;
    let entity_index = EntityIndexMap::from_mmcif(
        &entities,
        &entity_polymers,
        &pdbx_entity_branch,
        &atoms,
        &chemical_components,
        &struct_asym,
        &pdbx_molecule,
    );

    if atoms.is_empty() && coarse_spheres.is_empty() && coarse_gaussians.is_empty() {
        return Err("no _atom_site loop found in mmCIF input".to_string());
    }
    let (mut bonds, mut bond_metadata) = parse_cif_struct_conn_bonds(tables, &atoms);
    let (index_pair_raw_bonds, index_pair_metadata, index_pair_edges) =
        parse_molstar_bond_site_bonds(tables, &atoms);
    append_unique_bonds(
        &mut bonds,
        &mut bond_metadata,
        index_pair_raw_bonds.clone(),
        index_pair_metadata.clone(),
    );
    let index_pair_source_bonds = index_pair_raw_bonds
        .iter()
        .map(|raw_bond| {
            bonds
                .iter()
                .position(|bond| bond.a == raw_bond.a && bond.b == raw_bond.b)
                .unwrap_or(0)
        })
        .collect::<Vec<_>>();
    let (index_a, index_b): (Vec<_>, Vec<_>) = index_pair_edges.into_iter().unzip();
    let index_pair_bonds = IndexPairBonds::from_pairs(
        &index_a,
        &index_b,
        &index_pair_source_bonds,
        &index_pair_metadata,
        atoms.len(),
        f32::INFINITY,
        true,
    );
    let global_model_transform = parse_global_model_transform(tables);
    let quality_assessment = parse_cif_quality_assessment(tables, &atoms);
    let partial_charges = parse_cif_partial_charges(tables, &atoms);

    let (helices, struct_conf_sheets) = secondary_ranges_from_struct_conf(tables);
    let mut sheets = struct_conf_sheets;
    sheets.extend(secondary_ranges_from_tables(tables, "struct_sheet_range"));

    Ok(Molecule {
        source_data,
        atom_site_columns,
        global_model_transform,
        entries,
        experiments,
        bonds,
        bond_metadata,
        index_pair_bonds,
        atoms,
        atom_site_anisotrop,
        coarse_spheres,
        coarse_gaussians,
        assemblies,
        selected_assembly: None,
        helices,
        sheets,
        entities,
        entity_index,
        entity_polymers,
        entity_poly_seq,
        pdbx_entity_branch,
        pdbx_entity_branch_links,
        pdbx_branch_scheme,
        pdbx_nonpoly_scheme,
        pdbx_poly_seq_scheme,
        ihm_model_list,
        ihm_model_groups,
        ihm_model_group_links,
        ihm_cross_link_restraints,
        struct_asym,
        pdbx_molecule,
        chemical_components,
        chemical_component_atoms,
        chemical_component_bonds,
        chemical_component_angles,
        quality_assessment,
        partial_charges,
        rings: Vec::new(),
        resonance: Default::default(),
        derived_aromatic_bonds: Default::default(),
        derived_resonance_bonds: Default::default(),
    })
}

fn parse_cif_quality_assessment(tables: &[CifTable], atoms: &[Atom]) -> QualityAssessmentData {
    let Some(metric) = tables.iter().find(|table| table.name == "ma_qa_metric") else {
        return QualityAssessmentData::default();
    };
    let Some(local) = tables
        .iter()
        .find(|table| table.name == "ma_qa_metric_local")
    else {
        return QualityAssessmentData::default();
    };
    let (Some(metric_id), Some(metric_mode), Some(metric_name), Some(_ordinal_id)) = (
        metric.header_index("_ma_qa_metric.id"),
        metric.header_index("_ma_qa_metric.mode"),
        metric.header_index("_ma_qa_metric.name"),
        local.header_index("_ma_qa_metric_local.ordinal_id"),
    ) else {
        return QualityAssessmentData::default();
    };

    let mut plddt_metric = None;
    let mut qmean_metric = None;
    for row in metric.row_indices() {
        if !metric
            .clean_at(row, metric_mode)
            .eq_ignore_ascii_case("local")
        {
            continue;
        }
        let name = metric.clean_at(row, metric_name).to_ascii_lowercase();
        let id = metric.i32_at(row, metric_id);
        if plddt_metric.is_none() && name.contains("plddt") {
            plddt_metric = id;
        }
        if qmean_metric.is_none() && name.contains("qmean") {
            qmean_metric = id;
        }
    }

    let mut data = QualityAssessmentData {
        has_plddt_metric: plddt_metric.is_some(),
        has_qmean_metric: qmean_metric.is_some(),
        plddt: vec![None; atoms.len()],
        qmean: vec![None; atoms.len()],
    };
    let (Some(model_id), Some(asym_id), Some(seq_id), Some(local_metric_id), Some(metric_value)) = (
        local.header_index("_ma_qa_metric_local.model_id"),
        local.header_index("_ma_qa_metric_local.label_asym_id"),
        local.header_index("_ma_qa_metric_local.label_seq_id"),
        local.header_index("_ma_qa_metric_local.metric_id"),
        local.header_index("_ma_qa_metric_local.metric_value"),
    ) else {
        return data;
    };

    for row in local.row_indices() {
        let Some(metric_id) = local.i32_at(row, local_metric_id) else {
            continue;
        };
        let target = if Some(metric_id) == plddt_metric {
            &mut data.plddt
        } else if Some(metric_id) == qmean_metric {
            &mut data.qmean
        } else {
            continue;
        };
        let Some(value) = local.float_at(row, metric_value) else {
            continue;
        };
        let model = local.i32_at(row, model_id).unwrap_or(1);
        let chain = local.clean_at(row, asym_id);
        let sequence = local.clean_at(row, seq_id);
        for (atom_index, atom) in atoms.iter().enumerate() {
            if atom.model_num == model && atom.chain == chain && atom.residue_seq == sequence {
                target[atom_index] = Some(value);
            }
        }
    }
    data
}

fn parse_cif_partial_charges(tables: &[CifTable], atoms: &[Atom]) -> PartialChargeData {
    let has_atom_site = tables.iter().any(|table| table.name == "atom_site");
    let meta = tables
        .iter()
        .find(|table| table.name == "sb_ncbr_partial_atomic_charges_meta");
    let charges = tables
        .iter()
        .find(|table| table.name == "sb_ncbr_partial_atomic_charges");
    let is_applicable = has_atom_site && meta.is_some() && charges.is_some();
    let mut data = PartialChargeData {
        is_applicable,
        atom: vec![None; atoms.len()],
        residue: vec![None; atoms.len()],
        ..PartialChargeData::default()
    };
    let Some(charges) = charges else {
        return data;
    };
    let (Some(type_id), Some(atom_id), Some(charge)) = (
        charges.header_index("_sb_ncbr_partial_atomic_charges.type_id"),
        charges.header_index("_sb_ncbr_partial_atomic_charges.atom_id"),
        charges.header_index("_sb_ncbr_partial_atomic_charges.charge"),
    ) else {
        return data;
    };

    let mut by_atom_id = std::collections::BTreeMap::<usize, f32>::new();
    for row in charges.row_indices() {
        if charges.i32_at(row, type_id) != Some(1) {
            continue;
        }
        let (Some(id), Some(value)) = (
            charges.usize_at(row, atom_id),
            charges.float_at(row, charge),
        ) else {
            continue;
        };
        by_atom_id.insert(id, value);
        data.max_absolute_atom_charge = data.max_absolute_atom_charge.max(value.abs());
    }
    if by_atom_id.is_empty() {
        return data;
    }

    let mut residue_sums = std::collections::BTreeMap::<(i32, String, String, String), f32>::new();
    for atom in atoms {
        let key = (
            atom.model_num,
            atom.chain.clone(),
            atom.residue_seq.clone(),
            atom.insertion_code.clone(),
        );
        *residue_sums.entry(key).or_default() += by_atom_id.get(&atom.id).copied().unwrap_or(0.0);
    }
    for (atom_index, atom) in atoms.iter().enumerate() {
        data.atom[atom_index] = by_atom_id.get(&atom.id).copied();
        let key = (
            atom.model_num,
            atom.chain.clone(),
            atom.residue_seq.clone(),
            atom.insertion_code.clone(),
        );
        data.residue[atom_index] = residue_sums.get(&key).copied();
    }
    data.max_absolute_residue_charge = residue_sums
        .values()
        .fold(0.0f32, |max, value| max.max(value.abs()));
    data
}

fn cif_data_block_name(tokens: &[String]) -> String {
    tokens
        .iter()
        .find_map(|token| token.strip_prefix("data_"))
        .unwrap_or("")
        .to_string()
}

pub(in crate::parser) fn source_categories(tables: &[CifTable]) -> Vec<SourceCategory> {
    tables
        .iter()
        .map(|table| SourceCategory {
            name: table.name.clone(),
            row_count: table.row_count(),
            column_count: table.headers.len(),
        })
        .collect()
}

fn complete_cif_entities(
    mut entities: Vec<Entity>,
    atoms: &[Atom],
    coarse_spheres: &[CoarseSphere],
    coarse_gaussians: &[CoarseGaussian],
) -> Result<Vec<Entity>, String> {
    let mut referenced = Vec::<String>::new();
    for id in atoms
        .iter()
        .map(|atom| atom.entity_id.as_str())
        .chain(
            coarse_spheres
                .iter()
                .map(|sphere| sphere.entity_id.as_str()),
        )
        .chain(
            coarse_gaussians
                .iter()
                .map(|gaussian| gaussian.entity_id.as_str()),
        )
    {
        if !id.is_empty() && !referenced.iter().any(|existing| existing == id) {
            referenced.push(id.to_string());
        }
    }

    if entities.is_empty() {
        for id in &referenced {
            let type_name = atoms
                .iter()
                .find(|atom| atom.entity_id == *id)
                .map(|atom| entity_type_from_component(&atom.residue))
                .unwrap_or("polymer")
                .to_string();
            entities.push(Entity {
                id: id.clone(),
                type_name,
                description: String::new(),
            });
        }
        return Ok(entities);
    }

    for id in &referenced {
        if !entities.iter().any(|entity| entity.id == *id) {
            return Err(format!("missing _entity row for referenced entity id {id}"));
        }
    }
    Ok(entities)
}

fn parse_cif_chemical_components(tables: &[CifTable]) -> Vec<ChemicalComponent> {
    let Some(table) = tables.iter().find(|t| t.name == "chem_comp") else {
        return Vec::new();
    };
    let Some(id_col) = table.header_index("_chem_comp.id") else {
        return Vec::new();
    };
    let name_col = table.header_index("_chem_comp.name");
    let type_col = table.header_index("_chem_comp.type");
    let formula_col = table.header_index("_chem_comp.formula");
    let formula_weight_col = table.header_index("_chem_comp.formula_weight");
    let one_letter_code_col = table.header_index("_chem_comp.one_letter_code");
    let three_letter_code_col = table.header_index("_chem_comp.three_letter_code");
    let mon_nstd_flag_col = table.header_index("_chem_comp.mon_nstd_flag");
    let synonyms_col = table.header_index("_chem_comp.pdbx_synonyms");
    let formal_charge_col = table.header_index("_chem_comp.pdbx_formal_charge");
    let initial_date_col = table.header_index("_chem_comp.pdbx_initial_date");
    let modified_date_col = table.header_index("_chem_comp.pdbx_modified_date");
    let ambiguous_flag_col = table.header_index("_chem_comp.pdbx_ambiguous_flag");
    let release_status_col = table.header_index("_chem_comp.pdbx_release_status");
    let mut components = Vec::new();
    for row_index in table.row_indices() {
        let id = table.clean_at(row_index, id_col);
        if id.is_empty() {
            continue;
        }
        let type_name = type_col
            .map(|col| table.clean_at(row_index, col).to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .unwrap_or_else(|| "other".to_string());
        let component = ChemicalComponent {
            id,
            name: name_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            type_name,
            formula: formula_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            formula_weight: formula_weight_col.and_then(|col| table.float_at(row_index, col)),
            one_letter_code: one_letter_code_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            three_letter_code: three_letter_code_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            mon_nstd_flag: mon_nstd_flag_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            pdbx_synonyms: synonyms_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            pdbx_formal_charge: formal_charge_col.and_then(|col| table.i32_at(row_index, col)),
            pdbx_initial_date: initial_date_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            pdbx_modified_date: modified_date_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            pdbx_ambiguous_flag: ambiguous_flag_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
            pdbx_release_status: release_status_col
                .map(|col| table.clean_at(row_index, col))
                .unwrap_or_default(),
        };
        if let Some(existing) = components
            .iter_mut()
            .find(|existing: &&mut ChemicalComponent| existing.id == component.id)
        {
            *existing = component;
        } else {
            components.push(component);
        }
    }
    components
}

fn parse_cif_chemical_component_atoms(tables: &[CifTable]) -> Vec<ChemicalComponentAtom> {
    let Some(table) = tables.iter().find(|t| t.name == "chem_comp_atom") else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let (Some(comp_id), Some(atom_id)) = (
        idx("_chem_comp_atom.comp_id"),
        idx("_chem_comp_atom.atom_id"),
    ) else {
        return Vec::new();
    };
    let alt_atom_id = idx("_chem_comp_atom.alt_atom_id");
    let type_symbol = idx("_chem_comp_atom.type_symbol");
    let charge = idx("_chem_comp_atom.charge");
    let aromatic = idx("_chem_comp_atom.pdbx_aromatic_flag");
    let leaving_atom = idx("_chem_comp_atom.pdbx_leaving_atom_flag");
    let stereo_config = idx("_chem_comp_atom.pdbx_stereo_config");
    let model_x = idx("_chem_comp_atom.model_Cartn_x");
    let model_y = idx("_chem_comp_atom.model_Cartn_y");
    let model_z = idx("_chem_comp_atom.model_Cartn_z");
    let ideal_x = idx("_chem_comp_atom.pdbx_model_Cartn_x_ideal");
    let ideal_y = idx("_chem_comp_atom.pdbx_model_Cartn_y_ideal");
    let ideal_z = idx("_chem_comp_atom.pdbx_model_Cartn_z_ideal");

    table
        .row_indices()
        .filter_map(|row| {
            let comp_id = table.clean_at(row, comp_id);
            let atom_id = table.clean_at(row, atom_id);
            if comp_id.is_empty() || atom_id.is_empty() {
                return None;
            }
            Some(ChemicalComponentAtom {
                comp_id,
                atom_id,
                alt_atom_id: alt_atom_id
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                type_symbol: type_symbol
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                charge: charge.and_then(|col| table.i32_at(row, col)).unwrap_or(0),
                aromatic: flag_yes(aromatic.map(|col| table.clean_at(row, col)).as_deref()),
                leaving_atom: flag_yes(leaving_atom.map(|col| table.clean_at(row, col)).as_deref()),
                stereo_config: stereo_config
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                model_cartn: optional_vec3(table, row, model_x, model_y, model_z),
                ideal_cartn: optional_vec3(table, row, ideal_x, ideal_y, ideal_z),
            })
        })
        .collect()
}

fn parse_cif_chemical_component_angles(tables: &[CifTable]) -> Vec<ChemicalComponentAngle> {
    let Some(table) = tables.iter().find(|t| t.name == "chem_comp_angle") else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let (Some(comp_id), Some(atom_id_1), Some(atom_id_2), Some(atom_id_3)) = (
        idx("_chem_comp_angle.comp_id"),
        idx("_chem_comp_angle.atom_id_1"),
        idx("_chem_comp_angle.atom_id_2"),
        idx("_chem_comp_angle.atom_id_3"),
    ) else {
        return Vec::new();
    };
    let value_angle = idx("_chem_comp_angle.value_angle");
    let value_angle_esd = idx("_chem_comp_angle.value_angle_esd");

    table
        .row_indices()
        .filter_map(|row| {
            let comp_id = table.clean_at(row, comp_id);
            let atom_id_1 = table.clean_at(row, atom_id_1);
            let atom_id_2 = table.clean_at(row, atom_id_2);
            let atom_id_3 = table.clean_at(row, atom_id_3);
            if comp_id.is_empty()
                || atom_id_1.is_empty()
                || atom_id_2.is_empty()
                || atom_id_3.is_empty()
            {
                return None;
            }
            Some(ChemicalComponentAngle {
                comp_id,
                atom_id_1,
                atom_id_2,
                atom_id_3,
                value_angle: value_angle.and_then(|col| table.float_at(row, col)),
                value_angle_esd: value_angle_esd.and_then(|col| table.float_at(row, col)),
            })
        })
        .collect()
}

fn parse_cif_entries(tables: &[CifTable]) -> Vec<Entry> {
    let Some(table) = tables.iter().find(|t| t.name == "entry") else {
        return Vec::new();
    };
    let Some(id_col) = table.header_index("_entry.id") else {
        return Vec::new();
    };
    table
        .row_indices()
        .filter_map(|row| {
            let id = table.clean_at(row, id_col);
            (!id.is_empty()).then_some(Entry { id })
        })
        .collect()
}

fn flag_yes(value: Option<&str>) -> bool {
    value.is_some_and(|value| {
        value.eq_ignore_ascii_case("Y")
            || value.eq_ignore_ascii_case("YES")
            || value.eq_ignore_ascii_case("1")
            || value.eq_ignore_ascii_case("TRUE")
    })
}

fn optional_vec3(
    table: &CifTable,
    row: usize,
    x: Option<usize>,
    y: Option<usize>,
    z: Option<usize>,
) -> Option<Vec3> {
    Some(Vec3 {
        x: table.float_at(row, x?)?,
        y: table.float_at(row, y?)?,
        z: table.float_at(row, z?)?,
    })
}

fn parse_cif_experiments(tables: &[CifTable]) -> Vec<Experiment> {
    let Some(table) = tables.iter().find(|t| t.name == "exptl") else {
        return Vec::new();
    };
    let Some(method_col) = table.header_index("_exptl.method") else {
        return Vec::new();
    };
    table
        .row_indices()
        .filter_map(|row| {
            let method = table.clean_at(row, method_col);
            (!method.is_empty()).then_some(Experiment { method })
        })
        .collect()
}

fn parse_cif_entities(tables: &[CifTable]) -> Vec<Entity> {
    let Some(table) = tables.iter().find(|t| t.name == "entity") else {
        return Vec::new();
    };
    let Some(id_col) = table.header_index("_entity.id") else {
        return Vec::new();
    };
    let type_col = table.header_index("_entity.type");
    let desc_col = table.header_index("_entity.pdbx_description");
    table
        .row_indices()
        .filter_map(|row| {
            let id = table.clean_at(row, id_col);
            (!id.is_empty()).then(|| Entity {
                id,
                type_name: type_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                description: desc_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
            })
        })
        .collect()
}

fn parse_cif_entity_polymers(tables: &[CifTable]) -> Vec<EntityPoly> {
    let Some(table) = tables.iter().find(|t| t.name == "entity_poly") else {
        return Vec::new();
    };
    let Some(entity_id_col) = table.header_index("_entity_poly.entity_id") else {
        return Vec::new();
    };
    let type_col = table.header_index("_entity_poly.type");
    let sequence_col = table
        .header_index("_entity_poly.pdbx_seq_one_letter_code")
        .or_else(|| table.header_index("_entity_poly.pdbx_seq_one_letter_code_can"));
    let nstd_linkage_col = table.header_index("_entity_poly.nstd_linkage");
    let nstd_monomer_col = table.header_index("_entity_poly.nstd_monomer");
    table
        .row_indices()
        .filter_map(|row| {
            let entity_id = table.clean_at(row, entity_id_col);
            (!entity_id.is_empty()).then(|| EntityPoly {
                entity_id,
                polymer_type: type_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                sequence: sequence_col
                    .map(|col| table.clean_at(row, col).replace(['\n', ' '], ""))
                    .unwrap_or_default(),
                nstd_linkage: nstd_linkage_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                nstd_monomer: nstd_monomer_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
            })
        })
        .collect()
}

fn parse_cif_entity_poly_seq(tables: &[CifTable]) -> Vec<EntityPolySeq> {
    let Some(table) = tables.iter().find(|t| t.name == "entity_poly_seq") else {
        return Vec::new();
    };
    let (Some(entity_id_col), Some(num_col), Some(mon_id_col)) = (
        table.header_index("_entity_poly_seq.entity_id"),
        table.header_index("_entity_poly_seq.num"),
        table.header_index("_entity_poly_seq.mon_id"),
    ) else {
        return Vec::new();
    };
    let hetero_col = table.header_index("_entity_poly_seq.hetero");
    table
        .row_indices()
        .filter_map(|row| {
            let entity_id = table.clean_at(row, entity_id_col);
            let num = table.i32_at(row, num_col)?;
            let mon_id = table.clean_at(row, mon_id_col);
            (!entity_id.is_empty() && !mon_id.is_empty()).then(|| EntityPolySeq {
                entity_id,
                num,
                mon_id,
                hetero: hetero_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
            })
        })
        .collect()
}

fn parse_cif_pdbx_entity_branch(tables: &[CifTable]) -> Vec<PdbxEntityBranch> {
    let Some(table) = tables.iter().find(|t| t.name == "pdbx_entity_branch") else {
        return Vec::new();
    };
    let Some(entity_id_col) = table.header_index("_pdbx_entity_branch.entity_id") else {
        return Vec::new();
    };
    let type_col = table.header_index("_pdbx_entity_branch.type");
    table
        .row_indices()
        .filter_map(|row| {
            let entity_id = table.clean_at(row, entity_id_col);
            (!entity_id.is_empty()).then(|| PdbxEntityBranch {
                entity_id,
                type_name: type_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
            })
        })
        .collect()
}

fn parse_cif_pdbx_entity_branch_links(tables: &[CifTable]) -> Vec<PdbxEntityBranchLink> {
    let Some(table) = tables.iter().find(|t| t.name == "pdbx_entity_branch_link") else {
        return Vec::new();
    };
    let (Some(link_id_col), Some(entity_id_col)) = (
        table.header_index("_pdbx_entity_branch_link.link_id"),
        table.header_index("_pdbx_entity_branch_link.entity_id"),
    ) else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let details_col = idx("_pdbx_entity_branch_link.details");
    let num_1_col = idx("_pdbx_entity_branch_link.entity_branch_list_num_1");
    let num_2_col = idx("_pdbx_entity_branch_link.entity_branch_list_num_2");
    let comp_1_col = idx("_pdbx_entity_branch_link.comp_id_1");
    let comp_2_col = idx("_pdbx_entity_branch_link.comp_id_2");
    let atom_1_col = idx("_pdbx_entity_branch_link.atom_id_1");
    let leaving_1_col = idx("_pdbx_entity_branch_link.leaving_atom_id_1");
    let stereo_1_col = idx("_pdbx_entity_branch_link.atom_stereo_config_1");
    let atom_2_col = idx("_pdbx_entity_branch_link.atom_id_2");
    let leaving_2_col = idx("_pdbx_entity_branch_link.leaving_atom_id_2");
    let stereo_2_col = idx("_pdbx_entity_branch_link.atom_stereo_config_2");
    let order_col = idx("_pdbx_entity_branch_link.value_order");
    table
        .row_indices()
        .filter_map(|row| {
            let link_id = table.i32_at(row, link_id_col)?;
            let entity_id = table.clean_at(row, entity_id_col);
            (!entity_id.is_empty()).then(|| PdbxEntityBranchLink {
                link_id,
                details: details_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                entity_id,
                entity_branch_list_num_1: num_1_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
                entity_branch_list_num_2: num_2_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
                comp_id_1: comp_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                comp_id_2: comp_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                atom_id_1: atom_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                leaving_atom_id_1: leaving_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                atom_stereo_config_1: stereo_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                atom_id_2: atom_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                leaving_atom_id_2: leaving_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                atom_stereo_config_2: stereo_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                value_order: order_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
            })
        })
        .collect()
}

fn parse_cif_pdbx_branch_scheme(tables: &[CifTable]) -> Vec<PdbxBranchScheme> {
    let Some(table) = tables.iter().find(|t| t.name == "pdbx_branch_scheme") else {
        return Vec::new();
    };
    let (Some(entity_id_col), Some(asym_id_col), Some(mon_id_col), Some(num_col)) = (
        table.header_index("_pdbx_branch_scheme.entity_id"),
        table.header_index("_pdbx_branch_scheme.asym_id"),
        table.header_index("_pdbx_branch_scheme.mon_id"),
        table.header_index("_pdbx_branch_scheme.num"),
    ) else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let hetero_col = idx("_pdbx_branch_scheme.hetero");
    let pdb_asym_col = idx("_pdbx_branch_scheme.pdb_asym_id");
    let pdb_seq_col = idx("_pdbx_branch_scheme.pdb_seq_num");
    let pdb_mon_col = idx("_pdbx_branch_scheme.pdb_mon_id");
    let auth_asym_col = idx("_pdbx_branch_scheme.auth_asym_id");
    let auth_seq_col = idx("_pdbx_branch_scheme.auth_seq_num");
    let auth_mon_col = idx("_pdbx_branch_scheme.auth_mon_id");
    table
        .row_indices()
        .filter_map(|row| {
            let entity_id = table.clean_at(row, entity_id_col);
            let asym_id = table.clean_at(row, asym_id_col);
            let mon_id = table.clean_at(row, mon_id_col);
            let num = table.i32_at(row, num_col)?;
            (!entity_id.is_empty() && !asym_id.is_empty() && !mon_id.is_empty()).then(|| {
                PdbxBranchScheme {
                    entity_id,
                    hetero: hetero_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    asym_id,
                    mon_id,
                    num,
                    pdb_asym_id: pdb_asym_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_seq_num: pdb_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_mon_id: pdb_mon_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    auth_asym_id: auth_asym_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    auth_seq_num: auth_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    auth_mon_id: auth_mon_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                }
            })
        })
        .collect()
}

fn parse_cif_pdbx_nonpoly_scheme(tables: &[CifTable]) -> Vec<PdbxNonpolyScheme> {
    let Some(table) = tables.iter().find(|t| t.name == "pdbx_nonpoly_scheme") else {
        return Vec::new();
    };
    let (Some(asym_id_col), Some(entity_id_col), Some(mon_id_col)) = (
        table.header_index("_pdbx_nonpoly_scheme.asym_id"),
        table.header_index("_pdbx_nonpoly_scheme.entity_id"),
        table.header_index("_pdbx_nonpoly_scheme.mon_id"),
    ) else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let pdb_strand_col = idx("_pdbx_nonpoly_scheme.pdb_strand_id");
    let ndb_seq_col = idx("_pdbx_nonpoly_scheme.ndb_seq_num");
    let pdb_seq_col = idx("_pdbx_nonpoly_scheme.pdb_seq_num");
    let auth_seq_col = idx("_pdbx_nonpoly_scheme.auth_seq_num");
    let pdb_mon_col = idx("_pdbx_nonpoly_scheme.pdb_mon_id");
    let auth_mon_col = idx("_pdbx_nonpoly_scheme.auth_mon_id");
    let ins_col = idx("_pdbx_nonpoly_scheme.pdb_ins_code");
    table
        .row_indices()
        .filter_map(|row| {
            let asym_id = table.clean_at(row, asym_id_col);
            let entity_id = table.clean_at(row, entity_id_col);
            let mon_id = table.clean_at(row, mon_id_col);
            (!asym_id.is_empty() && !entity_id.is_empty() && !mon_id.is_empty()).then(|| {
                PdbxNonpolyScheme {
                    asym_id,
                    entity_id,
                    mon_id,
                    pdb_strand_id: pdb_strand_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    ndb_seq_num: ndb_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_seq_num: pdb_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    auth_seq_num: auth_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_mon_id: pdb_mon_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    auth_mon_id: auth_mon_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_ins_code: ins_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                }
            })
        })
        .collect()
}

fn parse_cif_pdbx_poly_seq_scheme(tables: &[CifTable]) -> Vec<PdbxPolySeqScheme> {
    let Some(table) = tables.iter().find(|t| t.name == "pdbx_poly_seq_scheme") else {
        return Vec::new();
    };
    let (Some(asym_id_col), Some(entity_id_col), Some(seq_id_col), Some(mon_id_col)) = (
        table.header_index("_pdbx_poly_seq_scheme.asym_id"),
        table.header_index("_pdbx_poly_seq_scheme.entity_id"),
        table.header_index("_pdbx_poly_seq_scheme.seq_id"),
        table.header_index("_pdbx_poly_seq_scheme.mon_id"),
    ) else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let ndb_seq_col = idx("_pdbx_poly_seq_scheme.ndb_seq_num");
    let pdb_seq_col = idx("_pdbx_poly_seq_scheme.pdb_seq_num");
    let auth_seq_col = idx("_pdbx_poly_seq_scheme.auth_seq_num");
    let pdb_mon_col = idx("_pdbx_poly_seq_scheme.pdb_mon_id");
    let auth_mon_col = idx("_pdbx_poly_seq_scheme.auth_mon_id");
    let pdb_strand_col = idx("_pdbx_poly_seq_scheme.pdb_strand_id");
    let ins_col = idx("_pdbx_poly_seq_scheme.pdb_ins_code");
    let hetero_col = idx("_pdbx_poly_seq_scheme.hetero");
    table
        .row_indices()
        .filter_map(|row| {
            let asym_id = table.clean_at(row, asym_id_col);
            let entity_id = table.clean_at(row, entity_id_col);
            let seq_id = table.i32_at(row, seq_id_col)?;
            let mon_id = table.clean_at(row, mon_id_col);
            (!asym_id.is_empty() && !entity_id.is_empty() && !mon_id.is_empty()).then(|| {
                PdbxPolySeqScheme {
                    asym_id,
                    entity_id,
                    seq_id,
                    mon_id,
                    ndb_seq_num: ndb_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_seq_num: pdb_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    auth_seq_num: auth_seq_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_mon_id: pdb_mon_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    auth_mon_id: auth_mon_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_strand_id: pdb_strand_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    pdb_ins_code: ins_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                    hetero: hetero_col
                        .map(|col| table.clean_at(row, col))
                        .unwrap_or_default(),
                }
            })
        })
        .collect()
}

fn parse_cif_ihm_model_list(tables: &[CifTable]) -> Vec<IhmModelList> {
    let Some(table) = tables.iter().find(|t| t.name == "ihm_model_list") else {
        return Vec::new();
    };
    let Some(model_id_col) = table.header_index("_ihm_model_list.model_id") else {
        return Vec::new();
    };
    let model_name_col = table.header_index("_ihm_model_list.model_name");
    let assembly_id_col = table.header_index("_ihm_model_list.assembly_id");
    let protocol_id_col = table.header_index("_ihm_model_list.protocol_id");
    let representation_id_col = table.header_index("_ihm_model_list.representation_id");
    table
        .row_indices()
        .filter_map(|row| {
            Some(IhmModelList {
                model_id: table.i32_at(row, model_id_col)?,
                model_name: model_name_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                assembly_id: assembly_id_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
                protocol_id: protocol_id_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
                representation_id: representation_id_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
            })
        })
        .collect()
}

fn parse_cif_ihm_model_groups(tables: &[CifTable]) -> Vec<IhmModelGroup> {
    let Some(table) = tables.iter().find(|t| t.name == "ihm_model_group") else {
        return Vec::new();
    };
    let Some(id_col) = table.header_index("_ihm_model_group.id") else {
        return Vec::new();
    };
    let name_col = table.header_index("_ihm_model_group.name");
    let details_col = table.header_index("_ihm_model_group.details");
    table
        .row_indices()
        .filter_map(|row| {
            Some(IhmModelGroup {
                id: table.i32_at(row, id_col)?,
                name: name_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                details: details_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
            })
        })
        .collect()
}

fn parse_cif_ihm_model_group_links(tables: &[CifTable]) -> Vec<IhmModelGroupLink> {
    let Some(table) = tables.iter().find(|t| t.name == "ihm_model_group_link") else {
        return Vec::new();
    };
    let (Some(model_id_col), Some(group_id_col)) = (
        table.header_index("_ihm_model_group_link.model_id"),
        table.header_index("_ihm_model_group_link.group_id"),
    ) else {
        return Vec::new();
    };
    table
        .row_indices()
        .filter_map(|row| {
            Some(IhmModelGroupLink {
                model_id: table.i32_at(row, model_id_col)?,
                group_id: table.i32_at(row, group_id_col)?,
            })
        })
        .collect()
}

fn parse_cif_ihm_cross_link_restraints(tables: &[CifTable]) -> Vec<IhmCrossLinkRestraint> {
    let Some(table) = tables.iter().find(|t| t.name == "ihm_cross_link_restraint") else {
        return Vec::new();
    };
    let Some(id_col) = table.header_index("_ihm_cross_link_restraint.id") else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let group_id_col = idx("_ihm_cross_link_restraint.group_id");
    let entity_1_col = idx("_ihm_cross_link_restraint.entity_id_1");
    let entity_2_col = idx("_ihm_cross_link_restraint.entity_id_2");
    let asym_1_col = idx("_ihm_cross_link_restraint.asym_id_1");
    let asym_2_col = idx("_ihm_cross_link_restraint.asym_id_2");
    let comp_1_col = idx("_ihm_cross_link_restraint.comp_id_1");
    let comp_2_col = idx("_ihm_cross_link_restraint.comp_id_2");
    let seq_1_col = idx("_ihm_cross_link_restraint.seq_id_1");
    let seq_2_col = idx("_ihm_cross_link_restraint.seq_id_2");
    let atom_1_col = idx("_ihm_cross_link_restraint.atom_id_1");
    let atom_2_col = idx("_ihm_cross_link_restraint.atom_id_2");
    let restraint_type_col = idx("_ihm_cross_link_restraint.restraint_type");
    let conditional_col = idx("_ihm_cross_link_restraint.conditional_crosslink_flag");
    let granularity_col = idx("_ihm_cross_link_restraint.model_granularity");
    let threshold_col = idx("_ihm_cross_link_restraint.distance_threshold");
    let psi_col = idx("_ihm_cross_link_restraint.psi");
    let sigma_1_col = idx("_ihm_cross_link_restraint.sigma_1");
    let sigma_2_col = idx("_ihm_cross_link_restraint.sigma_2");
    table
        .row_indices()
        .filter_map(|row| {
            Some(IhmCrossLinkRestraint {
                id: table.i32_at(row, id_col)?,
                group_id: group_id_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
                entity_id_1: entity_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                entity_id_2: entity_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                asym_id_1: asym_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                asym_id_2: asym_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                comp_id_1: comp_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                comp_id_2: comp_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                seq_id_1: seq_1_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
                seq_id_2: seq_2_col
                    .and_then(|col| table.i32_at(row, col))
                    .unwrap_or(0),
                atom_id_1: atom_1_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                atom_id_2: atom_2_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                restraint_type: restraint_type_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                conditional_crosslink_flag: conditional_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                model_granularity: granularity_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
                distance_threshold: threshold_col.and_then(|col| table.float_at(row, col)),
                psi: psi_col.and_then(|col| table.float_at(row, col)),
                sigma_1: sigma_1_col.and_then(|col| table.float_at(row, col)),
                sigma_2: sigma_2_col.and_then(|col| table.float_at(row, col)),
            })
        })
        .collect()
}

fn parse_cif_struct_asym(tables: &[CifTable]) -> Vec<StructAsym> {
    let Some(table) = tables.iter().find(|t| t.name == "struct_asym") else {
        return Vec::new();
    };
    let (Some(id_col), Some(entity_id_col)) = (
        table.header_index("_struct_asym.id"),
        table.header_index("_struct_asym.entity_id"),
    ) else {
        return Vec::new();
    };
    let details_col = table.header_index("_struct_asym.details");
    table
        .row_indices()
        .filter_map(|row| {
            let id = table.clean_at(row, id_col);
            let entity_id = table.clean_at(row, entity_id_col);
            (!id.is_empty() && !entity_id.is_empty()).then(|| StructAsym {
                id,
                entity_id,
                details: details_col
                    .map(|col| table.clean_at(row, col))
                    .unwrap_or_default(),
            })
        })
        .collect()
}

fn parse_cif_pdbx_molecule(tables: &[CifTable]) -> Vec<PdbxMolecule> {
    let Some(table) = tables.iter().find(|table| table.name == "pdbx_molecule") else {
        return Vec::new();
    };
    let (Some(asym_id_col), Some(prd_id_col)) = (
        table.header_index("_pdbx_molecule.asym_id"),
        table.header_index("_pdbx_molecule.prd_id"),
    ) else {
        return Vec::new();
    };
    table
        .row_indices()
        .filter_map(|row| {
            let asym_id = table.clean_at(row, asym_id_col);
            let prd_id = table.clean_at(row, prd_id_col);
            (!asym_id.is_empty() && !prd_id.is_empty()).then_some(PdbxMolecule { asym_id, prd_id })
        })
        .collect()
}

fn struct_asym_entity_id(struct_asym: &[StructAsym], asym_id: &str) -> Option<String> {
    struct_asym
        .iter()
        .find(|asym| asym.id == asym_id)
        .map(|asym| asym.entity_id.clone())
}

fn parse_cif_chemical_component_bonds(tables: &[CifTable]) -> Vec<ChemicalComponentBond> {
    let Some(table) = tables.iter().find(|t| t.name == "chem_comp_bond") else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let (Some(comp_id), Some(atom_id_1), Some(atom_id_2)) = (
        idx("_chem_comp_bond.comp_id"),
        idx("_chem_comp_bond.atom_id_1"),
        idx("_chem_comp_bond.atom_id_2"),
    ) else {
        return Vec::new();
    };
    let value_order = idx("_chem_comp_bond.value_order");
    let aromatic = idx("_chem_comp_bond.pdbx_aromatic_flag");
    let stereo_config = idx("_chem_comp_bond.pdbx_stereo_config");
    let ordinal = idx("_chem_comp_bond.pdbx_ordinal");
    table
        .row_indices()
        .filter_map(|row| {
            let comp_id = table.clean_at(row, comp_id);
            let atom_id_1 = table.clean_at(row, atom_id_1);
            let atom_id_2 = table.clean_at(row, atom_id_2);
            if comp_id.is_empty() || atom_id_1.is_empty() || atom_id_2.is_empty() {
                return None;
            }
            let value_order = value_order
                .map(|n| table.clean_at(row, n))
                .unwrap_or_else(|| "sing".to_string());
            let aromatic = aromatic
                .map(|n| table.clean_at(row, n))
                .is_some_and(|value| value.eq_ignore_ascii_case("Y"));
            let order_flags = molstar_bond_site_order_flags(&value_order);
            let flags = if aromatic || order_flags.contains(BondFlags::AROMATIC) {
                order_flags
                    .union(BondFlags::AROMATIC)
                    .union(BondFlags::RESONANCE)
            } else {
                order_flags
            };
            Some(ChemicalComponentBond {
                comp_id,
                atom_id_1,
                atom_id_2,
                order: molstar_bond_site_order(&value_order),
                flags,
                stereo_config: stereo_config
                    .map(|n| table.clean_at(row, n))
                    .unwrap_or_default(),
                ordinal: ordinal.and_then(|n| table.i32_at(row, n)),
            })
        })
        .collect()
}

fn parse_cif_atom_site_anisotrop(tables: &[CifTable]) -> Vec<AtomSiteAnisotrop> {
    let Some(table) = tables.iter().find(|t| t.name == "atom_site_anisotrop") else {
        return Vec::new();
    };
    let Some(id_col) = table.header_index("_atom_site_anisotrop.id") else {
        return Vec::new();
    };
    let indices = [
        [
            table.header_index("_atom_site_anisotrop.U[1][1]"),
            table.header_index("_atom_site_anisotrop.U[1][2]"),
            table.header_index("_atom_site_anisotrop.U[1][3]"),
        ],
        [
            table.header_index("_atom_site_anisotrop.U[2][1]"),
            table.header_index("_atom_site_anisotrop.U[2][2]"),
            table.header_index("_atom_site_anisotrop.U[2][3]"),
        ],
        [
            table.header_index("_atom_site_anisotrop.U[3][1]"),
            table.header_index("_atom_site_anisotrop.U[3][2]"),
            table.header_index("_atom_site_anisotrop.U[3][3]"),
        ],
    ];
    if indices.iter().flatten().any(Option::is_none) {
        return Vec::new();
    }
    table
        .row_indices()
        .filter_map(|row| {
            let atom_id = table.usize_at(row, id_col)?;
            let mut u = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    u[i][j] = table.float_at(row, indices[i][j]?)?;
                }
            }
            Some(AtomSiteAnisotrop { atom_id, u })
        })
        .collect()
}

fn parse_ihm_sphere_obj_site(tables: &[CifTable]) -> Vec<CoarseSphere> {
    let Some(table) = tables.iter().find(|t| t.name == "ihm_sphere_obj_site") else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let (Some(x), Some(y), Some(z)) = (
        idx("_ihm_sphere_obj_site.Cartn_x"),
        idx("_ihm_sphere_obj_site.Cartn_y"),
        idx("_ihm_sphere_obj_site.Cartn_z"),
    ) else {
        return Vec::new();
    };
    let id = idx("_ihm_sphere_obj_site.id");
    let model_num = idx("_ihm_sphere_obj_site.model_id");
    let entity_id = idx("_ihm_sphere_obj_site.entity_id");
    let asym_id = idx("_ihm_sphere_obj_site.asym_id");
    let seq_begin = idx("_ihm_sphere_obj_site.seq_id_begin");
    let seq_end = idx("_ihm_sphere_obj_site.seq_id_end");
    let radius =
        idx("_ihm_sphere_obj_site.object_radius").or_else(|| idx("_ihm_sphere_obj_site.radius"));
    let rmsf = idx("_ihm_sphere_obj_site.rmsf");
    table
        .row_indices()
        .filter_map(|row| {
            Some(CoarseSphere {
                id: id.and_then(|n| table.usize_at(row, n)).unwrap_or(row + 1),
                model_num: model_num.and_then(|n| table.i32_at(row, n)).unwrap_or(1),
                entity_id: entity_id
                    .map(|n| table.clean_at(row, n))
                    .unwrap_or_default(),
                asym_id: asym_id.map(|n| table.clean_at(row, n)).unwrap_or_default(),
                seq_id_begin: seq_begin.and_then(|n| table.i32_at(row, n)).unwrap_or(0),
                seq_id_end: seq_end.and_then(|n| table.i32_at(row, n)).unwrap_or(0),
                position: Vec3 {
                    x: table.float_at(row, x)?,
                    y: table.float_at(row, y)?,
                    z: table.float_at(row, z)?,
                },
                radius: radius
                    .and_then(|n| table.float_at(row, n))
                    .unwrap_or(1.0)
                    .max(0.01),
                rmsf: rmsf.and_then(|n| table.float_at(row, n)).unwrap_or(0.0),
            })
        })
        .collect()
}

fn parse_ihm_gaussian_obj_site(tables: &[CifTable]) -> Vec<CoarseGaussian> {
    let Some(table) = tables.iter().find(|t| t.name == "ihm_gaussian_obj_site") else {
        return Vec::new();
    };
    let idx = |name: &str| table.header_index(name);
    let x = idx("_ihm_gaussian_obj_site.mean_Cartn_x")
        .or_else(|| idx("_ihm_gaussian_obj_site.Cartn_x"));
    let y = idx("_ihm_gaussian_obj_site.mean_Cartn_y")
        .or_else(|| idx("_ihm_gaussian_obj_site.Cartn_y"));
    let z = idx("_ihm_gaussian_obj_site.mean_Cartn_z")
        .or_else(|| idx("_ihm_gaussian_obj_site.Cartn_z"));
    let (Some(x), Some(y), Some(z)) = (x, y, z) else {
        return Vec::new();
    };
    let id = idx("_ihm_gaussian_obj_site.id");
    let model_num = idx("_ihm_gaussian_obj_site.model_id");
    let entity_id = idx("_ihm_gaussian_obj_site.entity_id");
    let asym_id = idx("_ihm_gaussian_obj_site.asym_id");
    let seq_begin = idx("_ihm_gaussian_obj_site.seq_id_begin");
    let seq_end = idx("_ihm_gaussian_obj_site.seq_id_end");
    let weight = idx("_ihm_gaussian_obj_site.weight");
    let cov = [
        [
            idx("_ihm_gaussian_obj_site.covariance_matrix[1][1]"),
            idx("_ihm_gaussian_obj_site.covariance_matrix[1][2]"),
            idx("_ihm_gaussian_obj_site.covariance_matrix[1][3]"),
        ],
        [
            idx("_ihm_gaussian_obj_site.covariance_matrix[2][1]"),
            idx("_ihm_gaussian_obj_site.covariance_matrix[2][2]"),
            idx("_ihm_gaussian_obj_site.covariance_matrix[2][3]"),
        ],
        [
            idx("_ihm_gaussian_obj_site.covariance_matrix[3][1]"),
            idx("_ihm_gaussian_obj_site.covariance_matrix[3][2]"),
            idx("_ihm_gaussian_obj_site.covariance_matrix[3][3]"),
        ],
    ];
    table
        .row_indices()
        .filter_map(|row| {
            let mut covariance = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    covariance[i][j] = cov[i][j]
                        .and_then(|n| table.float_at(row, n))
                        .unwrap_or(if i == j { 1.0 } else { 0.0 });
                }
            }
            Some(CoarseGaussian {
                id: id.and_then(|n| table.usize_at(row, n)).unwrap_or(row + 1),
                model_num: model_num.and_then(|n| table.i32_at(row, n)).unwrap_or(1),
                entity_id: entity_id
                    .map(|n| table.clean_at(row, n))
                    .unwrap_or_default(),
                asym_id: asym_id.map(|n| table.clean_at(row, n)).unwrap_or_default(),
                seq_id_begin: seq_begin.and_then(|n| table.i32_at(row, n)).unwrap_or(0),
                seq_id_end: seq_end.and_then(|n| table.i32_at(row, n)).unwrap_or(0),
                position: Vec3 {
                    x: table.float_at(row, x)?,
                    y: table.float_at(row, y)?,
                    z: table.float_at(row, z)?,
                },
                weight: weight.and_then(|n| table.float_at(row, n)).unwrap_or(1.0),
                covariance,
            })
        })
        .collect()
}

fn parse_cif_struct_conn_bonds(
    tables: &[CifTable],
    atoms: &[Atom],
) -> (Vec<Bond>, Vec<BondMetadata>) {
    let Some(table) = tables.iter().find(|t| t.name == "struct_conn") else {
        return (Vec::new(), Vec::new());
    };
    let idx = |name: &str| table.header_index(name);
    let row_id = idx("_struct_conn.id");
    let conn_type = idx("_struct_conn.conn_type_id");
    let distance = idx("_struct_conn.pdbx_dist_value");
    let value_order = idx("_struct_conn.pdbx_value_order");
    let label_then_auth = |label: &str, auth: &str| {
        idx(label)
            .map(|index| CifAtomIdColumn::new(index, CifAtomIdKind::Label))
            .or_else(|| idx(auth).map(|index| CifAtomIdColumn::new(index, CifAtomIdKind::Auth)))
    };
    let auth_seq_then_label = |auth: &str, label: &str| {
        idx(auth)
            .map(|index| CifAtomIdColumn::new(index, CifAtomIdKind::Auth))
            .or_else(|| idx(label).map(|index| CifAtomIdColumn::new(index, CifAtomIdKind::Label)))
    };
    let p1_chain = label_then_auth(
        "_struct_conn.ptnr1_label_asym_id",
        "_struct_conn.ptnr1_auth_asym_id",
    );
    let p1_comp =
        idx("_struct_conn.ptnr1_label_comp_id").or_else(|| idx("_struct_conn.ptnr1_auth_comp_id"));
    let p1_seq = auth_seq_then_label(
        "_struct_conn.ptnr1_auth_seq_id",
        "_struct_conn.ptnr1_label_seq_id",
    );
    let p1_atom = label_then_auth(
        "_struct_conn.ptnr1_label_atom_id",
        "_struct_conn.ptnr1_auth_atom_id",
    );
    let p1_alt = idx("_struct_conn.pdbx_ptnr1_label_alt_id")
        .or_else(|| idx("_struct_conn.ptnr1_label_alt_id"));
    let p1_ins = idx("_struct_conn.pdbx_ptnr1_PDB_ins_code");
    let p1_symmetry = idx("_struct_conn.ptnr1_symmetry");
    let p2_chain = label_then_auth(
        "_struct_conn.ptnr2_label_asym_id",
        "_struct_conn.ptnr2_auth_asym_id",
    );
    let p2_comp =
        idx("_struct_conn.ptnr2_label_comp_id").or_else(|| idx("_struct_conn.ptnr2_auth_comp_id"));
    let p2_seq = auth_seq_then_label(
        "_struct_conn.ptnr2_auth_seq_id",
        "_struct_conn.ptnr2_label_seq_id",
    );
    let p2_atom = label_then_auth(
        "_struct_conn.ptnr2_label_atom_id",
        "_struct_conn.ptnr2_auth_atom_id",
    );
    let p2_alt = idx("_struct_conn.pdbx_ptnr2_label_alt_id")
        .or_else(|| idx("_struct_conn.ptnr2_label_alt_id"));
    let p2_ins = idx("_struct_conn.pdbx_ptnr2_PDB_ins_code");
    let p2_symmetry = idx("_struct_conn.ptnr2_symmetry");
    let (Some(p1_chain), Some(p1_seq), Some(p1_atom), Some(p2_chain), Some(p2_seq), Some(p2_atom)) =
        (p1_chain, p1_seq, p1_atom, p2_chain, p2_seq, p2_atom)
    else {
        return (Vec::new(), Vec::new());
    };

    let mut bonds = Vec::new();
    let mut metadata = Vec::new();
    for row_index in 0..table.row_count() {
        let a_candidates = find_cif_atoms(
            atoms,
            CifAtomIdValue::from_table(table, row_index, p1_chain),
            CifAtomIdValue::from_table(table, row_index, p1_seq),
            CifAtomIdValue::from_table(table, row_index, p1_atom),
            p1_alt.map(|n| table.clean_at(row_index, n)).as_deref(),
            p1_ins.map(|n| table.clean_at(row_index, n)).as_deref(),
        );
        let b_candidates = find_cif_atoms(
            atoms,
            CifAtomIdValue::from_table(table, row_index, p2_chain),
            CifAtomIdValue::from_table(table, row_index, p2_seq),
            CifAtomIdValue::from_table(table, row_index, p2_atom),
            p2_alt.map(|n| table.clean_at(row_index, n)).as_deref(),
            p2_ins.map(|n| table.clean_at(row_index, n)).as_deref(),
        );
        let conn_type_id = conn_type
            .map(|n| table.clean_at(row_index, n))
            .unwrap_or_default();
        let order_value = value_order
            .map(|n| table.clean_at(row_index, n))
            .unwrap_or_default();
        let flags = struct_conn_flags(&conn_type_id);
        let p1_comp_id = p1_comp
            .map(|n| table.clean_at(row_index, n))
            .unwrap_or_default();
        let p2_comp_id = p2_comp
            .map(|n| table.clean_at(row_index, n))
            .unwrap_or_default();
        let p1_atom_id = table.clean_at(row_index, p1_atom.index);
        let p2_atom_id = table.clean_at(row_index, p2_atom.index);
        let order = struct_conn_order(
            &order_value,
            &p1_comp_id,
            &p1_atom_id,
            &p2_comp_id,
            &p2_atom_id,
        );
        let dist = distance.and_then(|n| table.float_at(row_index, n));
        let struct_conn = StructConnMetadata {
            id: row_id
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_else(|| (row_index + 1).to_string()),
            row_index,
            partner_a_atom_index: usize::MAX,
            partner_b_atom_index: usize::MAX,
            conn_type_id: conn_type_id.clone(),
            value_order: order_value.clone(),
            partner_a_symmetry: p1_symmetry
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_default(),
            partner_b_symmetry: p2_symmetry
                .map(|n| table.clean_at(row_index, n))
                .unwrap_or_default(),
            partner_a_comp_id: p1_comp_id,
            partner_b_comp_id: p2_comp_id,
        };
        for &a in &a_candidates {
            for &b in &b_candidates {
                if a == b || !alt_locs_are_compatible(&atoms[a], &atoms[b]) {
                    continue;
                }
                let bond = if a < b {
                    Bond { a, b }
                } else {
                    Bond { a: b, b: a }
                };
                let mut struct_conn = struct_conn.clone();
                struct_conn.partner_a_atom_index = a;
                struct_conn.partner_b_atom_index = b;
                bonds.push(bond);
                metadata.push(BondMetadata {
                    source: BondSource::StructConn,
                    order,
                    flags,
                    key: row_index as i32,
                    distance: dist,
                    operator_a: -1,
                    operator_b: -1,
                    struct_conn: Some(struct_conn),
                });
            }
        }
    }
    (bonds, metadata)
}

fn parse_molstar_bond_site_bonds(
    tables: &[CifTable],
    atoms: &[Atom],
) -> (Vec<Bond>, Vec<BondMetadata>, Vec<(usize, usize)>) {
    let Some(table) = tables.iter().find(|t| t.name == "molstar_bond_site") else {
        return (Vec::new(), Vec::new(), Vec::new());
    };
    let idx = |name: &str| table.header_index(name);
    let (Some(atom_id_1), Some(atom_id_2)) = (
        idx("_molstar_bond_site.atom_id_1"),
        idx("_molstar_bond_site.atom_id_2"),
    ) else {
        return (Vec::new(), Vec::new(), Vec::new());
    };
    let value_order = idx("_molstar_bond_site.value_order");
    let type_id = idx("_molstar_bond_site.type_id");
    let mut bonds = Vec::new();
    let mut metadata = Vec::new();
    let mut edges = Vec::new();
    for row_index in table.row_indices() {
        let Some(a_id) = table.usize_at(row_index, atom_id_1) else {
            continue;
        };
        let Some(b_id) = table.usize_at(row_index, atom_id_2) else {
            continue;
        };
        let Some(a) = atoms.iter().position(|atom| atom.id == a_id) else {
            continue;
        };
        let Some(b) = atoms.iter().position(|atom| atom.id == b_id) else {
            continue;
        };
        if a == b {
            continue;
        }
        let value_order = value_order
            .map(|n| table.clean_at(row_index, n))
            .unwrap_or_default();
        let type_id = type_id
            .map(|n| table.clean_at(row_index, n))
            .unwrap_or_default();
        let mut flags = molstar_bond_site_import_order_flags(&value_order);
        flags = flags.union(struct_conn_flags(&type_id));
        let bond = if a < b {
            Bond { a, b }
        } else {
            Bond { a: b, b: a }
        };
        bonds.push(bond);
        edges.push((a, b));
        metadata.push(BondMetadata {
            source: BondSource::IndexPair,
            order: molstar_bond_site_order(&value_order),
            flags,
            key: -1,
            distance: None,
            operator_a: -1,
            operator_b: -1,
            struct_conn: None,
        });
    }
    (bonds, metadata, edges)
}

fn append_unique_bonds(
    bonds: &mut Vec<Bond>,
    metadata: &mut Vec<BondMetadata>,
    incoming_bonds: Vec<Bond>,
    incoming_metadata: Vec<BondMetadata>,
) {
    for (bond, meta) in incoming_bonds.into_iter().zip(incoming_metadata) {
        if bonds
            .iter()
            .any(|existing| existing.a == bond.a && existing.b == bond.b)
        {
            continue;
        }
        bonds.push(bond);
        metadata.push(meta);
    }
}

#[derive(Clone, Copy, Debug)]
enum CifAtomIdKind {
    Label,
    Auth,
}

#[derive(Clone, Copy, Debug)]
struct CifAtomIdColumn {
    index: usize,
    kind: CifAtomIdKind,
}

impl CifAtomIdColumn {
    fn new(index: usize, kind: CifAtomIdKind) -> Self {
        CifAtomIdColumn { index, kind }
    }
}

#[derive(Debug)]
struct CifAtomIdValue {
    value: String,
    kind: CifAtomIdKind,
}

impl CifAtomIdValue {
    fn from_table(table: &CifTable, row: usize, column: CifAtomIdColumn) -> Self {
        CifAtomIdValue {
            value: table.clean_at(row, column.index),
            kind: column.kind,
        }
    }
}

fn find_cif_atoms(
    atoms: &[Atom],
    chain: CifAtomIdValue,
    residue_seq: CifAtomIdValue,
    name: CifAtomIdValue,
    alt_id: Option<&str>,
    insertion_code: Option<&str>,
) -> Vec<usize> {
    let alt_id = alt_id.unwrap_or("");
    let insertion_code = insertion_code.unwrap_or("");
    atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| {
            let chain_matches =
                cif_atom_id_matches(chain.kind, &atom.chain, &atom.auth_chain, &chain.value);
            let seq_matches = cif_atom_id_matches(
                residue_seq.kind,
                &atom.residue_seq,
                &atom.auth_residue_seq,
                &residue_seq.value,
            );
            let name_matches =
                cif_atom_id_matches(name.kind, &atom.name, &atom.auth_name, &name.value);
            let insertion_matches = atom.insertion_code == insertion_code;
            (chain_matches
                && seq_matches
                && name_matches
                && insertion_matches
                && (alt_id.is_empty() || atom.alt_id.is_empty() || atom.alt_id == alt_id))
                .then_some(index)
        })
        .collect()
}

fn cif_atom_id_matches(kind: CifAtomIdKind, label: &str, auth: &str, expected: &str) -> bool {
    match kind {
        CifAtomIdKind::Label => label == expected,
        CifAtomIdKind::Auth => auth == expected,
    }
}

fn struct_conn_order(
    value_order: &str,
    comp_id_1: &str,
    atom_id_1: &str,
    comp_id_2: &str,
    atom_id_2: &str,
) -> i8 {
    match value_order.to_ascii_lowercase().as_str() {
        "sing" => 1,
        "doub" => 2,
        "trip" => 3,
        "quad" => 4,
        _ => inter_bond_order_from_table(comp_id_1, atom_id_1, comp_id_2, atom_id_2),
    }
}

pub(crate) fn inter_bond_order_from_table(
    comp_id_1: &str,
    atom_id_1: &str,
    comp_id_2: &str,
    atom_id_2: &str,
) -> i8 {
    let mut a = (
        comp_id_1.to_ascii_uppercase(),
        atom_id_1.to_ascii_uppercase(),
    );
    let mut b = (
        comp_id_2.to_ascii_uppercase(),
        atom_id_2.to_ascii_uppercase(),
    );
    if a.0 > b.0 {
        std::mem::swap(&mut a, &mut b);
    }
    match (a.0.as_str(), a.1.as_str(), b.0.as_str(), b.1.as_str()) {
        ("LYS", "NZ", "RET", "C15") => 2,
        _ => 1,
    }
}

fn molstar_bond_site_order(value_order: &str) -> i8 {
    match value_order.to_ascii_lowercase().as_str() {
        "doub" => 2,
        "trip" => 3,
        "quad" => 4,
        _ => 1,
    }
}

fn molstar_bond_site_order_flags(value_order: &str) -> BondFlags {
    match value_order.to_ascii_lowercase().as_str() {
        "arom" | "delo" => BondFlags::AROMATIC.union(BondFlags::RESONANCE),
        _ => BondFlags::NONE,
    }
}

fn molstar_bond_site_import_order_flags(value_order: &str) -> BondFlags {
    match value_order.to_ascii_lowercase().as_str() {
        "arom" => BondFlags::AROMATIC,
        _ => BondFlags::NONE,
    }
}

fn struct_conn_flags(conn_type_id: &str) -> BondFlags {
    match conn_type_id.to_ascii_lowercase().as_str() {
        "disulf" => BondFlags::COVALENT.union(BondFlags::DISULFIDE),
        "hydrog" => BondFlags::HYDROGEN_BOND,
        "metalc" => BondFlags::METALLIC_COORDINATION,
        "covale" => BondFlags::COVALENT,
        _ => BondFlags::NONE,
    }
}

fn alt_locs_are_compatible(a: &Atom, b: &Atom) -> bool {
    a.alt_id.is_empty() || b.alt_id.is_empty() || a.alt_id == b.alt_id
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum StructConfRangeKind {
    Helix,
    Sheet,
}

fn secondary_ranges_from_struct_conf(
    tables: &[CifTable],
) -> (Vec<SecondaryRange>, Vec<SecondaryRange>) {
    let Some(table) = tables.iter().find(|t| t.name == "struct_conf") else {
        return (Vec::new(), Vec::new());
    };
    let idx = |suffix: &str| table.header_index(&format!("_struct_conf.{suffix}"));
    let Some(conf_type) = idx("conf_type_id") else {
        return (
            table
                .row_indices()
                .filter_map(|row| secondary_range_from_row(table, &idx, row))
                .collect(),
            Vec::new(),
        );
    };
    let mut helices = Vec::new();
    let mut sheets = Vec::new();

    for row in table.row_indices() {
        let Some(kind) = struct_conf_range_kind(&table.clean_at(row, conf_type)) else {
            continue;
        };
        let Some(range) = secondary_range_from_row(table, &idx, row) else {
            continue;
        };
        match kind {
            StructConfRangeKind::Helix => helices.push(range),
            StructConfRangeKind::Sheet => sheets.push(range),
        }
    }

    (helices, sheets)
}

fn struct_conf_range_kind(conf_type_id: &str) -> Option<StructConfRangeKind> {
    let value = conf_type_id.trim().to_ascii_lowercase();
    if value == "strn" {
        Some(StructConfRangeKind::Sheet)
    } else if value.starts_with("helx") || value.contains("helix") {
        Some(StructConfRangeKind::Helix)
    } else {
        None
    }
}

fn secondary_ranges_from_tables(tables: &[CifTable], category: &str) -> Vec<SecondaryRange> {
    let Some(table) = tables.iter().find(|t| t.name == category) else {
        return Vec::new();
    };
    let prefix = format!("_{category}.");
    let idx = |suffix: &str| table.header_index(&format!("{prefix}{suffix}"));
    table
        .row_indices()
        .filter_map(|row| secondary_range_from_row(table, &idx, row))
        .collect()
}

fn secondary_range_from_row(
    table: &CifTable,
    idx: &dyn Fn(&str) -> Option<usize>,
    row: usize,
) -> Option<SecondaryRange> {
    let (chain, start, end) = secondary_coordinate_columns(table, &idx);
    let start_insertion_code = idx("pdbx_beg_PDB_ins_code");
    let end_insertion_code = idx("pdbx_end_PDB_ins_code");
    let (Some(chain), Some(start), Some(end)) = (chain, start, end) else {
        return None;
    };
    Some(SecondaryRange {
        chain: table.clean_at(row, chain),
        start: table.i32_at(row, start)?,
        start_insertion_code: start_insertion_code
            .map(|idx| table.clean_at(row, idx))
            .unwrap_or_default(),
        end: table.i32_at(row, end)?,
        end_insertion_code: end_insertion_code
            .map(|idx| table.clean_at(row, idx))
            .unwrap_or_default(),
    })
}

fn secondary_coordinate_columns(
    table: &CifTable,
    idx: &dyn Fn(&str) -> Option<usize>,
) -> (Option<usize>, Option<usize>, Option<usize>) {
    let label_start = idx("beg_label_seq_id");
    let label_end = idx("end_label_seq_id");
    let use_label = table.row_count() > 0
        && label_start.is_some_and(|column| first_row_value_is_present(table, column))
        && label_end.is_some_and(|column| first_row_value_is_present(table, column));

    let chain = idx("beg_label_asym_id")
        .filter(|&column| first_row_value_is_present(table, column))
        .or_else(|| idx("beg_auth_asym_id"));
    if use_label {
        (chain, label_start, label_end)
    } else {
        (chain, idx("beg_auth_seq_id"), idx("end_auth_seq_id"))
    }
}

fn first_row_value_is_present(table: &CifTable, column: usize) -> bool {
    table.row_count() > 0 && is_present_cif_value(&table.raw_at(0, column))
}

fn parse_global_model_transform(tables: &[CifTable]) -> Option<GlobalModelTransform> {
    let table = tables
        .iter()
        .find(|table| table.name == GlobalModelTransform::DESCRIPTOR)?;
    if table.row_count() == 0 {
        return None;
    }
    let prefix = format!("_{}.", GlobalModelTransform::DESCRIPTOR);
    let mut matrix = [[0.0; 4]; 4];
    for (row, matrix_row) in matrix.iter_mut().enumerate() {
        for (col, value) in matrix_row.iter_mut().enumerate() {
            let name = format!("{prefix}matrix[{}][{}]", row + 1, col + 1);
            let index = table.header_index(&name)?;
            *value = table.float_at(0, index)?;
        }
    }
    Some(GlobalModelTransform { matrix })
}

fn parse_cif_assemblies(tables: &[CifTable]) -> Vec<Assembly> {
    let Some(gen) = tables.iter().find(|t| t.name == "pdbx_struct_assembly_gen") else {
        return Vec::new();
    };
    let Some(ops) = tables.iter().find(|t| t.name == "pdbx_struct_oper_list") else {
        return Vec::new();
    };

    let op_idx = |name: &str| ops.header_index(name);
    let Some(op_id) = op_idx("_pdbx_struct_oper_list.id") else {
        return Vec::new();
    };
    let m = [
        op_idx("_pdbx_struct_oper_list.matrix[1][1]"),
        op_idx("_pdbx_struct_oper_list.matrix[1][2]"),
        op_idx("_pdbx_struct_oper_list.matrix[1][3]"),
        op_idx("_pdbx_struct_oper_list.vector[1]"),
        op_idx("_pdbx_struct_oper_list.matrix[2][1]"),
        op_idx("_pdbx_struct_oper_list.matrix[2][2]"),
        op_idx("_pdbx_struct_oper_list.matrix[2][3]"),
        op_idx("_pdbx_struct_oper_list.vector[2]"),
        op_idx("_pdbx_struct_oper_list.matrix[3][1]"),
        op_idx("_pdbx_struct_oper_list.matrix[3][2]"),
        op_idx("_pdbx_struct_oper_list.matrix[3][3]"),
        op_idx("_pdbx_struct_oper_list.vector[3]"),
    ];
    if m.iter().any(Option::is_none) {
        return Vec::new();
    }
    let m = m.map(Option::unwrap);
    let mut op_map = Vec::<(String, Transform)>::new();
    for row_index in ops.row_indices() {
        let values: Vec<f32> = m
            .iter()
            .map(|&idx| ops.float_at(row_index, idx).unwrap_or(0.0))
            .collect();
        op_map.push((
            ops.clean_at(row_index, op_id),
            Transform {
                m: [
                    [values[0], values[1], values[2], values[3]],
                    [values[4], values[5], values[6], values[7]],
                    [values[8], values[9], values[10], values[11]],
                ],
            },
        ));
    }

    let gen_idx = |name: &str| gen.header_index(name);
    let (Some(assembly_id), Some(oper_expression), Some(asym_id_list)) = (
        gen_idx("_pdbx_struct_assembly_gen.assembly_id"),
        gen_idx("_pdbx_struct_assembly_gen.oper_expression"),
        gen_idx("_pdbx_struct_assembly_gen.asym_id_list"),
    ) else {
        return Vec::new();
    };
    let assembly_metadata = parse_cif_assembly_metadata(tables);

    let mut assemblies = Vec::<Assembly>::new();
    for row_index in gen.row_indices() {
        let id = gen.clean_at(row_index, assembly_id);
        let op_id_groups = expand_oper_expression(&gen.clean_at(row_index, oper_expression));
        let mut transforms = Vec::new();
        let mut oper_list_ids = Vec::new();
        for ids in op_id_groups {
            if let Some(transform) = compose_operator_transforms(&ids, &op_map) {
                transforms.push(transform);
                oper_list_ids.push(ids);
            }
        }
        if transforms.is_empty() {
            continue;
        }
        let asym_ids: Vec<String> = gen
            .clean_at(row_index, asym_id_list)
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect();
        if let Some(existing) = assemblies.iter_mut().find(|a: &&mut Assembly| a.id == id) {
            let start_oper_id = existing.transforms.len();
            existing.transforms.extend(transforms.iter().copied());
            existing.generators.push(AssemblyGenerator::from_transforms(
                &id,
                asym_ids.clone(),
                start_oper_id,
                transforms,
                oper_list_ids,
            ));
            for asym_id in asym_ids {
                if !existing.asym_ids.iter().any(|id| id == &asym_id) {
                    existing.asym_ids.push(asym_id);
                }
            }
        } else {
            let metadata = assembly_metadata
                .iter()
                .find(|metadata| metadata.id == id)
                .cloned()
                .unwrap_or_default();
            assemblies.push(Assembly {
                id: id.clone(),
                details: metadata.details,
                oligomeric_details: metadata.oligomeric_details,
                oligomeric_count: metadata.oligomeric_count,
                asym_ids: asym_ids.clone(),
                transforms: transforms.clone(),
                generators: vec![AssemblyGenerator::from_transforms(
                    &id,
                    asym_ids,
                    0,
                    transforms,
                    oper_list_ids,
                )],
            });
        }
    }
    assemblies
}

#[derive(Clone, Debug, Default)]
struct AssemblyMetadata {
    id: String,
    details: String,
    oligomeric_details: String,
    oligomeric_count: Option<i32>,
}

fn parse_cif_assembly_metadata(tables: &[CifTable]) -> Vec<AssemblyMetadata> {
    let Some(table) = tables.iter().find(|t| t.name == "pdbx_struct_assembly") else {
        return Vec::new();
    };
    let Some(id_idx) = table.header_index("_pdbx_struct_assembly.id") else {
        return Vec::new();
    };
    let details_idx = table.header_index("_pdbx_struct_assembly.details");
    let oligomeric_details_idx = table.header_index("_pdbx_struct_assembly.oligomeric_details");
    let oligomeric_count_idx = table.header_index("_pdbx_struct_assembly.oligomeric_count");

    table
        .row_indices()
        .map(|row_index| AssemblyMetadata {
            id: table.clean_at(row_index, id_idx),
            details: details_idx
                .map(|idx| table.clean_at(row_index, idx))
                .unwrap_or_default(),
            oligomeric_details: oligomeric_details_idx
                .map(|idx| table.clean_at(row_index, idx))
                .unwrap_or_default(),
            oligomeric_count: oligomeric_count_idx
                .and_then(|idx| table.clean_at(row_index, idx).parse::<i32>().ok()),
        })
        .collect()
}

pub(crate) fn expand_oper_expression(value: &str) -> Vec<Vec<String>> {
    let groups = oper_expression_groups(value);
    let operator_names = groups
        .iter()
        .map(|group| expand_oper_group(group))
        .filter(|group| !group.is_empty())
        .collect::<Vec<_>>();
    if operator_names.is_empty() {
        return Vec::new();
    }
    let mut out = Vec::new();
    let mut current = vec![String::new(); operator_names.len()];
    expand_oper_expression_molstar(
        &operator_names,
        operator_names.len() - 1,
        &mut current,
        &mut out,
    );
    out
}

fn expand_oper_expression_molstar(
    operator_names: &[Vec<String>],
    group_index: usize,
    current: &mut [String],
    out: &mut Vec<Vec<String>>,
) {
    for value in &operator_names[group_index] {
        current[group_index] = value.clone();
        if group_index == 0 {
            out.push(current.to_vec());
        } else {
            expand_oper_expression_molstar(operator_names, group_index - 1, current, out);
        }
    }
}

fn oper_expression_groups(value: &str) -> Vec<String> {
    let trimmed = value.trim();
    if !trimmed.contains('(') {
        return vec![trimmed.to_string()];
    }
    let mut groups = Vec::new();
    let mut current = String::new();
    let mut in_group = false;
    for ch in trimmed.chars() {
        match ch {
            '(' => {
                if !current.trim().is_empty() {
                    groups.push(current.trim().to_string());
                    current.clear();
                }
                in_group = true;
            }
            ')' => {
                if in_group && !current.trim().is_empty() {
                    groups.push(current.trim().to_string());
                }
                current.clear();
                in_group = false;
            }
            ']' if !in_group && current.trim().is_empty() => {}
            _ => current.push(ch),
        }
    }
    if !current.trim().is_empty() {
        groups.push(current.trim().to_string());
    }
    groups
}

fn expand_oper_group(value: &str) -> Vec<String> {
    let mut out = Vec::new();
    for part in value.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        if let Some((a, b)) = part.split_once('-') {
            if let (Ok(a), Ok(b)) = (a.trim().parse::<usize>(), b.trim().parse::<usize>()) {
                out.extend((a..=b).map(|i| i.to_string()));
                continue;
            }
        }
        out.push(part.to_string());
    }
    out
}

pub(crate) fn compose_operator_transforms(
    ids: &[String],
    op_map: &[(String, Transform)],
) -> Option<Transform> {
    let mut combined = Transform::identity();
    for id in ids {
        let transform = op_map.iter().find(|(op_id, _)| op_id == id)?.1;
        combined = transform.then(combined);
    }
    Some(combined)
}

mod binary;

pub(crate) use binary::ColumnData;

fn normalize_type_symbol_molstar(value: &str) -> String {
    value.trim().to_ascii_uppercase()
}

fn select_alt_loc(mut molecule: Molecule, requested: &str) -> Molecule {
    let requested = requested.trim();
    if requested == "all" {
        return molecule;
    }
    let source_atoms = molecule.atoms;
    let mut index_map = vec![None; source_atoms.len()];
    let mut chosen = Vec::<Atom>::new();
    if requested.is_empty() || requested == "highest-occupancy" {
        let mut best_by_site = Vec::<((i32, String, String, String, String, String), usize)>::new();
        for (source_index, atom) in source_atoms.iter().enumerate() {
            let key = (
                atom.model_num,
                atom.chain.clone(),
                atom.residue_seq.clone(),
                atom.insertion_code.clone(),
                atom.residue.clone(),
                atom.name.clone(),
            );
            if let Some((_, best_index)) = best_by_site.iter_mut().find(|(site, _)| site == &key) {
                let best = &source_atoms[*best_index];
                let prefers_atom = atom.occupancy > best.occupancy
                    || (atom.occupancy == best.occupancy
                        && best.alt_id.is_empty()
                        && !atom.alt_id.is_empty())
                    || (atom.occupancy == best.occupancy
                        && atom.alt_id == "A"
                        && best.alt_id != "A");
                if prefers_atom {
                    *best_index = source_index;
                }
            } else {
                best_by_site.push((key, source_index));
            }
        }
        let mut keep = best_by_site
            .into_iter()
            .map(|(_, source_index)| source_index)
            .collect::<Vec<_>>();
        keep.sort_unstable();
        for source_index in keep {
            index_map[source_index] = Some(chosen.len());
            chosen.push(source_atoms[source_index].clone());
        }
    } else {
        for (source_index, atom) in source_atoms.into_iter().enumerate() {
            if atom.alt_id.is_empty() || atom.alt_id == requested {
                index_map[source_index] = Some(chosen.len());
                chosen.push(atom);
            }
        }
    }
    molecule.atoms = chosen;
    let (bonds, bond_metadata) = remap_bonds(&molecule.bonds, &molecule.bond_metadata, &index_map);
    molecule.bonds = bonds;
    molecule.bond_metadata = bond_metadata;
    if let Some(index_pair_bonds) = molecule.index_pair_bonds.as_ref() {
        molecule.index_pair_bonds = IndexPairBonds::from_bonds(
            &molecule.bonds,
            &molecule.bond_metadata,
            molecule.atoms.len(),
            index_pair_bonds.max_distance,
            index_pair_bonds.cacheable,
        );
    }
    molecule
}

pub(crate) fn unique_alt_locs(atoms: &[Atom]) -> Vec<String> {
    let mut out = Vec::new();
    for atom in atoms {
        if !atom.alt_id.is_empty() && !out.iter().any(|id| id == &atom.alt_id) {
            out.push(atom.alt_id.clone());
        }
    }
    out
}

fn remap_bonds(
    bonds: &[Bond],
    metadata: &[BondMetadata],
    index_map: &[Option<usize>],
) -> (Vec<Bond>, Vec<BondMetadata>) {
    let pairs = bonds
        .iter()
        .enumerate()
        .filter_map(|(index, bond)| {
            let a = index_map.get(bond.a).and_then(|i| *i)?;
            let b = index_map.get(bond.b).and_then(|i| *i)?;
            Some((
                Bond { a, b },
                metadata.get(index).cloned().unwrap_or_default(),
            ))
        })
        .collect::<Vec<_>>();
    pairs.into_iter().unzip()
}

fn apply_assembly(mut molecule: Molecule, id: &str) -> Molecule {
    let Some(assembly) = molecule.assemblies.iter().find(|a| a.id == id).cloned() else {
        return molecule;
    };
    molecule.selected_assembly = Some(assembly);
    molecule
}

#[cfg(test)]
mod tests;
