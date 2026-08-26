use crate::chemistry::{infer_element_from_name, normalize_element};
use crate::model::{
    Assembly, AssemblyGenerator, Atom, AtomSiteAnisotrop, AtomSiteColumnPresence, Bond, BondFlags,
    BondMetadata, BondSource, Entity, EntityIndexMap, EntityPolySeq, Molecule, SecondaryRange,
    SourceData, StructConnMetadata, Transform, Vec3,
};
use std::collections::BTreeMap;

use super::normalize_type_symbol_molstar;

pub(crate) fn parse_pdb(text: &str) -> Result<Molecule, String> {
    let mut atoms = Vec::new();
    let mut serial_to_index = Vec::<(usize, usize)>::new();
    let mut atom_site_anisotrop = Vec::new();
    let mut b_iso_defined = false;
    let assemblies = parse_pdb_assemblies(text);
    let (mut helices, mut sheets) = parse_pdb_secondary(text);
    let mut model_num = 1;
    let mut next_model_num = 1;

    for line in text.lines() {
        if line.starts_with("MODEL ") {
            model_num = next_model_num;
            next_model_num += 1;
        } else if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            let serial = field(line, 6, 11)
                .trim()
                .parse::<usize>()
                .unwrap_or(atoms.len() + 1);
            let name = field(line, 12, 16).trim().to_string();
            let alt_id = field(line, 16, 17).trim().to_string();
            let residue = field(line, 17, 20).trim().to_string();
            let chain = field(line, 21, 22).trim().to_string();
            let residue_seq = field(line, 22, 26).trim().to_string();
            let insertion_code = field(line, 26, 27).trim().to_string();
            let x64 = parse_f64(field(line, 30, 38))?;
            let y64 = parse_f64(field(line, 38, 46))?;
            let z64 = parse_f64(field(line, 46, 54))?;
            let (x, y, z) = (x64 as f32, y64 as f32, z64 as f32);
            let occupancy = parse_js_number_f32(field(line, 54, 60)).unwrap_or(1.0);
            b_iso_defined |= !field(line, 60, 66).trim().is_empty();
            let b_iso = parse_js_number_f32(field(line, 60, 66)).unwrap_or(0.0);
            let formal_charge = parse_pdb_charge(field(line, 78, 80));
            let element = normalize_element({
                let e = field(line, 76, 78).trim();
                if e.is_empty() {
                    infer_element_from_name(&name)
                } else {
                    e.to_string()
                }
            });
            serial_to_index.push((serial, atoms.len()));
            atoms.push(Atom {
                id: serial,
                source_index: atoms.len(),
                model_num,
                auth_name: name.clone(),
                type_symbol: normalize_type_symbol_molstar(&element),
                name,
                element,
                chain: chain.clone(),
                auth_chain: chain,
                entity_id: String::new(),
                residue: residue.clone(),
                auth_residue: residue,
                group_pdb: if line.starts_with("HETATM") {
                    "HETATM".to_string()
                } else {
                    "ATOM".to_string()
                },
                residue_seq: residue_seq.clone(),
                auth_residue_seq: residue_seq,
                insertion_code,
                alt_id,
                occupancy,
                b_iso,
                formal_charge,
                position: Vec3 { x, y, z },
                position64: [x as f64, y as f64, z as f64],
                het: line.starts_with("HETATM"),
                operator_name: String::new(),
            });
        } else if line.starts_with("ANISOU") {
            if let Some(row) = parse_pdb_anisou(line) {
                atom_site_anisotrop.push(row);
            }
        }
    }

    if atoms.is_empty() {
        return Err("no ATOM/HETATM records found in PDB input".to_string());
    }
    assign_pdb_label_atom_ids(&mut atoms);
    let (bonds, bond_metadata) = parse_pdb_conect(text, &atoms, &serial_to_index);
    let entities = assign_pdb_entities(text, &mut atoms);
    let seqres = parse_pdb_seqres(text);
    assign_pdb_label_seq_ids(&seqres, &mut atoms);
    normalize_pdb_secondary_ranges(&mut helices, &atoms);
    normalize_pdb_secondary_ranges(&mut sheets, &atoms);
    let entity_poly_seq = pdb_entity_poly_seq(&seqres, &atoms, &entities);
    let entity_index = EntityIndexMap::from_mmcif(&entities, &[], &[], &atoms, &[], &[], &[]);
    Ok(Molecule {
        source_data: SourceData::pdb(pdb_id(text)),
        atom_site_columns: AtomSiteColumnPresence {
            occupancy_defined: true,
            b_iso_defined,
            xyz_defined: true,
        },
        global_model_transform: None,
        entries: Vec::new(),
        experiments: Vec::new(),
        atoms,
        atom_site_anisotrop,
        bonds,
        bond_metadata,
        index_pair_bonds: None,
        coarse_spheres: Vec::new(),
        coarse_gaussians: Vec::new(),
        assemblies,
        selected_assembly: None,
        helices,
        sheets,
        entities,
        entity_index,
        entity_polymers: Vec::new(),
        entity_poly_seq,
        pdbx_entity_branch: Vec::new(),
        pdbx_entity_branch_links: Vec::new(),
        pdbx_branch_scheme: Vec::new(),
        pdbx_nonpoly_scheme: Vec::new(),
        pdbx_poly_seq_scheme: Vec::new(),
        ihm_model_list: Vec::new(),
        ihm_model_groups: Vec::new(),
        ihm_model_group_links: Vec::new(),
        ihm_cross_link_restraints: Vec::new(),
        struct_asym: Vec::new(),
        pdbx_molecule: Vec::new(),
        chemical_components: Vec::new(),
        chemical_component_atoms: Vec::new(),
        chemical_component_bonds: Vec::new(),
        chemical_component_angles: Vec::new(),
        quality_assessment: Default::default(),
        partial_charges: Default::default(),
        rings: Vec::new(),
        resonance: Default::default(),
        derived_aromatic_bonds: Default::default(),
        derived_resonance_bonds: Default::default(),
    })
}

fn assign_pdb_label_atom_ids(atoms: &mut [Atom]) {
    let mut current_residue = None::<(i32, String, String, String)>;
    let mut counts = BTreeMap::<String, usize>::new();
    for atom in atoms {
        let residue = (
            atom.model_num,
            atom.auth_chain.clone(),
            atom.auth_residue_seq.clone(),
            atom.insertion_code.clone(),
        );
        if current_residue.as_ref() != Some(&residue) {
            current_residue = Some(residue);
            counts.clear();
        }
        let count = counts.entry(atom.auth_name.clone()).or_default();
        atom.name = if *count == 0 {
            atom.auth_name.clone()
        } else {
            format!("{}_{}", atom.auth_name, *count)
        };
        *count += 1;
    }
}

fn parse_pdb_conect(
    text: &str,
    atoms: &[Atom],
    serial_to_index: &[(usize, usize)],
) -> (Vec<Bond>, Vec<BondMetadata>) {
    let mut bonds = Vec::<Bond>::new();
    let mut metadata = Vec::<BondMetadata>::new();
    let mut current_atom = None::<usize>;
    let mut bond_index = BTreeMap::<usize, usize>::new();

    for line in text.lines().filter(|line| line.starts_with("CONECT")) {
        let Ok(serial_a) = field(line, 6, 11).trim().parse::<usize>() else {
            continue;
        };
        let Some(atom_a) = lookup_serial(serial_to_index, serial_a) else {
            continue;
        };
        if current_atom != Some(atom_a) {
            current_atom = Some(atom_a);
            bond_index.clear();
        }

        for (start, end) in [(11, 16), (16, 21), (21, 26), (26, 31)] {
            let Ok(serial_b) = field(line, start, end).trim().parse::<usize>() else {
                continue;
            };
            let Some(atom_b) = lookup_serial(serial_to_index, serial_b) else {
                continue;
            };
            if atom_a > atom_b {
                continue;
            }
            if let Some(&index) = bond_index.get(&atom_b) {
                if let Some(entry) = metadata.get_mut(index) {
                    entry.order = entry.order.saturating_add(1);
                    if let Some(struct_conn) = &mut entry.struct_conn {
                        struct_conn.value_order = pdb_conect_value_order(entry.order).to_string();
                    }
                }
                continue;
            }

            let Some((a, b)) = atoms.get(atom_a).zip(atoms.get(atom_b)) else {
                continue;
            };
            let metallic = BondMetadata::computed_for_atoms(a, b)
                .flags
                .contains(BondFlags::METALLIC_COORDINATION);
            let conn_type_id = if metallic { "metalc" } else { "covale" };
            let row_index = metadata.len();
            bonds.push(Bond {
                a: atom_a,
                b: atom_b,
            });
            metadata.push(BondMetadata {
                source: BondSource::StructConn,
                order: 1,
                flags: if metallic {
                    BondFlags::METALLIC_COORDINATION
                } else {
                    BondFlags::COVALENT
                },
                key: row_index as i32,
                distance: None,
                operator_a: -1,
                operator_b: -1,
                struct_conn: Some(StructConnMetadata {
                    id: format!("{conn_type_id}{}", row_index + 1),
                    row_index,
                    partner_a_atom_index: atom_a,
                    partner_b_atom_index: atom_b,
                    conn_type_id: conn_type_id.to_string(),
                    value_order: "sing".to_string(),
                    partner_a_symmetry: String::new(),
                    partner_b_symmetry: String::new(),
                    partner_a_comp_id: a.residue.clone(),
                    partner_b_comp_id: b.residue.clone(),
                }),
            });
            bond_index.insert(atom_b, row_index);
        }
    }

    (bonds, metadata)
}

fn pdb_conect_value_order(order: i8) -> &'static str {
    match order {
        2 => "doub",
        3 => "trip",
        4 => "quad",
        _ => "sing",
    }
}

fn parse_pdb_seqres(text: &str) -> Vec<(String, Vec<String>)> {
    let mut chains = Vec::<(String, Vec<String>)>::new();
    for line in text.lines().filter(|line| line.starts_with("SEQRES")) {
        let chain = field(line, 11, 12).to_string();
        let residues = field(line, 19, line.len())
            .split_whitespace()
            .map(str::to_string)
            .collect::<Vec<_>>();
        if let Some((_, existing)) = chains.iter_mut().find(|(id, _)| id == &chain) {
            existing.extend(residues);
        } else {
            chains.push((chain, residues));
        }
    }
    chains
}

fn pdb_entity_poly_seq(
    seqres: &[(String, Vec<String>)],
    atoms: &[Atom],
    entities: &[Entity],
) -> Vec<EntityPolySeq> {
    let mut processed = Vec::<String>::new();
    let mut rows = Vec::new();
    for (chain, residues) in seqres {
        let Some(entity_id) = atoms.iter().find_map(|atom| {
            if atom.auth_chain != *chain {
                return None;
            }
            entities
                .iter()
                .find(|entity| entity.id == atom.entity_id && entity.type_name == "polymer")
                .map(|entity| entity.id.clone())
        }) else {
            continue;
        };
        if processed.iter().any(|id| id == &entity_id) {
            continue;
        }
        processed.push(entity_id.clone());
        rows.extend(
            residues
                .iter()
                .enumerate()
                .map(|(index, residue)| EntityPolySeq {
                    entity_id: entity_id.clone(),
                    num: index as i32 + 1,
                    mon_id: residue.clone(),
                    hetero: "no".to_string(),
                }),
        );
    }
    rows
}

fn assign_pdb_label_seq_ids(seqres: &[(String, Vec<String>)], atoms: &mut [Atom]) {
    if atoms.is_empty() {
        return;
    }
    let use_linear = !seqres.is_empty() || atoms.iter().any(|atom| !atom.insertion_code.is_empty());
    if !use_linear {
        for atom in atoms {
            atom.residue_seq.clear();
        }
        return;
    }

    let first_model = atoms[0].model_num;
    let mut observed = Vec::<(String, Vec<String>)>::new();
    let mut previous = None::<(String, String, String)>;
    for atom in atoms
        .iter()
        .take_while(|atom| atom.model_num == first_model)
    {
        let key = (
            atom.auth_chain.clone(),
            atom.auth_residue_seq.clone(),
            atom.insertion_code.clone(),
        );
        if previous.as_ref() == Some(&key) {
            continue;
        }
        if let Some((_, residues)) = observed
            .iter_mut()
            .find(|(chain, _)| chain == &atom.auth_chain)
        {
            residues.push(atom.auth_residue.clone());
        } else {
            observed.push((atom.auth_chain.clone(), vec![atom.auth_residue.clone()]));
        }
        previous = Some(key);
    }

    let alignments = observed
        .into_iter()
        .filter_map(|(chain, observed)| {
            let sequence = seqres.iter().find(|(id, _)| id == &chain)?.1.as_slice();
            Some((chain, pdb_seqres_alignment(&observed, sequence)))
        })
        .collect::<Vec<_>>();

    let mut current_model = atoms[0].model_num;
    let mut current_chain = atoms[0].auth_chain.clone();
    let mut current_auth_seq = pdb_integer(&atoms[0].auth_residue_seq);
    let mut current_insertion = atoms[0].insertion_code.clone();
    let mut residue_index = 0usize;
    let mut current_label_seq = initial_pdb_label_seq(
        alignment_for_chain(&alignments, &current_chain),
        residue_index,
        current_auth_seq,
    );

    for atom in atoms {
        let auth_seq = pdb_integer(&atom.auth_residue_seq);
        if atom.model_num != current_model {
            current_model = atom.model_num;
            current_chain = atom.auth_chain.clone();
            current_auth_seq = auth_seq;
            current_insertion = atom.insertion_code.clone();
            residue_index = 0;
            current_label_seq = initial_pdb_label_seq(
                alignment_for_chain(&alignments, &current_chain),
                residue_index,
                current_auth_seq,
            );
        } else if atom.auth_chain != current_chain {
            current_chain = atom.auth_chain.clone();
            current_auth_seq = auth_seq;
            current_insertion = atom.insertion_code.clone();
            residue_index = 0;
            current_label_seq = initial_pdb_label_seq(
                alignment_for_chain(&alignments, &current_chain),
                residue_index,
                current_auth_seq,
            );
        } else if auth_seq != current_auth_seq || atom.insertion_code != current_insertion {
            residue_index += 1;
            current_label_seq = alignment_for_chain(&alignments, &current_chain)
                .and_then(|alignment| alignment.get(residue_index).copied().flatten())
                .unwrap_or(current_label_seq + 1);
            current_auth_seq = auth_seq;
            current_insertion = atom.insertion_code.clone();
        }
        atom.residue_seq = current_label_seq.to_string();
    }
}

fn alignment_for_chain<'a>(
    alignments: &'a [(String, Vec<Option<i32>>)],
    chain: &str,
) -> Option<&'a [Option<i32>]> {
    alignments
        .iter()
        .find(|(id, _)| id == chain)
        .map(|(_, alignment)| alignment.as_slice())
}

fn initial_pdb_label_seq(
    alignment: Option<&[Option<i32>]>,
    residue_index: usize,
    _auth_seq: i32,
) -> i32 {
    alignment
        .and_then(|alignment| alignment.get(residue_index).copied().flatten())
        .unwrap_or(1)
}

fn pdb_integer(value: &str) -> i32 {
    value.trim().parse::<i32>().unwrap_or(0)
}

fn pdb_seqres_alignment(observed: &[String], sequence: &[String]) -> Vec<Option<i32>> {
    const GAP: i32 = -11;
    const EXTEND: i32 = -1;
    const NEG_INFINITY: i32 = i32::MIN / 4;

    let n = observed.len();
    let m = sequence.len();
    let mut score = vec![vec![0i32; m + 1]; n + 1];
    let mut vertical = vec![vec![0i32; m + 1]; n + 1];
    let mut horizontal = vec![vec![0i32; m + 1]; n + 1];
    for row in 0..=n {
        score[row][0] = GAP;
        horizontal[row][0] = NEG_INFINITY;
    }
    for column in 0..=m {
        score[0][column] = GAP;
        vertical[0][column] = NEG_INFINITY;
    }
    score[0][0] = 0;

    for row in 1..=n {
        for column in 1..=m {
            vertical[row][column] =
                (score[row - 1][column] + GAP).max(vertical[row - 1][column] + EXTEND);
            horizontal[row][column] =
                (score[row][column - 1] + GAP).max(horizontal[row][column - 1] + EXTEND);
            let substitution = if observed[row - 1] == sequence[column - 1] {
                5
            } else {
                -3
            };
            score[row][column] = (score[row - 1][column - 1] + substitution)
                .max(vertical[row][column])
                .max(horizontal[row][column]);
        }
    }

    #[derive(Clone, Copy)]
    enum Matrix {
        Score,
        Vertical,
        Horizontal,
    }
    let mut row = n;
    let mut column = m;
    let mut matrix = if score[row][column] >= vertical[row][column] {
        Matrix::Score
    } else if vertical[row][column] >= horizontal[row][column] {
        Matrix::Vertical
    } else {
        Matrix::Horizontal
    };
    let mut trace = Vec::<(Option<usize>, Option<usize>)>::new();
    while row > 0 && column > 0 {
        match matrix {
            Matrix::Score => {
                let substitution = if observed[row - 1] == sequence[column - 1] {
                    5
                } else {
                    -3
                };
                if score[row][column] == score[row - 1][column - 1] + substitution {
                    trace.push((Some(row - 1), Some(column - 1)));
                    row -= 1;
                    column -= 1;
                } else if score[row][column] == vertical[row][column] {
                    matrix = Matrix::Vertical;
                } else if score[row][column] == horizontal[row][column] {
                    matrix = Matrix::Horizontal;
                } else {
                    row -= 1;
                    column -= 1;
                }
            }
            Matrix::Vertical => {
                if vertical[row][column] == vertical[row - 1][column] + EXTEND {
                    trace.push((Some(row - 1), None));
                    row -= 1;
                } else if vertical[row][column] == score[row - 1][column] + GAP {
                    trace.push((Some(row - 1), None));
                    row -= 1;
                    matrix = Matrix::Score;
                } else {
                    row -= 1;
                }
            }
            Matrix::Horizontal => {
                if horizontal[row][column] == horizontal[row][column - 1] + EXTEND {
                    trace.push((None, Some(column - 1)));
                    column -= 1;
                } else if horizontal[row][column] == score[row][column - 1] + GAP {
                    trace.push((None, Some(column - 1)));
                    column -= 1;
                    matrix = Matrix::Score;
                } else {
                    column -= 1;
                }
            }
        }
    }
    while row > 0 {
        trace.push((Some(row - 1), None));
        row -= 1;
    }
    while column > 0 {
        trace.push((None, Some(column - 1)));
        column -= 1;
    }
    trace.reverse();

    let mut alignment = vec![None; n];
    for (observed_index, sequence_index) in trace {
        if let (Some(observed_index), Some(sequence_index)) = (observed_index, sequence_index) {
            alignment[observed_index] = Some(sequence_index as i32 + 1);
        }
    }
    alignment
}

fn parse_pdb_anisou(line: &str) -> Option<AtomSiteAnisotrop> {
    let atom_id = field(line, 6, 11).trim().parse::<usize>().ok()?;
    let u11 = parse_pdb_anisou_value(field(line, 28, 35))?;
    let u22 = parse_pdb_anisou_value(field(line, 35, 42))?;
    let u33 = parse_pdb_anisou_value(field(line, 42, 49))?;
    let u12 = parse_pdb_anisou_value(field(line, 49, 56))?;
    let u13 = parse_pdb_anisou_value(field(line, 56, 63))?;
    let u23 = parse_pdb_anisou_value(field(line, 63, 70))?;
    Some(AtomSiteAnisotrop {
        atom_id,
        u: [[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]],
    })
}

fn parse_pdb_anisou_value(value: &str) -> Option<f32> {
    value.trim().parse::<f32>().ok().map(|v| v / 10000.0)
}

fn pdb_id(text: &str) -> String {
    text.lines()
        .find(|line| line.starts_with("HEADER") && line.len() >= 66)
        .map(|line| line[62..66].trim().to_string())
        .filter(|id| !id.is_empty())
        .unwrap_or_default()
}

fn assign_pdb_entities(text: &str, atoms: &mut [Atom]) -> Vec<Entity> {
    let compounds = pdb_compound_chains(text);
    let hetero_names = pdb_hetero_names(text);
    let mut entities = Vec::<Entity>::new();
    let mut polymer_entities = Vec::<(String, String)>::new();
    let mut non_polymer_entities = Vec::<(String, String)>::new();
    let mut water_entity = None::<String>;

    for (chain, description) in compounds {
        let description = if description.is_empty() {
            format!("Polymer {}", entities.len() + 1)
        } else {
            description
        };
        let id = next_pdb_entity(&mut entities, "polymer", description);
        polymer_entities.push((chain, id));
    }

    for atom in atoms {
        let entity_type = crate::model::entity_type_from_component(&atom.residue);
        let id = if entity_type == "polymer" {
            if let Some((_, id)) = polymer_entities
                .iter()
                .find(|(chain, _)| chain == &atom.auth_chain)
            {
                id.clone()
            } else {
                let id = next_pdb_entity(
                    &mut entities,
                    "polymer",
                    format!("Polymer {}", polymer_entities.len() + 1),
                );
                polymer_entities.push((atom.auth_chain.clone(), id.clone()));
                id
            }
        } else if entity_type == "water" {
            if let Some(id) = &water_entity {
                id.clone()
            } else {
                let id = next_pdb_entity(&mut entities, "water", "Water".to_string());
                water_entity = Some(id.clone());
                id
            }
        } else if let Some((_, id)) = non_polymer_entities
            .iter()
            .find(|(component, _)| component == &atom.residue)
        {
            id.clone()
        } else {
            let description = hetero_names
                .iter()
                .find(|(component, _)| component == &atom.residue)
                .map(|(_, description)| description.clone())
                .unwrap_or_else(|| atom.residue.clone());
            let id = next_pdb_entity(&mut entities, entity_type, description);
            non_polymer_entities.push((atom.residue.clone(), id.clone()));
            id
        };
        atom.entity_id = id;
    }

    entities
}

fn next_pdb_entity(entities: &mut Vec<Entity>, type_name: &str, description: String) -> String {
    let id = (entities.len() + 1).to_string();
    entities.push(Entity {
        id: id.clone(),
        type_name: type_name.to_string(),
        description,
    });
    id
}

fn pdb_compound_chains(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut current_field = "";
    let mut description = String::new();
    for line in text.lines().filter(|line| line.starts_with("COMPND")) {
        let payload = field(line, 10, 80).trim();
        let (field_name, value) = payload
            .split_once(':')
            .map(|(name, value)| (name.trim(), value.trim()))
            .unwrap_or((current_field, payload));
        current_field = field_name;
        let value = value.trim_end_matches(';').trim();
        match current_field {
            "MOL_ID" => description.clear(),
            "MOLECULE" => {
                if !description.is_empty() && !value.is_empty() {
                    description.push(' ');
                }
                description.push_str(value);
            }
            "CHAIN" => {
                out.extend(
                    value
                        .split(',')
                        .map(str::trim)
                        .filter(|chain| !chain.is_empty())
                        .map(|chain| (chain.to_string(), description.clone())),
                );
            }
            _ => {}
        }
    }
    out
}

fn pdb_hetero_names(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::<(String, String)>::new();
    for line in text.lines().filter(|line| line.starts_with("HETNAM")) {
        let component = field(line, 11, 14).trim();
        let description = field(line, 15, 80).trim();
        if component.is_empty() {
            continue;
        }
        if let Some((_, existing)) = out.iter_mut().find(|(id, _)| id == component) {
            if !existing.is_empty() && !description.is_empty() {
                existing.push(' ');
            }
            existing.push_str(description);
        } else {
            out.push((component.to_string(), description.to_string()));
        }
    }
    out
}

fn parse_pdb_secondary(text: &str) -> (Vec<SecondaryRange>, Vec<SecondaryRange>) {
    let mut helices = Vec::new();
    let mut sheets = Vec::new();
    for line in text.lines() {
        if line.starts_with("HELIX") {
            if let Some(range) = parse_pdb_helix_range(line) {
                helices.push(range);
            }
        } else if line.starts_with("SHEET") {
            if let Some(range) = parse_pdb_sheet_range(line) {
                sheets.push(range);
            }
        }
    }
    (helices, sheets)
}

fn normalize_pdb_secondary_ranges(ranges: &mut [SecondaryRange], atoms: &[Atom]) {
    for range in ranges {
        let start = pdb_label_secondary_boundary(
            atoms,
            &range.chain,
            range.start,
            &range.start_insertion_code,
        );
        let end =
            pdb_label_secondary_boundary(atoms, &range.chain, range.end, &range.end_insertion_code);
        let (Some((start_chain, start_seq_id)), Some((end_chain, end_seq_id))) = (start, end)
        else {
            continue;
        };
        if start_chain != end_chain {
            continue;
        }
        range.chain = start_chain;
        range.start = start_seq_id;
        range.end = end_seq_id;
    }
}

fn pdb_label_secondary_boundary(
    atoms: &[Atom],
    auth_chain: &str,
    auth_seq_id: i32,
    insertion_code: &str,
) -> Option<(String, i32)> {
    atoms.iter().find_map(|atom| {
        (atom.auth_chain == auth_chain
            && pdb_integer(&atom.auth_residue_seq) == auth_seq_id
            && atom.insertion_code == insertion_code)
            .then(|| {
                atom.residue_seq
                    .trim()
                    .parse::<i32>()
                    .ok()
                    .map(|label_seq_id| (atom.chain.clone(), label_seq_id))
            })
            .flatten()
    })
}

fn parse_pdb_helix_range(line: &str) -> Option<SecondaryRange> {
    Some(SecondaryRange {
        chain: field(line, 19, 20).trim().to_string(),
        start: parse_pdb_i32_field(line, 21, 25)?,
        start_insertion_code: field(line, 25, 26).trim().to_string(),
        end: parse_pdb_i32_field(line, 33, 37)?,
        end_insertion_code: field(line, 37, 38).trim().to_string(),
    })
}

fn parse_pdb_sheet_range(line: &str) -> Option<SecondaryRange> {
    Some(SecondaryRange {
        chain: field(line, 21, 22).trim().to_string(),
        start: parse_pdb_i32_field(line, 22, 26)?,
        start_insertion_code: field(line, 26, 27).trim().to_string(),
        end: parse_pdb_i32_field(line, 33, 37)?,
        end_insertion_code: field(line, 37, 38).trim().to_string(),
    })
}

fn parse_pdb_i32_field(line: &str, start: usize, end: usize) -> Option<i32> {
    field(line, start, end).trim().parse::<i32>().ok()
}

fn parse_pdb_charge(value: &str) -> i32 {
    let value = value.trim();
    if value.len() != 2 {
        return 0;
    }
    let mut chars = value.chars();
    let Some(magnitude) = chars.next().and_then(|ch| ch.to_digit(10)) else {
        return 0;
    };
    match chars.next() {
        Some('+') => magnitude as i32,
        Some('-') => -(magnitude as i32),
        _ => 0,
    }
}

fn parse_pdb_assemblies(text: &str) -> Vec<Assembly> {
    let mut assemblies = Vec::new();
    let mut current_id = String::new();
    let mut current_chains = Vec::new();
    let mut current_rows = [[0.0; 4]; 3];
    let mut seen_rows = [false; 3];

    let flush = |assemblies: &mut Vec<Assembly>,
                 current_id: &str,
                 current_chains: &[String],
                 current_rows: [[f32; 4]; 3],
                 seen_rows: [bool; 3]| {
        if current_id.is_empty() || !seen_rows.iter().all(|v| *v) {
            return;
        }
        let transform = Transform { m: current_rows };
        if let Some(existing) = assemblies.iter_mut().find(|a| a.id == current_id) {
            existing.transforms.push(transform);
            existing.generators.push(AssemblyGenerator::from_transforms(
                current_id,
                current_chains.to_vec(),
                existing.transforms.len() - 1,
                vec![transform],
                vec![Vec::new()],
            ));
            for chain in current_chains {
                if !existing.asym_ids.iter().any(|id| id == chain) {
                    existing.asym_ids.push(chain.clone());
                }
            }
        } else {
            assemblies.push(Assembly {
                id: current_id.to_string(),
                details: String::new(),
                oligomeric_details: String::new(),
                oligomeric_count: None,
                asym_ids: current_chains.to_vec(),
                transforms: vec![transform],
                generators: vec![AssemblyGenerator::from_transforms(
                    current_id,
                    current_chains.to_vec(),
                    0,
                    vec![transform],
                    vec![Vec::new()],
                )],
            });
        }
    };

    for line in text.lines() {
        if line.starts_with("REMARK 350 BIOMOLECULE:") {
            flush(
                &mut assemblies,
                &current_id,
                &current_chains,
                current_rows,
                seen_rows,
            );
            current_id = line
                .split_once(':')
                .map(|(_, v)| v.split_whitespace().next().unwrap_or("1").to_string())
                .unwrap_or_else(|| "1".to_string());
            current_chains.clear();
            current_rows = [[0.0; 4]; 3];
            seen_rows = [false; 3];
        } else if line.starts_with("REMARK 350 APPLY THE FOLLOWING TO CHAINS:")
            || line.starts_with("REMARK 350                    AND CHAINS:")
        {
            if let Some((_, chains)) = line.split_once(':') {
                current_chains.extend(
                    chains
                        .split(',')
                        .map(|s| s.trim().trim_end_matches('.'))
                        .filter(|s| !s.is_empty())
                        .map(str::to_string),
                );
            }
        } else if line.starts_with("REMARK 350   BIOMT") {
            let row = field(line, 18, 19).trim().parse::<usize>().unwrap_or(0);
            if (1..=3).contains(&row) {
                let parts: Vec<f32> = line
                    .get(23..)
                    .unwrap_or("")
                    .split_whitespace()
                    .filter_map(|p| parse_js_number_f32(p).ok())
                    .collect();
                if parts.len() >= 4 {
                    current_rows[row - 1] = [parts[0], parts[1], parts[2], parts[3]];
                    seen_rows[row - 1] = true;
                    if seen_rows.iter().all(|v| *v) {
                        flush(
                            &mut assemblies,
                            &current_id,
                            &current_chains,
                            current_rows,
                            seen_rows,
                        );
                        current_rows = [[0.0; 4]; 3];
                        seen_rows = [false; 3];
                    }
                }
            }
        }
    }
    flush(
        &mut assemblies,
        &current_id,
        &current_chains,
        current_rows,
        seen_rows,
    );
    assemblies
}

fn field(line: &str, start: usize, end: usize) -> &str {
    if start >= line.len() {
        return "";
    }
    line.get(start..end.min(line.len())).unwrap_or("")
}

fn parse_f64(value: &str) -> Result<f64, String> {
    value
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("invalid coordinate: {}", value.trim()))
}

fn parse_js_number_f32(value: &str) -> Result<f32, std::num::ParseFloatError> {
    value.trim().parse::<f64>().map(|value| value as f32)
}

fn lookup_serial(pairs: &[(usize, usize)], serial: usize) -> Option<usize> {
    pairs.iter().find_map(|(s, i)| (*s == serial).then_some(*i))
}
