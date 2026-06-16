use crate::chemistry::{infer_element_from_name, normalize_element};
use crate::model::{
    Assembly, AssemblyGenerator, Atom, AtomSiteAnisotrop, AtomSiteColumnPresence, Bond,
    BondMetadata, EntityIndexMap, Molecule, SecondaryRange, SourceData, Transform, Vec3,
};

use super::normalize_type_symbol_molstar;

pub(crate) fn parse_pdb(text: &str) -> Result<Molecule, String> {
    let mut atoms = Vec::new();
    let mut serial_to_index = Vec::<(usize, usize)>::new();
    let mut bonds = Vec::new();
    let mut bond_metadata = Vec::new();
    let mut atom_site_anisotrop = Vec::new();
    let mut b_iso_defined = false;
    let assemblies = parse_pdb_assemblies(text);
    let (helices, sheets) = parse_pdb_secondary(text);
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
            let x = parse_f32(field(line, 30, 38))?;
            let y = parse_f32(field(line, 38, 46))?;
            let z = parse_f32(field(line, 46, 54))?;
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
                het: line.starts_with("HETATM"),
                operator_name: String::new(),
            });
        } else if line.starts_with("ANISOU") {
            if let Some(row) = parse_pdb_anisou(line) {
                atom_site_anisotrop.push(row);
            }
        } else if line.starts_with("CONECT") {
            let Some((a, rest)) = line
                .get(6..)
                .map(str::trim)
                .and_then(|s| s.split_once(char::is_whitespace))
            else {
                continue;
            };
            let Ok(serial_a) = a.trim().parse::<usize>() else {
                continue;
            };
            for part in rest.split_whitespace() {
                if let Ok(serial_b) = part.parse::<usize>() {
                    if let (Some(ia), Some(ib)) = (
                        lookup_serial(&serial_to_index, serial_a),
                        lookup_serial(&serial_to_index, serial_b),
                    ) {
                        if ia < ib {
                            bonds.push(Bond { a: ia, b: ib });
                            bond_metadata
                                .push(BondMetadata::pdb_conect(bond_metadata.len() as i32));
                        }
                    }
                }
            }
        }
    }

    if atoms.is_empty() {
        return Err("no ATOM/HETATM records found in PDB input".to_string());
    }
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
        entities: Vec::new(),
        entity_index: EntityIndexMap::default(),
        entity_polymers: Vec::new(),
        entity_poly_seq: Vec::new(),
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
        chemical_components: Vec::new(),
        chemical_component_atoms: Vec::new(),
        chemical_component_bonds: Vec::new(),
        chemical_component_angles: Vec::new(),
        rings: Vec::new(),
        resonance: Default::default(),
        derived_aromatic_bonds: Default::default(),
        derived_resonance_bonds: Default::default(),
    })
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
    line.get(start..end).unwrap_or("")
}

fn parse_f32(value: &str) -> Result<f32, String> {
    parse_js_number_f32(value).map_err(|_| format!("invalid coordinate: {}", value.trim()))
}

fn parse_js_number_f32(value: &str) -> Result<f32, std::num::ParseFloatError> {
    value.trim().parse::<f64>().map(|value| value as f32)
}

fn lookup_serial(pairs: &[(usize, usize)], serial: usize) -> Option<usize> {
    pairs.iter().find_map(|(s, i)| (*s == serial).then_some(*i))
}
