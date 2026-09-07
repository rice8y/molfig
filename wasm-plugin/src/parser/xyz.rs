/**
 * Copyright (c) 2021 mol* contributors, licensed under MIT, See LICENSE file for more info.
 *
 * Rust port of Mol*'s XYZ reader and model-format normalization.
 */
use crate::chemistry::normalize_element;
use crate::model::{
    Atom, AtomSiteColumnPresence, ChemicalComponent, Entity, EntityIndexMap, Molecule, SourceData,
    Vec3,
};

pub(crate) fn parse_xyz(text: &str) -> Result<Molecule, String> {
    let mut lines = text.lines();
    let mut atoms = Vec::new();
    let mut model_num = 0i32;

    loop {
        let Some(count_line) = lines.next() else {
            break;
        };
        let Ok(count) = count_line.trim().parse::<usize>() else {
            break;
        };
        if count == 0 {
            break;
        }

        lines
            .next()
            .ok_or_else(|| format!("XYZ model {} is missing its comment line", model_num + 1))?;

        for atom_index in 0..count {
            let line = lines.next().ok_or_else(|| {
                format!(
                    "XYZ model {} declares {count} atoms but contains only {atom_index}",
                    model_num + 1
                )
            })?;
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() < 4 {
                return Err(format!(
                    "XYZ model {} atom {} must contain an element symbol and three coordinates",
                    model_num + 1,
                    atom_index + 1
                ));
            }
            let x = parse_coordinate(fields[1], model_num, atom_index, "x")?;
            let y = parse_coordinate(fields[2], model_num, atom_index, "y")?;
            let z = parse_coordinate(fields[3], model_num, atom_index, "z")?;
            let type_symbol = fields[0].to_string();
            let element = normalize_element(type_symbol.clone());
            let source_index = atoms.len();

            atoms.push(Atom {
                id: atom_index,
                source_index,
                model_num,
                name: type_symbol.clone(),
                type_symbol,
                element,
                chain: "A".to_string(),
                auth_chain: "A".to_string(),
                entity_id: "1".to_string(),
                residue: "MOL".to_string(),
                auth_residue: "MOL".to_string(),
                group_pdb: String::new(),
                residue_seq: "1".to_string(),
                auth_residue_seq: "1".to_string(),
                insertion_code: String::new(),
                alt_id: String::new(),
                auth_name: fields[0].to_string(),
                occupancy: 1.0,
                b_iso: 0.0,
                formal_charge: 0,
                position: Vec3 { x, y, z },
                het: false,
                operator_name: String::new(),
            });
        }

        model_num += 1;
    }

    if atoms.is_empty() {
        return Err("XYZ input does not contain a molecule".to_string());
    }

    let entities = vec![Entity {
        id: "1".to_string(),
        type_name: "non-polymer".to_string(),
        description: "Unknown Entity".to_string(),
    }];
    let entity_index = EntityIndexMap::from_entities(&entities, &[], &[]);

    Ok(Molecule {
        source_data: SourceData::xyz(),
        atom_site_columns: AtomSiteColumnPresence {
            occupancy_defined: true,
            b_iso_defined: false,
            xyz_defined: true,
        },
        atoms,
        entities,
        entity_index,
        chemical_components: vec![ChemicalComponent {
            id: "MOL".to_string(),
            name: "Unknown Molecule".to_string(),
            type_name: "other".to_string(),
            mon_nstd_flag: "n".to_string(),
            ..ChemicalComponent::default()
        }],
        ..Molecule::default()
    })
}

fn parse_coordinate(
    value: &str,
    model_num: i32,
    atom_index: usize,
    axis: &str,
) -> Result<f32, String> {
    let value = value.parse::<f64>().map_err(|_| {
        format!(
            "invalid XYZ {axis} coordinate for model {} atom {}: {value}",
            model_num + 1,
            atom_index + 1
        )
    })? as f32;
    if !value.is_finite() {
        return Err(format!(
            "invalid XYZ {axis} coordinate for model {} atom {}: {value}",
            model_num + 1,
            atom_index + 1
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_xyz_like_molstar() {
        let molecule = parse_xyz("3\nwater\nO 0 0 0\nH 0.9572 0 0\nH -0.239 0.927 0\n").unwrap();

        assert_eq!(molecule.source_data.kind, "xyz");
        assert_eq!(molecule.atoms.len(), 3);
        assert_eq!(molecule.atoms[0].id, 0);
        assert_eq!(molecule.atoms[0].name, "O");
        assert_eq!(molecule.atoms[0].residue, "MOL");
        assert_eq!(molecule.atoms[0].chain, "A");
        assert_eq!(molecule.atoms[0].entity_id, "1");
        assert_eq!(molecule.entities[0].description, "Unknown Entity");
        assert_eq!(molecule.chemical_components[0].name, "Unknown Molecule");
    }

    #[test]
    fn keeps_xyz_frames_as_models() {
        let molecule = parse_xyz("1\nfirst\nH 0 0 0\n1\nsecond\nH 1 2 3\n").unwrap();

        assert_eq!(molecule.atoms.len(), 2);
        assert_eq!(molecule.atoms[0].model_num, 0);
        assert_eq!(molecule.atoms[1].model_num, 1);
        assert_eq!(molecule.atoms[0].id, 0);
        assert_eq!(molecule.atoms[1].id, 0);
    }
}
