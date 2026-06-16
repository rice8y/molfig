use super::{
    AssemblyGenerator, AtomicStructure, Bond, BondMetadata, CoarseGaussian, CoarseSphere, Molecule,
};

pub(super) fn expanded_bonds_from_units(
    expanded_indices: &[(usize, usize)],
    source_bonds: &[Bond],
    source_metadata: &[BondMetadata],
    structure: &AtomicStructure,
) -> (Vec<Bond>, Vec<BondMetadata>) {
    let mut bonds = Vec::new();
    let mut metadata = Vec::new();
    for (source_bond_index, source_bond) in source_bonds.iter().enumerate() {
        let a_positions = expanded_indices
            .iter()
            .enumerate()
            .filter(|(_, (_, source_index))| *source_index == source_bond.a)
            .collect::<Vec<_>>();
        let b_positions = expanded_indices
            .iter()
            .enumerate()
            .filter(|(_, (_, source_index))| *source_index == source_bond.b)
            .collect::<Vec<_>>();
        for (a_expanded, (a_unit, a_source)) in &a_positions {
            for (b_expanded, (b_unit, b_source)) in &b_positions {
                if a_unit == b_unit
                    || structure.inter_unit_bond_exists_exact(
                        *a_unit,
                        structure.unit_element_index(*a_unit, *a_source),
                        *b_unit,
                        structure.unit_element_index(*b_unit, *b_source),
                        source_bond_index,
                    )
                {
                    bonds.push(Bond {
                        a: *a_expanded,
                        b: *b_expanded,
                    });
                    metadata.push(
                        source_metadata
                            .get(source_bond_index)
                            .cloned()
                            .unwrap_or_default(),
                    );
                }
            }
        }
    }
    (bonds, metadata)
}

pub(super) fn expanded_coarse_for_geometry(
    molecule: &Molecule,
) -> (Vec<CoarseSphere>, Vec<CoarseGaussian>) {
    let Some(assembly) = molecule.selected_assembly.as_ref() else {
        return (
            molecule.coarse_spheres.clone(),
            molecule.coarse_gaussians.clone(),
        );
    };
    let generators = if assembly.generators.is_empty() {
        vec![AssemblyGenerator::from_transforms(
            &assembly.id,
            assembly.asym_ids.clone(),
            0,
            assembly.transforms.clone(),
            vec![Vec::new(); assembly.transforms.len()],
        )]
    } else {
        assembly.generators.clone()
    };
    let mut spheres = Vec::new();
    let mut gaussians = Vec::new();
    let mut operator_offset = 0usize;
    for generator in &generators {
        let operators = generator.operators_for_assembly(&assembly.id, operator_offset);
        operator_offset += operators.len();
        for operator in operators {
            for sphere in molecule.coarse_spheres.iter().filter(|sphere| {
                generator.asym_ids.is_empty() || generator.asym_ids.contains(&sphere.asym_id)
            }) {
                let mut expanded = sphere.clone();
                expanded.position = operator.transform.apply(sphere.position);
                spheres.push(expanded);
            }
            for gaussian in molecule.coarse_gaussians.iter().filter(|gaussian| {
                generator.asym_ids.is_empty() || generator.asym_ids.contains(&gaussian.asym_id)
            }) {
                let mut expanded = gaussian.clone();
                expanded.position = operator.transform.apply(gaussian.position);
                gaussians.push(expanded);
            }
        }
    }
    (spheres, gaussians)
}
