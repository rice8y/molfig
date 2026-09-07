use super::*;
use crate::model::PrincipalAxes;

#[test]
fn principal_axes_port_matches_molstar_axis_order_and_box_axes() {
    let Some(molstar_principal_axes) =
        read_molstar_source("mol-math/linear-algebra/matrix/principal-axes.ts")
    else {
        eprintln!("skipping pinned Mol* principal-axes source audit; artifacts is absent");
        return;
    };
    assert!(molstar_principal_axes.contains("svd(A, W, U, V);"));
    assert!(molstar_principal_axes.contains("Math.sqrt(W.data[0] / n3)"));
    assert!(molstar_principal_axes.contains("calculateBoxAxes(positions, momentsAxes)"));

    let axes = PrincipalAxes::of_positions(&[
        vec3(-1.0, 0.0, 0.0),
        vec3(1.0, 0.0, 0.0),
        vec3(0.0, -2.0, 0.0),
        vec3(0.0, 2.0, 0.0),
        vec3(0.0, 0.0, -3.0),
        vec3(0.0, 0.0, 3.0),
    ]);
    let _: &Axes3D = &axes.box_axes;

    assert_vec3_close(axes.moments_axes.origin, vec3(0.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(axes.moments_axes.dir_a, vec3(0.0, 0.0, 3.0), 0.000_001);
    assert_vec3_close(axes.moments_axes.dir_b, vec3(0.0, 2.0, 0.0), 0.000_001);
    assert_vec3_close(axes.moments_axes.dir_c, vec3(1.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(axes.box_axes.origin, vec3(0.0, 0.0, 0.0), 0.000_001);
    assert_vec3_close(axes.box_axes.dir_a, vec3(0.0, 0.0, 3.0), 0.000_001);
    assert_vec3_close(axes.box_axes.dir_b, vec3(0.0, 2.0, 0.0), 0.000_001);
    assert_vec3_close(axes.box_axes.dir_c, vec3(1.0, 0.0, 0.0), 0.000_001);
}

#[test]
fn principal_axes_degenerate_cases_keep_molstar_finite_origin_behavior() {
    let same_plane = PrincipalAxes::of_positions(&[
        vec3(0.1945, -0.0219, -0.0416),
        vec3(-0.0219, -0.0219, -0.0119),
    ]);
    assert!(same_plane.box_axes.origin.is_finite());
    assert_vec3_close(
        same_plane.box_axes.origin,
        same_plane.moments_axes.origin,
        0.000_001,
    );

    let same_point = PrincipalAxes::of_positions(&[
        vec3(0.1945, -0.0219, -0.0416),
        vec3(0.1945, -0.0219, -0.0416),
    ]);
    assert!(same_point.box_axes.origin.is_finite());
    assert_vec3_close(
        same_point.box_axes.origin,
        same_point.moments_axes.origin,
        0.000_001,
    );
}

#[test]
fn unit_and_structure_principal_axes_follow_molstar_position_semantics() {
    let cif = b"data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA GLY A 1 0.000 0.000 0.000
ATOM 2 C CB GLY A 1 1.000 0.000 0.000
#
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 2 A
#
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
2 1 0 0 10 0 1 0 0 0 0 1 0
#
";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            assembly: Some("1".to_string()),
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();
    assert_eq!(structure.units.len(), 1);
    assert_vec3_close(
        structure.units[0].props.principal_axes.moments_axes.origin,
        vec3(0.5, 0.0, 0.0),
        0.000_001,
    );
    assert_vec3_close(
        structure.principal_axes.moments_axes.origin,
        vec3(10.5, 0.0, 0.0),
        0.000_001,
    );
    assert_vec3_close(
        structure.principal_axes.box_axes.origin,
        vec3(10.5, 0.0, 0.0),
        0.000_001,
    );
}

#[test]
fn molstar_model_split_preserves_source_order_and_source_indices_are_frame_local() {
    let cif = b"data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.pdbx_PDB_model_num
ATOM 10 C CA CAX . GLY GLX A X 1 2 20 . 2.000 0.000 0.000 0.75 12.25 -1 1
ATOM 11 C CA CAY A ALA ALY A Y 1 1 10 A 1.000 0.000 0.000 0.50 10.50 1 1
ATOM 12 C CA CAZ B SER SRY A Z 1 1 10 . 9.000 0.000 0.000 0.25 99.00 0 2
#
";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();

    assert_eq!(structure.models.len(), 2);
    assert_eq!(structure.model.model_num, 1);
    assert_eq!(structure.model.hierarchy.atoms.len(), 2);
    assert_eq!(structure.model.hierarchy.chains.len(), 1);
    assert_eq!(structure.model.hierarchy.chains[0].auth_id, "X");
    assert_eq!(structure.model.hierarchy.chains[0].id, "A");
    assert_eq!(structure.model.hierarchy.chains[0].entity_id, "1");
    assert_eq!(
        structure
            .model
            .hierarchy
            .atoms
            .iter()
            .map(|atom| (
                atom.name.as_str(),
                atom.auth_name.as_str(),
                atom.alt_id.as_str()
            ))
            .collect::<Vec<_>>(),
        vec![("CA", "CAX", ""), ("CA", "CAY", "A")]
    );
    assert_eq!(
        structure
            .model
            .hierarchy
            .residues
            .iter()
            .map(|residue| {
                (
                    residue.comp_id.as_str(),
                    residue.auth_comp_id.as_str(),
                    residue.label_seq_id.as_str(),
                    residue.auth_seq_id.as_str(),
                    residue.insertion_code.as_str(),
                )
            })
            .collect::<Vec<_>>(),
        vec![
            ("GLY", "GLX", "2", "20", ""),
            ("ALA", "ALY", "1", "10", "A")
        ]
    );
    assert_eq!(structure.model.conformation.atom_ids, vec![10, 11]);
    assert_eq!(
        structure.model.conformation.positions,
        vec![vec3(2.0, 0.0, 0.0), vec3(1.0, 0.0, 0.0)]
    );
    assert_eq!(structure.model.conformation.x, vec![2.0, 1.0]);
    assert_eq!(structure.model.conformation.y, vec![0.0, 0.0]);
    assert_eq!(structure.model.conformation.z, vec![0.0, 0.0]);
    assert_eq!(structure.model.conformation.occupancies, vec![0.75, 0.5]);
    assert_eq!(structure.model.conformation.b_iso, vec![12.25, 10.5]);
    assert_eq!(structure.model.conformation.formal_charges, vec![-1, 1]);
    assert_eq!(structure.model.hierarchy.atom_source_index, vec![0, 1]);
    assert_eq!(structure.model.hierarchy.residue_source_index, vec![0, 1]);
    assert_eq!(
        structure.model.hierarchy.residue_atom_segments.offsets,
        vec![0, 1, 2]
    );
    assert_eq!(
        structure.model.hierarchy.residue_atom_segments.index,
        vec![0, 1]
    );
    assert_eq!(
        structure.model.hierarchy.chain_atom_segments.offsets,
        vec![0, 2]
    );
    assert_eq!(
        structure.model.hierarchy.chain_atom_segments.index,
        vec![0, 0]
    );
    assert_eq!(structure.properties.atom_source_index, vec![0, 1]);
    assert_eq!(structure.properties.atom_id, vec![10, 11]);
    assert_eq!(structure.properties.label_atom_id, vec!["CA", "CA"]);
    assert_eq!(structure.properties.auth_atom_id, vec!["CAX", "CAY"]);
    assert_eq!(structure.properties.label_alt_id, vec!["", "A"]);
    assert_eq!(structure.properties.label_comp_id, vec!["GLY", "ALA"]);
    assert_eq!(structure.properties.auth_comp_id, vec!["GLX", "ALY"]);
    assert_eq!(structure.properties.label_seq_id, vec!["2", "1"]);
    assert_eq!(structure.properties.auth_seq_id, vec!["20", "10"]);
    assert_eq!(structure.properties.pdbx_pdb_ins_code, vec!["", "A"]);
    assert_eq!(structure.properties.label_asym_id, vec!["A", "A"]);
    assert_eq!(structure.properties.auth_asym_id, vec!["X", "X"]);
    assert_eq!(structure.properties.label_entity_id, vec!["1", "1"]);

    assert_eq!(structure.models[0].hierarchy.atom_source_index, vec![0, 1]);
    assert_eq!(structure.models[1].model_num, 2);
    assert_eq!(structure.models[1].hierarchy.atom_source_index, vec![2]);
    assert_eq!(structure.models[1].conformation.atom_ids, vec![12]);
    assert_eq!(
        structure.models[1].conformation.positions,
        vec![vec3(9.0, 0.0, 0.0)]
    );
    assert_eq!(structure.models[1].conformation.x, vec![9.0]);
    assert_eq!(structure.models[1].conformation.y, vec![0.0]);
    assert_eq!(structure.models[1].conformation.z, vec![0.0]);
    assert_eq!(structure.models[1].conformation.occupancies, vec![0.25]);
    assert_eq!(structure.models[1].conformation.b_iso, vec![99.0]);
}

#[test]
fn multi_model_atomic_fixture_splits_models_from_atom_site_model_numbers() {
    let molecule = parse_molecule_with_options(
        include_bytes!("../../tests/fixtures/cif/multi-model-atomic.cif"),
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.atoms.len(), 4);
    assert!(molecule
        .source_data
        .categories
        .iter()
        .any(|category| category.name == "atom_site" && category.row_count == 4));

    let structure = molecule.atomic_structure();
    assert_eq!(structure.models.len(), 2);
    assert_eq!(
        structure
            .models
            .iter()
            .map(|model| model.model_num)
            .collect::<Vec<_>>(),
        vec![1, 2]
    );
    assert_eq!(structure.model.conformation.atom_ids, vec![1, 2]);
    assert_eq!(structure.model.hierarchy.atom_source_index, vec![0, 1]);
    assert_eq!(structure.models[1].conformation.atom_ids, vec![3, 4]);
    assert_eq!(structure.models[1].hierarchy.atom_source_index, vec![2, 3]);
    assert_eq!(
        structure.models[1].conformation.positions,
        vec![vec3(0.1, 0.2, 0.0), vec3(1.55, 0.2, 0.0)]
    );
    assert_eq!(structure.models[1].conformation.occupancies, vec![0.9, 0.9]);
    assert_eq!(structure.models[1].conformation.b_iso, vec![20.0, 21.0]);
}

#[test]
fn covalent_cross_link_fixture_uses_struct_conn_metadata_and_order_fallback() {
    let molecule = parse_molecule_with_options(
        include_bytes!("../../tests/fixtures/cif/covalent-cross-link.cif"),
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.atoms.len(), 2);
    assert_eq!(molecule.bonds.len(), 1);
    assert_eq!(molecule.bond_metadata.len(), 1);
    assert_eq!(molecule.atoms[molecule.bonds[0].a].name, "NZ");
    assert_eq!(molecule.atoms[molecule.bonds[0].b].name, "C15");

    let metadata = &molecule.bond_metadata[0];
    assert_eq!(metadata.source, BondSource::StructConn);
    assert_eq!(metadata.order, 2);
    assert_eq!(metadata.distance, Some(1.30));
    assert!(metadata.flags.contains(BondFlags::COVALENT));
    assert!(!metadata.flags.contains(BondFlags::COMPUTED));

    let struct_conn = metadata.struct_conn.as_ref().unwrap();
    assert_eq!(struct_conn.id, "lys-ret");
    assert_eq!(struct_conn.conn_type_id, "covale");
    assert_eq!(struct_conn.partner_a_comp_id, "LYS");
    assert_eq!(struct_conn.partner_b_comp_id, "RET");
    assert_eq!(struct_conn.partner_a_atom_index, 0);
    assert_eq!(struct_conn.partner_b_atom_index, 1);
}

#[test]
fn ihm_model_list_splits_atomic_and_coarse_models_per_molstar_read_integrative() {
    let cif = b"data_demo
loop_
_ihm_model_list.model_id
_ihm_model_list.model_name
_ihm_model_list.assembly_id
_ihm_model_list.protocol_id
_ihm_model_list.representation_id
101 one 1 1 1
102 two 1 1 1
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.ihm_model_id
ATOM 1 C CA GLY A 1 1 1.000 0.000 0.000 101
ATOM 2 C CA GLY A 1 1 2.000 0.000 0.000 102
#
loop_
_ihm_sphere_obj_site.id
_ihm_sphere_obj_site.model_id
_ihm_sphere_obj_site.entity_id
_ihm_sphere_obj_site.asym_id
_ihm_sphere_obj_site.seq_id_begin
_ihm_sphere_obj_site.seq_id_end
_ihm_sphere_obj_site.Cartn_x
_ihm_sphere_obj_site.Cartn_y
_ihm_sphere_obj_site.Cartn_z
_ihm_sphere_obj_site.object_radius
1 101 1 A 1 10 10.0 0.0 0.0 1.0
2 102 1 A 1 10 20.0 0.0 0.0 2.0
#
loop_
_ihm_gaussian_obj_site.id
_ihm_gaussian_obj_site.model_id
_ihm_gaussian_obj_site.entity_id
_ihm_gaussian_obj_site.asym_id
_ihm_gaussian_obj_site.seq_id_begin
_ihm_gaussian_obj_site.seq_id_end
_ihm_gaussian_obj_site.mean_Cartn_x
_ihm_gaussian_obj_site.mean_Cartn_y
_ihm_gaussian_obj_site.mean_Cartn_z
_ihm_gaussian_obj_site.weight
_ihm_gaussian_obj_site.covariance_matrix[1][1]
_ihm_gaussian_obj_site.covariance_matrix[1][2]
_ihm_gaussian_obj_site.covariance_matrix[1][3]
_ihm_gaussian_obj_site.covariance_matrix[2][1]
_ihm_gaussian_obj_site.covariance_matrix[2][2]
_ihm_gaussian_obj_site.covariance_matrix[2][3]
_ihm_gaussian_obj_site.covariance_matrix[3][1]
_ihm_gaussian_obj_site.covariance_matrix[3][2]
_ihm_gaussian_obj_site.covariance_matrix[3][3]
1 101 1 A 11 20 30.0 0.0 0.0 1.0 1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0
2 102 1 A 11 20 40.0 0.0 0.0 2.0 4.0 0.0 0.0 0.0 4.0 0.0 0.0 0.0 4.0
#
";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();

    assert_eq!(structure.models.len(), 2);
    assert_eq!(structure.coarse_models.len(), 2);
    assert_eq!(
        structure
            .models
            .iter()
            .map(|model| model.model_num)
            .collect::<Vec<_>>(),
        vec![101, 102]
    );
    assert_eq!(structure.models[0].conformation.atom_ids, vec![1]);
    assert_eq!(structure.models[1].conformation.atom_ids, vec![2]);
    assert_eq!(structure.coarse.conformation.spheres.len(), 1);
    assert_eq!(structure.coarse.conformation.gaussians.len(), 1);
    assert_eq!(
        structure.coarse_models[0].conformation.spheres.position(0),
        Some(vec3(10.0, 0.0, 0.0))
    );
    assert_eq!(
        structure.coarse_models[1].conformation.spheres.position(0),
        Some(vec3(20.0, 0.0, 0.0))
    );
    assert_eq!(
        structure.coarse_models[0]
            .conformation
            .gaussians
            .position(0),
        Some(vec3(30.0, 0.0, 0.0))
    );
    assert_eq!(
        structure.coarse_models[1]
            .conformation
            .gaussians
            .position(0),
        Some(vec3(40.0, 0.0, 0.0))
    );
    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| unit.kind)
            .collect::<Vec<_>>(),
        vec![UnitKind::Atomic, UnitKind::Spheres, UnitKind::Gaussians]
    );
    assert_eq!(structure.position(1, 0).unwrap(), vec3(10.0, 0.0, 0.0));
    assert_eq!(structure.position(2, 0).unwrap(), vec3(30.0, 0.0, 0.0));
}

#[test]
fn mixed_atomic_coarse_ihm_fixture_splits_each_kind_by_ihm_model_id() {
    let molecule = parse_molecule_with_options(
        include_bytes!("../../tests/fixtures/cif/mixed-atomic-coarse-ihm.cif"),
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    assert_eq!(molecule.atoms.len(), 2);
    assert_eq!(molecule.coarse_spheres.len(), 2);
    assert_eq!(molecule.coarse_gaussians.len(), 2);
    assert_eq!(
        molecule
            .source_data
            .categories
            .iter()
            .filter(|category| matches!(
                category.name.as_str(),
                "atom_site" | "ihm_model_list" | "ihm_sphere_obj_site" | "ihm_gaussian_obj_site"
            ))
            .map(|category| (category.name.as_str(), category.row_count))
            .collect::<Vec<_>>(),
        vec![
            ("ihm_model_list", 2),
            ("atom_site", 2),
            ("ihm_sphere_obj_site", 2),
            ("ihm_gaussian_obj_site", 2),
        ]
    );

    let structure = molecule.atomic_structure();
    assert_eq!(
        structure
            .models
            .iter()
            .map(|model| model.model_num)
            .collect::<Vec<_>>(),
        vec![101, 102]
    );
    assert_eq!(structure.coarse_models.len(), 2);
    assert_eq!(structure.models[0].conformation.atom_ids, vec![1]);
    assert_eq!(structure.models[1].conformation.atom_ids, vec![2]);
    assert_eq!(
        structure.coarse_models[0].conformation.spheres.position(0),
        Some(vec3(4.0, 0.0, 0.0))
    );
    assert_eq!(
        structure.coarse_models[1].conformation.spheres.position(0),
        Some(vec3(14.0, 0.0, 0.0))
    );
    assert_eq!(
        structure.coarse_models[0]
            .conformation
            .gaussians
            .position(0),
        Some(vec3(8.0, 0.0, 0.0))
    );
    assert_eq!(
        structure.coarse_models[1]
            .conformation
            .gaussians
            .position(0),
        Some(vec3(18.0, 0.0, 0.0))
    );
    assert_eq!(
        structure
            .units
            .iter()
            .map(|unit| unit.kind)
            .collect::<Vec<_>>(),
        vec![UnitKind::Atomic, UnitKind::Spheres, UnitKind::Gaussians]
    );
    assert_eq!(structure.position(0, 0).unwrap(), vec3(-2.0, 0.0, 0.0));
    assert_eq!(structure.position(1, 0).unwrap(), vec3(4.0, 0.0, 0.0));
    assert_eq!(structure.position(2, 0).unwrap(), vec3(8.0, 0.0, 0.0));
}

#[test]
fn molstar_atomic_numbers_treat_isotopes_and_late_elements_like_reference_table() {
    let cif = b"data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 D D1 HOH A 1 0.000 0.000 0.000
ATOM 2 T T1 HOH A 1 1.000 0.000 0.000
ATOM 3 MT MT1 UNK A 2 2.000 0.000 0.000
ATOM 4 DS DS1 UNK A 3 3.000 0.000 0.000
ATOM 5 OG OG1 UNK A 4 4.000 0.000 0.000
#
";
    let molecule = parse_molecule_with_options(
        cif,
        &MeshOptions {
            format: InputFormat::Cif,
            infer_bonds: false,
            ..MeshOptions::default()
        },
    )
    .unwrap();

    let structure = molecule.atomic_structure();

    assert_eq!(
        structure.model.hierarchy.derived.atom.atomic_number,
        vec![1, 1, 109, 0, 0]
    );
}
