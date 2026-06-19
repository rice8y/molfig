use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use crate::chemistry::{vdw_radius, vdw_radius64};
use crate::json::{json_escape, json_string_array};
use crate::model::{
    get_saccharide_shape, AtomicStructure, BondFlags, BondSource, Boundary, BoundingSphere, Face,
    Mesh, MeshMaterial, Molecule, MoleculeType, NucleotideAtoms, NucleotideBaseKind, PolymerType,
    PrincipalAxes, SaccharideShape, SecondaryRange, SecondaryStructureType, StructureUnit,
    TraceResidue, Transform, UnitKind, UnitOperator, Vec3,
};
use crate::options::{ColorTheme, MeshOptions, PolymerProfile, Representation};

mod geometry;

#[cfg(test)]
use crate::model::{CoarseElementKind, CoarseElements, CoarseRange};

#[allow(unused_imports)]
pub(crate) use geometry::add_oriented_ribbon;
pub(crate) use geometry::DVec3;
use geometry::{
    add_curve_segment_ribbon, add_curve_segment_sheet, add_curve_segment_tube,
    add_dashed_tube_path_cached, add_dashed_tube_samples_cached, add_ellipsoid,
    add_fixed_count_dashed_cylinder_cached, add_molstar_buffered_open_cylinder_cached,
    add_molstar_buffered_open_cylinder_with_radius64_cached, add_molstar_cylinder_caps_cached,
    add_open_cylinder_cached, add_oriented_ribbon_with_profile, add_ribbon, add_sheet, add_sphere,
    add_sphere_with_radius64, add_tube_path, fallback_side, helix_trace,
    molstar_cylinder_mesh_counts, molstar_sphere_mesh_counts, molstar_sphere_triangle_count,
    sample_path, sample_path_point_count, CurveSegmentScratch, CylinderPrimitiveCache,
    MolstarLocalTransform, MolstarPrimitiveTransform,
};
#[cfg(test)]
pub(crate) use geometry::{
    add_profile_tube_for_test, add_ribbon as add_ribbon_for_test, add_sheet as add_sheet_for_test,
    add_tube_path as add_tube_path_for_test, interpolate_curve_segment, interpolate_sizes,
    CurveSegmentControls, CurveSegmentState, TestTubeProfile,
};

pub(crate) fn build_mesh(molecule: &Molecule, options: &MeshOptions) -> Mesh {
    build_mesh_with_visible_bounding_sphere(molecule, options).0
}

pub(crate) fn render_materials(molecule: &Molecule, options: &MeshOptions) -> Vec<MeshMaterial> {
    let geometry = molecule.expanded_for_geometry();
    let options = resolved_mesh_options(&geometry, options);
    let structure = geometry.atomic_structure();
    let objects = build_semantic_render_objects_resolved_limited(
        &geometry,
        &options,
        None,
        Some(&structure),
        |_| {},
    );
    let cylinder_radial_segments = molstar_export_cylinder_radial_segments(
        objects
            .iter()
            .map(|object| render_object_export_cylinder_count(&object.object))
            .sum(),
    );
    let mut materials = Vec::new();
    for object in &objects {
        if object
            .object
            .mesh_estimate(&options, cylinder_radial_segments)
            .faces
            == 0
        {
            continue;
        }
        let Some(material) = object.material else {
            continue;
        };
        if !materials.contains(&material) {
            materials.push(material);
        }
    }
    materials
}

pub(crate) fn build_mesh_with_visible_bounding_sphere(
    molecule: &Molecule,
    options: &MeshOptions,
) -> (Mesh, Option<BoundingSphere>) {
    let (mesh, sphere, _) =
        build_mesh_with_visible_bounding_sphere_and_operator_snapshot(molecule, options, false);
    (mesh, sphere)
}

pub(crate) fn build_mesh_with_visible_bounding_sphere_and_operator_snapshot(
    molecule: &Molecule,
    options: &MeshOptions,
    capture_operators: bool,
) -> (Mesh, Option<BoundingSphere>, Vec<UnitOperator>) {
    let expansion = molecule.expanded_for_geometry_with_operator_snapshot(capture_operators);
    let geometry = expansion.molecule;
    let options = resolved_mesh_options(&geometry, options);
    let structure = geometry.atomic_structure();
    let objects = build_semantic_render_objects_resolved_limited(
        &geometry,
        &options,
        None,
        Some(&structure),
        |_| {},
    );
    let structure_sphere =
        molstar_visible_renderable_bounding_sphere_with_structure(&geometry, &options, &structure);
    let (mesh, mesh_slice_sphere, _) =
        flatten_semantic_render_objects_with_visible_bounding_sphere_and_stats(
            &objects,
            &geometry,
            &options,
            structure_sphere.is_none(),
        );
    let visible_bounding_sphere = structure_sphere.or(mesh_slice_sphere);
    (mesh, visible_bounding_sphere, expansion.assembly_operators)
}

pub(crate) struct RenderSummaries {
    pub(crate) render_objects_json: String,
    pub(crate) representation_json: String,
    pub(crate) structure: AtomicStructure,
    pub(crate) geometry: GeometryInfoSnapshot,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct GeometryInfoSnapshot {
    pub(crate) atom_count: usize,
    pub(crate) coarse_sphere_count: usize,
    pub(crate) coarse_gaussian_count: usize,
    pub(crate) bond_count: usize,
    pub(crate) bond_metadata: BondMetadataSnapshot,
    pub(crate) bounds: Option<(Vec3, Vec3)>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct BondMetadataSnapshot {
    pub(crate) count: usize,
    pub(crate) computed: usize,
    pub(crate) pdb_conect: usize,
    pub(crate) struct_conn: usize,
    pub(crate) index_pair: usize,
    pub(crate) chem_comp: usize,
    pub(crate) covalent: usize,
    pub(crate) metallic_coordination: usize,
    pub(crate) hydrogen_bond: usize,
    pub(crate) disulfide: usize,
    pub(crate) aromatic: usize,
    pub(crate) computed_flag: usize,
    pub(crate) resonance: usize,
    pub(crate) rings: usize,
    pub(crate) aromatic_rings: usize,
    pub(crate) delocalized_bonds: usize,
}

pub(crate) struct RenderScene {
    pub(crate) mesh: Mesh,
    pub(crate) visible_bounding_sphere: Option<BoundingSphere>,
    pub(crate) summaries: RenderSummaries,
    pub(crate) assembly_operators: Vec<UnitOperator>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct RenderObjectMeshStats {
    draw_count: usize,
    vertex_count: usize,
    group_count: usize,
}

pub(crate) fn render_summaries_json(molecule: &Molecule, options: &MeshOptions) -> RenderSummaries {
    let geometry = molecule.expanded_for_geometry();
    let options = resolved_mesh_options(&geometry, options);
    let structure = geometry.atomic_structure();
    let objects = build_semantic_render_objects_resolved_limited(
        &geometry,
        &options,
        None,
        Some(&structure),
        |_| {},
    );
    let object_stats = render_object_mesh_stats_from_estimates(&objects, &options);
    render_summaries_json_from_resolved(
        &options,
        structure,
        GeometryInfoSnapshot::from_molecule(&geometry),
        &objects,
        &object_stats,
    )
}

pub(crate) fn build_render_scene_with_summaries(
    molecule: &Molecule,
    options: &MeshOptions,
) -> RenderScene {
    let expansion =
        molecule.expanded_for_geometry_with_operator_snapshot(options.include_operator_metadata);
    let geometry = expansion.molecule;
    let options = resolved_mesh_options(&geometry, options);
    let structure = geometry.atomic_structure();
    let objects = build_semantic_render_objects_resolved_limited(
        &geometry,
        &options,
        None,
        Some(&structure),
        |_| {},
    );
    let structure_sphere =
        molstar_visible_renderable_bounding_sphere_with_structure(&geometry, &options, &structure);
    let (mesh, mesh_slice_sphere, object_stats) =
        flatten_semantic_render_objects_with_visible_bounding_sphere_and_stats(
            &objects,
            &geometry,
            &options,
            structure_sphere.is_none(),
        );
    let visible_bounding_sphere = structure_sphere.or(mesh_slice_sphere);
    let summaries = render_summaries_json_from_resolved(
        &options,
        structure,
        GeometryInfoSnapshot::from_molecule(&geometry),
        &objects,
        &object_stats,
    );
    RenderScene {
        mesh,
        visible_bounding_sphere,
        summaries,
        assembly_operators: expansion.assembly_operators,
    }
}

pub(crate) fn visible_renderable_bounding_sphere_for_export_with_structure(
    molecule: &Molecule,
    options: &MeshOptions,
    structure: &AtomicStructure,
) -> Option<BoundingSphere> {
    let options = resolved_mesh_options(molecule, options);
    molstar_visible_renderable_bounding_sphere_with_structure(molecule, &options, structure)
}

pub(crate) fn visible_renderable_bounding_sphere_report_for_export_with_structure(
    molecule: &Molecule,
    options: &MeshOptions,
    structure: &AtomicStructure,
) -> String {
    let options = resolved_mesh_options(molecule, options);
    let components =
        molstar_visible_renderable_component_spheres_with_structure(molecule, &options, structure);
    let scene = (!components.is_empty()).then(|| {
        let spheres = components
            .iter()
            .map(|(_, sphere)| sphere.clone())
            .collect::<Vec<_>>();
        Boundary::from_bounding_spheres(&spheres).sphere
    });
    let components_json = components
        .iter()
        .map(|(label, sphere)| {
            format!(
                "{{\"label\":\"{}\",\"center\":{},\"radius\":{:.17},\"center64\":{},\"radius64\":{:.17},\"extrema_count\":{}}}",
                json_escape(label),
                sphere_report_vec3_json(sphere.center),
                sphere.radius,
                sphere_report_vec3_64_json(sphere.center64()),
                sphere.radius64(),
                sphere.extrema.len().max(sphere.extrema64.len())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let scene_json = scene
        .as_ref()
        .map(|sphere| {
            format!(
                "{{\"center\":{},\"radius\":{:.17},\"center64\":{},\"radius64\":{:.17},\"extrema_count\":{}}}",
                sphere_report_vec3_json(sphere.center),
                sphere.radius,
                sphere_report_vec3_64_json(sphere.center64()),
                sphere.radius64(),
                sphere.extrema.len().max(sphere.extrema64.len())
            )
        })
        .unwrap_or_else(|| "null".to_string());
    format!(
        "{{\"component_count\":{},\"components\":[{}],\"scene\":{}}}",
        components.len(),
        components_json,
        scene_json
    )
}

fn resolved_mesh_options(molecule: &Molecule, options: &MeshOptions) -> MeshOptions {
    let coarse_count = molecule.coarse_spheres.len() + molecule.coarse_gaussians.len();
    let element_count = molecule.atoms.len() + coarse_count;
    options.resolved_for_quality(element_count, molecule.atoms.is_empty() && coarse_count > 0)
}

#[derive(Clone, Debug)]
pub(crate) enum RenderObject {
    Sphere {
        center: Vec3,
        radius: f64,
    },
    Cylinder {
        start: Vec3,
        end: Vec3,
        radius: f32,
    },
    LinkCylinder {
        start: Vec3,
        end: Vec3,
        radius: f32,
    },
    LinkCylinderWithSegments {
        start: Vec3,
        end: Vec3,
        radius: f64,
        radial_segments: usize,
    },
    Tube {
        points: Vec<Vec3>,
        radius: f32,
    },
    DashedTube {
        points: Vec<Vec3>,
        radius: f32,
    },
    FixedCountDashedCylinder {
        start: Vec3,
        end: Vec3,
        radius: f32,
        length_scale: f32,
        segment_count: usize,
    },
    #[allow(dead_code)]
    Ribbon {
        points: Vec<Vec3>,
        width: f32,
        thickness: f32,
    },
    Sheet {
        points: Vec<Vec3>,
        width: f32,
        thickness: f32,
        arrow_height: f32,
        start_cap: bool,
        end_cap: bool,
    },
    OrientedRibbon {
        centers: Vec<Vec3>,
        normals: Vec<Vec3>,
        width: f32,
        thickness: f32,
        profile: PolymerProfile,
        start_cap: bool,
        end_cap: bool,
        round_cap: bool,
    },
    PolymerTraceSegment {
        controls: geometry::CurveSegmentControls,
        widths: [f32; 3],
        heights: [f32; 3],
        tension: f64,
        shift: f64,
        overhang_width: f32,
        kind: PolymerTraceSegmentKind,
        start_cap: bool,
        end_cap: bool,
        initial: bool,
        final_residue: bool,
        swap_normal_binormal: bool,
    },
    NucleotideRing {
        center: Vec3,
        normal: Vec3,
        radius: f32,
        base: Option<NucleotideRingBase>,
        detail: usize,
        radial_segments: usize,
    },
    NucleotideBlock {
        geometry: NucleotideBlockGeometry,
        radius: f32,
        width: f32,
        depth: f32,
        radial_segments: usize,
    },
    DirectionWedge {
        center: Vec3,
        tangent: Vec3,
        up: Vec3,
        size: f32,
    },
    CarbohydrateSymbol {
        center: Vec3,
        normal: Vec3,
        direction: Vec3,
        shape: SaccharideShape,
        part: CarbohydrateSymbolPart,
    },
    Ellipsoid {
        center: Vec3,
        axes: [Vec3; 3],
    },
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
struct RenderObjectMeshEstimate {
    vertices: usize,
    faces: usize,
}

struct RenderObjectMeshPlan {
    estimate: RenderObjectMeshEstimate,
    dashed_samples: Option<Vec<Vec3>>,
}

impl RenderObjectMeshEstimate {
    fn from_counts((vertices, faces): (usize, usize)) -> Self {
        Self { vertices, faces }
    }

    fn add(self, other: Self) -> Self {
        Self {
            vertices: self.vertices.saturating_add(other.vertices),
            faces: self.faces.saturating_add(other.faces),
        }
    }

    fn scale(self, count: usize) -> Self {
        Self {
            vertices: self.vertices.saturating_mul(count),
            faces: self.faces.saturating_mul(count),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum CarbohydrateSymbolPart {
    Whole,
    Primary,
    Secondary,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum PolymerTraceSegmentKind {
    Ribbon {
        arrow_height: f32,
        swap_width_height: bool,
    },
    Tube {
        profile: PolymerProfile,
        round_cap: bool,
    },
    Sheet {
        arrow_height: f32,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum NucleotideRingBase {
    PurineConnector {
        trace: Vec3,
        n9: Vec3,
    },
    Purine {
        trace: Vec3,
        n1: Vec3,
        c2: Vec3,
        n3: Vec3,
        c4: Vec3,
        c5: Vec3,
        c6: Vec3,
        n7: Vec3,
        c8: Vec3,
        n9: Vec3,
    },
    PyrimidineConnector {
        trace: Vec3,
        n1: Vec3,
    },
    Pyrimidine {
        trace: Vec3,
        n1: Vec3,
        c2: Vec3,
        n3: Vec3,
        c4: Vec3,
        c5: Vec3,
        c6: Vec3,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct NucleotideBlockGeometry {
    trace: Vec3,
    anchor: Vec3,
    block: Option<NucleotideBlockBox>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct NucleotideBlockBox {
    p1: Vec3,
    p2: Vec3,
    p3: Vec3,
    p4: Vec3,
    height: f32,
}

#[derive(Clone, Debug)]
pub(crate) struct SemanticRenderObject {
    pub(crate) geometry_type: &'static str,
    pub(crate) visual: &'static str,
    pub(crate) representation: &'static str,
    pub(crate) secondary_type: &'static str,
    pub(crate) chain: Option<String>,
    pub(crate) residue_start: Option<i32>,
    pub(crate) residue_end: Option<i32>,
    pub(crate) group_id: usize,
    pub(crate) material: Option<MeshMaterial>,
    pub(crate) initial: bool,
    pub(crate) final_residue: bool,
    pub(crate) sec_struc_first: bool,
    pub(crate) sec_struc_last: bool,
    pub(crate) object: RenderObject,
}

#[derive(Clone, Copy, Debug, Default)]
struct TraceFlags {
    initial: bool,
    final_residue: bool,
    sec_struc_first: bool,
    sec_struc_last: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SecondaryTraceKind {
    Helix,
    Sheet,
}

struct SemanticMeta<'a> {
    representation: &'static str,
    secondary_type: &'static str,
    visual: Option<&'static str>,
    chain: Option<&'a str>,
    residue_start: Option<i32>,
    residue_end: Option<i32>,
    trace_flags: TraceFlags,
    material: Option<MeshMaterial>,
}

impl<'a> SemanticMeta<'a> {
    fn new(
        representation: &'static str,
        secondary_type: &'static str,
        chain: Option<&'a str>,
        residue_start: Option<i32>,
        residue_end: Option<i32>,
    ) -> Self {
        Self {
            representation,
            secondary_type,
            visual: None,
            chain,
            residue_start,
            residue_end,
            trace_flags: TraceFlags::default(),
            material: None,
        }
    }

    fn with_trace_flags(mut self, trace_flags: TraceFlags) -> Self {
        self.trace_flags = trace_flags;
        self
    }

    fn with_visual(mut self, visual: &'static str) -> Self {
        self.visual = Some(visual);
        self
    }

    fn with_material(mut self, material: MeshMaterial) -> Self {
        self.material = Some(material);
        self
    }
}

#[allow(dead_code)]
pub(crate) fn build_render_objects(
    molecule: &Molecule,
    options: &MeshOptions,
) -> Vec<RenderObject> {
    build_semantic_render_objects(molecule, options)
        .into_iter()
        .map(|object| object.object)
        .collect()
}

#[cfg(test)]
pub(crate) fn render_object_summary_json(molecule: &Molecule, options: &MeshOptions) -> String {
    let geometry = molecule.expanded_for_geometry();
    let options = resolved_mesh_options(&geometry, options);
    let structure = geometry.atomic_structure();
    let objects = build_semantic_render_objects_resolved_limited(
        &geometry,
        &options,
        None,
        Some(&structure),
        |_| {},
    );
    let (_, _, object_stats) =
        flatten_semantic_render_objects_with_visible_bounding_sphere_and_stats(
            &objects, &geometry, &options, false,
        );
    render_object_summary_json_from_resolved(&options, &objects, &object_stats)
}

fn render_summaries_json_from_resolved(
    options: &MeshOptions,
    structure: AtomicStructure,
    geometry: GeometryInfoSnapshot,
    objects: &[SemanticRenderObject],
    object_stats: &[RenderObjectMeshStats],
) -> RenderSummaries {
    let representation_json =
        representation_summary_json_from_resolved(options, &structure, objects);
    RenderSummaries {
        render_objects_json: render_object_summary_json_from_resolved(
            options,
            objects,
            object_stats,
        ),
        representation_json,
        structure,
        geometry,
    }
}

impl GeometryInfoSnapshot {
    fn from_molecule(molecule: &Molecule) -> Self {
        Self {
            atom_count: molecule.atoms.len(),
            coarse_sphere_count: molecule.coarse_spheres.len(),
            coarse_gaussian_count: molecule.coarse_gaussians.len(),
            bond_count: molecule.bonds.len(),
            bond_metadata: BondMetadataSnapshot::from_molecule(molecule),
            bounds: info_bounds_molecule(molecule),
        }
    }
}

impl BondMetadataSnapshot {
    fn from_molecule(molecule: &Molecule) -> Self {
        let count_source = |source: BondSource| {
            molecule
                .bond_metadata
                .iter()
                .filter(|metadata| metadata.source == source)
                .count()
        };
        let count_flag = |flag: BondFlags| {
            molecule
                .bond_metadata
                .iter()
                .filter(|metadata| metadata.flags.contains(flag))
                .count()
        };
        Self {
            count: molecule.bond_metadata.len(),
            computed: count_source(BondSource::Computed),
            pdb_conect: count_source(BondSource::PdbConect),
            struct_conn: count_source(BondSource::StructConn),
            index_pair: count_source(BondSource::IndexPair),
            chem_comp: count_source(BondSource::ChemComp),
            covalent: count_flag(BondFlags::COVALENT),
            metallic_coordination: count_flag(BondFlags::METALLIC_COORDINATION),
            hydrogen_bond: count_flag(BondFlags::HYDROGEN_BOND),
            disulfide: count_flag(BondFlags::DISULFIDE),
            aromatic: count_flag(BondFlags::AROMATIC),
            computed_flag: count_flag(BondFlags::COMPUTED),
            resonance: count_flag(BondFlags::RESONANCE),
            rings: molecule.resonance.ring_count,
            aromatic_rings: molecule.resonance.aromatic_ring_count,
            delocalized_bonds: molecule.resonance.delocalized_bond_count,
        }
    }
}

fn info_bounds_molecule(molecule: &Molecule) -> Option<(Vec3, Vec3)> {
    let mut points = molecule
        .atoms
        .iter()
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    for sphere in &molecule.coarse_spheres {
        let radius = Vec3::new(sphere.radius, sphere.radius, sphere.radius);
        points.push(sphere.position - radius);
        points.push(sphere.position + radius);
    }
    for gaussian in &molecule.coarse_gaussians {
        let extent = Vec3::new(
            gaussian.covariance[0][0].abs().sqrt().max(0.1),
            gaussian.covariance[1][1].abs().sqrt().max(0.1),
            gaussian.covariance[2][2].abs().sqrt().max(0.1),
        ) * gaussian.weight.abs().sqrt().max(0.1);
        points.push(gaussian.position - extent);
        points.push(gaussian.position + extent);
    }
    let first = points.first().copied()?;
    let mut min = first;
    let mut max = first;
    for point in &points[1..] {
        min = min.min(*point);
        max = max.max(*point);
    }
    Some((min, max))
}

fn render_object_summary_json_from_resolved(
    options: &MeshOptions,
    objects: &[SemanticRenderObject],
    object_stats: &[RenderObjectMeshStats],
) -> String {
    let semantic_group_count = objects
        .iter()
        .map(|object| object.group_id)
        .max()
        .map_or(0, |group_id| group_id + 1);
    let mut out = String::with_capacity(objects.len().saturating_mul(480).saturating_add(2));
    out.push('[');
    for (index, (object, stats)) in objects.iter().zip(object_stats).enumerate() {
        if index > 0 {
            out.push(',');
        }
        write!(
            out,
            "{{\"geometry_type\":\"{}\",\"visual\":\"{}\",\"representation\":\"{}\",\"secondary_type\":\"{}\",\"chain\":\"{}\",\"residue_start\":",
            json_escape(object.geometry_type),
            json_escape(object.visual),
            json_escape(object.representation),
            json_escape(object.secondary_type),
            json_escape(object.chain.as_deref().unwrap_or("")),
        )
        .expect("writing to String cannot fail");
        write_optional_i32(&mut out, object.residue_start);
        out.push_str(",\"residue_end\":");
        write_optional_i32(&mut out, object.residue_end);
        write!(
            out,
            ",\"group_id\":{},\"polymer_trace\":{{\"initial\":{},\"final\":{},\"sec_struc_first\":{},\"sec_struc_last\":{}}},\"value_cell\":{{\"group_id\":{},\"draw_count\":{},\"u_group_count\":{}}},\"valueCell\":{{\"drawCount\":{},\"uVertexCount\":{},\"uGroupCount\":{},\"instanceCount\":1,\"uInstanceCount\":1}}}}",
            object.group_id,
            bool_json(object.initial),
            bool_json(object.final_residue),
            bool_json(object.sec_struc_first),
            bool_json(object.sec_struc_last),
            object.group_id,
            object.object.face_estimate(options),
            semantic_group_count,
            stats.draw_count,
            stats.vertex_count,
            stats.group_count,
        )
        .expect("writing to String cannot fail");
    }
    out.push(']');
    out
}

fn write_optional_i32(out: &mut String, value: Option<i32>) {
    if let Some(value) = value {
        write!(out, "{value}").expect("writing to String cannot fail");
    } else {
        out.push_str("null");
    }
}

fn bool_json(value: bool) -> &'static str {
    if value {
        "true"
    } else {
        "false"
    }
}

#[cfg(test)]
pub(crate) fn representation_summary_json(molecule: &Molecule, options: &MeshOptions) -> String {
    let geometry = molecule.expanded_for_geometry();
    let options = resolved_mesh_options(&geometry, options);
    let structure = geometry.atomic_structure();
    let objects = build_semantic_render_objects_resolved_limited(
        &geometry,
        &options,
        None,
        Some(&structure),
        |_| {},
    );
    representation_summary_json_from_resolved(&options, &structure, &objects)
}

fn representation_summary_json_from_resolved(
    options: &MeshOptions,
    structure: &AtomicStructure,
    objects: &[SemanticRenderObject],
) -> String {
    let selected_visuals = selected_visuals(&structure, &options);
    let realized_visuals = realized_visuals(&structure, &options, &objects);
    format!(
        "{{\"name\":\"{}\",\"selected_visuals\":{},\"realized_visuals\":{}}}",
        json_escape(representation_name(options.representation)),
        json_string_array(&selected_visuals),
        json_string_array(&realized_visuals)
    )
}

#[cfg(test)]
pub(crate) fn render_object_span_summary_json(
    molecule: &Molecule,
    options: &MeshOptions,
) -> String {
    let geometry = molecule.expanded_for_geometry();
    let options = resolved_mesh_options(&geometry, options);
    let objects = build_semantic_render_objects_resolved(&geometry, &options);
    let cylinder_radial_segments = molstar_export_cylinder_radial_segments(
        objects
            .iter()
            .map(|object| render_object_export_cylinder_count(&object.object))
            .sum::<usize>(),
    );
    let (estimate, plans) = render_objects_mesh_plan(
        objects.iter().map(|object| &object.object),
        &options,
        cylinder_radial_segments,
    );
    let mut state = MeshBuilderState::with_capacity(estimate, objects.len(), false);
    let mut cylinder_cache = CylinderPrimitiveCache::default();
    let mut curve_scratch = CurveSegmentScratch::default();
    let spans = objects
        .iter()
        .zip(&plans)
        .enumerate()
        .map(|(index, (object, plan))| {
            let vertex_start = state.mesh.vertices.len();
            let face_start = state.mesh.faces.len();
            state.set_current_group(object.group_id);
            append_render_object_to_mesh(
                &mut state.mesh,
                &object.object,
                &options,
                cylinder_radial_segments,
                &mut cylinder_cache,
                &mut curve_scratch,
                Some(plan),
            );
            state.mark_appended(vertex_start, face_start);
            let vertex_end = state.mesh.vertices.len();
            let face_end = state.mesh.faces.len();
            let last_vertex = vertex_end
                .checked_sub(1)
                .and_then(|i| state.mesh.vertices.get(i).copied())
                .unwrap_or_default();
            let last_normal = vertex_end
                .checked_sub(1)
                .and_then(|i| state.mesh.normals.get(i).copied())
                .unwrap_or_default();
            let stl_facet_start = face_start * 3;
            let stl_facet_end = face_end * 3;
            let first_face = state
                .mesh
                .faces
                .get(face_start)
                .map(|face| render_span_face_json(&state.mesh, face))
                .unwrap_or_else(|| "null".to_string());
            let bool_json = |value| if value { "true" } else { "false" };
            format!(
                "{{\"index\":{},\"geometry_type\":\"{}\",\"visual\":\"{}\",\"representation\":\"{}\",\"secondary_type\":\"{}\",\"chain\":\"{}\",\"residue_start\":{},\"residue_end\":{},\"group_id\":{},\"polymer_trace\":{{\"initial\":{},\"final\":{},\"sec_struc_first\":{},\"sec_struc_last\":{}}},\"vertex_start\":{},\"vertex_end\":{},\"face_start\":{},\"face_end\":{},\"stl_facet_start\":{},\"stl_facet_end\":{},\"first_face\":{},\"last_vertex\":{},\"last_normal\":{},\"last_normal_length\":{:.6}}}",
                index,
                json_escape(object.geometry_type),
                json_escape(object.visual),
                json_escape(object.representation),
                json_escape(object.secondary_type),
                json_escape(object.chain.as_deref().unwrap_or("")),
                object
                    .residue_start
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "null".to_string()),
                object
                    .residue_end
                    .map(|value| value.to_string())
                    .unwrap_or_else(|| "null".to_string()),
                object.group_id,
                bool_json(object.initial),
                bool_json(object.final_residue),
                bool_json(object.sec_struc_first),
                bool_json(object.sec_struc_last),
                vertex_start,
                vertex_end,
                face_start,
                face_end,
                stl_facet_start,
                stl_facet_end,
                first_face,
                render_span_vec3_json(last_vertex),
                render_span_vec3_json(last_normal),
                last_normal.length(),
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"object_count\":{},\"vertex_count\":{},\"normal_count\":{},\"face_count\":{},\"spans\":[{}]}}",
        objects.len(),
        state.mesh.vertices.len(),
        state.mesh.normals.len(),
        state.mesh.faces.len(),
        spans
    )
}

pub(crate) fn render_object_stl_facet_context_json(
    molecule: &Molecule,
    options: &MeshOptions,
    stl_facet: usize,
    vertex_offset: [f64; 3],
) -> String {
    render_object_stl_facet_context_json_timed(molecule, options, stl_facet, vertex_offset, |_| {})
}

pub(crate) fn render_object_stl_facet_context_json_timed(
    molecule: &Molecule,
    options: &MeshOptions,
    stl_facet: usize,
    vertex_offset: [f64; 3],
    mut checkpoint: impl FnMut(&str),
) -> String {
    checkpoint("begin-render-stl-facet-context");
    let geometry = molecule
        .identity_assembly_trace_subset_for_geometry()
        .map(std::borrow::Cow::Owned)
        .unwrap_or_else(|| molecule.expanded_for_geometry());
    checkpoint("expand-geometry");
    let options = resolved_mesh_options(&geometry, options);
    checkpoint("resolve-geometry-options");
    render_object_stl_facet_context_from_resolved_geometry_json_timed(
        &geometry,
        &options,
        stl_facet,
        vertex_offset,
        None,
        checkpoint,
    )
}

pub(crate) fn render_object_stl_facet_context_for_geometry_json_timed(
    molecule: &Molecule,
    options: &MeshOptions,
    stl_facet: usize,
    vertex_offset: [f64; 3],
    structure: Option<&AtomicStructure>,
    mut checkpoint: impl FnMut(&str),
) -> String {
    checkpoint("begin-render-stl-facet-context");
    let options = resolved_mesh_options(molecule, options);
    checkpoint("resolve-geometry-options");
    render_object_stl_facet_context_from_resolved_geometry_json_timed(
        molecule,
        &options,
        stl_facet,
        vertex_offset,
        structure,
        checkpoint,
    )
}

fn render_object_stl_facet_context_from_resolved_geometry_json_timed(
    geometry: &Molecule,
    options: &MeshOptions,
    stl_facet: usize,
    vertex_offset: [f64; 3],
    structure: Option<&AtomicStructure>,
    mut checkpoint: impl FnMut(&str),
) -> String {
    let face_index = stl_facet / 3;
    let sparse_slot = stl_facet % 3;
    let objects = build_semantic_render_objects_resolved_until_face_timed(
        geometry,
        options,
        face_index,
        structure,
        |label| checkpoint(label),
    );
    checkpoint("build-semantic-render-objects");
    let cylinder_radial_segments = molstar_export_cylinder_radial_segments(
        objects
            .iter()
            .map(|object| render_object_export_cylinder_count(&object.object))
            .sum::<usize>(),
    );
    checkpoint("resolve-cylinder-segments");
    let mut face_start = 0usize;
    let mut cylinder_cache = CylinderPrimitiveCache::default();
    let mut curve_scratch = CurveSegmentScratch::default();

    for (index, object) in objects.iter().enumerate() {
        let estimated_face_count = object.object.face_estimate(&options);
        let estimated_face_end = face_start.saturating_add(estimated_face_count);
        if face_index < face_start || face_index >= estimated_face_end {
            face_start = estimated_face_end;
            continue;
        }

        let plan = render_objects_mesh_plan(
            std::iter::once(&object.object),
            &options,
            cylinder_radial_segments,
        )
        .1
        .pop()
        .expect("single render object plan");
        let mut mesh = mesh_with_capacity(plan.estimate);
        append_render_object_to_mesh(
            &mut mesh,
            &object.object,
            &options,
            cylinder_radial_segments,
            &mut cylinder_cache,
            &mut curve_scratch,
            Some(&plan),
        );
        checkpoint("append-target-render-object");
        let actual_face_count = mesh.faces.len();
        let face_end = face_start.saturating_add(actual_face_count);
        if face_index >= face_end {
            face_start = face_end;
            continue;
        }

        let local_face_index = face_index - face_start;
        let target_face = mesh
            .faces
            .get(local_face_index)
            .map(|face| render_stl_target_face_json(&mesh, face, vertex_offset))
            .unwrap_or_else(|| "null".to_string());
        let bool_json = |value| if value { "true" } else { "false" };
        return format!(
            "{{\"found\":true,\"stl_facet\":{},\"stl_sparse_slot\":{},\"face_index\":{},\"face_offset_in_span\":{},\"vertex_offset\":{},\"span\":{{\"index\":{},\"geometry_type\":\"{}\",\"visual\":\"{}\",\"representation\":\"{}\",\"secondary_type\":\"{}\",\"chain\":\"{}\",\"residue_start\":{},\"residue_end\":{},\"group_id\":{},\"polymer_trace\":{{\"initial\":{},\"final\":{},\"sec_struc_first\":{},\"sec_struc_last\":{}}},\"vertex_start\":null,\"vertex_end\":null,\"local_vertex_start\":0,\"local_vertex_end\":{},\"face_start\":{},\"face_end\":{},\"estimated_face_end\":{},\"estimated_face_count\":{},\"stl_facet_start\":{},\"stl_facet_end\":{}}},\"target_face\":{}}}",
            stl_facet,
            sparse_slot,
            face_index,
            local_face_index,
            render_f64_triplet_json(vertex_offset),
            index,
            json_escape(object.geometry_type),
            json_escape(object.visual),
            json_escape(object.representation),
            json_escape(object.secondary_type),
            json_escape(object.chain.as_deref().unwrap_or("")),
            object
                .residue_start
                .map(|value| value.to_string())
                .unwrap_or_else(|| "null".to_string()),
            object
                .residue_end
                .map(|value| value.to_string())
                .unwrap_or_else(|| "null".to_string()),
            object.group_id,
            bool_json(object.initial),
            bool_json(object.final_residue),
            bool_json(object.sec_struc_first),
            bool_json(object.sec_struc_last),
            mesh.vertices.len(),
            face_start,
            face_end,
            estimated_face_end,
            estimated_face_count,
            face_start * 3,
            face_end * 3,
            target_face,
        );
    }

    format!(
        "{{\"found\":false,\"stl_facet\":{},\"stl_sparse_slot\":{},\"face_index\":{},\"object_count\":{},\"face_count\":{},\"stl_facet_count\":{}}}",
        stl_facet,
        sparse_slot,
        face_index,
        objects.len(),
        face_start,
        face_start * 3,
    )
}

fn render_span_vec3_json(value: Vec3) -> String {
    format!("[{:.6},{:.6},{:.6}]", value.x, value.y, value.z)
}

fn sphere_report_vec3_json(value: Vec3) -> String {
    format!(
        "[{:.9},{:.9},{:.9}]",
        value.x as f64, value.y as f64, value.z as f64
    )
}

fn sphere_report_vec3_64_json(value: [f64; 3]) -> String {
    format!("[{:.17},{:.17},{:.17}]", value[0], value[1], value[2])
}

#[cfg(test)]
fn render_span_face_json(mesh: &Mesh, face: &crate::model::Face) -> String {
    let indices = [face.a, face.b, face.c];
    let vertices = indices
        .iter()
        .map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(render_span_vec3_json)
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    let normals = indices
        .iter()
        .map(|&index| {
            mesh.normals
                .get(index)
                .copied()
                .map(render_span_vec3_json)
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"indices\":[{},{},{}],\"vertices\":[{}],\"normals\":[{}]}}",
        indices[0], indices[1], indices[2], vertices, normals
    )
}

fn render_stl_target_face_json(
    mesh: &Mesh,
    face: &crate::model::Face,
    vertex_offset: [f64; 3],
) -> String {
    let indices = [face.a, face.b, face.c];
    let vertices = indices
        .iter()
        .map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(render_span_vec3_json)
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    let raw_vertices = indices
        .iter()
        .map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(render_precise_vec3_json)
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    let stl_vertices = indices
        .iter()
        .map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(|vertex| {
                    render_precise_vec3_json(stl_transformed_vertex(vertex, vertex_offset))
                })
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    let stl_vertex_bits = indices
        .iter()
        .map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(|vertex| render_vec3_bits_json(stl_transformed_vertex(vertex, vertex_offset)))
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    let transformed_vertices = indices
        .iter()
        .filter_map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(|vertex| stl_transformed_vertex(vertex, vertex_offset))
        })
        .collect::<Vec<_>>();
    let (stl_normal, stl_normal_bits) = if transformed_vertices.len() == 3 {
        let normal = stl_triangle_normal(
            transformed_vertices[0],
            transformed_vertices[1],
            transformed_vertices[2],
        );
        (
            render_precise_vec3_json(normal),
            render_vec3_bits_json(normal),
        )
    } else {
        ("null".to_string(), "null".to_string())
    };
    let normals = indices
        .iter()
        .map(|&index| {
            mesh.normals
                .get(index)
                .copied()
                .map(render_span_vec3_json)
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"indices\":[{},{},{}],\"vertices\":[{}],\"raw_vertices\":[{}],\"stl_vertices\":[{}],\"stl_vertex_bits\":[{}],\"stl_normal\":{},\"stl_normal_bits\":{},\"normals\":[{}]}}",
        indices[0],
        indices[1],
        indices[2],
        vertices,
        raw_vertices,
        stl_vertices,
        stl_vertex_bits,
        stl_normal,
        stl_normal_bits,
        normals
    )
}

fn stl_transformed_vertex(vertex: Vec3, offset: [f64; 3]) -> Vec3 {
    Vec3::new(
        (vertex.x as f64 + offset[0]) as f32,
        (vertex.y as f64 + offset[1]) as f32,
        (vertex.z as f64 + offset[2]) as f32,
    )
}

fn render_precise_vec3_json(value: Vec3) -> String {
    format!(
        "[{:.9},{:.9},{:.9}]",
        value.x as f64, value.y as f64, value.z as f64
    )
}

fn render_vec3_bits_json(value: Vec3) -> String {
    format!(
        "[\"0x{:08x}\",\"0x{:08x}\",\"0x{:08x}\"]",
        value.x.to_bits(),
        value.y.to_bits(),
        value.z.to_bits()
    )
}

fn stl_triangle_normal(a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let ab = [
        b.x as f64 - a.x as f64,
        b.y as f64 - a.y as f64,
        b.z as f64 - a.z as f64,
    ];
    let ac = [
        c.x as f64 - a.x as f64,
        c.y as f64 - a.y as f64,
        c.z as f64 - a.z as f64,
    ];
    let n = [
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    ];
    let len_sq = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if len_sq > 0.0 {
        let scale = 1.0 / len_sq.sqrt();
        Vec3::new(
            (n[0] * scale) as f32,
            (n[1] * scale) as f32,
            (n[2] * scale) as f32,
        )
    } else {
        Vec3::default()
    }
}

fn render_f64_triplet_json(value: [f64; 3]) -> String {
    format!("[{:.17},{:.17},{:.17}]", value[0], value[1], value[2])
}

#[cfg(test)]
pub(crate) fn polymer_trace_iterator_reference_json(name: &str, molecule: &Molecule) -> String {
    polymer_trace_iterator_reference_json_with_options(name, molecule, false)
}

#[cfg(test)]
const POLYMER_TRACE_ITERATOR_SOURCE_FIELDS: [&str; 20] = [
    "center",
    "centerPrev",
    "centerNext",
    "first",
    "last",
    "initial",
    "final",
    "secStrucFirst",
    "secStrucLast",
    "secStrucType",
    "moleculeType",
    "coarseBackboneFirst",
    "coarseBackboneLast",
    "p0",
    "p1",
    "p2",
    "p3",
    "p4",
    "d12",
    "d23",
];

#[cfg(test)]
pub(crate) fn polymer_trace_iterator_reference_json_with_helix_orientation(
    name: &str,
    molecule: &Molecule,
) -> String {
    polymer_trace_iterator_reference_json_with_options(name, molecule, true)
}

#[cfg(test)]
fn polymer_trace_iterator_reference_json_with_options(
    name: &str,
    molecule: &Molecule,
    use_helix_orientation: bool,
) -> String {
    let structure = molecule.atomic_structure();
    let mut trace = backbone_residues(molecule, &structure);
    apply_polymer_trace_terminal_flags(&structure, &mut trace);
    apply_cyclic_polymer_trace_flags(&structure, &mut trace);
    apply_polymer_trace_secondary_flags(&structure, &mut trace);
    let hierarchy = &structure.model.hierarchy;

    let polymer_ranges = paired_usize_json(&structure.ranges.polymer_ranges);
    let cyclic_polymer_map = structure
        .ranges
        .cyclic_polymer_map
        .iter()
        .map(|(from, to)| format!("[{},{}]", from, to))
        .collect::<Vec<_>>()
        .join(",");

    let mut records = Vec::new();
    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        let Some(start_residue) = hierarchy.residue_atom_segments.index.get(pair[0]).copied()
        else {
            continue;
        };
        let Some(end_residue) = hierarchy.residue_atom_segments.index.get(pair[1]).copied() else {
            continue;
        };
        for residue_index in start_residue..=end_residue {
            let Some(trace_index) =
                trace_residue_index_for_model_residue(hierarchy, &trace, residue_index)
            else {
                continue;
            };
            let Some(record) = polymer_trace_iterator_record_json(
                &structure,
                &trace,
                start_residue,
                end_residue,
                residue_index,
                trace_index,
                use_helix_orientation,
            ) else {
                continue;
            };
            records.push(record);
        }
    }

    format!(
        "{{\"molstar_reference_commit\":\"1b8117d3f10f7c978aabb5a0d3d47370635aefe4\",\"molstar_source\":\"artifacts/molstar/src/mol-repr/structure/visual/util/polymer/trace-iterator.ts\",\"case\":\"{}\"{},\"source_fields\":{},\"polymer_ranges\":{},\"cyclic_polymer_map\":[{}],\"records\":[{}]}}",
        json_escape(name),
        if use_helix_orientation {
            ",\"use_helix_orientation\":true"
        } else {
            ""
        },
        json_string_array(
            &POLYMER_TRACE_ITERATOR_SOURCE_FIELDS
                .iter()
                .map(|field| field.to_string())
                .collect::<Vec<_>>(),
        ),
        polymer_ranges,
        cyclic_polymer_map,
        records.join(",")
    )
}

#[cfg(test)]
pub(crate) fn coarse_polymer_trace_iterator_reference_json(
    name: &str,
    molecule: &Molecule,
    kind: CoarseElementKind,
) -> String {
    let structure = molecule.atomic_structure();
    let (unit_kind, unit_kind_name, elements, ranges) = match kind {
        CoarseElementKind::Spheres => (
            crate::model::AtomicUnitKind::Spheres,
            "spheres",
            &structure.coarse.hierarchy.spheres,
            &structure.coarse.hierarchy.spheres.polymer_ranges,
        ),
        CoarseElementKind::Gaussians => (
            crate::model::AtomicUnitKind::Gaussians,
            "gaussians",
            &structure.coarse.hierarchy.gaussians,
            &structure.coarse.hierarchy.gaussians.polymer_ranges,
        ),
    };
    let polymer_ranges = coarse_ranges_json(ranges);
    let records = structure
        .units
        .iter()
        .filter(|unit| unit.kind == unit_kind)
        .flat_map(|unit| {
            coarse_polymer_trace_iterator_unit_records(
                &structure,
                unit.id,
                &unit.elements,
                elements,
                ranges,
                kind,
            )
        })
        .collect::<Vec<_>>();

    format!(
        "{{\"molstar_reference_commit\":\"1b8117d3f10f7c978aabb5a0d3d47370635aefe4\",\"molstar_source\":\"artifacts/molstar/src/mol-repr/structure/visual/util/polymer/trace-iterator.ts\",\"case\":\"{}\",\"unit_kind\":\"{}\",\"source_fields\":{},\"polymer_ranges\":{},\"records\":[{}]}}",
        json_escape(name),
        unit_kind_name,
        json_string_array(
            &POLYMER_TRACE_ITERATOR_SOURCE_FIELDS
                .iter()
                .map(|field| field.to_string())
                .collect::<Vec<_>>(),
        ),
        polymer_ranges,
        records.join(",")
    )
}

#[cfg(test)]
fn coarse_ranges_json(ranges: &[CoarseRange]) -> String {
    format!(
        "[{}]",
        ranges
            .iter()
            .map(|range| format!("[{},{}]", range.start_element, range.end_element))
            .collect::<Vec<_>>()
            .join(",")
    )
}

#[cfg(test)]
fn coarse_polymer_trace_iterator_unit_records(
    structure: &AtomicStructure,
    unit_id: usize,
    unit_elements: &[usize],
    coarse_elements: &CoarseElements,
    ranges: &[CoarseRange],
    kind: CoarseElementKind,
) -> Vec<String> {
    let mut records = Vec::new();
    for range in ranges {
        let segment_start = lower_bound_usize(unit_elements, range.start_element);
        let segment_end = lower_bound_usize(unit_elements, range.end_element.saturating_add(1));
        if segment_start >= segment_end {
            continue;
        }
        for element_index in segment_start..segment_end {
            if let Some(record) = coarse_polymer_trace_iterator_record_json(
                structure,
                unit_id,
                unit_elements,
                coarse_elements,
                segment_start,
                segment_end,
                element_index,
                kind,
            ) {
                records.push(record);
            }
        }
    }
    records
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn coarse_polymer_trace_iterator_record_json(
    structure: &AtomicStructure,
    unit_id: usize,
    unit_elements: &[usize],
    coarse_elements: &CoarseElements,
    segment_start: usize,
    segment_end: usize,
    element_index: usize,
    kind: CoarseElementKind,
) -> Option<String> {
    let element_index_prev1 =
        coarse_polymer_trace_element_index(segment_start, segment_end, element_index as isize - 1);
    let element_index_next1 =
        coarse_polymer_trace_element_index(segment_start, segment_end, element_index as isize + 1);
    let source_prev = *unit_elements.get(element_index_prev1)?;
    let source = *unit_elements.get(element_index)?;
    let source_next = *unit_elements.get(element_index_next1)?;
    let element = coarse_elements.elements.get(source)?;
    let state = coarse_polymer_trace_iterator_state(
        structure,
        unit_elements,
        segment_start,
        segment_end,
        element_index,
        kind,
    )?;
    let bool_json = |value| if value { "true" } else { "false" };

    Some(format!(
        "{{\"unit_id\":{},\"element_index\":{},\"source_element\":{},\"chain\":\"{}\",\"seq_begin\":{},\"seq_end\":{},\"center_prev\":{},\"center\":{},\"center_next\":{},\"first\":{},\"last\":{},\"initial\":false,\"final\":false,\"sec_struc_first\":false,\"sec_struc_last\":false,\"sec_struc_type\":\"na\",\"molecule_type\":\"unknown\",\"coarse_backbone_first\":false,\"coarse_backbone_last\":false,\"p0\":{},\"p1\":{},\"p2\":{},\"p3\":{},\"p4\":{},\"d12\":{},\"d23\":{}}}",
        unit_id,
        element_index,
        source,
        json_escape(&element.asym_id),
        element.seq_id_begin,
        element.seq_id_end,
        source_prev,
        source,
        source_next,
        bool_json(element_index == segment_start),
        bool_json(element_index + 1 == segment_end),
        vec3_json(state.p0.to_vec3()),
        vec3_json(state.p1.to_vec3()),
        vec3_json(state.p2.to_vec3()),
        vec3_json(state.p3.to_vec3()),
        vec3_json(state.p4.to_vec3()),
        vec3_json(state.d12.to_vec3()),
        vec3_json(state.d23.to_vec3())
    ))
}

#[cfg(test)]
fn coarse_polymer_trace_iterator_state(
    structure: &AtomicStructure,
    unit_elements: &[usize],
    segment_start: usize,
    segment_end: usize,
    element_index: usize,
    kind: CoarseElementKind,
) -> Option<PolymerTraceIteratorStateSnapshot> {
    let element_index_prev2 =
        coarse_polymer_trace_element_index(segment_start, segment_end, element_index as isize - 2);
    let element_index_prev1 =
        coarse_polymer_trace_element_index(segment_start, segment_end, element_index as isize - 1);
    let element_index_next1 =
        coarse_polymer_trace_element_index(segment_start, segment_end, element_index as isize + 1);
    let element_index_next2 =
        coarse_polymer_trace_element_index(segment_start, segment_end, element_index as isize + 2);

    let mut p0 = DVec3::from_vec3(coarse_polymer_trace_position(
        structure,
        unit_elements,
        element_index_prev2,
        kind,
    )?);
    let mut p1 = DVec3::from_vec3(coarse_polymer_trace_position(
        structure,
        unit_elements,
        element_index_prev1,
        kind,
    )?);
    let p2 = DVec3::from_vec3(coarse_polymer_trace_position(
        structure,
        unit_elements,
        element_index,
        kind,
    )?);
    let mut p3 = DVec3::from_vec3(coarse_polymer_trace_position(
        structure,
        unit_elements,
        element_index_next1,
        kind,
    )?);
    let mut p4 = DVec3::from_vec3(coarse_polymer_trace_position(
        structure,
        unit_elements,
        element_index_next2,
        kind,
    )?);

    let f = 0.5;
    if element_index == element_index_prev1 {
        let dir = (p2 - p3) * f;
        p1 = p2 + dir;
        p0 = p1 + dir;
    } else if element_index_prev1 == element_index_prev2 {
        let dir = (p1 - p2) * f;
        p0 = p1 + dir;
    }
    if element_index == element_index_next1 {
        let dir = (p2 - p1) * f;
        p3 = p2 + dir;
        p4 = p3 + dir;
    } else if element_index_next1 == element_index_next2 {
        let dir = (p3 - p2) * f;
        p4 = p3 + dir;
    }

    Some(PolymerTraceIteratorStateSnapshot {
        p0,
        p1,
        p2,
        p3,
        p4,
        d12: DVec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        },
        d23: DVec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        },
    })
}

#[cfg(test)]
fn coarse_polymer_trace_position(
    structure: &AtomicStructure,
    unit_elements: &[usize],
    element_index: usize,
    kind: CoarseElementKind,
) -> Option<Vec3> {
    let source_index = *unit_elements.get(element_index)?;
    match kind {
        CoarseElementKind::Spheres => structure.coarse.conformation.spheres.position(source_index),
        CoarseElementKind::Gaussians => structure
            .coarse
            .conformation
            .gaussians
            .position(source_index),
    }
}

#[cfg(test)]
fn coarse_polymer_trace_element_index(
    segment_start: usize,
    segment_end: usize,
    element_index: isize,
) -> usize {
    (element_index.max(segment_start as isize) as usize).min(segment_end.saturating_sub(1))
}

#[cfg(test)]
fn lower_bound_usize(values: &[usize], target: usize) -> usize {
    values.partition_point(|value| *value < target)
}

#[cfg(test)]
fn paired_usize_json(values: &[usize]) -> String {
    format!(
        "[{}]",
        values
            .chunks_exact(2)
            .map(|pair| format!("[{},{}]", pair[0], pair[1]))
            .collect::<Vec<_>>()
            .join(",")
    )
}

#[cfg(test)]
fn polymer_trace_iterator_record_json(
    structure: &AtomicStructure,
    trace: &[TraceResidue],
    segment_min: usize,
    segment_max: usize,
    residue_index: usize,
    trace_index: usize,
    use_helix_orientation: bool,
) -> Option<String> {
    let hierarchy = &structure.model.hierarchy;
    let residue = hierarchy.residues.get(residue_index)?;
    let chain = hierarchy.chains.get(residue.chain_index)?;
    let _trace_residue = trace.get(trace_index)?;
    let prev_residue = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize - 1,
    );
    let next_residue = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize + 1,
    );
    let current_type = molstar_secondary_trace_type(structure, residue_index);
    let previous_type = molstar_secondary_trace_type(structure, prev_residue);
    let next_type = molstar_secondary_trace_type(structure, next_residue);
    let center_prev = hierarchy
        .derived
        .residue
        .trace_element_index
        .get(prev_residue)
        .and_then(|index| *index);
    let center = hierarchy
        .derived
        .residue
        .trace_element_index
        .get(residue_index)
        .and_then(|index| *index);
    let center_next = hierarchy
        .derived
        .residue
        .trace_element_index
        .get(next_residue)
        .and_then(|index| *index);
    let previous_coarse_backbone = polymer_trace_coarse_backbone(structure, prev_residue);
    let current_coarse_backbone = polymer_trace_coarse_backbone(structure, residue_index);
    let next_coarse_backbone = polymer_trace_coarse_backbone(structure, next_residue);
    let state = polymer_trace_iterator_state(
        structure,
        segment_min,
        segment_max,
        residue_index,
        current_type,
        use_helix_orientation,
    );
    let molecule_type = hierarchy
        .derived
        .residue
        .molecule_type
        .get(residue_index)
        .copied()
        .unwrap_or_default();
    let bool_json = |value| if value { "true" } else { "false" };
    Some(format!(
        "{{\"residue_index\":{},\"trace_index\":{},\"chain\":\"{}\",\"seq\":\"{}\",\"ins\":\"{}\",\"center_prev\":{},\"center\":{},\"center_next\":{},\"first\":{},\"last\":{},\"initial\":{},\"final\":{},\"sec_struc_first\":{},\"sec_struc_last\":{},\"sec_struc_type\":\"{}\",\"molecule_type\":\"{}\",\"coarse_backbone_first\":{},\"coarse_backbone_last\":{},\"p0\":{},\"p1\":{},\"p2\":{},\"p3\":{},\"p4\":{},\"d12\":{},\"d23\":{}}}",
        residue_index,
        trace_index,
        json_escape(&chain.id),
        json_escape(&residue.label_seq_id),
        json_escape(&residue.insertion_code),
        opt_usize_json(center_prev),
        opt_usize_json(center),
        opt_usize_json(center_next),
        bool_json(residue_index == segment_min),
        bool_json(residue_index == segment_max),
        bool_json(residue_index == prev_residue),
        bool_json(residue_index == next_residue),
        bool_json(previous_type != current_type),
        bool_json(current_type != next_type),
        secondary_trace_type_name(current_type),
        molecule_type_name(molecule_type),
        bool_json(previous_coarse_backbone != current_coarse_backbone),
        bool_json(current_coarse_backbone != next_coarse_backbone),
        vec3_json(state.p0.to_vec3()),
        vec3_json(state.p1.to_vec3()),
        vec3_json(state.p2.to_vec3()),
        vec3_json(state.p3.to_vec3()),
        vec3_json(state.p4.to_vec3()),
        vec3_json(state.d12.to_vec3()),
        vec3_json(state.d23.to_vec3())
    ))
}

#[derive(Clone, Copy, Debug, Default)]
struct PolymerTraceIteratorStateSnapshot {
    p0: DVec3,
    p1: DVec3,
    p2: DVec3,
    p3: DVec3,
    p4: DVec3,
    d12: DVec3,
    d23: DVec3,
}

fn polymer_trace_iterator_state(
    structure: &AtomicStructure,
    segment_min: usize,
    segment_max: usize,
    residue_index: usize,
    current_type: SecondaryStructureType,
    use_helix_orientation: bool,
) -> PolymerTraceIteratorStateSnapshot {
    let residue_prev3 = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize - 3,
    );
    let residue_prev2 = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize - 2,
    );
    let residue_prev1 = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize - 1,
    );
    let residue_next1 = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize + 1,
    );
    let residue_next2 = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize + 2,
    );
    let residue_next3 = polymer_trace_residue_index(
        structure,
        segment_min,
        segment_max,
        residue_index as isize + 3,
    );

    let ss_prev3 = molstar_secondary_trace_type(structure, residue_prev3);
    let ss_prev2 = molstar_secondary_trace_type(structure, residue_prev2);
    let ss_prev1 = molstar_secondary_trace_type(structure, residue_prev1);
    let ss = current_type;
    let ss_next1 = molstar_secondary_trace_type(structure, residue_next1);
    let ss_next2 = molstar_secondary_trace_type(structure, residue_next2);
    let ss_next3 = molstar_secondary_trace_type(structure, residue_next3);

    let helix_orientation_centers =
        use_helix_orientation.then(|| molstar_helix_orientation_centers(structure));
    let has_helix_orientation = helix_orientation_centers.is_some();
    let trace_position = |residue_index: usize, ss: SecondaryStructureType| {
        if let Some(centers) = &helix_orientation_centers {
            if is_helix_secondary(ss) {
                if let Some(center) = centers.get(residue_index).copied() {
                    if center.is_finite() {
                        return DVec3::from_vec3(center);
                    }
                }
            }
        }
        DVec3::from_vec3(polymer_trace_position(structure, residue_index))
    };

    let mut p0 = trace_position(residue_prev3, ss_prev3);
    let mut p1 = trace_position(residue_prev2, ss_prev2);
    let mut p2 = trace_position(residue_prev1, ss_prev1);
    let p3 = trace_position(residue_index, ss);
    let mut p4 = trace_position(residue_next1, ss_next1);
    let mut p5 = trace_position(residue_next2, ss_next2);
    let mut p6 = trace_position(residue_next3, ss_next3);

    let is_helix_prev3 = is_helix_secondary(ss_prev3);
    let is_helix_prev2 = is_helix_secondary(ss_prev2);
    let is_helix_prev1 = is_helix_secondary(ss_prev1);
    let is_helix = is_helix_secondary(ss);
    let is_helix_next1 = is_helix_secondary(ss_next1);
    let is_helix_next2 = is_helix_secondary(ss_next2);
    let is_helix_next3 = is_helix_secondary(ss_next3);

    let sec_struc_first = ss_prev1 != ss;
    let sec_struc_last = ss != ss_next1;
    if has_helix_orientation && !(is_helix && sec_struc_first && sec_struc_last) {
        if is_helix != is_helix_prev1 {
            if is_helix {
                p0 = p3;
                p1 = p3;
                p2 = p3;
            } else if is_helix_prev1 {
                let dir = (p2 - p3) * 2.0;
                p2 = p3 + dir;
                p1 = p2 + dir;
                p0 = p1 + dir;
            }
        } else if is_helix != is_helix_prev2 {
            if is_helix {
                p0 = p2;
                p1 = p2;
            } else if is_helix_prev2 {
                let dir = (p1 - p2) * 2.0;
                p1 = p2 + dir;
                p0 = p1 + dir;
            }
        } else if is_helix != is_helix_prev3 {
            if is_helix {
                p0 = p1;
            } else if is_helix_prev3 {
                let dir = (p0 - p1) * 2.0;
                p0 = p1 + dir;
            }
        }

        if is_helix != is_helix_next1 {
            if is_helix {
                p4 = p3;
                p5 = p3;
                p6 = p3;
            } else if is_helix_next1 {
                let dir = (p4 - p3) * 2.0;
                p4 = p3 + dir;
                p5 = p4 + dir;
                p6 = p5 + dir;
            }
        } else if is_helix != is_helix_next2 {
            if is_helix {
                p5 = p4;
                p6 = p4;
            } else if is_helix_next2 {
                let dir = (p5 - p4) * 2.0;
                p5 = p4 + dir;
                p6 = p5 + dir;
            }
        } else if is_helix != is_helix_next3 {
            if is_helix {
                p6 = p5;
            } else if is_helix_next3 {
                let dir = (p6 - p5) * 2.0;
                p6 = p5 + dir;
            }
        }
    }

    let (d01, d12, d23, d34) = if polymer_trace_coarse_backbone(structure, residue_index) {
        (
            triangle_normal(p1, p2, p3),
            triangle_normal(p2, p3, p4),
            triangle_normal(p3, p4, p5),
            triangle_normal(p4, p5, p6),
        )
    } else {
        (
            polymer_trace_from_to_vector(structure, residue_prev1, ss_prev1, has_helix_orientation),
            polymer_trace_from_to_vector(structure, residue_index, ss, has_helix_orientation),
            polymer_trace_from_to_vector(structure, residue_next1, ss_next1, has_helix_orientation),
            polymer_trace_from_to_vector(structure, residue_next2, ss_next2, has_helix_orientation),
        )
    };

    let helix_flag = is_helix && has_helix_orientation;
    let f = 1.5;
    if residue_index == residue_prev1 || (ss != ss_prev1 && helix_flag) {
        let dir = set_magnitude(p3 - p4, f);
        p2 = p3 + dir;
        p1 = p2 + dir;
        p0 = p1 + dir;
    } else if residue_prev1 == residue_prev2 || (ss != ss_prev2 && helix_flag) {
        let dir = set_magnitude(p2 - p3, f);
        p1 = p2 + dir;
        p0 = p1 + dir;
    } else if residue_prev2 == residue_prev3 || (ss != ss_prev3 && helix_flag) {
        let dir = set_magnitude(p1 - p2, f);
        p0 = p1 + dir;
    }
    if residue_index == residue_next1 || (ss != ss_next1 && helix_flag) {
        let dir = set_magnitude(p3 - p2, f);
        p4 = p3 + dir;
        p5 = p4 + dir;
        p6 = p5 + dir;
    } else if residue_next1 == residue_next2 || (ss != ss_next2 && helix_flag) {
        let dir = set_magnitude(p4 - p3, f);
        p5 = p4 + dir;
        p6 = p5 + dir;
    } else if residue_next2 == residue_next3 || (ss != ss_next3 && helix_flag) {
        let dir = set_magnitude(p5 - p4, f);
        p6 = p5 + dir;
    }

    PolymerTraceIteratorStateSnapshot {
        p0: polymer_trace_control_point(p0, p1, p2, ss_prev2, has_helix_orientation),
        p1: polymer_trace_control_point(p1, p2, p3, ss_prev1, has_helix_orientation),
        p2: polymer_trace_control_point(p2, p3, p4, ss, has_helix_orientation),
        p3: polymer_trace_control_point(p3, p4, p5, ss_next1, has_helix_orientation),
        p4: polymer_trace_control_point(p4, p5, p6, ss_next2, has_helix_orientation),
        d12: polymer_trace_direction(d01, d12, d23),
        d23: polymer_trace_direction(d12, d23, d34),
    }
}

fn polymer_trace_position(structure: &AtomicStructure, residue_index: usize) -> Vec3 {
    polymer_trace_atom_position(structure, residue_index).unwrap_or_default()
}

fn polymer_trace_radius(
    structure: &AtomicStructure,
    residue_index: usize,
    options: &MeshOptions,
) -> f32 {
    let _ = (structure, residue_index);
    molstar_cartoon_uniform_trace_radius(options)
}

fn polymer_trace_coarse_backbone(structure: &AtomicStructure, residue_index: usize) -> bool {
    structure
        .model
        .hierarchy
        .derived
        .residue
        .direction_from_element_index
        .get(residue_index)
        .and_then(|index| *index)
        .is_none()
        || structure
            .model
            .hierarchy
            .derived
            .residue
            .direction_to_element_index
            .get(residue_index)
            .and_then(|index| *index)
            .is_none()
}

fn polymer_trace_from_to_vector(
    structure: &AtomicStructure,
    residue_index: usize,
    ss: SecondaryStructureType,
    has_helix_orientation: bool,
) -> DVec3 {
    if has_helix_orientation && is_helix_secondary(ss) {
        return DVec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        };
    }
    let hierarchy = &structure.model.hierarchy;
    let position = |index: Option<usize>| {
        index
            .and_then(|atom_index| hierarchy.atoms.get(atom_index))
            .map(|atom| atom.position)
            .unwrap_or_default()
    };
    let from = hierarchy
        .derived
        .residue
        .direction_from_element_index
        .get(residue_index)
        .and_then(|index| *index);
    let to = hierarchy
        .derived
        .residue
        .direction_to_element_index
        .get(residue_index)
        .and_then(|index| *index);
    DVec3::from_vec3(position(to)) - DVec3::from_vec3(position(from))
}

fn polymer_trace_control_point(
    p1: DVec3,
    p2: DVec3,
    p3: DVec3,
    ss: SecondaryStructureType,
    has_helix_orientation: bool,
) -> DVec3 {
    if ss.contains(SecondaryStructureType::BETA)
        || (has_helix_orientation && is_helix_secondary(ss))
    {
        vec3_average4_f64(p1, p3, p2, p2)
    } else {
        p2
    }
}

fn polymer_trace_direction(v1: DVec3, v2: DVec3, v3: DVec3) -> DVec3 {
    vec3_average4_f64(match_direction(v1, v2), match_direction(v3, v2), v2, v2)
}

fn match_direction(a: DVec3, b: DVec3) -> DVec3 {
    if a.dot(b) > 0.0 {
        a
    } else {
        a * -1.0
    }
}

fn set_magnitude(v: DVec3, magnitude: f32) -> DVec3 {
    v.normalized() * magnitude as f64
}

fn triangle_normal(a: DVec3, b: DVec3, c: DVec3) -> DVec3 {
    (b - a).cross(c - a).normalized()
}

fn vec3_average4_f64(a: DVec3, b: DVec3, c: DVec3, d: DVec3) -> DVec3 {
    DVec3 {
        x: (a.x + (b.x + (c.x + d.x))) * 0.25,
        y: (a.y + (b.y + (c.y + d.y))) * 0.25,
        z: (a.z + (b.z + (c.z + d.z))) * 0.25,
    }
}

#[cfg(test)]
fn vec3_json(value: Vec3) -> String {
    format!(
        "[{},{},{}]",
        f32_json(value.x),
        f32_json(value.y),
        f32_json(value.z)
    )
}

#[cfg(test)]
fn f32_json(value: f32) -> String {
    let value = canonical_zero(value);
    if (value.round() - value).abs() < 0.000_05 {
        format!("{value:.0}")
    } else {
        let mut out = format!("{value:.4}");
        while out.contains('.') && out.ends_with('0') {
            out.pop();
        }
        if out.ends_with('.') {
            out.pop();
        }
        out
    }
}

#[cfg(test)]
fn canonical_zero(value: f32) -> f32 {
    if value.abs() < 0.000_05 {
        0.0
    } else {
        value
    }
}

#[cfg(test)]
fn opt_usize_json(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string())
}

#[cfg(test)]
fn secondary_trace_type_name(value: SecondaryStructureType) -> &'static str {
    if value.contains(SecondaryStructureType::HELIX) {
        "helix"
    } else if value.contains(SecondaryStructureType::BETA) {
        "beta"
    } else if value == SecondaryStructureType::NONE {
        "none"
    } else {
        "other"
    }
}

#[cfg(test)]
fn molecule_type_name(value: MoleculeType) -> &'static str {
    match value {
        MoleculeType::Unknown => "unknown",
        MoleculeType::Other => "other",
        MoleculeType::Water => "water",
        MoleculeType::Ion => "ion",
        MoleculeType::Lipid => "lipid",
        MoleculeType::Protein => "protein",
        MoleculeType::Rna => "rna",
        MoleculeType::Dna => "dna",
        MoleculeType::Pna => "pna",
        MoleculeType::Saccharide => "saccharide",
    }
}

impl RenderObject {
    fn mesh_estimate(
        &self,
        options: &MeshOptions,
        cylinder_radial_segments: usize,
    ) -> RenderObjectMeshEstimate {
        match self {
            RenderObject::Sphere { .. } | RenderObject::Ellipsoid { .. } => {
                RenderObjectMeshEstimate::from_counts(molstar_sphere_mesh_counts(
                    options.sphere_detail,
                ))
            }
            RenderObject::Cylinder { start, end, radius } => {
                let midpoint = molstar_link_midpoint_buffer(*start, *end);
                cylinder_mesh_estimate(
                    *start,
                    midpoint,
                    *radius as f64,
                    cylinder_radial_segments,
                    false,
                    false,
                )
                .add(cylinder_mesh_estimate(
                    midpoint,
                    *end,
                    *radius as f64,
                    cylinder_radial_segments,
                    false,
                    false,
                ))
            }
            RenderObject::LinkCylinder { start, end, radius } => cylinder_mesh_estimate(
                *start,
                molstar_link_midpoint_buffer(*start, *end),
                *radius as f64,
                options.radial_segments.max(3),
                false,
                false,
            ),
            RenderObject::LinkCylinderWithSegments {
                start,
                end,
                radius,
                radial_segments,
            } => cylinder_mesh_estimate(
                *start,
                molstar_link_midpoint_buffer(*start, *end),
                *radius,
                (*radial_segments).max(3),
                false,
                false,
            ),
            RenderObject::Tube { points, .. } => profile_tube_mesh_estimate(
                sample_path_point_count(points, 4),
                options.radial_segments.max(3),
                true,
                true,
            ),
            RenderObject::DashedTube { points, radius } => {
                dashed_tube_mesh_estimate(points, *radius, options.radial_segments.max(3))
            }
            RenderObject::FixedCountDashedCylinder {
                start,
                end,
                radius,
                length_scale,
                segment_count,
            } => fixed_count_dashed_cylinder_mesh_estimate(
                *start,
                *end,
                *radius,
                options.radial_segments.max(3),
                *length_scale,
                *segment_count,
            ),
            RenderObject::Ribbon { points, .. } => ribbon_mesh_estimate(sample_path_point_count(
                points,
                options.linear_segments.max(1),
            )),
            RenderObject::Sheet {
                points,
                arrow_height,
                start_cap,
                end_cap,
                ..
            } => sheet_mesh_estimate(
                sample_path_point_count(points, options.linear_segments.max(1)),
                *arrow_height,
                *start_cap,
                *end_cap,
            ),
            RenderObject::OrientedRibbon {
                centers,
                normals,
                profile,
                start_cap,
                end_cap,
                ..
            } => {
                let sample_count = if centers.len() < 2 || centers.len() != normals.len() {
                    0
                } else {
                    (centers.len() - 1)
                        .saturating_mul(options.linear_segments.max(1))
                        .saturating_add(1)
                };
                if options.radial_segments == 2 {
                    ribbon_mesh_estimate(sample_count)
                } else if options.radial_segments == 4 || *profile == PolymerProfile::Square {
                    sheet_mesh_estimate(sample_count, 0.0, *start_cap, *end_cap)
                } else {
                    profile_tube_mesh_estimate(
                        sample_count,
                        options.radial_segments.max(3),
                        *start_cap,
                        *end_cap,
                    )
                }
            }
            RenderObject::PolymerTraceSegment {
                shift,
                kind,
                start_cap,
                end_cap,
                initial,
                final_residue,
                ..
            } => {
                let segment_count = polymer_trace_segment_count(
                    options.linear_segments.max(1),
                    *shift,
                    *initial,
                    *final_residue,
                );
                let sample_count = segment_count + 1;
                match kind {
                    PolymerTraceSegmentKind::Ribbon { .. } => ribbon_mesh_estimate(sample_count),
                    PolymerTraceSegmentKind::Tube { .. } => profile_tube_mesh_estimate(
                        sample_count,
                        options.radial_segments.max(3),
                        *start_cap,
                        *end_cap,
                    ),
                    PolymerTraceSegmentKind::Sheet { arrow_height } => {
                        sheet_mesh_estimate(sample_count, *arrow_height, *start_cap, *end_cap)
                    }
                }
            }
            RenderObject::NucleotideRing {
                radius,
                base,
                detail,
                radial_segments,
                ..
            } => nucleotide_ring_mesh_estimate(*base, *radius, *detail, *radial_segments),
            RenderObject::NucleotideBlock {
                geometry,
                radius,
                radial_segments,
                ..
            } => {
                let mut estimate = cylinder_mesh_estimate(
                    geometry.anchor,
                    geometry.trace,
                    *radius as f64,
                    (*radial_segments).max(3),
                    false,
                    true,
                );
                if geometry.block.is_some() {
                    estimate = estimate.add(RenderObjectMeshEstimate {
                        vertices: 24,
                        faces: 12,
                    });
                }
                estimate
            }
            RenderObject::DirectionWedge { .. } => RenderObjectMeshEstimate {
                vertices: 24,
                faces: 8,
            },
            RenderObject::CarbohydrateSymbol { shape, part, .. } => {
                carbohydrate_symbol_mesh_estimate(*shape, *part)
            }
        }
    }

    fn face_estimate(&self, options: &MeshOptions) -> usize {
        match self {
            RenderObject::Sphere { .. } => molstar_sphere_triangle_count(options.sphere_detail),
            RenderObject::Cylinder { .. } => {
                let segments = molstar_export_cylinder_radial_segments(2);
                segments * 4
            }
            RenderObject::LinkCylinder { .. } => options.radial_segments.max(3) * 2,
            RenderObject::LinkCylinderWithSegments {
                radial_segments, ..
            } => (*radial_segments).max(3) * 2,
            RenderObject::Tube { points, .. } | RenderObject::DashedTube { points, .. } => points
                .len()
                .saturating_sub(1)
                .saturating_mul(options.radial_segments.max(3) * 2),
            RenderObject::FixedCountDashedCylinder { segment_count, .. } => {
                let dash_count = segment_count.div_ceil(2);
                dash_count * options.radial_segments.max(3) * 4
            }
            RenderObject::Ribbon { points, .. } => {
                sample_path_point_count(points, options.linear_segments.max(1)).saturating_sub(1)
                    * 4
            }
            RenderObject::Sheet {
                points,
                arrow_height,
                start_cap,
                end_cap,
                ..
            } => {
                let segments = sample_path_point_count(points, options.linear_segments.max(1))
                    .saturating_sub(1);
                let caps = usize::from(*start_cap || *arrow_height > 0.0)
                    + usize::from(*end_cap && *arrow_height == 0.0);
                segments * 8 + caps * 2
            }
            RenderObject::OrientedRibbon {
                centers,
                profile,
                start_cap,
                end_cap,
                ..
            } => {
                let sample_count = if centers.len() < 2 {
                    0
                } else {
                    (centers.len() - 1)
                        .saturating_mul(options.linear_segments.max(1))
                        .saturating_add(1)
                };
                if options.radial_segments == 2 {
                    sample_count.saturating_sub(1) * 4
                } else if options.radial_segments == 4 || *profile == PolymerProfile::Square {
                    let caps = usize::from(*start_cap) + usize::from(*end_cap);
                    sample_count.saturating_sub(1) * 8 + caps * 2
                } else {
                    let radial = options.radial_segments.max(3);
                    let caps = usize::from(*start_cap) + usize::from(*end_cap);
                    sample_count.saturating_sub(1) * radial * 2 + caps * radial
                }
            }
            RenderObject::PolymerTraceSegment {
                shift,
                kind,
                start_cap,
                end_cap,
                initial,
                final_residue,
                ..
            } => {
                let segment_count = polymer_trace_segment_count(
                    options.linear_segments.max(1),
                    *shift,
                    *initial,
                    *final_residue,
                );
                match kind {
                    PolymerTraceSegmentKind::Ribbon { .. } => segment_count * 4,
                    PolymerTraceSegmentKind::Tube { .. } => {
                        let radial = options.radial_segments.max(3);
                        let caps = usize::from(*start_cap) + usize::from(*end_cap);
                        segment_count * radial * 2 + caps * radial
                    }
                    PolymerTraceSegmentKind::Sheet { arrow_height } => {
                        let caps = if *start_cap {
                            2
                        } else if *arrow_height > 0.0 {
                            4
                        } else {
                            0
                        } + usize::from(*end_cap && *arrow_height == 0.0) * 2;
                        segment_count * 8 + caps
                    }
                }
            }
            RenderObject::NucleotideRing {
                base,
                detail,
                radial_segments,
                ..
            } => nucleotide_ring_face_count(*base, *detail, *radial_segments),
            RenderObject::NucleotideBlock {
                geometry,
                radial_segments,
                ..
            } => {
                let cylinder_faces = (*radial_segments).max(3) * 3;
                let box_faces = usize::from(geometry.block.is_some()) * 12;
                cylinder_faces + box_faces
            }
            RenderObject::DirectionWedge { .. } => 8,
            RenderObject::CarbohydrateSymbol { shape, part, .. } => {
                carbohydrate_symbol_face_count(*shape, *part)
            }
            RenderObject::Ellipsoid { .. } => molstar_sphere_triangle_count(options.sphere_detail),
        }
    }
}

fn cylinder_mesh_estimate(
    start: Vec3,
    end: Vec3,
    radius: f64,
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
) -> RenderObjectMeshEstimate {
    if DVec3::from_vec3(start).distance(DVec3::from_vec3(end)) <= 0.001 {
        return RenderObjectMeshEstimate::default();
    }
    RenderObjectMeshEstimate::from_counts(molstar_cylinder_mesh_counts(
        radial_segments,
        top_cap,
        bottom_cap,
        radius,
    ))
}

fn profile_tube_mesh_estimate(
    sample_count: usize,
    radial_segments: usize,
    start_cap: bool,
    end_cap: bool,
) -> RenderObjectMeshEstimate {
    if sample_count < 2 || radial_segments < 3 {
        return RenderObjectMeshEstimate::default();
    }
    let cap_count = usize::from(start_cap) + usize::from(end_cap);
    RenderObjectMeshEstimate {
        vertices: sample_count
            .saturating_mul(radial_segments)
            .saturating_add(cap_count.saturating_mul(radial_segments + 1)),
        faces: (sample_count - 1)
            .saturating_mul(radial_segments)
            .saturating_mul(2)
            .saturating_add(cap_count.saturating_mul(radial_segments)),
    }
}

fn ribbon_mesh_estimate(sample_count: usize) -> RenderObjectMeshEstimate {
    if sample_count < 2 {
        return RenderObjectMeshEstimate::default();
    }
    RenderObjectMeshEstimate {
        vertices: sample_count.saturating_mul(4),
        faces: (sample_count - 1).saturating_mul(4),
    }
}

fn sheet_mesh_estimate(
    sample_count: usize,
    arrow_height: f32,
    start_cap: bool,
    end_cap: bool,
) -> RenderObjectMeshEstimate {
    if sample_count < 2 {
        return RenderObjectMeshEstimate::default();
    }
    let arrow = arrow_height.max(0.0) > 0.0;
    let cap_count = if start_cap {
        1
    } else if arrow {
        2
    } else {
        0
    } + usize::from(end_cap && !arrow);
    RenderObjectMeshEstimate {
        vertices: sample_count
            .saturating_mul(8)
            .saturating_add(cap_count.saturating_mul(4)),
        faces: (sample_count - 1)
            .saturating_mul(8)
            .saturating_add(cap_count.saturating_mul(2)),
    }
}

fn dashed_tube_mesh_estimate(
    points: &[Vec3],
    radius: f32,
    radial_segments: usize,
) -> RenderObjectMeshEstimate {
    if points.len() < 2 {
        return RenderObjectMeshEstimate::default();
    }
    let samples = sample_path(points, 8);
    dashed_tube_mesh_estimate_from_samples(&samples, radius, radial_segments)
}

fn dashed_tube_mesh_estimate_from_samples(
    samples: &[Vec3],
    radius: f32,
    radial_segments: usize,
) -> RenderObjectMeshEstimate {
    if samples.len() < 2 {
        return RenderObjectMeshEstimate::default();
    }
    let dash_len = (radius * 3.8).max(0.55);
    let gap_len = (radius * 2.2).max(0.32);
    let period = dash_len + gap_len;
    let mut distance = 0.0;
    let mut dash_count = 0usize;

    for pair in samples.windows(2) {
        let length = pair[0].distance(pair[1]);
        if length <= 0.000_001 {
            continue;
        }
        let mut local = 0.0;
        while local < length {
            let phase = (distance + local) % period;
            let in_dash = phase < dash_len;
            let remaining_phase = if in_dash {
                dash_len - phase
            } else {
                period - phase
            };
            let step = remaining_phase.min(length - local);
            if in_dash && step > 0.02 {
                dash_count += 1;
            }
            local += step.max(0.000_001);
        }
        distance += length;
    }

    RenderObjectMeshEstimate::from_counts(molstar_cylinder_mesh_counts(
        radial_segments,
        true,
        true,
        radius as f64,
    ))
    .scale(dash_count)
}

fn fixed_count_dashed_cylinder_mesh_estimate(
    start: Vec3,
    end: Vec3,
    radius: f32,
    radial_segments: usize,
    length_scale: f32,
    segment_count: usize,
) -> RenderObjectMeshEstimate {
    let distance = DVec3::from_vec3(start).distance(DVec3::from_vec3(end)) * length_scale as f64;
    if distance <= 0.000_001 || segment_count == 0 {
        return RenderObjectMeshEstimate::default();
    }

    let dash_count = segment_count.div_ceil(2);
    let is_odd = !segment_count.is_multiple_of(2);
    let step = distance / (segment_count as f64 + 0.5);
    let full = RenderObjectMeshEstimate::from_counts(molstar_cylinder_mesh_counts(
        radial_segments,
        true,
        true,
        radius as f64,
    ));
    let half = RenderObjectMeshEstimate::from_counts(molstar_cylinder_mesh_counts(
        radial_segments,
        false,
        true,
        radius as f64,
    ));
    let mut estimate = RenderObjectMeshEstimate::default();
    for dash_index in 0..dash_count {
        let last_odd = is_odd && dash_index + 1 == dash_count;
        let dash_length = if last_odd { step * 0.5 } else { step };
        if dash_length > 0.001 {
            estimate = estimate.add(if last_odd { half } else { full });
        }
    }
    estimate
}

fn nucleotide_ring_mesh_estimate(
    base: Option<NucleotideRingBase>,
    radius: f32,
    detail: usize,
    radial_segments: usize,
) -> RenderObjectMeshEstimate {
    let Some(base) = base else {
        return RenderObjectMeshEstimate::default();
    };
    let (trace, anchor, ring_faces) = match base {
        NucleotideRingBase::PurineConnector { trace, n9 } => (trace, n9, 0),
        NucleotideRingBase::Purine { trace, n9, .. } => {
            (trace, n9, molstar_nucleotide_ring_5_6_face_count())
        }
        NucleotideRingBase::PyrimidineConnector { trace, n1 } => (trace, n1, 0),
        NucleotideRingBase::Pyrimidine { trace, n1, .. } => {
            (trace, n1, molstar_nucleotide_ring_6_face_count())
        }
    };
    cylinder_mesh_estimate(
        anchor,
        trace,
        radius as f64,
        radial_segments.max(3),
        false,
        false,
    )
    .add(RenderObjectMeshEstimate::from_counts(
        molstar_sphere_mesh_counts(detail),
    ))
    .add(RenderObjectMeshEstimate {
        vertices: ring_faces.saturating_mul(3),
        faces: ring_faces,
    })
}

fn carbohydrate_symbol_mesh_estimate(
    shape: SaccharideShape,
    part: CarbohydrateSymbolPart,
) -> RenderObjectMeshEstimate {
    if part == CarbohydrateSymbolPart::Secondary && !carbohydrate_symbol_has_secondary_part(shape) {
        return RenderObjectMeshEstimate::default();
    }
    match shape {
        SaccharideShape::FilledSphere => RenderObjectMeshEstimate::from_counts(
            molstar_sphere_mesh_counts(MOLSTAR_CARBOHYDRATE_SYMBOL_DETAIL),
        ),
        SaccharideShape::FilledCube | SaccharideShape::FlatBox => RenderObjectMeshEstimate {
            vertices: 24,
            faces: 12,
        },
        SaccharideShape::CrossedCube => RenderObjectMeshEstimate {
            vertices: 18,
            faces: 6,
        },
        SaccharideShape::FilledCone => RenderObjectMeshEstimate {
            vertices: 48,
            faces: 16,
        },
        SaccharideShape::DevidedCone => RenderObjectMeshEstimate {
            vertices: 24,
            faces: 8,
        },
        SaccharideShape::FilledStar => RenderObjectMeshEstimate {
            vertices: 60,
            faces: 20,
        },
        SaccharideShape::FilledDiamond => RenderObjectMeshEstimate {
            vertices: 24,
            faces: 8,
        },
        SaccharideShape::DividedDiamond => RenderObjectMeshEstimate {
            vertices: 12,
            faces: 4,
        },
        SaccharideShape::FlatDiamond | SaccharideShape::DiamondPrism => RenderObjectMeshEstimate {
            vertices: 24,
            faces: 12,
        },
        SaccharideShape::PentagonalPrism | SaccharideShape::Pentagon => RenderObjectMeshEstimate {
            vertices: 50,
            faces: 20,
        },
        SaccharideShape::HexagonalPrism | SaccharideShape::FlatHexagon => {
            RenderObjectMeshEstimate {
                vertices: 60,
                faces: 24,
            }
        }
        SaccharideShape::HeptagonalPrism => RenderObjectMeshEstimate {
            vertices: 70,
            faces: 28,
        },
    }
}

#[cfg(test)]
fn render_objects_mesh_estimate<'a>(
    objects: impl Iterator<Item = &'a RenderObject>,
    options: &MeshOptions,
    cylinder_radial_segments: usize,
) -> RenderObjectMeshEstimate {
    objects.fold(RenderObjectMeshEstimate::default(), |total, object| {
        total.add(object.mesh_estimate(options, cylinder_radial_segments))
    })
}

fn render_objects_mesh_plan<'a>(
    objects: impl Iterator<Item = &'a RenderObject>,
    options: &MeshOptions,
    cylinder_radial_segments: usize,
) -> (RenderObjectMeshEstimate, Vec<RenderObjectMeshPlan>) {
    let mut total = RenderObjectMeshEstimate::default();
    let mut plans = Vec::new();
    for object in objects {
        let (estimate, dashed_samples) = match object {
            RenderObject::DashedTube { points, radius } if points.len() >= 2 => {
                let samples = sample_path(points, 8);
                (
                    dashed_tube_mesh_estimate_from_samples(
                        &samples,
                        *radius,
                        options.radial_segments.max(3),
                    ),
                    Some(samples),
                )
            }
            _ => (
                object.mesh_estimate(options, cylinder_radial_segments),
                None,
            ),
        };
        total = total.add(estimate);
        plans.push(RenderObjectMeshPlan {
            estimate,
            dashed_samples,
        });
    }
    (total, plans)
}

fn mesh_with_capacity(estimate: RenderObjectMeshEstimate) -> Mesh {
    Mesh {
        vertices: Vec::with_capacity(estimate.vertices),
        normals: Vec::with_capacity(estimate.vertices),
        faces: Vec::with_capacity(estimate.faces),
        vertex_groups: Vec::with_capacity(estimate.vertices),
        face_groups: Vec::with_capacity(estimate.faces),
        face_materials: Vec::new(),
        sections: Vec::new(),
        group_count: 0,
    }
}

pub(crate) fn build_semantic_render_objects(
    molecule: &Molecule,
    options: &MeshOptions,
) -> Vec<SemanticRenderObject> {
    let options = resolved_mesh_options(molecule, options);
    build_semantic_render_objects_resolved(molecule, &options)
}

fn build_semantic_render_objects_resolved(
    molecule: &Molecule,
    options: &MeshOptions,
) -> Vec<SemanticRenderObject> {
    build_semantic_render_objects_resolved_limited(molecule, options, None, None, |_| {})
}

fn build_semantic_render_objects_resolved_until_face_timed(
    molecule: &Molecule,
    options: &MeshOptions,
    face_index: usize,
    structure: Option<&AtomicStructure>,
    checkpoint: impl FnMut(&str),
) -> Vec<SemanticRenderObject> {
    build_semantic_render_objects_resolved_limited(
        molecule,
        options,
        Some(face_index),
        structure,
        checkpoint,
    )
}

fn build_semantic_render_objects_resolved_limited(
    molecule: &Molecule,
    options: &MeshOptions,
    target_face_index: Option<usize>,
    prebuilt_structure: Option<&AtomicStructure>,
    mut checkpoint: impl FnMut(&str),
) -> Vec<SemanticRenderObject> {
    let center = if options.center {
        bounds_molecule(molecule)
            .map(|(min, max)| Vec3 {
                x: (min.x + max.x) * 0.5,
                y: (min.y + max.y) * 0.5,
                z: (min.z + max.z) * 0.5,
            })
            .unwrap_or_default()
    } else {
        Vec3::default()
    };
    let mut objects = Vec::new();
    let mut group_id = 0usize;
    let representation = representation_name(options.representation);
    let polymer_trace_visual = polymer_trace_visual_name(options.representation);

    match options.representation {
        Representation::Molstar | Representation::Cartoon | Representation::Ribbon => {
            let structure_storage;
            let structure = match prebuilt_structure {
                Some(structure) => structure,
                None => {
                    structure_storage = molecule.atomic_structure();
                    &structure_storage
                }
            };
            checkpoint("atomic-structure-for-representation");
            let mut trace = if target_face_index.is_some() {
                backbone_residues_from_atoms(molecule)
            } else {
                backbone_residues(molecule, structure)
            };
            checkpoint("backbone-residues");
            if target_face_index.is_none() {
                apply_polymer_trace_terminal_flags(structure, &mut trace);
                apply_cyclic_polymer_trace_flags(structure, &mut trace);
                apply_polymer_trace_secondary_flags(structure, &mut trace);
            }
            checkpoint("polymer-trace-flags");
            if options.representation == Representation::Ribbon {
                let backbone: Vec<(String, i32, String, Vec3)> = trace
                    .iter()
                    .map(|residue| {
                        (
                            residue.chain.clone(),
                            residue.seq,
                            residue.insertion_code.clone(),
                            residue.position,
                        )
                    })
                    .collect();
                let mut covered = Vec::<(String, i32, String)>::new();
                for range in &molecule.helices {
                    let residues: Vec<&TraceResidue> = trace
                        .iter()
                        .filter(|residue| residue_in_secondary_range(residue, range))
                        .collect();
                    let points: Vec<Vec3> = residues.iter().map(|r| r.position - center).collect();
                    let directions: Vec<Option<Vec3>> =
                        residues.iter().map(|residue| residue.direction).collect();
                    if points.len() == 1 {
                        covered.extend(
                            residues
                                .iter()
                                .map(|r| (range.chain.clone(), r.seq, r.insertion_code.clone())),
                        );
                        push_semantic_with_group(
                            &mut objects,
                            trace_group_id_for_residues(&trace, &residues),
                            SemanticMeta::new(
                                representation,
                                "helix",
                                Some(&range.chain),
                                Some(range.start),
                                Some(range.end),
                            )
                            .with_trace_flags(secondary_trace_flags(
                                &trace,
                                &residues,
                                molecule,
                                SecondaryTraceKind::Helix,
                            ))
                            .with_visual(polymer_trace_visual),
                            RenderObject::Sphere {
                                center: points[0],
                                radius: (molstar_trace_radius(options) * 2.0) as f64,
                            },
                        );
                    } else if points.len() >= 2 {
                        covered.extend(
                            residues
                                .iter()
                                .map(|r| (range.chain.clone(), r.seq, r.insertion_code.clone())),
                        );
                        let (centers, normals) = helix_trace(&points, &directions);
                        let trace_flags = secondary_trace_flags(
                            &trace,
                            &residues,
                            molecule,
                            SecondaryTraceKind::Helix,
                        );
                        let (start_cap, end_cap) =
                            secondary_trace_cap_flags(structure, &residues, trace_flags);
                        let (width, thickness) = if options.tubular_helices
                            && options.representation != Representation::Ribbon
                        {
                            let radius =
                                molstar_trace_height(options) * MOLSTAR_TUBULAR_HELIX_FACTOR;
                            (radius, radius)
                        } else {
                            (
                                molstar_trace_radius(options),
                                if options.representation == Representation::Ribbon {
                                    molstar_trace_radius(options)
                                } else {
                                    molstar_trace_height(options)
                                },
                            )
                        };
                        push_semantic_with_group(
                            &mut objects,
                            trace_group_id_for_residues(&trace, &residues),
                            SemanticMeta::new(
                                representation,
                                "helix",
                                Some(&range.chain),
                                Some(range.start),
                                Some(range.end),
                            )
                            .with_trace_flags(trace_flags)
                            .with_visual(polymer_trace_visual),
                            RenderObject::OrientedRibbon {
                                centers,
                                normals,
                                width,
                                thickness,
                                profile: if options.tubular_helices
                                    && options.representation != Representation::Ribbon
                                {
                                    PolymerProfile::Elliptical
                                } else {
                                    options.helix_profile
                                },
                                start_cap,
                                end_cap,
                                round_cap: options.round_cap
                                    && options.tubular_helices
                                    && options.representation != Representation::Ribbon,
                            },
                        );
                    }
                }

                for range in &molecule.sheets {
                    let residues: Vec<&TraceResidue> = trace
                        .iter()
                        .filter(|residue| residue_in_secondary_range(residue, range))
                        .collect();
                    let points: Vec<Vec3> = residues
                        .iter()
                        .map(|residue| residue.position - center)
                        .collect();
                    if points.len() == 1 {
                        covered.extend(residues.iter().map(|residue| {
                            (
                                range.chain.clone(),
                                residue.seq,
                                residue.insertion_code.clone(),
                            )
                        }));
                        push_semantic_with_group(
                            &mut objects,
                            trace_group_id_for_residues(&trace, &residues),
                            SemanticMeta::new(
                                representation,
                                "sheet",
                                Some(&range.chain),
                                Some(range.start),
                                Some(range.end),
                            )
                            .with_trace_flags(secondary_trace_flags(
                                &trace,
                                &residues,
                                molecule,
                                SecondaryTraceKind::Sheet,
                            ))
                            .with_visual(polymer_trace_visual),
                            RenderObject::Sphere {
                                center: points[0],
                                radius: (molstar_trace_radius(options) * 2.0) as f64,
                            },
                        );
                    } else if points.len() >= 2 {
                        let width = molstar_trace_radius(options);
                        let height = molstar_trace_height(options);
                        let trace_flags = secondary_trace_flags(
                            &trace,
                            &residues,
                            molecule,
                            SecondaryTraceKind::Sheet,
                        );
                        let (start_cap, end_cap) =
                            secondary_trace_cap_flags(structure, &residues, trace_flags);
                        covered.extend(residues.iter().map(|residue| {
                            (
                                range.chain.clone(),
                                residue.seq,
                                residue.insertion_code.clone(),
                            )
                        }));
                        push_semantic_with_group(
                            &mut objects,
                            trace_group_id_for_residues(&trace, &residues),
                            SemanticMeta::new(
                                representation,
                                "sheet",
                                Some(&range.chain),
                                Some(range.start),
                                Some(range.end),
                            )
                            .with_trace_flags(trace_flags)
                            .with_visual(polymer_trace_visual),
                            RenderObject::Sheet {
                                points,
                                width,
                                thickness: height,
                                arrow_height: height * options.sheet_arrow_factor,
                                start_cap,
                                end_cap,
                            },
                        );
                    }
                }

                for segment in uncovered_backbone_segments(&backbone, &covered) {
                    let points: Vec<Vec3> =
                        segment.points.into_iter().map(|p| p - center).collect();
                    if points.len() == 1 {
                        push_semantic_with_group(
                            &mut objects,
                            trace_group_id_for_segment(
                                &trace,
                                &segment.chain,
                                segment.start,
                                &segment.start_insertion_code,
                            ),
                            SemanticMeta::new(
                                representation,
                                "coil",
                                Some(&segment.chain),
                                Some(segment.start),
                                Some(segment.end),
                            )
                            .with_trace_flags(trace_flags_for_segment(
                                &trace,
                                &segment.chain,
                                segment.start,
                                &segment.start_insertion_code,
                                segment.end,
                                &segment.end_insertion_code,
                                false,
                            ))
                            .with_visual(polymer_trace_visual),
                            RenderObject::Sphere {
                                center: points[0],
                                radius: (molstar_trace_radius(options) * 2.0) as f64,
                            },
                        );
                    } else if points.len() >= 2 {
                        let object = if options.representation == Representation::Ribbon {
                            RenderObject::Tube {
                                points,
                                radius: molstar_trace_radius(options),
                            }
                        } else if options.representation == Representation::Molstar {
                            RenderObject::DashedTube {
                                points,
                                radius: molstar_trace_radius(options),
                            }
                        } else {
                            RenderObject::Tube {
                                points,
                                radius: molstar_trace_radius(options),
                            }
                        };
                        push_semantic_with_group(
                            &mut objects,
                            trace_group_id_for_segment(
                                &trace,
                                &segment.chain,
                                segment.start,
                                &segment.start_insertion_code,
                            ),
                            SemanticMeta::new(
                                representation,
                                "coil",
                                Some(&segment.chain),
                                Some(segment.start),
                                Some(segment.end),
                            )
                            .with_trace_flags(trace_flags_for_segment(
                                &trace,
                                &segment.chain,
                                segment.start,
                                &segment.start_insertion_code,
                                segment.end,
                                &segment.end_insertion_code,
                                false,
                            ))
                            .with_visual(polymer_trace_visual),
                            object,
                        );
                    }
                }
            } else {
                add_polymer_trace_segment_semantic_objects(
                    options,
                    center,
                    representation,
                    polymer_trace_visual,
                    &trace,
                    &structure,
                    &mut objects,
                );
            }
            checkpoint("add-polymer-trace");
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
            let selected = selected_visuals(&structure, options);
            checkpoint("selected-visuals");

            add_nucleotide_semantic_objects(
                &trace,
                options,
                center,
                representation,
                &mut group_id,
                &mut objects,
                &selected,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
            add_polymer_gap_semantic_objects(
                molecule,
                structure,
                options,
                center,
                representation,
                &mut objects,
                &selected,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
            add_direction_wedge_semantic_objects(
                &trace,
                options,
                center,
                representation,
                &mut group_id,
                &mut objects,
                &selected,
                &structure,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
            add_carbohydrate_symbol_semantic_objects(
                molecule,
                structure,
                center,
                representation,
                &mut group_id,
                &mut objects,
                &selected,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
            add_carbohydrate_link_semantic_objects(
                molecule,
                structure,
                options,
                center,
                representation,
                &mut group_id,
                &mut objects,
                &selected,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
            add_carbohydrate_terminal_link_semantic_objects(
                molecule,
                structure,
                options,
                center,
                representation,
                &mut group_id,
                &mut objects,
                &selected,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }

            if options.representation == Representation::Molstar
                && selected_has_molstar_component_visuals(&selected)
            {
                let mut branched_mask = None::<Vec<bool>>;
                if has_ligand_component(structure) {
                    let branched = branched_mask
                        .get_or_insert_with(|| molstar_branched_atom_mask(molecule, structure));
                    let ligand_mask = molstar_ligand_atom_mask(molecule, structure, branched);
                    add_molstar_component_ball_and_stick_semantic_objects(
                        molecule,
                        options,
                        center,
                        representation,
                        "ligand",
                        &ligand_mask,
                        &selected,
                        &mut objects,
                    );
                    if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                        return objects;
                    }
                }
                if has_branched_component(structure) {
                    let branched = branched_mask
                        .get_or_insert_with(|| molstar_branched_atom_mask(molecule, structure));
                    add_molstar_component_ball_and_stick_semantic_objects(
                        molecule,
                        options,
                        center,
                        representation,
                        "branched",
                        branched,
                        &selected,
                        &mut objects,
                    );
                    if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                        return objects;
                    }
                }
                if has_ion_component(structure) {
                    let ion_mask = molstar_ion_atom_mask(structure);
                    add_molstar_component_ball_and_stick_semantic_objects(
                        molecule,
                        options,
                        center,
                        representation,
                        "ion",
                        &ion_mask,
                        &selected,
                        &mut objects,
                    );
                    if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                        return objects;
                    }
                }
                if has_water_component(structure) {
                    let water_mask = molstar_water_atom_mask(structure);
                    add_molstar_component_ball_and_stick_semantic_objects(
                        molecule,
                        options,
                        center,
                        representation,
                        "water",
                        &water_mask,
                        &selected,
                        &mut objects,
                    );
                    if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                        return objects;
                    }
                }
            } else if trace.is_empty() {
                add_ball_and_stick_semantic_objects(
                    molecule,
                    options,
                    center,
                    representation,
                    &mut group_id,
                    &mut objects,
                );
                if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                    return objects;
                }
            }
        }
        Representation::Backbone => {
            let structure_storage;
            let structure = match prebuilt_structure {
                Some(structure) => structure,
                None => {
                    structure_storage = molecule.atomic_structure();
                    &structure_storage
                }
            };
            let selected = selected_visuals(structure, options);
            checkpoint("selected-visuals");
            add_polymer_backbone_semantic_objects(
                molecule,
                structure,
                options,
                center,
                representation,
                &mut group_id,
                &mut objects,
                &selected,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
            add_polymer_gap_semantic_objects(
                molecule,
                structure,
                options,
                center,
                representation,
                &mut objects,
                &selected,
            );
            if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                return objects;
            }
        }
        Representation::Spacefill | Representation::BallAndStick => {
            if options.representation == Representation::Spacefill {
                group_id = 0;
                let chain_materials = molstar_chain_materials(molecule);
                for atom in &molecule.atoms {
                    push_semantic(
                        &mut objects,
                        &mut group_id,
                        SemanticMeta::new(
                            representation,
                            "atom",
                            Some(&atom.chain),
                            atom.residue_seq.parse::<i32>().ok(),
                            atom.residue_seq.parse::<i32>().ok(),
                        )
                        .with_visual("element-sphere")
                        .with_material(molstar_atom_material(atom, &chain_materials, "atom")),
                        RenderObject::Sphere {
                            center: atom.position - center,
                            radius: molstar_spacefill_atom_radius(atom, options),
                        },
                    );
                    if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                        return objects;
                    }
                }
            } else {
                add_ball_and_stick_semantic_objects(
                    molecule,
                    options,
                    center,
                    representation,
                    &mut group_id,
                    &mut objects,
                );
                if semantic_objects_cover_target_face(&objects, options, target_face_index) {
                    return objects;
                }
            }
        }
    }
    add_coarse_semantic_objects(
        molecule,
        center,
        representation,
        &mut group_id,
        &mut objects,
    );
    if semantic_objects_cover_target_face(&objects, options, target_face_index) {
        return objects;
    }
    apply_molstar_default_materials(&mut objects, molecule, options.color_theme);
    objects
}

fn semantic_objects_cover_target_face(
    objects: &[SemanticRenderObject],
    options: &MeshOptions,
    target_face_index: Option<usize>,
) -> bool {
    let Some(target_face_index) = target_face_index else {
        return false;
    };
    target_face_index < semantic_objects_face_count(objects, options)
}

fn semantic_objects_face_count(objects: &[SemanticRenderObject], options: &MeshOptions) -> usize {
    objects
        .iter()
        .map(|object| object.object.face_estimate(options))
        .sum()
}

fn push_semantic(
    objects: &mut Vec<SemanticRenderObject>,
    group_id: &mut usize,
    meta: SemanticMeta<'_>,
    object: RenderObject,
) {
    push_semantic_with_group(objects, *group_id, meta, object);
    *group_id += 1;
}

fn push_semantic_with_group(
    objects: &mut Vec<SemanticRenderObject>,
    group_id: usize,
    meta: SemanticMeta<'_>,
    object: RenderObject,
) {
    let geometry_type = geometry_type(&object);
    objects.push(SemanticRenderObject {
        geometry_type,
        visual: meta.visual.unwrap_or(geometry_type),
        representation: meta.representation,
        secondary_type: meta.secondary_type,
        chain: meta.chain.map(str::to_string),
        residue_start: meta.residue_start,
        residue_end: meta.residue_end,
        group_id,
        initial: meta.trace_flags.initial,
        final_residue: meta.trace_flags.final_residue,
        sec_struc_first: meta.trace_flags.sec_struc_first,
        sec_struc_last: meta.trace_flags.sec_struc_last,
        material: meta.material,
        object,
    });
}

const MOLSTAR_MANY_DISTINCT_COLORS: [u32; 25] = [
    0x1b9e77, 0xd95f02, 0x7570b3, 0xe7298a, 0x66a61e, 0xe6ab02, 0xa6761d, 0x666666, 0xe41a1c,
    0x377eb8, 0x4daf4a, 0x984ea3, 0xff7f00, 0xffff33, 0xa65628, 0xf781bf, 0x999999, 0x66c2a5,
    0xfc8d62, 0x8da0cb, 0xe78ac3, 0xa6d854, 0xffd92f, 0xe5c494, 0xb3b3b3,
];

fn apply_molstar_default_materials(
    objects: &mut [SemanticRenderObject],
    molecule: &Molecule,
    color_theme: ColorTheme,
) {
    match color_theme {
        ColorTheme::ChainId => {
            let chain_materials = molstar_chain_materials(molecule);
            for object in objects {
                if object.material.is_none() {
                    object.material = Some(molstar_chain_material_for_key(
                        object.chain.as_deref(),
                        &chain_materials,
                    ));
                }
            }
        }
    }
}

fn molstar_chain_materials(molecule: &Molecule) -> BTreeMap<String, MeshMaterial> {
    let mut keys = Vec::<String>::new();
    for atom in &molecule.atoms {
        let key = molstar_atom_chain_key(atom);
        if !keys.iter().any(|existing| existing == &key) {
            keys.push(key);
        }
    }
    for sphere in &molecule.coarse_spheres {
        if !keys.iter().any(|existing| existing == &sphere.asym_id) {
            keys.push(sphere.asym_id.clone());
        }
    }
    for gaussian in &molecule.coarse_gaussians {
        if !keys.iter().any(|existing| existing == &gaussian.asym_id) {
            keys.push(gaussian.asym_id.clone());
        }
    }

    keys.into_iter()
        .enumerate()
        .map(|(index, key)| {
            (
                key,
                MeshMaterial::opaque(
                    MOLSTAR_MANY_DISTINCT_COLORS[index % MOLSTAR_MANY_DISTINCT_COLORS.len()],
                ),
            )
        })
        .collect()
}

fn molstar_atom_chain_key(atom: &crate::model::Atom) -> String {
    if atom.auth_chain.trim().is_empty() {
        atom.chain.clone()
    } else {
        atom.auth_chain.clone()
    }
}

fn molstar_chain_material_for_key(
    key: Option<&str>,
    chain_materials: &BTreeMap<String, MeshMaterial>,
) -> MeshMaterial {
    key.and_then(|key| chain_materials.get(key))
        .copied()
        .unwrap_or_else(|| {
            chain_materials
                .values()
                .next()
                .copied()
                .unwrap_or_else(|| MeshMaterial::opaque(0xfafafa))
        })
}

fn molstar_atom_material(
    atom: &crate::model::Atom,
    chain_materials: &BTreeMap<String, MeshMaterial>,
    component: &str,
) -> MeshMaterial {
    let alpha_tenths = match component {
        "branched" => 3,
        "water" | "lipid" => 6,
        _ => 10,
    };
    let color = if molstar_atom_uses_chain_carbon_color(atom, component) {
        molstar_chain_material_for_key(Some(&molstar_atom_chain_key(atom)), chain_materials).color
    } else {
        molstar_element_symbol_color(atom)
    };
    MeshMaterial::with_alpha_tenths(color, alpha_tenths)
}

fn molstar_atom_uses_chain_carbon_color(atom: &crate::model::Atom, component: &str) -> bool {
    molstar_atom_element_symbol(atom) == "C" && !matches!(component, "water" | "ion" | "lipid")
}

fn molstar_element_symbol_color(atom: &crate::model::Atom) -> u32 {
    let base = match molstar_atom_element_symbol(atom).as_str() {
        "H" => 0xffffff,
        "D" => 0xffffc0,
        "T" => 0xffffa0,
        "C" => 0x909090,
        "N" => 0x3050f8,
        "O" => 0xff0d0d,
        "F" => 0x90e050,
        "P" => 0xff8000,
        "S" => 0xffff30,
        "CL" => 0x1ff01f,
        "BR" => 0xa62929,
        "I" => 0x940094,
        "NA" => 0xab5cf2,
        "MG" => 0x8aff00,
        "K" => 0x8f40d4,
        "CA" => 0x3dff00,
        "MN" => 0x9c7ac7,
        "FE" => 0xe06633,
        "CO" => 0xf090a0,
        "NI" => 0x50d050,
        "CU" => 0xc88033,
        "ZN" => 0x7d80b0,
        "SE" => 0xffa100,
        _ => 0xffffff,
    };
    molstar_adjust_element_symbol_color(base)
}

fn molstar_atom_element_symbol(atom: &crate::model::Atom) -> String {
    let symbol = if atom.type_symbol.trim().is_empty() {
        atom.element.trim()
    } else {
        atom.type_symbol.trim()
    };
    symbol.to_ascii_uppercase()
}

fn molstar_adjust_element_symbol_color(color: u32) -> u32 {
    let mut lab = molstar_color_to_lab(color);
    // ElementSymbolColorTheme default lightness is 0.2; Color.darken(c, -0.2)
    // raises CIE Lab L by 18 * 0.2.
    lab[0] += 18.0 * 0.2;
    molstar_lab_to_color(lab)
}

fn molstar_color_to_lab(color: u32) -> [f64; 3] {
    const XN: f64 = 0.950470;
    const YN: f64 = 1.0;
    const ZN: f64 = 1.088830;
    let r = molstar_rgb_xyz(((color >> 16) & 0xff) as f64);
    let g = molstar_rgb_xyz(((color >> 8) & 0xff) as f64);
    let b = molstar_rgb_xyz((color & 0xff) as f64);
    let x = molstar_xyz_lab((0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / XN);
    let y = molstar_xyz_lab((0.2126729 * r + 0.7151522 * g + 0.0721750 * b) / YN);
    let z = molstar_xyz_lab((0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / ZN);
    let l = 116.0 * y - 16.0;
    [
        if l < 0.0 { 0.0 } else { l },
        500.0 * (x - y),
        200.0 * (y - z),
    ]
}

fn molstar_lab_to_color(lab: [f64; 3]) -> u32 {
    const XN: f64 = 0.950470;
    const YN: f64 = 1.0;
    const ZN: f64 = 1.088830;
    let mut y = (lab[0] + 16.0) / 116.0;
    let mut x = if lab[1].is_nan() {
        y
    } else {
        y + lab[1] / 500.0
    };
    let mut z = if lab[2].is_nan() {
        y
    } else {
        y - lab[2] / 200.0
    };
    y = YN * molstar_lab_xyz(y);
    x = XN * molstar_lab_xyz(x);
    z = ZN * molstar_lab_xyz(z);
    let r = molstar_xyz_rgb(3.2404542 * x - 1.5371385 * y - 0.4985314 * z);
    let g = molstar_xyz_rgb(-0.9692660 * x + 1.8760108 * y + 0.0415560 * z);
    let b = molstar_xyz_rgb(0.0556434 * x - 0.2040259 * y + 1.0572252 * z);
    (molstar_round_u8(r) << 16) | (molstar_round_u8(g) << 8) | molstar_round_u8(b)
}

fn molstar_rgb_xyz(mut c: f64) -> f64 {
    c /= 255.0;
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

fn molstar_xyz_lab(t: f64) -> f64 {
    const T0: f64 = 0.137931034;
    const T2: f64 = 0.12841855;
    const T3: f64 = 0.008856452;
    if t > T3 {
        t.powf(1.0 / 3.0)
    } else {
        t / T2 + T0
    }
}

fn molstar_lab_xyz(t: f64) -> f64 {
    const T0: f64 = 0.137931034;
    const T1: f64 = 0.206896552;
    const T2: f64 = 0.12841855;
    if t > T1 {
        t * t * t
    } else {
        T2 * (t - T0)
    }
}

fn molstar_xyz_rgb(c: f64) -> f64 {
    255.0
        * if c <= 0.00304 {
            12.92 * c
        } else {
            1.055 * c.powf(1.0 / 2.4) - 0.055
        }
}

fn molstar_round_u8(value: f64) -> u32 {
    (value.clamp(0.0, 255.0) + 0.5).floor() as u32
}

fn trace_group_id_for_residues(trace: &[TraceResidue], residues: &[&TraceResidue]) -> usize {
    residues
        .first()
        .map(|residue| {
            trace_group_id_for_segment(trace, &residue.chain, residue.seq, &residue.insertion_code)
        })
        .unwrap_or(0)
}

fn trace_group_id_for_segment(
    trace: &[TraceResidue],
    chain: &str,
    seq: i32,
    insertion_code: &str,
) -> usize {
    trace
        .iter()
        .position(|residue| {
            residue.chain == chain && residue.seq == seq && residue.insertion_code == insertion_code
        })
        .unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
fn add_polymer_trace_segment_semantic_objects(
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    polymer_trace_visual: &'static str,
    trace: &[TraceResidue],
    structure: &AtomicStructure,
    objects: &mut Vec<SemanticRenderObject>,
) {
    let hierarchy = &structure.model.hierarchy;
    let use_helix_orientation =
        options.tubular_helices && options.representation != Representation::Ribbon;
    let object_start = objects.len();

    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        let Some(start_residue) = hierarchy.residue_atom_segments.index.get(pair[0]).copied()
        else {
            continue;
        };
        let Some(end_residue) = hierarchy.residue_atom_segments.index.get(pair[1]).copied() else {
            continue;
        };

        for residue_index in start_residue..=end_residue {
            let Some(trace_index) =
                trace_residue_index_for_model_residue(hierarchy, trace, residue_index)
            else {
                continue;
            };
            let residue = &trace[trace_index];
            let current_type = molstar_secondary_trace_type(structure, residue_index);
            let previous_residue = polymer_trace_residue_index(
                structure,
                start_residue,
                end_residue,
                residue_index as isize - 1,
            );
            let next_residue = polymer_trace_residue_index(
                structure,
                start_residue,
                end_residue,
                residue_index as isize + 1,
            );
            let w0 = polymer_trace_radius(structure, previous_residue, options);
            let w1 = polymer_trace_radius(structure, residue_index, options);
            let w2 = polymer_trace_radius(structure, next_residue, options);
            let previous_type = molstar_secondary_trace_type(structure, previous_residue);
            let next_type = molstar_secondary_trace_type(structure, next_residue);
            let sec_struc_first = previous_type != current_type;
            let sec_struc_last = current_type != next_type;
            let current_coarse_backbone = polymer_trace_coarse_backbone(structure, residue_index);
            let start_cap = sec_struc_first
                || polymer_trace_coarse_backbone(structure, previous_residue)
                    != current_coarse_backbone
                || residue.initial;
            let end_cap = sec_struc_last
                || current_coarse_backbone
                    != polymer_trace_coarse_backbone(structure, next_residue)
                || residue.final_residue;
            let state = polymer_trace_iterator_state(
                structure,
                start_residue,
                end_residue,
                residue_index,
                current_type,
                use_helix_orientation,
            );
            let center = DVec3::from_vec3(center);
            let controls = geometry::CurveSegmentControls {
                sec_struc_first,
                sec_struc_last,
                p0: state.p0 - center,
                p1: state.p1 - center,
                p2: state.p2 - center,
                p3: state.p3 - center,
                p4: state.p4 - center,
                d12: state.d12,
                d23: state.d23,
            };
            let is_sheet = current_type.contains(SecondaryStructureType::BETA);
            let is_helix = is_helix_secondary(current_type);
            let secondary_type = if is_sheet {
                "sheet"
            } else if is_helix {
                "helix"
            } else {
                "coil"
            };
            let tension = if is_helix && !options.tubular_helices {
                MOLSTAR_HELIX_TENSION
            } else {
                MOLSTAR_STANDARD_TENSION
            };
            let shift = if residue.is_nucleotide {
                MOLSTAR_NUCLEIC_BACKBONE_SHIFT
            } else {
                MOLSTAR_STANDARD_BACKBONE_SHIFT
            };
            let trace_flags = TraceFlags {
                initial: residue.initial,
                final_residue: residue.final_residue,
                sec_struc_first,
                sec_struc_last,
            };
            let meta = SemanticMeta::new(
                representation,
                secondary_type,
                Some(&residue.chain),
                Some(residue.seq),
                Some(residue.seq),
            )
            .with_trace_flags(trace_flags)
            .with_visual(polymer_trace_visual);

            if residue.initial && residue.final_residue {
                push_semantic_with_group(
                    objects,
                    trace_index,
                    meta,
                    RenderObject::Sphere {
                        center: controls.p2.to_vec3(),
                        radius: (w1 * 2.0) as f64,
                    },
                );
                continue;
            }

            let (widths, heights, kind, swap_normal_binormal) = if is_sheet {
                let h0 = w0 * MOLSTAR_TRACE_ASPECT_RATIO;
                let h1 = w1 * MOLSTAR_TRACE_ASPECT_RATIO;
                let h2 = w2 * MOLSTAR_TRACE_ASPECT_RATIO;
                (
                    [w0, w1, w2],
                    [h0, h1, h2],
                    if options.radial_segments == 2 {
                        PolymerTraceSegmentKind::Ribbon {
                            arrow_height: if sec_struc_last {
                                h1 * options.sheet_arrow_factor
                            } else {
                                0.0
                            },
                            swap_width_height: false,
                        }
                    } else {
                        PolymerTraceSegmentKind::Sheet {
                            arrow_height: if sec_struc_last {
                                h1 * options.sheet_arrow_factor
                            } else {
                                0.0
                            },
                        }
                    },
                    false,
                )
            } else {
                let (widths, heights) = if is_helix {
                    if options.tubular_helices && options.representation != Representation::Ribbon {
                        let factor = MOLSTAR_TRACE_ASPECT_RATIO * MOLSTAR_TUBULAR_HELIX_FACTOR;
                        let widths = [w0 * factor, w1 * factor, w2 * factor];
                        (widths, widths)
                    } else {
                        (
                            [w0, w1, w2],
                            [
                                w0 * MOLSTAR_TRACE_ASPECT_RATIO,
                                w1 * MOLSTAR_TRACE_ASPECT_RATIO,
                                w2 * MOLSTAR_TRACE_ASPECT_RATIO,
                            ],
                        )
                    }
                } else if residue.is_nucleotide {
                    (
                        [w0, w1, w2],
                        [
                            w0 * MOLSTAR_TRACE_ASPECT_RATIO,
                            w1 * MOLSTAR_TRACE_ASPECT_RATIO,
                            w2 * MOLSTAR_TRACE_ASPECT_RATIO,
                        ],
                    )
                } else {
                    ([w0, w1, w2], [w0, w1, w2])
                };
                let profile = if is_helix
                    && options.tubular_helices
                    && options.representation != Representation::Ribbon
                {
                    PolymerProfile::Elliptical
                } else if residue.is_nucleotide {
                    PolymerProfile::Square
                } else {
                    options.helix_profile
                };
                let kind = if options.radial_segments == 2 {
                    PolymerTraceSegmentKind::Ribbon {
                        arrow_height: 0.0,
                        swap_width_height: residue.is_nucleotide,
                    }
                } else if options.radial_segments == 4 {
                    PolymerTraceSegmentKind::Sheet { arrow_height: 0.0 }
                } else if widths[1] == heights[1] {
                    PolymerTraceSegmentKind::Tube {
                        profile: PolymerProfile::Elliptical,
                        round_cap: options.round_cap
                            && is_helix
                            && options.tubular_helices
                            && options.representation != Representation::Ribbon,
                    }
                } else if profile == PolymerProfile::Square {
                    PolymerTraceSegmentKind::Sheet { arrow_height: 0.0 }
                } else {
                    PolymerTraceSegmentKind::Tube {
                        profile,
                        round_cap: options.round_cap
                            && is_helix
                            && options.tubular_helices
                            && options.representation != Representation::Ribbon,
                    }
                };
                (widths, heights, kind, residue.is_nucleotide)
            };

            push_semantic_with_group(
                objects,
                trace_index,
                meta,
                RenderObject::PolymerTraceSegment {
                    controls,
                    widths,
                    heights,
                    tension,
                    shift,
                    overhang_width: w1,
                    kind,
                    start_cap,
                    end_cap,
                    initial: residue.initial,
                    final_residue: residue.final_residue,
                    swap_normal_binormal,
                },
            );
        }
    }

    if objects.len() == object_start && trace.len() == 1 {
        let residue = &trace[0];
        let radius = model_residue_index_for_trace_residue(hierarchy, residue)
            .map(|residue_index| polymer_trace_radius(structure, residue_index, options))
            .unwrap_or_else(|| molstar_trace_radius(options));
        push_semantic_with_group(
            objects,
            0,
            SemanticMeta::new(
                representation,
                "coil",
                Some(&residue.chain),
                Some(residue.seq),
                Some(residue.seq),
            )
            .with_trace_flags(trace_flags_from_residues(&[residue]))
            .with_visual(polymer_trace_visual),
            RenderObject::Sphere {
                center: residue.position - center,
                radius: (radius * 2.0) as f64,
            },
        );
    }
}

fn add_polymer_gap_semantic_objects(
    molecule: &Molecule,
    structure: &AtomicStructure,
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    objects: &mut Vec<SemanticRenderObject>,
    selected: &[String],
) {
    if !selected.iter().any(|visual| visual == "polymer-gap") {
        return;
    }

    let mut group_id = 0usize;
    for unit in &structure.units {
        if unit.kind != crate::model::UnitKind::Atomic {
            continue;
        }
        for gap in unit.props.gap_elements.chunks_exact(2) {
            let Some(atom_a) = molecule.atoms.get(gap[0]) else {
                continue;
            };
            let Some(atom_b) = molecule.atoms.get(gap[1]) else {
                continue;
            };
            if atom_a.position.distance(atom_b.position) <= 0.000_001 {
                continue;
            }

            let start = atom_a.position - center;
            let end = atom_b.position - center;
            let radius_a = polymer_gap_radius(atom_a, options);
            let radius_b = polymer_gap_radius(atom_b, options);
            let residue_start = atom_a.residue_seq.parse::<i32>().ok();
            let residue_end = atom_b.residue_seq.parse::<i32>().ok();

            push_semantic_with_group(
                objects,
                group_id,
                SemanticMeta::new(
                    representation,
                    "gap",
                    Some(&atom_a.chain),
                    residue_start,
                    residue_end,
                )
                .with_visual("polymer-gap"),
                RenderObject::FixedCountDashedCylinder {
                    start,
                    end,
                    radius: radius_a,
                    length_scale: 0.5,
                    segment_count: 10,
                },
            );
            push_semantic_with_group(
                objects,
                group_id + 1,
                SemanticMeta::new(
                    representation,
                    "gap",
                    Some(&atom_b.chain),
                    residue_end,
                    residue_start,
                )
                .with_visual("polymer-gap"),
                RenderObject::FixedCountDashedCylinder {
                    start: end,
                    end: start,
                    radius: radius_b,
                    length_scale: 0.5,
                    segment_count: 10,
                },
            );
            group_id += 2;
        }
    }
}

fn polymer_gap_radius(atom: &crate::model::Atom, options: &MeshOptions) -> f32 {
    let _ = atom;
    molstar_cartoon_uniform_trace_radius(options)
}

fn bounds_molecule(molecule: &Molecule) -> Option<(Vec3, Vec3)> {
    let mut points = molecule
        .atoms
        .iter()
        .map(|atom| atom.position)
        .collect::<Vec<_>>();
    for sphere in &molecule.coarse_spheres {
        let radius = Vec3::new(sphere.radius, sphere.radius, sphere.radius);
        points.push(sphere.position - radius);
        points.push(sphere.position + radius);
    }
    for gaussian in &molecule.coarse_gaussians {
        let axes = gaussian_axes(gaussian.covariance, gaussian.weight);
        let extent = axes.iter().fold(Vec3::default(), |acc, axis| {
            acc + Vec3::new(axis.x.abs(), axis.y.abs(), axis.z.abs())
        });
        points.push(gaussian.position - extent);
        points.push(gaussian.position + extent);
    }
    let first = points.first().copied()?;
    let mut min = first;
    let mut max = first;
    for point in &points[1..] {
        min = min.min(*point);
        max = max.max(*point);
    }
    Some((min, max))
}

fn representation_name(representation: Representation) -> &'static str {
    match representation {
        Representation::Molstar => "molstar",
        Representation::Spacefill => "spacefill",
        Representation::BallAndStick => "ball-and-stick",
        Representation::Cartoon => "cartoon",
        Representation::Ribbon => "ribbon",
        Representation::Backbone => "backbone",
    }
}

fn polymer_trace_visual_name(representation: Representation) -> &'static str {
    match representation {
        Representation::Ribbon => "polymer-tube",
        Representation::Backbone => "polymer-backbone-cylinder",
        _ => "polymer-trace",
    }
}

const SMALL_STRUCTURE_RESIDUE_COUNT: usize = 10;
const MEDIUM_STRUCTURE_RESIDUE_COUNT: usize = 5_000;
const LARGE_STRUCTURE_RESIDUE_COUNT: usize = 30_000;
const HIGH_SYMMETRY_UNIT_COUNT: usize = 10;

fn selected_visuals(structure: &AtomicStructure, options: &MeshOptions) -> Vec<String> {
    if !options.visuals.is_empty() {
        let allowed = match options.representation {
            Representation::Spacefill => &["element-sphere", "structure-element-sphere"][..],
            Representation::BallAndStick => &[
                "element-sphere",
                "intra-bond",
                "inter-bond",
                "structure-element-sphere",
                "structure-intra-bond",
            ][..],
            Representation::Molstar | Representation::Cartoon => &[
                "polymer-trace",
                "polymer-gap",
                "element-sphere",
                "intra-bond",
                "inter-bond",
                "structure-element-sphere",
                "structure-intra-bond",
                "nucleotide-block",
                "nucleotide-ring",
                "nucleotide-atomic-ring-fill",
                "nucleotide-atomic-bond",
                "nucleotide-atomic-element",
                "direction-wedge",
                "carbohydrate-symbol",
                "carbohydrate-link",
                "carbohydrate-terminal-link",
            ][..],
            Representation::Ribbon => &["polymer-tube", "polymer-gap"][..],
            Representation::Backbone => &[
                "polymer-backbone-cylinder",
                "polymer-backbone-sphere",
                "polymer-gap",
            ][..],
        };
        let mut selected = Vec::new();
        for visual in &options.visuals {
            if allowed.iter().any(|allowed| allowed == visual)
                && !selected.iter().any(|selected| selected == visual)
            {
                selected.push(visual.clone());
            }
        }
        return selected;
    }
    match options.representation {
        Representation::Spacefill => {
            if structure.symmetry_groups.len() > 5_000 {
                vec!["structure-element-sphere".to_string()]
            } else {
                vec!["element-sphere".to_string()]
            }
        }
        Representation::BallAndStick => {
            if molstar_structure_size(structure) >= MolstarStructureSize::Huge {
                vec!["element-sphere".to_string(), "intra-bond".to_string()]
            } else if structure.symmetry_groups.len() > 5_000 {
                vec![
                    "structure-element-sphere".to_string(),
                    "structure-intra-bond".to_string(),
                ]
            } else {
                vec![
                    "element-sphere".to_string(),
                    "intra-bond".to_string(),
                    "inter-bond".to_string(),
                ]
            }
        }
        Representation::Molstar => molstar_selected_visuals(structure),
        Representation::Cartoon => cartoon_selected_visuals(structure),
        Representation::Ribbon => ribbon_selected_visuals(structure),
        Representation::Backbone => backbone_selected_visuals(structure),
    }
}

fn selected_has_molstar_component_visuals(selected: &[String]) -> bool {
    selected.iter().any(|visual| {
        matches!(
            visual.as_str(),
            "element-sphere" | "structure-element-sphere" | "intra-bond" | "structure-intra-bond"
        )
    })
}

fn molstar_selected_visuals(structure: &AtomicStructure) -> Vec<String> {
    let mut visuals = cartoon_selected_visuals(structure);
    if has_ligand_component(structure) {
        append_visuals(&mut visuals, &ball_and_stick_default_visuals(structure));
    }
    if has_branched_component(structure) {
        append_visuals(&mut visuals, &ball_and_stick_default_visuals(structure));
    }
    if has_ion_component(structure) {
        append_visuals(&mut visuals, &ball_and_stick_default_visuals(structure));
    }
    if has_water_component(structure) {
        append_visuals(&mut visuals, &ball_and_stick_default_visuals(structure));
    }
    visuals
}

fn append_visuals(visuals: &mut Vec<String>, additions: &[String]) {
    for visual in additions {
        if !visuals.iter().any(|existing| existing == visual) {
            visuals.push(visual.clone());
        }
    }
}

fn ball_and_stick_default_visuals(structure: &AtomicStructure) -> Vec<String> {
    if molstar_structure_size(structure) >= MolstarStructureSize::Huge {
        vec!["element-sphere".to_string(), "intra-bond".to_string()]
    } else if structure.symmetry_groups.len() > 5_000 {
        vec![
            "structure-element-sphere".to_string(),
            "structure-intra-bond".to_string(),
        ]
    } else {
        vec![
            "element-sphere".to_string(),
            "intra-bond".to_string(),
            "inter-bond".to_string(),
        ]
    }
}

fn cartoon_selected_visuals(structure: &AtomicStructure) -> Vec<String> {
    let mut visuals = Vec::new();
    if structure.polymer_residue_count > 0 {
        visuals.push("polymer-trace".to_string());
    }
    if has_nucleotides(structure) {
        visuals.push("nucleotide-ring".to_string());
    }
    if structure.polymer_gap_count > 0 {
        visuals.push("polymer-gap".to_string());
    }
    visuals
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum MolstarStructureSize {
    Small,
    Medium,
    Large,
    Huge,
    Gigantic,
}

fn molstar_structure_size(structure: &AtomicStructure) -> MolstarStructureSize {
    let residue_count = structure.polymer_residue_count;
    if residue_count >= LARGE_STRUCTURE_RESIDUE_COUNT {
        if structure
            .symmetry_groups
            .iter()
            .any(|group| group.unit_ids.len() > HIGH_SYMMETRY_UNIT_COUNT)
        {
            MolstarStructureSize::Huge
        } else {
            MolstarStructureSize::Gigantic
        }
    } else if residue_count < SMALL_STRUCTURE_RESIDUE_COUNT {
        MolstarStructureSize::Small
    } else if residue_count < MEDIUM_STRUCTURE_RESIDUE_COUNT {
        MolstarStructureSize::Medium
    } else {
        MolstarStructureSize::Large
    }
}

fn ribbon_selected_visuals(structure: &AtomicStructure) -> Vec<String> {
    let mut visuals = Vec::new();
    if structure.polymer_residue_count > 0 {
        visuals.push("polymer-tube".to_string());
    }
    if structure.polymer_gap_count > 0 {
        visuals.push("polymer-gap".to_string());
    }
    visuals
}

fn backbone_selected_visuals(structure: &AtomicStructure) -> Vec<String> {
    let mut visuals = Vec::new();
    if structure.polymer_residue_count > 0 {
        visuals.push("polymer-backbone-cylinder".to_string());
        visuals.push("polymer-backbone-sphere".to_string());
    }
    if structure.polymer_gap_count > 0 {
        visuals.push("polymer-gap".to_string());
    }
    visuals
}

fn has_nucleotides(structure: &AtomicStructure) -> bool {
    structure
        .units
        .iter()
        .any(|unit| !unit.props.nucleotide_elements.is_empty())
}

fn has_ligand_component(structure: &AtomicStructure) -> bool {
    structure
        .model
        .hierarchy
        .atoms
        .iter()
        .enumerate()
        .any(|(atom_index, _)| atom_is_molstar_ligand(structure, atom_index))
}

fn has_branched_component(structure: &AtomicStructure) -> bool {
    structure
        .model
        .hierarchy
        .atoms
        .iter()
        .enumerate()
        .any(|(atom_index, _)| atom_is_molstar_branched(structure, atom_index))
}

fn has_ion_component(structure: &AtomicStructure) -> bool {
    structure
        .model
        .hierarchy
        .derived
        .residue
        .molecule_type
        .contains(&MoleculeType::Ion)
}

fn has_water_component(structure: &AtomicStructure) -> bool {
    structure
        .model
        .hierarchy
        .derived
        .residue
        .molecule_type
        .contains(&MoleculeType::Water)
}

fn molstar_ligand_atom_mask(
    molecule: &Molecule,
    structure: &AtomicStructure,
    branched_mask: &[bool],
) -> Vec<bool> {
    let direct = (0..molecule.atoms.len())
        .map(|atom_index| atom_is_molstar_ligand(&structure, atom_index))
        .collect::<Vec<_>>();
    expand_mask_to_connected_whole_residues(molecule, structure, &direct, Some(branched_mask))
}

fn molstar_branched_atom_mask(molecule: &Molecule, structure: &AtomicStructure) -> Vec<bool> {
    let direct = (0..molecule.atoms.len())
        .map(|atom_index| atom_is_molstar_branched(&structure, atom_index))
        .collect::<Vec<_>>();
    expand_mask_to_connected_whole_residues(molecule, structure, &direct, None)
}

fn molstar_ion_atom_mask(structure: &AtomicStructure) -> Vec<bool> {
    (0..structure.model.hierarchy.atoms.len())
        .map(|atom_index| {
            atom_molecule_type(structure, atom_index)
                .is_some_and(|molecule_type| molecule_type == MoleculeType::Ion)
        })
        .collect()
}

fn molstar_water_atom_mask(structure: &AtomicStructure) -> Vec<bool> {
    (0..structure.model.hierarchy.atoms.len())
        .map(|atom_index| {
            atom_molecule_type(structure, atom_index)
                .is_some_and(|molecule_type| molecule_type == MoleculeType::Water)
        })
        .collect()
}

fn atom_is_molstar_ligand(structure: &AtomicStructure, atom_index: usize) -> bool {
    if atom_is_molstar_branched(structure, atom_index) {
        return false;
    }
    let Some(atom) = structure.model.hierarchy.atoms.get(atom_index) else {
        return false;
    };
    let Some(residue) = structure.model.hierarchy.residues.get(atom.residue_index) else {
        return false;
    };
    let entity_type = structure
        .model
        .hierarchy
        .index
        .entity_type_from_chain(residue.chain_index)
        .unwrap_or_default();
    let molecule_type = atom_molecule_type(structure, atom_index).unwrap_or_default();
    if matches!(
        molecule_type,
        MoleculeType::Water
            | MoleculeType::Ion
            | MoleculeType::Lipid
            | MoleculeType::Protein
            | MoleculeType::Rna
            | MoleculeType::Dna
            | MoleculeType::Pna
            | MoleculeType::Saccharide
    ) {
        return false;
    }
    residue.is_het || entity_type == "non-polymer" || molecule_type == MoleculeType::Other
}

fn atom_is_molstar_branched(structure: &AtomicStructure, atom_index: usize) -> bool {
    let Some(atom) = structure.model.hierarchy.atoms.get(atom_index) else {
        return false;
    };
    let Some(residue) = structure.model.hierarchy.residues.get(atom.residue_index) else {
        return false;
    };
    let chain_entity_type = structure
        .model
        .hierarchy
        .index
        .entity_type_from_chain(residue.chain_index)
        .unwrap_or_default();
    chain_entity_type == "branched"
        || atom_molecule_type(structure, atom_index)
            .is_some_and(|molecule_type| molecule_type == MoleculeType::Saccharide)
}

fn atom_molecule_type(structure: &AtomicStructure, atom_index: usize) -> Option<MoleculeType> {
    let residue_index = structure
        .model
        .hierarchy
        .atoms
        .get(atom_index)?
        .residue_index;
    structure
        .model
        .hierarchy
        .derived
        .residue
        .molecule_type
        .get(residue_index)
        .copied()
}

fn expand_mask_to_connected_whole_residues(
    molecule: &Molecule,
    structure: &AtomicStructure,
    mask: &[bool],
    excluded: Option<&[bool]>,
) -> Vec<bool> {
    let mut expanded = mask.to_vec();
    let mut connected_residues = vec![false; structure.model.hierarchy.residues.len()];

    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        let a_selected = mask.get(bond.a).copied().unwrap_or(false);
        let b_selected = mask.get(bond.b).copied().unwrap_or(false);
        if a_selected == b_selected || !bond_allows_connected_component(molecule, bond_index) {
            continue;
        }
        let connected_atom_index = if a_selected { bond.b } else { bond.a };
        if excluded
            .and_then(|mask| mask.get(connected_atom_index))
            .copied()
            .unwrap_or(false)
        {
            continue;
        }
        if let Some(residue_index) = structure
            .model
            .hierarchy
            .atoms
            .get(connected_atom_index)
            .map(|atom| atom.residue_index)
        {
            connected_residues[residue_index] = true;
        }
    }

    for (atom_index, atom) in structure.model.hierarchy.atoms.iter().enumerate() {
        if !connected_residues
            .get(atom.residue_index)
            .copied()
            .unwrap_or(false)
        {
            continue;
        }
        if excluded
            .and_then(|mask| mask.get(atom_index))
            .copied()
            .unwrap_or(false)
        {
            continue;
        }
        if let Some(slot) = expanded.get_mut(atom_index) {
            *slot = true;
        }
    }
    expanded
}

fn bond_allows_connected_component(molecule: &Molecule, bond_index: usize) -> bool {
    let flags = molecule
        .bond_metadata
        .get(bond_index)
        .map(|metadata| metadata.flags)
        .unwrap_or_else(|| {
            let Some(bond) = molecule.bonds.get(bond_index) else {
                return BondFlags::NONE;
            };
            molecule
                .atoms
                .get(bond.a)
                .zip(molecule.atoms.get(bond.b))
                .map(|(a, b)| crate::model::BondMetadata::computed_for_atoms(a, b).flags)
                .unwrap_or(BondFlags::NONE)
        });
    flags.contains(BondFlags::COVALENT) || flags.contains(BondFlags::METALLIC_COORDINATION)
}

const MOLSTAR_BACKBONE_SIZE_FACTOR: f32 = 0.3;
const MOLSTAR_STANDARD_BACKBONE_SHIFT: f64 = 0.5;
const MOLSTAR_NUCLEIC_BACKBONE_SHIFT: f64 = 0.3;
const MOLSTAR_STANDARD_TENSION: f64 = 0.5;
const MOLSTAR_DIRECTION_WEDGE_TENSION: f64 = 0.9;
const MOLSTAR_HELIX_TENSION: f64 = 0.9;
const MOLSTAR_TRACE_SIZE_FACTOR: f32 = 0.2;
const MOLSTAR_TRACE_SIZE_FACTOR64: f64 = 0.2;
const MOLSTAR_TRACE_ASPECT_RATIO: f32 = 5.0;
const MOLSTAR_TUBULAR_HELIX_FACTOR: f32 = 1.5;
const MOLSTAR_BALL_AND_STICK_SIZE_FACTOR64: f64 = 0.15;
const MOLSTAR_BACKBONE_SIZE_FACTOR64: f64 = 0.3;
const MOLSTAR_BOND_SIZE_ASPECT_RATIO64: f64 = 2.0 / 3.0;
const MOLSTAR_COMPONENT_BOND_RADIAL_SEGMENTS: usize = 36;

fn molstar_ball_and_stick_atom_radius(atom: &crate::model::Atom, options: &MeshOptions) -> f64 {
    vdw_radius64(&atom.type_symbol)
        * MOLSTAR_BALL_AND_STICK_SIZE_FACTOR64
        * molstar_radius_scale64(options)
}

fn molstar_spacefill_atom_radius(atom: &crate::model::Atom, options: &MeshOptions) -> f64 {
    vdw_radius64(&atom.type_symbol) * molstar_radius_scale64(options)
}

fn molstar_option_atom_radius64(options: &MeshOptions) -> f64 {
    let atom_radius = if options.atom_radius.to_bits() == 0.28f32.to_bits() {
        0.28
    } else {
        options.atom_radius as f64
    };
    atom_radius * molstar_radius_scale64(options)
}

fn molstar_ball_and_stick_bond_radius64(
    a: &crate::model::Atom,
    b: &crate::model::Atom,
    options: &MeshOptions,
) -> f64 {
    vdw_radius64(&a.type_symbol).min(vdw_radius64(&b.type_symbol))
        * MOLSTAR_BALL_AND_STICK_SIZE_FACTOR64
        * MOLSTAR_BOND_SIZE_ASPECT_RATIO64
        * molstar_radius_scale64(options)
}

fn polymer_trace_segment_count(
    linear_segments: usize,
    shift: f64,
    initial: bool,
    final_residue: bool,
) -> usize {
    let linear_segments = linear_segments.max(1);
    if initial {
        ((linear_segments as f64 * shift).round() as usize).max(1)
    } else if final_residue {
        ((linear_segments as f64 * (1.0 - shift)).round() as usize).max(1)
    } else {
        linear_segments
    }
}

fn molstar_trace_radius(options: &MeshOptions) -> f32 {
    options.ribbon_radius * options.radius_scale
}

fn molstar_cartoon_uniform_trace_radius(options: &MeshOptions) -> f32 {
    MOLSTAR_TRACE_SIZE_FACTOR * options.radius_scale
}

fn molstar_radius_scale64(options: &MeshOptions) -> f64 {
    if options.radius_scale.to_bits() == 1.0f32.to_bits() {
        1.0
    } else {
        options.radius_scale as f64
    }
}

fn molstar_trace_radius64(options: &MeshOptions) -> f64 {
    let radius = if options.ribbon_radius.to_bits() == MOLSTAR_TRACE_SIZE_FACTOR.to_bits() {
        MOLSTAR_TRACE_SIZE_FACTOR64
    } else {
        options.ribbon_radius as f64
    };
    radius * molstar_radius_scale64(options)
}

fn molstar_cartoon_uniform_trace_radius64(options: &MeshOptions) -> f64 {
    MOLSTAR_TRACE_SIZE_FACTOR64 * molstar_radius_scale64(options)
}

fn molstar_trace_height(options: &MeshOptions) -> f32 {
    molstar_trace_radius(options) * MOLSTAR_TRACE_ASPECT_RATIO
}

#[derive(Clone, Copy, Debug)]
struct PolymerBackboneLink {
    from_group: usize,
    to_group: usize,
    shift: f64,
}

fn add_polymer_backbone_semantic_objects(
    molecule: &Molecule,
    structure: &AtomicStructure,
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
    selected: &[String],
) {
    let trace = backbone_residues(molecule, structure);
    if trace.is_empty() {
        return;
    }

    let radius = MOLSTAR_BACKBONE_SIZE_FACTOR * options.radius_scale;
    if selected
        .iter()
        .any(|visual| visual == "polymer-backbone-cylinder")
    {
        for link in polymer_backbone_links(structure, &trace) {
            let from = &trace[link.from_group];
            let to = &trace[link.to_group];
            let middle = from.position + (to.position - from.position) * link.shift as f32;
            push_semantic_with_group(
                objects,
                link.from_group,
                SemanticMeta::new(
                    representation,
                    "backbone",
                    Some(&from.chain),
                    Some(from.seq),
                    Some(to.seq),
                )
                .with_visual("polymer-backbone-cylinder"),
                RenderObject::Cylinder {
                    start: from.position - center,
                    end: middle - center,
                    radius,
                },
            );
            push_semantic_with_group(
                objects,
                link.to_group,
                SemanticMeta::new(
                    representation,
                    "backbone",
                    Some(&to.chain),
                    Some(to.seq),
                    Some(from.seq),
                )
                .with_visual("polymer-backbone-cylinder"),
                RenderObject::Cylinder {
                    start: to.position - center,
                    end: middle - center,
                    radius,
                },
            );
        }
    }

    if selected
        .iter()
        .any(|visual| visual == "polymer-backbone-sphere")
    {
        for (group, residue) in trace.iter().enumerate() {
            push_semantic_with_group(
                objects,
                group,
                SemanticMeta::new(
                    representation,
                    "backbone",
                    Some(&residue.chain),
                    Some(residue.seq),
                    Some(residue.seq),
                )
                .with_visual("polymer-backbone-sphere"),
                RenderObject::Sphere {
                    center: residue.position - center,
                    radius: radius as f64,
                },
            );
        }
    }

    *group_id = (*group_id).max(trace.len());
}

fn polymer_backbone_links(
    structure: &AtomicStructure,
    trace: &[TraceResidue],
) -> Vec<PolymerBackboneLink> {
    let hierarchy = &structure.model.hierarchy;
    let mut links = Vec::new();

    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        let Some(start_residue) = hierarchy.residue_atom_segments.index.get(pair[0]).copied()
        else {
            continue;
        };
        let Some(end_residue) = hierarchy.residue_atom_segments.index.get(pair[1]).copied() else {
            continue;
        };

        let mut first_group = None;
        let mut previous_group = None;
        for residue_index in start_residue..=end_residue {
            let Some(group) =
                trace_residue_index_for_model_residue(hierarchy, trace, residue_index)
            else {
                continue;
            };
            first_group.get_or_insert(group);
            if let Some(from_group) = previous_group {
                links.push(PolymerBackboneLink {
                    from_group,
                    to_group: group,
                    shift: polymer_backbone_shift(trace[from_group].is_nucleotide),
                });
            }
            previous_group = Some(group);
        }

        if let (Some(from_group), Some(to_group)) = (previous_group, first_group) {
            if structure
                .ranges
                .cyclic_polymer_map
                .get(&end_residue)
                .copied()
                == Some(start_residue)
            {
                links.push(PolymerBackboneLink {
                    from_group,
                    to_group,
                    shift: polymer_backbone_shift(trace[from_group].is_nucleotide),
                });
            }
        }
    }

    links
}

fn polymer_backbone_shift(is_nucleotide: bool) -> f64 {
    if is_nucleotide {
        MOLSTAR_NUCLEIC_BACKBONE_SHIFT
    } else {
        MOLSTAR_STANDARD_BACKBONE_SHIFT
    }
}

fn realized_visuals(
    structure: &AtomicStructure,
    options: &MeshOptions,
    objects: &[SemanticRenderObject],
) -> Vec<String> {
    let selected = selected_visuals(structure, options);
    let mut realized = Vec::new();
    match options.representation {
        Representation::Spacefill => {
            if objects
                .iter()
                .any(|object| object.visual == "element-sphere")
            {
                push_visual_if_selected(&mut realized, &selected, "element-sphere");
            }
            push_visual_if_present(
                &mut realized,
                &selected,
                objects,
                "structure-element-sphere",
            );
        }
        Representation::BallAndStick => {
            if objects
                .iter()
                .any(|object| object.secondary_type == "atom" && object.geometry_type == "sphere")
            {
                push_visual_if_selected(&mut realized, &selected, "element-sphere");
            }
            if objects.iter().any(|object| object.secondary_type == "bond") {
                push_visual_if_selected(&mut realized, &selected, "intra-bond");
            }
            if !structure.inter_unit_bonds.is_empty() {
                push_visual_if_selected(&mut realized, &selected, "inter-bond");
            }
        }
        Representation::Molstar | Representation::Cartoon => {
            if objects
                .iter()
                .any(|object| matches!(object.secondary_type, "helix" | "sheet" | "coil"))
            {
                push_visual_if_selected(&mut realized, &selected, "polymer-trace");
            }
            if objects
                .iter()
                .any(|object| object.geometry_type == "nucleotide-ring")
            {
                push_visual_if_selected(&mut realized, &selected, "nucleotide-ring");
            }
            push_visual_if_present(&mut realized, &selected, objects, "nucleotide-block");
            push_visual_if_present(&mut realized, &selected, objects, "direction-wedge");
            push_visual_if_present(&mut realized, &selected, objects, "polymer-gap");
            push_visual_if_present(&mut realized, &selected, objects, "carbohydrate-symbol");
            push_visual_if_present(&mut realized, &selected, objects, "carbohydrate-link");
            push_visual_if_present(
                &mut realized,
                &selected,
                objects,
                "carbohydrate-terminal-link",
            );
            push_visual_if_present(&mut realized, &selected, objects, "element-sphere");
            push_visual_if_present(
                &mut realized,
                &selected,
                objects,
                "structure-element-sphere",
            );
            push_visual_if_present(&mut realized, &selected, objects, "intra-bond");
            push_visual_if_present(&mut realized, &selected, objects, "structure-intra-bond");
            push_visual_if_present(&mut realized, &selected, objects, "inter-bond");
        }
        Representation::Ribbon => {
            if objects
                .iter()
                .any(|object| matches!(object.secondary_type, "helix" | "sheet" | "coil"))
            {
                push_visual_if_selected(&mut realized, &selected, "polymer-tube");
            }
            push_visual_if_present(&mut realized, &selected, objects, "polymer-gap");
        }
        Representation::Backbone => {
            push_visual_if_present(
                &mut realized,
                &selected,
                objects,
                "polymer-backbone-cylinder",
            );
            push_visual_if_present(&mut realized, &selected, objects, "polymer-backbone-sphere");
            push_visual_if_present(&mut realized, &selected, objects, "polymer-gap");
        }
    }
    realized
}

fn push_visual_if_selected(realized: &mut Vec<String>, selected: &[String], visual: &str) {
    if selected.iter().any(|selected| selected == visual)
        && !realized.iter().any(|realized| realized == visual)
    {
        realized.push(visual.to_string());
    }
}

fn push_visual_if_present(
    realized: &mut Vec<String>,
    selected: &[String],
    objects: &[SemanticRenderObject],
    visual: &str,
) {
    if objects.iter().any(|object| object.visual == visual) {
        push_visual_if_selected(realized, selected, visual);
    }
}

fn geometry_type(object: &RenderObject) -> &'static str {
    match object {
        RenderObject::Sphere { .. } => "sphere",
        RenderObject::Cylinder { .. }
        | RenderObject::LinkCylinder { .. }
        | RenderObject::LinkCylinderWithSegments { .. } => "cylinder",
        RenderObject::Tube { .. } => "tube",
        RenderObject::DashedTube { .. } => "dashed-tube",
        RenderObject::FixedCountDashedCylinder { .. } => "dashed-cylinder",
        RenderObject::Ribbon { .. } => "ribbon",
        RenderObject::Sheet { .. } => "sheet",
        RenderObject::OrientedRibbon { .. } => "oriented-ribbon",
        RenderObject::PolymerTraceSegment { kind, .. } => match kind {
            PolymerTraceSegmentKind::Ribbon { .. } => "ribbon",
            PolymerTraceSegmentKind::Tube { .. } => "tube",
            PolymerTraceSegmentKind::Sheet { .. } => "sheet",
        },
        RenderObject::NucleotideRing { .. } => "nucleotide-ring",
        RenderObject::NucleotideBlock { .. } => "nucleotide-block",
        RenderObject::DirectionWedge { .. } => "direction-wedge",
        RenderObject::CarbohydrateSymbol { .. } => "carbohydrate-symbol",
        RenderObject::Ellipsoid { .. } => "ellipsoid",
    }
}

fn add_coarse_semantic_objects(
    molecule: &Molecule,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
) {
    *group_id = 0;
    for sphere in &molecule.coarse_spheres {
        push_semantic(
            objects,
            group_id,
            SemanticMeta::new(
                representation,
                "coarse-sphere",
                Some(&sphere.asym_id),
                Some(sphere.seq_id_begin),
                Some(sphere.seq_id_end),
            ),
            RenderObject::Sphere {
                center: sphere.position - center,
                radius: sphere.radius as f64,
            },
        );
    }
    *group_id = 0;
    for gaussian in &molecule.coarse_gaussians {
        push_semantic(
            objects,
            group_id,
            SemanticMeta::new(
                representation,
                "coarse-gaussian",
                Some(&gaussian.asym_id),
                Some(gaussian.seq_id_begin),
                Some(gaussian.seq_id_end),
            ),
            RenderObject::Ellipsoid {
                center: gaussian.position - center,
                axes: gaussian_axes(gaussian.covariance, gaussian.weight),
            },
        );
    }
}

fn gaussian_axes(covariance: [[f32; 3]; 3], weight: f32) -> [Vec3; 3] {
    let scale = weight.abs().sqrt().max(0.1);
    [
        Vec3::new(covariance[0][0].abs().sqrt().max(0.1) * scale, 0.0, 0.0),
        Vec3::new(0.0, covariance[1][1].abs().sqrt().max(0.1) * scale, 0.0),
        Vec3::new(0.0, 0.0, covariance[2][2].abs().sqrt().max(0.1) * scale),
    ]
}

fn add_nucleotide_semantic_objects(
    trace: &[TraceResidue],
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
    selected: &[String],
) {
    *group_id = 0;
    let include_ring = selected.iter().any(|visual| visual == "nucleotide-ring");
    let include_block = selected.iter().any(|visual| visual == "nucleotide-block");
    if !include_ring && !include_block {
        return;
    }
    let nucleotides = trace
        .iter()
        .filter(|residue| residue.is_nucleotide)
        .collect::<Vec<_>>();
    for residue in &nucleotides {
        let object_group_id = *group_id;
        if include_ring {
            let normal = residue.base_normal.unwrap_or_else(|| {
                fallback_side(
                    Vec3 {
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    },
                    None,
                )
            });
            push_semantic_with_group(
                objects,
                object_group_id,
                SemanticMeta::new(
                    representation,
                    "nucleotide",
                    Some(&residue.chain),
                    Some(residue.seq),
                    Some(residue.seq),
                )
                .with_trace_flags(trace_flags_from_residues(&[*residue]))
                .with_visual("nucleotide-ring"),
                RenderObject::NucleotideRing {
                    center: residue.base_center.unwrap_or(residue.position) - center,
                    normal,
                    radius: 0.2,
                    base: residue.nucleotide_atoms.and_then(|atoms| {
                        nucleotide_ring_base(atoms, residue.nucleotide_base_kind, center)
                    }),
                    detail: options.sphere_detail.min(3),
                    radial_segments: options.radial_segments,
                },
            );
        }
        if include_block {
            if let Some(geometry) = residue.nucleotide_atoms.and_then(|atoms| {
                nucleotide_block_geometry(atoms, residue.nucleotide_base_kind, center)
            }) {
                push_semantic_with_group(
                    objects,
                    object_group_id,
                    SemanticMeta::new(
                        representation,
                        "nucleotide",
                        Some(&residue.chain),
                        Some(residue.seq),
                        Some(residue.seq),
                    )
                    .with_trace_flags(trace_flags_from_residues(&[*residue]))
                    .with_visual("nucleotide-block"),
                    RenderObject::NucleotideBlock {
                        geometry,
                        radius: 0.2,
                        width: 4.5,
                        depth: 0.4,
                        radial_segments: options.radial_segments,
                    },
                );
            }
        }
        *group_id += 1;
    }
}

#[allow(clippy::too_many_arguments)]
fn add_direction_wedge_semantic_objects(
    trace: &[TraceResidue],
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
    selected: &[String],
    structure: &AtomicStructure,
) {
    if !selected.iter().any(|visual| visual == "direction-wedge") {
        return;
    }
    *group_id = 0;
    let mut had_polymer_range = false;
    let hierarchy = &structure.model.hierarchy;
    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        had_polymer_range = true;
        let Some(start_residue) = hierarchy.residue_atom_segments.index.get(pair[0]).copied()
        else {
            continue;
        };
        let Some(end_residue) = hierarchy.residue_atom_segments.index.get(pair[1]).copied() else {
            continue;
        };

        for residue_index in start_residue..=end_residue {
            let Some(trace_index) =
                trace_residue_index_for_model_residue(hierarchy, trace, residue_index)
            else {
                *group_id += 1;
                continue;
            };
            let residue = &trace[trace_index];
            let current_type = molstar_secondary_trace_type(structure, residue_index);
            let previous_residue = polymer_trace_residue_index(
                structure,
                start_residue,
                end_residue,
                residue_index as isize - 1,
            );
            let next_residue = polymer_trace_residue_index(
                structure,
                start_residue,
                end_residue,
                residue_index as isize + 1,
            );
            let previous_type = if previous_residue == residue_index {
                SecondaryStructureType::NONE
            } else {
                molstar_secondary_trace_type(structure, previous_residue)
            };
            let next_type = if next_residue == residue_index {
                SecondaryStructureType::NONE
            } else {
                molstar_secondary_trace_type(structure, next_residue)
            };
            let sec_struc_first = previous_type != current_type;
            let sec_struc_last = current_type != next_type;
            let is_sheet = current_type.contains(SecondaryStructureType::BETA);
            if is_sheet && sec_struc_last {
                *group_id += 1;
                continue;
            }

            let iterator_state = polymer_trace_iterator_state(
                structure,
                start_residue,
                end_residue,
                residue_index,
                current_type,
                false,
            );
            let center = DVec3::from_vec3(center);
            let controls = geometry::CurveSegmentControls {
                sec_struc_first,
                sec_struc_last,
                p0: iterator_state.p0 - center,
                p1: iterator_state.p1 - center,
                p2: iterator_state.p2 - center,
                p3: iterator_state.p3 - center,
                p4: iterator_state.p4 - center,
                d12: iterator_state.d12,
                d23: iterator_state.d23,
            };
            let tension = if residue.is_nucleotide || is_sheet {
                MOLSTAR_STANDARD_TENSION
            } else {
                MOLSTAR_DIRECTION_WEDGE_TENSION
            };
            let shift = if residue.is_nucleotide {
                MOLSTAR_NUCLEIC_BACKBONE_SHIFT
            } else {
                MOLSTAR_STANDARD_BACKBONE_SHIFT
            };
            let mut curve_state = geometry::CurveSegmentState::new(1);
            geometry::interpolate_curve_segment(&mut curve_state, &controls, tension, shift);
            let vectors = if residue.is_nucleotide {
                &curve_state.binormal_vectors
            } else {
                &curve_state.normal_vectors
            };
            let up = (vectors[0] + vectors[1]).normalized();
            let tangent = (controls.p3 - controls.p1).normalized().to_vec3();

            push_semantic(
                objects,
                group_id,
                SemanticMeta::new(
                    representation,
                    "direction",
                    Some(&residue.chain),
                    Some(residue.seq),
                    Some(residue.seq),
                )
                .with_trace_flags(TraceFlags {
                    initial: residue.initial,
                    final_residue: residue.final_residue,
                    sec_struc_first,
                    sec_struc_last,
                })
                .with_visual("direction-wedge"),
                RenderObject::DirectionWedge {
                    center: controls.p2.to_vec3(),
                    tangent,
                    up,
                    size: polymer_trace_radius(structure, residue_index, options),
                },
            );
        }
    }
    if !had_polymer_range {
        for residue in trace {
            let Some(tangent) = residue.direction else {
                *group_id += 1;
                continue;
            };
            push_semantic(
                objects,
                group_id,
                SemanticMeta::new(
                    representation,
                    "direction",
                    Some(&residue.chain),
                    Some(residue.seq),
                    Some(residue.seq),
                )
                .with_trace_flags(trace_flags_from_residues(&[residue]))
                .with_visual("direction-wedge"),
                RenderObject::DirectionWedge {
                    center: residue.position - center,
                    tangent,
                    up: residue
                        .base_normal
                        .unwrap_or_else(|| fallback_side(tangent, None)),
                    size: 0.2,
                },
            );
        }
    }
}

fn add_carbohydrate_symbol_semantic_objects(
    _molecule: &Molecule,
    structure: &AtomicStructure,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
    selected: &[String],
) {
    if !selected
        .iter()
        .any(|visual| visual == "carbohydrate-symbol")
    {
        return;
    }
    *group_id = 0;

    for (carbohydrate_index, carb) in structure.carbohydrates.elements.iter().enumerate() {
        let shape = get_saccharide_shape(
            carb.component.component_type,
            carb.ring_element_indices.len(),
        );
        let (chain, seq) = carbohydrate_residue_metadata(&structure, carb.residue_index);
        let base_group = carbohydrate_index * 2;
        let meta = SemanticMeta::new(representation, "carbohydrate", chain.as_deref(), seq, seq)
            .with_visual("carbohydrate-symbol");

        let part = if carbohydrate_symbol_has_secondary_part(shape) {
            CarbohydrateSymbolPart::Primary
        } else {
            CarbohydrateSymbolPart::Whole
        };

        push_semantic_with_group(
            objects,
            base_group,
            meta,
            RenderObject::CarbohydrateSymbol {
                center: carb.geometry.center - center,
                normal: carb.geometry.normal,
                direction: carb.geometry.direction,
                shape,
                part,
            },
        );

        if carbohydrate_symbol_has_secondary_part(shape) {
            let meta =
                SemanticMeta::new(representation, "carbohydrate", chain.as_deref(), seq, seq)
                    .with_visual("carbohydrate-symbol");
            push_semantic_with_group(
                objects,
                base_group + 1,
                meta,
                RenderObject::CarbohydrateSymbol {
                    center: carb.geometry.center - center,
                    normal: carb.geometry.normal,
                    direction: carb.geometry.direction,
                    shape,
                    part: CarbohydrateSymbolPart::Secondary,
                },
            );
        }

        *group_id = base_group + 1 + usize::from(carbohydrate_symbol_has_secondary_part(shape));
    }
}

fn add_carbohydrate_link_semantic_objects(
    molecule: &Molecule,
    structure: &AtomicStructure,
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
    selected: &[String],
) {
    if !selected.iter().any(|visual| visual == "carbohydrate-link") {
        return;
    }
    let carbohydrates = &structure.carbohydrates;
    *group_id = 0;

    for (link_index, link) in carbohydrates.links.iter().enumerate() {
        let Some(carb_a) = carbohydrates.elements.get(link.carbohydrate_index_a) else {
            continue;
        };
        let Some(carb_b) = carbohydrates.elements.get(link.carbohydrate_index_b) else {
            continue;
        };
        let (chain, start_seq) = carbohydrate_residue_metadata(&structure, carb_a.residue_index);
        let (_, end_seq) = carbohydrate_residue_metadata(&structure, carb_b.residue_index);
        push_semantic_with_group(
            objects,
            link_index,
            SemanticMeta::new(
                representation,
                "carbohydrate",
                chain.as_deref(),
                start_seq,
                end_seq,
            )
            .with_visual("carbohydrate-link"),
            RenderObject::LinkCylinder {
                start: carb_a.geometry.center - center,
                end: carb_b.geometry.center - center,
                radius: carbohydrate_link_radius(molecule, &structure, carb_a, options),
            },
        );
        *group_id = link_index + 1;
    }
}

fn add_carbohydrate_terminal_link_semantic_objects(
    molecule: &Molecule,
    structure: &AtomicStructure,
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
    selected: &[String],
) {
    if !selected
        .iter()
        .any(|visual| visual == "carbohydrate-terminal-link")
    {
        return;
    }
    let carbohydrates = &structure.carbohydrates;
    *group_id = 0;

    for (link_index, link) in carbohydrates.terminal_links.iter().enumerate() {
        let Some(carb) = carbohydrates.elements.get(link.carbohydrate_index) else {
            continue;
        };
        let Some(element_position) =
            carbohydrate_terminal_element_position(molecule, &structure, link)
        else {
            continue;
        };

        let (carb_chain, carb_seq) = carbohydrate_residue_metadata(&structure, carb.residue_index);
        let (element_chain, element_seq) = carbohydrate_terminal_element_metadata(&structure, link);
        let (start, end, chain, start_seq, end_seq) = if link.from_carbohydrate {
            (
                carb.geometry.center,
                element_position,
                carb_chain,
                carb_seq,
                element_seq,
            )
        } else {
            (
                element_position,
                carb.geometry.center,
                element_chain,
                element_seq,
                carb_seq,
            )
        };

        push_semantic_with_group(
            objects,
            link_index,
            SemanticMeta::new(
                representation,
                "carbohydrate",
                chain.as_deref(),
                start_seq,
                end_seq,
            )
            .with_visual("carbohydrate-terminal-link"),
            RenderObject::LinkCylinder {
                start: start - center,
                end: end - center,
                radius: carbohydrate_terminal_link_radius(
                    molecule, &structure, carb, link, options,
                ),
            },
        );
        *group_id = link_index + 1;
    }
}

fn carbohydrate_residue_metadata(
    structure: &AtomicStructure,
    residue_index: usize,
) -> (Option<String>, Option<i32>) {
    let Some(residue) = structure.model.hierarchy.residues.get(residue_index) else {
        return (None, None);
    };
    let chain = structure
        .model
        .hierarchy
        .chains
        .get(residue.chain_index)
        .map(|chain| chain.id.clone());
    let seq = residue.label_seq_id.trim().parse::<i32>().ok();
    (chain, seq)
}

fn carbohydrate_terminal_element_metadata(
    structure: &AtomicStructure,
    link: &crate::model::CarbohydrateTerminalLink,
) -> (Option<String>, Option<i32>) {
    let Some(unit) = structure.units.get(link.element_unit_id) else {
        return (None, None);
    };
    let Some(&element) = unit.elements.get(link.element_index) else {
        return (None, None);
    };
    let Some(&residue_index) = unit.residue_index_by_element.get(element) else {
        return (None, None);
    };
    carbohydrate_residue_metadata(structure, residue_index)
}

fn carbohydrate_terminal_element_position(
    molecule: &Molecule,
    structure: &AtomicStructure,
    link: &crate::model::CarbohydrateTerminalLink,
) -> Option<Vec3> {
    structure
        .position(link.element_unit_id, link.element_index)
        .or_else(|| {
            structure
                .units
                .get(link.element_unit_id)
                .and_then(|unit| unit.atom_indices.get(link.element_index))
                .and_then(|&source_atom| molecule.atoms.get(source_atom))
                .map(|atom| atom.position)
        })
}

fn carbohydrate_link_radius(
    molecule: &Molecule,
    structure: &AtomicStructure,
    carb: &crate::model::CarbohydrateElement,
    options: &MeshOptions,
) -> f32 {
    let radius = structure
        .units
        .get(carb.unit_id)
        .and_then(|unit| {
            carb.ring_element_indices
                .first()
                .and_then(|&i| unit.atom_indices.get(i))
        })
        .and_then(|&source_atom| molecule.atoms.get(source_atom))
        .map(|atom| vdw_radius(&atom.type_symbol))
        .unwrap_or(1.0);
    radius * 0.3 * options.radius_scale
}

fn carbohydrate_terminal_link_radius(
    molecule: &Molecule,
    structure: &AtomicStructure,
    carb: &crate::model::CarbohydrateElement,
    link: &crate::model::CarbohydrateTerminalLink,
    options: &MeshOptions,
) -> f32 {
    let radius = if link.from_carbohydrate {
        structure
            .units
            .get(carb.unit_id)
            .and_then(|unit| {
                carb.ring_element_indices
                    .first()
                    .and_then(|&i| unit.atom_indices.get(i))
            })
            .and_then(|&source_atom| molecule.atoms.get(source_atom))
    } else {
        structure
            .units
            .get(link.element_unit_id)
            .and_then(|unit| unit.atom_indices.get(link.element_index))
            .and_then(|&source_atom| molecule.atoms.get(source_atom))
    }
    .map(|atom| vdw_radius(&atom.type_symbol))
    .unwrap_or(1.0);
    radius * 0.2 * options.radius_scale
}

fn nucleotide_ring_base(
    atoms: NucleotideAtoms,
    kind: Option<NucleotideBaseKind>,
    center_offset: Vec3,
) -> Option<NucleotideRingBase> {
    let kind = kind.unwrap_or_else(|| match (atoms.c4, atoms.n9) {
        (Some(c4), Some(n9)) if c4.distance(n9) < 1.6 => NucleotideBaseKind::Purine,
        _ => NucleotideBaseKind::Pyrimidine,
    });
    let translate = |position: Vec3| position - center_offset;
    match kind {
        NucleotideBaseKind::Purine => {
            let trace = translate(atoms.trace?);
            let n9 = translate(atoms.n9?);
            match (
                atoms.n1,
                atoms.c2,
                atoms.n3,
                atoms.c4,
                atoms.c5.or(atoms.n5),
                atoms.c6,
                atoms.n7.or(atoms.c7),
                atoms.c8,
            ) {
                (
                    Some(n1),
                    Some(c2),
                    Some(n3),
                    Some(c4),
                    Some(c5),
                    Some(c6),
                    Some(n7),
                    Some(c8),
                ) => Some(NucleotideRingBase::Purine {
                    trace,
                    n1: translate(n1),
                    c2: translate(c2),
                    n3: translate(n3),
                    c4: translate(c4),
                    c5: translate(c5),
                    c6: translate(c6),
                    n7: translate(n7),
                    c8: translate(c8),
                    n9,
                }),
                _ => Some(NucleotideRingBase::PurineConnector { trace, n9 }),
            }
        }
        NucleotideBaseKind::Pyrimidine => {
            let trace = translate(atoms.trace?);
            let n1 = translate(atoms.n1.or(atoms.c1)?);
            match (atoms.c2, atoms.n3, atoms.c4, atoms.c5, atoms.c6) {
                (Some(c2), Some(n3), Some(c4), Some(c5), Some(c6)) => {
                    Some(NucleotideRingBase::Pyrimidine {
                        trace,
                        n1,
                        c2: translate(c2),
                        n3: translate(n3),
                        c4: translate(c4),
                        c5: translate(c5),
                        c6: translate(c6),
                    })
                }
                _ => Some(NucleotideRingBase::PyrimidineConnector { trace, n1 }),
            }
        }
    }
}

fn nucleotide_block_geometry(
    atoms: NucleotideAtoms,
    kind: Option<NucleotideBaseKind>,
    center_offset: Vec3,
) -> Option<NucleotideBlockGeometry> {
    let kind = kind.unwrap_or_else(|| match (atoms.c4, atoms.n9) {
        (Some(c4), Some(n9)) if c4.distance(n9) < 1.6 => NucleotideBaseKind::Purine,
        _ => NucleotideBaseKind::Pyrimidine,
    });
    let translate = |position: Vec3| position - center_offset;
    let trace = translate(atoms.trace?);
    let (anchor, block) = match kind {
        NucleotideBaseKind::Purine => (
            translate(atoms.n9?),
            match (atoms.n1, atoms.c4, atoms.c6, atoms.c2) {
                (Some(p1), Some(p2), Some(p3), Some(p4)) => Some(NucleotideBlockBox {
                    p1: translate(p1),
                    p2: translate(p2),
                    p3: translate(p3),
                    p4: translate(p4),
                    height: 4.5,
                }),
                _ => None,
            },
        ),
        NucleotideBaseKind::Pyrimidine => (
            translate(atoms.n1.or(atoms.c1)?),
            match (atoms.n3, atoms.c6, atoms.c4, atoms.c2) {
                (Some(p1), Some(p2), Some(p3), Some(p4)) => Some(NucleotideBlockBox {
                    p1: translate(p1),
                    p2: translate(p2),
                    p3: translate(p3),
                    p4: translate(p4),
                    height: 3.0,
                }),
                _ => None,
            },
        ),
    };
    Some(NucleotideBlockGeometry {
        trace,
        anchor,
        block,
    })
}

fn trace_flags_from_residues(residues: &[&TraceResidue]) -> TraceFlags {
    let Some(first) = residues.first() else {
        return TraceFlags::default();
    };
    let last = residues.last().copied().unwrap_or(*first);
    TraceFlags {
        initial: first.initial,
        final_residue: last.final_residue,
        sec_struc_first: false,
        sec_struc_last: false,
    }
}

fn secondary_trace_flags(
    _trace: &[TraceResidue],
    residues: &[&TraceResidue],
    _molecule: &Molecule,
    _kind: SecondaryTraceKind,
) -> TraceFlags {
    let mut flags = trace_flags_from_residues(residues);
    let Some(first) = residues.first().copied() else {
        return flags;
    };
    let last = residues.last().copied().unwrap_or(first);

    flags.sec_struc_first = first.sec_struc_first;
    flags.sec_struc_last = last.sec_struc_last;
    flags
}

fn secondary_trace_cap_flags(
    structure: &AtomicStructure,
    residues: &[&TraceResidue],
    flags: TraceFlags,
) -> (bool, bool) {
    let Some(first) = residues.first().copied() else {
        return (false, false);
    };
    let last = residues.last().copied().unwrap_or(first);
    (
        flags.sec_struc_first || trace_residue_is_polymer_range_boundary(structure, first, true),
        flags.sec_struc_last || trace_residue_is_polymer_range_boundary(structure, last, false),
    )
}

fn trace_residue_is_polymer_range_boundary(
    structure: &AtomicStructure,
    trace_residue: &TraceResidue,
    start: bool,
) -> bool {
    let hierarchy = &structure.model.hierarchy;
    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        let element_index = if start { pair[0] } else { pair[1] };
        let Some(residue_index) = hierarchy
            .residue_atom_segments
            .index
            .get(element_index)
            .copied()
        else {
            continue;
        };
        if trace_residue_matches_model_residue(hierarchy, trace_residue, residue_index) {
            return true;
        }
    }
    false
}

fn trace_residue_matches_model_residue(
    hierarchy: &crate::model::AtomicHierarchy,
    trace_residue: &TraceResidue,
    residue_index: usize,
) -> bool {
    let Some(residue) = hierarchy.residues.get(residue_index) else {
        return false;
    };
    let Some(chain) = hierarchy.chains.get(residue.chain_index) else {
        return false;
    };
    let Some(seq) = residue.label_seq_id.trim().parse::<i32>().ok() else {
        return false;
    };
    trace_residue.chain == chain.id
        && trace_residue.seq == seq
        && trace_residue.insertion_code == residue.insertion_code
}

fn trace_flags_for_segment(
    trace: &[TraceResidue],
    chain: &str,
    start: i32,
    start_insertion_code: &str,
    end: i32,
    end_insertion_code: &str,
    secondary: bool,
) -> TraceFlags {
    let residues = trace
        .iter()
        .filter(|residue| {
            residue.chain == chain
                && residue_position_cmp(
                    residue.seq,
                    &residue.insertion_code,
                    start,
                    start_insertion_code,
                )
                .is_ge()
                && residue_position_cmp(
                    residue.seq,
                    &residue.insertion_code,
                    end,
                    end_insertion_code,
                )
                .is_le()
        })
        .collect::<Vec<_>>();
    let mut flags = trace_flags_from_residues(&residues);
    if secondary {
        flags.sec_struc_first = true;
        flags.sec_struc_last = true;
    }
    flags
}

fn apply_polymer_trace_terminal_flags(structure: &AtomicStructure, trace: &mut [TraceResidue]) {
    if trace.is_empty() {
        return;
    }

    let hierarchy = &structure.model.hierarchy;
    if structure.ranges.polymer_ranges.is_empty() {
        return;
    }

    for residue in trace.iter_mut() {
        residue.initial = false;
        residue.final_residue = false;
    }

    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        let Some(start_residue) = hierarchy.residue_atom_segments.index.get(pair[0]).copied()
        else {
            continue;
        };
        let Some(end_residue) = hierarchy.residue_atom_segments.index.get(pair[1]).copied() else {
            continue;
        };

        if !structure
            .ranges
            .cyclic_polymer_map
            .contains_key(&start_residue)
        {
            set_trace_terminal_flag(hierarchy, trace, start_residue, true);
        }
        if !structure
            .ranges
            .cyclic_polymer_map
            .contains_key(&end_residue)
        {
            set_trace_terminal_flag(hierarchy, trace, end_residue, false);
        }
    }
}

fn set_trace_terminal_flag(
    hierarchy: &crate::model::AtomicHierarchy,
    trace: &mut [TraceResidue],
    residue_index: usize,
    initial: bool,
) {
    let Some(residue) = hierarchy.residues.get(residue_index) else {
        return;
    };
    let Some(chain) = hierarchy.chains.get(residue.chain_index) else {
        return;
    };
    let Some(seq) = residue.label_seq_id.trim().parse::<i32>().ok() else {
        return;
    };
    let Some(trace_residue) = trace.iter_mut().find(|trace_residue| {
        trace_residue.chain == chain.id
            && trace_residue.seq == seq
            && trace_residue.insertion_code == residue.insertion_code
    }) else {
        return;
    };
    if initial {
        trace_residue.initial = true;
    } else {
        trace_residue.final_residue = true;
    }
}

fn apply_cyclic_polymer_trace_flags(structure: &AtomicStructure, trace: &mut [TraceResidue]) {
    if trace.is_empty() {
        return;
    }
    if structure.ranges.cyclic_polymer_map.is_empty() {
        return;
    }

    for residue_index in structure.ranges.cyclic_polymer_map.keys().copied() {
        let Some(residue) = structure.model.hierarchy.residues.get(residue_index) else {
            continue;
        };
        let Some(chain) = structure
            .model
            .hierarchy
            .chains
            .get(residue.chain_index)
            .map(|chain| chain.id.as_str())
        else {
            continue;
        };
        let Some(seq) = residue.label_seq_id.trim().parse::<i32>().ok() else {
            continue;
        };
        if let Some(trace_residue) = trace.iter_mut().find(|trace_residue| {
            trace_residue.chain == chain
                && trace_residue.seq == seq
                && trace_residue.insertion_code == residue.insertion_code
        }) {
            trace_residue.initial = false;
            trace_residue.final_residue = false;
        }
    }
}

fn apply_polymer_trace_secondary_flags(structure: &AtomicStructure, trace: &mut [TraceResidue]) {
    if trace.is_empty() {
        return;
    }
    for residue in trace.iter_mut() {
        residue.sec_struc_first = false;
        residue.sec_struc_last = false;
    }

    if structure.ranges.polymer_ranges.is_empty() {
        return;
    }
    let hierarchy = &structure.model.hierarchy;
    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        let Some(start_residue) = hierarchy.residue_atom_segments.index.get(pair[0]).copied()
        else {
            continue;
        };
        let Some(end_residue) = hierarchy.residue_atom_segments.index.get(pair[1]).copied() else {
            continue;
        };
        for residue_index in start_residue..=end_residue {
            let Some(trace_index) =
                trace_residue_index_for_model_residue(hierarchy, trace, residue_index)
            else {
                continue;
            };
            let previous_residue = polymer_trace_residue_index(
                &structure,
                start_residue,
                end_residue,
                residue_index as isize - 1,
            );
            let next_residue = polymer_trace_residue_index(
                &structure,
                start_residue,
                end_residue,
                residue_index as isize + 1,
            );
            let current_type = molstar_secondary_trace_type(&structure, residue_index);
            let previous_type = molstar_secondary_trace_type(&structure, previous_residue);
            let next_type = molstar_secondary_trace_type(&structure, next_residue);
            trace[trace_index].sec_struc_first = previous_type != current_type;
            trace[trace_index].sec_struc_last = current_type != next_type;
        }
    }
}

fn polymer_trace_residue_index(
    structure: &AtomicStructure,
    segment_min: usize,
    segment_max: usize,
    residue_index: isize,
) -> usize {
    if residue_index < segment_min as isize {
        if let Some(&cyclic_index) = structure.ranges.cyclic_polymer_map.get(&segment_min) {
            return (cyclic_index as isize - (segment_min as isize - residue_index - 1)).max(0)
                as usize;
        }
        segment_min
    } else if residue_index > segment_max as isize {
        if let Some(&cyclic_index) = structure.ranges.cyclic_polymer_map.get(&segment_max) {
            return (cyclic_index as isize + (residue_index - segment_max as isize - 1)).max(0)
                as usize;
        }
        segment_max
    } else {
        residue_index as usize
    }
}

fn molstar_helix_orientation_centers(structure: &AtomicStructure) -> Vec<Vec3> {
    let hierarchy = &structure.model.hierarchy;
    let residue_count = hierarchy.derived.residue.polymer_type.len();
    let mut centers = vec![Vec3::new(f32::NAN, f32::NAN, f32::NAN); residue_count];

    for pair in structure.ranges.polymer_ranges.chunks_exact(2) {
        let Some(start_residue) = hierarchy.residue_atom_segments.index.get(pair[0]).copied()
        else {
            continue;
        };
        let Some(end_residue) = hierarchy.residue_atom_segments.index.get(pair[1]).copied() else {
            continue;
        };
        if end_residue.saturating_sub(start_residue) + 1 < 4 {
            continue;
        }

        let mut trace_window = [DVec3::default(); 4];

        for (i, residue_index) in (start_residue..=end_residue).enumerate() {
            trace_window[0] = trace_window[1];
            trace_window[1] = trace_window[2];
            trace_window[2] = trace_window[3];
            let Some(trace_position) = polymer_trace_atom_position(structure, residue_index) else {
                continue;
            };
            trace_window[3] = DVec3::from_vec3(trace_position);

            if i < 3 {
                continue;
            }

            let r12 = trace_window[1] - trace_window[0];
            let r23 = trace_window[2] - trace_window[1];
            let r34 = trace_window[3] - trace_window[2];

            let diff13 = r12 - r23;
            let diff24 = r23 - r34;
            let diff13_len = diff13.squared_length().sqrt();
            let diff24_len = diff24.squared_length().sqrt();
            if diff13_len == 0.0 || diff24_len == 0.0 {
                continue;
            }

            let tmp = molstar_vec3_angle(diff13, diff24).cos();
            let radius = (diff24_len * diff13_len).sqrt() / (2.0 * (1.0 - tmp)).max(2.0);
            let first_center = trace_window[1] - diff13 * (radius / diff13_len);
            let second_center = trace_window[2] - diff24 * (radius / diff24_len);
            let center_index = residue_index - 2;
            if let Some(center) = centers.get_mut(center_index) {
                *center = first_center.to_vec3();
            }
            if let Some(center) = centers.get_mut(center_index + 1) {
                *center = second_center.to_vec3();
            }
        }

        if let (Some(first_axis_a), Some(first_axis_b), Some(first_trace)) = (
            centers.get(start_residue + 1).copied(),
            centers.get(start_residue + 2).copied(),
            polymer_trace_atom_position(structure, start_residue),
        ) {
            let first_axis_a = DVec3::from_vec3(first_axis_a);
            let axis = (first_axis_a - DVec3::from_vec3(first_axis_b)).normalized();
            if axis.squared_length() > 0.0 {
                centers[start_residue] = molstar_project_point_on_vector_d(
                    DVec3::from_vec3(first_trace),
                    axis,
                    first_axis_a,
                )
                .to_vec3();
            }
        }

        if end_residue >= 2 {
            if let (Some(last_axis_a), Some(last_axis_b), Some(last_trace)) = (
                centers.get(end_residue - 1).copied(),
                centers.get(end_residue - 2).copied(),
                polymer_trace_atom_position(structure, end_residue),
            ) {
                let last_axis_a = DVec3::from_vec3(last_axis_a);
                let axis = (last_axis_a - DVec3::from_vec3(last_axis_b)).normalized();
                if axis.squared_length() > 0.0 {
                    centers[end_residue] = molstar_project_point_on_vector_d(
                        DVec3::from_vec3(last_trace),
                        axis,
                        last_axis_a,
                    )
                    .to_vec3();
                }
            }
        }
    }

    centers
}

fn polymer_trace_atom_position(structure: &AtomicStructure, residue_index: usize) -> Option<Vec3> {
    structure
        .model
        .hierarchy
        .derived
        .residue
        .trace_element_index
        .get(residue_index)
        .and_then(|index| *index)
        .and_then(|atom_index| structure.model.hierarchy.atoms.get(atom_index))
        .map(|atom| atom.position)
}

fn molstar_vec3_angle(a: DVec3, b: DVec3) -> f64 {
    let denominator = (a.squared_length() * b.squared_length()).sqrt();
    if denominator == 0.0 {
        return std::f64::consts::PI / 2.0;
    }
    (a.dot(b) / denominator).clamp(-1.0, 1.0).acos()
}

fn molstar_project_point_on_vector_d(point: DVec3, vector: DVec3, origin: DVec3) -> DVec3 {
    let out = point - origin;
    origin + vector * (vector.dot(out) / vector.dot(vector))
}

fn molstar_secondary_trace_type(
    structure: &AtomicStructure,
    residue_index: usize,
) -> SecondaryStructureType {
    let secondary_type = structure
        .model
        .secondary_structure
        .residue_type(residue_index);
    if secondary_type.contains(SecondaryStructureType::HELIX) {
        SecondaryStructureType::HELIX
    } else {
        secondary_type
    }
}

fn is_helix_secondary(value: SecondaryStructureType) -> bool {
    value.contains(SecondaryStructureType::HELIX)
}

fn trace_residue_index_for_model_residue(
    hierarchy: &crate::model::AtomicHierarchy,
    trace: &[TraceResidue],
    residue_index: usize,
) -> Option<usize> {
    let residue = hierarchy.residues.get(residue_index)?;
    let chain = hierarchy.chains.get(residue.chain_index)?;
    let seq = residue.label_seq_id.trim().parse::<i32>().ok()?;
    trace.iter().position(|trace_residue| {
        trace_residue.chain == chain.id
            && trace_residue.seq == seq
            && trace_residue.insertion_code == residue.insertion_code
    })
}

fn model_residue_index_for_trace_residue(
    hierarchy: &crate::model::AtomicHierarchy,
    trace_residue: &TraceResidue,
) -> Option<usize> {
    hierarchy
        .residues
        .iter()
        .enumerate()
        .find_map(|(residue_index, residue)| {
            let chain = hierarchy.chains.get(residue.chain_index)?;
            let seq = residue.label_seq_id.trim().parse::<i32>().ok()?;
            (chain.id == trace_residue.chain
                && seq == trace_residue.seq
                && residue.insertion_code == trace_residue.insertion_code)
                .then_some(residue_index)
        })
}

fn residue_in_secondary_range(residue: &TraceResidue, range: &SecondaryRange) -> bool {
    residue.chain == range.chain
        && residue_position_cmp(
            residue.seq,
            &residue.insertion_code,
            range.start,
            &range.start_insertion_code,
        )
        .is_ge()
        && residue_position_cmp(
            residue.seq,
            &residue.insertion_code,
            range.end,
            &range.end_insertion_code,
        )
        .is_le()
}

fn residue_position_cmp(
    seq: i32,
    insertion_code: &str,
    other_seq: i32,
    other_insertion_code: &str,
) -> std::cmp::Ordering {
    seq.cmp(&other_seq)
        .then_with(|| insertion_code.cmp(other_insertion_code))
}

fn add_ball_and_stick_semantic_objects(
    molecule: &Molecule,
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    group_id: &mut usize,
    objects: &mut Vec<SemanticRenderObject>,
) {
    *group_id = 0;
    let chain_materials = molstar_chain_materials(molecule);
    for atom in &molecule.atoms {
        push_semantic(
            objects,
            group_id,
            SemanticMeta::new(
                representation,
                "atom",
                Some(&atom.chain),
                atom.residue_seq.parse::<i32>().ok(),
                atom.residue_seq.parse::<i32>().ok(),
            )
            .with_visual("element-sphere")
            .with_material(molstar_atom_material(atom, &chain_materials, "atom")),
            RenderObject::Sphere {
                center: atom.position - center,
                radius: molstar_option_atom_radius64(options),
            },
        );
    }

    *group_id = 0;
    for bond in &molecule.bonds {
        let atom_a = &molecule.atoms[bond.a];
        push_semantic(
            objects,
            group_id,
            SemanticMeta::new(
                representation,
                "bond",
                Some(&molecule.atoms[bond.a].chain),
                molecule.atoms[bond.a].residue_seq.parse::<i32>().ok(),
                molecule.atoms[bond.b].residue_seq.parse::<i32>().ok(),
            )
            .with_visual("intra-bond")
            .with_material(molstar_atom_material(atom_a, &chain_materials, "bond")),
            RenderObject::Cylinder {
                start: molecule.atoms[bond.a].position - center,
                end: molecule.atoms[bond.b].position - center,
                radius: options.bond_radius,
            },
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn add_molstar_component_ball_and_stick_semantic_objects(
    molecule: &Molecule,
    options: &MeshOptions,
    center: Vec3,
    representation: &'static str,
    component: &'static str,
    atom_mask: &[bool],
    selected: &[String],
    objects: &mut Vec<SemanticRenderObject>,
) {
    let sphere_visual = selected
        .iter()
        .any(|visual| visual == "element-sphere")
        .then_some("element-sphere")
        .or_else(|| {
            selected
                .iter()
                .any(|visual| visual == "structure-element-sphere")
                .then_some("structure-element-sphere")
        });
    let bond_visual = selected
        .iter()
        .any(|visual| visual == "intra-bond")
        .then_some("intra-bond")
        .or_else(|| {
            selected
                .iter()
                .any(|visual| visual == "structure-intra-bond")
                .then_some("structure-intra-bond")
        });

    if let Some(visual) = sphere_visual {
        let mut group_id = 0usize;
        let chain_materials = molstar_chain_materials(molecule);
        for (atom_index, atom) in molecule.atoms.iter().enumerate() {
            if !atom_mask.get(atom_index).copied().unwrap_or(false) {
                continue;
            }
            push_semantic(
                objects,
                &mut group_id,
                SemanticMeta::new(
                    representation,
                    component,
                    Some(&atom.chain),
                    atom.residue_seq.parse::<i32>().ok(),
                    atom.residue_seq.parse::<i32>().ok(),
                )
                .with_visual(visual)
                .with_material(molstar_atom_material(
                    atom,
                    &chain_materials,
                    component,
                )),
                RenderObject::Sphere {
                    center: atom.position - center,
                    radius: molstar_ball_and_stick_atom_radius(atom, options),
                },
            );
        }
    }

    if let Some(visual) = bond_visual {
        let mut group_id = 0usize;
        let chain_materials = molstar_chain_materials(molecule);
        let mut unit_slot_half_links = BTreeSet::<(usize, usize)>::new();
        for (atom_a, atom_b) in molstar_component_unit_slot_half_links(molecule, atom_mask) {
            let Some(a) = molecule.atoms.get(atom_a) else {
                continue;
            };
            let Some(b) = molecule.atoms.get(atom_b) else {
                continue;
            };
            unit_slot_half_links.insert((atom_a, atom_b));
            push_molstar_component_bond_semantic_object(
                objects,
                &mut group_id,
                representation,
                component,
                visual,
                a,
                b,
                center,
                options,
                &chain_materials,
            );
        }
        for bond in &molecule.bonds {
            if !atom_mask.get(bond.a).copied().unwrap_or(false)
                || !atom_mask.get(bond.b).copied().unwrap_or(false)
            {
                continue;
            }
            let Some(a) = molecule.atoms.get(bond.a) else {
                continue;
            };
            let Some(b) = molecule.atoms.get(bond.b) else {
                continue;
            };
            if !unit_slot_half_links.contains(&(bond.a, bond.b)) {
                push_molstar_component_bond_semantic_object(
                    objects,
                    &mut group_id,
                    representation,
                    component,
                    visual,
                    a,
                    b,
                    center,
                    options,
                    &chain_materials,
                );
            }
            if !unit_slot_half_links.contains(&(bond.b, bond.a)) {
                push_molstar_component_bond_semantic_object(
                    objects,
                    &mut group_id,
                    representation,
                    component,
                    visual,
                    b,
                    a,
                    center,
                    options,
                    &chain_materials,
                );
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn push_molstar_component_bond_semantic_object(
    objects: &mut Vec<SemanticRenderObject>,
    group_id: &mut usize,
    representation: &'static str,
    component: &'static str,
    visual: &'static str,
    a: &crate::model::Atom,
    b: &crate::model::Atom,
    center: Vec3,
    options: &MeshOptions,
    chain_materials: &BTreeMap<String, MeshMaterial>,
) {
    let a_position = a.position - center;
    let b_position = b.position - center;
    let radius = molstar_ball_and_stick_bond_radius64(a, b, options);
    push_semantic(
        objects,
        group_id,
        SemanticMeta::new(
            representation,
            component,
            Some(&a.chain),
            a.residue_seq.parse::<i32>().ok(),
            b.residue_seq.parse::<i32>().ok(),
        )
        .with_visual(visual)
        .with_material(molstar_atom_material(a, chain_materials, component)),
        RenderObject::LinkCylinderWithSegments {
            start: a_position,
            end: b_position,
            radius,
            radial_segments: MOLSTAR_COMPONENT_BOND_RADIAL_SEGMENTS,
        },
    );
}

#[derive(Clone, Debug, Default)]
struct ComponentBondGroup {
    atoms: Vec<usize>,
    edges: Vec<(usize, usize)>,
}

fn molstar_component_unit_slot_half_links(
    molecule: &Molecule,
    atom_mask: &[bool],
) -> Vec<(usize, usize)> {
    let mut groups = Vec::<ComponentBondGroup>::new();
    let mut group_by_key = BTreeMap::<(String, String), usize>::new();
    let mut atom_local = vec![None; molecule.atoms.len()];

    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        if !atom_mask.get(atom_index).copied().unwrap_or(false) {
            continue;
        }
        let key = (atom.chain.clone(), atom.operator_name.clone());
        let group_index = if let Some(group_index) = group_by_key.get(&key).copied() {
            group_index
        } else {
            let group_index = groups.len();
            groups.push(ComponentBondGroup::default());
            group_by_key.insert(key, group_index);
            group_index
        };
        let local_index = groups[group_index].atoms.len();
        groups[group_index].atoms.push(atom_index);
        atom_local[atom_index] = Some((group_index, local_index));
    }

    for bond in &molecule.bonds {
        let Some((group_a, local_a)) = atom_local.get(bond.a).copied().flatten() else {
            continue;
        };
        let Some((group_b, local_b)) = atom_local.get(bond.b).copied().flatten() else {
            continue;
        };
        if group_a != group_b {
            continue;
        }
        groups[group_a].edges.push((local_a, local_b));
    }

    groups
        .iter()
        .flat_map(|group| molstar_edge_builder_directed_slots(&group.atoms, &group.edges))
        .collect()
}

fn molstar_edge_builder_directed_slots(
    atoms: &[usize],
    edges: &[(usize, usize)],
) -> Vec<(usize, usize)> {
    let mut bucket_sizes = vec![0usize; atoms.len()];
    for &(a, b) in edges {
        if a < atoms.len() {
            bucket_sizes[a] += 1;
        }
        if b < atoms.len() {
            bucket_sizes[b] += 1;
        }
    }

    let mut offsets = vec![0usize; atoms.len() + 1];
    let mut cursor = 0usize;
    for (index, size) in bucket_sizes.iter().enumerate() {
        offsets[index] = cursor;
        cursor += *size;
    }
    offsets[atoms.len()] = cursor;

    let mut bucket_fill = vec![0usize; atoms.len()];
    let mut slots = vec![None; cursor];
    for &(a, b) in edges {
        if a >= atoms.len() || b >= atoms.len() {
            continue;
        }
        let slot_ab = offsets[a] + bucket_fill[a];
        bucket_fill[a] += 1;
        let slot_ba = offsets[b] + bucket_fill[b];
        bucket_fill[b] += 1;
        slots[slot_ab] = Some((atoms[a], atoms[b]));
        slots[slot_ba] = Some((atoms[b], atoms[a]));
    }

    slots.into_iter().flatten().collect()
}

fn molstar_visible_renderable_bounding_sphere_with_structure(
    molecule: &Molecule,
    options: &MeshOptions,
    structure: &AtomicStructure,
) -> Option<BoundingSphere> {
    let spheres =
        molstar_visible_renderable_component_spheres_with_structure(molecule, options, structure)
            .into_iter()
            .map(|(_, sphere)| sphere)
            .collect::<Vec<_>>();
    (!spheres.is_empty()).then(|| Boundary::from_bounding_spheres(&spheres).sphere)
}

fn molstar_visible_renderable_component_spheres_with_structure(
    molecule: &Molecule,
    options: &MeshOptions,
    structure: &AtomicStructure,
) -> Vec<(&'static str, BoundingSphere)> {
    let selected = selected_visuals(&structure, options);
    let has_visual = |name: &str| selected.iter().any(|visual| visual == name);
    let mut spheres = Vec::new();

    let trace_padding = molstar_cartoon_uniform_trace_radius64(options);
    let tube_padding = molstar_trace_radius64(options);
    let bond_padding = MOLSTAR_BACKBONE_SIZE_FACTOR64 * molstar_radius_scale64(options);
    for group in &structure.symmetry_groups {
        let units = group
            .unit_ids
            .iter()
            .filter_map(|unit_id| structure.unit_by_id(*unit_id))
            .collect::<Vec<_>>();
        let Some(unit) = units.first().copied() else {
            continue;
        };
        let Some(unit_sphere) = molstar_unit_invariant_bounding_sphere(molecule, unit) else {
            continue;
        };
        if unit_sphere.radius <= 0.0 {
            continue;
        }
        let has_polymer = !unit.props.polymer_elements.is_empty();
        let has_nucleotide = !unit.props.nucleotide_elements.is_empty();
        let has_bonds =
            unit.props.intra_unit_bond_count > 0 || unit.props.inter_unit_bond_count > 0;

        if has_visual("polymer-trace") && has_polymer {
            spheres.push((
                "polymer-trace",
                molstar_units_transform_bounding_sphere(
                    &molstar_expand_bounding_sphere(&unit_sphere, trace_padding),
                    &units,
                ),
            ));
        }
        if has_visual("polymer-tube") && has_polymer {
            spheres.push((
                "polymer-tube",
                molstar_units_transform_bounding_sphere(
                    &molstar_expand_bounding_sphere(&unit_sphere, tube_padding),
                    &units,
                ),
            ));
        }
        if has_visual("polymer-gap") && !unit.props.gap_elements.is_empty() {
            spheres.push((
                "polymer-gap",
                molstar_units_transform_bounding_sphere(
                    &molstar_expand_bounding_sphere(&unit_sphere, trace_padding),
                    &units,
                ),
            ));
        }
        if (has_visual("nucleotide-ring")
            || has_visual("nucleotide-block")
            || has_visual("direction-wedge"))
            && has_nucleotide
        {
            spheres.push((
                "nucleotide",
                molstar_units_transform_bounding_sphere(
                    &molstar_expand_bounding_sphere(&unit_sphere, trace_padding),
                    &units,
                ),
            ));
        }
        if (has_visual("polymer-backbone-cylinder") || has_visual("polymer-backbone-sphere"))
            && has_polymer
        {
            spheres.push((
                "polymer-backbone",
                molstar_units_transform_bounding_sphere(
                    &molstar_expand_bounding_sphere(&unit_sphere, bond_padding),
                    &units,
                ),
            ));
        }
        if (has_visual("intra-bond")
            || has_visual("inter-bond")
            || has_visual("structure-intra-bond"))
            && has_bonds
        {
            spheres.push((
                "bond",
                molstar_units_transform_bounding_sphere(
                    &molstar_expand_bounding_sphere(&unit_sphere, bond_padding),
                    &units,
                ),
            ));
        }
        if (has_visual("element-sphere") || has_visual("structure-element-sphere"))
            && (!has_polymer
                || matches!(
                    options.representation,
                    Representation::Spacefill | Representation::BallAndStick
                ))
        {
            spheres.push((
                "element-sphere",
                molstar_units_transform_bounding_sphere(
                    &molstar_expand_bounding_sphere(
                        &unit_sphere,
                        unit_max_theme_size(molecule, unit) * molstar_radius_scale64(options)
                            + 0.05,
                    ),
                    &units,
                ),
            ));
        }
    }

    if (has_visual("carbohydrate-symbol")
        || has_visual("carbohydrate-link")
        || has_visual("carbohydrate-terminal-link"))
        && structure.boundary.sphere.radius > 0.0
    {
        spheres.push((
            "carbohydrate",
            molstar_expand_bounding_sphere(
                &structure.boundary.sphere,
                MOLSTAR_CARBOHYDRATE_SYMBOL_SIZE_FACTOR as f64,
            ),
        ));
    }

    spheres
}

fn molstar_unit_invariant_bounding_sphere(
    molecule: &Molecule,
    unit: &StructureUnit,
) -> Option<BoundingSphere> {
    match unit.kind {
        UnitKind::Atomic => {
            let mut positions = Vec::with_capacity(unit.elements.len());
            let mut radii = Vec::with_capacity(unit.elements.len());
            for &element in &unit.elements {
                let atom = molecule.atoms.get(element)?;
                positions.push(atom.position);
                radii.push(vdw_radius(&atom.type_symbol));
            }
            Some(Boundary::from_positions_and_radii(&positions, &radii).sphere)
        }
        UnitKind::Spheres | UnitKind::Gaussians => Some(unit.props.boundary.sphere.clone()),
    }
}

fn molstar_units_transform_bounding_sphere(
    invariant_sphere: &BoundingSphere,
    units: &[&StructureUnit],
) -> BoundingSphere {
    if units.is_empty() {
        return invariant_sphere.clone();
    }
    if units.len() == 1 {
        return molstar_transform_bounding_sphere(invariant_sphere, units[0].operator.transform);
    }

    if invariant_sphere.extrema.len() > 1 && units.len() <= 14 {
        let mut positions = Vec::with_capacity(invariant_sphere.extrema.len() * units.len());
        for unit in units {
            for &point in &invariant_sphere.extrema {
                positions.push(unit.operator.transform.apply(point));
            }
        }
        Boundary::from_positions(&positions).sphere
    } else {
        let centers_and_radii = units
            .iter()
            .map(|unit| {
                (
                    unit.operator.transform.apply(invariant_sphere.center),
                    invariant_sphere.radius
                        * molstar_transform_max_scale_on_axis(unit.operator.transform),
                )
            })
            .collect::<Vec<_>>();
        let centers = centers_and_radii
            .iter()
            .map(|(center, _)| *center)
            .collect::<Vec<_>>();
        let radii = centers_and_radii
            .iter()
            .map(|(_, radius)| *radius)
            .collect::<Vec<_>>();
        Boundary::from_positions_and_radii(&centers, &radii).sphere
    }
}

fn molstar_transform_bounding_sphere(
    sphere: &BoundingSphere,
    transform: Transform,
) -> BoundingSphere {
    if transform.is_identity() {
        return sphere.clone();
    }
    let scale = molstar_transform_max_scale_on_axis(transform);
    BoundingSphere {
        center: transform.apply(sphere.center),
        radius: sphere.radius * scale,
        extrema: sphere
            .extrema
            .iter()
            .map(|&point| transform.apply(point))
            .collect(),
        extrema64: sphere
            .extrema64
            .iter()
            .map(|&[x, y, z]| {
                let point = transform.apply(Vec3::new(x as f32, y as f32, z as f32));
                [point.x as f64, point.y as f64, point.z as f64]
            })
            .collect(),
        center64: None,
        radius64: sphere
            .radius64()
            .is_finite()
            .then_some(sphere.radius64() * scale as f64),
    }
}

fn molstar_transform_max_scale_on_axis(transform: Transform) -> f32 {
    let sx = transform.m[0][0] * transform.m[0][0]
        + transform.m[1][0] * transform.m[1][0]
        + transform.m[2][0] * transform.m[2][0];
    let sy = transform.m[0][1] * transform.m[0][1]
        + transform.m[1][1] * transform.m[1][1]
        + transform.m[2][1] * transform.m[2][1];
    let sz = transform.m[0][2] * transform.m[0][2]
        + transform.m[1][2] * transform.m[1][2]
        + transform.m[2][2] * transform.m[2][2];
    sx.max(sy).max(sz).sqrt()
}

fn unit_max_theme_size(molecule: &Molecule, unit: &crate::model::StructureUnit) -> f64 {
    match unit.kind {
        UnitKind::Atomic => unit
            .elements
            .iter()
            .filter_map(|&atom_index| molecule.atoms.get(atom_index))
            .map(|atom| vdw_radius64(&atom.type_symbol))
            .fold(0.0, f64::max),
        UnitKind::Spheres | UnitKind::Gaussians => 1.0,
    }
}

fn molstar_expand_bounding_sphere(sphere: &BoundingSphere, delta: f64) -> BoundingSphere {
    let delta32 = delta as f32;
    let mut out = BoundingSphere {
        center: sphere.center,
        radius: sphere.radius + delta32,
        extrema: Vec::new(),
        center64: sphere.center64,
        radius64: Some(sphere.radius64() + delta),
        extrema64: Vec::new(),
    };
    if sphere.radius < 1e-12 || sphere.extrema.len() <= 1 {
        return out;
    }

    let moments_axes = PrincipalAxes::calculate_moments_axes64(&sphere.extrema);
    let axes = PrincipalAxes::calculate_normalized_axes64(&moments_axes);
    let center64 = sphere.center64();
    let radius64 = out.radius64();
    let delta64 = delta;
    let dir_a64 = [
        axes.dir_a[0] * delta64,
        axes.dir_a[1] * delta64,
        axes.dir_a[2] * delta64,
    ];
    let dir_b64 = [
        axes.dir_b[0] * delta64,
        axes.dir_b[1] * delta64,
        axes.dir_b[2] * delta64,
    ];
    let dir_c64 = [
        axes.dir_c[0] * delta64,
        axes.dir_c[1] * delta64,
        axes.dir_c[2] * delta64,
    ];

    let normalize64 = |v: [f64; 3]| {
        let length = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if length > 0.0 {
            [v[0] / length, v[1] / length, v[2] / length]
        } else {
            v
        }
    };
    let dot64 = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let add_scaled64 = |mut point: [f64; 3], direction: [f64; 3], sign: f64| {
        point[0] += direction[0] * sign;
        point[1] += direction[1] * sign;
        point[2] += direction[2] * sign;
        point
    };
    let distance64 = |a: [f64; 3], b: [f64; 3]| {
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    };
    let source_extrema64 = if sphere.extrema64.len() == sphere.extrema.len() {
        sphere.extrema64.clone()
    } else {
        sphere
            .extrema
            .iter()
            .map(|point| [point.x as f64, point.y as f64, point.z as f64])
            .collect()
    };

    out.extrema64 = source_extrema64
        .iter()
        .copied()
        .map(|extreme64| {
            let mut direction = normalize64([
                extreme64[0] - center64[0],
                extreme64[1] - center64[1],
                extreme64[2] - center64[2],
            ]);
            let mut point = extreme64;

            point = add_scaled64(
                point,
                dir_a64,
                if dot64(direction, dir_a64) < 0.0 {
                    -1.0
                } else {
                    1.0
                },
            );
            point = add_scaled64(
                point,
                dir_b64,
                if dot64(direction, dir_b64) < 0.0 {
                    -1.0
                } else {
                    1.0
                },
            );
            point = add_scaled64(
                point,
                dir_c64,
                if dot64(direction, dir_c64) < 0.0 {
                    -1.0
                } else {
                    1.0
                },
            );

            if distance64(center64, point) > radius64 {
                if sphere.extrema.len() >= 14 {
                    direction = normalize64([
                        point[0] - center64[0],
                        point[1] - center64[1],
                        point[2] - center64[2],
                    ]);
                }
                point = [
                    center64[0] + direction[0] * radius64,
                    center64[1] + direction[1] * radius64,
                    center64[2] + direction[2] * radius64,
                ];
            }

            point
        })
        .collect();
    out.extrema = out
        .extrema64
        .iter()
        .map(|point| Vec3::new(point[0] as f32, point[1] as f32, point[2] as f32))
        .collect();
    out
}

fn flatten_semantic_render_objects_with_visible_bounding_sphere_and_stats(
    objects: &[SemanticRenderObject],
    _molecule: &Molecule,
    options: &MeshOptions,
    collect_visible_bounding_sphere: bool,
) -> (Mesh, Option<BoundingSphere>, Vec<RenderObjectMeshStats>) {
    let cylinder_radial_segments = molstar_export_cylinder_radial_segments(
        objects
            .iter()
            .map(|object| render_object_export_cylinder_count(&object.object))
            .sum::<usize>(),
    );
    let (estimate, plans) = render_objects_mesh_plan(
        objects.iter().map(|object| &object.object),
        options,
        cylinder_radial_segments,
    );
    let mut state = MeshBuilderState::with_capacity(estimate, objects.len(), true);
    let mut object_spheres = Vec::with_capacity(if collect_visible_bounding_sphere {
        objects.len()
    } else {
        0
    });
    let mut cylinder_cache = CylinderPrimitiveCache::default();
    let mut curve_scratch = CurveSegmentScratch::default();
    let mut object_stats = Vec::with_capacity(objects.len());
    for (object, plan) in objects.iter().zip(&plans) {
        state.set_current_group(object.group_id);
        state.set_current_material(object.material);
        state.set_current_section(Some(object.visual));
        let vertex_start = state.mesh.vertices.len();
        let face_start = state.mesh.faces.len();
        append_render_object_to_mesh(
            &mut state.mesh,
            &object.object,
            options,
            cylinder_radial_segments,
            &mut cylinder_cache,
            &mut curve_scratch,
            Some(plan),
        );
        state.mark_appended(vertex_start, face_start);
        let vertex_count = state.mesh.vertices.len().saturating_sub(vertex_start);
        let face_count = state.mesh.faces.len().saturating_sub(face_start);
        object_stats.push(render_object_mesh_stats(
            &object.object,
            vertex_count,
            face_count,
        ));
        if collect_visible_bounding_sphere && face_count > 0 && vertex_count > 0 {
            let boundary = Boundary::from_positions(&state.mesh.vertices[vertex_start..]);
            if boundary.sphere.radius > 0.0 {
                object_spheres.push(boundary.sphere);
            }
        }
    }
    let visible_sphere = if object_spheres.is_empty() {
        None
    } else {
        Some(Boundary::from_bounding_spheres(&object_spheres).sphere)
    };
    (state.into_mesh(), visible_sphere, object_stats)
}

fn render_object_mesh_stats(
    object: &RenderObject,
    vertex_count: usize,
    face_count: usize,
) -> RenderObjectMeshStats {
    if matches!(object, RenderObject::Cylinder { .. }) && face_count > 0 {
        let radial_segments =
            molstar_export_cylinder_radial_segments(render_object_export_cylinder_count(object));
        return RenderObjectMeshStats {
            draw_count: radial_segments * 4 * 3,
            vertex_count: (radial_segments + 1) * 4,
            group_count: 1,
        };
    }
    RenderObjectMeshStats {
        draw_count: face_count.saturating_mul(3),
        vertex_count,
        group_count: usize::from(face_count > 0),
    }
}

fn render_object_mesh_stats_from_estimates(
    objects: &[SemanticRenderObject],
    options: &MeshOptions,
) -> Vec<RenderObjectMeshStats> {
    let cylinder_radial_segments = molstar_export_cylinder_radial_segments(
        objects
            .iter()
            .map(|object| render_object_export_cylinder_count(&object.object))
            .sum(),
    );
    objects
        .iter()
        .map(|object| {
            let estimate = object
                .object
                .mesh_estimate(options, cylinder_radial_segments);
            render_object_mesh_stats(&object.object, estimate.vertices, estimate.faces)
        })
        .collect()
}

#[allow(dead_code)]
fn flatten_render_objects(
    objects: &[RenderObject],
    _molecule: &Molecule,
    options: &MeshOptions,
) -> Mesh {
    let groups = (0..objects.len()).collect::<Vec<_>>();
    flatten_render_objects_with_groups(objects, &groups, _molecule, options)
}

fn flatten_render_objects_with_groups(
    objects: &[RenderObject],
    groups: &[usize],
    _molecule: &Molecule,
    options: &MeshOptions,
) -> Mesh {
    flatten_render_objects_with_groups_and_visible_bounding_sphere(
        objects, groups, None, None, _molecule, options, false,
    )
    .0
}

fn flatten_render_objects_with_groups_and_visible_bounding_sphere(
    objects: &[RenderObject],
    groups: &[usize],
    materials: Option<&[Option<MeshMaterial>]>,
    section_keys: Option<&[String]>,
    _molecule: &Molecule,
    options: &MeshOptions,
    collect_visible_bounding_sphere: bool,
) -> (Mesh, Option<BoundingSphere>) {
    let cylinder_radial_segments = molstar_export_cylinder_radial_segments(
        objects
            .iter()
            .map(render_object_export_cylinder_count)
            .sum::<usize>(),
    );
    let (estimate, plans) =
        render_objects_mesh_plan(objects.iter(), options, cylinder_radial_segments);
    let mut state = MeshBuilderState::with_capacity(
        estimate,
        section_keys.map_or(0, <[String]>::len),
        materials.is_some(),
    );
    let mut object_spheres = Vec::with_capacity(if collect_visible_bounding_sphere {
        objects.len()
    } else {
        0
    });
    let mut cylinder_cache = CylinderPrimitiveCache::default();
    let mut curve_scratch = CurveSegmentScratch::default();
    for (index, (object, plan)) in objects.iter().zip(&plans).enumerate() {
        let group = groups.get(index).copied().unwrap_or(index);
        state.set_current_group(group);
        state.set_current_material(
            materials.and_then(|materials| materials.get(index).copied().flatten()),
        );
        state
            .set_current_section(section_keys.and_then(|keys| keys.get(index).map(String::as_str)));
        let vertex_start = state.mesh.vertices.len();
        let face_start = state.mesh.faces.len();
        append_render_object_to_mesh(
            &mut state.mesh,
            object,
            options,
            cylinder_radial_segments,
            &mut cylinder_cache,
            &mut curve_scratch,
            Some(plan),
        );
        state.mark_appended(vertex_start, face_start);
        if collect_visible_bounding_sphere
            && state.mesh.faces.len() > face_start
            && state.mesh.vertices.len() > vertex_start
        {
            let boundary = Boundary::from_positions(&state.mesh.vertices[vertex_start..]);
            if boundary.sphere.radius > 0.0 {
                object_spheres.push(boundary.sphere);
            }
        }
    }
    let visible_sphere = if object_spheres.is_empty() {
        None
    } else {
        Some(Boundary::from_bounding_spheres(&object_spheres).sphere)
    };
    (state.into_mesh(), visible_sphere)
}

fn append_render_object_to_mesh(
    mesh: &mut Mesh,
    object: &RenderObject,
    options: &MeshOptions,
    cylinder_radial_segments: usize,
    cylinder_cache: &mut CylinderPrimitiveCache,
    curve_scratch: &mut CurveSegmentScratch,
    plan: Option<&RenderObjectMeshPlan>,
) {
    #[cfg(debug_assertions)]
    let (vertex_start, face_start, expected) = (
        mesh.vertices.len(),
        mesh.faces.len(),
        plan.map_or_else(
            || object.mesh_estimate(options, cylinder_radial_segments),
            |plan| plan.estimate,
        ),
    );
    match object {
        RenderObject::Sphere { center, radius } => {
            add_sphere_with_radius64(mesh, *center, *radius, options.sphere_detail);
        }
        RenderObject::Cylinder { start, end, radius } => {
            add_molstar_export_bond_cylinder(
                mesh,
                *start,
                *end,
                *radius,
                cylinder_radial_segments,
                cylinder_cache,
            );
        }
        RenderObject::LinkCylinder { start, end, radius } => {
            let midpoint = molstar_link_midpoint_buffer(*start, *end);
            add_molstar_buffered_open_cylinder_cached(
                mesh,
                *start,
                midpoint,
                *radius,
                options.radial_segments.max(3),
                cylinder_cache,
            );
        }
        RenderObject::LinkCylinderWithSegments {
            start,
            end,
            radius,
            radial_segments,
        } => {
            let midpoint = molstar_link_midpoint_buffer(*start, *end);
            add_molstar_buffered_open_cylinder_with_radius64_cached(
                mesh,
                *start,
                midpoint,
                *radius,
                (*radial_segments).max(3),
                cylinder_cache,
            );
        }
        RenderObject::Tube { points, radius } => {
            add_tube_path(mesh, points, *radius, options.radial_segments.max(3));
        }
        RenderObject::DashedTube { points, radius } => {
            if let Some(samples) = plan.and_then(|plan| plan.dashed_samples.as_deref()) {
                add_dashed_tube_samples_cached(
                    mesh,
                    samples,
                    *radius,
                    options.radial_segments.max(3),
                    cylinder_cache,
                );
            } else {
                add_dashed_tube_path_cached(
                    mesh,
                    points,
                    *radius,
                    options.radial_segments.max(3),
                    cylinder_cache,
                );
            }
        }
        RenderObject::FixedCountDashedCylinder {
            start,
            end,
            radius,
            length_scale,
            segment_count,
        } => add_fixed_count_dashed_cylinder_cached(
            mesh,
            *start,
            *end,
            *radius,
            options.radial_segments.max(3),
            *length_scale,
            *segment_count,
            cylinder_cache,
        ),
        RenderObject::Ribbon {
            points,
            width,
            thickness,
        } => add_ribbon(mesh, points, *width, *thickness, options.linear_segments),
        RenderObject::Sheet {
            points,
            width,
            thickness,
            arrow_height,
            start_cap,
            end_cap,
        } => add_sheet(
            mesh,
            points,
            *width,
            *thickness,
            *arrow_height,
            *start_cap,
            *end_cap,
            options.linear_segments,
        ),
        RenderObject::OrientedRibbon {
            centers,
            normals,
            width,
            thickness,
            profile,
            start_cap,
            end_cap,
            round_cap,
        } => add_oriented_ribbon_with_profile(
            mesh,
            centers,
            normals,
            *width,
            *thickness,
            *profile,
            *start_cap,
            *end_cap,
            *round_cap,
            options.linear_segments,
            options.radial_segments,
        ),
        RenderObject::PolymerTraceSegment {
            controls,
            widths,
            heights,
            tension,
            shift,
            overhang_width,
            kind,
            start_cap,
            end_cap,
            initial,
            final_residue,
            swap_normal_binormal,
        } => match kind {
            PolymerTraceSegmentKind::Ribbon {
                arrow_height,
                swap_width_height,
            } => add_curve_segment_ribbon(
                mesh,
                controls,
                *widths,
                *heights,
                *tension,
                *shift,
                *overhang_width,
                *arrow_height,
                *initial,
                *final_residue,
                *swap_normal_binormal,
                *swap_width_height,
                options.linear_segments,
                curve_scratch,
            ),
            PolymerTraceSegmentKind::Tube { profile, round_cap } => add_curve_segment_tube(
                mesh,
                controls,
                *widths,
                *heights,
                *tension,
                *shift,
                *overhang_width,
                *start_cap,
                *end_cap,
                *round_cap,
                *initial,
                *final_residue,
                *swap_normal_binormal,
                options.linear_segments,
                options.radial_segments,
                *profile,
                curve_scratch,
            ),
            PolymerTraceSegmentKind::Sheet { arrow_height } => add_curve_segment_sheet(
                mesh,
                controls,
                *widths,
                *heights,
                *tension,
                *shift,
                *overhang_width,
                *arrow_height,
                *start_cap,
                *end_cap,
                *initial,
                *final_residue,
                *swap_normal_binormal,
                options.linear_segments,
                curve_scratch,
            ),
        },
        RenderObject::NucleotideRing {
            center,
            normal,
            radius,
            base,
            detail,
            radial_segments,
        } => add_nucleotide_ring(
            mesh,
            *center,
            *normal,
            *radius,
            *base,
            *detail,
            *radial_segments,
            cylinder_cache,
        ),
        RenderObject::NucleotideBlock {
            geometry,
            radius,
            width,
            depth,
            radial_segments,
        } => add_nucleotide_block(
            mesh,
            *geometry,
            *radius,
            *width,
            *depth,
            *radial_segments,
            cylinder_cache,
        ),
        RenderObject::DirectionWedge {
            center,
            tangent,
            up,
            size,
        } => add_direction_wedge(mesh, *center, *tangent, *up, *size),
        RenderObject::CarbohydrateSymbol {
            center,
            normal,
            direction,
            shape,
            part,
        } => add_carbohydrate_symbol(mesh, *center, *normal, *direction, *shape, *part),
        RenderObject::Ellipsoid { center, axes } => {
            add_ellipsoid(mesh, *center, *axes, options.sphere_detail)
        }
    }
    #[cfg(debug_assertions)]
    {
        debug_assert_eq!(
            mesh.vertices.len() - vertex_start,
            expected.vertices,
            "render object vertex estimate drifted for {object:?}"
        );
        debug_assert_eq!(
            mesh.faces.len() - face_start,
            expected.faces,
            "render object face estimate drifted for {object:?}"
        );
    }
}

fn render_object_export_cylinder_count(object: &RenderObject) -> usize {
    match object {
        // Mol* solid bond cylinder impostors store two half-cylinders per link.
        RenderObject::Cylinder { .. } => 2,
        _ => 0,
    }
}

fn molstar_export_cylinder_radial_segments(cylinder_count: usize) -> usize {
    if cylinder_count < 2_000 {
        36
    } else if cylinder_count < 20_000 {
        24
    } else {
        12
    }
}

fn add_molstar_export_bond_cylinder(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    radial_segments: usize,
    cylinder_cache: &mut CylinderPrimitiveCache,
) {
    let midpoint = molstar_link_midpoint_buffer(start, end);
    add_molstar_buffered_open_cylinder_cached(
        mesh,
        start,
        midpoint,
        radius,
        radial_segments,
        cylinder_cache,
    );
    add_molstar_buffered_open_cylinder_cached(
        mesh,
        midpoint,
        end,
        radius,
        radial_segments,
        cylinder_cache,
    );
}

fn molstar_link_midpoint_buffer(start: Vec3, end: Vec3) -> Vec3 {
    Vec3::new(
        ((start.x as f64 + end.x as f64) * 0.5) as f32,
        ((start.y as f64 + end.y as f64) * 0.5) as f32,
        ((start.z as f64 + end.z as f64) * 0.5) as f32,
    )
}

#[derive(Default)]
struct MeshBuilderState {
    current_group: Option<usize>,
    current_material: Option<MeshMaterial>,
    current_section: Option<String>,
    mesh: Mesh,
}

impl MeshBuilderState {
    fn with_capacity(
        estimate: RenderObjectMeshEstimate,
        section_count: usize,
        include_materials: bool,
    ) -> Self {
        let mut mesh = mesh_with_capacity(estimate);
        mesh.sections = Vec::with_capacity(section_count);
        if include_materials {
            mesh.face_materials = Vec::with_capacity(estimate.faces);
        }
        Self {
            current_group: None,
            current_material: None,
            current_section: None,
            mesh,
        }
    }

    fn set_current_group(&mut self, group: usize) {
        self.current_group = Some(group);
    }

    fn set_current_material(&mut self, material: Option<MeshMaterial>) {
        self.current_material = material;
    }

    fn set_current_section(&mut self, section: Option<&str>) {
        self.current_section = section.map(str::to_string);
    }

    fn mark_appended(&mut self, vertex_start: usize, face_start: usize) {
        let group = self
            .current_group
            .expect("MeshBuilderState current group must be set before append");
        debug_assert_eq!(self.mesh.vertex_groups.len(), vertex_start);
        debug_assert_eq!(self.mesh.face_groups.len(), face_start);
        debug_assert_eq!(self.mesh.vertices.len(), self.mesh.normals.len());
        let new_vertices = self.mesh.vertices.len().saturating_sub(vertex_start);
        let new_faces = self.mesh.faces.len().saturating_sub(face_start);
        self.mesh
            .vertex_groups
            .extend(std::iter::repeat_n(group, new_vertices));
        self.mesh
            .face_groups
            .extend(std::iter::repeat_n(group, new_faces));
        if let Some(material) = self.current_material {
            self.mesh
                .face_materials
                .extend(std::iter::repeat_n(material, new_faces));
        } else if !self.mesh.face_materials.is_empty() {
            self.mesh.face_materials.extend(std::iter::repeat_n(
                MeshMaterial::opaque(0xfafafa),
                new_faces,
            ));
        }
        if new_faces > 0 {
            self.mesh.group_count = self.mesh.group_count.max(group + 1);
            self.mark_section(vertex_start, face_start);
        }
    }

    fn mark_section(&mut self, vertex_start: usize, face_start: usize) {
        let Some(key) = &self.current_section else {
            return;
        };
        let vertex_end = self.mesh.vertices.len();
        let face_end = self.mesh.faces.len();
        if let Some(section) = self.mesh.sections.last_mut() {
            if section.key == *key
                && section.vertex_end == vertex_start
                && section.face_end == face_start
            {
                section.vertex_end = vertex_end;
                section.face_end = face_end;
                return;
            }
        }
        self.mesh.sections.push(crate::model::MeshSection {
            key: key.clone(),
            vertex_start,
            vertex_end,
            face_start,
            face_end,
        });
    }

    fn into_mesh(self) -> Mesh {
        self.mesh
    }
}

fn backbone_residues(molecule: &Molecule, structure: &AtomicStructure) -> Vec<TraceResidue> {
    #[derive(Clone, Copy, Debug)]
    struct DerivedTraceAtoms {
        molecule_type: MoleculeType,
        polymer_type: PolymerType,
        trace: Option<Vec3>,
        direction: Option<Vec3>,
    }

    #[derive(Clone, Debug)]
    struct ResidueAtoms {
        chain: String,
        residue: String,
        seq: i32,
        insertion_code: String,
        trace: Option<Vec3>,
        carbonyl_c: Option<Vec3>,
        carbonyl_o: Option<Vec3>,
        c1: Option<Vec3>,
        c3: Option<Vec3>,
        c4: Option<Vec3>,
        o3: Option<Vec3>,
        nucleotide_atoms: NucleotideAtoms,
    }

    let hierarchy = &structure.model.hierarchy;
    let derived_residues = hierarchy
        .residues
        .iter()
        .enumerate()
        .filter_map(|(residue_index, residue)| {
            let chain = hierarchy.chains.get(residue.chain_index)?;
            let seq = residue.label_seq_id.trim().parse::<i32>().ok()?;
            let atom_position = |index: Option<usize>| {
                index
                    .and_then(|atom_index| hierarchy.atoms.get(atom_index))
                    .map(|atom| atom.position)
            };
            let trace = hierarchy
                .derived
                .residue
                .trace_element_index
                .get(residue_index)
                .and_then(|index| atom_position(*index));
            let from = hierarchy
                .derived
                .residue
                .direction_from_element_index
                .get(residue_index)
                .and_then(|index| atom_position(*index));
            let to = hierarchy
                .derived
                .residue
                .direction_to_element_index
                .get(residue_index)
                .and_then(|index| atom_position(*index));
            let direction = match (from, to) {
                (Some(from), Some(to)) => {
                    let direction = (to - from).normalized();
                    (direction.length() > 0.000_001).then_some(direction)
                }
                _ => None,
            };
            Some((
                chain.id.clone(),
                seq,
                residue.insertion_code.clone(),
                DerivedTraceAtoms {
                    molecule_type: hierarchy.derived.residue.molecule_type[residue_index],
                    polymer_type: hierarchy.derived.residue.polymer_type[residue_index],
                    trace,
                    direction,
                },
            ))
        })
        .collect::<Vec<_>>();

    let mut residues = Vec::<ResidueAtoms>::new();
    for atom in &molecule.atoms {
        let seq = atom
            .residue_seq
            .trim()
            .parse::<i32>()
            .unwrap_or(atom.id as i32);
        let index = residues
            .iter()
            .position(|residue| {
                residue.chain == atom.chain
                    && residue.seq == seq
                    && residue.insertion_code == atom.insertion_code
            })
            .unwrap_or_else(|| {
                residues.push(ResidueAtoms {
                    chain: atom.chain.clone(),
                    residue: atom.residue.clone(),
                    seq,
                    insertion_code: atom.insertion_code.clone(),
                    trace: None,
                    carbonyl_c: None,
                    carbonyl_o: None,
                    c1: None,
                    c3: None,
                    c4: None,
                    o3: None,
                    nucleotide_atoms: NucleotideAtoms::default(),
                });
                residues.len() - 1
            });
        let residue = &mut residues[index];
        let atom_name = atom.name.trim();
        let is_nucleotide = is_nucleotide_residue(&residue.residue);
        residue
            .nucleotide_atoms
            .record_atom(atom_name, atom.position);
        match atom_name {
            "CA" if !is_nucleotide => residue.trace = Some(atom.position),
            "O3'" | "O3*" if is_nucleotide => {
                residue.trace = Some(atom.position);
                residue.o3 = Some(atom.position);
                residue.nucleotide_atoms.set_trace(atom.position);
            }
            "P" if residue.trace.is_none() => residue.trace = Some(atom.position),
            "C4'" | "C4*" => {
                residue.c4 = Some(atom.position);
                if residue.trace.is_none() {
                    residue.trace = Some(atom.position);
                }
            }
            "C1'" | "C1*" => residue.c1 = Some(atom.position),
            "C3'" | "C3*" => residue.c3 = Some(atom.position),
            "C" => residue.carbonyl_c = Some(atom.position),
            "O" | "OXT" => residue.carbonyl_o = Some(atom.position),
            _ => {}
        }
    }

    let mut out = residues
        .into_iter()
        .filter_map(|residue| {
            let derived = derived_residues
                .iter()
                .find(|(chain, seq, insertion_code, _)| {
                    *chain == residue.chain
                        && *seq == residue.seq
                        && *insertion_code == residue.insertion_code
                })
                .map(|(_, _, _, derived)| *derived);
            let position = derived
                .and_then(|derived| derived.trace)
                .or(residue.trace)?;
            let is_nucleotide = derived.is_some_and(|derived| {
                matches!(
                    derived.molecule_type,
                    MoleculeType::Rna | MoleculeType::Dna | MoleculeType::Pna
                )
            }) || is_nucleotide_residue(&residue.residue);
            let direction = if let Some(derived) = derived {
                if derived.polymer_type == PolymerType::None {
                    None
                } else {
                    derived.direction
                }
            } else if is_nucleotide {
                let (from, to) = if is_dna_residue(&residue.residue) {
                    (residue.c3, residue.c1)
                } else {
                    (residue.c4, residue.c3)
                };
                match (from, to) {
                    (Some(from), Some(to)) => {
                        let direction = (to - from).normalized();
                        (direction.length() > 0.000_001).then_some(direction)
                    }
                    _ => None,
                }
            } else {
                match (residue.carbonyl_c, residue.carbonyl_o) {
                    (Some(from), Some(to)) => {
                        let direction = (to - from).normalized();
                        (direction.length() > 0.000_001).then_some(direction)
                    }
                    _ => None,
                }
            };
            let mut nucleotide_atoms = residue.nucleotide_atoms;
            if is_nucleotide && nucleotide_atoms.trace.is_none() {
                nucleotide_atoms.set_trace(position);
            }
            let nucleotide_base_kind = nucleotide_base_kind(&residue.residue, nucleotide_atoms);
            let base_atoms = nucleotide_base_atoms(nucleotide_atoms, nucleotide_base_kind);
            let base_center = centroid(&base_atoms);
            let base_normal = nucleotide_base_normal(nucleotide_atoms, nucleotide_base_kind);
            Some(TraceResidue {
                chain: residue.chain,
                residue: residue.residue,
                seq: residue.seq,
                insertion_code: residue.insertion_code,
                position,
                direction,
                initial: false,
                final_residue: false,
                sec_struc_first: false,
                sec_struc_last: false,
                is_nucleotide,
                base_center,
                base_normal,
                nucleotide_atoms: is_nucleotide.then_some(nucleotide_atoms),
                nucleotide_base_kind,
            })
        })
        .collect::<Vec<_>>();
    out.sort_by(|a, b| {
        a.chain
            .cmp(&b.chain)
            .then_with(|| a.seq.cmp(&b.seq))
            .then_with(|| a.insertion_code.cmp(&b.insertion_code))
    });
    let len = out.len();
    for i in 0..len {
        out[i].initial = i == 0 || out[i - 1].chain != out[i].chain;
        out[i].final_residue = i + 1 == len || out[i + 1].chain != out[i].chain;
    }
    out
}

fn backbone_residues_from_atoms(molecule: &Molecule) -> Vec<TraceResidue> {
    #[derive(Clone, Debug)]
    struct ResidueAtoms {
        chain: String,
        residue: String,
        seq: i32,
        insertion_code: String,
        trace: Option<Vec3>,
        carbonyl_c: Option<Vec3>,
        carbonyl_o: Option<Vec3>,
        c1: Option<Vec3>,
        c3: Option<Vec3>,
        c4: Option<Vec3>,
        o3: Option<Vec3>,
        nucleotide_atoms: NucleotideAtoms,
    }

    let mut residues = Vec::<ResidueAtoms>::new();
    for atom in &molecule.atoms {
        let seq = atom
            .residue_seq
            .trim()
            .parse::<i32>()
            .unwrap_or(atom.id as i32);
        let index = residues
            .iter()
            .position(|residue| {
                residue.chain == atom.chain
                    && residue.seq == seq
                    && residue.insertion_code == atom.insertion_code
            })
            .unwrap_or_else(|| {
                residues.push(ResidueAtoms {
                    chain: atom.chain.clone(),
                    residue: atom.residue.clone(),
                    seq,
                    insertion_code: atom.insertion_code.clone(),
                    trace: None,
                    carbonyl_c: None,
                    carbonyl_o: None,
                    c1: None,
                    c3: None,
                    c4: None,
                    o3: None,
                    nucleotide_atoms: NucleotideAtoms::default(),
                });
                residues.len() - 1
            });
        let residue = &mut residues[index];
        let atom_name = atom.name.trim();
        let is_nucleotide = is_nucleotide_residue(&residue.residue);
        residue
            .nucleotide_atoms
            .record_atom(atom_name, atom.position);
        match atom_name {
            "CA" if !is_nucleotide => residue.trace = Some(atom.position),
            "O3'" | "O3*" if is_nucleotide => {
                residue.trace = Some(atom.position);
                residue.o3 = Some(atom.position);
                residue.nucleotide_atoms.set_trace(atom.position);
            }
            "P" if residue.trace.is_none() => residue.trace = Some(atom.position),
            "C4'" | "C4*" => {
                residue.c4 = Some(atom.position);
                if residue.trace.is_none() {
                    residue.trace = Some(atom.position);
                }
            }
            "C1'" | "C1*" => residue.c1 = Some(atom.position),
            "C3'" | "C3*" => residue.c3 = Some(atom.position),
            "C" => residue.carbonyl_c = Some(atom.position),
            "O" | "OXT" => residue.carbonyl_o = Some(atom.position),
            _ => {}
        }
    }

    let mut out = residues
        .into_iter()
        .filter_map(|residue| {
            let position = residue.trace?;
            let is_nucleotide = is_nucleotide_residue(&residue.residue);
            let direction = if is_nucleotide {
                let (from, to) = if is_dna_residue(&residue.residue) {
                    (residue.c3, residue.c1)
                } else {
                    (residue.c4, residue.c3)
                };
                match (from, to) {
                    (Some(from), Some(to)) => {
                        let direction = (to - from).normalized();
                        (direction.length() > 0.000_001).then_some(direction)
                    }
                    _ => None,
                }
            } else {
                match (residue.carbonyl_c, residue.carbonyl_o) {
                    (Some(from), Some(to)) => {
                        let direction = (to - from).normalized();
                        (direction.length() > 0.000_001).then_some(direction)
                    }
                    _ => None,
                }
            };
            let mut nucleotide_atoms = residue.nucleotide_atoms;
            if is_nucleotide && nucleotide_atoms.trace.is_none() {
                nucleotide_atoms.set_trace(position);
            }
            let nucleotide_base_kind = nucleotide_base_kind(&residue.residue, nucleotide_atoms);
            let base_atoms = nucleotide_base_atoms(nucleotide_atoms, nucleotide_base_kind);
            let base_center = centroid(&base_atoms);
            let base_normal = nucleotide_base_normal(nucleotide_atoms, nucleotide_base_kind);
            Some(TraceResidue {
                chain: residue.chain,
                residue: residue.residue,
                seq: residue.seq,
                insertion_code: residue.insertion_code,
                position,
                direction,
                initial: false,
                final_residue: false,
                sec_struc_first: false,
                sec_struc_last: false,
                is_nucleotide,
                base_center,
                base_normal,
                nucleotide_atoms: is_nucleotide.then_some(nucleotide_atoms),
                nucleotide_base_kind,
            })
        })
        .collect::<Vec<_>>();
    out.sort_by(|a, b| {
        a.chain
            .cmp(&b.chain)
            .then_with(|| a.seq.cmp(&b.seq))
            .then_with(|| a.insertion_code.cmp(&b.insertion_code))
    });
    let len = out.len();
    for i in 0..len {
        out[i].initial = i == 0 || out[i - 1].chain != out[i].chain;
        out[i].final_residue = i + 1 == len || out[i + 1].chain != out[i].chain;
    }
    out
}

fn is_nucleotide_residue(residue: &str) -> bool {
    matches!(
        residue.trim().to_ascii_uppercase().as_str(),
        "A" | "C"
            | "G"
            | "U"
            | "T"
            | "I"
            | "DA"
            | "DC"
            | "DG"
            | "DT"
            | "DU"
            | "DI"
            | "APN"
            | "GPN"
            | "CPN"
            | "TPN"
    )
}

fn is_dna_residue(residue: &str) -> bool {
    matches!(
        residue.trim().to_ascii_uppercase().as_str(),
        "DA" | "DC" | "DG" | "DT" | "DU" | "DI"
    )
}

fn nucleotide_base_kind(residue: &str, atoms: NucleotideAtoms) -> Option<NucleotideBaseKind> {
    let residue = residue.trim().to_ascii_uppercase();
    if matches!(
        residue.as_str(),
        "A" | "G" | "I" | "DA" | "DG" | "DI" | "APN" | "GPN"
    ) {
        return Some(NucleotideBaseKind::Purine);
    }
    if matches!(
        residue.as_str(),
        "C" | "T" | "U" | "DC" | "DT" | "DU" | "CPN" | "TPN"
    ) {
        return Some(NucleotideBaseKind::Pyrimidine);
    }
    match (atoms.c4, atoms.n9) {
        (Some(c4), Some(n9)) if c4.distance(n9) < 1.6 => Some(NucleotideBaseKind::Purine),
        _ => Some(NucleotideBaseKind::Pyrimidine),
    }
}

fn centroid(points: &[Vec3]) -> Option<Vec3> {
    if points.is_empty() {
        return None;
    }
    Some(
        points
            .iter()
            .copied()
            .fold(Vec3::default(), |sum, point| sum + point)
            / points.len() as f32,
    )
}

fn nucleotide_base_atoms(atoms: NucleotideAtoms, kind: Option<NucleotideBaseKind>) -> Vec<Vec3> {
    match kind.unwrap_or_else(|| {
        nucleotide_base_kind("", atoms).unwrap_or(NucleotideBaseKind::Pyrimidine)
    }) {
        NucleotideBaseKind::Purine => [
            atoms.n1,
            atoms.c2,
            atoms.n3,
            atoms.c4,
            atoms.c5.or(atoms.n5),
            atoms.c6,
            atoms.n7.or(atoms.c7),
            atoms.c8,
            atoms.n9,
        ]
        .into_iter()
        .flatten()
        .collect(),
        NucleotideBaseKind::Pyrimidine => [
            atoms.n1.or(atoms.c1),
            atoms.c2,
            atoms.n3,
            atoms.c4,
            atoms.c5,
            atoms.c6,
        ]
        .into_iter()
        .flatten()
        .collect(),
    }
}

fn nucleotide_base_normal(
    atoms: NucleotideAtoms,
    kind: Option<NucleotideBaseKind>,
) -> Option<Vec3> {
    let kind = kind.unwrap_or_else(|| {
        nucleotide_base_kind("", atoms).unwrap_or(NucleotideBaseKind::Pyrimidine)
    });
    let n1 = match kind {
        NucleotideBaseKind::Purine => atoms.n1,
        NucleotideBaseKind::Pyrimidine => atoms.n1.or(atoms.c1),
    }?;
    let c4 = atoms.c4?;
    let c5 = atoms.c5.or_else(|| {
        (kind == NucleotideBaseKind::Purine)
            .then_some(atoms.n5)
            .flatten()
    })?;
    let normal = (c4 - n1).cross(c5 - n1).normalized();
    (normal.length() > 0.000_001).then_some(normal)
}

#[derive(Clone, Debug)]
struct BackboneSegment {
    chain: String,
    start: i32,
    start_insertion_code: String,
    end: i32,
    end_insertion_code: String,
    points: Vec<Vec3>,
}

fn uncovered_backbone_segments(
    backbone: &[(String, i32, String, Vec3)],
    covered: &[(String, i32, String)],
) -> Vec<BackboneSegment> {
    let mut segments = Vec::<BackboneSegment>::new();
    let mut current = Vec::<(String, i32, String, Vec3)>::new();
    let mut previous: Option<(String, i32, String, Vec3)> = None;
    let keep_singletons = backbone.len() == 1;

    for (chain, seq, insertion_code, position) in backbone {
        let adjacent_to_previous =
            previous
                .as_ref()
                .is_some_and(|(last_chain, last_seq, _, last_position)| {
                    trace_points_connect(
                        last_chain,
                        *last_seq,
                        *last_position,
                        chain,
                        *seq,
                        *position,
                    )
                });
        let is_covered = covered
            .iter()
            .any(|(c, s, i)| c == chain && s == seq && i == insertion_code);

        if !adjacent_to_previous && !current.is_empty() {
            if current.len() >= 2 || keep_singletons {
                segments.push(backbone_segment(std::mem::take(&mut current)));
            } else {
                current.clear();
            }
        }

        if is_covered {
            if !current.is_empty() {
                current.push((chain.clone(), *seq, insertion_code.clone(), *position));
                if current.len() >= 2 || keep_singletons {
                    segments.push(backbone_segment(std::mem::take(&mut current)));
                } else {
                    current.clear();
                }
            }
            previous = Some((chain.clone(), *seq, insertion_code.clone(), *position));
            continue;
        }

        if current.is_empty() {
            if let Some((last_chain, last_seq, last_insertion_code, last_position)) = &previous {
                if trace_points_connect(
                    last_chain,
                    *last_seq,
                    *last_position,
                    chain,
                    *seq,
                    *position,
                ) {
                    current.push((
                        last_chain.clone(),
                        *last_seq,
                        last_insertion_code.clone(),
                        *last_position,
                    ));
                }
            }
        }
        current.push((chain.clone(), *seq, insertion_code.clone(), *position));
        previous = Some((chain.clone(), *seq, insertion_code.clone(), *position));
    }

    if current.len() >= 2 || (keep_singletons && !current.is_empty()) {
        segments.push(backbone_segment(current));
    }
    segments
}

fn trace_points_connect(
    last_chain: &str,
    last_seq: i32,
    last_position: Vec3,
    chain: &str,
    seq: i32,
    position: Vec3,
) -> bool {
    if last_chain != chain || seq <= last_seq {
        return false;
    }

    let ca_distance = last_position.distance(position);
    if seq == last_seq + 1 {
        return ca_distance <= 5.0;
    }

    let missing_residues = seq - last_seq - 1;
    missing_residues <= 24 && ca_distance <= 18.0
}

fn backbone_segment(points: Vec<(String, i32, String, Vec3)>) -> BackboneSegment {
    let chain = points.first().map(|p| p.0.clone()).unwrap_or_default();
    let start = points.first().map(|p| p.1).unwrap_or_default();
    let start_insertion_code = points.first().map(|p| p.2.clone()).unwrap_or_default();
    let end = points.last().map(|p| p.1).unwrap_or(start);
    let end_insertion_code = points.last().map(|p| p.2.clone()).unwrap_or_default();
    let points = points.into_iter().map(|(_, _, _, p)| p).collect();
    BackboneSegment {
        chain,
        start,
        start_insertion_code,
        end,
        end_insertion_code,
        points,
    }
}

fn add_nucleotide_ring(
    mesh: &mut Mesh,
    _center: Vec3,
    _normal: Vec3,
    radius: f32,
    base: Option<NucleotideRingBase>,
    detail: usize,
    radial_segments: usize,
    cylinder_cache: &mut CylinderPrimitiveCache,
) {
    if let Some(base) = base {
        add_nucleotide_named_atom_ring(mesh, base, radius, detail, radial_segments, cylinder_cache);
    }
}

#[allow(clippy::too_many_arguments)]
fn add_nucleotide_block(
    mesh: &mut Mesh,
    geometry: NucleotideBlockGeometry,
    radius: f32,
    width: f32,
    depth: f32,
    radial_segments: usize,
    cylinder_cache: &mut CylinderPrimitiveCache,
) {
    add_molstar_cylinder_caps_cached(
        mesh,
        geometry.anchor,
        geometry.trace,
        radius,
        radial_segments,
        false,
        true,
        cylinder_cache,
    );
    if let Some(block) = geometry.block {
        add_molstar_nucleotide_block_box(mesh, block, width, depth);
    }
}

fn add_molstar_nucleotide_block_box(
    mesh: &mut Mesh,
    block: NucleotideBlockBox,
    width: f32,
    depth: f32,
) {
    let v12 = (block.p2 - block.p1).normalized();
    let v34 = (block.p4 - block.p3).normalized();
    let up = v12.cross(v34).normalized();
    let center = block.p1 + v12 * (block.height / 2.0 - 0.2);
    let (x_axis, y_axis, z_axis) = molstar_target_to_axes(block.p1, block.p2, up);
    add_molstar_box_primitive(
        mesh,
        center,
        x_axis * width,
        y_axis * depth,
        z_axis * block.height,
    );
}

fn molstar_target_to_axes(eye: Vec3, target: Vec3, up: Vec3) -> (Vec3, Vec3, Vec3) {
    let z = (eye - target).normalized();
    let x = up.cross(z).normalized();
    let y = z.cross(x);
    (x, y, z)
}

fn add_molstar_box_primitive(
    mesh: &mut Mesh,
    center: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
    z_axis: Vec3,
) {
    geometry::add_molstar_box_primitive(
        mesh,
        MolstarPrimitiveTransform::from_axes(center, x_axis, y_axis, z_axis),
        false,
    );
}

fn add_direction_wedge(mesh: &mut Mesh, center: Vec3, tangent: Vec3, up: Vec3, size: f32) {
    let (axis, side, up) = oriented_basis(tangent, up);
    add_molstar_wedge_primitive(
        mesh,
        center,
        axis * (size * 6.0),
        up * (size * -4.0),
        side * (size * -4.0),
    );
}

const MOLSTAR_CARBOHYDRATE_SYMBOL_DETAIL: usize = 0;
const MOLSTAR_CARBOHYDRATE_SYMBOL_SIZE_FACTOR: f32 = 1.75;
const MOLSTAR_CARBOHYDRATE_SYMBOL_SIDE_FACTOR: f32 = 2.0 * 0.806;

fn add_carbohydrate_symbol(
    mesh: &mut Mesh,
    center: Vec3,
    normal: Vec3,
    direction: Vec3,
    shape: SaccharideShape,
    part: CarbohydrateSymbolPart,
) {
    let radius = MOLSTAR_CARBOHYDRATE_SYMBOL_SIZE_FACTOR;
    let side = MOLSTAR_CARBOHYDRATE_SYMBOL_SIZE_FACTOR * MOLSTAR_CARBOHYDRATE_SYMBOL_SIDE_FACTOR;
    let transform = carbohydrate_symbol_transform(center, normal, direction);
    let primary = part != CarbohydrateSymbolPart::Secondary;
    let secondary = part == CarbohydrateSymbolPart::Secondary;

    match shape {
        SaccharideShape::FilledSphere => {
            if primary {
                add_sphere(mesh, center, radius, MOLSTAR_CARBOHYDRATE_SYMBOL_DETAIL);
            }
        }
        SaccharideShape::FilledCube => {
            if primary {
                geometry::add_molstar_box_primitive(mesh, transform.scale_uniformly(side), false);
            }
        }
        SaccharideShape::CrossedCube => {
            let mut transform = transform.scale_uniformly(side);
            if secondary {
                transform = transform.mul_local(MolstarLocalTransform::rot_z90x180());
            }
            geometry::add_molstar_box_primitive(mesh, transform, true);
        }
        SaccharideShape::FilledCone => {
            if primary {
                geometry::add_molstar_pyramid_primitive(
                    mesh,
                    transform.scale_uniformly(side * 1.2),
                    8,
                    true,
                );
            }
        }
        SaccharideShape::DevidedCone => {
            let mut transform = transform.scale_uniformly(side * 1.2);
            if secondary {
                transform = transform.mul_local(MolstarLocalTransform::rot_z90());
            }
            geometry::add_molstar_perforated_octagonal_pyramid_primitive(mesh, transform);
        }
        SaccharideShape::FlatBox => {
            if primary {
                geometry::add_molstar_box_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zy90())
                        .scale(Vec3::new(side, side, side / 2.0)),
                    false,
                );
            }
        }
        SaccharideShape::FilledStar => {
            if primary {
                geometry::add_molstar_star_primitive(
                    mesh,
                    transform
                        .scale_uniformly(side)
                        .mul_local(MolstarLocalTransform::rot_zy90()),
                );
            }
        }
        SaccharideShape::FilledDiamond => {
            if primary {
                geometry::add_molstar_octahedron_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zy90())
                        .scale(Vec3::new(side * 1.4, side * 1.4, side * 1.4)),
                    false,
                );
            }
        }
        SaccharideShape::DividedDiamond => {
            let mut transform = transform
                .mul_local(MolstarLocalTransform::rot_zy90())
                .scale(Vec3::new(side * 1.4, side * 1.4, side * 1.4));
            if secondary {
                transform = transform.mul_local(MolstarLocalTransform::rot_y90());
            }
            geometry::add_molstar_octahedron_primitive(mesh, transform, true);
        }
        SaccharideShape::FlatDiamond => {
            if primary {
                geometry::add_molstar_prism_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zy90())
                        .scale(Vec3::new(side, side / 2.0, side / 2.0)),
                    4,
                    false,
                );
            }
        }
        SaccharideShape::DiamondPrism => {
            if primary {
                geometry::add_molstar_prism_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zy90())
                        .scale(Vec3::new(side, side, side / 2.0)),
                    4,
                    false,
                );
            }
        }
        SaccharideShape::PentagonalPrism | SaccharideShape::Pentagon => {
            if primary {
                geometry::add_molstar_prism_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zy90())
                        .scale(Vec3::new(side, side, side / 2.0)),
                    5,
                    false,
                );
            }
        }
        SaccharideShape::HexagonalPrism => {
            if primary {
                geometry::add_molstar_prism_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zy90())
                        .scale(Vec3::new(side, side, side / 2.0)),
                    6,
                    false,
                );
            }
        }
        SaccharideShape::HeptagonalPrism => {
            if primary {
                geometry::add_molstar_prism_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zy90())
                        .scale(Vec3::new(side, side, side / 2.0)),
                    7,
                    false,
                );
            }
        }
        SaccharideShape::FlatHexagon => {
            if primary {
                geometry::add_molstar_prism_primitive(
                    mesh,
                    transform
                        .mul_local(MolstarLocalTransform::rot_zyz90())
                        .scale(Vec3::new(side / 1.5, side, side / 2.0)),
                    6,
                    true,
                );
            }
        }
    }
}

fn carbohydrate_symbol_transform(
    center: Vec3,
    normal: Vec3,
    direction: Vec3,
) -> MolstarPrimitiveTransform {
    let target = if direction.length() > 0.000_001 {
        center + direction
    } else {
        center + Vec3::new(0.0, 0.0, 1.0)
    };
    MolstarPrimitiveTransform::from_target_to(center, target, normal)
}

fn carbohydrate_symbol_face_count(shape: SaccharideShape, part: CarbohydrateSymbolPart) -> usize {
    if part == CarbohydrateSymbolPart::Secondary && !carbohydrate_symbol_has_secondary_part(shape) {
        return 0;
    }
    match shape {
        SaccharideShape::FilledSphere => {
            molstar_sphere_triangle_count(MOLSTAR_CARBOHYDRATE_SYMBOL_DETAIL)
        }
        SaccharideShape::FilledCube | SaccharideShape::FlatBox => 12,
        SaccharideShape::CrossedCube => 6,
        SaccharideShape::FilledCone => 16,
        SaccharideShape::DevidedCone => 8,
        SaccharideShape::FilledStar => 20,
        SaccharideShape::FilledDiamond => 8,
        SaccharideShape::DividedDiamond => 4,
        SaccharideShape::FlatDiamond | SaccharideShape::DiamondPrism => 12,
        SaccharideShape::PentagonalPrism | SaccharideShape::Pentagon => 20,
        SaccharideShape::HexagonalPrism | SaccharideShape::FlatHexagon => 24,
        SaccharideShape::HeptagonalPrism => 28,
    }
}

fn carbohydrate_symbol_has_secondary_part(shape: SaccharideShape) -> bool {
    matches!(
        shape,
        SaccharideShape::CrossedCube
            | SaccharideShape::DevidedCone
            | SaccharideShape::DividedDiamond
    )
}

fn add_molstar_wedge_primitive(
    mesh: &mut Mesh,
    center: Vec3,
    x_axis: Vec3,
    y_axis: Vec3,
    z_axis: Vec3,
) {
    geometry::add_molstar_wedge_primitive(
        mesh,
        MolstarPrimitiveTransform::from_axes(center, x_axis, y_axis, z_axis),
    );
}

fn oriented_basis(axis: Vec3, up_hint: Vec3) -> (Vec3, Vec3, Vec3) {
    let axis = if axis.length() > 0.000_001 {
        axis.normalized()
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    };
    let mut up = up_hint - axis * axis.dot(up_hint);
    if up.length() <= 0.000_001 {
        up = fallback_side(axis, None);
    }
    up = up.normalized();
    let side = axis.cross(up).normalized();
    (axis, side, up)
}

fn add_nucleotide_named_atom_ring(
    mesh: &mut Mesh,
    base: NucleotideRingBase,
    radius: f32,
    detail: usize,
    radial_segments: usize,
    cylinder_cache: &mut CylinderPrimitiveCache,
) {
    let radial_segments = radial_segments.max(3);
    match base {
        NucleotideRingBase::PurineConnector { trace, n9 } => {
            add_open_cylinder_cached(mesh, n9, trace, radius, radial_segments, cylinder_cache);
            add_sphere(mesh, n9, radius, detail);
        }
        NucleotideRingBase::Purine {
            trace,
            n1,
            c2,
            n3,
            c4,
            c5,
            c6,
            n7,
            c8,
            n9,
        } => {
            add_open_cylinder_cached(mesh, n9, trace, radius, radial_segments, cylinder_cache);
            add_sphere(mesh, n9, radius, detail);
            add_molstar_nucleotide_ring_5_6_faces(
                mesh,
                radius,
                [n1, c2, n3, c4, c5, c6, n7, c8, n9],
            );
        }
        NucleotideRingBase::PyrimidineConnector { trace, n1 } => {
            add_open_cylinder_cached(mesh, n1, trace, radius, radial_segments, cylinder_cache);
            add_sphere(mesh, n1, radius, detail);
        }
        NucleotideRingBase::Pyrimidine {
            trace,
            n1,
            c2,
            n3,
            c4,
            c5,
            c6,
        } => {
            add_open_cylinder_cached(mesh, n1, trace, radius, radial_segments, cylinder_cache);
            add_sphere(mesh, n1, radius, detail);
            add_molstar_nucleotide_ring_6_faces(mesh, radius, [n1, c2, n3, c4, c5, c6]);
        }
    }
}

const MOLSTAR_RING_5_6_STRIP_INDICES: [usize; 20] = [
    0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 14, 15, 12, 13, 8, 9, 10, 11, 0, 1,
];
const MOLSTAR_RING_5_6_TOP_FAN_INDICES: [usize; 9] = [8, 12, 14, 16, 6, 4, 2, 0, 10];
const MOLSTAR_RING_5_6_BOTTOM_FAN_INDICES: [usize; 9] = [9, 11, 1, 3, 5, 7, 17, 15, 13];
const MOLSTAR_RING_6_STRIP_INDICES: [usize; 14] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1];
const MOLSTAR_RING_6_TOP_FAN_INDICES: [usize; 6] = [0, 10, 8, 6, 4, 2];
const MOLSTAR_RING_6_BOTTOM_FAN_INDICES: [usize; 6] = [1, 3, 5, 7, 9, 11];

fn nucleotide_ring_face_count(
    base: Option<NucleotideRingBase>,
    detail: usize,
    radial_segments: usize,
) -> usize {
    let radial = radial_segments.max(3);
    match base {
        Some(NucleotideRingBase::PurineConnector { .. })
        | Some(NucleotideRingBase::PyrimidineConnector { .. }) => {
            radial * 2 + molstar_sphere_triangle_count(detail)
        }
        Some(NucleotideRingBase::Purine { .. }) => {
            radial * 2
                + molstar_sphere_triangle_count(detail)
                + molstar_nucleotide_ring_5_6_face_count()
        }
        Some(NucleotideRingBase::Pyrimidine { .. }) => {
            radial * 2
                + molstar_sphere_triangle_count(detail)
                + molstar_nucleotide_ring_6_face_count()
        }
        None => 0,
    }
}

fn molstar_nucleotide_ring_5_6_face_count() -> usize {
    molstar_triangle_strip_face_count(&MOLSTAR_RING_5_6_STRIP_INDICES)
        + molstar_triangle_fan_face_count(&MOLSTAR_RING_5_6_TOP_FAN_INDICES)
        + molstar_triangle_fan_face_count(&MOLSTAR_RING_5_6_BOTTOM_FAN_INDICES)
}

fn molstar_nucleotide_ring_6_face_count() -> usize {
    molstar_triangle_strip_face_count(&MOLSTAR_RING_6_STRIP_INDICES)
        + molstar_triangle_fan_face_count(&MOLSTAR_RING_6_TOP_FAN_INDICES)
        + molstar_triangle_fan_face_count(&MOLSTAR_RING_6_BOTTOM_FAN_INDICES)
}

fn molstar_triangle_strip_face_count(indices: &[usize]) -> usize {
    if indices.len() < 4 {
        return 0;
    }
    indices[2..].chunks_exact(2).count() * 2
}

fn molstar_triangle_fan_face_count(indices: &[usize]) -> usize {
    indices.len().saturating_sub(2)
}

fn add_molstar_nucleotide_ring_5_6_faces(mesh: &mut Mesh, thickness: f32, points: [Vec3; 9]) {
    let positions = molstar_shifted_ring_positions(
        (points[3] - points[0])
            .cross(points[4] - points[0])
            .normalized()
            * thickness,
        &points,
    );
    add_molstar_triangle_strip(mesh, &positions, &MOLSTAR_RING_5_6_STRIP_INDICES);
    add_molstar_triangle_fan(mesh, &positions, &MOLSTAR_RING_5_6_TOP_FAN_INDICES);
    add_molstar_triangle_fan(mesh, &positions, &MOLSTAR_RING_5_6_BOTTOM_FAN_INDICES);
}

fn add_molstar_nucleotide_ring_6_faces(mesh: &mut Mesh, thickness: f32, points: [Vec3; 6]) {
    let positions = molstar_shifted_ring_positions(
        (points[3] - points[0])
            .cross(points[4] - points[0])
            .normalized()
            * thickness,
        &points,
    );
    add_molstar_triangle_strip(mesh, &positions, &MOLSTAR_RING_6_STRIP_INDICES);
    add_molstar_triangle_fan(mesh, &positions, &MOLSTAR_RING_6_TOP_FAN_INDICES);
    add_molstar_triangle_fan(mesh, &positions, &MOLSTAR_RING_6_BOTTOM_FAN_INDICES);
}

fn molstar_shifted_ring_positions(shift: Vec3, points: &[Vec3]) -> Vec<Vec3> {
    let mut positions = Vec::with_capacity(points.len() * 2);
    for &point in points {
        positions.push(point + shift);
        positions.push(point - shift);
    }
    positions
}

fn add_molstar_triangle_strip(mesh: &mut Mesh, positions: &[Vec3], indices: &[usize]) {
    if indices.len() < 4 {
        return;
    }
    let mut c = positions[indices[0]];
    let mut d = positions[indices[1]];
    for pair in indices[2..].chunks_exact(2) {
        let a = c;
        let b = d;
        c = positions[pair[0]];
        d = positions[pair[1]];
        add_molstar_triangle(mesh, a, b, c);
        add_molstar_triangle(mesh, b, d, c);
    }
}

fn add_molstar_triangle_fan(mesh: &mut Mesh, positions: &[Vec3], indices: &[usize]) {
    if indices.len() < 3 {
        return;
    }
    let a = positions[indices[0]];
    for i in 2..indices.len() {
        let b = positions[indices[i - 1]];
        let c = positions[indices[i]];
        add_molstar_triangle(mesh, a, c, b);
    }
}

fn add_molstar_triangle(mesh: &mut Mesh, a: Vec3, b: Vec3, c: Vec3) {
    let normal = (b - a).cross(c - a).normalized();
    let base = mesh.vertices.len();
    mesh.vertices.push(a);
    mesh.vertices.push(b);
    mesh.vertices.push(c);
    mesh.normals.push(normal);
    mesh.normals.push(normal);
    mesh.normals.push(normal);
    mesh.faces.push(Face {
        a: base,
        b: base + 1,
        c: base + 2,
    });
}

#[cfg(test)]
mod tests {
    use std::f32::consts::PI;

    use super::*;
    use crate::model::{StructureUnit, UnitKind, UnitProps, UnitSymmetryGroup, UnitTraits};

    fn symmetry_group(unit_count: usize) -> UnitSymmetryGroup {
        UnitSymmetryGroup {
            kind: UnitKind::Atomic,
            model_id: 0,
            invariant_id: 0,
            elements: Vec::new(),
            unit_ids: (0..unit_count).collect(),
            operator_names: Vec::new(),
            operator_instance_ids: Vec::new(),
            unit_index_map: Vec::new(),
            hash_code: 0,
            transform_hash: 0,
        }
    }

    fn atomic_unit_with_props(props: UnitProps) -> StructureUnit {
        StructureUnit {
            id: 0,
            invariant_id: 0,
            chain_group_id: 0,
            kind: UnitKind::Atomic,
            traits: UnitTraits::NONE,
            model_index: 0,
            chain_index: 0,
            chain_indices: Vec::new(),
            elements: Vec::new(),
            atom_indices: Vec::new(),
            residue_indices: Vec::new(),
            residue_index_by_element: Vec::new(),
            chain_index_by_element: Vec::new(),
            props,
            operator: Default::default(),
        }
    }

    fn assert_render_object_mesh_estimate(
        object: RenderObject,
        options: &MeshOptions,
        cylinder_radial_segments: usize,
    ) {
        let estimate = object.mesh_estimate(options, cylinder_radial_segments);
        let mut mesh = Mesh::default();
        let mut cylinder_cache = CylinderPrimitiveCache::default();
        let mut curve_scratch = CurveSegmentScratch::default();
        append_render_object_to_mesh(
            &mut mesh,
            &object,
            options,
            cylinder_radial_segments,
            &mut cylinder_cache,
            &mut curve_scratch,
            None,
        );
        assert_eq!(
            estimate.vertices,
            mesh.vertices.len(),
            "vertex estimate mismatch for {object:?}"
        );
        assert_eq!(
            estimate.faces,
            mesh.faces.len(),
            "face estimate mismatch for {object:?}"
        );
        assert_eq!(
            render_object_mesh_stats(&object, estimate.vertices, estimate.faces),
            render_object_mesh_stats(&object, mesh.vertices.len(), mesh.faces.len()),
            "value-cell estimate mismatch for {object:?}"
        );
    }

    #[test]
    fn render_object_mesh_estimates_match_every_geometry_builder() {
        let options = MeshOptions {
            sphere_detail: 1,
            linear_segments: 6,
            radial_segments: 12,
            ..MeshOptions::default()
        };
        let cylinder_radial_segments = 24;
        let points = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.3, 0.1),
            Vec3::new(2.0, -0.2, 0.2),
        ];
        let controls = CurveSegmentControls {
            sec_struc_first: false,
            sec_struc_last: false,
            p0: DVec3::new(-2.0, 0.0, 0.0),
            p1: DVec3::new(-1.0, 0.2, 0.0),
            p2: DVec3::new(0.0, 0.0, 0.0),
            p3: DVec3::new(1.0, 0.2, 0.0),
            p4: DVec3::new(2.0, 0.0, 0.0),
            d12: DVec3::new(0.0, 1.0, 0.0),
            d23: DVec3::new(0.0, 1.0, 0.0),
        };
        let block = NucleotideBlockBox {
            p1: Vec3::new(0.0, 0.0, 0.0),
            p2: Vec3::new(1.0, 0.0, 0.0),
            p3: Vec3::new(0.0, 1.0, 0.0),
            p4: Vec3::new(1.0, 1.0, 0.0),
            height: 1.2,
        };
        let mut objects = vec![
            RenderObject::Sphere {
                center: Vec3::default(),
                radius: 1.0,
            },
            RenderObject::Cylinder {
                start: Vec3::default(),
                end: Vec3::new(2.0, 0.0, 0.0),
                radius: 0.2,
            },
            RenderObject::LinkCylinder {
                start: Vec3::default(),
                end: Vec3::new(2.0, 0.0, 0.0),
                radius: 0.2,
            },
            RenderObject::LinkCylinderWithSegments {
                start: Vec3::default(),
                end: Vec3::new(2.0, 0.0, 0.0),
                radius: 0.2,
                radial_segments: 4,
            },
            RenderObject::Tube {
                points: points.clone(),
                radius: 0.3,
            },
            RenderObject::DashedTube {
                points: points.clone(),
                radius: 0.2,
            },
            RenderObject::FixedCountDashedCylinder {
                start: Vec3::default(),
                end: Vec3::new(4.0, 0.0, 0.0),
                radius: 0.2,
                length_scale: 1.0,
                segment_count: 5,
            },
            RenderObject::Ribbon {
                points: points.clone(),
                width: 0.8,
                thickness: 0.2,
            },
            RenderObject::Sheet {
                points: points.clone(),
                width: 0.8,
                thickness: 0.2,
                arrow_height: 0.6,
                start_cap: false,
                end_cap: true,
            },
            RenderObject::OrientedRibbon {
                centers: points.clone(),
                normals: vec![Vec3::new(0.0, 1.0, 0.0); points.len()],
                width: 0.8,
                thickness: 0.2,
                profile: PolymerProfile::Rounded,
                start_cap: true,
                end_cap: true,
                round_cap: true,
            },
            RenderObject::PolymerTraceSegment {
                controls: controls.clone(),
                widths: [0.8; 3],
                heights: [0.2; 3],
                tension: 0.5,
                shift: 0.5,
                overhang_width: 0.8,
                kind: PolymerTraceSegmentKind::Ribbon {
                    arrow_height: 0.0,
                    swap_width_height: false,
                },
                start_cap: false,
                end_cap: false,
                initial: true,
                final_residue: false,
                swap_normal_binormal: false,
            },
            RenderObject::PolymerTraceSegment {
                controls: controls.clone(),
                widths: [0.8; 3],
                heights: [0.2; 3],
                tension: 0.5,
                shift: 0.5,
                overhang_width: 0.8,
                kind: PolymerTraceSegmentKind::Tube {
                    profile: PolymerProfile::Rounded,
                    round_cap: true,
                },
                start_cap: true,
                end_cap: true,
                initial: false,
                final_residue: false,
                swap_normal_binormal: false,
            },
            RenderObject::PolymerTraceSegment {
                controls,
                widths: [0.8; 3],
                heights: [0.2; 3],
                tension: 0.5,
                shift: 0.5,
                overhang_width: 0.8,
                kind: PolymerTraceSegmentKind::Sheet { arrow_height: 0.6 },
                start_cap: false,
                end_cap: true,
                initial: false,
                final_residue: true,
                swap_normal_binormal: false,
            },
            RenderObject::NucleotideRing {
                center: Vec3::default(),
                normal: Vec3::new(0.0, 0.0, 1.0),
                radius: 0.2,
                base: Some(NucleotideRingBase::Pyrimidine {
                    trace: Vec3::new(-1.0, 0.0, 0.0),
                    n1: Vec3::new(0.0, 0.0, 0.0),
                    c2: Vec3::new(0.5, 0.8, 0.0),
                    n3: Vec3::new(1.5, 0.8, 0.0),
                    c4: Vec3::new(2.0, 0.0, 0.0),
                    c5: Vec3::new(1.5, -0.8, 0.0),
                    c6: Vec3::new(0.5, -0.8, 0.0),
                }),
                detail: 1,
                radial_segments: 12,
            },
            RenderObject::NucleotideBlock {
                geometry: NucleotideBlockGeometry {
                    trace: Vec3::new(-1.0, 0.0, 0.0),
                    anchor: Vec3::default(),
                    block: Some(block),
                },
                radius: 0.2,
                width: 1.0,
                depth: 0.4,
                radial_segments: 12,
            },
            RenderObject::DirectionWedge {
                center: Vec3::default(),
                tangent: Vec3::new(1.0, 0.0, 0.0),
                up: Vec3::new(0.0, 1.0, 0.0),
                size: 0.5,
            },
            RenderObject::Ellipsoid {
                center: Vec3::default(),
                axes: [
                    Vec3::new(1.0, 0.0, 0.0),
                    Vec3::new(0.0, 2.0, 0.0),
                    Vec3::new(0.0, 0.0, 3.0),
                ],
            },
        ];

        for shape in [
            SaccharideShape::FilledSphere,
            SaccharideShape::FilledCube,
            SaccharideShape::CrossedCube,
            SaccharideShape::FilledCone,
            SaccharideShape::DevidedCone,
            SaccharideShape::FlatBox,
            SaccharideShape::FilledStar,
            SaccharideShape::FilledDiamond,
            SaccharideShape::DividedDiamond,
            SaccharideShape::FlatDiamond,
            SaccharideShape::DiamondPrism,
            SaccharideShape::PentagonalPrism,
            SaccharideShape::Pentagon,
            SaccharideShape::HexagonalPrism,
            SaccharideShape::HeptagonalPrism,
            SaccharideShape::FlatHexagon,
        ] {
            objects.push(RenderObject::CarbohydrateSymbol {
                center: Vec3::default(),
                normal: Vec3::new(0.0, 1.0, 0.0),
                direction: Vec3::new(1.0, 0.0, 0.0),
                shape,
                part: CarbohydrateSymbolPart::Primary,
            });
            if carbohydrate_symbol_has_secondary_part(shape) {
                objects.push(RenderObject::CarbohydrateSymbol {
                    center: Vec3::default(),
                    normal: Vec3::new(0.0, 1.0, 0.0),
                    direction: Vec3::new(1.0, 0.0, 0.0),
                    shape,
                    part: CarbohydrateSymbolPart::Secondary,
                });
            }
        }

        let flatten_cylinder_radial_segments = molstar_export_cylinder_radial_segments(
            objects
                .iter()
                .map(render_object_export_cylinder_count)
                .sum(),
        );
        let total = render_objects_mesh_estimate(
            objects.iter(),
            &options,
            flatten_cylinder_radial_segments,
        );
        let mesh = flatten_render_objects(&objects, &Molecule::default(), &options);
        assert_eq!(mesh.vertices.len(), total.vertices);
        assert_eq!(mesh.normals.len(), total.vertices);
        assert_eq!(mesh.faces.len(), total.faces);
        assert!(mesh.vertices.capacity() >= total.vertices);
        assert!(mesh.normals.capacity() >= total.vertices);
        assert!(mesh.faces.capacity() >= total.faces);
        assert!(mesh.vertex_groups.capacity() >= total.vertices);
        assert!(mesh.face_groups.capacity() >= total.faces);

        for object in objects {
            assert_render_object_mesh_estimate(object, &options, cylinder_radial_segments);
        }
    }

    #[test]
    fn polymer_trace_average4_uses_molstar_nested_vec3_add_order() {
        let actual = vec3_average4_f64(
            DVec3::from_vec3(Vec3::new(1.0e20, 0.0, 0.0)),
            DVec3::from_vec3(Vec3::new(-1.0e20, 0.0, 0.0)),
            DVec3::from_vec3(Vec3::new(1.0, 0.0, 0.0)),
            DVec3::from_vec3(Vec3::new(1.0, 0.0, 0.0)),
        );

        assert_eq!(
            actual.x, 0.0,
            "Mol* evaluates setControlPoint/setDirection as a + (b + (c + d)) before scaling"
        );
    }

    #[test]
    fn molstar_spacefill_visual_selection_switches_at_symmetry_group_threshold() {
        let mut structure = AtomicStructure {
            element_count: 1,
            symmetry_groups: vec![symmetry_group(1); 5_001],
            ..AtomicStructure::default()
        };
        let options = MeshOptions {
            representation: Representation::Spacefill,
            ..MeshOptions::default()
        };

        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["structure-element-sphere".to_string()]
        );

        structure.symmetry_groups.truncate(5_000);
        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["element-sphere".to_string()]
        );
    }

    #[test]
    fn molstar_ball_and_stick_visual_selection_matches_huge_and_symmetry_branches() {
        let options = MeshOptions {
            representation: Representation::BallAndStick,
            ..MeshOptions::default()
        };
        let mut structure = AtomicStructure {
            polymer_residue_count: 29_999,
            symmetry_groups: vec![symmetry_group(1); 5_001],
            ..AtomicStructure::default()
        };

        assert_eq!(
            selected_visuals(&structure, &options),
            vec![
                "structure-element-sphere".to_string(),
                "structure-intra-bond".to_string()
            ]
        );

        structure.polymer_residue_count = 30_000;
        structure.symmetry_groups = vec![symmetry_group(11)];
        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["element-sphere".to_string(), "intra-bond".to_string()]
        );

        structure.polymer_residue_count = 9;
        structure.symmetry_groups = vec![symmetry_group(1)];
        assert_eq!(
            selected_visuals(&structure, &options),
            vec![
                "element-sphere".to_string(),
                "intra-bond".to_string(),
                "inter-bond".to_string()
            ]
        );
    }

    #[test]
    fn molstar_cartoon_visual_selection_matches_get_cartoon_params() {
        let options = MeshOptions {
            representation: Representation::Cartoon,
            ..MeshOptions::default()
        };
        let molstar_options = MeshOptions {
            representation: Representation::Molstar,
            ..MeshOptions::default()
        };
        let mut structure = AtomicStructure {
            polymer_residue_count: 3,
            ..AtomicStructure::default()
        };

        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["polymer-trace".to_string()]
        );
        assert_eq!(
            selected_visuals(&structure, &molstar_options),
            vec!["polymer-trace".to_string()]
        );

        structure.units = vec![atomic_unit_with_props(UnitProps {
            nucleotide_elements: vec![0],
            ..UnitProps::default()
        })];
        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["polymer-trace".to_string(), "nucleotide-ring".to_string()]
        );

        structure.polymer_gap_count = 1;
        structure.units[0].props.gap_elements = vec![0, 2];
        assert_eq!(
            selected_visuals(&structure, &options),
            vec![
                "polymer-trace".to_string(),
                "nucleotide-ring".to_string(),
                "polymer-gap".to_string()
            ]
        );

        let selected = selected_visuals(&structure, &options);
        assert!(!selected.iter().any(|visual| matches!(
            visual.as_str(),
            "nucleotide-block"
                | "direction-wedge"
                | "nucleotide-atomic-ring-fill"
                | "nucleotide-atomic-bond"
                | "nucleotide-atomic-element"
        )));
    }

    #[test]
    fn molstar_ribbon_visual_selection_matches_putty_params() {
        let options = MeshOptions {
            representation: Representation::Ribbon,
            ..MeshOptions::default()
        };
        let mut structure = AtomicStructure {
            polymer_residue_count: 3,
            ..AtomicStructure::default()
        };

        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["polymer-tube".to_string()]
        );

        structure.units = vec![atomic_unit_with_props(UnitProps {
            nucleotide_elements: vec![0],
            ..UnitProps::default()
        })];
        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["polymer-tube".to_string()]
        );

        structure.polymer_gap_count = 1;
        structure.units[0].props.gap_elements = vec![0, 2];
        assert_eq!(
            selected_visuals(&structure, &options),
            vec!["polymer-tube".to_string(), "polymer-gap".to_string()]
        );
    }

    #[test]
    fn oriented_ribbon_value_cell_estimate_tracks_radial_dispatch() {
        for (radial_segments, profile) in [
            (2, PolymerProfile::Elliptical),
            (4, PolymerProfile::Elliptical),
            (12, PolymerProfile::Square),
        ] {
            let object = RenderObject::OrientedRibbon {
                centers: vec![Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
                normals: vec![Vec3::new(0.0, 0.0, 1.0); 2],
                width: 2.0,
                thickness: 0.5,
                profile,
                start_cap: true,
                end_cap: true,
                round_cap: false,
            };
            let options = MeshOptions {
                linear_segments: 1,
                radial_segments,
                ..MeshOptions::default()
            };
            let mesh = flatten_render_objects(
                std::slice::from_ref(&object),
                &Molecule::default(),
                &options,
            );

            assert_eq!(
                object.face_estimate(&options),
                mesh.faces.len(),
                "radial_segments={radial_segments} profile={profile:?}"
            );
        }
    }

    #[test]
    fn nucleotide_ring_faces_follow_molstar_strip_and_fan_indices() {
        // Derived from artifacts/molstar/src/mol-repr/structure/visual/nucleotide-ring-mesh.ts:
        // positionsRing6, stripIndicesRing6, fanIndicesTopRing6, fanIndicesBottomRing6.
        let mut mesh = Mesh::default();
        add_molstar_nucleotide_ring_6_faces(
            &mut mesh,
            0.2,
            [
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(1.0, 1.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                Vec3::new(-1.0, 1.0, 0.0),
                Vec3::new(-1.0, 0.0, 0.0),
            ],
        );

        assert_eq!(mesh.faces.len(), 20);
        assert_eq!(mesh.vertices.len(), 60);
        assert_eq!(mesh.normals.len(), 60);
        assert_eq!(
            (mesh.faces[0].a, mesh.faces[0].b, mesh.faces[0].c),
            (0, 1, 2)
        );
        assert_eq!(
            (mesh.faces[1].a, mesh.faces[1].b, mesh.faces[1].c),
            (3, 4, 5)
        );
        assert_vec3_close(mesh.vertices[0], Vec3::new(0.0, 0.0, 0.2));
        assert_vec3_close(mesh.vertices[1], Vec3::new(0.0, 0.0, -0.2));
        assert_vec3_close(mesh.vertices[2], Vec3::new(1.0, 0.0, 0.2));
        assert_vec3_close(mesh.vertices[3], Vec3::new(0.0, 0.0, -0.2));
        assert_vec3_close(mesh.vertices[4], Vec3::new(1.0, 0.0, -0.2));
        assert_vec3_close(mesh.vertices[5], Vec3::new(1.0, 0.0, 0.2));
        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, -1.0, 0.0));
        assert_vec3_close(mesh.normals[3], Vec3::new(0.0, -1.0, 0.0));
    }

    #[test]
    fn nucleotide_ring_mesh_uses_molstar_default_size_and_detail() {
        let mut mesh = Mesh::default();
        let mut cylinder_cache = CylinderPrimitiveCache::default();
        add_nucleotide_named_atom_ring(
            &mut mesh,
            NucleotideRingBase::Pyrimidine {
                trace: Vec3::new(0.0, -1.5, 0.0),
                n1: Vec3::new(0.0, 0.0, 0.0),
                c2: Vec3::new(1.0, 0.0, 0.0),
                n3: Vec3::new(1.0, 1.0, 0.0),
                c4: Vec3::new(0.0, 1.0, 0.0),
                c5: Vec3::new(-1.0, 1.0, 0.0),
                c6: Vec3::new(-1.0, 0.0, 0.0),
            },
            0.2,
            0,
            16,
            &mut cylinder_cache,
        );

        assert_eq!(mesh.faces.len(), 72);
        assert_eq!(
            mesh.faces.len(),
            16 * 2 + molstar_sphere_triangle_count(0) + molstar_nucleotide_ring_6_face_count()
        );
        assert_eq!(mesh.vertices.len(), 106);
        assert!(
            mesh.vertices
                .iter()
                .any(|vertex| vertex.distance(Vec3::new(0.0, 0.0, 0.2)) < 0.000_001),
            "ring face must use sizeFactor/thicknessFactor default thickness"
        );
    }

    #[test]
    fn nucleotide_ring_mesh_uses_resolved_detail_and_radial_segments() {
        let mut mesh = Mesh::default();
        let mut cylinder_cache = CylinderPrimitiveCache::default();
        let base = NucleotideRingBase::Pyrimidine {
            trace: Vec3::new(0.0, -1.5, 0.0),
            n1: Vec3::new(0.0, 0.0, 0.0),
            c2: Vec3::new(1.0, 0.0, 0.0),
            n3: Vec3::new(1.0, 1.0, 0.0),
            c4: Vec3::new(0.0, 1.0, 0.0),
            c5: Vec3::new(-1.0, 1.0, 0.0),
            c6: Vec3::new(-1.0, 0.0, 0.0),
        };
        add_nucleotide_named_atom_ring(&mut mesh, base, 0.2, 1, 12, &mut cylinder_cache);

        let expected_faces =
            12 * 2 + molstar_sphere_triangle_count(1) + molstar_nucleotide_ring_6_face_count();
        assert_eq!(mesh.faces.len(), expected_faces);
        assert_eq!(
            nucleotide_ring_face_count(Some(base), 1, 12),
            expected_faces
        );

        let mut connector = Mesh::default();
        let mut connector_cache = CylinderPrimitiveCache::default();
        add_nucleotide_named_atom_ring(
            &mut connector,
            NucleotideRingBase::PyrimidineConnector {
                trace: Vec3::new(0.0, -1.5, 0.0),
                n1: Vec3::new(0.0, 0.0, 0.0),
            },
            0.2,
            1,
            8,
            &mut connector_cache,
        );
        assert_eq!(
            connector.faces.len(),
            8 * 2 + molstar_sphere_triangle_count(1)
        );
        assert_eq!(
            nucleotide_ring_face_count(
                Some(NucleotideRingBase::PyrimidineConnector {
                    trace: Vec3::new(0.0, -1.5, 0.0),
                    n1: Vec3::new(0.0, 0.0, 0.0),
                }),
                1,
                8
            ),
            connector.faces.len()
        );

        let mut fallback = Mesh::default();
        let mut fallback_cache = CylinderPrimitiveCache::default();
        add_nucleotide_ring(
            &mut fallback,
            Vec3::default(),
            Vec3::new(0.0, 0.0, 1.0),
            0.2,
            None,
            1,
            8,
            &mut fallback_cache,
        );
        assert!(
            fallback.faces.is_empty(),
            "Mol* does not emit a generic nucleotide annulus when named base atoms are missing"
        );
        assert_eq!(nucleotide_ring_face_count(None, 1, 8), 0);
    }

    #[test]
    fn nucleotide_block_box_primitive_matches_molstar_box_order() {
        // Derived from artifacts/molstar/src/mol-geo/primitive/box.ts:
        // Box() uses polygon(4, true), four side quads, then bottom and top quads.
        let mut mesh = Mesh::default();
        add_molstar_box_primitive(
            &mut mesh,
            Vec3::default(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );

        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.normals.len(), 24);
        assert_eq!(mesh.faces.len(), 12);
        assert_faces_eq(
            &mesh.faces,
            &[
                [0, 1, 2],
                [2, 3, 0],
                [4, 5, 6],
                [6, 7, 4],
                [8, 9, 10],
                [10, 11, 8],
                [12, 13, 14],
                [14, 15, 12],
                [16, 17, 18],
                [18, 19, 16],
                [20, 21, 22],
                [22, 23, 20],
            ],
        );
        assert_vec3_close(mesh.vertices[0], Vec3::new(0.5, 0.5, -0.5));
        assert_vec3_close(mesh.vertices[1], Vec3::new(-0.5, 0.5, -0.5));
        assert_vec3_close(mesh.vertices[2], Vec3::new(-0.5, 0.5, 0.5));
        assert_vec3_close(mesh.vertices[16], Vec3::new(0.5, -0.5, -0.5));
        assert_vec3_close(mesh.vertices[20], Vec3::new(0.5, 0.5, 0.5));
        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, 1.0, 0.0));
        assert_vec3_close(mesh.normals[16], Vec3::new(0.0, 0.0, -1.0));
        assert_vec3_close(mesh.normals[20], Vec3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn nucleotide_block_mesh_uses_molstar_connector_and_box_transform() {
        // Derived from artifacts/molstar/src/mol-repr/structure/visual/nucleotide-block-mesh.ts:
        // addCylinder(p5, trace, bottomCap=true), then targetTo(p1, p2, vC)
        // scaled by width/depth/height and translated to p1 + v12 * (height / 2 - 0.2).
        let mut mesh = Mesh::default();
        let mut cylinder_cache = CylinderPrimitiveCache::default();
        add_nucleotide_block(
            &mut mesh,
            NucleotideBlockGeometry {
                trace: Vec3::new(0.0, 0.0, -1.0),
                anchor: Vec3::new(0.0, 0.0, 0.0),
                block: Some(NucleotideBlockBox {
                    p1: Vec3::new(0.0, 0.0, 0.0),
                    p2: Vec3::new(0.0, 0.0, 1.0),
                    p3: Vec3::new(1.0, 0.0, 0.0),
                    p4: Vec3::new(1.0, 1.0, 0.0),
                    height: 0.9,
                }),
            },
            0.2,
            0.9,
            0.4,
            16,
            &mut cylinder_cache,
        );

        let box_start = 67;
        assert_eq!(mesh.faces.len(), 60);
        assert_eq!(mesh.vertices.len(), box_start + 24);
        assert_eq!(mesh.normals.len(), box_start + 24);
        assert_vec3_close(mesh.vertices[box_start], Vec3::new(-0.2, -0.45, 0.7));
        assert_vec3_close(mesh.vertices[box_start + 1], Vec3::new(-0.2, 0.45, 0.7));
        assert_vec3_close(mesh.vertices[box_start + 2], Vec3::new(-0.2, 0.45, -0.2));
        assert_vec3_close(mesh.normals[box_start], Vec3::new(-2.5, 0.0, 0.0));
    }

    #[test]
    fn direction_wedge_primitive_matches_molstar_triangle_prism_order() {
        // Derived from artifacts/molstar/src/mol-geo/primitive/wedge.ts:
        // Wedge() is PrimitiveBuilder over polygon(3, false) with six side
        // triangles followed by bottom and top base triangles.
        let mut mesh = Mesh::default();
        add_molstar_wedge_primitive(
            &mut mesh,
            Vec3::default(),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        );

        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.normals.len(), 24);
        assert_eq!(mesh.faces.len(), 8);
        assert_faces_eq(
            &mesh.faces,
            &[
                [0, 1, 2],
                [3, 4, 5],
                [6, 7, 8],
                [9, 10, 11],
                [12, 13, 14],
                [15, 16, 17],
                [18, 19, 20],
                [21, 22, 23],
            ],
        );

        assert_vec3_close(mesh.vertices[0], Vec3::new(0.707_106_77, 0.0, -0.5));
        assert_vec3_close(
            mesh.vertices[1],
            Vec3::new(-0.353_553_35, 0.612_372_46, -0.5),
        );
        assert_vec3_close(
            mesh.vertices[2],
            Vec3::new(-0.353_553_35, 0.612_372_46, 0.5),
        );
        assert_vec3_close(
            mesh.vertices[18],
            Vec3::new(-0.353_553_47, -0.612_372_4, -0.5),
        );
        assert_vec3_close(
            mesh.vertices[19],
            Vec3::new(-0.353_553_35, 0.612_372_46, -0.5),
        );
        assert_vec3_close(mesh.vertices[20], Vec3::new(0.707_106_77, 0.0, -0.5));
        assert_vec3_close(mesh.vertices[21], Vec3::new(0.707_106_77, 0.0, 0.5));
        assert_vec3_close(
            mesh.vertices[22],
            Vec3::new(-0.353_553_35, 0.612_372_46, 0.5),
        );
        assert_vec3_close(
            mesh.vertices[23],
            Vec3::new(-0.353_553_47, -0.612_372_4, 0.5),
        );

        assert_vec3_close(mesh.normals[0], Vec3::new(0.5, 0.866_025_4, 0.0));
        assert_vec3_close(mesh.normals[18], Vec3::new(0.0, 0.0, -1.0));
        assert_vec3_close(mesh.normals[21], Vec3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn direction_wedge_mesh_applies_molstar_target_to_rot_scale_axes() {
        // Derived from artifacts/molstar/src/mol-repr/structure/visual/polymer-direction-wedge.ts:
        // targetTo(p3, p1, up), rotY90Z180, scale(height, width, depth),
        // then setTranslation(p2).
        let mut mesh = Mesh::default();
        let center = Vec3::new(1.0, 2.0, 3.0);
        let size = 0.5;
        add_direction_wedge(
            &mut mesh,
            center,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            size,
        );

        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.normals.len(), 24);
        assert_eq!(mesh.faces.len(), 8);

        let radius = std::f32::consts::FRAC_1_SQRT_2;
        let points = [
            Vec3::new(radius, 0.0, 0.0),
            Vec3::new(
                (2.0 * PI / 3.0).cos() * radius,
                (2.0 * PI / 3.0).sin() * radius,
                0.0,
            ),
            Vec3::new(
                (4.0 * PI / 3.0).cos() * radius,
                (4.0 * PI / 3.0).sin() * radius,
                0.0,
            ),
        ];
        let x_axis = Vec3::new(size * 6.0, 0.0, 0.0);
        let y_axis = Vec3::new(0.0, size * -4.0, 0.0);
        let z_axis = Vec3::new(0.0, 0.0, size * -4.0);
        let transform =
            |point: Vec3| center + x_axis * point.x + y_axis * point.y + z_axis * point.z;

        assert_vec3_close(
            mesh.vertices[0],
            transform(Vec3::new(points[0].x, points[0].y, -0.5)),
        );
        assert_vec3_close(
            mesh.vertices[1],
            transform(Vec3::new(points[1].x, points[1].y, -0.5)),
        );
        assert_vec3_close(
            mesh.vertices[2],
            transform(Vec3::new(points[1].x, points[1].y, 0.5)),
        );
        assert_vec3_close(
            mesh.vertices[18],
            transform(Vec3::new(points[2].x, points[2].y, -0.5)),
        );
        assert_vec3_close(
            mesh.vertices[19],
            transform(Vec3::new(points[1].x, points[1].y, -0.5)),
        );
        assert_vec3_close(
            mesh.vertices[20],
            transform(Vec3::new(points[0].x, points[0].y, -0.5)),
        );
        assert_vec3_close(
            mesh.vertices[21],
            transform(Vec3::new(points[0].x, points[0].y, 0.5)),
        );
        assert_vec3_close(
            mesh.vertices[22],
            transform(Vec3::new(points[1].x, points[1].y, 0.5)),
        );
        assert_vec3_close(
            mesh.vertices[23],
            transform(Vec3::new(points[2].x, points[2].y, 0.5)),
        );

        assert_vec3_close(mesh.normals[0], Vec3::new(0.166_666_67, -0.433_012_7, 0.0));
        assert_vec3_close(mesh.normals[18], Vec3::new(0.0, 0.0, 0.5));
        assert_vec3_close(mesh.normals[21], Vec3::new(0.0, 0.0, -0.5));
    }

    #[test]
    fn carbohydrate_secondary_symbol_parts_use_molstar_primitive_halves() {
        // Derived from artifacts/molstar/src/mol-repr/structure/visual/carbohydrate-symbol-mesh.ts:
        // crossed/divided symbols append one primitive per picking group, with
        // the second primitive receiving the Mol* local rotation.
        let center = Vec3::default();
        let normal = Vec3::new(0.0, 1.0, 0.0);
        let direction = Vec3::new(1.0, 0.0, 0.0);

        let mut primary = Mesh::default();
        add_carbohydrate_symbol(
            &mut primary,
            center,
            normal,
            direction,
            SaccharideShape::CrossedCube,
            CarbohydrateSymbolPart::Primary,
        );
        let mut secondary = Mesh::default();
        add_carbohydrate_symbol(
            &mut secondary,
            center,
            normal,
            direction,
            SaccharideShape::CrossedCube,
            CarbohydrateSymbolPart::Secondary,
        );
        assert_eq!(primary.faces.len(), 6);
        assert_eq!(secondary.faces.len(), 6);
        assert_eq!(primary.vertices.len(), 18);
        assert_eq!(secondary.vertices.len(), 18);
        assert!(
            primary.vertices[0].distance(secondary.vertices[0]) > 0.001,
            "secondary crossed-cube half must be locally rotated"
        );

        let mut cone = Mesh::default();
        add_carbohydrate_symbol(
            &mut cone,
            center,
            normal,
            direction,
            SaccharideShape::DevidedCone,
            CarbohydrateSymbolPart::Primary,
        );
        assert_eq!(cone.faces.len(), 8);
        assert_eq!(cone.vertices.len(), 24);
        assert_eq!(
            carbohydrate_symbol_face_count(
                SaccharideShape::DevidedCone,
                CarbohydrateSymbolPart::Primary
            ),
            8
        );
        assert_eq!(
            carbohydrate_symbol_face_count(
                SaccharideShape::DevidedCone,
                CarbohydrateSymbolPart::Secondary
            ),
            8
        );
    }

    #[test]
    fn meshbuilder_state_groups_follow_molstar_current_group_append_order() {
        // Mol* MeshBuilder appends currentGroup alongside each added vertex.
        // molfig stores the export metadata per face, so the equivalent is one
        // group entry per newly appended face, in the same face append order,
        // while vertex_groups preserve the Mol* group-buffer state layout.
        let options = MeshOptions {
            sphere_detail: 0,
            center: false,
            assembly: None,
            ..MeshOptions::default()
        };
        let objects = vec![
            RenderObject::Sphere {
                center: Vec3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
            RenderObject::Cylinder {
                start: Vec3::new(1.0, 0.0, 0.0),
                end: Vec3::new(1.0, 0.0, 0.0),
                radius: 0.25,
            },
            RenderObject::Cylinder {
                start: Vec3::new(0.0, 0.0, 0.0),
                end: Vec3::new(0.0, 1.0, 0.0),
                radius: 0.25,
            },
            RenderObject::Sphere {
                center: Vec3::new(2.0, 0.0, 0.0),
                radius: 0.5,
            },
        ];
        let groups = [2, 9, 0, 2];

        let mesh =
            flatten_render_objects_with_groups(&objects, &groups, &Molecule::default(), &options);

        let sphere_faces = molstar_sphere_triangle_count(0);
        let cylinder_segments = molstar_export_cylinder_radial_segments(4);
        let cylinder_faces = cylinder_segments * 4;
        let sphere_vertices = 12;
        let cylinder_vertices = (cylinder_segments + 1) * 4;
        let mut expected_vertex_groups = Vec::new();
        expected_vertex_groups.extend(std::iter::repeat_n(2, sphere_vertices));
        expected_vertex_groups.extend(std::iter::repeat_n(0, cylinder_vertices));
        expected_vertex_groups.extend(std::iter::repeat_n(2, sphere_vertices));
        let mut expected = Vec::new();
        expected.extend(std::iter::repeat_n(2, sphere_faces));
        expected.extend(std::iter::repeat_n(0, cylinder_faces));
        expected.extend(std::iter::repeat_n(2, sphere_faces));

        assert_eq!(
            mesh.faces.len(),
            sphere_faces + cylinder_faces + sphere_faces
        );
        assert_eq!(mesh.vertices.len(), mesh.normals.len());
        assert_eq!(mesh.vertex_groups.len(), mesh.vertices.len());
        assert_eq!(mesh.vertex_groups, expected_vertex_groups);
        assert_eq!(mesh.face_groups, expected);
        assert_eq!(mesh.group_count, 3);
    }

    fn assert_vec3_close(actual: Vec3, expected: Vec3) {
        assert!(
            actual.distance(expected) < 0.000_001,
            "actual={actual:?} expected={expected:?}"
        );
    }

    fn assert_faces_eq(actual: &[Face], expected: &[[usize; 3]]) {
        assert_eq!(actual.len(), expected.len());
        for (i, (actual, expected)) in actual.iter().zip(expected).enumerate() {
            assert_eq!(
                [actual.a, actual.b, actual.c],
                *expected,
                "face {i} does not match Mol* primitive order"
            );
        }
    }
}
