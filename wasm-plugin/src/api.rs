use std::borrow::Cow;

use crate::export::{
    export_maquette_material_map_json, export_maquette_material_map_json_from_obj,
    export_mtl_from_materials, export_obj_with_metadata, export_ply_with_metadata,
    export_stl_facet_context_json, export_stl_with_metadata, ExportMetadata,
    ExportOperatorMetadata, ExportVec3,
};
use crate::json::{json_escape, json_string_array};
use crate::mesh::{
    build_mesh, build_mesh_with_visible_bounding_sphere,
    build_mesh_with_visible_bounding_sphere_and_operator_snapshot,
    build_render_scene_with_summaries, render_materials,
    render_object_stl_facet_context_for_geometry_json_timed, render_object_stl_facet_context_json,
    render_object_stl_facet_context_json_timed, render_summaries_json,
    visible_renderable_bounding_sphere_for_export_with_structure,
    visible_renderable_bounding_sphere_report_for_export_with_structure, BondMetadataSnapshot,
    RenderScene, RenderSummaries,
};
use crate::model::{
    Assembly, AtomicStructure, BoundingSphere, Mesh, Molecule, MoleculeType, SecondaryRange,
    SourceData, UnitKind, UnitOperator, Vec3,
};
use crate::options::MeshOptions;
use crate::parser::{
    parse_molecule_with_options, parse_molecule_with_options_and_metadata, ParsedMolecule,
};

pub fn convert_to_obj(data: &[u8], options_json: &[u8]) -> Result<Vec<u8>, String> {
    let options = MeshOptions::from_json(options_json)?;
    let molecule = parse_molecule_with_options(data, &options)?;
    let (mesh, export_center, assembly_operators) =
        build_mesh_for_obj_export_with_operator_snapshot(&molecule, &options);
    validate_mesh_for_export(&mesh)?;
    let mut metadata = export_metadata_for_molecule(&molecule, &options, &assembly_operators);
    metadata.vertex_offset = export_center.negated();
    Ok(export_obj_with_metadata(&mesh, &metadata).into_bytes())
}

pub fn convert_to_obj_bundle(data: &[u8], options_json: &[u8]) -> Result<Vec<u8>, String> {
    let options = MeshOptions::from_json(options_json)?;
    let molecule = parse_molecule_with_options(data, &options)?;
    let (mesh, export_center, assembly_operators) =
        build_mesh_for_obj_export_with_operator_snapshot(&molecule, &options);
    validate_mesh_for_export(&mesh)?;
    let mut metadata = export_metadata_for_molecule(&molecule, &options, &assembly_operators);
    metadata.vertex_offset = export_center.negated();
    let materials_json = export_maquette_material_map_json(&mesh);
    let obj = export_obj_with_metadata(&mesh, &metadata);
    obj_bundle(&materials_json, obj)
}

pub fn convert_to_mtl(data: &[u8], options_json: &[u8]) -> Result<Vec<u8>, String> {
    let options = MeshOptions::from_json(options_json)?;
    let molecule = parse_molecule_with_options(data, &options)?;
    Ok(export_mtl_from_materials(&render_materials(&molecule, &options)).into_bytes())
}

pub fn maquette_material_map(obj: &[u8]) -> Result<Vec<u8>, String> {
    Ok(export_maquette_material_map_json_from_obj(obj)?.into_bytes())
}

fn obj_bundle(materials_json: &str, obj: String) -> Result<Vec<u8>, String> {
    const MAX_HEADER_LEN: usize = 99_999_999;
    let materials_len = materials_json.len();
    if materials_len > MAX_HEADER_LEN {
        return Err("OBJ material map is too large to bundle".to_string());
    }
    let mut out = Vec::with_capacity(8 + materials_json.len() + obj.len());
    out.extend_from_slice(format!("{materials_len:08}").as_bytes());
    out.extend_from_slice(materials_json.as_bytes());
    out.extend_from_slice(obj.as_bytes());
    Ok(out)
}

fn render_object_bundle(
    materials_json: &str,
    info_json: &str,
    mesh: Vec<u8>,
) -> Result<Vec<u8>, String> {
    const MAX_HEADER_LEN: usize = 99_999_999;
    let materials_len = materials_json.len();
    let info_len = info_json.len();
    if materials_len > MAX_HEADER_LEN {
        return Err("render-object material map is too large to bundle".to_string());
    }
    if info_len > MAX_HEADER_LEN {
        return Err("render-object info is too large to bundle".to_string());
    }
    let mut out = Vec::with_capacity(16 + materials_len + info_len + mesh.len());
    out.extend_from_slice(format!("{materials_len:08}").as_bytes());
    out.extend_from_slice(format!("{info_len:08}").as_bytes());
    out.extend_from_slice(materials_json.as_bytes());
    out.extend_from_slice(info_json.as_bytes());
    out.extend_from_slice(&mesh);
    Ok(out)
}

pub fn convert_to_stl(data: &[u8], options_json: &[u8]) -> Result<Vec<u8>, String> {
    let options = MeshOptions::from_json(options_json)?;
    let molecule = parse_molecule_with_options(data, &options)?;
    let (mesh, export_center) = build_mesh_for_obj_export(&molecule, &options);
    validate_mesh_for_export(&mesh)?;
    let metadata = ExportMetadata {
        vertex_offset: export_center.negated(),
        ..ExportMetadata::default()
    };
    Ok(export_stl_with_metadata(&mesh, &metadata))
}

pub fn convert_to_ply(data: &[u8], options_json: &[u8]) -> Result<Vec<u8>, String> {
    let options = MeshOptions::from_json(options_json)?;
    let molecule = parse_molecule_with_options(data, &options)?;
    let (mesh, assembly_operators) =
        build_mesh_for_export_with_operator_snapshot(&molecule, &options);
    validate_mesh_for_export(&mesh)?;
    let metadata = export_metadata_for_molecule(&molecule, &options, &assembly_operators);
    Ok(export_ply_with_metadata(&mesh, &metadata).into_bytes())
}

pub fn convert_to_render_object_bundle(
    data: &[u8],
    options_json: &[u8],
) -> Result<Vec<u8>, String> {
    let options = MeshOptions::from_json(options_json)?;
    let mesh_format = render_object_mesh_format(options_json)?;
    let ParsedMolecule {
        molecule,
        available_alt_locs,
    } = parse_molecule_with_options_and_metadata(data, &options)?;
    let (mut scene, export_center) = build_render_scene_for_export(&molecule, &options);
    let info_json = molecule_info_json_with_summaries(
        &options,
        &available_alt_locs,
        &molecule,
        Some(&scene.summaries),
    );

    match mesh_format {
        RenderObjectMeshFormat::Obj => {
            validate_mesh_for_export(&scene.mesh)?;
            let mut metadata =
                export_metadata_for_molecule(&molecule, &options, &scene.assembly_operators);
            metadata.vertex_offset = export_center.negated();
            let materials_json = export_maquette_material_map_json(&scene.mesh);
            let obj = export_obj_with_metadata(&scene.mesh, &metadata).into_bytes();
            render_object_bundle(&materials_json, &info_json, obj)
        }
        RenderObjectMeshFormat::Stl => {
            validate_mesh_for_export(&scene.mesh)?;
            let metadata = ExportMetadata {
                vertex_offset: export_center.negated(),
                ..ExportMetadata::default()
            };
            render_object_bundle(
                "{}",
                &info_json,
                export_stl_with_metadata(&scene.mesh, &metadata),
            )
        }
        RenderObjectMeshFormat::Ply => {
            if options.center {
                center_mesh_on_export_center(&mut scene.mesh, export_center.to_vec3());
            }
            validate_mesh_for_export(&scene.mesh)?;
            let metadata =
                export_metadata_for_molecule(&molecule, &options, &scene.assembly_operators);
            render_object_bundle(
                "{}",
                &info_json,
                export_ply_with_metadata(&scene.mesh, &metadata).into_bytes(),
            )
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RenderObjectMeshFormat {
    Obj,
    Stl,
    Ply,
}

fn render_object_mesh_format(options_json: &[u8]) -> Result<RenderObjectMeshFormat, String> {
    let text =
        std::str::from_utf8(options_json).map_err(|_| "options must be UTF-8 JSON".to_string())?;
    let value = json_string_field(text, "mesh-format")
        .or_else(|| json_string_field(text, "mesh_format"))
        .unwrap_or_else(|| "obj".to_string());
    match value.to_ascii_lowercase().as_str() {
        "obj" => Ok(RenderObjectMeshFormat::Obj),
        "stl" => Ok(RenderObjectMeshFormat::Stl),
        "ply" => Ok(RenderObjectMeshFormat::Ply),
        other => Err(format!(
            "mesh-format must be one of \"obj\", \"stl\", or \"ply\"; got {other}"
        )),
    }
}

fn json_string_field(text: &str, key: &str) -> Option<String> {
    let marker = format!("\"{key}\"");
    let mut tail = text.split_once(&marker)?.1;
    tail = tail.split_once(':')?.1.trim_start();
    if !tail.starts_with('"') {
        return None;
    }
    let mut value = String::new();
    let mut escape = false;
    for ch in tail[1..].chars() {
        if escape {
            value.push(ch);
            escape = false;
        } else if ch == '\\' {
            escape = true;
        } else if ch == '"' {
            return Some(value);
        } else {
            value.push(ch);
        }
    }
    None
}

pub fn stl_facet_semantic_context(
    data: &[u8],
    options_json: &[u8],
    facet_index: usize,
) -> Result<String, String> {
    let options = MeshOptions::from_json(options_json)?;
    let molecule = parse_molecule_with_options(data, &options)?;
    let mut geometry_options = options.clone();
    let export_center = if options.center {
        geometry_options.center = false;
        let (mesh, visible_sphere) =
            build_mesh_with_visible_bounding_sphere(&molecule, &geometry_options);
        visible_sphere
            .map(|sphere| export_box_center_from_sphere(&sphere))
            .unwrap_or_else(|| export_mesh_boundary_sphere_center(&mesh))
    } else {
        ExportVec3::default()
    };
    Ok(render_object_stl_facet_context_json(
        &molecule,
        &geometry_options,
        facet_index,
        [-export_center.x, -export_center.y, -export_center.z],
    ))
}

pub fn stl_export_facet_context(
    data: &[u8],
    options_json: &[u8],
    facet_index: usize,
) -> Result<String, String> {
    stl_export_facet_context_timed(data, options_json, facet_index, |_| {})
}

pub fn stl_export_facet_context_timed(
    data: &[u8],
    options_json: &[u8],
    facet_index: usize,
    mut checkpoint: impl FnMut(&str),
) -> Result<String, String> {
    let options = MeshOptions::from_json(options_json)?;
    checkpoint("parse-options");
    let molecule = parse_molecule_with_options(data, &options)?;
    checkpoint("parse-molecule");
    let mut geometry_options = options.clone();
    geometry_options.center = false;
    let geometry = export_context_geometry(&molecule);
    checkpoint("export-context-geometry");
    let structure = geometry.atomic_structure();
    checkpoint("atomic-structure-for-export-context");
    let visible_sphere_report = options.center.then(|| {
        visible_renderable_bounding_sphere_report_for_export_with_structure(
            &geometry,
            &geometry_options,
            &structure,
        )
    });
    if visible_sphere_report.is_some() {
        checkpoint("visible-renderable-bounding-sphere-report");
    }
    let (export_center, export_box) = if options.center {
        if let Some(sphere) = visible_renderable_bounding_sphere_for_export_with_structure(
            &geometry,
            &geometry_options,
            &structure,
        ) {
            checkpoint("visible-renderable-bounding-sphere");
            let export_box = export_box_context_from_sphere(&sphere);
            (
                export_box_center(export_box.min, export_box.max),
                Some(export_box),
            )
        } else {
            let (mesh, export_center) = build_mesh_for_obj_export(&geometry, &options);
            checkpoint("build-export-mesh");
            validate_mesh_for_export(&mesh)?;
            checkpoint("validate-export-mesh");
            let context = export_stl_facet_context_json(&mesh, export_center, facet_index);
            checkpoint("render-export-facet-context");
            return Ok(context);
        }
    } else {
        checkpoint("skip-export-center");
        (ExportVec3::default(), None)
    };
    let vertex_offset = export_center.negated();
    let context = render_object_stl_facet_context_for_geometry_json_timed(
        &geometry,
        &geometry_options,
        facet_index,
        [vertex_offset.x, vertex_offset.y, vertex_offset.z],
        Some(&structure),
        |label| checkpoint(label),
    );
    checkpoint("build-flat-export-mesh-context");
    let (flat_mesh, flat_export_center) = build_mesh_for_obj_export(&geometry, &options);
    let flat_context = export_stl_facet_context_json(&flat_mesh, flat_export_center, facet_index);
    checkpoint("render-flat-export-facet-context");
    let center_source = if options.center {
        "visible-renderable-bounding-sphere"
    } else {
        "disabled"
    };
    Ok(prefix_export_facet_context_json(
        &context,
        &flat_context,
        facet_index,
        export_center,
        export_box,
        options.center,
        center_source,
        visible_sphere_report.as_deref(),
    ))
}

pub fn stl_facet_semantic_context_with_vertex_offset(
    data: &[u8],
    options_json: &[u8],
    facet_index: usize,
    vertex_offset: [f64; 3],
) -> Result<String, String> {
    stl_facet_semantic_context_with_vertex_offset_timed(
        data,
        options_json,
        facet_index,
        vertex_offset,
        |_| {},
    )
}

pub fn stl_facet_semantic_context_with_vertex_offset_timed(
    data: &[u8],
    options_json: &[u8],
    facet_index: usize,
    vertex_offset: [f64; 3],
    mut checkpoint: impl FnMut(&str),
) -> Result<String, String> {
    let options = MeshOptions::from_json(options_json)?;
    checkpoint("parse-options");
    let molecule = parse_molecule_with_options(data, &options)?;
    checkpoint("parse-molecule");
    let mut geometry_options = options.clone();
    geometry_options.center = false;
    let context = render_object_stl_facet_context_json_timed(
        &molecule,
        &geometry_options,
        facet_index,
        vertex_offset,
        checkpoint,
    );
    Ok(context)
}

fn export_context_geometry(molecule: &Molecule) -> Cow<'_, Molecule> {
    if let Some(molecule) = molecule.identity_assembly_subset_for_geometry() {
        return Cow::Owned(molecule);
    }
    if molecule.selected_assembly.is_some() {
        molecule.expanded_for_geometry()
    } else {
        Cow::Borrowed(molecule)
    }
}

fn prefix_export_facet_context_json(
    context: &str,
    flat_context: &str,
    stl_facet: usize,
    export_center: ExportVec3,
    export_box: Option<ExportBoxContext>,
    centered: bool,
    center_source: &str,
    visible_sphere_report: Option<&str>,
) -> String {
    let sparse_slot_has_face = context.contains("\"found\":true") && stl_facet % 3 == 0;
    let box_json = export_box
        .map(|context| context.json_fields())
        .unwrap_or_else(|| {
            "\"export_box_min\":null,\"export_box_max\":null,\"export_box_extrema_count\":0,\"export_box_min_indices\":[null,null,null],\"export_box_max_indices\":[null,null,null],\"export_box_min_points\":[null,null,null],\"export_box_max_points\":[null,null,null]".to_string()
        });
    let visible_sphere_report = visible_sphere_report.unwrap_or("null");
    let prefix = format!(
        "\"export_center\":{},{},\"export_centered\":{},\"export_center_source\":\"{}\",\"visible_sphere_report\":{},\"sparse_slot_has_face\":{},",
        export_vec3_json(export_center),
        box_json,
        centered,
        json_escape(center_source),
        visible_sphere_report,
        sparse_slot_has_face,
    );
    if let Some(rest) = context.strip_prefix('{') {
        let mut json = format!("{{{prefix}{rest}");
        json.pop();
        json.push_str(",\"flat_export_context\":");
        json.push_str(flat_context);
        json.push('}');
        json
    } else {
        format!(
            "{{{prefix}\"context_parse_error\":\"expected object\",\"context\":\"{}\"}}",
            json_escape(context)
        )
    }
}

fn export_vec3_json(value: ExportVec3) -> String {
    format!("[{:.17},{:.17},{:.17}]", value.x, value.y, value.z)
}

fn build_mesh_for_export_with_operator_snapshot(
    molecule: &Molecule,
    options: &MeshOptions,
) -> (Mesh, Vec<UnitOperator>) {
    if !options.center {
        let (mesh, _, operators) = build_mesh_with_visible_bounding_sphere_and_operator_snapshot(
            molecule,
            options,
            options.include_operator_metadata,
        );
        return (mesh, operators);
    }

    let mut geometry_options = options.clone();
    geometry_options.center = false;
    let (mut mesh, visible_sphere, operators) =
        build_mesh_with_visible_bounding_sphere_and_operator_snapshot(
            molecule,
            &geometry_options,
            options.include_operator_metadata,
        );
    let center = visible_sphere
        .map(|sphere| export_box_center_from_sphere(&sphere))
        .unwrap_or_else(|| export_mesh_boundary_sphere_center(&mesh));
    center_mesh_on_export_center(&mut mesh, center.to_vec3());
    (mesh, operators)
}

fn build_mesh_for_obj_export(molecule: &Molecule, options: &MeshOptions) -> (Mesh, ExportVec3) {
    if !options.center {
        return (build_mesh(molecule, options), ExportVec3::default());
    }

    let mut geometry_options = options.clone();
    geometry_options.center = false;
    let (mesh, visible_sphere) =
        build_mesh_with_visible_bounding_sphere(molecule, &geometry_options);
    let center = visible_sphere
        .map(|sphere| export_box_center_from_sphere(&sphere))
        .unwrap_or_else(|| export_mesh_boundary_sphere_center(&mesh));
    (mesh, center)
}

fn build_mesh_for_obj_export_with_operator_snapshot(
    molecule: &Molecule,
    options: &MeshOptions,
) -> (Mesh, ExportVec3, Vec<UnitOperator>) {
    if !options.center {
        let (mesh, _, operators) = build_mesh_with_visible_bounding_sphere_and_operator_snapshot(
            molecule,
            options,
            options.include_operator_metadata,
        );
        return (mesh, ExportVec3::default(), operators);
    }

    let mut geometry_options = options.clone();
    geometry_options.center = false;
    let (mesh, visible_sphere, operators) =
        build_mesh_with_visible_bounding_sphere_and_operator_snapshot(
            molecule,
            &geometry_options,
            options.include_operator_metadata,
        );
    let center = visible_sphere
        .map(|sphere| export_box_center_from_sphere(&sphere))
        .unwrap_or_else(|| export_mesh_boundary_sphere_center(&mesh));
    (mesh, center, operators)
}

fn build_render_scene_for_export(
    molecule: &Molecule,
    options: &MeshOptions,
) -> (RenderScene, ExportVec3) {
    if !options.center {
        return (
            build_render_scene_with_summaries(molecule, options),
            ExportVec3::default(),
        );
    }

    let mut geometry_options = options.clone();
    geometry_options.center = false;
    let scene = build_render_scene_with_summaries(molecule, &geometry_options);
    let center = scene
        .visible_bounding_sphere
        .as_ref()
        .map(export_box_center_from_sphere)
        .unwrap_or_else(|| export_mesh_boundary_sphere_center(&scene.mesh));
    (scene, center)
}

fn center_mesh_on_export_center(mesh: &mut Mesh, center: Vec3) {
    if mesh.vertices.is_empty() {
        return;
    }

    if center == Vec3::default() {
        return;
    }

    for vertex in &mut mesh.vertices {
        *vertex = *vertex - center;
    }
}

fn export_mesh_boundary_sphere_center(mesh: &Mesh) -> ExportVec3 {
    let sphere = crate::model::Boundary::from_positions(&mesh.vertices).sphere;
    export_box_center_from_sphere(&sphere)
}

fn export_box_center_from_sphere(sphere: &BoundingSphere) -> ExportVec3 {
    let context = export_box_context_from_sphere(sphere);
    export_box_center(context.min, context.max)
}

fn export_box_center(min: ExportVec3, max: ExportVec3) -> ExportVec3 {
    ExportVec3::new(
        (min.x + max.x) * 0.5,
        (min.y + max.y) * 0.5,
        (min.z + max.z) * 0.5,
    )
}

#[derive(Clone, Copy, Debug)]
struct ExportBoxContext {
    min: ExportVec3,
    max: ExportVec3,
    extrema_count: usize,
    min_indices: [Option<usize>; 3],
    max_indices: [Option<usize>; 3],
    min_points: [Option<ExportVec3>; 3],
    max_points: [Option<ExportVec3>; 3],
}

impl ExportBoxContext {
    fn json_fields(self) -> String {
        format!(
            "\"export_box_min\":{},\"export_box_max\":{},\"export_box_extrema_count\":{},\"export_box_min_indices\":{},\"export_box_max_indices\":{},\"export_box_min_points\":{},\"export_box_max_points\":{}",
            export_vec3_json(self.min),
            export_vec3_json(self.max),
            self.extrema_count,
            option_usize_array_json(self.min_indices),
            option_usize_array_json(self.max_indices),
            option_vec3_array_json(self.min_points),
            option_vec3_array_json(self.max_points),
        )
    }
}

fn export_box_context_from_sphere(sphere: &BoundingSphere) -> ExportBoxContext {
    if sphere.extrema64.len() >= 14 {
        export_box_context_from_points(
            sphere
                .extrema64
                .iter()
                .copied()
                .map(|[x, y, z]| ExportVec3::new(x, y, z)),
        )
    } else if sphere.extrema.len() >= 14 {
        export_box_context_from_points(
            sphere
                .extrema
                .iter()
                .map(|point| ExportVec3::new(point.x as f64, point.y as f64, point.z as f64)),
        )
    } else {
        let center = ExportVec3::new(
            sphere.center.x as f64,
            sphere.center.y as f64,
            sphere.center.z as f64,
        );
        let radius = sphere.radius as f64;
        ExportBoxContext {
            min: ExportVec3::new(center.x - radius, center.y - radius, center.z - radius),
            max: ExportVec3::new(center.x + radius, center.y + radius, center.z + radius),
            extrema_count: 0,
            min_indices: [None, None, None],
            max_indices: [None, None, None],
            min_points: [None, None, None],
            max_points: [None, None, None],
        }
    }
}

fn export_box_context_from_points(
    points: impl IntoIterator<Item = ExportVec3>,
) -> ExportBoxContext {
    let mut min = [f64::INFINITY, f64::INFINITY, f64::INFINITY];
    let mut max = [f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
    let mut min_indices = [None, None, None];
    let mut max_indices = [None, None, None];
    let mut min_points = [None, None, None];
    let mut max_points = [None, None, None];
    let mut count = 0usize;
    for point in points {
        let values = [point.x, point.y, point.z];
        for axis in 0..3 {
            if values[axis] < min[axis] {
                min[axis] = values[axis];
                min_indices[axis] = Some(count);
                min_points[axis] = Some(point);
            }
            if values[axis] > max[axis] {
                max[axis] = values[axis];
                max_indices[axis] = Some(count);
                max_points[axis] = Some(point);
            }
        }
        count += 1;
    }
    ExportBoxContext {
        min: ExportVec3::new(min[0], min[1], min[2]),
        max: ExportVec3::new(max[0], max[1], max[2]),
        extrema_count: count,
        min_indices,
        max_indices,
        min_points,
        max_points,
    }
}

fn option_usize_array_json(values: [Option<usize>; 3]) -> String {
    let items = values
        .into_iter()
        .map(|value| {
            value
                .map(|value| value.to_string())
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{items}]")
}

fn option_vec3_array_json(values: [Option<ExportVec3>; 3]) -> String {
    let items = values
        .into_iter()
        .map(|value| {
            value
                .map(export_vec3_json)
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{items}]")
}

pub(crate) fn validate_mesh_for_export(mesh: &Mesh) -> Result<(), String> {
    for (index, vertex) in mesh.vertices.iter().enumerate() {
        if !vec3_is_finite(*vertex) {
            return Err(format!("mesh vertex {index} contains NaN or infinity"));
        }
    }
    for (index, normal) in mesh.normals.iter().enumerate() {
        if !vec3_is_finite(*normal) {
            return Err(format!("mesh normal {index} contains NaN or infinity"));
        }
    }
    if mesh.normals.len() != mesh.vertices.len() {
        return Err(format!(
            "mesh normal count {} does not match vertex count {}",
            mesh.normals.len(),
            mesh.vertices.len()
        ));
    }
    for (index, face) in mesh.faces.iter().enumerate() {
        if face.a >= mesh.vertices.len()
            || face.b >= mesh.vertices.len()
            || face.c >= mesh.vertices.len()
        {
            return Err(format!(
                "mesh face {index} references an out-of-range vertex"
            ));
        }
    }
    if !mesh.vertex_groups.is_empty() && mesh.vertex_groups.len() != mesh.vertices.len() {
        return Err(format!(
            "mesh vertex group count {} does not match vertex count {}",
            mesh.vertex_groups.len(),
            mesh.vertices.len()
        ));
    }
    if mesh.face_groups.len() != mesh.faces.len() {
        return Err(format!(
            "mesh face group count {} does not match face count {}",
            mesh.face_groups.len(),
            mesh.faces.len()
        ));
    }
    if !mesh.face_materials.is_empty() && mesh.face_materials.len() != mesh.faces.len() {
        return Err(format!(
            "mesh face material count {} does not match face count {}",
            mesh.face_materials.len(),
            mesh.faces.len()
        ));
    }
    validate_mesh_sections_for_export(mesh)?;
    Ok(())
}

fn vec3_is_finite(value: Vec3) -> bool {
    value.x.is_finite() && value.y.is_finite() && value.z.is_finite()
}

fn validate_mesh_sections_for_export(mesh: &Mesh) -> Result<(), String> {
    let mut vertex_cursor = 0;
    let mut face_cursor = 0;
    for (index, section) in mesh.sections.iter().enumerate() {
        if section.key.is_empty() {
            return Err(format!("mesh section {index} has an empty key"));
        }
        if section.vertex_start != vertex_cursor {
            return Err(format!(
                "mesh section {index} starts at vertex {} but expected {vertex_cursor}",
                section.vertex_start
            ));
        }
        if section.vertex_end < section.vertex_start || section.vertex_end > mesh.vertices.len() {
            return Err(format!("mesh section {index} has an invalid vertex range"));
        }
        if section.face_start != face_cursor {
            return Err(format!(
                "mesh section {index} starts at face {} but expected {face_cursor}",
                section.face_start
            ));
        }
        if section.face_end < section.face_start || section.face_end > mesh.faces.len() {
            return Err(format!("mesh section {index} has an invalid face range"));
        }
        vertex_cursor = section.vertex_end;
        face_cursor = section.face_end;
    }
    if !mesh.sections.is_empty()
        && (vertex_cursor != mesh.vertices.len() || face_cursor != mesh.faces.len())
    {
        return Err("mesh sections do not cover the full exported mesh".to_string());
    }
    Ok(())
}

fn export_metadata_for_molecule(
    molecule: &Molecule,
    options: &MeshOptions,
    assembly_operators: &[UnitOperator],
) -> ExportMetadata {
    let mut metadata = ExportMetadata {
        obj_basename: options.obj_basename.clone(),
        include_operator_metadata: options.include_operator_metadata,
        include_face_groups: options.obj_groups,
        ..ExportMetadata::default()
    };
    let Some(assembly) = &molecule.selected_assembly else {
        return metadata;
    };

    metadata.assembly_id = Some(assembly.id.clone());
    if !options.include_operator_metadata {
        return metadata;
    }

    for operator in assembly_operators {
        if metadata
            .operators
            .iter()
            .any(|existing: &ExportOperatorMetadata| {
                existing.name == operator.name
                    && existing.instance_id == operator.instance_id
                    && existing.assembly_id == operator.assembly_id
                    && existing.oper_id == operator.oper_id
                    && existing.oper_list_ids == operator.oper_list_ids
            })
        {
            continue;
        }
        metadata.operators.push(ExportOperatorMetadata {
            name: operator.name.clone(),
            instance_id: operator.instance_id.clone(),
            assembly_id: operator.assembly_id.clone(),
            oper_id: operator.oper_id,
            oper_list_ids: operator.oper_list_ids.clone(),
            is_identity: operator.is_identity,
        });
    }

    metadata
}

pub fn molecule_info(data: &[u8], options_json: &[u8]) -> Result<Vec<u8>, String> {
    let options = MeshOptions::from_json(options_json)?;
    let ParsedMolecule {
        molecule,
        available_alt_locs,
    } = parse_molecule_with_options_and_metadata(data, &options)?;
    Ok(
        molecule_info_json_with_summaries(&options, &available_alt_locs, &molecule, None)
            .into_bytes(),
    )
}

fn molecule_info_json_with_summaries(
    options: &MeshOptions,
    available_alt_locs: &[String],
    molecule: &Molecule,
    summaries: Option<&RenderSummaries>,
) -> String {
    let assemblies = molecule
        .assemblies
        .iter()
        .map(assembly_json)
        .collect::<Vec<_>>()
        .join(",");
    let helices = molecule
        .helices
        .iter()
        .map(|range| secondary_range_json("helix", range))
        .collect::<Vec<_>>()
        .join(",");
    let sheets = molecule
        .sheets
        .iter()
        .map(|range| secondary_range_json("sheet", range))
        .collect::<Vec<_>>()
        .join(",");
    let summaries_storage;
    let summaries = if let Some(summaries) = summaries {
        summaries
    } else {
        summaries_storage = render_summaries_json(molecule, options);
        &summaries_storage
    };
    let geometry = &summaries.geometry;
    let (min, max) = geometry.bounds.unwrap_or_default();
    let assembly_id = options.assembly.as_deref().unwrap_or("asymmetric-unit");
    format!(
        "{{\"entry_count\":{},\"experiment_count\":{},\"atom_count\":{},\"anisotrop_count\":{},\"coarse_sphere_count\":{},\"coarse_gaussian_count\":{},\"entity_count\":{},\"entity_poly_count\":{},\"entity_poly_seq_count\":{},\"pdbx_entity_branch_count\":{},\"pdbx_entity_branch_link_count\":{},\"pdbx_branch_scheme_count\":{},\"pdbx_nonpoly_scheme_count\":{},\"pdbx_poly_seq_scheme_count\":{},\"ihm_model_count\":{},\"ihm_model_group_count\":{},\"ihm_model_group_link_count\":{},\"ihm_cross_link_restraint_count\":{},\"struct_asym_count\":{},\"chem_comp_count\":{},\"chem_comp_atom_count\":{},\"chem_comp_bond_count\":{},\"chem_comp_angle_count\":{},\"bond_count\":{},\"bond_metadata\":{},\"assembly_count\":{},\"alt_locs\":{},\"alt_locs_info\":{{\"policy\":\"{}\",\"available\":{}}},\"assemblies\":[{}],\"assembly\":{{\"id\":\"{}\"}},\"source_data\":{},\"structure\":{},\"secondary_structure\":{{\"helices\":[{}],\"sheets\":[{}]}},\"representation\":{},\"render_objects\":{},\"bounds\":{{\"min\":[{:.4},{:.4},{:.4}],\"max\":[{:.4},{:.4},{:.4}]}}}}",
        molecule.entries.len(),
        molecule.experiments.len(),
        geometry.atom_count,
        molecule.atom_site_anisotrop.len(),
        geometry.coarse_sphere_count,
        geometry.coarse_gaussian_count,
        molecule.entities.len(),
        molecule.entity_polymers.len(),
        molecule.entity_poly_seq.len(),
        molecule.pdbx_entity_branch.len(),
        molecule.pdbx_entity_branch_links.len(),
        molecule.pdbx_branch_scheme.len(),
        molecule.pdbx_nonpoly_scheme.len(),
        molecule.pdbx_poly_seq_scheme.len(),
        molecule.ihm_model_list.len(),
        molecule.ihm_model_groups.len(),
        molecule.ihm_model_group_links.len(),
        molecule.ihm_cross_link_restraints.len(),
        molecule.struct_asym.len(),
        molecule.chemical_components.len(),
        molecule.chemical_component_atoms.len(),
        molecule.chemical_component_bonds.len(),
        molecule.chemical_component_angles.len(),
        geometry.bond_count,
        bond_metadata_json(&geometry.bond_metadata),
        molecule.assemblies.len(),
        json_string_array(available_alt_locs),
        json_escape(if options.alt_loc.is_empty() { "default" } else { &options.alt_loc }),
        json_string_array(available_alt_locs),
        assemblies,
        json_escape(assembly_id),
        source_data_json(&molecule.source_data),
        atomic_structure_json(&summaries.structure),
        helices,
        sheets,
        summaries.representation_json,
        summaries.render_objects_json,
        min.x,
        min.y,
        min.z,
        max.x,
        max.y,
        max.z
    )
}

fn source_data_json(source_data: &SourceData) -> String {
    let categories = source_categories_json(&source_data.categories);
    let db_categories = source_categories_json(&source_data.db_categories);
    let frame_categories = source_categories_json(&source_data.frame_categories);
    format!(
        "{{\"kind\":\"{}\",\"name\":\"{}\",\"original_kind\":\"{}\",\"category_count\":{},\"categories\":[{}],\"db\":{{\"category_count\":{},\"categories\":[{}]}},\"frame\":{{\"category_count\":{},\"categories\":[{}]}}}}",
        json_escape(&source_data.kind),
        json_escape(&source_data.name),
        json_escape(&source_data.original_kind),
        source_data.categories.len(),
        categories,
        source_data.db_categories.len(),
        db_categories,
        source_data.frame_categories.len(),
        frame_categories
    )
}

fn source_categories_json(categories: &[crate::model::SourceCategory]) -> String {
    categories
        .iter()
        .map(|category| {
            format!(
                "{{\"name\":\"{}\",\"row_count\":{},\"column_count\":{}}}",
                json_escape(&category.name),
                category.row_count,
                category.column_count
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn bond_metadata_json(metadata: &BondMetadataSnapshot) -> String {
    format!(
        "{{\"count\":{},\"computed\":{},\"pdb_conect\":{},\"struct_conn\":{},\"index_pair\":{},\"chem_comp\":{},\"covalent\":{},\"metallic_coordination\":{},\"hydrogen_bond\":{},\"disulfide\":{},\"aromatic\":{},\"computed_flag\":{},\"resonance\":{},\"rings\":{},\"aromatic_rings\":{},\"delocalized_bonds\":{}}}",
        metadata.count,
        metadata.computed,
        metadata.pdb_conect,
        metadata.struct_conn,
        metadata.index_pair,
        metadata.chem_comp,
        metadata.covalent,
        metadata.metallic_coordination,
        metadata.hydrogen_bond,
        metadata.disulfide,
        metadata.aromatic,
        metadata.computed_flag,
        metadata.resonance,
        metadata.rings,
        metadata.aromatic_rings,
        metadata.delocalized_bonds
    )
}

fn atomic_structure_json(structure: &AtomicStructure) -> String {
    let count_unit_kind = |kind: UnitKind| {
        structure
            .units
            .iter()
            .filter(|unit| unit.kind == kind)
            .count()
    };
    format!(
        "{{\"model_count\":{},\"unit_count\":{},\"unit_kind_counts\":{{\"atomic\":{},\"spheres\":{},\"gaussians\":{}}},\"symmetry_group_count\":{},\"element_count\":{},\"source_element_count\":{},\"polymer_residue_count\":{},\"polymer_gap_count\":{},\"bonding\":{{\"intra_unit_count\":{},\"inter_unit_count\":{}}},\"derived\":{{\"protein_residue_count\":{},\"nucleotide_residue_count\":{},\"water_residue_count\":{},\"trace_element_count\":{},\"atomic_number_count\":{}}},\"boundary\":{{\"sphere_center\":[{:.4},{:.4},{:.4}],\"sphere_radius\":{:.4},\"box_min\":[{:.4},{:.4},{:.4}],\"box_max\":[{:.4},{:.4},{:.4}]}},\"lookup3d\":{{\"unit_count\":{}}},\"coordinate_system\":{{\"name\":\"{}\",\"assembly_id\":\"{}\",\"oper_id\":{},\"is_identity\":{}}},\"chain_count\":{},\"residue_count\":{},\"atom_count\":{},\"conformation\":{{\"atom_id_count\":{},\"position_count\":{},\"occupancy_count\":{},\"b_iso_count\":{},\"formal_charge_count\":{},\"occupancy_defined\":{},\"b_iso_defined\":{},\"xyz_defined\":{},\"element_to_anisotrop_count\":{},\"anisotropic_displacement_count\":{}}},\"segments\":{{\"residue_count\":{},\"chain_count\":{}}},\"ranges\":{{\"polymer_count\":{},\"gap_count\":{},\"cyclic_count\":{}}},\"alt_loc_count\":{}}}",
        structure.models.len(),
        structure.units.len(),
        count_unit_kind(UnitKind::Atomic),
        count_unit_kind(UnitKind::Spheres),
        count_unit_kind(UnitKind::Gaussians),
        structure.symmetry_groups.len(),
        structure.element_count,
        structure.model.hierarchy.atoms.len(),
        structure.polymer_residue_count,
        structure.polymer_gap_count,
        structure.intra_unit_bond_count,
        structure.inter_unit_bonds.len(),
        structure
            .model
            .hierarchy
            .derived
            .molecule_type_count(MoleculeType::Protein),
        structure
            .model
            .hierarchy
            .derived
            .molecule_type_count(MoleculeType::Rna)
            + structure
                .model
                .hierarchy
                .derived
                .molecule_type_count(MoleculeType::Dna),
        structure
            .model
            .hierarchy
            .derived
            .molecule_type_count(MoleculeType::Water),
        structure.model.hierarchy.derived.trace_element_count(),
        structure.model.hierarchy.derived.atom.atomic_number.len(),
        structure.boundary.sphere.center.x,
        structure.boundary.sphere.center.y,
        structure.boundary.sphere.center.z,
        structure.boundary.sphere.radius,
        structure.boundary.box_min.x,
        structure.boundary.box_min.y,
        structure.boundary.box_min.z,
        structure.boundary.box_max.x,
        structure.boundary.box_max.y,
        structure.boundary.box_max.z,
        structure.units.len(),
        json_escape(&structure.coordinate_system.name),
        json_escape(&structure.coordinate_system.assembly_id),
        structure.coordinate_system.oper_id,
        if structure.coordinate_system.is_identity { "true" } else { "false" },
        structure.model.hierarchy.chains.len(),
        structure.model.hierarchy.residues.len(),
        structure.model.hierarchy.atoms.len(),
        structure.model.conformation.atom_ids.len(),
        structure.model.conformation.positions.len(),
        structure.model.conformation.occupancies.len(),
        structure.model.conformation.b_iso.len(),
        structure.model.conformation.formal_charges.len(),
        if structure.model.conformation.occupancy_defined {
            "true"
        } else {
            "false"
        },
        if structure.model.conformation.b_iso_defined {
            "true"
        } else {
            "false"
        },
        if structure.model.conformation.xyz_defined {
            "true"
        } else {
            "false"
        },
        structure.model.conformation.element_to_anisotrop.len(),
        structure
            .model
            .conformation
            .anisotropic_displacement
            .iter()
            .filter(|value| value.is_some())
            .count(),
        structure.model.hierarchy.residue_atom_segments.count,
        structure.model.hierarchy.chain_atom_segments.count,
        structure.ranges.polymer_ranges.len() / 2,
        structure.ranges.gap_ranges.len() / 2,
        structure.ranges.cyclic_polymer_map.len(),
        structure.alt_loc_count()
    )
}

fn assembly_json(assembly: &Assembly) -> String {
    let transform_count = assembly
        .generators
        .iter()
        .map(|g| g.transforms.len())
        .sum::<usize>();
    let generator_count = assembly.generators.len();
    let oligomeric_count = assembly
        .oligomeric_count
        .map(|value| value.to_string())
        .unwrap_or_else(|| "null".to_string());
    format!(
        "{{\"id\":\"{}\",\"details\":\"{}\",\"oligomeric_details\":\"{}\",\"oligomeric_count\":{},\"asym_ids\":{},\"operator_count\":{},\"generator_count\":{}}}",
        json_escape(&assembly.id),
        json_escape(&assembly.details),
        json_escape(&assembly.oligomeric_details),
        oligomeric_count,
        json_string_array(&assembly.asym_ids),
        transform_count,
        generator_count
    )
}

fn secondary_range_json(kind: &str, range: &SecondaryRange) -> String {
    format!(
        "{{\"type\":\"{}\",\"chain\":\"{}\",\"start\":{},\"start_insertion_code\":\"{}\",\"end\":{},\"end_insertion_code\":\"{}\"}}",
        json_escape(kind),
        json_escape(&range.chain),
        range.start,
        json_escape(&range.start_insertion_code),
        range.end,
        json_escape(&range.end_insertion_code)
    )
}
