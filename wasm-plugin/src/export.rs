use std::fmt::Write as _;

use crate::json::{json_escape, json_string_array};
use crate::model::{Mesh, MeshMaterial, Vec3};

#[cfg(test)]
mod tests;

#[derive(Clone, Debug)]
pub(crate) struct ExportMetadata {
    pub(crate) assembly_id: Option<String>,
    pub(crate) operators: Vec<ExportOperatorMetadata>,
    pub(crate) obj_basename: Option<String>,
    pub(crate) include_operator_metadata: bool,
    pub(crate) include_face_groups: bool,
    pub(crate) vertex_offset: ExportVec3,
}

impl Default for ExportMetadata {
    fn default() -> Self {
        Self {
            assembly_id: None,
            operators: Vec::new(),
            obj_basename: None,
            include_operator_metadata: true,
            include_face_groups: true,
            vertex_offset: ExportVec3::default(),
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct ExportVec3 {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

impl ExportVec3 {
    pub(crate) const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub(crate) fn negated(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }

    pub(crate) fn to_vec3(self) -> Vec3 {
        Vec3::new(self.x as f32, self.y as f32, self.z as f32)
    }
}

#[derive(Clone, Debug)]
pub(crate) struct ExportOperatorMetadata {
    pub(crate) name: String,
    pub(crate) instance_id: String,
    pub(crate) assembly_id: String,
    pub(crate) oper_id: i32,
    pub(crate) oper_list_ids: Vec<String>,
    pub(crate) is_identity: bool,
}

#[allow(dead_code)]
pub(crate) fn export_obj(mesh: &Mesh) -> String {
    export_obj_with_metadata(mesh, &ExportMetadata::default())
}

pub(crate) fn export_obj_with_metadata(mesh: &Mesh, metadata: &ExportMetadata) -> String {
    let basename = metadata.obj_basename.as_deref().unwrap_or("molfig");
    let mut out = String::with_capacity(obj_output_capacity(mesh));
    writeln!(out, "mtllib {basename}.mtl").expect("writing to String cannot fail");
    if metadata.include_operator_metadata {
        if let Some(metadata_json) = export_metadata_json(metadata) {
            out.push_str("# molfig_operator_metadata ");
            out.push_str(&metadata_json);
            out.push('\n');
        }
    }
    if let Some(sections) = molstar_obj_sections(mesh) {
        export_obj_sections(mesh, metadata, sections, &mut out);
    } else {
        export_obj_unsectioned(mesh, metadata, &mut out);
    }
    out
}

fn export_obj_sections(
    mesh: &Mesh,
    metadata: &ExportMetadata,
    sections: &[crate::model::MeshSection],
    out: &mut String,
) {
    let mut current_group = None;
    let mut current_material = None;
    for section in sections {
        for v in &mesh.vertices[section.vertex_start..section.vertex_end] {
            write_obj_vertex(out, *v, metadata.vertex_offset);
        }
        for n in &mesh.normals[section.vertex_start..section.vertex_end] {
            write_obj_normal(out, *n);
        }
        write_obj_faces(
            mesh,
            metadata,
            section.face_start,
            section.face_end,
            out,
            &mut current_group,
            &mut current_material,
        );
    }
}

fn export_obj_unsectioned(mesh: &Mesh, metadata: &ExportMetadata, out: &mut String) {
    for v in &mesh.vertices {
        write_obj_vertex(out, *v, metadata.vertex_offset);
    }
    for n in &mesh.normals {
        write_obj_normal(out, *n);
    }
    let mut current_group = None;
    let mut current_material = None;
    write_obj_faces(
        mesh,
        metadata,
        0,
        mesh.faces.len(),
        out,
        &mut current_group,
        &mut current_material,
    );
}

fn write_obj_vertex(out: &mut String, v: crate::model::Vec3, offset: ExportVec3) {
    let v = molstar_obj_vertex_transform(v, offset);
    writeln!(
        out,
        "v {} {} {}",
        molstar_rounded_float64(v[0], 1000.0),
        molstar_rounded_float64(v[1], 1000.0),
        molstar_rounded_float64(v[2], 1000.0)
    )
    .expect("writing to String cannot fail");
}

fn molstar_obj_vertex_transform(v: crate::model::Vec3, offset: ExportVec3) -> [f64; 3] {
    let x = v.x as f64;
    let y = v.y as f64;
    let z = v.z as f64;
    let m = [
        1.0, 0.0, 0.0, 0.0, //
        0.0, 1.0, 0.0, 0.0, //
        0.0, 0.0, 1.0, 0.0, //
        offset.x, offset.y, offset.z, 1.0,
    ];
    let denominator = m[3] * x + m[7] * y + m[11] * z + m[15];
    let w = 1.0 / if denominator != 0.0 { denominator } else { 1.0 };
    [
        (m[0] * x + m[4] * y + m[8] * z + m[12]) * w,
        (m[1] * x + m[5] * y + m[9] * z + m[13]) * w,
        (m[2] * x + m[6] * y + m[10] * z + m[14]) * w,
    ]
}

fn write_obj_normal(out: &mut String, n: crate::model::Vec3) {
    writeln!(
        out,
        "vn {} {} {}",
        molstar_rounded_float64(n.x as f64, 100.0),
        molstar_rounded_float64(n.y as f64, 100.0),
        molstar_rounded_float64(n.z as f64, 100.0)
    )
    .expect("writing to String cannot fail");
}

fn write_obj_faces(
    mesh: &Mesh,
    metadata: &ExportMetadata,
    face_start: usize,
    face_end: usize,
    out: &mut String,
    current_group: &mut Option<usize>,
    current_material: &mut Option<MeshMaterial>,
) {
    for face_index in face_start..face_end {
        let f = &mesh.faces[face_index];
        let group = mesh.face_group(face_index);
        if metadata.include_face_groups && *current_group != Some(group) {
            writeln!(out, "g molfig_group_{group}").expect("writing to String cannot fail");
            *current_group = Some(group);
        }
        if let Some(material) = mesh.face_material(face_index) {
            if *current_material != Some(material) {
                out.push_str("usemtl ");
                write_molstar_material_id(out, material);
                out.push('\n');
                *current_material = Some(material);
            }
        }
        let [a, b, c] = molstar_obj_face_indices(f);
        writeln!(out, "f {0}//{0} {1}//{1} {2}//{2}", a, b, c)
            .expect("writing to String cannot fail");
    }
}

fn obj_output_capacity(mesh: &Mesh) -> usize {
    let index_digits = decimal_digits(mesh.vertices.len().max(1));
    64usize
        .saturating_add(mesh.vertices.len().saturating_mul(48))
        .saturating_add(mesh.normals.len().saturating_mul(40))
        .saturating_add(
            mesh.faces
                .len()
                .saturating_mul(16usize.saturating_add(index_digits.saturating_mul(6))),
        )
        .saturating_add(mesh.sections.len().saturating_mul(32))
}

fn decimal_digits(mut value: usize) -> usize {
    let mut digits = 1;
    while value >= 10 {
        value /= 10;
        digits += 1;
    }
    digits
}

fn molstar_obj_sections(mesh: &Mesh) -> Option<&[crate::model::MeshSection]> {
    if mesh.sections.is_empty() {
        return None;
    }
    let mut vertex_cursor = 0;
    let mut face_cursor = 0;
    for section in &mesh.sections {
        if section.vertex_start != vertex_cursor
            || section.vertex_end < section.vertex_start
            || section.vertex_end > mesh.vertices.len()
            || section.vertex_end > mesh.normals.len()
            || section.face_start != face_cursor
            || section.face_end < section.face_start
            || section.face_end > mesh.faces.len()
        {
            return None;
        }
        vertex_cursor = section.vertex_end;
        face_cursor = section.face_end;
    }
    if vertex_cursor == mesh.vertices.len() && face_cursor == mesh.faces.len() {
        Some(&mesh.sections)
    } else {
        None
    }
}

#[cfg(test)]
pub(crate) fn export_mtl(mesh: &Mesh) -> String {
    export_mtl_from_materials(&mesh_materials_in_first_use_order(mesh))
}

pub(crate) fn export_maquette_material_map_json(mesh: &Mesh) -> String {
    let entries = mesh_materials_in_first_use_order(mesh)
        .iter()
        .map(|material| {
            format!(
                "\"{}\":\"#{:06x}\"",
                json_escape(&molstar_material_id(*material)),
                material.color & 0x00ff_ffff
            )
        })
        .collect::<Vec<_>>();
    format!("{{{}}}", entries.join(","))
}

pub(crate) fn export_maquette_material_map_json_from_obj(obj: &[u8]) -> Result<String, String> {
    let text =
        std::str::from_utf8(obj).map_err(|_| "generated OBJ must be UTF-8 text".to_string())?;
    let mut materials = Vec::<(&str, u32)>::new();
    for line in text.lines() {
        let Some(material_id) = line.strip_prefix("usemtl ") else {
            continue;
        };
        let material_id = material_id.trim();
        if materials
            .iter()
            .any(|(existing_id, _)| *existing_id == material_id)
        {
            continue;
        }
        let color = maquette_color_from_molstar_material_id(material_id).ok_or_else(|| {
            format!("generated OBJ contains unsupported material id: {material_id}")
        })?;
        materials.push((material_id, color));
    }
    let entries = materials
        .iter()
        .map(|(material_id, color)| format!("\"{}\":\"#{color:06x}\"", json_escape(material_id)))
        .collect::<Vec<_>>();
    Ok(format!("{{{}}}", entries.join(",")))
}

fn maquette_color_from_molstar_material_id(material_id: &str) -> Option<u32> {
    let hex = material_id.strip_prefix("0x")?;
    let color = hex.get(..6)?;
    (hex.len() > 6)
        .then(|| u32::from_str_radix(color, 16).ok())
        .flatten()
}

pub(crate) fn export_mtl_from_materials(materials: &[MeshMaterial]) -> String {
    let mut out = String::new();
    for material in materials {
        out.push_str("newmtl ");
        out.push_str(&molstar_material_id(*material));
        out.push('\n');
        out.push_str("illum 2\n");
        out.push_str("Ns 163\n");
        out.push_str("Ni 0.001\n");
        out.push_str("Ka 0 0 0\n");
        out.push_str("Kd ");
        out.push_str(&molstar_color_component(material.color >> 16));
        out.push(' ');
        out.push_str(&molstar_color_component(material.color >> 8));
        out.push(' ');
        out.push_str(&molstar_color_component(material.color));
        out.push('\n');
        out.push_str("Ks 0.25 0.25 0.25\n");
        out.push_str("d ");
        out.push_str(&molstar_alpha(*material));
        out.push('\n');
    }
    out
}

fn mesh_materials_in_first_use_order(mesh: &Mesh) -> Vec<MeshMaterial> {
    let mut materials = Vec::new();
    for material in &mesh.face_materials {
        if !materials.contains(material) {
            materials.push(*material);
        }
    }
    materials
}

#[allow(dead_code)]
pub(crate) fn export_ply(mesh: &Mesh) -> String {
    export_ply_with_metadata(mesh, &ExportMetadata::default())
}

pub(crate) fn export_ply_with_metadata(mesh: &Mesh, metadata: &ExportMetadata) -> String {
    let mut out = String::with_capacity(ply_output_capacity(mesh));
    writeln!(out, "ply").expect("writing to String cannot fail");
    writeln!(out, "format ascii 1.0").expect("writing to String cannot fail");
    writeln!(out, "comment Exported by molfig").expect("writing to String cannot fail");
    writeln!(
        out,
        "comment molfig_group_count {}",
        mesh.effective_group_count()
    )
    .expect("writing to String cannot fail");
    if metadata.include_operator_metadata {
        if let Some(metadata_json) = export_metadata_json(metadata) {
            out.push_str("comment molfig_operator_metadata ");
            out.push_str(&metadata_json);
            out.push('\n');
        }
    }
    write!(
        out,
        "comment molfig_face_group_property molfig_group\nelement vertex {}\nproperty float x\nproperty float y\nproperty float z\nelement face {}\nproperty list uchar int vertex_indices\nproperty int molfig_group\nend_header\n",
        mesh.vertices.len(),
        mesh.faces.len()
    )
    .expect("writing to String cannot fail");
    for v in &mesh.vertices {
        writeln!(out, "{:.5} {:.5} {:.5}", v.x, v.y, v.z).expect("writing to String cannot fail");
    }
    for (face_index, f) in mesh.faces.iter().enumerate() {
        writeln!(
            out,
            "3 {} {} {} {}",
            f.a,
            f.b,
            f.c,
            mesh.face_group(face_index)
        )
        .expect("writing to String cannot fail");
    }
    out
}

fn ply_output_capacity(mesh: &Mesh) -> usize {
    let index_digits = decimal_digits(mesh.vertices.len().max(1));
    let group_digits = decimal_digits(mesh.effective_group_count().max(1));
    256usize
        .saturating_add(mesh.vertices.len().saturating_mul(48))
        .saturating_add(
            mesh.faces.len().saturating_mul(
                8usize
                    .saturating_add(index_digits.saturating_mul(3))
                    .saturating_add(group_digits),
            ),
        )
}

#[allow(dead_code)]
pub(crate) fn export_stl(mesh: &Mesh) -> Vec<u8> {
    export_stl_with_metadata(mesh, &ExportMetadata::default())
}

pub(crate) fn export_stl_with_metadata(mesh: &Mesh, metadata: &ExportMetadata) -> Vec<u8> {
    let draw_count = mesh.faces.len().saturating_mul(3);
    let size = 84usize
        .checked_add(
            draw_count
                .checked_mul(50)
                .expect("STL facet record count overflow"),
        )
        .expect("STL output size overflow");
    let mut out = vec![0u8; size];
    let header = b"Exported from Mol* 5.9.0";
    out[..header.len()].copy_from_slice(header);
    let draw_count = u32::try_from(draw_count).expect("STL draw count exceeds u32");
    out[80..84].copy_from_slice(&draw_count.to_le_bytes());
    for (i, face) in mesh.faces.iter().enumerate() {
        let offset = 84 + i * 3 * 50;
        let a = molstar_stl_vertex_transform(mesh.vertices[face.a], metadata.vertex_offset);
        let b = molstar_stl_vertex_transform(mesh.vertices[face.b], metadata.vertex_offset);
        let c = molstar_stl_vertex_transform(mesh.vertices[face.c], metadata.vertex_offset);
        let n = molstar_triangle_normal(a, b, c);
        let values = [n.x, n.y, n.z, a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z];
        for (j, value) in values.iter().enumerate() {
            out[offset + j * 4..offset + j * 4 + 4].copy_from_slice(&value.to_le_bytes());
        }
    }
    out
}

pub(crate) fn export_stl_facet_context_json(
    mesh: &Mesh,
    export_center: ExportVec3,
    stl_facet: usize,
) -> String {
    let stl_facet_count = mesh.faces.len().saturating_mul(3);
    let face_index = stl_facet / 3;
    let sparse_slot = stl_facet % 3;
    let vertex_offset = export_center.negated();
    if stl_facet >= stl_facet_count {
        return format!(
            "{{\"found\":false,\"stl_facet\":{},\"stl_sparse_slot\":{},\"face_index\":{},\"mesh_vertex_count\":{},\"mesh_face_count\":{},\"stl_facet_count\":{},\"export_center\":{},\"vertex_offset\":{}}}",
            stl_facet,
            sparse_slot,
            face_index,
            mesh.vertices.len(),
            mesh.faces.len(),
            stl_facet_count,
            export_vec3_json(export_center),
            export_vec3_json(vertex_offset),
        );
    }

    if sparse_slot != 0 {
        return format!(
            "{{\"found\":true,\"stl_facet\":{},\"stl_sparse_slot\":{},\"face_index\":{},\"mesh_vertex_count\":{},\"mesh_face_count\":{},\"stl_facet_count\":{},\"export_center\":{},\"vertex_offset\":{},\"sparse_slot_has_face\":false,\"stl_normal\":[0.000000000,0.000000000,0.000000000],\"stl_normal_bits\":[\"0x00000000\",\"0x00000000\",\"0x00000000\"],\"stl_vertices\":[[0.000000000,0.000000000,0.000000000],[0.000000000,0.000000000,0.000000000],[0.000000000,0.000000000,0.000000000]],\"stl_vertex_bits\":[[\"0x00000000\",\"0x00000000\",\"0x00000000\"],[\"0x00000000\",\"0x00000000\",\"0x00000000\"],[\"0x00000000\",\"0x00000000\",\"0x00000000\"]],\"target_face\":null}}",
            stl_facet,
            sparse_slot,
            face_index,
            mesh.vertices.len(),
            mesh.faces.len(),
            stl_facet_count,
            export_vec3_json(export_center),
            export_vec3_json(vertex_offset),
        );
    }

    let Some(face) = mesh.faces.get(face_index) else {
        return format!(
            "{{\"found\":false,\"stl_facet\":{},\"stl_sparse_slot\":{},\"face_index\":{},\"mesh_vertex_count\":{},\"mesh_face_count\":{},\"stl_facet_count\":{},\"export_center\":{},\"vertex_offset\":{}}}",
            stl_facet,
            sparse_slot,
            face_index,
            mesh.vertices.len(),
            mesh.faces.len(),
            stl_facet_count,
            export_vec3_json(export_center),
            export_vec3_json(vertex_offset),
        );
    };
    let face_group = mesh.face_group(face_index);
    let section = mesh.sections.iter().find(|section| {
        section.face_start <= face_index
            && face_index < section.face_end
            && section.vertex_start <= face.a
            && face.a < section.vertex_end
    });
    render_export_stl_face_context_json(
        mesh,
        face,
        stl_facet,
        sparse_slot,
        face_index,
        face_group,
        section,
        stl_facet_count,
        export_center,
        vertex_offset,
    )
}

fn molstar_stl_vertex_transform(v: crate::model::Vec3, offset: ExportVec3) -> crate::model::Vec3 {
    let v = molstar_obj_vertex_transform(v, offset);
    crate::model::Vec3::new(v[0] as f32, v[1] as f32, v[2] as f32)
}

fn molstar_triangle_normal(
    a: crate::model::Vec3,
    b: crate::model::Vec3,
    c: crate::model::Vec3,
) -> crate::model::Vec3 {
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
        crate::model::Vec3::new(
            (n[0] * scale) as f32,
            (n[1] * scale) as f32,
            (n[2] * scale) as f32,
        )
    } else {
        crate::model::Vec3::default()
    }
}

#[allow(clippy::too_many_arguments)]
fn render_export_stl_face_context_json(
    mesh: &Mesh,
    face: &crate::model::Face,
    stl_facet: usize,
    sparse_slot: usize,
    face_index: usize,
    face_group: usize,
    section: Option<&crate::model::MeshSection>,
    stl_facet_count: usize,
    export_center: ExportVec3,
    vertex_offset: ExportVec3,
) -> String {
    let indices = [face.a, face.b, face.c];
    let raw_vertices = indices
        .iter()
        .map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(precise_vec3_json)
                .unwrap_or_else(|| "null".to_string())
        })
        .collect::<Vec<_>>()
        .join(",");
    let stl_vertices_vec = indices
        .iter()
        .map(|&index| {
            mesh.vertices
                .get(index)
                .copied()
                .map(|vertex| molstar_stl_vertex_transform(vertex, vertex_offset))
        })
        .collect::<Option<Vec<_>>>();
    let (stl_vertices, stl_vertex_bits, stl_normal, stl_normal_bits) =
        if let Some(vertices) = stl_vertices_vec {
            let normal = molstar_triangle_normal(vertices[0], vertices[1], vertices[2]);
            (
                vertices
                    .iter()
                    .copied()
                    .map(precise_vec3_json)
                    .collect::<Vec<_>>()
                    .join(","),
                vertices
                    .iter()
                    .copied()
                    .map(vec3_bits_json)
                    .collect::<Vec<_>>()
                    .join(","),
                precise_vec3_json(normal),
                vec3_bits_json(normal),
            )
        } else {
            (
                "null,null,null".to_string(),
                "null,null,null".to_string(),
                "null".to_string(),
                "null".to_string(),
            )
        };

    format!(
        "{{\"found\":true,\"stl_facet\":{},\"stl_sparse_slot\":{},\"face_index\":{},\"face_group\":{},\"section\":{},\"mesh_vertex_count\":{},\"mesh_face_count\":{},\"stl_facet_count\":{},\"export_center\":{},\"vertex_offset\":{},\"sparse_slot_has_face\":true,\"stl_normal\":{},\"stl_normal_bits\":{},\"stl_vertices\":[{}],\"stl_vertex_bits\":[{}],\"target_face\":{{\"indices\":[{},{},{}],\"raw_vertices\":[{}]}}}}",
        stl_facet,
        sparse_slot,
        face_index,
        face_group,
        export_section_json(section),
        mesh.vertices.len(),
        mesh.faces.len(),
        stl_facet_count,
        export_vec3_json(export_center),
        export_vec3_json(vertex_offset),
        stl_normal,
        stl_normal_bits,
        stl_vertices,
        stl_vertex_bits,
        indices[0],
        indices[1],
        indices[2],
        raw_vertices,
    )
}

fn export_section_json(section: Option<&crate::model::MeshSection>) -> String {
    if let Some(section) = section {
        format!(
            "{{\"key\":\"{}\",\"vertex_start\":{},\"vertex_end\":{},\"face_start\":{},\"face_end\":{},\"face_offset\":{}}}",
            json_escape(&section.key),
            section.vertex_start,
            section.vertex_end,
            section.face_start,
            section.face_end,
            section.face_end.saturating_sub(section.face_start)
        )
    } else {
        "null".to_string()
    }
}

fn export_vec3_json(value: ExportVec3) -> String {
    format!("[{:.17},{:.17},{:.17}]", value.x, value.y, value.z)
}

fn precise_vec3_json(value: crate::model::Vec3) -> String {
    format!(
        "[{:.9},{:.9},{:.9}]",
        value.x as f64, value.y as f64, value.z as f64
    )
}

fn vec3_bits_json(value: crate::model::Vec3) -> String {
    format!(
        "[\"0x{:08x}\",\"0x{:08x}\",\"0x{:08x}\"]",
        value.x.to_bits(),
        value.y.to_bits(),
        value.z.to_bits()
    )
}

fn molstar_float(value: f32, precision_multiplier: f64) -> String {
    molstar_float64(value as f64, precision_multiplier)
}

fn molstar_float64(value: f64, precision_multiplier: f64) -> String {
    molstar_rounded_float64(value, precision_multiplier).to_string()
}

fn molstar_rounded_float64(value: f64, precision_multiplier: f64) -> f64 {
    let rounded = js_round(value * precision_multiplier) / precision_multiplier;
    if rounded == 0.0 {
        0.0
    } else {
        rounded
    }
}

fn js_round(value: f64) -> f64 {
    (value + 0.5).floor()
}

fn molstar_obj_face_indices(face: &crate::model::Face) -> [usize; 3] {
    [face.a + 1, face.b + 1, face.c + 1]
}

fn molstar_material_id(material: MeshMaterial) -> String {
    let mut out = String::with_capacity(12);
    write_molstar_material_id(&mut out, material);
    out
}

fn write_molstar_material_id(out: &mut String, material: MeshMaterial) {
    write!(out, "0x{:06x}", material.color & 0x00ff_ffff).expect("writing to String cannot fail");
    let alpha_tenths = material.alpha_tenths.min(10);
    if alpha_tenths == 10 {
        out.push('1');
    } else {
        write!(out, "0.{alpha_tenths}").expect("writing to String cannot fail");
    }
}

fn molstar_alpha(material: MeshMaterial) -> String {
    let alpha_tenths = material.alpha_tenths.min(10);
    if alpha_tenths == 10 {
        "1".to_string()
    } else {
        format!("0.{alpha_tenths}")
    }
}

fn molstar_color_component(component_source: u32) -> String {
    let component = (component_source & 0xff) as f32 / 255.0;
    molstar_float(component, 1000.0)
}

fn export_metadata_json(metadata: &ExportMetadata) -> Option<String> {
    if metadata.assembly_id.is_none() && metadata.operators.is_empty() {
        return None;
    }

    let assembly_id = metadata
        .assembly_id
        .as_ref()
        .map(|id| format!("\"{}\"", json_escape(id)))
        .unwrap_or_else(|| "null".to_string());
    let operators = metadata
        .operators
        .iter()
        .map(|operator| {
            format!(
                "{{\"name\":\"{}\",\"instance_id\":\"{}\",\"assembly_id\":\"{}\",\"oper_id\":{},\"oper_list_ids\":{},\"is_identity\":{}}}",
                json_escape(&operator.name),
                json_escape(&operator.instance_id),
                json_escape(&operator.assembly_id),
                operator.oper_id,
                json_string_array(&operator.oper_list_ids),
                if operator.is_identity { "true" } else { "false" }
            )
        })
        .collect::<Vec<_>>()
        .join(",");

    Some(format!(
        "{{\"assembly_id\":{},\"operator_count\":{},\"operators\":[{}]}}",
        assembly_id,
        metadata.operators.len(),
        operators
    ))
}
