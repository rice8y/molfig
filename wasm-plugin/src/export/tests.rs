use super::{
    export_mtl, export_mtl_from_materials, export_obj, export_obj_with_metadata,
    export_stl_with_metadata, molstar_float, molstar_float64, molstar_material_id,
    molstar_obj_vertex_transform, molstar_triangle_normal, ExportMetadata, ExportVec3,
};
use crate::model::{Face, Mesh, MeshMaterial, MeshSection, Vec3};

#[test]
fn molstar_float_matches_javascript_half_rounding() {
    assert_eq!(molstar_float(0.0625, 1000.0), "0.063");
    assert_eq!(molstar_float(-0.0625, 1000.0), "-0.062");
    assert_eq!(molstar_float(0.125, 100.0), "0.13");
    assert_eq!(molstar_float(-0.125, 100.0), "-0.12");
    assert_eq!(molstar_float(-0.005, 100.0), "0");
}

#[test]
fn export_obj_uses_molstar_rounding_for_vertices_and_normals() {
    let mesh = Mesh {
        vertices: vec![
            Vec3::new(-0.0625, 0.0625, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        normals: vec![
            Vec3::new(-0.125, 0.125, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
        ],
        faces: vec![Face { a: 0, b: 1, c: 2 }],
        vertex_groups: vec![0, 0, 0],
        face_groups: vec![0],
        face_materials: Vec::new(),
        sections: Vec::new(),
        group_count: 1,
    };

    let obj = export_obj(&mesh);

    assert!(obj.contains("\nv -0.062 0.063 0\n"));
    assert!(obj.contains("\nvn -0.12 0.13 1\n"));
    assert!(obj.contains("\nf 1//1 2//2 3//3\n"));
}

#[test]
fn export_obj_applies_vertex_offset_before_js_float_formatting() {
    let mesh = Mesh {
        vertices: vec![
            Vec3::new(160.369_13, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        normals: vec![
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
        ],
        faces: vec![Face { a: 0, b: 1, c: 2 }],
        vertex_groups: vec![0, 0, 0],
        face_groups: vec![0],
        face_materials: Vec::new(),
        sections: Vec::new(),
        group_count: 1,
    };
    let metadata = ExportMetadata {
        vertex_offset: ExportVec3::new(366.066_38, 0.0, 0.0),
        ..ExportMetadata::default()
    };

    let obj = export_obj_with_metadata(&mesh, &metadata);

    assert!(obj.contains("\nv 526.436 0 0\n"));
    assert_ne!(
        molstar_float(
            mesh.vertices[0].x + metadata.vertex_offset.to_vec3().x,
            1000.0
        ),
        "526.436"
    );
    assert_eq!(
        molstar_float64(mesh.vertices[0].x as f64 + metadata.vertex_offset.x, 1000.0),
        "526.436"
    );
}

#[test]
fn export_obj_vertex_offset_uses_molstar_center_transform_mat4_order() {
    let vertex = Vec3::new(6.484_824_7, -2.972_500_6, 16.397_552);
    let offset = ExportVec3::new(
        0.847_906_646_900_810_3,
        2.941_000_461_578_369,
        -18.151_618_698_611_856,
    );

    let transformed = molstar_obj_vertex_transform(vertex, offset);

    assert_eq!(transformed[0], 7.332_731_304_340_996);
    assert_eq!(transformed[1], -0.031_500_101_089_477_54);
    assert_eq!(transformed[2], -1.754_066_208_377_480_5);
    assert_eq!(molstar_float64(transformed[1], 1000.0), "-0.032");
}

#[test]
fn export_stl_vertex_offset_uses_molstar_center_transform_before_float32_buffer() {
    let mesh = Mesh {
        vertices: vec![
            Vec3::new(f32::from_bits(0xc47f_ffdf), 0.0, 0.0),
            Vec3::new(f32::from_bits(0xc47f_bfdf), 0.0, 0.0),
            Vec3::new(f32::from_bits(0xc47f_ffdf), 1.0, 0.0),
        ],
        normals: vec![Vec3::new(0.0, 0.0, 1.0); 3],
        faces: vec![Face { a: 0, b: 1, c: 2 }],
        vertex_groups: vec![0, 0, 0],
        face_groups: vec![0],
        face_materials: Vec::new(),
        sections: Vec::new(),
        group_count: 1,
    };
    let metadata = ExportMetadata {
        vertex_offset: ExportVec3::new(366.066_38, 0.0, 0.0),
        ..ExportMetadata::default()
    };

    let stl = export_stl_with_metadata(&mesh, &metadata);

    assert_eq!(stl_f32(&stl, 96).to_bits(), 0xc424_7b9f);
    assert_ne!(stl_f32(&stl, 96).to_bits(), 0xc424_7ba0);
    assert_eq!(stl_f32(&stl, 100), 0.0);
    assert_eq!(stl_f32(&stl, 104), 0.0);
    assert_eq!(stl_f32(&stl, 84), 0.0);
    assert_eq!(stl_f32(&stl, 88), 0.0);
    assert_eq!(stl_f32(&stl, 92), 1.0);
}

#[test]
fn export_stl_facet_normal_uses_molstar_vec3_triangle_normal_staging() {
    let a = Vec3::new(7.332_731_2, -0.031_500_1, -1.754_066_2);
    let b = Vec3::new(7.341_112, -0.016_123_772, -1.745_331_5);
    let c = Vec3::new(7.325_991, -0.018_999_817, -1.760_224_9);
    let normal = molstar_triangle_normal(a, b, c);

    assert_eq!(normal.x.to_bits(), 0xbf32_f852);
    assert_eq!(normal.y.to_bits(), 0xbccb_da6f);
    assert_eq!(normal.z.to_bits(), 0x3f36_ef52);

    let mesh = Mesh {
        vertices: vec![a, b, c],
        normals: vec![Vec3::new(1.0, 0.0, 0.0); 3],
        faces: vec![Face { a: 0, b: 1, c: 2 }],
        vertex_groups: vec![0, 0, 0],
        face_groups: vec![0],
        face_materials: Vec::new(),
        sections: Vec::new(),
        group_count: 1,
    };

    let stl = export_stl_with_metadata(&mesh, &ExportMetadata::default());

    assert_eq!(stl_f32(&stl, 84).to_bits(), normal.x.to_bits());
    assert_eq!(stl_f32(&stl, 88).to_bits(), normal.y.to_bits());
    assert_eq!(stl_f32(&stl, 92).to_bits(), normal.z.to_bits());
}

#[test]
fn export_obj_emits_molstar_material_switches_on_face_color_changes() {
    let mesh = Mesh {
        vertices: vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        normals: vec![
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
        ],
        faces: vec![
            Face { a: 0, b: 1, c: 2 },
            Face { a: 0, b: 2, c: 3 },
            Face { a: 1, b: 2, c: 3 },
        ],
        vertex_groups: vec![0, 0, 0, 0],
        face_groups: vec![0, 0, 0],
        face_materials: vec![
            MeshMaterial::opaque(0x1b9e77),
            MeshMaterial::opaque(0x1b9e77),
            MeshMaterial::with_alpha_tenths(0xff2618, 6),
        ],
        sections: Vec::new(),
        group_count: 1,
    };

    let obj = export_obj(&mesh);

    assert_eq!(
        molstar_material_id(MeshMaterial::opaque(0x1b9e77)),
        "0x1b9e771"
    );
    assert_eq!(
        molstar_material_id(MeshMaterial::with_alpha_tenths(0xff2618, 6)),
        "0xff26180.6"
    );
    let switches = obj
        .lines()
        .filter(|line| line.starts_with("usemtl "))
        .collect::<Vec<_>>();
    assert_eq!(switches, vec!["usemtl 0x1b9e771", "usemtl 0xff26180.6"]);
    assert!(obj.contains("\nusemtl 0x1b9e771\nf 1//1 2//2 3//3\n"));
    assert!(obj.contains("\nf 1//1 3//3 4//4\nusemtl 0xff26180.6\n"));
}

fn stl_f32(stl: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes(stl[offset..offset + 4].try_into().unwrap())
}

#[test]
fn export_obj_uses_molstar_render_object_section_order_when_available() {
    let mesh = Mesh {
        vertices: vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(2.0, 1.0, 0.0),
        ],
        normals: vec![
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ],
        faces: vec![Face { a: 0, b: 1, c: 2 }, Face { a: 3, b: 4, c: 5 }],
        vertex_groups: vec![0, 0, 0, 1, 1, 1],
        face_groups: vec![0, 1],
        face_materials: vec![
            MeshMaterial::opaque(0x1b9e77),
            MeshMaterial::opaque(0xff2618),
        ],
        sections: vec![
            MeshSection {
                key: "first-visual".to_string(),
                vertex_start: 0,
                vertex_end: 3,
                face_start: 0,
                face_end: 1,
            },
            MeshSection {
                key: "second-visual".to_string(),
                vertex_start: 3,
                vertex_end: 6,
                face_start: 1,
                face_end: 2,
            },
        ],
        group_count: 2,
    };

    let obj = export_obj(&mesh);
    let lines = obj.lines().collect::<Vec<_>>();
    assert_eq!(
        lines,
        vec![
            "mtllib molfig.mtl",
            "v 0 0 0",
            "v 1 0 0",
            "v 0 1 0",
            "vn 0 0 1",
            "vn 0 0 1",
            "vn 0 0 1",
            "g molfig_group_0",
            "usemtl 0x1b9e771",
            "f 1//1 2//2 3//3",
            "v 2 0 0",
            "v 3 0 0",
            "v 2 1 0",
            "vn 0 1 0",
            "vn 0 1 0",
            "vn 0 1 0",
            "g molfig_group_1",
            "usemtl 0xff26181",
            "f 4//4 5//5 6//6",
        ]
    );
}

#[test]
fn export_mtl_matches_molstar_material_library_rows() {
    let mtl = export_mtl_from_materials(&[
        MeshMaterial::opaque(0x1b9e77),
        MeshMaterial::with_alpha_tenths(0xff2618, 6),
    ]);

    assert_eq!(
        mtl,
        concat!(
            "newmtl 0x1b9e771\n",
            "illum 2\n",
            "Ns 163\n",
            "Ni 0.001\n",
            "Ka 0 0 0\n",
            "Kd 0.106 0.62 0.467\n",
            "Ks 0.25 0.25 0.25\n",
            "d 1\n",
            "newmtl 0xff26180.6\n",
            "illum 2\n",
            "Ns 163\n",
            "Ni 0.001\n",
            "Ka 0 0 0\n",
            "Kd 1 0.149 0.094\n",
            "Ks 0.25 0.25 0.25\n",
            "d 0.6\n",
        )
    );
}

#[test]
fn export_mtl_uses_first_face_material_order_without_duplicates() {
    let mesh = Mesh {
        vertices: Vec::new(),
        normals: Vec::new(),
        faces: vec![
            Face { a: 0, b: 1, c: 2 },
            Face { a: 0, b: 2, c: 3 },
            Face { a: 1, b: 2, c: 3 },
        ],
        vertex_groups: Vec::new(),
        face_groups: vec![0, 0, 0],
        face_materials: vec![
            MeshMaterial::opaque(0x1b9e77),
            MeshMaterial::opaque(0x1b9e77),
            MeshMaterial::opaque(0xff2618),
        ],
        sections: Vec::new(),
        group_count: 1,
    };

    let mtl = export_mtl(&mesh);
    let materials = mtl
        .lines()
        .filter_map(|line| line.strip_prefix("newmtl "))
        .collect::<Vec<_>>();
    assert_eq!(materials, vec!["0x1b9e771", "0xff26181"]);
}
