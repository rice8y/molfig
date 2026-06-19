use std::f32::consts::PI;
use std::sync::OnceLock;

use crate::model::{Face, Mesh, Vec3};
use crate::options::PolymerProfile;

const MOLSTAR_NUMBER_EPSILON: f32 = f64::EPSILON as f32;
const SPHERE_PRIMITIVE_DETAIL_COUNT: usize = 6;
static SPHERE_PRIMITIVES: [OnceLock<Primitive>; SPHERE_PRIMITIVE_DETAIL_COUNT] =
    [const { OnceLock::new() }; SPHERE_PRIMITIVE_DETAIL_COUNT];
const MAX_TUBE_RADIAL_SEGMENTS: usize = 56;
const TUBE_TRIG_TABLE_COUNT: usize = (MAX_TUBE_RADIAL_SEGMENTS + 1) * 2;
static TUBE_TRIG_TABLES: [OnceLock<TubeTrigTable>; TUBE_TRIG_TABLE_COUNT] =
    [const { OnceLock::new() }; TUBE_TRIG_TABLE_COUNT];

pub(super) fn add_sphere(mesh: &mut Mesh, center: Vec3, radius: f32, detail: usize) {
    add_sphere_with_radius64(
        mesh,
        center,
        molstar_js_number_from_common_f32(radius),
        detail,
    );
}

pub(super) fn add_sphere_with_radius64(mesh: &mut Mesh, center: Vec3, radius: f64, detail: usize) {
    let primitive = molstar_sphere_primitive(detail);
    let normal_scale = (1.0 / radius) as f32;
    let center = [center.x as f64, center.y as f64, center.z as f64];
    let base = mesh.vertices.len();
    for (&vertex, &normal) in primitive.vertices.iter().zip(&primitive.normals) {
        mesh.vertices.push(Vec3::new(
            (radius * vertex.x as f64 + center[0]) as f32,
            (radius * vertex.y as f64 + center[1]) as f32,
            (radius * vertex.z as f64 + center[2]) as f32,
        ));
        mesh.normals.push(normal * normal_scale);
    }
    for face in &primitive.faces {
        mesh.faces.push(Face {
            a: base + face.a,
            b: base + face.b,
            c: base + face.c,
        });
    }
}

pub(super) fn molstar_sphere_triangle_count(detail: usize) -> usize {
    20 * 4usize.pow(detail as u32)
}

pub(super) fn molstar_sphere_mesh_counts(detail: usize) -> (usize, usize) {
    let primitive = molstar_sphere_primitive(detail);
    (primitive.vertices.len(), primitive.faces.len())
}

pub(super) fn molstar_cylinder_mesh_counts(
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
    radius: f64,
) -> (usize, usize) {
    let radial_segments = radial_segments.max(3);
    let cap_count = usize::from(top_cap) + usize::from(bottom_cap);
    if radial_segments <= 4 {
        let (cap_vertices, cap_faces) = if radial_segments == 3 { (3, 1) } else { (4, 2) };
        return (
            radial_segments * 4 + cap_count * cap_vertices,
            radial_segments * 2 + cap_count * cap_faces,
        );
    }

    let cap_count = if radius > 0.0 { cap_count } else { 0 };
    (
        (radial_segments + 1) * 2 + cap_count * (radial_segments * 2 + 1),
        radial_segments * 2 + cap_count * radial_segments,
    )
}

#[derive(Clone, Debug)]
struct Primitive {
    vertices: Vec<Vec3>,
    normals: Vec<Vec3>,
    faces: Vec<Face>,
}

struct TubeTrigTable {
    cos: Vec<f64>,
    sin: Vec<f64>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct CylinderPrimitiveKey {
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
    radius_bits: u64,
}

#[derive(Default)]
pub(super) struct CylinderPrimitiveCache {
    entries: Vec<(CylinderPrimitiveKey, Primitive)>,
}

impl CylinderPrimitiveCache {
    fn get(
        &mut self,
        radial_segments: usize,
        top_cap: bool,
        bottom_cap: bool,
        radius: f64,
    ) -> &Primitive {
        let key = CylinderPrimitiveKey {
            radial_segments: radial_segments.max(3),
            top_cap,
            bottom_cap,
            radius_bits: radius.to_bits(),
        };
        let index = self
            .entries
            .iter()
            .position(|(existing, _)| *existing == key)
            .unwrap_or_else(|| {
                let primitive = build_molstar_cylinder_primitive_with_radius64(
                    key.radial_segments,
                    top_cap,
                    bottom_cap,
                    radius,
                );
                self.entries.push((key, primitive));
                self.entries.len() - 1
            });
        &self.entries[index].1
    }

    #[cfg(test)]
    fn len(&self) -> usize {
        self.entries.len()
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MolstarPrimitiveTransform {
    center: Vec3,
    axes: [Vec3; 3],
}

impl MolstarPrimitiveTransform {
    pub(super) fn from_axes(center: Vec3, x_axis: Vec3, y_axis: Vec3, z_axis: Vec3) -> Self {
        Self {
            center,
            axes: [x_axis, y_axis, z_axis],
        }
    }

    pub(super) fn from_target_to(eye: Vec3, target: Vec3, up: Vec3) -> Self {
        let z = (eye - target).normalized();
        let mut x = up.cross(z).normalized();
        if x.length() <= 0.000_001 {
            x = fallback_side(z, None);
        }
        let y = z.cross(x);
        Self::from_axes(eye, x, y, z)
    }

    pub(super) fn scale(self, scale: Vec3) -> Self {
        Self::from_axes(
            self.center,
            self.axes[0] * scale.x,
            self.axes[1] * scale.y,
            self.axes[2] * scale.z,
        )
    }

    pub(super) fn scale_uniformly(self, scale: f32) -> Self {
        self.scale(Vec3::new(scale, scale, scale))
    }

    pub(super) fn mul_local(self, rhs: MolstarLocalTransform) -> Self {
        let transform_direction = |direction: Vec3| {
            self.axes[0] * direction.x + self.axes[1] * direction.y + self.axes[2] * direction.z
        };
        Self::from_axes(
            self.center,
            transform_direction(rhs.axes[0]),
            transform_direction(rhs.axes[1]),
            transform_direction(rhs.axes[2]),
        )
    }

    fn transform_position(self, point: Vec3) -> Vec3 {
        molstar_transform_position(self.center, self.axes, point)
    }

    fn direction_transform(self) -> [f32; 9] {
        molstar_direction_transform(self.axes)
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MolstarLocalTransform {
    axes: [Vec3; 3],
}

impl MolstarLocalTransform {
    fn from_axes(x_axis: Vec3, y_axis: Vec3, z_axis: Vec3) -> Self {
        Self {
            axes: [x_axis, y_axis, z_axis],
        }
    }

    fn mul(self, rhs: Self) -> Self {
        let transform_direction = |direction: Vec3| {
            self.axes[0] * direction.x + self.axes[1] * direction.y + self.axes[2] * direction.z
        };
        Self::from_axes(
            transform_direction(rhs.axes[0]),
            transform_direction(rhs.axes[1]),
            transform_direction(rhs.axes[2]),
        )
    }

    pub(super) fn rot_y90() -> Self {
        Self::from_axes(
            Vec3::new(0.0, 0.0, -1.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        )
    }

    pub(super) fn rot_z90() -> Self {
        Self::from_axes(
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        )
    }

    fn rot_x180() -> Self {
        Self::from_axes(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, -1.0),
        )
    }

    pub(super) fn rot_zy90() -> Self {
        Self::rot_z90().mul(Self::rot_y90())
    }

    pub(super) fn rot_zyz90() -> Self {
        Self::rot_zy90().mul(Self::rot_z90())
    }

    pub(super) fn rot_z90x180() -> Self {
        Self::rot_z90().mul(Self::rot_x180())
    }
}

fn add_molstar_primitive(
    mesh: &mut Mesh,
    transform: MolstarPrimitiveTransform,
    primitive: &Primitive,
) {
    let base = mesh.vertices.len();
    let normal_transform = transform.direction_transform();
    for (&vertex, &normal) in primitive.vertices.iter().zip(&primitive.normals) {
        mesh.vertices.push(transform.transform_position(vertex));
        mesh.normals
            .push(molstar_transform_direction(normal_transform, normal));
    }
    for face in &primitive.faces {
        mesh.faces.push(Face {
            a: base + face.a,
            b: base + face.b,
            c: base + face.c,
        });
    }
}

fn molstar_sphere_primitive(detail: usize) -> &'static Primitive {
    SPHERE_PRIMITIVES
        .get(detail)
        .unwrap_or_else(|| panic!("sphere detail {detail} exceeds the supported maximum 5"))
        .get_or_init(|| build_molstar_sphere_primitive(detail))
}

fn build_molstar_sphere_primitive(detail: usize) -> Primitive {
    let t = (1.0 + 5.0f64.sqrt()) / 2.0;
    let icosahedron_vertices = [
        DVec3::new(-1.0, t, 0.0),
        DVec3::new(1.0, t, 0.0),
        DVec3::new(-1.0, -t, 0.0),
        DVec3::new(1.0, -t, 0.0),
        DVec3::new(0.0, -1.0, t),
        DVec3::new(0.0, 1.0, t),
        DVec3::new(0.0, -1.0, -t),
        DVec3::new(0.0, 1.0, -t),
        DVec3::new(t, 0.0, -1.0),
        DVec3::new(t, 0.0, 1.0),
        DVec3::new(-t, 0.0, -1.0),
        DVec3::new(-t, 0.0, 1.0),
    ];
    const ICOSAHEDRON_INDICES: [usize; 60] = [
        0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7,
        1, 8, 3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9,
        8, 1,
    ];

    let source_vertices = ICOSAHEDRON_INDICES
        .iter()
        .map(|&index| {
            let vertex = icosahedron_vertices[index];
            DVec3::new(
                vertex.x as f32 as f64,
                vertex.y as f32 as f64,
                vertex.z as f32 as f64,
            )
        })
        .collect::<Vec<_>>();
    let source_indices = (0..source_vertices.len()).collect::<Vec<_>>();

    molstar_polyhedron_primitive_d(&source_vertices, &source_indices, 1.0, detail)
}

#[cfg(test)]
fn molstar_polyhedron_primitive(
    source_vertices: &[Vec3],
    source_indices: &[usize],
    radius: f32,
    detail: usize,
) -> Primitive {
    let mut builder = PolyhedronBuilder::default();
    for face in source_indices.chunks_exact(3) {
        let a = source_vertices[face[0]];
        let b = source_vertices[face[1]];
        let c = source_vertices[face[2]];
        subdivide_molstar_polyhedron_face(&mut builder, a, b, c, detail);
    }
    for vertex in &mut builder.vertices {
        *vertex = vertex.normalized() * radius;
    }
    let normals = molstar_indexed_vertex_normals(&builder.vertices, &builder.faces);
    Primitive {
        vertices: builder.vertices,
        normals,
        faces: builder.faces,
    }
}

fn molstar_polyhedron_primitive_d(
    source_vertices: &[DVec3],
    source_indices: &[usize],
    radius: f64,
    detail: usize,
) -> Primitive {
    let mut builder = PolyhedronBuilderD::default();
    for face in source_indices.chunks_exact(3) {
        let a = source_vertices[face[0]];
        let b = source_vertices[face[1]];
        let c = source_vertices[face[2]];
        subdivide_molstar_polyhedron_face_d(&mut builder, a, b, c, detail);
    }
    let d_vertices = builder
        .vertices
        .iter()
        .map(|vertex| vertex.normalized() * radius)
        .collect::<Vec<_>>();
    let vertices = d_vertices.iter().copied().map(DVec3::to_vec3).collect();
    let normals = molstar_indexed_vertex_normals_d(&d_vertices, &builder.faces);
    Primitive {
        vertices,
        normals,
        faces: builder.faces,
    }
}

#[derive(Default)]
struct MolstarPrimitiveBuilder {
    vertices: Vec<Vec3>,
    normals: Vec<Vec3>,
    faces: Vec<Face>,
}

impl MolstarPrimitiveBuilder {
    fn add(&mut self, a: Vec3, b: Vec3, c: Vec3) {
        let base = self.vertices.len();
        let normal = (b - a).cross(c - a).normalized();
        self.vertices.extend([a, b, c]);
        self.normals.extend([normal; 3]);
        self.faces.push(Face {
            a: base,
            b: base + 1,
            c: base + 2,
        });
    }

    fn add_quad(&mut self, a: Vec3, b: Vec3, c: Vec3, d: Vec3) {
        let base = self.vertices.len();
        let normal = (b - a).cross(c - a).normalized();
        self.vertices.extend([a, b, c, d]);
        self.normals.extend([normal; 4]);
        self.faces.push(Face {
            a: base,
            b: base + 1,
            c: base + 2,
        });
        self.faces.push(Face {
            a: base + 2,
            b: base + 3,
            c: base,
        });
    }

    fn into_primitive(self) -> Primitive {
        Primitive {
            vertices: self.vertices,
            normals: self.normals,
            faces: self.faces,
        }
    }
}

fn molstar_create_primitive(source_vertices: &[Vec3], source_indices: &[usize]) -> Primitive {
    let mut builder = MolstarPrimitiveBuilder::default();
    for face in source_indices.chunks_exact(3) {
        builder.add(
            source_vertices[face[0]],
            source_vertices[face[1]],
            source_vertices[face[2]],
        );
    }
    builder.into_primitive()
}

fn molstar_polygon(side_count: usize, shifted: bool) -> Vec<Vec3> {
    let radius = if side_count <= 4 {
        std::f32::consts::FRAC_1_SQRT_2
    } else {
        0.6
    };
    molstar_polygon_with_radius(side_count, shifted, radius)
}

fn molstar_polygon_with_radius(side_count: usize, shifted: bool, radius: f32) -> Vec<Vec3> {
    let offset = usize::from(shifted);
    (0..side_count)
        .map(|i| {
            let angle = (i * 2 + offset) as f32 / side_count as f32 * PI;
            Vec3::new(angle.cos() * radius, angle.sin() * radius, 0.0)
        })
        .collect()
}

fn molstar_box_primitive(perforated: bool) -> Primitive {
    let points = molstar_polygon(4, true);
    let mut builder = MolstarPrimitiveBuilder::default();

    for i in 0..4 {
        let ni = (i + 1) % 4;
        let a = Vec3::new(points[i].x, points[i].y, -0.5);
        let b = Vec3::new(points[ni].x, points[ni].y, -0.5);
        let c = Vec3::new(points[ni].x, points[ni].y, 0.5);
        let d = Vec3::new(points[i].x, points[i].y, 0.5);
        if perforated {
            builder.add(a, b, c);
        } else {
            builder.add_quad(a, b, c, d);
        }
    }

    let a = Vec3::new(points[0].x, points[0].y, -0.5);
    let b = Vec3::new(points[1].x, points[1].y, -0.5);
    let c = Vec3::new(points[2].x, points[2].y, -0.5);
    let d = Vec3::new(points[3].x, points[3].y, -0.5);
    if perforated {
        builder.add(c, b, a);
    } else {
        builder.add_quad(d, c, b, a);
    }

    let a = Vec3::new(points[0].x, points[0].y, 0.5);
    let b = Vec3::new(points[1].x, points[1].y, 0.5);
    let c = Vec3::new(points[2].x, points[2].y, 0.5);
    let d = Vec3::new(points[3].x, points[3].y, 0.5);
    if perforated {
        builder.add(a, b, c);
    } else {
        builder.add_quad(a, b, c, d);
    }

    builder.into_primitive()
}

fn molstar_wedge_primitive() -> Primitive {
    let points = molstar_polygon(3, false);
    let mut builder = MolstarPrimitiveBuilder::default();

    for i in 0..3 {
        let ni = (i + 1) % 3;
        let a = Vec3::new(points[i].x, points[i].y, -0.5);
        let b = Vec3::new(points[ni].x, points[ni].y, -0.5);
        let c = Vec3::new(points[ni].x, points[ni].y, 0.5);
        let d = Vec3::new(points[i].x, points[i].y, 0.5);
        builder.add(a, b, c);
        builder.add(c, d, a);
    }

    let a = Vec3::new(points[0].x, points[0].y, -0.5);
    let b = Vec3::new(points[1].x, points[1].y, -0.5);
    let c = Vec3::new(points[2].x, points[2].y, -0.5);
    builder.add(c, b, a);

    let a = Vec3::new(points[0].x, points[0].y, 0.5);
    let b = Vec3::new(points[1].x, points[1].y, 0.5);
    let c = Vec3::new(points[2].x, points[2].y, 0.5);
    builder.add(a, b, c);

    builder.into_primitive()
}

fn molstar_prism_primitive(
    side_count: usize,
    shifted: bool,
    top_cap: bool,
    bottom_cap: bool,
) -> Primitive {
    molstar_prism_primitive_with_radius(side_count, shifted, top_cap, bottom_cap, {
        if side_count <= 4 {
            std::f32::consts::FRAC_1_SQRT_2
        } else {
            0.6
        }
    })
}

fn molstar_prism_primitive_with_radius(
    side_count: usize,
    shifted: bool,
    top_cap: bool,
    bottom_cap: bool,
    radius: f32,
) -> Primitive {
    assert!(side_count >= 3, "need at least 3 points to build a prism");
    let points = molstar_polygon_with_radius(side_count, shifted, radius);
    let mut builder = MolstarPrimitiveBuilder::default();
    let half_height = 0.5;

    for i in 0..side_count {
        let ni = (i + 1) % side_count;
        let a = Vec3::new(points[i].x, points[i].y, -half_height);
        let b = Vec3::new(points[ni].x, points[ni].y, -half_height);
        let c = Vec3::new(points[ni].x, points[ni].y, half_height);
        let d = Vec3::new(points[i].x, points[i].y, half_height);
        builder.add_quad(a, b, c, d);
    }

    if side_count == 3 {
        if top_cap {
            let a = Vec3::new(points[0].x, points[0].y, -half_height);
            let b = Vec3::new(points[1].x, points[1].y, -half_height);
            let c = Vec3::new(points[2].x, points[2].y, -half_height);
            builder.add(c, b, a);
        }
        if bottom_cap {
            let a = Vec3::new(points[0].x, points[0].y, half_height);
            let b = Vec3::new(points[1].x, points[1].y, half_height);
            let c = Vec3::new(points[2].x, points[2].y, half_height);
            builder.add(a, b, c);
        }
    } else if side_count == 4 {
        if top_cap {
            let a = Vec3::new(points[0].x, points[0].y, -half_height);
            let b = Vec3::new(points[1].x, points[1].y, -half_height);
            let c = Vec3::new(points[2].x, points[2].y, -half_height);
            let d = Vec3::new(points[3].x, points[3].y, -half_height);
            builder.add_quad(d, c, b, a);
        }
        if bottom_cap {
            let a = Vec3::new(points[0].x, points[0].y, half_height);
            let b = Vec3::new(points[1].x, points[1].y, half_height);
            let c = Vec3::new(points[2].x, points[2].y, half_height);
            let d = Vec3::new(points[3].x, points[3].y, half_height);
            builder.add_quad(a, b, c, d);
        }
    } else {
        let on = Vec3::new(0.0, 0.0, -half_height);
        let op = Vec3::new(0.0, 0.0, half_height);
        for i in 0..side_count {
            let ni = (i + 1) % side_count;
            if top_cap {
                let a = Vec3::new(points[i].x, points[i].y, -half_height);
                let b = Vec3::new(points[ni].x, points[ni].y, -half_height);
                builder.add(on, b, a);
            }
            if bottom_cap {
                let a = Vec3::new(points[i].x, points[i].y, half_height);
                let b = Vec3::new(points[ni].x, points[ni].y, half_height);
                builder.add(a, b, op);
            }
        }
    }

    builder.into_primitive()
}

fn molstar_transform_primitive_rot_x90(mut primitive: Primitive) -> Primitive {
    for vertex in &mut primitive.vertices {
        *vertex = Vec3::new(vertex.x, -vertex.z, vertex.y);
    }
    for normal in &mut primitive.normals {
        *normal = Vec3::new(normal.x, -normal.z, normal.y);
    }
    primitive
}

fn molstar_pyramid_primitive(side_count: usize, shifted: bool) -> Primitive {
    let points = molstar_polygon(side_count, shifted);
    let mut builder = MolstarPrimitiveBuilder::default();
    let on = Vec3::new(0.0, 0.0, -0.5);
    let op = Vec3::new(0.0, 0.0, 0.5);

    for i in 0..side_count {
        let ni = (i + 1) % side_count;
        let a = Vec3::new(points[i].x, points[i].y, -0.5);
        let b = Vec3::new(points[ni].x, points[ni].y, -0.5);
        builder.add(a, b, op);
    }

    if side_count == 3 {
        let a = Vec3::new(points[0].x, points[0].y, -0.5);
        let b = Vec3::new(points[1].x, points[1].y, -0.5);
        let c = Vec3::new(points[2].x, points[2].y, -0.5);
        builder.add(c, b, a);
    } else if side_count == 4 {
        let a = Vec3::new(points[0].x, points[0].y, -0.5);
        let b = Vec3::new(points[1].x, points[1].y, -0.5);
        let c = Vec3::new(points[2].x, points[2].y, -0.5);
        let d = Vec3::new(points[3].x, points[3].y, -0.5);
        builder.add_quad(d, c, b, a);
    } else {
        for i in 0..side_count {
            let ni = (i + 1) % side_count;
            let a = Vec3::new(points[i].x, points[i].y, -0.5);
            let b = Vec3::new(points[ni].x, points[ni].y, -0.5);
            builder.add(on, b, a);
        }
    }

    builder.into_primitive()
}

fn molstar_perforated_octagonal_pyramid_primitive() -> Primitive {
    let points = molstar_polygon(8, true);
    let mut vertices = Vec::with_capacity(10);
    for point in points {
        vertices.push(Vec3::new(point.x, point.y, -0.5));
    }
    vertices.push(Vec3::new(0.0, 0.0, -0.5));
    vertices.push(Vec3::new(0.0, 0.0, 0.5));
    molstar_create_primitive(
        &vertices,
        &[
            0, 1, 8, 1, 2, 8, 4, 5, 8, 5, 6, 8, 2, 3, 9, 3, 4, 9, 6, 7, 9, 7, 0, 9,
        ],
    )
}

fn molstar_star_primitive(
    point_count: usize,
    outer_radius: f32,
    inner_radius: f32,
    thickness: f32,
) -> Primitive {
    let mut builder = MolstarPrimitiveBuilder::default();
    let op = Vec3::new(0.0, 0.0, thickness / 2.0);
    let on = Vec3::new(0.0, 0.0, -thickness / 2.0);

    let mut inner_points = Vec::with_capacity(point_count);
    let mut outer_points = Vec::with_capacity(point_count);
    for i in 0..point_count {
        let co = (i * 2 + 1) as f32 / point_count as f32 * PI;
        let ci = (i * 2 + 2) as f32 / point_count as f32 * PI;
        outer_points.push(Vec3::new(
            co.cos() * outer_radius,
            co.sin() * outer_radius,
            0.0,
        ));
        inner_points.push(Vec3::new(
            ci.cos() * inner_radius,
            ci.sin() * inner_radius,
            0.0,
        ));
    }

    for i in 0..point_count {
        let ni = (i + 1) % point_count;
        let a = outer_points[i];
        let b = inner_points[i];
        let c = outer_points[ni];
        builder.add(op, a, b);
        builder.add(b, a, on);
        builder.add(op, b, c);
        builder.add(c, b, on);
    }

    builder.into_primitive()
}

fn molstar_octahedron_primitive(perforated: bool) -> Primitive {
    let vertices = [
        Vec3::new(0.5, 0.0, 0.0),
        Vec3::new(-0.5, 0.0, 0.0),
        Vec3::new(0.0, 0.5, 0.0),
        Vec3::new(0.0, -0.5, 0.0),
        Vec3::new(0.0, 0.0, 0.5),
        Vec3::new(0.0, 0.0, -0.5),
    ];
    let indices: &[usize] = if perforated {
        &[0, 2, 4, 0, 4, 3, 1, 2, 5, 1, 5, 3]
    } else {
        &[
            0, 2, 4, 0, 4, 3, 0, 3, 5, 0, 5, 2, 1, 2, 5, 1, 5, 3, 1, 3, 4, 1, 4, 2,
        ]
    };
    molstar_create_primitive(&vertices, indices)
}

pub(super) fn add_molstar_box_primitive(
    mesh: &mut Mesh,
    transform: MolstarPrimitiveTransform,
    perforated: bool,
) {
    add_molstar_primitive(mesh, transform, &molstar_box_primitive(perforated));
}

pub(super) fn add_molstar_wedge_primitive(mesh: &mut Mesh, transform: MolstarPrimitiveTransform) {
    add_molstar_primitive(mesh, transform, &molstar_wedge_primitive());
}

pub(super) fn add_molstar_prism_primitive(
    mesh: &mut Mesh,
    transform: MolstarPrimitiveTransform,
    side_count: usize,
    shifted: bool,
) {
    add_molstar_primitive(
        mesh,
        transform,
        &molstar_prism_primitive(side_count, shifted, true, true),
    );
}

pub(super) fn add_molstar_pyramid_primitive(
    mesh: &mut Mesh,
    transform: MolstarPrimitiveTransform,
    side_count: usize,
    shifted: bool,
) {
    add_molstar_primitive(
        mesh,
        transform,
        &molstar_pyramid_primitive(side_count, shifted),
    );
}

pub(super) fn add_molstar_perforated_octagonal_pyramid_primitive(
    mesh: &mut Mesh,
    transform: MolstarPrimitiveTransform,
) {
    add_molstar_primitive(
        mesh,
        transform,
        &molstar_perforated_octagonal_pyramid_primitive(),
    );
}

pub(super) fn add_molstar_star_primitive(mesh: &mut Mesh, transform: MolstarPrimitiveTransform) {
    add_molstar_primitive(mesh, transform, &molstar_star_primitive(5, 1.0, 0.5, 0.5));
}

pub(super) fn add_molstar_octahedron_primitive(
    mesh: &mut Mesh,
    transform: MolstarPrimitiveTransform,
    perforated: bool,
) {
    add_molstar_primitive(mesh, transform, &molstar_octahedron_primitive(perforated));
}

#[derive(Default)]
#[cfg(test)]
struct PolyhedronBuilder {
    vertices: Vec<Vec3>,
    faces: Vec<Face>,
    vertex_keys: Vec<[i32; 3]>,
}

#[derive(Default)]
struct PolyhedronBuilderD {
    vertices: Vec<DVec3>,
    faces: Vec<Face>,
    vertex_keys: Vec<String>,
}

impl PolyhedronBuilderD {
    fn add_triangle(&mut self, a: DVec3, b: DVec3, c: DVec3) {
        let a = self.add_vertex(a);
        let b = self.add_vertex(b);
        let c = self.add_vertex(c);
        self.faces.push(Face { a, b, c });
    }

    fn add_vertex(&mut self, vertex: DVec3) -> usize {
        let key = molstar_polyhedron_vertex_key_d(vertex);
        if let Some(index) = self
            .vertex_keys
            .iter()
            .position(|candidate| *candidate == key)
        {
            index
        } else {
            let index = self.vertices.len();
            self.vertices.push(vertex);
            self.vertex_keys.push(key);
            index
        }
    }
}

#[cfg(test)]
impl PolyhedronBuilder {
    fn add_triangle(&mut self, a: Vec3, b: Vec3, c: Vec3) {
        let a = self.add_vertex(a);
        let b = self.add_vertex(b);
        let c = self.add_vertex(c);
        self.faces.push(Face { a, b, c });
    }

    fn add_vertex(&mut self, vertex: Vec3) -> usize {
        let key = molstar_polyhedron_vertex_key(vertex);
        if let Some(index) = self
            .vertex_keys
            .iter()
            .position(|candidate| *candidate == key)
        {
            index
        } else {
            let index = self.vertices.len();
            self.vertices.push(vertex);
            self.vertex_keys.push(key);
            index
        }
    }
}

#[cfg(test)]
fn subdivide_molstar_polyhedron_face(
    builder: &mut PolyhedronBuilder,
    a: Vec3,
    b: Vec3,
    c: Vec3,
    detail: usize,
) {
    let cols = 1usize << detail;
    let mut vertices = Vec::with_capacity(cols + 1);

    for i in 0..=cols {
        let mut row = Vec::with_capacity(cols - i + 1);
        let aj = lerp_vec3(a, c, i as f32 / cols as f32);
        let bj = lerp_vec3(b, c, i as f32 / cols as f32);
        let rows = cols - i;
        for j in 0..=rows {
            if j == 0 && i == cols {
                row.push(aj);
            } else {
                row.push(lerp_vec3(aj, bj, j as f32 / rows as f32));
            }
        }
        vertices.push(row);
    }

    for i in 0..cols {
        for j in 0..(2 * (cols - i) - 1) {
            let k = j / 2;
            if j % 2 == 0 {
                builder.add_triangle(vertices[i][k + 1], vertices[i + 1][k], vertices[i][k]);
            } else {
                builder.add_triangle(
                    vertices[i][k + 1],
                    vertices[i + 1][k + 1],
                    vertices[i + 1][k],
                );
            }
        }
    }
}

fn subdivide_molstar_polyhedron_face_d(
    builder: &mut PolyhedronBuilderD,
    a: DVec3,
    b: DVec3,
    c: DVec3,
    detail: usize,
) {
    let cols = 1usize << detail;
    let mut vertices = Vec::with_capacity(cols + 1);

    for i in 0..=cols {
        let mut row = Vec::with_capacity(cols - i + 1);
        let aj = lerp_dvec3(a, c, i as f64 / cols as f64);
        let bj = lerp_dvec3(b, c, i as f64 / cols as f64);
        let rows = cols - i;
        for j in 0..=rows {
            if j == 0 && i == cols {
                row.push(aj);
            } else {
                row.push(lerp_dvec3(aj, bj, j as f64 / rows as f64));
            }
        }
        vertices.push(row);
    }

    for i in 0..cols {
        for j in 0..(2 * (cols - i) - 1) {
            let k = j / 2;
            if j % 2 == 0 {
                builder.add_triangle(vertices[i][k + 1], vertices[i + 1][k], vertices[i][k]);
            } else {
                builder.add_triangle(
                    vertices[i][k + 1],
                    vertices[i + 1][k + 1],
                    vertices[i + 1][k],
                );
            }
        }
    }
}

#[cfg(test)]
fn molstar_indexed_vertex_normals(vertices: &[Vec3], faces: &[Face]) -> Vec<Vec3> {
    let mut normals = vec![Vec3::default(); vertices.len()];
    for face in faces {
        let a = vertices[face.a];
        let b = vertices[face.b];
        let c = vertices[face.c];
        let normal = (c - b).cross(a - b);
        normals[face.a] = normals[face.a] + normal;
        normals[face.b] = normals[face.b] + normal;
        normals[face.c] = normals[face.c] + normal;
    }
    normals.into_iter().map(Vec3::normalized).collect()
}

fn molstar_indexed_vertex_normals_d(vertices: &[DVec3], faces: &[Face]) -> Vec<Vec3> {
    let mut normals = vec![Vec3::default(); vertices.len()];
    for face in faces {
        let a = vertices[face.a];
        let b = vertices[face.b];
        let c = vertices[face.c];
        let normal = (c - b).cross(a - b);
        for index in [face.a, face.b, face.c] {
            normals[index].x = (normals[index].x as f64 + normal.x) as f32;
            normals[index].y = (normals[index].y as f64 + normal.y) as f32;
            normals[index].z = (normals[index].z as f64 + normal.z) as f32;
        }
    }
    normals
        .into_iter()
        .map(|normal| {
            let x = normal.x as f64;
            let y = normal.y as f64;
            let z = normal.z as f64;
            let len_sq = x * x + y * y + z * z;
            if len_sq > 0.0 {
                let scale = 1.0 / len_sq.sqrt();
                Vec3::new((x * scale) as f32, (y * scale) as f32, (z * scale) as f32)
            } else {
                Vec3::default()
            }
        })
        .collect()
}

#[cfg(test)]
fn molstar_polyhedron_vertex_key(vertex: Vec3) -> [i32; 3] {
    [
        (vertex.x * 100_000.0).round() as i32,
        (vertex.y * 100_000.0).round() as i32,
        (vertex.z * 100_000.0).round() as i32,
    ]
}

fn molstar_polyhedron_vertex_key_d(vertex: DVec3) -> String {
    format!(
        "{}|{}|{}",
        molstar_to_fixed_5(vertex.x),
        molstar_to_fixed_5(vertex.y),
        molstar_to_fixed_5(vertex.z)
    )
}

fn molstar_to_fixed_5(value: f64) -> String {
    if value == 0.0 {
        "0.00000".to_string()
    } else {
        format!("{value:.5}")
    }
}

#[cfg(test)]
fn lerp_vec3(a: Vec3, b: Vec3, t: f32) -> Vec3 {
    a * (1.0 - t) + b * t
}

fn lerp_dvec3(a: DVec3, b: DVec3, t: f64) -> DVec3 {
    a * (1.0 - t) + b * t
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f32 = 0.000_01;

    #[test]
    fn sphere_primitive_detail_zero_matches_molstar_icosahedron_order() {
        // Derived from artifacts/molstar/src/mol-geo/primitive/sphere.ts and
        // polyhedron.ts at the pinned Mol* commit.
        let mut mesh = Mesh::default();
        add_sphere(&mut mesh, Vec3::new(1.0, 2.0, 3.0), 2.0, 0);

        assert_eq!(mesh.vertices.len(), 12);
        assert_eq!(mesh.normals.len(), 12);
        assert_eq!(mesh.faces.len(), 20);

        let expected_unit_vertices = [
            Vec3::new(-0.850_650_8, 0.0, 0.525_731_1),
            Vec3::new(0.0, 0.525_731_1, 0.850_650_8),
            Vec3::new(-0.525_731_1, 0.850_650_8, 0.0),
            Vec3::new(0.525_731_1, 0.850_650_8, 0.0),
            Vec3::new(0.0, 0.525_731_1, -0.850_650_8),
            Vec3::new(-0.850_650_8, 0.0, -0.525_731_1),
            Vec3::new(0.850_650_8, 0.0, 0.525_731_1),
            Vec3::new(0.0, -0.525_731_1, 0.850_650_8),
            Vec3::new(-0.525_731_1, -0.850_650_8, 0.0),
            Vec3::new(0.0, -0.525_731_1, -0.850_650_8),
            Vec3::new(0.850_650_8, 0.0, -0.525_731_1),
            Vec3::new(0.525_731_1, -0.850_650_8, 0.0),
        ];
        for (i, expected) in expected_unit_vertices.iter().enumerate() {
            assert_vec3_close(mesh.vertices[i], Vec3::new(1.0, 2.0, 3.0) + *expected * 2.0);
            assert_vec3_close(mesh.normals[i], *expected * 0.5);
        }
        assert_faces_eq(
            &mesh.faces,
            &[
                [0, 1, 2],
                [1, 3, 2],
                [3, 4, 2],
                [4, 5, 2],
                [5, 0, 2],
                [1, 6, 3],
                [0, 7, 1],
                [5, 8, 0],
                [4, 9, 5],
                [3, 10, 4],
                [6, 7, 11],
                [7, 8, 11],
                [8, 9, 11],
                [9, 10, 11],
                [10, 6, 11],
                [6, 1, 7],
                [7, 0, 8],
                [8, 5, 9],
                [9, 4, 10],
                [10, 3, 6],
            ],
        );
    }

    #[test]
    fn sphere_primitive_counts_follow_molstar_subdivided_icosahedron() {
        for detail in 0..=3 {
            let mut mesh = Mesh::default();
            add_sphere(&mut mesh, Vec3::default(), 1.0, detail);

            assert_eq!(mesh.vertices.len(), 10 * 4usize.pow(detail as u32) + 2);
            assert_eq!(mesh.normals.len(), mesh.vertices.len());
            assert_eq!(mesh.faces.len(), molstar_sphere_triangle_count(detail));
        }
    }

    #[test]
    fn sphere_primitive_cache_matches_uncached_builder_for_all_supported_details() {
        for detail in 0..SPHERE_PRIMITIVE_DETAIL_COUNT {
            let cached = molstar_sphere_primitive(detail);
            assert!(std::ptr::eq(cached, molstar_sphere_primitive(detail)));

            let uncached = build_molstar_sphere_primitive(detail);
            assert_eq!(cached.vertices.len(), uncached.vertices.len());
            assert_eq!(cached.normals.len(), uncached.normals.len());
            assert_eq!(cached.faces.len(), uncached.faces.len());
            for (index, (actual, expected)) in
                cached.vertices.iter().zip(&uncached.vertices).enumerate()
            {
                assert_eq!(actual.x.to_bits(), expected.x.to_bits(), "vertex {index} x");
                assert_eq!(actual.y.to_bits(), expected.y.to_bits(), "vertex {index} y");
                assert_eq!(actual.z.to_bits(), expected.z.to_bits(), "vertex {index} z");
            }
            for (index, (actual, expected)) in
                cached.normals.iter().zip(&uncached.normals).enumerate()
            {
                assert_eq!(actual.x.to_bits(), expected.x.to_bits(), "normal {index} x");
                assert_eq!(actual.y.to_bits(), expected.y.to_bits(), "normal {index} y");
                assert_eq!(actual.z.to_bits(), expected.z.to_bits(), "normal {index} z");
            }
            for (index, (actual, expected)) in cached.faces.iter().zip(&uncached.faces).enumerate()
            {
                assert_eq!(
                    [actual.a, actual.b, actual.c],
                    [expected.a, expected.b, expected.c],
                    "face {index}"
                );
            }
        }
    }

    #[test]
    fn sphere_primitive_detail_three_matches_molstar_reference_vertices() {
        let primitive = molstar_sphere_primitive(3);
        assert_faces_eq(&[primitive.faces[113]], &[[67, 72, 71]]);
        let expected = [
            (
                67,
                Vec3::new(0.080_178_73, 0.988_302_05, 0.129_731_91),
                Vec3::new(0.078_267_95, 0.987_694_14, 0.135_404_74),
            ),
            (
                72,
                Vec3::new(0.234_579_70, 0.963_828_44, 0.126_519_31),
                Vec3::new(0.229_469_79, 0.964_327_16, 0.131_972_58),
            ),
            (
                71,
                Vec3::new(0.152_696_59, 0.988_273_14, 0.0),
                Vec3::new(0.151_797_31, 0.988_411_66, 3.339_272_7e-10),
            ),
        ];
        for (index, vertex, normal) in expected {
            assert_eq!(
                primitive.vertices[index].x.to_bits(),
                vertex.x.to_bits(),
                "vertex {index} x"
            );
            assert_eq!(
                primitive.vertices[index].y.to_bits(),
                vertex.y.to_bits(),
                "vertex {index} y"
            );
            assert_eq!(
                primitive.vertices[index].z.to_bits(),
                vertex.z.to_bits(),
                "vertex {index} z"
            );
            assert_eq!(
                primitive.normals[index].x.to_bits(),
                normal.x.to_bits(),
                "normal {index} x"
            );
            assert_eq!(
                primitive.normals[index].y.to_bits(),
                normal.y.to_bits(),
                "normal {index} y"
            );
            assert_eq!(
                primitive.normals[index].z.to_bits(),
                normal.z.to_bits(),
                "normal {index} z"
            );
        }
    }

    #[test]
    fn ellipsoid_primitive_counts_follow_molstar_sphere_detail() {
        for detail in 0..=3 {
            let mut mesh = Mesh::default();
            add_ellipsoid(
                &mut mesh,
                Vec3::default(),
                [
                    Vec3::new(2.0, 0.0, 0.0),
                    Vec3::new(0.0, 3.0, 0.0),
                    Vec3::new(0.0, 0.0, 4.0),
                ],
                detail,
            );

            assert_eq!(mesh.vertices.len(), 10 * 4usize.pow(detail as u32) + 2);
            assert_eq!(mesh.normals.len(), mesh.vertices.len());
            assert_eq!(mesh.faces.len(), molstar_sphere_triangle_count(detail));
        }
    }

    #[test]
    fn ellipsoid_primitive_uses_molstar_sphere_transform() {
        // Equivalent to Mol* setEllipsoidMat(center, +X major, +Y minor,
        // Vec3(2, 3, 4)): targetTo gives scaled axes +Z, +Y, -X.
        let mut mesh = Mesh::default();
        add_ellipsoid(
            &mut mesh,
            Vec3::new(1.0, 2.0, 3.0),
            [
                Vec3::new(0.0, 0.0, 2.0),
                Vec3::new(0.0, 3.0, 0.0),
                Vec3::new(-4.0, 0.0, 0.0),
            ],
            0,
        );

        assert_eq!(mesh.vertices.len(), 12);
        assert_eq!(mesh.normals.len(), 12);
        assert_eq!(mesh.faces.len(), 20);
        assert_eq!(mesh.faces[0].a, 0);
        assert_eq!(mesh.faces[0].b, 1);
        assert_eq!(mesh.faces[0].c, 2);

        assert_vec3_close(mesh.vertices[0], Vec3::new(-1.102_924, 2.0, 1.298_698));
        assert_vec3_close(mesh.vertices[1], Vec3::new(-2.402_604, 3.577_193, 3.0));
        assert_vec3_close(mesh.vertices[2], Vec3::new(1.0, 4.551_953, 1.948_538));
        assert_vec3_close(mesh.normals[0], Vec3::new(-0.131_433, 0.0, -0.425_325));
        assert_vec3_close(mesh.normals[1], Vec3::new(-0.212_663, 0.175_244, 0.0));
    }

    #[test]
    fn cylinder_primitive_matches_molstar_vertex_and_cap_order() {
        // Derived from artifacts/molstar/src/mol-geo/primitive/cylinder.ts at
        // the pinned Mol* commit. Mol* keeps a duplicate seam column and
        // duplicates cap centers/surrounding vertices per segment.
        let primitive = molstar_cylinder_primitive(8, true, true);

        assert_eq!(primitive.vertices.len(), 52);
        assert_eq!(primitive.normals.len(), 52);
        assert_eq!(primitive.faces.len(), 32);

        assert_vec3_close(primitive.vertices[0], Vec3::new(0.0, 0.5, 1.0));
        assert_vec3_close(
            primitive.vertices[1],
            Vec3::new(0.707_106_77, 0.5, 0.707_106_77),
        );
        assert_vec3_close(primitive.vertices[8], Vec3::new(0.0, 0.5, 1.0));
        assert_vec3_close(primitive.vertices[9], Vec3::new(0.0, -0.5, 1.0));
        assert_vec3_close(
            primitive.normals[1],
            Vec3::new(0.707_106_77, 0.0, 0.707_106_77),
        );

        assert_eq!(primitive.faces[0].a, 0);
        assert_eq!(primitive.faces[0].b, 9);
        assert_eq!(primitive.faces[0].c, 1);
        assert_eq!(primitive.faces[1].a, 9);
        assert_eq!(primitive.faces[1].b, 10);
        assert_eq!(primitive.faces[1].c, 1);

        assert_vec3_close(primitive.vertices[18], Vec3::new(0.0, 0.5, 0.0));
        assert_vec3_close(primitive.normals[18], Vec3::new(0.0, 1.0, 0.0));
        assert_eq!(primitive.faces[16].a, 26);
        assert_eq!(primitive.faces[16].b, 27);
        assert_eq!(primitive.faces[16].c, 18);

        assert_vec3_close(primitive.vertices[35], Vec3::new(0.0, -0.5, 0.0));
        assert_vec3_close(primitive.normals[35], Vec3::new(0.0, -1.0, 0.0));
        assert_eq!(primitive.faces[24].a, 44);
        assert_eq!(primitive.faces[24].b, 43);
        assert_eq!(primitive.faces[24].c, 35);
    }

    #[test]
    fn cylinder_primitive_segments_three_matches_molstar_golden_order() {
        // Small full-order golden derived from
        // artifacts/molstar/src/mol-geo/primitive/cylinder.ts.
        let primitive = molstar_cylinder_primitive(3, true, true);

        assert_vec3s_close(
            &primitive.vertices,
            &[
                Vec3::new(0.0, 0.5, 1.0),
                Vec3::new(0.866_025_4, 0.5, -0.5),
                Vec3::new(-0.866_025_4, 0.5, -0.5),
                Vec3::new(0.0, 0.5, 1.0),
                Vec3::new(0.0, -0.5, 1.0),
                Vec3::new(0.866_025_4, -0.5, -0.5),
                Vec3::new(-0.866_025_4, -0.5, -0.5),
                Vec3::new(0.0, -0.5, 1.0),
                Vec3::new(0.0, 0.5, 0.0),
                Vec3::new(0.0, 0.5, 0.0),
                Vec3::new(0.0, 0.5, 0.0),
                Vec3::new(0.0, 0.5, 1.0),
                Vec3::new(0.866_025_4, 0.5, -0.5),
                Vec3::new(-0.866_025_4, 0.5, -0.5),
                Vec3::new(0.0, 0.5, 1.0),
                Vec3::new(0.0, -0.5, 0.0),
                Vec3::new(0.0, -0.5, 0.0),
                Vec3::new(0.0, -0.5, 0.0),
                Vec3::new(0.0, -0.5, 1.0),
                Vec3::new(0.866_025_4, -0.5, -0.5),
                Vec3::new(-0.866_025_4, -0.5, -0.5),
                Vec3::new(0.0, -0.5, 1.0),
            ],
        );
        assert_faces_eq(
            &primitive.faces,
            &[
                [0, 4, 1],
                [4, 5, 1],
                [1, 5, 2],
                [5, 6, 2],
                [2, 6, 3],
                [6, 7, 3],
                [11, 12, 8],
                [12, 13, 9],
                [13, 14, 10],
                [19, 18, 15],
                [20, 19, 16],
                [21, 20, 17],
            ],
        );
    }

    #[test]
    fn cylinder_primitive_cache_reuses_exact_radius_and_cap_keys() {
        let mut cache = CylinderPrimitiveCache::default();
        let cases = [
            (3, true, true, 0.2_f64),
            (4, false, true, 0.2_f64),
            (8, false, false, 0.2_f64),
            (12, true, false, 0.125_f64),
            (36, true, true, 0.1_f32 as f64),
        ];

        for (segments, top_cap, bottom_cap, radius) in cases {
            let cached_ptr = cache.get(segments, top_cap, bottom_cap, radius) as *const Primitive;
            let cached_again_ptr =
                cache.get(segments, top_cap, bottom_cap, radius) as *const Primitive;
            assert_eq!(cached_ptr, cached_again_ptr);

            let cached = cache.get(segments, top_cap, bottom_cap, radius);
            let uncached = build_molstar_cylinder_primitive_with_radius64(
                segments, top_cap, bottom_cap, radius,
            );
            assert_eq!(cached.vertices.len(), uncached.vertices.len());
            assert_eq!(cached.normals.len(), uncached.normals.len());
            assert_eq!(cached.faces.len(), uncached.faces.len());
            for (index, (actual, expected)) in
                cached.vertices.iter().zip(&uncached.vertices).enumerate()
            {
                assert_eq!(actual.x.to_bits(), expected.x.to_bits(), "vertex {index} x");
                assert_eq!(actual.y.to_bits(), expected.y.to_bits(), "vertex {index} y");
                assert_eq!(actual.z.to_bits(), expected.z.to_bits(), "vertex {index} z");
            }
            for (index, (actual, expected)) in
                cached.normals.iter().zip(&uncached.normals).enumerate()
            {
                assert_eq!(actual.x.to_bits(), expected.x.to_bits(), "normal {index} x");
                assert_eq!(actual.y.to_bits(), expected.y.to_bits(), "normal {index} y");
                assert_eq!(actual.z.to_bits(), expected.z.to_bits(), "normal {index} z");
            }
            for (index, (actual, expected)) in cached.faces.iter().zip(&uncached.faces).enumerate()
            {
                assert_eq!(
                    [actual.a, actual.b, actual.c],
                    [expected.a, expected.b, expected.c],
                    "face {index}"
                );
            }
        }

        assert_eq!(cache.len(), cases.len());
        cache.get(8, false, false, f64::from_bits(0.2_f64.to_bits() + 1));
        assert_eq!(cache.len(), cases.len() + 1);
    }

    #[test]
    fn cylinder_primitive_props_support_taper_and_height_segments_like_molstar() {
        let primitive = molstar_cylinder_primitive_with_props(CylinderPrimitiveProps {
            radial_segments: 8,
            height_segments: 2,
            radius_top: 0.5,
            radius_bottom: 1.0,
            height: 2.0,
            top_cap: true,
            bottom_cap: true,
            theta_start: 0.0,
            theta_length: std::f64::consts::PI * 2.0,
        });

        assert_eq!(primitive.vertices.len(), 61);
        assert_eq!(primitive.normals.len(), 61);
        assert_eq!(primitive.faces.len(), 48);
        assert_vec3_close(primitive.vertices[0], Vec3::new(0.0, 1.0, 0.5));
        assert_vec3_close(primitive.vertices[9], Vec3::new(0.0, 0.0, 0.75));
        assert_vec3_close(primitive.vertices[18], Vec3::new(0.0, -1.0, 1.0));
        assert_vec3_close(primitive.normals[0], Vec3::new(0.0, 0.242_536, 0.970_143));
        assert_eq!(
            (
                primitive.faces[0].a,
                primitive.faces[0].b,
                primitive.faces[0].c
            ),
            (0, 9, 1)
        );
        assert_eq!(
            (
                primitive.faces[2].a,
                primitive.faces[2].b,
                primitive.faces[2].c
            ),
            (9, 18, 10)
        );
        assert_vec3_close(primitive.vertices[27], Vec3::new(0.0, 1.0, 0.0));
        assert_vec3_close(primitive.vertices[35], Vec3::new(0.0, 1.0, 0.5));
        assert_vec3_close(primitive.vertices[44], Vec3::new(0.0, -1.0, 0.0));
        assert_vec3_close(primitive.vertices[52], Vec3::new(0.0, -1.0, 1.0));
    }

    #[test]
    fn add_cylinder_uses_molstar_primitive_transform() {
        let mut mesh = Mesh::default();
        add_cylinder(
            &mut mesh,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            0.5,
            8,
        );

        assert_eq!(mesh.vertices.len(), 52);
        assert_eq!(mesh.normals.len(), 52);
        assert_eq!(mesh.faces.len(), 32);
        assert_vec3_close(mesh.vertices[0], Vec3::new(0.0, 2.0, 0.5));
        assert_vec3_close(mesh.vertices[9], Vec3::new(0.0, 0.0, 0.5));
        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, 0.0, 1.0));
        assert_vec3_close(mesh.normals[18], Vec3::new(0.0, 0.5, 0.0));
        assert_vec3_close(mesh.normals[35], Vec3::new(0.0, -0.5, 0.0));
    }

    #[test]
    fn fixed_count_dashed_cylinder_matches_molstar_gap_spacing() {
        let mut mesh = Mesh::default();
        let primitive = molstar_get_cylinder_primitive(8, true, true);

        add_fixed_count_dashed_cylinder(
            &mut mesh,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 21.0, 0.0),
            0.25,
            8,
            0.5,
            10,
        );

        assert_eq!(mesh.vertices.len(), primitive.vertices.len() * 5);
        assert_eq!(mesh.normals.len(), primitive.normals.len() * 5);
        assert_eq!(mesh.faces.len(), primitive.faces.len() * 5);
        for dash in 0..5 {
            let offset = dash * primitive.vertices.len();
            let start_y = 1.0 + dash as f32 * 2.0;
            let end_y = start_y + 1.0;
            assert_vec3_close(mesh.vertices[offset], Vec3::new(0.0, end_y, 0.25));
            assert_vec3_close(mesh.vertices[offset + 9], Vec3::new(0.0, start_y, 0.25));
        }
    }

    #[test]
    fn fixed_count_dashed_cylinder_halves_odd_final_dash_like_molstar() {
        let mut mesh = Mesh::default();
        let full = molstar_get_cylinder_primitive(8, true, true);
        let final_without_top = molstar_get_cylinder_primitive(8, false, true);

        add_fixed_count_dashed_cylinder(
            &mut mesh,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 11.0, 0.0),
            0.25,
            8,
            1.0,
            5,
        );

        assert_eq!(
            mesh.vertices.len(),
            full.vertices.len() * 2 + final_without_top.vertices.len()
        );
        assert_eq!(
            mesh.faces.len(),
            full.faces.len() * 2 + final_without_top.faces.len()
        );
        let (min_y, max_y) = mesh
            .vertices
            .iter()
            .fold((f32::INFINITY, f32::NEG_INFINITY), |(min_y, max_y), v| {
                (min_y.min(v.y), max_y.max(v.y))
            });
        assert!((min_y - 2.0).abs() < 0.000_01);
        assert!((max_y - 11.0).abs() < 0.000_01);
    }

    #[test]
    fn fixed_count_dashed_cylinder_stub_cap_keeps_odd_final_top_cap_like_molstar() {
        let mut mesh = Mesh::default();
        let full = molstar_get_cylinder_primitive(8, true, true);

        add_fixed_count_dashed_cylinder_with_stub_cap(
            &mut mesh,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 11.0, 0.0),
            0.25,
            8,
            1.0,
            5,
            true,
        );

        assert_eq!(mesh.vertices.len(), full.vertices.len() * 3);
        assert_eq!(mesh.faces.len(), full.faces.len() * 3);
        let (min_y, max_y) = mesh
            .vertices
            .iter()
            .fold((f32::INFINITY, f32::NEG_INFINITY), |(min_y, max_y), v| {
                (min_y.min(v.y), max_y.max(v.y))
            });
        assert!((min_y - 2.0).abs() < 0.000_01);
        assert!((max_y - 11.0).abs() < 0.000_01);
    }

    #[test]
    fn fixed_count_dashed_cylinder_uses_molstar_non_matching_direction_rotation() {
        let mut mesh = Mesh::default();

        add_fixed_count_dashed_cylinder(
            &mut mesh,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, -2.0, 0.0),
            0.25,
            8,
            1.0,
            2,
        );

        assert_vec3_close(mesh.vertices[0], Vec3::new(0.0, -1.6, -0.25));
    }

    #[test]
    fn meshbuilder_cylinder_uses_molstar_prism_path_for_low_radial_segments() {
        let primitive = molstar_get_cylinder_primitive(3, true, true);
        assert_eq!(primitive.vertices.len(), 18);
        assert_eq!(primitive.faces.len(), 8);
        assert_vec3_close(primitive.vertices[0], Vec3::new(0.5, 0.5, 0.866_025_4));

        let mut mesh = Mesh::default();
        add_cylinder(
            &mut mesh,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            1.0,
            3,
        );

        assert_eq!(mesh.vertices.len(), primitive.vertices.len());
        assert_eq!(mesh.faces.len(), primitive.faces.len());
        assert_vec3_close(mesh.vertices[0], Vec3::new(0.5, 1.0, 0.866_025_4));
    }

    #[test]
    fn sphere_meshbuilder_appends_vertex_normal_and_index_buffers_in_molstar_order() {
        // Mol* MeshBuilder.addPrimitive appends all transformed positions and
        // inverse-direction transformed normals in primitive vertex order, then
        // appends primitive indices with the starting vertex offset.
        // See artifacts/molstar/src/mol-geo/geometry/mesh/mesh-builder.ts.
        let mut mesh = Mesh {
            vertices: vec![Vec3::new(-1.0, 0.0, 0.0)],
            normals: vec![Vec3::new(0.0, -1.0, 0.0)],
            faces: vec![Face { a: 0, b: 0, c: 0 }],
            vertex_groups: Vec::new(),
            face_groups: Vec::new(),
            face_materials: Vec::new(),
            sections: Vec::new(),
            group_count: 0,
        };
        let vertex_offset = mesh.vertices.len();
        let face_offset = mesh.faces.len();
        let center = Vec3::new(1.0, 2.0, 3.0);
        let radius = 2.0;
        let primitive = molstar_sphere_primitive(0);

        add_sphere(&mut mesh, center, radius, 0);

        assert_eq!(
            mesh.vertices.len(),
            vertex_offset + primitive.vertices.len()
        );
        assert_eq!(mesh.normals.len(), vertex_offset + primitive.normals.len());
        for (i, (&vertex, &normal)) in primitive
            .vertices
            .iter()
            .zip(&primitive.normals)
            .enumerate()
        {
            assert_vec3_close(mesh.vertices[vertex_offset + i], center + vertex * radius);
            assert_vec3_close(mesh.normals[vertex_offset + i], normal * (1.0 / radius));
        }
        assert_eq!(mesh.faces.len(), face_offset + primitive.faces.len());
        for (i, expected) in primitive.faces.iter().enumerate() {
            let actual = mesh.faces[face_offset + i];
            assert_eq!(
                (actual.a, actual.b, actual.c),
                (
                    vertex_offset + expected.a,
                    vertex_offset + expected.b,
                    vertex_offset + expected.c,
                ),
                "primitive face {i}"
            );
        }
    }

    #[test]
    fn cylinder_meshbuilder_appends_vertex_normal_and_index_buffers_in_molstar_order() {
        // Mol* MeshBuilder.addPrimitive appends transformed primitive
        // vertices/normals first, then primitive indices as primitive_index +
        // current vertex offset, preserving primitive order.
        // See artifacts/molstar/src/mol-geo/geometry/mesh/mesh-builder.ts.
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(-1.0, 0.0, 0.0),
                Vec3::new(0.0, -1.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ],
            normals: vec![Vec3::new(0.0, 0.0, 1.0); 4],
            faces: vec![Face { a: 0, b: 1, c: 2 }],
            vertex_groups: Vec::new(),
            face_groups: Vec::new(),
            face_materials: Vec::new(),
            sections: Vec::new(),
            group_count: 0,
        };
        let vertex_offset = mesh.vertices.len();
        let face_offset = mesh.faces.len();
        let primitive = molstar_get_cylinder_primitive(3, true, true);

        add_cylinder(
            &mut mesh,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            1.0,
            3,
        );

        assert_eq!(
            mesh.vertices.len(),
            vertex_offset + primitive.vertices.len()
        );
        assert_eq!(mesh.normals.len(), vertex_offset + primitive.normals.len());
        for (i, (&vertex, &normal)) in primitive
            .vertices
            .iter()
            .zip(&primitive.normals)
            .enumerate()
        {
            assert_vec3_close(
                mesh.vertices[vertex_offset + i],
                vertex + Vec3::new(0.0, 0.5, 0.0),
            );
            assert_vec3_close(mesh.normals[vertex_offset + i], normal);
        }
        assert_eq!(mesh.faces.len(), face_offset + primitive.faces.len());
        for (i, expected) in primitive.faces.iter().enumerate() {
            let actual = mesh.faces[face_offset + i];
            assert_eq!(
                (actual.a, actual.b, actual.c),
                (
                    vertex_offset + expected.a,
                    vertex_offset + expected.b,
                    vertex_offset + expected.c,
                ),
                "primitive face {i}"
            );
        }
    }

    #[test]
    fn prism_primitive_caps_follow_molstar_polygon_cap_rules() {
        // Derived from artifacts/molstar/src/mol-geo/primitive/prism.ts:
        // pentagonal and larger caps use one center triangle per side.
        let primitive = molstar_prism_primitive(5, false, true, true);

        assert_eq!(primitive.vertices.len(), 50);
        assert_eq!(primitive.normals.len(), 50);
        assert_eq!(primitive.faces.len(), 20);
        assert_eq!(
            (
                primitive.faces[10].a,
                primitive.faces[10].b,
                primitive.faces[10].c
            ),
            (20, 21, 22)
        );
        assert_eq!(
            (
                primitive.faces[11].a,
                primitive.faces[11].b,
                primitive.faces[11].c
            ),
            (23, 24, 25)
        );
        assert_vec3_close(primitive.vertices[20], Vec3::new(0.0, 0.0, -0.5));
        assert_vec3_close(
            primitive.vertices[21],
            Vec3::new(0.185_410_2, 0.570_633_9, -0.5),
        );
        assert_vec3_close(primitive.vertices[22], Vec3::new(0.6, 0.0, -0.5));
        assert_vec3_close(primitive.vertices[23], Vec3::new(0.6, 0.0, 0.5));
        assert_vec3_close(
            primitive.vertices[24],
            Vec3::new(0.185_410_2, 0.570_633_9, 0.5),
        );
        assert_vec3_close(primitive.vertices[25], Vec3::new(0.0, 0.0, 0.5));
        assert_vec3_close(primitive.normals[20], Vec3::new(0.0, 0.0, -1.0));
        assert_vec3_close(primitive.normals[23], Vec3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn pyramid_primitives_follow_molstar_full_and_perforated_order() {
        // Derived from artifacts/molstar/src/mol-geo/primitive/pyramid.ts.
        let pyramid = molstar_pyramid_primitive(8, true);
        assert_eq!(pyramid.vertices.len(), 48);
        assert_eq!(pyramid.normals.len(), 48);
        assert_eq!(pyramid.faces.len(), 16);
        assert_vec3_close(
            pyramid.vertices[0],
            Vec3::new(0.554_327_7, 0.229_610_07, -0.5),
        );
        assert_vec3_close(
            pyramid.vertices[1],
            Vec3::new(0.229_610_07, 0.554_327_7, -0.5),
        );
        assert_vec3_close(pyramid.vertices[2], Vec3::new(0.0, 0.0, 0.5));

        let perforated = molstar_perforated_octagonal_pyramid_primitive();
        assert_eq!(perforated.vertices.len(), 24);
        assert_eq!(perforated.normals.len(), 24);
        assert_eq!(perforated.faces.len(), 8);
        assert_vec3_close(
            perforated.vertices[0],
            Vec3::new(0.554_327_7, 0.229_610_07, -0.5),
        );
        assert_vec3_close(
            perforated.vertices[1],
            Vec3::new(0.229_610_07, 0.554_327_7, -0.5),
        );
        assert_vec3_close(perforated.vertices[2], Vec3::new(0.0, 0.0, -0.5));
    }

    #[test]
    fn primitive_transform_uses_molstar_direction_transform_for_normals() {
        let mut mesh = Mesh::default();
        add_molstar_box_primitive(
            &mut mesh,
            MolstarPrimitiveTransform::from_axes(
                Vec3::default(),
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::new(0.0, 3.0, 0.0),
                Vec3::new(0.0, 0.0, 4.0),
            ),
            false,
        );

        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.normals.len(), 24);
        assert_eq!(mesh.faces.len(), 12);
        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, 1.0 / 3.0, 0.0));
        assert_vec3_close(mesh.normals[16], Vec3::new(0.0, 0.0, -0.25));
        assert_vec3_close(mesh.normals[20], Vec3::new(0.0, 0.0, 0.25));
    }

    #[test]
    fn polyhedron_primitive_uses_indexed_area_weighted_normal_smoothing() {
        let vertices = [
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(-0.5, -0.5, 0.5),
            Vec3::new(-0.5, 0.5, -0.5),
            Vec3::new(0.5, -0.5, -0.5),
        ];
        let indices = [2, 1, 0, 0, 3, 2, 1, 3, 0, 2, 3, 1];
        let primitive = molstar_polyhedron_primitive(&vertices, &indices, 2.0, 1);

        assert_eq!(primitive.vertices.len(), 10);
        assert_eq!(primitive.normals.len(), 10);
        assert_eq!(primitive.faces.len(), 16);
        assert!(
            primitive.vertices.len() < primitive.faces.len() * 3,
            "Mol* polyhedron shares vertices before smoothing normals"
        );
        for (vertex, normal) in primitive.vertices.iter().zip(&primitive.normals) {
            assert!(
                (vertex.length() - 2.0).abs() <= EPS,
                "polyhedron radius must be applied before normal smoothing"
            );
            assert!(
                (normal.length() - 1.0).abs() <= EPS,
                "indexed normals must be normalized after area-weight accumulation"
            );
        }
    }

    #[test]
    fn initial_curve_segment_can_copy_molstar_polymer_tube_size_window() {
        // Mol* polymer-tube-mesh copies the initial size window together
        // with the shortened curve frame vectors.
        let controls = CurveSegmentControls {
            sec_struc_first: false,
            sec_struc_last: false,
            p0: DVec3::from_vec3(Vec3::new(-2.0, 0.0, 0.0)),
            p1: DVec3::from_vec3(Vec3::new(-1.0, 0.2, 0.0)),
            p2: DVec3::from_vec3(Vec3::new(0.0, 0.0, 0.0)),
            p3: DVec3::from_vec3(Vec3::new(1.0, 0.2, 0.0)),
            p4: DVec3::from_vec3(Vec3::new(2.0, 0.0, 0.0)),
            d12: DVec3::from_vec3(Vec3::new(0.0, 1.0, 0.0)),
            d23: DVec3::from_vec3(Vec3::new(0.0, 1.0, 0.0)),
        };

        let samples = curve_segment_samples(
            &controls,
            [1.0, 2.0, 4.0],
            [10.0, 20.0, 40.0],
            0.5,
            0.5,
            2.0,
            true,
            false,
            true,
            false,
            8,
        );

        assert_eq!(samples.centers.len(), 5);
        assert_f32s_close(&samples.widths, &[2.0, 2.25, 2.5, 2.75, 3.0]);
        assert_f32s_close(&samples.heights, &[20.0, 22.5, 25.0, 27.5, 30.0]);
    }

    #[test]
    fn initial_curve_segment_keeps_molstar_polymer_trace_size_window() {
        // Mol* polymer-trace-mesh adjusts the initial curve frame before
        // interpolateSizes, so width/height stay at the leading size window.
        let controls = CurveSegmentControls {
            sec_struc_first: false,
            sec_struc_last: false,
            p0: DVec3::from_vec3(Vec3::new(-2.0, 0.0, 0.0)),
            p1: DVec3::from_vec3(Vec3::new(-1.0, 0.2, 0.0)),
            p2: DVec3::from_vec3(Vec3::new(0.0, 0.0, 0.0)),
            p3: DVec3::from_vec3(Vec3::new(1.0, 0.2, 0.0)),
            p4: DVec3::from_vec3(Vec3::new(2.0, 0.0, 0.0)),
            d12: DVec3::from_vec3(Vec3::new(0.0, 1.0, 0.0)),
            d23: DVec3::from_vec3(Vec3::new(0.0, 1.0, 0.0)),
        };

        let samples = curve_segment_samples(
            &controls,
            [1.0, 2.0, 4.0],
            [10.0, 20.0, 40.0],
            0.5,
            0.5,
            2.0,
            true,
            false,
            false,
            false,
            8,
        );

        assert_eq!(samples.centers.len(), 5);
        let expected_start = (controls.p2
            + (controls.p2 - DVec3::from_vec3(samples.centers[1])).normalized() * 4.0)
            .to_vec3();
        assert_vec3_close(samples.centers[0], expected_start);
        assert_f32s_close(&samples.widths, &[1.5, 1.625, 1.75, 1.875, 2.0]);
        assert_f32s_close(&samples.heights, &[15.0, 16.25, 17.5, 18.75, 20.0]);
    }

    #[test]
    fn curve_segment_ribbon_can_swap_molstar_nucleic_width_height_values() {
        let controls = CurveSegmentControls {
            sec_struc_first: false,
            sec_struc_last: false,
            p0: DVec3::from_vec3(Vec3::new(-2.0, 0.0, 0.0)),
            p1: DVec3::from_vec3(Vec3::new(-1.0, 0.0, 0.0)),
            p2: DVec3::from_vec3(Vec3::new(0.0, 0.0, 0.0)),
            p3: DVec3::from_vec3(Vec3::new(1.0, 0.0, 0.0)),
            p4: DVec3::from_vec3(Vec3::new(2.0, 0.0, 0.0)),
            d12: DVec3::from_vec3(Vec3::new(0.0, 1.0, 0.0)),
            d23: DVec3::from_vec3(Vec3::new(0.0, 1.0, 0.0)),
        };
        let mut scratch = CurveSegmentScratch::default();
        let mut standard = Mesh::default();
        add_curve_segment_ribbon(
            &mut standard,
            &controls,
            [1.0; 3],
            [5.0; 3],
            0.5,
            0.5,
            1.0,
            0.0,
            false,
            false,
            false,
            false,
            4,
            &mut scratch,
        );
        let mut swapped = Mesh::default();
        add_curve_segment_ribbon(
            &mut swapped,
            &controls,
            [1.0; 3],
            [5.0; 3],
            0.5,
            0.5,
            1.0,
            0.0,
            false,
            false,
            false,
            true,
            4,
            &mut scratch,
        );

        assert_eq!(standard.vertices.len(), 20);
        assert_eq!(swapped.vertices.len(), 20);
        assert!(
            vertex_axis_span(&standard.vertices, Axis::Y)
                > vertex_axis_span(&swapped.vertices, Axis::Y) * 4.5,
            "Mol* nucleic radialSegments=2 swaps width/height arrays before addRibbon"
        );
    }

    #[test]
    fn curve_segment_scratch_reuses_state_and_sample_buffers() {
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
        let mut scratch = CurveSegmentScratch::default();

        curve_segment_samples_into(
            &mut scratch,
            &controls,
            [1.0, 2.0, 4.0],
            [10.0, 20.0, 40.0],
            0.5,
            0.5,
            2.0,
            false,
            false,
            true,
            false,
            16,
        );
        let state_ptrs = [
            scratch.state.curve_points.as_ptr() as usize,
            scratch.state.tangent_vectors.as_ptr() as usize,
            scratch.state.normal_vectors.as_ptr() as usize,
            scratch.state.binormal_vectors.as_ptr() as usize,
            scratch.state.width_values.as_ptr() as usize,
            scratch.state.height_values.as_ptr() as usize,
        ];
        let sample_ptrs = [
            scratch.samples.centers.as_ptr() as usize,
            scratch.samples.normals.as_ptr() as usize,
            scratch.samples.binormals.as_ptr() as usize,
            scratch.samples.widths.as_ptr() as usize,
            scratch.samples.heights.as_ptr() as usize,
        ];

        curve_segment_samples_into(
            &mut scratch,
            &controls,
            [1.0, 2.0, 4.0],
            [10.0, 20.0, 40.0],
            0.5,
            0.5,
            2.0,
            true,
            false,
            true,
            true,
            8,
        );

        assert_eq!(
            state_ptrs,
            [
                scratch.state.curve_points.as_ptr() as usize,
                scratch.state.tangent_vectors.as_ptr() as usize,
                scratch.state.normal_vectors.as_ptr() as usize,
                scratch.state.binormal_vectors.as_ptr() as usize,
                scratch.state.width_values.as_ptr() as usize,
                scratch.state.height_values.as_ptr() as usize,
            ]
        );
        assert_eq!(
            sample_ptrs,
            [
                scratch.samples.centers.as_ptr() as usize,
                scratch.samples.normals.as_ptr() as usize,
                scratch.samples.binormals.as_ptr() as usize,
                scratch.samples.widths.as_ptr() as usize,
                scratch.samples.heights.as_ptr() as usize,
            ]
        );
    }

    #[test]
    fn sheet_builder_uses_raw_molstar_frame_vectors() {
        let mut mesh = Mesh::default();
        let samples = CurveSamples {
            centers: vec![Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
            normals: vec![Vec3::new(0.0, 2.0, 0.0); 2],
            binormals: vec![Vec3::new(0.0, 0.0, 3.0); 2],
            widths: vec![1.0; 2],
            heights: vec![1.0; 2],
        };

        add_sheet_samples(&mut mesh, &samples, 0.0, true, false);

        assert_vec3_close(mesh.vertices[0], Vec3::new(0.0, 2.0, 3.0));
        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, 2.0, 0.0));
        assert_vec3_close(mesh.normals[2], Vec3::new(0.0, 0.0, -3.0));
        assert_vec3_close(mesh.normals[16], Vec3::new(-6.0, 0.0, 0.0));
    }

    #[test]
    fn tube_caps_use_raw_molstar_cross_normals() {
        let mut mesh = Mesh::default();
        add_profile_tube_for_test(
            &mut mesh,
            vec![Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
            vec![Vec3::new(0.0, 2.0, 0.0); 2],
            vec![Vec3::new(0.0, 0.0, 3.0); 2],
            vec![1.0; 2],
            vec![1.0; 2],
            4,
            TestTubeProfile::Elliptical,
            true,
            false,
            false,
        );

        assert_vec3_close(mesh.normals[8], Vec3::new(-6.0, 0.0, 0.0));
    }

    #[test]
    fn oriented_profile_radial_two_uses_molstar_ribbon_builder() {
        let mut mesh = Mesh::default();
        add_oriented_ribbon_with_profile(
            &mut mesh,
            &[Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
            &[Vec3::new(0.0, 0.0, 1.0); 2],
            2.0,
            0.5,
            PolymerProfile::Elliptical,
            true,
            true,
            false,
            1,
            2,
        );

        assert_eq!(mesh.vertices.len(), 8);
        assert_eq!(mesh.faces.len(), 4);
        assert_vec3_close(mesh.vertices[0], Vec3::new(0.0, 0.0, 0.5));
        assert_vec3_close(mesh.vertices[1], Vec3::new(0.0, 0.0, -0.5));
    }

    #[test]
    fn oriented_profile_radial_four_uses_molstar_sheet_builder() {
        let mut mesh = Mesh::default();
        add_oriented_ribbon_with_profile(
            &mut mesh,
            &[Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
            &[Vec3::new(0.0, 0.0, 1.0); 2],
            2.0,
            0.5,
            PolymerProfile::Elliptical,
            true,
            true,
            false,
            1,
            4,
        );

        assert_eq!(mesh.vertices.len(), 24);
        assert_eq!(mesh.faces.len(), 12);
        assert_vec3_close(mesh.vertices[0], Vec3::new(0.0, -2.0, 0.5));
        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn oriented_square_profile_keeps_oriented_frame_samples() {
        let mut mesh = Mesh::default();
        add_oriented_ribbon_with_profile(
            &mut mesh,
            &[Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
            &[Vec3::new(0.0, 0.0, 1.0); 2],
            2.0,
            0.5,
            PolymerProfile::Square,
            false,
            false,
            false,
            1,
            12,
        );

        assert_eq!(mesh.vertices.len(), 16);
        assert_eq!(mesh.faces.len(), 8);
        assert_vec3_close(mesh.vertices[0], Vec3::new(0.0, -2.0, 0.5));
        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn round_tube_cap_smoothing_uses_molstar_epsilon_after_sqrt() {
        let mut mesh = Mesh::default();
        add_profile_tube_for_test(
            &mut mesh,
            vec![Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
            vec![Vec3::new(0.0, 1.0, 0.0); 2],
            vec![Vec3::new(0.0, 0.0, 1.0); 2],
            vec![1.0; 2],
            vec![1.0; 2],
            4,
            TestTubeProfile::Elliptical,
            true,
            false,
            true,
        );

        assert!(
            mesh.vertices[0].length() <= MOLSTAR_NUMBER_EPSILON * 2.0,
            "round-cap endpoint must use Number.EPSILON, not sqrt(Number.EPSILON)"
        );
    }

    #[test]
    fn rounded_tube_corner_normals_are_normalized_like_molstar() {
        let mut mesh = Mesh::default();
        add_profile_tube_for_test(
            &mut mesh,
            vec![Vec3::default(), Vec3::new(1.0, 0.0, 0.0)],
            vec![Vec3::new(0.0, 2.0, 0.0); 2],
            vec![Vec3::new(0.0, 0.0, 3.0); 2],
            vec![1.0; 2],
            vec![3.0; 2],
            4,
            TestTubeProfile::Rounded,
            false,
            false,
            false,
        );

        assert_vec3_close(mesh.normals[0], Vec3::new(0.0, 0.0, 1.0));
        assert_vec3_close(mesh.normals[1], Vec3::new(0.0, 0.0, 1.0));
        assert_vec3_close(mesh.normals[2], Vec3::new(0.0, 0.0, -1.0));
        assert_vec3_close(mesh.normals[3], Vec3::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn normalize_molstar_uses_js_double_reciprocal_scale_before_float32_write() {
        let value = Vec3::new(
            f32::from_bits(0xc504_a74b),
            f32::from_bits(0x4295_f4a0),
            f32::from_bits(0x4501_2dcb),
        );

        let actual = normalize_molstar(value);
        assert_eq!(actual.x.to_bits(), 0xbf37_58bf);
        assert_eq!(actual.y.to_bits(), 0x3ccf_42bd);
        assert_eq!(actual.z.to_bits(), 0x3f32_8b54);

        let rust_f32_divide = value / value.squared_length().sqrt();
        assert_ne!(actual.x.to_bits(), rust_f32_divide.x.to_bits());
    }

    #[test]
    fn rounded_tube_surface_point_keeps_molstar_number_staging_until_buffer_write() {
        let mut mesh = Mesh::default();
        add_profile_tube_for_test(
            &mut mesh,
            vec![
                Vec3::new(570.862_3, 0.0, 0.0),
                Vec3::new(571.862_3, 0.0, 0.0),
            ],
            vec![Vec3::new(6.937_694, 1.0, 0.0); 2],
            vec![Vec3::new(-8.312_337, 0.0, 1.0); 2],
            vec![1.035_237_6; 2],
            vec![3.283_583_6; 2],
            8,
            TestTubeProfile::Rounded,
            false,
            false,
            false,
        );

        assert_eq!(mesh.vertices[0].x.to_bits(), 0x4413_7365);
        assert_ne!(
            mesh.vertices[0].x.to_bits(),
            0x4413_7364,
            "Mol* keeps surfacePoint as a JS-number Vec3 until ChunkedArray.add3 writes Float32"
        );
    }

    #[test]
    fn tube_cos_sin_values_follow_molstar_cache_sequence() {
        let (cos, sin) = molstar_tube_cos_sin(28, false);
        assert_eq!(cos[8].to_bits(), 0xbfcc_7b90_e302_4580);
        assert_eq!(sin[8].to_bits(), 0x3fef_329c_0558_e969);

        let (cos, sin) = molstar_tube_cos_sin(28, true);
        assert_eq!(cos[8].to_bits(), 0xbfd5_234a_ca69_a9f9);
        assert_eq!(sin[8].to_bits(), 0x3fee_344a_d05d_3f87);
    }

    #[test]
    fn tube_cos_sin_cache_matches_uncached_builder_for_all_supported_keys() {
        for radial_segments in 2..=MAX_TUBE_RADIAL_SEGMENTS {
            for shifted in [false, true] {
                let (cached_cos, cached_sin) = molstar_tube_cos_sin(radial_segments, shifted);
                let (cached_cos_again, cached_sin_again) =
                    molstar_tube_cos_sin(radial_segments, shifted);
                assert!(std::ptr::eq(cached_cos, cached_cos_again));
                assert!(std::ptr::eq(cached_sin, cached_sin_again));

                let uncached = build_molstar_tube_cos_sin(radial_segments, shifted);
                assert_eq!(cached_cos.len(), uncached.cos.len());
                assert_eq!(cached_sin.len(), uncached.sin.len());
                for (index, (actual, expected)) in cached_cos.iter().zip(&uncached.cos).enumerate()
                {
                    assert_eq!(
                        actual.to_bits(),
                        expected.to_bits(),
                        "cos radial_segments={radial_segments} shifted={shifted} index={index}"
                    );
                }
                for (index, (actual, expected)) in cached_sin.iter().zip(&uncached.sin).enumerate()
                {
                    assert_eq!(
                        actual.to_bits(),
                        expected.to_bits(),
                        "sin radial_segments={radial_segments} shifted={shifted} index={index}"
                    );
                }
            }
        }
    }

    #[test]
    fn curve_segment_tension_keeps_molstar_js_number_precision() {
        let controls = CurveSegmentControls {
            sec_struc_first: false,
            sec_struc_last: false,
            p0: DVec3::new(9.725000381469727, -5.230999946594238, 12.611000061035156),
            p1: DVec3::new(10.869000434875488, -2.7249999046325684, 15.244000434875488),
            p2: DVec3::new(7.271999835968018, -2.3489999771118164, 16.45199966430664),
            p3: DVec3::new(7.002999782562256, -6.13700008392334, 16.729000091552734),
            p4: DVec3::new(10.32699966430664, -6.208000183105469, 18.59000015258789),
            d12: DVec3::new(
                -0.1995002031326294,
                -0.28675001859664917,
                1.1344995498657227,
            ),
            d23: DVec3::new(
                -0.16650032997131348,
                -0.5555000305175781,
                1.0339999198913574,
            ),
        };
        let mut state = CurveSegmentState::new(14);
        interpolate_curve_segment(&mut state, &controls, 0.9, 0.5);

        assert_eq!(state.tangent_vectors[0].z.to_bits(), 0x3e1a_22ed);
        assert_eq!(state.normal_vectors[10].y.to_bits(), 0xbd3a_e613);
        assert_eq!(state.binormal_vectors[10].y.to_bits(), 0x3ef8_41d1);
    }

    #[derive(Clone, Copy)]
    enum Axis {
        Y,
    }

    fn vertex_axis_span(vertices: &[Vec3], axis: Axis) -> f32 {
        let values = vertices.iter().map(|v| match axis {
            Axis::Y => v.y,
        });
        let mut min = f32::INFINITY;
        let mut max = f32::NEG_INFINITY;
        for value in values {
            min = min.min(value);
            max = max.max(value);
        }
        max - min
    }

    fn assert_vec3_close(actual: Vec3, expected: Vec3) {
        assert!(
            (actual.x - expected.x).abs() <= EPS
                && (actual.y - expected.y).abs() <= EPS
                && (actual.z - expected.z).abs() <= EPS,
            "actual={actual:?} expected={expected:?}"
        );
    }

    fn assert_vec3s_close(actual: &[Vec3], expected: &[Vec3]) {
        assert_eq!(actual.len(), expected.len());
        for (i, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (actual.x - expected.x).abs() <= EPS
                    && (actual.y - expected.y).abs() <= EPS
                    && (actual.z - expected.z).abs() <= EPS,
                "vertex {i}: actual={actual:?} expected={expected:?}"
            );
        }
    }

    fn assert_f32s_close(actual: &[f32], expected: &[f32]) {
        assert_eq!(actual.len(), expected.len());
        for (i, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (actual - expected).abs() <= EPS,
                "value {i}: actual={actual:?} expected={expected:?}"
            );
        }
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

pub(super) fn add_ellipsoid(mesh: &mut Mesh, center: Vec3, axes: [Vec3; 3], detail: usize) {
    let primitive = molstar_sphere_primitive(detail);
    let normal_transform = molstar_direction_transform(axes);
    let base = mesh.vertices.len();

    for (&vertex, &normal) in primitive.vertices.iter().zip(&primitive.normals) {
        mesh.vertices
            .push(molstar_transform_position(center, axes, vertex));
        mesh.normals
            .push(molstar_transform_direction(normal_transform, normal));
    }
    for face in &primitive.faces {
        mesh.faces.push(Face {
            a: base + face.a,
            b: base + face.b,
            c: base + face.c,
        });
    }
}

fn molstar_transform_position(center: Vec3, axes: [Vec3; 3], point: Vec3) -> Vec3 {
    center + axes[0] * point.x + axes[1] * point.y + axes[2] * point.z
}

fn molstar_transform_direction(transform: [f32; 9], direction: Vec3) -> Vec3 {
    Vec3::new(
        direction.x * transform[0] + direction.y * transform[3] + direction.z * transform[6],
        direction.x * transform[1] + direction.y * transform[4] + direction.z * transform[7],
        direction.x * transform[2] + direction.y * transform[5] + direction.z * transform[8],
    )
}

fn molstar_direction_transform(axes: [Vec3; 3]) -> [f32; 9] {
    let mut out = [
        axes[0].x, axes[0].y, axes[0].z, axes[1].x, axes[1].y, axes[1].z, axes[2].x, axes[2].y,
        axes[2].z,
    ];
    let a = out;
    let a00 = a[0];
    let a01 = a[1];
    let a02 = a[2];
    let a10 = a[3];
    let a11 = a[4];
    let a12 = a[5];
    let a20 = a[6];
    let a21 = a[7];
    let a22 = a[8];

    let b01 = a22 * a11 - a12 * a21;
    let b11 = -a22 * a10 + a12 * a20;
    let b21 = a21 * a10 - a11 * a20;
    let det = a00 * b01 + a01 * b11 + a02 * b21;

    if det != 0.0 {
        let det = 1.0 / det;
        out = [
            b01 * det,
            (-a22 * a01 + a02 * a21) * det,
            (a12 * a01 - a02 * a11) * det,
            b11 * det,
            (a22 * a00 - a02 * a20) * det,
            (-a12 * a00 + a02 * a10) * det,
            b21 * det,
            (-a21 * a00 + a01 * a20) * det,
            (a11 * a00 - a01 * a10) * det,
        ];
    }

    [
        out[0], out[3], out[6], out[1], out[4], out[7], out[2], out[5], out[8],
    ]
}

#[cfg(test)]
pub(super) fn add_cylinder(mesh: &mut Mesh, start: Vec3, end: Vec3, radius: f32, segments: usize) {
    let mut cache = CylinderPrimitiveCache::default();
    add_cylinder_with_caps_cached(mesh, start, end, radius, segments, true, true, &mut cache);
}

pub(super) fn add_open_cylinder_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    cache: &mut CylinderPrimitiveCache,
) {
    add_cylinder_with_caps_cached(mesh, start, end, radius, segments, false, false, cache);
}

pub(super) fn add_molstar_buffered_open_cylinder_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    cache: &mut CylinderPrimitiveCache,
) {
    add_molstar_buffered_open_cylinder_with_radius64_cached(
        mesh,
        start,
        end,
        molstar_js_number_from_common_f32(radius),
        segments,
        cache,
    );
}

pub(super) fn add_molstar_buffered_open_cylinder_with_radius64_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f64,
    segments: usize,
    cache: &mut CylinderPrimitiveCache,
) {
    let start = DVec3::from_vec3(start);
    let end = DVec3::from_vec3(end);
    let dir = end - start;
    let length = start.distance(end);
    add_cylinder_dvec3_from_dir_with_caps_and_match_dir_cached(
        mesh,
        start,
        dir,
        length,
        radius,
        segments,
        CylinderBuildMode {
            start_cap: false,
            end_cap: false,
            match_dir: true,
        },
        cache,
    );
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_molstar_cylinder_caps_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    top_cap: bool,
    bottom_cap: bool,
    cache: &mut CylinderPrimitiveCache,
) {
    add_cylinder_with_caps_cached(
        mesh, start, end, radius, segments, bottom_cap, top_cap, cache,
    );
}

#[allow(clippy::too_many_arguments)]
fn add_cylinder_with_caps_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    start_cap: bool,
    end_cap: bool,
    cache: &mut CylinderPrimitiveCache,
) {
    add_cylinder_with_caps_and_match_dir_cached(
        mesh,
        start,
        end,
        radius,
        segments,
        CylinderBuildMode {
            start_cap,
            end_cap,
            match_dir: true,
        },
        cache,
    );
}

#[derive(Clone, Copy)]
struct CylinderBuildMode {
    start_cap: bool,
    end_cap: bool,
    match_dir: bool,
}

fn add_cylinder_with_caps_and_match_dir_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    mode: CylinderBuildMode,
    cache: &mut CylinderPrimitiveCache,
) {
    let axis = end - start;
    let length = axis.length();
    if length <= 0.001 {
        return;
    }
    let primitive = cache.get(
        segments,
        mode.end_cap,
        mode.start_cap,
        molstar_js_number_from_common_f32(radius),
    );
    let rotation = MolstarRotation::from_up_to_axis(axis, mode.match_dir);
    let center = start + axis * 0.5;
    let base = mesh.vertices.len();

    for (&vertex, &normal) in primitive.vertices.iter().zip(&primitive.normals) {
        let scaled_vertex = Vec3::new(vertex.x, vertex.y * length, vertex.z);
        let scaled_normal = Vec3::new(normal.x, normal.y / length, normal.z);
        mesh.vertices.push(center + rotation.apply(scaled_vertex));
        mesh.normals.push(rotation.apply(scaled_normal));
    }
    for face in &primitive.faces {
        mesh.faces.push(Face {
            a: base + face.a,
            b: base + face.b,
            c: base + face.c,
        });
    }
}

#[cfg(test)]
fn molstar_get_cylinder_primitive(
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
) -> Primitive {
    molstar_get_cylinder_primitive_with_radius(radial_segments, top_cap, bottom_cap, 1.0)
}

#[cfg(test)]
fn molstar_get_cylinder_primitive_with_radius(
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
    radius: f32,
) -> Primitive {
    molstar_get_cylinder_primitive_with_radius64(
        radial_segments,
        top_cap,
        bottom_cap,
        molstar_js_number_from_common_f32(radius),
    )
}

#[cfg(test)]
fn molstar_get_cylinder_primitive_with_radius64(
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
    radius: f64,
) -> Primitive {
    build_molstar_cylinder_primitive_with_radius64(radial_segments, top_cap, bottom_cap, radius)
}

fn build_molstar_cylinder_primitive_with_radius64(
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
    radius: f64,
) -> Primitive {
    let radial_segments = radial_segments.max(3);
    if radial_segments <= 4 {
        let prism = molstar_prism_primitive_with_radius(
            radial_segments,
            true,
            top_cap,
            bottom_cap,
            radius as f32,
        );
        molstar_transform_primitive_rot_x90(prism)
    } else {
        molstar_cylinder_primitive_with_props(CylinderPrimitiveProps {
            radial_segments,
            height_segments: 1,
            radius_top: radius,
            radius_bottom: radius,
            height: 1.0,
            top_cap,
            bottom_cap,
            theta_start: 0.0,
            theta_length: std::f64::consts::PI * 2.0,
        })
    }
}

#[cfg(test)]
fn molstar_cylinder_primitive(
    radial_segments: usize,
    top_cap: bool,
    bottom_cap: bool,
) -> Primitive {
    molstar_cylinder_primitive_with_props(CylinderPrimitiveProps {
        radial_segments,
        height_segments: 1,
        radius_top: 1.0,
        radius_bottom: 1.0,
        height: 1.0,
        top_cap,
        bottom_cap,
        theta_start: 0.0,
        theta_length: std::f64::consts::PI * 2.0,
    })
}

#[derive(Clone, Copy)]
struct CylinderPrimitiveProps {
    radial_segments: usize,
    height_segments: usize,
    radius_top: f64,
    radius_bottom: f64,
    height: f64,
    top_cap: bool,
    bottom_cap: bool,
    theta_start: f64,
    theta_length: f64,
}

fn molstar_cylinder_primitive_with_props(props: CylinderPrimitiveProps) -> Primitive {
    let CylinderPrimitiveProps {
        radial_segments,
        height_segments,
        radius_top,
        radius_bottom,
        height,
        top_cap,
        bottom_cap,
        theta_start,
        theta_length,
    } = props;
    let radial_segments = radial_segments.max(3);
    let height_segments = height_segments.max(1);
    let height = height.max(0.000_001);
    let half_height = height * 0.5;
    let slope = (radius_bottom - radius_top) / height;
    let mut vertices = Vec::new();
    let mut normals = Vec::new();
    let mut faces = Vec::new();
    let mut index = 0usize;
    let mut index_array = Vec::with_capacity(height_segments + 1);

    for y in 0..=height_segments {
        let mut index_row = Vec::with_capacity(radial_segments + 1);
        let v = y as f64 / height_segments as f64;
        let radius = v * (radius_bottom - radius_top) + radius_top;
        for x in 0..=radial_segments {
            let u = x as f64 / radial_segments as f64;
            let theta = u * theta_length + theta_start;
            let sin_theta = theta.sin();
            let cos_theta = theta.cos();
            vertices.push(Vec3::new(
                (radius * sin_theta) as f32,
                (-v * height + half_height) as f32,
                (radius * cos_theta) as f32,
            ));
            normals.push(normalize_cylinder_normal(sin_theta, slope, cos_theta));
            index_row.push(index);
            index += 1;
        }
        index_array.push(index_row);
    }

    for x in 0..radial_segments {
        for y in 0..height_segments {
            let a = index_array[y][x];
            let b = index_array[y + 1][x];
            let c = index_array[y + 1][x + 1];
            let d = index_array[y][x + 1];
            faces.push(Face { a, b, c: d });
            faces.push(Face { a: b, b: c, c: d });
        }
    }

    if top_cap && radius_top > 0.0 {
        generate_molstar_cylinder_cap(
            &mut vertices,
            &mut normals,
            &mut faces,
            &mut index,
            CylinderCapProps {
                radial_segments,
                radius: radius_top,
                half_height,
                theta_start,
                theta_length,
                top: true,
            },
        );
    }
    if bottom_cap && radius_bottom > 0.0 {
        generate_molstar_cylinder_cap(
            &mut vertices,
            &mut normals,
            &mut faces,
            &mut index,
            CylinderCapProps {
                radial_segments,
                radius: radius_bottom,
                half_height,
                theta_start,
                theta_length,
                top: false,
            },
        );
    }

    Primitive {
        vertices,
        normals,
        faces,
    }
}

#[derive(Clone, Copy)]
struct CylinderCapProps {
    radial_segments: usize,
    radius: f64,
    half_height: f64,
    theta_start: f64,
    theta_length: f64,
    top: bool,
}

fn generate_molstar_cylinder_cap(
    vertices: &mut Vec<Vec3>,
    normals: &mut Vec<Vec3>,
    faces: &mut Vec<Face>,
    index: &mut usize,
    props: CylinderCapProps,
) {
    let CylinderCapProps {
        radial_segments,
        radius,
        half_height,
        theta_start,
        theta_length,
        top,
    } = props;
    let sign = if top { 1.0 } else { -1.0 };
    let center_index_start = *index;

    for _ in 1..=radial_segments {
        vertices.push(Vec3::new(0.0, (half_height * sign) as f32, 0.0));
        normals.push(Vec3::new(0.0, sign as f32, 0.0));
        *index += 1;
    }

    let center_index_end = *index;
    for x in 0..=radial_segments {
        let u = x as f64 / radial_segments as f64;
        let theta = u * theta_length + theta_start;
        vertices.push(Vec3::new(
            (radius * theta.sin()) as f32,
            (half_height * sign) as f32,
            (radius * theta.cos()) as f32,
        ));
        normals.push(Vec3::new(0.0, sign as f32, 0.0));
        *index += 1;
    }

    for x in 0..radial_segments {
        let c = center_index_start + x;
        let i = center_index_end + x;
        if top {
            faces.push(Face { a: i, b: i + 1, c });
        } else {
            faces.push(Face { a: i + 1, b: i, c });
        }
    }
}

fn normalize_cylinder_normal(x: f64, y: f64, z: f64) -> Vec3 {
    let len_sq = x * x + y * y + z * z;
    if len_sq > 0.0 {
        let scale = 1.0 / len_sq.sqrt();
        Vec3::new((x * scale) as f32, (y * scale) as f32, (z * scale) as f32)
    } else {
        Vec3::default()
    }
}

fn molstar_js_number_from_common_f32(value: f32) -> f64 {
    if value.to_bits() == 0.2f32.to_bits() {
        0.2
    } else if value.to_bits() == 0.3f32.to_bits() {
        0.3
    } else if value.to_bits() == 0.5f32.to_bits() {
        0.5
    } else if value.to_bits() == 1.0f32.to_bits() {
        1.0
    } else {
        value as f64
    }
}

#[derive(Clone, Copy, Debug)]
struct MolstarRotation {
    m: [[f32; 3]; 3],
}

impl MolstarRotation {
    fn from_up_to_axis(axis: Vec3, match_dir: bool) -> Self {
        let target = axis.normalized();
        let source = if !match_dir || target.y > 0.0 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, -1.0, 0.0)
        };
        let dot = source.dot(target).clamp(-1.0, 1.0);
        let angle = dot.acos();
        if angle.abs() < 0.0001 {
            return Self::identity();
        }
        if (angle - PI).abs() < f32::EPSILON {
            let rotation_axis = if source.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 0.0, 1.0)
            };
            return Self::from_axis_angle(rotation_axis, PI);
        }
        Self::from_axis_angle(source.cross(target), angle)
    }

    fn identity() -> Self {
        Self {
            m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    fn from_axis_angle(axis: Vec3, angle: f32) -> Self {
        let axis = axis.normalized();
        let s = angle.sin();
        let c = angle.cos();
        let t = 1.0 - c;
        let x = axis.x;
        let y = axis.y;
        let z = axis.z;
        Self {
            m: [
                [x * x * t + c, x * y * t - z * s, x * z * t + y * s],
                [y * x * t + z * s, y * y * t + c, y * z * t - x * s],
                [z * x * t - y * s, z * y * t + x * s, z * z * t + c],
            ],
        }
    }

    fn apply(self, v: Vec3) -> Vec3 {
        Vec3::new(
            self.m[0][0] * v.x + self.m[0][1] * v.y + self.m[0][2] * v.z,
            self.m[1][0] * v.x + self.m[1][1] * v.y + self.m[1][2] * v.z,
            self.m[2][0] * v.x + self.m[2][1] * v.y + self.m[2][2] * v.z,
        )
    }
}

pub(crate) fn add_tube_path(mesh: &mut Mesh, points: &[Vec3], radius: f32, segments: usize) {
    if points.len() < 2 {
        return;
    }
    let samples = sample_sheet_path(points, 4, radius, radius);
    if samples.centers.len() < 2 {
        return;
    }
    add_profile_tube(
        mesh,
        &samples,
        segments.max(3),
        TubeProfile::Elliptical,
        true,
        true,
        false,
    );
}

pub(super) fn add_dashed_tube_path_cached(
    mesh: &mut Mesh,
    points: &[Vec3],
    radius: f32,
    segments: usize,
    cache: &mut CylinderPrimitiveCache,
) {
    if points.len() < 2 {
        return;
    }
    let samples = sample_path(points, 8);
    add_dashed_tube_samples_cached(mesh, &samples, radius, segments, cache);
}

pub(super) fn add_dashed_tube_samples_cached(
    mesh: &mut Mesh,
    samples: &[Vec3],
    radius: f32,
    segments: usize,
    cache: &mut CylinderPrimitiveCache,
) {
    if samples.len() < 2 {
        return;
    }
    let dash_len = (radius * 3.8).max(0.55);
    let gap_len = (radius * 2.2).max(0.32);
    let period = dash_len + gap_len;
    let mut distance = 0.0;

    for pair in samples.windows(2) {
        let start = pair[0];
        let end = pair[1];
        let delta = end - start;
        let length = delta.length();
        if length <= 0.000_001 {
            continue;
        }
        let direction = delta / length;
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
                add_cylinder_with_caps_cached(
                    mesh,
                    start + direction * local,
                    start + direction * (local + step),
                    radius,
                    segments,
                    true,
                    true,
                    cache,
                );
            }
            local += step.max(0.000_001);
        }
        distance += length;
    }
}

#[cfg(test)]
pub(super) fn add_fixed_count_dashed_cylinder(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    length_scale: f32,
    segment_count: usize,
) {
    let mut cache = CylinderPrimitiveCache::default();
    add_fixed_count_dashed_cylinder_cached(
        mesh,
        start,
        end,
        radius,
        segments,
        length_scale,
        segment_count,
        &mut cache,
    );
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_fixed_count_dashed_cylinder_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    length_scale: f32,
    segment_count: usize,
    cache: &mut CylinderPrimitiveCache,
) {
    add_fixed_count_dashed_cylinder_with_stub_cap_cached(
        mesh,
        start,
        end,
        radius,
        segments,
        length_scale,
        segment_count,
        false,
        cache,
    );
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn add_fixed_count_dashed_cylinder_with_stub_cap(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    length_scale: f32,
    segment_count: usize,
    stub_cap: bool,
) {
    let mut cache = CylinderPrimitiveCache::default();
    add_fixed_count_dashed_cylinder_with_stub_cap_cached(
        mesh,
        start,
        end,
        radius,
        segments,
        length_scale,
        segment_count,
        stub_cap,
        &mut cache,
    );
}

#[allow(clippy::too_many_arguments)]
fn add_fixed_count_dashed_cylinder_with_stub_cap_cached(
    mesh: &mut Mesh,
    start: Vec3,
    end: Vec3,
    radius: f32,
    segments: usize,
    length_scale: f32,
    segment_count: usize,
    stub_cap: bool,
    cache: &mut CylinderPrimitiveCache,
) {
    let start = DVec3::from_vec3(start);
    let end = DVec3::from_vec3(end);
    let distance = start.distance(end) * length_scale as f64;
    if distance <= 0.000_001 || segment_count == 0 {
        return;
    }

    let dash_count = segment_count.div_ceil(2);
    let is_odd = !segment_count.is_multiple_of(2);
    let step = distance / (segment_count as f64 + 0.5);
    let direction = (end - start).set_magnitude(step);
    let mut cursor = start;

    for dash_index in 0..dash_count {
        cursor = cursor + direction;
        let is_last_odd_dash = is_odd && dash_index == dash_count - 1;
        let dash_length = if is_last_odd_dash { step * 0.5 } else { step };
        let end_cap = !is_last_odd_dash || stub_cap;
        add_cylinder_dvec3_from_dir_with_caps_and_match_dir_cached(
            mesh,
            cursor,
            direction,
            dash_length,
            molstar_js_number_from_common_f32(radius),
            segments,
            CylinderBuildMode {
                start_cap: true,
                end_cap,
                match_dir: false,
            },
            cache,
        );
        cursor = cursor + direction;
    }
}

fn add_cylinder_dvec3_from_dir_with_caps_and_match_dir_cached(
    mesh: &mut Mesh,
    start: DVec3,
    dir: DVec3,
    length: f64,
    radius: f64,
    segments: usize,
    mode: CylinderBuildMode,
    cache: &mut CylinderPrimitiveCache,
) {
    if length <= 0.001 {
        return;
    }
    let primitive = cache.get(segments, mode.end_cap, mode.start_cap, radius);
    let mat_dir = dir.set_magnitude(length * 0.5);
    let rotation = MolstarRotationD::from_up_to_mat_dir(mat_dir, mode.match_dir);
    let center = start + mat_dir;
    let base = mesh.vertices.len();

    for (&vertex, &normal) in primitive.vertices.iter().zip(&primitive.normals) {
        let transformed_vertex =
            molstar_transform_cylinder_position_d(center, rotation, length, vertex);
        let scaled_normal = DVec3::new(normal.x as f64, normal.y as f64 / length, normal.z as f64);
        mesh.vertices.push(transformed_vertex.to_vec3());
        mesh.normals.push(rotation.apply(scaled_normal).to_vec3());
    }
    for face in &primitive.faces {
        mesh.faces.push(Face {
            a: base + face.a,
            b: base + face.b,
            c: base + face.c,
        });
    }
}

pub(super) fn sample_path(points: &[Vec3], subdivisions: usize) -> Vec<Vec3> {
    let points = filtered_path_points(points);
    let points = points.as_slice();
    if points.len() < 2 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for i in 0..points.len() - 1 {
        let p0 = if i == 0 { points[i] } else { points[i - 1] };
        let p1 = points[i];
        let p2 = points[i + 1];
        let p3 = if i + 2 < points.len() {
            points[i + 2]
        } else {
            points[i + 1]
        };
        for step in 0..subdivisions {
            let t = step as f32 / subdivisions as f32;
            out.push(catmull_rom(p0, p1, p2, p3, t));
        }
    }
    out.push(*points.last().unwrap());
    out
}

pub(super) fn sample_path_point_count(points: &[Vec3], subdivisions: usize) -> usize {
    if subdivisions == 0 {
        return 0;
    }
    let mut filtered_count = 0usize;
    let mut previous = None;
    for &point in points {
        if previous.is_none_or(|previous: Vec3| previous.distance(point) > 0.000_1) {
            filtered_count += 1;
            previous = Some(point);
        }
    }
    if filtered_count < 2 {
        0
    } else {
        (filtered_count - 1)
            .saturating_mul(subdivisions)
            .saturating_add(1)
    }
}

fn filtered_path_points(points: &[Vec3]) -> Vec<Vec3> {
    let mut out = Vec::with_capacity(points.len());
    for point in points {
        if out
            .last()
            .is_none_or(|previous: &Vec3| previous.distance(*point) > 0.000_1)
        {
            out.push(*point);
        }
    }
    out
}

fn catmull_rom(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
    let t2 = t * t;
    let t3 = t2 * t;
    (p1 * 2.0
        + (p2 - p0) * t
        + (p0 * 2.0 - p1 * 5.0 + p2 * 4.0 - p3) * t2
        + (p3 - p0 + (p1 - p2) * 3.0) * t3)
        * 0.5
}

fn molstar_spline(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32, tension: f32) -> Vec3 {
    molstar_spline_f64(
        DVec3::from_vec3(p0),
        DVec3::from_vec3(p1),
        DVec3::from_vec3(p2),
        DVec3::from_vec3(p3),
        t as f64,
        tension as f64,
    )
    .to_vec3()
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct DVec3 {
    pub(crate) x: f64,
    pub(crate) y: f64,
    pub(crate) z: f64,
}

impl DVec3 {
    pub(crate) const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub(crate) fn from_vec3(value: Vec3) -> Self {
        Self {
            x: value.x as f64,
            y: value.y as f64,
            z: value.z as f64,
        }
    }

    pub(crate) fn to_vec3(self) -> Vec3 {
        Vec3::new(self.x as f32, self.y as f32, self.z as f32)
    }

    pub(crate) fn squared_length(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub(crate) fn length(self) -> f64 {
        self.squared_length().sqrt()
    }

    pub(crate) fn distance(self, other: Self) -> f64 {
        (self - other).length()
    }

    pub(crate) fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub(crate) fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    fn is_zero(self) -> bool {
        self.x == 0.0 && self.y == 0.0 && self.z == 0.0
    }

    pub(crate) fn normalized(self) -> Self {
        let len_sq = self.squared_length();
        if len_sq > 0.0 {
            self * (1.0 / len_sq.sqrt())
        } else {
            self
        }
    }

    fn set_magnitude(self, length: f64) -> Self {
        self.normalized() * length
    }
}

impl std::ops::Add for DVec3 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::Sub for DVec3 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Self {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl std::ops::Mul<f64> for DVec3 {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Self {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl std::ops::Div<f64> for DVec3 {
    type Output = Self;

    fn div(self, rhs: f64) -> Self {
        Self {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

fn molstar_spline_f64(p0: DVec3, p1: DVec3, p2: DVec3, p3: DVec3, t: f64, tension: f64) -> DVec3 {
    DVec3 {
        x: molstar_spline_scalar_f64(p0.x, p1.x, p2.x, p3.x, t, tension),
        y: molstar_spline_scalar_f64(p0.y, p1.y, p2.y, p3.y, t, tension),
        z: molstar_spline_scalar_f64(p0.z, p1.z, p2.z, p3.z, t, tension),
    }
}

fn molstar_spline_scalar_f64(p0: f64, p1: f64, p2: f64, p3: f64, t: f64, tension: f64) -> f64 {
    let v0 = (p2 - p0) * tension;
    let v1 = (p3 - p1) * tension;
    let t2 = t * t;
    let t3 = t * t2;
    (2.0 * p1 - 2.0 * p2 + v0 + v1) * t3 + (-3.0 * p1 + 3.0 * p2 - 2.0 * v0 - v1) * t2 + v0 * t + p1
}

pub(super) fn helix_trace(points: &[Vec3], directions: &[Option<Vec3>]) -> (Vec<Vec3>, Vec<Vec3>) {
    let mut centers = Vec::with_capacity(points.len());
    let mut normals = Vec::with_capacity(points.len());
    let window = if points.len() >= 7 { 3 } else { 2 };
    let mut previous_normal: Option<Vec3> = None;

    for i in 0..points.len() {
        let start = i.saturating_sub(window);
        let end = (i + window).min(points.len() - 1);
        let mut center = Vec3::default();
        let count = end - start + 1;
        for point in &points[start..=end] {
            center = center + *point;
        }
        center = center / count as f32;

        let mut normal = directions
            .get(i)
            .and_then(|direction| *direction)
            .unwrap_or_else(|| (points[i] - center).normalized());
        if normal.length() <= 0.000_001 {
            normal = previous_normal.unwrap_or(Vec3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            });
        }
        if let Some(previous) = previous_normal {
            if normal.dot(previous) < 0.0 {
                normal = normal * -1.0;
            }
        }

        centers.push(points[i]);
        normals.push(normal);
        previous_normal = Some(normal);
    }

    (centers, normals)
}

#[allow(dead_code)]
fn helix_orientation_centers(points: &[Vec3]) -> Option<Vec<Vec3>> {
    if points.len() < 4 {
        return None;
    }
    let mut centers = vec![Vec3::default(); points.len()];
    let mut have_center = vec![false; points.len()];

    for i in 3..points.len() {
        let a1 = points[i - 3];
        let a2 = points[i - 2];
        let a3 = points[i - 1];
        let a4 = points[i];

        let r12 = a2 - a1;
        let r23 = a3 - a2;
        let r34 = a4 - a3;
        let diff13 = r12 - r23;
        let diff24 = r23 - r34;
        let diff13_len = diff13.length();
        let diff24_len = diff24.length();
        if diff13_len <= 0.000_001 || diff24_len <= 0.000_001 {
            continue;
        }

        let axis = diff13.cross(diff24).normalized();
        if axis.length() <= 0.000_001 {
            continue;
        }
        let cos_angle = diff13
            .normalized()
            .dot(diff24.normalized())
            .clamp(-1.0, 1.0);
        let radius = (diff24_len * diff13_len).sqrt() / (2.0 * (1.0 - cos_angle)).max(2.0);

        let c1 = a2 - diff13 * (radius / diff13_len);
        let c2 = a3 - diff24 * (radius / diff24_len);
        centers[i - 2] = c1;
        centers[i - 1] = c2;
        have_center[i - 2] = true;
        have_center[i - 1] = true;
    }

    if points.len() >= 3 && have_center[1] && have_center[2] {
        let axis = (centers[1] - centers[2]).normalized();
        if axis.length() > 0.000_001 {
            centers[0] = project_point_on_vector(points[0], axis, centers[1]);
            have_center[0] = true;
        }
    }

    let last = points.len() - 1;
    if last >= 2 && have_center[last - 1] && have_center[last - 2] {
        let axis = (centers[last - 1] - centers[last - 2]).normalized();
        if axis.length() > 0.000_001 {
            centers[last] = project_point_on_vector(points[last], axis, centers[last - 1]);
            have_center[last] = true;
        }
    }

    if have_center.iter().all(|have| *have)
        && centers
            .windows(2)
            .all(|pair| pair[0].distance(pair[1]).is_finite() && pair[0].distance(pair[1]) > 0.01)
    {
        Some(centers)
    } else {
        None
    }
}

#[allow(dead_code)]
fn helix_orientation_normals(
    centers: &[Vec3],
    trace_points: &[Vec3],
    directions: &[Option<Vec3>],
) -> Vec<Vec3> {
    let mut normals = Vec::with_capacity(centers.len());
    let mut previous_normal = None;
    for i in 0..centers.len() {
        let tangent = path_tangent(centers, i);
        let carbonyl_direction = directions
            .get(i)
            .and_then(|direction| *direction)
            .unwrap_or(Vec3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            });
        let mut normal = orthogonalize(
            tangent,
            match previous_normal {
                Some(previous) => previous,
                None => carbonyl_direction,
            },
        );
        if normal.length() <= 0.000_001 {
            normal = orthogonalize(
                tangent,
                Vec3 {
                    x: 1.0,
                    y: 0.0,
                    z: 0.0,
                },
            );
        }
        if normal.length() <= 0.000_001 {
            normal = (trace_points
                .get(i)
                .map(|point| *point - centers[i])
                .unwrap_or_default()
                - tangent
                    * trace_points
                        .get(i)
                        .map(|point| (*point - centers[i]).dot(tangent))
                        .unwrap_or_default())
            .normalized();
        }
        if let Some(previous) = previous_normal {
            normal = match_direction(normal, previous);
        }
        normals.push(normal);
        previous_normal = Some(normal);
    }
    normals
}

#[allow(dead_code)]
fn project_point_on_vector(point: Vec3, direction: Vec3, origin: Vec3) -> Vec3 {
    let direction = direction.normalized();
    if direction.length() <= 0.000_001 {
        return point;
    }
    origin + direction * (point - origin).dot(direction)
}

#[allow(dead_code)]
pub(crate) fn add_oriented_ribbon(
    mesh: &mut Mesh,
    centers: &[Vec3],
    normals: &[Vec3],
    width: f32,
    thickness: f32,
) {
    add_oriented_ribbon_with_profile(
        mesh,
        centers,
        normals,
        width,
        thickness,
        PolymerProfile::Elliptical,
        true,
        true,
        false,
        12,
        32,
    );
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_oriented_ribbon_with_profile(
    mesh: &mut Mesh,
    centers: &[Vec3],
    normals: &[Vec3],
    width: f32,
    thickness: f32,
    profile: PolymerProfile,
    start_cap: bool,
    end_cap: bool,
    round_cap: bool,
    linear_segments: usize,
    radial_segments: usize,
) {
    if centers.len() < 2 || centers.len() != normals.len() {
        return;
    }
    let samples =
        sample_oriented_profile_path(centers, normals, linear_segments.max(1), width, thickness);
    if samples.centers.len() < 2 {
        return;
    }

    if radial_segments == 2 {
        add_ribbon_samples(mesh, &samples, 0.0);
        return;
    }
    if radial_segments == 4 || profile == PolymerProfile::Square {
        add_sheet_samples(mesh, &samples, 0.0, start_cap, end_cap);
        return;
    }

    match profile {
        PolymerProfile::Square => unreachable!("square profile routed to sheet samples above"),
        PolymerProfile::Rounded => add_profile_tube(
            mesh,
            &samples,
            radial_segments.max(3),
            TubeProfile::Rounded,
            start_cap,
            end_cap,
            round_cap,
        ),
        PolymerProfile::Elliptical => add_profile_tube(
            mesh,
            &samples,
            radial_segments.max(3),
            TubeProfile::Elliptical,
            start_cap,
            end_cap,
            round_cap,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_curve_segment_tube(
    mesh: &mut Mesh,
    controls: &CurveSegmentControls,
    widths: [f32; 3],
    heights: [f32; 3],
    tension: f64,
    shift: f64,
    overhang_width: f32,
    start_cap: bool,
    end_cap: bool,
    round_cap: bool,
    initial: bool,
    final_residue: bool,
    swap_normal_binormal: bool,
    linear_segments: usize,
    radial_segments: usize,
    profile: PolymerProfile,
    scratch: &mut CurveSegmentScratch,
) {
    let samples = curve_segment_samples_into(
        scratch,
        controls,
        widths,
        heights,
        tension,
        shift,
        overhang_width,
        initial,
        final_residue,
        false,
        swap_normal_binormal,
        linear_segments,
    );
    let profile = match profile {
        PolymerProfile::Rounded => TubeProfile::Rounded,
        _ => TubeProfile::Elliptical,
    };
    add_profile_tube(
        mesh,
        &samples,
        radial_segments.max(3),
        profile,
        start_cap,
        end_cap,
        round_cap,
    );
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_curve_segment_ribbon(
    mesh: &mut Mesh,
    controls: &CurveSegmentControls,
    widths: [f32; 3],
    heights: [f32; 3],
    tension: f64,
    shift: f64,
    overhang_width: f32,
    arrow_height: f32,
    initial: bool,
    final_residue: bool,
    swap_normal_binormal: bool,
    swap_width_height: bool,
    linear_segments: usize,
    scratch: &mut CurveSegmentScratch,
) {
    let samples = curve_segment_samples_into(
        scratch,
        controls,
        widths,
        heights,
        tension,
        shift,
        overhang_width,
        initial,
        final_residue,
        false,
        swap_normal_binormal,
        linear_segments,
    );
    if swap_width_height {
        std::mem::swap(&mut samples.widths, &mut samples.heights);
    }
    add_ribbon_samples(mesh, &samples, arrow_height);
}

#[allow(clippy::too_many_arguments)]
pub(super) fn add_curve_segment_sheet(
    mesh: &mut Mesh,
    controls: &CurveSegmentControls,
    widths: [f32; 3],
    heights: [f32; 3],
    tension: f64,
    shift: f64,
    overhang_width: f32,
    arrow_height: f32,
    start_cap: bool,
    end_cap: bool,
    initial: bool,
    final_residue: bool,
    swap_normal_binormal: bool,
    linear_segments: usize,
    scratch: &mut CurveSegmentScratch,
) {
    let samples = curve_segment_samples_into(
        scratch,
        controls,
        widths,
        heights,
        tension,
        shift,
        overhang_width,
        initial,
        final_residue,
        false,
        swap_normal_binormal,
        linear_segments,
    );
    add_sheet_samples(mesh, &samples, arrow_height, start_cap, end_cap);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TubeProfile {
    Elliptical,
    Rounded,
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum TestTubeProfile {
    Elliptical,
    Rounded,
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
pub(crate) fn add_profile_tube_for_test(
    mesh: &mut Mesh,
    centers: Vec<Vec3>,
    normals: Vec<Vec3>,
    binormals: Vec<Vec3>,
    widths: Vec<f32>,
    heights: Vec<f32>,
    radial_segments: usize,
    profile: TestTubeProfile,
    start_cap: bool,
    end_cap: bool,
    round_cap: bool,
) {
    let samples = CurveSamples {
        centers,
        normals,
        binormals,
        widths,
        heights,
    };
    let profile = match profile {
        TestTubeProfile::Elliptical => TubeProfile::Elliptical,
        TestTubeProfile::Rounded => TubeProfile::Rounded,
    };
    add_profile_tube(
        mesh,
        &samples,
        radial_segments,
        profile,
        start_cap,
        end_cap,
        round_cap,
    );
}

#[derive(Clone, Debug, Default)]
pub(crate) struct CurveSamples {
    pub(crate) centers: Vec<Vec3>,
    pub(crate) normals: Vec<Vec3>,
    pub(crate) binormals: Vec<Vec3>,
    pub(crate) widths: Vec<f32>,
    pub(crate) heights: Vec<f32>,
}

#[derive(Debug)]
pub(super) struct CurveSegmentScratch {
    state: CurveSegmentState,
    samples: CurveSamples,
}

impl Default for CurveSegmentScratch {
    fn default() -> Self {
        Self {
            state: CurveSegmentState::new(1),
            samples: CurveSamples::default(),
        }
    }
}

fn molstar_tube_cos_sin(radial_segments: usize, shifted: bool) -> (&'static [f64], &'static [f64]) {
    assert!(
        radial_segments <= MAX_TUBE_RADIAL_SEGMENTS,
        "tube radial segments {radial_segments} exceeds the supported maximum {MAX_TUBE_RADIAL_SEGMENTS}"
    );
    let index = radial_segments * 2 + usize::from(shifted);
    let table = TUBE_TRIG_TABLES[index]
        .get_or_init(|| build_molstar_tube_cos_sin(radial_segments, shifted));
    (&table.cos, &table.sin)
}

fn build_molstar_tube_cos_sin(radial_segments: usize, shifted: bool) -> TubeTrigTable {
    let offset = if shifted { 1.0 } else { 0.0 };
    let mut cos = Vec::with_capacity(radial_segments);
    let mut sin = Vec::with_capacity(radial_segments);
    for j in 0..radial_segments {
        let phi = (j as f64 * 2.0 + offset) / radial_segments as f64 * std::f64::consts::PI;
        cos.push(phi.cos());
        sin.push(phi.sin());
    }
    TubeTrigTable { cos, sin }
}

fn add_profile_tube(
    mesh: &mut Mesh,
    samples: &CurveSamples,
    radial_segments: usize,
    profile: TubeProfile,
    start_cap: bool,
    end_cap: bool,
    round_cap: bool,
) {
    let point_count = samples.centers.len();
    if point_count < 2
        || samples.normals.len() != point_count
        || samples.binormals.len() != point_count
        || samples.widths.len() != point_count
        || samples.heights.len() != point_count
        || radial_segments < 3
    {
        return;
    }

    let base = mesh.vertices.len();
    let q1 = (radial_segments as f32 / 4.0).round() as usize;
    let q3 = q1 * 3;
    let round_cap = round_cap && point_count > 1 && (start_cap || end_cap);
    let double_round_cap = round_cap && start_cap && end_cap;
    let half_linear_segments = (point_count - 1) as f32 / 2.0;
    let (cos_values, sin_values) =
        molstar_tube_cos_sin(radial_segments, profile == TubeProfile::Rounded);

    for i in 0..point_count {
        let u = samples.normals[i];
        let v = samples.binormals[i];
        let center = samples.centers[i];
        let mut width = samples.widths[i];
        let mut height = samples.heights[i];
        let mut cap_smoothing_factor = 1.0;
        let mut cap_normal_smoothing_vector = Vec3::default();
        if round_cap {
            let i_f = i as f32;
            let linear_segments = (point_count - 1) as f32;
            let smooth_start = if double_round_cap {
                i_f <= half_linear_segments
            } else {
                start_cap
            };
            let denom = if double_round_cap {
                half_linear_segments
            } else {
                linear_segments
            };
            let distance = if double_round_cap {
                if smooth_start {
                    half_linear_segments - i_f
                } else {
                    i_f - half_linear_segments
                }
            } else if smooth_start {
                linear_segments - i_f
            } else {
                i_f
            };
            let smoothing = (1.0 - (distance / denom).powi(2)).sqrt();
            cap_smoothing_factor = if smoothing < MOLSTAR_NUMBER_EPSILON {
                MOLSTAR_NUMBER_EPSILON
            } else {
                smoothing
            };
            width *= cap_smoothing_factor;
            height *= cap_smoothing_factor;
            cap_normal_smoothing_vector = if smooth_start {
                normalize_molstar(v.cross(u))
            } else {
                normalize_molstar(u.cross(v))
            };
        }
        let rounded = profile == TubeProfile::Rounded && height > width;

        for j in 0..radial_segments {
            let cos = cos_values[j];
            let sin = sin_values[j];
            let (vertex, mut normal, normal_is_normalized) = if rounded {
                let h = if v.dot(Vec3::new(1.0, 0.0, 0.0)) < 0.0 {
                    if j < q1 || j >= q3 {
                        height - width
                    } else {
                        -height + width
                    }
                } else if j >= q1 && j < q3 {
                    -height + width
                } else {
                    height - width
                };
                let vertex = molstar_vec3_add3_scaled_then_scale_and_add(
                    center,
                    u,
                    v,
                    width as f64 * cos,
                    width as f64 * sin,
                    h as f64,
                );
                let normal = if j == q1 || j + 1 == q1 {
                    v
                } else if j == q3 || j + 1 == q3 {
                    v * -1.0
                } else {
                    molstar_vec3_normalize_add2_scaled(u, v, cos, sin)
                };
                let normal_is_normalized = !(j == q1 || j + 1 == q1 || j == q3 || j + 1 == q3);
                (vertex, normal, normal_is_normalized)
            } else {
                (
                    molstar_vec3_add3_scaled(center, u, v, height as f64 * cos, width as f64 * sin),
                    molstar_vec3_normalize_add2_scaled(
                        u,
                        v,
                        width as f64 * cos,
                        height as f64 * sin,
                    ),
                    true,
                )
            };
            if !normal_is_normalized {
                normal = normalize_molstar_f64(normal);
            }
            if round_cap {
                normal = slerp_unit(cap_normal_smoothing_vector, normal, cap_smoothing_factor);
            }
            mesh.vertices.push(vertex);
            mesh.normals.push(normal);
        }
    }

    let half = (radial_segments as f32 / 2.0).round() as usize;
    for i in 0..point_count - 1 {
        for j in 0..half {
            let next = (j + 1) % radial_segments;
            mesh.faces.push(Face {
                a: base + i * radial_segments + next,
                b: base + (i + 1) * radial_segments + next,
                c: base + i * radial_segments + j,
            });
            mesh.faces.push(Face {
                a: base + (i + 1) * radial_segments + next,
                b: base + (i + 1) * radial_segments + j,
                c: base + i * radial_segments + j,
            });
        }
        for j in half..radial_segments {
            let next = (j + 1) % radial_segments;
            mesh.faces.push(Face {
                a: base + i * radial_segments + next,
                b: base + (i + 1) * radial_segments + j,
                c: base + i * radial_segments + j,
            });
            mesh.faces.push(Face {
                a: base + (i + 1) * radial_segments + next,
                b: base + (i + 1) * radial_segments + j,
                c: base + i * radial_segments + next,
            });
        }
    }

    if start_cap {
        add_profile_tube_cap(
            mesh,
            &samples.centers,
            &samples.normals,
            &samples.binormals,
            &samples.widths,
            &samples.heights,
            radial_segments,
            profile,
            0,
            true,
            round_cap,
        );
    }
    if end_cap {
        add_profile_tube_cap(
            mesh,
            &samples.centers,
            &samples.normals,
            &samples.binormals,
            &samples.widths,
            &samples.heights,
            radial_segments,
            profile,
            point_count - 1,
            false,
            round_cap,
        );
    }
}

#[allow(clippy::too_many_arguments)]
fn curve_segment_samples_into<'a>(
    scratch: &'a mut CurveSegmentScratch,
    controls: &CurveSegmentControls,
    widths: [f32; 3],
    heights: [f32; 3],
    tension: f64,
    shift: f64,
    overhang_width: f32,
    initial: bool,
    final_residue: bool,
    copy_initial_size_window: bool,
    swap_normal_binormal: bool,
    linear_segments: usize,
) -> &'a mut CurveSamples {
    const OVERHANG_FACTOR: f32 = 2.0;

    let linear_segments = linear_segments.max(1);
    let CurveSegmentScratch { state, samples } = scratch;
    state.prepare(linear_segments);
    interpolate_curve_segment(state, controls, tension, shift);
    let mut sizes_interpolated = false;
    if copy_initial_size_window {
        interpolate_sizes(
            state, widths[0], widths[1], widths[2], heights[0], heights[1], heights[2], shift,
        );
        sizes_interpolated = true;
    }

    let mut segment_count = linear_segments;
    if initial {
        segment_count = ((linear_segments as f64 * shift).round() as usize).max(1);
        let offset = linear_segments - segment_count;
        let source = offset..offset + segment_count + 1;
        state.curve_points.copy_within(source.clone(), 0);
        state.normal_vectors.copy_within(source.clone(), 0);
        state.binormal_vectors.copy_within(source.clone(), 0);
        if copy_initial_size_window {
            state.width_values.copy_within(source.clone(), 0);
            state.height_values.copy_within(source, 0);
        }

        let next = state.curve_points[1.min(segment_count)];
        let direction = (controls.p2 - DVec3::from_vec3(next)).normalized();
        state.curve_points[0] =
            (controls.p2 + direction * (overhang_width as f64 * OVERHANG_FACTOR as f64)).to_vec3();
    } else if final_residue {
        segment_count = ((linear_segments as f64 * (1.0 - shift)).round() as usize).max(1);
        let previous = state.curve_points[segment_count.saturating_sub(1)];
        let direction = (controls.p2 - DVec3::from_vec3(previous)).normalized();
        state.curve_points[segment_count] =
            (controls.p2 + direction * (overhang_width as f64 * OVERHANG_FACTOR as f64)).to_vec3();
    }
    if !sizes_interpolated {
        interpolate_sizes(
            state, widths[0], widths[1], widths[2], heights[0], heights[1], heights[2], shift,
        );
    }

    samples.centers.clear();
    samples.normals.clear();
    samples.binormals.clear();
    samples.widths.clear();
    samples.heights.clear();
    let sample_count = segment_count + 1;
    samples.centers.reserve(sample_count);
    samples.normals.reserve(sample_count);
    samples.binormals.reserve(sample_count);
    samples.widths.reserve(sample_count);
    samples.heights.reserve(sample_count);

    for i in 0..=segment_count {
        samples.centers.push(state.curve_points[i]);
        if swap_normal_binormal {
            samples.normals.push(state.binormal_vectors[i] * -1.0);
            samples.binormals.push(state.normal_vectors[i]);
        } else {
            samples.normals.push(state.normal_vectors[i]);
            samples.binormals.push(state.binormal_vectors[i]);
        }
        samples.widths.push(state.width_values[i]);
        samples.heights.push(state.height_values[i]);
    }

    samples
}

#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn curve_segment_samples(
    controls: &CurveSegmentControls,
    widths: [f32; 3],
    heights: [f32; 3],
    tension: f64,
    shift: f64,
    overhang_width: f32,
    initial: bool,
    final_residue: bool,
    copy_initial_size_window: bool,
    swap_normal_binormal: bool,
    linear_segments: usize,
) -> CurveSamples {
    let mut scratch = CurveSegmentScratch::default();
    curve_segment_samples_into(
        &mut scratch,
        controls,
        widths,
        heights,
        tension,
        shift,
        overhang_width,
        initial,
        final_residue,
        copy_initial_size_window,
        swap_normal_binormal,
        linear_segments,
    )
    .clone()
}

#[allow(clippy::too_many_arguments)]
fn add_profile_tube_cap(
    mesh: &mut Mesh,
    centers: &[Vec3],
    normals: &[Vec3],
    binormals: &[Vec3],
    widths: &[f32],
    heights: &[f32],
    radial_segments: usize,
    profile: TubeProfile,
    sample_index: usize,
    start: bool,
    round_cap: bool,
) {
    let u = normals[sample_index];
    let v = binormals[sample_index];
    let center = centers[sample_index];
    let width = if round_cap { 0.0 } else { widths[sample_index] };
    let mut height = if round_cap {
        0.0
    } else {
        heights[sample_index]
    };
    let q1 = (radial_segments as f32 / 4.0).round() as usize;
    let q3 = q1 * 3;
    let rounded = profile == TubeProfile::Rounded && height > width;
    if rounded {
        height -= width;
    }
    let (cos_values, sin_values) =
        molstar_tube_cos_sin(radial_segments, profile == TubeProfile::Rounded);
    let normal = if start { v.cross(u) } else { u.cross(v) };
    let center_index = mesh.vertices.len();
    mesh.vertices.push(center);
    mesh.normals.push(normal);
    let ring_start = mesh.vertices.len();
    for j in 0..radial_segments {
        let cos = cos_values[j];
        let sin = sin_values[j];
        let vertex = if rounded {
            let h = if j < q1 || j >= q3 { height } else { -height };
            molstar_vec3_add3_scaled_then_scale_and_add(
                center,
                u,
                v,
                width as f64 * cos,
                width as f64 * sin,
                h as f64,
            )
        } else {
            molstar_vec3_add3_scaled(center, u, v, height as f64 * cos, width as f64 * sin)
        };
        mesh.vertices.push(vertex);
        mesh.normals.push(normal);
        let next = (j + 1) % radial_segments;
        if start {
            mesh.faces.push(Face {
                a: ring_start + next,
                b: ring_start + j,
                c: center_index,
            });
        } else {
            mesh.faces.push(Face {
                a: ring_start + j,
                b: ring_start + next,
                c: center_index,
            });
        }
    }
}

fn sample_oriented_profile_path(
    centers: &[Vec3],
    normals: &[Vec3],
    subdivisions: usize,
    width: f32,
    height: f32,
) -> CurveSamples {
    let mut out = CurveSamples {
        centers: Vec::new(),
        normals: Vec::new(),
        binormals: Vec::new(),
        widths: Vec::new(),
        heights: Vec::new(),
    };
    if centers.len() < 2 || centers.len() != normals.len() || subdivisions == 0 {
        return out;
    }

    let mut previous_normal: Option<Vec3> = None;
    let tension = 0.9;
    for i in 0..centers.len() - 1 {
        let p0 = if i == 0 { centers[i] } else { centers[i - 1] };
        let p1 = centers[i];
        let p2 = centers[i + 1];
        let p3 = if i + 2 < centers.len() {
            centers[i + 2]
        } else {
            centers[i + 1]
        };
        for step in 0..subdivisions {
            let t = step as f32 / subdivisions as f32;
            let center = molstar_spline(p0, p1, p2, p3, t, tension);
            let tangent = spline_tangent(p0, p1, p2, p3, t, tension);
            let target = slerp_unit(normals[i], match_direction(normals[i + 1], normals[i]), t);
            let mut normal = orthogonalize(tangent, target);
            if let Some(previous) = previous_normal {
                normal = match_direction(normal, previous);
            }
            out.centers.push(center);
            out.normals.push(normal);
            out.binormals.push(tangent.cross(normal).normalized());
            out.widths.push(width);
            out.heights.push(height);
            previous_normal = Some(normal);
        }
    }

    let last = centers.len() - 1;
    let p0 = if centers.len() > 2 {
        centers[last - 2]
    } else {
        centers[last - 1]
    };
    let tangent = spline_tangent(
        p0,
        centers[last - 1],
        centers[last],
        centers[last],
        1.0,
        tension,
    );
    let mut normal = orthogonalize(tangent, normals[last]);
    if let Some(previous) = previous_normal {
        normal = match_direction(normal, previous);
    }
    out.centers.push(centers[last]);
    out.normals.push(normal);
    out.binormals.push(tangent.cross(normal).normalized());
    out.widths.push(width);
    out.heights.push(height);
    smooth_profile_normals(&mut out);
    out
}

fn smooth_profile_normals(samples: &mut CurveSamples) {
    if samples.centers.len() < 3 || samples.centers.len() != samples.normals.len() {
        return;
    }
    let original = samples.normals.clone();
    for i in 1..samples.normals.len() - 1 {
        let tangent = path_tangent(&samples.centers, i);
        let blended = (original[i - 1] + original[i] + original[i + 1]) / 3.0;
        let mut normal = orthogonalize(tangent, blended);
        normal = match_direction(normal, samples.normals[i - 1]);
        samples.normals[i] = normal;
        samples.binormals[i] = tangent.cross(normal).normalized();
    }
}

fn spline_tangent(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32, tension: f32) -> Vec3 {
    let eps = 0.01;
    let before = molstar_spline(p0, p1, p2, p3, (t - eps).clamp(0.0, 1.0), tension);
    let after = molstar_spline(p0, p1, p2, p3, (t + eps).clamp(0.0, 1.0), tension);
    (after - before).normalized()
}

#[allow(dead_code)]
fn sample_curve_segment_path(
    centers: &[Vec3],
    normals: &[Vec3],
    linear_segments: usize,
    tension: f64,
    shift: f64,
    width: f32,
    height: f32,
) -> CurveSamples {
    let mut out = CurveSamples {
        centers: Vec::new(),
        normals: Vec::new(),
        binormals: Vec::new(),
        widths: Vec::new(),
        heights: Vec::new(),
    };
    if centers.len() < 2 || centers.len() != normals.len() || linear_segments == 0 {
        return out;
    }

    for i in 0..centers.len() {
        let controls = CurveSegmentControls {
            sec_struc_first: i == 0,
            sec_struc_last: i + 1 == centers.len(),
            p0: DVec3::from_vec3(curve_control_point(centers, i as isize - 2)),
            p1: DVec3::from_vec3(curve_control_point(centers, i as isize - 1)),
            p2: DVec3::from_vec3(centers[i]),
            p3: DVec3::from_vec3(curve_control_point(centers, i as isize + 1)),
            p4: DVec3::from_vec3(curve_control_point(centers, i as isize + 2)),
            d12: DVec3::from_vec3(normals[i]),
            d23: DVec3::from_vec3(normals[(i + 1).min(normals.len() - 1)]),
        };
        let mut state = CurveSegmentState::new(linear_segments);
        interpolate_curve_segment(&mut state, &controls, tension, shift);
        interpolate_sizes(
            &mut state, width, width, width, height, height, height, shift,
        );

        let (start, end) = if centers.len() == 1 {
            (0, linear_segments)
        } else if i == 0 {
            let segment_count = ((linear_segments as f64 * shift).round() as usize).max(1);
            let start = linear_segments - segment_count;
            let away = (controls.p2
                - DVec3::from_vec3(state.curve_points[(start + 1).min(linear_segments)]))
            .normalized();
            state.curve_points[start] = (controls.p2 + away * (width as f64 * 2.0)).to_vec3();
            (start, linear_segments)
        } else if i + 1 == centers.len() {
            let segment_count = ((linear_segments as f64 * (1.0 - shift)).round() as usize).max(1);
            let away = (controls.p2
                - DVec3::from_vec3(state.curve_points[segment_count.saturating_sub(1)]))
            .normalized();
            state.curve_points[segment_count] =
                (controls.p2 + away * (width as f64 * 2.0)).to_vec3();
            (0, segment_count)
        } else {
            (0, linear_segments)
        };

        for j in start..=end {
            if !out.centers.is_empty() && j == start {
                continue;
            }
            out.centers.push(state.curve_points[j]);
            out.normals.push(state.normal_vectors[j]);
            out.binormals.push(state.binormal_vectors[j]);
            out.widths.push(state.width_values[j]);
            out.heights.push(state.height_values[j]);
        }
    }

    out
}

fn curve_control_point(points: &[Vec3], index: isize) -> Vec3 {
    if points.is_empty() {
        return Vec3::default();
    }
    if index >= 0 && (index as usize) < points.len() {
        return points[index as usize];
    }
    if points.len() == 1 {
        return points[0];
    }
    if index < 0 {
        let delta = points[0] - points[1];
        points[0] + delta * (-index as f32)
    } else {
        let last = points.len() - 1;
        let delta = points[last] - points[last - 1];
        points[last] + delta * (index as usize - last) as f32
    }
}

#[derive(Clone, Debug)]
pub(crate) struct CurveSegmentState {
    pub(crate) curve_points: Vec<Vec3>,
    pub(crate) tangent_vectors: Vec<Vec3>,
    pub(crate) normal_vectors: Vec<Vec3>,
    pub(crate) binormal_vectors: Vec<Vec3>,
    pub(crate) width_values: Vec<f32>,
    pub(crate) height_values: Vec<f32>,
    pub(crate) linear_segments: usize,
}

impl CurveSegmentState {
    pub(crate) fn new(linear_segments: usize) -> Self {
        let n = linear_segments + 1;
        Self {
            curve_points: vec![Vec3::default(); n],
            tangent_vectors: vec![Vec3::default(); n],
            normal_vectors: vec![Vec3::default(); n],
            binormal_vectors: vec![Vec3::default(); n],
            width_values: vec![0.0; n],
            height_values: vec![0.0; n],
            linear_segments,
        }
    }

    fn prepare(&mut self, linear_segments: usize) {
        let n = linear_segments + 1;
        self.curve_points.resize(n, Vec3::default());
        self.tangent_vectors.resize(n, Vec3::default());
        self.normal_vectors.resize(n, Vec3::default());
        self.binormal_vectors.resize(n, Vec3::default());
        self.width_values.resize(n, 0.0);
        self.height_values.resize(n, 0.0);
        self.linear_segments = linear_segments;
    }
}

#[derive(Clone, Debug)]
pub(crate) struct CurveSegmentControls {
    pub(crate) sec_struc_first: bool,
    pub(crate) sec_struc_last: bool,
    pub(crate) p0: DVec3,
    pub(crate) p1: DVec3,
    pub(crate) p2: DVec3,
    pub(crate) p3: DVec3,
    pub(crate) p4: DVec3,
    pub(crate) d12: DVec3,
    pub(crate) d23: DVec3,
}

pub(crate) fn interpolate_curve_segment(
    state: &mut CurveSegmentState,
    controls: &CurveSegmentControls,
    tension: f64,
    shift: f64,
) {
    interpolate_points_and_tangents(state, controls, tension, shift);
    interpolate_curve_normals(state, controls);
}

pub(crate) fn interpolate_points_and_tangents(
    state: &mut CurveSegmentState,
    controls: &CurveSegmentControls,
    tension: f64,
    shift: f64,
) {
    let shift1 = 1.0 - shift;
    let tension_beg = if controls.sec_struc_first {
        0.5
    } else {
        tension
    };
    let tension_end = if controls.sec_struc_last {
        0.5
    } else {
        tension
    };

    for j in 0..=state.linear_segments {
        let t = j as f64 / state.linear_segments as f64;
        let (point, tan_a, tan_b) = if t < shift1 {
            let te = lerp_f64(tension_beg, tension, t);
            (
                molstar_spline_f64(
                    controls.p0,
                    controls.p1,
                    controls.p2,
                    controls.p3,
                    t + shift,
                    te,
                ),
                molstar_spline_f64(
                    controls.p0,
                    controls.p1,
                    controls.p2,
                    controls.p3,
                    t + shift + 0.01,
                    tension_beg,
                ),
                molstar_spline_f64(
                    controls.p0,
                    controls.p1,
                    controls.p2,
                    controls.p3,
                    t + shift - 0.01,
                    tension_beg,
                ),
            )
        } else {
            let te = lerp_f64(tension, tension_end, t);
            (
                molstar_spline_f64(
                    controls.p1,
                    controls.p2,
                    controls.p3,
                    controls.p4,
                    t - shift1,
                    te,
                ),
                molstar_spline_f64(
                    controls.p1,
                    controls.p2,
                    controls.p3,
                    controls.p4,
                    t - shift1 + 0.01,
                    te,
                ),
                molstar_spline_f64(
                    controls.p1,
                    controls.p2,
                    controls.p3,
                    controls.p4,
                    t - shift1 - 0.01,
                    te,
                ),
            )
        };
        state.curve_points[j] = point.to_vec3();
        state.tangent_vectors[j] = normalize_dvec3_molstar(tan_a - tan_b).to_vec3();
    }
}

pub(crate) fn interpolate_curve_normals(
    state: &mut CurveSegmentState,
    controls: &CurveSegmentControls,
) {
    let n = state.curve_points.len();
    if n == 0 {
        return;
    }
    let first_tangent = DVec3::from_vec3(state.tangent_vectors[0]);
    let last_tangent = DVec3::from_vec3(state.tangent_vectors[n - 1]);
    let first_normal = orthogonalize_f64(first_tangent, controls.d12);
    let last_normal =
        match_direction_f64(orthogonalize_f64(last_tangent, controls.d23), first_normal);

    let mut previous_normal = first_normal;
    let n1 = n - 1;
    for i in 0..n {
        let j = smoothstep_f64(0.0, n1 as f64, i as f64) * n1 as f64;
        let t = if i == 0 { 0.0 } else { 1.0 / (n as f64 - j) };
        let tangent = DVec3::from_vec3(state.tangent_vectors[i]);
        let normal = orthogonalize_f64(tangent, slerp_unit_f64(previous_normal, last_normal, t));
        state.normal_vectors[i] = normal.to_vec3();
        previous_normal = normal;
        state.binormal_vectors[i] =
            normalize_dvec3_molstar(cross_dvec3_molstar(tangent, normal)).to_vec3();
    }

    for i in 1..n1 {
        let tangent = DVec3::from_vec3(state.tangent_vectors[i]);
        let normal = (DVec3::from_vec3(state.normal_vectors[i - 1])
            + (DVec3::from_vec3(state.normal_vectors[i + 1])
                + DVec3::from_vec3(state.normal_vectors[i])))
            * (1.0 / 3.0);
        state.normal_vectors[i] = normal.to_vec3();
        state.binormal_vectors[i] =
            normalize_dvec3_molstar(cross_dvec3_molstar(tangent, normal)).to_vec3();
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn interpolate_sizes(
    state: &mut CurveSegmentState,
    w0: f32,
    w1: f32,
    w2: f32,
    h0: f32,
    h1: f32,
    h2: f32,
    shift: f64,
) {
    let shift1 = 1.0 - shift;
    for i in 0..=state.linear_segments {
        let t = i as f64 / state.linear_segments as f64;
        if t < shift1 {
            state.width_values[i] = lerp_f64(w0 as f64, w1 as f64, t + shift) as f32;
            state.height_values[i] = lerp_f64(h0 as f64, h1 as f64, t + shift) as f32;
        } else {
            state.width_values[i] = lerp_f64(w1 as f64, w2 as f64, t - shift1) as f32;
            state.height_values[i] = lerp_f64(h1 as f64, h2 as f64, t - shift1) as f32;
        }
    }
}

fn lerp_f64(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

fn smoothstep_f64(edge0: f64, edge1: f64, x: f64) -> f64 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

fn slerp_unit_f64(a: DVec3, b: DVec3, t: f64) -> DVec3 {
    let dot = a.dot(b).clamp(-1.0, 1.0);
    let theta = dot.acos() * t;
    let rel = normalize_dvec3_molstar(DVec3::new(
        b.x + a.x * -dot,
        b.y + a.y * -dot,
        b.z + a.z * -dot,
    ));
    let cos = theta.cos();
    let sin = theta.sin();
    DVec3::new(
        a.x * cos + rel.x * sin,
        a.y * cos + rel.y * sin,
        a.z * cos + rel.z * sin,
    )
}

fn match_direction_f64(value: DVec3, reference: DVec3) -> DVec3 {
    if value.dot(reference) > 0.0 {
        value
    } else {
        value * -1.0
    }
}

fn orthogonalize_f64(axis: DVec3, direction: DVec3) -> DVec3 {
    let mut out = normalize_dvec3_molstar(cross_dvec3_molstar(
        cross_dvec3_molstar(axis, direction),
        axis,
    ));
    if !out.is_zero() {
        return out;
    }
    out = normalize_dvec3_molstar(cross_dvec3_molstar(
        cross_dvec3_molstar(
            axis,
            DVec3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
        ),
        axis,
    ));
    if !out.is_zero() {
        return out;
    }
    out = normalize_dvec3_molstar(cross_dvec3_molstar(
        cross_dvec3_molstar(
            axis,
            DVec3 {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
        ),
        axis,
    ));
    if !out.is_zero() {
        return out;
    }
    let fallback = normalize_dvec3_molstar(direction);
    if !fallback.is_zero() {
        fallback
    } else {
        DVec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

fn cross_dvec3_molstar(a: DVec3, b: DVec3) -> DVec3 {
    let ax = a.x;
    let ay = a.y;
    let az = a.z;
    let bx = b.x;
    let by = b.y;
    let bz = b.z;
    DVec3::new(ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx)
}

fn normalize_dvec3_molstar(value: DVec3) -> DVec3 {
    let x = value.x;
    let y = value.y;
    let z = value.z;
    let len = x * x + y * y + z * z;
    if len > 0.0 {
        let scale = 1.0 / len.sqrt();
        DVec3::new(value.x * scale, value.y * scale, value.z * scale)
    } else {
        value
    }
}

fn orthogonalize(axis: Vec3, direction: Vec3) -> Vec3 {
    let mut out = normalize_molstar(axis.cross(direction).cross(axis));
    if out != Vec3::default() {
        return out;
    }
    out = normalize_molstar(
        axis.cross(Vec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        })
        .cross(axis),
    );
    if out != Vec3::default() {
        return out;
    }
    out = normalize_molstar(
        axis.cross(Vec3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        })
        .cross(axis),
    );
    if out != Vec3::default() {
        return out;
    }
    let fallback = normalize_molstar(direction);
    if fallback != Vec3::default() {
        fallback
    } else {
        Vec3 {
            x: 1.0,
            y: 0.0,
            z: 0.0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct MolstarRotationD {
    m: [[f64; 3]; 3],
}

impl MolstarRotationD {
    fn from_up_to_mat_dir(target: DVec3, match_dir: bool) -> Self {
        let source = if !match_dir || target.y > 0.0 {
            DVec3::new(0.0, 1.0, 0.0)
        } else {
            DVec3::new(0.0, -1.0, 0.0)
        };
        let denominator = (source.squared_length() * target.squared_length()).sqrt();
        let angle = if denominator == 0.0 {
            std::f64::consts::FRAC_PI_2
        } else {
            (source.dot(target) / denominator).clamp(-1.0, 1.0).acos()
        };
        if angle.abs() < 0.0001 {
            return Self::identity();
        }
        if (angle - std::f64::consts::PI).abs() < f64::EPSILON {
            let rotation_axis = if source.x.abs() < 0.9 {
                DVec3::new(1.0, 0.0, 0.0)
            } else {
                DVec3::new(0.0, 0.0, 1.0)
            };
            return Self::from_axis_angle(rotation_axis, std::f64::consts::PI);
        }
        Self::from_axis_angle(source.cross(target), angle)
    }

    fn identity() -> Self {
        Self {
            m: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    fn from_axis_angle(axis: DVec3, angle: f64) -> Self {
        let axis = axis.normalized();
        let (s, c) = angle.sin_cos();
        let oc = 1.0 - c;
        Self {
            m: [
                [
                    oc * axis.x * axis.x + c,
                    oc * axis.x * axis.y - axis.z * s,
                    oc * axis.x * axis.z + axis.y * s,
                ],
                [
                    oc * axis.y * axis.x + axis.z * s,
                    oc * axis.y * axis.y + c,
                    oc * axis.y * axis.z - axis.x * s,
                ],
                [
                    oc * axis.z * axis.x - axis.y * s,
                    oc * axis.z * axis.y + axis.x * s,
                    oc * axis.z * axis.z + c,
                ],
            ],
        }
    }

    fn apply(self, v: DVec3) -> DVec3 {
        DVec3::new(
            self.m[0][0] * v.x + self.m[0][1] * v.y + self.m[0][2] * v.z,
            self.m[1][0] * v.x + self.m[1][1] * v.y + self.m[1][2] * v.z,
            self.m[2][0] * v.x + self.m[2][1] * v.y + self.m[2][2] * v.z,
        )
    }
}

fn molstar_transform_cylinder_position_d(
    center: DVec3,
    rotation: MolstarRotationD,
    length: f64,
    vertex: Vec3,
) -> DVec3 {
    let x = vertex.x as f64;
    let y = vertex.y as f64;
    let z = vertex.z as f64;
    DVec3::new(
        rotation.m[0][0] * x + (rotation.m[0][1] * length) * y + rotation.m[0][2] * z + center.x,
        rotation.m[1][0] * x + (rotation.m[1][1] * length) * y + rotation.m[1][2] * z + center.y,
        rotation.m[2][0] * x + (rotation.m[2][1] * length) * y + rotation.m[2][2] * z + center.z,
    )
}

fn match_direction(value: Vec3, reference: Vec3) -> Vec3 {
    if value.dot(reference) > 0.0 {
        value
    } else {
        value * -1.0
    }
}

fn slerp_unit(a: Vec3, b: Vec3, t: f32) -> Vec3 {
    let dot = a.dot(b).clamp(-1.0, 1.0);
    let theta = dot.acos() * t;
    let rel = normalize_molstar(b - a * dot);
    a * theta.cos() + rel * theta.sin()
}

fn normalize_molstar(value: Vec3) -> Vec3 {
    let x = value.x as f64;
    let y = value.y as f64;
    let z = value.z as f64;
    let len_sq = x * x + y * y + z * z;
    if len_sq > 0.0 {
        let scale = 1.0 / len_sq.sqrt();
        Vec3::new((x * scale) as f32, (y * scale) as f32, (z * scale) as f32)
    } else {
        value
    }
}

fn normalize_molstar_f64(value: Vec3) -> Vec3 {
    let x = value.x as f64;
    let y = value.y as f64;
    let z = value.z as f64;
    let len_sq = x * x + y * y + z * z;
    if len_sq > 0.0 {
        let scale = 1.0 / len_sq.sqrt();
        Vec3::new((x * scale) as f32, (y * scale) as f32, (z * scale) as f32)
    } else {
        value
    }
}

fn molstar_vec3_add3_scaled(center: Vec3, a: Vec3, b: Vec3, sa: f64, sb: f64) -> Vec3 {
    Vec3::new(
        (a.x as f64 * sa + b.x as f64 * sb + center.x as f64) as f32,
        (a.y as f64 * sa + b.y as f64 * sb + center.y as f64) as f32,
        (a.z as f64 * sa + b.z as f64 * sb + center.z as f64) as f32,
    )
}

fn molstar_vec3_add3_scaled_then_scale_and_add(
    center: Vec3,
    a: Vec3,
    b: Vec3,
    sa: f64,
    sb: f64,
    scale: f64,
) -> Vec3 {
    Vec3::new(
        (a.x as f64 * sa + b.x as f64 * sb + center.x as f64 + a.x as f64 * scale) as f32,
        (a.y as f64 * sa + b.y as f64 * sb + center.y as f64 + a.y as f64 * scale) as f32,
        (a.z as f64 * sa + b.z as f64 * sb + center.z as f64 + a.z as f64 * scale) as f32,
    )
}

fn molstar_vec3_normalize_add2_scaled(a: Vec3, b: Vec3, sa: f64, sb: f64) -> Vec3 {
    let x = a.x as f64 * sa + b.x as f64 * sb;
    let y = a.y as f64 * sa + b.y as f64 * sb;
    let z = a.z as f64 * sa + b.z as f64 * sb;
    let len_sq = x * x + y * y + z * z;
    if len_sq > 0.0 {
        let scale = 1.0 / len_sq.sqrt();
        Vec3::new((x * scale) as f32, (y * scale) as f32, (z * scale) as f32)
    } else {
        Vec3::new(x as f32, y as f32, z as f32)
    }
}

fn path_tangent(points: &[Vec3], i: usize) -> Vec3 {
    if points.len() < 2 {
        return Vec3::default();
    }
    let tangent = if i == 0 {
        points[1] - points[0]
    } else if i + 1 == points.len() {
        points[i] - points[i - 1]
    } else {
        points[i + 1] - points[i - 1]
    };
    tangent.normalized()
}

pub(super) fn fallback_side(tangent: Vec3, previous_side: Option<Vec3>) -> Vec3 {
    if let Some(previous) = previous_side {
        let side = (previous - tangent * previous.dot(tangent)).normalized();
        if side.length() > 0.000_001 {
            return side;
        }
    }
    let helper = if tangent.z.abs() < 0.9 {
        Vec3 {
            x: 0.0,
            y: 0.0,
            z: 1.0,
        }
    } else {
        Vec3 {
            x: 0.0,
            y: 1.0,
            z: 0.0,
        }
    };
    tangent.cross(helper).normalized()
}

pub(crate) fn add_ribbon(
    mesh: &mut Mesh,
    points: &[Vec3],
    width: f32,
    thickness: f32,
    linear_segments: usize,
) {
    if points.len() < 2 {
        return;
    }
    let samples = sample_sheet_path(points, linear_segments.max(1), width, thickness);
    add_ribbon_samples(mesh, &samples, 0.0);
}

fn add_ribbon_samples(mesh: &mut Mesh, samples: &CurveSamples, arrow_height: f32) {
    let linear_segments = samples.centers.len().saturating_sub(1);
    if linear_segments == 0 {
        return;
    }
    let base = mesh.vertices.len();

    for i in 0..samples.centers.len() {
        let center = samples.centers[i];
        let actual_height = if arrow_height == 0.0 {
            samples.heights[i]
        } else {
            arrow_height * (1.0 - i as f32 / linear_segments as f32)
        };
        let vertical = samples.normals[i] * actual_height;
        let torsion = samples.binormals[i];

        mesh.vertices.push(center + vertical);
        mesh.normals.push(torsion * -1.0);
        mesh.vertices.push(center - vertical);
        mesh.normals.push(torsion * -1.0);
        mesh.vertices.push(center + vertical);
        mesh.normals.push(torsion);
        mesh.vertices.push(center - vertical);
        mesh.normals.push(torsion);
    }

    for i in 0..linear_segments {
        let a = base + i * 4;
        let b = a + 4;
        mesh.faces.push(Face {
            a,
            b: b + 1,
            c: a + 1,
        });
        mesh.faces.push(Face { a, b, c: b + 1 });
        mesh.faces.push(Face {
            a: a + 3,
            b: b + 3,
            c: a + 2,
        });
        mesh.faces.push(Face {
            a: a + 2,
            b: b + 3,
            c: b + 2,
        });
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn add_sheet(
    mesh: &mut Mesh,
    points: &[Vec3],
    width: f32,
    thickness: f32,
    arrow_height: f32,
    start_cap: bool,
    end_cap: bool,
    linear_segments: usize,
) {
    if points.len() < 2 {
        return;
    }
    let samples = sample_sheet_path(points, linear_segments.max(1), width, thickness);
    add_sheet_samples(mesh, &samples, arrow_height, start_cap, end_cap);
}

fn add_sheet_samples(
    mesh: &mut Mesh,
    samples: &CurveSamples,
    arrow_height: f32,
    start_cap: bool,
    end_cap: bool,
) {
    let linear_segments = samples.centers.len().saturating_sub(1);
    if linear_segments == 0 {
        return;
    }

    let base = mesh.vertices.len();
    let arrow_height = arrow_height.max(0.0);
    let offset_length = if arrow_height > 0.0 {
        let length = samples.centers[0].distance(*samples.centers.last().unwrap());
        if length > 0.000_001 {
            arrow_height / length
        } else {
            0.0
        }
    } else {
        0.0
    };

    for i in 0..samples.centers.len() {
        let center = samples.centers[i];
        let normal = samples.normals[i];
        let binormal = samples.binormals[i];
        let width = samples.widths[i];
        let height = samples.heights[i];
        let actual_height = if arrow_height == 0.0 {
            height
        } else {
            arrow_height * (1.0 - i as f32 / linear_segments as f32)
        };
        let vertical = normal * actual_height;
        let horizontal = binormal * width;
        let normal_offset = if arrow_height > 0.0 {
            normal.cross(binormal) * offset_length
        } else {
            Vec3::default()
        };
        let torsion = binormal;

        let p_top_right = center + horizontal + vertical;
        let p_top_left = center - horizontal + vertical;
        let p_bottom_left = center - horizontal - vertical;
        let p_bottom_right = center + horizontal - vertical;

        mesh.vertices.push(p_top_right);
        mesh.normals.push(normal + normal_offset);
        mesh.vertices.push(p_top_left);
        mesh.normals.push(normal + normal_offset);
        mesh.vertices.push(p_top_left);
        mesh.normals.push(torsion * -1.0);
        mesh.vertices.push(p_bottom_left);
        mesh.normals.push(torsion * -1.0);
        mesh.vertices.push(p_bottom_left);
        mesh.normals.push(normal * -1.0 + normal_offset);
        mesh.vertices.push(p_bottom_right);
        mesh.normals.push(normal * -1.0 + normal_offset);
        mesh.vertices.push(p_bottom_right);
        mesh.normals.push(torsion);
        mesh.vertices.push(p_top_right);
        mesh.normals.push(torsion);
    }

    for i in 0..linear_segments {
        for j in 0..2 {
            let a = base + i * 8 + 2 * j;
            let b = base + i * 8 + 2 * j + 1;
            let c = base + (i + 1) * 8 + 2 * j + 1;
            let d = base + (i + 1) * 8 + 2 * j;
            mesh.faces.push(Face { a, b: c, c: b });
            mesh.faces.push(Face { a, b: d, c });
        }
        for j in 2..4 {
            let a = base + i * 8 + 2 * j;
            let b = base + i * 8 + 2 * j + 1;
            let c = base + (i + 1) * 8 + 2 * j + 1;
            let d = base + (i + 1) * 8 + 2 * j;
            mesh.faces.push(Face { a, b: d, c: b });
            mesh.faces.push(Face { a: d, b: c, c: b });
        }
    }

    if start_cap {
        let width = samples.widths[0];
        let height = samples.heights[0];
        let h = if arrow_height == 0.0 {
            height
        } else {
            arrow_height
        };
        add_sheet_cap(mesh, samples, 0, width, h, h, false);
    } else if arrow_height > 0.0 {
        let width = samples.widths[0];
        let height = samples.heights[0];
        add_sheet_cap(mesh, samples, 0, width, arrow_height, -height, false);
        add_sheet_cap(mesh, samples, 0, width, -arrow_height, height, false);
    }
    if end_cap && arrow_height == 0.0 {
        let width = samples.widths[linear_segments];
        let height = samples.heights[linear_segments];
        add_sheet_cap(mesh, samples, linear_segments, width, height, height, true);
    }
}

fn sample_sheet_path(
    points: &[Vec3],
    subdivisions: usize,
    width: f32,
    height: f32,
) -> CurveSamples {
    let centers = sample_path(points, subdivisions);
    let mut out = CurveSamples {
        centers,
        normals: Vec::new(),
        binormals: Vec::new(),
        widths: Vec::new(),
        heights: Vec::new(),
    };
    let mut previous_side = None;
    for i in 0..out.centers.len() {
        let tangent = path_tangent(&out.centers, i);
        let mut side = fallback_side(tangent, previous_side);
        if let Some(previous) = previous_side {
            side = match_direction(side, previous);
        }
        let normal = side.cross(tangent).normalized();
        out.normals.push(normal);
        out.binormals.push(side);
        out.widths.push(width);
        out.heights.push(height);
        previous_side = Some(side);
    }
    out
}

fn add_sheet_cap(
    mesh: &mut Mesh,
    samples: &CurveSamples,
    sample_index: usize,
    width: f32,
    left_height: f32,
    right_height: f32,
    flip: bool,
) {
    let center = samples.centers[sample_index];
    let normal = samples.normals[sample_index];
    let binormal = samples.binormals[sample_index];
    let vertical_left = normal * left_height;
    let vertical_right = normal * right_height;
    let horizontal = binormal * width;
    let cap_normal = binormal.cross(normal);

    let p1 = center + horizontal + vertical_right;
    let p2 = center + horizontal - vertical_left;
    let p3 = center - horizontal - vertical_left;
    let p4 = center - horizontal + vertical_right;
    let points = if left_height < right_height {
        [p4, p3, p2, p1]
    } else {
        [p1, p2, p3, p4]
    };
    let base = mesh.vertices.len();
    for point in points {
        mesh.vertices.push(point);
        mesh.normals
            .push(if flip { cap_normal * -1.0 } else { cap_normal });
    }
    if flip {
        mesh.faces.push(Face {
            a: base,
            b: base + 1,
            c: base + 2,
        });
        mesh.faces.push(Face {
            a: base + 2,
            b: base + 3,
            c: base,
        });
    } else {
        mesh.faces.push(Face {
            a: base + 2,
            b: base + 1,
            c: base,
        });
        mesh.faces.push(Face {
            a: base,
            b: base + 3,
            c: base + 2,
        });
    }
}
