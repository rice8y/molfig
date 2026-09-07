use std::ops::Range;

use crate::mesh::{NativePrimitive, NativeRenderGeometry, NativeRenderScene};
use crate::model::{Mesh, MeshMaterial, Vec3};

use super::camera::{CameraState, ProjectedVertex};
use super::color::{color_f32, quantize};
use super::framebuffer::Framebuffer;
use super::options::RendererOptions;
use super::postprocessing::angle_metal_math::{divide, fused_multiply_add, inverse_sqrt};
use super::postprocessing::fog_opaque_fragment;
use super::shading::{shade_material_base_color_linear, shade_unlit_material_linear};
use super::style::ResolvedStyle;
use super::transparency::WboitTarget;

enum RasterTarget<'a> {
    Opaque,
    Wboit(&'a mut WboitTarget),
}

pub(super) struct RasterPass<'a> {
    target: RasterTarget<'a>,
    object_alpha: Option<f64>,
}

impl<'a> RasterPass<'a> {
    pub(super) fn opaque() -> Self {
        Self {
            target: RasterTarget::Opaque,
            object_alpha: None,
        }
    }

    pub(super) fn wboit(target: &'a mut WboitTarget) -> Self {
        Self {
            target: RasterTarget::Wboit(target),
            object_alpha: None,
        }
    }

    fn alpha(&self, material: MeshMaterial) -> f64 {
        self.object_alpha
            .unwrap_or(f64::from(material.alpha_tenths) / 10.0)
    }

    fn accepts(&self, material: MeshMaterial) -> bool {
        let alpha = self.alpha(material);
        match self.target {
            RasterTarget::Opaque => alpha >= 1.0,
            RasterTarget::Wboit(_) => alpha > 0.0 && alpha < 1.0,
        }
    }

    fn is_opaque(&self) -> bool {
        matches!(self.target, RasterTarget::Opaque)
    }
}

pub(super) fn rasterize_scene(
    scene: &NativeRenderScene,
    camera: &CameraState,
    renderer: &RendererOptions,
    style: ResolvedStyle,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    for object in &scene.render_objects {
        if !object.visible {
            continue;
        }
        // Mol* selects render passes in JavaScript-number precision, then
        // uploads the same clamped alpha * alphaFactor as a float32 uAlpha.
        // Export materials are not the renderer's uniform state.
        pass.object_alpha = Some((object.alpha * object.alpha_factor).clamp(0.0, 1.0));
        match &object.geometry {
            NativeRenderGeometry::Mesh { ranges } => {
                for range in ranges {
                    rasterize_mesh_faces(
                        &scene.mesh,
                        range.faces.clone(),
                        object.double_sided,
                        object.flip_sided,
                        camera,
                        renderer,
                        style,
                        framebuffer,
                        pass,
                    );
                }
            }
            NativeRenderGeometry::Primitives { indices, .. } => {
                for &index in indices {
                    if let Some(primitive) = scene.primitives.get(index) {
                        rasterize_primitive(primitive, camera, renderer, style, framebuffer, pass);
                    }
                }
            }
        }
    }
    pass.object_alpha = None;
}

#[allow(clippy::too_many_arguments)]
fn rasterize_mesh_faces(
    mesh: &Mesh,
    face_range: Range<usize>,
    double_sided: bool,
    flip_sided: bool,
    camera: &CameraState,
    renderer: &RendererOptions,
    style: ResolvedStyle,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    for face_index in face_range {
        let Some(face) = mesh.faces.get(face_index) else {
            continue;
        };
        let Some(&a_world) = mesh.vertices.get(face.a) else {
            continue;
        };
        let Some(&b_world) = mesh.vertices.get(face.b) else {
            continue;
        };
        let Some(&c_world) = mesh.vertices.get(face.c) else {
            continue;
        };
        let fallback_normal = (b_world - a_world).cross(c_world - a_world).normalized();
        let na = mesh.normals.get(face.a).copied().unwrap_or(fallback_normal);
        let nb = mesh.normals.get(face.b).copied().unwrap_or(fallback_normal);
        let nc = mesh.normals.get(face.c).copied().unwrap_or(fallback_normal);
        // `mesh.vert` normalizes the attribute before applying the normal
        // matrix, then normalizes the transformed value once more. Preserve
        // those per-vertex varyings; normalizing only after barycentric
        // interpolation incorrectly weights non-unit source normals.
        let world_normals = [na, nb, nc].map(angle_normalize);
        let view_normals = world_normals.map(|normal| mesh_normal_to_view(camera, normal));
        let view_positions = [a_world, b_world, c_world].map(|world| camera.view_position(world));
        let material = mesh
            .face_material(face_index)
            .unwrap_or_else(|| MeshMaterial::opaque(0xcccccc));
        if !pass.accepts(material) {
            continue;
        }
        let vertex_colors = [face.a, face.b, face.c].map(|vertex| {
            mesh.vertex_color(vertex)
                .unwrap_or_else(|| color_f32(material.color))
        });

        let direct = [
            camera.project(a_world, framebuffer.width, framebuffer.height),
            camera.project(b_world, framebuffer.width, framebuffer.height),
            camera.project(c_world, framebuffer.width, framebuffer.height),
        ];
        if let [Some(a), Some(b), Some(c)] = direct {
            rasterize_mesh_screen_triangle(
                [a, b, c],
                view_positions,
                world_normals,
                view_normals,
                vertex_colors,
                if camera.orthographic {
                    [1.0; 3]
                } else {
                    [a.depth, b.depth, c.depth]
                },
                double_sided,
                flip_sided,
                material,
                camera,
                renderer,
                style,
                framebuffer,
                pass,
            );
            continue;
        }

        // `CameraState::project` deliberately rejects vertices outside the
        // near/far interval. WebGL does not reject the whole primitive in that
        // case: fixed-function clipping intersects it with all six homogeneous
        // clip planes, and the generated vertices carry the shader varyings.
        let polygon = clip_mesh_triangle_to_view_volume(
            [a_world, b_world, c_world],
            view_positions,
            world_normals,
            view_normals,
            vertex_colors,
            camera,
        );
        if polygon.len() < 3 {
            continue;
        }
        for index in 1..polygon.len() - 1 {
            let vertices = [polygon[0], polygon[index], polygon[index + 1]];
            let (Some(a), Some(b), Some(c)) = (
                project_line_clip_vertex(
                    vertices[0].clip,
                    camera,
                    framebuffer.width,
                    framebuffer.height,
                ),
                project_line_clip_vertex(
                    vertices[1].clip,
                    camera,
                    framebuffer.width,
                    framebuffer.height,
                ),
                project_line_clip_vertex(
                    vertices[2].clip,
                    camera,
                    framebuffer.width,
                    framebuffer.height,
                ),
            ) else {
                continue;
            };
            rasterize_mesh_screen_triangle(
                [a, b, c],
                [
                    vertices[0].view_position,
                    vertices[1].view_position,
                    vertices[2].view_position,
                ],
                [
                    vertices[0].world_normal,
                    vertices[1].world_normal,
                    vertices[2].world_normal,
                ],
                [
                    vertices[0].view_normal,
                    vertices[1].view_normal,
                    vertices[2].view_normal,
                ],
                [vertices[0].color, vertices[1].color, vertices[2].color],
                [vertices[0].clip.w, vertices[1].clip.w, vertices[2].clip.w],
                double_sided,
                flip_sided,
                material,
                camera,
                renderer,
                style,
                framebuffer,
                pass,
            );
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct MeshClipVertex {
    clip: LineClipVertex,
    world: Vec3,
    view_position: Vec3,
    world_normal: Vec3,
    view_normal: Vec3,
    color: [f32; 3],
}

fn clip_mesh_triangle_to_view_volume(
    world: [Vec3; 3],
    view_positions: [Vec3; 3],
    world_normals: [Vec3; 3],
    view_normals: [Vec3; 3],
    colors: [[f32; 3]; 3],
    camera: &CameraState,
) -> Vec<MeshClipVertex> {
    let projection = camera.projection_matrix();
    let mut input = world
        .into_iter()
        .zip(view_positions)
        .zip(world_normals)
        .zip(view_normals)
        .zip(colors)
        .map(
            |((((world, view_position), world_normal), view_normal), color)| MeshClipVertex {
                clip: project_view_to_line_clip(view_position, projection, camera.orthographic),
                world,
                view_position,
                world_normal,
                view_normal,
                color,
            },
        )
        .collect::<Vec<_>>();
    if input
        .iter()
        .any(|vertex| !line_clip_vertex_is_finite(vertex.clip))
    {
        return Vec::new();
    }

    let mut output = Vec::with_capacity(9);
    for plane in 0..6 {
        if input.is_empty() {
            break;
        }
        output.clear();
        let mut previous = *input.last().expect("non-empty clipped mesh polygon");
        let mut previous_distance = line_clip_plane_distance(previous.clip, plane);
        for &current in &input {
            let current_distance = line_clip_plane_distance(current.clip, plane);
            let previous_inside = previous_distance >= 0.0;
            let current_inside = current_distance >= 0.0;
            if previous_inside != current_inside {
                let alpha = previous_distance / (previous_distance - current_distance);
                output.push(interpolate_mesh_clip_vertex(previous, current, alpha));
            }
            if current_inside {
                output.push(current);
            }
            previous = current;
            previous_distance = current_distance;
        }
        std::mem::swap(&mut input, &mut output);
    }
    input
}

fn interpolate_mesh_clip_vertex(
    start: MeshClipVertex,
    end: MeshClipVertex,
    alpha: f32,
) -> MeshClipVertex {
    MeshClipVertex {
        clip: interpolate_line_clip_vertex(start.clip, end.clip, alpha),
        world: interpolate_vec3(start.world, end.world, alpha),
        view_position: interpolate_vec3(start.view_position, end.view_position, alpha),
        world_normal: interpolate_vec3(start.world_normal, end.world_normal, alpha),
        view_normal: interpolate_vec3(start.view_normal, end.view_normal, alpha),
        color: interpolate_color(start.color, end.color, alpha),
    }
}

fn interpolate_vec3(start: Vec3, end: Vec3, alpha: f32) -> Vec3 {
    Vec3::new(
        fused_multiply_add(alpha, end.x - start.x, start.x),
        fused_multiply_add(alpha, end.y - start.y, start.y),
        fused_multiply_add(alpha, end.z - start.z, start.z),
    )
}

fn interpolate_color(start: [f32; 3], end: [f32; 3], alpha: f32) -> [f32; 3] {
    [
        fused_multiply_add(alpha, end[0] - start[0], start[0]),
        fused_multiply_add(alpha, end[1] - start[1], start[1]),
        fused_multiply_add(alpha, end[2] - start[2], start[2]),
    ]
}

#[allow(clippy::too_many_arguments)]
fn rasterize_mesh_screen_triangle(
    [mut a, mut b, mut c]: [ProjectedVertex; 3],
    [a_view, b_view, c_view]: [Vec3; 3],
    [na_world, nb_world, nc_world]: [Vec3; 3],
    [na_view, nb_view, nc_view]: [Vec3; 3],
    [ca, cb, cc]: [[f32; 3]; 3],
    interpolation_w: [f32; 3],
    double_sided: bool,
    flip_sided: bool,
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
    style: ResolvedStyle,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    // Apple GPU rasterization snaps post-viewport vertex coordinates to an
    // 8-bit subpixel grid before coverage and depth interpolation. This is
    // observable in the pinned Mol*/ANGLE path on narrow cartoon triangles.
    for vertex in [&mut a, &mut b, &mut c] {
        vertex.x = snap_to_raster_subpixel(vertex.x);
        vertex.y = snap_to_raster_subpixel(vertex.y);
    }
    let area = edge(a.x, a.y, b.x, b.y, c.x, c.y);
    if !area.is_finite() || area.abs() < 0.000_001 {
        return;
    }
    // Framebuffer rows are stored top-down, so the viewport Y flip maps a
    // WebGL CCW front face to a positive value under this module's negated
    // cross-product edge convention. Mol* disables culling for quality
    // levels whose resolved `uDoubleSided` is true and flips the shaded
    // normal on interior fragments.
    let front_facing = if flip_sided { area < 0.0 } else { area > 0.0 };
    // All current scene visuals use Mol* `transparentBackfaces = off`.
    // Its transparency shader discards interior fragments even when the
    // geometry's uDoubleSided flag disables fixed-function face culling.
    // Applying only uDoubleSided would count each transparent shell twice.
    if !front_facing && (!double_sided || !pass.is_opaque()) {
        return;
    }
    let min_x = a.x.min(b.x).min(c.x).floor().max(0.0) as usize;
    let max_x =
        a.x.max(b.x)
            .max(c.x)
            .ceil()
            .min(framebuffer.width as f32 - 1.0) as usize;
    let min_y = a.y.min(b.y).min(c.y).floor().max(0.0) as usize;
    let max_y =
        a.y.max(b.y)
            .max(c.y)
            .ceil()
            .min(framebuffer.height as f32 - 1.0) as usize;
    if min_x > max_x || min_y > max_y {
        return;
    }
    let inv_area = 1.0 / area;
    let area64 = edge64(a.x, a.y, b.x, b.y, c.x, c.y);
    let inverse_area64 = 1.0 / area64;
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let px = x as f32 + 0.5;
            let py = y as f32 + 0.5;
            let e0 = edge(b.x, b.y, c.x, c.y, px, py);
            let e1 = edge(c.x, c.y, a.x, a.y, px, py);
            let e2 = edge(a.x, a.y, b.x, b.y, px, py);
            if !edge_is_owned(e0, b.x, b.y, c.x, c.y, area)
                || !edge_is_owned(e1, c.x, c.y, a.x, a.y, area)
                || !edge_is_owned(e2, a.x, a.y, b.x, b.y, area)
            {
                continue;
            }
            let w0 = e0 * inv_area;
            let w1 = e1 * inv_area;
            let w2 = e2 * inv_area;
            let inverse_w =
                w0 / interpolation_w[0] + w1 / interpolation_w[1] + w2 / interpolation_w[2];
            if inverse_w <= 0.0 {
                continue;
            }
            let reciprocal_w = 1.0 / inverse_w;
            let depth01 = interpolate_depth01(&a, &b, &c, px, py, inverse_area64);
            let index = y * framebuffer.width + x;
            if !fragment_visible(pass, framebuffer, index, depth01) {
                continue;
            }
            let pa = w0 / interpolation_w[0] * reciprocal_w;
            let pb = w1 / interpolation_w[1] * reciprocal_w;
            let pc = w2 / interpolation_w[2] * reciprocal_w;
            let view_position = a_view * pa + b_view * pb + c_view * pc;
            let base_color = [
                ca[0] * pa + cb[0] * pb + cc[0] * pc,
                ca[1] * pa + cb[1] * pb + cc[1] * pc,
                ca[2] * pa + cb[2] * pb + cc[2] * pc,
            ];
            let mut world_normal = angle_normalize(na_world * pa + nb_world * pb + nc_world * pc);
            let mut view_normal = angle_normalize(na_view * pa + nb_view * pb + nc_view * pc);
            if double_sided && !front_facing {
                world_normal = world_normal * -1.0;
                view_normal = view_normal * -1.0;
            }
            if flip_sided {
                world_normal = world_normal * -1.0;
                view_normal = view_normal * -1.0;
            }
            let depth = if camera.orthographic {
                -view_position.z
            } else {
                reciprocal_w
            };
            write_fragment(
                pass,
                framebuffer,
                index,
                depth,
                depth01,
                world_normal,
                view_normal,
                view_position,
                material,
                base_color,
                camera,
                renderer,
                style,
            );
        }
    }
}

fn rasterize_primitive(
    primitive: &NativePrimitive,
    camera: &CameraState,
    renderer: &RendererOptions,
    style: ResolvedStyle,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    let material = match primitive {
        NativePrimitive::Point { material, .. }
        | NativePrimitive::Line { material, .. }
        | NativePrimitive::Sphere { material, .. }
        | NativePrimitive::Cylinder { material, .. } => *material,
    };
    if !pass.accepts(material) {
        return;
    }
    match primitive {
        NativePrimitive::Point {
            center,
            size,
            material,
        } => rasterize_point(
            *center,
            *size,
            *material,
            camera,
            renderer,
            framebuffer,
            pass,
        ),
        NativePrimitive::Line {
            start,
            end,
            size,
            material,
        } => rasterize_line(
            *start,
            *end,
            *size,
            *material,
            camera,
            renderer,
            framebuffer,
            pass,
        ),
        NativePrimitive::Sphere {
            center,
            radius,
            material,
        } => rasterize_sphere_impostor(
            *center,
            *radius,
            *material,
            camera,
            renderer,
            style,
            framebuffer,
            pass,
        ),
        NativePrimitive::Cylinder {
            start,
            end,
            radius,
            top_cap,
            bottom_cap,
            material,
        } => rasterize_cylinder_impostor(
            *start,
            *end,
            *radius,
            *top_cap,
            *bottom_cap,
            *material,
            camera,
            renderer,
            style,
            framebuffer,
            pass,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn rasterize_point(
    center: Vec3,
    size: f32,
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    let Some(projected) = camera.project(center, framebuffer.width, framebuffer.height) else {
        return;
    };
    let point_size = (size * renderer.viewport.pixel_ratio as f32).max(1.0);
    let half = point_size * 0.5;
    let min_x = (projected.x - half).floor().max(0.0) as usize;
    let max_x = (projected.x + half)
        .ceil()
        .min(framebuffer.width as f32 - 1.0) as usize;
    let min_y = (projected.y - half).floor().max(0.0) as usize;
    let max_y = (projected.y + half)
        .ceil()
        .min(framebuffer.height as f32 - 1.0) as usize;
    if min_x > max_x || min_y > max_y {
        return;
    }
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            // LineRepresentation selects Points.pointStyle = circle. WebGL's
            // point coordinate spans the point-size square at pixel centers.
            let point_x = (x as f32 + 0.5 - (projected.x - half)) / point_size;
            let point_y = (y as f32 + 0.5 - (projected.y - half)) / point_size;
            let dx = point_x - 0.5;
            let dy = point_y - 0.5;
            if dx * dx + dy * dy > 0.25 {
                continue;
            }
            let index = y * framebuffer.width + x;
            if !fragment_visible(pass, framebuffer, index, projected.depth01) {
                continue;
            }
            write_unlit_fragment(
                pass,
                framebuffer,
                index,
                projected.depth,
                projected.depth01,
                material,
                camera,
                renderer,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn rasterize_line(
    start: Vec3,
    end: Vec3,
    size: f32,
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    let Some(clip_vertices) = molstar_line_clip_vertices(
        start,
        end,
        size,
        camera,
        renderer,
        framebuffer.width,
        framebuffer.height,
    ) else {
        return;
    };
    for triangle in [[0usize, 1usize, 2usize], [1usize, 3usize, 2usize]] {
        rasterize_clipped_line_triangle(
            [
                clip_vertices[triangle[0]],
                clip_vertices[triangle[1]],
                clip_vertices[triangle[2]],
            ],
            material,
            camera,
            renderer,
            framebuffer,
            pass,
        );
    }
}

/// Clip-space output of Mol*'s `lines.vert` shader. The six fixed-function
/// clip distances are kept in OpenGL convention; ANGLE's Metal Z remap turns
/// `-w <= z <= w` into Metal's equivalent `0 <= z' <= w` interval.
#[derive(Clone, Copy, Debug, Default)]
struct LineClipVertex {
    x: f32,
    y: f32,
    z: f32,
    w: f32,
}

#[allow(clippy::too_many_arguments)]
fn molstar_line_clip_vertices(
    start: Vec3,
    end: Vec3,
    size: f32,
    camera: &CameraState,
    renderer: &RendererOptions,
    width: usize,
    height: usize,
) -> Option<[LineClipVertex; 4]> {
    if width == 0 || height == 0 {
        return None;
    }
    let mut view_start = camera.view_position(start);
    let mut view_end = camera.view_position(end);
    if !view_start.is_finite() || !view_end.is_finite() {
        return None;
    }

    let projection = camera.projection_matrix();
    if !camera.orthographic {
        // Port `trimSegment` from the pinned Mol* lines vertex shader. It is
        // intentionally applied only when the camera plane (z = 0), rather
        // than the near plane, is crossed. Ordinary near-plane crossing is
        // left to homogeneous fixed-function clipping below.
        if view_start.z < 0.0 && view_end.z >= 0.0 {
            trim_molstar_line_segment(view_start, &mut view_end, projection);
        } else if view_end.z < 0.0 && view_start.z >= 0.0 {
            trim_molstar_line_segment(view_end, &mut view_start, projection);
        }
    }

    let clip_start = project_view_to_line_clip(view_start, projection, camera.orthographic);
    let clip_end = project_view_to_line_clip(view_end, projection, camera.orthographic);
    if !line_clip_vertex_is_finite(clip_start) || !line_clip_vertex_is_finite(clip_end) {
        return None;
    }

    // Preserve the shader's NDC/aspect/normalize staging. In particular, the
    // perpendicular is formed after multiplying dir.x by the viewport aspect,
    // then its X component is divided by that aspect again.
    let ndc_start_x = divide(clip_start.x, clip_start.w);
    let ndc_start_y = divide(clip_start.y, clip_start.w);
    let ndc_end_x = divide(clip_end.x, clip_end.w);
    let ndc_end_y = divide(clip_end.y, clip_end.w);
    let aspect = divide(width as f32, height as f32);
    let dir_x = (ndc_end_x - ndc_start_x) * aspect;
    let dir_y = ndc_end_y - ndc_start_y;
    let length_squared = fused_multiply_add(dir_y, dir_y, dir_x * dir_x);
    if !length_squared.is_finite() || length_squared <= f32::MIN_POSITIVE {
        return None;
    }
    let inverse_length = inverse_sqrt(length_squared);
    let normalized_x = dir_x * inverse_length;
    let normalized_y = dir_y * inverse_length;
    let linewidth = (size * renderer.viewport.pixel_ratio as f32).max(1.0);
    let offset_x = divide(normalized_y, aspect) * linewidth;
    let offset_y = -normalized_x * linewidth;
    let offset_x = divide(offset_x, height as f32);
    let offset_y = divide(offset_y, height as f32);

    Some([
        offset_line_clip_vertex(clip_start, -offset_x, -offset_y),
        offset_line_clip_vertex(clip_start, offset_x, offset_y),
        offset_line_clip_vertex(clip_end, -offset_x, -offset_y),
        offset_line_clip_vertex(clip_end, offset_x, offset_y),
    ])
}

fn trim_molstar_line_segment(start: Vec3, end: &mut Vec3, projection: [f32; 16]) {
    let near_estimate = divide(-0.5 * projection[14], projection[10]);
    let alpha = divide(near_estimate - start.z, end.z - start.z);
    end.x = fused_multiply_add(alpha, end.x - start.x, start.x);
    end.y = fused_multiply_add(alpha, end.y - start.y, start.y);
    end.z = fused_multiply_add(alpha, end.z - start.z, start.z);
}

fn project_view_to_line_clip(
    view: Vec3,
    projection: [f32; 16],
    orthographic: bool,
) -> LineClipVertex {
    LineClipVertex {
        x: projection[0] * view.x + projection[8] * view.z + projection[12],
        y: projection[5] * view.y + projection[9] * view.z + projection[13],
        z: projection[10] * view.z + projection[14],
        w: if orthographic { 1.0 } else { -view.z },
    }
}

fn offset_line_clip_vertex(
    mut vertex: LineClipVertex,
    offset_x: f32,
    offset_y: f32,
) -> LineClipVertex {
    vertex.x = fused_multiply_add(offset_x, vertex.w, vertex.x);
    vertex.y = fused_multiply_add(offset_y, vertex.w, vertex.y);
    vertex
}

fn line_clip_vertex_is_finite(vertex: LineClipVertex) -> bool {
    vertex.x.is_finite() && vertex.y.is_finite() && vertex.z.is_finite() && vertex.w.is_finite()
}

#[allow(clippy::too_many_arguments)]
fn rasterize_clipped_line_triangle(
    triangle: [LineClipVertex; 3],
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    if triangle
        .iter()
        .all(|&vertex| line_vertex_inside_all_clip_planes(vertex))
    {
        let (Some(a), Some(b), Some(c)) = (
            project_line_clip_vertex(triangle[0], camera, framebuffer.width, framebuffer.height),
            project_line_clip_vertex(triangle[1], camera, framebuffer.width, framebuffer.height),
            project_line_clip_vertex(triangle[2], camera, framebuffer.width, framebuffer.height),
        ) else {
            return;
        };
        rasterize_unlit_screen_triangle(a, b, c, material, camera, renderer, framebuffer, pass);
        return;
    }

    let polygon = clip_line_triangle_to_view_volume(triangle);
    if polygon.len() < 3 {
        return;
    }
    let Some(anchor) =
        project_line_clip_vertex(polygon[0], camera, framebuffer.width, framebuffer.height)
    else {
        return;
    };
    for index in 1..polygon.len() - 1 {
        let (Some(b), Some(c)) = (
            project_line_clip_vertex(
                polygon[index],
                camera,
                framebuffer.width,
                framebuffer.height,
            ),
            project_line_clip_vertex(
                polygon[index + 1],
                camera,
                framebuffer.width,
                framebuffer.height,
            ),
        ) else {
            continue;
        };
        rasterize_unlit_screen_triangle(
            anchor,
            b,
            c,
            material,
            camera,
            renderer,
            framebuffer,
            pass,
        );
    }
}

fn clip_line_triangle_to_view_volume(triangle: [LineClipVertex; 3]) -> Vec<LineClipVertex> {
    let mut input = triangle.to_vec();
    let mut output = Vec::with_capacity(9);
    for plane in 0..6 {
        if input.is_empty() {
            break;
        }
        output.clear();
        let mut previous = *input.last().expect("non-empty clipped line polygon");
        let mut previous_distance = line_clip_plane_distance(previous, plane);
        for &current in &input {
            let current_distance = line_clip_plane_distance(current, plane);
            let previous_inside = previous_distance >= 0.0;
            let current_inside = current_distance >= 0.0;
            if previous_inside != current_inside {
                let alpha = previous_distance / (previous_distance - current_distance);
                output.push(interpolate_line_clip_vertex(previous, current, alpha));
            }
            if current_inside {
                output.push(current);
            }
            previous = current;
            previous_distance = current_distance;
        }
        std::mem::swap(&mut input, &mut output);
    }
    input
}

fn line_vertex_inside_all_clip_planes(vertex: LineClipVertex) -> bool {
    (0..6).all(|plane| line_clip_plane_distance(vertex, plane) >= 0.0)
}

fn line_clip_plane_distance(vertex: LineClipVertex, plane: usize) -> f32 {
    match plane {
        0 => vertex.x + vertex.w,
        1 => vertex.w - vertex.x,
        2 => vertex.y + vertex.w,
        3 => vertex.w - vertex.y,
        4 => vertex.z + vertex.w,
        5 => vertex.w - vertex.z,
        _ => unreachable!("six homogeneous clip planes"),
    }
}

fn interpolate_line_clip_vertex(
    start: LineClipVertex,
    end: LineClipVertex,
    alpha: f32,
) -> LineClipVertex {
    LineClipVertex {
        x: fused_multiply_add(alpha, end.x - start.x, start.x),
        y: fused_multiply_add(alpha, end.y - start.y, start.y),
        z: fused_multiply_add(alpha, end.z - start.z, start.z),
        w: fused_multiply_add(alpha, end.w - start.w, start.w),
    }
}

fn project_line_clip_vertex(
    vertex: LineClipVertex,
    camera: &CameraState,
    width: usize,
    height: usize,
) -> Option<ProjectedVertex> {
    if !line_clip_vertex_is_finite(vertex) || vertex.w <= 0.0 {
        return None;
    }
    let ndc_x = f64::from(vertex.x) / f64::from(vertex.w);
    let ndc_y = f64::from(vertex.y) / f64::from(vertex.w);
    let x = ((ndc_x * 0.5 + 0.5) * width as f64) as f32;
    let y = ((0.5 - ndc_y * 0.5) * height as f64) as f32;
    let depth01 = ((vertex.z + vertex.w) * 0.5 / vertex.w).clamp(0.0, 1.0);
    let depth = if camera.orthographic {
        camera.depth_from_depth01(depth01)
    } else {
        vertex.w
    };
    Some(ProjectedVertex {
        x,
        y,
        depth,
        depth01,
    })
}

#[allow(clippy::too_many_arguments)]
fn rasterize_unlit_screen_triangle(
    mut a: ProjectedVertex,
    mut b: ProjectedVertex,
    mut c: ProjectedVertex,
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    for vertex in [&mut a, &mut b, &mut c] {
        vertex.x = snap_to_raster_subpixel(vertex.x);
        vertex.y = snap_to_raster_subpixel(vertex.y);
    }
    let area = edge(a.x, a.y, b.x, b.y, c.x, c.y);
    if !area.is_finite() || area.abs() < 0.000_001 {
        return;
    }
    let min_x = a.x.min(b.x).min(c.x).floor().max(0.0) as usize;
    let max_x =
        a.x.max(b.x)
            .max(c.x)
            .ceil()
            .min(framebuffer.width as f32 - 1.0) as usize;
    let min_y = a.y.min(b.y).min(c.y).floor().max(0.0) as usize;
    let max_y =
        a.y.max(b.y)
            .max(c.y)
            .ceil()
            .min(framebuffer.height as f32 - 1.0) as usize;
    let inverse_area64 = 1.0 / edge64(a.x, a.y, b.x, b.y, c.x, c.y);
    for y in min_y..=max_y {
        for x in min_x..=max_x {
            let px = x as f32 + 0.5;
            let py = y as f32 + 0.5;
            let e0 = edge(b.x, b.y, c.x, c.y, px, py);
            let e1 = edge(c.x, c.y, a.x, a.y, px, py);
            let e2 = edge(a.x, a.y, b.x, b.y, px, py);
            if !edge_is_owned(e0, b.x, b.y, c.x, c.y, area)
                || !edge_is_owned(e1, c.x, c.y, a.x, a.y, area)
                || !edge_is_owned(e2, a.x, a.y, b.x, b.y, area)
            {
                continue;
            }
            let depth01 = interpolate_depth01(&a, &b, &c, px, py, inverse_area64);
            let index = y * framebuffer.width + x;
            if !fragment_visible(pass, framebuffer, index, depth01) {
                continue;
            }
            let depth = if (a.depth - b.depth).abs() <= f32::EPSILON
                && (a.depth - c.depth).abs() <= f32::EPSILON
            {
                a.depth
            } else {
                let inverse_area = 1.0 / area;
                let w0 = e0 * inverse_area;
                let w1 = e1 * inverse_area;
                let w2 = e2 * inverse_area;
                1.0 / (w0 / a.depth + w1 / b.depth + w2 / c.depth)
            };
            write_unlit_fragment(
                pass,
                framebuffer,
                index,
                depth,
                depth01,
                material,
                camera,
                renderer,
            );
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn write_unlit_fragment(
    pass: &mut RasterPass<'_>,
    framebuffer: &mut Framebuffer,
    index: usize,
    depth: f32,
    depth01: f32,
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
) {
    let mut color = shade_unlit_material_linear(material, renderer);
    color[3] = pass.alpha(material) as f32;
    match &mut pass.target {
        RasterTarget::Opaque => {
            framebuffer.color[index] =
                fog_opaque_fragment(color, camera, renderer, depth01).map(quantize);
            framebuffer.depth[index] = depth;
            framebuffer.depth01[index] = depth01;
            framebuffer.normal[index] = Vec3::default();
        }
        RasterTarget::Wboit(target) => {
            target.accumulate(index, depth01, color, camera, renderer);
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn rasterize_sphere_impostor(
    center: Vec3,
    radius: f32,
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
    style: ResolvedStyle,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    const MAPPING: [[f32; 2]; 6] = [
        [-1.0, 1.0],
        [-1.0, -1.0],
        [1.0, 1.0],
        [-1.0, -1.0],
        [1.0, -1.0],
        [1.0, 1.0],
    ];
    const TRIANGLES: [[usize; 3]; 2] = [[0, 1, 2], [3, 4, 5]];

    let view_center = camera.view_position(center);
    let camera_depth = -view_center.z;
    if camera_depth <= camera.near || camera_depth >= camera.far || radius <= 0.0 {
        return;
    }
    let projection = camera.projection_matrix();
    let clip_z = fused_multiply_add(projection[10], view_center.z, projection[14]);
    let clip_w = if camera.orthographic {
        1.0
    } else {
        camera_depth
    };
    let center_depth01 = ((clip_z + clip_w) * 0.5 / clip_w).clamp(0.0, 1.0);
    let mut projected = [None; 6];
    let mut view_points = [Vec3::default(); 6];

    if camera.orthographic {
        for (index, mapping) in MAPPING.iter().enumerate() {
            let world =
                center + camera.right * (mapping[0] * radius) + camera.up * (mapping[1] * radius);
            let Some(vertex) = camera.project(world, framebuffer.width, framebuffer.height) else {
                continue;
            };
            view_points[index] = camera.screen_space_to_view_space(
                vertex.x / framebuffer.width as f32,
                1.0 - vertex.y / framebuffer.height as f32,
                center_depth01,
            );
            projected[index] = Some(vertex);
        }
    } else {
        let radius_squared = radius * radius;
        let pzr2 = fused_multiply_add(view_center.z, view_center.z, -radius_squared);
        let scaled_center = view_center * radius;
        let vx = fused_multiply_add(view_center.x, view_center.x, pzr2).sqrt();
        let vy = fused_multiply_add(view_center.y, view_center.y, pzr2).sqrt();
        let min_x = divide(
            fused_multiply_add(vx, view_center.x, -scaled_center.z),
            fused_multiply_add(vx, view_center.z, scaled_center.x),
        ) * projection[0];
        let max_x = divide(
            fused_multiply_add(vx, view_center.x, scaled_center.z),
            fused_multiply_add(vx, view_center.z, -scaled_center.x),
        ) * projection[0];
        let min_y = divide(
            fused_multiply_add(vy, view_center.y, -scaled_center.z),
            fused_multiply_add(vy, view_center.z, scaled_center.y),
        ) * projection[5];
        let max_y = divide(
            fused_multiply_add(vy, view_center.y, scaled_center.z),
            fused_multiply_add(vy, view_center.z, -scaled_center.y),
        ) * projection[5];
        for (index, mapping) in MAPPING.iter().enumerate() {
            let mut ndc_x = (max_x + min_x) * -0.5;
            let mut ndc_y = (max_y + min_y) * -0.5;
            ndc_x -= mapping[0] * (max_x - min_x) * 0.5;
            ndc_y -= mapping[1] * (max_y - min_y) * 0.5;
            let coords_x = ndc_x * 0.5 + 0.5;
            let coords_y = ndc_y * 0.5 + 0.5;
            view_points[index] =
                camera.clip_position_to_view_space(ndc_x * clip_w, ndc_y * clip_w, clip_z, clip_w);
            projected[index] = Some(ProjectedVertex {
                x: coords_x * framebuffer.width as f32,
                y: (1.0 - coords_y) * framebuffer.height as f32,
                depth: camera_depth,
                depth01: center_depth01,
            });
        }
    }

    for triangle in TRIANGLES {
        let (Some(mut a), Some(mut b), Some(mut c)) = (
            projected[triangle[0]],
            projected[triangle[1]],
            projected[triangle[2]],
        ) else {
            continue;
        };
        for vertex in [&mut a, &mut b, &mut c] {
            vertex.x = snap_to_raster_subpixel(vertex.x);
            vertex.y = snap_to_raster_subpixel(vertex.y);
        }
        let area = edge(a.x, a.y, b.x, b.y, c.x, c.y);
        if !area.is_finite() || area.abs() < 0.000_001 {
            continue;
        }
        let min_x = a.x.min(b.x).min(c.x).floor().max(0.0) as usize;
        let max_x =
            a.x.max(b.x)
                .max(c.x)
                .ceil()
                .min(framebuffer.width as f32 - 1.0) as usize;
        let min_y = a.y.min(b.y).min(c.y).floor().max(0.0) as usize;
        let max_y =
            a.y.max(b.y)
                .max(c.y)
                .ceil()
                .min(framebuffer.height as f32 - 1.0) as usize;
        let inverse_area = 1.0 / area;
        for y in min_y..=max_y {
            for x in min_x..=max_x {
                let px = x as f32 + 0.5;
                let py = y as f32 + 0.5;
                let e0 = edge(b.x, b.y, c.x, c.y, px, py);
                let e1 = edge(c.x, c.y, a.x, a.y, px, py);
                let e2 = edge(a.x, a.y, b.x, b.y, px, py);
                if !edge_is_owned(e0, b.x, b.y, c.x, c.y, area)
                    || !edge_is_owned(e1, c.x, c.y, a.x, a.y, area)
                    || !edge_is_owned(e2, a.x, a.y, b.x, b.y, area)
                {
                    continue;
                }
                let w0 = e0 * inverse_area;
                let w1 = e1 * inverse_area;
                let w2 = e2 * inverse_area;
                let view_point = interpolate_sphere_varying(&view_points, triangle, w0, w1, w2);
                let ray_origin = if camera.orthographic {
                    view_point
                } else {
                    Vec3::default()
                };
                let ray_direction = if camera.orthographic {
                    Vec3::new(0.0, 0.0, 1.0)
                } else {
                    angle_normalize(view_point)
                };
                let Some(hit) = intersect_sphere_shader(
                    ray_origin,
                    ray_direction,
                    view_center,
                    radius,
                    camera.orthographic,
                ) else {
                    continue;
                };
                let depth = -hit.position.z;
                if depth <= camera.near || depth >= camera.far {
                    continue;
                }
                let index = y * framebuffer.width + x;
                let depth01 = camera.impostor_depth01(hit.position.z);
                if !fragment_visible(pass, framebuffer, index, depth01) {
                    continue;
                }
                let world_normal = camera.right * hit.normal.x + camera.up * hit.normal.y
                    - camera.forward * hit.normal.z;
                write_fragment(
                    pass,
                    framebuffer,
                    index,
                    depth,
                    depth01,
                    world_normal,
                    hit.normal,
                    hit.position,
                    material,
                    color_f32(material.color),
                    camera,
                    renderer,
                    style,
                );
            }
        }
    }
}

fn interpolate_sphere_varying(
    view_points: &[Vec3; 6],
    triangle: [usize; 3],
    w0: f32,
    w1: f32,
    w2: f32,
) -> Vec3 {
    let mut point = view_points[triangle[0]] * w0
        + view_points[triangle[1]] * w1
        + view_points[triangle[2]] * w2;
    // `vPoint.z` is identical for every sphere proxy vertex. The fixed-function
    // interpolator preserves that constant bitwise; recombining it through
    // barycentric weights introduces a spurious ULP at analytic silhouettes.
    point.z = view_points[triangle[0]].z;
    point
}

#[allow(clippy::too_many_arguments)]
fn rasterize_cylinder_impostor(
    start: Vec3,
    end: Vec3,
    radius: f32,
    top_cap: bool,
    bottom_cap: bool,
    material: MeshMaterial,
    camera: &CameraState,
    renderer: &RendererOptions,
    style: ResolvedStyle,
    framebuffer: &mut Framebuffer,
    pass: &mut RasterPass<'_>,
) {
    const MAPPING: [[f32; 3]; 6] = [
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [1.0, -1.0, 1.0],
    ];
    const TRIANGLES: [[usize; 3]; 4] = [[0, 1, 2], [1, 4, 2], [2, 4, 3], [4, 5, 3]];

    let midpoint = (start + end) * 0.5;
    let camera_direction = if camera.orthographic {
        camera.forward * -1.0
    } else {
        angle_normalize(camera.position - midpoint)
    };
    let mut axis = end - start;
    if angle_dot(camera_direction, axis) < 0.0 {
        axis = axis * -1.0;
    }
    let raw_left = angle_cross(camera_direction, axis);
    let raw_up = angle_cross(raw_left, axis);
    let left = angle_normalize(raw_left) * radius;
    let up = angle_normalize(raw_up) * radius;
    let mut world_vertices = [Vec3::default(); 6];
    let mut projected_vertices = [None; 6];
    for (index, mapping) in MAPPING.iter().enumerate() {
        let position = midpoint + axis * mapping[0] + left * mapping[1] + up * mapping[2];
        world_vertices[index] = position;
        projected_vertices[index] = camera.project(position, framebuffer.width, framebuffer.height);
    }

    for triangle in TRIANGLES {
        let (Some(mut a), Some(mut b), Some(mut c)) = (
            projected_vertices[triangle[0]],
            projected_vertices[triangle[1]],
            projected_vertices[triangle[2]],
        ) else {
            continue;
        };
        for vertex in [&mut a, &mut b, &mut c] {
            vertex.x = snap_to_raster_subpixel(vertex.x);
            vertex.y = snap_to_raster_subpixel(vertex.y);
        }
        let area = edge(a.x, a.y, b.x, b.y, c.x, c.y);
        if !area.is_finite() || area.abs() < 0.000_001 {
            continue;
        }
        let min_x = a.x.min(b.x).min(c.x).floor().max(0.0) as usize;
        let max_x =
            a.x.max(b.x)
                .max(c.x)
                .ceil()
                .min(framebuffer.width as f32 - 1.0) as usize;
        let min_y = a.y.min(b.y).min(c.y).floor().max(0.0) as usize;
        let max_y =
            a.y.max(b.y)
                .max(c.y)
                .ceil()
                .min(framebuffer.height as f32 - 1.0) as usize;
        let inverse_area = 1.0 / area;
        for y in min_y..=max_y {
            for x in min_x..=max_x {
                let px = x as f32 + 0.5;
                let py = y as f32 + 0.5;
                let e0 = edge(b.x, b.y, c.x, c.y, px, py);
                let e1 = edge(c.x, c.y, a.x, a.y, px, py);
                let e2 = edge(a.x, a.y, b.x, b.y, px, py);
                if !edge_is_owned(e0, b.x, b.y, c.x, c.y, area)
                    || !edge_is_owned(e1, c.x, c.y, a.x, a.y, area)
                    || !edge_is_owned(e2, a.x, a.y, b.x, b.y, area)
                {
                    continue;
                }
                let w0 = e0 * inverse_area;
                let w1 = e1 * inverse_area;
                let w2 = e2 * inverse_area;
                let inverse_depth = w0 / a.depth + w1 / b.depth + w2 / c.depth;
                if inverse_depth <= 0.0 {
                    continue;
                }
                let interpolated_depth = 1.0 / inverse_depth;
                let pa = w0 / a.depth * interpolated_depth;
                let pb = w1 / b.depth * interpolated_depth;
                let pc = w2 / c.depth * interpolated_depth;
                let ray_origin = world_vertices[triangle[0]] * pa
                    + world_vertices[triangle[1]] * pb
                    + world_vertices[triangle[2]] * pc;
                let ray_direction = if camera.orthographic {
                    camera.forward
                } else {
                    angle_normalize(ray_origin - camera.position)
                };
                let Some(hit) = intersect_cylinder_shader(
                    ray_origin,
                    ray_direction,
                    start,
                    end,
                    radius,
                    top_cap,
                    bottom_cap,
                ) else {
                    continue;
                };
                let view_position = camera.view_position(hit.position);
                let depth = -view_position.z;
                if depth <= camera.near || depth >= camera.far {
                    continue;
                }
                let index = y * framebuffer.width + x;
                let depth01 = camera.impostor_depth01(view_position.z);
                if !fragment_visible(pass, framebuffer, index, depth01) {
                    continue;
                }
                write_fragment(
                    pass,
                    framebuffer,
                    index,
                    depth,
                    depth01,
                    hit.normal,
                    camera.normal_to_view(hit.normal),
                    view_position,
                    material,
                    color_f32(material.color),
                    camera,
                    renderer,
                    style,
                );
            }
        }
    }
}

fn fragment_visible(
    _pass: &RasterPass<'_>,
    framebuffer: &Framebuffer,
    index: usize,
    depth01: f32,
) -> bool {
    // The WBOIT pass shares the completed opaque depth buffer and deliberately
    // never updates it. Consequently every transparent fragment in front of
    // opaque geometry contributes, independent of transparent draw order.
    depth01 < framebuffer.depth01[index]
}

#[allow(clippy::too_many_arguments)]
fn write_fragment(
    pass: &mut RasterPass<'_>,
    framebuffer: &mut Framebuffer,
    index: usize,
    depth: f32,
    depth01: f32,
    world_normal: Vec3,
    view_normal: Vec3,
    view_position: Vec3,
    material: MeshMaterial,
    base_color: [f32; 3],
    camera: &CameraState,
    renderer: &RendererOptions,
    style: ResolvedStyle,
) {
    let mut color = shade_material_base_color_linear(
        material,
        base_color,
        view_normal,
        view_position,
        renderer,
        style,
    );
    color[3] = pass.alpha(material) as f32;
    match &mut pass.target {
        RasterTarget::Opaque => {
            framebuffer.color[index] =
                fog_opaque_fragment(color, camera, renderer, depth01).map(quantize);
            framebuffer.depth[index] = depth;
            framebuffer.depth01[index] = depth01;
            framebuffer.normal[index] = world_normal;
        }
        RasterTarget::Wboit(target) => target.accumulate(index, depth01, color, camera, renderer),
    }
}

#[derive(Clone, Copy)]
struct RayHit {
    position: Vec3,
    normal: Vec3,
}

fn intersect_sphere_shader(
    ray_origin: Vec3,
    ray_direction: Vec3,
    center: Vec3,
    radius: f32,
    orthographic: bool,
) -> Option<RayHit> {
    let sphere_direction = if orthographic {
        ray_origin - center
    } else {
        center
    };
    let b = angle_dot(ray_direction, sphere_direction);
    let determinant = fused_multiply_add(
        b,
        b,
        radius * radius - angle_dot(sphere_direction, sphere_direction),
    );
    if determinant < 0.0 || radius <= 0.0 {
        return None;
    }
    let root = determinant.sqrt();
    let t = if orthographic { b + root } else { b - root };
    let position = Vec3::new(
        fused_multiply_add(ray_direction.x, t, ray_origin.x),
        fused_multiply_add(ray_direction.y, t, ray_origin.y),
        fused_multiply_add(ray_direction.z, t, ray_origin.z),
    );
    Some(RayHit {
        position,
        normal: angle_normalize(position - center),
    })
}

fn intersect_cylinder_shader(
    origin: Vec3,
    direction: Vec3,
    start: Vec3,
    end: Vec3,
    radius: f32,
    top_cap: bool,
    bottom_cap: bool,
) -> Option<RayHit> {
    let axis = end - start;
    let offset = origin - start;
    let axis_squared = angle_dot(axis, axis);
    let axis_ray = angle_dot(axis, direction);
    let axis_offset = angle_dot(axis, offset);
    let k2 = fused_multiply_add(-axis_ray, axis_ray, axis_squared);
    let k1 = fused_multiply_add(
        -axis_offset,
        axis_ray,
        axis_squared * angle_dot(offset, direction),
    );
    let k0 = fused_multiply_add(
        -(radius * radius),
        axis_squared,
        fused_multiply_add(
            -axis_offset,
            axis_offset,
            axis_squared * angle_dot(offset, offset),
        ),
    );
    let determinant = fused_multiply_add(-k2, k0, k1 * k1);
    if determinant < 0.0 || k2.abs() <= f32::MIN_POSITIVE || radius <= 0.0 {
        return None;
    }
    let root = determinant.sqrt();
    let mut t = divide(-k1 - root, k2);
    let y = fused_multiply_add(t, axis_ray, axis_offset);
    if y > 0.0 && y < axis_squared {
        let axis_scale = divide(y, axis_squared);
        let position = origin + direction * t;
        let normal = Vec3::new(
            fused_multiply_add(
                -axis.x,
                axis_scale,
                fused_multiply_add(direction.x, t, offset.x),
            ),
            fused_multiply_add(
                -axis.y,
                axis_scale,
                fused_multiply_add(direction.y, t, offset.y),
            ),
            fused_multiply_add(
                -axis.z,
                axis_scale,
                fused_multiply_add(direction.z, t, offset.z),
            ),
        ) * divide(1.0, radius);
        return Some(RayHit { position, normal });
    }

    if top_cap && y < 0.0 && axis_ray.abs() > f32::MIN_POSITIVE {
        t = divide(-axis_offset, axis_ray);
        if fused_multiply_add(k2, t, k1).abs() < root {
            return Some(RayHit {
                position: origin + direction * t,
                normal: (axis * -1.0).normalized(),
            });
        }
    } else if bottom_cap && y >= 0.0 && axis_ray.abs() > f32::MIN_POSITIVE {
        t = divide(axis_squared - axis_offset, axis_ray);
        if fused_multiply_add(k2, t, k1).abs() < root {
            return Some(RayHit {
                position: origin + direction * t,
                normal: axis.normalized(),
            });
        }
    }
    None
}

fn angle_dot(a: Vec3, b: Vec3) -> f32 {
    fused_multiply_add(b.z, a.z, fused_multiply_add(b.y, a.y, a.x * b.x))
}

fn angle_cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        fused_multiply_add(a.y, b.z, -(a.z * b.y)),
        fused_multiply_add(a.z, b.x, -(a.x * b.z)),
        fused_multiply_add(a.x, b.y, -(a.y * b.x)),
    )
}

fn angle_normalize(value: Vec3) -> Vec3 {
    let length_squared = angle_dot(value, value);
    if length_squared <= 0.000_000_000_001 {
        Vec3::default()
    } else {
        value * inverse_sqrt(length_squared)
    }
}

fn mesh_normal_to_view(camera: &CameraState, normal: Vec3) -> Vec3 {
    let view = camera.view_matrix();
    angle_normalize(Vec3::new(
        fused_multiply_add(
            view[8],
            normal.z,
            fused_multiply_add(view[4], normal.y, view[0] * normal.x),
        ),
        fused_multiply_add(
            view[9],
            normal.z,
            fused_multiply_add(view[5], normal.y, view[1] * normal.x),
        ),
        fused_multiply_add(
            view[10],
            normal.z,
            fused_multiply_add(view[6], normal.y, view[2] * normal.x),
        ),
    ))
}

fn edge(ax: f32, ay: f32, bx: f32, by: f32, px: f32, py: f32) -> f32 {
    (px - ax) * (by - ay) - (py - ay) * (bx - ax)
}

fn edge_is_owned(value: f32, ax: f32, ay: f32, bx: f32, by: f32, area: f32) -> bool {
    if area > 0.0 {
        value > 0.0 || (value == 0.0 && ((by > ay) || (by == ay && bx < ax)))
    } else {
        value < 0.0 || (value == 0.0 && ((by < ay) || (by == ay && bx > ax)))
    }
}

fn edge64(ax: f32, ay: f32, bx: f32, by: f32, px: f32, py: f32) -> f64 {
    (px as f64 - ax as f64) * (by as f64 - ay as f64)
        - (py as f64 - ay as f64) * (bx as f64 - ax as f64)
}

fn interpolate_depth01(
    a: &ProjectedVertex,
    b: &ProjectedVertex,
    c: &ProjectedVertex,
    px: f32,
    py: f32,
    inverse_area: f64,
) -> f32 {
    let w0 = edge64(b.x, b.y, c.x, c.y, px, py) * inverse_area;
    let w1 = edge64(c.x, c.y, a.x, a.y, px, py) * inverse_area;
    let w2 = 1.0 - w0 - w1;
    truncate_positive_f64_to_f32(
        w0 * a.depth01 as f64 + w1 * b.depth01 as f64 + w2 * c.depth01 as f64,
    )
}

/// Apple depth interpolation truncates positive depth32Float plane results
/// toward zero at the final f32 boundary. Keep the plane calculation in f64
/// (matching the measured interpolation coefficients) and reproduce only the
/// target conversion here.
fn truncate_positive_f64_to_f32(value: f64) -> f32 {
    let rounded = value as f32;
    if rounded > 0.0 && rounded as f64 > value {
        f32::from_bits(rounded.to_bits() - 1)
    } else {
        rounded
    }
}

fn snap_to_raster_subpixel(value: f32) -> f32 {
    (value * 256.0).round() * (1.0 / 256.0)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::excessive_precision)]

    use super::*;
    use crate::model::{BoundingSphere, Face};
    use crate::render::camera::resolve_camera;

    #[test]
    fn object_alpha_controls_partition_before_float32_uniform_conversion() {
        let material = MeshMaterial::opaque(0x1b9e77);
        let mut target = WboitTarget::new(1, 1);
        for alpha in [-0.1f64, 0.0, 0.371, 0.7, 1.0 - f64::EPSILON, 1.0, 1.4] {
            let clamped = alpha.clamp(0.0, 1.0);
            let mut opaque = RasterPass::opaque();
            let mut transparent = RasterPass::wboit(&mut target);
            opaque.object_alpha = Some(clamped);
            transparent.object_alpha = Some(clamped);
            assert_eq!(opaque.accepts(material), clamped == 1.0);
            assert_eq!(
                transparent.accepts(material),
                clamped > 0.0 && clamped < 1.0
            );
            assert_eq!(transparent.alpha(material), clamped);
        }
    }

    #[test]
    fn scene_alpha_factor_reaches_mesh_and_all_analytic_geometry_without_material_edits() {
        use crate::mesh::build_native_render_scene_with_summaries;
        use crate::options::MeshOptions;
        use crate::parser::parse_molecule_with_options;
        use crate::render::style::resolve_style;

        let (mut renderer, camera) = line_test_camera("perspective");
        renderer.camera.fog.name = "off".into();
        let mut kinds = std::collections::BTreeSet::new();
        for representation in [
            "default",
            "spacefill",
            "lines",
            "points",
            "molecular-surface",
        ] {
            let provider = if matches!(representation, "lines" | "points") {
                "default"
            } else {
                representation
            };
            let options = MeshOptions::from_json(format!(
                r#"{{"format":"xyz","representation":"{provider}","quality":"high","style":"illustrative"}}"#
            ).as_bytes()).unwrap();
            let molecule = parse_molecule_with_options(
                b"2\npaired atoms\nC -0.65 0 0\nC 0.65 0 0\n",
                &options,
            )
            .unwrap();
            let mut scene = build_native_render_scene_with_summaries(&molecule, &options);
            // Lines/points are selected by visual providers, not exposed as
            // public representation names. Exercise their native ownership
            // and shader paths directly with opaque source materials.
            if matches!(representation, "lines" | "points") {
                for primitive in &mut scene.primitives {
                    *primitive = if representation == "lines" {
                        NativePrimitive::Line {
                            start: Vec3::new(-0.65, 0.0, 0.0),
                            end: Vec3::new(0.65, 0.0, 0.0),
                            size: 2.0,
                            material: MeshMaterial::opaque(0x1b9e77),
                        }
                    } else {
                        NativePrimitive::Point {
                            center: Vec3::default(),
                            size: 3.0,
                            material: MeshMaterial::opaque(0x1b9e77),
                        }
                    };
                }
                for object in &mut scene.render_objects {
                    object.geometry_kind = representation;
                    if let NativeRenderGeometry::Primitives { kind, .. } = &mut object.geometry {
                        *kind = representation;
                    }
                }
            }
            let style = resolve_style(&options, &renderer);
            for object in &scene.render_objects {
                kinds.insert(object.geometry_kind);
            }
            for factor in [0.0, 0.371, 0.7, 1.0, 1.4] {
                for object in &mut scene.render_objects {
                    object.alpha_factor = factor;
                }
                let mut opaque = Framebuffer::new(128, 128, [255; 4]);
                rasterize_scene(
                    &scene,
                    &camera,
                    &renderer,
                    style,
                    &mut opaque,
                    &mut RasterPass::opaque(),
                );
                assert_eq!(
                    opaque.depth01.iter().any(|&depth| depth < 1.0),
                    factor >= 1.0,
                    "{representation}, {factor}"
                );
                let mut target = WboitTarget::new(128, 128);
                rasterize_scene(
                    &scene,
                    &camera,
                    &renderer,
                    style,
                    &mut opaque,
                    &mut RasterPass::wboit(&mut target),
                );
                let packed = target.packed_nearest_depth_alpha_rgba8();
                let expected = super::super::postprocessing::pack_depth_alpha_rgba8(
                    0.5,
                    factor.clamp(0.0, 1.0) as f32,
                )[3];
                let mut covered = 0;
                for pixel in packed {
                    if pixel != [255; 4] {
                        covered += 1;
                        assert_eq!(pixel[3], expected, "{representation}, {factor}");
                    }
                }
                assert_eq!(
                    covered > 0,
                    factor > 0.0 && factor < 1.0,
                    "{representation}, {factor}"
                );
            }
        }
        assert_eq!(
            kinds.into_iter().collect::<Vec<_>>(),
            ["cylinders", "lines", "mesh", "points", "spheres"]
        );
    }

    fn line_test_camera(mode: &str) -> (RendererOptions, CameraState) {
        let options = RendererOptions::from_json(
            format!(
                r#"{{"viewport":{{"width":128,"height":128}},"camera":{{"mode":"{mode}","view":{{"name":"snapshot","params":{{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":20}}}}}}}}"#
            )
            .as_bytes(),
        )
        .unwrap();
        let sphere = BoundingSphere {
            center: Vec3::default(),
            radius: 2.0,
            ..BoundingSphere::default()
        };
        let camera = resolve_camera(&options, &sphere, 128, 128).unwrap();
        (options, camera)
    }

    fn clipped_line_polygons(vertices: [LineClipVertex; 4]) -> Vec<Vec<LineClipVertex>> {
        [[0usize, 1usize, 2usize], [1usize, 3usize, 2usize]]
            .into_iter()
            .map(|triangle| {
                clip_line_triangle_to_view_volume([
                    vertices[triangle[0]],
                    vertices[triangle[1]],
                    vertices[triangle[2]],
                ])
            })
            .collect()
    }

    #[test]
    fn perspective_line_crossing_camera_plane_uses_molstar_trim_then_fixed_clipping() {
        let (renderer, camera) = line_test_camera("perspective");
        let vertices = molstar_line_clip_vertices(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 11.0),
            2.0,
            &camera,
            &renderer,
            128,
            128,
        )
        .unwrap();
        let projection = camera.projection_matrix();
        let near_estimate = divide(-0.5 * projection[14], projection[10]);
        assert_eq!(vertices[2].w.to_bits(), (-near_estimate).to_bits());
        assert!(vertices[2].w > 0.0 && vertices[2].w < camera.near);

        let polygons = clipped_line_polygons(vertices);
        assert!(polygons.iter().any(|polygon| polygon.len() >= 3));
        assert!(polygons
            .iter()
            .flatten()
            .copied()
            .all(line_vertex_inside_all_clip_planes));

        let mut framebuffer = Framebuffer::new(128, 128, [0, 0, 0, 0]);
        rasterize_line(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 11.0),
            2.0,
            MeshMaterial::opaque(0xff0000),
            &camera,
            &renderer,
            &mut framebuffer,
            &mut RasterPass::opaque(),
        );
        assert!(framebuffer.depth01.iter().any(|&depth| depth < 1.0));
    }

    #[test]
    fn perspective_line_fully_behind_camera_is_rejected_by_homogeneous_clipping() {
        let (renderer, camera) = line_test_camera("perspective");
        let vertices = molstar_line_clip_vertices(
            Vec3::new(-1.0, 0.0, 11.0),
            Vec3::new(1.0, 1.0, 12.0),
            2.0,
            &camera,
            &renderer,
            128,
            128,
        )
        .unwrap();
        assert!(clipped_line_polygons(vertices).iter().all(Vec::is_empty));
    }

    #[test]
    fn orthographic_line_skips_perspective_camera_plane_trim() {
        let (renderer, camera) = line_test_camera("orthographic");
        let end = Vec3::new(1.0, 1.0, 11.0);
        let vertices = molstar_line_clip_vertices(
            Vec3::new(-1.0, 0.0, 0.0),
            end,
            2.0,
            &camera,
            &renderer,
            128,
            128,
        )
        .unwrap();
        let direct =
            project_view_to_line_clip(camera.view_position(end), camera.projection_matrix(), true);
        assert_eq!(vertices[2].z.to_bits(), direct.z.to_bits());
        assert_eq!(vertices[2].w, 1.0);
    }

    #[test]
    fn ordinary_near_plane_crossing_is_left_to_fixed_function_clipping() {
        let (renderer, camera) = line_test_camera("perspective");
        let end = Vec3::new(1.0, 1.0, 3.0);
        let vertices = molstar_line_clip_vertices(
            Vec3::new(-1.0, 0.0, 0.0),
            end,
            2.0,
            &camera,
            &renderer,
            128,
            128,
        )
        .unwrap();
        assert_eq!(
            vertices[2].w.to_bits(),
            (-camera.view_position(end).z).to_bits()
        );
        assert!(vertices[2].w < camera.near);
        let polygons = clipped_line_polygons(vertices);
        assert!(polygons.iter().any(|polygon| polygon.len() >= 3));
        assert!(polygons
            .iter()
            .flatten()
            .any(|&vertex| { line_clip_plane_distance(vertex, 4).abs() <= 0.000_01 }));
    }

    #[test]
    fn transparent_mesh_discards_interior_even_when_double_sided() {
        let (renderer, camera) = line_test_camera("perspective");
        for (reverse, alpha_tenths) in [(false, 7), (true, 7), (true, 10)] {
            let mesh = Mesh {
                vertices: vec![
                    Vec3::new(-1.0, -1.0, 0.0),
                    Vec3::new(1.0, -1.0, 0.0),
                    Vec3::new(0.0, 1.0, 0.0),
                ],
                faces: vec![if reverse {
                    Face { a: 0, b: 2, c: 1 }
                } else {
                    Face { a: 0, b: 1, c: 2 }
                }],
                face_materials: vec![MeshMaterial::with_alpha_tenths(0x1b9e77, alpha_tenths)],
                ..Mesh::default()
            };
            let mut framebuffer = Framebuffer::new(128, 128, [0; 4]);
            let mut target = WboitTarget::new(128, 128);
            let mut pass = if alpha_tenths == 10 {
                RasterPass::opaque()
            } else {
                RasterPass::wboit(&mut target)
            };
            rasterize_mesh_faces(
                &mesh,
                0..1,
                true,
                false,
                &camera,
                &renderer,
                ResolvedStyle {
                    ignore_light: true,
                    occlusion: false,
                    outline: false,
                },
                &mut framebuffer,
                &mut pass,
            );
            if alpha_tenths == 10 {
                assert!(framebuffer.depth01.iter().any(|&depth| depth < 1.0));
            } else {
                let layer = target.evaluated_framebuffer();
                assert_eq!(layer.color.iter().any(|color| color[3] > 0), !reverse);
                assert_eq!(
                    target
                        .packed_nearest_depth_alpha_rgba8()
                        .iter()
                        .any(|&depth| depth != [255; 4]),
                    !reverse
                );
            }
        }
    }

    #[test]
    fn mesh_near_plane_crossing_preserves_clipped_varyings_and_renders() {
        let (renderer, camera) = line_test_camera("perspective");
        let world = [
            Vec3::new(-1.5, -1.0, 0.0),
            Vec3::new(1.5, -1.0, 0.0),
            Vec3::new(0.0, 1.5, 3.0),
        ];
        let normals = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        assert!(camera.project(world[2], 128, 128).is_none());

        let world_normals = normals.map(angle_normalize);
        let view_normals = world_normals.map(|normal| mesh_normal_to_view(&camera, normal));
        let view_positions = world.map(|position| camera.view_position(position));
        let polygon = clip_mesh_triangle_to_view_volume(
            world,
            view_positions,
            world_normals,
            view_normals,
            [[0.0; 3]; 3],
            &camera,
        );
        assert_eq!(polygon.len(), 4);
        assert!(polygon
            .iter()
            .all(|vertex| line_vertex_inside_all_clip_planes(vertex.clip)));
        assert!(polygon
            .iter()
            .any(|vertex| line_clip_plane_distance(vertex.clip, 4).abs() <= 0.000_01));
        assert!(polygon.iter().any(|vertex| {
            vertex.world_normal.x > 0.0
                && vertex.world_normal.z > 0.0
                && vertex.world_normal.y == 0.0
        }));

        let mesh = Mesh {
            vertices: world.to_vec(),
            normals: normals.to_vec(),
            faces: vec![Face { a: 0, b: 1, c: 2 }],
            face_materials: vec![MeshMaterial::opaque(0x1b9e77)],
            ..Mesh::default()
        };
        let mut framebuffer = Framebuffer::new(128, 128, [0, 0, 0, 0]);
        rasterize_mesh_faces(
            &mesh,
            0..1,
            true,
            false,
            &camera,
            &renderer,
            ResolvedStyle {
                ignore_light: true,
                occlusion: false,
                outline: false,
            },
            &mut framebuffer,
            &mut RasterPass::opaque(),
        );
        assert!(framebuffer.depth01.iter().any(|&depth| depth < 1.0));
        assert!(framebuffer.color.iter().any(|&color| color[3] == 255));
    }

    #[test]
    fn mesh_fully_outside_far_plane_is_rejected_by_homogeneous_clipping() {
        let (_, camera) = line_test_camera("perspective");
        let world = [
            Vec3::new(-1.0, -1.0, -20.0),
            Vec3::new(1.0, -1.0, -20.0),
            Vec3::new(0.0, 1.0, -20.0),
        ];
        let world_normals = [Vec3::new(0.0, 0.0, 1.0); 3];
        let view_normals = world_normals.map(|normal| mesh_normal_to_view(&camera, normal));
        let view_positions = world.map(|position| camera.view_position(position));
        let polygon = clip_mesh_triangle_to_view_volume(
            world,
            view_positions,
            world_normals,
            view_normals,
            [[0.0; 3]; 3],
            &camera,
        );
        assert!(polygon.is_empty());
    }

    #[test]
    fn mesh_far_and_screen_plane_crossings_are_clipped_without_rejection() {
        let (_, camera) = line_test_camera("perspective");
        let world_normals = [Vec3::new(0.0, 0.0, 1.0); 3];
        let view_normals = world_normals.map(|normal| mesh_normal_to_view(&camera, normal));
        for (world, plane) in [
            (
                [
                    Vec3::new(-1.0, -1.0, 0.0),
                    Vec3::new(1.0, -1.0, 0.0),
                    Vec3::new(0.0, 1.0, -25.0),
                ],
                5,
            ),
            (
                [
                    Vec3::new(-1.0, -1.0, 0.0),
                    Vec3::new(20.0, -1.0, 0.0),
                    Vec3::new(-1.0, 1.0, 0.0),
                ],
                1,
            ),
        ] {
            let view_positions = world.map(|position| camera.view_position(position));
            let polygon = clip_mesh_triangle_to_view_volume(
                world,
                view_positions,
                world_normals,
                view_normals,
                [[0.0; 3]; 3],
                &camera,
            );
            assert!(polygon.len() >= 3);
            assert!(polygon
                .iter()
                .all(|vertex| line_vertex_inside_all_clip_planes(vertex.clip)));
            assert!(polygon
                .iter()
                .any(|vertex| { line_clip_plane_distance(vertex.clip, plane).abs() <= 0.000_01 }));
        }
    }

    #[test]
    fn mesh_vertex_shader_normalization_makes_normal_magnitude_irrelevant() {
        let (renderer, camera) = line_test_camera("perspective");
        let vertices = vec![
            Vec3::new(-1.5, -1.0, 0.0),
            Vec3::new(1.5, -1.0, 0.0),
            Vec3::new(0.0, 1.5, 0.0),
        ];
        let render = |normals: Vec<Vec3>| {
            let mesh = Mesh {
                vertices: vertices.clone(),
                normals,
                faces: vec![Face { a: 0, b: 1, c: 2 }],
                face_materials: vec![MeshMaterial::opaque(0x1b9e77)],
                ..Mesh::default()
            };
            let mut framebuffer = Framebuffer::new(128, 128, [0, 0, 0, 0]);
            rasterize_mesh_faces(
                &mesh,
                0..1,
                true,
                false,
                &camera,
                &renderer,
                ResolvedStyle {
                    ignore_light: false,
                    occlusion: false,
                    outline: false,
                },
                &mut framebuffer,
                &mut RasterPass::opaque(),
            );
            framebuffer
        };
        let unit = render(vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]);
        let scaled = render(vec![
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 3.0, 0.0),
            Vec3::new(0.0, 0.0, 4.0),
        ]);
        assert_eq!(unit.color, scaled.color);
        assert_eq!(unit.depth01, scaled.depth01);
        assert_eq!(unit.normal, scaled.normal);
        assert!(unit.color.iter().any(|&color| color[3] == 255));
    }

    #[test]
    fn mesh_vertex_colors_are_interpolated_before_fragment_shading() {
        let (renderer, camera) = line_test_camera("perspective");
        let mesh = Mesh {
            vertices: vec![
                Vec3::new(-1.5, -1.0, 0.0),
                Vec3::new(1.5, -1.0, 0.0),
                Vec3::new(0.0, 1.5, 0.0),
            ],
            normals: vec![Vec3::new(0.0, 0.0, 1.0); 3],
            faces: vec![Face { a: 0, b: 1, c: 2 }],
            vertex_colors: vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            face_materials: vec![MeshMaterial::opaque(0xff0000)],
            ..Mesh::default()
        };
        let mut framebuffer = Framebuffer::new(128, 128, [0, 0, 0, 0]);
        rasterize_mesh_faces(
            &mesh,
            0..1,
            true,
            false,
            &camera,
            &renderer,
            ResolvedStyle {
                ignore_light: true,
                occlusion: false,
                outline: false,
            },
            &mut framebuffer,
            &mut RasterPass::opaque(),
        );

        assert!(framebuffer.color.iter().any(|color| {
            color[3] == 255 && color[..3].iter().filter(|&&channel| channel > 0).count() >= 2
        }));
    }

    #[test]
    fn mesh_view_position_keeps_vertex_varying_rounding() {
        let (_, camera) = line_test_camera("perspective");
        let world = [
            Vec3::new(-1.234_567, 2.345_678, -0.456_789),
            Vec3::new(3.210_987, -1.876_543, 1.234_567),
            Vec3::new(0.765_432, 1.111_111, -2.222_222),
        ];
        let view = world.map(|position| camera.view_position(position));
        let mut difference = None;
        for numerator_a in 1..16 {
            for numerator_b in 1..16 - numerator_a {
                let pa = numerator_a as f32 / 17.0;
                let pb = numerator_b as f32 / 17.0;
                let pc = 1.0 - pa - pb;
                let interpolated = view[0] * pa + view[1] * pb + view[2] * pc;
                let interpolated_world = world[0] * pa + world[1] * pb + world[2] * pc;
                let recomputed = camera.view_position(interpolated_world);
                if [interpolated.x, interpolated.y, interpolated.z].map(f32::to_bits)
                    != [recomputed.x, recomputed.y, recomputed.z].map(f32::to_bits)
                {
                    difference = Some((interpolated, recomputed));
                    break;
                }
            }
            if difference.is_some() {
                break;
            }
        }
        let (interpolated, recomputed) =
            difference.expect("vertex varying and fragment recomputation must round differently");
        assert_ne!(
            [interpolated.x, interpolated.y, interpolated.z].map(f32::to_bits),
            [recomputed.x, recomputed.y, recomputed.z].map(f32::to_bits)
        );
    }

    #[test]
    fn positive_depth_conversion_truncates_instead_of_rounding_up() {
        let lower = f32::from_bits(0x3f00_0000);
        let upper = f32::from_bits(lower.to_bits() + 1);
        let value = lower as f64 + (upper as f64 - lower as f64) * 0.75;
        assert_eq!(value as f32, upper);
        assert_eq!(truncate_positive_f64_to_f32(value), lower);
        assert_eq!(truncate_positive_f64_to_f32(lower as f64), lower);
    }

    #[test]
    fn apple_raster_subpixel_grid_matches_reference_metal_vertices() {
        assert_eq!(snap_to_raster_subpixel(590.5085), 590.5078125);
        assert_eq!(snap_to_raster_subpixel(229.28821), 229.2890625);
    }

    #[test]
    fn sphere_varying_interpolation_preserves_constant_vertex_z() {
        let z = f32::from_bits(0xc1ab_2454);
        let points = [
            Vec3::new(1.0, 2.0, z),
            Vec3::new(3.0, 4.0, z),
            Vec3::new(5.0, 6.0, z),
            Vec3::default(),
            Vec3::default(),
            Vec3::default(),
        ];
        let weights = [0.793_707_4, 0.169_582_2, 0.036_710_426];
        let naive_z = (points[0] * weights[0] + points[1] * weights[1] + points[2] * weights[2]).z;
        assert_eq!(naive_z.to_bits(), 0xc1ab_2455);
        assert_eq!(
            interpolate_sphere_varying(&points, [0, 1, 2], weights[0], weights[1], weights[2],)
                .z
                .to_bits(),
            0xc1ab_2454,
        );
    }

    #[test]
    fn snapped_depth_plane_matches_reference_metal_fragment() {
        let a = ProjectedVertex {
            x: 590.5078125,
            y: 229.48046875,
            depth: 132.82796,
            depth01: f32::from_bits(1055606440),
        };
        let b = ProjectedVertex {
            x: 587.1640625,
            y: 229.2890625,
            depth: 132.99644,
            depth01: f32::from_bits(1055667670),
        };
        let c = ProjectedVertex {
            x: 590.6328125,
            y: 229.875,
            depth: 132.87575,
            depth01: f32::from_bits(1055623821),
        };
        let inverse_area = 1.0 / edge64(a.x, a.y, b.x, b.y, c.x, c.y);
        assert_eq!(
            interpolate_depth01(&a, &b, &c, 588.5, 229.5, inverse_area).to_bits(),
            1055650034
        );
    }

    #[test]
    fn top_left_edge_ownership_matches_reference_metal_coverage() {
        let a = (1.5, 1.5);
        let b = (5.5, 1.5);
        let c = (1.5, 5.5);
        let area = edge(a.0, a.1, b.0, b.1, c.0, c.1);
        let covered = |x: f32, y: f32| {
            edge_is_owned(edge(b.0, b.1, c.0, c.1, x, y), b.0, b.1, c.0, c.1, area)
                && edge_is_owned(edge(c.0, c.1, a.0, a.1, x, y), c.0, c.1, a.0, a.1, area)
                && edge_is_owned(edge(a.0, a.1, b.0, b.1, x, y), a.0, a.1, b.0, b.1, area)
        };
        assert!(covered(1.5, 1.5));
        assert!(covered(4.5, 1.5));
        assert!(covered(1.5, 4.5));
        assert!(!covered(5.5, 1.5));
        assert!(!covered(4.5, 2.5));
        assert!(!covered(2.5, 4.5));
    }

    #[test]
    fn top_down_viewport_maps_webgl_ccw_front_face_to_positive_edge_area() {
        // NDC (-1,-1), (1,-1), (-1,1) is counter-clockwise. After WebGL's
        // bottom-up coordinates are stored as top-down rows it becomes the
        // bottom-left, bottom-right, top-left triangle below.
        let bottom_left = (1.5, 5.5);
        let bottom_right = (5.5, 5.5);
        let top_left = (1.5, 1.5);
        assert!(
            edge(
                bottom_left.0,
                bottom_left.1,
                bottom_right.0,
                bottom_right.1,
                top_left.0,
                top_left.1,
            ) > 0.0
        );
    }
}
