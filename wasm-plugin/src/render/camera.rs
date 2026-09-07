use crate::model::{BoundingSphere, Vec3};

use super::options::RendererOptions;
use super::postprocessing::angle_metal_math::{
    divide, fused_multiply_add, inverse_sqrt, reciprocal,
};

#[derive(Clone, Copy)]
pub(super) struct CameraState {
    pub(super) position: Vec3,
    pub(super) target: Vec3,
    pub(super) up: Vec3,
    pub(super) right: Vec3,
    pub(super) forward: Vec3,
    pub(super) radius: f32,
    pub(super) radius_max: f32,
    pub(super) distance: f32,
    pub(super) near: f32,
    pub(super) far: f32,
    pub(super) fov: f32,
    pub(super) aspect: f32,
    pub(super) orthographic: bool,
    fov64: f64,
    aspect64: f64,
    near64: f64,
    far64: f64,
    normalized_far64: f64,
    distance64: f64,
    viewport_width64: f64,
    viewport_height64: f64,
    view_matrix64: [f64; 16],
    view_offset_x: f32,
    view_offset_y: f32,
}

impl CameraState {
    /// Mol*'s ordinary Canvas3D camera uses scale 1. XR may override this
    /// value, but XR is outside the static renderer contract. Keep the value
    /// explicit because SSAO radius and blur depth bias are multiplied by it
    /// before their float32 uniform upload.
    pub(super) const fn scale(&self) -> f64 {
        1.0
    }

    /// Column-major world-to-view matrix using Mol*'s right-handed camera
    /// basis, in the same flat-array layout passed to a WebGL `mat4` uniform.
    pub(super) fn view_matrix(&self) -> [f32; 16] {
        self.view_matrix64.map(|value| value as f32)
    }

    /// Column-major OpenGL clip-space projection matrix, constructed with the
    /// same double-precision camera arithmetic Mol* uses before uniform upload.
    pub(super) fn projection_matrix(&self) -> [f32; 16] {
        self.projection_matrix64().map(|value| value as f32)
    }

    pub(super) fn inverse_projection_matrix(&self) -> [f32; 16] {
        // Mol* stores matrices in double-backed JavaScript arrays, runs its
        // generic `Mat4.invert`, and converts to float32 only while uploading
        // the uniform. Using the rounded projection uniform as the inverse
        // input loses bits, especially in orthographic and jittered views.
        invert_mat4_molstar64(self.projection_matrix64()).map(|value| value as f32)
    }

    /// Mol* multiplies its JavaScript Number projection and view matrices on
    /// the CPU, then uploads the resulting global uniform as float32.
    pub(super) fn projection_view_matrix(&self) -> [f32; 16] {
        multiply_mat4_column_major64(self.projection_matrix64(), self.view_matrix64)
            .map(|value| value as f32)
    }

    fn projection_matrix64(&self) -> [f64; 16] {
        let offset_x = self.view_offset_x as f64;
        let offset_y = self.view_offset_y as f64;
        if self.orthographic {
            // Preserve Camera.update/updateOrtho staging rather than reducing
            // it algebraically. Although the expressions are equivalent, the
            // intermediate JavaScript Number rounding is observable after the
            // float32 uniform upload.
            let height = 2.0 * js_tan64(self.fov64 / 2.0) * self.distance64;
            let zoom = self.viewport_height64 / height;
            let full_left = -self.viewport_width64 / 2.0;
            let full_right = self.viewport_width64 / 2.0;
            let full_top = self.viewport_height64 / 2.0;
            let full_bottom = -self.viewport_height64 / 2.0;
            let dx = (full_right - full_left) / (2.0 * zoom);
            let dy = (full_top - full_bottom) / (2.0 * zoom);
            let cx = (full_right + full_left) / 2.0;
            let cy = (full_top + full_bottom) / 2.0;
            let mut left = cx - dx;
            let mut right = cx + dx;
            let mut top = cy + dy;
            let mut bottom = cy - dy;
            if offset_x != 0.0 || offset_y != 0.0 {
                // Multisample jitter uses a full-size view offset, so both
                // Mol* width ratios are exactly one.
                let zoom_w = zoom;
                let zoom_h = zoom;
                let scale_w = (full_right - full_left) / self.viewport_width64;
                let scale_h = (full_top - full_bottom) / self.viewport_height64;
                left += scale_w * (offset_x / zoom_w);
                right = left + scale_w * (self.viewport_width64 / zoom_w);
                top -= scale_h * (offset_y / zoom_h);
                bottom = top - scale_h * (self.viewport_height64 / zoom_h);
            }
            let w = 1.0 / (right - left);
            let h = 1.0 / (top - bottom);
            let p = 1.0 / (self.far64 - self.near64);
            let x = (right + left) * w;
            let y = (top + bottom) * h;
            let z = (self.far64 + self.near64) * p;
            [
                2.0 * w,
                0.0,
                0.0,
                0.0,
                0.0,
                2.0 * h,
                0.0,
                0.0,
                0.0,
                0.0,
                -2.0 * p,
                0.0,
                -x,
                -y,
                -z,
                1.0,
            ]
        } else {
            let top = self.near64 * js_tan64(self.fov64 * 0.5);
            let height = 2.0 * top;
            let width = self.aspect64 * height;
            let left = -0.5 * width + offset_x * width / self.aspect_pixel_width();
            let top = top - offset_y * height / self.aspect_pixel_height();
            let right = left + width;
            let bottom = top - height;
            [
                2.0 * self.near64 / (right - left),
                0.0,
                0.0,
                0.0,
                0.0,
                2.0 * self.near64 / (top - bottom),
                0.0,
                0.0,
                (right + left) / (right - left),
                (top + bottom) / (top - bottom),
                -(self.far64 + self.near64) / (self.far64 - self.near64),
                -1.0,
                0.0,
                0.0,
                -2.0 * self.far64 * self.near64 / (self.far64 - self.near64),
                0.0,
            ]
        }
    }

    fn aspect_pixel_width(&self) -> f64 {
        self.viewport_width64
    }

    fn aspect_pixel_height(&self) -> f64 {
        self.viewport_height64
    }

    pub(super) fn project(
        &self,
        point: Vec3,
        width: usize,
        height: usize,
    ) -> Option<ProjectedVertex> {
        let view = self.view_position(point);
        let camera_z = -view.z;
        if !camera_z.is_finite() || camera_z <= self.near || camera_z >= self.far {
            return None;
        }
        let projection = self.projection_matrix();
        // ANGLE's translated Metal vertex path preserves the separate
        // multiply/add staging of the projection-matrix product.
        let clip_x = projection[0] * view.x + projection[8] * view.z + projection[12];
        let clip_y = projection[5] * view.y + projection[9] * view.z + projection[13];
        let clip_w = if self.orthographic { 1.0 } else { camera_z };
        let clip_z = projection[10] * view.z + projection[14];
        let metal_clip_z = (clip_z + clip_w) * 0.5;
        Some(ProjectedVertex {
            // Metal's fixed-function viewport transform does not first round
            // clip / w to a shader-visible float32 NDC value. Preserve the
            // ratio through the viewport mapping before returning to the
            // float32 screen-coordinate domain used by the rasterizer.
            x: clip_to_viewport_coordinate(clip_x, clip_w, width, false),
            y: clip_to_viewport_coordinate(clip_y, clip_w, height, true),
            depth: camera_z,
            depth01: (metal_clip_z / clip_w).clamp(0.0, 1.0),
        })
    }

    pub(super) fn normal_to_view(&self, normal: Vec3) -> Vec3 {
        let view = self.view_matrix();
        Vec3::new(
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
        )
        .normalized()
    }

    pub(super) fn view_position(&self, point: Vec3) -> Vec3 {
        let view = self.view_matrix();
        Vec3::new(
            fused_multiply_add(
                view[8],
                point.z,
                fused_multiply_add(
                    view[4],
                    point.y,
                    fused_multiply_add(view[0], point.x, view[12]),
                ),
            ),
            fused_multiply_add(
                view[9],
                point.z,
                fused_multiply_add(
                    view[5],
                    point.y,
                    fused_multiply_add(view[1], point.x, view[13]),
                ),
            ),
            fused_multiply_add(
                view[10],
                point.z,
                fused_multiply_add(
                    view[6],
                    point.y,
                    fused_multiply_add(view[2], point.x, view[14]),
                ),
            ),
        )
    }
}

#[inline]
fn clip_to_viewport_coordinate(clip: f32, clip_w: f32, extent: usize, inverted: bool) -> f32 {
    let ndc = f64::from(clip) / f64::from(clip_w);
    let normalized = if inverted {
        0.5 - ndc * 0.5
    } else {
        ndc * 0.5 + 0.5
    };
    (normalized * extent as f64) as f32
}

#[derive(Clone, Copy)]
pub(super) struct ProjectedVertex {
    pub(super) x: f32,
    pub(super) y: f32,
    pub(super) depth: f32,
    pub(super) depth01: f32,
}

pub(super) fn resolve_camera(
    renderer: &RendererOptions,
    sphere: &BoundingSphere,
    width: u32,
    height: u32,
) -> Result<CameraState, String> {
    let fov64 = renderer.camera.fov.to_radians();
    let fov = fov64 as f32;
    let aspect64 = width as f64 / height as f64;
    let aspect = aspect64 as f32;
    let target;
    let position;
    let up;
    let radius;
    let radius_max;
    let distance64;
    let radius64;
    let radius_max64;
    let view_position64;
    let view_target64;
    let mut view_up64;
    if renderer.camera.view.name == "snapshot" {
        let params = &renderer.camera.view.params;
        let snapshot_position = params.position.expect("validated snapshot camera position");
        let snapshot_target = params.target.expect("validated snapshot camera target");
        target = array_vec3(snapshot_target);
        position = array_vec3(snapshot_position);
        up = array_vec3(params.up.expect("validated snapshot camera up")).normalized();
        radius64 = params.radius.expect("validated snapshot camera radius");
        radius_max64 = params
            .radius_max
            .expect("validated snapshot camera radius-max");
        radius = radius64 as f32;
        radius_max = radius_max64 as f32;
        distance64 = vector_distance64(snapshot_position, snapshot_target);
        view_position64 = snapshot_position;
        view_target64 = snapshot_target;
        view_up64 = params.up.expect("validated snapshot camera up");
    } else {
        view_target64 = sphere.center64();
        target = array_vec3(view_target64);
        radius64 = sphere.radius64().max(0.01);
        radius = radius64 as f32;
        radius_max = radius;
        radius_max64 = radius64;
        let aspect_factor64 = if height < width { 1.0 } else { aspect64 };
        let distance = if renderer.camera.mode == "orthographic" {
            (radius64 / aspect_factor64) / js_tan64(fov64 * 0.5)
        } else {
            (radius64 / aspect_factor64) / js_sin64(fov64 * 0.5)
        }
        .abs();
        distance64 = distance;
        if renderer.camera.view.name == "orbit" {
            let params = &renderer.camera.view.params;
            let azimuth = params.azimuth.unwrap_or(0.0).to_radians();
            let elevation = params.elevation.unwrap_or(0.0).to_radians();
            let direction = normalize64([
                js_sin64(azimuth) * js_cos64(elevation),
                js_sin64(elevation),
                js_cos64(azimuth) * js_cos64(elevation),
            ]);
            view_position64 = add_scaled64(view_target64, direction, distance64);
            let forward = normalize64(subtract64(view_target64, view_position64));
            let right = normalize64(cross64(forward, [0.0, 1.0, 0.0]));
            view_up64 = normalize64(cross64(right, forward));
            let roll = params.roll.unwrap_or(0.0);
            if roll != 0.0 {
                let roll = roll.to_radians();
                view_up64 = normalize64(add_scaled64(
                    scale64(view_up64, js_cos64(roll)),
                    right,
                    js_sin64(roll),
                ));
            }
            position = array_vec3(view_position64);
            up = array_vec3(view_up64);
        } else {
            view_position64 = [
                view_target64[0],
                view_target64[1],
                view_target64[2] + distance64,
            ];
            view_up64 = [0.0, 1.0, 0.0];
            position = array_vec3(view_position64);
            up = Vec3::new(0.0, 1.0, 0.0);
        }
    }
    let forward = (target - position).normalized();
    if forward.squared_length() <= 0.000_001 {
        return Err("renderer camera position and target must differ".into());
    }
    let right = forward.cross(up).normalized();
    if right.squared_length() <= 0.000_001 {
        return Err("renderer camera up must not be parallel to its viewing direction".into());
    }
    let distance = distance64 as f32;
    let min_near64 = renderer.camera.clipping.min_near;
    let normalized_far64 = if renderer.camera.clipping.far {
        radius64
    } else {
        radius_max64
    }
    .max(renderer.camera.clipping.min_far);
    let near64 = radius_max64
        .min(min_near64)
        .max(distance64 - radius64)
        .max(0.01);
    let mut far64 = (distance64 + normalized_far64).max(min_near64);
    if near64 == far64 {
        far64 = near64 + 0.01;
    }
    let near = near64 as f32;
    let far = far64 as f32;
    let view_matrix64 = molstar_look_at64(view_position64, view_target64, view_up64);
    let view_matrix = view_matrix64.map(|value| value as f32);
    // These vectors are the rows of the actual float32 view uniform. Keeping
    // the analytic proxy basis tied to that boundary avoids a second,
    // independently rounded normalization path for non-axis-aligned views.
    let right = Vec3::new(view_matrix[0], view_matrix[4], view_matrix[8]);
    let up = Vec3::new(view_matrix[1], view_matrix[5], view_matrix[9]);
    let forward = Vec3::new(-view_matrix[2], -view_matrix[6], -view_matrix[10]);
    Ok(CameraState {
        position,
        target,
        up,
        right,
        forward,
        radius,
        radius_max,
        distance,
        near,
        far,
        fov,
        aspect,
        orthographic: renderer.camera.mode == "orthographic",
        fov64,
        aspect64,
        near64,
        far64,
        normalized_far64,
        distance64,
        viewport_width64: width as f64,
        viewport_height64: height as f64,
        view_matrix64,
        view_offset_x: 0.0,
        view_offset_y: 0.0,
    })
}

impl CameraState {
    pub(super) fn with_view_offset(self, x: f32, y: f32) -> Self {
        Self {
            view_offset_x: x,
            view_offset_y: y,
            ..self
        }
    }

    /// Resolve Mol*'s camera fog distances in JavaScript Number precision,
    /// then apply the float32 WebGL-uniform boundary once.
    pub(in crate::render) fn fog_range(&self, intensity: f64) -> (f32, f32) {
        let fog_near = self.distance64 + self.normalized_far64 * ((50.0 - intensity) / 50.0);
        (fog_near as f32, self.far64 as f32)
    }

    #[allow(dead_code)]
    pub(super) fn screen_ray(&self, x: f32, y: f32, width: usize, height: usize) -> (Vec3, Vec3) {
        let ndc_x = (x + self.view_offset_x) / width as f32 * 2.0 - 1.0;
        let ndc_y = 1.0 - (y + self.view_offset_y) / height as f32 * 2.0;
        let half_height = self.distance * (self.fov * 0.5).tan();
        if self.orthographic {
            let origin = self.position
                + self.right * (ndc_x * half_height * self.aspect)
                + self.up * (ndc_y * half_height);
            (origin, self.forward)
        } else {
            let direction = self.forward
                + self.right * (ndc_x * (self.fov * 0.5).tan() * self.aspect)
                + self.up * (ndc_y * (self.fov * 0.5).tan());
            let length_squared = fused_multiply_add(
                direction.z,
                direction.z,
                fused_multiply_add(direction.y, direction.y, direction.x * direction.x),
            );
            let direction = direction * inverse_sqrt(length_squared);
            (self.position, direction)
        }
    }

    #[allow(dead_code)]
    pub(super) fn projected_radius_pixels(&self, point: Vec3, radius: f32, height: usize) -> f32 {
        let depth = (point - self.position).dot(self.forward).max(self.near);
        let reference_depth = if self.orthographic {
            self.distance
        } else {
            depth
        };
        radius / (reference_depth * (self.fov * 0.5).tan()) * height as f32 * 0.5
    }

    pub(super) fn view_position_at_pixel(
        &self,
        x: usize,
        y: usize,
        depth01: f32,
        width: usize,
        height: usize,
    ) -> Vec3 {
        self.view_position_at_fragment(
            x as f32 + 0.5,
            (height - y) as f32 - 0.5,
            depth01,
            width,
            height,
        )
    }

    pub(super) fn view_position_at_fragment(
        &self,
        fragment_x: f32,
        fragment_y: f32,
        depth01: f32,
        width: usize,
        height: usize,
    ) -> Vec3 {
        let inverse_width = reciprocal(width as f32);
        let inverse_height = reciprocal(height as f32);
        let ndc_x = fused_multiply_add(fragment_x, inverse_width * 2.0, -1.0);
        let ndc_y = fused_multiply_add(fragment_y, inverse_height * 2.0, -1.0);
        self.ndc_to_view_space(ndc_x, ndc_y, depth01)
    }

    pub(super) fn view_position_at_neighbor_fragment(
        &self,
        fragment_x: f32,
        fragment_y: f32,
        depth01: f32,
        width: usize,
        height: usize,
    ) -> Vec3 {
        let inverse_width = reciprocal(width as f32);
        let inverse_height = reciprocal(height as f32);
        let ndc_x = fused_multiply_add(fragment_x, inverse_width * 2.0, -1.0);
        let ndc_y = fused_multiply_add(fragment_y, inverse_height * 2.0, -1.0);
        self.ndc_to_view_space_with_fused_w(ndc_x, ndc_y, depth01)
    }

    pub(super) fn project_view_position(
        &self,
        position: Vec3,
        width: usize,
        height: usize,
        bounds: [f32; 4],
    ) -> Option<(usize, usize)> {
        let [offset_x, offset_y] = self.project_view_position_offset(position, bounds)?;
        let gl_x = (offset_x * width as f32).floor() as usize;
        let gl_y = (offset_y * height as f32).floor() as usize;
        let x = gl_x.min(width.saturating_sub(1));
        let y = height - 1 - gl_y.min(height.saturating_sub(1));
        Some((x, y))
    }

    pub(super) fn project_view_position_offset(
        &self,
        position: Vec3,
        bounds: [f32; 4],
    ) -> Option<[f32; 2]> {
        if !position.x.is_finite() || !position.y.is_finite() || !position.z.is_finite() {
            return None;
        }
        let projection = self.projection_matrix();
        let clip_x = fused_multiply_add(projection[8], position.z, projection[0] * position.x)
            + projection[12];
        let clip_y = fused_multiply_add(projection[9], position.z, projection[5] * position.y)
            + projection[13];
        let clip_w = if self.orthographic { 1.0 } else { -position.z };
        let offset_x = fused_multiply_add(divide(clip_x, clip_w), 0.5, 0.5);
        let offset_y = fused_multiply_add(divide(clip_y, clip_w), 0.5, 0.5);
        if !offset_x.is_finite()
            || !offset_y.is_finite()
            || offset_x < bounds[0]
            || offset_x > bounds[2]
            || offset_y < bounds[1]
            || offset_y > bounds[3]
        {
            return None;
        }
        Some([offset_x, offset_y])
    }

    pub(super) fn screen_space_to_view_space(
        &self,
        coords_x: f32,
        coords_y: f32,
        depth01: f32,
    ) -> Vec3 {
        let x = fused_multiply_add(coords_x, 2.0, -1.0);
        let y = fused_multiply_add(coords_y, 2.0, -1.0);
        self.ndc_to_view_space(x, y, depth01)
    }

    /// Apply Mol*'s `uInvProjection * gl_Position` vertex-shader path without
    /// round-tripping through normalized depth. Sphere impostors interpolate
    /// this varying directly; reconstructing it from depth01 loses one or more
    /// float32 ULPs at analytic silhouette boundaries.
    pub(super) fn clip_position_to_view_space(&self, x: f32, y: f32, z: f32, w: f32) -> Vec3 {
        let inverse = self.inverse_projection_matrix();
        let px = fused_multiply_add(
            inverse[12],
            w,
            fused_multiply_add(
                inverse[8],
                z,
                fused_multiply_add(inverse[4], y, inverse[0] * x),
            ),
        );
        let py = fused_multiply_add(
            inverse[13],
            w,
            fused_multiply_add(
                inverse[9],
                z,
                fused_multiply_add(inverse[5], y, inverse[1] * x),
            ),
        );
        let pz = fused_multiply_add(
            inverse[14],
            w,
            fused_multiply_add(
                inverse[10],
                z,
                fused_multiply_add(inverse[6], y, inverse[2] * x),
            ),
        );
        let pw = fused_multiply_add(
            inverse[15],
            w,
            fused_multiply_add(
                inverse[11],
                z,
                fused_multiply_add(inverse[7], y, inverse[3] * x),
            ),
        );
        Vec3::new(divide(px, pw), divide(py, pw), divide(pz, pw))
    }

    pub(super) fn screen_space_to_view_space_with_fused_w(
        &self,
        coords_x: f32,
        coords_y: f32,
        depth01: f32,
    ) -> Vec3 {
        let x = fused_multiply_add(coords_x, 2.0, -1.0);
        let y = fused_multiply_add(coords_y, 2.0, -1.0);
        self.ndc_to_view_space_with_fused_w(x, y, depth01)
    }

    fn ndc_to_view_space(&self, x: f32, y: f32, depth01: f32) -> Vec3 {
        let inverse = self.inverse_projection_matrix();
        let z = depth01 * 2.0 - 1.0;
        let px = fused_multiply_add(inverse[8], z, inverse[0] * x) + inverse[12];
        let py = fused_multiply_add(inverse[9], z, inverse[5] * y) + inverse[13];
        let pz = fused_multiply_add(inverse[10], z, inverse[14]);
        let pw = inverse[11] * z + inverse[15];
        Vec3::new(divide(px, pw), divide(py, pw), divide(pz, pw))
    }

    fn ndc_to_view_space_with_fused_w(&self, x: f32, y: f32, depth01: f32) -> Vec3 {
        let inverse = self.inverse_projection_matrix();
        let z = depth01 * 2.0 - 1.0;
        let px = fused_multiply_add(inverse[8], z, inverse[0] * x) + inverse[12];
        let py = fused_multiply_add(inverse[9], z, inverse[5] * y) + inverse[13];
        let pz = fused_multiply_add(inverse[10], z, inverse[14]);
        let pw = fused_multiply_add(inverse[11], z, inverse[15]);
        Vec3::new(divide(px, pw), divide(py, pw), divide(pz, pw))
    }

    pub(super) fn view_z_from_depth01(&self, depth01: f32) -> f32 {
        let inverse = self.inverse_projection_matrix();
        let z = depth01 * 2.0 - 1.0;
        let pz = fused_multiply_add(inverse[10], z, inverse[14]);
        let pw = inverse[11] * z + inverse[15];
        divide(pz, pw)
    }

    /// Reproduce the shared Mol* sphere/cylinder fragment-shader `calcDepth`
    /// helper from camera-space Z and the uploaded float32 projection uniform.
    pub(super) fn impostor_depth01(&self, view_z: f32) -> f32 {
        if !view_z.is_finite() {
            return 1.0;
        }
        let projection = self.projection_matrix();
        let clip_z = fused_multiply_add(view_z, projection[10], projection[14]);
        let clip_w = fused_multiply_add(view_z, projection[11], projection[15]);
        fused_multiply_add(divide(clip_z, clip_w), 0.5, 0.5).clamp(0.0, 1.0)
    }

    pub(super) fn depth_from_depth01(&self, depth: f32) -> f32 {
        if depth >= 1.0 {
            return self.far;
        }
        if self.orthographic {
            -fused_multiply_add(depth, self.near - self.far, -self.near)
        } else {
            let denominator = fused_multiply_add(self.far - self.near, depth, -self.far);
            -divide(self.near * self.far, denominator)
        }
    }
}

fn array_vec3(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0] as f32, value[1] as f32, value[2] as f32)
}

#[inline]
fn js_sin64(value: f64) -> f64 {
    libm::sin(value)
}

#[inline]
fn js_cos64(value: f64) -> f64 {
    libm::cos(value)
}

#[inline]
fn js_tan64(value: f64) -> f64 {
    libm::tan(value)
}

fn molstar_look_at64(eye: [f64; 3], center: [f64; 3], up: [f64; 3]) -> [f64; 16] {
    let mut z = normalize64([eye[0] - center[0], eye[1] - center[1], eye[2] - center[2]]);
    let mut x = normalize64(cross64(up, z));
    let mut y = normalize64(cross64(z, x));
    // Preserve the zero-vector behavior of Mol*'s gl-matrix-derived lookAt.
    if !x.iter().all(|value| value.is_finite()) {
        x = [0.0; 3];
    }
    if !y.iter().all(|value| value.is_finite()) {
        y = [0.0; 3];
    }
    if !z.iter().all(|value| value.is_finite()) {
        z = [0.0; 3];
    }
    [
        x[0],
        y[0],
        z[0],
        0.0,
        x[1],
        y[1],
        z[1],
        0.0,
        x[2],
        y[2],
        z[2],
        0.0,
        -dot64(x, eye),
        -dot64(y, eye),
        -dot64(z, eye),
        1.0,
    ]
}

fn normalize64(value: [f64; 3]) -> [f64; 3] {
    let length = dot64(value, value).sqrt();
    if length == 0.0 {
        [0.0; 3]
    } else {
        [value[0] / length, value[1] / length, value[2] / length]
    }
}

fn cross64(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn subtract64(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale64(value: [f64; 3], scale: f64) -> [f64; 3] {
    [value[0] * scale, value[1] * scale, value[2] * scale]
}

fn add_scaled64(base: [f64; 3], direction: [f64; 3], scale: f64) -> [f64; 3] {
    [
        base[0] + direction[0] * scale,
        base[1] + direction[1] * scale,
        base[2] + direction[2] * scale,
    ]
}

fn dot64(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Mol* `Mat4.tryInvert`, kept in the same statement order so JavaScript
/// Number intermediates round at the same places before uniform upload.
fn invert_mat4_molstar64(a: [f64; 16]) -> [f64; 16] {
    let a00 = a[0];
    let a01 = a[1];
    let a02 = a[2];
    let a03 = a[3];
    let a10 = a[4];
    let a11 = a[5];
    let a12 = a[6];
    let a13 = a[7];
    let a20 = a[8];
    let a21 = a[9];
    let a22 = a[10];
    let a23 = a[11];
    let a30 = a[12];
    let a31 = a[13];
    let a32 = a[14];
    let a33 = a[15];

    let b00 = a00 * a11 - a01 * a10;
    let b01 = a00 * a12 - a02 * a10;
    let b02 = a00 * a13 - a03 * a10;
    let b03 = a01 * a12 - a02 * a11;
    let b04 = a01 * a13 - a03 * a11;
    let b05 = a02 * a13 - a03 * a12;
    let b06 = a20 * a31 - a21 * a30;
    let b07 = a20 * a32 - a22 * a30;
    let b08 = a20 * a33 - a23 * a30;
    let b09 = a21 * a32 - a22 * a31;
    let b10 = a21 * a33 - a23 * a31;
    let b11 = a22 * a33 - a23 * a32;

    let det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
    debug_assert!(det != 0.0, "camera projection matrix must be invertible");
    let det = 1.0 / det;

    [
        (a11 * b11 - a12 * b10 + a13 * b09) * det,
        (a02 * b10 - a01 * b11 - a03 * b09) * det,
        (a31 * b05 - a32 * b04 + a33 * b03) * det,
        (a22 * b04 - a21 * b05 - a23 * b03) * det,
        (a12 * b08 - a10 * b11 - a13 * b07) * det,
        (a00 * b11 - a02 * b08 + a03 * b07) * det,
        (a32 * b02 - a30 * b05 - a33 * b01) * det,
        (a20 * b05 - a22 * b02 + a23 * b01) * det,
        (a10 * b10 - a11 * b08 + a13 * b06) * det,
        (a01 * b08 - a00 * b10 - a03 * b06) * det,
        (a30 * b04 - a31 * b02 + a33 * b00) * det,
        (a21 * b02 - a20 * b04 - a23 * b00) * det,
        (a11 * b07 - a10 * b09 - a12 * b06) * det,
        (a00 * b09 - a01 * b07 + a02 * b06) * det,
        (a31 * b01 - a30 * b03 - a32 * b00) * det,
        (a20 * b03 - a21 * b01 + a22 * b00) * det,
    ]
}

fn multiply_mat4_column_major64(a: [f64; 16], b: [f64; 16]) -> [f64; 16] {
    let mut out = [0.0; 16];
    for column in 0..4 {
        for row in 0..4 {
            out[column * 4 + row] = (0..4)
                .map(|index| a[index * 4 + row] * b[column * 4 + index])
                .sum();
        }
    }
    out
}

fn vector_distance64(a: [f64; 3], b: [f64; 3]) -> f64 {
    let x = a[0] - b[0];
    let y = a[1] - b[1];
    let z = a[2] - b[2];
    (x * x + y * y + z * z).sqrt()
}

#[cfg(test)]
mod tests {
    use super::{clip_to_viewport_coordinate, js_cos64, js_sin64, js_tan64, resolve_camera};
    use crate::model::{BoundingSphere, Vec3};
    use crate::render::options::RendererOptions;
    use crate::render::postprocessing::angle_metal_math::{divide, fused_multiply_add};

    #[test]
    fn viewport_transform_keeps_fixed_function_clip_ratio_precision() {
        let clip_y = f32::from_bits(0xc17a_e766);
        let clip_w = f32::from_bits(0x430c_a276);
        let screen_y = clip_to_viewport_coordinate(clip_y, clip_w, 937, true);
        assert_eq!(screen_y.to_bits(), 0x4402_2f5f);
        assert_eq!((screen_y * 256.0).round() / 256.0, 520.738_3);
    }

    #[test]
    fn portable_number_trigonometry_matches_pinned_v8_bits() {
        let angle = std::f64::consts::PI / 8.0;
        assert_eq!(js_tan64(angle).to_bits(), 0x3fda_8279_99fc_ef32);
        assert_eq!(js_sin64(angle).to_bits(), 0x3fd8_7de2_a6ae_a963);
        assert_eq!(js_cos64(angle).to_bits(), 0x3fed_906b_cf32_8d46);
        for (degrees, sin, cos, tan) in [
            (
                -61.75_f64,
                0xbfec_3041_c5fe_202f,
                0x3fde_4ade_92c7_3350,
                0xbffd_c706_d63d_b1af,
            ),
            (
                17.125,
                0x3fd2_d863_99fd_e8cd,
                0x3fee_94cd_fa88_930c,
                0x3fd3_b833_fad2_c808,
            ),
            (
                35.0,
                0x3fe2_5abc_f87c_4978,
                0x3fea_367e_5915_8747,
                0x3fe6_6819_a3a0_bf7a,
            ),
            (
                89.5,
                0x3fef_ffb0_2599_c9cd,
                0x3f81_df37_c495_4c0b,
                0x405c_a5ac_7197_8b14,
            ),
        ] {
            let angle = degrees.to_radians();
            assert_eq!(js_sin64(angle).to_bits(), sin);
            assert_eq!(js_cos64(angle).to_bits(), cos);
            assert_eq!(js_tan64(angle).to_bits(), tan);
        }
    }

    #[test]
    fn sphere_inverse_projection_matches_reference_metal_vpoint_z() {
        let clip_z = f32::from_bits(0x4085_fa44);
        let clip_w = f32::from_bits(0x41ab_2455);
        let inverse_z = f32::from_bits(0xbc1b_da40);
        let inverse_w = f32::from_bits(0x3d47_17a5);
        let projected_w = fused_multiply_add(inverse_w, clip_w, inverse_z * clip_z);
        assert_eq!(projected_w.to_bits(), 0x3f80_0001);
        assert_eq!(divide(-clip_w, projected_w).to_bits(), 0xc1ab_2454);
    }

    #[test]
    fn clipping_matches_molstar_min_far_and_clip_far_equations() {
        let sphere = BoundingSphere {
            center: Vec3::new(0.0, 0.0, 0.0),
            radius: 2.0,
            ..BoundingSphere::default()
        };
        let clipped = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":20}},"clipping":{"far":true,"min-near":0.3,"min-far":5}}}"#,
        )
        .unwrap();
        let full = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":20}},"clipping":{"far":false,"min-near":0.3,"min-far":5}}}"#,
        )
        .unwrap();

        let clipped = resolve_camera(&clipped, &sphere, 800, 600).unwrap();
        let full = resolve_camera(&full, &sphere, 800, 600).unwrap();
        assert_eq!(clipped.near, 8.0);
        assert_eq!(clipped.far, 15.0);
        assert_eq!(full.near, 8.0);
        assert_eq!(full.far, 30.0);

        let (fog_near, fog_far) = clipped.fog_range(15.3);
        assert_eq!(fog_near.to_bits(), 0x4157_851f);
        assert_eq!(fog_far, 15.0);
    }

    #[test]
    fn clipping_preserves_snapshot_radius_max_and_normalized_far_for_fog() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,0.01],"target":[0,0,0],"up":[0,1,0],"radius":0.01,"radius-max":0.01}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 800, 600).unwrap();
        assert_eq!(camera.radius_max.to_bits(), 0.01f32.to_bits());
        assert_eq!(camera.near64.to_bits(), 0.01f64.to_bits());
        assert_eq!(camera.far64.to_bits(), 1.0f64.to_bits());
        let (fog_near, fog_far) = camera.fog_range(15.0);
        assert_eq!(fog_near.to_bits(), 0.017f32.to_bits());
        assert_eq!(fog_far.to_bits(), 1.0f32.to_bits());

        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":1}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 800, 600).unwrap();
        assert_eq!(camera.radius, 2.0);
        assert_eq!(camera.radius_max, 1.0);
    }

    #[test]
    fn pinned_xyz_projection_view_matches_browser_uniform_bits() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[-0.03077494221759732,0.03827120753606447,21.39274083202126],"target":[-0.03077494221759732,0.03827120753606447,0.00003476357507388535],"up":[0,1,0],"radius":4.186634185850473,"radius-max":4.186634185850473}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 1024, 937).unwrap();
        let expected = [
            0x400d_61e4,
            0,
            0,
            0,
            0,
            0x401a_827a,
            0,
            0,
            0,
            0,
            0xc0a3_832c,
            0xbf80_0000,
            0x3d8b_3bad,
            0xbdbd_3985,
            0x4085_fa5d,
            0x41ab_2455,
        ];
        assert_eq!(camera.projection_view_matrix().map(f32::to_bits), expected);
    }

    #[test]
    fn pinned_xyz_snapshot_camera_matches_browser_number_state() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[-0.03077494221759732,0.03827120753606447,21.39274083202126],"target":[-0.03077494221759732,0.03827120753606447,0.00003476357507388535],"up":[0,1,0],"radius":4.186634185850473,"radius-max":4.186634185850473}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 1024, 937).unwrap();
        assert_eq!(
            camera.near64.to_bits(),
            17.206_071_882_595_715_f64.to_bits()
        );
        assert_eq!(camera.far64.to_bits(), 25.579_340_254_296_66_f64.to_bits());
        let expected_view = [
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            0.0,
            0.030_774_942_217_597_32,
            -0.038_271_207_536_064_47,
            -21.392_740_832_021_26,
            1.0,
        ];
        let expected_projection = [
            2.209_099_714_788_662,
            0.0,
            0.0,
            0.0,
            0.0,
            2.414_213_562_373_095,
            0.0,
            0.0,
            0.0,
            0.0,
            -5.109_762_429_387_051,
            -1.0,
            0.0,
            0.0,
            -105.125_011_545_616_22,
            0.0,
        ];
        assert_eq!(
            camera.view_matrix64.map(f64::to_bits),
            expected_view.map(f64::to_bits)
        );
        assert_eq!(
            camera.projection_matrix64().map(f64::to_bits),
            expected_projection.map(f64::to_bits)
        );
        let (fog_near, fog_far) = camera.fog_range(15.0);
        assert_eq!(
            fog_near.to_bits(),
            (24.323_349_998_541_52_f64 as f32).to_bits()
        );
        assert_eq!(
            fog_far.to_bits(),
            (25.579_340_254_296_66_f64 as f32).to_bits()
        );
    }

    #[test]
    fn pinned_pdb_snapshot_camera_matches_browser_number_state() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[131.63666556103,125.63926427158606,292.5377510373967],"target":[131.63666556103,125.63926427158606,135.88791972158091],"up":[0,1,0],"radius":55.94729512734871,"radius-max":55.94729512734871}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 1024, 937).unwrap();
        assert_eq!(
            camera.near64.to_bits(),
            100.702_536_188_467_05_f64.to_bits()
        );
        assert_eq!(camera.far64.to_bits(), 212.597_126_443_164_48_f64.to_bits());
        let expected_projection = [
            2.209_099_714_788_661_7,
            0.0,
            0.0,
            0.0,
            0.0,
            2.414_213_562_373_095,
            0.0,
            0.0,
            0.0,
            0.0,
            -2.799_953_616_332_036_4,
            -1.0,
            0.0,
            0.0,
            -382.664_966_563_173_15,
            0.0,
        ];
        assert_eq!(
            camera.projection_matrix64().map(f64::to_bits),
            expected_projection.map(f64::to_bits)
        );
        let expected_projection_view = [
            2.209_099_714_788_661_7,
            0.0,
            0.0,
            0.0,
            0.0,
            2.414_213_562_373_095,
            0.0,
            0.0,
            0.0,
            0.0,
            -2.799_953_616_332_036_4,
            -1.0,
            -290.798_520_346_601_8,
            -303.320_015_771_040_46,
            436.427_167_367_626_57,
            292.537_751_037_396_7,
        ];
        assert_eq!(
            camera.projection_view_matrix().map(f32::to_bits),
            expected_projection_view.map(|value| (value as f32).to_bits())
        );
    }

    #[test]
    fn oblique_view_uses_one_number_precision_look_at_uniform_boundary() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[12.345678901,-3.210987654,9.876543219],"target":[1.234567891,2.345678912,-0.456789123],"up":[0.125,1,0.25],"radius":3,"radius-max":5}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 701, 997).unwrap();
        let expected = [
            0x3f33_a9ac,
            0x3e3e_9eaa,
            0x3f30_077c,
            0,
            0x3db6_3683,
            0x3f6f_4e70,
            0xbeb0_107f,
            0,
            0xbf34_f018,
            0x3e9a_e2ed,
            0x3f23_b509,
            0,
            0xbfb2_f106,
            0xc012_32d1,
            0xc17e_8bb5,
            0x3f80_0000,
        ];
        assert_eq!(camera.view_matrix().map(f32::to_bits), expected);
        assert_eq!(camera.right.x.to_bits(), expected[0]);
        assert_eq!(camera.right.y.to_bits(), expected[4]);
        assert_eq!(camera.right.z.to_bits(), expected[8]);
        assert_eq!(camera.up.x.to_bits(), expected[1]);
        assert_eq!(camera.up.y.to_bits(), expected[5]);
        assert_eq!(camera.up.z.to_bits(), expected[9]);
        assert_eq!(
            camera.forward.x.to_bits(),
            (-f32::from_bits(expected[2])).to_bits()
        );
        assert_eq!(
            camera.forward.y.to_bits(),
            (-f32::from_bits(expected[6])).to_bits()
        );
        assert_eq!(
            camera.forward.z.to_bits(),
            (-f32::from_bits(expected[10])).to_bits()
        );
    }

    #[test]
    fn auto_fit_preserves_molstar_number_precision_and_aspect_rules() {
        let sphere = BoundingSphere {
            center: Vec3::new(1.0, 2.0, 3.0),
            radius: 2.0,
            center64: Some([1.000_000_03, 2.000_000_06, 3.000_000_09]),
            radius64: Some(2.000_000_12),
            ..BoundingSphere::default()
        };
        let perspective = RendererOptions::default();
        let landscape = resolve_camera(&perspective, &sphere, 128, 112).unwrap();
        let portrait = resolve_camera(&perspective, &sphere, 112, 128).unwrap();
        let expected_landscape =
            sphere.radius64() / js_sin64(perspective.camera.fov.to_radians() * 0.5);
        let expected_portrait = sphere.radius64()
            / (112.0 / 128.0)
            / js_sin64(perspective.camera.fov.to_radians() * 0.5);
        assert_eq!(landscape.distance64.to_bits(), expected_landscape.to_bits());
        assert_eq!(portrait.distance64.to_bits(), expected_portrait.to_bits());
        assert_eq!(
            landscape.view_matrix64[12].to_bits(),
            (-1.000_000_03f64).to_bits()
        );

        let orthographic =
            RendererOptions::from_json(br#"{"camera":{"mode":"orthographic"}}"#).unwrap();
        let camera = resolve_camera(&orthographic, &sphere, 128, 112).unwrap();
        let expected = sphere.radius64() / js_tan64(orthographic.camera.fov.to_radians() * 0.5);
        assert_eq!(camera.distance64.to_bits(), expected.to_bits());
    }

    #[test]
    fn orthographic_projection_uniform_has_the_molstar_zero_structure() {
        let renderer =
            RendererOptions::from_json(br#"{"camera":{"mode":"orthographic"}}"#).unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 701, 997).unwrap();
        assert_eq!(
            camera
                .projection_matrix64()
                .map(|value| (value as f32).to_bits()),
            camera.projection_matrix().map(f32::to_bits)
        );
        assert_eq!(camera.projection_matrix64()[6].to_bits(), 0);
    }

    #[test]
    fn orthographic_projection_and_inverse_match_pinned_v8_staging() {
        let renderer =
            RendererOptions::from_json(br#"{"camera":{"mode":"orthographic"}}"#).unwrap();
        let sphere = BoundingSphere {
            center: Vec3::new(1.0, 2.0, 3.0),
            radius: 2.0,
            center64: Some([1.000_000_03, 2.000_000_06, 3.000_000_09]),
            radius64: Some(2.000_000_12),
            ..BoundingSphere::default()
        };
        let camera = resolve_camera(&renderer, &sphere, 701, 997)
            .unwrap()
            .with_view_offset(0.375, -0.125);
        assert_eq!(camera.distance64.to_bits(), 0x401b_7810_5707_0775);
        assert_eq!(camera.near64.to_bits(), 0x4013_7810_4ef9_71e0);
        assert_eq!(camera.far64.to_bits(), 0x4021_bc08_2f8a_4e85);
        assert_eq!(
            camera.projection_matrix().map(f32::to_bits),
            [
                0x3eff_ffff,
                0,
                0,
                0,
                0,
                0x3eb3_fef8,
                0,
                0,
                0,
                0,
                0xbeff_ffff,
                0,
                0xba8c_3be4,
                0xb983_7766,
                0xc05b_c082,
                0x3f80_0000,
            ]
        );
        assert_eq!(
            camera.inverse_projection_matrix().map(f32::to_bits),
            [
                0x4000_0001,
                0x8000_0000,
                0x8000_0000,
                0,
                0x8000_0000,
                0x4036_0c6b,
                0x8000_0000,
                0x8000_0000,
                0x8000_0000,
                0x8000_0000,
                0xc000_0001,
                0x8000_0000,
                0x3b0c_3be5,
                0x3a3a_fa86,
                0xc0db_c083,
                0x3f80_0000,
            ]
        );
    }

    #[test]
    fn impostor_depth_uses_the_uploaded_projection_uniform() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":5}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 701, 997).unwrap();
        let projection = camera.projection_matrix();
        let view_z = -9.125_f32;
        let clip_z = fused_multiply_add(view_z, projection[10], projection[14]);
        let clip_w = fused_multiply_add(view_z, projection[11], projection[15]);
        let expected = fused_multiply_add(divide(clip_z, clip_w), 0.5, 0.5);
        assert_eq!(
            camera.impostor_depth01(view_z).to_bits(),
            expected.to_bits()
        );
    }

    #[test]
    fn ssao_sample_projection_uses_shader_bounds_without_near_plane_rejection() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":5}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 10, 10).unwrap();

        // The SSAO fragment shader applies only projected uBounds. It does not
        // discard a sample merely because its view-space Z crossed the camera.
        assert_eq!(
            camera.project_view_position(Vec3::new(0.0, 0.0, 1.0), 10, 10, [0.0, 0.0, 1.0, 1.0]),
            Some((5, 4))
        );

        let projection = camera.projection_matrix();
        let x = 1.1 / projection[0];
        let sample = Vec3::new(x, 0.0, -1.0);
        assert_eq!(
            camera.project_view_position(sample, 10, 10, [0.0, 0.0, 1.0, 1.0]),
            None
        );
        assert_eq!(
            camera.project_view_position(sample, 10, 10, [0.0, 0.0, 1.1, 1.0]),
            Some((9, 4))
        );
    }

    #[test]
    fn depth_decode_uses_molstar_shader_operation_order() {
        for mode in ["perspective", "orthographic"] {
            let source = format!(r#"{{"camera":{{"mode":"{mode}"}}}}"#);
            let renderer = RendererOptions::from_json(source.as_bytes()).unwrap();
            let camera = resolve_camera(&renderer, &BoundingSphere::default(), 701, 997).unwrap();
            let depth = f32::from_bits(0x3f42_1357);
            let expected_view_z = if camera.orthographic {
                fused_multiply_add(depth, camera.near - camera.far, -camera.near)
            } else {
                divide(
                    camera.near * camera.far,
                    fused_multiply_add(camera.far - camera.near, depth, -camera.far),
                )
            };
            assert_eq!(
                camera.depth_from_depth01(depth).to_bits(),
                (-expected_view_z).to_bits()
            );
        }
    }
}
