use crate::model::Vec3;

use super::super::camera::CameraState;
use super::super::color::{color_f32, parse_color, quantize};
use super::super::framebuffer::Framebuffer;
use super::super::options::OcclusionParams;
use super::angle_metal_math::{divide, divide_blur_sum, fused_multiply_add, inverse_sqrt};
use super::fog::fog_factor;
#[cfg(test)]
use super::packing::pack_depth_alpha_rgba8;
use super::packing::{
    pack_unit_interval_rgba8, packed_unit_interval_roundtrip, unpack_depth_alpha_rgba,
    unpack_unit_interval_rg,
};
use super::ssao_sinc::sinc_approximation;

fn framebuffer_depth01_at(framebuffer: &Framebuffer, x: isize, y: isize) -> f32 {
    let sx = x.clamp(0, framebuffer.width.saturating_sub(1) as isize) as usize;
    let sy = y.clamp(0, framebuffer.height.saturating_sub(1) as isize) as usize;
    framebuffer.depth01[sy * framebuffer.width + sx]
}

fn angle_metal_cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        fused_multiply_add(a.y, b.z, -(a.z * b.y)),
        fused_multiply_add(a.z, b.x, -(a.x * b.z)),
        fused_multiply_add(a.x, b.y, -(a.y * b.x)),
    )
}

fn angle_metal_dot3(a: Vec3, b: Vec3) -> f32 {
    fused_multiply_add(b.z, a.z, fused_multiply_add(b.y, a.y, a.x * b.x))
}

fn angle_metal_normalize(value: Vec3) -> Vec3 {
    let length_squared = angle_metal_dot3(value, value);
    if length_squared <= 0.000_000_000_001 {
        Vec3::default()
    } else {
        value * inverse_sqrt(length_squared)
    }
}

fn angle_metal_tbn_transform(tangent: Vec3, bitangent: Vec3, normal: Vec3, sample: Vec3) -> Vec3 {
    let component = |tangent: f32, bitangent: f32, normal: f32| {
        fused_multiply_add(
            normal,
            sample.z,
            fused_multiply_add(bitangent, sample.y, tangent * sample.x),
        )
    };
    Vec3::new(
        component(tangent.x, bitangent.x, normal.x),
        component(tangent.y, bitangent.y, normal.y),
        component(tangent.z, bitangent.z, normal.z),
    )
}

/// Reconstruct the view-space normal from the depth buffer exactly as Mol*'s
/// SSAO pass does. In particular, derivative selection uses non-linear WebGL
/// depth while the cross product uses reconstructed view-space positions.
fn reconstructed_view_normal_basis(
    framebuffer: &Framebuffer,
    camera: &CameraState,
    x: usize,
    y: usize,
) -> Vec3 {
    let x = x as isize;
    let y = y as isize;
    let inverse_width = super::angle_metal_math::reciprocal(framebuffer.width as f32);
    let inverse_height = super::angle_metal_math::reciprocal(framebuffer.height as f32);
    let fragment_x = x as f32 + 0.5;
    let fragment_y = framebuffer.height as f32 - y as f32 - 0.5;
    let coords_x = fragment_x * inverse_width;
    let coords_y = fragment_y * inverse_height;
    let position_at_coords = |sx: isize, sy: isize, cx: f32, cy: f32| {
        let sx = sx.clamp(0, framebuffer.width.saturating_sub(1) as isize) as usize;
        let sy = sy.clamp(0, framebuffer.height.saturating_sub(1) as isize) as usize;
        camera.screen_space_to_view_space_with_fused_w(
            cx,
            cy,
            framebuffer.depth01[sy * framebuffer.width + sx],
        )
    };
    let position_at_fragment = |sx: isize, sy: isize, fx: f32, fy: f32| {
        let sx = sx.clamp(0, framebuffer.width.saturating_sub(1) as isize) as usize;
        let sy = sy.clamp(0, framebuffer.height.saturating_sub(1) as isize) as usize;
        camera.view_position_at_neighbor_fragment(
            fx,
            fy,
            framebuffer.depth01[sy * framebuffer.width + sx],
            framebuffer.width,
            framebuffer.height,
        )
    };
    let fragment_center_position = camera.view_position_at_pixel(
        x as usize,
        y as usize,
        framebuffer.depth01[y as usize * framebuffer.width + x as usize],
        framebuffer.width,
        framebuffer.height,
    );
    let normalized_center_position = position_at_coords(x, y, coords_x, coords_y);
    // ANGLE's Metal lowering preserves the normalized texture-coordinate
    // arithmetic for the view-space Y component, while X/Z follow the
    // reassociated fragment-coordinate path on the pinned reference GPU.
    let center_position = Vec3::new(
        fragment_center_position.x,
        normalized_center_position.y,
        fragment_center_position.z,
    );
    let fragment_left_position = position_at_fragment(x - 1, y, fragment_x - 1.0, fragment_y);
    let fragment_right_position = position_at_fragment(x + 1, y, fragment_x + 1.0, fragment_y);
    let normalized_left_position = position_at_coords(x - 1, y, coords_x - inverse_width, coords_y);
    let normalized_right_position =
        position_at_coords(x + 1, y, coords_x + inverse_width, coords_y);
    let left_position = Vec3::new(
        fragment_left_position.x,
        normalized_left_position.y,
        fragment_left_position.z,
    );
    let right_position = Vec3::new(
        fragment_right_position.x,
        normalized_right_position.y,
        fragment_right_position.z,
    );
    // GLSL texture coordinates grow upward, whereas framebuffer rows grow
    // downward. Therefore shader "down" is y + 1 in this buffer.
    let down_position = position_at_coords(x, y + 1, coords_x, coords_y - inverse_height);
    let up_position = position_at_coords(x, y - 1, coords_x, coords_y + inverse_height);

    let left = center_position - left_position;
    let right = right_position - center_position;
    let down = center_position - down_position;
    let up = up_position - center_position;

    let depth01 = |dx: isize, dy: isize| framebuffer_depth01_at(framebuffer, x + dx, y + dy);
    let center_depth = depth01(0, 0);
    let horizontal_error_left =
        fused_multiply_add(depth01(-1, 0), 2.0, -depth01(-2, 0) - center_depth).abs();
    let horizontal_error_right =
        fused_multiply_add(depth01(1, 0), 2.0, -depth01(2, 0) - center_depth).abs();
    let vertical_error_down =
        fused_multiply_add(depth01(0, 1), 2.0, -depth01(0, 2) - center_depth).abs();
    let vertical_error_up =
        fused_multiply_add(depth01(0, -1), 2.0, -depth01(0, -2) - center_depth).abs();

    let horizontal = if horizontal_error_left < horizontal_error_right {
        left
    } else {
        right
    };
    let vertical = if vertical_error_down < vertical_error_up {
        down
    } else {
        up
    };
    angle_metal_cross(horizontal, vertical)
}

fn ssao_background_factor() -> f32 {
    packed_unit_interval_roundtrip(1.0)
}

fn blur_sample_coordinates(x: usize, y: usize, offset: isize, horizontal: bool) -> (isize, isize) {
    if horizontal {
        (x as isize + offset, y as isize)
    } else {
        // GLSL texture coordinates grow upward, while this framebuffer is
        // stored top-down. Reversing the row offset also preserves the
        // shader's accumulation order.
        (x as isize, y as isize - offset)
    }
}

fn accumulate_blur_sample(sum: f32, sample: f32, weight: f32) -> f32 {
    fused_multiply_add(sample, weight, sum)
}

#[derive(Clone)]
pub(in crate::render) struct OcclusionFactors {
    width: usize,
    height: usize,
    values: Vec<f32>,
}

pub(in crate::render) fn ssao_target_dimensions(
    width: usize,
    height: usize,
    pixel_ratio: f64,
    resolution_scale: f64,
) -> (usize, usize, f64) {
    // SsaoPass.calcSsaoScale downscales high-DPI drawing buffers to CSS-pixel
    // density before applying the author-selected resolution scale.
    let scale = (1.0 / pixel_ratio).min(1.0) * resolution_scale;
    let target_width = ((width as f64 * scale).floor() as usize).max(1);
    let target_height = ((height as f64 * scale).floor() as usize).max(1);
    (target_width, target_height, scale)
}

fn ssao_full_viewport_bounds(width: usize, height: usize, scale: f64) -> [f32; 4] {
    // `SsaoPass.update` derives uBounds from the unscaled Canvas3D viewport,
    // retaining the fractional product in the denominator even though the
    // SSAO target itself was floor-sized. For a full viewport this can make
    // the upper bound slightly larger than one.
    let scaled_width = width as f64 * scale;
    let scaled_height = height as f64 * scale;
    [
        0.0,
        0.0,
        (scaled_width.ceil() / scaled_width) as f32,
        (scaled_height.ceil() / scaled_height) as f32,
    ]
}

fn nearest_scaled_values(
    source: &[f32],
    source_width: usize,
    source_height: usize,
    target_width: usize,
    target_height: usize,
) -> Vec<f32> {
    assert_eq!(source.len(), source_width * source_height);
    if source_width == target_width && source_height == target_height {
        return source.to_vec();
    }
    let mut target = Vec::with_capacity(target_width * target_height);
    for y in 0..target_height {
        let source_y = ((2 * y + 1) * source_height / (2 * target_height))
            .min(source_height.saturating_sub(1));
        for x in 0..target_width {
            let source_x = ((2 * x + 1) * source_width / (2 * target_width))
                .min(source_width.saturating_sub(1));
            target.push(source[source_y * source_width + source_x]);
        }
    }
    target
}

fn scaled_depth_framebuffer(
    framebuffer: &Framebuffer,
    target_width: usize,
    target_height: usize,
) -> Framebuffer {
    let depth01 = nearest_scaled_values(
        &framebuffer.depth01,
        framebuffer.width,
        framebuffer.height,
        target_width,
        target_height,
    );
    Framebuffer {
        width: target_width,
        height: target_height,
        color: vec![[0; 4]; target_width * target_height],
        depth: vec![f32::INFINITY; target_width * target_height],
        depth01,
        normal: vec![Vec3::default(); target_width * target_height],
    }
}

fn nearest_depth_at_offset<T: Copy>(
    values: &[T],
    width: usize,
    height: usize,
    [offset_x, offset_y]: [f32; 2],
) -> T {
    let gl_x = (offset_x * width as f32).floor() as usize;
    let gl_y = (offset_y * height as f32).floor() as usize;
    let x = gl_x.min(width.saturating_sub(1));
    let y = height - 1 - gl_y.min(height.saturating_sub(1));
    values[y * width + x]
}

fn linear_depth_at_offset(
    values: &[f32],
    width: usize,
    height: usize,
    [offset_x, offset_y]: [f32; 2],
) -> f32 {
    let sample_x = fused_multiply_add(offset_x, width as f32, -0.5);
    let sample_y = fused_multiply_add(1.0 - offset_y, height as f32, -0.5);
    let x0 = sample_x.floor() as isize;
    let y0 = sample_y.floor() as isize;
    let tx = sample_x - x0 as f32;
    let ty = sample_y - y0 as f32;
    let at = |x: isize, y: isize| {
        let x = x.clamp(0, width.saturating_sub(1) as isize) as usize;
        let y = y.clamp(0, height.saturating_sub(1) as isize) as usize;
        values[y * width + x]
    };
    let upper = fused_multiply_add(at(x0 + 1, y0) - at(x0, y0), tx, at(x0, y0));
    let lower = fused_multiply_add(at(x0 + 1, y0 + 1) - at(x0, y0 + 1), tx, at(x0, y0 + 1));
    fused_multiply_add(lower - upper, ty, upper)
}

fn copy_depth_target(
    source: &[f32],
    source_width: usize,
    source_height: usize,
    target_width: usize,
    target_height: usize,
    source_linear: bool,
) -> Vec<f32> {
    let sample = if source_linear {
        linear_depth_at_offset
    } else {
        nearest_depth_at_offset
    };
    let mut target = Vec::with_capacity(target_width * target_height);
    for y in 0..target_height {
        let offset_y = 1.0 - divide(y as f32 + 0.5, target_height as f32);
        for x in 0..target_width {
            let offset_x = divide(x as f32 + 0.5, target_width as f32);
            target.push(sample(
                source,
                source_width,
                source_height,
                [offset_x, offset_y],
            ));
        }
    }
    target
}

struct DepthPyramid<'a> {
    full: &'a [f32],
    full_width: usize,
    full_height: usize,
    full_linear: bool,
    half: Vec<f32>,
    half_width: usize,
    half_height: usize,
    quarter: Vec<f32>,
    quarter_width: usize,
    quarter_height: usize,
}

impl<'a> DepthPyramid<'a> {
    fn new(full: &'a [f32], width: usize, height: usize, full_linear: bool) -> Self {
        let half_width = (width / 2).max(1);
        let half_height = (height / 2).max(1);
        let half = copy_depth_target(full, width, height, half_width, half_height, full_linear);
        let quarter_width = (width / 4).max(1);
        let quarter_height = (height / 4).max(1);
        let quarter = copy_depth_target(
            &half,
            half_width,
            half_height,
            quarter_width,
            quarter_height,
            true,
        );
        Self {
            full,
            full_width: width,
            full_height: height,
            full_linear,
            half,
            half_width,
            half_height,
            quarter,
            quarter_width,
            quarter_height,
        }
    }

    fn sample_full(&self, offset: [f32; 2]) -> f32 {
        if self.full_linear {
            linear_depth_at_offset(self.full, self.full_width, self.full_height, offset)
        } else {
            nearest_depth_at_offset(self.full, self.full_width, self.full_height, offset)
        }
    }

    fn sample_mapped(&self, offset: [f32; 2], self_offset: [f32; 2]) -> f32 {
        let dx = offset[0] - self_offset[0];
        let dy = offset[1] - self_offset[1];
        let distance = fused_multiply_add(dy, dy, dx * dx).sqrt();
        if distance > 0.1 {
            linear_depth_at_offset(
                &self.quarter,
                self.quarter_width,
                self.quarter_height,
                offset,
            )
        } else if distance > 0.05 {
            linear_depth_at_offset(&self.half, self.half_width, self.half_height, offset)
        } else {
            self.sample_full(offset)
        }
    }
}

fn nearest_rgba8_at_offset(
    values: &[[u8; 4]],
    width: usize,
    height: usize,
    offset: [f32; 2],
) -> [f32; 4] {
    let packed = nearest_depth_at_offset(values, width, height, offset);
    packed.map(|value| value as f32 / 255.0)
}

fn linear_rgba8_at_offset(
    values: &[[u8; 4]],
    width: usize,
    height: usize,
    [offset_x, offset_y]: [f32; 2],
) -> [f32; 4] {
    let sample_x = fused_multiply_add(offset_x, width as f32, -0.5);
    let sample_y = fused_multiply_add(1.0 - offset_y, height as f32, -0.5);
    let x0 = sample_x.floor() as isize;
    let y0 = sample_y.floor() as isize;
    let tx = sample_x - x0 as f32;
    let ty = sample_y - y0 as f32;
    let at = |x: isize, y: isize| {
        let x = x.clamp(0, width.saturating_sub(1) as isize) as usize;
        let y = y.clamp(0, height.saturating_sub(1) as isize) as usize;
        values[y * width + x]
    };
    let upper_left = at(x0, y0);
    let upper_right = at(x0 + 1, y0);
    let lower_left = at(x0, y0 + 1);
    let lower_right = at(x0 + 1, y0 + 1);
    std::array::from_fn(|channel| {
        let upper_left = upper_left[channel] as f32;
        let lower_left = lower_left[channel] as f32;
        let upper = fused_multiply_add(upper_right[channel] as f32 - upper_left, tx, upper_left);
        let lower = fused_multiply_add(lower_right[channel] as f32 - lower_left, tx, lower_left);
        let byte_value = fused_multiply_add(lower - upper, ty, upper);
        (byte_value * 16.0).round() / (16.0 * 255.0)
    })
}

fn copy_rgba8_target(
    source: &[[u8; 4]],
    source_width: usize,
    source_height: usize,
    target_width: usize,
    target_height: usize,
    source_linear: bool,
) -> Vec<[u8; 4]> {
    let sample = if source_linear {
        linear_rgba8_at_offset
    } else {
        nearest_rgba8_at_offset
    };
    let mut target = Vec::with_capacity(target_width * target_height);
    for y in 0..target_height {
        let offset_y = 1.0 - divide(y as f32 + 0.5, target_height as f32);
        for x in 0..target_width {
            let offset_x = divide(x as f32 + 0.5, target_width as f32);
            target.push(
                sample(source, source_width, source_height, [offset_x, offset_y]).map(quantize),
            );
        }
    }
    target
}

struct PackedDepthAlphaPyramid {
    full: Vec<[u8; 4]>,
    full_width: usize,
    full_height: usize,
    full_linear: bool,
    half: Vec<[u8; 4]>,
    half_width: usize,
    half_height: usize,
    quarter: Vec<[u8; 4]>,
    quarter_width: usize,
    quarter_height: usize,
}

impl PackedDepthAlphaPyramid {
    fn new(
        packed: &[[u8; 4]],
        width: usize,
        height: usize,
        full_linear: bool,
        mapped_levels: bool,
    ) -> Self {
        assert_eq!(packed.len(), width * height);
        let full = packed.to_vec();
        let (half, half_width, half_height, quarter, quarter_width, quarter_height) =
            if mapped_levels {
                let half_width = (width / 2).max(1);
                let half_height = (height / 2).max(1);
                let half =
                    copy_rgba8_target(&full, width, height, half_width, half_height, full_linear);
                let quarter_width = (width / 4).max(1);
                let quarter_height = (height / 4).max(1);
                let quarter = copy_rgba8_target(
                    &half,
                    half_width,
                    half_height,
                    quarter_width,
                    quarter_height,
                    true,
                );
                (
                    half,
                    half_width,
                    half_height,
                    quarter,
                    quarter_width,
                    quarter_height,
                )
            } else {
                (Vec::new(), 0, 0, Vec::new(), 0, 0)
            };
        Self {
            full,
            full_width: width,
            full_height: height,
            full_linear,
            half,
            half_width,
            half_height,
            quarter,
            quarter_width,
            quarter_height,
        }
    }

    fn sample_full(&self, offset: [f32; 2]) -> (f32, f32) {
        let rgba = if self.full_linear {
            linear_rgba8_at_offset(&self.full, self.full_width, self.full_height, offset)
        } else {
            nearest_rgba8_at_offset(&self.full, self.full_width, self.full_height, offset)
        };
        unpack_depth_alpha_rgba(rgba)
    }

    fn sample_mapped(&self, offset: [f32; 2], self_offset: [f32; 2]) -> (f32, f32) {
        let dx = offset[0] - self_offset[0];
        let dy = offset[1] - self_offset[1];
        let distance = fused_multiply_add(dy, dy, dx * dx).sqrt();
        let rgba = if distance > 0.1 {
            linear_rgba8_at_offset(
                &self.quarter,
                self.quarter_width,
                self.quarter_height,
                offset,
            )
        } else if distance > 0.05 {
            linear_rgba8_at_offset(&self.half, self.half_width, self.half_height, offset)
        } else {
            return self.sample_full(offset);
        };
        unpack_depth_alpha_rgba(rgba)
    }
}

pub(in crate::render) fn compute_occlusion_factors(
    framebuffer: &Framebuffer,
    camera: &CameraState,
    params: &OcclusionParams,
    pixel_ratio: f64,
) -> OcclusionFactors {
    let (width, height, scale) = ssao_target_dimensions(
        framebuffer.width,
        framebuffer.height,
        pixel_ratio,
        params.resolution_scale,
    );
    let bounds = ssao_full_viewport_bounds(framebuffer.width, framebuffer.height, scale);
    let scaled = scaled_depth_framebuffer(framebuffer, width, height);
    OcclusionFactors {
        width,
        height,
        values: compute_occlusion_factors_internal(
            &scaled,
            camera,
            params,
            bounds,
            scale != 1.0,
            None,
            None,
        ),
    }
}

pub(in crate::render) fn compute_occlusion_factors_including_transparency(
    opaque: &Framebuffer,
    transparent_depth_alpha_rgba8: &[[u8; 4]],
    camera: &CameraState,
    params: &OcclusionParams,
    pixel_ratio: f64,
) -> (OcclusionFactors, OcclusionFactors) {
    assert_eq!(
        transparent_depth_alpha_rgba8.len(),
        opaque.width * opaque.height
    );
    let (width, height, scale) = ssao_target_dimensions(
        opaque.width,
        opaque.height,
        pixel_ratio,
        params.resolution_scale,
    );
    let bounds = ssao_full_viewport_bounds(opaque.width, opaque.height, scale);
    let scaled_opaque = scaled_depth_framebuffer(opaque, width, height);
    // The source transparent-depth texture is nearest-filtered and copied
    // into an RGBA8 linear target when SSAO is downscaled. Keep the packed
    // clear/background bytes intact across that copy and unpack only after
    // the target attachment conversion.
    let scaled_transparent_packed = if width == opaque.width && height == opaque.height {
        transparent_depth_alpha_rgba8.to_vec()
    } else {
        copy_rgba8_target(
            transparent_depth_alpha_rgba8,
            opaque.width,
            opaque.height,
            width,
            height,
            false,
        )
    };
    let scaled_transparent_depth = scaled_transparent_packed
        .iter()
        .map(|packed| unpack_depth_alpha_rgba(packed.map(|channel| f32::from(channel) / 255.0)).0)
        .collect::<Vec<_>>();
    let transparent = Framebuffer {
        width,
        height,
        color: vec![[0; 4]; width * height],
        depth: vec![f32::INFINITY; width * height],
        depth01: scaled_transparent_depth.clone(),
        normal: vec![Vec3::default(); width * height],
    };
    let transparent_occluder = Some(scaled_transparent_packed.as_slice());
    (
        OcclusionFactors {
            width,
            height,
            values: compute_occlusion_factors_internal(
                &scaled_opaque,
                camera,
                params,
                bounds,
                scale != 1.0,
                Some(&scaled_opaque.depth01),
                transparent_occluder,
            ),
        },
        OcclusionFactors {
            width,
            height,
            values: compute_occlusion_factors_internal(
                &transparent,
                camera,
                params,
                bounds,
                scale != 1.0,
                Some(&scaled_opaque.depth01),
                transparent_occluder,
            ),
        },
    )
}

fn compute_occlusion_factors_internal(
    framebuffer: &Framebuffer,
    camera: &CameraState,
    params: &OcclusionParams,
    bounds: [f32; 4],
    full_depth_linear: bool,
    opaque_depth01: Option<&[f32]>,
    transparent_depth_alpha: Option<&[[u8; 4]]>,
) -> Vec<f32> {
    let opaque_depth01 = opaque_depth01.unwrap_or(&framebuffer.depth01);
    let width = framebuffer.width;
    let height = framebuffer.height;
    // The SSAO shader writes every pixel, including background, through an
    // RG8 target before either blur pass samples it. Preserve that render-
    // target roundtrip for the background value as well.
    let mut factors = vec![ssao_background_factor(); width * height];
    let samples = molstar_ssao_samples(params.samples.clamp(1, 256));
    let radius = ssao_uniform_radius(params, camera);
    let bias = params.bias as f32;
    let multi_scale = params.multi_scale.name == "on";
    let opaque_pyramid =
        multi_scale.then(|| DepthPyramid::new(opaque_depth01, width, height, full_depth_linear));
    let transparent_pyramid = transparent_depth_alpha.map(|packed| {
        PackedDepthAlphaPyramid::new(packed, width, height, full_depth_linear, multi_scale)
    });
    let mut levels = params.multi_scale.params.levels.clone();
    levels.sort_by(|a, b| a.radius.partial_cmp(&b.radius).unwrap());
    let levels = levels
        .iter()
        .map(|level| {
            (
                (2.0f64.powf(level.radius) * camera.scale()) as f32,
                level.bias as f32,
            )
        })
        .collect::<Vec<_>>();
    let near_threshold = params.multi_scale.params.near_threshold as f32;
    let far_threshold = params.multi_scale.params.far_threshold as f32;
    for y in 0..height {
        for x in 0..width {
            let index = y * width + x;
            let center_depth01 = framebuffer.depth01[index];
            // Mol*'s `isBackground` helper tests the depth sentinel exactly;
            // values merely greater than one are not silently reclassified.
            if center_depth01 == 1.0 {
                continue;
            }
            let self_view_normal_cross = reconstructed_view_normal_basis(framebuffer, camera, x, y);
            let self_view_position =
                camera.view_position_at_pixel(x, y, center_depth01, width, height);
            let self_view_normal = angle_metal_normalize(self_view_normal_cross);
            if self_view_normal.squared_length() <= 0.000_001 {
                continue;
            }
            let random_raw = ssao_noise_raw_vector(x, y, width, height);
            let random_length_squared = angle_metal_dot3(random_raw, random_raw);
            let random_inverse_length = inverse_sqrt(random_length_squared);
            let random = random_raw * random_inverse_length;
            let tangent_projection = angle_metal_dot3(random, self_view_normal);
            let tangent_raw = Vec3::new(
                fused_multiply_add(-self_view_normal.x, tangent_projection, random.x),
                fused_multiply_add(-self_view_normal.y, tangent_projection, random.y),
                random.z - self_view_normal.z * tangent_projection,
            );
            let tangent_length_squared = angle_metal_dot3(tangent_raw, tangent_raw);
            let tangent_inverse_length = inverse_sqrt(tangent_length_squared);
            let mut tangent = tangent_raw * tangent_inverse_length;
            if tangent.squared_length() <= 0.000_001 {
                tangent = if self_view_normal.x.abs() < 0.9 {
                    self_view_normal
                        .cross(Vec3::new(1.0, 0.0, 0.0))
                        .normalized()
                } else {
                    self_view_normal
                        .cross(Vec3::new(0.0, 1.0, 0.0))
                        .normalized()
                };
            }
            let bitangent = angle_metal_cross(self_view_normal, tangent);
            if multi_scale {
                let self_offset = [
                    divide(x as f32 + 0.5, width as f32),
                    divide((height - y) as f32 - 0.5, height as f32),
                ];
                let neighbor_offset = [self_offset[0] + divide(1.0, width as f32), self_offset[1]];
                let pixel_delta = camera.screen_space_to_view_space(
                    neighbor_offset[0],
                    neighbor_offset[1],
                    center_depth01,
                ) - camera.screen_space_to_view_space(
                    self_offset[0],
                    self_offset[1],
                    center_depth01,
                );
                let pixel_size = angle_metal_dot3(pixel_delta, pixel_delta).sqrt();
                let opaque_pyramid = opaque_pyramid.as_ref().expect("multi-scale depth pyramid");
                let mut occluded = 0.0f32;
                for &(level_radius, level_bias) in &levels {
                    if pixel_size * near_threshold > level_radius
                        || pixel_size * far_threshold < level_radius
                    {
                        continue;
                    }
                    let mut level_occlusion = 0.0;
                    let mut valid = samples.len() as f32;
                    for sample in &samples {
                        let view_offset = angle_metal_tbn_transform(
                            tangent,
                            bitangent,
                            self_view_normal,
                            *sample,
                        );
                        let sample_view_position = Vec3::new(
                            fused_multiply_add(view_offset.x, level_radius, self_view_position.x),
                            fused_multiply_add(view_offset.y, level_radius, self_view_position.y),
                            fused_multiply_add(view_offset.z, level_radius, self_view_position.z),
                        );
                        let Some(sample_offset) =
                            camera.project_view_position_offset(sample_view_position, bounds)
                        else {
                            valid -= 1.0;
                            continue;
                        };
                        let opaque_depth = opaque_pyramid.sample_mapped(sample_offset, self_offset);
                        let mut sample_factor = depth_sample_occlusion(
                            camera,
                            opaque_depth,
                            sample_view_position,
                            self_view_position,
                            level_radius,
                        ) * level_bias;
                        if let Some(transparent) = transparent_pyramid.as_ref() {
                            let (transparent_depth, transparent_alpha) =
                                transparent.sample_mapped(sample_offset, self_offset);
                            let transparent_factor = depth_sample_occlusion(
                                camera,
                                transparent_depth,
                                sample_view_position,
                                self_view_position,
                                level_radius,
                            ) * level_bias
                                * transparent_alpha;
                            sample_factor = sample_factor.max(transparent_factor);
                        }
                        level_occlusion += sample_factor;
                    }
                    if valid > 0.0 {
                        occluded = occluded.max(divide(level_occlusion, valid));
                    }
                }
                factors[index] = packed_unit_interval_roundtrip(
                    fused_multiply_add(-bias, occluded, 1.0).clamp(0.01, 1.0),
                );
                continue;
            }
            let mut occluded = 0.0;
            let mut valid = samples.len() as f32;
            for sample in &samples {
                let view_offset =
                    angle_metal_tbn_transform(tangent, bitangent, self_view_normal, *sample);
                let sample_view_position = Vec3::new(
                    fused_multiply_add(view_offset.x, radius, self_view_position.x),
                    fused_multiply_add(view_offset.y, radius, self_view_position.y),
                    fused_multiply_add(view_offset.z, radius, self_view_position.z),
                );
                let sample_offset = if full_depth_linear {
                    let Some(offset) =
                        camera.project_view_position_offset(sample_view_position, bounds)
                    else {
                        valid -= 1.0;
                        continue;
                    };
                    offset
                } else {
                    let Some((sx, sy)) =
                        camera.project_view_position(sample_view_position, width, height, bounds)
                    else {
                        valid -= 1.0;
                        continue;
                    };
                    [
                        divide(sx as f32 + 0.5, width as f32),
                        divide((height - sy) as f32 - 0.5, height as f32),
                    ]
                };
                let sample_depth = if full_depth_linear {
                    linear_depth_at_offset(opaque_depth01, width, height, sample_offset)
                } else {
                    nearest_depth_at_offset(opaque_depth01, width, height, sample_offset)
                };
                let mut sample_factor = depth_sample_occlusion(
                    camera,
                    sample_depth,
                    sample_view_position,
                    self_view_position,
                    radius,
                );
                if transparent_depth_alpha.is_some() {
                    let (transparent_depth, transparent_alpha) = transparent_pyramid
                        .as_ref()
                        .expect("transparent depth/alpha pyramid")
                        .sample_full(sample_offset);
                    sample_factor = combine_opaque_and_transparent_occlusion(
                        sample_factor,
                        depth_sample_occlusion(
                            camera,
                            transparent_depth,
                            sample_view_position,
                            self_view_position,
                            radius,
                        ),
                        transparent_alpha,
                    );
                }
                occluded += sample_factor;
            }
            if valid > 0.0 {
                let average = divide(occluded, valid);
                factors[index] = packed_unit_interval_roundtrip(
                    fused_multiply_add(-bias, average, 1.0).clamp(0.01, 1.0),
                );
            }
        }
    }
    bilateral_blur_occlusion(
        &factors,
        &framebuffer.depth01,
        width,
        height,
        camera,
        params,
        bounds,
    )
}

/// Evaluate Mol*'s `uRadius = 2^radius * camera.scale` as JavaScript Number
/// arithmetic, then round once at the highp-float uniform boundary.
fn ssao_uniform_radius(params: &OcclusionParams, camera: &CameraState) -> f32 {
    (2.0f64.powf(params.radius) * camera.scale()) as f32
}

fn sample_range_attenuation(radius: f32, depth_delta: f32) -> f32 {
    // Keep the shader expression literal. IEEE division by zero produces
    // infinity, which `smootherstep` clamps to one; Mol* has no epsilon term.
    smootherstep(0.0, 1.0, divide(radius, depth_delta.abs()))
}

fn depth_sample_occlusion(
    camera: &CameraState,
    sample_depth01: f32,
    sample_view_position: Vec3,
    self_view_position: Vec3,
    radius: f32,
) -> f32 {
    if sample_depth01 == 1.0 {
        return 0.0;
    }
    let sample_view_z = camera.view_z_from_depth01(sample_depth01);
    if sample_view_z < sample_view_position.z + 0.025 {
        return 0.0;
    }
    sample_range_attenuation(radius, self_view_position.z - sample_view_z)
}

fn combine_opaque_and_transparent_occlusion(
    opaque: f32,
    transparent: f32,
    transparent_alpha: f32,
) -> f32 {
    opaque.max(transparent * transparent_alpha)
}

#[allow(clippy::too_many_arguments)]
pub(in crate::render) fn apply_occlusion_factors(
    framebuffer: &Framebuffer,
    color: &mut [[f32; 4]],
    factors: &OcclusionFactors,
    params: &OcclusionParams,
    camera: &CameraState,
    fog_intensity: f64,
    background: u32,
    transparent_background: bool,
    offset: (f32, f32),
) -> Result<(), String> {
    if color.len() != framebuffer.color.len() {
        return Err("occlusion color dimensions do not match the framebuffer".into());
    }
    if factors.values.len() != factors.width * factors.height {
        return Err("occlusion factor storage does not match its dimensions".into());
    }
    let occlusion_color = color_f32(parse_color(&params.color)?);
    let fog_color = color_f32(background);
    for (index, pixel) in color.iter_mut().enumerate() {
        if framebuffer.depth01[index] < 1.0 {
            let x = index % framebuffer.width;
            let y = index / framebuffer.width;
            let factor = sample_occlusion_factor(
                factors,
                (framebuffer.width, framebuffer.height),
                (x, y),
                offset,
            );
            let factor = valid_occlusion_factor(factor);
            let fog = fog_factor(camera, fog_intensity, framebuffer.depth01[index]);
            for (channel, occlusion) in occlusion_color.iter().enumerate() {
                let source = pixel[channel];
                let lower = if transparent_background {
                    occlusion * (1.0 - fog)
                } else {
                    occlusion * (1.0 - fog) + fog_color[channel] * fog
                };
                let mixed = lower * (1.0 - factor) + source * factor;
                pixel[channel] = mixed;
            }
        }
    }
    Ok(())
}

fn valid_occlusion_factor(raw: f32) -> f32 {
    // postprocessing.frag getSsao/getSsaoTransparent apply this only after
    // texture sampling and unpacking, not in the SSAO or blur attachments.
    if raw > 0.001 && raw <= 0.999 {
        raw
    } else {
        1.0
    }
}

fn sample_occlusion_factor(
    factors: &OcclusionFactors,
    (output_width, output_height): (usize, usize),
    (x, y): (usize, usize),
    (offset_x, offset_y): (f32, f32),
) -> f32 {
    let factor_width = factors.width;
    let factor_height = factors.height;
    if factor_width == output_width
        && factor_height == output_height
        && offset_x == 0.0
        && offset_y == 0.0
    {
        return factors.values[y * factor_width + x];
    }

    // `ssaoDepthTexture` is linearly sampled in normalized coordinates during
    // composition. Convert the full-resolution fragment center plus temporal
    // jitter to the lower-resolution texture's texel-center coordinate.
    let sample_x = ((x as f32 + 0.5 + offset_x) * factor_width as f32 / output_width as f32) - 0.5;
    let sample_y =
        ((y as f32 + 0.5 - offset_y) * factor_height as f32 / output_height as f32) - 0.5;
    let x0 = sample_x.floor() as isize;
    let y0 = sample_y.floor() as isize;
    let tx = sample_x - x0 as f32;
    let ty = sample_y - y0 as f32;
    let at = |sx: isize, sy: isize| {
        let sx = sx.clamp(0, factor_width.saturating_sub(1) as isize) as usize;
        let sy = sy.clamp(0, factor_height.saturating_sub(1) as isize) as usize;
        pack_unit_interval_rgba8(factors.values[sy * factor_width + sx])
    };
    let upper_left = at(x0, y0);
    let upper_right = at(x0 + 1, y0);
    let lower_left = at(x0, y0 + 1);
    let lower_right = at(x0 + 1, y0 + 1);
    let sample_channel = |channel: usize| {
        let upper_left = upper_left[channel] as f32;
        let lower_left = lower_left[channel] as f32;
        let upper = fused_multiply_add(upper_right[channel] as f32 - upper_left, tx, upper_left);
        let lower = fused_multiply_add(lower_right[channel] as f32 - lower_left, tx, lower_left);
        let byte_value = fused_multiply_add(lower - upper, ty, upper);
        // Apple GPU's normalized RGBA8 sampler rounds linear interpolation to
        // 1/16 of a byte before normalizing the channel to [0, 1]. Direct
        // Metal probes confirm positive half-way cases round upward.
        (byte_value * 16.0).round() / (16.0 * 255.0)
    };
    unpack_unit_interval_rg([sample_channel(0), sample_channel(1)])
}

fn bilateral_blur_occlusion(
    factors: &[f32],
    depth01: &[f32],
    width: usize,
    height: usize,
    camera: &CameraState,
    params: &OcclusionParams,
    bounds: [f32; 4],
) -> Vec<f32> {
    let kernel_size = params.blur_kernel_size.clamp(1, 25) | 1;
    // Like uRadius, uBlurDepthBias is multiplied by Camera.scale in
    // JavaScript-number precision before WebGL converts it to float32.
    let blur_depth_bias = (params.blur_depth_bias * camera.scale()) as f32;
    let half = kernel_size as isize / 2;
    let kernel = molstar_blur_kernel(kernel_size);
    let decoded_depth01 = depth01
        .iter()
        .map(|&value| packed_unit_interval_roundtrip(value))
        .collect::<Vec<_>>();
    let packed_depth = decoded_depth01
        .iter()
        .map(|&value| camera.depth_from_depth01(value))
        .collect::<Vec<_>>();
    let inverse_size = [divide(1.0, width as f32), divide(1.0, height as f32)];
    let pass = |source: &[f32], horizontal: bool| {
        let mut target = source.to_vec();
        for y in 0..height {
            for x in 0..width {
                let index = y * width + x;
                let coords = [
                    divide(x as f32 + 0.5, width as f32),
                    divide((height - y) as f32 - 0.5, height as f32),
                ];
                let self_depth = packed_depth[index];
                if outside_blur_bounds(coords, bounds)
                    || decoded_depth01[index] == 0.0
                    || (!horizontal && decoded_depth01[index] == 1.0)
                {
                    target[index] = ssao_background_factor();
                    continue;
                }
                let pixel_size =
                    blur_pixel_size(camera, coords, inverse_size[0], decoded_depth01[index]);
                let mut sum = 0.0;
                let mut kernel_sum = 0.0;
                for offset in -half..=half {
                    if offset.abs() > 1 && offset.abs() as f32 * pixel_size > 0.8 {
                        continue;
                    }
                    let sample_coords = if horizontal {
                        [
                            fused_multiply_add(offset as f32, inverse_size[0], coords[0]),
                            coords[1],
                        ]
                    } else {
                        [
                            coords[0],
                            fused_multiply_add(offset as f32, inverse_size[1], coords[1]),
                        ]
                    };
                    if outside_blur_bounds(sample_coords, bounds) {
                        continue;
                    }
                    // uBounds can extend beyond [0, 1] after SSAO target
                    // downscaling. Accepted samples still use CLAMP_TO_EDGE;
                    // they must not be discarded at the integer image edge.
                    let (sx, sy) = blur_sample_coordinates(x, y, offset, horizontal);
                    let sx = sx.clamp(0, width as isize - 1) as usize;
                    let sy = sy.clamp(0, height as isize - 1) as usize;
                    let sample_index = sy * width + sx;
                    let sample_depth = packed_depth[sample_index];
                    if decoded_depth01[sample_index] == 0.0
                        || decoded_depth01[sample_index] == 1.0
                        || (self_depth - sample_depth).abs() >= blur_depth_bias
                    {
                        continue;
                    }
                    let weight = kernel[offset.unsigned_abs()];
                    sum = accumulate_blur_sample(sum, source[sample_index], weight);
                    kernel_sum += weight;
                }
                // The shader divides even when every sample is rejected
                // (notably at zero depth bias). Its RG8 target stores the
                // resulting non-finite factor as zero, not the input factor.
                target[index] =
                    packed_unit_interval_roundtrip(divide_blur_sum(sum, kernel_sum, kernel_size));
            }
        }
        target
    };
    let horizontal = pass(factors, true);
    pass(&horizontal, false)
}

fn outside_blur_bounds(coords: [f32; 2], bounds: [f32; 4]) -> bool {
    coords[0] < bounds[0] || coords[1] < bounds[1] || coords[0] > bounds[2] || coords[1] > bounds[3]
}

fn blur_pixel_size(
    camera: &CameraState,
    coords: [f32; 2],
    inverse_width: f32,
    depth01: f32,
) -> f32 {
    let center = camera.screen_space_to_view_space(coords[0], coords[1], depth01);
    let right = camera.screen_space_to_view_space(coords[0] + inverse_width, coords[1], depth01);
    let delta = right - center;
    angle_metal_dot3(delta, delta).sqrt()
}

/// Port Mol*'s `getBlurKernel` JavaScript evaluation order and cast the
/// resulting Number values only at the WebGL float-uniform boundary.
fn molstar_blur_kernel(kernel_size: usize) -> Vec<f32> {
    let sigma = kernel_size as f64 / 3.0;
    let normalization = 1.0 / ((2.0 * std::f64::consts::PI).sqrt() * sigma);
    let half_kernel_size = kernel_size.div_ceil(2);
    (0..half_kernel_size)
        .map(|offset| {
            let x = offset as f64;
            (normalization * (-x * x / (2.0 * sigma * sigma)).exp()) as f32
        })
        .collect()
}

fn smootherstep(edge0: f32, edge1: f32, value: f32) -> f32 {
    let x = divide(value - edge0, edge1 - edge0).clamp(0.0, 1.0);
    let polynomial = fused_multiply_add(x, fused_multiply_add(x, 6.0, -15.0), 10.0);
    x * x * x * polynomial
}

fn ssao_noise_raw_vector(x: usize, y: usize, width: usize, height: usize) -> Vec3 {
    let [first, second] = ssao_noise_pair(x, y, width, height);
    Vec3::new(first * 2.0 - 1.0, second * 2.0 - 1.0, 0.0)
}

fn ssao_noise_pair(x: usize, y: usize, width: usize, height: usize) -> [f32; 2] {
    ssao_noise_coordinates(x, y, width, height).map(|(x, y)| ssao_noise(x, y))
}

#[allow(clippy::approx_constant)]
fn ssao_noise_coordinates(x: usize, y: usize, width: usize, height: usize) -> [(f32, f32); 2] {
    let inverse_width = (width as f32).recip();
    let inverse_height = (height as f32).recip();
    let fragment_x = x as f32 + 0.5;
    let fragment_y = (height - y) as f32 - 0.5;
    let coords_x = fragment_x * inverse_width;
    let coords_y = fragment_y * inverse_height;
    // Both hashes consume this rounded coordinate. Fusing the second hash's
    // offset with fragment scaling changes the shared-expression boundary;
    // a probe that emits only the second hash can be optimized differently.
    let second_x = coords_x + std::f32::consts::PI;
    let second_y = coords_y + 2.71828;
    [(coords_x, coords_y), (second_x, second_y)]
}

fn ssao_noise_phase(x: f32, y: f32) -> f32 {
    let dot = fused_multiply_add(y, 78.233, x * 12.9898);
    let quotient = (dot / std::f32::consts::PI).floor();
    fused_multiply_add(-std::f32::consts::PI, quotient, dot)
}

#[inline]
fn portable_sine(value: f32) -> f32 {
    let magnitude = value.abs();
    if !((1.0 / 4096.0)..=std::f32::consts::PI).contains(&magnitude) {
        // The GPU returns x below 2^-12; sinf also rounds to x there. Keep
        // libm's exceptional-value handling and general-domain fallback.
        return libm::sinf(value);
    }

    let (quadrant, residual) = sine_reduced_quadrant(magnitude);
    let folded = if quadrant == 1.0 {
        1.0 - residual.abs()
    } else {
        residual.abs()
    };
    if folded == 0.0 {
        return 0.0_f32.copysign(value);
    }

    // The reference evaluates a sinc-like special-function instruction and
    // then multiplies by its argument, with a binary32 rounding at each
    // boundary. A single correctly rounded sinf skips both this reduction
    // and the intermediate rounding. Retain the instruction's fixed-point
    // input boundary too; its remaining approximation error is documented
    // separately from these argument-reduction and multiplication steps.
    let sinc = sinc_approximation(folded);
    let result = sinc * folded;
    let sign = if quadrant == 2.0 && residual > 0.0 {
        -1.0
    } else {
        1.0
    };
    result * sign * value.signum()
}

fn sine_reduced_quadrant(magnitude: f32) -> (f32, f32) {
    // Split 2/pi as in the reference's compiled argument reduction. The
    // integer quadrant uses round-to-nearest-even, followed by three fused
    // operations. Retain their order instead of first combining the three
    // constants or rounding the converted angle before subtracting.
    let high = f32::from_bits(0x3f22_f983);
    let middle = f32::from_bits(0x32dc_9c88);
    let low = f32::from_bits(0x25a9_4fe1);
    let quadrant = fused_multiply_add(magnitude, high, 12_582_912.0) - 12_582_912.0;
    let residual = fused_multiply_add(magnitude, high, -quadrant);
    let residual = fused_multiply_add(magnitude, middle, residual);
    let residual = fused_multiply_add(magnitude, low, residual);
    (quadrant, residual)
}

#[allow(clippy::excessive_precision)]
fn ssao_noise(x: f32, y: f32) -> f32 {
    let sn = ssao_noise_phase(x, y);
    // Rust's primitive `sin` delegates to the host libm on native targets but
    // to compiler-builtins in WASM. Use one pure-Rust implementation so the
    // released renderer and the native parity harness consume identical SSAO
    // noise on every target, including the reference's argument-reduction
    // and special-function multiplication boundaries.
    let value = (portable_sine(sn) * 43_758.545_3).fract();
    if value < 0.0 {
        value + 1.0
    } else {
        value
    }
}

#[derive(Clone, Copy)]
struct Pcg32 {
    state: u32,
}

impl Pcg32 {
    fn new() -> Self {
        Self { state: 26_699 }
    }

    fn int(&mut self) -> u32 {
        let old_state = self.state;
        self.state = self
            .state
            .wrapping_mul(1_664_525)
            .wrapping_add(1_013_904_223);
        let xor_shifted = ((old_state >> 18) ^ old_state) >> 5;
        xor_shifted.rotate_right(old_state >> 27)
    }

    fn float(&mut self) -> f64 {
        self.int() as f64 / 4_294_967_296.0
    }
}

fn random_hemisphere_vector(pcg: &mut Pcg32) -> [f64; 3] {
    loop {
        let x = pcg.float() * 2.0 - 1.0;
        let y = pcg.float() * 2.0 - 1.0;
        let xy = x * x + y * y;
        if xy < 1.0 {
            let sign = if pcg.float() < 0.5 { -1.0 } else { 1.0 };
            let z = 2.0 * (1.0 - xy).sqrt() * sign;
            let inverse_length = 1.0 / (x * x + y * y + z * z).sqrt();
            let scale = pcg.float();
            let mut vector = [
                x * inverse_length * scale,
                y * inverse_length * scale,
                z * inverse_length * scale,
            ];
            if vector[2] < 0.0 {
                vector[2] = -vector[2];
            }
            return vector;
        }
    }
}

fn extend_blue_noise_vectors(
    blue_noise: &mut Vec<[f64; 3]>,
    pcg: &mut Pcg32,
    count: usize,
    candidate_count: usize,
) {
    if blue_noise.is_empty() && count > 0 {
        blue_noise.push(random_hemisphere_vector(pcg));
    }
    while blue_noise.len() < count {
        let mut best = [0.0; 3];
        let mut best_distance = -1.0f64;
        for _ in 0..candidate_count {
            let candidate = random_hemisphere_vector(pcg);
            let min_distance = blue_noise
                .iter()
                .map(|existing| {
                    let dx = candidate[0] - existing[0];
                    let dy = candidate[1] - existing[1];
                    let dz = candidate[2] - existing[2];
                    (dx * dx + dy * dy + dz * dz).sqrt()
                })
                .fold(f64::INFINITY, f64::min);
            if min_distance > best_distance {
                best_distance = min_distance;
                best = candidate;
            }
        }
        blue_noise.push(best);
    }
}

fn molstar_ssao_samples(count: usize) -> Vec<Vec3> {
    let mut pcg = Pcg32::new();
    // Creating Mol*'s SSAO renderable eagerly calls `getSamples(32)`. A later
    // non-default sample count reuses that module-level 32-vector prefix and
    // only extends the blue-noise cache. Reproduce that observable staging;
    // generating a large count from an empty cache chooses a different set as
    // soon as candidateCount exceeds ten.
    let initial_count = 32;
    let mut blue_noise = Vec::with_capacity(count.max(initial_count));
    extend_blue_noise_vectors(&mut blue_noise, &mut pcg, initial_count, 10);
    if count > initial_count {
        let candidate_count = (count / 10).clamp(10, 30);
        extend_blue_noise_vectors(&mut blue_noise, &mut pcg, count, candidate_count);
    }
    blue_noise
        .into_iter()
        .take(count)
        .enumerate()
        .map(|(index, sample)| {
            let n = count as f64;
            let i = index as f64;
            let scale = 0.1 + ((i * i + 2.0 * i + 1.0) / (n * n)) * 0.9;
            Vec3::new(
                (sample[0] * scale) as f32,
                (sample[1] * scale) as f32,
                (sample[2] * scale) as f32,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {
    #![allow(clippy::approx_constant, clippy::excessive_precision)]

    use super::*;
    use crate::render::postprocessing::packing::unpack_unit_interval_rgba8;

    #[test]
    fn compose_rejects_occlusion_texture_values_only_after_sampling() {
        for raw in [0.0, 0.001, f32::from_bits(0.999f32.to_bits() + 1), 1.0] {
            assert_eq!(valid_occlusion_factor(raw), 1.0);
        }
        for raw in [f32::from_bits(0.001f32.to_bits() + 1), 0.5, 0.999] {
            assert_eq!(valid_occlusion_factor(raw).to_bits(), raw.to_bits());
        }
    }

    #[test]
    fn background_factor_uses_the_rg8_render_target_roundtrip() {
        assert_eq!(ssao_background_factor(), 255.0 / 256.0);
    }

    fn blur_test_camera() -> CameraState {
        let renderer = crate::render::options::RendererOptions::default();
        let sphere = crate::model::BoundingSphere {
            radius: 1.0,
            ..crate::model::BoundingSphere::default()
        };
        crate::render::camera::resolve_camera(&renderer, &sphere, 64, 64).unwrap()
    }

    #[test]
    fn blur_near_clip_writes_packed_no_occlusion() {
        let factors = [0.25; 9];
        let mut depth = [0.5; 9];
        depth[4] = 0.0;
        let params = OcclusionParams {
            blur_kernel_size: 3,
            ..OcclusionParams::default()
        };
        let result = bilateral_blur_occlusion(
            &factors,
            &depth,
            3,
            3,
            &blur_test_camera(),
            &params,
            [0.0, 0.0, 1.0, 1.0],
        );
        assert_eq!(result[4], ssao_background_factor());
        // The near-clipped sample is excluded from both neighbors' sums.
        assert_eq!(result[3], packed_unit_interval_roundtrip(0.25));
        assert_eq!(result[5], packed_unit_interval_roundtrip(0.25));
    }

    #[test]
    fn blur_zero_depth_bias_does_not_preserve_input_factors() {
        let mut depth = [0.5; 9];
        depth[4] = 0.0;
        let params = OcclusionParams {
            blur_depth_bias: 0.0,
            ..OcclusionParams::default()
        };
        let result = bilateral_blur_occlusion(
            &[0.25; 9],
            &depth,
            3,
            3,
            &blur_test_camera(),
            &params,
            [0.0, 0.0, 1.0, 1.0],
        );
        for (index, value) in result.into_iter().enumerate() {
            assert_eq!(
                value,
                if index == 4 {
                    ssao_background_factor()
                } else {
                    0.0
                }
            );
        }
    }

    #[test]
    fn blur_bounds_are_inclusive_normalized_coordinates() {
        let bounds = [0.25, 0.25, 0.75, 0.75];
        for point in [[0.25, 0.25], [0.75, 0.75], [0.5, 0.5]] {
            assert!(!outside_blur_bounds(point, bounds));
        }
        for point in [[0.0, 0.5], [1.0, 0.5], [0.5, 0.0], [0.5, 1.0]] {
            assert!(outside_blur_bounds(point, bounds));
        }
        let params = OcclusionParams {
            blur_kernel_size: 1,
            ..OcclusionParams::default()
        };
        let result = bilateral_blur_occlusion(
            &[0.5; 9],
            &[0.5; 9],
            3,
            3,
            &blur_test_camera(),
            &params,
            bounds,
        );
        for (index, value) in result.into_iter().enumerate() {
            assert_eq!(
                value,
                if index == 4 {
                    0.5
                } else {
                    ssao_background_factor()
                }
            );
        }
    }

    #[test]
    fn blur_clamps_edge_samples_accepted_by_expanded_bounds() {
        let factors = [0.25, 0.25, 0.25, 0.5, 0.5, 0.5, 0.75, 0.75, 0.75];
        let camera = blur_test_camera();
        let params = OcclusionParams {
            blur_kernel_size: 3,
            ..OcclusionParams::default()
        };
        let clipped = bilateral_blur_occlusion(
            &factors,
            &[0.5; 9],
            3,
            3,
            &camera,
            &params,
            [0.0, 0.0, 1.0, 1.0],
        );
        let extended = bilateral_blur_occlusion(
            &factors,
            &[0.5; 9],
            3,
            3,
            &camera,
            &params,
            [0.0, 0.0, 1.0, 1.2],
        );
        // GLSL's positive Y points to the top row. The extra sample repeats
        // that row, so its lower factor gains weight only at the upper edge.
        assert!(extended[1] < clipped[1]);
        assert_eq!(extended[4], clipped[4]);
        assert_eq!(extended[7], clipped[7]);
    }

    #[test]
    fn vertical_blur_preserves_glsl_offset_iteration_order() {
        assert_eq!(blur_sample_coordinates(7, 11, -1, false), (7, 12));
        assert_eq!(blur_sample_coordinates(7, 11, 0, false), (7, 11));
        assert_eq!(blur_sample_coordinates(7, 11, 1, false), (7, 10));
        assert_eq!(blur_sample_coordinates(7, 11, -1, true), (6, 11));
        assert_eq!(blur_sample_coordinates(7, 11, 1, true), (8, 11));
    }

    #[test]
    fn blur_accumulation_uses_one_fused_rounding() {
        let sample = 0.074_821_845f32;
        let weight = 0.227_215_02f32;
        let sum = 0.037_286_54f32;
        assert_eq!(
            accumulate_blur_sample(sum, sample, weight).to_bits(),
            0x3d5e_5c3d
        );
        assert_eq!((sample * weight + sum).to_bits(), 0x3d5e_5c3e);
    }

    #[test]
    fn every_supported_blur_kernel_matches_molstar_number_staging() {
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        for kernel_size in (1..=25).step_by(2) {
            let kernel = molstar_blur_kernel(kernel_size);
            assert_eq!(kernel.len(), kernel_size.div_ceil(2));
            for byte in kernel.into_iter().flat_map(f32::to_le_bytes) {
                hash = (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        // Generated by the pinned TypeScript `getBlurKernel` implementation,
        // after Float32Array/WebGL uniform conversion.
        assert_eq!(hash, 0x58f8_5c42_c857_e68f);
    }

    #[test]
    fn noise_pair_reuses_the_rounded_coordinate_before_adding_its_offset() {
        let [first, second] = ssao_noise_coordinates(0, 4, 1024, 937);
        assert_eq!(
            second.0.to_bits(),
            (first.0 + std::f32::consts::PI).to_bits()
        );
        assert_eq!(second.1.to_bits(), (first.1 + 2.71828).to_bits());
        assert_ne!(
            fused_multiply_add(937.0 - 4.0 - 0.5, (937.0f32).recip(), 2.71828).to_bits(),
            second.1.to_bits(),
        );
    }

    #[test]
    fn portable_sine_has_stable_render_domain_bits() {
        for (input, expected) in [
            (0x4033_28d8, 0x3eab_d19f),
            (0x3fe6_2c5a, 0x3f79_685a),
            (0x4033_f8b8, 0x3ea5_af8f),
            (0x3fe7_cc5a, 0x3f78_a79e),
            (0x3ff3_2a5a, 0x3f72_4671),
            (0x401b_ae2d, 0x3f26_b12d),
            (0x3fff_a6ee, 0x3f69_11a1),
            (0x4040_b5a2, 0x3e05_4282),
            (0x4026_5aa4, 0x3f04_2032),
        ] {
            assert_eq!(portable_sine(f32::from_bits(input)).to_bits(), expected);
        }
    }

    #[test]
    fn sine_argument_reduction_retains_split_constant_fma_boundaries() {
        for (input, quadrant, residual) in [
            (0x39a4_b5bd, 0x0000_0000, 0x3951_b716),
            (0x3ec9_0fda, 0x0000_0000, 0x3e80_0000),
            (0x3f49_0fda, 0x0000_0000, 0x3f00_0000),
            (0x3f86_0a88, 0x3f80_0000, 0xbeaa_aac3),
            (0x3fc9_0fda, 0x3f80_0000, 0xb34e_6e54),
            (0x4006_0a88, 0x3f80_0000, 0x3eaa_aa79),
            (0x4016_cbe4, 0x3f80_0000, 0x3f00_0000),
            (0x4049_0fcc, 0x4000_0000, 0xb615_0dc6),
            (0x4049_0fda, 0x4000_0000, 0xb3ce_6e54),
        ] {
            let actual = sine_reduced_quadrant(f32::from_bits(input));
            assert_eq!(
                (actual.0.to_bits(), actual.1.to_bits()),
                (quadrant, residual)
            );
        }
    }

    #[test]
    fn portable_sine_keeps_tiny_inputs_and_general_domain_fallback() {
        for bits in [0, 1, 0x007f_ffff, 0x0080_0000, 0x397f_ffff] {
            let value = f32::from_bits(bits);
            assert_eq!(portable_sine(value).to_bits(), bits);
            assert_eq!(portable_sine(-value).to_bits(), bits | 0x8000_0000);
        }
        for value in [4.0, -4.0, 1000.0, -1000.0, f32::MAX, f32::MIN] {
            assert_eq!(portable_sine(value).to_bits(), libm::sinf(value).to_bits());
        }
        for value in [f32::NAN, f32::INFINITY, f32::NEG_INFINITY] {
            assert!(portable_sine(value).is_nan());
        }
    }

    #[test]
    fn portable_sine_is_odd_and_bounded_across_noise_phase_domain() {
        for index in 1..=65_536 {
            let value = index as f32 * (std::f32::consts::PI / 65_536.0);
            let positive = portable_sine(value);
            let negative = portable_sine(-value);
            assert_eq!(negative.to_bits(), (-positive).to_bits());
            assert!(positive.abs() <= 1.0);
            assert!((positive - libm::sinf(value)).abs() <= 3.0e-7);
        }
        assert!(portable_sine(std::f32::consts::PI) < 0.0);
    }

    #[test]
    fn arbitrary_sample_counts_reuse_molstar_eager_32_vector_prefix() {
        fn hash_samples(samples: &[Vec3]) -> u64 {
            let mut hash = 0xcbf2_9ce4_8422_2325u64;
            for byte in samples.iter().flat_map(|sample| {
                [sample.x, sample.y, sample.z]
                    .into_iter()
                    .flat_map(f32::to_le_bytes)
            }) {
                hash = (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3);
            }
            hash
        }

        // These hashes come from the pinned Mol* TypeScript implementation
        // after its renderable constructor has evaluated getSamples(32).
        assert_eq!(
            hash_samples(&molstar_ssao_samples(32)),
            0x8e44_4a14_2762_0fab
        );
        assert_eq!(
            hash_samples(&molstar_ssao_samples(128)),
            0x8e40_677e_f504_b647
        );
        assert_eq!(
            hash_samples(&molstar_ssao_samples(256)),
            0xf8cb_fd9d_a430_78fd
        );
    }

    #[test]
    fn fractional_radius_is_rounded_only_at_the_uniform_boundary() {
        let params = OcclusionParams {
            radius: 5.123_456_789,
            ..OcclusionParams::default()
        };
        let renderer = super::super::super::options::RendererOptions::default();
        let sphere = crate::model::BoundingSphere {
            center: Vec3::default(),
            radius: 1.0,
            ..crate::model::BoundingSphere::default()
        };
        let camera = crate::render::camera::resolve_camera(&renderer, &sphere, 64, 64).unwrap();
        let staged = ssao_uniform_radius(&params, &camera);
        let prematurely_rounded = 2.0f32.powf(params.radius as f32);
        assert_eq!(staged.to_bits(), 0x420b_6f8e);
        assert_ne!(staged.to_bits(), prematurely_rounded.to_bits());
    }

    #[test]
    fn default_radius_matches_the_pinned_molstar_uniform_capture() {
        let params = OcclusionParams::default();
        let renderer = super::super::super::options::RendererOptions::default();
        let sphere = crate::model::BoundingSphere {
            center: Vec3::default(),
            radius: 1.0,
            ..crate::model::BoundingSphere::default()
        };
        let camera = crate::render::camera::resolve_camera(&renderer, &sphere, 64, 64).unwrap();

        // Captured from `uRadius` in the pinned Mol* PDB SSAO pass. Ordinary
        // static Canvas3D rendering keeps Camera.scale at its default 1.0.
        assert_eq!(
            ssao_uniform_radius(&params, &camera).to_bits(),
            32.0f32.to_bits()
        );
    }

    #[test]
    fn range_attenuation_keeps_the_shader_zero_delta_semantics() {
        assert_eq!(
            sample_range_attenuation(32.0, 0.0).to_bits(),
            1.0f32.to_bits()
        );
        assert_eq!(
            sample_range_attenuation(1.0, 2.0).to_bits(),
            0.5f32.to_bits()
        );
        assert_eq!(
            sample_range_attenuation(1.0, f32::INFINITY).to_bits(),
            0.0f32.to_bits()
        );
    }

    #[test]
    fn apple_rg8_linear_filter_matches_metal_fractional_byte_quantization() {
        let factors = OcclusionFactors {
            width: 2,
            height: 2,
            values: [[0, 255], [0, 255], [0, 255], [255, 254]]
                .map(unpack_unit_interval_rgba8)
                .to_vec(),
        };
        assert_eq!(
            sample_occlusion_factor(&factors, (2, 2), (0, 0), (0.375, -0.125)).to_bits(),
            0.996_031_76f32.to_bits(),
        );
    }

    #[test]
    fn ssao_target_size_combines_resolution_scale_and_high_dpi_downscaling() {
        assert_eq!(
            ssao_target_dimensions(120, 90, 1.5, 1.0),
            (80, 60, 2.0 / 3.0)
        );
        assert_eq!(ssao_target_dimensions(800, 600, 2.0, 0.5), (200, 150, 0.25));
        assert_eq!(ssao_target_dimensions(1, 1, 4.0, 0.1), (1, 1, 0.025));
    }

    #[test]
    fn full_viewport_bounds_retain_the_unfloored_scaled_extent() {
        let bounds = ssao_full_viewport_bounds(7, 5, 0.6);
        assert_eq!(bounds.map(f32::to_bits), [0, 0, 0x3f98_6186, 0x3f80_0000]);
    }

    #[test]
    fn depth_downsampling_uses_normalized_nearest_texture_coordinates() {
        let source = (0..16).map(|value| value as f32).collect::<Vec<_>>();
        assert_eq!(
            nearest_scaled_values(&source, 4, 4, 2, 2),
            vec![5.0, 7.0, 13.0, 15.0]
        );
        assert_eq!(nearest_scaled_values(&source, 4, 4, 4, 4), source);
    }

    #[test]
    fn multiscale_depth_pyramid_follows_nearest_then_linear_copy_targets() {
        let source = (0..8)
            .flat_map(|y| (0..8).map(move |x| (y * 10 + x) as f32))
            .collect::<Vec<_>>();
        let pyramid = DepthPyramid::new(&source, 8, 8, false);
        assert_eq!(pyramid.half_width, 4);
        assert_eq!(pyramid.half_height, 4);
        assert_eq!(
            pyramid.half,
            vec![
                1.0, 3.0, 5.0, 7.0, 21.0, 23.0, 25.0, 27.0, 41.0, 43.0, 45.0, 47.0, 61.0, 63.0,
                65.0, 67.0
            ]
        );
        assert_eq!(pyramid.quarter_width, 2);
        assert_eq!(pyramid.quarter_height, 2);
        assert_eq!(pyramid.quarter, vec![12.0, 16.0, 52.0, 56.0]);
    }

    #[test]
    fn multiscale_depth_mapping_uses_strict_half_and_quarter_thresholds() {
        let full = [0.1];
        let pyramid = DepthPyramid {
            full: &full,
            full_width: 1,
            full_height: 1,
            full_linear: false,
            half: vec![0.2],
            half_width: 1,
            half_height: 1,
            quarter: vec![0.3],
            quarter_width: 1,
            quarter_height: 1,
        };
        let self_offset = [0.0, 0.0];
        assert_eq!(
            pyramid.sample_mapped([0.05, 0.0], self_offset).to_bits(),
            0.1f32.to_bits()
        );
        assert_eq!(
            pyramid.sample_mapped([0.0501, 0.0], self_offset).to_bits(),
            0.2f32.to_bits()
        );
        assert_eq!(
            pyramid.sample_mapped([0.1, 0.0], self_offset).to_bits(),
            0.2f32.to_bits()
        );
        assert_eq!(
            pyramid.sample_mapped([0.1001, 0.0], self_offset).to_bits(),
            0.3f32.to_bits()
        );
    }

    #[test]
    fn transparent_depth_filters_packed_rgba8_before_decoding() {
        let depth = [0.123_456_7, 0.654_321];
        let alpha = [0.3, 0.6];
        let packed = depth
            .into_iter()
            .zip(alpha)
            .map(|(depth, alpha)| pack_depth_alpha_rgba8(depth, alpha))
            .collect::<Vec<_>>();
        let pyramid = PackedDepthAlphaPyramid::new(&packed, 2, 1, true, false);
        let (filtered_depth, filtered_alpha) = pyramid.sample_full([0.5, 0.5]);
        assert_eq!(filtered_depth.to_bits(), 0x3ec7_1c71);
        assert_eq!(filtered_alpha.to_bits(), 0x3ee6_e6e7);
        assert_ne!(
            filtered_alpha.to_bits(),
            ((alpha[0] + alpha[1]) * 0.5).to_bits()
        );
    }

    #[test]
    fn transparent_depth_filters_the_actual_white_clear_target() {
        let fragment = pack_depth_alpha_rgba8(0.25, 0.5);
        let actual = PackedDepthAlphaPyramid::new(&[fragment, [255; 4]], 2, 1, true, false);
        let legacy_background = pack_depth_alpha_rgba8(1.0, 0.0);
        let legacy =
            PackedDepthAlphaPyramid::new(&[fragment, legacy_background], 2, 1, true, false);

        let actual_filtered = actual.sample_full([0.5, 0.5]);
        let legacy_filtered = legacy.sample_full([0.5, 0.5]);
        assert_eq!(actual_filtered.0.to_bits(), 0x3f1f_ffff);
        assert_eq!(actual_filtered.1.to_bits(), 0x3f40_4040);
        assert_ne!(actual_filtered.0.to_bits(), legacy_filtered.0.to_bits());
        assert_ne!(actual_filtered.1.to_bits(), legacy_filtered.1.to_bits());
    }

    #[test]
    fn multiscale_occlusion_executes_every_level_on_a_finite_depth_target() {
        let mut framebuffer = Framebuffer::new(8, 8, [255, 255, 255, 255]);
        framebuffer.depth01.fill(0.5);
        let renderer = super::super::super::options::RendererOptions::default();
        let sphere = crate::model::BoundingSphere {
            center: Vec3::default(),
            radius: 1.0,
            ..crate::model::BoundingSphere::default()
        };
        let camera = crate::render::camera::resolve_camera(&renderer, &sphere, 8, 8).unwrap();
        let params = OcclusionParams {
            samples: 4,
            multi_scale: super::super::super::options::OcclusionMultiScale {
                name: "on".into(),
                ..super::super::super::options::OcclusionMultiScale::default()
            },
            blur_kernel_size: 1,
            ..OcclusionParams::default()
        };

        let factors = compute_occlusion_factors(&framebuffer, &camera, &params, 1.0);
        assert_eq!((factors.width, factors.height), (8, 8));
        assert!(factors.values.iter().all(|value| value.is_finite()));
    }

    #[test]
    fn lower_resolution_occlusion_is_sampled_at_full_resolution_fragment_centers() {
        let left = unpack_unit_interval_rgba8([0, 100]);
        let right = unpack_unit_interval_rgba8([0, 200]);
        let factors = OcclusionFactors {
            width: 2,
            height: 1,
            values: vec![left, right],
        };
        assert_eq!(
            sample_occlusion_factor(&factors, (4, 1), (0, 0), (0.0, 0.0)).to_bits(),
            left.to_bits()
        );
        assert_eq!(
            sample_occlusion_factor(&factors, (4, 1), (3, 0), (0.0, 0.0)).to_bits(),
            right.to_bits()
        );
        let inner = sample_occlusion_factor(&factors, (4, 1), (1, 0), (0.0, 0.0));
        assert!(inner > left && inner < right);
    }

    #[test]
    fn transparent_ssao_uses_alpha_weighted_max_from_molstar_shader() {
        assert_eq!(
            combine_opaque_and_transparent_occlusion(0.2, 0.9, 0.5).to_bits(),
            0.45f32.to_bits()
        );
        assert_eq!(
            combine_opaque_and_transparent_occlusion(0.8, 0.9, 0.5).to_bits(),
            0.8f32.to_bits()
        );
        assert_eq!(
            combine_opaque_and_transparent_occlusion(0.2, 0.9, 0.0).to_bits(),
            0.2f32.to_bits()
        );
    }
}
