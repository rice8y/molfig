use super::super::camera::CameraState;
use super::super::color::{color_f32, parse_color, quantize};
use super::super::framebuffer::Framebuffer;
use super::super::options::OutlineParams;
use super::fog::fog_factor;
use super::packing::packed_unit_interval_roundtrip;

#[allow(clippy::too_many_arguments)]
pub(in crate::render) fn apply_outline(
    framebuffer: &Framebuffer,
    color: &mut [[f32; 4]],
    camera: &CameraState,
    params: &OutlineParams,
    pixel_ratio: f64,
    fog_intensity: f64,
    background: u32,
    transparent_background: bool,
) -> Result<Vec<f32>, String> {
    if color.len() != framebuffer.color.len() {
        return Err("outline color dimensions do not match the framebuffer".into());
    }
    let outline_color = parse_color(&params.color)?;
    let outline_rgb = color_f32(outline_color);
    let fog_rgb = color_f32(background);
    let width = framebuffer.width;
    let height = framebuffer.height;
    let mut outline_depth01 = vec![1.0; width * height];
    for y in 0..height {
        for x in 0..width {
            let index = y * width + x;
            let center_depth01 = framebuffer.depth01[index];
            let center_depth = camera.depth_from_depth01(center_depth01);
            let center_view_z = if center_depth01 < 1.0 {
                -center_depth
            } else {
                2.0 * camera.far
            };
            // OutlinePass.update multiplies the public threshold by 50 and by
            // the WebGL pixel ratio before passing uOutlineThreshold.
            let coords_x = (x as f32 + 0.5) / width as f32;
            let coords_y = (height as f32 - y as f32 - 0.5) / height as f32;
            let view_position =
                camera.screen_space_to_view_space(coords_x, coords_y, center_depth01);
            let view_position_right = camera.screen_space_to_view_space(
                coords_x + 1.0 / width as f32,
                coords_y,
                center_depth01,
            );
            let threshold = (view_position_right - view_position).length()
                * outline_threshold(params, pixel_ratio);
            let mut best_depth01 = 1.0;
            for oy in -1isize..=1 {
                for ox in -1isize..=1 {
                    let nx = clamp_texture_coord(x as isize + ox, width);
                    // GLSL iterates bottom-up texture Y while the CPU target
                    // is stored top-down. Reverse only the storage mapping so
                    // equal-depth tie updates retain shader loop order.
                    let ny = clamp_texture_coord(y as isize - oy, height);
                    let other_depth01 = framebuffer.depth01[ny * width + nx];
                    let other = camera.depth_from_depth01(other_depth01);
                    let other_view_z = if other_depth01 < 1.0 {
                        -other
                    } else {
                        2.0 * camera.far
                    };
                    if (center_view_z - other_view_z).abs() > threshold
                        && center_depth01 > other_depth01
                        && other_depth01 <= best_depth01
                    {
                        best_depth01 = other_depth01;
                    }
                }
            }
            let mut outlined = best_depth01 < 1.0;
            if outlined && center_depth01 < 1.0 {
                let depth_at = |sx: isize, sy: isize| -> f32 {
                    let sx = clamp_texture_coord(sx, width);
                    let sy = clamp_texture_coord(sy, height);
                    let depth01 = framebuffer.depth01[sy * width + sx];
                    if depth01 < 1.0 {
                        -camera.depth_from_depth01(depth01)
                    } else {
                        2.0 * camera.far
                    }
                };
                let left = depth_at(x as isize - 1, y as isize);
                let right = depth_at(x as isize + 1, y as isize);
                let up = depth_at(x as isize, y as isize - 1);
                let down = depth_at(x as isize, y as isize + 1);
                let curvature = (left + right - 2.0 * center_view_z)
                    .abs()
                    .max((up + down - 2.0 * center_view_z).abs());
                if curvature < threshold * 0.75 {
                    outlined = false;
                }
            }
            if outlined {
                outline_depth01[index] = packed_unit_interval_roundtrip(best_depth01);
            }
        }
    }
    let dilation = outline_scale(params, pixel_ratio);
    if dilation > 0 {
        let original = outline_depth01.clone();
        for y in 0..height {
            for x in 0..width {
                let mut best_depth01: f32 = 1.0;
                for oy in -dilation..=dilation {
                    for ox in -dilation..=dilation {
                        if ox * ox + oy * oy > dilation * dilation {
                            continue;
                        }
                        let nx = clamp_texture_coord(x as isize + ox, width);
                        let ny = clamp_texture_coord(y as isize - oy, height);
                        best_depth01 = best_depth01.min(original[ny * width + nx]);
                    }
                }
                outline_depth01[y * width + x] = best_depth01;
            }
        }
    }
    for (index, depth01) in outline_depth01.iter().copied().enumerate() {
        if depth01 < 1.0 {
            let fog = fog_factor(camera, fog_intensity, depth01);
            let pixel = &mut color[index];
            if transparent_background {
                let alpha = 1.0 - fog;
                for channel in 0..3 {
                    pixel[channel] = outline_rgb[channel] * alpha;
                }
                pixel[3] = alpha;
            } else {
                for channel in 0..3 {
                    pixel[channel] = outline_rgb[channel] * (1.0 - fog) + fog_rgb[channel] * fog;
                }
            }
        }
    }
    Ok(outline_depth01)
}

#[allow(clippy::too_many_arguments)]
pub(in crate::render) fn apply_transparent_outline(
    opaque: &Framebuffer,
    transparent: &Framebuffer,
    color: &mut [[f32; 4]],
    transparent_alpha: &[f32],
    opaque_outline_depth01: &[f32],
    camera: &CameraState,
    params: &OutlineParams,
    pixel_ratio: f64,
    fog_intensity: f64,
) -> Result<(), String> {
    if !params.include_transparent {
        return Ok(());
    }
    let len = opaque.width * opaque.height;
    if transparent.depth01.len() != len
        || color.len() != len
        || transparent_alpha.len() != len
        || opaque_outline_depth01.len() != len
    {
        return Err("transparent outline dimensions do not match the framebuffer".into());
    }
    let outline_rgb = color_f32(parse_color(&params.color)?);
    let width = opaque.width;
    let height = opaque.height;
    let view_z = |depth01: f32| {
        if depth01 < 1.0 {
            -camera.depth_from_depth01(depth01)
        } else {
            2.0 * camera.far
        }
    };
    let mut outline_depth01 = vec![1.0; len];
    let mut outline_alpha = vec![0.0; len];
    for y in 0..height {
        for x in 0..width {
            let index = y * width + x;
            let self_depth = transparent.depth01[index];
            let self_view_z = view_z(self_depth);
            let coords_x = (x as f32 + 0.5) / width as f32;
            let coords_y = (height as f32 - y as f32 - 0.5) / height as f32;
            let view_position = camera.screen_space_to_view_space(coords_x, coords_y, self_depth);
            let view_position_right = camera.screen_space_to_view_space(
                coords_x + 1.0 / width as f32,
                coords_y,
                self_depth,
            );
            let threshold = (view_position_right - view_position).length()
                * outline_threshold(params, pixel_ratio);
            let mut best_depth = 1.0;
            let mut best_alpha = 0.0;
            for oy in -1isize..=1 {
                for ox in -1isize..=1 {
                    let nx = clamp_texture_coord(x as isize + ox, width);
                    let ny = clamp_texture_coord(y as isize - oy, height);
                    let sample_index = ny * width + nx;
                    let sample_depth = transparent.depth01[sample_index];
                    if (self_view_z - view_z(sample_depth)).abs() > threshold
                        && self_depth > sample_depth
                        && sample_depth <= best_depth
                    {
                        best_depth = sample_depth;
                        best_alpha = transparent_alpha[sample_index];
                    }
                }
            }
            if best_depth >= 1.0 {
                continue;
            }

            let opaque_depth = opaque.depth01[index];
            let opaque_coords = camera.screen_space_to_view_space(coords_x, coords_y, opaque_depth);
            let opaque_right = camera.screen_space_to_view_space(
                coords_x + 1.0 / width as f32,
                coords_y,
                opaque_depth,
            );
            let opaque_pixel_size =
                (opaque_right - opaque_coords).length() * outline_threshold(params, pixel_ratio);
            if (view_z(opaque_depth) - view_z(best_depth)).abs() < opaque_pixel_size {
                continue;
            }

            if self_depth < 1.0 {
                let depth_at = |sx: isize, sy: isize| -> f32 {
                    let sx = clamp_texture_coord(sx, width);
                    let sy = clamp_texture_coord(sy, height);
                    view_z(transparent.depth01[sy * width + sx])
                };
                let curvature = (depth_at(x as isize - 1, y as isize)
                    + depth_at(x as isize + 1, y as isize)
                    - 2.0 * self_view_z)
                    .abs()
                    .max(
                        (depth_at(x as isize, y as isize - 1)
                            + depth_at(x as isize, y as isize + 1)
                            - 2.0 * self_view_z)
                            .abs(),
                    );
                if curvature < threshold * 0.75 {
                    continue;
                }
            }
            // tOutlines stores transparent depth directly in its 8-bit alpha
            // channel and packs transparent alpha through pack2x4. Preserve
            // both intermediate render-target quantizations before dilation.
            outline_depth01[index] = f32::from(quantize(best_depth)) / 255.0;
            outline_alpha[index] = packed_outline_alpha(best_alpha);
        }
    }

    let dilation = outline_scale(params, pixel_ratio);
    if dilation > 0 {
        let original_depth = outline_depth01.clone();
        let original_alpha = outline_alpha.clone();
        for y in 0..height {
            for x in 0..width {
                let index = y * width + x;
                let mut best_depth = 1.0;
                let mut best_alpha = 0.0;
                for oy in -dilation..=dilation {
                    for ox in -dilation..=dilation {
                        if ox * ox + oy * oy > dilation * dilation {
                            continue;
                        }
                        let nx = clamp_texture_coord(x as isize + ox, width);
                        let ny = clamp_texture_coord(y as isize - oy, height);
                        let sample_index = ny * width + nx;
                        if original_depth[sample_index] < best_depth {
                            best_depth = original_depth[sample_index];
                            best_alpha = original_alpha[sample_index];
                        }
                    }
                }
                outline_depth01[index] = best_depth;
                outline_alpha[index] = best_alpha;
            }
        }
    }

    for index in 0..len {
        let depth = outline_depth01[index];
        if depth >= 1.0 {
            continue;
        }
        let opaque_outline_depth = opaque_outline_depth01[index];
        if opaque_outline_depth < 1.0 && opaque_outline_depth < depth {
            color[index] = [0.0; 4];
            continue;
        }
        let outline_alpha = outline_alpha[index];
        let fog = fog_factor(camera, fog_intensity, depth);
        let current_alpha = color[index][3];
        let final_alpha = current_alpha.max(outline_alpha * (1.0 - fog));
        for (channel, value) in outline_rgb.iter().enumerate() {
            color[index][channel] = *value * final_alpha;
        }
        color[index][3] = final_alpha;
    }
    Ok(())
}

fn outline_scale(params: &OutlineParams, pixel_ratio: f64) -> isize {
    ((params.scale * pixel_ratio).round().max(1.0) - 1.0) as isize
}

fn outline_threshold(params: &OutlineParams, pixel_ratio: f64) -> f32 {
    (50.0 * params.threshold * pixel_ratio) as f32
}

fn clamp_texture_coord(value: isize, size: usize) -> usize {
    value.clamp(0, size.saturating_sub(1) as isize) as usize
}

fn packed_outline_alpha(alpha: f32) -> f32 {
    let packed_high_nibble = ((alpha.clamp(0.0, 0.5) * 2.0 * 15.0) + 0.5).floor();
    packed_high_nibble / 15.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::render::options::RendererOptions;
    use crate::render::postprocessing::packed_depth_alpha_roundtrip;

    #[test]
    fn scale_and_threshold_include_webgl_pixel_ratio() {
        let params = OutlineParams {
            scale: 2.0,
            threshold: 0.33,
            ..OutlineParams::default()
        };
        assert_eq!(outline_scale(&params, 1.5), 2);
        assert_eq!(
            outline_threshold(&params, 1.5).to_bits(),
            24.75f32.to_bits()
        );
    }

    #[test]
    fn transparent_outline_alpha_matches_pack2x4_roundtrip() {
        assert_eq!(packed_outline_alpha(0.0).to_bits(), 0.0f32.to_bits());
        assert_eq!(
            packed_outline_alpha(0.3).to_bits(),
            (9.0f32 / 15.0).to_bits()
        );
        assert_eq!(packed_outline_alpha(0.5).to_bits(), 1.0f32.to_bits());
        assert_eq!(packed_outline_alpha(0.8).to_bits(), 1.0f32.to_bits());
    }

    #[test]
    fn outline_sampling_clamps_to_texture_edges() {
        assert_eq!(clamp_texture_coord(-1, 4), 0);
        assert_eq!(clamp_texture_coord(2, 4), 2);
        assert_eq!(clamp_texture_coord(4, 4), 3);
    }

    #[test]
    fn equal_transparent_depth_ties_follow_bottom_up_glsl_loop_order() {
        let renderer = RendererOptions::default();
        let sphere = crate::model::Boundary::from_positions(&[
            crate::model::Vec3::new(-1.0, 0.0, 0.0),
            crate::model::Vec3::new(1.0, 0.0, 0.0),
        ])
        .sphere;
        let camera = crate::render::camera::resolve_camera(&renderer, &sphere, 5, 5).unwrap();
        let opaque = Framebuffer::new(5, 5, [255; 4]);
        let mut transparent = Framebuffer::new(5, 5, [0; 4]);
        let mut transparent_alpha = vec![0.0; 25];
        let top = 7;
        let center = 12;
        let bottom = 17;
        transparent.depth01[top] = 0.5;
        transparent.depth01[bottom] = 0.5;
        transparent_alpha[top] = 0.4;
        transparent_alpha[bottom] = 0.2;
        let mut color = crate::render::postprocessing::read_color(&transparent);

        apply_transparent_outline(
            &opaque,
            &transparent,
            &mut color,
            &transparent_alpha,
            &[1.0; 25],
            &camera,
            &OutlineParams {
                threshold: 0.01,
                ..OutlineParams::default()
            },
            1.0,
            0.0,
        )
        .unwrap();

        // GLSL visits screen-bottom first and screen-top last. Its `<=` tie
        // therefore retains the top sample's alpha: 0.4 is doubled, packed to
        // the high nibble, and restored as 0.8 in the compose pass.
        assert_eq!(color[center].map(quantize), [0, 0, 0, 204]);
    }

    #[test]
    fn transparent_depth_edge_writes_premultiplied_four_bit_outline() {
        let renderer = RendererOptions::default();
        let sphere = crate::model::Boundary::from_positions(&[
            crate::model::Vec3::new(-1.0, 0.0, 0.0),
            crate::model::Vec3::new(1.0, 0.0, 0.0),
        ])
        .sphere;
        let camera = crate::render::camera::resolve_camera(&renderer, &sphere, 9, 9).unwrap();
        let opaque = Framebuffer::new(9, 9, [255; 4]);
        let (depth, alpha) = packed_depth_alpha_roundtrip(0.5, 0.3);
        let mut transparent = Framebuffer::new(9, 9, [0; 4]);
        let center = 4 * 9 + 4;
        transparent.depth01[center] = depth;
        transparent.color[center] = [15, 30, 45, quantize(alpha)];
        let mut transparent_alpha = vec![0.0; 81];
        transparent_alpha[center] = alpha;
        let mut color = crate::render::postprocessing::read_color(&transparent);

        apply_transparent_outline(
            &opaque,
            &transparent,
            &mut color,
            &transparent_alpha,
            &[1.0; 81],
            &camera,
            &OutlineParams::default(),
            1.0,
            0.0,
        )
        .unwrap();

        assert_eq!(color[center].map(quantize), [15, 30, 45, quantize(alpha)]);
        for (index, pixel) in color.iter().enumerate() {
            let x = index % 9;
            let y = index / 9;
            if x.abs_diff(4) <= 1 && y.abs_diff(4) <= 1 && index != center {
                assert_eq!(pixel.map(quantize), [0, 0, 0, 153]);
            } else if index != center {
                assert_eq!(pixel.map(quantize), [0; 4]);
            }
        }
    }
}
