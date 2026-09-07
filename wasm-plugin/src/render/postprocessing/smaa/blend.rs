//! SMAA gamma-correct neighborhood blending pass.

use super::texture::{quantize, sample_render_texel};
use crate::render::framebuffer::Framebuffer;

pub(super) fn apply(framebuffer: &mut Framebuffer, weights: &[[u8; 4]]) {
    let width = framebuffer.width;
    let height = framebuffer.height;
    let source = framebuffer.color.clone();
    let mut output = source.clone();

    for y in 0..height {
        for x in 0..width {
            let index = y * width + x;
            let center_weights =
                sample_render_texel(weights, width, height, x as isize, y as isize);
            let bottom_weights =
                sample_render_texel(weights, width, height, x as isize, y as isize + 1);
            let right_weights =
                sample_render_texel(weights, width, height, x as isize + 1, y as isize);
            let crossing = [
                center_weights[0],
                bottom_weights[1],
                center_weights[2],
                right_weights[3],
            ];
            if crossing.iter().sum::<f32>() < 1e-5 {
                continue;
            }

            let mut offset = choose_offset(crossing);
            if offset[0].abs() > offset[1].abs() {
                offset[1] = 0.0;
            } else {
                offset[0] = 0.0;
            }
            let opposite_x = x as isize + glsl_sign(offset[0]);
            // Positive texture-space Y points upward, while framebuffer rows
            // are stored top-down.
            let opposite_y = y as isize - glsl_sign(offset[1]);
            let center = sample_render_texel(&source, width, height, x as isize, y as isize);
            let opposite = sample_render_texel(&source, width, height, opposite_x, opposite_y);
            let amount = offset[0].abs().max(offset[1].abs());
            let mut mixed = [0.0; 4];
            for channel in 0..3 {
                mixed[channel] = mix(
                    center[channel].powf(2.2),
                    opposite[channel].powf(2.2),
                    amount,
                )
                .powf(1.0 / 2.2);
            }
            mixed[3] = mix(center[3], opposite[3], amount);
            output[index] = mixed.map(quantize);
        }
    }
    framebuffer.color = output;
}

#[inline]
fn mix(a: f32, b: f32, amount: f32) -> f32 {
    a * (1.0 - amount) + b * amount
}

#[inline]
fn glsl_sign(value: f32) -> isize {
    if value > 0.0 {
        1
    } else if value < 0.0 {
        -1
    } else {
        0
    }
}

#[inline]
fn choose_offset(weights: [f32; 4]) -> [f32; 2] {
    [
        if weights[3] > weights[2] {
            weights[3]
        } else {
            -weights[2]
        },
        if weights[1] > weights[0] {
            -weights[1]
        } else {
            weights[0]
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::{choose_offset, glsl_sign};

    #[test]
    fn shader_swizzles_horizontal_and_vertical_weight_pairs() {
        assert_eq!(choose_offset([0.1, 0.2, 0.3, 0.4]), [0.4, -0.2]);
        assert_eq!(choose_offset([0.5, 0.2, 0.7, 0.4]), [-0.7, 0.5]);
    }

    #[test]
    fn shader_sign_preserves_zero() {
        assert_eq!(glsl_sign(-0.5), -1);
        assert_eq!(glsl_sign(0.0), 0);
        assert_eq!(glsl_sign(0.5), 1);
    }
}
