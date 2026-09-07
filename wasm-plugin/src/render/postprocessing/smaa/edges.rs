//! SMAA color edge detection pass.

use super::texture::sample_render_texel;
use crate::render::framebuffer::Framebuffer;

pub(super) fn detect(framebuffer: &Framebuffer, threshold: f32) -> Vec<[u8; 4]> {
    let width = framebuffer.width;
    let height = framebuffer.height;
    let mut output = vec![[0, 0, 0, 255]; width * height];

    for y in 0..height {
        for x in 0..width {
            let center = rgb(sample_render_texel(
                &framebuffer.color,
                width,
                height,
                x as isize,
                y as isize,
            ));
            let left = rgb(sample_render_texel(
                &framebuffer.color,
                width,
                height,
                x as isize - 1,
                y as isize,
            ));
            let top = rgb(sample_render_texel(
                &framebuffer.color,
                width,
                height,
                x as isize,
                y as isize - 1,
            ));
            let delta_left = color_delta(center, left);
            let delta_top = color_delta(center, top);
            let mut edge_left = delta_left >= threshold;
            let mut edge_top = delta_top >= threshold;
            if !edge_left && !edge_top {
                continue;
            }

            let right = rgb(sample_render_texel(
                &framebuffer.color,
                width,
                height,
                x as isize + 1,
                y as isize,
            ));
            let bottom = rgb(sample_render_texel(
                &framebuffer.color,
                width,
                height,
                x as isize,
                y as isize + 1,
            ));
            let left_left = rgb(sample_render_texel(
                &framebuffer.color,
                width,
                height,
                x as isize - 2,
                y as isize,
            ));
            let top_top = rgb(sample_render_texel(
                &framebuffer.color,
                width,
                height,
                x as isize,
                y as isize - 2,
            ));
            let max_delta = delta_left
                .max(delta_top)
                .max(color_delta(center, right))
                .max(color_delta(center, bottom))
                .max(color_delta(center, left_left))
                .max(color_delta(center, top_top));
            edge_left &= delta_left >= 0.5 * max_delta;
            edge_top &= delta_top >= 0.5 * max_delta;
            output[y * width + x] = [
                if edge_left { 255 } else { 0 },
                if edge_top { 255 } else { 0 },
                0,
                0,
            ];
        }
    }
    output
}

#[inline]
fn rgb(color: [f32; 4]) -> [f32; 3] {
    [color[0], color[1], color[2]]
}

#[inline]
fn color_delta(a: [f32; 3], b: [f32; 3]) -> f32 {
    (a[0] - b[0])
        .abs()
        .max((a[1] - b[1]).abs())
        .max((a[2] - b[2]).abs())
}
