//! Texture sampling rules used by the pinned Mol* SMAA shaders.

pub(super) const AREA_WIDTH: usize = 160;
pub(super) const AREA_HEIGHT: usize = 560;
pub(super) const SEARCH_WIDTH: usize = 66;
pub(super) const SEARCH_HEIGHT: usize = 33;

pub(super) static AREA_RGB: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/smaa-area.rgb"));
pub(super) static SEARCH_R: &[u8] = include_bytes!(concat!(env!("OUT_DIR"), "/smaa-search.r"));

#[inline]
pub(super) fn sample_render_texel(
    texture: &[[u8; 4]],
    width: usize,
    height: usize,
    x: isize,
    y: isize,
) -> [f32; 4] {
    let x = x.clamp(0, width as isize - 1) as usize;
    let y = y.clamp(0, height as isize - 1) as usize;
    texture[y * width + x].map(|channel| channel as f32 / 255.0)
}

#[inline]
pub(super) fn sample_render_linear(
    texture: &[[u8; 4]],
    width: usize,
    height: usize,
    uv: [f32; 2],
) -> [f32; 4] {
    let x = uv[0] * width as f32 - 0.5;
    let y_gl = uv[1] * height as f32 - 0.5;
    let x0 = x.floor() as isize;
    let y0 = y_gl.floor() as isize;
    let tx = x - x0 as f32;
    let ty = y_gl - y0 as f32;
    let a = render_texel(texture, width, height, x0, y0);
    let b = render_texel(texture, width, height, x0 + 1, y0);
    let c = render_texel(texture, width, height, x0, y0 + 1);
    let d = render_texel(texture, width, height, x0 + 1, y0 + 1);
    let mut result = [0.0; 4];
    for channel in 0..4 {
        let lower = mix(a[channel], b[channel], tx);
        let upper = mix(c[channel], d[channel], tx);
        result[channel] = mix(lower, upper, ty);
    }
    result
}

#[inline]
pub(super) fn sample_area_linear(uv: [f32; 2]) -> [f32; 2] {
    let x = uv[0] * AREA_WIDTH as f32 - 0.5;
    let y = uv[1] * AREA_HEIGHT as f32 - 0.5;
    let x0 = x.floor() as isize;
    let y0 = y.floor() as isize;
    let tx = quantize_area_fraction(x - x0 as f32);
    let ty = quantize_area_fraction(y - y0 as f32);
    let a = area_texel(x0, y0);
    let b = area_texel(x0 + 1, y0);
    let c = area_texel(x0, y0 + 1);
    let d = area_texel(x0 + 1, y0 + 1);
    let filtered = [0, 1].map(|channel| {
        let upper = (b[channel] - a[channel]).mul_add(tx, a[channel]);
        let lower = (d[channel] - c[channel]).mul_add(tx, c[channel]);
        (lower - upper).mul_add(ty, upper)
    });
    filtered.map(quantize_linear_sample)
}

/// The Apple/ANGLE RGBA8 reference resolves normalized linear-filter results
/// to 1/16-byte precision before exposing the channel to the shader.
#[inline]
fn quantize_linear_sample(value: f32) -> f32 {
    (value * 16.0).round() / (255.0 * 16.0)
}

#[inline]
fn quantize_area_fraction(value: f32) -> f32 {
    (value * 256.0).round() / 256.0
}

#[inline]
pub(super) fn sample_search_nearest(uv: [f32; 2]) -> f32 {
    let x = ((uv[0] * SEARCH_WIDTH as f32).floor() as isize).clamp(0, SEARCH_WIDTH as isize - 1)
        as usize;
    let y = ((uv[1] * SEARCH_HEIGHT as f32).floor() as isize).clamp(0, SEARCH_HEIGHT as isize - 1)
        as usize;
    SEARCH_R[y * SEARCH_WIDTH + x] as f32 / 255.0
}

#[inline]
fn render_texel(
    texture: &[[u8; 4]],
    width: usize,
    height: usize,
    x: isize,
    y_gl: isize,
) -> [f32; 4] {
    let x = x.clamp(0, width as isize - 1) as usize;
    let y_gl = y_gl.clamp(0, height as isize - 1) as usize;
    let y = height - 1 - y_gl;
    let pixel = texture[y * width + x];
    pixel.map(|channel| channel as f32 / 255.0)
}

#[inline]
fn area_texel(x: isize, y: isize) -> [f32; 2] {
    let x = x.clamp(0, AREA_WIDTH as isize - 1) as usize;
    let y = y.clamp(0, AREA_HEIGHT as isize - 1) as usize;
    let offset = (y * AREA_WIDTH + x) * 3;
    [AREA_RGB[offset] as f32, AREA_RGB[offset + 1] as f32]
}

#[inline]
fn mix(a: f32, b: f32, amount: f32) -> f32 {
    a * (1.0 - amount) + b * amount
}

pub(super) use crate::render::color::quantize;

#[cfg(test)]
mod tests {
    use super::{quantize_area_fraction, quantize_linear_sample, sample_area_linear};

    #[test]
    fn area_lookup_uses_eight_bit_subtexel_coordinates() {
        assert_eq!(quantize_area_fraction(0.5009), 0.5);
        assert_eq!(quantize_area_fraction(0.5021), 0.503_906_25);
    }

    #[test]
    fn normalized_linear_result_uses_sixteen_steps_per_byte() {
        assert_eq!(quantize_linear_sample(91.48), 91.5 / 255.0);
        assert_eq!(quantize_linear_sample(91.53), 91.5 / 255.0);
    }

    #[test]
    fn area_lookup_matches_reference_sampler_boundary() {
        let uv = [(16.0 + 0.5) / 160.0, (49.732_05 + 0.5) / 560.0];
        let sampled = sample_area_linear(uv);
        assert_eq!(sampled[0], 0.0);
        assert_eq!(sampled[1], 91.5 / 255.0);
    }
}
