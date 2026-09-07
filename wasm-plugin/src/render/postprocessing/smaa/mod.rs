//! CPU port of Mol*'s WebGL SMAA 1x Medium pipeline.
//!
//! The pass order, lookup textures, RGBA8 intermediate targets, thresholds,
//! search limits, and gamma-correct neighborhood blend follow the pinned Mol*
//! `SmaaPass` and its SMAA v2.8 shaders.

mod blend;
mod edges;
mod texture;
mod varying;
mod weights;

use super::super::framebuffer::Framebuffer;
use super::super::options::SmaaParams;

pub(super) fn apply(framebuffer: &mut Framebuffer, params: &SmaaParams) {
    if framebuffer.width == 0 || framebuffer.height == 0 {
        return;
    }
    let edges = edges::detect(framebuffer, params.edge_threshold);
    let weights = weights::calculate(
        &edges,
        framebuffer.width,
        framebuffer.height,
        params.max_search_steps,
    );
    blend::apply(framebuffer, &weights);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_lookup_textures_have_expected_layout() {
        assert_eq!(texture::AREA_RGB.len(), 160 * 560 * 3);
        assert_eq!(texture::SEARCH_R.len(), 66 * 33);
        assert_eq!(fnv1a64(texture::AREA_RGB), 0xaf03_9577_0c5c_33d3);
        assert_eq!(fnv1a64(texture::SEARCH_R), 0x0af4_9891_dcb1_560d);
    }

    #[test]
    fn uniform_image_is_unchanged() {
        let mut framebuffer = Framebuffer::new(8, 6, [27, 158, 119, 255]);
        let before = framebuffer.color.clone();
        apply(&mut framebuffer, &SmaaParams::default());
        assert_eq!(framebuffer.color, before);
    }

    fn fnv1a64(bytes: &[u8]) -> u64 {
        bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
            (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3)
        })
    }
}
