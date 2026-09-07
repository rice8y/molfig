//! Floating-point values within Mol*'s postprocessing fragment shader.
//!
//! Geometry and WBOIT evaluation have already written RGBA8 textures. SSAO,
//! outlines, and transparency blending do not write another texture between
//! operations: only the completed fragment is converted back to RGBA8.

use crate::render::color::quantize;
use crate::render::framebuffer::Framebuffer;

pub(in crate::render) fn read_color(framebuffer: &Framebuffer) -> Vec<[f32; 4]> {
    framebuffer
        .color
        .iter()
        .map(|pixel| pixel.map(|channel| f32::from(channel) / 255.0))
        .collect()
}

pub(in crate::render) fn write_composed_color(
    framebuffer: &mut Framebuffer,
    opaque: &[[f32; 4]],
    transparent: Option<&[[f32; 4]]>,
) {
    assert_eq!(framebuffer.color.len(), opaque.len());
    if let Some(layer) = transparent {
        assert_eq!(opaque.len(), layer.len());
    }
    for (index, destination) in framebuffer.color.iter_mut().enumerate() {
        let mut color = opaque[index];
        if let Some(layer) = transparent {
            let source = layer[index];
            if source[3] != 0.0 {
                for channel in 0..4 {
                    color[channel] = source[channel] + color[channel] * (1.0 - source[3]);
                }
            }
        }
        *destination = color.map(quantize);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn composition_rounds_only_after_transparent_occlusion_and_blending() {
        let mut framebuffer = Framebuffer::new(1, 1, [252, 251, 250, 255]);
        let opaque = read_color(&framebuffer);
        let source = [18, 110, 83, 178].map(|channel| channel as f32 / 255.0);
        let factor = 0.713;
        let layer = [[
            source[0] * factor,
            source[1] * factor,
            source[2] * factor,
            source[3],
        ]];
        write_composed_color(&mut framebuffer, &opaque, Some(&layer));
        let expected = [89, 154, 135, 255];
        assert_eq!(framebuffer.color[0], expected);

        let prematurely_rounded = layer.map(|pixel| pixel.map(|v| f32::from(quantize(v)) / 255.0));
        write_composed_color(&mut framebuffer, &opaque, Some(&prematurely_rounded));
        assert_ne!(framebuffer.color[0], expected);
    }

    #[test]
    fn outline_alpha_and_transparent_background_remain_float_until_output() {
        let mut framebuffer = Framebuffer::new(1, 1, [0; 4]);
        let opaque = [[0.13, 0.07, 0.04, 0.37]];
        let transparent = [[0.11, 0.23, 0.19, 0.63]];
        write_composed_color(&mut framebuffer, &opaque, Some(&transparent));
        assert_eq!(framebuffer.color[0], [40, 65, 52, 196]);
        // A disabled transparent blend must preserve all opaque channels.
        write_composed_color(&mut framebuffer, &opaque, Some(&[[1.0, 1.0, 1.0, 0.0]]));
        assert_eq!(framebuffer.color[0], opaque[0].map(quantize));
    }
}
