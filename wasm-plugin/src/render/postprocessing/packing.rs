use super::super::color::quantize;
use super::angle_metal_math::fused_multiply_add;

/// Simulate Mol*'s `packUnitIntervalToRG` write to an RGBA8 render target and
/// the following `unpackRGToUnitInterval` texture read.
pub(super) fn packed_unit_interval_roundtrip(value: f32) -> f32 {
    unpack_unit_interval_rgba8(pack_unit_interval_rgba8(value))
}

pub(super) fn unpack_unit_interval_rgba8([x, y]: [u8; 2]) -> f32 {
    let x = x as f32 / 255.0;
    let y = y as f32 / 255.0;
    unpack_unit_interval_rg([x, y])
}

pub(super) fn unpack_unit_interval_rg([x, y]: [f32; 2]) -> f32 {
    fused_multiply_add(y, 255.0 / 256.0, x * (255.0 / (256.0 * 256.0)))
}

pub(super) fn pack_unit_interval_rgba8(value: f32) -> [u8; 2] {
    let value = value.clamp(0.0, 1.0);
    let mut x = (value * 256.0).fract();
    let mut y = value - x * (1.0 / 256.0);
    x *= 256.0 / 255.0;
    y *= 256.0 / 255.0;
    [quantize(x), quantize(y)]
}

/// Simulate `packDepthWithAlphaToRGBA` followed by an RGBA8 attachment write
/// and `unpackRGBAToDepthWithAlpha`. Mol* uses this path for the nearest
/// transparent depth texture consumed by SSAO and outlines.
#[cfg(test)]
pub(in crate::render) fn packed_depth_alpha_roundtrip(depth: f32, alpha: f32) -> (f32, f32) {
    let [x, y, z, a] = pack_depth_alpha_rgba8(depth, alpha);
    unpack_depth_alpha_rgba([
        x as f32 / 255.0,
        y as f32 / 255.0,
        z as f32 / 255.0,
        a as f32 / 255.0,
    ])
}

pub(in crate::render) fn unpack_depth_alpha_rgba([x, y, z, alpha]: [f32; 4]) -> (f32, f32) {
    let downscale = 255.0 / 256.0;
    let unpacked_depth = fused_multiply_add(
        z,
        downscale,
        fused_multiply_add(y, downscale / 256.0, x * (downscale / 65_536.0)),
    );
    (unpacked_depth, alpha)
}

pub(in crate::render) fn pack_depth_alpha_rgba8(depth: f32, alpha: f32) -> [u8; 4] {
    let depth = depth.clamp(0.0, 1.0);
    let x = (depth * 65_536.0).fract();
    let y = (depth * 256.0).fract();
    let z = depth;
    let upscale = 256.0 / 255.0;
    [
        quantize(x * upscale),
        quantize((y - x / 256.0) * upscale),
        quantize((z - y / 256.0) * upscale),
        quantize(alpha),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rg8_unit_interval_roundtrip_matches_shader_normalization() {
        assert_eq!(pack_unit_interval_rgba8(1.0), [0, 255]);
        assert_eq!(unpack_unit_interval_rgba8([0, 255]), 255.0 / 256.0);
        assert_eq!(unpack_unit_interval_rgba8([17, 203]).to_bits(), 0x3f4b_1100);
    }

    #[test]
    fn valid_rg8_unit_interval_values_repack_without_loss() {
        for y in 0..255 {
            for x in 0..=255 {
                let packed = [x, y];
                assert_eq!(
                    pack_unit_interval_rgba8(unpack_unit_interval_rgba8(packed)),
                    packed
                );
            }
        }
        assert_eq!(
            pack_unit_interval_rgba8(unpack_unit_interval_rgba8([0, 255])),
            [0, 255]
        );
    }

    #[test]
    fn transparent_depth_roundtrip_preserves_packed_alpha_and_depth() {
        assert_eq!(pack_depth_alpha_rgba8(0.5, 0.3), [0, 0, 128, 77]);
        let (depth, alpha) = packed_depth_alpha_roundtrip(0.5, 0.3);
        assert_eq!(depth, 0.5);
        assert_eq!(alpha, 77.0 / 255.0);

        let (depth, alpha) = packed_depth_alpha_roundtrip(0.123_456_7, 0.6);
        assert!((depth - 0.123_456_7).abs() < 1.0e-7);
        assert_eq!(alpha, 153.0 / 255.0);
    }

    #[test]
    fn valid_packed_depth_alpha_values_repack_without_loss() {
        for index in 0..=1024 {
            let depth = index as f32 / 1024.0;
            let alpha = ((index * 73) % 1025) as f32 / 1024.0;
            let packed = pack_depth_alpha_rgba8(depth, alpha);
            let (unpacked_depth, unpacked_alpha) =
                unpack_depth_alpha_rgba(packed.map(|channel| channel as f32 / 255.0));
            assert_eq!(
                pack_depth_alpha_rgba8(unpacked_depth, unpacked_alpha),
                packed
            );
        }
    }
}
