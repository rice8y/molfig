use super::framebuffer::Framebuffer;
use super::options::MultiSampleOptions;

/// Mol* `JitterVectors`, after its fixed 1/16 scaling.
const JITTER_LEVEL_0: &[(f32, f32)] = &[(0.0, 0.0)];
const JITTER_LEVEL_1: &[(f32, f32)] = &[(0.0, 0.0), (-0.25, -0.25)];
const JITTER_LEVEL_2: &[(f32, f32)] =
    &[(0.0, 0.0), (0.375, -0.125), (-0.375, 0.125), (0.125, 0.375)];
const JITTER_LEVEL_3: &[(f32, f32)] = &[
    (0.0, 0.0),
    (-0.0625, 0.1875),
    (0.3125, 0.0625),
    (-0.1875, -0.3125),
    (-0.3125, 0.3125),
    (-0.4375, -0.0625),
    (0.1875, 0.4375),
    (0.4375, -0.4375),
];
const JITTER_LEVEL_4: &[(f32, f32)] = &[
    (0.0, 0.0),
    (-0.0625, -0.1875),
    (-0.1875, 0.125),
    (0.25, -0.0625),
    (-0.3125, -0.125),
    (0.125, 0.3125),
    (0.3125, 0.1875),
    (0.1875, -0.3125),
    (-0.125, 0.375),
    (0.0, -0.4375),
    (-0.25, -0.375),
    (-0.375, 0.25),
    (-0.5, 0.0),
    (0.4375, -0.25),
    (0.375, 0.4375),
    (-0.4375, -0.5),
];
const JITTER_LEVEL_5: &[(f32, f32)] = &[
    (0.0, 0.0),
    (-0.4375, -0.3125),
    (-0.1875, -0.3125),
    (-0.3125, -0.25),
    (-0.0625, -0.25),
    (-0.125, -0.125),
    (-0.375, -0.0625),
    (-0.25, 0.0),
    (-0.4375, 0.0625),
    (-0.0625, 0.125),
    (-0.375, 0.1875),
    (-0.1875, 0.1875),
    (-0.4375, 0.375),
    (-0.1875, 0.375),
    (-0.3125, 0.4375),
    (-0.0625, 0.4375),
    (0.3125, -0.4375),
    (0.0625, -0.375),
    (0.375, -0.3125),
    (0.25, -0.25),
    (0.125, -0.1875),
    (0.4375, -0.125),
    (0.0625, -0.0625),
    (0.25, -0.0625),
    (0.125, 0.0625),
    (0.375, 0.125),
    (0.0, 0.25),
    (0.25, 0.25),
    (0.125, 0.3125),
    (0.4375, 0.3125),
    (0.3125, 0.375),
    (0.1875, 0.4375),
];

pub(super) fn jitter_offsets(options: &MultiSampleOptions) -> &'static [(f32, f32)] {
    if options.mode == "off" {
        return JITTER_LEVEL_0;
    }
    match options.sample_level {
        0 => JITTER_LEVEL_0,
        1 => JITTER_LEVEL_1,
        2 => JITTER_LEVEL_2,
        3 => JITTER_LEVEL_3,
        4 => JITTER_LEVEL_4,
        _ => JITTER_LEVEL_5,
    }
}

pub(super) fn accumulate(
    accumulation: &mut [[f32; 4]],
    framebuffer: &Framebuffer,
    mode: &str,
    sample_index: usize,
    sample_count: usize,
) {
    let weight = sample_weight(mode, sample_index, sample_count);
    for (sum, color) in accumulation.iter_mut().zip(&framebuffer.color) {
        for channel in 0..4 {
            let source = f32::from(color[channel]) / 255.0 * weight;
            sum[channel] = truncate_to_f16_blend_target(sum[channel] + source);
        }
    }
}

/// Mol* uses equal weights for temporal accumulation. In its synchronous
/// `on` mode only, it varies the weights across a symmetric 1/32 range to
/// cancel systematic accumulation error in the fp16 compose target.
fn sample_weight(mode: &str, sample_index: usize, sample_count: usize) -> f32 {
    let sample_count = sample_count.max(1) as f64;
    let base = 1.0 / sample_count;
    if mode != "on" {
        return base as f32;
    }
    let centered = -0.5 + (sample_index as f64 + 0.5) / sample_count;
    (base + (1.0 / 32.0) * centered) as f32
}

pub(super) fn resolve(accumulation: Vec<[f32; 4]>) -> Vec<u8> {
    let mut rgba = Vec::with_capacity(accumulation.len() * 4);
    for sum in accumulation {
        for channel in sum {
            rgba.push(super::color::quantize(channel));
        }
    }
    rgba
}

/// Convert a finite non-negative blend result to the binary16 value stored by
/// the pinned ANGLE/Metal fp16 color attachment. Direct render-pass probes on
/// the reference Apple GPU show truncation toward zero at this attachment
/// boundary, rather than IEEE round-to-nearest-even.
fn truncate_to_f16_blend_target(value: f32) -> f32 {
    let bits = value.to_bits();
    let sign = ((bits >> 16) & 0x8000) as u16;
    let exponent = ((bits >> 23) & 0xff) as i32 - 127 + 15;
    let mantissa = bits & 0x007f_ffff;
    let half = if exponent <= 0 {
        if exponent < -10 {
            sign
        } else {
            let mantissa = mantissa | 0x0080_0000;
            let shift = 14 - exponent;
            sign | (mantissa >> shift) as u16
        }
    } else if exponent >= 31 {
        sign | 0x7c00
    } else {
        sign | ((exponent as u16) << 10) | (mantissa >> 13) as u16
    };
    f16_bits_to_f32(half)
}

fn f16_bits_to_f32(bits: u16) -> f32 {
    let sign = (u32::from(bits & 0x8000)) << 16;
    let exponent = u32::from((bits >> 10) & 0x1f);
    let mantissa = u32::from(bits & 0x03ff);
    let converted = if exponent == 0 {
        if mantissa == 0 {
            sign
        } else {
            let mut mantissa = mantissa;
            let mut exponent = -14i32;
            while mantissa & 0x0400 == 0 {
                mantissa <<= 1;
                exponent -= 1;
            }
            mantissa &= 0x03ff;
            sign | (((exponent + 127) as u32) << 23) | (mantissa << 13)
        }
    } else if exponent == 31 {
        sign | 0x7f80_0000 | (mantissa << 13)
    } else {
        sign | ((exponent + 112) << 23) | (mantissa << 13)
    };
    f32::from_bits(converted)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn molstar_default_jitter_vectors_are_level_two() {
        let options = MultiSampleOptions::default();
        assert_eq!(jitter_offsets(&options), JITTER_LEVEL_2);
    }

    #[test]
    fn temporal_fp16_blend_accumulation_matches_reference_metal_boundary_colors() {
        let colors = [[27, 158, 119, 255], [252, 251, 250, 255]];
        for (foreground_count, expected) in [
            (1, [195, 227, 217, 255]),
            (2, [139, 204, 184, 255]),
            (3, [83, 181, 152, 255]),
        ] {
            let mut accumulation = vec![[0.0; 4]];
            for sample_index in 0..4 {
                let framebuffer = Framebuffer {
                    width: 1,
                    height: 1,
                    color: vec![colors[usize::from(sample_index >= foreground_count)]],
                    depth: vec![0.0],
                    depth01: vec![0.0],
                    normal: vec![Default::default()],
                };
                accumulate(&mut accumulation, &framebuffer, "temporal", sample_index, 4);
            }
            assert_eq!(resolve(accumulation), expected);
        }
    }

    #[test]
    fn sample_weights_match_molstar_rounding_error_distribution() {
        let weights = (0..4)
            .map(|index| sample_weight("on", index, 4))
            .collect::<Vec<_>>();
        assert_eq!(
            weights,
            [0.238_281_25, 0.246_093_75, 0.253_906_25, 0.261_718_75]
        );
        assert_eq!(weights.iter().sum::<f32>(), 1.0);
        assert_eq!(sample_weight("temporal", 0, 4), 0.25);
        assert_eq!(sample_weight("temporal", 3, 4), 0.25);
    }
}
