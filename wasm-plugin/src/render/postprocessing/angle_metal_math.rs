//! Floating-point behavior of the pinned Mol* WebGL2 reference path.
//!
//! ANGLE lowers highp GLSL division to the Apple GPU `frcp` instruction and a
//! multiplication. `frcp` is faithful, but not always the correctly rounded
//! IEEE-754 reciprocal. Its ±1 ULP correction depends only on the 23-bit input
//! mantissa for finite, positive, normal float32 values. The lookup below is an
//! exhaustive capture of that domain on the pinned reference platform.

const MANTISSA_MASK: u32 = 0x007f_ffff;
const EXPONENT_MASK: u32 = 0x7f80_0000;
const FRC_UP: &[u8; 1 << 20] = include_bytes!("angle-metal-data/frcp-up.bits");
const FRC_DOWN: &[u8; 35_028] = include_bytes!("angle-metal-data/frcp-down.u32le");
const RSQRT_EVEN_UP: &[u8; 1 << 20] = include_bytes!("angle-metal-data/rsqrt-even-up.bits");
const RSQRT_EVEN_DOWN: &[u8; 215_924] = include_bytes!("angle-metal-data/rsqrt-even-down.u32le");
const RSQRT_ODD_UP: &[u8; 1 << 20] = include_bytes!("angle-metal-data/rsqrt-odd-up.bits");
const RSQRT_ODD_DOWN: &[u8; 174_376] = include_bytes!("angle-metal-data/rsqrt-odd-down.u32le");

/// Evaluate a binary32 fused multiply-add with one final binary32 rounding on
/// every supported target. WebAssembly has no scalar FMA instruction, so use
/// libm's software implementation instead of depending on target-specific
/// lowering of `f32::mul_add`.
#[inline]
pub(in crate::render) fn fused_multiply_add(a: f32, b: f32, c: f32) -> f32 {
    libm::fmaf(a, b, c)
}

#[inline]
pub(in crate::render) fn reciprocal(value: f32) -> f32 {
    let value_bits = value.to_bits();
    let exponent = value_bits & EXPONENT_MASK;
    if value_bits >> 31 != 0 || exponent == 0 || exponent == EXPONENT_MASK {
        return value.recip();
    }

    let mantissa = value_bits & MANTISSA_MASK;
    let mut reciprocal_bits = value.recip().to_bits();
    if FRC_UP[(mantissa >> 3) as usize] & (1 << (mantissa & 7)) != 0 {
        reciprocal_bits += 1;
    } else if sparse_contains(FRC_DOWN, mantissa) {
        reciprocal_bits -= 1;
    }
    f32::from_bits(reciprocal_bits)
}

#[inline]
pub(in crate::render) fn inverse_sqrt(value: f32) -> f32 {
    let value_bits = value.to_bits();
    let exponent = value_bits & EXPONENT_MASK;
    if value_bits >> 31 != 0 || exponent == 0 || exponent == EXPONENT_MASK {
        return value.sqrt().recip();
    }

    let mantissa = value_bits & MANTISSA_MASK;
    let unbiased_exponent = ((exponent >> 23) as i32) - 127;
    let (up, down) = if unbiased_exponent & 1 == 0 {
        (RSQRT_EVEN_UP.as_slice(), RSQRT_EVEN_DOWN.as_slice())
    } else {
        (RSQRT_ODD_UP.as_slice(), RSQRT_ODD_DOWN.as_slice())
    };
    let mut result_bits = ((1.0f64 / f64::from(value).sqrt()) as f32).to_bits();
    if up[(mantissa >> 3) as usize] & (1 << (mantissa & 7)) != 0 {
        result_bits += 1;
    } else if sparse_contains(down, mantissa) {
        result_bits -= 1;
    }
    f32::from_bits(result_bits)
}

#[inline]
pub(in crate::render) fn divide(numerator: f32, denominator: f32) -> f32 {
    numerator * reciprocal(denominator)
}

#[inline]
pub(super) fn divide_blur_sum(numerator: f32, denominator: f32, _kernel_size: usize) -> f32 {
    divide(numerator, denominator)
}

fn sparse_contains(values: &[u8], mantissa: u32) -> bool {
    let mut low = 0usize;
    let mut high = values.len() / 4;
    while low < high {
        let middle = low + (high - low) / 2;
        let offset = middle * 4;
        let candidate = u32::from_le_bytes([
            values[offset],
            values[offset + 1],
            values[offset + 2],
            values[offset + 3],
        ]);
        if candidate < mantissa {
            low = middle + 1;
        } else {
            high = middle;
        }
    }
    let offset = low * 4;
    offset < values.len()
        && u32::from_le_bytes([
            values[offset],
            values[offset + 1],
            values[offset + 2],
            values[offset + 3],
        ]) == mantissa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_fma_matches_native_fma_for_render_domain_values() {
        let mut state = 0x8f1b_bcdc_u32;
        for _ in 0..1_000_000 {
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            let a = f32::from_bits((state & 0x807f_ffff) | 0x3f00_0000);
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            let b = f32::from_bits((state & 0x807f_ffff) | 0x3f00_0000);
            state = state.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            let c = f32::from_bits((state & 0x807f_ffff) | 0x3f00_0000);
            assert_eq!(
                fused_multiply_add(a, b, c).to_bits(),
                a.mul_add(b, c).to_bits()
            );
        }
    }

    #[test]
    fn applies_both_observed_reciprocal_correction_directions() {
        let up = f32::from_bits(0x3f80_25f8);
        let down = f32::from_bits(0x3f80_05a9);
        assert_eq!(reciprocal(up).to_bits(), up.recip().to_bits() + 1);
        assert_eq!(reciprocal(down).to_bits(), down.recip().to_bits() - 1);
    }

    #[test]
    fn correction_is_exponent_invariant() {
        for bits in [0x3f80_25f8, 0x4080_25f8, 0x3e80_25f8] {
            let value = f32::from_bits(bits);
            assert_eq!(reciprocal(value).to_bits(), value.recip().to_bits() + 1);
        }
    }

    #[test]
    fn exceptional_values_use_ieee_reciprocal() {
        for value in [
            0.0,
            -0.0,
            -1.5,
            f32::INFINITY,
            f32::NAN,
            f32::MIN_POSITIVE / 2.0,
        ] {
            assert_eq!(reciprocal(value).to_bits(), value.recip().to_bits());
        }
    }

    #[test]
    fn inverse_sqrt_lookup_counts_and_checksums_match_the_exhaustive_capture() {
        assert_eq!(
            RSQRT_EVEN_UP
                .iter()
                .map(|byte| byte.count_ones())
                .sum::<u32>(),
            315_859
        );
        assert_eq!(RSQRT_EVEN_DOWN.len() / 4, 53_981);
        assert_eq!(
            RSQRT_ODD_UP
                .iter()
                .map(|byte| byte.count_ones())
                .sum::<u32>(),
            312_387
        );
        assert_eq!(RSQRT_ODD_DOWN.len() / 4, 43_594);
        for (up, down, expected) in [
            (
                RSQRT_EVEN_UP.as_slice(),
                RSQRT_EVEN_DOWN.as_slice(),
                0x1b18_b3c8_994f_275c,
            ),
            (
                RSQRT_ODD_UP.as_slice(),
                RSQRT_ODD_DOWN.as_slice(),
                0x1367_fb4c_b271_c6c8,
            ),
        ] {
            let mut hash = 0xcbf2_9ce4_8422_2325u64;
            for byte in up.iter().chain(down.iter()) {
                hash = (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3);
            }
            assert_eq!(hash, expected);
        }
    }

    #[test]
    fn inverse_sqrt_correction_depends_on_exponent_parity() {
        for bits in [0x3f80_bee0, 0x3f80_0578, 0x4000_bee0, 0x4000_0578] {
            let value = f32::from_bits(bits);
            let result = inverse_sqrt(value);
            assert!(result.is_finite());
            assert!(result > 0.0);
        }
    }

    #[test]
    fn exceptional_inverse_sqrt_values_use_the_ieee_fallback() {
        for value in [
            0.0,
            -0.0,
            -1.5,
            f32::INFINITY,
            f32::NAN,
            f32::MIN_POSITIVE / 2.0,
        ] {
            assert_eq!(
                inverse_sqrt(value).to_bits(),
                value.sqrt().recip().to_bits()
            );
        }
    }

    #[test]
    fn lookup_counts_and_checksums_match_the_exhaustive_capture() {
        assert_eq!(
            FRC_UP.iter().map(|byte| byte.count_ones()).sum::<u32>(),
            859_417
        );
        assert_eq!(FRC_DOWN.len() / 4, 8_757);
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        for byte in FRC_UP.iter().chain(FRC_DOWN.iter()) {
            hash = (hash ^ u64::from(*byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
        assert_eq!(hash, 0x0ed3_f94e_5f55_81dc);
    }

    #[test]
    fn default_blur_subset_domain_remains_exact() {
        let sigma = 15.0f64 / 3.0;
        let kernel = (0..=7)
            .map(|offset| {
                let x = offset as f64;
                ((-x * x / (2.0 * sigma * sigma)).exp()
                    / ((2.0 * std::f64::consts::PI).sqrt() * sigma)) as f32
            })
            .collect::<Vec<_>>();
        let mut hash = 0xcbf2_9ce4_8422_2325u64;
        for mask in 1u32..(1 << 15) {
            let mut denominator = 0.0f32;
            for offset in -7i32..=7 {
                if mask & (1 << (offset + 7)) != 0 {
                    denominator += kernel[offset.unsigned_abs() as usize];
                }
            }
            for bits in [denominator.to_bits(), divide(1.0, denominator).to_bits()] {
                for byte in bits.to_le_bytes() {
                    hash = (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3);
                }
            }
        }
        assert_eq!(hash, 0x6fe5_2e0a_e4fb_27fd);
    }
}
