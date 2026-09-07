//! Interpolation of the full-screen SMAA vertex shader's search offsets.

use super::super::angle_metal_math::fused_multiply_add;

const TILE_SIZE: usize = 32;

/// Evaluate a normalized axis varying from its vertex endpoints. The reference
/// interpolator rounds a tile origin before evaluating pixel centers within
/// that tile, then truncates the positive or negative result toward zero.
/// Adding an offset to an already rounded fragment UV loses these boundaries.
pub(super) fn axis_offsets(extent: usize, offset: f32) -> Vec<f32> {
    let inverse = 1.0 / extent as f32;
    let scale = extent as f32 * inverse;
    let low = offset * inverse;
    let high = fused_multiply_add(offset, inverse, scale);
    let slope = (high - low) / extent as f32;
    let mut values = Vec::with_capacity(extent);
    for origin in (0..extent).step_by(TILE_SIZE) {
        let base = fused_multiply_add(origin as f32, slope, low);
        for within in 0..TILE_SIZE.min(extent - origin) {
            let value = f64::from(base) + (within as f64 + 0.5) * f64::from(slope);
            values.push(truncate_to_f32(value));
        }
    }
    values
}

#[inline]
fn truncate_to_f32(value: f64) -> f32 {
    let rounded = value as f32;
    if f64::from(rounded.abs()) > value.abs() {
        f32::from_bits(rounded.to_bits() - 1)
    } else {
        rounded
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn power_of_two_axis_preserves_exact_vertex_offsets() {
        for extent in [1, 16, 64, 256] {
            for offset in [-1.25, -0.25, -0.125, 0.125, 0.25, 1.25] {
                let values = axis_offsets(extent, offset);
                for (coordinate, value) in values.into_iter().enumerate() {
                    assert_eq!(value, (coordinate as f32 + 0.5 + offset) / extent as f32);
                }
            }
        }
    }

    #[test]
    fn non_power_of_two_axis_rounds_each_tile_origin_once() {
        let extent = 157;
        let inverse = 1.0 / extent as f32;
        let values = axis_offsets(extent, 0.0);
        let mut differs_from_fragment_division = false;
        for (coordinate, &value) in values.iter().enumerate() {
            let origin = coordinate / TILE_SIZE * TILE_SIZE;
            let base = fused_multiply_add(origin as f32, inverse, 0.0);
            let exact = f64::from(base) + ((coordinate - origin) as f64 + 0.5) * f64::from(inverse);
            assert_eq!(value, truncate_to_f32(exact));
            differs_from_fragment_division |= value != (coordinate as f32 + 0.5) / extent as f32;
        }
        assert!(differs_from_fragment_division);
    }

    #[test]
    fn interpolation_truncation_is_sign_symmetric_and_preserves_zero() {
        let halfway = 1.0 + 3.0 * f64::from(f32::EPSILON) / 4.0;
        assert_eq!(truncate_to_f32(halfway), 1.0);
        assert_eq!(truncate_to_f32(-halfway), -1.0);
        assert_eq!(truncate_to_f32(0.0).to_bits(), 0);
        assert_eq!(truncate_to_f32(-0.0).to_bits(), (-0.0f32).to_bits());
    }
}
