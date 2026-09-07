//! SMAA blending-weight calculation pass.

use super::super::angle_metal_math::fused_multiply_add;
use super::texture::{
    quantize, sample_area_linear, sample_render_linear, sample_render_texel, sample_search_nearest,
};
use super::varying::axis_offsets;

const AREA_MAX_DISTANCE: f32 = 16.0;
const AREA_PIXEL_SIZE: [f32; 2] = [1.0 / 160.0, 1.0 / 560.0];

struct SearchContext<'a> {
    edges: &'a [[u8; 4]],
    width: usize,
    height: usize,
    inverse: [f32; 2],
    max_steps: usize,
}

pub(super) fn calculate(
    edges: &[[u8; 4]],
    width: usize,
    height: usize,
    max_steps: usize,
) -> Vec<[u8; 4]> {
    let inverse = [1.0 / width as f32, 1.0 / height as f32];
    let context = SearchContext {
        edges,
        width,
        height,
        inverse,
        max_steps,
    };
    let mut output = vec![[0; 4]; width * height];
    let x_left = axis_offsets(width, -0.25);
    let x_right = axis_offsets(width, 1.25);
    let x_crossing = axis_offsets(width, -0.125);
    let y_crossing = axis_offsets(height, 0.125);
    let y_up = axis_offsets(height, 0.25);
    let y_down = axis_offsets(height, -1.25);

    for y in 0..height {
        let y_gl = height - 1 - y;
        for x in 0..width {
            // This varying is defined in pixel units by the pinned vertex
            // shader. At a fragment center its mathematical value is exact;
            // deriving it again through reciprocal texture sizes introduces
            // avoidable f32 error for non-power-of-two dimensions.
            let pixel = [x as f32 + 0.5, height as f32 - y as f32 - 0.5];
            let offset0 = [x_left[x], y_crossing[y_gl], x_right[x], y_crossing[y_gl]];
            let offset1 = [x_crossing[x], y_up[y_gl], x_crossing[x], y_down[y_gl]];
            let search_distance = 2.0 * max_steps as f32;
            let offset2 = [
                offset0[0] - search_distance * inverse[0],
                offset0[2] + search_distance * inverse[0],
                offset1[1] - search_distance * inverse[1],
                offset1[3] + search_distance * inverse[1],
            ];
            let mut weights = [0.0; 4];
            // This lookup is at the current fragment center. Reading the
            // corresponding texel directly avoids introducing CPU reciprocal
            // error into a sample that raster texture filtering resolves to
            // the exact texel center.
            let edge = sample_render_texel(edges, width, height, x as isize, y as isize);

            if edge[1] > 0.0 {
                let mut coordinates = [
                    search_x_left(&context, [offset0[0], offset0[1]], offset2[0]),
                    offset1[1],
                ];
                let left = context.sample(coordinates)[0];
                let left_distance = coordinates[0];
                coordinates[0] = search_x_right(&context, [offset0[2], offset0[3]], offset2[1]);
                let distances = [
                    pixel_distance(left_distance, inverse[0], pixel[0]),
                    pixel_distance(coordinates[0], inverse[0], pixel[0]),
                ];
                coordinates[1] -= inverse[1];
                // The pinned GLSL macro converts `ivec2(1, 0)` with
                // `float(offset)`. GLSL scalar construction takes the first
                // vector component, so the resulting 1.0 is applied to both
                // axes through multiplication by `uTexSizeInv`.
                let right = sample_level_zero_offset(&context, coordinates, [1, 0])[0];
                let area = area(
                    [distances[0].abs().sqrt(), distances[1].abs().sqrt()],
                    left,
                    right,
                );
                weights[0] = area[0];
                weights[1] = area[1];
            }

            if edge[0] > 0.0 {
                let mut coordinates = [
                    offset0[0],
                    search_y_up(&context, [offset1[0], offset1[1]], offset2[2]),
                ];
                let top = context.sample(coordinates)[1];
                let top_distance = coordinates[1];
                coordinates[1] = search_y_down(&context, [offset1[2], offset1[3]], offset2[3]);
                let distances = [
                    pixel_distance(top_distance, inverse[1], pixel[1]),
                    pixel_distance(coordinates[1], inverse[1], pixel[1]),
                ];
                coordinates[1] -= inverse[1];
                // Likewise, `float(ivec2(0, 1))` becomes 0.0 in the pinned
                // shader, so this sample receives no additional offset.
                let bottom = sample_level_zero_offset(&context, coordinates, [0, 1])[1];
                let area = area(
                    [distances[0].abs().sqrt(), distances[1].abs().sqrt()],
                    top,
                    bottom,
                );
                weights[2] = area[0];
                weights[3] = area[1];
            }
            output[y * width + x] = weights.map(quantize);
        }
    }
    output
}

/// Division and subtraction contract into a reciprocal multiply-add in the
/// reference shader. Rounding the quotient first can erase a small distance;
/// its subsequent square root then selects a different area-texture sample.
#[inline]
fn pixel_distance(coordinate: f32, inverse_extent: f32, pixel: f32) -> f32 {
    fused_multiply_add(coordinate, inverse_extent.recip(), -pixel)
}

impl SearchContext<'_> {
    #[inline]
    fn sample(&self, uv: [f32; 2]) -> [f32; 4] {
        sample_render_linear(self.edges, self.width, self.height, uv)
    }
}

/// Mirrors the pinned GLSL `float(offset)` conversion in
/// `SMAASampleLevelZeroOffset`. Constructing a scalar from an `ivec2` selects
/// its first component, and that scalar is then multiplied by both inverse
/// texture dimensions.
#[inline]
fn sample_level_zero_offset(
    context: &SearchContext<'_>,
    coordinate: [f32; 2],
    offset: [i32; 2],
) -> [f32; 4] {
    let scalar_offset = offset[0] as f32;
    context.sample([
        coordinate[0] + scalar_offset * context.inverse[0],
        coordinate[1] + scalar_offset * context.inverse[1],
    ])
}

fn search_x_left(context: &SearchContext<'_>, mut coordinate: [f32; 2], end: f32) -> f32 {
    let mut edge = [0.0, 1.0];
    for _ in 0..context.max_steps {
        let sampled = context.sample(coordinate);
        edge = [sampled[0], sampled[1]];
        coordinate[0] -= 2.0 * context.inverse[0];
        if !(coordinate[0] > end && edge[1] > 0.8281 && edge[0] == 0.0) {
            break;
        }
    }
    // The vertex offset, one-pixel search bias, and final two-pixel step
    // combine before the multiply-add in the reference shader.
    coordinate[0] = fused_multiply_add(3.25, context.inverse[0], coordinate[0]);
    coordinate[0] -= context.inverse[0] * search_length(edge, 0.0, 0.5);
    coordinate[0]
}

fn search_x_right(context: &SearchContext<'_>, mut coordinate: [f32; 2], end: f32) -> f32 {
    let mut edge = [0.0, 1.0];
    for _ in 0..context.max_steps {
        let sampled = context.sample(coordinate);
        edge = [sampled[0], sampled[1]];
        coordinate[0] += 2.0 * context.inverse[0];
        if !(coordinate[0] < end && edge[1] > 0.8281 && edge[0] == 0.0) {
            break;
        }
    }
    coordinate[0] = fused_multiply_add(-3.25, context.inverse[0], coordinate[0]);
    coordinate[0] += context.inverse[0] * search_length(edge, 0.5, 0.5);
    coordinate[0]
}

fn search_y_up(context: &SearchContext<'_>, mut coordinate: [f32; 2], end: f32) -> f32 {
    let mut edge = [1.0, 0.0];
    for _ in 0..context.max_steps {
        let sampled = context.sample(coordinate);
        edge = [sampled[0], sampled[1]];
        coordinate[1] += 2.0 * context.inverse[1];
        if !(coordinate[1] > end && edge[0] > 0.8281 && edge[1] == 0.0) {
            break;
        }
    }
    coordinate[1] = fused_multiply_add(-3.25, context.inverse[1], coordinate[1]);
    coordinate[1] += context.inverse[1] * search_length([edge[1], edge[0]], 0.0, 0.5);
    coordinate[1]
}

fn search_y_down(context: &SearchContext<'_>, mut coordinate: [f32; 2], end: f32) -> f32 {
    let mut edge = [1.0, 0.0];
    for _ in 0..context.max_steps {
        let sampled = context.sample(coordinate);
        edge = [sampled[0], sampled[1]];
        coordinate[1] -= 2.0 * context.inverse[1];
        if !(coordinate[1] < end && edge[0] > 0.8281 && edge[1] == 0.0) {
            break;
        }
    }
    coordinate[1] = fused_multiply_add(3.25, context.inverse[1], coordinate[1]);
    coordinate[1] -= context.inverse[1] * search_length([edge[1], edge[0]], 0.5, 0.5);
    coordinate[1]
}

#[inline]
fn search_length(mut edge: [f32; 2], bias: f32, scale: f32) -> f32 {
    edge[0] = bias + edge[0] * scale;
    255.0 * sample_search_nearest(edge)
}

#[inline]
fn area(distance: [f32; 2], first: f32, second: f32) -> [f32; 2] {
    let coordinate = [
        AREA_MAX_DISTANCE * glsl_round(4.0 * first) + distance[0],
        AREA_MAX_DISTANCE * glsl_round(4.0 * second) + distance[1],
    ];
    sample_area_linear([
        AREA_PIXEL_SIZE[0] * coordinate[0] + 0.5 * AREA_PIXEL_SIZE[0],
        AREA_PIXEL_SIZE[1] * coordinate[1] + 0.5 * AREA_PIXEL_SIZE[1],
    ])
}

#[inline]
fn glsl_round(value: f32) -> f32 {
    value.signum() * (value.abs() + 0.5).floor()
}

#[cfg(test)]
mod tests {
    use super::{pixel_distance, sample_level_zero_offset, SearchContext};

    #[test]
    fn pixel_distance_preserves_the_residual_before_square_root() {
        assert_eq!(0.6f32 / 0.2 - 3.0, 0.0);
        assert_eq!(pixel_distance(0.6, 0.2, 3.0), f32::EPSILON);
    }

    #[test]
    fn pinned_vector_to_scalar_offset_uses_first_component_for_both_axes() {
        let edges = [[0, 0, 0, 0], [64, 0, 0, 0], [128, 0, 0, 0], [255, 0, 0, 0]];
        let context = SearchContext {
            edges: &edges,
            width: 2,
            height: 2,
            inverse: [0.5, 0.5],
            max_steps: 16,
        };

        assert_eq!(
            sample_level_zero_offset(&context, [0.25, 0.25], [1, 0]),
            context.sample([0.75, 0.75])
        );
        assert_eq!(
            sample_level_zero_offset(&context, [0.25, 0.25], [0, 1]),
            context.sample([0.25, 0.25])
        );
    }
}
