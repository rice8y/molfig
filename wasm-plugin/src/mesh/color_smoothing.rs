// Strict CPU port of Mol* 5.9.0 mesh color smoothing and OBJ color
// quantization:
// - src/mol-geo/geometry/mesh/color-smoothing.ts
// - src/extensions/geo-export/mesh-exporter.ts

use crate::model::{Mesh, MeshMaterial, Vec3};

#[derive(Clone, Copy, Debug)]
pub(super) struct ColorSmoothingParams {
    pub(super) resolution: f64,
    pub(super) stride: usize,
    pub(super) box_min: [f64; 3],
    pub(super) box_max: [f64; 3],
    pub(super) alpha_tenths: u8,
}

pub(super) fn apply_mesh_color_smoothing(
    mesh: &mut Mesh,
    group_colors: &[u32],
    params: ColorSmoothingParams,
) {
    if mesh.vertices.is_empty()
        || mesh.faces.is_empty()
        || group_colors.is_empty()
        || !params.resolution.is_finite()
        || params.resolution <= 0.0
        || params.box_min.iter().any(|value| !value.is_finite())
        || params.box_max.iter().any(|value| !value.is_finite())
    {
        return;
    }

    let (vertex_colors, render_vertex_colors) = smooth_vertex_colors(mesh, group_colors, params);
    // The native renderer consumes the same per-vertex `vColor` values as the
    // Mol* mesh shader. Keep these before the exporter-only face selection and
    // median-cut quantization below.
    mesh.vertex_colors = render_vertex_colors;
    let mut face_colors = mesh
        .faces
        .iter()
        .map(|face| vertex_colors.get(face.a).copied().unwrap_or(0))
        .collect::<Vec<_>>();
    quantize_obj_face_colors(&mut face_colors, mesh.vertices.len());
    mesh.face_materials = face_colors
        .into_iter()
        .map(|color| MeshMaterial::with_alpha_tenths(color, params.alpha_tenths))
        .collect();
}

fn smooth_vertex_colors(
    mesh: &Mesh,
    group_colors: &[u32],
    params: ColorSmoothingParams,
) -> (Vec<u32>, Vec<[f32; 3]>) {
    let resolution = params.resolution;
    let pad = 1.0 + resolution;
    let min = [
        params.box_min[0] - pad,
        params.box_min[1] - pad,
        params.box_min[2] - pad,
    ];
    let max = [
        params.box_max[0] + pad,
        params.box_max[1] + pad,
        params.box_max[2] + pad,
    ];
    let scale_factor = 1.0 / resolution;
    let grid_dim = [
        ((max[0] - min[0]) * scale_factor).ceil() as usize + 2,
        ((max[1] - min[1]) * scale_factor).ceil() as usize + 2,
        ((max[2] - min[2]) * scale_factor).ceil() as usize + 2,
    ];
    let [xn, yn, _] = grid_dim;
    let (width, height) = volume_texture_2d_layout(grid_dim);
    let pixel_count = width.saturating_mul(height);
    let mut data = vec![0.0_f32; pixel_count.saturating_mul(3)];
    let mut count = vec![0.0_f32; pixel_count];
    let stride = params.stride.max(1);

    for vertex_index in (0..mesh.vertices.len()).step_by(stride) {
        let vertex = mesh.vertices[vertex_index];
        let v = [
            (vertex.x as f64 - min[0]) * scale_factor,
            (vertex.y as f64 - min[1]) * scale_factor,
            (vertex.z as f64 - min[2]) * scale_factor,
        ];
        let cell = [
            v[0].floor() as isize,
            v[1].floor() as isize,
            v[2].floor() as isize,
        ];
        let group = mesh.vertex_groups.get(vertex_index).copied().unwrap_or(0);
        let color = group_colors.get(group).copied().unwrap_or(0);
        let rgb = [
            ((color >> 16) & 0xff) as f64,
            ((color >> 8) & 0xff) as f64,
            (color & 0xff) as f64,
        ];

        let p = 2_isize;
        let beg = [
            (cell[0] - p).max(0) as usize,
            (cell[1] - p).max(0) as usize,
            (cell[2] - p).max(0) as usize,
        ];
        let end = [
            (cell[0] + p + 2).clamp(0, grid_dim[0] as isize) as usize,
            (cell[1] + p + 2).clamp(0, grid_dim[1] as isize) as usize,
            (cell[2] + p + 2).clamp(0, grid_dim[2] as isize) as usize,
        ];

        for x in beg[0]..end[0] {
            let dx = x as f64 - v[0];
            for y in beg[1]..end[1] {
                let dy = y as f64 - v[1];
                for z in beg[2]..end[2] {
                    let dz = z as f64 - v[2];
                    let distance = (dx * dx + dy * dy + dz * dz).sqrt();
                    if distance > 2.0 {
                        continue;
                    }
                    let weight = 2.0 - distance;
                    let index = volume_texture_index(x, y, z, xn, yn, width);
                    let color_offset = index * 3;
                    for channel in 0..3 {
                        // Float32Array compound assignment rounds after every
                        // accumulation in JavaScript.
                        data[color_offset + channel] =
                            (data[color_offset + channel] as f64 + rgb[channel] * weight) as f32;
                    }
                    count[index] = (count[index] as f64 + weight) as f32;
                }
            }
        }
    }

    let mut grid = vec![0_u8; pixel_count.saturating_mul(3)];
    for (index, &weight) in count.iter().enumerate() {
        if weight == 0.0 {
            continue;
        }
        let offset = index * 3;
        for channel in 0..3 {
            grid[offset + channel] =
                js_math_round_u8(data[offset + channel] as f64 / weight as f64);
        }
    }

    let export_colors = mesh
        .vertices
        .iter()
        .map(|&vertex| interpolate_grid_color(vertex, min, scale_factor, grid_dim, width, &grid))
        .collect();
    let render_colors = mesh
        .vertices
        .iter()
        .map(|&vertex| {
            sample_grid_color_shader(vertex, min, scale_factor, grid_dim, width, height, &grid)
        })
        .collect();
    (export_colors, render_colors)
}

/// Port of `texture3dFrom2dLinear` for the normalized uint8 color grid used by
/// Mol*'s WebGL mesh path. Unlike `getTrilinearlyInterpolated`, this preserves
/// the filtered sub-byte value visible to the vertex shader instead of writing
/// the result through a `Uint8Array` for geometry export.
fn sample_grid_color_shader(
    vertex: Vec3,
    min: [f64; 3],
    scale_factor: f64,
    grid_dim: [usize; 3],
    width: usize,
    height: usize,
    grid: &[u8],
) -> [f32; 3] {
    let min = min.map(|value| value as f32);
    let scale_factor = scale_factor as f32;
    let dim = grid_dim.map(|value| value as f32);
    let position = [
        ((vertex.x - min[0]) * scale_factor) / dim[0],
        ((vertex.y - min[1]) * scale_factor) / dim[1],
        ((vertex.z - min[2]) * scale_factor) / dim[2],
    ];
    let scaled_z = position[2] * dim[2];
    let slice0 = scaled_z.floor() as isize;
    let slice1 = slice0 + 1;
    let color0 = sample_grid_slice_shader(position, slice0, grid_dim, width, height, grid);
    let color1 = sample_grid_slice_shader(position, slice1, grid_dim, width, height, grid);
    let delta = (scaled_z - slice0 as f32).abs();
    [0, 1, 2].map(|channel| (color1[channel] - color0[channel]).mul_add(delta, color0[channel]))
}

fn sample_grid_slice_shader(
    position: [f32; 3],
    slice: isize,
    grid_dim: [usize; 3],
    width: usize,
    height: usize,
    grid: &[u8],
) -> [f32; 3] {
    let xn = grid_dim[0];
    let yn = grid_dim[1];
    let slice_offset = slice.saturating_mul(xn as isize);
    let column = slice_offset.rem_euclid(width as isize) / xn as isize;
    let row = slice_offset.div_euclid(width as isize);
    let uv = [
        (column as f32 * xn as f32 + position[0] * xn as f32) / width as f32,
        (row as f32 * yn as f32 + position[1] * yn as f32) / height as f32,
    ];
    sample_normalized_rgb8_linear(uv, width, height, grid)
}

fn sample_normalized_rgb8_linear(
    uv: [f32; 2],
    width: usize,
    height: usize,
    grid: &[u8],
) -> [f32; 3] {
    let x = uv[0] * width as f32 - 0.5;
    let y = uv[1] * height as f32 - 0.5;
    let x0 = x.floor() as isize;
    let y0 = y.floor() as isize;
    let tx = x - x0 as f32;
    let ty = y - y0 as f32;
    let texel = |sx: isize, sy: isize, channel: usize| {
        let sx = sx.clamp(0, width.saturating_sub(1) as isize) as usize;
        let sy = sy.clamp(0, height.saturating_sub(1) as isize) as usize;
        grid.get((sy * width + sx) * 3 + channel)
            .copied()
            .unwrap_or(0) as f32
    };
    [0, 1, 2].map(|channel| {
        let upper_left = texel(x0, y0, channel);
        let upper_right = texel(x0 + 1, y0, channel);
        let lower_left = texel(x0, y0 + 1, channel);
        let lower_right = texel(x0 + 1, y0 + 1, channel);
        let upper = (upper_right - upper_left).mul_add(tx, upper_left);
        let lower = (lower_right - lower_left).mul_add(tx, lower_left);
        let byte = (lower - upper).mul_add(ty, upper);
        (byte * 16.0).round() / (16.0 * 255.0)
    })
}

#[allow(clippy::needless_range_loop)]
fn interpolate_grid_color(
    vertex: Vec3,
    min: [f64; 3],
    scale_factor: f64,
    grid_dim: [usize; 3],
    width: usize,
    grid: &[u8],
) -> u32 {
    let v = [
        (vertex.x as f64 - min[0]) * scale_factor,
        (vertex.y as f64 - min[1]) * scale_factor,
        (vertex.z as f64 - min[2]) * scale_factor,
    ];
    let floor = [v[0].floor(), v[1].floor(), v[2].floor()];
    let ceil = [v[0].ceil(), v[1].ceil(), v[2].ceil()];
    let fraction = [
        js_divide(v[0] - floor[0], ceil[0] - floor[0]),
        js_divide(v[1] - floor[1], ceil[1] - floor[1]),
        js_divide(v[2] - floor[2], ceil[2] - floor[2]),
    ];
    let x0 = floor[0].max(0.0) as usize;
    let y0 = floor[1].max(0.0) as usize;
    let z0 = floor[2].max(0.0) as usize;
    let x1 = ceil[0].max(0.0) as usize;
    let y1 = ceil[1].max(0.0) as usize;
    let z1 = ceil[2].max(0.0) as usize;
    let [xn, yn, _] = grid_dim;
    let corners = [
        volume_texture_index(x0, y0, z0, xn, yn, width),
        volume_texture_index(x1, y0, z0, xn, yn, width),
        volume_texture_index(x0, y0, z1, xn, yn, width),
        volume_texture_index(x1, y0, z1, xn, yn, width),
        volume_texture_index(x0, y1, z0, xn, yn, width),
        volume_texture_index(x1, y1, z0, xn, yn, width),
        volume_texture_index(x0, y1, z1, xn, yn, width),
        volume_texture_index(x1, y1, z1, xn, yn, width),
    ];

    let mut rgb = [0_u8; 3];
    for channel in 0..3 {
        let sample = |corner: usize| grid.get(corner * 3 + channel).copied().unwrap_or(0) as f64;
        let s00 = lerp(sample(corners[0]), sample(corners[1]), fraction[0]);
        let s01 = lerp(sample(corners[2]), sample(corners[3]), fraction[0]);
        let s10 = lerp(sample(corners[4]), sample(corners[5]), fraction[0]);
        let s11 = lerp(sample(corners[6]), sample(corners[7]), fraction[0]);
        let s0 = lerp(s00, s10, fraction[1]);
        let s1 = lerp(s01, s11, fraction[1]);
        rgb[channel] = js_uint8(lerp(s0, s1, fraction[2]));
    }
    ((rgb[0] as u32) << 16) | ((rgb[1] as u32) << 8) | rgb[2] as u32
}

fn volume_texture_2d_layout(dim: [usize; 3]) -> (usize, usize) {
    let area = dim[0].saturating_mul(dim[1]).saturating_mul(dim[2]);
    let square_dim = (area as f64).sqrt();
    let power_of_two_size = 2_f64.powf(square_dim.log2().ceil()) as usize;
    let mut width = dim[0];
    let mut height = dim[1];
    if power_of_two_size < width.saturating_mul(dim[2]) {
        let columns = (power_of_two_size / width).max(1);
        let rows = dim[2].div_ceil(columns);
        width *= columns;
        height *= rows;
    } else {
        width *= dim[2];
    }
    (width, height)
}

fn volume_texture_index(x: usize, y: usize, z: usize, xn: usize, yn: usize, width: usize) -> usize {
    let column = ((z * xn) % width) / xn;
    let row = (z * xn) / width;
    row * yn * width + y * width + column * xn + x
}

fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

fn js_divide(numerator: f64, denominator: f64) -> f64 {
    if denominator == 0.0 {
        f64::NAN
    } else {
        numerator / denominator
    }
}

fn js_math_round_u8(value: f64) -> u8 {
    if !value.is_finite() || value <= 0.0 {
        0
    } else if value >= 255.0 {
        255
    } else {
        (value + 0.5).floor() as u8
    }
}

fn js_uint8(value: f64) -> u8 {
    if !value.is_finite() {
        return 0;
    }
    value.trunc().rem_euclid(256.0) as u8
}

fn quantize_obj_face_colors(face_colors: &mut [u32], vertex_count: usize) {
    if vertex_count <= 1024 {
        return;
    }
    let quantized_count = vertex_count;
    let mut prefix = Vec::with_capacity(quantized_count);
    for index in 0..quantized_count {
        prefix.push(face_colors.get(index).copied().unwrap_or(0));
    }
    let mut colors = Vec::new();
    for &color in &prefix {
        if !colors.contains(&color) {
            colors.push(color);
        }
    }
    let mut color_map = Vec::with_capacity(colors.len());
    median_cut(&mut colors, 0, 0, &mut color_map);

    for (target, source) in face_colors.iter_mut().zip(prefix) {
        if let Some((_, mapped)) = color_map.iter().rev().find(|(color, _)| *color == source) {
            *target = *mapped;
        }
    }
}

fn median_cut(colors: &mut [u32], start: usize, depth: usize, color_map: &mut Vec<(u32, u32)>) {
    if colors.is_empty() {
        return;
    }
    if colors.len() == 1 || depth >= 10 {
        let mut sum = [0.0_f64; 3];
        for &color in colors.iter() {
            sum[0] += ((color >> 16) & 0xff) as f64;
            sum[1] += ((color >> 8) & 0xff) as f64;
            sum[2] += (color & 0xff) as f64;
        }
        let count = colors.len() as f64;
        let average = ((js_math_round_u8(sum[0] / count) as u32) << 16)
            | ((js_math_round_u8(sum[1] / count) as u32) << 8)
            | js_math_round_u8(sum[2] / count) as u32;
        for &color in colors.iter() {
            color_map.push((color, average));
        }
        return;
    }

    let mut min = [255_u8; 3];
    let mut max = [0_u8; 3];
    for &color in colors.iter() {
        let rgb = [
            ((color >> 16) & 0xff) as u8,
            ((color >> 8) & 0xff) as u8,
            (color & 0xff) as u8,
        ];
        for channel in 0..3 {
            min[channel] = min[channel].min(rgb[channel]);
            max[channel] = max[channel].max(rgb[channel]);
        }
    }
    let mut channel = 0;
    if max[1] - min[1] > max[channel] - min[channel] {
        channel = 1;
    }
    if max[2] - min[2] > max[channel] - min[channel] {
        channel = 2;
    }
    molstar_sort_colors(colors, channel);
    let middle = ((start + start + colors.len() - 1) >> 1) - start;
    let (left, right) = colors.split_at_mut(middle + 1);
    median_cut(left, start, depth + 1, color_map);
    median_cut(right, start + middle + 1, depth + 1, color_map);
}

fn molstar_sort_colors(colors: &mut [u32], channel: usize) {
    if colors.len() < 2 {
        return;
    }
    quick_sort(colors, 0, colors.len() - 1, channel);
}

fn quick_sort(colors: &mut [u32], mut low: usize, mut high: usize, channel: usize) {
    while low < high {
        if high - low < 16 {
            insertion_sort(colors, low, high, channel);
            return;
        }
        let (left_end, right_end) = partition(colors, low, high, channel);
        if left_end.saturating_sub(low) < high.saturating_sub(right_end) {
            if left_end > low {
                quick_sort(colors, low, left_end - 1, channel);
            }
            low = right_end + 1;
        } else {
            if right_end < high {
                quick_sort(colors, right_end + 1, high, channel);
            }
            if left_end == 0 {
                return;
            }
            high = left_end - 1;
        }
    }
}

fn partition(colors: &mut [u32], left: usize, right: usize, channel: usize) -> (usize, usize) {
    let pivot = median_pivot_index(colors, left, right, channel);
    colors.swap(left, pivot);
    let mut equals = left + 1;
    let mut tail = right;
    while compare_color(colors[tail], colors[left], channel).is_gt() {
        tail -= 1;
    }
    let mut index = left + 1;
    while index <= tail {
        match compare_color(colors[index], colors[left], channel) {
            std::cmp::Ordering::Greater => {
                colors.swap(index, tail);
                tail -= 1;
                while compare_color(colors[tail], colors[left], channel).is_gt() {
                    tail -= 1;
                }
            }
            std::cmp::Ordering::Equal => {
                colors.swap(index, equals);
                equals += 1;
                index += 1;
            }
            std::cmp::Ordering::Less => index += 1,
        }
    }
    for index in left..equals {
        colors.swap(index, left + tail - index);
    }
    (tail - equals + left + 1, tail)
}

fn median_pivot_index(colors: &[u32], left: usize, right: usize, channel: usize) -> usize {
    let middle = (left + right) >> 1;
    if compare_color(colors[left], colors[right], channel).is_gt() {
        if compare_color(colors[left], colors[middle], channel).is_gt() {
            if compare_color(colors[middle], colors[right], channel).is_gt() {
                middle
            } else {
                right
            }
        } else {
            left
        }
    } else if compare_color(colors[right], colors[middle], channel).is_gt() {
        if compare_color(colors[middle], colors[left], channel).is_gt() {
            middle
        } else {
            left
        }
    } else {
        right
    }
}

fn insertion_sort(colors: &mut [u32], start: usize, end: usize, channel: usize) {
    for index in start + 1..=end {
        let mut cursor = index;
        while cursor > start && compare_color(colors[cursor - 1], colors[cursor], channel).is_gt() {
            colors.swap(cursor - 1, cursor);
            cursor -= 1;
        }
    }
}

fn compare_color(a: u32, b: u32, channel: usize) -> std::cmp::Ordering {
    color_channel(a, channel).cmp(&color_channel(b, channel))
}

fn color_channel(color: u32, channel: usize) -> u8 {
    match channel {
        0 => ((color >> 16) & 0xff) as u8,
        1 => ((color >> 8) & 0xff) as u8,
        _ => (color & 0xff) as u8,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Face;

    #[test]
    fn smoothing_blends_group_colors_at_shared_surface_vertices() {
        let mut mesh = Mesh {
            vertices: vec![
                Vec3::new(-0.25, 0.0, 0.0),
                Vec3::new(0.25, 0.0, 0.0),
                Vec3::new(0.0, 0.25, 0.0),
            ],
            normals: vec![Vec3::default(); 3],
            faces: vec![Face { a: 2, b: 0, c: 1 }],
            vertex_groups: vec![0, 1, 0],
            vertex_colors: Vec::new(),
            face_groups: vec![0],
            face_materials: Vec::new(),
            sections: Vec::new(),
            group_count: 2,
        };
        apply_mesh_color_smoothing(
            &mut mesh,
            &[0xff0000, 0x0000ff],
            ColorSmoothingParams {
                resolution: 0.5,
                stride: 1,
                box_min: [-1.0; 3],
                box_max: [1.0; 3],
                alpha_tenths: 10,
            },
        );
        let color = mesh.face_materials[0].color;
        assert_eq!(mesh.vertex_colors.len(), mesh.vertices.len());
        assert!(mesh
            .vertex_colors
            .iter()
            .flatten()
            .all(|channel| channel.is_finite() && (0.0..=1.0).contains(channel)));
        assert_ne!(color, 0xff0000);
        assert_ne!(color, 0x0000ff);
        assert_eq!(mesh.face_materials[0].alpha_tenths, 10);
    }

    #[test]
    fn obj_quantization_leaves_small_mesh_colors_unchanged() {
        let mut colors = vec![0x112233, 0x445566];
        quantize_obj_face_colors(&mut colors, 1024);
        assert_eq!(colors, vec![0x112233, 0x445566]);
    }

    #[test]
    fn normalized_rgb8_sampler_keeps_sub_byte_vertex_shader_values() {
        let grid = [
            0, 0, 0, // upper left
            255, 0, 0, // upper right
            0, 255, 0, // lower left
            0, 0, 255, // lower right
        ];
        assert_eq!(
            sample_normalized_rgb8_linear([0.5, 0.5], 2, 2, &grid),
            [0.25; 3]
        );
    }
}
