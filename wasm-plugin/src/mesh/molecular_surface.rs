// Molecular surface field generation is a strict CPU port of Mol* 5.9.0:
// - src/mol-math/geometry/molecular-surface.ts
// - src/mol-math/geometry/lookup3d/grid.ts
//
// Copyright (c) 2018-2025 Mol* contributors.
// Licensed under the MIT License. See the Mol* LICENSE file for details.

use crate::model::{Mesh, Vec3};

use super::surface::{marching_cubes_mesh, GaussianDensityGrid};

#[derive(Clone, Copy, Debug)]
pub(super) struct MolecularSurfacePoint {
    pub(super) position: Vec3,
    pub(super) radius: f64,
    pub(super) group_id: usize,
}

impl MolecularSurfacePoint {
    pub(super) const fn new(position: Vec3, radius: f64, group_id: usize) -> Self {
        Self {
            position,
            radius,
            group_id,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct MolecularSurfaceParams {
    pub(super) resolution: f64,
    pub(super) probe_radius: f64,
    pub(super) probe_positions: usize,
}

impl Default for MolecularSurfaceParams {
    fn default() -> Self {
        Self {
            resolution: 0.5,
            probe_radius: 1.4,
            probe_positions: 36,
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct MolecularSurfaceGrid {
    pub(super) grid: GaussianDensityGrid,
    pub(super) max_radius: f64,
}

pub(super) fn build_molecular_surface_mesh_in_box(
    points: &[MolecularSurfacePoint],
    params: MolecularSurfaceParams,
    box_min: Vec3,
    box_max: Vec3,
) -> Mesh {
    build_molecular_surface_mesh_in_box64(
        points,
        params,
        [box_min.x as f64, box_min.y as f64, box_min.z as f64],
        [box_max.x as f64, box_max.y as f64, box_max.z as f64],
    )
}

pub(super) fn build_molecular_surface_mesh_in_box64(
    points: &[MolecularSurfacePoint],
    params: MolecularSurfaceParams,
    box_min: [f64; 3],
    box_max: [f64; 3],
) -> Mesh {
    let result = build_molecular_surface_grid_in_box64(points, params, box_min, box_max, false);
    debug_assert!(result.max_radius.is_finite());
    marching_cubes_mesh(&result.grid, params.probe_radius)
}

pub(super) fn build_structure_molecular_surface_mesh_in_box64(
    points: &[MolecularSurfacePoint],
    params: MolecularSurfaceParams,
    box_min: [f64; 3],
    box_max: [f64; 3],
) -> Mesh {
    let result = build_molecular_surface_grid_in_box64(points, params, box_min, box_max, true);
    debug_assert!(result.max_radius.is_finite());
    marching_cubes_mesh(&result.grid, params.probe_radius)
}

#[cfg(test)]
pub(super) fn build_molecular_surface_grid_in_box(
    points: &[MolecularSurfacePoint],
    params: MolecularSurfaceParams,
    box_min: Vec3,
    box_max: Vec3,
) -> MolecularSurfaceGrid {
    build_molecular_surface_grid_in_box64(
        points,
        params,
        [box_min.x as f64, box_min.y as f64, box_min.z as f64],
        [box_max.x as f64, box_max.y as f64, box_max.z as f64],
        false,
    )
}

pub(super) fn build_molecular_surface_grid_in_box64(
    points: &[MolecularSurfacePoint],
    params: MolecularSurfaceParams,
    box_min: [f64; 3],
    box_max: [f64; 3],
    include_structure_sentinel: bool,
) -> MolecularSurfaceGrid {
    let resolution = params.resolution;
    let probe_radius = params.probe_radius;
    assert!(
        resolution.is_finite() && resolution > 0.0,
        "Molecular surface resolution must be finite and positive"
    );
    assert!(
        probe_radius.is_finite() && probe_radius >= 0.0,
        "Molecular surface probe radius must be finite and non-negative"
    );
    assert!(
        params.probe_positions > 0,
        "Molecular surface probe position count must be positive"
    );
    if points.is_empty() {
        return MolecularSurfaceGrid {
            grid: GaussianDensityGrid::empty(resolution),
            max_radius: 0.0,
        };
    }

    let mut positions = Vec::with_capacity(points.len());
    let mut radii = Vec::with_capacity(points.len());
    let mut ids = Vec::with_capacity(points.len());
    let mut max_radius = 0.0_f64;
    for point in points {
        assert!(
            point.position.is_finite() && point.radius.is_finite() && point.radius >= 0.0,
            "Molecular surface input must be finite"
        );
        positions.push(point.position);
        max_radius = max_radius.max(point.radius);
        // Mol* writes the physical radius plus probe radius to Float32Array.
        radii.push((point.radius + probe_radius) as f32);
        ids.push(point.group_id as f32);
    }

    let pad = max_radius + resolution;
    let origin = [box_min[0] - pad, box_min[1] - pad, box_min[2] - pad];
    let expanded_max = [box_max[0] + pad, box_max[1] + pad, box_max[2] + pad];
    let scale_factor = 1.0 / resolution;
    let dimensions = [
        (expanded_max[0] * scale_factor - origin[0] * scale_factor).ceil() as usize,
        (expanded_max[1] * scale_factor - origin[1] * scale_factor).ceil() as usize,
        (expanded_max[2] * scale_factor - origin[2] * scale_factor).ceil() as usize,
    ];
    let cell_count = dimensions[0]
        .checked_mul(dimensions[1])
        .and_then(|value| value.checked_mul(dimensions[2]))
        .expect("Molecular surface grid is too large");

    let grid_x = fill_grid_dimension(dimensions[0], origin[0], resolution);
    let grid_y = fill_grid_dimension(dimensions[1], origin[1], resolution);
    let grid_z = fill_grid_dimension(dimensions[2], origin[2], resolution);
    let mut field = vec![-1001.0_f32; cell_count];
    let mut id_field = vec![-1.0_f32; cell_count];

    let lookup = RadiusLookupGrid::new(
        &positions,
        &radii,
        box_min,
        box_max,
        max_radius * 2.0,
        include_structure_sentinel,
    );
    let angle_tables = angle_tables(params.probe_positions);
    let mut state = MolecularSurfaceState {
        positions: &positions,
        radii: &radii,
        ids: &ids,
        lookup: &lookup,
        dimensions,
        origin,
        scale_factor,
        resolution,
        probe_radius,
        probe_positions: params.probe_positions,
        cos_table: &angle_tables.0,
        sin_table: &angle_tables.1,
        grid_x: &grid_x,
        grid_y: &grid_y,
        grid_z: &grid_z,
        field: &mut field,
        id_field: &mut id_field,
        neighbours: Vec::new(),
        last_clip: None,
    };
    state.project_points();
    state.project_torii();

    MolecularSurfaceGrid {
        grid: GaussianDensityGrid {
            dimensions,
            field,
            id_field,
            origin,
            resolution,
        },
        max_radius,
    }
}

struct MolecularSurfaceState<'a> {
    positions: &'a [Vec3],
    radii: &'a [f32],
    ids: &'a [f32],
    lookup: &'a RadiusLookupGrid,
    dimensions: [usize; 3],
    origin: [f64; 3],
    scale_factor: f64,
    resolution: f64,
    probe_radius: f64,
    probe_positions: usize,
    cos_table: &'a [f32],
    sin_table: &'a [f32],
    grid_x: &'a [f32],
    grid_y: &'a [f32],
    grid_z: &'a [f32],
    field: &'a mut [f32],
    id_field: &'a mut [f32],
    neighbours: Vec<usize>,
    last_clip: Option<usize>,
}

impl MolecularSurfaceState<'_> {
    fn obscured(&mut self, x: f64, y: f64, z: f64, a: usize, b: Option<usize>) -> Option<usize> {
        if let Some(index) = self.last_clip {
            if index != a && Some(index) != b && self.single_atom_obscures(index, x, y, z) {
                return Some(index);
            }
            self.last_clip = None;
        }

        for position in 0..self.neighbours.len() {
            let index = self.neighbours[position];
            if index != a && Some(index) != b && self.single_atom_obscures(index, x, y, z) {
                self.last_clip = Some(index);
                return Some(index);
            }
        }
        None
    }

    fn single_atom_obscures(&self, index: usize, x: f64, y: f64, z: f64) -> bool {
        let radius = self.radii[index] as f64;
        let point = self.positions[index];
        let dx = point.x as f64 - x;
        let dy = point.y as f64 - y;
        let dz = point.z as f64 - z;
        dx * dx + dy * dy + dz * dz < radius * radius
    }

    fn project_points(&mut self) {
        let [dim_x, dim_y, dim_z] = self.dimensions;
        let yz = dim_y * dim_z;
        let grid_extension = usize::from(self.probe_radius < self.resolution * 2.0) as isize;

        for point_index in 0..self.positions.len() {
            let point = self.positions[point_index];
            let vx = point.x as f64;
            let vy = point.y as f64;
            let vz = point.z as f64;
            let radius = self.radii[point_index] as f64;
            let radius_sq = radius * radius;
            let extended = grid_extension > 0;
            self.neighbours = self.lookup.find(vx, vy, vz, radius);

            let grid_radius = (radius * self.scale_factor).ceil() as isize + grid_extension;
            let atom_x = (self.scale_factor * (vx - self.origin[0])).floor() as isize;
            let atom_y = (self.scale_factor * (vy - self.origin[1])).floor() as isize;
            let atom_z = (self.scale_factor * (vz - self.origin[2])).floor() as isize;
            let begin_x = (atom_x - grid_radius).max(0) as usize;
            let begin_y = (atom_y - grid_radius).max(0) as usize;
            let begin_z = (atom_z - grid_radius).max(0) as usize;
            let end_x = (atom_x + grid_radius + 2).min(dim_x as isize) as usize;
            let end_y = (atom_y + grid_radius + 2).min(dim_y as isize) as usize;
            let end_z = (atom_z + grid_radius + 2).min(dim_z as isize) as usize;

            for x in begin_x..end_x {
                let dx = self.grid_x[x] as f64 - vx;
                let x_offset = x * yz;
                for y in begin_y..end_y {
                    let dy = self.grid_y[y] as f64 - vy;
                    let distance_xy_sq = dx * dx + dy * dy;
                    let xy_offset = y * dim_z + x_offset;
                    for z in begin_z..end_z {
                        let dz = self.grid_z[z] as f64 - vz;
                        let distance_sq = distance_xy_sq + dz * dz;
                        if !extended && distance_sq >= radius_sq {
                            continue;
                        }

                        let distance = distance_sq.sqrt();
                        let projection = radius / distance;
                        let surface_x = dx * projection + vx;
                        let surface_y = dy * projection + vy;
                        let surface_z = dz * projection + vz;
                        let obscured = self
                            .obscured(surface_x, surface_y, surface_z, point_index, None)
                            .is_some();
                        let offset = z + xy_offset;

                        if distance_sq < radius_sq {
                            if self.field[offset] < 0.0 {
                                self.field[offset] *= -1.0;
                            }
                            if !obscured {
                                let delta = radius - distance;
                                if delta < self.field[offset] as f64 {
                                    self.field[offset] = delta as f32;
                                    self.id_field[offset] = self.ids[point_index];
                                }
                            }
                        } else if extended && !obscured {
                            let delta = radius - distance;
                            if delta > self.field[offset] as f64 {
                                self.field[offset] = delta as f32;
                            }
                        }
                    }
                }
            }
        }
    }

    fn project_torii(&mut self) {
        for a in 0..self.positions.len() {
            let point = self.positions[a];
            self.neighbours = self.lookup.find(
                point.x as f64,
                point.y as f64,
                point.z as f64,
                self.radii[a] as f64,
            );
            let neighbours = self.neighbours.clone();
            for b in neighbours {
                if a < b {
                    self.project_torus(a, b);
                }
            }
        }
    }

    fn project_torus(&mut self, a: usize, b: usize) {
        let radius_a = self.radii[a] as f64;
        let radius_b = self.radii[b] as f64;
        let pa = self.positions[a];
        let pb = self.positions[b];
        let mut atob = [
            pb.x as f64 - pa.x as f64,
            pb.y as f64 - pa.y as f64,
            pb.z as f64 - pa.z as f64,
        ];
        let distance_sq = dot(atob, atob);
        let distance = distance_sq.sqrt();
        let cosine_a = (radius_a * radius_a + distance * distance - radius_b * radius_b)
            / (2.0 * radius_a * distance);
        let midpoint_distance = radius_a * cosine_a;
        normalize(&mut atob);

        let mut normal_a = normal_to_line(atob);
        normalize(&mut normal_a);
        let mut normal_b = cross(atob, normal_a);
        normalize(&mut normal_b);
        let intersection_radius =
            (radius_a * radius_a - midpoint_distance * midpoint_distance).sqrt();
        scale(&mut normal_a, intersection_radius);
        scale(&mut normal_b, intersection_radius);
        scale(&mut atob, midpoint_distance);
        let midpoint = [
            atob[0] + pa.x as f64,
            atob[1] + pa.y as f64,
            atob[2] + pa.z as f64,
        ];
        self.last_clip = None;

        let grid_radius = 2 + (self.probe_radius * self.scale_factor).floor() as isize;
        let [dim_x, dim_y, dim_z] = self.dimensions;
        let yz = dim_y * dim_z;

        for index in 0..self.probe_positions {
            let cosine = self.cos_table[index] as f64;
            let sine = self.sin_table[index] as f64;
            let probe = [
                midpoint[0] + cosine * normal_a[0] + sine * normal_b[0],
                midpoint[1] + cosine * normal_a[1] + sine * normal_b[1],
                midpoint[2] + cosine * normal_a[2] + sine * normal_b[2],
            ];
            if self
                .obscured(probe[0], probe[1], probe[2], a, Some(b))
                .is_some()
            {
                continue;
            }

            let atom_x = (self.scale_factor * (probe[0] - self.origin[0])).floor() as isize;
            let atom_y = (self.scale_factor * (probe[1] - self.origin[1])).floor() as isize;
            let atom_z = (self.scale_factor * (probe[2] - self.origin[2])).floor() as isize;
            let begin_x = (atom_x - grid_radius).max(0) as usize;
            let begin_y = (atom_y - grid_radius).max(0) as usize;
            let begin_z = (atom_z - grid_radius).max(0) as usize;
            let end_x = (atom_x + grid_radius + 2).min(dim_x as isize) as usize;
            let end_y = (atom_y + grid_radius + 2).min(dim_y as isize) as usize;
            let end_z = (atom_z + grid_radius + 2).min(dim_z as isize) as usize;

            for x in begin_x..end_x {
                let dx = probe[0] - self.grid_x[x] as f64;
                let x_offset = x * yz;
                for y in begin_y..end_y {
                    let dy = probe[1] - self.grid_y[y] as f64;
                    let distance_xy_sq = dx * dx + dy * dy;
                    let xy_offset = y * dim_z + x_offset;
                    for z in begin_z..end_z {
                        let dz = probe[2] - self.grid_z[z] as f64;
                        let distance_sq = distance_xy_sq + dz * dz;
                        let offset = z + xy_offset;
                        let current = self.field[offset] as f64;
                        if current > 0.0 && distance_sq < current * current {
                            self.field[offset] = distance_sq.sqrt() as f32;
                            let projection = dx * atob[0] + dy * atob[1] + dz * atob[2];
                            self.id_field[offset] = self.ids[if projection < 0.0 { b } else { a }];
                        }
                    }
                }
            }
        }
    }
}

#[derive(Clone, Debug)]
struct RadiusLookupGrid {
    positions: Vec<Vec3>,
    radii: Vec<f32>,
    size: [usize; 3],
    min: [f64; 3],
    delta: [f64; 3],
    grid: Vec<u32>,
    bucket_offsets: Vec<usize>,
    bucket_counts: Vec<usize>,
    bucket_array: Vec<usize>,
    max_radius: f64,
}

impl RadiusLookupGrid {
    fn new(
        positions: &[Vec3],
        radii: &[f32],
        box_min: [f64; 3],
        box_max: [f64; 3],
        cell_size: f64,
        include_structure_sentinel: bool,
    ) -> Self {
        let min = [box_min[0] - 0.5, box_min[1] - 0.5, box_min[2] - 0.5];
        let max = [box_max[0] + 0.5, box_max[1] + 0.5, box_max[2] + 0.5];
        let extent = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];
        let (mut size, mut delta) = if cell_size != 0.0 {
            (
                [
                    (extent[0] / cell_size).ceil().max(1.0) as usize,
                    (extent[1] / cell_size).ceil().max(1.0) as usize,
                    (extent[2] / cell_size).ceil().max(1.0) as usize,
                ],
                [cell_size; 3],
            )
        } else {
            ([1, 1, 1], extent)
        };
        let volume = size[0] * size[1] * size[2];
        if volume > (1 << 24) {
            let factor = (volume as f64 / (1 << 24) as f64).cbrt();
            size = [
                (size[0] as f64 / factor).ceil() as usize,
                (size[1] as f64 / factor).ceil() as usize,
                (size[2] as f64 / factor).ceil() as usize,
            ];
            delta = [
                extent[0] / size[0] as f64,
                extent[1] / size[1] as f64,
                extent[2] / size[2] as f64,
            ];
        }

        let mut grid = vec![0_u32; size[0] * size[1] * size[2]];
        let mut bucket_indices = Vec::with_capacity(positions.len());
        let mut bucket_count = 0_usize;
        for position in positions {
            let x = (((position.x as f64 - min[0]) / delta[0]).floor() as isize)
                .clamp(0, size[0] as isize - 1) as usize;
            let y = (((position.y as f64 - min[1]) / delta[1]).floor() as isize)
                .clamp(0, size[1] as isize - 1) as usize;
            let z = (((position.z as f64 - min[2]) / delta[2]).floor() as isize)
                .clamp(0, size[2] as isize - 1) as usize;
            let offset = ((x * size[1]) + y) * size[2] + z;
            grid[offset] += 1;
            if grid[offset] == 1 {
                bucket_count += 1;
            }
            bucket_indices.push(offset);
        }
        if include_structure_sentinel {
            // OrderedSet.ofRange(0, id.length) includes its end. The
            // undefined position becomes bucket index zero in Int32Array.
            bucket_indices.push(0);
        }

        let mut bucket_counts = vec![0_usize; bucket_count];
        let mut next = 0_usize;
        for cell in &mut grid {
            if *cell > 0 {
                let count = *cell as usize;
                *cell = (next + 1) as u32;
                bucket_counts[next] = count;
                next += 1;
            }
        }
        let mut bucket_offsets = vec![0_usize; bucket_count];
        for index in 1..bucket_count {
            bucket_offsets[index] = bucket_offsets[index - 1] + bucket_counts[index - 1];
        }
        let mut bucket_fill = vec![0_usize; bucket_count];
        let mut bucket_array = vec![0_usize; bucket_indices.len()];
        for (position_index, grid_index) in bucket_indices.into_iter().enumerate() {
            let bucket = grid[grid_index] as usize - 1;
            let destination = bucket_offsets[bucket] + bucket_fill[bucket];
            bucket_array[destination] = position_index;
            bucket_fill[bucket] += 1;
        }

        Self {
            positions: positions.to_vec(),
            radii: radii.to_vec(),
            size,
            min,
            delta,
            grid,
            bucket_offsets,
            bucket_counts,
            bucket_array,
            max_radius: radii.iter().copied().fold(0.0_f32, f32::max) as f64,
        }
    }

    fn find(&self, x: f64, y: f64, z: f64, input_radius: f64) -> Vec<usize> {
        let radius = input_radius + self.max_radius;
        let radius_sq = radius * radius;
        let low = [
            (((x - radius - self.min[0]) / self.delta[0]).floor() as isize).max(0),
            (((y - radius - self.min[1]) / self.delta[1]).floor() as isize).max(0),
            (((z - radius - self.min[2]) / self.delta[2]).floor() as isize).max(0),
        ];
        let high = [
            (((x + radius - self.min[0]) / self.delta[0]).floor() as isize)
                .min(self.size[0] as isize - 1),
            (((y + radius - self.min[1]) / self.delta[1]).floor() as isize)
                .min(self.size[1] as isize - 1),
            (((z + radius - self.min[2]) / self.delta[2]).floor() as isize)
                .min(self.size[2] as isize - 1),
        ];
        if low[0] > high[0] || low[1] > high[1] || low[2] > high[2] {
            return Vec::new();
        }

        let mut result = Vec::new();
        for grid_x in low[0]..=high[0] {
            for grid_y in low[1]..=high[1] {
                for grid_z in low[2]..=high[2] {
                    let grid_index = ((grid_x as usize * self.size[1]) + grid_y as usize)
                        * self.size[2]
                        + grid_z as usize;
                    let bucket = self.grid[grid_index];
                    if bucket == 0 {
                        continue;
                    }
                    let bucket = bucket as usize - 1;
                    let begin = self.bucket_offsets[bucket];
                    let end = begin + self.bucket_counts[bucket];
                    for position in begin..end {
                        let index = self.bucket_array[position];
                        if index >= self.positions.len() {
                            continue;
                        }
                        let point = self.positions[index];
                        let dx = point.x as f64 - x;
                        let dy = point.y as f64 - y;
                        let dz = point.z as f64 - z;
                        let distance_sq = dx * dx + dy * dy + dz * dz;
                        if distance_sq <= radius_sq
                            && distance_sq.sqrt() - self.radii[index] as f64 <= input_radius
                        {
                            result.push(index);
                        }
                    }
                }
            }
        }
        result
    }
}

fn fill_grid_dimension(length: usize, start: f64, step: f64) -> Vec<f32> {
    (0..length)
        .map(|index| (start + step * index as f64) as f32)
        .collect()
}

fn angle_tables(count: usize) -> (Vec<f32>, Vec<f32>) {
    let step = std::f64::consts::TAU / count as f64;
    let mut theta = 0.0_f64;
    let mut cosine = Vec::with_capacity(count);
    let mut sine = Vec::with_capacity(count);
    for _ in 0..count {
        cosine.push(theta.cos() as f32);
        sine.push(theta.sin() as f32);
        theta += step;
    }
    (cosine, sine)
}

#[cfg(test)]
pub(super) fn molecular_surface_lookup_contract(
    points: &[MolecularSurfacePoint],
    probe_radius: f64,
    box_min: [f64; 3],
    box_max: [f64; 3],
) -> ([usize; 3], [f64; 3], [f64; 3], [u64; 3]) {
    let positions = points
        .iter()
        .map(|point| point.position)
        .collect::<Vec<_>>();
    let radii = points
        .iter()
        .map(|point| (point.radius + probe_radius) as f32)
        .collect::<Vec<_>>();
    let max_radius = points
        .iter()
        .map(|point| point.radius)
        .fold(0.0_f64, f64::max);
    let lookup =
        RadiusLookupGrid::new(&positions, &radii, box_min, box_max, max_radius * 2.0, true);
    let hash = |values: &[usize]| {
        values
            .iter()
            .fold(0xcbf29ce484222325_u64, |mut hash, &value| {
                for byte in (value as i32).to_le_bytes() {
                    hash ^= byte as u64;
                    hash = hash.wrapping_mul(0x100000001b3);
                }
                hash
            })
    };
    (
        lookup.size,
        lookup.min,
        lookup.delta,
        [
            hash(&lookup.bucket_offsets),
            hash(&lookup.bucket_counts),
            hash(&lookup.bucket_array),
        ],
    )
}

fn normal_to_line(point: [f64; 3]) -> [f64; 3] {
    let mut out = [1.0_f64; 3];
    if point[0] != 0.0 {
        out[0] = (point[1] + point[2]) / -point[0];
    } else if point[1] != 0.0 {
        out[1] = (point[0] + point[2]) / -point[1];
    } else if point[2] != 0.0 {
        out[2] = (point[0] + point[1]) / -point[2];
    }
    out
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(value: &mut [f64; 3]) {
    let length = dot(*value, *value).sqrt();
    let scale = if length > 0.0 { 1.0 / length } else { 0.0 };
    value[0] *= scale;
    value[1] *= scale;
    value[2] *= scale;
}

fn scale(value: &mut [f64; 3], factor: f64) {
    value[0] *= factor;
    value[1] *= factor;
    value[2] *= factor;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(x: f32, y: f32, z: f32, radius: f64, group_id: usize) -> MolecularSurfacePoint {
        MolecularSurfacePoint::new(Vec3::new(x, y, z), radius, group_id)
    }

    #[test]
    fn default_parameters_match_molstar() {
        let params = MolecularSurfaceParams::default();
        assert_eq!(params.resolution, 0.5);
        assert_eq!(params.probe_radius, 1.4);
        assert_eq!(params.probe_positions, 36);
    }

    #[test]
    fn radius_lookup_uses_molstar_bucket_and_radius_order() {
        let positions = [
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let radii = [2.0, 2.0, 2.0];
        let lookup = RadiusLookupGrid::new(
            &positions,
            &radii,
            [-1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            3.0,
            false,
        );
        assert_eq!(lookup.find(-1.0, 0.0, 0.0, 2.0), vec![0, 1, 2]);
    }

    #[test]
    fn molecular_surface_grid_and_mesh_are_deterministic() {
        let points = [point(-1.0, 0.0, 0.0, 1.5, 3), point(1.0, 0.0, 0.0, 1.5, 7)];
        let params = MolecularSurfaceParams {
            resolution: 0.5,
            ..MolecularSurfaceParams::default()
        };
        let min = Vec3::new(-1.0, 0.0, 0.0);
        let max = Vec3::new(1.0, 0.0, 0.0);
        let first = build_molecular_surface_grid_in_box(&points, params, min, max);
        let second = build_molecular_surface_grid_in_box(&points, params, min, max);
        assert_eq!(first.grid.dimensions, second.grid.dimensions);
        assert_eq!(first.grid.field, second.grid.field);
        assert_eq!(first.grid.id_field, second.grid.id_field);
        assert_eq!(first.max_radius, 1.5);

        let mesh = build_molecular_surface_mesh_in_box(&points, params, min, max);
        assert!(!mesh.vertices.is_empty());
        assert!(!mesh.faces.is_empty());
        assert!(mesh.vertices.iter().all(|vertex| vertex.is_finite()));
        assert!(mesh
            .vertex_groups
            .iter()
            .all(|&group| group == 3 || group == 7));
    }
}
