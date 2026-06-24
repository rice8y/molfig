// Gaussian density and marching cubes are strict Rust ports of Mol* 5.9.0:
// - src/mol-math/geometry/gaussian-density/cpu.ts
// - src/mol-geo/util/marching-cubes/algorithm.ts
// - src/mol-geo/util/marching-cubes/builder.ts
//
// Copyright (c) 2018-2022 Mol* contributors.
// Licensed under the MIT License. See the Mol* LICENSE file for details.

use crate::model::{Face, Mesh, Vec3};

use super::surface_tables::{CUBE_EDGES, EDGE_ID_INFO, EDGE_TABLE, TRI_TABLE};

#[derive(Clone, Copy, Debug)]
pub(super) struct GaussianDensityPoint {
    pub(super) position: Vec3,
    pub(super) radius: f64,
    pub(super) group_id: usize,
}

impl GaussianDensityPoint {
    pub(super) const fn new(position: Vec3, radius: f64, group_id: usize) -> Self {
        Self {
            position,
            radius,
            group_id,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub(super) struct GaussianDensityParams {
    pub(super) resolution: f32,
    pub(super) radius_offset: f32,
    pub(super) smoothness: f32,
}

impl Default for GaussianDensityParams {
    fn default() -> Self {
        Self {
            resolution: 1.0,
            radius_offset: 0.0,
            smoothness: 1.5,
        }
    }
}

#[derive(Clone, Debug)]
pub(super) struct GaussianDensityGrid {
    pub(super) dimensions: [usize; 3],
    pub(super) field: Vec<f32>,
    pub(super) id_field: Vec<f32>,
    pub(super) origin: [f64; 3],
    pub(super) resolution: f64,
}

impl GaussianDensityGrid {
    pub(super) fn empty(resolution: f64) -> Self {
        Self {
            dimensions: [0, 0, 0],
            field: Vec::new(),
            id_field: Vec::new(),
            origin: [0.0; 3],
            resolution,
        }
    }

    #[inline]
    fn offset(&self, x: usize, y: usize, z: usize) -> usize {
        z + self.dimensions[2] * (y + self.dimensions[1] * x)
    }

    #[inline]
    fn scalar(&self, x: usize, y: usize, z: usize) -> f32 {
        self.field[self.offset(x, y, z)]
    }

    #[inline]
    fn group(&self, x: usize, y: usize, z: usize) -> f32 {
        self.id_field[self.offset(x, y, z)]
    }
}

#[cfg(test)]
pub(super) fn build_gaussian_density_grid(
    points: &[GaussianDensityPoint],
    params: GaussianDensityParams,
) -> GaussianDensityGrid {
    let boundary = points.first().map(|first| {
        points[1..]
            .iter()
            .fold((first.position, first.position), |(min, max), point| {
                (min.min(point.position), max.max(point.position))
            })
    });
    build_gaussian_density_grid_with_boundary(points, params, boundary)
}

pub(super) fn build_gaussian_density_grid_in_box(
    points: &[GaussianDensityPoint],
    params: GaussianDensityParams,
    box_min: Vec3,
    box_max: Vec3,
) -> GaussianDensityGrid {
    build_gaussian_density_grid_with_boundary(points, params, Some((box_min, box_max)))
}

fn build_gaussian_density_grid_with_boundary(
    points: &[GaussianDensityPoint],
    params: GaussianDensityParams,
    boundary: Option<(Vec3, Vec3)>,
) -> GaussianDensityGrid {
    let resolution = params.resolution as f64;
    assert!(
        resolution.is_finite() && resolution > 0.0,
        "Gaussian density resolution must be finite and positive"
    );
    assert!(
        params.radius_offset.is_finite() && params.smoothness.is_finite(),
        "Gaussian density parameters must be finite"
    );
    if points.is_empty() {
        return GaussianDensityGrid::empty(resolution);
    }

    let (box_min, box_max) = boundary.expect("non-empty Gaussian density input needs a box");
    let min = [box_min.x as f64, box_min.y as f64, box_min.z as f64];
    let max = [box_max.x as f64, box_max.y as f64, box_max.z as f64];
    let mut radii = Vec::with_capacity(points.len());
    let mut max_radius = 0.0_f64;

    for point in points {
        assert!(
            point.position.is_finite() && point.radius.is_finite(),
            "Gaussian density input must be finite"
        );
        // Mol* keeps maxRadius as a JS number, but writes every per-point
        // radius to Float32Array before the accumulation loop.
        let radius = point.radius + params.radius_offset as f64;
        assert!(
            radius.is_finite() && radius > 0.0,
            "Gaussian density radii plus radius_offset must be positive"
        );
        max_radius = max_radius.max(radius);
        radii.push(radius as f32);
    }

    let pad = max_radius * 2.0 + resolution;
    let origin = [min[0] - pad, min[1] - pad, min[2] - pad];
    let expanded_max = [max[0] + pad, max[1] + pad, max[2] + pad];
    let scale_factor = 1.0 / resolution;
    let dimensions = [
        (expanded_max[0] * scale_factor - origin[0] * scale_factor).ceil() as usize,
        (expanded_max[1] * scale_factor - origin[1] * scale_factor).ceil() as usize,
        (expanded_max[2] * scale_factor - origin[2] * scale_factor).ceil() as usize,
    ];
    let cell_count = dimensions[0]
        .checked_mul(dimensions[1])
        .and_then(|value| value.checked_mul(dimensions[2]))
        .expect("Gaussian density grid is too large");

    let grid_x = fill_grid_dimension(dimensions[0], origin[0], resolution);
    let grid_y = fill_grid_dimension(dimensions[1], origin[1], resolution);
    let grid_z = fill_grid_dimension(dimensions[2], origin[2], resolution);
    let mut field = vec![0.0_f32; cell_count];
    let mut id_field = vec![-1.0_f32; cell_count];
    let mut dominant_density = vec![0.0_f32; cell_count];
    let alpha = params.smoothness as f64;
    let dim_x = dimensions[0];
    let dim_y = dimensions[1];
    let dim_z = dimensions[2];
    let yz = dim_z * dim_y;

    for (point_index, point) in points.iter().enumerate() {
        let vx = point.position.x as f64;
        let vy = point.position.y as f64;
        let vz = point.position.z as f64;
        let radius = radii[point_index] as f64;
        let radius_sq = radius * radius;
        let radius_sq_inv = 1.0 / radius_sq;
        let radius2 = radius * 2.0;
        let radius2_sq = radius2 * radius2;
        let grid_radius = (radius2 * scale_factor).ceil() as isize;
        let atom_x = (scale_factor * (vx - origin[0])).floor() as isize;
        let atom_y = (scale_factor * (vy - origin[1])).floor() as isize;
        let atom_z = (scale_factor * (vz - origin[2])).floor() as isize;

        let begin_x = (atom_x - grid_radius).max(0) as usize;
        let begin_y = (atom_y - grid_radius).max(0) as usize;
        let begin_z = (atom_z - grid_radius).max(0) as usize;
        let end_x = (atom_x + grid_radius + 2).min(dim_x as isize) as usize;
        let end_y = (atom_y + grid_radius + 2).min(dim_y as isize) as usize;
        let end_z = (atom_z + grid_radius + 2).min(dim_z as isize) as usize;

        for x in begin_x..end_x {
            let dx = grid_x[x] as f64 - vx;
            let x_offset = x * yz;
            for y in begin_y..end_y {
                let dy = grid_y[y] as f64 - vy;
                let dxy_sq = dx * dx + dy * dy;
                let xy_offset = y * dim_z + x_offset;
                for z in begin_z..end_z {
                    let dz = grid_z[z] as f64 - vz;
                    let distance_sq = dxy_sq + dz * dz;
                    if distance_sq <= radius2_sq {
                        let density = faster_exp(-alpha * (distance_sq * radius_sq_inv));
                        let offset = z + xy_offset;

                        // Typed-array compound assignment reads Float32, performs
                        // the addition as a JS number, then writes Float32 again.
                        field[offset] = (field[offset] as f64 + density as f64) as f32;
                        if density > dominant_density[offset] {
                            dominant_density[offset] = density;
                            id_field[offset] = point.group_id as f32;
                        }
                    }
                }
            }
        }
    }

    GaussianDensityGrid {
        dimensions,
        field,
        id_field,
        origin,
        resolution,
    }
}

pub(super) fn marching_cubes_mesh(grid: &GaussianDensityGrid, iso_level: f64) -> Mesh {
    assert!(
        iso_level.is_finite(),
        "Marching cubes iso level must be finite"
    );
    let [nx, ny, nz] = grid.dimensions;
    let expected_len = nx
        .checked_mul(ny)
        .and_then(|value| value.checked_mul(nz))
        .expect("Marching cubes grid is too large");
    assert_eq!(
        grid.field.len(),
        expected_len,
        "Marching cubes scalar field length does not match dimensions"
    );
    assert_eq!(
        grid.id_field.len(),
        expected_len,
        "Marching cubes id field length does not match dimensions"
    );
    if nx < 2 || ny < 2 || nz < 2 {
        return Mesh::default();
    }

    let mut state = MarchingCubesState::new(grid, iso_level);
    for z in 0..(nz - 1) {
        for y in 0..(ny - 1) {
            for x in 0..(nx - 1) {
                state.process_cell(x, y, z);
            }
        }
        state.clear_edge_vertex_index_slice(z);
    }

    state.finish()
}

#[cfg(test)]
pub(super) fn build_gaussian_surface_mesh(
    points: &[GaussianDensityPoint],
    params: GaussianDensityParams,
) -> Mesh {
    let grid = build_gaussian_density_grid(points, params);
    marching_cubes_mesh(&grid, (-(params.smoothness as f64)).exp())
}

pub(super) fn build_gaussian_surface_mesh_in_box(
    points: &[GaussianDensityPoint],
    params: GaussianDensityParams,
    box_min: Vec3,
    box_max: Vec3,
) -> Mesh {
    let grid = build_gaussian_density_grid_in_box(points, params, box_min, box_max);
    marching_cubes_mesh(&grid, (-(params.smoothness as f64)).exp())
}

fn fill_grid_dimension(length: usize, start: f64, step: f64) -> Vec<f32> {
    (0..length)
        .map(|index| (start + step * index as f64) as f32)
        .collect()
}

fn faster_exp(value: f64) -> f32 {
    faster_pow2(1.442_695_040 * value)
}

fn faster_pow2(value: f64) -> f32 {
    let clipped = if value < -126.0 { -126.0 } else { value };
    let int_bits = ((1_u32 << 23) as f64 * (clipped + 126.942_695_04)).trunc() as i32;
    f32::from_bits(int_bits as u32)
}

struct MarchingCubesState<'a> {
    grid: &'a GaussianDensityGrid,
    iso_level: f64,
    vertices_on_edges: Vec<i32>,
    vertex_list: [i32; 12],
    cell: [usize; 3],
    vertices: Vec<Vec3>,
    normals: Vec<Vec3>,
    groups: Vec<f32>,
    faces: Vec<Face>,
}

impl<'a> MarchingCubesState<'a> {
    fn new(grid: &'a GaussianDensityGrid, iso_level: f64) -> Self {
        let [nx, ny, _] = grid.dimensions;
        Self {
            grid,
            iso_level,
            vertices_on_edges: vec![0; 3 * nx * ny * 2],
            vertex_list: [0; 12],
            cell: [0; 3],
            vertices: Vec::new(),
            normals: Vec::new(),
            groups: Vec::new(),
            faces: Vec::new(),
        }
    }

    fn clear_edge_vertex_index_slice(&mut self, z: usize) {
        let [nx, ny, _] = self.grid.dimensions;
        let half = 3 * nx * ny;
        let range = if z % 2 == 0 {
            0..half
        } else {
            half..self.vertices_on_edges.len()
        };
        self.vertices_on_edges[range].fill(0);
    }

    fn edge_offset(&self, edge_number: usize) -> usize {
        let [nx, ny, _] = self.grid.dimensions;
        let info = EDGE_ID_INFO[edge_number];
        let x = self.cell[0] + info.i;
        let y = self.cell[1] + info.j;
        let z = (self.cell[2] + info.k) % 2;
        3 * (nx * (z * ny + y) + x) + info.edge
    }

    fn interpolate(&mut self, edge_number: usize) -> i32 {
        let edge_offset = self.edge_offset(edge_number);
        let cached = self.vertices_on_edges[edge_offset];
        if cached > 0 {
            return cached - 1;
        }

        let edge = CUBE_EDGES[edge_number];
        let low = [
            self.cell[0] + edge.a.i,
            self.cell[1] + edge.a.j,
            self.cell[2] + edge.a.k,
        ];
        let high = [
            self.cell[0] + edge.b.i,
            self.cell[1] + edge.b.j,
            self.cell[2] + edge.b.k,
        ];
        let value0 = self.grid.scalar(low[0], low[1], low[2]) as f64;
        let value1 = self.grid.scalar(high[0], high[1], high[2]) as f64;
        let t = (self.iso_level - value0) / (value0 - value1);

        let u = self.grid.group(low[0], low[1], low[2]);
        let v = self.grid.group(high[0], high[1], high[2]);
        let mut group = if t < 0.5 { u } else { v };
        if group == -1.0 {
            group = if t < 0.5 { v } else { u };
        }
        if group == -2.0 {
            return -1;
        }
        self.groups.push(group);

        let vertex_id = self.vertices.len() as i32;
        self.vertices.push(Vec3::new(
            (low[0] as f64 + t * (low[0] as f64 - high[0] as f64)) as f32,
            (low[1] as f64 + t * (low[1] as f64 - high[1] as f64)) as f32,
            (low[2] as f64 + t * (low[2] as f64 - high[2] as f64)) as f32,
        ));
        self.vertices_on_edges[edge_offset] = vertex_id + 1;

        let [nx, ny, nz] = self.grid.dimensions;
        let normal0 = [
            self.grid.scalar(low[0].saturating_sub(1), low[1], low[2]) as f64
                - self.grid.scalar((low[0] + 1).min(nx - 1), low[1], low[2]) as f64,
            self.grid.scalar(low[0], low[1].saturating_sub(1), low[2]) as f64
                - self.grid.scalar(low[0], (low[1] + 1).min(ny - 1), low[2]) as f64,
            self.grid.scalar(low[0], low[1], low[2].saturating_sub(1)) as f64
                - self.grid.scalar(low[0], low[1], (low[2] + 1).min(nz - 1)) as f64,
        ];
        let normal1 = [
            self.grid
                .scalar(high[0].saturating_sub(1), high[1], high[2]) as f64
                - self
                    .grid
                    .scalar((high[0] + 1).min(nx - 1), high[1], high[2]) as f64,
            self.grid
                .scalar(high[0], high[1].saturating_sub(1), high[2]) as f64
                - self
                    .grid
                    .scalar(high[0], (high[1] + 1).min(ny - 1), high[2]) as f64,
            self.grid
                .scalar(high[0], high[1], high[2].saturating_sub(1)) as f64
                - self
                    .grid
                    .scalar(high[0], high[1], (high[2] + 1).min(nz - 1)) as f64,
        ];
        let normal = [
            normal0[0] + t * (normal0[0] - normal1[0]),
            normal0[1] + t * (normal0[1] - normal1[1]),
            normal0[2] + t * (normal0[2] - normal1[2]),
        ];
        let sign = if self.iso_level >= 0.0 { 1.0 } else { -1.0 };
        self.normals.push(Vec3::new(
            (sign * normal[0]) as f32,
            (sign * normal[1]) as f32,
            (sign * normal[2]) as f32,
        ));

        vertex_id
    }

    fn process_cell(&mut self, x: usize, y: usize, z: usize) {
        let mut table_index = 0_usize;
        if (self.grid.scalar(x, y, z) as f64) < self.iso_level {
            table_index |= 1;
        }
        if (self.grid.scalar(x + 1, y, z) as f64) < self.iso_level {
            table_index |= 2;
        }
        if (self.grid.scalar(x + 1, y + 1, z) as f64) < self.iso_level {
            table_index |= 4;
        }
        if (self.grid.scalar(x, y + 1, z) as f64) < self.iso_level {
            table_index |= 8;
        }
        if (self.grid.scalar(x, y, z + 1) as f64) < self.iso_level {
            table_index |= 16;
        }
        if (self.grid.scalar(x + 1, y, z + 1) as f64) < self.iso_level {
            table_index |= 32;
        }
        if (self.grid.scalar(x + 1, y + 1, z + 1) as f64) < self.iso_level {
            table_index |= 64;
        }
        if (self.grid.scalar(x, y + 1, z + 1) as f64) < self.iso_level {
            table_index |= 128;
        }
        if table_index == 0 || table_index == 255 {
            return;
        }

        self.cell = [x, y, z];
        let edge_info = EDGE_TABLE[table_index];
        for edge_number in 0..12 {
            if edge_info & (1 << edge_number) != 0 {
                self.vertex_list[edge_number] = self.interpolate(edge_number);
            }
        }

        for triangle in TRI_TABLE[table_index].chunks_exact(3) {
            let mut a = self.vertex_list[triangle[0] as usize];
            let b = self.vertex_list[triangle[1] as usize];
            let mut c = self.vertex_list[triangle[2] as usize];
            if self.iso_level < 0.0 {
                std::mem::swap(&mut a, &mut c);
            }
            if a >= 0 && b >= 0 && c >= 0 {
                self.faces.push(Face {
                    a: a as usize,
                    b: b as usize,
                    c: c as usize,
                });
            }
        }
    }

    fn finish(mut self) -> Mesh {
        for vertex in &mut self.vertices {
            vertex.x = (vertex.x as f64 * self.grid.resolution + self.grid.origin[0]) as f32;
            vertex.y = (vertex.y as f64 * self.grid.resolution + self.grid.origin[1]) as f32;
            vertex.z = (vertex.z as f64 * self.grid.resolution + self.grid.origin[2]) as f32;
        }

        let vertex_groups = self
            .groups
            .iter()
            .map(|&group| {
                if group.is_finite() && group >= 0.0 {
                    group as usize
                } else {
                    0
                }
            })
            .collect::<Vec<_>>();
        let face_groups = self
            .faces
            .iter()
            .map(|face| vertex_groups[face.a])
            .collect::<Vec<_>>();
        let group_count = vertex_groups
            .iter()
            .copied()
            .max()
            .map_or(0, |group| group + 1);

        Mesh {
            vertices: self.vertices,
            normals: self.normals,
            faces: self.faces,
            vertex_groups,
            face_groups,
            face_materials: Vec::new(),
            sections: Vec::new(),
            group_count,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> GaussianDensityParams {
        GaussianDensityParams {
            resolution: 1.0,
            radius_offset: 0.0,
            smoothness: 1.5,
        }
    }

    fn point(x: f32, y: f32, z: f32, radius: f64, group_id: usize) -> GaussianDensityPoint {
        GaussianDensityPoint::new(Vec3::new(x, y, z), radius, group_id)
    }

    #[test]
    fn single_point_density_uses_molstar_float32_layout() {
        let grid = build_gaussian_density_grid(&[point(0.0, 0.0, 0.0, 1.0, 7)], params());

        assert_eq!(grid.dimensions, [6, 6, 6]);
        assert_eq!(grid.origin, [-3.0, -3.0, -3.0]);
        let center = grid.offset(3, 3, 3);
        assert_eq!(grid.field[center].to_bits(), faster_exp(0.0).to_bits());
        assert_eq!(grid.id_field[center], 7.0);
        assert_eq!(grid.offset(1, 2, 3), 3 + 6 * (2 + 6));
    }

    #[test]
    fn density_expands_the_supplied_component_boundary_box() {
        let grid = build_gaussian_density_grid_in_box(
            &[point(0.0, 0.0, 0.0, 1.0, 7)],
            params(),
            Vec3::new(-2.0, -1.0, -1.0),
            Vec3::new(4.0, 1.0, 1.0),
        );

        assert_eq!(grid.origin, [-5.0, -4.0, -4.0]);
        assert_eq!(grid.dimensions, [12, 8, 8]);
    }

    #[test]
    fn multiple_points_accumulate_in_input_order_and_keep_first_equal_group() {
        let points = [point(-1.0, 0.0, 0.0, 1.0, 3), point(1.0, 0.0, 0.0, 1.0, 9)];
        let grid = build_gaussian_density_grid(&points, params());
        let midpoint = grid.offset(4, 3, 3);
        let contribution = faster_exp(-(params().smoothness as f64));
        let expected = (contribution as f64 + contribution as f64) as f32;

        assert_eq!(grid.field[midpoint].to_bits(), expected.to_bits());
        assert_eq!(grid.id_field[midpoint], 3.0);
    }

    #[test]
    fn marching_cubes_preserves_vertex_group_and_face_order() {
        let mut grid = GaussianDensityGrid {
            dimensions: [2, 2, 2],
            field: vec![1.0; 8],
            id_field: vec![9.0; 8],
            origin: [10.0, 20.0, 30.0],
            resolution: 2.0,
        };
        let corner = grid.offset(0, 0, 0);
        grid.field[corner] = 0.0;
        grid.id_field[corner] = 2.0;

        let mesh = marching_cubes_mesh(&grid, 0.5);

        assert_eq!(mesh.vertices.len(), 3);
        assert_eq!(mesh.normals.len(), 3);
        assert_eq!(mesh.vertex_groups, [2, 9, 2]);
        assert_eq!(mesh.faces.len(), 1);
        assert_eq!(
            (mesh.faces[0].a, mesh.faces[0].b, mesh.faces[0].c),
            (0, 2, 1)
        );
        assert_eq!(mesh.vertices[0], Vec3::new(11.0, 20.0, 30.0));
        assert_eq!(mesh.vertices[1], Vec3::new(10.0, 21.0, 30.0));
        assert_eq!(mesh.vertices[2], Vec3::new(10.0, 20.0, 31.0));
    }

    #[test]
    fn gaussian_surface_is_deterministic_and_finite() {
        let points = [
            point(-0.75, 0.0, 0.0, 1.4, 4),
            point(0.75, 0.0, 0.0, 1.4, 8),
            point(0.0, 1.0, 0.25, 1.1, 12),
        ];
        let params = GaussianDensityParams {
            resolution: 0.5,
            radius_offset: 0.2,
            smoothness: 1.5,
        };
        let a = build_gaussian_surface_mesh(&points, params);
        let b = build_gaussian_surface_mesh(&points, params);

        assert!(!a.vertices.is_empty());
        assert!(!a.faces.is_empty());
        assert_eq!(a.vertices.len(), a.normals.len());
        assert_eq!(a.vertices.len(), a.vertex_groups.len());
        assert_eq!(a.faces.len(), a.face_groups.len());
        assert!(a.vertices.iter().all(|value| value.is_finite()));
        assert!(a.normals.iter().all(|value| value.is_finite()));
        assert!(a.faces.iter().all(|face| face.a < a.vertices.len()
            && face.b < a.vertices.len()
            && face.c < a.vertices.len()));
        assert!(a.vertex_groups.contains(&4));
        assert!(a.vertex_groups.contains(&8));
        assert!(a.vertex_groups.contains(&12));

        assert_eq!(
            a.vertices
                .iter()
                .flat_map(|value| [value.x.to_bits(), value.y.to_bits(), value.z.to_bits()])
                .collect::<Vec<_>>(),
            b.vertices
                .iter()
                .flat_map(|value| [value.x.to_bits(), value.y.to_bits(), value.z.to_bits()])
                .collect::<Vec<_>>()
        );
        assert_eq!(
            a.normals
                .iter()
                .flat_map(|value| [value.x.to_bits(), value.y.to_bits(), value.z.to_bits()])
                .collect::<Vec<_>>(),
            b.normals
                .iter()
                .flat_map(|value| [value.x.to_bits(), value.y.to_bits(), value.z.to_bits()])
                .collect::<Vec<_>>()
        );
        assert_eq!(a.vertex_groups, b.vertex_groups);
        assert_eq!(
            a.faces
                .iter()
                .map(|face| (face.a, face.b, face.c))
                .collect::<Vec<_>>(),
            b.faces
                .iter()
                .map(|face| (face.a, face.b, face.c))
                .collect::<Vec<_>>()
        );
    }
}
