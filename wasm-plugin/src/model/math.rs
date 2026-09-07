#[derive(Clone, Copy, Debug)]
pub struct Transform {
    pub m: [[f32; 4]; 3],
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Axes3D {
    pub origin: Vec3,
    pub dir_a: Vec3,
    pub dir_b: Vec3,
    pub dir_c: Vec3,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub(crate) struct Axes3D64 {
    pub origin: [f64; 3],
    pub dir_a: [f64; 3],
    pub dir_b: [f64; 3],
    pub dir_c: [f64; 3],
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PrincipalAxes {
    pub moments_axes: Axes3D,
    pub box_axes: Axes3D,
}

impl Vec3 {
    pub const fn new(x: f32, y: f32, z: f32) -> Vec3 {
        Vec3 { x, y, z }
    }

    pub(crate) fn length(self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub(crate) fn squared_distance(self, other: Vec3) -> f32 {
        let d = self - other;
        d.x * d.x + d.y * d.y + d.z * d.z
    }

    pub(crate) fn distance(self, other: Vec3) -> f32 {
        self.squared_distance(other).sqrt()
    }

    pub(crate) fn normalized(self) -> Vec3 {
        let len = self.length();
        if len <= 0.000_001 {
            Vec3::default()
        } else {
            self / len
        }
    }

    pub(crate) fn cross(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }

    pub(crate) fn dot(self, other: Vec3) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub(crate) fn squared_length(self) -> f32 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub(crate) fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub(crate) fn min(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x.min(other.x),
            y: self.y.min(other.y),
            z: self.z.min(other.z),
        }
    }

    pub(crate) fn max(self, other: Vec3) -> Vec3 {
        Vec3 {
            x: self.x.max(other.x),
            y: self.y.max(other.y),
            z: self.z.max(other.z),
        }
    }
}

impl Axes3D {
    fn create(origin: Vec3, dir_a: Vec3, dir_b: Vec3, dir_c: Vec3) -> Self {
        Self {
            origin,
            dir_a,
            dir_b,
            dir_c,
        }
    }

    fn normalize(&self) -> Self {
        Self {
            origin: self.origin,
            dir_a: molstar_vec3_normalize(self.dir_a),
            dir_b: molstar_vec3_normalize(self.dir_b),
            dir_c: molstar_vec3_normalize(self.dir_c),
        }
    }
}

impl Axes3D64 {
    fn create(origin: [f64; 3], dir_a: [f64; 3], dir_b: [f64; 3], dir_c: [f64; 3]) -> Self {
        Self {
            origin,
            dir_a,
            dir_b,
            dir_c,
        }
    }

    fn normalize(&self) -> Self {
        Self {
            origin: self.origin,
            dir_a: molstar_vec3_normalize64(self.dir_a),
            dir_b: molstar_vec3_normalize64(self.dir_b),
            dir_c: molstar_vec3_normalize64(self.dir_c),
        }
    }

    fn to_f32(&self) -> Axes3D {
        Axes3D::create(
            Vec3::new(
                self.origin[0] as f32,
                self.origin[1] as f32,
                self.origin[2] as f32,
            ),
            Vec3::new(
                self.dir_a[0] as f32,
                self.dir_a[1] as f32,
                self.dir_a[2] as f32,
            ),
            Vec3::new(
                self.dir_b[0] as f32,
                self.dir_b[1] as f32,
                self.dir_b[2] as f32,
            ),
            Vec3::new(
                self.dir_c[0] as f32,
                self.dir_c[1] as f32,
                self.dir_c[2] as f32,
            ),
        )
    }
}

impl PrincipalAxes {
    pub(crate) fn of_positions(positions: &[Vec3]) -> Self {
        let moments_axes = Self::calculate_moments_axes(positions);
        let box_axes = Self::calculate_box_axes(positions, &moments_axes);
        Self {
            moments_axes,
            box_axes,
        }
    }

    pub(crate) fn calculate_moments_axes(positions: &[Vec3]) -> Axes3D {
        Self::calculate_moments_axes64(positions).to_f32()
    }

    pub(crate) fn calculate_moments_axes64(positions: &[Vec3]) -> Axes3D64 {
        if positions.is_empty() {
            return Axes3D64::default();
        }
        if positions.len() == 1 {
            return Axes3D64::create(
                [
                    positions[0].x as f64,
                    positions[0].y as f64,
                    positions[0].z as f64,
                ],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            );
        }

        let n = positions.len();
        let n3 = n as f64 / 3.0;
        let mut mean = [0.0f64; 3];
        for position in positions {
            mean[0] += position.x as f64;
            mean[1] += position.y as f64;
            mean[2] += position.z as f64;
        }
        mean[0] /= n as f64;
        mean[1] /= n as f64;
        mean[2] /= n as f64;

        let centered = positions
            .iter()
            .map(|position| {
                [
                    (position.x as f64 - mean[0]) as f32,
                    (position.y as f64 - mean[1]) as f32,
                    (position.z as f64 - mean[2]) as f32,
                ]
            })
            .collect::<Vec<_>>();

        let mut a = [0.0f32; 9];
        for row in 0..3 {
            for col in 0..3 {
                let mut sum = 0.0f64;
                for point in &centered {
                    sum += point[row] as f64 * point[col] as f64;
                }
                a[row * 3 + col] = sum as f32;
            }
        }

        let (w, u) = molstar_svd_3x3(a);
        let scale_a = (w[0] as f64 / n3).sqrt();
        let scale_b = (w[1] as f64 / n3).sqrt();
        let scale_c = (w[2] as f64 / n3).sqrt();
        Axes3D64::create(
            mean,
            [
                u[0] as f64 * scale_a,
                u[3] as f64 * scale_a,
                u[6] as f64 * scale_a,
            ],
            [
                u[1] as f64 * scale_b,
                u[4] as f64 * scale_b,
                u[7] as f64 * scale_b,
            ],
            [
                u[2] as f64 * scale_c,
                u[5] as f64 * scale_c,
                u[8] as f64 * scale_c,
            ],
        )
    }

    pub(crate) fn calculate_normalized_axes(moments_axes: &Axes3D) -> Axes3D {
        let mut axes = moments_axes.clone();
        if axes.dir_c.length() < 0.000_001 {
            axes.dir_c = axes.dir_a.cross(axes.dir_b);
        }
        axes.normalize()
    }

    pub(crate) fn calculate_normalized_axes64(moments_axes: &Axes3D64) -> Axes3D64 {
        let mut axes = moments_axes.clone();
        if molstar_vec3_magnitude64(axes.dir_c) < 0.000_001 {
            axes.dir_c = molstar_vec3_cross64(axes.dir_a, axes.dir_b);
        }
        axes.normalize()
    }

    pub(crate) fn calculate_box_axes(positions: &[Vec3], moments_axes: &Axes3D) -> Axes3D {
        if positions.is_empty() {
            return Axes3D::default();
        }
        if positions.len() == 1 {
            return moments_axes.clone();
        }

        let mut d1a = f32::NEG_INFINITY;
        let mut d1b = f32::NEG_INFINITY;
        let mut d2a = f32::NEG_INFINITY;
        let mut d2b = f32::NEG_INFINITY;
        let mut d3a = f32::NEG_INFINITY;
        let mut d3b = f32::NEG_INFINITY;

        let center = moments_axes.origin;
        let axes = Self::calculate_normalized_axes(moments_axes);

        for &position in positions {
            let projected = molstar_project_point_on_vector(position, axes.dir_a, center);
            let dp = axes.dir_a.dot(molstar_vec3_normalize(projected - center));
            let dt = projected.distance(center);
            if dp > 0.0 {
                if dt > d1a {
                    d1a = dt;
                }
            } else if dt > d1b {
                d1b = dt;
            }

            let projected = molstar_project_point_on_vector(position, axes.dir_b, center);
            let dp = axes.dir_b.dot(molstar_vec3_normalize(projected - center));
            let dt = projected.distance(center);
            if dp > 0.0 {
                if dt > d2a {
                    d2a = dt;
                }
            } else if dt > d2b {
                d2b = dt;
            }

            let projected = molstar_project_point_on_vector(position, axes.dir_c, center);
            let dp = axes.dir_c.dot(molstar_vec3_normalize(projected - center));
            let dt = projected.distance(center);
            if dp > 0.0 {
                if dt > d3a {
                    d3a = dt;
                }
            } else if dt > d3b {
                d3b = dt;
            }
        }

        let dir_a = molstar_vec3_set_magnitude(axes.dir_a, (d1a + d1b) / 2.0);
        let dir_b = molstar_vec3_set_magnitude(axes.dir_b, (d2a + d2b) / 2.0);
        let dir_c = molstar_vec3_set_magnitude(axes.dir_c, (d3a + d3b) / 2.0);

        let ok_a = dir_a.is_finite();
        let ok_b = dir_b.is_finite();
        let ok_c = dir_c.is_finite();
        let mut origin = Vec3::default();
        for (a, b, c) in [
            (d1a, d2a, d3a),
            (d1a, d2a, -d3b),
            (d1a, -d2b, -d3b),
            (d1a, -d2b, d3a),
            (-d1b, -d2b, -d3b),
            (-d1b, -d2b, d3a),
            (-d1b, d2a, d3a),
            (-d1b, d2a, -d3b),
        ] {
            let mut corner = center;
            if ok_a {
                corner = corner + axes.dir_a * a;
            }
            if ok_b {
                corner = corner + axes.dir_b * b;
            }
            if ok_c {
                corner = corner + axes.dir_c * c;
            }
            origin = origin + corner;
        }
        origin = origin * (1.0 / 8.0);

        Axes3D::create(origin, dir_a, dir_b, dir_c)
    }
}

fn molstar_vec3_normalize(v: Vec3) -> Vec3 {
    let mut out = v;
    let mut len = v.squared_length();
    if len > 0.0 {
        len = 1.0 / len.sqrt();
        out = v * len;
    }
    out
}

fn molstar_vec3_magnitude64(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn molstar_vec3_normalize64(v: [f64; 3]) -> [f64; 3] {
    let mut out = v;
    let mut len = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if len > 0.0 {
        len = 1.0 / len.sqrt();
        out = [v[0] * len, v[1] * len, v[2] * len];
    }
    out
}

fn molstar_vec3_cross64(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn molstar_vec3_set_magnitude(v: Vec3, length: f32) -> Vec3 {
    molstar_vec3_normalize(v) * length
}

fn molstar_project_point_on_vector(point: Vec3, vector: Vec3, origin: Vec3) -> Vec3 {
    let out = point - origin;
    let scalar = vector.dot(out) / vector.squared_length();
    origin + vector * scalar
}

fn molstar_svd_3x3(a: [f32; 9]) -> ([f32; 3], [f32; 9]) {
    let mut amt = [0.0f32; 9];
    for i in 0..3 {
        for j in 0..3 {
            amt[i * 3 + j] = a[j * 3 + i];
        }
    }
    let mut w = [0.0f32; 3];
    let mut vt = [0.0f32; 9];
    molstar_jacobi_svd_impl(&mut amt, 3, &mut w, &mut vt, 3, 3, 3, 3);

    let mut u = [0.0f32; 9];
    for i in 0..3 {
        for j in 0..3 {
            u[i * 3 + j] = amt[j * 3 + i];
        }
    }
    (w, u)
}

#[allow(clippy::too_many_arguments)]
fn molstar_jacobi_svd_impl(
    at: &mut [f32; 9],
    astep: usize,
    out_w: &mut [f32; 3],
    vt: &mut [f32; 9],
    vstep: usize,
    m: usize,
    n: usize,
    n1: usize,
) {
    const EPSILON: f64 = 0.0000001192092896;
    const FLT_MIN: f64 = 1E-37;
    let eps = EPSILON * 2.0;
    let minval = FLT_MIN;
    let max_iter = m.max(30);
    let mut w = [0.0f64; 24];

    for i in 0..n {
        let mut sd = 0.0;
        for k in 0..m {
            let t = at[i * astep + k] as f64;
            sd += t * t;
        }
        w[i] = sd;
        for k in 0..n {
            vt[i * vstep + k] = 0.0;
        }
        vt[i * vstep + i] = 1.0;
    }

    for _ in 0..max_iter {
        let mut changed = 0;
        for i in 0..n - 1 {
            for j in i + 1..n {
                let ai = i * astep;
                let aj = j * astep;
                let mut a = w[i];
                let mut p = 0.0;
                let mut b = w[j];

                p += at[ai] as f64 * at[aj] as f64;
                p += at[ai + 1] as f64 * at[aj + 1] as f64;
                for k in 2..m {
                    p += at[ai + k] as f64 * at[aj + k] as f64;
                }
                if p.abs() <= eps * (a * b).sqrt() {
                    continue;
                }

                p *= 2.0;
                let beta = a - b;
                let gamma = molstar_hypot(p, beta);
                let (c, s) = if beta < 0.0 {
                    let delta = (gamma - beta) * 0.5;
                    let s = (delta / gamma).sqrt();
                    let c = p / (gamma * s * 2.0);
                    (c, s)
                } else {
                    let c = ((gamma + beta) / (gamma * 2.0)).sqrt();
                    let s = p / (gamma * c * 2.0);
                    (c, s)
                };

                a = 0.0;
                b = 0.0;
                for k in 0..m {
                    let t0 = c * at[ai + k] as f64 + s * at[aj + k] as f64;
                    let t1 = -s * at[ai + k] as f64 + c * at[aj + k] as f64;
                    at[ai + k] = t0 as f32;
                    at[aj + k] = t1 as f32;
                    a += t0 * t0;
                    b += t1 * t1;
                }
                w[i] = a;
                w[j] = b;
                changed = 1;

                let vi = i * vstep;
                let vj = j * vstep;
                for k in 0..n {
                    let t0 = c * vt[vi + k] as f64 + s * vt[vj + k] as f64;
                    let t1 = -s * vt[vi + k] as f64 + c * vt[vj + k] as f64;
                    vt[vi + k] = t0 as f32;
                    vt[vj + k] = t1 as f32;
                }
            }
        }
        if changed == 0 {
            break;
        }
    }

    for i in 0..n {
        let mut sd = 0.0;
        for k in 0..m {
            let t = at[i * astep + k] as f64;
            sd += t * t;
        }
        w[i] = sd.sqrt();
    }

    for i in 0..n - 1 {
        let mut j = i;
        for k in i + 1..n {
            if w[j] < w[k] {
                j = k;
            }
        }
        if i != j {
            w.swap(i, j);
            for k in 0..m {
                at.swap(i * astep + k, j * astep + k);
            }
            for k in 0..n {
                vt.swap(i * vstep + k, j * vstep + k);
            }
        }
    }

    for i in 0..n {
        out_w[i] = w[i] as f32;
    }

    let mut seed = 0x1234 as f64;
    for i in 0..n1 {
        let mut sd = if i < n { w[i] } else { 0.0 };
        while sd <= minval {
            let val0 = 1.0 / m as f64;
            for k in 0..m {
                seed = seed * 214013.0 + 2531011.0;
                let shifted = js_signed_right_shift(seed, 16);
                let val = if (shifted & 0x7fff) & 256 != 0 {
                    val0
                } else {
                    -val0
                };
                at[i * astep + k] = val as f32;
            }
            for _ in 0..2 {
                for j in 0..i {
                    sd = 0.0;
                    for k in 0..m {
                        sd += at[i * astep + k] as f64 * at[j * astep + k] as f64;
                    }
                    let mut asum = 0.0;
                    for k in 0..m {
                        let t = at[i * astep + k] as f64 - sd * at[j * astep + k] as f64;
                        at[i * astep + k] = t as f32;
                        asum += t.abs();
                    }
                    let scale = if asum != 0.0 { 1.0 / asum } else { 0.0 };
                    for k in 0..m {
                        at[i * astep + k] = (at[i * astep + k] as f64 * scale) as f32;
                    }
                }
            }
            sd = 0.0;
            for k in 0..m {
                let t = at[i * astep + k] as f64;
                sd += t * t;
            }
            sd = sd.sqrt();
        }

        let s = 1.0 / sd;
        for k in 0..m {
            at[i * astep + k] = (at[i * astep + k] as f64 * s) as f32;
        }
    }
}

fn molstar_hypot(mut a: f64, mut b: f64) -> f64 {
    a = a.abs();
    b = b.abs();
    if a > b {
        b /= a;
        return a * (1.0 + b * b).sqrt();
    }
    if b > 0.0 {
        a /= b;
        return b * (1.0 + a * a).sqrt();
    }
    0.0
}

fn js_signed_right_shift(value: f64, bits: u32) -> i32 {
    let two32 = 4_294_967_296.0;
    let mut int = value.trunc() % two32;
    if int < 0.0 {
        int += two32;
    }
    let unsigned = int as u32;
    (unsigned as i32) >> bits
}

impl Transform {
    pub(crate) fn identity() -> Self {
        Self {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
            ],
        }
    }

    pub(crate) fn apply(self, v: Vec3) -> Vec3 {
        Vec3 {
            x: self.m[0][0] * v.x + self.m[0][1] * v.y + self.m[0][2] * v.z + self.m[0][3],
            y: self.m[1][0] * v.x + self.m[1][1] * v.y + self.m[1][2] * v.z + self.m[1][3],
            z: self.m[2][0] * v.x + self.m[2][1] * v.y + self.m[2][2] * v.z + self.m[2][3],
        }
    }

    pub(crate) fn is_identity(self) -> bool {
        let identity = Transform::identity();
        self.m.iter().enumerate().all(|(row, values)| {
            values
                .iter()
                .enumerate()
                .all(|(col, value)| (*value - identity.m[row][col]).abs() < 0.000_001)
        })
    }

    pub(crate) fn then(self, next: Transform) -> Transform {
        let mut m = [[0.0; 4]; 3];
        for (row, m_row) in m.iter_mut().enumerate() {
            for (col, value) in m_row.iter_mut().take(3).enumerate() {
                *value = next.m[row][0] * self.m[0][col]
                    + next.m[row][1] * self.m[1][col]
                    + next.m[row][2] * self.m[2][col];
            }
            m_row[3] = next.m[row][0] * self.m[0][3]
                + next.m[row][1] * self.m[1][3]
                + next.m[row][2] * self.m[2][3]
                + next.m[row][3];
        }
        Transform { m }
    }
}

impl std::ops::Add for Vec3 {
    type Output = Vec3;
    fn add(self, rhs: Vec3) -> Vec3 {
        Vec3 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, rhs: Vec3) -> Vec3 {
        Vec3 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl std::ops::Mul<f32> for Vec3 {
    type Output = Vec3;
    fn mul(self, rhs: f32) -> Vec3 {
        Vec3 {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl std::ops::Div<f32> for Vec3 {
    type Output = Vec3;
    fn div(self, rhs: f32) -> Vec3 {
        Vec3 {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}
