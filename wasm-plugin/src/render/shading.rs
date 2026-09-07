use crate::model::{MeshMaterial, Vec3};

#[cfg(test)]
use super::color::quantize;
use super::color::{color_f32, parse_color};
use super::options::{DirectionalLight, RendererOptions};
use super::style::ResolvedStyle;

const MOLSTAR_LIGHT: Vec3 = Vec3 {
    x: 0.383_022_22,
    y: -0.321_393_82,
    z: -0.866_025_4,
};

#[cfg(test)]
pub(super) fn shade_material(
    material: MeshMaterial,
    geometric_normal: Vec3,
    view_position: Vec3,
    renderer: &RendererOptions,
    style: ResolvedStyle,
) -> [u8; 4] {
    shade_material_linear(material, geometric_normal, view_position, renderer, style).map(quantize)
}

/// Evaluate the Mol* material shader before the destination attachment
/// conversion. Opaque rendering stores this value in RGBA8, while WBOIT
/// consumes it in float32 accumulation attachments and must not inherit an
/// early eight-bit quantization.
#[cfg(test)]
pub(super) fn shade_material_linear(
    material: MeshMaterial,
    geometric_normal: Vec3,
    view_position: Vec3,
    renderer: &RendererOptions,
    style: ResolvedStyle,
) -> [f32; 4] {
    shade_material_base_color_linear(
        material,
        color_f32(material.color),
        geometric_normal,
        view_position,
        renderer,
        style,
    )
}

/// Evaluate a mesh fragment whose RGB value was produced by Mol*'s vertex
/// color varying. Alpha remains the render-object uniform represented by
/// `material`, rather than becoming a varying.
pub(super) fn shade_material_base_color_linear(
    material: MeshMaterial,
    base: [f32; 3],
    geometric_normal: Vec3,
    view_position: Vec3,
    renderer: &RendererOptions,
    style: ResolvedStyle,
) -> [f32; 4] {
    let alpha = material.alpha_tenths.min(10) as f32 / 10.0;
    let mut out = base;
    if !style.ignore_light {
        let metalness = renderer.shading.material.metalness.clamp(0.0, 0.99);
        let roughness = renderer.shading.material.roughness.clamp(0.0525, 1.0);
        let diffuse_color = scale3(base, 1.0 - metalness);
        let specular_color = [
            0.04 * (1.0 - metalness) + base[0] * metalness,
            0.04 * (1.0 - metalness) + base[1] * metalness,
            0.04 * (1.0 - metalness) + base[2] * metalness,
        ];
        let ambient = molstar_scaled_uniform_color(
            parse_color(&renderer.lighting.ambient.color).unwrap_or(0xffffff),
            renderer.lighting.ambient.intensity,
        );
        out = mul3(diffuse_color, ambient);
        // Mol*'s mesh, sphere, and cylinder fragment shaders negate their
        // geometric surface normal before evaluating apply_light_color.glsl.
        // `vViewPosition` remains the camera-space fragment position, so its
        // normalized value points in the same (camera-to-fragment) hemisphere
        // as that shader normal. Keeping the conventional outward normal and
        // fragment-to-camera view vector here reverses both vectors and makes
        // the default Mol* light illuminate the back hemisphere instead.
        let normal = (geometric_normal * -1.0).normalized();
        let view_dir = view_position.normalized();
        for light in &renderer.lighting.directional {
            let direction = light_direction(light);
            let light_color = molstar_scaled_uniform_color(
                parse_color(&light.color).unwrap_or(0xffffff),
                light.intensity,
            );
            let n_dot_l = normal.dot(direction).max(0.0);
            let diffuse = scale3(mul3(diffuse_color, light_color), n_dot_l);
            out = add3(out, diffuse);
            let brdf = brdf_ggx(direction, view_dir, normal, specular_color, roughness);
            // Mol* multiplies punctual light color by PI before evaluating the
            // BRDF; Lambert's reciprocal PI above cancels for the diffuse term.
            out = add3(
                out,
                scale3(mul3(light_color, brdf), n_dot_l * std::f32::consts::PI),
            );
        }
        if metalness > 0.0 {
            let (single_scatter, multi_scatter) =
                ggx_multiscattering(normal, view_dir, specular_color, roughness);
            let radiance = scale3(ambient, metalness);
            let cosine_irradiance = scale3(radiance, 1.0 / std::f32::consts::PI);
            let energy = [
                1.0 - single_scatter[0] - multi_scatter[0],
                1.0 - single_scatter[1] - multi_scatter[1],
                1.0 - single_scatter[2] - multi_scatter[2],
            ];
            out = add3(out, mul3(radiance, single_scatter));
            out = add3(out, mul3(multi_scatter, cosine_irradiance));
            out = add3(out, mul3(mul3(diffuse_color, energy), cosine_irradiance));
        }
        out = [
            out[0].clamp(0.01, 0.99),
            out[1].clamp(0.01, 0.99),
            out[2].clamp(0.01, 0.99),
        ];
    }
    out = scale3(out, renderer.lighting.exposure);
    [out[0], out[1], out[2], alpha]
}

/// Mol* points and wide lines use `assign_material_color` directly and never
/// include the mesh/sphere/cylinder lighting chunk. Exposure still applies at
/// the renderer level and WBOIT consumes the same pre-quantized float value.
pub(super) fn shade_unlit_material_linear(
    material: MeshMaterial,
    renderer: &RendererOptions,
) -> [f32; 4] {
    let base = color_f32(material.color);
    let exposure = renderer.lighting.exposure;
    [
        base[0] * exposure,
        base[1] * exposure,
        base[2] * exposure,
        material.alpha_tenths.min(10) as f32 / 10.0,
    ]
}

fn brdf_ggx(
    light_dir: Vec3,
    view_dir: Vec3,
    normal: Vec3,
    f0: [f32; 3],
    roughness: f32,
) -> [f32; 3] {
    let alpha = roughness * roughness;
    let alpha2 = alpha * alpha;
    let half_dir = (light_dir + view_dir).normalized();
    let dot_nl = normal.dot(light_dir).clamp(0.0, 1.0);
    let dot_nv = normal.dot(view_dir).clamp(0.0, 1.0);
    let dot_nh = normal.dot(half_dir).clamp(0.0, 1.0);
    let dot_vh = view_dir.dot(half_dir).clamp(0.0, 1.0);
    let fresnel = 2.0f32.powf((-5.55473 * dot_vh - 6.98316) * dot_vh);
    let f = [
        f0[0] * (1.0 - fresnel) + fresnel,
        f0[1] * (1.0 - fresnel) + fresnel,
        f0[2] * (1.0 - fresnel) + fresnel,
    ];
    let gv = dot_nl * (alpha2 + (1.0 - alpha2) * dot_nv * dot_nv).sqrt();
    let gl = dot_nv * (alpha2 + (1.0 - alpha2) * dot_nl * dot_nl).sqrt();
    let visibility = 0.5 / (gv + gl).max(1.0e-6);
    let denominator = dot_nh * dot_nh * (alpha2 - 1.0) + 1.0;
    let distribution = alpha2 / (std::f32::consts::PI * denominator * denominator).max(1.0e-6);
    scale3(f, visibility * distribution)
}

fn ggx_multiscattering(
    normal: Vec3,
    view_dir: Vec3,
    specular_color: [f32; 3],
    roughness: f32,
) -> ([f32; 3], [f32; 3]) {
    let dot_nv = normal.dot(view_dir).clamp(0.0, 1.0);
    let rx = 1.0 - roughness;
    let ry = 0.0425 - 0.0275 * roughness;
    let rz = 1.04 - 0.572 * roughness;
    let rw = -0.04 + 0.022 * roughness;
    let a004 = (rx * rx).min(2.0f32.powf(-9.28 * dot_nv)) * rx + ry;
    let fab_x = -1.04 * a004 + rz;
    let fab_y = 1.04 * a004 + rw;
    let single = [
        specular_color[0] * fab_x + fab_y,
        specular_color[1] * fab_x + fab_y,
        specular_color[2] * fab_x + fab_y,
    ];
    let energy_multiple = 1.0 - fab_x - fab_y;
    let f_avg = [
        specular_color[0] + (1.0 - specular_color[0]) * 0.047619,
        specular_color[1] + (1.0 - specular_color[1]) * 0.047619,
        specular_color[2] + (1.0 - specular_color[2]) * 0.047619,
    ];
    let multi = [
        single[0] * f_avg[0] / (1.0 - energy_multiple * f_avg[0]) * energy_multiple,
        single[1] * f_avg[1] / (1.0 - energy_multiple * f_avg[1]) * energy_multiple,
        single[2] * f_avg[2] / (1.0 - energy_multiple * f_avg[2]) * energy_multiple,
    ];
    (single, multi)
}

pub(super) fn light_direction(light: &DirectionalLight) -> Vec3 {
    if (light.inclination - 150.0).abs() < f64::EPSILON
        && (light.azimuth - 320.0).abs() < f64::EPSILON
    {
        return MOLSTAR_LIGHT;
    }
    // Mol* computes these values in JavaScript Number precision and converts
    // them to float32 only when `uniform3fv` uploads the regular number array.
    // Parsing and evaluating the angles as float32 first creates visible ULP
    // errors, especially at the cardinal directions.
    let radians_per_degree = std::f64::consts::PI / 180.0;
    let inclination = light.inclination * radians_per_degree;
    let azimuth = light.azimuth * radians_per_degree;
    Vec3::new(
        (azimuth.cos() * inclination.sin()) as f32,
        (azimuth.sin() * inclination.sin()) as f32,
        inclination.cos() as f32,
    )
}

pub(super) fn molstar_scaled_uniform_color(color: u32, intensity: f64) -> [f32; 3] {
    // Color.toVec3Normalized and Vec3.scale run in JavaScript Number
    // precision. WebGL's uniform upload performs the sole float32 conversion.
    [
        ((((color >> 16) & 0xff) as f64 / 255.0) * intensity) as f32,
        ((((color >> 8) & 0xff) as f64 / 255.0) * intensity) as f32,
        (((color & 0xff) as f64 / 255.0) * intensity) as f32,
    ]
}

fn add3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn mul3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] * b[0], a[1] * b[1], a[2] * b[2]]
}

fn scale3(a: [f32; 3], value: f32) -> [f32; 3] {
    [a[0] * value, a[1] * value, a[2] * value]
}

#[cfg(test)]
mod tests {
    #![allow(clippy::field_reassign_with_default)]

    use super::*;

    fn style(ignore_light: bool) -> ResolvedStyle {
        ResolvedStyle {
            ignore_light,
            occlusion: false,
            outline: false,
        }
    }

    #[test]
    fn default_light_matches_molstar_spherical_direction() {
        let direction = light_direction(&DirectionalLight::default());
        assert_eq!(direction.x.to_bits(), MOLSTAR_LIGHT.x.to_bits());
        assert_eq!(direction.y.to_bits(), MOLSTAR_LIGHT.y.to_bits());
        assert_eq!(direction.z.to_bits(), MOLSTAR_LIGHT.z.to_bits());
    }

    #[test]
    fn custom_light_directions_match_javascript_number_then_uniform3fv_staging() {
        let cases = [
            (180.0, 0.0, [0x250d_3132, 0x0000_0000, 0xbf80_0000]),
            (90.0, 90.0, [0x248d_3132, 0x3f80_0000, 0x248d_3132]),
            (
                33.333_333_333_333,
                271.234_567_890_123,
                [0x3c41_fa7b, 0xbf0c_a443, 0x3f55_e287],
            ),
            (
                12.345_678_901_234,
                -359.876_543_210_987,
                [0x3e5a_f0b0, 0x39f1_8a27, 0x3f7a_1482],
            ),
            (67.125, 123.875, [0xbf03_77db, 0x3f43_d46f, 0x3ec7_0691]),
        ];
        for (inclination, azimuth, expected) in cases {
            let mut light = DirectionalLight::default();
            light.inclination = inclination;
            light.azimuth = azimuth;
            let direction = light_direction(&light);
            assert_eq!(
                [
                    direction.x.to_bits(),
                    direction.y.to_bits(),
                    direction.z.to_bits(),
                ],
                expected,
                "inclination={inclination}, azimuth={azimuth}"
            );
        }
    }

    #[test]
    fn light_color_matches_javascript_scale_then_uniform3fv_staging() {
        let color = molstar_scaled_uniform_color(0x123456, 0.123_456_789_012_3);
        assert_eq!(
            color.map(f32::to_bits),
            [0x3c0e_c7ab, 0x3cce_3cdb, 0x3d2a_8af0]
        );
    }

    #[test]
    fn default_light_illuminates_the_molstar_front_fragment_hemisphere() {
        let material = MeshMaterial::opaque(0x1b9e77);
        let renderer = RendererOptions::default();
        let lit = shade_material(
            material,
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, -10.0),
            &renderer,
            style(false),
        );
        // This is the pinned apply-light-color.glsl result after RGBA8
        // quantization. Before reproducing the fragment-normal sign, the same
        // front fragment received only 40% ambient light: [11, 63, 48, 255].
        assert_eq!(lit, [26, 147, 111, 255]);
    }

    #[test]
    fn lighting_matches_metal_rgba8_across_normals_and_materials() {
        let cases = [
            (
                "front-green-r1-m0",
                0x1b9e77,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, -10.0),
                1.0,
                0.0,
                [26, 147, 111, 255],
            ),
            (
                "tilt-a-green-r1-m0",
                0x1b9e77,
                Vec3::new(0.45, 0.10, 0.89),
                Vec3::new(-0.25, 0.10, -8.0),
                1.0,
                0.0,
                [22, 124, 94, 255],
            ),
            (
                "tilt-b-green-r075-m0",
                0x1b9e77,
                Vec3::new(-0.55, 0.20, 0.81),
                Vec3::new(0.40, -0.15, -12.0),
                0.75,
                0.0,
                [30, 159, 121, 255],
            ),
            (
                "tilt-c-orange-r05-m0",
                0xd55e00,
                Vec3::new(0.15, -0.70, 0.70),
                Vec3::new(-0.20, 0.30, -6.0),
                0.50,
                0.0,
                [127, 56, 3, 255],
            ),
            (
                "tilt-d-blue-r025-m0",
                0x0072b2,
                Vec3::new(-0.75, -0.20, 0.63),
                Vec3::new(0.10, 0.20, -9.0),
                0.25,
                0.0,
                [3, 98, 153, 255],
            ),
            (
                "tilt-e-yellow-rmin-m0",
                0xf0e442,
                Vec3::new(0.30, 0.80, 0.52),
                Vec3::new(-0.50, -0.25, -11.0),
                0.0525,
                0.0,
                [181, 172, 50, 255],
            ),
            (
                "tilt-f-purple-r09-m025",
                0xcc79a7,
                Vec3::new(-0.15, 0.35, 0.92),
                Vec3::new(0.35, -0.45, -7.0),
                0.90,
                0.25,
                [171, 102, 140, 255],
            ),
            (
                "tilt-g-sky-r06-m05",
                0x56b4e9,
                Vec3::new(0.60, -0.10, 0.79),
                Vec3::new(-0.15, 0.25, -13.0),
                0.60,
                0.50,
                [40, 81, 104, 255],
            ),
            (
                "tilt-h-green-r035-m075",
                0x009e73,
                Vec3::new(-0.40, -0.45, 0.80),
                Vec3::new(0.20, 0.15, -5.0),
                0.35,
                0.75,
                [3, 67, 49, 255],
            ),
            (
                "front-orange-r015-m099",
                0xe69f00,
                Vec3::new(0.05, 0.05, 1.0),
                Vec3::new(-0.05, -0.10, -10.0),
                0.15,
                0.99,
                [90, 62, 3, 255],
            ),
            (
                "front-black-r1-m0",
                0x000000,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, -10.0),
                1.0,
                0.0,
                [3, 3, 3, 255],
            ),
            (
                "back-white-r1-m0",
                0xffffff,
                Vec3::new(0.0, 0.0, -1.0),
                Vec3::new(0.0, 0.0, -10.0),
                1.0,
                0.0,
                [102, 102, 102, 255],
            ),
        ];

        for (name, color, normal, view_position, roughness, metalness, expected) in cases {
            let mut renderer = RendererOptions::default();
            renderer.shading.material.roughness = roughness;
            renderer.shading.material.metalness = metalness;
            assert_eq!(
                shade_material(
                    MeshMaterial::opaque(color),
                    normal,
                    view_position,
                    &renderer,
                    style(false),
                ),
                expected,
                "{name}"
            );
        }
    }

    #[test]
    fn ignore_light_remains_the_exact_material_rgba() {
        let material = MeshMaterial::with_alpha_tenths(0x1b9e77, 6);
        assert_eq!(
            shade_material(
                material,
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(0.0, 0.0, -10.0),
                &RendererOptions::default(),
                style(true),
            ),
            [0x1b, 0x9e, 0x77, 153]
        );
    }
}
