use super::super::camera::CameraState;
use super::super::color::color_f32;
use super::super::options::RendererOptions;

/// Apply the opaque branch of Mol*'s `apply_fog` fragment chunk while the
/// shaded color is still float32. The geometry shader performs this mix before
/// writing the RGBA8 color target; quantizing lighting first is observably
/// different for lit fragments.
pub(in crate::render) fn fog_opaque_fragment(
    mut color: [f32; 4],
    camera: &CameraState,
    renderer: &RendererOptions,
    depth01: f32,
) -> [f32; 4] {
    if renderer.camera.fog.name != "on" {
        return color;
    }
    let factor = fog_factor(camera, renderer.camera.fog.params.intensity, depth01);
    if renderer.background.transparent {
        // Opaque fragments enter this branch with alpha one. Mol* nevertheless
        // evaluates the general premultiplied-alpha expression literally.
        let fog_alpha = (1.0 - factor) * color[3];
        for channel in &mut color[..3] {
            *channel *= fog_alpha;
        }
        color[3] = fog_alpha;
    } else {
        let fog_color = color_f32(renderer.background.color_value);
        for (channel, fog) in fog_color.iter().enumerate() {
            color[channel] = color[channel] * (1.0 - factor) + fog * factor;
        }
    }
    color
}

pub(in crate::render) fn fog_factor(camera: &CameraState, intensity: f64, depth01: f32) -> f32 {
    let (fog_near, fog_far) = fog_range(camera, intensity);
    let denom = (fog_far - fog_near).max(0.0001);
    let depth = camera.depth_from_depth01(depth01);
    let normalized = ((depth - fog_near) / denom).clamp(0.0, 1.0);
    normalized * normalized * (3.0 - 2.0 * normalized)
}

pub(in crate::render) fn fog_range(camera: &CameraState, intensity: f64) -> (f32, f32) {
    camera.fog_range(intensity)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::BoundingSphere;
    use crate::render::camera::resolve_camera;
    use crate::render::color::quantize;

    #[test]
    fn opaque_fog_mixes_float_lighting_before_rgba8_attachment_conversion() {
        let renderer = RendererOptions::from_json(
            br##"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":2}},"fog":{"name":"on","params":{"intensity":50}}},"background":{"color":"#4080c0"}}"##,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 64, 64).unwrap();
        let depth01 = camera.impostor_depth01(-11.0);
        let factor = fog_factor(&camera, 50.0, depth01);
        assert_eq!(factor.to_bits(), 0.5f32.to_bits());

        let source = [0.003, 0.654_321, 0.333_333_34, 1.0];
        let fogged = fog_opaque_fragment(source, &camera, &renderer, depth01);
        let fog = color_f32(0x4080c0);
        let expected = [
            source[0] * 0.5 + fog[0] * 0.5,
            source[1] * 0.5 + fog[1] * 0.5,
            source[2] * 0.5 + fog[2] * 0.5,
            1.0,
        ];
        assert_eq!(fogged.map(f32::to_bits), expected.map(f32::to_bits));

        let attachment = fogged.map(quantize);
        let prematurely_quantized = source.map(quantize).map(|value| f32::from(value) / 255.0);
        let old_attachment = [
            prematurely_quantized[0] * 0.5 + fog[0] * 0.5,
            prematurely_quantized[1] * 0.5 + fog[1] * 0.5,
            prematurely_quantized[2] * 0.5 + fog[2] * 0.5,
            1.0,
        ]
        .map(quantize);
        assert_eq!(attachment, [32, 147, 138, 255]);
        assert_eq!(old_attachment, [33, 147, 138, 255]);
    }

    #[test]
    fn transparent_background_premultiplies_opaque_fog_before_quantization() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"view":{"name":"snapshot","params":{"position":[0,0,10],"target":[0,0,0],"up":[0,1,0],"radius":2,"radius-max":2}},"fog":{"name":"on","params":{"intensity":50}}},"background":{"transparent":true}}"#,
        )
        .unwrap();
        let camera = resolve_camera(&renderer, &BoundingSphere::default(), 64, 64).unwrap();
        let depth01 = camera.impostor_depth01(-11.0);
        let source = [0.2, 0.4, 0.8, 1.0];
        assert_eq!(
            fog_opaque_fragment(source, &camera, &renderer, depth01).map(f32::to_bits),
            [0.1f32, 0.2, 0.4, 0.5].map(f32::to_bits)
        );
    }
}
