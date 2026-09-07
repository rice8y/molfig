//! Mol* weighted blended order-independent transparency.
//!
//! The pinned renderer accumulates transparent fragments in two float32
//! attachments and evaluates them into an RGBA8 color target. Keeping this as
//! a separate target is important: transparent fragments test only against
//! opaque depth and never hide one another through depth writes.

use super::camera::CameraState;
use super::color::quantize;
use super::framebuffer::Framebuffer;
use super::options::RendererOptions;
use super::postprocessing::{fog_factor, pack_depth_alpha_rgba8, unpack_depth_alpha_rgba};

pub(super) struct WboitTarget {
    width: usize,
    height: usize,
    accumulated_rgb: Vec<[f32; 3]>,
    accumulated_weight: Vec<f32>,
    revealage: Vec<f32>,
    nearest_depth01: Vec<f32>,
    nearest_alpha: Vec<f32>,
}

impl WboitTarget {
    pub(super) fn new(width: usize, height: usize) -> Self {
        let len = width.saturating_mul(height);
        Self {
            width,
            height,
            accumulated_rgb: vec![[0.0; 3]; len],
            accumulated_weight: vec![0.0; len],
            revealage: vec![1.0; len],
            nearest_depth01: vec![1.0; len],
            nearest_alpha: vec![0.0; len],
        }
    }

    pub(super) fn accumulate(
        &mut self,
        index: usize,
        fragment_depth: f32,
        mut color: [f32; 4],
        camera: &CameraState,
        renderer: &RendererOptions,
    ) {
        let pre_fog_alpha = color[3];
        if pre_fog_alpha >= 1.0 {
            return;
        }
        if fragment_depth < self.nearest_depth01[index] {
            self.nearest_depth01[index] = fragment_depth;
            self.nearest_alpha[index] = pre_fog_alpha;
        }

        if renderer.camera.fog.name == "on" {
            let fog = fog_factor(camera, renderer.camera.fog.params.intensity, fragment_depth);
            let fog_alpha = (1.0 - fog) * color[3];
            if renderer.background.transparent {
                for channel in &mut color[..3] {
                    *channel *= fog_alpha;
                }
            }
            color[3] = fog_alpha;
        } else if renderer.background.transparent {
            let alpha = color[3];
            for channel in &mut color[..3] {
                *channel *= alpha;
            }
        }

        let alpha = color[3];
        let depth_weight = ((1.0 - fragment_depth) * (1.0 - fragment_depth)).clamp(0.01, 1.0);
        let weight = alpha * depth_weight;
        for (channel, accumulated) in self.accumulated_rgb[index].iter_mut().enumerate() {
            *accumulated += color[channel] * alpha * weight;
        }
        self.accumulated_weight[index] += if renderer.background.transparent {
            alpha * alpha * weight
        } else {
            alpha * weight
        };
        self.revealage[index] *= 1.0 - alpha;
    }

    /// Evaluate WBOIT and compose it over the postprocessed opaque color.
    /// When Mol* postprocessing is enabled, evaluation first lands in its
    /// RGBA8 transparent color target; otherwise evaluation blends directly
    /// into the RGBA8 main color target.
    pub(super) fn compose(&self, framebuffer: &mut Framebuffer, intermediate_rgba8: bool) {
        if intermediate_rgba8 {
            let layer = self.evaluated_framebuffer();
            Self::compose_evaluated(framebuffer, &layer);
            return;
        }
        for index in 0..framebuffer.color.len() {
            let alpha = 1.0 - self.revealage[index];
            if alpha <= 0.0 {
                continue;
            }
            let denominator = self.accumulated_weight[index].clamp(0.000_000_01, 50_000.0);
            let evaluated = [
                self.accumulated_rgb[index][0] / denominator,
                self.accumulated_rgb[index][1] / denominator,
                self.accumulated_rgb[index][2] / denominator,
            ];
            let destination = framebuffer.color[index].map(|value| f32::from(value) / 255.0);
            let source_rgb = evaluated.map(|value| value * alpha);
            let source_alpha = alpha;
            for channel in 0..3 {
                framebuffer.color[index][channel] =
                    quantize(source_rgb[channel] + destination[channel] * (1.0 - source_alpha));
            }
            framebuffer.color[index][3] =
                quantize(source_alpha + destination[3] * (1.0 - source_alpha));
        }
    }

    pub(super) fn evaluated_framebuffer(&self) -> Framebuffer {
        let (depth01, _) = self.packed_nearest_depth_and_alpha();
        let color = (0..self.revealage.len())
            .map(|index| {
                let alpha = 1.0 - self.revealage[index];
                if alpha <= 0.0 {
                    return [0; 4];
                }
                let denominator = self.accumulated_weight[index].clamp(0.000_000_01, 50_000.0);
                [
                    quantize(self.accumulated_rgb[index][0] / denominator * alpha),
                    quantize(self.accumulated_rgb[index][1] / denominator * alpha),
                    quantize(self.accumulated_rgb[index][2] / denominator * alpha),
                    quantize(alpha),
                ]
            })
            .collect::<Vec<_>>();
        Framebuffer {
            width: self.width,
            height: self.height,
            color,
            depth: vec![f32::INFINITY; depth01.len()],
            depth01,
            normal: vec![Default::default(); self.revealage.len()],
        }
    }

    pub(super) fn compose_evaluated(destination: &mut Framebuffer, layer: &Framebuffer) {
        assert_eq!(destination.width, layer.width);
        assert_eq!(destination.height, layer.height);
        for index in 0..destination.color.len() {
            let alpha = f32::from(layer.color[index][3]) / 255.0;
            if alpha <= 0.0 {
                continue;
            }
            for channel in 0..3 {
                let source = f32::from(layer.color[index][channel]) / 255.0;
                let background = f32::from(destination.color[index][channel]) / 255.0;
                destination.color[index][channel] = quantize(source + background * (1.0 - alpha));
            }
            let background_alpha = f32::from(destination.color[index][3]) / 255.0;
            destination.color[index][3] = quantize(alpha + background_alpha * (1.0 - alpha));
        }
    }

    pub(super) fn packed_nearest_depth_and_alpha(&self) -> (Vec<f32>, Vec<f32>) {
        self.packed_nearest_depth_alpha_rgba8()
            .into_iter()
            .map(|packed| unpack_depth_alpha_rgba(packed.map(|channel| f32::from(channel) / 255.0)))
            .unzip()
    }

    /// Return the actual transparent-depth render-target bytes. Mol* clears
    /// this RGBA8 target to white, so background texels must remain `[255; 4]`
    /// through later linear filtering rather than being synthesized by
    /// packing a decoded `(depth: 1, alpha: 0)` pair.
    pub(super) fn packed_nearest_depth_alpha_rgba8(&self) -> Vec<[u8; 4]> {
        self.nearest_depth01
            .iter()
            .copied()
            .zip(self.nearest_alpha.iter().copied())
            .map(|(depth, alpha)| {
                if depth >= 1.0 {
                    [255; 4]
                } else {
                    pack_depth_alpha_rgba8(depth, alpha)
                }
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlapping_fragments_follow_molstar_wboit_equations() {
        let mut target = WboitTarget::new(1, 1);
        let mut renderer = RendererOptions::default();
        renderer.camera.fog.name = "off".into();
        let camera = test_camera(&renderer);
        target.accumulate(0, 0.25, [1.0, 0.0, 0.0, 0.5], &camera, &renderer);
        target.accumulate(0, 0.75, [0.0, 0.0, 1.0, 0.5], &camera, &renderer);

        assert_eq!(target.revealage[0], 0.25);
        assert_eq!(target.accumulated_rgb[0], [0.140625, 0.0, 0.015625]);
        assert_eq!(target.accumulated_weight[0], 0.15625);
        assert_eq!(target.nearest_depth01[0], 0.25);
        assert_eq!(target.nearest_alpha[0], 0.5);
        let mut framebuffer = Framebuffer::new(1, 1, [255; 4]);
        target.compose(&mut framebuffer, false);
        assert_eq!(framebuffer.color[0], [236, 64, 83, 255]);
    }

    #[test]
    fn transparent_background_uses_premultiplied_wboit_denominator() {
        let mut target = WboitTarget::new(1, 1);
        let mut renderer = RendererOptions::default();
        renderer.camera.fog.name = "off".into();
        renderer.background.transparent = true;
        let camera = test_camera(&renderer);
        target.accumulate(0, 0.5, [0.2, 0.4, 0.8, 0.5], &camera, &renderer);

        let mut framebuffer = Framebuffer::new(1, 1, [0; 4]);
        target.compose(&mut framebuffer, false);
        assert_eq!(framebuffer.color[0], [26, 51, 102, 128]);
    }

    #[test]
    fn transparent_depth_target_preserves_white_clear_texels() {
        let target = WboitTarget::new(2, 1);
        assert_eq!(target.packed_nearest_depth_alpha_rgba8(), vec![[255; 4]; 2]);
        let (depth, alpha) = target.packed_nearest_depth_and_alpha();
        // The shader's packed-depth dot product maps white to the largest
        // float below one. Do not normalize it to the ordinary depth sentinel.
        assert_eq!(depth, vec![f32::from_bits(0x3f7f_ffff); 2]);
        assert_eq!(alpha, vec![1.0; 2]);
    }

    fn test_camera(renderer: &RendererOptions) -> CameraState {
        let sphere = crate::model::Boundary::from_positions(&[
            crate::model::Vec3::new(-1.0, 0.0, 0.0),
            crate::model::Vec3::new(1.0, 0.0, 0.0),
        ])
        .sphere;
        super::super::camera::resolve_camera(renderer, &sphere, 1, 1).unwrap()
    }
}
