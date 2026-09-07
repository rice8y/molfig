//! Native renderer pipeline orchestration and temporary pass implementations.
//!
//! The renderer consumes Molfig's semantic scene result directly. Geometry is
//! never serialized to OBJ, STL, or PLY on the rendering path. The current
//! indexed pass is deliberately self-contained so camera, depth, material,
//! postprocessing, and byte quantization stay under one reproducible contract.

use super::camera::resolve_camera;
use super::color::parse_color;
use super::framebuffer::Framebuffer;
use super::multisample::{accumulate, jitter_offsets, resolve};
use super::options::{OcclusionParams, OutlineParams, RendererOptions};
use super::postprocessing::{
    apply_occlusion_factors, apply_outline, apply_smaa, apply_transparent_outline,
    compute_occlusion_factors, compute_occlusion_factors_including_transparency, read_color,
    write_composed_color,
};
use super::raster::{rasterize_scene, RasterPass};
use super::report::render_info_json;
use super::style::resolve_style;
use super::transparency::WboitTarget;
use crate::mesh::NativeRenderScene;
use crate::model::{BoundingSphere, Mesh};
use crate::options::MeshOptions;

#[derive(Clone, Debug)]
pub(crate) struct RenderedImage {
    pub(crate) width: u32,
    pub(crate) height: u32,
    pub(crate) rgba: Vec<u8>,
    pub(crate) render_info_json: String,
}

struct OcclusionCache {
    opaque: super::postprocessing::OcclusionFactors,
    transparent: Option<super::postprocessing::OcclusionFactors>,
}

pub(crate) fn render_scene(
    scene: &NativeRenderScene,
    geometry_options: &MeshOptions,
    renderer: &RendererOptions,
) -> Result<RenderedImage, String> {
    let width = ((f64::from(renderer.viewport.width) * renderer.viewport.pixel_ratio).round()
        as u32)
        .max(1);
    let height = ((f64::from(renderer.viewport.height) * renderer.viewport.pixel_ratio).round()
        as u32)
        .max(1);
    let sphere = scene
        .visible_bounding_sphere
        .clone()
        .unwrap_or_else(|| mesh_bounding_sphere(&scene.mesh));
    let camera = resolve_camera(renderer, &sphere, width, height)?;
    let style = resolve_style(geometry_options, renderer);
    let background_rgb = parse_color(&renderer.background.color)?;
    let background = [
        (background_rgb >> 16) as u8,
        (background_rgb >> 8) as u8,
        background_rgb as u8,
        if renderer.background.transparent {
            0
        } else {
            255
        },
    ];
    let offsets = jitter_offsets(&renderer.multi_sample);
    let mut accumulation = vec![[0.0f32; 4]; width as usize * height as usize];
    let mut reused_occlusion = None;
    for (sample_index, &(offset_x, offset_y)) in offsets.iter().enumerate() {
        let sample_camera = camera.with_view_offset(offset_x, offset_y);
        let framebuffer = render_sample(
            scene,
            renderer,
            style,
            background,
            background_rgb,
            &sample_camera,
            width as usize,
            height as usize,
            (offset_x, offset_y),
            &mut reused_occlusion,
        )?;
        accumulate(
            &mut accumulation,
            &framebuffer,
            &renderer.multi_sample.mode,
            sample_index,
            offsets.len(),
        );
    }
    let rgba = resolve(accumulation);
    let render_info_json = render_info_json(
        scene,
        geometry_options,
        renderer,
        style,
        &camera,
        width,
        height,
    );
    Ok(RenderedImage {
        width,
        height,
        rgba,
        render_info_json,
    })
}

#[allow(clippy::too_many_arguments)]
fn render_sample(
    scene: &NativeRenderScene,
    renderer: &RendererOptions,
    style: super::style::ResolvedStyle,
    background: [u8; 4],
    background_rgb: u32,
    camera: &super::camera::CameraState,
    width: usize,
    height: usize,
    occlusion_offset: (f32, f32),
    reused_occlusion: &mut Option<OcclusionCache>,
) -> Result<Framebuffer, String> {
    let mut framebuffer = Framebuffer::new(width, height, background);
    let mut opaque_pass = RasterPass::opaque();
    rasterize_scene(
        scene,
        camera,
        renderer,
        style,
        &mut framebuffer,
        &mut opaque_pass,
    );
    let mut transparency = scene_has_transparency(scene).then(|| WboitTarget::new(width, height));
    if let Some(target) = transparency.as_mut() {
        let mut transparent_pass = RasterPass::wboit(target);
        rasterize_scene(
            scene,
            camera,
            renderer,
            style,
            &mut framebuffer,
            &mut transparent_pass,
        );
    }
    let mut transparent_layer = if style.occlusion || style.outline {
        transparency
            .as_ref()
            .map(WboitTarget::evaluated_framebuffer)
    } else {
        None
    };
    let mut post_color = if style.occlusion || style.outline {
        read_color(&framebuffer)
    } else {
        Vec::new()
    };
    let mut transparent_color = transparent_layer.as_ref().map(read_color);
    if style.occlusion {
        let occlusion_defaults = OcclusionParams::default();
        let occlusion_params = renderer
            .postprocessing
            .occlusion
            .as_ref()
            .map(|pass| &pass.params)
            .unwrap_or(&occlusion_defaults);
        let include_transparency = transparency.is_some()
            && scene_transparency_min(scene) < occlusion_params.transparent_threshold;
        if renderer.multi_sample.reuse_occlusion && reused_occlusion.is_some() {
            let cached = reused_occlusion.as_ref().expect("checked occlusion cache");
            apply_occlusion_factors(
                &framebuffer,
                &mut post_color,
                &cached.opaque,
                occlusion_params,
                camera,
                renderer.camera.fog.params.intensity,
                background_rgb,
                renderer.background.transparent,
                occlusion_offset,
            )?;
            if let (Some(layer), Some(factors)) =
                (transparent_layer.as_mut(), cached.transparent.as_ref())
            {
                apply_occlusion_factors(
                    layer,
                    transparent_color
                        .as_mut()
                        .expect("transparent color target"),
                    factors,
                    occlusion_params,
                    camera,
                    renderer.camera.fog.params.intensity,
                    background_rgb,
                    true,
                    occlusion_offset,
                )?;
            }
        } else {
            let (factors, transparent_factors) = if include_transparency {
                let transparent_depth_alpha = transparency
                    .as_ref()
                    .expect("checked transparent target")
                    .packed_nearest_depth_alpha_rgba8();
                let (opaque, transparent) = compute_occlusion_factors_including_transparency(
                    &framebuffer,
                    &transparent_depth_alpha,
                    camera,
                    occlusion_params,
                    renderer.viewport.pixel_ratio,
                );
                (opaque, Some(transparent))
            } else {
                (
                    compute_occlusion_factors(
                        &framebuffer,
                        camera,
                        occlusion_params,
                        renderer.viewport.pixel_ratio,
                    ),
                    None,
                )
            };
            apply_occlusion_factors(
                &framebuffer,
                &mut post_color,
                &factors,
                occlusion_params,
                camera,
                renderer.camera.fog.params.intensity,
                background_rgb,
                renderer.background.transparent,
                (0.0, 0.0),
            )?;
            if let (Some(layer), Some(transparent_factors)) =
                (transparent_layer.as_mut(), transparent_factors.as_ref())
            {
                apply_occlusion_factors(
                    layer,
                    transparent_color
                        .as_mut()
                        .expect("transparent color target"),
                    transparent_factors,
                    occlusion_params,
                    camera,
                    renderer.camera.fog.params.intensity,
                    background_rgb,
                    true,
                    (0.0, 0.0),
                )?;
            }
            if renderer.multi_sample.reuse_occlusion {
                *reused_occlusion = Some(OcclusionCache {
                    opaque: factors,
                    transparent: transparent_factors,
                });
            }
        }
    }
    if style.outline {
        let outline_defaults = OutlineParams::default();
        let outline_params = renderer
            .postprocessing
            .outline
            .as_ref()
            .map(|pass| &pass.params)
            .unwrap_or(&outline_defaults);
        let opaque_outline_depth01 = apply_outline(
            &framebuffer,
            &mut post_color,
            camera,
            outline_params,
            renderer.viewport.pixel_ratio,
            renderer.camera.fog.params.intensity,
            background_rgb,
            renderer.background.transparent,
        )?;
        if let (Some(target), Some(layer)) = (transparency.as_ref(), transparent_layer.as_mut()) {
            let (_, transparent_alpha) = target.packed_nearest_depth_and_alpha();
            apply_transparent_outline(
                &framebuffer,
                layer,
                transparent_color
                    .as_mut()
                    .expect("transparent color target"),
                &transparent_alpha,
                &opaque_outline_depth01,
                camera,
                outline_params,
                renderer.viewport.pixel_ratio,
                renderer.camera.fog.params.intensity,
            )?;
        }
    }
    if style.occlusion || style.outline {
        write_composed_color(&mut framebuffer, &post_color, transparent_color.as_deref());
    } else if let Some(transparency) = transparency {
        transparency.compose(&mut framebuffer, false);
    }
    let antialiasing = renderer
        .postprocessing
        .antialiasing
        .clone()
        .unwrap_or_default();
    if antialiasing.name == "smaa" {
        apply_smaa(&mut framebuffer, &antialiasing.params);
    }
    Ok(framebuffer)
}

fn scene_has_transparency(scene: &NativeRenderScene) -> bool {
    scene
        .render_objects
        .iter()
        .filter(|object| object.visible)
        .any(|object| {
            let alpha = (object.alpha * object.alpha_factor).clamp(0.0, 1.0);
            alpha > 0.0 && (alpha < 1.0 || object.transparency_average > 0.0)
        })
}

fn scene_transparency_min(scene: &NativeRenderScene) -> f64 {
    let mut minimum = 1.0f64;
    for object in scene.render_objects.iter().filter(|object| object.visible) {
        let alpha = (object.alpha * object.alpha_factor).clamp(0.0, 1.0);
        if alpha < 1.0 {
            minimum = minimum.min(1.0 - alpha);
        }
        if object.transparency_min > 0.0 {
            minimum = minimum.min(object.transparency_min);
        }
    }
    minimum
}

#[cfg(test)]
fn transparency_min_for_alpha_tenths(values: impl IntoIterator<Item = u8>) -> f64 {
    values
        .into_iter()
        .filter(|&alpha| alpha < 10)
        .map(|alpha| 1.0 - f64::from(alpha) / 10.0)
        .fold(1.0, f64::min)
}

fn mesh_bounding_sphere(mesh: &Mesh) -> BoundingSphere {
    crate::model::Boundary::from_positions(&mesh.vertices).sphere
}

#[cfg(test)]
mod tests {
    #![allow(clippy::neg_cmp_op_on_partial_ord)]

    use super::transparency_min_for_alpha_tenths;

    #[test]
    fn transparent_ssao_threshold_uses_molstar_javascript_number_comparison() {
        let threshold = 0.4f64;
        let transparency = transparency_min_for_alpha_tenths([10, 3, 6]);
        assert_eq!(transparency, 1.0 - 0.6f64);
        assert!(!(transparency < threshold));
        assert!(transparency_min_for_alpha_tenths([7]) < threshold);
        assert_eq!(transparency_min_for_alpha_tenths([10]), 1.0);
    }
}
