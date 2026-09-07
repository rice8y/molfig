use crate::mesh::{NativePrimitive, NativeRenderGeometry, NativeRenderScene};
use crate::options::{MeshOptions, RenderStyle};

use super::camera::CameraState;
use super::color::parse_color;
use super::multisample::jitter_offsets;
use super::options::{OcclusionParams, OutlineParams, RendererOptions};
use super::postprocessing::{fog_range, ssao_target_dimensions};
use super::shading::{light_direction, molstar_scaled_uniform_color};
use super::style::ResolvedStyle;

pub(super) fn render_info_json(
    scene: &NativeRenderScene,
    options: &MeshOptions,
    renderer: &RendererOptions,
    style: ResolvedStyle,
    camera: &CameraState,
    width: u32,
    height: u32,
) -> String {
    let style_name = if options.style == RenderStyle::Illustrative {
        "illustrative"
    } else {
        "default"
    };
    let view_matrix = camera.view_matrix();
    let projection_matrix = camera.projection_matrix();
    let inverse_projection_matrix = camera.inverse_projection_matrix();
    let projection_view_matrix = camera.projection_view_matrix();
    let staged_projection_view_matrix = multiply_mat4_column_major(projection_matrix, view_matrix);
    let directional = renderer
        .lighting
        .directional
        .iter()
        .map(|light| {
            let direction = light_direction(light);
            let color = parse_color(&light.color).unwrap_or(0xffffff);
            serde_json::json!({
                "color": light.color,
                "intensity": light.intensity,
                "inclination": light.inclination,
                "azimuth": light.azimuth,
                "uniform_direction": [direction.x, direction.y, direction.z],
                "uniform_color_times_intensity": molstar_scaled_uniform_color(color, light.intensity),
            })
        })
        .collect::<Vec<_>>();
    let ambient_color = parse_color(&renderer.lighting.ambient.color).unwrap_or(0xffffff);
    let occlusion_defaults = OcclusionParams::default();
    let occlusion = renderer
        .postprocessing
        .occlusion
        .as_ref()
        .map(|pass| &pass.params)
        .unwrap_or(&occlusion_defaults);
    let outline_defaults = OutlineParams::default();
    let outline = renderer
        .postprocessing
        .outline
        .as_ref()
        .map(|pass| &pass.params)
        .unwrap_or(&outline_defaults);
    let antialiasing = renderer
        .postprocessing
        .antialiasing
        .clone()
        .unwrap_or_default();
    let (fog_near, fog_far) = fog_range(camera, renderer.camera.fog.params.intensity);
    let scene_transparency_min = transparency_min(scene);
    let transparent_ssao =
        style.occlusion && scene_transparency_min < occlusion.transparent_threshold;
    let jitter_offsets = jitter_offsets(&renderer.multi_sample)
        .iter()
        .map(|&(x, y)| [x, y])
        .collect::<Vec<_>>();
    let render_objects = native_render_objects_value(scene);

    let viewport = serde_json::json!({
        "width": renderer.viewport.width,
        "height": renderer.viewport.height,
        "pixel_ratio": renderer.viewport.pixel_ratio,
        "drawing_buffer_width": width,
        "drawing_buffer_height": height,
    });
    let resolved_style = serde_json::json!({
        "ignore_light": style.ignore_light,
        "occlusion": style.occlusion,
        "outline": style.outline,
        "shadow": false,
        "antialiasing": antialiasing.name,
    });
    let background = serde_json::json!({
        "color": renderer.background.color,
        "transparent": renderer.background.transparent,
    });
    let material = serde_json::json!({
        "metalness": renderer.shading.material.metalness,
        "roughness": renderer.shading.material.roughness,
        "bumpiness": renderer.shading.material.bumpiness,
    });
    let shading = serde_json::json!({
        "ignore_light": style.ignore_light,
        "material": material,
    });
    let ambient = serde_json::json!({
        "color": renderer.lighting.ambient.color,
        "intensity": renderer.lighting.ambient.intensity,
        "uniform_color_times_intensity": molstar_scaled_uniform_color(
            ambient_color,
            renderer.lighting.ambient.intensity,
        ),
    });
    let lighting = serde_json::json!({
        "exposure": renderer.lighting.exposure,
        "ambient": ambient,
        "directional": directional,
    });
    let transparency = serde_json::json!({
        "mode": renderer.transparency.mode,
        "scene_transparency_min": scene_transparency_min,
        "transparent_ssao": transparent_ssao,
    });
    let multi_sample = serde_json::json!({
        "mode": renderer.multi_sample.mode,
        "sample_level": renderer.multi_sample.sample_level,
        "sample_count": jitter_offsets.len(),
        "reuse_occlusion": renderer.multi_sample.reuse_occlusion,
        "jitter_offsets": jitter_offsets,
    });
    let (ssao_target_width, ssao_target_height, ssao_target_scale) = ssao_target_dimensions(
        width as usize,
        height as usize,
        renderer.viewport.pixel_ratio,
        occlusion.resolution_scale,
    );
    let occlusion_params = serde_json::json!({
        "samples": occlusion.samples,
        "multi_scale": {
            "name": occlusion.multi_scale.name,
            "params": {
                "levels": occlusion.multi_scale.params.levels.iter().map(|level| {
                    serde_json::json!({
                        "radius": level.radius,
                        "bias": level.bias,
                    })
                }).collect::<Vec<_>>(),
                "near_threshold": occlusion.multi_scale.params.near_threshold,
                "far_threshold": occlusion.multi_scale.params.far_threshold,
            },
        },
        "radius": occlusion.radius,
        "bias": occlusion.bias,
        "blur_kernel_size": occlusion.blur_kernel_size,
        "blur_depth_bias": occlusion.blur_depth_bias,
        "resolution_scale": occlusion.resolution_scale,
        "color": occlusion.color,
        "transparent_threshold": occlusion.transparent_threshold,
        "target": {
            "width": ssao_target_width,
            "height": ssao_target_height,
            "scale": ssao_target_scale,
        },
    });
    let outline_params = serde_json::json!({
        "color": outline.color,
        "scale": outline.scale,
        "threshold": outline.threshold,
        "include_transparent": outline.include_transparent,
        "pixel_threshold": 50.0 * outline.threshold * renderer.viewport.pixel_ratio,
        "pixel_scale": (outline.scale * renderer.viewport.pixel_ratio).round().max(1.0) - 1.0,
    });
    let antialiasing_params = serde_json::json!({
        "edge_threshold": antialiasing.params.edge_threshold,
        "max_search_steps": antialiasing.params.max_search_steps,
    });
    let postprocessing = serde_json::json!({
        "occlusion": {
            "name": if style.occlusion { "on" } else { "off" },
            "params": occlusion_params,
        },
        "outline": {
            "name": if style.outline { "on" } else { "off" },
            "params": outline_params,
        },
        "shadow": { "name": "off" },
        "antialiasing": {
            "name": antialiasing.name,
            "params": antialiasing_params,
        },
    });
    let clipping = serde_json::json!({
        "far": renderer.camera.clipping.far,
        "min_near": renderer.camera.clipping.min_near,
        "min_far": renderer.camera.clipping.min_far,
        "force_full": false,
    });
    let fog = serde_json::json!({
        "name": renderer.camera.fog.name,
        "intensity": renderer.camera.fog.params.intensity,
        "near": fog_near,
        "far": fog_far,
    });
    let camera_info = serde_json::json!({
        "mode": renderer.camera.mode,
        "view": renderer.camera.view.name,
        "fov": renderer.camera.fov,
        "fov_degrees": renderer.camera.fov,
        "fov_radians_uniform": camera.fov,
        "aspect": camera.aspect,
        "position": [camera.position.x, camera.position.y, camera.position.z],
        "target": [camera.target.x, camera.target.y, camera.target.z],
        "up": [camera.up.x, camera.up.y, camera.up.z],
        "right": [camera.right.x, camera.right.y, camera.right.z],
        "forward": [camera.forward.x, camera.forward.y, camera.forward.z],
        "distance": camera.distance,
        "radius": camera.radius,
        "radius_max": camera.radius_max,
        "near": camera.near,
        "far": camera.far,
        "clipping": clipping,
        "fog": fog,
        "matrix_layout": "column-major-webgl-uniform",
        "matrix_number_format": "f32-shortest-roundtrip-decimal",
        "view_matrix": view_matrix,
        "projection_matrix": projection_matrix,
        "inverse_projection_matrix": inverse_projection_matrix,
        "projection_view_matrix": projection_view_matrix,
        "staged_projection_view_matrix": staged_projection_view_matrix,
    });

    serde_json::json!({
        "renderer": "molstar-native-cpu",
        "bundle_version": 2,
        "render_info_version": 2,
        "molstar_commit": crate::MOLSTAR_REFERENCE_COMMIT,
        "pixel_width": width,
        "pixel_height": height,
        "viewport": viewport,
        "triangle_count": scene.mesh.faces.len(),
        "analytic_primitive_count": scene.primitives.len(),
        "render_object_count": scene.render_objects.len(),
        "render_objects": render_objects,
        "style": style_name,
        "resolved_style": resolved_style,
        "background": background,
        "shading": shading,
        "lighting": lighting,
        "transparency": transparency,
        "multi_sample": multi_sample,
        "postprocessing": postprocessing,
        "camera": camera_info,
    })
    .to_string()
}

pub(crate) fn native_render_objects_json(scene: &NativeRenderScene) -> String {
    serde_json::Value::Array(native_render_objects_value(scene)).to_string()
}

fn native_render_objects_value(scene: &NativeRenderScene) -> Vec<serde_json::Value> {
    scene
        .render_objects
        .iter()
        .map(|object| {
            let geometry = match &object.geometry {
                NativeRenderGeometry::Mesh { ranges } => serde_json::json!({
                    "kind": "mesh",
                    "ranges": ranges.iter().map(|range| serde_json::json!({
                        "vertex_start": range.vertices.start,
                        "vertex_end": range.vertices.end,
                        "face_start": range.faces.start,
                        "face_end": range.faces.end,
                    })).collect::<Vec<_>>(),
                }),
                NativeRenderGeometry::Primitives { kind, indices } => {
                    let spans = contiguous_index_spans(indices);
                    serde_json::json!({
                        "kind": kind,
                        "primitive_count": indices.len(),
                        "primitive_spans": spans.iter().map(|&(start, end)| serde_json::json!({
                            "start": start,
                            "end": end,
                        })).collect::<Vec<_>>(),
                    })
                }
            };
            let material = object.material.map(|material| {
                serde_json::json!({
                    "color": format!("#{:06x}", material.color),
                    "alpha": f64::from(material.alpha_tenths) / 10.0,
                })
            });
            let bounding_sphere = object.bounding_sphere.as_ref().map(|sphere| {
                serde_json::json!({
                    "center": [sphere.center.x, sphere.center.y, sphere.center.z],
                    "center64": sphere.center64(),
                    "radius": sphere.radius,
                    "radius64": sphere.radius64(),
                    "extrema_count": sphere.extrema.len().max(sphere.extrema64.len()),
                })
            });
            serde_json::json!({
                "order": object.order,
                "geometry_type": object.geometry_kind,
                "semantic_geometry_types": object.semantic_geometry_kinds,
                "visual": object.visual,
                "representation": object.representation,
                "representation_order": object.representation_order,
                "color_theme": object.color_theme,
                "component": object.component,
                "tag": object.tag,
                "geometry": geometry,
                "draw_count": object.draw_count,
                "vertex_count": object.vertex_count,
                "group_count": object.group_count,
                "instance_count": object.instance_count,
                "material": material,
                "bounding_sphere": bounding_sphere,
                "alpha": object.alpha,
                "alpha_factor": object.alpha_factor,
                "alpha_min": object.alpha_min,
                "alpha_average": object.alpha_average,
                "transparency_average": object.transparency_average,
                "transparency_min": object.transparency_min,
                "state": {
                    "visible": object.visible,
                    "opaque": object.opaque,
                    "write_depth": object.write_depth,
                    "double_sided": object.double_sided,
                    "flip_sided": object.flip_sided,
                    "backface_transparency": object.backface_transparency,
                },
            })
        })
        .collect()
}

fn contiguous_index_spans(indices: &[usize]) -> Vec<(usize, usize)> {
    let Some((&first, rest)) = indices.split_first() else {
        return Vec::new();
    };
    let mut spans = Vec::new();
    let mut start = first;
    let mut end = first.saturating_add(1);
    for &index in rest {
        if index == end {
            end = end.saturating_add(1);
        } else {
            spans.push((start, end));
            start = index;
            end = index.saturating_add(1);
        }
    }
    spans.push((start, end));
    spans
}

fn transparency_min(scene: &NativeRenderScene) -> f64 {
    scene
        .mesh
        .face_materials
        .iter()
        .map(|material| material.alpha_tenths)
        .chain(scene.primitives.iter().map(|primitive| match primitive {
            NativePrimitive::Point { material, .. }
            | NativePrimitive::Line { material, .. }
            | NativePrimitive::Sphere { material, .. }
            | NativePrimitive::Cylinder { material, .. } => material.alpha_tenths,
        }))
        .filter(|&alpha| alpha < 10)
        .map(|alpha| 1.0 - f64::from(alpha) / 10.0)
        .fold(1.0, f64::min)
}

fn multiply_mat4_column_major(a: [f32; 16], b: [f32; 16]) -> [f32; 16] {
    let mut out = [0.0; 16];
    for column in 0..4 {
        for row in 0..4 {
            // Mol* matrices use ordinary JavaScript arrays. Matrix products
            // therefore run as Number arithmetic and are converted to f32
            // only when WebGL uploads the uniform.
            out[column * 4 + row] = (0..4)
                .map(|index| f64::from(a[index * 4 + row]) * f64::from(b[column * 4 + index]))
                .sum::<f64>() as f32;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{BoundingSphere, Vec3};
    use crate::render::camera::resolve_camera;

    #[test]
    fn reported_matrices_share_webgl_column_major_layout() {
        let renderer = RendererOptions::default();
        let camera = resolve_camera(
            &renderer,
            &BoundingSphere {
                center: Vec3::new(1.0, 2.0, 3.0),
                radius: 4.0,
                ..BoundingSphere::default()
            },
            800,
            600,
        )
        .unwrap();
        let view = camera.view_matrix();
        let projection = camera.projection_matrix();
        let combined = multiply_mat4_column_major(projection, view);
        let point = [2.0, -1.0, 0.5, 1.0];
        let combined_point = transform_column_major(combined, point);
        let staged_point = transform_column_major(projection, transform_column_major(view, point));
        for (combined, staged) in combined_point.into_iter().zip(staged_point) {
            assert!((combined - staged).abs() <= 0.000_004);
        }
        assert_eq!(view[15], 1.0);
        assert_eq!(projection[11], -1.0);
    }

    #[test]
    fn column_major_product_preserves_identity_on_both_sides() {
        let identity = [
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ];
        let matrix = [
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
        ];
        assert_eq!(multiply_mat4_column_major(identity, matrix), matrix);
        assert_eq!(multiply_mat4_column_major(matrix, identity), matrix);
    }

    #[test]
    fn json_matrix_decimals_round_trip_every_float32_bit() {
        let renderer = RendererOptions::from_json(
            br#"{"camera":{"mode":"perspective","fov":37.125,"view":{"name":"orbit","params":{"azimuth":123.456789,"elevation":-27.75,"roll":8.5}}}}"#,
        )
        .unwrap();
        let camera = resolve_camera(
            &renderer,
            &BoundingSphere {
                center: Vec3::new(-12.25, 0.03125, 87.5),
                radius: 33.125,
                ..BoundingSphere::default()
            },
            1003,
            719,
        )
        .unwrap();
        for matrix in [
            camera.view_matrix(),
            camera.projection_matrix(),
            camera.inverse_projection_matrix(),
            camera.projection_view_matrix(),
            multiply_mat4_column_major(camera.projection_matrix(), camera.view_matrix()),
        ] {
            let encoded = serde_json::to_string(&matrix).unwrap();
            let decoded: [f32; 16] = serde_json::from_str(&encoded).unwrap();
            assert_eq!(decoded.map(f32::to_bits), matrix.map(f32::to_bits));
        }
    }

    #[test]
    fn primitive_report_spans_are_half_open_and_compact() {
        assert_eq!(contiguous_index_spans(&[]), Vec::<(usize, usize)>::new());
        assert_eq!(contiguous_index_spans(&[7]), vec![(7, 8)]);
        assert_eq!(
            contiguous_index_spans(&[0, 1, 2, 8, 9, 4, 11]),
            vec![(0, 3), (8, 10), (4, 5), (11, 12)]
        );
    }

    fn transform_column_major(matrix: [f32; 16], value: [f32; 4]) -> [f32; 4] {
        [0, 1, 2, 3].map(|row| {
            (0..4)
                .map(|column| matrix[column * 4 + row] * value[column])
                .sum()
        })
    }
}
