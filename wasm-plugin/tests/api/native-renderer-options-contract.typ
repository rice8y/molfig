// Contract test for strict Mol*-oriented native renderer options.

#import "../../../package/lib.typ" as molfig

#let water = read("../fixtures/xyz/water.xyz", encoding: none)
#let result = molfig.render-result(
  water,
  format: "xyz",
  representation: "ball-and-stick",
  color-theme: "element-symbol",
  style: "illustrative",
  output-format: "png",
  renderer: (
    viewport: (width: 80, height: 60, pixel-ratio: 1.5),
    camera: (
      mode: "orthographic",
      view: (name: "orbit", params: (azimuth: 35, elevation: 24, roll: 4)),
      fog: (name: "off"),
      clipping: (far: false, min-near: 0.3, min-far: 12.75),
    ),
    background: (color: "#f8f8f8", transparent: false),
    shading: (ignore-light: false, material: (metalness: 0, roughness: 0.4, bumpiness: 0)),
    lighting: (
      exposure: 1,
      ambient: (color: "#ffffff", intensity: 0.4),
      directional: ((color: "#ffffff", intensity: 0.6, inclination: 150, azimuth: 320),),
    ),
    transparency: (mode: "wboit"),
    multi-sample: (mode: "on", sample-level: 1, reuse-occlusion: false),
    postprocessing: (
      occlusion: (
        name: "on",
        params: (
          samples: 4,
          multi-scale: (
            name: "on",
            params: (
              levels: (
                (radius: 2, bias: 1),
                (radius: 5, bias: 0.75),
              ),
              near-threshold: 10,
              far-threshold: 1500,
            ),
          ),
        ),
      ),
      outline: (name: "on", params: (color: "#102030", scale: 2, threshold: 0.25)),
      shadow: (name: "off"),
      antialiasing: (name: "off"),
    ),
  ),
)

#assert.eq(result.pixel-width, 120)
#assert.eq(result.pixel-height, 90)
#assert.eq(result.pixels.len(), 120 * 90 * 4)
#assert.eq(result.output-format, "png")
#assert.eq(result.image.at(0), 137)
#assert.eq(str(result.image.slice(1, 4)), "PNG")
#assert.eq(result.render-info.camera.mode, "orthographic")
#assert.eq(result.render-info.render_info_version, 2)
#assert.eq(result.render-info.viewport.pixel_ratio, 1.5)
#assert.eq(result.render-info.viewport.drawing_buffer_width, 120)
#assert.eq(result.render-info.viewport.drawing_buffer_height, 90)
#assert(not result.render-info.resolved_style.ignore_light)
#assert(result.render-info.resolved_style.occlusion)
#assert(result.render-info.resolved_style.outline)
#assert.eq(result.render-info.resolved_style.antialiasing, "off")
#assert.eq(result.render-info.multi_sample.mode, "on")
#assert.eq(result.render-info.multi_sample.sample_level, 1)
#assert.eq(result.render-info.multi_sample.sample_count, 2)
#assert.eq(result.render-info.multi_sample.jitter_offsets, ((0, 0), (-0.25, -0.25)))
#assert(not result.render-info.multi_sample.reuse_occlusion)
#assert.eq(result.render-info.postprocessing.occlusion.name, "on")
#assert.eq(result.render-info.postprocessing.occlusion.params.multi_scale.name, "on")
#assert.eq(result.render-info.postprocessing.occlusion.params.multi_scale.params.levels.len(), 2)
#assert.eq(result.render-info.postprocessing.occlusion.params.multi_scale.params.levels.at(1).radius, 5)
#assert.eq(result.render-info.postprocessing.occlusion.params.multi_scale.params.levels.at(1).bias, 0.75)
#assert.eq(result.render-info.postprocessing.occlusion.params.multi_scale.params.near_threshold, 10)
#assert.eq(result.render-info.postprocessing.occlusion.params.multi_scale.params.far_threshold, 1500)
#assert.eq(result.render-info.postprocessing.occlusion.params.target.width, 80)
#assert.eq(result.render-info.postprocessing.occlusion.params.target.height, 60)
#assert.eq(result.render-info.postprocessing.outline.name, "on")
#assert.eq(result.render-info.postprocessing.outline.params.color, "#102030")
#assert.eq(result.render-info.postprocessing.outline.params.pixel_threshold, 18.75)
#assert.eq(result.render-info.postprocessing.outline.params.pixel_scale, 2)
#assert.eq(result.render-info.postprocessing.antialiasing.name, "off")
#assert.eq(result.render-info.camera.view, "orbit")
#assert.eq(result.render-info.camera.fog.name, "off")
#assert(not result.render-info.camera.clipping.far)
#assert.eq(result.render-info.camera.clipping.min_near, 0.3)
#assert.eq(result.render-info.camera.clipping.min_far, 12.75)
#assert(not result.render-info.camera.clipping.force_full)
#assert.eq(result.render-info.camera.matrix_layout, "column-major-webgl-uniform")
#assert.eq(result.render-info.camera.view_matrix.len(), 16)
#assert.eq(result.render-info.camera.projection_matrix.len(), 16)
#assert.eq(result.render-info.camera.inverse_projection_matrix.len(), 16)
#assert.eq(result.render-info.camera.projection_view_matrix.len(), 16)
#assert.eq(result.render-info.camera.staged_projection_view_matrix.len(), 16)
#assert.eq(result.render-info.lighting.directional.at(0).uniform_direction.len(), 3)
