// Contract test for passing maquette config dictionaries through unchanged.

#import "../../../package/lib.typ" as molfig
#import "@preview/maquette:0.2.0": get-obj-info, get-stl-info, get-ply-info, render-obj, render-stl, render-ply

#let water-cif = read("../fixtures/cif/water.cif", encoding: none)
#let mesh-options = (
  format: "cif",
  representation: "spacefill",
  sphere-detail: 1,
  center: false,
  assembly: "asymmetric-unit",
)
#let passthrough-config = (
  azimuth: 17,
  elevation: 23,
  background: "",
  lights: (
    ambient: 0.31,
    directional: 0.79,
  ),
)

#for (mesh-format, export-mesh, get-info, render-mesh) in (
  ("obj", molfig.to-obj, get-obj-info, render-obj),
  ("stl", molfig.to-stl, get-stl-info, render-stl),
  ("ply", molfig.to-ply, get-ply-info, render-ply),
) {
  let render-options = mesh-options + (mesh-format: mesh-format)
  let raw-mesh = export-mesh(water-cif, ..mesh-options)
  let direct-info = get-info(raw-mesh, passthrough-config)
  let changed-info = get-info(raw-mesh, passthrough-config + (azimuth: 71))
  let wrapped-info = molfig.mesh-info(water-cif, config: passthrough-config, ..render-options)

  assert(direct-info != changed-info)
  assert.eq(wrapped-info, direct-info)

  for output-format in ("svg", "png") {
    let object = molfig.render-object(
      water-cif,
      output-format: output-format,
      config: passthrough-config,
      width: 32mm,
      height: 24mm,
      ..render-options,
    )

    assert.eq(object.kind, "render-object")
    assert.eq(object.mesh, raw-mesh)
    assert.eq(object.info.atom_count, molfig.info(water-cif, ..mesh-options).atom_count)
    assert(object.content != none)

    // SVG is eager, so compare the actual image including generated OBJ colors.
    if output-format == "svg" {
      let config = if mesh-format == "obj" {
        passthrough-config + (materials: object.materials)
      } else {
        passthrough-config
      }
      assert.eq(object.content, render-mesh(raw-mesh, config, format: "svg", width: 32mm, height: 24mm))
    }

    // maquette 0.2.0 defers PNG rendering until layout. Place both public
    // render paths in the document so compilation executes their rasterizers.
    object.content
    molfig.render(
      water-cif,
      output-format: output-format,
      config: passthrough-config,
      width: 32mm,
      height: 24mm,
      ..render-options,
    )
  }
}
