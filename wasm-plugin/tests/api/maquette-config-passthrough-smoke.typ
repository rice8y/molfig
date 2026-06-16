// Contract smoke for passing maquette config dictionaries through unchanged.

#import "../../../package/lib.typ" as molfig
#import "@preview/maquette:0.1.0": get-ply-info

#let water-cif = read("../fixtures/cif/water.cif", encoding: none)
#let mesh-options = (
  format: "cif",
  representation: "spacefill",
  sphere-detail: 1,
  center: false,
  assembly: "asymmetric-unit",
)
#let render-options = mesh-options + (mesh-format: "ply")
#let passthrough-config = (
  azimuth: 17,
  elevation: 23,
  background: "",
  lights: (
    ambient: 0.31,
    directional: 0.79,
  ),
)

#let raw-ply = molfig.to-ply(water-cif, ..mesh-options)
#let direct-info = get-ply-info(raw-ply, json.encode(passthrough-config))
#let changed-info = get-ply-info(raw-ply, json.encode(passthrough-config + (azimuth: 71)))
#let wrapped-info = molfig.mesh-info(water-cif, config: passthrough-config, ..render-options)

#assert(direct-info != changed-info)
#assert.eq(wrapped-info, direct-info)

#let render-object = molfig.render-object(
  water-cif,
  output-format: "svg",
  config: passthrough-config,
  width: 32mm,
  height: 24mm,
  ..render-options,
)

#assert.eq(render-object.kind, "render-object")
#assert.eq(render-object.mesh, raw-ply)
#assert.eq(render-object.info.atom_count, molfig.info(water-cif, ..mesh-options).atom_count)
#assert(render-object.content != none)
