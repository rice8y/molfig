// Compile-time smoke test for the public Typst API.

#import "../../../package/lib.typ" as molfig

#let water-pdb = read("../fixtures/pdb/water.pdb", encoding: none)
#let water-cif = read("../fixtures/cif/water.cif", encoding: none)
#let sheet-pdb = "SHEET    1 S1  1 SER A   1  THR A   4\nATOM      1  CA  SER A   1       0.000   0.000   0.000  1.00 10.00           C\nATOM      2  CA  THR A   2       1.100   0.250   0.200  1.00 10.00           C\nATOM      3  CA  SER A   3       2.200  -0.200   0.100  1.00 10.00           C\nATOM      4  CA  THR A   4       3.300   0.150   0.000  1.00 10.00           C\nEND\n"

#let pdb-info = molfig.info(water-pdb, format: "pdb")
#let cif-info = molfig.info(water-cif, format: "cif")
#let ball-info = molfig.info(water-pdb, format: "pdb", representation: "ball-and-stick")
#let spacefill-info = molfig.info(water-cif, format: "cif", representation: "spacefill")
#let sheet-info = molfig.info(
  sheet-pdb,
  format: "pdb",
  representation: "molstar",
  infer-bonds: false,
  center: false,
  sheet-arrow-factor: 0.75,
)

#assert.eq(pdb-info.atom_count, 3)
#assert.eq(pdb-info.bond_count, 2)
#assert.eq(cif-info.atom_count, pdb-info.atom_count)
#assert.eq(ball-info.representation.name, "ball-and-stick")
#assert.eq(ball-info.representation.selected_visuals, ("element-sphere", "intra-bond", "inter-bond"))
#assert(ball-info.representation.realized_visuals.any(visual => visual == "element-sphere"))
#assert(ball-info.representation.realized_visuals.any(visual => visual == "intra-bond"))
#assert(not ball-info.representation.realized_visuals.any(visual => visual == "inter-bond"))
#assert.eq(spacefill-info.representation.name, "spacefill")
#assert.eq(spacefill-info.representation.selected_visuals, ("element-sphere",))
#assert.eq(spacefill-info.representation.realized_visuals, ("element-sphere",))
#assert(sheet-info.render_objects.any(object => object.geometry_type == "sheet"))

#if molfig.v15-or-later() {
  let water-pdb-path = path("../fixtures/pdb/water.pdb")
  let path-info = molfig.info(water-pdb-path, format: "pdb")
  let path-obj = molfig.to-obj(
    water-pdb-path,
    format: "pdb",
    representation: "ball-and-stick",
    sphere-detail: 1,
  )
  assert.eq(path-info.atom_count, pdb-info.atom_count)
  assert.eq(path-info.bond_count, pdb-info.bond_count)
  assert(str(path-obj).contains("\nv "))
}

#let obj = molfig.to-obj(
  water-pdb,
  format: "pdb",
  representation: "ball-and-stick",
  sphere-detail: 1,
  helix-profile: "rounded",
  round-cap: true,
  sheet-arrow-factor: 1.0,
  tubular-helices: true,
  linear-segments: 6,
  radial-segments: 12,
)

#assert(str(obj).contains("\nv "))
#assert(str(obj).contains("\nf "))

#let colored-object = molfig.render-object(
  water-pdb,
  format: "pdb",
  representation: "spacefill",
  color-theme: "chain-id",
  mesh-format: "obj",
  sphere-detail: 1,
  width: 24mm,
  height: 24mm,
  config: (
    background: "",
    materials: (
      "0xff0d0d6": "#123456",
    ),
  ),
)

#assert(colored-object.materials.len() > 0)
#assert(colored-object.materials.values().all(color => color.starts-with("#")))

#let stl = molfig.to-stl(water-pdb, format: "pdb", sphere-detail: 1, round-cap: true)
#let ply = molfig.to-ply(water-cif, format: "cif", sphere-detail: 1, sheet-arrow-factor: 0.5)

#assert(stl.len() > 84)
#assert(str(ply).starts-with("ply\n"))

#let mesh-meta = molfig.mesh-info(water-cif, format: "cif", mesh-format: "ply", sphere-detail: 1)
#assert(mesh-meta != none)

#let rendered = molfig.render(
  water-cif,
  format: "cif",
  representation: "spacefill",
  mesh-format: "ply",
  helix-profile: "square",
  round-cap: true,
  sheet-arrow-factor: 0.6,
  tubular-helices: false,
  linear-segments: 8,
  radial-segments: 16,
  width: 42mm,
  height: 34mm,
  config: (
    azimuth: 25,
    elevation: 18,
    background: "",
  ),
)

#assert(rendered != none)
