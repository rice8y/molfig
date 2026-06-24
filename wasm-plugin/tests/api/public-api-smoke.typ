// Compile-time smoke test for the public Typst API.

#import "../../../package/lib.typ" as molfig

#let water-pdb = read("../fixtures/pdb/water.pdb", encoding: none)
#let water-cif = read("../fixtures/cif/water.cif", encoding: none)
#let annotation-cif = read("../fixtures/cif/viewer-default-annotations.cif", encoding: none)
#let annotation-bcif = read("../fixtures/bcif/viewer-default-annotations.bcif", encoding: none)
#let sheet-pdb = "SHEET    1 S1  1 SER A   1  THR A   4\nATOM      1  CA  SER A   1       0.000   0.000   0.000  1.00 10.00           C\nATOM      2  CA  THR A   2       1.100   0.250   0.200  1.00 10.00           C\nATOM      3  CA  SER A   3       2.200  -0.200   0.100  1.00 10.00           C\nATOM      4  CA  THR A   4       3.300   0.150   0.000  1.00 10.00           C\nEND\n"

#let pdb-info = molfig.info(water-pdb, format: "pdb")
#let cif-info = molfig.info(water-cif, format: "cif")
#let auto-info = molfig.info(water-pdb, format: "pdb", representation: "auto")
#let ball-info = molfig.info(water-pdb, format: "pdb", representation: "ball-and-stick")
#let spacefill-info = molfig.info(water-cif, format: "cif", representation: "spacefill")
#let surface-info = molfig.info(water-cif, format: "cif", representation: "surface")
#let annotation-info = molfig.info(annotation-cif, format: "cif", representation: "default")
#let annotation-bcif-info = molfig.info(annotation-bcif, format: "bcif", representation: "default")
#let qmean-info = molfig.info(
  annotation-cif,
  format: "cif",
  representation: "auto",
  color-theme: "qmean-score",
)
#let sheet-info = molfig.info(
  sheet-pdb,
  format: "pdb",
  representation: "cartoon",
  infer-bonds: false,
  center: false,
  sheet-arrow-factor: 0.75,
)
#let polymer-cartoon-info = molfig.info(
  sheet-pdb,
  format: "pdb",
  representation: "polymer-cartoon",
  infer-bonds: false,
  center: false,
)
#assert.eq(pdb-info.atom_count, 3)
#assert.eq(pdb-info.bond_count, 2)
#assert.eq(cif-info.atom_count, pdb-info.atom_count)
#assert.eq(pdb-info.representation.name, "default")
#assert.eq(auto-info.representation.name, "auto")
#assert.eq(pdb-info.representation.selected_visuals, ("element-sphere", "intra-bond", "inter-bond"))
#assert.eq(auto-info.representation.selected_visuals, pdb-info.representation.selected_visuals)
#assert.eq(ball-info.representation.name, "ball-and-stick")
#assert.eq(ball-info.representation.selected_visuals, ("element-sphere", "intra-bond", "inter-bond"))
#assert(ball-info.representation.realized_visuals.any(visual => visual == "element-sphere"))
#assert(ball-info.representation.realized_visuals.any(visual => visual == "intra-bond"))
#assert(not ball-info.representation.realized_visuals.any(visual => visual == "inter-bond"))
#assert.eq(spacefill-info.representation.name, "spacefill")
#assert.eq(spacefill-info.representation.selected_visuals, ("element-sphere",))
#assert.eq(spacefill-info.representation.realized_visuals, ("element-sphere",))
#assert.eq(surface-info.representation.name, "surface")
#assert.eq(surface-info.representation.selected_visuals, ("molecular-surface-mesh",))
#assert.eq(surface-info.representation.realized_visuals, ("molecular-surface-mesh",))
#assert(surface-info.render_objects.all(object => object.color_theme == "entity-id"))
#assert.eq(annotation-info.atom_count, 4)
#assert.eq(annotation-bcif-info.atom_count, annotation-info.atom_count)
#assert(annotation-info.render_objects.any(object => object.color_theme == "plddt-confidence"))
#assert(annotation-bcif-info.render_objects.any(object => object.color_theme == "plddt-confidence"))
#assert(qmean-info.render_objects.any(object => object.color_theme == "qmean-score"))
#assert.eq(sheet-info.representation.name, "cartoon")
#assert.eq(polymer-cartoon-info.representation.name, "polymer-cartoon")
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

#let themed-info = molfig.info(
  water-pdb,
  format: "pdb",
  representation: "cartoon",
  theme: (
    globalName: "element-symbol",
    carbonColor: "chain-id",
    symmetryColor: "operator-name",
  ),
  sphere-detail: 1,
)
#assert(themed-info.render_objects.all(object => object.color_theme == "element-symbol"))
#assert(themed-info.render_objects.all(object => object.carbon_color_theme == "element-symbol"))

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
#assert.eq(colored-object.info, molfig.info(
  water-pdb,
  format: "pdb",
  representation: "spacefill",
  color-theme: "chain-id",
  sphere-detail: 1,
))
#assert.eq(colored-object.mesh, molfig.to-obj(
  water-pdb,
  format: "pdb",
  representation: "spacefill",
  color-theme: "chain-id",
  sphere-detail: 1,
))

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
