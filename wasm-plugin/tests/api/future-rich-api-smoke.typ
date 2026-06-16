// Smoke test for the richer public Typst API: BinaryCIF, assembly selection,
// altLoc policy, cartoon/ribbon representations, and render-object output.

#import "../../../package/lib.typ" as molfig

#let rich-cif = read("../fixtures/cif/assembly-altloc-secondary.cif", encoding: none)
#let rich-pdb = read("../fixtures/pdb/assembly-altloc-secondary.pdb", encoding: none)
#let water-bcif = read("../fixtures/bcif/water.bcif", encoding: none)

#let bcif-info = molfig.info(water-bcif, format: "bcif")
#assert.eq(bcif-info.atom_count, 3)
#assert.eq(bcif-info.bond_count, 2)

#let rich-info = molfig.info(
  rich-cif,
  format: "mmcif",
  assembly: "1",
  alt-loc: "highest-occupancy",
)

#assert.eq(rich-info.atom_count, 11)
#assert.eq(rich-info.assembly.id, "1")
#assert.eq(rich-info.alt_locs_info.policy, "highest-occupancy")
#assert.eq(rich-info.secondary_structure.helices.len(), 1)
#assert.eq(rich-info.secondary_structure.sheets.len(), 1)
#assert(rich-info.render_objects.any(object => object.secondary_type == "helix"))
#assert(rich-info.render_objects.any(object => object.polymer_trace.sec_struc_first))
#assert(rich-info.render_objects.any(object => object.geometry_type == "sheet" and object.secondary_type == "sheet"))

#let cartoon = molfig.to-ply(
  rich-cif,
  format: "mmcif",
  representation: "cartoon",
  assembly: "1",
  alt-loc: "highest-occupancy",
  sphere-detail: 1,
  helix-profile: "rounded",
  round-cap: true,
  sheet-arrow-factor: 0.8,
  tubular-helices: true,
  linear-segments: 6,
  radial-segments: 12,
)

#let ribbon = molfig.to-obj(
  rich-pdb,
  format: "pdb",
  representation: "ribbon",
  assembly: "1",
  alt-loc: "A",
  sphere-detail: 1,
)

#assert(str(cartoon).starts-with("ply\n"))
#assert(str(ribbon).contains("\nv "))

#let object = molfig.render-object(
  rich-cif,
  format: "mmcif",
  representation: "cartoon",
  mesh-format: "ply",
  assembly: "1",
  alt-loc: "highest-occupancy",
  helix-profile: "rounded",
  round-cap: true,
  sheet-arrow-factor: 0.8,
  tubular-helices: true,
  linear-segments: 6,
  radial-segments: 12,
  config: (
    azimuth: 30,
    elevation: 18,
    background: "",
  ),
  width: 54mm,
  height: 42mm,
)

#assert.eq(object.kind, "render-object")
#assert.eq(object.format, "ply")
#assert.eq(object.mesh_format, "ply")
#assert.eq(object.mesh, molfig.to-ply(
  rich-cif,
  format: "mmcif",
  representation: "cartoon",
  assembly: "1",
  alt-loc: "highest-occupancy",
  helix-profile: "rounded",
  round-cap: true,
  sheet-arrow-factor: 0.8,
  tubular-helices: true,
  linear-segments: 6,
  radial-segments: 12,
))
#assert(object.info.render_objects.any(item => item.geometry_type == "tube"))
#assert(object.info.render_objects.any(item => item.visual == "polymer-trace"))
#assert(object.info.render_objects.any(item => item.value_cell.u_group_count >= 1))
#assert(object.info.render_objects.any(item => item.valueCell.uGroupCount >= 1))
#assert(object.info.render_objects.any(item => item.valueCell.drawCount > 0))
#assert.eq(object.info.representation.name, "cartoon")
#assert.eq(object.info.representation.selected_visuals, ("polymer-trace",))
#assert.eq(object.info.representation.realized_visuals, ("polymer-trace",))
#assert(object.info.representation.selected_visuals.any(visual => visual == "polymer-trace"))
#assert(object.info.representation.realized_visuals.any(visual => visual == "polymer-trace"))
#assert(object.content != none)
