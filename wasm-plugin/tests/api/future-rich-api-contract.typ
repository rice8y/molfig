// Contract test for the richer public Typst API: BinaryCIF, assembly selection,
// altLoc policy, polymer-cartoon/ribbon representations, and native rendering.

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
  representation: "cartoon",
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

#let polymer-cartoon = molfig.to-ply(
  rich-cif,
  format: "mmcif",
  representation: "polymer-cartoon",
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

#assert(str(polymer-cartoon).starts-with("ply\n"))
#assert(str(ribbon).contains("\nv "))

#let result = molfig.render-result(
  rich-cif,
  format: "mmcif",
  representation: "polymer-cartoon",
  assembly: "1",
  alt-loc: "highest-occupancy",
  helix-profile: "rounded",
  round-cap: true,
  sheet-arrow-factor: 0.8,
  tubular-helices: true,
  linear-segments: 6,
  radial-segments: 12,
  renderer: (
    viewport: (width: 112, height: 84),
    camera: (view: (name: "orbit", params: (azimuth: 30, elevation: 18)),),
    background: (color: "#ffffff", transparent: true),
  ),
  width: 54mm,
  height: 42mm,
)

#assert.eq(result.kind, "render-result")
#assert.eq(result.pixel-width, 112)
#assert.eq(result.info.render_objects, result.render-info.render_objects)
#assert(result.info.render_objects.any(item => item.semantic_geometry_types.contains("tube")))
#assert(result.info.render_objects.any(item => item.visual == "polymer-trace"))
#assert(result.info.render_objects.any(item => item.group_count >= 1))
#assert(result.info.render_objects.any(item => item.instance_count >= 1))
#assert(result.info.render_objects.any(item => item.draw_count > 0))
#assert.eq(result.info.representation.name, "polymer-cartoon")
#assert.eq(result.info.representation.selected_visuals, ("polymer-trace",))
#assert.eq(result.info.representation.realized_visuals, ("polymer-trace",))
#assert(result.render-info.triangle_count > 0)
#assert(result.content != none)
