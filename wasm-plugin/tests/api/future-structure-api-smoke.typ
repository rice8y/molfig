// Compile-time contract for BinaryCIF, assembly, altLoc, cartoon/ribbon, and
// render-object support.

#import "../../../package/lib.typ" as molfig

#let pdb = read("../fixtures/pdb/assembly-altloc-helix.pdb", encoding: none)
#let cif = read("../fixtures/cif/assembly-altloc-helix.cif", encoding: none)
#let bcif = read("../fixtures/bcif/assembly-altloc-helix.bcif", encoding: none)

#let base = molfig.info(cif, format: "cif", alt-loc: "A", assembly: "asymmetric-unit")
#let biological = molfig.info(cif, format: "cif", alt-loc: "A", assembly: "1")
#let binary = molfig.info(bcif, format: "bcif", alt-loc: "A", assembly: "1")

#assert.eq(base.atom_count, 13)
#assert.eq(base.alt_locs, ("A", "B"))
#assert.eq(base.assemblies.first().id, "1")
#assert.eq(biological.atom_count, 26)
#assert.eq(binary.atom_count, biological.atom_count)
#assert(binary.render_objects.any(object => object.secondary_type == "helix"))
#assert(binary.render_objects.any(object => object.geometry_type == "tube"))

#let cartoon-obj = molfig.to-obj(
  cif,
  format: "cif",
  representation: "cartoon",
  alt-loc: "highest-occupancy",
  assembly: "1",
)

#let ribbon-ply = molfig.to-ply(
  pdb,
  format: "pdb",
  representation: "ribbon",
  alt-loc: "A",
  assembly: "asymmetric-unit",
)

#assert(str(cartoon-obj).contains("\nv "))
#assert(str(ribbon-ply).starts-with("ply\n"))

#let object = molfig.render-object(
  bcif,
  format: "bcif",
  representation: "cartoon",
  alt-loc: "A",
  assembly: "1",
  mesh-format: "obj",
  config: (
    azimuth: 30,
    elevation: 18,
    background: "",
  ),
)

#assert.eq(object.format, "obj")
#assert(object.info.render_objects.any(item => item.representation == "cartoon"))
#assert(object.mesh.len() > 0)
#assert(object.content != none)
