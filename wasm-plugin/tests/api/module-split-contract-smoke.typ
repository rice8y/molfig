// Compile-time regression smoke for a Rust module split.
// The goal is to keep parser normalization, mesh export, metadata, and
// maquette-facing render-object behavior stable while implementation files move.

#import "../../../package/lib.typ" as molfig

#let peptide-pdb = read("../fixtures/pdb/tiny-peptide.pdb", encoding: none)
#let peptide-cif = read("../fixtures/cif/tiny-peptide.cif", encoding: none)
#let water-bcif = read("../fixtures/bcif/water.bcif", encoding: none)

#let pdb-info = molfig.info(peptide-pdb, format: "pdb", assembly: "asymmetric-unit")
#let cif-info = molfig.info(peptide-cif, format: "mmcif", assembly: "asymmetric-unit")

#assert.eq(pdb-info.atom_count, 6)
#assert.eq(cif-info.atom_count, pdb-info.atom_count)
#assert.eq(cif-info.bond_count, pdb-info.bond_count)
#assert.eq(cif-info.bounds.min, pdb-info.bounds.min)
#assert.eq(cif-info.bounds.max, pdb-info.bounds.max)

#let bcif-info = molfig.info(water-bcif, format: "binarycif", assembly: "asymmetric-unit")
#assert.eq(bcif-info.atom_count, 3)
#assert.eq(bcif-info.bond_count, 2)

#let options = (
  representation: "ball-and-stick",
  sphere-detail: 1,
  center: true,
  assembly: "asymmetric-unit",
)

#let obj = molfig.to-obj(peptide-cif, format: "mmcif", ..options)
#let stl = molfig.to-stl(peptide-cif, format: "mmcif", ..options)
#let ply = molfig.to-ply(peptide-cif, format: "mmcif", ..options)

#assert(str(obj).contains("\nv "))
#assert(str(obj).contains("\nf "))
#assert(stl.len() > 84)
#assert(str(ply).starts-with("ply\n"))
#assert(str(ply).contains("element vertex"))
#assert(not str(obj).contains("nan"))
#assert(not str(ply).contains("nan"))

#let object = molfig.render-object(
  peptide-cif,
  format: "mmcif",
  mesh-format: "obj",
  config: (
    azimuth: 20,
    elevation: 15,
    background: "",
  ),
  width: 40mm,
  height: 32mm,
  ..options,
)

#assert.eq(object.kind, "render-object")
#assert.eq(object.mesh_format, "obj")
#assert.eq(object.info.atom_count, cif-info.atom_count)
#assert(object.mesh.len() > 0)
#assert(object.content != none)
