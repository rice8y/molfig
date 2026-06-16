// Compile-time smoke test for mixed atomic and coarse IHM CIF input.

#import "../../../package/lib.typ" as molfig

#let mixed-cif = read("../fixtures/cif/mixed-atomic-coarse-ihm.cif", encoding: none)

#let info = molfig.info(
  mixed-cif,
  format: "cif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert.eq(info.atom_count, 2)
#assert.eq(info.ihm_model_count, 2)
#assert.eq(info.coarse_sphere_count, 2)
#assert.eq(info.coarse_gaussian_count, 2)
#assert.eq(info.structure.model_count, 2)
#assert.eq(info.structure.unit_kind_counts.atomic, 1)
#assert.eq(info.structure.unit_kind_counts.spheres, 1)
#assert.eq(info.structure.unit_kind_counts.gaussians, 1)
#assert(info.render_objects.any(object => object.secondary_type == "coarse-sphere"))
#assert(info.render_objects.any(object => object.secondary_type == "coarse-gaussian"))

#let obj = molfig.to-obj(
  mixed-cif,
  format: "cif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert(str(obj).contains("\nv "))
#assert(str(obj).contains("\nf "))
