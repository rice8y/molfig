// Compile-time smoke test for IHM-only coarse BinaryCIF input.

#import "../../../package/lib.typ" as molfig

#let ihm-bcif = read("../fixtures/bcif/ihm-only.bcif", encoding: none)

#let info = molfig.info(
  ihm-bcif,
  format: "bcif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert.eq(info.atom_count, 0)
#assert.eq(info.ihm_model_count, 1)
#assert.eq(info.ihm_model_group_count, 1)
#assert.eq(info.ihm_model_group_link_count, 1)
#assert.eq(info.ihm_cross_link_restraint_count, 1)
#assert.eq(info.coarse_sphere_count, 1)
#assert.eq(info.coarse_gaussian_count, 1)
#assert.eq(info.structure.unit_kind_counts.atomic, 0)
#assert.eq(info.structure.unit_kind_counts.spheres, 1)
#assert.eq(info.structure.unit_kind_counts.gaussians, 1)

#let obj = molfig.to-obj(
  ihm-bcif,
  format: "bcif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert(str(obj).contains("\nv "))
#assert(str(obj).contains("\nf "))
