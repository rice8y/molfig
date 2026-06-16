// Compile-time smoke test for IHM-only coarse CIF input.

#import "../../../package/lib.typ" as molfig

#let sphere-cif = read("../fixtures/cif/ihm-sphere-only.cif", encoding: none)
#let gaussian-cif = read("../fixtures/cif/ihm-gaussian-only.cif", encoding: none)

#let sphere-info = molfig.info(
  sphere-cif,
  format: "cif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert.eq(sphere-info.atom_count, 0)
#assert.eq(sphere-info.ihm_model_count, 1)
#assert.eq(sphere-info.coarse_sphere_count, 1)
#assert.eq(sphere-info.coarse_gaussian_count, 0)
#assert.eq(sphere-info.structure.unit_kind_counts.atomic, 0)
#assert.eq(sphere-info.structure.unit_kind_counts.spheres, 1)
#assert.eq(sphere-info.structure.unit_kind_counts.gaussians, 0)
#assert(sphere-info.render_objects.any(object => object.secondary_type == "coarse-sphere"))

#let sphere-obj = molfig.to-obj(
  sphere-cif,
  format: "cif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert(str(sphere-obj).contains("\nv "))
#assert(str(sphere-obj).contains("\nf "))

#let gaussian-info = molfig.info(
  gaussian-cif,
  format: "cif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert.eq(gaussian-info.atom_count, 0)
#assert.eq(gaussian-info.ihm_model_count, 1)
#assert.eq(gaussian-info.coarse_sphere_count, 0)
#assert.eq(gaussian-info.coarse_gaussian_count, 1)
#assert.eq(gaussian-info.structure.unit_kind_counts.atomic, 0)
#assert.eq(gaussian-info.structure.unit_kind_counts.spheres, 0)
#assert.eq(gaussian-info.structure.unit_kind_counts.gaussians, 1)
#assert(gaussian-info.render_objects.any(object => object.secondary_type == "coarse-gaussian"))

#let gaussian-obj = molfig.to-obj(
  gaussian-cif,
  format: "cif",
  center: false,
  assembly: "asymmetric-unit",
)

#assert(str(gaussian-obj).contains("\nv "))
#assert(str(gaussian-obj).contains("\nf "))
