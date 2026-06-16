# BinaryCIF Fixtures

`assembly-altloc-helix.bcif` is the BinaryCIF counterpart of
`../cif/assembly-altloc-helix.cif`. It covers the same small contract surface:

- `atom_site` with alternate locations, occupancies, auth/entity ids, and
  insertion codes;
- `struct_conf` helix metadata for cartoon/ribbon tests, including boundary
  insertion codes;
- `struct_sheet_range` sheet metadata with boundary insertion codes;
- `pdbx_struct_assembly*` categories with two operators for assembly expansion;
- `struct_conn` covalent links for deterministic smoke tests.

Regenerate it from the repository root with:

```sh
node wasm-plugin/tests/fixtures/bcif/generate-bcif-fixture.mjs
```

The `.sha256` file is checked by `wasm-plugin/tests/validate-fixtures.mjs` so
accidental binary fixture drift is visible even before a BinaryCIF parser lands
in the Rust/WASM implementation.

Every text mmCIF fixture in `../cif` has a same-stem BinaryCIF counterpart in
this directory. Most are generated directly from their text fixture by
`generate-bcif-fixture.mjs`, while `assembly-altloc-helix.bcif` stays
hand-expanded to keep extra Mol* typed-column coverage. The Rust regression
suite checks that each counterpart preserves atoms, bonds, coarse spheres,
coarse gaussians, and all text-fixture categories.

`water.bcif` is the smallest BinaryCIF fixture that mirrors `../cif/water.cif`.
It is useful for minimal `format: "bcif"` smoke tests before the
assembly/altLoc path is enabled.

`ihm-only.bcif` is a coarse-only BinaryCIF fixture with one
`_ihm_sphere_obj_site` row and one `_ihm_gaussian_obj_site` row. It is used by
the Typst API smoke tests to verify IHM parsing without any `_atom_site` rows.
