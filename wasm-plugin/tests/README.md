# Molfig Test Design

This test area covers the Typst package surface for PDB/CIF to OBJ/STL/PLY to
maquette integration. It is scoped to fixtures, examples, docs, and API
contracts; implementation code should live outside this directory.

## Fixtures

The fixtures are intentionally small:

- `tests/fixtures/pdb/water.pdb` and `tests/fixtures/cif/water.cif` check the
  smallest useful molecule with explicit bonds.
- `tests/fixtures/pdb/tiny-peptide.pdb` and
  `tests/fixtures/cif/tiny-peptide.cif` check chain, residue, and backbone-like
  metadata.
- `tests/fixtures/pdb/atomic-protein-no-altloc.pdb` and
  `tests/fixtures/cif/atomic-protein-no-altloc.cif` check an atomic protein
  fragment with no alternate-location ids.
- `tests/fixtures/pdb/atomic-protein-altloc-tie.pdb` and
  `tests/fixtures/cif/atomic-protein-altloc-tie.cif` check deterministic A/B
  alternate-location selection when occupancies tie.
- `tests/fixtures/pdb/mixed-protein-nucleic.pdb` and
  `tests/fixtures/cif/mixed-protein-nucleic.cif` check that protein and
  nucleic-acid residues derive correctly in one structure.
- `tests/fixtures/cif/mixed-atomic-coarse-ihm.cif` checks that one IHM entry can
  carry atomic atoms, coarse spheres, and coarse gaussians across the same model
  list.
- `tests/fixtures/pdb/assembly-altloc-helix.pdb`,
  `tests/fixtures/cif/assembly-altloc-helix.cif`, and
  `tests/fixtures/bcif/assembly-altloc-helix.bcif` check the future BinaryCIF,
  biological assembly, alternate-location, and helix/cartoon/ribbon contract.
- `tests/fixtures/pdb/assembly-altloc-secondary.pdb` and
  `tests/fixtures/cif/assembly-altloc-secondary.cif` add a sheet range beside
  helix metadata so cartoon/ribbon smoke tests can cover both secondary
  structure families.
- `tests/fixtures/bcif/water.bcif` is a small BinaryCIF counterpart to the
  existing water fixture for minimal decoder smoke coverage.

Each PDB/CIF pair should produce the same parse summary. The exact mesh can
vary by tessellation quality, but the output must be deterministic for a fixed
input and option set.

## Expected Outputs

`tests/expected/mesh-contract.json` defines stable expectations:

- exact atom and bond counts;
- exact chain, residue, element, and bounding-box metadata;
- required OBJ/STL/PLY structural markers;
- minimum mesh primitive counts;
- forbidden invalid numeric strings.

## Runnable Checks

Run fixture validation and the public API smoke tests from the repository root:

```sh
node wasm-plugin/tests/validate-fixtures.mjs
typst compile --root . wasm-plugin/tests/api/public-api-smoke.typ /tmp/molfig-public-api-smoke.pdf
typst compile --root . wasm-plugin/tests/api/module-split-contract-smoke.typ /tmp/molfig-module-split-contract-smoke.pdf
typst compile --root package package/examples/1CRN.typ /tmp/molfig-example-1crn.pdf
typst compile --root package package/examples/1FYY.typ /tmp/molfig-example-1fyy.pdf
```

The richer structure-aware API smoke tests are also runnable:

```sh
node wasm-plugin/tests/fixtures/bcif/generate-bcif-fixture.mjs
node wasm-plugin/tests/validate-fixtures.mjs
typst compile --root . wasm-plugin/tests/api/future-structure-api-smoke.typ /tmp/molfig-future-structure-api-smoke.pdf
typst compile --root . wasm-plugin/tests/api/future-rich-api-smoke.typ /tmp/molfig-future-rich-api-smoke.pdf
typst compile --root package package/examples/9Z4O.typ /tmp/molfig-example-9z4o.pdf
```

Negative tests should fail compilation and match stderr:

```sh
typst compile --root . wasm-plugin/tests/api/negative-bad-input-format.typ /tmp/negative.pdf
typst compile --root . wasm-plugin/tests/api/negative-bad-output-format.typ /tmp/negative.pdf
```

## Review Focus

Reviewers should check that:

- PDB and CIF paths share one normalized structure model before meshing;
- OBJ, STL, and PLY exporters use the same mesh data;
- maquette integration receives generated mesh assets, not raw molecular text;
- cache keys include source path/content hash, representation, output format,
  quality, scale, centering, and color options;
- a Rust module split preserves the public Typst wrappers and the shared
  normalized structure model used by `info`, exporters, and `render-object`;
- assembly id, alternate-location selection, BinaryCIF/text CIF equivalence,
  cartoon/ribbon secondary-structure handling, and cartoon tuning options are
  observable in fixtures or smoke tests;
- error messages name the bad option and the accepted values.
