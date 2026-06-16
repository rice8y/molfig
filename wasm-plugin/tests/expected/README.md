# Expected Output Contract

This directory defines the review contract for generated OBJ/STL/PLY files.
It intentionally avoids committing large golden meshes because sphere and
cylinder tessellation can change while the public behavior remains correct.
Small `.contract` files under `obj/`, `stl/`, and `ply/` may pin exact
package-owned bytes or counts when a fixture needs a stricter regression
surface without committing a full mesh.

`ply/reference-fixtures.txt` is the manifest for package-owned PLY contract
diffs. Each listed `.contract` file names its input fixture and exporter
options so the Rust focused test can replay molfig export bytes and compare
counts, byte length, and the stable test hash.

Use `mesh-contract.json` as the source of truth:

- parse counts must match exactly;
- bounding boxes are in Angstrom and compare with `numericTolerance`;
- exported files must have deterministic names and no `nan`, `undefined`, or
  infinite numeric values;
- OBJ, STL, and PLY outputs must expose enough vertices/faces to prove that the
  package emitted real surface geometry, not only atom-center markers;
- PDB, CIF, and BinaryCIF versions of the same fixture must produce equivalent
  parse summaries and materially equivalent mesh statistics where the parser
  supports the format.
- future structure options such as assembly expansion, alternate-location
  selection, cartoon/ribbon representation, and `render-object` should be
  covered by `futureOptions` entries before implementation starts.

Recommended future test flow:

1. Convert each fixture to OBJ, STL, and PLY with default settings.
2. Parse the exported mesh headers and primitive counts.
3. Compare counts and parse summaries to `mesh-contract.json`.
4. Compile the Typst smoke files in `tests/api`.
5. Render one example document and inspect that maquette receives a mesh file
   path rather than raw molecular text.
