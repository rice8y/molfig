# Typst API Contract and Integration Tests

Run from the repository root:

```sh
typst compile --root . wasm-plugin/tests/api/public-api-contract.typ /tmp/molfig-public-api-contract.pdf
typst compile --root . wasm-plugin/tests/api/module-split-contract.typ /tmp/molfig-module-split-contract.pdf
typst compile --root . wasm-plugin/tests/api/native-renderer-options-contract.typ /tmp/molfig-native-renderer-options-contract.pdf
typst compile --root . wasm-plugin/tests/api/future-structure-api-contract.typ /tmp/molfig-future-structure-api-contract.pdf
typst compile --root . wasm-plugin/tests/api/future-rich-api-contract.typ /tmp/molfig-future-rich-api-contract.pdf
typst compile --root . wasm-plugin/tests/api/9r1o-reference-integration.typ /tmp/molfig-9r1o-reference-integration.pdf
node wasm-plugin/scripts/validate-illustrative-rendering.mjs
```

The contract and integration tests cover:

- bytes input from `read(..., encoding: none)`;
- Typst 0.15+ path input from `path(...)`;
- PDB, mmCIF, BinaryCIF, and XYZ parsing;
- OBJ/STL/PLY export;
- equivalent normalized metadata across PDB and mmCIF fixtures;
- stable native `render-result` shape for module-split work;
- assembly selection;
- alternate-location selection;
- Viewer `default`, `auto`, and `cartoon` presets, plus `spacefill`,
  `polymer-cartoon`, `ball-and-stick`, `ribbon`, and `backbone` representations;
- orthogonal `illustrative` renderer parameters and metadata, with invariant
  OBJ/MTL color-theme materials;
- ViewerAuto pLDDT annotation dispatch from text CIF and BinaryCIF, plus an
  explicitly selected QMEAN color theme;
- Mol*-style `selected_visuals` and `realized_visuals` representation
  metadata;
- cartoon tuning options: `helix-profile`, `round-cap`, and
  `sheet-arrow-factor`;
- semantic render metadata for `dashed-tube` and `sheet` geometry;
- native `render(...)` and `render-result(...)` output;
- SVG-by-default and explicit PNG render output, both backed by the same RGBA8
  pixels;
- strict native renderer option validation and resolved state metadata;
- self-contained 9R1O PDB export and native rendering without a checked-in
  reference OBJ.

Negative compile-failure tests are intentionally expected to fail compilation:

```sh
typst compile --root . wasm-plugin/tests/api/negative-bad-input-format.typ /tmp/negative.pdf
typst compile --root . wasm-plugin/tests/api/negative-bad-renderer-option.typ /tmp/negative.pdf
typst compile --root . wasm-plugin/tests/api/negative-bad-output-format.typ /tmp/negative.pdf
```

Their stderr should mention the invalid option and accepted values.
