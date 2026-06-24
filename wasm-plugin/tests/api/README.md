# Typst API Smoke Tests

Run from the repository root:

```sh
typst compile --root . wasm-plugin/tests/api/public-api-smoke.typ /tmp/molfig-public-api-smoke.pdf
typst compile --root . wasm-plugin/tests/api/module-split-contract-smoke.typ /tmp/molfig-module-split-contract-smoke.pdf
typst compile --root . wasm-plugin/tests/api/maquette-config-passthrough-smoke.typ /tmp/molfig-maquette-config-passthrough-smoke.pdf
typst compile --root . wasm-plugin/tests/api/future-structure-api-smoke.typ /tmp/molfig-future-structure-api-smoke.pdf
typst compile --root . wasm-plugin/tests/api/future-rich-api-smoke.typ /tmp/molfig-future-rich-api-smoke.pdf
typst compile --root . wasm-plugin/tests/api/9r1o-reference-smoke.typ /tmp/molfig-9r1o-reference-smoke.pdf
```

The smoke tests cover:

- bytes input from `read(..., encoding: none)`;
- Typst 0.15+ path input from `path(...)`;
- PDB, mmCIF, and BinaryCIF parsing;
- OBJ/STL/PLY export;
- equivalent normalized metadata across PDB and mmCIF fixtures;
- stable render-object shape for module-split work;
- assembly selection;
- alternate-location selection;
- Viewer `default`, `auto`, `cartoon`, and illustrative `spacefill` presets,
  plus `polymer-cartoon`, `ball-and-stick`, `ribbon`, and `backbone` representations;
- ViewerAuto pLDDT annotation dispatch from text CIF and BinaryCIF, plus an
  explicitly selected QMEAN color theme;
- Mol*-style `selected_visuals` and `realized_visuals` representation
  metadata;
- cartoon tuning options: `helix-profile`, `round-cap`, and
  `sheet-arrow-factor`;
- semantic render-object metadata for `dashed-tube` and `sheet` geometry;
- `render(...)` and `render-object(...)` maquette integration;
- maquette config passthrough for mesh metadata and render-object content.
- self-contained 9R1O PDB to OBJ generation and maquette rendering without a
  checked-in reference OBJ.

Negative smoke tests are intentionally expected to fail compilation:

```sh
typst compile --root . wasm-plugin/tests/api/negative-bad-input-format.typ /tmp/negative.pdf
typst compile --root . wasm-plugin/tests/api/negative-bad-output-format.typ /tmp/negative.pdf
```

Their stderr should mention the invalid option and accepted values.
