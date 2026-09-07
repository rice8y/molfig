// Compile from the repository root:
// typst compile --root . package/examples/theme.typ package/examples/theme.pdf
// Uses the local package and fails on theme/material contract regressions.
#import "theme-samples.typ" as samples
#let manifest = toml(read("../typst.toml", encoding: none))
#let package-import = "@preview/" + manifest.package.name + ":" + manifest.package.version

#set page(paper: "a4", margin: 16mm, numbering: "1")
#set text(font: "Helvetica Neue", size: 10pt)
#set heading(numbering: "1.")
#set par(leading: 0.65em)

= Theme gallery

Every panel pairs executable settings with the rendering produced by this
checkout. Swatches show the unlit OBJ material colors. Compilation also checks
theme metadata, expected colors, and the no-op controls.

The samples use published structural data: RCSB PDB 1HHO and 1CRN,
AlphaFold human insulin precursor AF-P01308-F1 (model v6), the official
SWISS-MODEL QMEAN example model_001, and Mol\*'s charge-annotated 7QPD.
The prediction scores and computed charges are genuine published annotations,
not experimental measurements or manually assigned test values.
See #raw("data/README.md") for source URLs, licenses, and the reproducible
PDB + QMEAN JSON to CIF conversion. Source coordinates and scores are preserved.

The QMEAN data and its rendered panel are adapted from SWISS-MODEL,
Computational Structural Biology Group, SIB / University of Basel, under
#link("https://creativecommons.org/licenses/by-sa/4.0/")[CC BY-SA 4.0].
QMEAN: Benkert et al. (2011),
#link("https://doi.org/10.1093/bioinformatics/btq662")[doi:10.1093/bioinformatics/btq662].
AlphaFold: Jumper et al. (2021),
#link("https://doi.org/10.1038/s41586-021-03819-2")[doi:10.1038/s41586-021-03819-2],
data under #link("https://creativecommons.org/licenses/by/4.0/")[CC BY 4.0].

Use the following setup once, then append any panel's render call.
For checkout-local debugging, replace the package import with
#raw("#import \"../lib.typ\" as molfig").

#text(size: 8pt, raw(samples.common-code(package-import), lang: "typ", block: true))

#pagebreak()
= All seven global themes

An explicit #raw("globalName") replaces every component theme. In particular,
#raw("globalName: \"chain-id\"") forces chain coloring; the default
#raw("color-theme: \"chain-id\"") instead preserves the representation preset.
Structural themes use ball-and-stick; per-residue annotation themes use cartoon.
Both are selected explicitly to avoid ViewerAuto annotation selection.

#for item in samples.global-samples {
  samples.panel(item)
  v(6pt)
}

#pagebreak()
= All three carbon subthemes

These calls set #raw("carbonColor") alone. Carbon changes; nitrogen and oxygen
must not. No #raw("globalName") or #raw("color-theme") is supplied.

#for item in samples.carbon-samples {
  samples.panel(item)
  v(6pt)
}

#pagebreak()
= Spacefill preset control

Carbon subthemes apply only within element-symbol coloring. The first panel
must match spacefill with an empty theme dictionary; the second explicitly
replaces its illustrative preset.

#for item in samples.preset-samples {
  samples.panel(item)
  v(6pt)
}

#pagebreak()
= Assembly versus crystal symmetry

All three panels must match. The public API currently cannot generate crystal
contacts or crystallographic symmetry expansions; these are negative controls
for biological assemblies, not demonstrations of a crystal-symmetry scene.

These panels reuse the bundled RCSB PDB 1CRN biological assembly 1.
The build checks identical geometry and material assignments, not only similar
images. No crystal-symmetry scene is implied by this control.

#for item in samples.symmetry-samples {
  samples.panel(item)
  v(6pt)
}
