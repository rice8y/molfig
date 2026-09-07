// Executable source of truth for theme.typ and docs/documentation.typ.
// Always load this checkout's WASM, not a previously installed Molfig package.
#import "../lib.typ" as molfig

#let view = (
  mesh-format: "obj",
  quality: "medium",
  center: true,
  output-format: "png",
  config: (azimuth: 35, elevation: 24, projection: "orthographic", width: 768, height: 648, background: ""),
  width: 64mm,
  height: 54mm,
)

#let chain-colors = ("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
#let element-colors = ("#999999", "#4259ff", "#ff2618")
#let noncarbon-colors = ("#4259ff", "#ff2618")

#let sample(id, title, theme, expected, note, representation: "ball-and-stick", assembly: "1", colors: (), file: "1HHO.cif", format: "cif") = (
  id: id, title: title, theme: theme, expected: expected, note: note,
  representation: representation, assembly: assembly, colors: colors, file: file, format: format,
)

#let global-samples = (
  sample("chain-id", "Chain ID", (globalName: "chain-id"), "chain-id",
    [1HHO oxyhemoglobin: alpha and beta chains have distinct colors, repeated in the second assembly copy.],
    colors: chain-colors.slice(0, 2)),
  sample("entity-id", "Entity ID", (globalName: "entity-id"), "entity-id",
    [1HHO: alpha/beta proteins, heme, oxygen, phosphate, and water are separate entities. Copies retain entity colors.],
    colors: chain-colors.slice(0, 2)),
  sample("operator-name", "Operator name", (globalName: "operator-name"), "operator-name",
    [1HHO assembly 1: the identity and twofold-related copies have different colors, regardless of chain or entity.],
    colors: chain-colors.slice(0, 2)),
  sample("element-symbol", "Element symbol", (globalName: "element-symbol", carbonColor: "element-symbol"), "element-symbol",
    [1HHO: carbon is grey, nitrogen blue, and oxygen red. Heme iron, sulfur, and phosphorus keep their element colors.],
    colors: element-colors),
  sample("plddt-confidence", "pLDDT confidence", (globalName: "plddt-confidence"), "plddt-confidence",
    [AlphaFold human insulin precursor, model v6: genuine predicted pLDDT spans three confidence bands; no residue exceeds 90.],
    file: "AF-P01308-F1-model_v6.cif", assembly: "asymmetric-unit", representation: "cartoon",
    colors: ("#ff7d45", "#ffdb13", "#65cbf3")),
  sample("qmean-score", "QMEAN score", (globalName: "qmean-score"), "qmean-score",
    [Official SWISS-MODEL QMEAN example model_001. Published local scores are mapped to their original residues, without rounding or rescaling.],
    file: "qmean-model_001.cif", assembly: "asymmetric-unit", representation: "cartoon", colors: ("#ff5000",)),
  sample("sb-ncbr-partial-charges", "SB-NCBR partial charges", (globalName: "sb-ncbr-partial-charges"), "sb-ncbr-partial-charges",
    [Mol\*'s 7QPD example: computed SQE+qp/Schindler 2021 atomic charges, summed per residue. Negative is red; positive is blue.],
    file: "7qpd.fw2.cif", assembly: "asymmetric-unit", representation: "cartoon"),
)

#let carbon-samples = (
  sample("carbon-chain", "Carbon: chain ID", (carbonColor: "chain-id"), "element-symbol",
    [Only carbon follows its chain. Nitrogen and oxygen keep their element colors. No global theme is set.],
    colors: chain-colors.slice(0, 2) + noncarbon-colors),
  sample("carbon-operator", "Carbon: operator name", (carbonColor: "operator-name"), "element-symbol",
    [Only carbon distinguishes the two assembly copies. This exercises carbonColor without globalName.],
    colors: chain-colors.slice(0, 2) + noncarbon-colors),
  sample("carbon-element", "Carbon: element symbol", (carbonColor: "element-symbol"), "element-symbol",
    [Carbon becomes grey using only carbonColor; nitrogen and oxygen are unchanged.],
    colors: element-colors),
)

#let preset-samples = (
  sample("spacefill-preset", "Spacefill: preset retained", (carbonColor: "element-symbol"), "illustrative",
    [A carbon subtheme alone does not replace the illustrative spacefill preset. Carbon stays lightened by entity.],
    representation: "spacefill"),
  sample("spacefill-element", "Spacefill: explicit element theme", (globalName: "element-symbol", carbonColor: "element-symbol"), "element-symbol",
    [Select element-symbol explicitly to obtain grey carbon in spacefill. Geometry and camera are unchanged.],
    representation: "spacefill", colors: element-colors),
)

#let symmetry-samples = (
  sample("symmetry-base", "Assembly: baseline", (:), "chain-id",
    [RCSB PDB 1CRN, assembly 1: cartoon geometry with its preset chain color.],
    representation: "cartoon", file: "1crn.bcif", format: "bcif"),
  sample("symmetry-operator", "Assembly: symmetryColor only", (symmetryColor: "operator-name"), "chain-id",
    [The rendering must equal the baseline. Biological assembly operators are not crystal-symmetry operators.],
    representation: "cartoon", file: "1crn.bcif", format: "bcif"),
  sample("symmetry-qmean", "Assembly: global + symmetry", (globalName: "chain-id", symmetryColor: "qmean-score"), "chain-id",
    [The rendering must still equal the baseline. An inapplicable symmetry override must not replace the chain theme with missing-score grey.],
    representation: "cartoon", file: "1crn.bcif", format: "bcif"),
)

#let code-value(value) = if type(value) == length { str(value / 1mm) + "mm" } else { repr(value) }

#let common-code(package-import) = "#import \"" + package-import + "\"\n\n// Source files and attribution: examples/data/README.md.\n#let view = (\n" + view.pairs().map(((key, value)) => {
  if type(value) == dictionary {
    "  " + key + ": (\n" + value.pairs().map(((k, v)) => "    " + k + ": " + code-value(v) + ",").join("\n") + "\n  ),"
  } else { "  " + key + ": " + code-value(value) + "," }
}).join("\n") + "\n)"

#let sample-code(item) = {
  let theme = if item.theme.len() == 0 { "  theme: (:)," } else {
    "  theme: (\n" + item.theme.pairs().map(((key, value)) => "    " + key + ": " + repr(value) + ",").join("\n") + "\n  ),"
  }
  "#molfig.render(\n  read(\"data/" + item.file + "\", encoding: none),\n  ..view,\n  format: " + repr(item.format) + ",\n  representation: " + repr(item.representation) + ",\n  assembly: " + repr(item.assembly) + ",\n" + theme + "\n)"
}

#let result(item) = {
  let data = read("data/" + item.file, encoding: none)
  let object = molfig.render-object(data, ..view, format: item.format,
    representation: item.representation, assembly: item.assembly, theme: item.theme)
  // Compile-time contracts accompany the visible results. A stale WASM or a
  // broken override fails the gallery build instead of silently printing it.
  assert(object.info.render_objects.len() > 0, message: item.id)
  assert(object.info.atom_count > 0, message: item.id)
  assert(object.info.structure.unit_count > 0, message: item.id)
  assert(object.info.render_objects.all(o => o.color_theme == item.expected), message: item.id)
  let colors = object.materials.values().dedup().sorted()
  for color in item.colors {
    assert(colors.contains(color), message: item.id + ": missing " + color)
  }
  if item.expected == "element-symbol" {
    let carbon = item.theme.at("carbonColor", default: "chain-id")
    assert(object.info.render_objects.all(o => o.carbon_color_theme == carbon), message: item.id)
  }
  if item.id == "plddt-confidence" {
    assert.eq(colors, item.colors.sorted(), message: item.id)
    assert.eq(object.info.atom_count, 839, message: item.id)
  }
  if item.id == "qmean-score" {
    assert(colors.len() > 10 and not colors.contains("#aaaaaa"), message: item.id)
    assert.eq(object.info.atom_count, 2306, message: item.id)
  }
  if item.id == "sb-ncbr-partial-charges" {
    assert(colors.len() > 10 and not colors.contains("#66ff00"), message: item.id)
  }
  if item.id.starts-with("symmetry-") {
    let base = molfig.render-object(data, ..view, format: item.format, representation: "cartoon", assembly: item.assembly)
    assert(object.info.render_objects.any(o => o.tag == "polymer"), message: item.id)
    assert.eq(object.mesh, base.mesh, message: item.id)
    assert.eq(object.materials, base.materials, message: item.id)
  }
  if item.id == "spacefill-preset" {
    let base = molfig.render-object(data, ..view, format: item.format, representation: "spacefill", assembly: item.assembly)
    assert.eq(object.mesh, base.mesh, message: item.id)
    assert.eq(object.materials, base.materials, message: item.id)
  }
  (
    image: object.content,
    colors: colors,
  )
}

#let panel(item) = {
  let output = result(item)
  // Annotation themes may have hundreds of materials. Show a bounded sample.
  let swatches = if output.colors.len() <= 12 { output.colors } else {
    range(12).map(i => output.colors.at(int(calc.round(i * (output.colors.len() - 1) / 11))))
  }
  block(breakable: false, width: 100%, inset: 7pt, stroke: 0.5pt + luma(82%), radius: 2pt)[
    #text(font: "Helvetica Neue", size: 10pt, weight: "bold", item.title)
    #v(4pt)
    #grid(columns: (1fr, 64mm), column-gutter: 3mm,
      [#set text(size: 8pt)
       #raw(sample-code(item), lang: "typ", block: true)
       #v(3pt)
       #item.note],
      [#output.image
       #align(center, grid(columns: swatches.map(_ => 3.5mm), column-gutter: 1mm,
         ..swatches.map(color => box(width: 3.5mm, height: 2mm, fill: rgb(color), stroke: 0.3pt + luma(65%)))))
       #align(center, text(size: 7pt, str(output.colors.len()) + if output.colors.len() == 1 { " material color" } else { " material colors" } + if output.colors.len() > 12 { " (12 shown)" } else { "" }))],
    )
  ]
}
