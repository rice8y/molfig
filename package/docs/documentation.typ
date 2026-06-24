#import "@preview/mantys:1.0.2": *

#let manifest = toml(read("../typst.toml", encoding: none))
#let package-id = manifest.package.name
#let package-version = manifest.package.version
#let product-name = "Molfig"
#let package-import = "@preview/" + package-id + ":" + package-version
#let rendered-9r1o-pdf = "../examples/9R1O.pdf"
#let rendered-representations-pdf = "../examples/representations.pdf"
#let rendered-render-object-pdf = "../examples/render-object.pdf"
#let rendered-theme-pdf = "../examples/theme.pdf"
#let rendered-exports-pdf = "../examples/exports.pdf"
#let rendered-metadata-pdf = "../examples/metadata.pdf"
#let rendered-draft-pdf = "../examples/9M1U.pdf"
#let example-9r1o-code = "#import \"" + package-import + "\"\n\n#set page(width: auto, height: auto, margin: 0mm)\n\n// Uses structural data from RCSB PDB / wwPDB.\n// PDB ID: 9R1O\n// PDB DOI: https://doi.org/10.2210/pdb9R1O/pdb\n// PDB archive data files are available under CC0 1.0.\n#let pdb = path(\"9R1O.pdb\")\n\n#molfig.render(\n  pdb,\n  format: \"pdb\",\n  representation: \"cartoon\",\n  assembly: \"1\",\n  mesh-format: \"obj\",\n  quality: \"high\",\n  center: true,\n  output-format: \"svg\",\n  config: (\n    azimuth: 35,\n    elevation: 24,\n    background: \"\",\n  ),\n)"

#let ic(value) = raw(str(value))

#let code-lines(parts) = text(
  size: 0.9em,
  grid(
    columns: 1,
    row-gutter: 0pt,
    ..parts.map(ic),
  ),
)

#let compact-ic(value) = text(size: 0.9em, ic(value))

#let color-code(value) = [
  #box(
    width: 0.72em,
    height: 0.72em,
    fill: rgb(value),
    stroke: 0.4pt + luma(60%),
    radius: 1pt,
  )#h(3pt)#ic(value)
]

#let element-color(symbol, value) = [#symbol #color-code(value)]

#let code(lang, body, title: none, file: none) = sourcecode(
  title: title,
  file: file,
  raw(body, block: true, lang: lang),
)

#let shell(body) = codesnippet(raw(body, block: true, lang: "bash"))

#let row2(a, b) = ([#a], [#b])
#let row3(a, b, c) = ([#a], [#b], [#c])
#let row4(a, b, c, d) = ([#a], [#b], [#c], [#d])

#let table2(head-a, head-b, rows) = table(
  columns: (28%, 72%),
  inset: 6pt,
  stroke: luma(82%),
  fill: (_, y) => if y == 0 { luma(94%) } else { none },
  [*#head-a*], [*#head-b*],
  ..rows.flatten(),
)

#let table3(head-a, head-b, head-c, rows) = table(
  columns: (26%, 22%, 52%),
  inset: 6pt,
  stroke: luma(82%),
  fill: (_, y) => if y == 0 { luma(94%) } else { none },
  [*#head-a*], [*#head-b*], [*#head-c*],
  ..rows.flatten(),
)

#let table4(head-a, head-b, head-c, head-d, rows) = table(
  columns: (20%, 16%, 24%, 40%),
  inset: 5pt,
  stroke: luma(82%),
  fill: (_, y) => if y == 0 { luma(94%) } else { none },
  [*#head-a*], [*#head-b*], [*#head-c*], [*#head-d*],
  ..rows.flatten(),
)

#let option-row(name, default, description) = row3(arg(name), default, description)

#let option-table(rows) = table3([Option], [Default], [Meaning], rows)

#let pipeline-row(stage, control, description) = row3(strong(stage), control, description)

#let pipeline-table(rows) = table3([Stage], [Control], [Role], rows)

#let example-result(
  pdf,
  caption,
  pdb-id,
  doi,
  width: 100%,
  source-note: none,
) = [
  #figure(
    block(
      width: width,
      inset: 4pt,
      stroke: luma(82%),
      radius: 2pt,
      image(pdf, width: 100%),
    ),
    caption: caption,
  )

  #info-alert[
    Structural data source: RCSB PDB / wwPDB, PDB ID #ic(pdb-id),
    #link(doi)[#doi]. PDB archive data files are available under CC0 1.0.
    #if source-note != none { source-note }
  ]
]

#show: mantys(
  ..manifest,
  title: [#product-name],
  subtitle: [Static molecular structure rendering for Typst],
  date: datetime.today(),
  abstract: [
    #product-name turns PDB, mmCIF, and BinaryCIF structures into static OBJ,
    STL, or PLY meshes inside a Rust WebAssembly plugin, then hands those meshes to
    #std.link("https://typst.app/universe/package/maquette", [maquette]) for
    document rendering.
  ],
  show-index: true,
  wrap-snippets: true,
  theme: create-theme(
    fonts: (
      serif: ("Times New Roman", "Georgia"),
      sans: ("Helvetica Neue", "Arial"),
      mono: ("Menlo", "Courier New"),
    ),
    text: (
      size: 11pt,
      font: ("Times New Roman", "Georgia"),
      fill: rgb(35, 31, 32),
    ),
    heading: (
      font: ("Helvetica Neue", "Arial"),
      fill: rgb(35, 31, 32),
    ),
    emph: (
      link: rgb("#1f4f73"),
    ),
    code: (
      size: 9pt,
      font: ("Menlo", "Courier New"),
      fill: rgb("#555555"),
    ),
  ),
)

= Overview <sec:overview>

#product-name is a Typst package for molecular figures in static documents. It
accepts PDB, mmCIF, and BinaryCIF input, converts the structure into a CPU-side
Mol\*-style Model/Structure/Unit representation, builds static geometry, exports
the mesh as OBJ, STL, or PLY, and delegates final rendering to maquette.

Mol\* is an open-source molecular visualization toolkit and web viewer. Molfig
uses Mol\* as the compatibility target for structure interpretation, static
geometry, and OBJ/STL export behavior; see
#link("https://molstar.org/")[molstar.org] and
#link("https://github.com/molstar/molstar")[molstar/molstar].

== Design Goals <sec:design-goals>

#product-name is designed around four constraints:

1. Typst packages cannot invent paths into the caller's project. On Typst 0.15.0
   or later, the caller may pass a project-side #ic("path(...)") value and
   #product-name will read it as bytes. For Typst 0.14-compatible documents, the
   caller should pass bytes produced by #ic("read(..., encoding: none)").
2. The plugin should run in Typst's WebAssembly environment without WebGL.
3. The molecular path should preserve Mol\*-style structure concepts rather
   than flattening everything at parse time.
4. The final image should be a document-native static rendering through
   maquette, not an interactive viewer snapshot.

#warning-alert[
  #product-name generates static meshes and is not a Mol\* WebGL renderer.
  #ic("representation: \"surface\"") computes the Mol\* molecular-surface
  field and marching-cubes mesh on the CPU. Gaussian volume and density-map
  visuals remain outside the static export contract; ViewerAuto's Huge and
  Gigantic Gaussian surfaces are likewise generated on the CPU.
]

== Rendering Pipeline <sec:pipeline>

#pipeline-table((
  pipeline-row("Read", ic("path(...) or read(..., encoding: none)"), [The user document supplies either a Typst 0.15+ project path or already-read structure bytes.]),
  pipeline-row("Parse", ic("format"), [The plugin parses PDB, text CIF/mmCIF, or BinaryCIF MessagePack data.]),
  pipeline-row("Model", [#ic("assembly"), #ic("alt-loc")], [Model, Structure, Unit, unit operators, altLoc policy, bonds, secondary structure, and coarse data are selected.]),
  pipeline-row("Geometry", [#ic("representation"), #ic("quality")], [The static geometry layer builds spheres, cylinders, tubes, sheets, ribbons, nucleotide visuals, carbohydrate visuals, and related mesh spans.]),
  pipeline-row("Export", ic("mesh-format"), [The mesh is serialized as OBJ, binary STL, or ASCII PLY bytes.]),
  pipeline-row("Render", ic("config"), [The exported bytes and maquette configuration are passed to maquette.]),
))

= Quickstart <sec:quickstart>

The following example renders PDB entry 9R1O. Put #ic("9R1O.typ") and
#ic("9R1O.pdb") in the same directory, then compile the Typst file.

#code("typ", example-9r1o-code, title: "Complete 9R1O example", file: "9R1O.typ")

#info-alert[
  The #ic("path(\"9R1O.pdb\")") input in this example requires Typst 0.15.0 or
  later. On an older Typst version, replace it with
  #ic("read(\"9R1O.pdb\", encoding: none)") and pass the resulting bytes to
  #ic("molfig.render").
]

#example-result(
  rendered-9r1o-pdf,
  [RCSB PDB entry 9R1O rendered by #product-name with the complete Quickstart settings.],
  "9R1O",
  "https://doi.org/10.2210/pdb9R1O/pdb",
  source-note: [Deposition authors: Petrenas, Ozga, Chubb, and Woolfson.],
)

== Reading Structure Files <sec:reading-files>

On Typst 0.15.0 or later, prefer passing a #ic("path(...)") value for
file-backed structure data. The path is resolved by Typst relative to the caller's
file, and #product-name reads it internally with #ic("encoding: none"). This keeps
package calls concise while preserving the caller-side project boundary.

#table3([Format], [Typst 0.15+ path], [Real structure used in this manual], (
  row3([PDB], [#ic("path(\"9R1O.pdb\")")], [RCSB PDB entry 9R1O.]),
  row3([mmCIF], [#ic("path(\"1FYY.cif\")")], [RCSB PDB entry 1FYY.]),
  row3([BinaryCIF], [#ic("path(\"1CRN.bcif\")")], [RCSB PDB entry 1CRN.]),
))

For documents that must also compile on Typst 0.14, read the file in the caller
document and pass the resulting bytes. Always use #ic("encoding: none") for
external structure files; it preserves BinaryCIF bytes and avoids Unicode
decoding loss for fixed-width PDB columns.

For the same files, use #ic("read(\"9R1O.pdb\", encoding: none)"),
#ic("read(\"1FYY.cif\", encoding: none)"), or
#ic("read(\"1CRN.bcif\", encoding: none)"). Complete, rendered examples for
all three archive formats appear below.

#info-alert[
  These files are RCSB PDB / wwPDB entries 9R1O, 1FYY, and 1CRN. Their PDB
  DOIs are #link("https://doi.org/10.2210/pdb9R1O/pdb")[9R1O],
  #link("https://doi.org/10.2210/pdb1FYY/pdb")[1FYY], and
  #link("https://doi.org/10.2210/pdb1CRN/pdb")[1CRN]. PDB archive data files
  are available under CC0 1.0.
]

Passing a Typst string is accepted for small inline examples. A string is treated
as inline molecular text, not as a file path.

== Choosing The First Options <sec:first-options>

#table3([Question], [Recommended setting], [Why], (
  row3([I want the Mol\* Viewer Cartoon style.], [#ic("representation: \"cartoon\"")], [Uses the Viewer Quick Styles Cartoon preset: polymer cartoon plus atomic detail for ligands and other non-polymer components.]),
  row3([I want only the polymer cartoon.], [#ic("representation: \"polymer-cartoon\"")], [Uses the Mol\* Cartoon provider defaults: polymer trace, nucleotide ring, and polymer gap visuals.]),
  row3([I need readable geometry diffs.], [#ic("mesh-format: \"obj\"")], [OBJ is text and preserves face groups.]),
  row3([I need compact triangle bytes.], [#ic("mesh-format: \"stl\"")], [STL is binary and simple for downstream mesh tools.]),
  row3([I want text mesh plus group metadata.], [#ic("mesh-format: \"ply\"")], [PLY is ASCII and carries package-owned group comments/properties.]),
  row3([Compilation is slow.], [#ic("quality: \"auto\"")], [Lets Molfig pick lower tessellation for large structures.]),
))

= Public API <sec:api>

Every public command takes #arg[data] as its first positional argument. It can be
bytes from #ic("read(..., encoding: none)"), inline molecular text as a string,
or a Typst 0.15+ #ic("path(...)") value. The same mesh options are shared by
render, render-object, export, and metadata commands unless noted otherwise.

#table3([Command], [Returns], [Use when], (
  row3([@cmd:render[-]], [content], [You want the figure directly in the document.]),
  row3([@cmd:render-object[-]], [dictionary], [You want the rendered content plus raw mesh bytes and metadata.]),
  row3([@cmd:to-obj[-]], [bytes], [You need OBJ mesh bytes for maquette or an external tool.]),
  row3([@cmd:to-mtl[-]], [bytes], [You need the companion material library for OBJ output.]),
  row3([@cmd:to-stl[-]], [bytes], [You need binary STL triangle data.]),
  row3([@cmd:to-ply[-]], [bytes], [You need ASCII PLY output with Molfig metadata comments/properties.]),
  row3([@cmd:info[-]], [dictionary], [You want structure and render planning metadata without rendering.]),
  row3([@cmd:mesh-info[-]], [dictionary], [You want maquette's mesh metadata for the generated OBJ/STL/PLY bytes.]),
))

== Rendering Commands <sec:rendering-commands>

#command(
  "render",
  arg("data"),
  arg("format", _value: values.value.with("auto")),
  arg("mesh-format", _value: values.value.with("obj")),
  arg("representation", _value: values.value.with("default")),
  arg("config", _value: values.value.with((:))),
  arg("width", _value: values.value.with("auto")),
  arg("height", _value: values.value.with("auto")),
  arg("output-format", _value: values.value.with("png")),
)[
  Converts molecular bytes to a static mesh and delegates rendering to
  maquette.

  #argument("data", types: (bytes, str), command: "render")[
    Molecular structure input. Pass #ic("path(\"file.pdb\")") on Typst 0.15+,
    or #ic("read(\"file.pdb\", encoding: none)") for Typst 0.14-compatible
    documents.
  ]

  #argument("format", choices: ("auto", "pdb", "cif", "mmcif", "bcif", "binarycif", "binary-cif"), default: "auto", command: "render")[
    Input parser selection. Use an explicit format for reproducible documents.
  ]

  #argument("mesh-format", choices: ("obj", "stl", "ply"), default: "obj", command: "render")[
    Intermediate mesh format sent to maquette.
  ]

  #argument("config", types: ("dictionary",), default: (:), command: "render")[
    JSON-encoded and passed to maquette. For OBJ rendering, Molfig adds the
    active color theme as a maquette material map; explicit
    #ic("config.materials") entries override generated colors with the same
    material id.
  ]

  #argument("width", default: "auto", command: "render")[
    Width forwarded to the maquette rendering command.
  ]

  #argument("height", default: "auto", command: "render")[
    Height forwarded to the maquette rendering command.
  ]

  #argument("output-format", choices: ("png", "svg"), default: "png", command: "render")[
    Output format forwarded to maquette. PNG uses maquette's Z-buffer raster
    output and is recommended for high-poly meshes, large assemblies, and
    spacefill representations. SVG preserves vector geometry, but is intended
    for small to moderately sized meshes because each visible mesh face
    contributes SVG content that Typst must parse.
  ]
]

#command("render-object", arg("data"), arg("mesh-format", _value: values.value.with("obj")))[
  Returns a dictionary with:

  #table2([Key], [Value], (
    row2([#ic("kind")], [The string #ic("\"render-object\"").]),
    row2([#ic("format")], [The selected mesh format.]),
    row2([#ic("mesh_format")], [The selected mesh format with snake_case naming.]),
    row2([#ic("mesh")], [Generated OBJ/STL/PLY bytes.]),
    row2([#ic("materials")], [Generated maquette material map for OBJ, or an empty dictionary for STL/PLY.]),
    row2([#ic("info")], [The same dictionary returned by @cmd:info[-] for the same mesh options.]),
    row2([#ic("content")], [The Typst content returned by @cmd:render[-].]),
  ))
]

#code("typ", "#import \"" + package-import + "\"\n\n// Structural data: RCSB PDB / wwPDB entry 1FYY.\n// https://doi.org/10.2210/pdb1FYY/pdb (CC0 1.0)\n#let object = molfig.render-object(\n  read(\"1FYY.cif\", encoding: none),\n  format: \"cif\",\n  representation: \"cartoon\",\n  mesh-format: \"obj\",\n  assembly: \"1\",\n  config: (azimuth: 30, elevation: 18, background: \"\"),\n  width: 92mm,\n  height: 64mm,\n  output-format: \"svg\",\n)\n\n#object.content\n#align(center)[\n  Atoms: #object.info.atom_count \\\n  Mesh format: #object.mesh_format\n]", title: "Render and inspect RCSB PDB entry 1FYY", file: "render-object.typ")

#example-result(
  rendered-render-object-pdf,
  [The complete #ic("render-object") example for 1FYY, including values read from its returned metadata.],
  "1FYY",
  "https://doi.org/10.2210/pdb1FYY/pdb",
  width: 74%,
  source-note: [Primary citation: Volk et al., #emph[Biochemistry] 39,
    14040--14053 (2000),
    #link("https://doi.org/10.1021/bi001669l")[doi:10.1021/bi001669l].],
)

== Export Commands <sec:export-commands>

#command("to-obj", arg("data"))[
  Returns OBJ bytes. OBJ is the most useful interchange format when a downstream
  tool wants text geometry, material references, and group labels.
]

#command("to-mtl", arg("data"))[
  Returns the companion MTL material library for OBJ output. Use it alongside
  @cmd:to-obj[-] when exporting outside Typst.
]

#command("to-stl", arg("data"))[
  Returns binary STL bytes. STL follows Mol\* static exporter behavior: the
  header starts with #ic("\"Exported from Mol*\"") and the two-byte facet
  attribute field is kept at zero.
]

#command("to-ply", arg("data"))[
  Returns ASCII PLY bytes. PLY is a Molfig-owned static mesh format because the
  pinned Mol\* geo-export extension does not provide PLY output.
]

#code("typ", "#import \"" + package-import + "\"\n\n// Structural data: RCSB PDB / wwPDB entry 1CRN.\n// https://doi.org/10.2210/pdb1CRN/pdb (CC0 1.0)\n#let data = read(\"1CRN.bcif\", encoding: none)\n#let obj = molfig.to-obj(data, format: \"bcif\", representation: \"cartoon\")\n#let mtl = molfig.to-mtl(data, format: \"bcif\", representation: \"cartoon\")\n#let stl = molfig.to-stl(data, format: \"bcif\", representation: \"cartoon\")\n#let ply = molfig.to-ply(data, format: \"bcif\", representation: \"cartoon\")\n\n#table(\n  columns: (1fr, 1fr),\n  table.header([*Export*], [*Generated bytes*]),\n  [OBJ], [#obj.len()],\n  [MTL], [#mtl.len()],\n  [STL], [#stl.len()],\n  [PLY], [#ply.len()],\n)", title: "Export RCSB PDB entry 1CRN to every mesh format", file: "exports.typ")

#example-result(
  rendered-exports-pdf,
  [Byte lengths produced by the complete 1CRN export example.],
  "1CRN",
  "https://doi.org/10.2210/pdb1CRN/pdb",
  width: 64%,
  source-note: [Primary citation: Teeter, #emph[Proceedings of the National
    Academy of Sciences] 81, 6014--6018 (1984),
    #link("https://doi.org/10.1073/pnas.81.19.6014")[doi:10.1073/pnas.81.19.6014].],
)

== Metadata Commands <sec:metadata-commands>

#command("info", arg("data"))[
  Parses the structure and returns molecular and mesh-planning metadata without
  maquette rendering. This is the fastest public command for debugging parser,
  assembly, altLoc, bond, secondary structure, and representation behavior.
]

#command("mesh-info", arg("data"), arg("mesh-format", _value: values.value.with("obj")))[
  Generates a mesh, then delegates mesh metadata extraction to maquette's
  #ic("get-obj-info"), #ic("get-stl-info"), or #ic("get-ply-info") helper.
]

= Shared Mesh Options <sec:shared-options>

The following options are accepted by @cmd:render[-], @cmd:render-object[-],
@cmd:to-obj[-], @cmd:to-mtl[-], @cmd:to-stl[-], @cmd:to-ply[-], and
@cmd:info[-].

== Input Selection <sec:input-selection>

#option-table((
  option-row("format", ic("\"auto\""), [Input parser. Supported values: #ic("\"auto\""), #ic("\"pdb\""), #ic("\"cif\""), #ic("\"mmcif\""), #ic("\"bcif\""), #ic("\"binarycif\""), #ic("\"binary-cif\"").]),
  option-row("block-index", ic("none"), [Zero-based BinaryCIF data block index. Useful when one BinaryCIF file contains multiple data blocks.]),
  option-row("block-header", ic("\"\""), [Exact BinaryCIF block header. When set, it takes precedence over #arg[block-index].]),
))

#table3([Format], [Typical extension], [Notes], (
  row3([#ic("\"pdb\"")], [#ic(".pdb")], [Fixed-width text. Keep bytes to preserve columns.]),
  row3([#ic("\"cif\"") or #ic("\"mmcif\"")], [#ic(".cif"), #ic(".mmcif")], [Text CIF categories and loops.]),
  row3([#ic("\"bcif\"")], [#ic(".bcif")], [BinaryCIF MessagePack data with encoded columns and optional masks.]),
  row3([#ic("\"auto\"")], [any], [Convenient for exploration, but explicit formats are better for release documents.]),
))

BinaryCIF numeric atom-site columns are decoded into typed accessors before
falling back to string cells. Molfig does not round-trip BinaryCIF through
synthetic CIF text.

== Structure Selection <sec:structure-options>

#option-table((
  option-row("assembly", ic("\"1\""), [Biological assembly id. Use #ic("\"asymmetric-unit\""), #ic("\"none\""), or an empty string for source coordinates without assembly operators.]),
  option-row("alt-loc", ic("\"\""), [Alternate location id or policy. Supports concrete ids such as #ic("\"A\""), plus #ic("\"highest-occupancy\"") and #ic("\"all\"").]),
  option-row("infer-bonds", ic("true"), [Infer covalent bonds when explicit bond data are absent.]),
  option-row("center", ic("true"), [Recenters exported vertices around Mol\*-style visible export bounds.]),
))

Assembly expansion keeps the source model and unit operators as the semantic
structure boundary, then materializes mesh coordinates for export. This means
metadata can still report selected assembly operators even though OBJ/STL/PLY
are static mesh formats.

== Representation Selection <sec:representation-options>

#option-table((
  option-row("representation", ic("\"default\""), [Uses the Mol\* Viewer configured default preset, falling back to its automatic size-dependent selection. Use #ic("\"cartoon\"") to pin Viewer Quick Styles Cartoon or #ic("\"polymer-cartoon\"") for the standalone Cartoon provider.]),
  option-row("color-theme", ic("\"chain-id\""), [Selects #ic("\"chain-id\""), #ic("\"element-symbol\""), #ic("\"entity-id\""), #ic("\"operator-name\""), #ic("\"plddt-confidence\""), #ic("\"qmean-score\""), or #ic("\"sb-ncbr-partial-charges\""). Color is serialized through OBJ material references and the companion MTL output. STL and the current PLY output do not carry these colors, so themed document rendering requires #ic("mesh-format: \"obj\"").]),
  option-row("theme", ic("(:)"), [Applies Mol\* Viewer preset overrides. Supported keys are #ic("globalName"), #ic("carbonColor"), and #ic("symmetryColor"). An empty dictionary preserves the preset's normal component-specific themes.]),
))

The representation names follow Mol\* Viewer Quick Styles rather than only the
low-level representation registry. In the Viewer, Cartoon applies the
#ic("polymer-and-ligand") preset; it is not limited to the standalone Cartoon
provider. Spacefill is available as #ic("\"spacefill\""). Surface corresponds
to Mol\*'s #ic("molecular-surface") preset and is selected with
#ic("representation: \"surface\""). It uses the #ic("all") component,
#ic("entity-id") colors with the water override, a physical size theme,
#ic("probeRadius: 1.4"), #ic("probePositions: 36"), and the Mol\* CPU
marching-cubes path.

#info-alert[
  #ic("\"default\"") and #ic("\"auto\"") preserve their distinct public
  names and follow Mol\* size routing for Small, Medium, and Large structures.
  #ic("\"default\"") first applies the Viewer annotation priority described
  below, while #ic("\"auto\"") requests size routing directly. Huge and
  Gigantic structures use Mol\*'s CPU Gaussian-surface branches, including
  trace-only and structure-level routing where prescribed by ViewerAuto.
]

=== Cartoon, Spacefill, And Surface <sec:representation-comparison>

The same structure and camera make the geometry difference explicit. Cartoon
emphasizes polymer topology and secondary structure; Spacefill exposes atomic
packing with van der Waals spheres; Surface shows the continuous
solvent-excluded molecular envelope computed by the CPU molecular-surface path.

#code("typ", "#import \"" + package-import + "\"\n\n// Structural data: RCSB PDB / wwPDB entry 1CRN.\n// https://doi.org/10.2210/pdb1CRN/pdb (CC0 1.0)\n#let data = read(\"1CRN.bcif\", encoding: none)\n#let view(rep) = molfig.render(\n  data, format: \"bcif\", representation: rep,\n  mesh-format: \"obj\", quality: \"high\", center: true,\n  output-format: \"svg\",\n  config: (azimuth: 35, elevation: 24, background: \"\"),\n  width: 55mm, height: 50mm,\n)\n#let panel(label, rep) = [\n  #align(center, strong(label))\n  #block(\n    width: 55mm, height: 50mm, clip: true,\n    align(center + horizon, scale(\n      x: 165%, y: 165%, origin: center + horizon, view(rep),\n    )),\n  )\n]\n#grid(\n  columns: (1fr, 1fr, 1fr), column-gutter: 2mm,\n  panel([Cartoon], \"cartoon\"),\n  panel([Spacefill], \"spacefill\"),\n  panel([Surface], \"surface\"),\n)", title: "Compare three representations of RCSB PDB entry 1CRN", file: "representations.typ")

#example-result(
  rendered-representations-pdf,
  [RCSB PDB entry 1CRN rendered with identical camera and quality settings as Cartoon, Spacefill, and Surface.],
  "1CRN",
  "https://doi.org/10.2210/pdb1CRN/pdb",
  source-note: [Primary citation: Teeter, #emph[Proceedings of the National
    Academy of Sciences] 81, 6014--6018 (1984),
    #link("https://doi.org/10.1073/pnas.81.19.6014")[doi:10.1073/pnas.81.19.6014].],
)

=== Viewer Theme Overrides <sec:viewer-theme-overrides>

Mol\* Viewer presets accept a #ic("theme") dictionary through their common
representation parameters. #ic("globalName") replaces the provider theme for each component,
#ic("carbonColor") controls carbon atoms in ball-and-stick ligand,
non-standard, and branched components, and #ic("symmetryColor") replaces the
polymer theme only for non-assembly crystallographic symmetry units. Water,
ion, and lipid components keep element-symbol carbon coloring as in the Viewer
preset.

#code("typ", "#import \"" + package-import + "\"\n\n// Structural data: RCSB PDB / wwPDB entry 1CRN.\n// https://doi.org/10.2210/pdb1CRN/pdb (CC0 1.0)\n#molfig.render(\n  read(\"1CRN.bcif\", encoding: none),\n  format: \"bcif\",\n  representation: \"cartoon\",\n  theme: (\n    globalName: \"element-symbol\",\n    carbonColor: \"chain-id\",\n    symmetryColor: \"operator-name\",\n  ),\n  mesh-format: \"obj\",\n  quality: \"high\",\n  center: true,\n  output-format: \"svg\",\n  config: (azimuth: 35, elevation: 24, background: \"\"),\n  width: 92mm,\n  height: 68mm,\n)", title: "Apply Mol* Viewer theme overrides to RCSB PDB entry 1CRN", file: "theme.typ")

#example-result(
  rendered-theme-pdf,
  [RCSB PDB entry 1CRN with the Viewer theme override example applied.],
  "1CRN",
  "https://doi.org/10.2210/pdb1CRN/pdb",
  width: 66%,
  source-note: [Primary citation: Teeter (1984),
    #link("https://doi.org/10.1073/pnas.81.19.6014")[doi:10.1073/pnas.81.19.6014].],
)

=== Annotation Color Rules <sec:annotation-colors>

With #ic("representation: \"default\""), Molfig follows the ViewerAuto preset
and selects the first applicable embedded annotation in this order:
#ic("plddt-confidence"), #ic("qmean-score"), then
#ic("sb-ncbr-partial-charges"). These names may also be selected explicitly
with #arg[color-theme]. #ic("representation: \"auto\"") skips annotation
selection and performs only structure-size routing.

#table3([Theme], [Value rule], [Colors], (
  row3([#ic("plddt-confidence")], [#ic("score <= 50"), #ic("<= 70"), #ic("<= 90"), then #ic("> 90")], [#color-code("#ff7d45"), #color-code("#ffdb13"), #color-code("#65cbf3"), #color-code("#0053d6")]),
  row3([#ic("qmean-score")], [Values through #ic("0.5") use orange; #ic("0.5..1.0") interpolates to blue.], [#color-code("#ff5000") to #color-code("#025afd")]),
  row3([#ic("sb-ncbr-partial-charges")], [Residue charge: negative through zero to positive.], [#color-code("#ff0000") through #color-code("#ffffff") to #color-code("#0000ff")]),
))

Missing pLDDT values fall back to #ic("B_iso_or_equiv"). Unavailable or
negative quality scores use #color-code("#aaaaaa"). Missing partial-charge
locations use #color-code("#66ff00"), matching the Mol\* provider. Annotation
colors require OBJ when rendered through maquette because STL and the current
PLY schema have no material color channel.

=== Chain ID Color Rules <sec:chain-id-colors>

The #ic("\"chain-id\"") theme follows these rules:

The Viewer #ic("\"spacefill\"") preset intentionally overrides this option
with Mol\* illustrative entity-id coloring. It lightens carbon colors in CIE
Lab and forces water entities to #color-code("#ff0d0d").

- For atomic models, the author chain id (#ic("auth_asym_id")) is used when it
  is present; otherwise the label chain id (#ic("label_asym_id")) is used. PDB
  input normally has the same value for both. Coarse IHM spheres and Gaussians
  use their asym id.
- Unique chain ids are collected in model order. Atomic chains are collected
  first, followed by previously unseen coarse-sphere and coarse-Gaussian chain
  ids. Assembly copies retain the source chain id and therefore retain its
  color.
- Colors are assigned from the following 25-color Mol\* many-distinct palette.
  If a structure contains more than 25 chain ids, the palette repeats from the
  first color.

#table2([Palette positions], [Colors], (
  row2([1--5], [#color-code("#1b9e77"), #color-code("#d95f02"), #color-code("#7570b3"), #color-code("#e7298a"), #color-code("#66a61e")]),
  row2([6--10], [#color-code("#e6ab02"), #color-code("#a6761d"), #color-code("#666666"), #color-code("#e41a1c"), #color-code("#377eb8")]),
  row2([11--15], [#color-code("#4daf4a"), #color-code("#984ea3"), #color-code("#ff7f00"), #color-code("#ffff33"), #color-code("#a65628")]),
  row2([16--20], [#color-code("#f781bf"), #color-code("#999999"), #color-code("#66c2a5"), #color-code("#fc8d62"), #color-code("#8da0cb")]),
  row2([21--25], [#color-code("#e78ac3"), #color-code("#a6d854"), #color-code("#ffd92f"), #color-code("#e5c494"), #color-code("#b3b3b3")]),
))

Chain-associated geometry, including cartoon, ribbon, backbone, nucleotide,
carbohydrate, gap, and coarse visuals, uses the assigned chain color unless the
visual has an explicit atomic material. Atomic visuals apply a
chain-and-element rule:

- Carbon uses its chain color for ordinary atoms, polymers, ligands, and
  branched components.
- Water, ion, and lipid components use element-symbol colors for every atom,
  including carbon.
- Non-carbon atoms use their element-symbol color. The symbol is read from
  #ic("type_symbol") when available, otherwise from the parsed element field.
  Unknown symbols fall back to white.
- Bond geometry inherits the material of its source atom. Component visuals
  that emit two directed half-bonds consequently color each half from its
  adjacent atom.

Element-symbol colors include Mol\*'s default lightness adjustment. The final
RGB values written to OBJ materials are:

#table2([Elements], [Final colors], (
  row2([Hydrogen isotopes], [#element-color("H", "#ffffff"), #element-color("D", "#ffffca"), #element-color("T", "#ffffaa")]),
  row2([Common organic elements], [#element-color("C", "#999999"), #element-color("N", "#4259ff"), #element-color("O", "#ff2618"), #element-color("P", "#ff8a14"), #element-color("S", "#ffff3e"), #element-color("Se", "#ffab17")]),
  row2([Halogens], [#element-color("F", "#9aea5a"), #element-color("Cl", "#37fb2e"), #element-color("Br", "#b13431"), #element-color("I", "#9e179e")]),
  row2([Alkali and alkaline-earth metals], [#element-color("Na", "#b566fd"), #element-color("Mg", "#95ff1f"), #element-color("K", "#994ade"), #element-color("Ca", "#4eff1e")]),
  row2([Transition metals], [#element-color("Mn", "#a683d1"), #element-color("Fe", "#eb703c"), #element-color("Co", "#fb9aaa"), #element-color("Ni", "#5cda5a"), #element-color("Cu", "#d3893c"), #element-color("Zn", "#8689ba")]),
  row2([Other or unknown], [#color-code("#ffffff")]),
))

The companion MTL records opacity #ic("0.3") for branched components,
#ic("0.6") for water and lipids, and #ic("1.0") otherwise. Molfig's automatic
maquette material map currently forwards RGB colors only; it does not forward
these MTL opacity values.

#table4([Value], [Best for], [Main geometry], [Notes], (
  row4([#ic("\"default\"")], [Viewer-compatible automatic figures], [ViewerAuto annotation theme, then size-selected geometry], [Uses pLDDT, QMEAN, or SB-NCBR partial charges when applicable, with Gaussian surfaces for Huge/Gigantic structures.]),
  row4([#ic("\"auto\"")], [Explicit automatic selection], [Atomic detail, polymer Cartoon, or a Gaussian surface according to structure size], [Huge/Gigantic routing follows the ViewerAuto size thresholds.]),
  row4([#ic("\"cartoon\"")], [General figures], [Polymer Cartoon plus atomic detail for non-polymers and carbohydrate visuals], [Pins the Mol\* Viewer Quick Styles Cartoon preset.]),
  row4([#ic("\"polymer-cartoon\"")], [Polymer-only figures], [Polymer trace, nucleotide ring, and polymer gap visuals], [Standalone Mol\* Cartoon provider defaults.]),
  row4([#ic("\"spacefill\"")], [Atomic packing], [Atom spheres], [Viewer Quick Styles Spacefill with illustrative entity-id colors and water override.]),
  row4([#ic("\"surface\"")], [Solvent-accessible shape], [CPU molecular-surface field and marching-cubes mesh], [Viewer Quick Styles Surface with entity-id colors and red water override.]),
  row4([#ic("\"ball-and-stick\"")], [Ligands and small molecules], [Atom spheres plus bond cylinders], [Good for local chemistry and explicit bond inspection.]),
  row4([#ic("\"ribbon\"")], [Backbone shape], [Ribbon-oriented polymer geometry], [Good for broad fold visibility.]),
  row4([#ic("\"backbone\"")], [Trace-only views], [Backbone cylinders and spheres], [Lower-detail alternative for polymer paths.]),
))

== Geometry And Quality <sec:geometry-options>

#option-table((
  option-row("quality", ic("\"custom\""), [#ic("\"custom\"") preserves explicit numeric controls. #ic("\"auto\"") chooses a preset from structure size. Also accepts #ic("\"highest\""), #ic("\"higher\""), #ic("\"high\""), #ic("\"medium\""), #ic("\"low\""), #ic("\"lower\""), #ic("\"lowest\"").]),
  option-row("sphere-detail", ic("2"), [Sphere tessellation detail. Public input is clamped by the plugin.]),
  option-row("linear-segments", ic("8"), [Curve subdivisions for polymer and tube paths.]),
  option-row("radial-segments", ic("16"), [Radial subdivisions for cylinders, tubes, and ribbon-like profiles.]),
  option-row("radius-scale", ic("1.0"), [Global radius multiplier.]),
  option-row("atom-radius", ic("0.28"), [Explicit atom radius for simple representations.]),
  option-row("bond-radius", ic("0.12"), [Explicit bond radius for cylinder-based bond visuals.]),
  option-row("ribbon-radius", ic("0.2"), [Tube/ribbon radius control.]),
  option-row("ribbon-width", ic("0.55"), [Ribbon width control.]),
  option-row("helix-profile", ic("\"elliptical\""), [#ic("\"elliptical\""), #ic("\"rounded\""), or #ic("\"square\""). #ic("\"sheet\"") is accepted as an alias for square.]),
  option-row("round-cap", ic("false"), [Enable rounded caps where supported by the selected cartoon/ribbon path.]),
  option-row("sheet-arrow-factor", ic("1.5"), [Scales sheet arrow geometry.]),
  option-row("tubular-helices", ic("false"), [Prefer tubular helix geometry in supported cartoon paths.]),
))

Quality presets override the explicit numeric tessellation controls. Use
#ic("quality: \"custom\"") when exact reproducibility matters, and use
#ic("quality: \"auto\"") or a lower preset when large structures compile too
slowly.

#table4([Preset], [Sphere detail], [Radial segments], [Linear segments], (
  row4([#ic("\"highest\"")], [3], [36], [18]),
  row4([#ic("\"higher\"")], [3], [28], [14]),
  row4([#ic("\"high\"")], [2], [20], [10]),
  row4([#ic("\"medium\"")], [1], [12], [8]),
  row4([#ic("\"low\"")], [0], [8], [3]),
  row4([#ic("\"lower\"")], [0], [4], [2]),
  row4([#ic("\"lowest\"")], [0], [2], [1]),
))

#info-alert[
  The public default is #ic("quality: \"custom\""), so the numeric defaults in
  this manual are honored unless you explicitly choose a preset or #ic("\"auto\"").
]

= Maquette Rendering <sec:maquette>

@cmd:render[-] calls one of maquette's mesh renderers based on
#arg[mesh-format]: #ic("render-obj"), #ic("render-stl"), or #ic("render-ply").
The #arg[config] dictionary is JSON-encoded and passed through. For OBJ,
Molfig derives maquette's #ic("materials") dictionary from the exported
material ids so #ic("color-theme: \"chain-id\"") is visible in the document.
User-supplied material entries are merged last and therefore take precedence.
STL has no material channel, and Molfig's current PLY schema contains no vertex
or face color properties. Consequently, maquette renders STL and PLY without
the selected color theme.

The complete 1FYY #ic("render-object") example above demonstrates camera,
background, dimensions, and raster-output passthrough, and its rendered result
shows the exact content returned in #ic("object.content").

Molfig does not validate maquette-specific keys. Unknown keys are passed to
maquette, so maquette remains the source of truth for camera and renderer
settings.

== Choosing A Mesh Format <sec:choosing-format>

#table4([Format], [Command], [Strength], [Tradeoff], (
  row4([OBJ], [@cmd:to-obj[-]], [Readable text, group labels, material references], [Larger than binary STL and needs MTL when material rows are consumed externally.]),
  row4([MTL], [@cmd:to-mtl[-]], [Companion material rows for OBJ], [Not a geometry format by itself.]),
  row4([STL], [@cmd:to-stl[-]], [Compact binary triangle stream], [No stable text metadata channel; facet attribute bytes are zero.]),
  row4([PLY], [@cmd:to-ply[-]], [Text mesh with Molfig group metadata], [Package-owned output, not a pinned Mol\* exporter format.]),
))

For visual documents, start with OBJ. Switch to STL only when the downstream
consumer specifically wants binary STL, and use PLY when the group metadata is
more useful than OBJ/MTL compatibility.

= Metadata Reference <sec:metadata>

@cmd:info[-] returns a dictionary. It is intended for document logic,
regression tests, and debugging. The exact set of keys may grow as Molfig's
Mol\*-parity layer grows, but the current major groups are:

#table3([Key], [Type], [Meaning], (
  row3([#ic("atom_count")], [integer], [Visible atom count after selected assembly/altLoc handling and geometry expansion.]),
  row3([#ic("bond_count")], [integer], [Visible bond count used by static geometry.]),
  row3([#ic("alt_locs")], [array], [Available alternate location ids from the source coordinates.]),
  row3([#ic("alt_locs_info")], [dictionary], [Selected altLoc policy and available ids.]),
  row3([#ic("assemblies")], [array], [Available biological assembly summaries.]),
  row3([#ic("assembly")], [dictionary], [Selected assembly id.]),
  row3([#ic("source_data")], [dictionary], [Input source categories, row counts, and original source kind.]),
  row3([#ic("structure")], [dictionary], [Model/Structure/Unit counts, boundary, conformation, segments, ranges, and lookup3d summary.]),
  row3([#ic("secondary_structure")], [dictionary], [Helix and sheet ranges.]),
  row3([#ic("representation")], [dictionary], [Selected and realized visual names.]),
  row3([#ic("render_objects")], [array], [Semantic render-object spans with geometry type, visual, chain, residue range, group id, component, representation tag/order, provider color-theme metadata, and value-cell style counts.]),
  row3([#ic("bond_metadata")], [dictionary], [Counts by bond source and flags such as struct_conn, index_pair, chem_comp, aromatic, and resonance.]),
  row3([#ic("bounds")], [dictionary], [Expanded coordinate bounds before export centering.]),
))

#code("typ", "#import \"" + package-import + "\"\n\n// Structural data: RCSB PDB / wwPDB entry 9R1O.\n// https://doi.org/10.2210/pdb9R1O/pdb (CC0 1.0)\n#let meta = molfig.info(\n  read(\"9R1O.pdb\", encoding: none),\n  format: \"pdb\",\n  representation: \"cartoon\",\n  assembly: \"1\",\n  alt-loc: \"highest-occupancy\",\n)\n\n#table(\n  columns: (1fr, 1fr),\n  table.header([*Property*], [*Value*]),\n  [Atoms], [#meta.atom_count],\n  [Bonds], [#meta.bond_count],\n  [Assembly], [#meta.assembly.id],\n  [Units], [#meta.structure.unit_count],\n  [Realized visuals], [#meta.representation.realized_visuals.join(\", \")],\n)", title: "Inspect the Model/Structure/Unit result for RCSB PDB entry 9R1O", file: "metadata.typ")

#example-result(
  rendered-metadata-pdf,
  [Selected structure and representation metadata produced from PDB entry 9R1O.],
  "9R1O",
  "https://doi.org/10.2210/pdb9R1O/pdb",
  width: 68%,
  source-note: [Deposition authors: Petrenas, Ozga, Chubb, and Woolfson.],
)

The #ic("render_objects") array can be used for stronger assertions, for
example checking that at least one object has #ic("secondary_type == \"helix\"")
or that #ic("realized_visuals") contains #ic("polymer-trace").

== Assembly And Unit Metadata <sec:assembly-metadata>

Assembly selection is visible through both #ic("info(...).structure") and mesh
export metadata. OBJ and PLY include #ic("molfig_operator_metadata") comments
when an assembly is selected. STL has no equivalent text channel, so assembly
operator metadata should be inspected through @cmd:info[-] or
@cmd:render-object[-].

#table(
  columns: (25%, 35%, 40%),
  inset: 3pt,
  stroke: luma(82%),
  fill: (_, y) => if y == 0 { luma(94%) } else { none },
  [*Path*], [*Example*], [*Meaning*],
  ..(
    row3([#compact-ic("structure.unit_count")], [#compact-ic("meta.structure.unit_count")], [Number of units after structure construction.]),
    row3([#code-lines(("structure.", "symmetry_group_count"))], [#code-lines(("meta.structure.", "symmetry_group_count"))], [Mol\*-style symmetry grouping count.]),
    row3([#code-lines(("structure.", "coordinate_system"))], [#code-lines(("meta.structure.", "coordinate_system"))], [Selected coordinate operator summary.]),
    row3([#compact-ic("structure.boundary")], [#code-lines(("meta.structure.boundary.", "sphere_radius"))], [Boundary sphere and box used by structure-level geometry.]),
    row3([#compact-ic("structure.lookup3d")], [#code-lines(("meta.structure.lookup3d.", "unit_count"))], [Lookup3D summary for constructed units.]),
  ).flatten(),
)

= Practical Workflows <sec:workflows>

== A Publication Figure <sec:publication-figure>

The complete 9R1O Quickstart is the publication-oriented example: it pins the
PDB entry, biological assembly, representation, tessellation quality, camera,
centering, mesh format, and raster output. Its code, rendered PDF, and source
attribution are shown together in the Quickstart section.

== A Lightweight Draft Figure <sec:draft-figure>

#code("typ", "#import \"" + package-import + "\"\n\n// Structural data: RCSB PDB / wwPDB entry 9M1U.\n// https://doi.org/10.2210/pdb9M1U/pdb (CC0 1.0)\n#molfig.render(\n  read(\"9M1U.pdb\", encoding: none),\n  format: \"pdb\",\n  representation: \"cartoon\",\n  assembly: \"1\",\n  quality: \"auto\",\n  mesh-format: \"obj\",\n  output-format: \"png\",\n  config: (elevation: 45, background: \"\"),\n)", title: "Draft a large real assembly from RCSB PDB entry 9M1U", file: "9M1U.typ")

#example-result(
  rendered-draft-pdf,
  [RCSB PDB entry 9M1U rendered as a Cartoon with automatic quality and PNG output.],
  "9M1U",
  "https://doi.org/10.2210/pdb9M1U/pdb",
  width: 68%,
  source-note: [Structure authors: Liu, Zhang, and Xu. Primary citation:
    Zhang et al. (2026), #emph[The EMBO Journal],
    #link("https://doi.org/10.1038/s44318-026-00823-y")[doi:10.1038/s44318-026-00823-y].],
)

== Export For External Tools <sec:external-export>

The complete 1CRN export example in the Export Commands section generates OBJ,
MTL, STL, and PLY bytes from a real BinaryCIF entry and typesets their measured
sizes. Typst does not write arbitrary files from a package call; use those byte
values inside document logic or expose them through a workflow that may write
artifacts outside Typst.

= Troubleshooting <sec:troubleshooting>

#table3([Symptom], [Likely cause], [Action], (
  row3([Plugin says Molfig expects bytes.], [The input was not bytes, inline string data, or a Typst 0.15+ path value.], [For example, pass #ic("path(\"9R1O.pdb\")") on Typst 0.15+, or #ic("read(\"9R1O.pdb\", encoding: none)") for Typst 0.14-compatible documents.]),
  row3([BinaryCIF reports a missing encoding.], [The file is not valid BinaryCIF MessagePack or a column lacks required BinaryCIF encoding metadata.], [Check the source file and use #ic("format: \"cif\"") only for text CIF/mmCIF.]),
  row3([The figure is sparse or missing chains.], [Assembly or altLoc selection filtered the structure.], [Inspect #ic("molfig.info(...).assemblies") and #ic("alt_locs_info"), then set #arg[assembly] and #arg[alt-loc] explicitly.]),
  row3([A large assembly compiles slowly.], [High tessellation or too much assembly geometry.], [Use #ic("quality: \"auto\""), lower #arg[radial-segments] and #arg[linear-segments], or choose a lighter representation.]),
  row3([SVG parsing fails with #ic("\"nodes limit reached\"").], [A high-poly mesh, commonly a large spacefill representation, expands to more SVG nodes than Typst accepts.], [Use #ic("output-format: \"png\""). If vector output is required, reduce #arg[quality], #arg[sphere-detail], #arg[linear-segments], or #arg[radial-segments].]),
  row3([Bonds are missing.], [The file lacks explicit bonds and inference is disabled.], [Leave #arg[infer-bonds] as #value(true), or inspect #ic("bond_metadata") to see available bond sources.]),
  row3([Rendered view differs from external OBJ inspection.], [Different camera, centering, or mesh format.], [Pin #arg[center], #arg[mesh-format], maquette #arg[config], and representation quality options.]),
))

For reproducible documents, follow the 9R1O Quickstart pattern: use a stable
archive identifier, explicit input and mesh formats, a pinned representation,
assembly, quality, centering policy, output format, and camera.

= License And Notices <sec:license-notices>

Molfig package code is distributed under the MIT License. The package also
contains third-party material that is covered separately:

#table3([Material], [Scope], [License / notice], (
  row3([#ic("Mol*")], [#ic("molfig.wasm") contains Rust behavior and generated reference data derived from #ic("Mol*") parity work.], [MIT License; copyright (c) 2017 - now, #ic("Mol*") contributors.]),
  row3([PDB archive data], [#ic("examples/data/*.pdb"), #ic("examples/data/*.cif"), #ic("examples/data/*.bcif"), and generated example PDFs.], [CC0 1.0 Universal Public Domain Dedication; attribution to PDB IDs and structure authors is encouraged.]),
  row3([Typst helper packages], [#ic("maquette") for rendering and #ic("mantys") for this manual.], [Imported by Typst; not vendored by Molfig.]),
))

The package includes #ic("NOTICE.md") and #ic("THIRD_PARTY_NOTICES.md") with
the full distribution notice. Example data attributions are also listed in
#ic("examples/data/README.md") and #ic("examples/data/ATTRIBUTION.tsv").

For figures made from PDB archive data, cite the PDB ID and DOI. For example:

#code("text", "Structural data source: RCSB PDB / wwPDB, PDB ID 9R1O,\nhttps://doi.org/10.2210/pdb9R1O/pdb.\nPDB archive data files are available under CC0 1.0.", title: "Suggested data attribution")

When describing #ic("Mol*") parity or comparing against #ic("Mol*") output, cite #ic("Mol*"):

#code("text", "Sehnal et al. Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures.\nNucleic Acids Research 49:W431-W437 (2021).\nhttps://doi.org/10.1093/nar/gkab314", title: "Suggested Mol* citation")

No endorsement by #ic("Mol*"), RCSB PDB, wwPDB, structure authors, or data
providers is implied.

= Development <sec:development>

The package consists of a Typst wrapper, a Rust WebAssembly plugin, tests, and
documentation:

#table2([Path], [Role], (
  row2([#ic("package/lib.typ")], [Public Typst API and maquette bridge.]),
  row2([#ic("wasm-plugin/src/")], [Rust parser, Model/Structure/Unit layer, geometry builders, exporters, and Typst plugin ABI.]),
  row2([#ic("wasm-plugin/tests/")], [Rust and Typst smoke/regression tests.]),
  row2([#ic("package/docs/documentation.typ")], [This Mantys manual.]),
  row2([#ic("package/docs/documentation.pdf")], [Compiled reference manual.]),
  row2([#ic("package/molfig.wasm")], [Checked-in WebAssembly plugin consumed by Typst.]),
))

#shell("cd wasm-plugin\ncargo fmt --check\ncargo test\ncargo build --release --target wasm32-unknown-unknown\ncp target/wasm32-unknown-unknown/release/molfig.wasm ../package/molfig.wasm\ncd ..\ntypst compile --root package package/docs/documentation.typ package/docs/documentation.pdf")

Regenerate #ic("package/molfig.wasm") after Rust changes that affect the Typst plugin.
Regenerate this PDF after documentation or public API changes.
