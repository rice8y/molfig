#import "@preview/mantys:1.0.2": *

#let manifest = toml(read("../typst.toml", encoding: none))
#let package-id = manifest.package.name
#let package-version = manifest.package.version
#let product-name = "Molfig"
#let package-import = "@preview/" + package-id + ":" + package-version
#let rendered-9r1o-pdf = "../examples/9R1O.pdf"
#let example-9r1o-code = "#import \"" + package-import + "\" as molfig\n\n#set page(width: auto, height: auto, margin: 0mm)\n\n// Uses structural data from RCSB PDB / wwPDB.\n// PDB ID: 9R1O\n// PDB DOI: https://doi.org/10.2210/pdb9R1O/pdb\n// PDB archive data files are available under CC0 1.0.\n#let pdb = if molfig.v15-or-later() {\n  path(\"9R1O.pdb\")\n} else {\n  read(\"9R1O.pdb\", encoding: none)\n}\n\n#molfig.render(\n  pdb,\n  format: \"pdb\",\n  representation: \"molstar\",\n  assembly: \"1\",\n  mesh-format: \"obj\",\n  quality: \"high\",\n  center: true,\n  output-format: \"png\",\n  config: (\n    azimuth: 35,\n    elevation: 24,\n    zoom: 1.0,\n    background: \"\",\n  ),\n)"

#let ic(value) = raw(str(value))

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
  #product-name is not a Mol\* WebGL renderer. It intentionally excludes
  interactive WebGL surface and volume visuals such as molecular surface,
  gaussian surface, gaussian volume, and density maps. Its public contract is
  structure parsing, static mesh generation, OBJ/STL/PLY export, and maquette
  rendering.
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

Compile the example first when you want to refresh the figure embedded below:

#shell("typst compile 9R1O.typ")

#figure(
  block(
    width: 100%,
    inset: 4pt,
    stroke: luma(82%),
    radius: 2pt,
    image(rendered-9r1o-pdf, width: 100%),
  ),
  caption: [
    9R1O rendered by #product-name from PDB entry 9R1O.
  ],
) <fig:9r1o-example>

#info-alert[
  Structural data source: RCSB PDB / wwPDB, PDB ID #ic("9R1O"),
  DOI #ic("10.2210/pdb9R1O/pdb"). PDB archive data files are distributed under
  CC0 1.0.
]

== Reading Structure Files <sec:reading-files>

On Typst 0.15.0 or later, prefer passing a #ic("path(...)") value for
file-backed structure data. The path is resolved by Typst relative to the caller's
file, and #product-name reads it internally with #ic("encoding: none"). This keeps
package calls concise while preserving the caller-side project boundary.

#code("typ", "#let pdb = path(\"entry.pdb\")\n#let cif = path(\"entry.cif\")\n#let bcif = path(\"entry.bcif\")", title: "Typst 0.15+ path inputs")

For documents that must also compile on Typst 0.14, read the file in the caller
document and pass the resulting bytes. Always use #ic("encoding: none") for
external structure files; it preserves BinaryCIF bytes and avoids Unicode
decoding loss for fixed-width PDB columns.

#code("typ", "#let pdb = read(\"entry.pdb\", encoding: none)\n#let cif = read(\"entry.cif\", encoding: none)\n#let bcif = read(\"entry.bcif\", encoding: none)", title: "Typst 0.14-compatible byte inputs")

Passing a Typst string is accepted for small inline examples. A string is treated
as inline molecular text, not as a file path.

== Choosing The First Options <sec:first-options>

#table3([Question], [Recommended setting], [Why], (
  row3([I want the closest default molecular figure.], [#ic("representation: \"molstar\"")], [Uses the package's Mol\*-style static visual set.]),
  row3([I want a protein cartoon.], [#ic("representation: \"cartoon\"")], [Uses polymer trace, secondary structure, nucleotide, and gap visuals.]),
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
  row3([@cmd:v15-or-later[-]], [bool], [You want one source file that uses #ic("path(...)") on Typst 0.15+ and byte input on older Typst.]),
))

== Rendering Commands <sec:rendering-commands>

#command(
  "render",
  arg("data"),
  arg("format", _value: values.value.with("auto")),
  arg("mesh-format", _value: values.value.with("obj")),
  arg("representation", _value: values.value.with("molstar")),
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
    JSON-encoded and passed to maquette unchanged. Use it for camera,
    background, image, and render settings supported by maquette.
  ]

  #argument("width", default: "auto", command: "render")[
    Width forwarded to the maquette rendering command.
  ]

  #argument("height", default: "auto", command: "render")[
    Height forwarded to the maquette rendering command.
  ]

  #argument("output-format", choices: ("png", "svg"), default: "png", command: "render")[
    Output format forwarded to maquette. PNG is the safest default for complex
    meshes.
  ]
]

#command("render-object", arg("data"), arg("mesh-format", _value: values.value.with("obj")))[
  Returns a dictionary with:

  #table2([Key], [Value], (
    row2([#ic("kind")], [The string #ic("\"render-object\"").]),
    row2([#ic("format")], [The selected mesh format.]),
    row2([#ic("mesh_format")], [The selected mesh format with snake_case naming.]),
    row2([#ic("mesh")], [Generated OBJ/STL/PLY bytes.]),
    row2([#ic("info")], [The same dictionary returned by @cmd:info[-] for the same mesh options.]),
    row2([#ic("content")], [The Typst content returned by @cmd:render[-].]),
  ))
]

#code("typ", "#let object = molfig.render-object(\n  data,\n  format: \"mmcif\",\n  representation: \"cartoon\",\n  mesh-format: \"ply\",\n  assembly: \"1\",\n  config: (azimuth: 30, elevation: 18, background: \"\"),\n  width: 54mm,\n)\n\n#object.content\n\n#metadata((\n  mesh-format: object.mesh_format,\n  atom-count: object.info.atom_count,\n))", title: "Render and inspect in one call")

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

#code("typ", "#let data = read(\"structure.bcif\", encoding: none)\n\n#let obj = molfig.to-obj(data, format: \"bcif\", representation: \"molstar\")\n#let mtl = molfig.to-mtl(data, format: \"bcif\", representation: \"molstar\")\n#let stl = molfig.to-stl(data, format: \"bcif\", representation: \"molstar\")\n#let ply = molfig.to-ply(data, format: \"bcif\", representation: \"molstar\")", title: "Direct mesh export")

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

== Compatibility Helper <sec:compatibility-helper>

#command("v15-or-later")[
  Returns #value(true) when #ic("sys.version >= version(0, 15, 0)"). Use it in
  source files that should prefer Typst 0.15+ project paths while remaining
  compilable on Typst 0.14.
]

#code("typ", "#let pdb = if molfig.v15-or-later() {\n  path(\"entry.pdb\")\n} else {\n  read(\"entry.pdb\", encoding: none)\n}", title: "Version-gated file input")

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
  option-row("representation", ic("\"molstar\""), [Static Mol\*-style default visual set. Also accepts #ic("\"default\"") and #ic("\"auto\"").]),
))

#table4([Value], [Best for], [Main geometry], [Notes], (
  row4([#ic("\"molstar\"")], [General figures], [Polymer traces, element spheres, bonds, carbohydrates, nucleotides, gaps], [Closest default for static Molfig output.]),
  row4([#ic("\"spacefill\"")], [Atomic packing], [Atom spheres], [Uses atom radii and #arg[radius-scale].]),
  row4([#ic("\"ball-and-stick\"")], [Ligands and small molecules], [Atom spheres plus bond cylinders], [Good for local chemistry and explicit bond inspection.]),
  row4([#ic("\"cartoon\"")], [Proteins and nucleic acids], [Polymer trace tubes/sheets/ribbons, nucleotide visuals, gaps], [Uses secondary structure and polymer trace metadata.]),
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
The #arg[config] dictionary is JSON-encoded and passed through unchanged.

#code("typ", "#let cif = read(\"structure.cif\", encoding: none)\n\n#molfig.render(\n  cif,\n  format: \"cif\",\n  representation: \"ball-and-stick\",\n  mesh-format: \"ply\",\n  config: (\n    azimuth: 45,\n    elevation: 20,\n    zoom: 1.15,\n    background: \"\",\n  ),\n  width: 75mm,\n)", title: "Camera and background passthrough")

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
  row3([#ic("render_objects")], [array], [Semantic render-object spans with geometry type, visual, chain, residue range, group id, and value-cell style counts.]),
  row3([#ic("bond_metadata")], [dictionary], [Counts by bond source and flags such as struct_conn, index_pair, chem_comp, aromatic, and resonance.]),
  row3([#ic("bounds")], [dictionary], [Expanded coordinate bounds before export centering.]),
))

#code("typ", "#let pdb = read(\"9R1O.pdb\", encoding: none)\n#let meta = molfig.info(pdb, format: \"pdb\", assembly: \"1\")\n\nAtoms: #meta.atom_count\nBonds: #meta.bond_count\nAssembly: #meta.assembly.id", title: "Basic 9R1O metadata query")

Render-object metadata is useful when validating what Molfig actually built:

#code("typ", "#let meta = molfig.info(\n  data,\n  format: \"mmcif\",\n  representation: \"cartoon\",\n  assembly: \"1\",\n  alt-loc: \"highest-occupancy\",\n)\n\n#assert(meta.render_objects.any(object => object.secondary_type == \"helix\"))\n#assert(meta.representation.realized_visuals.contains(\"polymer-trace\"))", title: "Representation assertions")

== Assembly And Unit Metadata <sec:assembly-metadata>

Assembly selection is visible through both #ic("info(...).structure") and mesh
export metadata. OBJ and PLY include #ic("molfig_operator_metadata") comments
when an assembly is selected. STL has no equivalent text channel, so assembly
operator metadata should be inspected through @cmd:info[-] or
@cmd:render-object[-].

#table3([Path], [Example], [Meaning], (
  row3([#ic("structure.unit_count")], [#ic("meta.structure.unit_count")], [Number of units after structure construction.]),
  row3([#ic("structure.symmetry_group_count")], [#ic("meta.structure.symmetry_group_count")], [Mol\*-style symmetry grouping count.]),
  row3([#ic("structure.coordinate_system")], [#ic("meta.structure.coordinate_system")], [Selected coordinate operator summary.]),
  row3([#ic("structure.boundary")], [#ic("meta.structure.boundary.sphere_radius")], [Boundary sphere and box used by structure-level geometry.]),
  row3([#ic("structure.lookup3d")], [#ic("meta.structure.lookup3d.unit_count")], [Lookup3D summary for constructed units.]),
))

= Practical Workflows <sec:workflows>

== A Publication Figure <sec:publication-figure>

#code("typ", "#import \"" + package-import + "\" as molfig\n\n#let pdb = read(\"9R1O.pdb\", encoding: none)\n\n#figure(\n  molfig.render(\n    pdb,\n    format: \"pdb\",\n    representation: \"molstar\",\n    assembly: \"1\",\n    mesh-format: \"obj\",\n    quality: \"high\",\n    center: true,\n    config: (\n      azimuth: 35,\n      elevation: 24,\n      background: \"\",\n    ),\n    width: 90mm,\n  ),\n  caption: [9R1O rendered with #product-name.],\n)", title: "Pinned 9R1O figure")

== A Lightweight Draft Figure <sec:draft-figure>

#code("typ", "#molfig.render(\n  read(\"data/large-assembly.bcif\", encoding: none),\n  format: \"bcif\",\n  representation: \"molstar\",\n  assembly: \"1\",\n  quality: \"auto\",\n  mesh-format: \"stl\",\n  config: (background: \"\"),\n  width: 60mm,\n)", title: "Large assembly draft")

== Export For External Tools <sec:external-export>

#code("typ", "#let data = read(\"data/ligand.cif\", encoding: none)\n#let mesh = molfig.to-obj(\n  data,\n  format: \"cif\",\n  representation: \"ball-and-stick\",\n  assembly: \"asymmetric-unit\",\n  center: false,\n)\n#let material = molfig.to-mtl(data, format: \"cif\", representation: \"ball-and-stick\")", title: "OBJ and MTL bytes")

Typst does not write arbitrary files from a package call. Use these byte values
inside Typst document logic, or expose them through a workflow that is allowed
to write artifacts outside Typst.

= Troubleshooting <sec:troubleshooting>

#table3([Symptom], [Likely cause], [Action], (
  row3([Plugin says Molfig expects bytes.], [The input was not bytes, inline string data, or a Typst 0.15+ path value.], [Pass #ic("path(\"entry.pdb\")") on Typst 0.15+, or pass #ic("read(\"entry.pdb\", encoding: none)") for Typst 0.14-compatible documents.]),
  row3([BinaryCIF reports a missing encoding.], [The file is not valid BinaryCIF MessagePack or a column lacks required BinaryCIF encoding metadata.], [Check the source file and use #ic("format: \"cif\"") only for text CIF/mmCIF.]),
  row3([The figure is sparse or missing chains.], [Assembly or altLoc selection filtered the structure.], [Inspect #ic("molfig.info(...).assemblies") and #ic("alt_locs_info"), then set #arg[assembly] and #arg[alt-loc] explicitly.]),
  row3([A large assembly compiles slowly.], [High tessellation or too much assembly geometry.], [Use #ic("quality: \"auto\""), lower #arg[radial-segments] and #arg[linear-segments], or choose a lighter representation.]),
  row3([Bonds are missing.], [The file lacks explicit bonds and inference is disabled.], [Leave #arg[infer-bonds] as #value(true), or inspect #ic("bond_metadata") to see available bond sources.]),
  row3([Rendered view differs from external OBJ inspection.], [Different camera, centering, or mesh format.], [Pin #arg[center], #arg[mesh-format], maquette #arg[config], and representation quality options.]),
))

For reproducible documents, avoid implicit choices:

#code("typ", "#let pdb = read(\"9R1O.pdb\", encoding: none)\n\n#molfig.render(\n  pdb,\n  format: \"pdb\",\n  representation: \"molstar\",\n  assembly: \"1\",\n  mesh-format: \"obj\",\n  quality: \"custom\",\n  sphere-detail: 2,\n  linear-segments: 8,\n  radial-segments: 16,\n  center: true,\n  config: (azimuth: 35, elevation: 24, background: \"\"),\n)", title: "Reproducible 9R1O option block")

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
