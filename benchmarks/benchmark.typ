#let manifest = toml(read("../package/typst.toml", encoding: none))
#let molfig-version = sys.inputs.at("molfig-version", default: manifest.package.version)
#let molfig-namespace = sys.inputs.at("molfig-namespace", default: "preview")
#let package-spec = "@" + molfig-namespace + "/molfig:" + molfig-version
#import package-spec as molfig

#let cases = (
  (
    id: "1crn-bcif-spacefill",
    file: "../package/examples/data/1crn.bcif",
    format: "bcif",
    representation: "spacefill",
    quality: "high",
  ),
  (
    id: "1fyy-cif-surface",
    file: "../package/examples/data/1FYY.cif",
    format: "cif",
    representation: "surface",
    quality: "high",
  ),
  (
    id: "9r1o-pdb-cartoon",
    file: "../package/examples/data/9R1O.pdb",
    format: "pdb",
    representation: "cartoon",
    quality: "high",
  ),
  (
    id: "9z4o-pdb-spacefill",
    file: "../package/examples/data/9Z4O.pdb",
    format: "pdb",
    representation: "spacefill",
    quality: "high",
  ),
  (
    id: "9m1u-pdb-cartoon-auto",
    file: "../package/examples/data/9M1U.pdb",
    format: "pdb",
    representation: "cartoon",
    quality: "auto",
  ),
  (
    id: "9q12-pdb-cartoon",
    file: "../package/examples/data/9q12.pdb",
    format: "pdb",
    representation: "cartoon",
    quality: "high",
  ),
)

#let case-id = sys.inputs.at("case", default: "9r1o-pdb-cartoon")
#let mode = sys.inputs.at("mode", default: "export")
#let selected = cases.find(case => case.id == case-id)

#if selected == none {
  panic("unknown benchmark case: " + case-id)
}

#if mode != "export" and mode != "render" {
  panic("benchmark mode must be either export or render")
}

#let source = read(selected.file, encoding: none)

#set page(width: 80mm, height: 80mm, margin: 0mm)
#set text(size: 1pt, fill: white)

#if mode == "export" {
  let mesh = molfig.to-obj(
    source,
    format: selected.format,
    representation: selected.representation,
    quality: selected.quality,
    assembly: "1",
    center: true,
  )

  // Materialize the result so the exporter remains part of the measured work.
  [#mesh.len()]
} else {
  molfig.render(
    source,
    format: selected.format,
    representation: selected.representation,
    quality: selected.quality,
    assembly: "1",
    center: true,
    output-format: "png",
    width: 80mm,
    height: 80mm,
    config: (
      azimuth: 35,
      elevation: 24,
      background: "",
    ),
  )
}
