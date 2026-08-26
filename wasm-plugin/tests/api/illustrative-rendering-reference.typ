// Pixel-level Mol* Quick Style comparison fixture.
// Page 1: XYZ Ball & Stick. Page 2: PDB Cartoon.

#import "../../../package/lib.typ" as molfig

#set page(width: 1024pt, height: 937pt, margin: 0pt, fill: rgb("fcfbfa"))

#let render-reference(data, format, representation, assembly, distance) = align(
  center + horizon,
  molfig.render(
    data,
    format: format,
    representation: representation,
    assembly: assembly,
    color-theme: if format == "xyz" { "element-symbol" } else { "chain-id" },
    style: "illustrative",
    mesh-format: "obj",
    quality: "high",
    center: true,
    output-format: "png",
    config: (
      width: 1024,
      height: 937,
      // Mol* starts at +Z looking toward the origin, with +Y screen-up.
      // Use that camera explicitly instead of approximating it with angles.
      // With +Y as up, Maquette's azimuth 180° points the camera along +Z.
      // Keeping the spherical camera also preserves its radius-scaled distance.
      azimuth: 180,
      elevation: 0,
      up: (0, 1, 0),
      distance: distance,
      auto_center: false,
      center: if format == "xyz" {
        (-0.0253810993447476, 0.0171381860213884, -0.00097060000330877)
      } else {
        (2.65456386259908, -3.08523388114679, -1.16618677096460)
      },
      background: "#fcfbfa",
      antialias: 4,
    ),
    width: 1024pt,
    height: 937pt,
  ),
)

#render-reference(
  read("../../../package/examples/data/benzene.xyz", encoding: none),
  "xyz",
  "default",
  "asymmetric-unit",
  21.3927060684462,
)

#pagebreak()

#render-reference(
  read("../../../package/examples/data/9Z4O.pdb", encoding: none),
  "pdb",
  "cartoon",
  "1",
  156.649831315816,
)
