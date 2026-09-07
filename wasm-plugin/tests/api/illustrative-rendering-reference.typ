// Mol* Quick Style comparison fixture for the native renderer.
// Page 1: XYZ Ball & Stick. Page 2: PDB Cartoon.

#import "../../../package/lib.typ" as molfig

#set page(width: 1024pt, height: 937pt, margin: 0pt, fill: rgb("fcfbfa"))

#let render-reference(data, format, representation, assembly, position, target, radius) = align(
  center + horizon,
  molfig.render(
    data,
    format: format,
    representation: representation,
    assembly: assembly,
    color-theme: if format == "xyz" { "element-symbol" } else { "chain-id" },
    style: "illustrative",
    quality: "high",
    renderer: (
      viewport: (width: 1024, height: 937),
      camera: (
        view: (
          name: "snapshot",
          params: (
            position: position,
            target: target,
            up: (0, 1, 0),
            radius: radius,
            radius-max: radius,
          ),
        ),
      ),
      background: (color: "#fcfbfa"),
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
  (-0.03077494221759732, 0.03827120753606447, 21.39274083202126),
  (-0.03077494221759732, 0.03827120753606447, 0.00003476357507388535),
  4.186634185850473,
)

#pagebreak()

#render-reference(
  read("../../../package/examples/data/9Z4O.pdb", encoding: none),
  "pdb",
  "cartoon",
  "1",
  (131.63666556103, 125.63926427158606, 292.5377510373967),
  (131.63666556103, 125.63926427158606, 135.88791972158091),
  55.94729512734871,
)
