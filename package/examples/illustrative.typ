// Structural data: RCSB PDB / wwPDB entry 9Z4O (CC0 1.0).
// https://doi.org/10.2210/pdb9Z4O/pdb
#import "@preview/molfig:0.1.4"
#set page(width: 124mm, height: auto, margin: 4mm)
#set text(font: "New Computer Modern", size: 9pt)
#let data = path("data/9Z4O.pdb")
#let view(style) = molfig.render(
  data,
  format: "pdb",
  representation: "cartoon",
  assembly: "1",
  color-theme: "chain-id",
  style: style,
  mesh-format: "obj",
  quality: "high",
  center: true,
  output-format: "png",
  config: (up: (0, 1, 0), azimuth: 35, elevation: 24, background: "", antialias: 4),
  width: 55mm,
  height: 50mm,
)
#let panel(label, style) = [
  #align(center, strong(label))
  #block(
    width: 55mm,
    height: 50mm,
    clip: true,
    align(center + horizon, scale(
      x: 145%, y: 145%, origin: center + horizon, view(style),
    )),
  )
]
#grid(
  columns: (1fr, 1fr),
  column-gutter: 4mm,
  panel([Default], "default"),
  panel([Illustrative], "illustrative"),
)
