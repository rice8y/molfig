// PubChem CID 241 (benzene), PubChem3D conformer 000000F100000001.
// Record: https://pubchem.ncbi.nlm.nih.gov/compound/241
// 3D SDF retrieved through PubChem PUG REST on 2026-08-24, then
// transcribed to XYZ without changing atom order or coordinates.
#import "@preview/molfig:0.2.0"
#set page(width: 124mm, height: auto, margin: 4mm)
#set text(font: "New Computer Modern", size: 9pt)
#let xyz = path("data/benzene.xyz")
#let view(style) = molfig.render(
  xyz,
  format: "xyz",
  representation: "default",
  color-theme: "element-symbol",
  style: style,
  quality: "high",
  renderer: (
    viewport: (width: 550, height: 500),
    camera: (view: (name: "orbit", params: (azimuth: 35, elevation: 80)),),
    background: (color: "#ffffff", transparent: true),
  ),
  width: 55mm,
  height: 50mm,
)
#let panel(label, style) = [
  #align(center, strong(label))
  #block(
    width: 55mm,
    height: 50mm,
    clip: true,
    align(center + horizon, view(style)),
  )
]
#grid(
  columns: (1fr, 1fr),
  column-gutter: 4mm,
  panel([Default], "default"),
  panel([Illustrative], "illustrative"),
)
