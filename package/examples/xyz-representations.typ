// PubChem CID 241 (benzene), PubChem3D conformer 000000F100000001.
// Record: https://pubchem.ncbi.nlm.nih.gov/compound/241
// 3D SDF retrieved through PubChem PUG REST on 2026-08-24, then
// transcribed to XYZ without changing atom order or coordinates.
#import "@preview/molfig:0.1.4"
#set page(width: 180mm, height: auto, margin: 4mm)
#set text(font: "New Computer Modern", size: 9pt)
#let xyz = path("data/benzene.xyz")
#let view(representation) = molfig.render(
  xyz,
  format: "xyz",
  representation: representation,
  style: "default",
  mesh-format: "obj",
  quality: "high",
  center: true,
  output-format: "png",
  config: (azimuth: 35, elevation: 80, background: "", antialias: 4),
  width: 54mm,
  height: 50mm,
)
#let panel(label, representation) = [
  #align(center, strong(label))
  #block(
    width: 54mm,
    height: 50mm,
    clip: true,
    align(center + horizon, view(representation)),
  )
]
#grid(
  columns: (1fr, 1fr, 1fr),
  column-gutter: 5mm,
  // Mol* Viewer Default resolves small XYZ structures to Ball & Stick.
  panel([Ball-and-stick], "default"),
  panel([Spacefill], "spacefill"),
  panel([Surface], "surface"),
)
