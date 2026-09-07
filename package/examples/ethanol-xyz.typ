#import "@preview/molfig:0.2.0"

#set page(width: 82mm, height: 82mm, margin: 4mm)

// PubChem CID 702 (ethanol), PubChem3D conformer 000002BE00000001.
// Record: https://pubchem.ncbi.nlm.nih.gov/compound/702
// 3D SDF retrieved through PubChem PUG REST on 2026-08-24, then
// transcribed to XYZ without changing atom order or coordinates.
#let xyz = path("data/ethanol.xyz")

#molfig.render(
  xyz,
  format: "xyz",
  representation: "default",
  color-theme: "element-symbol",
  quality: "high",
  renderer: (
    viewport: (width: 740, height: 740),
    camera: (view: (name: "orbit", params: (azimuth: 125, elevation: 40)),),
    background: (color: "#ffffff", transparent: true),
  ),
  width: 74mm,
)
