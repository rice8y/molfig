// Uses structural data from RCSB PDB / wwPDB.
// PDB ID: 1FYY
// PDB DOI: https://doi.org/10.2210/pdb1FYY/pdb
// Deposition authors: Volk, D.E.; Rice, J.S.; Luxon, B.A.; Yeh, H.J.C.;
// Liang, C.; Xie, G.; Sayer, J.M.; Jerina, D.M.; Gorenstein, D.G.
// PDB archive data files are available under CC0 1.0.
// Primary citation: Volk, D.E. et al. (2000) Biochemistry 39: 14040-14053.
// Article DOI: https://doi.org/10.1021/bi001669l

#import "@preview/molfig:0.2.0"

#set page(width: 100mm, height: auto, margin: 4mm)
#set text(font: "New Computer Modern", size: 9pt)

#let result = molfig.render-result(
  read("data/1FYY.cif", encoding: none),
  format: "cif",
  representation: "cartoon",
  assembly: "1",
  renderer: (
    viewport: (width: 920, height: 640),
    camera: (view: (name: "orbit", params: (azimuth: 30, elevation: 18)),),
    background: (color: "#ffffff", transparent: true),
  ),
  width: 92mm,
  height: 64mm,
)

#result.content
#align(center)[
  #text(size: 8pt)[Atoms: #result.info.atom_count \ 
    Pixels: #result.pixel-width × #result.pixel-height]
]
