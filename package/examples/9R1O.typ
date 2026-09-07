#import "@preview/molfig:0.2.0"

#set page(width: auto, height: auto, margin: 0mm)

// Uses structural data from RCSB PDB / wwPDB.
// PDB ID: 9R1O
// PDB DOI: https://doi.org/10.2210/pdb9R1O/pdb
// Deposition authors: Petrenas, R.; Ozga, K.; Chubb, J.J.; Woolfson, D.N.
// PDB archive data files are available under CC0 1.0.
// Literature status: To be published.
#let pdb = read("data/9R1O.pdb", encoding: none)

#molfig.render(
  pdb,
  format: "pdb",
  representation: "cartoon",
  assembly: "1",
  quality: "high",
  renderer: (
    viewport: (width: 960, height: 720),
    camera: (view: (name: "orbit", params: (azimuth: 35, elevation: 24)),),
    background: (color: "#ffffff", transparent: true),
  ),
)
