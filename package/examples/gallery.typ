// Gallery of Molfig example renders.
// Structural data source: RCSB PDB / wwPDB.
// PDB IDs: 1CRN, 1FYY, 9M1U, 9Q12, 9R1O, and 9Z4O.
// PDB archive data files are available under CC0 1.0.

#import "@preview/molfig:0.2.0"

#set page(
  paper: "a4",
  flipped: true,
  margin: 4mm,
)
#set text(font: "New Computer Modern", size: 10pt)

#let entry(data, format) = block(
  width: 93mm,
  height: 99mm,
  clip: true,
  align(
    center + horizon,
    molfig.render(
      data,
      format: format,
      representation: "cartoon",
      quality: "high",
      renderer: (
        viewport: (width: 560, height: 594),
        camera: (view: (name: "orbit", params: (elevation: 45)),),
        background: (color: "#ffffff", transparent: true),
      ),
      width: 93mm,
      height: 99mm,
    ),
  ),
)

#align(
  center + horizon,
  grid(
    columns: (93mm, 93mm, 93mm),
    rows: (99mm, 99mm),
    column-gutter: 5mm,
    row-gutter: 4mm,
    entry(read("data/1crn.bcif", encoding: none), "bcif"),
    entry(read("data/1FYY.cif", encoding: none), "cif"),
    entry(read("data/9M1U.pdb", encoding: none), "pdb"),
    entry(read("data/9q12.pdb", encoding: none), "pdb"),
    entry(read("data/9R1O.pdb", encoding: none), "pdb"),
    entry(read("data/9Z4O.pdb", encoding: none), "pdb"),
  ),
)
