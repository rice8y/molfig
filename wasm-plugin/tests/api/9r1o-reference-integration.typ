// Compile-time integration test for export and native rendering of 9R1O.

#import "../../../package/lib.typ" as molfig

#let pdb = read("../../../package/examples/data/9R1O.pdb", encoding: none)
#let obj = molfig.to-obj(
  pdb,
  format: "pdb",
  representation: "cartoon",
  assembly: "1",
  sphere-detail: 1,
  quality: "auto",
)

#assert(obj.len() > 10000000)

#molfig.render(
  pdb,
  format: "pdb",
  representation: "cartoon",
  assembly: "1",
  quality: "auto",
  renderer: (
    viewport: (width: 320, height: 240),
    camera: (view: (name: "orbit", params: (azimuth: 25, elevation: 18)),),
    background: (color: "#ffffff", transparent: true),
  ),
)
