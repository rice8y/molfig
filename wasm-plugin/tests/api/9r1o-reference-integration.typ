// Compile-time integration test for generating and rendering Molfig's 9R1O OBJ.

#import "../../../package/lib.typ" as molfig
#import "@preview/maquette:0.2.0": render-obj

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

#render-obj(
  obj,
  (
    azimuth: 25,
    elevation: 18,
    background: "",
  ),
  format: "svg",
)
