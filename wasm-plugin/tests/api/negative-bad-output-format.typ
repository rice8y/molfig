// Negative smoke test.
// Expected result: Typst compilation fails with an actionable message that
// mentions the accepted output formats: obj, stl, ply.

#import "../../../package/lib.typ" as molfig

#let water = read("../fixtures/pdb/water.pdb", encoding: none)
#molfig.render(water, format: "pdb", mesh-format: "glb")
