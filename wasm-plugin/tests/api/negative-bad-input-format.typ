// Negative smoke test.
// Expected result: Typst compilation fails with an actionable message that
// mentions the accepted input formats: pdb, cif.

#import "../../../package/lib.typ" as molfig

#let water = read("../fixtures/pdb/water.pdb", encoding: none)
#molfig.to-obj(water, format: "xyz")
