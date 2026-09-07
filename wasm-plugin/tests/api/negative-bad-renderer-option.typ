// Negative compile-failure test.
// Expected result: Typst compilation fails and identifies the unknown renderer key.

#import "../../../package/lib.typ" as molfig

#let water = read("../fixtures/pdb/water.pdb", encoding: none)
#molfig.render(water, format: "pdb", renderer: (unknown-pass: true))
