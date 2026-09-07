// Negative compile-failure test.
// Expected result: Typst compilation fails and identifies SVG and PNG as the
// accepted rendering output formats.

#import "../../../package/lib.typ" as molfig

#let water = read("../fixtures/pdb/water.pdb", encoding: none)
#molfig.render(water, format: "pdb", output-format: "obj")
