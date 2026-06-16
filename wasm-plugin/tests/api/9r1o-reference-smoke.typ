// Compile-time smoke test for rendering the Mol* reference OBJ with maquette.

#import "@preview/maquette:0.1.0": render-obj

#let obj = read("../../../package/examples/data/9R1O.obj", encoding: none)

#assert(obj.len() > 1000000)

#render-obj(
  obj,
  json.encode((
    azimuth: 25,
    elevation: 18,
    background: "",
  )),
  format: "svg",
)
