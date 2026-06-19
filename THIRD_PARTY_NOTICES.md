# Third-Party Notices

This file records third-party material that is included in, derived by, or redistributed with Molfig. It is intended to travel with both source releases and package releases.

Molfig itself is licensed under the MIT License. See `LICENSE`.

## Mol*

Project: Mol* (`molstar/molstar`)

Source: https://github.com/molstar/molstar

Reference commit used for parity work: `1b8117d3f10f7c978aabb5a0d3d47370635aefe4`

License: MIT License

Copyright notice: Copyright (c) 2017 - now, Mol* contributors

Scope in this repository:

- `wasm-plugin/src/` contains Rust code that ports, adapts, or behaviorally matches Mol* Model/Structure/Unit, representation, geometry, and geo-exporter behavior.
- `wasm-plugin/src/model/reference_data/*.ts` contains Mol* generated reference data files and constants used by the Rust model classification layer.
- `package/molfig.wasm` is built from the Rust implementation and therefore includes Mol*-derived behavior and reference data.
- `wasm-plugin/artifacts/molstar/` is a local pinned Mol* checkout for parity development. It is ignored by git and is not intended to be published as part of the Molfig Typst package.

Mol* MIT License text:

```text
The MIT License

Copyright (c) 2017 - now, Mol* contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

Recommended Mol* citation when describing the parity target or generated figures:

David Sehnal, Sebastian Bittrich, Mandar Deshpande, Radka Svobodova, Karel Berka, Vaclav Bazgier, Sameer Velankar, Stephen K. Burley, Jaroslav Koca, and Alexander S. Rose. Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures. Nucleic Acids Research 49:W431-W437 (2021). https://doi.org/10.1093/nar/gkab314

## PDB Archive Example Data

Source: RCSB PDB / wwPDB PDB archive

RCSB PDB policies: https://www.rcsb.org/pages/policies

CC0 1.0: https://creativecommons.org/publicdomain/zero/1.0/

License/dedication: CC0 1.0 Universal Public Domain Dedication

Scope in this repository:

- `package/examples/data/1crn.bcif`
- `package/examples/data/1FYY.cif`
- `package/examples/data/9M1U.pdb`
- `package/examples/data/9q12.pdb`
- `package/examples/data/9R1O.pdb`
- `package/examples/data/9Z4O.pdb`
- `package/examples/*.typ` files that load these data files
- `package/examples/*.pdf` files generated from these data files
- Documentation examples that display or cite these data files

RCSB PDB states that data files in the PDB archive are available under CC0 1.0.
RCSB PDB also encourages attribution to the original structure-data authors where possible. CC0 does not imply endorsement by the authors, RCSB PDB, wwPDB, or Creative Commons, and the data are provided without warranty.

Example data attribution:

| File | PDB ID | Format | PDB DOI | Structure authors / status |
| --- | --- | --- | --- | --- |
| `package/examples/data/1crn.bcif` | 1CRN | BinaryCIF | https://doi.org/10.2210/pdb1CRN/pdb | Hendrickson, W.A.; Teeter, M.M. Primary citation: Teeter, M.M. (1984) Proc Natl Acad Sci U S A 81:6014-6018. Article DOI: https://doi.org/10.1073/pnas.81.19.6014 |
| `package/examples/data/1FYY.cif` | 1FYY | mmCIF | https://doi.org/10.2210/pdb1FYY/pdb | Volk, D.E.; Rice, J.S.; Luxon, B.A.; Yeh, H.J.C.; Liang, C.; Xie, G.; Sayer, J.M.; Jerina, D.M.; Gorenstein, D.G. Primary citation: Biochemistry 39:14040-14053 (2000). Article DOI: https://doi.org/10.1021/bi001669l |
| `package/examples/data/9M1U.pdb` | 9M1U | PDB | https://doi.org/10.2210/pdb9M1U/pdb | Liu, H.; Zhang, X.; Xu, H.E. Primary citation: Zhang, X. et al. (2026), EMBO J. Article DOI: https://doi.org/10.1038/s44318-026-00823-y |
| `package/examples/data/9q12.pdb` | 9Q12 | PDB | https://doi.org/10.2210/pdb9Q12/pdb | Wang, Y.; Liu, B.; He, Y.; Feigon, J. Literature status in the included PDB file: to be published. |
| `package/examples/data/9R1O.pdb` | 9R1O | PDB | https://doi.org/10.2210/pdb9R1O/pdb | Petrenas, R.; Ozga, K.; Chubb, J.J.; Woolfson, D.N. Literature status in the included PDB file: to be published. |
| `package/examples/data/9Z4O.pdb` | 9Z4O | PDB | https://doi.org/10.2210/pdb9Z4O/pdb | Ge, Y.; de Almeida Magalhaes, T.; Wu, H.; Yadav, G.P.; Wang, Z.; Salic, A.; Jiang, J.; Huang, P. Literature status in the included PDB file: to be published. |

Suggested wording for documents that use the bundled data:

```text
Structural data source: RCSB PDB / wwPDB, PDB ID <ID>,
https://doi.org/10.2210/pdb<ID>/pdb. PDB archive data files are available
under CC0 1.0.
```

## Runtime And Documentation Package References

Molfig delegates Typst-side mesh rendering to `maquette` through Typst package imports. The maquette package is not vendored in this repository.

The public manual source imports `mantys` and its dependencies for documentation layout. Those documentation packages are not vendored in this repository.
