# Example Structure Data

The files in this directory are example molecular structure data used by Molfig examples, documentation, and regression checks.

## RCSB PDB / wwPDB examples

Source: RCSB PDB / wwPDB PDB archive

License/dedication: CC0 1.0 Universal Public Domain Dedication

RCSB PDB policies: https://www.rcsb.org/pages/policies

CC0 1.0: https://creativecommons.org/publicdomain/zero/1.0/

RCSB PDB states that data files in the PDB archive are available under CC0 1.0.
RCSB PDB also encourages attribution to the original structure-data authors where possible. No endorsement by the authors, RCSB PDB, wwPDB, or Creative Commons is implied.

| File | PDB ID | Format | PDB DOI | Structure authors / status |
| --- | --- | --- | --- | --- |
| `1crn.bcif` | 1CRN | BinaryCIF | https://doi.org/10.2210/pdb1CRN/pdb | Hendrickson, W.A.; Teeter, M.M. Primary citation: Teeter, M.M. (1984) Proc Natl Acad Sci U S A 81:6014-6018. Article DOI: https://doi.org/10.1073/pnas.81.19.6014 |
| `1FYY.cif` | 1FYY | mmCIF | https://doi.org/10.2210/pdb1FYY/pdb | Volk, D.E.; Rice, J.S.; Luxon, B.A.; Yeh, H.J.C.; Liang, C.; Xie, G.; Sayer, J.M.; Jerina, D.M.; Gorenstein, D.G. Primary citation: Biochemistry 39:14040-14053 (2000). Article DOI: https://doi.org/10.1021/bi001669l |
| `9M1U.pdb` | 9M1U | PDB | https://doi.org/10.2210/pdb9M1U/pdb | Liu, H.; Zhang, X.; Xu, H.E. Primary citation: Zhang, X. et al. (2026), EMBO J. Article DOI: https://doi.org/10.1038/s44318-026-00823-y |
| `9q12.pdb` | 9Q12 | PDB | https://doi.org/10.2210/pdb9Q12/pdb | Wang, Y.; Liu, B.; He, Y.; Feigon, J. Literature status in the included PDB file: to be published. |
| `9R1O.pdb` | 9R1O | PDB | https://doi.org/10.2210/pdb9R1O/pdb | Petrenas, R.; Ozga, K.; Chubb, J.J.; Woolfson, D.N. Literature status in the included PDB file: to be published. |
| `9Z4O.pdb` | 9Z4O | PDB | https://doi.org/10.2210/pdb9Z4O/pdb | Ge, Y.; de Almeida Magalhaes, T.; Wu, H.; Yadav, G.P.; Wang, Z.; Salic, A.; Jiang, J.; Huang, P. Literature status in the included PDB file: to be published. |

## PubChem XYZ validation corpus

The XYZ corpus contains PubChem3D conformer coordinates retrieved as 3D SDF
through PubChem PUG REST on 2026-08-24. Each file preserves the source atom
order and four-decimal coordinates.

| File | Compound | PubChem CID | Formula | Conformer ID | Atoms | Source SDF bonds |
| --- | --- | ---: | --- | --- | ---: | ---: |
| `ethanol.xyz` | ethanol | 702 | C2H6O | `000002BE00000001` | 9 | 8 |
| `benzene.xyz` | benzene | 241 | C6H6 | `000000F100000001` | 12 | 12 |
| `aspirin.xyz` | aspirin | 2244 | C9H8O4 | `000008C400000001` | 21 | 21 |
| `caffeine.xyz` | caffeine | 2519 | C8H10N4O2 | `000009D700000001` | 24 | 25 |

`XYZ_VALIDATION.json` pins each source URL, conformer ID, molecular formula,
atom and source-bond counts, exact source bond endpoints and orders, element
composition, coordinate bounds, and SHA-256 digest. The offline validator at
`wasm-plugin/tests/validate-pubchem-xyz.mjs` additionally checks canonical XYZ
syntax, finite and non-coincident coordinates, and exact agreement between
Mol*-style XYZ bond inference and the source SDF connectivity.

- Compound records: https://pubchem.ncbi.nlm.nih.gov/compound/702, https://pubchem.ncbi.nlm.nih.gov/compound/241, https://pubchem.ncbi.nlm.nih.gov/compound/2244, https://pubchem.ncbi.nlm.nih.gov/compound/2519
- Retrieval service: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- General PubChem citation: Kim S, Chen J, Cheng T, et al. PubChem 2025 update. Nucleic Acids Res. 2025;53(D1):D1516-D1525. https://doi.org/10.1093/nar/gkae1059
- PubChem data submission policy: https://pubchem.ncbi.nlm.nih.gov/docs/data-submission-policy
- NCBI molecular data usage policy: https://www.ncbi.nlm.nih.gov/home/about/policies/

The PubChem data submission policy states that PubChem-generated information is
made available without cost and without restriction. NCBI also notes that some
submitters may claim rights in contributed molecular data. This corpus contains
only PubChem-generated conformer coordinates and records its provenance here.

Suggested wording:

```text
Structural data source: RCSB PDB / wwPDB, PDB ID <ID>,
https://doi.org/10.2210/pdb<ID>/pdb. PDB archive data files are available
under CC0 1.0.
```

For an XYZ corpus record:

```text
Coordinate source: PubChem CID <CID> (<compound>), PubChem3D conformer
<conformer ID>, retrieved through PubChem PUG REST on 2026-08-24.
https://pubchem.ncbi.nlm.nih.gov/compound/<CID>
https://doi.org/10.1093/nar/gkae1059
```
