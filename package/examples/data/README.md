# Example Structure Data

The files in this directory are example molecular structure data used by Molfig examples, documentation, and regression checks.

## Published theme-gallery data

`../theme.typ` and the documentation use published structures and annotations, not project-authored molecular specimens. Prediction scores and calculated partial charges are real computational results, not experimental measurements. The additional source files below were obtained on 2026-09-07; gallery builds work offline. The existing `1crn.bcif` supplies the biological-assembly symmetry no-op controls.

### Oxyhemoglobin: chain, entity, operator, element, and carbon themes

`1HHO.cif` is an unchanged download of [RCSB PDB entry 1HHO](https://doi.org/10.2210/pdb1HHO/pdb), human oxyhemoglobin, from [the PDB mmCIF archive](https://files.rcsb.org/download/1HHO.cif). Biological assembly 1 contains identity and twofold-related copies of the alpha/beta pair, plus the deposited nonpolymer components. Deposition author: Shaanan, B.; primary citation: Shaanan (1983), *Journal of Molecular Biology* 171, 31-59. PDB data: [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).

### AlphaFold: pLDDT confidence

`AF-P01308-F1-model_v6.cif` is the unchanged [AlphaFold DB model for human insulin precursor P01308](https://alphafold.ebi.ac.uk/entry/P01308), downloaded from [the versioned model URL](https://alphafold.ebi.ac.uk/files/AF-P01308-F1-model_v6.cif). Its 110 residues and 839 atoms carry the published per-residue pLDDT values in `_ma_qa_metric_local`; the example has three confidence bands and no score above 90. No scores were edited to fill a missing band. Source: Google DeepMind / EMBL-EBI, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), with the [AlphaFold DB license and disclaimer](https://alphafold.ebi.ac.uk/assets/License-Disclaimer.pdf). Citation: Jumper et al. (2021), *Nature* 596, 583-589, [doi:10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2).

### SWISS-MODEL: QMEAN local quality

`QMEAN.json` and `qmean-model_001_processed.pdb` are unchanged downloads from the [official QMEAN example](https://swissmodel.expasy.org/qmean/QMEAN): [JSON results](https://swissmodel.expasy.org/qmean/QMEAN.json) and [processed model_001 PDB](https://swissmodel.expasy.org/qmean/QMEAN/model_001_processed.pdb). The source result was created on 2019-03-04. `qmean-model_001.cif` is a format conversion of this published model and its full-precision local quality scores, not a newly constructed model or an independently calculated QMEAN result.

Run `node package/examples/data/convert-qmean.mjs` from the repository root to reproduce the CIF, or append `--check` to verify it without writing. The converter preserves all 2,306 atom coordinates, atom order, and author residue identifiers. It maps JSON `models.model_001.scores.local_scores.A[seq_id - 1]` to `_ma_qa_metric_local` for the 295 observed residues; the six unmodeled N-terminal residues stay unmodeled. It checks residue names against the published sequence and checks every local score against the rounded B-factor field of the official processed PDB. The [QMEAN documentation](https://swissmodel.expasy.org/qmean/help) explains this processed-PDB numbering and B-factor convention. No score is inferred from a crystallographic B-factor, rounded to the PDB precision, or rescaled.

Source and attribution: SWISS-MODEL, Computational Structural Biology Group, SIB Swiss Institute of Bioinformatics, Biozentrum, University of Basel. The source data, converted CIF, and adapted QMEAN renderings are distributed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), as specified in the [SWISS-MODEL terms](https://swissmodel.expasy.org/docs/terms_of_use); these data are not relicensed under Molfig's software license. Changes: PDB-to-CIF conversion with the existing local QMEAN annotations, and molecular rendering. QMEAN citation: Benkert, Biasini, and Schwede (2011), *Bioinformatics* 27, 343-350, [doi:10.1093/bioinformatics/btq662](https://doi.org/10.1093/bioinformatics/btq662). SWISS-MODEL citation: Waterhouse et al. (2018), *Nucleic Acids Research* 46(W1), W296-W303, [doi:10.1093/nar/gky427](https://doi.org/10.1093/nar/gky427). No endorsement is implied.

### Mol* / SB-NCBR: computed partial charges

`7qpd.fw2.cif` is copied unchanged from the [Mol* example at commit `1b8117d3f10f7c978aabb5a0d3d47370635aefe4`](https://github.com/molstar/molstar/blob/1b8117d3f10f7c978aabb5a0d3d47370635aefe4/examples/7qpd.fw2.cif), also present in this checkout's Mol* reference repository. The underlying experimental structure is [PDB 7QPD](https://doi.org/10.2210/pdb7QPD/pdb), a human MHC class I complex. The file contains five computed charge sets in `_sb_ncbr_partial_atomic_charges`; Molfig selects the first set, `SQE+qp/Schindler 2021 (CCD_gen)`, and the theme sums atomic charges per residue. These are the source's calculated charges, not hand-assigned values. Underlying PDB coordinates: CC0 1.0. Mol* example distribution: MIT, copyright Mol* contributors; the original license is retained in `7qpd-LICENSE.txt`.

### SHA-256 pins

| File | SHA-256 |
| --- | --- |
| `1HHO.cif` | `4d7b3ff09044c5d0a55f7adbfcfac71b0e33cbedd6c63896ecac5d865a05c434` |
| `AF-P01308-F1-model_v6.cif` | `0d71181346ce6d218d681417e849d3cc6f7b04dfd11fa85702b2984da2b81b14` |
| `QMEAN.json` | `4b77017e858a5dd284ff90e52f057281403a7083e83f061e2c58a10b4e641e1d` |
| `qmean-model_001_processed.pdb` | `ca974cfa3a579c955040f3b3a549512ab9cb632668994c09eec4780f1292ec96` |
| `qmean-model_001.cif` | `a13822595cb6458d55f2113ea39a2f7444eb9c99b3f169a5d23bb85dccb59fce` |
| `7qpd.fw2.cif` | `96df5e5515fddf3db5e390211a45784a016e3efd25be5adaaeed6afbd99501e9` |

## RCSB PDB / wwPDB examples

Source: RCSB PDB / wwPDB PDB archive

License/dedication: CC0 1.0 Universal Public Domain Dedication

RCSB PDB policies: https://www.rcsb.org/pages/policies

CC0 1.0: https://creativecommons.org/publicdomain/zero/1.0/

RCSB PDB states that data files in the PDB archive are available under CC0 1.0. RCSB PDB also encourages attribution to the original structure-data authors where possible. No endorsement by the authors, RCSB PDB, wwPDB, or Creative Commons is implied.

| File | PDB ID | Format | PDB DOI | Structure authors / status |
| --- | --- | --- | --- | --- |
| `1crn.bcif` | 1CRN | BinaryCIF | https://doi.org/10.2210/pdb1CRN/pdb | Hendrickson, W.A.; Teeter, M.M. Primary citation: Teeter, M.M. (1984) Proc Natl Acad Sci U S A 81:6014-6018. Article DOI: https://doi.org/10.1073/pnas.81.19.6014 |
| `1FYY.cif` | 1FYY | mmCIF | https://doi.org/10.2210/pdb1FYY/pdb | Volk, D.E.; Rice, J.S.; Luxon, B.A.; Yeh, H.J.C.; Liang, C.; Xie, G.; Sayer, J.M.; Jerina, D.M.; Gorenstein, D.G. Primary citation: Biochemistry 39:14040-14053 (2000). Article DOI: https://doi.org/10.1021/bi001669l |
| `9M1U.pdb` | 9M1U | PDB | https://doi.org/10.2210/pdb9M1U/pdb | Liu, H.; Zhang, X.; Xu, H.E. Primary citation: Zhang, X. et al. (2026), EMBO J. Article DOI: https://doi.org/10.1038/s44318-026-00823-y |
| `9q12.pdb` | 9Q12 | PDB | https://doi.org/10.2210/pdb9Q12/pdb | Wang, Y.; Liu, B.; He, Y.; Feigon, J. Literature status in the included PDB file: to be published. |
| `9R1O.pdb` | 9R1O | PDB | https://doi.org/10.2210/pdb9R1O/pdb | Petrenas, R.; Ozga, K.; Chubb, J.J.; Woolfson, D.N. Literature status in the included PDB file: to be published. |
| `9Z4O.pdb` | 9Z4O | PDB | https://doi.org/10.2210/pdb9Z4O/pdb | Ge, Y.; de Almeida Magalhaes, T.; Wu, H.; Yadav, G.P.; Wang, Z.; Salic, A.; Jiang, J.; Huang, P. Literature status in the included PDB file: to be published. |

## PubChem XYZ validation corpus

The XYZ corpus contains PubChem3D conformer coordinates retrieved as 3D SDF through PubChem PUG REST on 2026-08-24. Each file preserves the source atom order and four-decimal coordinates.

| File | Compound | PubChem CID | Formula | Conformer ID | Atoms | Source SDF bonds |
| --- | --- | ---: | --- | --- | ---: | ---: |
| `ethanol.xyz` | ethanol | 702 | C2H6O | `000002BE00000001` | 9 | 8 |
| `benzene.xyz` | benzene | 241 | C6H6 | `000000F100000001` | 12 | 12 |
| `aspirin.xyz` | aspirin | 2244 | C9H8O4 | `000008C400000001` | 21 | 21 |
| `caffeine.xyz` | caffeine | 2519 | C8H10N4O2 | `000009D700000001` | 24 | 25 |

`XYZ_VALIDATION.json` pins each source URL, conformer ID, molecular formula, atom and source-bond counts, exact source bond endpoints and orders, element composition, coordinate bounds, and SHA-256 digest. The offline validator at `wasm-plugin/tests/validate-pubchem-xyz.mjs` additionally checks canonical XYZ syntax, finite and non-coincident coordinates, and exact agreement between Mol*-style XYZ bond inference and the source SDF connectivity.

- Compound records: https://pubchem.ncbi.nlm.nih.gov/compound/702, https://pubchem.ncbi.nlm.nih.gov/compound/241, https://pubchem.ncbi.nlm.nih.gov/compound/2244, https://pubchem.ncbi.nlm.nih.gov/compound/2519
- Retrieval service: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- General PubChem citation: Kim S, Chen J, Cheng T, et al. PubChem 2025 update. Nucleic Acids Res. 2025;53(D1):D1516-D1525. https://doi.org/10.1093/nar/gkae1059
- PubChem data submission policy: https://pubchem.ncbi.nlm.nih.gov/docs/data-submission-policy
- NCBI molecular data usage policy: https://www.ncbi.nlm.nih.gov/home/about/policies/

The PubChem data submission policy states that PubChem-generated information is made available without cost and without restriction. NCBI also notes that some submitters may claim rights in contributed molecular data. This corpus contains only PubChem-generated conformer coordinates and records its provenance here.

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
