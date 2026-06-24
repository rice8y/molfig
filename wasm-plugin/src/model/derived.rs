use std::collections::BTreeMap;
use std::sync::OnceLock;

use crate::chemistry::atomic_number;

use super::{AtomicChain, AtomicHierarchy, AtomicResidue, ChemicalComponent};

const MOLSTAR_ION_NAMES: &str = include_str!("reference_data/ions.ts");
const MOLSTAR_LIPID_NAMES: &str = include_str!("reference_data/lipids.ts");
const MOLSTAR_SACCHARIDE_NAMES: &str = include_str!("reference_data/saccharides.ts");
const MOLSTAR_CARBOHYDRATE_CONSTANTS: &str =
    include_str!("reference_data/carbohydrate_constants.ts");

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MoleculeType {
    #[default]
    Unknown,
    Other,
    Water,
    Ion,
    Lipid,
    Protein,
    Rna,
    Dna,
    Pna,
    Saccharide,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PolymerType {
    #[default]
    None,
    PeptideL,
    GammaPeptide,
    BetaPeptide,
    Rna,
    Dna,
    Pna,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SaccharideShape {
    #[default]
    FilledSphere,
    FilledCube,
    CrossedCube,
    DividedDiamond,
    FilledCone,
    DevidedCone,
    FlatBox,
    FilledStar,
    FilledDiamond,
    FlatDiamond,
    FlatHexagon,
    Pentagon,
    DiamondPrism,
    PentagonalPrism,
    HexagonalPrism,
    HeptagonalPrism,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SaccharideType {
    #[default]
    Hexose,
    HexNAc,
    Hexosamine,
    Hexuronate,
    Deoxyhexose,
    DeoxyhexNAc,
    DiDeoxyhexose,
    Pentose,
    Deoxynonulosonate,
    DiDeoxynonulosonate,
    Unknown,
    Assigned,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SaccharideCompIdMapType {
    #[default]
    Default,
    Glycam,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SaccharideComponent {
    pub abbr: String,
    pub name: String,
    pub color: u32,
    pub component_type: SaccharideType,
}

impl Default for SaccharideComponent {
    fn default() -> Self {
        unknown_saccharide_component()
    }
}

pub fn get_saccharide_name(saccharide_type: SaccharideType) -> &'static str {
    match saccharide_type {
        SaccharideType::Hexose => "Hexose",
        SaccharideType::HexNAc => "HexNAc",
        SaccharideType::Hexosamine => "Hexosamine",
        SaccharideType::Hexuronate => "Hexuronate",
        SaccharideType::Deoxyhexose => "Deoxyhexose",
        SaccharideType::DeoxyhexNAc => "DeoxyhexNAc",
        SaccharideType::DiDeoxyhexose => "Di-deoxyhexose",
        SaccharideType::Pentose => "Pentose",
        SaccharideType::Deoxynonulosonate => "Deoxynonulosonate",
        SaccharideType::DiDeoxynonulosonate => "Di-deoxynonulosonate",
        SaccharideType::Unknown => "Unknown",
        SaccharideType::Assigned => "Assigned",
    }
}

pub fn get_saccharide_shape(
    saccharide_type: SaccharideType,
    ring_member_count: usize,
) -> SaccharideShape {
    if saccharide_type == SaccharideType::Unknown {
        match ring_member_count {
            4 => SaccharideShape::DiamondPrism,
            5 => SaccharideShape::PentagonalPrism,
            6 => SaccharideShape::HexagonalPrism,
            7 => SaccharideShape::HeptagonalPrism,
            _ => SaccharideShape::FlatHexagon,
        }
    } else {
        match saccharide_type {
            SaccharideType::Hexose => SaccharideShape::FilledSphere,
            SaccharideType::HexNAc => SaccharideShape::FilledCube,
            SaccharideType::Hexosamine => SaccharideShape::CrossedCube,
            SaccharideType::Hexuronate => SaccharideShape::DividedDiamond,
            SaccharideType::Deoxyhexose => SaccharideShape::FilledCone,
            SaccharideType::DeoxyhexNAc => SaccharideShape::DevidedCone,
            SaccharideType::DiDeoxyhexose => SaccharideShape::FlatBox,
            SaccharideType::Pentose => SaccharideShape::FilledStar,
            SaccharideType::Deoxynonulosonate => SaccharideShape::FilledDiamond,
            SaccharideType::DiDeoxynonulosonate => SaccharideShape::FlatDiamond,
            SaccharideType::Unknown => SaccharideShape::FlatHexagon,
            SaccharideType::Assigned => SaccharideShape::Pentagon,
        }
    }
}

pub fn saccharide_component(comp_id: &str) -> Option<SaccharideComponent> {
    saccharide_component_with_map(comp_id, SaccharideCompIdMapType::Default)
}

pub fn saccharide_component_with_map(
    comp_id: &str,
    map_type: SaccharideCompIdMapType,
) -> Option<SaccharideComponent> {
    let comp_id = comp_id.to_ascii_uppercase();
    if comp_id.is_empty() || comp_id.contains('\'') {
        return None;
    }

    let maps = saccharide_component_maps();
    maps.default
        .get(&comp_id)
        .cloned()
        .or_else(|| {
            (map_type == SaccharideCompIdMapType::Glycam)
                .then(|| maps.glycam.get(&comp_id).cloned())
                .flatten()
        })
        .or_else(|| {
            molstar_generated_name_set_contains(MOLSTAR_SACCHARIDE_NAMES, &comp_id)
                .then(unknown_saccharide_component)
        })
}

struct SaccharideComponentMaps {
    default: BTreeMap<String, SaccharideComponent>,
    glycam: BTreeMap<String, SaccharideComponent>,
}

fn saccharide_component_maps() -> &'static SaccharideComponentMaps {
    static MAPS: OnceLock<SaccharideComponentMaps> = OnceLock::new();
    MAPS.get_or_init(|| {
        let mut default = BTreeMap::new();
        let mut glycam = BTreeMap::new();
        for mono in monosaccharides() {
            for name in saccharide_names_for_map("CommonSaccharideNames", &mono.abbr) {
                default.insert(name, mono.clone());
            }
            for name in saccharide_names_for_map("CharmmSaccharideNames", &mono.abbr) {
                default.insert(name, mono.clone());
            }
            for name in saccharide_names_for_map("GlycamSaccharideNames", &mono.abbr) {
                glycam.entry(name).or_insert_with(|| mono.clone());
            }
        }
        SaccharideComponentMaps { default, glycam }
    })
}

#[derive(Clone, Debug, Default)]
pub struct AtomicDerivedData {
    pub residue: ResidueDerivedData,
    pub atom: AtomDerivedData,
}

impl AtomicDerivedData {
    pub(super) fn from_hierarchy(
        hierarchy: &AtomicHierarchy,
        chemical_components: &[ChemicalComponent],
    ) -> Self {
        let residue = ResidueDerivedData::from_hierarchy(hierarchy, chemical_components);
        let atom = AtomDerivedData::from_hierarchy(hierarchy, &residue);
        AtomicDerivedData { residue, atom }
    }

    pub fn molecule_type_count(&self, molecule_type: MoleculeType) -> usize {
        self.residue
            .molecule_type
            .iter()
            .filter(|ty| **ty == molecule_type)
            .count()
    }

    pub fn trace_element_count(&self) -> usize {
        self.residue
            .trace_element_index
            .iter()
            .filter(|index| index.is_some())
            .count()
    }
}

#[derive(Clone, Debug, Default)]
pub struct ResidueDerivedData {
    pub component_type: Vec<String>,
    pub molecule_type: Vec<MoleculeType>,
    pub polymer_type: Vec<PolymerType>,
    pub is_non_standard: Vec<bool>,
    pub trace_element_index: Vec<Option<usize>>,
    pub direction_from_element_index: Vec<Option<usize>>,
    pub direction_to_element_index: Vec<Option<usize>>,
}

impl ResidueDerivedData {
    pub(super) fn from_hierarchy(
        hierarchy: &AtomicHierarchy,
        chemical_components: &[ChemicalComponent],
    ) -> Self {
        let mut molecule_type = Vec::with_capacity(hierarchy.residues.len());
        let mut component_type = Vec::with_capacity(hierarchy.residues.len());
        let mut polymer_type = Vec::with_capacity(hierarchy.residues.len());
        let is_non_standard = vec![false; hierarchy.residues.len()];
        let mut trace_element_index = Vec::with_capacity(hierarchy.residues.len());
        let mut direction_from_element_index = Vec::with_capacity(hierarchy.residues.len());
        let mut direction_to_element_index = Vec::with_capacity(hierarchy.residues.len());

        for residue in &hierarchy.residues {
            let comp_type = chemical_component_type(chemical_components, &residue.comp_id);
            let mol_type = residue_molecule_type(residue, &comp_type);
            let poly_type = residue_polymer_type(&comp_type, mol_type);
            component_type.push(comp_type);
            molecule_type.push(mol_type);
            polymer_type.push(poly_type);
            trace_element_index.push(residue_trace_atom(hierarchy, residue, poly_type, mol_type));
            let (from, to) = residue_direction_atoms(hierarchy, residue, poly_type);
            direction_from_element_index.push(from);
            direction_to_element_index.push(to);
        }

        ResidueDerivedData {
            component_type,
            molecule_type,
            polymer_type,
            is_non_standard,
            trace_element_index,
            direction_from_element_index,
            direction_to_element_index,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct AtomDerivedData {
    pub atomic_number: Vec<u8>,
    pub is_protein: Vec<bool>,
    pub is_nucleotide: Vec<bool>,
    pub is_water: Vec<bool>,
}

impl AtomDerivedData {
    fn from_hierarchy(hierarchy: &AtomicHierarchy, residue: &ResidueDerivedData) -> Self {
        let mut is_protein = Vec::with_capacity(hierarchy.atoms.len());
        let mut is_nucleotide = Vec::with_capacity(hierarchy.atoms.len());
        let mut is_water = Vec::with_capacity(hierarchy.atoms.len());
        let mut atomic_numbers = Vec::with_capacity(hierarchy.atoms.len());
        for atom in &hierarchy.atoms {
            let molecule_type = residue
                .molecule_type
                .get(atom.residue_index)
                .copied()
                .unwrap_or_default();
            atomic_numbers.push(atomic_number(&atom.element));
            is_protein.push(molecule_type == MoleculeType::Protein);
            is_nucleotide.push(matches!(
                molecule_type,
                MoleculeType::Rna | MoleculeType::Dna | MoleculeType::Pna
            ));
            is_water.push(molecule_type == MoleculeType::Water);
        }
        AtomDerivedData {
            atomic_number: atomic_numbers,
            is_protein,
            is_nucleotide,
            is_water,
        }
    }
}

fn residue_molecule_type(residue: &AtomicResidue, comp_type: &str) -> MoleculeType {
    let comp_id = residue.comp_id.to_ascii_uppercase();
    if is_pna_residue(&comp_id) {
        return MoleculeType::Pna;
    }
    if is_protein_component_type(comp_type) {
        return MoleculeType::Protein;
    }
    if is_rna_component_type(comp_type) {
        return MoleculeType::Rna;
    }
    if is_dna_component_type(comp_type) {
        return MoleculeType::Dna;
    }
    if is_saccharide_component_type(comp_type) {
        return MoleculeType::Saccharide;
    }
    if is_water_residue(&comp_id) {
        return MoleculeType::Water;
    }
    if is_ion_residue(&comp_id) {
        return MoleculeType::Ion;
    }
    if is_lipid_residue(&comp_id) {
        return MoleculeType::Lipid;
    }
    if is_other_component_type(comp_type) {
        if is_saccharide_residue(&comp_id) {
            return MoleculeType::Saccharide;
        }
        if is_protein_residue(&comp_id) {
            return MoleculeType::Protein;
        }
        if is_rna_residue(&comp_id) {
            return MoleculeType::Rna;
        }
        if is_dna_residue(&comp_id) {
            return MoleculeType::Dna;
        }
        return MoleculeType::Other;
    }
    MoleculeType::Unknown
}

fn residue_polymer_type(comp_type: &str, molecule_type: MoleculeType) -> PolymerType {
    match molecule_type {
        MoleculeType::Protein if is_gamma_peptide_component_type(comp_type) => {
            PolymerType::GammaPeptide
        }
        MoleculeType::Protein if is_beta_peptide_component_type(comp_type) => {
            PolymerType::BetaPeptide
        }
        MoleculeType::Protein if is_peptide_terminus_component_type(comp_type) => PolymerType::None,
        MoleculeType::Protein => PolymerType::PeptideL,
        MoleculeType::Rna => PolymerType::Rna,
        MoleculeType::Dna => PolymerType::Dna,
        MoleculeType::Pna => PolymerType::Pna,
        _ => PolymerType::None,
    }
}

fn residue_trace_atom(
    hierarchy: &AtomicHierarchy,
    residue: &AtomicResidue,
    polymer_type: PolymerType,
    molecule_type: MoleculeType,
) -> Option<usize> {
    let trace = match polymer_type {
        PolymerType::PeptideL | PolymerType::GammaPeptide | PolymerType::BetaPeptide => {
            residue_atom_by_names(hierarchy, residue, &["CA"])
        }
        PolymerType::Rna | PolymerType::Dna => {
            residue_atom_by_names(hierarchy, residue, &["O3'", "O3*"])
        }
        PolymerType::Pna => residue_atom_by_names(hierarchy, residue, &["N4'", "N4*"]),
        PolymerType::None => None,
    };
    trace
        .or_else(|| match polymer_type {
            PolymerType::PeptideL => {
                residue_atom_by_names(hierarchy, residue, &["CA", "CA1", "BB", "BAS"])
            }
            PolymerType::GammaPeptide | PolymerType::BetaPeptide => {
                residue_atom_by_names(hierarchy, residue, &["CA"])
            }
            PolymerType::Rna | PolymerType::Dna | PolymerType::Pna => {
                residue_atom_by_names(hierarchy, residue, &["P"])
            }
            PolymerType::None => None,
        })
        .or_else(|| {
            matches!(
                molecule_type,
                MoleculeType::Protein | MoleculeType::Rna | MoleculeType::Dna | MoleculeType::Pna
            )
            .then(|| residue_atom_by_element(hierarchy, residue, "C"))
            .flatten()
        })
}

fn residue_direction_atoms(
    hierarchy: &AtomicHierarchy,
    residue: &AtomicResidue,
    polymer_type: PolymerType,
) -> (Option<usize>, Option<usize>) {
    match polymer_type {
        PolymerType::PeptideL => (
            residue_atom_by_names(hierarchy, residue, &["C"]),
            residue_atom_by_names(hierarchy, residue, &["O", "OC1", "O1", "OX1", "OXT", "OT1"]),
        ),
        PolymerType::GammaPeptide => (
            residue_atom_by_names(hierarchy, residue, &["C"]),
            residue_atom_by_names(hierarchy, residue, &["O"]),
        ),
        PolymerType::BetaPeptide => (
            residue_atom_by_names(hierarchy, residue, &["C"]),
            residue_atom_by_names(hierarchy, residue, &["O"]),
        ),
        PolymerType::Rna => (
            residue_atom_by_names(hierarchy, residue, &["C4'", "C4*"]),
            residue_atom_by_names(hierarchy, residue, &["C3'", "C3*"]),
        ),
        PolymerType::Dna => (
            residue_atom_by_names(hierarchy, residue, &["C3'", "C3*"]),
            residue_atom_by_names(hierarchy, residue, &["C1'", "C1*"]),
        ),
        PolymerType::Pna => (
            residue_atom_by_names(hierarchy, residue, &["N4'", "N4*"]),
            residue_atom_by_names(hierarchy, residue, &["C7'", "C7*"]),
        ),
        PolymerType::None => (None, None),
    }
}

fn residue_atom_by_names(
    hierarchy: &AtomicHierarchy,
    residue: &AtomicResidue,
    names: &[&str],
) -> Option<usize> {
    (residue.start_atom..residue.end_atom).find(|atom_index| {
        hierarchy
            .atoms
            .get(*atom_index)
            .is_some_and(|atom| names.iter().any(|name| atom.name == *name))
    })
}

fn residue_atom_by_element(
    hierarchy: &AtomicHierarchy,
    residue: &AtomicResidue,
    element: &str,
) -> Option<usize> {
    (residue.start_atom..residue.end_atom).find(|atom_index| {
        hierarchy
            .atoms
            .get(*atom_index)
            .is_some_and(|atom| atom.element == element)
    })
}

pub(super) fn chemical_component_type(components: &[ChemicalComponent], comp_id: &str) -> String {
    components
        .iter()
        .rev()
        .find(|component| component.id == comp_id)
        .map(|component| component.type_name.to_ascii_lowercase())
        .unwrap_or_else(|| default_component_type(comp_id))
}

pub(super) fn entity_subtype_from_component(comp_id: &str, comp_type: &str) -> String {
    let comp_id = comp_id.to_ascii_uppercase();
    let comp_type = comp_type.to_ascii_lowercase();
    if matches!(
        comp_type.as_str(),
        "l-peptide linking"
            | "l-peptide nh3 amino terminus"
            | "l-peptide cooh carboxy terminus"
            | "l-gamma-peptide, c-delta linking"
            | "l-beta-peptide, c-gamma linking"
    ) {
        "polypeptide(L)".to_string()
    } else if matches!(
        comp_type.as_str(),
        "d-peptide linking"
            | "d-peptide nh3 amino terminus"
            | "d-peptide cooh carboxy terminus"
            | "d-gamma-peptide, c-delta linking"
            | "d-beta-peptide, c-gamma linking"
    ) {
        "polypeptide(D)".to_string()
    } else if is_rna_component_type(&comp_type) {
        "polyribonucleotide".to_string()
    } else if is_dna_component_type(&comp_type) {
        "polydeoxyribonucleotide".to_string()
    } else if is_saccharide_component_type(&comp_type) || is_saccharide_residue(&comp_id) {
        "oligosaccharide".to_string()
    } else if is_pna_residue(&comp_id) {
        "peptide nucleic acid".to_string()
    } else if is_protein_residue(&comp_id) {
        "polypeptide(L)".to_string()
    } else if is_rna_residue(&comp_id) {
        "polyribonucleotide".to_string()
    } else if is_dna_residue(&comp_id) {
        "polydeoxyribonucleotide".to_string()
    } else if comp_type == "ion" || is_ion_residue(&comp_id) {
        "ion".to_string()
    } else if comp_type == "lipid" || is_lipid_residue(&comp_id) {
        "lipid".to_string()
    } else if matches!(comp_type.as_str(), "peptide linking" | "peptide-like") {
        "peptide-like".to_string()
    } else {
        "other".to_string()
    }
}

pub(crate) fn is_polymer_name(comp_id: &str) -> bool {
    let comp_id = comp_id.to_ascii_uppercase();
    is_protein_residue(&comp_id)
        || is_rna_residue(&comp_id)
        || is_dna_residue(&comp_id)
        || is_pna_residue(&comp_id)
}

pub(crate) fn entity_type_from_component(comp_id: &str) -> &'static str {
    let comp_id = comp_id.to_ascii_uppercase();
    if is_water_residue(&comp_id) {
        "water"
    } else if is_polymer_name(&comp_id) {
        "polymer"
    } else if is_saccharide_residue(&comp_id) {
        "branched"
    } else {
        "non-polymer"
    }
}

pub(crate) fn is_common_protein_cap(comp_id: &str) -> bool {
    matches!(
        comp_id.to_ascii_uppercase().as_str(),
        "NME" | "ACE" | "NH2" | "FOR" | "FMT"
    )
}

pub(crate) fn is_non_polymer_residue_component_type(comp_type: &str) -> bool {
    let comp_type = comp_type.to_ascii_lowercase();
    comp_type.contains("non-polymer")
        || comp_type.contains("amino terminus")
        || comp_type.contains("carboxy terminus")
        || comp_type.contains("peptide-like")
}

pub(crate) fn is_saccharide_component_type_name(comp_type: &str) -> bool {
    is_saccharide_component_type(&comp_type.to_ascii_lowercase())
}

pub(super) fn component_is_non_standard(components: &[ChemicalComponent], comp_id: &str) -> bool {
    components
        .iter()
        .rev()
        .find(|component| component.id == comp_id)
        .map_or_else(
            || !is_polymer_name(comp_id),
            |component| component.mon_nstd_flag.starts_with('n'),
        )
}

fn default_component_type(comp_id: &str) -> String {
    let comp_id = comp_id.to_ascii_uppercase();
    if is_protein_residue(&comp_id) {
        "peptide linking".to_string()
    } else if is_rna_residue(&comp_id) {
        "rna linking".to_string()
    } else if is_dna_residue(&comp_id) {
        "dna linking".to_string()
    } else if is_saccharide_residue(&comp_id) {
        "saccharide".to_string()
    } else {
        "other".to_string()
    }
}

fn is_protein_component_type(comp_type: &str) -> bool {
    matches!(
        comp_type,
        "d-peptide linking"
            | "d-peptide nh3 amino terminus"
            | "d-peptide cooh carboxy terminus"
            | "d-gamma-peptide, c-delta linking"
            | "d-beta-peptide, c-gamma linking"
            | "l-peptide linking"
            | "l-peptide nh3 amino terminus"
            | "l-peptide cooh carboxy terminus"
            | "l-gamma-peptide, c-delta linking"
            | "l-beta-peptide, c-gamma linking"
            | "peptide linking"
            | "peptide-like"
    )
}

fn is_gamma_peptide_component_type(comp_type: &str) -> bool {
    matches!(
        comp_type,
        "d-gamma-peptide, c-delta linking" | "l-gamma-peptide, c-delta linking"
    )
}

fn is_beta_peptide_component_type(comp_type: &str) -> bool {
    matches!(
        comp_type,
        "d-beta-peptide, c-gamma linking" | "l-beta-peptide, c-gamma linking"
    )
}

fn is_peptide_terminus_component_type(comp_type: &str) -> bool {
    matches!(
        comp_type,
        "d-peptide nh3 amino terminus"
            | "d-peptide cooh carboxy terminus"
            | "l-peptide nh3 amino terminus"
            | "l-peptide cooh carboxy terminus"
    )
}

fn is_rna_component_type(comp_type: &str) -> bool {
    matches!(
        comp_type,
        "rna linking" | "l-rna linking" | "rna oh 5 prime terminus" | "rna oh 3 prime terminus"
    )
}

fn is_dna_component_type(comp_type: &str) -> bool {
    matches!(
        comp_type,
        "dna linking" | "l-dna linking" | "dna oh 5 prime terminus" | "dna oh 3 prime terminus"
    )
}

fn is_saccharide_component_type(comp_type: &str) -> bool {
    matches!(
        comp_type,
        "d-saccharide, beta linking"
            | "l-saccharide, beta linking"
            | "d-saccharide, alpha linking"
            | "l-saccharide, alpha linking"
            | "l-saccharide"
            | "d-saccharide"
            | "saccharide"
            | "d-saccharide 1,4 and 1,4 linking"
            | "l-saccharide 1,4 and 1,4 linking"
            | "d-saccharide 1,4 and 1,6 linking"
            | "l-saccharide 1,4 and 1,6 linking"
    )
}

fn is_other_component_type(comp_type: &str) -> bool {
    matches!(comp_type, "non-polymer" | "other")
}

pub(super) fn is_protein_residue(comp_id: &str) -> bool {
    matches!(
        comp_id,
        "ALA"
            | "ARG"
            | "ASN"
            | "ASP"
            | "CYS"
            | "GLN"
            | "GLU"
            | "GLY"
            | "HIS"
            | "ILE"
            | "LEU"
            | "LYS"
            | "MET"
            | "PHE"
            | "PRO"
            | "SER"
            | "THR"
            | "TRP"
            | "TYR"
            | "VAL"
            | "SEC"
            | "PYL"
            | "UNK"
            | "MSE"
            | "SEP"
            | "TPO"
            | "PTR"
            | "PCA"
            | "HYP"
            | "HSD"
            | "HSE"
            | "HSP"
            | "LSN"
            | "ASPP"
            | "GLUP"
            | "HID"
            | "HIE"
            | "HIP"
            | "LYN"
            | "ASH"
            | "GLH"
            | "DAL"
            | "DAR"
            | "DSG"
            | "DAS"
            | "DCY"
            | "DGL"
            | "DGN"
            | "DHI"
            | "DIL"
            | "DLE"
            | "DLY"
            | "MED"
            | "DPN"
            | "DPR"
            | "DSN"
            | "DTH"
            | "DTR"
            | "DTY"
            | "DVA"
            | "DNE"
    )
}

pub(super) fn is_rna_residue(comp_id: &str) -> bool {
    matches!(comp_id, "A" | "C" | "T" | "G" | "I" | "U" | "N")
}

fn is_pna_residue(comp_id: &str) -> bool {
    matches!(comp_id, "APN" | "CPN" | "TPN" | "GPN")
}

fn is_water_residue(comp_id: &str) -> bool {
    matches!(
        comp_id,
        "SOL" | "WAT" | "HOH" | "H2O" | "W" | "DOD" | "D3O" | "TIP" | "TIP3" | "TIP4" | "SPC"
    )
}

fn is_ion_residue(comp_id: &str) -> bool {
    molstar_generated_name_set_contains(MOLSTAR_ION_NAMES, comp_id)
}

fn is_lipid_residue(comp_id: &str) -> bool {
    molstar_generated_name_set_contains(MOLSTAR_LIPID_NAMES, comp_id)
}

pub(super) fn is_dna_residue(comp_id: &str) -> bool {
    matches!(comp_id, "DA" | "DC" | "DT" | "DG" | "DI" | "DU" | "DN")
}

fn is_saccharide_residue(comp_id: &str) -> bool {
    saccharide_component(comp_id).is_some()
}

fn molstar_generated_name_set_contains(source: &str, comp_id: &str) -> bool {
    if comp_id.is_empty() || comp_id.contains('\'') {
        return false;
    }
    let mut quoted = String::with_capacity(comp_id.len() + 2);
    quoted.push('\'');
    quoted.push_str(comp_id);
    quoted.push('\'');
    source.contains(&quoted)
}

fn unknown_saccharide_component() -> SaccharideComponent {
    SaccharideComponent {
        abbr: "Unk".to_string(),
        name: "Unknown".to_string(),
        color: saccharide_color("Secondary"),
        component_type: SaccharideType::Unknown,
    }
}

fn monosaccharides() -> &'static [SaccharideComponent] {
    static MONOSACCHARIDES: OnceLock<Vec<SaccharideComponent>> = OnceLock::new();
    MONOSACCHARIDES
        .get_or_init(parse_monosaccharides)
        .as_slice()
}

fn parse_monosaccharides() -> Vec<SaccharideComponent> {
    let Some(section) = source_section(MOLSTAR_CARBOHYDRATE_CONSTANTS, "const Monosaccharides")
    else {
        return Vec::new();
    };
    section
        .lines()
        .filter_map(|line| {
            let abbr = quoted_field(line, "abbr")?;
            let name = quoted_field(line, "name")?;
            let color = saccharide_color(symbol_field(line, "color", "SaccharideColors.")?);
            let component_type = saccharide_type(symbol_field(line, "type", "SaccharideType.")?)?;
            Some(SaccharideComponent {
                abbr,
                name,
                color,
                component_type,
            })
        })
        .collect()
}

fn saccharide_names_for_map(map_name: &str, abbr: &str) -> Vec<String> {
    let Some(section) =
        source_section(MOLSTAR_CARBOHYDRATE_CONSTANTS, &format!("const {map_name}"))
    else {
        return Vec::new();
    };
    let mut in_entry = false;
    let mut entry = String::new();
    for line in section.lines() {
        let trimmed = line.trim_start();
        if !in_entry && !line_has_object_key(trimmed, abbr) {
            continue;
        }
        in_entry = true;
        entry.push_str(line);
        entry.push('\n');
        if trimmed.contains("],") || trimmed.ends_with(']') {
            return quoted_list_values(&entry);
        }
    }
    if in_entry {
        quoted_list_values(&entry)
    } else {
        Vec::new()
    }
}

fn quoted_list_values(entry: &str) -> Vec<String> {
    let list = entry.split_once('[').map(|(_, list)| list).unwrap_or(entry);
    let mut values = Vec::new();
    let mut rest = list;
    while let Some((_, after_open)) = rest.split_once('\'') {
        let Some((value, after_close)) = after_open.split_once('\'') else {
            break;
        };
        values.push(value.to_string());
        rest = after_close;
    }
    values
}

fn source_section<'a>(source: &'a str, marker: &str) -> Option<&'a str> {
    let start = source.find(marker)?;
    let after_start = &source[start..];
    let end = [after_start.find("\n};"), after_start.find("\n];")]
        .into_iter()
        .flatten()
        .min()?;
    Some(&after_start[..end])
}

fn line_has_object_key(line: &str, key: &str) -> bool {
    line.strip_prefix(key)
        .is_some_and(|rest| rest.trim_start().starts_with(':'))
        || line
            .strip_prefix('\'')
            .and_then(|rest| rest.strip_prefix(key))
            .and_then(|rest| rest.strip_prefix('\''))
            .is_some_and(|rest| rest.trim_start().starts_with(':'))
}

fn quoted_field(line: &str, field: &str) -> Option<String> {
    let pattern = format!("{field}: '");
    let after = line.split_once(&pattern)?.1;
    Some(after.split_once('\'')?.0.to_string())
}

fn symbol_field<'a>(line: &'a str, field: &str, prefix: &str) -> Option<&'a str> {
    let pattern = format!("{field}: {prefix}");
    let after = line.split_once(&pattern)?.1;
    after
        .split(|ch: char| !(ch.is_ascii_alphanumeric() || ch == '_'))
        .next()
}

fn saccharide_color(name: &str) -> u32 {
    match name {
        "Blue" => 0x0090bc,
        "Green" => 0x00a651,
        "Yellow" => 0xffd400,
        "Orange" => 0xf47920,
        "Pink" => 0xf69ea1,
        "Purple" => 0xa54399,
        "LightBlue" => 0x8fcce9,
        "Brown" => 0xa17a4d,
        "Red" => 0xed1c24,
        "Secondary" => 0xf1ece1,
        _ => 0xf1ece1,
    }
}

fn saccharide_type(name: &str) -> Option<SaccharideType> {
    Some(match name {
        "Hexose" => SaccharideType::Hexose,
        "HexNAc" => SaccharideType::HexNAc,
        "Hexosamine" => SaccharideType::Hexosamine,
        "Hexuronate" => SaccharideType::Hexuronate,
        "Deoxyhexose" => SaccharideType::Deoxyhexose,
        "DeoxyhexNAc" => SaccharideType::DeoxyhexNAc,
        "DiDeoxyhexose" => SaccharideType::DiDeoxyhexose,
        "Pentose" => SaccharideType::Pentose,
        "Deoxynonulosonate" => SaccharideType::Deoxynonulosonate,
        "DiDeoxynonulosonate" => SaccharideType::DiDeoxynonulosonate,
        "Unknown" => SaccharideType::Unknown,
        "Assigned" => SaccharideType::Assigned,
        _ => return None,
    })
}

pub(super) fn is_polymer_residue(hierarchy: &AtomicHierarchy, residue_index: usize) -> bool {
    let molecule_type = hierarchy
        .derived
        .residue
        .molecule_type
        .get(residue_index)
        .copied()
        .unwrap_or_default();
    matches!(
        molecule_type,
        MoleculeType::Protein | MoleculeType::Rna | MoleculeType::Dna | MoleculeType::Pna
    ) && hierarchy
        .derived
        .residue
        .trace_element_index
        .get(residue_index)
        .is_some_and(|index| index.is_some())
}

pub(super) fn chain_residue_indices(
    hierarchy: &AtomicHierarchy,
    chain_index: usize,
    chain: &AtomicChain,
) -> Vec<usize> {
    (chain.start_residue..chain.end_residue)
        .filter(|residue_index| {
            hierarchy
                .residues
                .get(*residue_index)
                .is_some_and(|residue| residue.chain_index == chain_index)
        })
        .collect()
}

pub(super) fn are_backbone_connected(
    hierarchy: &AtomicHierarchy,
    start: usize,
    end: usize,
) -> bool {
    let polymer_start = hierarchy
        .derived
        .residue
        .polymer_type
        .get(start)
        .copied()
        .unwrap_or_default();
    let polymer_end = hierarchy
        .derived
        .residue
        .polymer_type
        .get(end)
        .copied()
        .unwrap_or_default();
    if polymer_start == PolymerType::None || polymer_end == PolymerType::None {
        return false;
    }
    if hierarchy
        .derived
        .residue
        .trace_element_index
        .get(start)
        .and_then(|index| *index)
        .is_none()
        || hierarchy
            .derived
            .residue
            .trace_element_index
            .get(end)
            .and_then(|index| *index)
            .is_none()
    {
        return false;
    }

    let mut start_element = residue_backbone_start_atom(hierarchy, start, polymer_start);
    let mut end_element = residue_backbone_end_atom(hierarchy, end, polymer_end);
    if start_element.is_none() || end_element.is_none() {
        start_element = residue_coarse_backbone_atom(hierarchy, start, polymer_start);
        end_element = residue_coarse_backbone_atom(hierarchy, end, polymer_end);
    }
    let (Some(start_element), Some(end_element)) = (start_element, end_element) else {
        return false;
    };
    let is_coarse = hierarchy
        .derived
        .residue
        .direction_from_element_index
        .get(start)
        .and_then(|index| *index)
        .is_none()
        || hierarchy
            .derived
            .residue
            .direction_to_element_index
            .get(start)
            .and_then(|index| *index)
            .is_none()
        || hierarchy
            .derived
            .residue
            .direction_from_element_index
            .get(end)
            .and_then(|index| *index)
            .is_none()
        || hierarchy
            .derived
            .residue
            .direction_to_element_index
            .get(end)
            .and_then(|index| *index)
            .is_none();
    let Some(start_atom) = hierarchy.atoms.get(start_element) else {
        return false;
    };
    let Some(end_atom) = hierarchy.atoms.get(end_element) else {
        return false;
    };
    start_atom.position.distance(end_atom.position) < if is_coarse { 10.0 } else { 3.0 }
}

fn residue_backbone_start_atom(
    hierarchy: &AtomicHierarchy,
    residue_index: usize,
    polymer_type: PolymerType,
) -> Option<usize> {
    let residue = hierarchy.residues.get(residue_index)?;
    match polymer_type {
        PolymerType::PeptideL => residue_atom_by_names(hierarchy, residue, &["N"]),
        PolymerType::GammaPeptide => residue_atom_by_names(hierarchy, residue, &["N"]),
        PolymerType::BetaPeptide => residue_atom_by_names(hierarchy, residue, &["N"]),
        PolymerType::Rna | PolymerType::Dna => residue_atom_by_names(hierarchy, residue, &["P"]),
        PolymerType::Pna => residue_atom_by_names(hierarchy, residue, &["N1'", "N1*"]),
        PolymerType::None => None,
    }
}

fn residue_backbone_end_atom(
    hierarchy: &AtomicHierarchy,
    residue_index: usize,
    polymer_type: PolymerType,
) -> Option<usize> {
    let residue = hierarchy.residues.get(residue_index)?;
    match polymer_type {
        PolymerType::PeptideL => residue_atom_by_names(hierarchy, residue, &["C"]),
        PolymerType::GammaPeptide => residue_atom_by_names(hierarchy, residue, &["CD"]),
        PolymerType::BetaPeptide => residue_atom_by_names(hierarchy, residue, &["CG"]),
        PolymerType::Rna | PolymerType::Dna => {
            residue_atom_by_names(hierarchy, residue, &["O3'", "O3*"])
        }
        PolymerType::Pna => residue_atom_by_names(hierarchy, residue, &["C'", "C*"]),
        PolymerType::None => None,
    }
}

fn residue_coarse_backbone_atom(
    hierarchy: &AtomicHierarchy,
    residue_index: usize,
    polymer_type: PolymerType,
) -> Option<usize> {
    let residue = hierarchy.residues.get(residue_index)?;
    match polymer_type {
        PolymerType::PeptideL => {
            residue_atom_by_names(hierarchy, residue, &["CA", "CA1", "BB", "BAS"])
        }
        PolymerType::GammaPeptide | PolymerType::BetaPeptide => {
            residue_atom_by_names(hierarchy, residue, &["CA"])
        }
        PolymerType::Rna | PolymerType::Dna | PolymerType::Pna => {
            residue_atom_by_names(hierarchy, residue, &["P"])
        }
        PolymerType::None => None,
    }
}

pub(super) fn residue_seq_id(residue: &AtomicResidue) -> Option<i32> {
    residue.label_seq_id.trim().parse::<i32>().ok()
}
