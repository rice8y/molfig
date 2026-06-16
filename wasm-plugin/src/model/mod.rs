use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::sync::atomic::{AtomicU64, Ordering};

#[derive(Clone, Debug)]
pub struct Atom {
    pub id: usize,
    pub source_index: usize,
    pub model_num: i32,
    pub name: String,
    pub type_symbol: String,
    pub element: String,
    pub chain: String,
    pub auth_chain: String,
    pub entity_id: String,
    pub residue: String,
    pub auth_residue: String,
    pub group_pdb: String,
    pub residue_seq: String,
    pub auth_residue_seq: String,
    pub insertion_code: String,
    pub alt_id: String,
    pub auth_name: String,
    pub occupancy: f32,
    pub b_iso: f32,
    pub formal_charge: i32,
    pub position: Vec3,
    pub het: bool,
    pub operator_name: String,
}

#[derive(Clone, Copy, Debug)]
pub struct Bond {
    pub a: usize,
    pub b: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MolstarBondSiteEntry {
    pub atom_id_1: usize,
    pub atom_id_2: usize,
    pub value_order: Option<&'static str>,
    pub type_id: Option<&'static str>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IndexPairBonds {
    pub bonds: IndexPairGraph,
    pub max_distance: f32,
    pub cacheable: bool,
    pub has_operators: bool,
    pub by_same_operator: BTreeMap<i32, Vec<usize>>,
    source_bond_by_slot: Vec<usize>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct IntraUnitBonds {
    pub vertex_count: usize,
    pub offset: Vec<usize>,
    pub a: Vec<usize>,
    pub b: Vec<usize>,
    pub edge_count: usize,
    pub props: IntraUnitBondProps,
    pub can_remap: bool,
    pub cacheable: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct IntraUnitBondProps {
    pub key: Vec<i32>,
    pub order: Vec<i8>,
    pub flags: Vec<BondFlags>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct IndexPairGraph {
    pub vertex_count: usize,
    pub offset: Vec<usize>,
    pub a: Vec<usize>,
    pub b: Vec<usize>,
    pub edge_count: usize,
    pub props: IndexPairProps,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct IndexPairProps {
    pub key: Vec<i32>,
    pub operator_a: Vec<i32>,
    pub operator_b: Vec<i32>,
    pub order: Vec<i8>,
    pub distance: Vec<f32>,
    pub flag: Vec<BondFlags>,
}

impl IntraUnitBonds {
    fn from_edges(
        atom_a: &[usize],
        atom_b: &[usize],
        metadata: &[BondMetadata],
        atom_count: usize,
        can_remap: bool,
        cacheable: bool,
    ) -> Self {
        let edge_count = atom_a.len();
        let mut bucket_sizes = vec![0usize; atom_count];
        for (&a, &b) in atom_a.iter().zip(atom_b) {
            if a < atom_count {
                bucket_sizes[a] += 1;
            }
            if b < atom_count {
                bucket_sizes[b] += 1;
            }
        }
        let mut offset = vec![0usize; atom_count + 1];
        let mut cursor = 0usize;
        for (index, size) in bucket_sizes.iter().enumerate() {
            offset[index] = cursor;
            cursor += *size;
        }
        offset[atom_count] = cursor;

        let mut graph = IntraUnitBonds {
            vertex_count: atom_count,
            offset,
            a: vec![0; cursor],
            b: vec![0; cursor],
            edge_count,
            props: IntraUnitBondProps {
                key: vec![0; cursor],
                order: vec![0; cursor],
                flags: vec![BondFlags::NONE; cursor],
            },
            can_remap,
            cacheable,
        };
        let mut bucket_fill = vec![0usize; atom_count];
        for edge_index in 0..edge_count {
            let x = atom_a[edge_index];
            let y = atom_b[edge_index];
            if x >= atom_count || y >= atom_count {
                continue;
            }
            let slot_xy = graph.offset[x] + bucket_fill[x];
            bucket_fill[x] += 1;
            let slot_yx = graph.offset[y] + bucket_fill[y];
            bucket_fill[y] += 1;
            let metadata = metadata.get(edge_index).cloned().unwrap_or_default();

            graph.a[slot_xy] = x;
            graph.b[slot_xy] = y;
            graph.assign_slot(slot_xy, &metadata);

            graph.a[slot_yx] = y;
            graph.b[slot_yx] = x;
            graph.assign_slot(slot_yx, &metadata);
        }
        graph
    }

    fn assign_slot(&mut self, slot: usize, metadata: &BondMetadata) {
        self.props.key[slot] = metadata.key;
        self.props.order[slot] = metadata.order;
        self.props.flags[slot] = metadata.flags;
    }
}

impl IndexPairBonds {
    pub fn from_bonds(
        bonds: &[Bond],
        metadata: &[BondMetadata],
        atom_count: usize,
        max_distance: f32,
        cacheable: bool,
    ) -> Option<Self> {
        let mut index_a = Vec::new();
        let mut index_b = Vec::new();
        let mut source_bonds = Vec::new();
        let mut pair_metadata = Vec::new();
        for (source_bond, metadata) in metadata.iter().enumerate() {
            if metadata.source != BondSource::IndexPair {
                continue;
            }
            let Some(bond) = bonds.get(source_bond) else {
                continue;
            };
            index_a.push(bond.a);
            index_b.push(bond.b);
            source_bonds.push(source_bond);
            pair_metadata.push(metadata.clone());
        }
        Self::from_pairs(
            &index_a,
            &index_b,
            &source_bonds,
            &pair_metadata,
            atom_count,
            max_distance,
            cacheable,
        )
    }

    pub fn from_pairs(
        index_a: &[usize],
        index_b: &[usize],
        source_bonds: &[usize],
        metadata: &[BondMetadata],
        atom_count: usize,
        max_distance: f32,
        cacheable: bool,
    ) -> Option<Self> {
        if index_a.is_empty() || index_a.len() != index_b.len() || index_a.len() != metadata.len() {
            return None;
        }
        if index_a
            .iter()
            .chain(index_b.iter())
            .any(|&index| index >= atom_count)
        {
            return None;
        }
        let source_bonds = if source_bonds.len() == index_a.len() {
            source_bonds.to_vec()
        } else {
            (0..index_a.len()).collect()
        };
        let (bonds, source_bond_by_slot) =
            IndexPairGraph::from_edges(index_a, index_b, &source_bonds, metadata, atom_count);
        let has_operators = metadata
            .iter()
            .all(|metadata| metadata.operator_a >= 0 && metadata.operator_b >= 0);
        let by_same_operator = if has_operators {
            bonds.by_same_operator()
        } else {
            BTreeMap::new()
        };
        Some(IndexPairBonds {
            bonds,
            max_distance,
            cacheable,
            has_operators,
            by_same_operator,
            source_bond_by_slot,
        })
    }

    pub fn contains_bond(&self, bond_index: usize) -> bool {
        self.source_bond_by_slot.contains(&bond_index)
    }

    pub fn has_operators(&self, metadata: &[BondMetadata]) -> bool {
        let _ = metadata;
        self.has_operators
    }

    pub fn by_same_operator(&self, metadata: &[BondMetadata]) -> BTreeMap<i32, Vec<usize>> {
        let _ = metadata;
        self.by_same_operator.clone()
    }

    pub fn get_edge_index_for_operators(
        &self,
        i: usize,
        j: usize,
        op_i: i32,
        op_j: i32,
    ) -> Option<usize> {
        self.bonds.get_edge_index_for_operators(i, j, op_i, op_j)
    }
}

impl IndexPairGraph {
    fn from_edges(
        index_a: &[usize],
        index_b: &[usize],
        source_bonds: &[usize],
        metadata: &[BondMetadata],
        atom_count: usize,
    ) -> (Self, Vec<usize>) {
        let edge_count = index_a.len();
        let mut bucket_sizes = vec![0usize; atom_count];
        for (&a, &b) in index_a.iter().zip(index_b) {
            if a < atom_count {
                bucket_sizes[a] += 1;
            }
            if b < atom_count {
                bucket_sizes[b] += 1;
            }
        }
        let mut offset = vec![0usize; atom_count + 1];
        let mut cursor = 0usize;
        for (index, size) in bucket_sizes.iter().enumerate() {
            offset[index] = cursor;
            cursor += *size;
        }
        offset[atom_count] = cursor;

        let mut graph = IndexPairGraph {
            vertex_count: atom_count,
            offset,
            a: vec![0; cursor],
            b: vec![0; cursor],
            edge_count,
            props: IndexPairProps {
                key: vec![-1; cursor],
                operator_a: vec![-1; cursor],
                operator_b: vec![-1; cursor],
                order: vec![1; cursor],
                distance: vec![-1.0; cursor],
                flag: vec![BondFlags::COVALENT; cursor],
            },
        };
        let mut source_bond_by_slot = vec![0; cursor];
        let mut bucket_fill = vec![0usize; atom_count];
        for edge_index in 0..edge_count {
            let x = index_a[edge_index];
            let y = index_b[edge_index];
            if x >= atom_count || y >= atom_count {
                continue;
            }
            let slot_xy = graph.offset[x] + bucket_fill[x];
            bucket_fill[x] += 1;
            let slot_yx = graph.offset[y] + bucket_fill[y];
            bucket_fill[y] += 1;
            let metadata = metadata.get(edge_index).cloned().unwrap_or_default();
            let source_bond = source_bonds.get(edge_index).copied().unwrap_or(edge_index);

            graph.a[slot_xy] = x;
            graph.b[slot_xy] = y;
            graph.assign_slot(slot_xy, &metadata, metadata.operator_a, metadata.operator_b);
            source_bond_by_slot[slot_xy] = source_bond;

            graph.a[slot_yx] = y;
            graph.b[slot_yx] = x;
            graph.assign_slot(slot_yx, &metadata, metadata.operator_b, metadata.operator_a);
            source_bond_by_slot[slot_yx] = source_bond;
        }
        (graph, source_bond_by_slot)
    }

    fn assign_slot(
        &mut self,
        slot: usize,
        metadata: &BondMetadata,
        operator_a: i32,
        operator_b: i32,
    ) {
        self.props.key[slot] = metadata.key;
        self.props.operator_a[slot] = operator_a;
        self.props.operator_b[slot] = operator_b;
        self.props.order[slot] = metadata.order;
        self.props.distance[slot] = metadata.distance.unwrap_or(-1.0);
        self.props.flag[slot] = metadata.flags;
    }

    fn by_same_operator(&self) -> BTreeMap<i32, Vec<usize>> {
        let mut map = BTreeMap::<i32, Vec<usize>>::new();
        for slot in 0..self.props.operator_a.len() {
            if self.props.operator_a[slot] == self.props.operator_b[slot] {
                map.entry(self.props.operator_a[slot])
                    .or_default()
                    .push(slot);
            }
        }
        map
    }

    fn get_edge_index_for_operators(
        &self,
        i: usize,
        j: usize,
        op_i: i32,
        op_j: i32,
    ) -> Option<usize> {
        let (a, b, op_a, op_b) = if i < j {
            (i, j, op_i, op_j)
        } else {
            (j, i, op_j, op_i)
        };
        let start = *self.offset.get(a)?;
        let end = *self.offset.get(a + 1)?;
        (start..end).find(|&slot| {
            self.b[slot] == b
                && self.props.operator_a[slot] == op_a
                && self.props.operator_b[slot] == op_b
        })
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct BondFlags {
    pub bits: u16,
}

impl BondFlags {
    pub const NONE: BondFlags = BondFlags { bits: 0x0 };
    pub const COVALENT: BondFlags = BondFlags { bits: 0x1 };
    pub const METALLIC_COORDINATION: BondFlags = BondFlags { bits: 0x2 };
    pub const HYDROGEN_BOND: BondFlags = BondFlags { bits: 0x4 };
    pub const DISULFIDE: BondFlags = BondFlags { bits: 0x8 };
    pub const AROMATIC: BondFlags = BondFlags { bits: 0x10 };
    pub const COMPUTED: BondFlags = BondFlags { bits: 0x20 };
    pub const RESONANCE: BondFlags = BondFlags { bits: 0x40 };

    pub const fn contains(self, other: BondFlags) -> bool {
        self.bits & other.bits == other.bits
    }

    pub const fn union(self, other: BondFlags) -> BondFlags {
        BondFlags {
            bits: self.bits | other.bits,
        }
    }

    pub const fn without(self, other: BondFlags) -> BondFlags {
        BondFlags {
            bits: self.bits & !other.bits,
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub enum BondSource {
    #[default]
    Computed,
    PdbConect,
    StructConn,
    IndexPair,
    ChemComp,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct StructConnMetadata {
    pub id: String,
    pub row_index: usize,
    pub partner_a_atom_index: usize,
    pub partner_b_atom_index: usize,
    pub conn_type_id: String,
    pub value_order: String,
    pub partner_a_symmetry: String,
    pub partner_b_symmetry: String,
    pub partner_a_comp_id: String,
    pub partner_b_comp_id: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct BondMetadata {
    pub source: BondSource,
    pub order: i8,
    pub flags: BondFlags,
    pub key: i32,
    pub distance: Option<f32>,
    pub operator_a: i32,
    pub operator_b: i32,
    pub struct_conn: Option<StructConnMetadata>,
}

impl Default for BondMetadata {
    fn default() -> Self {
        Self::computed()
    }
}

impl BondMetadata {
    pub fn computed() -> Self {
        BondMetadata {
            source: BondSource::Computed,
            order: 1,
            flags: BondFlags::COVALENT.union(BondFlags::COMPUTED),
            key: -1,
            distance: None,
            operator_a: -1,
            operator_b: -1,
            struct_conn: None,
        }
    }

    pub fn computed_for_atoms(a: &Atom, b: &Atom) -> Self {
        let is_metal = (is_metal_element(&a.element) || is_metal_element(&b.element))
            && !(is_hydrogen_element(&a.element) || is_hydrogen_element(&b.element));
        let flags = if is_metal {
            BondFlags::METALLIC_COORDINATION
        } else {
            BondFlags::COVALENT
        }
        .union(BondFlags::COMPUTED);
        BondMetadata {
            flags,
            ..BondMetadata::computed()
        }
    }

    pub fn pdb_conect(key: i32) -> Self {
        BondMetadata {
            source: BondSource::PdbConect,
            order: 1,
            flags: BondFlags::COVALENT,
            key,
            distance: None,
            operator_a: -1,
            operator_b: -1,
            struct_conn: None,
        }
    }
}

fn is_hydrogen_element(element: &str) -> bool {
    matches!(
        element.trim().to_ascii_uppercase().as_str(),
        "H" | "D" | "T"
    )
}

fn is_metal_element(element: &str) -> bool {
    matches!(
        element.trim().to_ascii_uppercase().as_str(),
        "LI" | "NA"
            | "K"
            | "RB"
            | "CS"
            | "FR"
            | "BE"
            | "MG"
            | "CA"
            | "SR"
            | "BA"
            | "RA"
            | "AL"
            | "GA"
            | "IN"
            | "SN"
            | "TL"
            | "PB"
            | "BI"
            | "SC"
            | "TI"
            | "V"
            | "CR"
            | "MN"
            | "FE"
            | "CO"
            | "NI"
            | "CU"
            | "ZN"
            | "Y"
            | "ZR"
            | "NB"
            | "MO"
            | "TC"
            | "RU"
            | "RH"
            | "PD"
            | "AG"
            | "CD"
            | "LA"
            | "HF"
            | "TA"
            | "W"
            | "RE"
            | "OS"
            | "IR"
            | "PT"
            | "AU"
            | "HG"
            | "AC"
            | "RF"
            | "DB"
            | "SG"
            | "BH"
            | "HS"
            | "MT"
            | "CE"
            | "PR"
            | "ND"
            | "PM"
            | "SM"
            | "EU"
            | "GD"
            | "TB"
            | "DY"
            | "HO"
            | "ER"
            | "TM"
            | "YB"
            | "LU"
            | "TH"
            | "PA"
            | "U"
            | "NP"
            | "PU"
            | "AM"
            | "CM"
            | "BK"
            | "CF"
            | "ES"
            | "FM"
            | "MD"
            | "NO"
            | "LR"
    )
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChemicalComponent {
    pub id: String,
    pub name: String,
    pub type_name: String,
    pub formula: String,
    pub formula_weight: Option<f32>,
    pub one_letter_code: String,
    pub three_letter_code: String,
    pub mon_nstd_flag: String,
    pub pdbx_synonyms: String,
    pub pdbx_formal_charge: Option<i32>,
    pub pdbx_initial_date: String,
    pub pdbx_modified_date: String,
    pub pdbx_ambiguous_flag: String,
    pub pdbx_release_status: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChemicalComponentBond {
    pub comp_id: String,
    pub atom_id_1: String,
    pub atom_id_2: String,
    pub order: i8,
    pub flags: BondFlags,
    pub stereo_config: String,
    pub ordinal: Option<i32>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChemicalComponentAtom {
    pub comp_id: String,
    pub atom_id: String,
    pub alt_atom_id: String,
    pub type_symbol: String,
    pub charge: i32,
    pub aromatic: bool,
    pub leaving_atom: bool,
    pub stereo_config: String,
    pub model_cartn: Option<Vec3>,
    pub ideal_cartn: Option<Vec3>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ChemicalComponentAngle {
    pub comp_id: String,
    pub atom_id_1: String,
    pub atom_id_2: String,
    pub atom_id_3: String,
    pub value_angle: Option<f32>,
    pub value_angle_esd: Option<f32>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Entity {
    pub id: String,
    pub type_name: String,
    pub description: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EntityIndexMap {
    pub ids: Vec<String>,
    pub subtype: Vec<String>,
    pub id_to_index: BTreeMap<String, usize>,
}

impl EntityIndexMap {
    pub fn from_entities(
        entities: &[Entity],
        entity_polymers: &[EntityPoly],
        entity_branches: &[PdbxEntityBranch],
    ) -> Self {
        let mut ids = Vec::with_capacity(entities.len());
        let mut subtype = Vec::with_capacity(entities.len());
        let mut id_to_index = BTreeMap::new();
        for (index, entity) in entities.iter().enumerate() {
            ids.push(entity.id.clone());
            id_to_index.insert(entity.id.clone(), index);
            subtype.push(entity_subtype(entity, entity_polymers, entity_branches));
        }
        EntityIndexMap {
            ids,
            subtype,
            id_to_index,
        }
    }

    pub fn get_entity_index(&self, id: &str) -> Option<usize> {
        self.id_to_index.get(id).copied()
    }
}

fn entity_subtype(
    entity: &Entity,
    entity_polymers: &[EntityPoly],
    entity_branches: &[PdbxEntityBranch],
) -> String {
    entity_polymers
        .iter()
        .find(|polymer| polymer.entity_id == entity.id)
        .map(|polymer| polymer.polymer_type.clone())
        .or_else(|| {
            entity_branches
                .iter()
                .find(|branch| branch.entity_id == entity.id)
                .map(|branch| branch.type_name.clone())
        })
        .filter(|subtype| !subtype.is_empty())
        .unwrap_or_else(|| entity.type_name.clone())
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EntityPoly {
    pub entity_id: String,
    pub polymer_type: String,
    pub sequence: String,
    pub nstd_linkage: String,
    pub nstd_monomer: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EntityPolySeq {
    pub entity_id: String,
    pub num: i32,
    pub mon_id: String,
    pub hetero: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StructureSequence {
    pub sequences: Vec<SequenceEntity>,
    pub by_entity_key: BTreeMap<usize, usize>,
}

impl StructureSequence {
    fn from_model_parts(
        entity_index: &EntityIndexMap,
        entities: &[Entity],
        entity_poly_seq: &[EntityPolySeq],
        hierarchy: &AtomicHierarchy,
        coarse: Option<&CoarseHierarchy>,
    ) -> Self {
        if !entity_poly_seq.is_empty() {
            return StructureSequence::from_entity_poly_seq(entity_index, entity_poly_seq);
        }
        let mut sequence = StructureSequence::from_atomic_hierarchy(entities, hierarchy);
        if let Some(coarse) = coarse.filter(|coarse| coarse.is_defined) {
            sequence.merge_molstar(StructureSequence::from_coarse_elements(
                entity_index,
                &coarse.spheres,
            ));
            sequence.merge_molstar(StructureSequence::from_coarse_elements(
                entity_index,
                &coarse.gaussians,
            ));
        }
        sequence
    }

    fn from_entity_poly_seq(
        entity_index: &EntityIndexMap,
        entity_poly_seq: &[EntityPolySeq],
    ) -> Self {
        let mut sequence = StructureSequence::default();
        let mut row = 0usize;
        while row < entity_poly_seq.len() {
            let start = row;
            let entity_id = &entity_poly_seq[start].entity_id;
            while row + 1 < entity_poly_seq.len()
                && entity_poly_seq[row + 1].entity_id == *entity_id
            {
                row += 1;
            }
            row += 1;
            if let Some(entity_key) = entity_index.get_entity_index(entity_id) {
                let mut residues = Vec::<SequenceResidue>::new();
                let mut index_by_seq_id = BTreeMap::<i32, usize>::new();
                let mut micro_het = BTreeMap::<i32, Vec<String>>::new();
                for entry in &entity_poly_seq[start..row] {
                    if let Some(&residue_index) = index_by_seq_id.get(&entry.num) {
                        let variants = micro_het
                            .entry(entry.num)
                            .or_insert_with(|| vec![residues[residue_index].comp_id.clone()]);
                        if !variants.iter().any(|comp_id| comp_id == &entry.mon_id) {
                            variants.push(entry.mon_id.clone());
                        }
                    } else {
                        index_by_seq_id.insert(entry.num, residues.len());
                        residues.push(SequenceResidue {
                            comp_id: entry.mon_id.clone(),
                            seq_id: entry.num,
                        });
                    }
                }
                sequence.push_entity(
                    entity_key,
                    entity_id.clone(),
                    residues,
                    Vec::new(),
                    micro_het,
                );
            }
        }
        sequence
    }

    fn from_atomic_hierarchy(entities: &[Entity], hierarchy: &AtomicHierarchy) -> Self {
        let mut sequence = StructureSequence::default();
        for (chain_index, chain) in hierarchy.chains.iter().enumerate() {
            let Some(entity_key) = hierarchy
                .index
                .chain_entity_index
                .get(chain_index)
                .copied()
                .flatten()
            else {
                continue;
            };
            if sequence.by_entity_key.contains_key(&entity_key)
                || entities
                    .get(entity_key)
                    .is_none_or(|entity| entity.type_name != "polymer")
            {
                continue;
            }
            let residues = hierarchy.residues[chain.start_residue..chain.end_residue]
                .iter()
                .map(|residue| SequenceResidue {
                    comp_id: residue.comp_id.clone(),
                    seq_id: residue.label_seq_id.parse::<i32>().unwrap_or(0),
                })
                .collect();
            sequence.push_entity(
                entity_key,
                chain.entity_id.clone(),
                residues,
                Vec::new(),
                BTreeMap::new(),
            );
        }
        sequence
    }

    fn from_coarse_elements(
        entity_index: &EntityIndexMap,
        elements: &CoarseElements,
    ) -> StructureSequence {
        let mut sequence = StructureSequence::default();
        for segment in 0..elements.chain_element_segments.count {
            let start = elements.chain_element_segments.offsets[segment];
            let end = elements.chain_element_segments.offsets[segment + 1];
            let Some(first) = elements.elements.get(start) else {
                continue;
            };
            let Some(entity_key) = entity_index.get_entity_index(&first.entity_id) else {
                continue;
            };
            if sequence.by_entity_key.contains_key(&entity_key) {
                continue;
            }
            let ranges = elements.elements[start..end]
                .iter()
                .map(|element| SequenceRange {
                    seq_id_begin: element.seq_id_begin,
                    seq_id_end: element.seq_id_end,
                })
                .collect();
            sequence.push_entity(
                entity_key,
                first.entity_id.clone(),
                Vec::new(),
                ranges,
                BTreeMap::new(),
            );
        }
        sequence
    }

    fn push_entity(
        &mut self,
        entity_key: usize,
        entity_id: String,
        residues: Vec<SequenceResidue>,
        ranges: Vec<SequenceRange>,
        micro_het: BTreeMap<i32, Vec<String>>,
    ) {
        if self.by_entity_key.contains_key(&entity_key) {
            return;
        }
        let index = self.sequences.len();
        let index_by_seq_id = residues
            .iter()
            .enumerate()
            .map(|(index, residue)| (residue.seq_id, index))
            .collect();
        self.sequences.push(SequenceEntity {
            entity_id,
            residues,
            ranges,
            index_by_seq_id,
            micro_het,
        });
        self.by_entity_key.insert(entity_key, index);
    }

    fn merge_molstar(&mut self, other: StructureSequence) {
        let offset = self.sequences.len();
        self.sequences.extend(other.sequences);
        for (entity_key, sequence_index) in other.by_entity_key {
            self.by_entity_key
                .insert(entity_key, sequence_index + offset);
        }
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SequenceEntity {
    pub entity_id: String,
    pub residues: Vec<SequenceResidue>,
    pub ranges: Vec<SequenceRange>,
    pub index_by_seq_id: BTreeMap<i32, usize>,
    pub micro_het: BTreeMap<i32, Vec<String>>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SequenceResidue {
    pub comp_id: String,
    pub seq_id: i32,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SequenceRange {
    pub seq_id_begin: i32,
    pub seq_id_end: i32,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PdbxEntityBranch {
    pub entity_id: String,
    pub type_name: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PdbxEntityBranchLink {
    pub link_id: i32,
    pub details: String,
    pub entity_id: String,
    pub entity_branch_list_num_1: i32,
    pub entity_branch_list_num_2: i32,
    pub comp_id_1: String,
    pub comp_id_2: String,
    pub atom_id_1: String,
    pub leaving_atom_id_1: String,
    pub atom_stereo_config_1: String,
    pub atom_id_2: String,
    pub leaving_atom_id_2: String,
    pub atom_stereo_config_2: String,
    pub value_order: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PdbxBranchScheme {
    pub entity_id: String,
    pub hetero: String,
    pub asym_id: String,
    pub mon_id: String,
    pub num: i32,
    pub pdb_asym_id: String,
    pub pdb_seq_num: String,
    pub pdb_mon_id: String,
    pub auth_asym_id: String,
    pub auth_seq_num: String,
    pub auth_mon_id: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PdbxNonpolyScheme {
    pub asym_id: String,
    pub entity_id: String,
    pub mon_id: String,
    pub pdb_strand_id: String,
    pub ndb_seq_num: String,
    pub pdb_seq_num: String,
    pub auth_seq_num: String,
    pub pdb_mon_id: String,
    pub auth_mon_id: String,
    pub pdb_ins_code: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PdbxPolySeqScheme {
    pub asym_id: String,
    pub entity_id: String,
    pub seq_id: i32,
    pub mon_id: String,
    pub ndb_seq_num: String,
    pub pdb_seq_num: String,
    pub auth_seq_num: String,
    pub pdb_mon_id: String,
    pub auth_mon_id: String,
    pub pdb_strand_id: String,
    pub pdb_ins_code: String,
    pub hetero: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct IhmModelList {
    pub model_id: i32,
    pub model_name: String,
    pub assembly_id: i32,
    pub protocol_id: i32,
    pub representation_id: i32,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IhmModelGroup {
    pub id: i32,
    pub name: String,
    pub details: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IhmModelGroupLink {
    pub model_id: i32,
    pub group_id: i32,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct IhmCrossLinkRestraint {
    pub id: i32,
    pub group_id: i32,
    pub entity_id_1: String,
    pub entity_id_2: String,
    pub asym_id_1: String,
    pub asym_id_2: String,
    pub comp_id_1: String,
    pub comp_id_2: String,
    pub seq_id_1: i32,
    pub seq_id_2: i32,
    pub atom_id_1: String,
    pub atom_id_2: String,
    pub restraint_type: String,
    pub conditional_crosslink_flag: String,
    pub model_granularity: String,
    pub distance_threshold: Option<f32>,
    pub psi: Option<f32>,
    pub sigma_1: Option<f32>,
    pub sigma_2: Option<f32>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StructAsym {
    pub id: String,
    pub entity_id: String,
    pub details: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Entry {
    pub id: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Experiment {
    pub method: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AtomSiteAnisotrop {
    pub atom_id: usize,
    pub u: AnisotropicDisplacement,
}

pub type AnisotropicDisplacement = [[f32; 3]; 3];

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AtomSiteColumnPresence {
    pub occupancy_defined: bool,
    pub b_iso_defined: bool,
    pub xyz_defined: bool,
}

impl Default for AtomSiteColumnPresence {
    fn default() -> Self {
        Self {
            occupancy_defined: true,
            b_iso_defined: true,
            xyz_defined: true,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Molecule {
    pub source_data: SourceData,
    pub atom_site_columns: AtomSiteColumnPresence,
    pub global_model_transform: Option<GlobalModelTransform>,
    pub entries: Vec<Entry>,
    pub experiments: Vec<Experiment>,
    pub atoms: Vec<Atom>,
    pub atom_site_anisotrop: Vec<AtomSiteAnisotrop>,
    pub bonds: Vec<Bond>,
    pub bond_metadata: Vec<BondMetadata>,
    pub index_pair_bonds: Option<IndexPairBonds>,
    pub coarse_spheres: Vec<CoarseSphere>,
    pub coarse_gaussians: Vec<CoarseGaussian>,
    pub assemblies: Vec<Assembly>,
    pub selected_assembly: Option<Assembly>,
    pub helices: Vec<SecondaryRange>,
    pub sheets: Vec<SecondaryRange>,
    pub entities: Vec<Entity>,
    pub entity_index: EntityIndexMap,
    pub entity_polymers: Vec<EntityPoly>,
    pub entity_poly_seq: Vec<EntityPolySeq>,
    pub pdbx_entity_branch: Vec<PdbxEntityBranch>,
    pub pdbx_entity_branch_links: Vec<PdbxEntityBranchLink>,
    pub pdbx_branch_scheme: Vec<PdbxBranchScheme>,
    pub pdbx_nonpoly_scheme: Vec<PdbxNonpolyScheme>,
    pub pdbx_poly_seq_scheme: Vec<PdbxPolySeqScheme>,
    pub ihm_model_list: Vec<IhmModelList>,
    pub ihm_model_groups: Vec<IhmModelGroup>,
    pub ihm_model_group_links: Vec<IhmModelGroupLink>,
    pub ihm_cross_link_restraints: Vec<IhmCrossLinkRestraint>,
    pub struct_asym: Vec<StructAsym>,
    pub chemical_components: Vec<ChemicalComponent>,
    pub chemical_component_atoms: Vec<ChemicalComponentAtom>,
    pub chemical_component_bonds: Vec<ChemicalComponentBond>,
    pub chemical_component_angles: Vec<ChemicalComponentAngle>,
    pub rings: Vec<Ring>,
    pub resonance: Resonance,
    pub(crate) derived_aromatic_bonds: BTreeSet<usize>,
    pub(crate) derived_resonance_bonds: BTreeSet<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct GlobalModelTransform {
    pub matrix: [[f32; 4]; 4],
}

impl GlobalModelTransform {
    pub const DESCRIPTOR: &'static str = "molstar_global_model_transform_info";

    pub fn to_property_string(&self) -> String {
        self.matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect::<Vec<_>>()
            .join(";")
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SourceData {
    pub kind: String,
    pub name: String,
    pub original_kind: String,
    pub categories: Vec<SourceCategory>,
    pub db_categories: Vec<SourceCategory>,
    pub frame_categories: Vec<SourceCategory>,
}

impl SourceData {
    pub fn pdb(name: impl Into<String>) -> Self {
        SourceData {
            kind: "pdb".to_string(),
            name: name.into(),
            original_kind: String::new(),
            categories: Vec::new(),
            db_categories: Vec::new(),
            frame_categories: Vec::new(),
        }
    }

    pub fn mmcif(name: impl Into<String>, categories: Vec<SourceCategory>) -> Self {
        let db_categories = categories
            .iter()
            .filter(|category| is_mmcif_db_category(&category.name))
            .cloned()
            .collect();
        SourceData {
            kind: "mmCIF".to_string(),
            name: name.into(),
            original_kind: String::new(),
            categories: categories.clone(),
            db_categories,
            frame_categories: categories,
        }
    }

    pub fn from_pdb_as_mmcif(name: impl Into<String>, categories: Vec<SourceCategory>) -> Self {
        let db_categories = categories
            .iter()
            .filter(|category| is_mmcif_db_category(&category.name))
            .cloned()
            .collect();
        SourceData {
            kind: "mmCIF".to_string(),
            name: name.into(),
            original_kind: "pdb".to_string(),
            categories: categories.clone(),
            db_categories,
            frame_categories: categories,
        }
    }
}

fn is_mmcif_db_category(name: &str) -> bool {
    matches!(
        name,
        "entry"
            | "exptl"
            | "entity"
            | "entity_poly"
            | "entity_poly_seq"
            | "struct_asym"
            | "atom_site"
            | "atom_site_anisotrop"
            | "chem_comp"
            | "chem_comp_atom"
            | "chem_comp_bond"
            | "chem_comp_angle"
            | "struct_conn"
            | "struct_conf"
            | "struct_sheet_range"
            | "pdbx_struct_assembly"
            | "pdbx_struct_assembly_gen"
            | "pdbx_struct_oper_list"
            | "pdbx_entity_branch"
            | "pdbx_entity_branch_link"
            | "pdbx_branch_scheme"
            | "pdbx_nonpoly_scheme"
            | "pdbx_poly_seq_scheme"
            | "ihm_model_list"
            | "ihm_model_group"
            | "ihm_model_group_link"
            | "ihm_sphere_obj_site"
            | "ihm_gaussian_obj_site"
            | "ihm_cross_link_restraint"
            | "molstar_bond_site"
    )
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SourceCategory {
    pub name: String,
    pub row_count: usize,
    pub column_count: usize,
}

impl Molecule {
    pub fn atomic_structure(&self) -> AtomicStructure {
        AtomicStructure::from_molecule(self)
    }

    pub fn carbohydrates(&self) -> Carbohydrates {
        self.atomic_structure().carbohydrates
    }

    pub(crate) fn expanded_for_geometry(&self) -> Molecule {
        let structure = self.atomic_structure();
        if self.selected_assembly.is_none() {
            return self.clone();
        }

        let mut atoms = Vec::new();
        let mut expanded_indices = Vec::<(usize, usize)>::new();
        for unit in &structure.units {
            if unit.kind != UnitKind::Atomic {
                continue;
            }
            for &atom_index in &unit.elements {
                let Some(source_atom) = self.atoms.get(atom_index) else {
                    continue;
                };
                let mut atom = source_atom.clone();
                atom.position = unit.operator.transform.apply(source_atom.position);
                atom.operator_name = unit.operator.name.clone();
                expanded_indices.push((unit.id, atom_index));
                atoms.push(atom);
            }
        }
        let (bonds, bond_metadata) = expanded_bonds_from_units(
            &expanded_indices,
            &self.bonds,
            &self.bond_metadata,
            &structure,
        );
        let (coarse_spheres, coarse_gaussians) = expanded_coarse_for_geometry(self);

        let mut molecule = self.clone();
        molecule.atoms = atoms;
        molecule.bonds = bonds;
        molecule.bond_metadata = bond_metadata;
        molecule.index_pair_bonds = None;
        molecule.coarse_spheres = coarse_spheres;
        molecule.coarse_gaussians = coarse_gaussians;
        molecule.selected_assembly = None;
        molecule
    }

    pub(crate) fn identity_assembly_subset_for_geometry(&self) -> Option<Molecule> {
        let assembly = self.selected_assembly.as_ref()?;
        let generators = if assembly.generators.is_empty() {
            vec![AssemblyGenerator::from_transforms(
                &assembly.id,
                assembly.asym_ids.clone(),
                0,
                assembly.transforms.clone(),
                vec![Vec::new(); assembly.transforms.len()],
            )]
        } else {
            assembly.generators.clone()
        };
        let mut operator_offset = 0usize;
        let mut asym_ids = BTreeSet::new();
        for generator in generators {
            let operators = generator.operators_for_assembly(&assembly.id, operator_offset);
            operator_offset += operators.len();
            if operators.len() != 1 || !operators[0].transform.is_identity() {
                return None;
            }
            asym_ids.extend(generator.asym_ids);
        }
        if asym_ids.is_empty() {
            return None;
        }

        let mut atom_index_map = vec![None; self.atoms.len()];
        let mut atoms = Vec::new();
        for (source_index, atom) in self.atoms.iter().enumerate() {
            if asym_ids.contains(&atom.chain) || asym_ids.contains(&atom.auth_chain) {
                atom_index_map[source_index] = Some(atoms.len());
                atoms.push(atom.clone());
            }
        }
        if atoms.len() == self.atoms.len() {
            let mut molecule = self.clone();
            molecule.selected_assembly = None;
            return Some(molecule);
        }

        let mut bonds = Vec::new();
        let mut bond_metadata = Vec::new();
        for (bond, metadata) in self.bonds.iter().zip(&self.bond_metadata) {
            let (Some(a), Some(b)) = (
                atom_index_map.get(bond.a).copied().flatten(),
                atom_index_map.get(bond.b).copied().flatten(),
            ) else {
                continue;
            };
            bonds.push(Bond { a, b });
            bond_metadata.push(metadata.clone());
        }

        let mut molecule = self.clone();
        molecule.atoms = atoms;
        molecule.bonds = bonds;
        molecule.bond_metadata = bond_metadata;
        molecule.index_pair_bonds = None;
        molecule.coarse_spheres = self
            .coarse_spheres
            .iter()
            .filter(|sphere| asym_ids.contains(&sphere.asym_id))
            .cloned()
            .collect();
        molecule.coarse_gaussians = self
            .coarse_gaussians
            .iter()
            .filter(|gaussian| asym_ids.contains(&gaussian.asym_id))
            .cloned()
            .collect();
        molecule.helices = self
            .helices
            .iter()
            .filter(|range| asym_ids.contains(&range.chain))
            .cloned()
            .collect();
        molecule.sheets = self
            .sheets
            .iter()
            .filter(|range| asym_ids.contains(&range.chain))
            .cloned()
            .collect();
        molecule.selected_assembly = None;
        Some(molecule)
    }

    pub(crate) fn identity_assembly_trace_subset_for_geometry(&self) -> Option<Molecule> {
        let mut molecule = self.identity_assembly_subset_for_geometry()?;
        molecule.bonds.clear();
        molecule.bond_metadata.clear();
        molecule.index_pair_bonds = None;
        Some(molecule)
    }

    pub(crate) fn refresh_topology_metadata(&mut self) {
        for bond_index in std::mem::take(&mut self.derived_aromatic_bonds) {
            if let Some(metadata) = self.bond_metadata.get_mut(bond_index) {
                metadata.flags = metadata.flags.without(BondFlags::AROMATIC);
            }
        }
        for bond_index in std::mem::take(&mut self.derived_resonance_bonds) {
            if let Some(metadata) = self.bond_metadata.get_mut(bond_index) {
                metadata.flags = metadata.flags.without(BondFlags::RESONANCE);
            }
        }
        apply_chemical_component_bonds(self);
        if self.bond_metadata.len() < self.bonds.len() {
            self.bond_metadata
                .extend((self.bond_metadata.len()..self.bonds.len()).map(|index| {
                    let Some(bond) = self.bonds.get(index) else {
                        return BondMetadata::computed();
                    };
                    self.atoms
                        .get(bond.a)
                        .zip(self.atoms.get(bond.b))
                        .map(|(a, b)| BondMetadata::computed_for_atoms(a, b))
                        .unwrap_or_else(BondMetadata::computed)
                }));
        }
        assign_intra_bond_orders(self);
        self.rings = detect_rings(self);
        assign_ring_resonance(self);
        self.resonance = build_resonance(self);
    }

    pub fn molstar_bond_site_entries(&self) -> Vec<MolstarBondSiteEntry> {
        let structure = self.atomic_structure();
        let mut entries = Vec::new();
        let mut added = HashSet::<(usize, usize)>::new();

        let mut add = |atom_a: usize, atom_b: usize, order: i8, flags: BondFlags| {
            let Some(a) = self.atoms.get(atom_a).map(|atom| atom.id) else {
                return;
            };
            let Some(b) = self.atoms.get(atom_b).map(|atom| atom.id) else {
                return;
            };
            let (atom_id_1, atom_id_2) = if a > b { (b, a) } else { (a, b) };
            if !added.insert((atom_id_1, atom_id_2)) {
                return;
            }
            entries.push(MolstarBondSiteEntry {
                atom_id_1,
                atom_id_2,
                value_order: molstar_bond_site_export_value_order(order, flags),
                type_id: molstar_bond_site_export_type_id(flags),
            });
        };

        for unit in &structure.units {
            if unit.kind != UnitKind::Atomic {
                continue;
            }
            for (source_bond, bond) in self.bonds.iter().enumerate() {
                let Some(index_a) = unit
                    .atom_indices
                    .iter()
                    .position(|&source| source == bond.a)
                else {
                    continue;
                };
                let Some(index_b) = unit
                    .atom_indices
                    .iter()
                    .position(|&source| source == bond.b)
                else {
                    continue;
                };
                let metadata = self
                    .bond_metadata
                    .get(source_bond)
                    .cloned()
                    .unwrap_or_default();
                if metadata_allows_intra_bond(&metadata, unit)
                    && metadata_allows_bond_distance(
                        &metadata,
                        source_bond,
                        self.index_pair_bonds.as_ref(),
                        unit,
                        index_a,
                        unit,
                        index_b,
                    )
                {
                    add(bond.a, bond.b, metadata.order, metadata.flags);
                }
            }
        }

        for bond in &structure.inter_unit_bonds {
            let Some(unit_a) = structure.unit_by_id(bond.unit_a) else {
                continue;
            };
            let Some(unit_b) = structure.unit_by_id(bond.unit_b) else {
                continue;
            };
            let Some(&atom_a) = unit_a.atom_indices.get(bond.index_a) else {
                continue;
            };
            let Some(&atom_b) = unit_b.atom_indices.get(bond.index_b) else {
                continue;
            };
            add(atom_a, atom_b, bond.order, bond.flags);
        }

        entries.sort_by(|a, b| {
            a.atom_id_1
                .cmp(&b.atom_id_1)
                .then_with(|| a.atom_id_2.cmp(&b.atom_id_2))
        });
        entries
    }
}

fn molstar_bond_site_export_value_order(order: i8, flags: BondFlags) -> Option<&'static str> {
    let mut value_order = match order {
        1 => Some("sing"),
        2 => Some("doub"),
        3 => Some("trip"),
        4 => Some("quad"),
        _ => None,
    };
    if flags.contains(BondFlags::AROMATIC) {
        value_order = Some("arom");
    }
    value_order
}

fn molstar_bond_site_export_type_id(flags: BondFlags) -> Option<&'static str> {
    if flags.contains(BondFlags::DISULFIDE) {
        Some("disulf")
    } else if flags.contains(BondFlags::COVALENT) {
        Some("covale")
    } else if flags.contains(BondFlags::METALLIC_COORDINATION) {
        Some("metalc")
    } else if flags.contains(BondFlags::HYDROGEN_BOND) {
        Some("hydrog")
    } else {
        None
    }
}

mod coarse;

pub use coarse::{
    CoarseConformation, CoarseElement, CoarseElementKey, CoarseElementKind, CoarseElementReference,
    CoarseElements, CoarseGaussian, CoarseGaussianConformation, CoarseHierarchy, CoarseIndex,
    CoarseModel, CoarseRange, CoarseSegmentation, CoarseSphere, CoarseSphereConformation,
};

mod branched;
pub use branched::{BranchedEntityLinkMap, BranchedEntityLinkPlacement, BranchedSequenceMap};

mod carbohydrates;
pub use carbohydrates::{
    CarbohydrateElement, CarbohydrateLink, CarbohydrateSymbolGeometry, CarbohydrateTerminalLink,
    Carbohydrates, PartialCarbohydrateElement,
};

mod expansion;
mod topology;

use expansion::{expanded_bonds_from_units, expanded_coarse_for_geometry};
#[allow(unused_imports)]
pub(crate) use topology::intra_bond_order_from_table;
use topology::{
    apply_chemical_component_bonds, assign_intra_bond_orders, assign_ring_resonance,
    build_resonance, detect_rings,
};
pub use topology::{DelocalizedTriplets, Resonance, Ring};

mod assembly;
mod geometry_data;
mod math;

pub use assembly::{Assembly, AssemblyGenerator, AssemblyOperator, SecondaryRange};
pub(crate) use geometry_data::{
    Face, Mesh, MeshMaterial, MeshSection, NucleotideAtoms, NucleotideBaseKind, TraceResidue,
};
pub use math::{Axes3D, PrincipalAxes, Transform, Vec3};

#[derive(Clone, Debug, Default)]
pub struct AtomicStructure {
    pub model: AtomicModel,
    pub models: Vec<AtomicModel>,
    pub coarse: CoarseModel,
    pub coarse_models: Vec<CoarseModel>,
    pub units: Vec<StructureUnit>,
    pub coordinate_system: UnitOperator,
    pub boundary: Boundary,
    pub principal_axes: PrincipalAxes,
    pub lookup3d: StructureLookup3D,
    pub ranges: AtomicRanges,
    pub properties: StructureProperties,
    pub element_count: usize,
    pub polymer_residue_count: usize,
    pub polymer_gap_count: usize,
    pub symmetry_groups: Vec<UnitSymmetryGroup>,
    pub intra_unit_bond_count: usize,
    pub inter_unit_bonds: Vec<InterUnitBond>,
    pub inter_unit_bond_graph: InterUnitBonds,
    pub carbohydrates: Carbohydrates,
}

impl AtomicStructure {
    fn from_molecule(molecule: &Molecule) -> Self {
        let mut models = atomic_models(&molecule.atoms, molecule);
        let mut coarse_models = coarse_models(molecule);
        if coarse_models.is_empty() {
            coarse_models.push(CoarseModel::default());
        }
        let empty_coarse = CoarseModel::default();
        for (model_index, model) in models.iter_mut().enumerate() {
            let model_coarse = coarse_models.get(model_index).unwrap_or(&empty_coarse);
            model.sequence = StructureSequence::from_model_parts(
                &molecule.entity_index,
                &molecule.entities,
                &molecule.entity_poly_seq,
                &model.hierarchy,
                Some(&model_coarse.hierarchy),
            );
        }
        let coarse = coarse_models.first().cloned().unwrap_or_default();
        let mut model = models
            .first()
            .cloned()
            .unwrap_or_else(|| AtomicModel::from_atoms(0, 1, &[], molecule));
        model.sequence = StructureSequence::from_model_parts(
            &molecule.entity_index,
            &molecule.entities,
            &molecule.entity_poly_seq,
            &model.hierarchy,
            Some(&coarse.hierarchy),
        );
        let properties = StructureProperties::from_hierarchy(&model.hierarchy);
        let ranges = AtomicRanges::from_hierarchy(
            &model.hierarchy,
            &model.sequence,
            model.conformation.xyz_defined,
        );
        let mut units = structure_units(molecule, &model, &coarse, &ranges);
        let (intra_unit_bond_count, inter_unit_bonds, inter_unit_bond_graph) = assign_unit_bonds(
            &mut units,
            &molecule.bonds,
            &molecule.bond_metadata,
            molecule.index_pair_bonds.as_ref(),
        );
        let symmetry_groups = unit_symmetry_groups(&units);
        let boundary = structure_boundary(&units);
        let principal_axes = structure_principal_axes(&units);
        let lookup3d = StructureLookup3D::from_units(&units);
        let polymer_residue_count = units
            .iter()
            .map(|unit| unit.props.polymer_elements.len())
            .sum();
        let polymer_gap_count = units
            .iter()
            .map(|unit| unit.props.gap_elements.len() / 2)
            .sum();
        let element_count = units.iter().map(|unit| unit.elements.len()).sum();
        let mut structure = AtomicStructure {
            model,
            models,
            coarse,
            coarse_models,
            units,
            coordinate_system: coordinate_system_operator(molecule.selected_assembly.as_ref()),
            boundary,
            principal_axes,
            lookup3d,
            ranges,
            properties,
            element_count,
            polymer_residue_count,
            polymer_gap_count,
            symmetry_groups,
            intra_unit_bond_count,
            inter_unit_bonds,
            inter_unit_bond_graph,
            carbohydrates: Carbohydrates::default(),
        };
        structure.carbohydrates = Carbohydrates::from_structure(molecule, &structure);
        structure
    }

    pub fn alt_loc_count(&self) -> usize {
        self.model
            .hierarchy
            .atoms
            .iter()
            .filter(|atom| !atom.alt_id.is_empty())
            .map(|atom| atom.alt_id.as_str())
            .collect::<HashSet<_>>()
            .len()
    }

    pub fn position(&self, unit_id: usize, element_index: usize) -> Option<Vec3> {
        let unit = self.unit_by_id(unit_id)?;
        let source_index = *unit.elements.get(element_index)?;
        let position = match unit.kind {
            UnitKind::Atomic => *self.model.conformation.positions.get(source_index)?,
            UnitKind::Spheres => self.coarse.conformation.spheres.position(source_index)?,
            UnitKind::Gaussians => self.coarse.conformation.gaussians.position(source_index)?,
        };
        Some(unit.operator.transform.apply(position))
    }

    pub fn carbohydrate_element_indices(
        &self,
        unit_id: usize,
        source_atom_index: usize,
    ) -> &[usize] {
        self.carbohydrates
            .get_element_indices(unit_id, source_atom_index)
    }

    pub fn carbohydrate_link_indices(&self, unit_id: usize, source_atom_index: usize) -> &[usize] {
        self.carbohydrates
            .get_link_indices(unit_id, source_atom_index)
    }

    pub fn carbohydrate_terminal_link_indices(
        &self,
        unit_id: usize,
        source_atom_index: usize,
    ) -> &[usize] {
        self.carbohydrates
            .get_terminal_link_indices(unit_id, source_atom_index)
    }

    pub(super) fn unit_element_index(
        &self,
        unit_id: usize,
        source_atom_index: usize,
    ) -> Option<usize> {
        self.unit_by_id(unit_id)?
            .elements
            .iter()
            .position(|element| *element == source_atom_index)
    }

    pub(super) fn inter_unit_bond_exists_exact(
        &self,
        unit_a: usize,
        index_a: Option<usize>,
        unit_b: usize,
        index_b: Option<usize>,
        source_bond: usize,
    ) -> bool {
        let (Some(index_a), Some(index_b)) = (index_a, index_b) else {
            return false;
        };
        self.inter_unit_bonds.iter().any(|bond| {
            bond.source_bond == source_bond
                && ((bond.unit_a == unit_a
                    && bond.index_a == index_a
                    && bond.unit_b == unit_b
                    && bond.index_b == index_b)
                    || (bond.unit_a == unit_b
                        && bond.index_a == index_b
                        && bond.unit_b == unit_a
                        && bond.index_b == index_a))
        })
    }

    pub fn units_are_molstar_sorted(&self) -> bool {
        units_are_molstar_sorted(&self.units)
    }

    pub fn unit_index_map(&self) -> BTreeMap<usize, usize> {
        self.units
            .iter()
            .enumerate()
            .map(|(index, unit)| (unit.id, index))
            .collect()
    }

    pub fn unit_index_by_id(&self, unit_id: usize) -> Option<usize> {
        if self.units_are_molstar_sorted() {
            self.units
                .binary_search_by_key(&unit_id, |unit| unit.id)
                .ok()
        } else {
            self.units.iter().position(|unit| unit.id == unit_id)
        }
    }

    pub fn unit_by_id(&self, unit_id: usize) -> Option<&StructureUnit> {
        self.unit_index_by_id(unit_id)
            .and_then(|index| self.units.get(index))
    }

    pub fn structure_element_loci(&self) -> Vec<(usize, Vec<usize>)> {
        self.units
            .iter()
            .map(|unit| (unit.id, (0..unit.elements.len()).collect()))
            .collect()
    }

    pub fn remap_element_loci_to(
        &self,
        loci: &[(usize, Vec<usize>)],
        target: &AtomicStructure,
    ) -> Vec<(usize, Vec<usize>)> {
        loci.iter()
            .filter_map(|(unit_id, indices)| {
                let source_unit = self.unit_by_id(*unit_id)?;
                let target_unit = target.unit_by_id(*unit_id)?;
                let remapped = remap_unit_indices(source_unit, indices, target_unit);
                (!remapped.is_empty()).then_some((*unit_id, remapped))
            })
            .collect()
    }

    pub fn remap_bond_loci_to(
        &self,
        bonds: &[(usize, usize, usize, usize)],
        target: &AtomicStructure,
    ) -> Vec<(usize, usize, usize, usize)> {
        bonds
            .iter()
            .filter_map(|&(unit_a_id, index_a, unit_b_id, index_b)| {
                let source_unit_a = self.unit_by_id(unit_a_id)?;
                let source_unit_b = self.unit_by_id(unit_b_id)?;
                let target_unit_a = target.unit_by_id(unit_a_id)?;
                let target_unit_b = target.unit_by_id(unit_b_id)?;
                let element_a = *source_unit_a.elements.get(index_a)?;
                let element_b = *source_unit_b.elements.get(index_b)?;
                let target_index_a = target_unit_a
                    .elements
                    .iter()
                    .position(|element| *element == element_a)?;
                let target_index_b = target_unit_b
                    .elements
                    .iter()
                    .position(|element| *element == element_b)?;
                Some((unit_a_id, target_index_a, unit_b_id, target_index_b))
            })
            .collect()
    }

    pub fn extend_element_loci_to_whole_chains(
        &self,
        loci: &[(usize, Vec<usize>)],
    ) -> Vec<(usize, Vec<usize>)> {
        let loci = normalized_element_loci(self, loci);
        let mut elements = Vec::new();
        let mut index = 0usize;
        while index < loci.len() {
            let (unit_id, _) = &loci[index];
            let Some(unit) = self.unit_by_id(*unit_id) else {
                index += 1;
                continue;
            };
            if unit.traits.contains(UnitTraits::PARTITIONED) {
                let start = index;
                index += 1;
                while index < loci.len() {
                    let Some(next) = self.unit_by_id(loci[index].0) else {
                        break;
                    };
                    if !next.are_same_chain_operator_group(unit) {
                        break;
                    }
                    index += 1;
                }
                let chain_indices = selected_chain_indices(self, &loci[start..index]);
                for candidate in &self.units {
                    if candidate.are_same_chain_operator_group(unit) {
                        collect_unit_chain_indices(candidate, &chain_indices, &mut elements);
                    }
                }
            } else {
                let chain_indices = selected_chain_indices(self, &loci[index..index + 1]);
                collect_unit_chain_indices(unit, &chain_indices, &mut elements);
                index += 1;
            }
        }
        elements
    }

    pub fn serial_mapping(&self) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
        let mut cumulative_unit_element_count = Vec::with_capacity(self.units.len());
        let mut unit_indices = Vec::with_capacity(self.element_count);
        let mut element_indices = Vec::with_capacity(self.element_count);
        let mut offset = 0usize;
        for (unit_index, unit) in self.units.iter().enumerate() {
            cumulative_unit_element_count.push(offset);
            for &element in &unit.elements {
                unit_indices.push(unit_index);
                element_indices.push(element);
                offset += 1;
            }
        }
        (cumulative_unit_element_count, unit_indices, element_indices)
    }

    pub fn serial_index(&self, unit_id: usize, element: usize) -> Option<usize> {
        let unit_index = self.unit_index_by_id(unit_id)?;
        let unit = self.units.get(unit_index)?;
        let in_unit_index = unit.elements.iter().position(|value| *value == element)?;
        Some(
            self.units
                .iter()
                .take(unit_index)
                .map(|unit| unit.elements.len())
                .sum::<usize>()
                + in_unit_index,
        )
    }

    pub fn are_unit_ids_equal(&self, other: &AtomicStructure) -> bool {
        self.element_count == other.element_count
            && self.units.len() == other.units.len()
            && self
                .units
                .iter()
                .zip(&other.units)
                .all(|(a, b)| a.id == b.id)
    }

    pub fn are_unit_ids_and_indices_equal(&self, other: &AtomicStructure) -> bool {
        self.are_unit_ids_equal(other)
            && self
                .units
                .iter()
                .zip(&other.units)
                .all(|(a, b)| a.elements == b.elements)
    }
}

fn normalized_element_loci(
    structure: &AtomicStructure,
    loci: &[(usize, Vec<usize>)],
) -> Vec<(usize, Vec<usize>)> {
    let mut by_unit = BTreeMap::<usize, BTreeSet<usize>>::new();
    for (unit_id, indices) in loci {
        let Some(unit) = structure.unit_by_id(*unit_id) else {
            continue;
        };
        let entry = by_unit.entry(*unit_id).or_default();
        for &index in indices {
            if index < unit.elements.len() {
                entry.insert(index);
            }
        }
    }
    structure
        .units
        .iter()
        .filter_map(|unit| {
            by_unit
                .get(&unit.id)
                .map(|indices| (unit.id, indices.iter().copied().collect::<Vec<_>>()))
                .filter(|(_, indices)| !indices.is_empty())
        })
        .collect()
}

fn remap_unit_indices(
    source_unit: &StructureUnit,
    indices: &[usize],
    target_unit: &StructureUnit,
) -> Vec<usize> {
    let selected = indices
        .iter()
        .filter_map(|&index| source_unit.elements.get(index).copied())
        .collect::<BTreeSet<_>>();
    target_unit
        .elements
        .iter()
        .enumerate()
        .filter_map(|(index, element)| selected.contains(element).then_some(index))
        .collect()
}

fn selected_chain_indices(
    structure: &AtomicStructure,
    loci: &[(usize, Vec<usize>)],
) -> BTreeSet<usize> {
    let mut chain_indices = BTreeSet::new();
    for (unit_id, indices) in loci {
        let Some(unit) = structure.unit_by_id(*unit_id) else {
            continue;
        };
        for &index in indices {
            if let Some(chain_index) = unit.chain_index_at(index) {
                chain_indices.insert(chain_index);
            }
        }
    }
    chain_indices
}

fn collect_unit_chain_indices(
    unit: &StructureUnit,
    chain_indices: &BTreeSet<usize>,
    elements: &mut Vec<(usize, Vec<usize>)>,
) {
    let indices = (0..unit.elements.len())
        .filter(|&index| {
            unit.chain_index_at(index)
                .is_some_and(|chain_index| chain_indices.contains(&chain_index))
        })
        .collect::<Vec<_>>();
    if !indices.is_empty() {
        elements.push((unit.id, indices));
    }
}

#[derive(Clone, Debug, Default)]
pub struct AtomicModel {
    pub id: String,
    pub model_num: i32,
    pub source_data: SourceData,
    pub global_model_transform: Option<GlobalModelTransform>,
    pub sequence: StructureSequence,
    pub hierarchy: AtomicHierarchy,
    pub secondary_structure: SecondaryStructure,
    pub conformation: AtomicConformation,
    pub custom_properties: CustomPropertyContainer,
    pub static_property_data: ModelPropertyData,
    pub dynamic_property_data: ModelPropertyData,
}

fn sorted_model_atoms(atoms: &[Atom]) -> Vec<&Atom> {
    atoms.iter().collect()
}

static MOLSTAR_UUID_COUNTER: AtomicU64 = AtomicU64::new(1);

fn molstar_uuid_create22() -> String {
    let counter = MOLSTAR_UUID_COUNTER.fetch_add(1, Ordering::Relaxed);
    let mut state = counter.wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ 0xd1b5_4a32_d192_ed03;
    let mut bytes = [0u8; 16];
    for byte in &mut bytes {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        state = state.wrapping_mul(0x2545_f491_4f6c_dd1d);
        *byte = state as u8;
    }
    base64_url_22(&bytes)
}

fn base64_url_22(bytes: &[u8; 16]) -> String {
    const CHARS: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_";
    let mut out = String::with_capacity(22);
    let mut index = 0;
    while index + 3 <= bytes.len() {
        let value = ((bytes[index] as u32) << 16)
            | ((bytes[index + 1] as u32) << 8)
            | bytes[index + 2] as u32;
        out.push(CHARS[((value >> 18) & 0x3f) as usize] as char);
        out.push(CHARS[((value >> 12) & 0x3f) as usize] as char);
        out.push(CHARS[((value >> 6) & 0x3f) as usize] as char);
        out.push(CHARS[(value & 0x3f) as usize] as char);
        index += 3;
    }
    if index < bytes.len() {
        let value = (bytes[index] as u32) << 16;
        out.push(CHARS[((value >> 18) & 0x3f) as usize] as char);
        out.push(CHARS[((value >> 12) & 0x3f) as usize] as char);
    }
    out
}

impl AtomicModel {
    fn from_atoms(
        _model_index: usize,
        model_num: i32,
        atoms: &[Atom],
        molecule: &Molecule,
    ) -> Self {
        let mut hierarchy = AtomicHierarchy::default();
        let mut atom_ids = Vec::with_capacity(atoms.len());
        let mut positions = Vec::with_capacity(atoms.len());
        let mut x = Vec::with_capacity(atoms.len());
        let mut y = Vec::with_capacity(atoms.len());
        let mut z = Vec::with_capacity(atoms.len());
        let mut occupancies = Vec::with_capacity(atoms.len());
        let mut b_iso = Vec::with_capacity(atoms.len());
        let mut formal_charges = Vec::with_capacity(atoms.len());
        let mut residue_keys = Vec::<ResidueKey>::new();
        let mut chain_keys = Vec::<ChainKey>::new();

        let sorted_atoms = sorted_model_atoms(atoms);
        for (element_index, atom) in sorted_atoms.iter().copied().enumerate() {
            let chain_key = ChainKey {
                id: atom.chain.clone(),
                entity_id: atom.entity_id.clone(),
            };
            let chain_index = chain_keys
                .iter()
                .position(|key| key == &chain_key)
                .unwrap_or_else(|| {
                    chain_keys.push(chain_key.clone());
                    hierarchy.chains.push(AtomicChain {
                        id: chain_key.id,
                        auth_id: atom.auth_chain.clone(),
                        entity_id: chain_key.entity_id,
                        start_residue: usize::MAX,
                        end_residue: 0,
                    });
                    chain_keys.len() - 1
                });

            let residue_key = ResidueKey {
                chain_index,
                label_seq_id: atom.residue_seq.clone(),
                auth_seq_id: atom.auth_residue_seq.clone(),
                insertion_code: atom.insertion_code.clone(),
            };
            let residue_index = residue_keys
                .iter()
                .position(|key| key == &residue_key)
                .unwrap_or_else(|| {
                    residue_keys.push(residue_key.clone());
                    hierarchy.residues.push(AtomicResidue {
                        chain_index,
                        comp_id: atom.residue.clone(),
                        auth_comp_id: atom.auth_residue.clone(),
                        group_pdb: atom.group_pdb.clone(),
                        label_seq_id: residue_key.label_seq_id,
                        auth_seq_id: residue_key.auth_seq_id,
                        insertion_code: residue_key.insertion_code,
                        start_atom: element_index,
                        end_atom: element_index + 1,
                        is_het: atom.het,
                    });
                    residue_keys.len() - 1
                });

            if let Some(residue) = hierarchy.residues.get_mut(residue_index) {
                residue.end_atom = element_index + 1;
                residue.is_het &= atom.het;
            }
            if let Some(chain) = hierarchy.chains.get_mut(chain_index) {
                chain.start_residue = chain.start_residue.min(residue_index);
                chain.end_residue = chain.end_residue.max(residue_index + 1);
            }

            hierarchy.atoms.push(AtomicAtom {
                source_index: atom.source_index,
                id: atom.id,
                model_num: atom.model_num,
                name: atom.name.clone(),
                auth_name: atom.auth_name.clone(),
                type_symbol: atom.type_symbol.clone(),
                element: atom.element.clone(),
                label_comp_id: atom.residue.clone(),
                auth_comp_id: atom.auth_residue.clone(),
                formal_charge: atom.formal_charge,
                alt_id: atom.alt_id.clone(),
                operator_name: atom_operator_name(atom),
                position: atom.position,
                residue_index,
                chain_index,
            });
            atom_ids.push(atom.id);
            positions.push(atom.position);
            x.push(atom.position.x);
            y.push(atom.position.y);
            z.push(atom.position.z);
            occupancies.push(atom.occupancy);
            b_iso.push(atom.b_iso);
            formal_charges.push(atom.formal_charge);
        }

        for chain in &mut hierarchy.chains {
            if chain.start_residue == usize::MAX {
                chain.start_residue = 0;
            }
        }

        hierarchy.build_segments();
        hierarchy.index =
            AtomicIndex::from_hierarchy(&hierarchy, &molecule.entity_index, &molecule.entities);
        hierarchy.derived =
            AtomicDerivedData::from_hierarchy(&hierarchy, &molecule.chemical_components);
        let sequence = StructureSequence::from_model_parts(
            &molecule.entity_index,
            &molecule.entities,
            &molecule.entity_poly_seq,
            &hierarchy,
            None,
        );
        let secondary_structure =
            SecondaryStructure::from_hierarchy(&hierarchy, &molecule.helices, &molecule.sheets);
        let (element_to_anisotrop, anisotropic_displacement) =
            atom_site_anisotrop_mapping(&sorted_atoms, &molecule.atom_site_anisotrop);

        let global_model_transform = molecule.global_model_transform.clone();
        let mut custom_properties = CustomPropertyContainer::default();
        let mut static_property_data = ModelPropertyData::default();
        if let Some(transform) = &global_model_transform {
            custom_properties.register(GlobalModelTransform::DESCRIPTOR);
            static_property_data.insert(
                GlobalModelTransform::DESCRIPTOR,
                transform.to_property_string(),
            );
        }

        let conformation_id = molstar_uuid_create22();
        let model_id = molstar_uuid_create22();

        AtomicModel {
            id: model_id,
            model_num,
            source_data: molecule.source_data.clone(),
            global_model_transform,
            sequence,
            hierarchy,
            secondary_structure,
            conformation: AtomicConformation {
                id: conformation_id,
                atom_ids,
                positions,
                x,
                y,
                z,
                occupancies,
                b_iso,
                formal_charges,
                occupancy_defined: molecule.atom_site_columns.occupancy_defined,
                b_iso_defined: molecule.atom_site_columns.b_iso_defined,
                xyz_defined: molecule.atom_site_columns.xyz_defined,
                element_to_anisotrop,
                anisotropic_displacement,
            },
            custom_properties,
            static_property_data,
            dynamic_property_data: ModelPropertyData::default(),
        }
    }
}

fn atom_site_anisotrop_mapping(
    atoms: &[&Atom],
    anisotrop: &[AtomSiteAnisotrop],
) -> (Vec<i32>, Vec<Option<AnisotropicDisplacement>>) {
    let mut atom_id_to_element = BTreeMap::new();
    for (element_index, atom) in atoms.iter().enumerate() {
        atom_id_to_element.insert(atom.id, element_index);
    }

    let mut element_to_anisotrop = vec![-1; atoms.len()];
    let mut anisotropic_displacement = vec![None; atoms.len()];
    for (anisotrop_index, row) in anisotrop.iter().enumerate() {
        let Some(&element_index) = atom_id_to_element.get(&row.atom_id) else {
            continue;
        };
        element_to_anisotrop[element_index] = anisotrop_index as i32;
        anisotropic_displacement[element_index] = Some(row.u);
    }

    (element_to_anisotrop, anisotropic_displacement)
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CustomPropertyContainer {
    pub descriptors: Vec<String>,
}

impl CustomPropertyContainer {
    pub fn register(&mut self, descriptor: impl Into<String>) {
        let descriptor = descriptor.into();
        if !self.descriptors.iter().any(|name| name == &descriptor) {
            self.descriptors.push(descriptor);
        }
    }

    pub fn contains(&self, descriptor: &str) -> bool {
        self.descriptors.iter().any(|name| name == descriptor)
    }

    pub fn is_empty(&self) -> bool {
        self.descriptors.is_empty()
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ModelPropertyData {
    pub values: BTreeMap<String, String>,
}

impl ModelPropertyData {
    pub fn insert(&mut self, name: impl Into<String>, value: impl Into<String>) {
        self.values.insert(name.into(), value.into());
    }

    pub fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name).map(String::as_str)
    }

    pub fn contains_key(&self, name: &str) -> bool {
        self.values.contains_key(name)
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}

#[derive(Clone, Debug, Default)]
pub struct AtomicHierarchy {
    pub atoms: Vec<AtomicAtom>,
    pub residues: Vec<AtomicResidue>,
    pub chains: Vec<AtomicChain>,
    pub derived: AtomicDerivedData,
    pub atom_source_index: Vec<usize>,
    pub residue_source_index: Vec<usize>,
    pub residue_atom_segments: AtomicSegmentation,
    pub chain_atom_segments: AtomicSegmentation,
    pub index: AtomicIndex,
}

impl AtomicHierarchy {
    pub fn tables(&self) -> AtomicHierarchyTables {
        AtomicHierarchyTables {
            atoms: self.atoms_table(),
            residues: self.residues_table(),
            chains: self.chains_table(),
        }
    }

    pub fn atoms_table(&self) -> AtomicAtomsTable {
        let mut table = AtomicAtomsTable::with_capacity(self.atoms.len());
        for atom in &self.atoms {
            table.type_symbol.push(atom.type_symbol.clone());
            table.label_atom_id.push(atom.name.clone());
            table.auth_atom_id.push(atom.auth_name.clone());
            table.label_alt_id.push(atom.alt_id.clone());
            table.label_comp_id.push(atom.label_comp_id.clone());
            table.auth_comp_id.push(atom.auth_comp_id.clone());
            table.pdbx_formal_charge.push(atom.formal_charge);
        }
        table
    }

    pub fn residues_table(&self) -> AtomicResiduesTable {
        let mut table = AtomicResiduesTable::with_capacity(self.residues.len());
        for residue in &self.residues {
            table.group_pdb.push(residue.group_pdb.clone());
            table.label_seq_id.push(residue.label_seq_id.clone());
            table.auth_seq_id.push(residue.auth_seq_id.clone());
            table.pdbx_pdb_ins_code.push(residue.insertion_code.clone());
        }
        table
    }

    pub fn chains_table(&self) -> AtomicChainsTable {
        let mut table = AtomicChainsTable::with_capacity(self.chains.len());
        for chain in &self.chains {
            table.label_asym_id.push(chain.id.clone());
            table.auth_asym_id.push(chain.auth_id.clone());
            table.label_entity_id.push(chain.entity_id.clone());
        }
        table
    }

    fn build_segments(&mut self) {
        let atom_count = self.atoms.len();
        self.atom_source_index = self.atoms.iter().map(|atom| atom.source_index).collect();
        self.residue_source_index = residue_source_index(&self.atom_source_index, &self.residues);
        self.residue_atom_segments =
            AtomicSegmentation::from_residue_ranges(atom_count, &self.residues);
        self.chain_atom_segments =
            AtomicSegmentation::from_chain_ranges(atom_count, &self.chains, &self.residues);
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomicHierarchyTables {
    pub atoms: AtomicAtomsTable,
    pub residues: AtomicResiduesTable,
    pub chains: AtomicChainsTable,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomicAtomsTable {
    pub type_symbol: Vec<String>,
    pub label_atom_id: Vec<String>,
    pub auth_atom_id: Vec<String>,
    pub label_alt_id: Vec<String>,
    pub label_comp_id: Vec<String>,
    pub auth_comp_id: Vec<String>,
    pub pdbx_formal_charge: Vec<i32>,
}

impl AtomicAtomsTable {
    pub const COLUMN_NAMES: &'static [&'static str] = &[
        "type_symbol",
        "label_atom_id",
        "auth_atom_id",
        "label_alt_id",
        "label_comp_id",
        "auth_comp_id",
        "pdbx_formal_charge",
    ];

    fn with_capacity(capacity: usize) -> Self {
        AtomicAtomsTable {
            type_symbol: Vec::with_capacity(capacity),
            label_atom_id: Vec::with_capacity(capacity),
            auth_atom_id: Vec::with_capacity(capacity),
            label_alt_id: Vec::with_capacity(capacity),
            label_comp_id: Vec::with_capacity(capacity),
            auth_comp_id: Vec::with_capacity(capacity),
            pdbx_formal_charge: Vec::with_capacity(capacity),
        }
    }

    pub fn row_count(&self) -> usize {
        self.type_symbol.len()
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomicResiduesTable {
    pub group_pdb: Vec<String>,
    pub label_seq_id: Vec<String>,
    pub auth_seq_id: Vec<String>,
    pub pdbx_pdb_ins_code: Vec<String>,
}

impl AtomicResiduesTable {
    pub const COLUMN_NAMES: &'static [&'static str] = &[
        "group_PDB",
        "label_seq_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
    ];

    fn with_capacity(capacity: usize) -> Self {
        AtomicResiduesTable {
            group_pdb: Vec::with_capacity(capacity),
            label_seq_id: Vec::with_capacity(capacity),
            auth_seq_id: Vec::with_capacity(capacity),
            pdbx_pdb_ins_code: Vec::with_capacity(capacity),
        }
    }

    pub fn row_count(&self) -> usize {
        self.group_pdb.len()
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomicChainsTable {
    pub label_asym_id: Vec<String>,
    pub auth_asym_id: Vec<String>,
    pub label_entity_id: Vec<String>,
}

impl AtomicChainsTable {
    pub const COLUMN_NAMES: &'static [&'static str] =
        &["label_asym_id", "auth_asym_id", "label_entity_id"];

    fn with_capacity(capacity: usize) -> Self {
        AtomicChainsTable {
            label_asym_id: Vec::with_capacity(capacity),
            auth_asym_id: Vec::with_capacity(capacity),
            label_entity_id: Vec::with_capacity(capacity),
        }
    }

    pub fn row_count(&self) -> usize {
        self.label_asym_id.len()
    }
}

fn residue_source_index(atom_source_index: &[usize], residues: &[AtomicResidue]) -> Vec<usize> {
    let mut first_source_by_residue = residues
        .iter()
        .enumerate()
        .map(|(residue_index, residue)| {
            let first_source = atom_source_index
                .iter()
                .take(residue.end_atom.min(atom_source_index.len()))
                .skip(residue.start_atom)
                .copied()
                .min()
                .unwrap_or(residue_index);
            (residue_index, first_source)
        })
        .collect::<Vec<_>>();
    first_source_by_residue.sort_by_key(|(_, first_source)| *first_source);
    let mut residue_source_index = vec![0; residues.len()];
    for (rank, (residue_index, _)) in first_source_by_residue.into_iter().enumerate() {
        residue_source_index[residue_index] = rank;
    }
    residue_source_index
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomicIndex {
    pub chain_entity_index: Vec<Option<usize>>,
    pub chain_entity_type: Vec<String>,
    pub label_asym_to_entity_index: BTreeMap<String, usize>,
    pub label_asym_to_chain_index: BTreeMap<String, usize>,
    pub entity_label_asym_to_chain_index: BTreeMap<String, usize>,
    pub auth_asym_auth_seq_to_chain_index: BTreeMap<String, usize>,
    pub chain_label_seq_to_residue_index: BTreeMap<String, usize>,
    pub chain_auth_seq_to_residue_index: BTreeMap<String, usize>,
}

impl AtomicIndex {
    fn from_hierarchy(
        hierarchy: &AtomicHierarchy,
        entity_index: &EntityIndexMap,
        entities: &[Entity],
    ) -> Self {
        let mut chain_entity_index = Vec::with_capacity(hierarchy.chains.len());
        let mut chain_entity_type = Vec::with_capacity(hierarchy.chains.len());
        let mut label_asym_to_entity_index = BTreeMap::new();
        let mut label_asym_to_chain_index = BTreeMap::new();
        let mut entity_label_asym_to_chain_index = BTreeMap::new();
        let mut auth_asym_auth_seq_to_chain_index = BTreeMap::new();
        let mut chain_label_seq_to_residue_index = BTreeMap::new();
        let mut chain_auth_seq_to_residue_index = BTreeMap::new();
        for (chain_index, chain) in hierarchy.chains.iter().enumerate() {
            let entity_key = entity_index.get_entity_index(&chain.entity_id);
            chain_entity_index.push(entity_key);
            chain_entity_type.push(
                entity_key
                    .and_then(|index| entities.get(index))
                    .map(|entity| entity.type_name.clone())
                    .unwrap_or_default(),
            );
            label_asym_to_chain_index
                .entry(chain.id.clone())
                .or_insert(chain_index);
            if let Some(entity_key) = entity_key {
                label_asym_to_entity_index
                    .entry(chain.id.clone())
                    .or_insert(entity_key);
                entity_label_asym_to_chain_index
                    .insert(atomic_index_key(entity_key, &chain.id), chain_index);
            }
        }
        for (residue_index, residue) in hierarchy.residues.iter().enumerate() {
            let Some(chain) = hierarchy.chains.get(residue.chain_index) else {
                continue;
            };
            auth_asym_auth_seq_to_chain_index
                .entry(sequence_key(&chain.auth_id, &residue.auth_seq_id))
                .or_insert(residue.chain_index);
            chain_label_seq_to_residue_index
                .entry(residue_index_key(
                    residue.chain_index,
                    &residue.label_seq_id,
                    &residue.insertion_code,
                ))
                .or_insert(residue_index);
            chain_auth_seq_to_residue_index
                .entry(residue_index_key(
                    residue.chain_index,
                    &residue.auth_seq_id,
                    &residue.insertion_code,
                ))
                .or_insert(residue_index);
        }
        AtomicIndex {
            chain_entity_index,
            chain_entity_type,
            label_asym_to_entity_index,
            label_asym_to_chain_index,
            entity_label_asym_to_chain_index,
            auth_asym_auth_seq_to_chain_index,
            chain_label_seq_to_residue_index,
            chain_auth_seq_to_residue_index,
        }
    }

    pub fn entity_from_chain(&self, chain_index: usize) -> Option<usize> {
        self.chain_entity_index.get(chain_index).copied().flatten()
    }

    pub fn entity_type_from_chain(&self, chain_index: usize) -> Option<&str> {
        self.chain_entity_type
            .get(chain_index)
            .map(|entity_type| entity_type.as_str())
    }

    pub fn entity_by_label_asym_id(&self, label_asym_id: &str) -> Option<usize> {
        self.label_asym_to_entity_index.get(label_asym_id).copied()
    }

    pub fn chain_by_label_asym_id(&self, label_asym_id: &str) -> Option<usize> {
        self.label_asym_to_chain_index.get(label_asym_id).copied()
    }

    pub fn chain_by_entity_and_label_asym_id(
        &self,
        entity_index: usize,
        label_asym_id: &str,
    ) -> Option<usize> {
        self.entity_label_asym_to_chain_index
            .get(&atomic_index_key(entity_index, label_asym_id))
            .copied()
    }

    pub fn chain_by_auth_asym_and_seq_id(
        &self,
        auth_asym_id: &str,
        auth_seq_id: &str,
    ) -> Option<usize> {
        self.auth_asym_auth_seq_to_chain_index
            .get(&sequence_key(auth_asym_id, auth_seq_id))
            .copied()
    }

    pub fn residue_by_label_key(
        &self,
        chain_index: usize,
        label_seq_id: &str,
        insertion_code: &str,
    ) -> Option<usize> {
        self.chain_label_seq_to_residue_index
            .get(&residue_index_key(
                chain_index,
                label_seq_id,
                insertion_code,
            ))
            .copied()
    }

    pub fn residue_by_auth_key(
        &self,
        chain_index: usize,
        auth_seq_id: &str,
        insertion_code: &str,
    ) -> Option<usize> {
        self.chain_auth_seq_to_residue_index
            .get(&residue_index_key(chain_index, auth_seq_id, insertion_code))
            .copied()
    }

    pub fn residue_by_entity_label_asym_and_auth_seq(
        &self,
        entity_index: usize,
        label_asym_id: &str,
        auth_seq_id: &str,
        insertion_code: &str,
    ) -> Option<usize> {
        let chain_index = self.chain_by_entity_and_label_asym_id(entity_index, label_asym_id)?;
        self.residue_by_auth_key(chain_index, auth_seq_id, insertion_code)
    }

    pub fn atom_by_label_key(
        &self,
        hierarchy: &AtomicHierarchy,
        chain_index: usize,
        label_seq_id: &str,
        insertion_code: &str,
        label_atom_id: &str,
        label_alt_id: Option<&str>,
    ) -> Option<usize> {
        let residue_index = self.residue_by_label_key(chain_index, label_seq_id, insertion_code)?;
        self.atom_on_residue(hierarchy, residue_index, label_atom_id, false, label_alt_id)
    }

    pub fn atom_by_auth_key(
        &self,
        hierarchy: &AtomicHierarchy,
        auth_asym_id: &str,
        auth_seq_id: &str,
        insertion_code: &str,
        auth_atom_id: &str,
        label_alt_id: Option<&str>,
    ) -> Option<usize> {
        let chain_index = self.chain_by_auth_asym_and_seq_id(auth_asym_id, auth_seq_id)?;
        let residue_index = self.residue_by_auth_key(chain_index, auth_seq_id, insertion_code)?;
        self.atom_on_residue(hierarchy, residue_index, auth_atom_id, true, label_alt_id)
    }

    pub fn atom_on_residue(
        &self,
        hierarchy: &AtomicHierarchy,
        residue_index: usize,
        atom_id: &str,
        use_auth_id: bool,
        label_alt_id: Option<&str>,
    ) -> Option<usize> {
        let residue = hierarchy.residues.get(residue_index)?;
        (residue.start_atom..residue.end_atom).find(|atom_index| {
            hierarchy.atoms.get(*atom_index).is_some_and(|atom| {
                let name = if use_auth_id {
                    atom.auth_name.as_str()
                } else {
                    atom.name.as_str()
                };
                name == atom_id
                    && label_alt_id
                        .filter(|alt| !alt.is_empty())
                        .is_none_or(|alt| atom.alt_id == alt)
            })
        })
    }
}

fn atomic_index_key(entity_index: usize, label_asym_id: &str) -> String {
    format!("{entity_index}\u{1f}{label_asym_id}")
}

fn sequence_key(asym_id: &str, seq_id: &str) -> String {
    format!("{asym_id}\u{1f}{seq_id}")
}

fn residue_index_key(chain_index: usize, seq_id: &str, insertion_code: &str) -> String {
    format!("{chain_index}\u{1f}{seq_id}\u{1f}{insertion_code}")
}

mod derived;

use derived::{are_backbone_connected, chain_residue_indices, is_polymer_residue, residue_seq_id};
pub use derived::{
    get_saccharide_name, get_saccharide_shape, saccharide_component, saccharide_component_with_map,
    AtomDerivedData, AtomicDerivedData, MoleculeType, PolymerType, ResidueDerivedData,
    SaccharideCompIdMapType, SaccharideComponent, SaccharideShape, SaccharideType,
};

#[derive(Clone, Debug)]
pub struct AtomicAtom {
    pub source_index: usize,
    pub id: usize,
    pub model_num: i32,
    pub name: String,
    pub auth_name: String,
    pub type_symbol: String,
    pub element: String,
    pub label_comp_id: String,
    pub auth_comp_id: String,
    pub formal_charge: i32,
    pub alt_id: String,
    pub operator_name: String,
    pub position: Vec3,
    pub residue_index: usize,
    pub chain_index: usize,
}

#[derive(Clone, Debug)]
pub struct AtomicResidue {
    pub chain_index: usize,
    pub comp_id: String,
    pub auth_comp_id: String,
    pub group_pdb: String,
    pub label_seq_id: String,
    pub auth_seq_id: String,
    pub insertion_code: String,
    pub start_atom: usize,
    pub end_atom: usize,
    pub is_het: bool,
}

#[derive(Clone, Debug)]
pub struct AtomicChain {
    pub id: String,
    pub auth_id: String,
    pub entity_id: String,
    pub start_residue: usize,
    pub end_residue: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SecondaryStructure {
    pub residue_type: Vec<SecondaryStructureType>,
    pub key: Vec<usize>,
    pub elements: Vec<SecondaryStructureElement>,
}

impl SecondaryStructure {
    fn from_hierarchy(
        hierarchy: &AtomicHierarchy,
        helices: &[SecondaryRange],
        sheets: &[SecondaryRange],
    ) -> Self {
        let residue_count = hierarchy.residues.len();
        let mut secondary_structure = SecondaryStructure {
            residue_type: vec![SecondaryStructureType::NONE; residue_count],
            key: vec![0; residue_count],
            elements: vec![SecondaryStructureElement::None],
        };

        for range in helices {
            let key = secondary_structure.elements.len();
            secondary_structure
                .elements
                .push(SecondaryStructureElement::Helix);
            secondary_structure.assign_range(range, SecondaryStructureType::HELIX, key, hierarchy);
        }
        for range in sheets {
            let key = secondary_structure.elements.len();
            secondary_structure
                .elements
                .push(SecondaryStructureElement::Sheet);
            secondary_structure.assign_range(
                range,
                SecondaryStructureType::BETA_SHEET,
                key,
                hierarchy,
            );
        }

        secondary_structure
    }

    fn assign_range(
        &mut self,
        range: &SecondaryRange,
        secondary_type: SecondaryStructureType,
        key: usize,
        hierarchy: &AtomicHierarchy,
    ) {
        let Some(start_residue) = secondary_range_start_residue(hierarchy, range) else {
            return;
        };
        let Some(chain) = hierarchy
            .chains
            .get(hierarchy.residues[start_residue].chain_index)
        else {
            return;
        };

        let mut residue_index = start_residue;
        while residue_index < chain.end_residue && residue_index < hierarchy.residues.len() {
            let residue = &hierarchy.residues[residue_index];
            self.residue_type[residue_index] = secondary_type;
            self.key[residue_index] = key;

            if secondary_range_reaches_end(residue, range) {
                break;
            }

            residue_index += 1;
        }
    }

    pub fn residue_type(&self, residue_index: usize) -> SecondaryStructureType {
        self.residue_type
            .get(residue_index)
            .copied()
            .unwrap_or(SecondaryStructureType::NONE)
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SecondaryStructureType {
    pub bits: u32,
}

impl SecondaryStructureType {
    pub const NONE: SecondaryStructureType = SecondaryStructureType { bits: 0x0 };
    pub const HELIX: SecondaryStructureType = SecondaryStructureType { bits: 0x2 };
    pub const BETA: SecondaryStructureType = SecondaryStructureType { bits: 0x4 };
    pub const BETA_SHEET: SecondaryStructureType = SecondaryStructureType {
        bits: Self::BETA.bits | 0x800000,
    };

    pub const fn contains(self, other: SecondaryStructureType) -> bool {
        self.bits & other.bits == other.bits
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SecondaryStructureElement {
    None,
    Helix,
    Sheet,
}

fn secondary_range_start_residue(
    hierarchy: &AtomicHierarchy,
    range: &SecondaryRange,
) -> Option<usize> {
    let chain_index = hierarchy.index.chain_by_label_asym_id(&range.chain)?;
    let chain = hierarchy.chains.get(chain_index)?;
    (chain.start_residue..chain.end_residue.min(hierarchy.residues.len())).find(|&residue_index| {
        let residue = &hierarchy.residues[residue_index];
        secondary_residue_seq_id(residue) == Some(range.start)
            && residue.insertion_code == range.start_insertion_code
    })
}

fn secondary_range_reaches_end(residue: &AtomicResidue, range: &SecondaryRange) -> bool {
    let Some(seq_id) = secondary_residue_seq_id(residue) else {
        return false;
    };
    seq_id > range.end
        || (seq_id == range.end && residue.insertion_code == range.end_insertion_code)
}

fn secondary_residue_seq_id(residue: &AtomicResidue) -> Option<i32> {
    residue.label_seq_id.trim().parse().ok()
}

#[derive(Clone, Debug, Default)]
pub struct AtomicConformation {
    pub id: String,
    pub atom_ids: Vec<usize>,
    pub positions: Vec<Vec3>,
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub occupancies: Vec<f32>,
    pub b_iso: Vec<f32>,
    pub formal_charges: Vec<i32>,
    pub occupancy_defined: bool,
    pub b_iso_defined: bool,
    pub xyz_defined: bool,
    pub element_to_anisotrop: Vec<i32>,
    pub anisotropic_displacement: Vec<Option<AnisotropicDisplacement>>,
}

#[derive(Clone, Debug, Default)]
pub struct AtomicSegmentation {
    pub offsets: Vec<usize>,
    pub index: Vec<usize>,
    pub count: usize,
}

impl AtomicSegmentation {
    fn from_residue_ranges(atom_count: usize, residues: &[AtomicResidue]) -> Self {
        let mut offsets = vec![0; residues.len() + 1];
        let mut index = vec![0; atom_count];
        for (residue_index, residue) in residues.iter().enumerate() {
            offsets[residue_index] = residue.start_atom;
            for value in index
                .iter_mut()
                .take(residue.end_atom.min(atom_count))
                .skip(residue.start_atom)
            {
                *value = residue_index;
            }
        }
        if let Some(last) = offsets.last_mut() {
            *last = atom_count;
        }
        AtomicSegmentation {
            offsets,
            index,
            count: residues.len(),
        }
    }

    fn from_chain_ranges(
        atom_count: usize,
        chains: &[AtomicChain],
        residues: &[AtomicResidue],
    ) -> Self {
        let mut offsets = vec![0; chains.len() + 1];
        let mut index = vec![0; atom_count];
        for (chain_index, chain) in chains.iter().enumerate() {
            let start_atom = residues
                .get(chain.start_residue)
                .map(|residue| residue.start_atom)
                .unwrap_or(0);
            let end_atom = chain
                .end_residue
                .checked_sub(1)
                .and_then(|residue_index| residues.get(residue_index))
                .map(|residue| residue.end_atom)
                .unwrap_or(start_atom);
            offsets[chain_index] = start_atom;
            for value in index
                .iter_mut()
                .take(end_atom.min(atom_count))
                .skip(start_atom)
            {
                *value = chain_index;
            }
        }
        if let Some(last) = offsets.last_mut() {
            *last = atom_count;
        }
        AtomicSegmentation {
            offsets,
            index,
            count: chains.len(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct StructureProperties {
    pub atom_source_index: Vec<usize>,
    pub atom_id: Vec<usize>,
    pub residue_index: Vec<usize>,
    pub chain_index: Vec<usize>,
    pub type_symbol: Vec<String>,
    pub label_atom_id: Vec<String>,
    pub auth_atom_id: Vec<String>,
    pub label_alt_id: Vec<String>,
    pub pdbx_formal_charge: Vec<i32>,
    pub group_pdb: Vec<String>,
    pub label_comp_id: Vec<String>,
    pub auth_comp_id: Vec<String>,
    pub label_seq_id: Vec<String>,
    pub auth_seq_id: Vec<String>,
    pub pdbx_pdb_ins_code: Vec<String>,
    pub label_asym_id: Vec<String>,
    pub auth_asym_id: Vec<String>,
    pub label_entity_id: Vec<String>,
}

impl StructureProperties {
    fn from_hierarchy(hierarchy: &AtomicHierarchy) -> Self {
        let mut properties = StructureProperties {
            atom_source_index: Vec::with_capacity(hierarchy.atoms.len()),
            atom_id: Vec::with_capacity(hierarchy.atoms.len()),
            residue_index: Vec::with_capacity(hierarchy.atoms.len()),
            chain_index: Vec::with_capacity(hierarchy.atoms.len()),
            type_symbol: Vec::with_capacity(hierarchy.atoms.len()),
            label_atom_id: Vec::with_capacity(hierarchy.atoms.len()),
            auth_atom_id: Vec::with_capacity(hierarchy.atoms.len()),
            label_alt_id: Vec::with_capacity(hierarchy.atoms.len()),
            pdbx_formal_charge: Vec::with_capacity(hierarchy.atoms.len()),
            group_pdb: Vec::with_capacity(hierarchy.atoms.len()),
            label_comp_id: Vec::with_capacity(hierarchy.atoms.len()),
            auth_comp_id: Vec::with_capacity(hierarchy.atoms.len()),
            label_seq_id: Vec::with_capacity(hierarchy.atoms.len()),
            auth_seq_id: Vec::with_capacity(hierarchy.atoms.len()),
            pdbx_pdb_ins_code: Vec::with_capacity(hierarchy.atoms.len()),
            label_asym_id: Vec::with_capacity(hierarchy.atoms.len()),
            auth_asym_id: Vec::with_capacity(hierarchy.atoms.len()),
            label_entity_id: Vec::with_capacity(hierarchy.atoms.len()),
        };
        for atom in &hierarchy.atoms {
            let residue = &hierarchy.residues[atom.residue_index];
            let chain = &hierarchy.chains[atom.chain_index];
            properties.atom_source_index.push(atom.source_index);
            properties.atom_id.push(atom.id);
            properties.residue_index.push(atom.residue_index);
            properties.chain_index.push(atom.chain_index);
            properties.type_symbol.push(atom.type_symbol.clone());
            properties.label_atom_id.push(atom.name.clone());
            properties.auth_atom_id.push(atom.auth_name.clone());
            properties.label_alt_id.push(atom.alt_id.clone());
            properties.pdbx_formal_charge.push(atom.formal_charge);
            properties.group_pdb.push(residue.group_pdb.clone());
            properties.label_comp_id.push(atom.label_comp_id.clone());
            properties.auth_comp_id.push(atom.auth_comp_id.clone());
            properties.label_seq_id.push(residue.label_seq_id.clone());
            properties.auth_seq_id.push(residue.auth_seq_id.clone());
            properties
                .pdbx_pdb_ins_code
                .push(residue.insertion_code.clone());
            properties.label_asym_id.push(chain.id.clone());
            properties.auth_asym_id.push(chain.auth_id.clone());
            properties.label_entity_id.push(chain.entity_id.clone());
        }
        properties
    }

    pub fn len(&self) -> usize {
        self.atom_id.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atom_id.is_empty()
    }

    pub fn atom_key(&self, element: usize) -> Option<usize> {
        (element < self.len()).then_some(element)
    }

    pub fn atom_id(&self, element: usize) -> Option<usize> {
        self.atom_id.get(element).copied()
    }

    pub fn atom_source_index(&self, element: usize) -> Option<usize> {
        self.atom_source_index.get(element).copied()
    }

    pub fn atom_type_symbol(&self, element: usize) -> Option<&str> {
        self.type_symbol.get(element).map(String::as_str)
    }

    pub fn atom_label_atom_id(&self, element: usize) -> Option<&str> {
        self.label_atom_id.get(element).map(String::as_str)
    }

    pub fn atom_auth_atom_id(&self, element: usize) -> Option<&str> {
        self.auth_atom_id.get(element).map(String::as_str)
    }

    pub fn atom_label_alt_id(&self, element: usize) -> Option<&str> {
        self.label_alt_id.get(element).map(String::as_str)
    }

    pub fn atom_label_comp_id(&self, element: usize) -> Option<&str> {
        self.label_comp_id.get(element).map(String::as_str)
    }

    pub fn atom_auth_comp_id(&self, element: usize) -> Option<&str> {
        self.auth_comp_id.get(element).map(String::as_str)
    }

    pub fn atom_formal_charge(&self, element: usize) -> Option<i32> {
        self.pdbx_formal_charge.get(element).copied()
    }

    pub fn residue_key(&self, element: usize) -> Option<usize> {
        self.residue_index.get(element).copied()
    }

    pub fn residue_group_pdb(&self, element: usize) -> Option<&str> {
        self.group_pdb.get(element).map(String::as_str)
    }

    pub fn residue_label_comp_id(&self, element: usize) -> Option<&str> {
        self.label_comp_id.get(element).map(String::as_str)
    }

    pub fn residue_auth_comp_id(&self, element: usize) -> Option<&str> {
        self.auth_comp_id.get(element).map(String::as_str)
    }

    pub fn residue_label_seq_id(&self, element: usize) -> Option<&str> {
        self.label_seq_id.get(element).map(String::as_str)
    }

    pub fn residue_auth_seq_id(&self, element: usize) -> Option<&str> {
        self.auth_seq_id.get(element).map(String::as_str)
    }

    pub fn residue_pdb_ins_code(&self, element: usize) -> Option<&str> {
        self.pdbx_pdb_ins_code.get(element).map(String::as_str)
    }

    pub fn chain_key(&self, element: usize) -> Option<usize> {
        self.chain_index.get(element).copied()
    }

    pub fn chain_label_asym_id(&self, element: usize) -> Option<&str> {
        self.label_asym_id.get(element).map(String::as_str)
    }

    pub fn chain_auth_asym_id(&self, element: usize) -> Option<&str> {
        self.auth_asym_id.get(element).map(String::as_str)
    }

    pub fn chain_label_entity_id(&self, element: usize) -> Option<&str> {
        self.label_entity_id.get(element).map(String::as_str)
    }
}

#[derive(Clone, Debug, Default)]
pub struct AtomicRanges {
    pub polymer_ranges: Vec<usize>,
    pub gap_ranges: Vec<usize>,
    pub cyclic_polymer_map: BTreeMap<usize, usize>,
}

impl AtomicRanges {
    fn from_hierarchy(
        hierarchy: &AtomicHierarchy,
        sequence: &StructureSequence,
        xyz_defined: bool,
    ) -> Self {
        let mut polymer_ranges = Vec::new();
        let mut gap_ranges = Vec::new();
        let mut cyclic_polymer_map = BTreeMap::new();
        for (chain_index, chain) in hierarchy.chains.iter().enumerate() {
            let residues = chain_residue_indices(hierarchy, chain_index, chain);
            if residues.is_empty() {
                continue;
            }
            if let (Some(&first), Some(&last)) = (residues.first(), residues.last()) {
                let first_seq = residue_seq_id(&hierarchy.residues[first]);
                let last_seq = residue_seq_id(&hierarchy.residues[last]);
                let max_seq = sequence_max_seq_id_for_chain(hierarchy, sequence, chain_index);
                if first_seq == Some(1)
                    && last_seq == max_seq
                    && xyz_defined
                    && are_backbone_connected(hierarchy, first, last)
                {
                    cyclic_polymer_map.insert(first, last);
                    cyclic_polymer_map.insert(last, first);
                }
            }

            let mut start_index: Option<usize> = None;
            let mut start_residue: Option<usize> = None;
            let mut previous_residue: Option<usize> = None;
            let mut previous_start = 0usize;
            let mut previous_end = 0usize;
            let mut previous_seq_id: Option<i32> = None;
            for (residue_position, residue_index) in residues.iter().copied().enumerate() {
                let Some(residue) = hierarchy.residues.get(residue_index) else {
                    continue;
                };
                let seq_id = residue_seq_id(residue);
                if is_polymer_residue(hierarchy, residue_index)
                    && hierarchy
                        .derived
                        .residue
                        .trace_element_index
                        .get(residue_index)
                        .and_then(|index| *index)
                        .is_some()
                {
                    if let Some(start_element) = start_index {
                        let sequential =
                            matches!((previous_seq_id, seq_id), (Some(a), Some(b)) if b == a + 1);
                        if !sequential {
                            if start_residue.is_some() && previous_residue.is_some() {
                                push_atomic_range(
                                    &mut polymer_ranges,
                                    AtomicRange {
                                        start_element,
                                        end_element: previous_end,
                                    },
                                );
                                push_atomic_range(
                                    &mut gap_ranges,
                                    AtomicRange {
                                        start_element: previous_start,
                                        end_element: residue.end_atom,
                                    },
                                );
                            }
                            start_index = Some(residue.start_atom);
                            start_residue = Some(residue_index);
                        } else if residue_position + 1 == residues.len() {
                            if start_residue.is_some() {
                                push_atomic_range(
                                    &mut polymer_ranges,
                                    AtomicRange {
                                        start_element,
                                        end_element: residue.end_atom,
                                    },
                                );
                            }
                        } else if xyz_defined
                            && previous_residue.is_some_and(|previous| {
                                !are_backbone_connected(hierarchy, residue_index, previous)
                            })
                        {
                            if start_residue.is_some() && previous_residue.is_some() {
                                push_atomic_range(
                                    &mut polymer_ranges,
                                    AtomicRange {
                                        start_element,
                                        end_element: previous_end,
                                    },
                                );
                                push_atomic_range(
                                    &mut gap_ranges,
                                    AtomicRange {
                                        start_element: previous_start,
                                        end_element: residue.end_atom,
                                    },
                                );
                            }
                            start_index = Some(residue.start_atom);
                            start_residue = Some(residue_index);
                        }
                    } else {
                        start_index = Some(residue.start_atom);
                        start_residue = Some(residue_index);
                    }
                } else if let (Some(start_element), Some(_start), Some(_previous)) =
                    (start_index.take(), start_residue.take(), previous_residue)
                {
                    push_atomic_range(
                        &mut polymer_ranges,
                        AtomicRange {
                            start_element,
                            end_element: previous_end,
                        },
                    );
                }

                previous_start = residue.start_atom;
                previous_end = residue.end_atom;
                previous_seq_id = seq_id;
                previous_residue = Some(residue_index);
            }
        }
        AtomicRanges {
            polymer_ranges,
            gap_ranges,
            cyclic_polymer_map,
        }
    }
}

fn push_atomic_range(ranges: &mut Vec<usize>, range: AtomicRange) {
    if range.start_element < range.end_element {
        ranges.push(range.start_element);
        ranges.push(range.end_element - 1);
    }
}

fn sequence_max_seq_id_for_chain(
    hierarchy: &AtomicHierarchy,
    sequence: &StructureSequence,
    chain_index: usize,
) -> Option<i32> {
    let entity_key = hierarchy.index.entity_from_chain(chain_index)?;
    let sequence_index = sequence.by_entity_key.get(&entity_key).copied()?;
    let entity = sequence.sequences.get(sequence_index)?;
    if let Some(residue) = entity.residues.last() {
        return Some(residue.seq_id);
    }
    entity
        .ranges
        .iter()
        .map(|range| range.seq_id_end)
        .max()
        .map(|max_end| max_end + 1)
}

#[derive(Clone, Debug)]
struct AtomicRange {
    start_element: usize,
    end_element: usize,
}

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum UnitKind {
    Atomic,
    Spheres,
    Gaussians,
}

pub type AtomicUnitKind = UnitKind;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UnitTraits {
    pub bits: u32,
}

impl UnitTraits {
    pub const NONE: UnitTraits = UnitTraits { bits: 0 };
    pub const MULTI_CHAIN: UnitTraits = UnitTraits { bits: 0x1 };
    pub const PARTITIONED: UnitTraits = UnitTraits { bits: 0x2 };
    pub const FAST_BOUNDARY: UnitTraits = UnitTraits { bits: 0x4 };
    pub const WATER: UnitTraits = UnitTraits { bits: 0x8 };
    pub const COARSE_GRAINED: UnitTraits = UnitTraits { bits: 0x10 };

    pub const fn contains(self, other: UnitTraits) -> bool {
        self.bits & other.bits == other.bits
    }

    pub const fn union(self, other: UnitTraits) -> UnitTraits {
        UnitTraits {
            bits: self.bits | other.bits,
        }
    }
}

#[derive(Clone, Debug)]
pub struct UnitOperator {
    pub name: String,
    pub instance_id: String,
    pub assembly_id: String,
    pub oper_id: i32,
    pub oper_list_ids: Vec<String>,
    pub transform: Transform,
    pub is_identity: bool,
    pub suffix: String,
}

impl Default for UnitOperator {
    fn default() -> Self {
        UnitOperator {
            name: "1_555".to_string(),
            instance_id: "1_555".to_string(),
            assembly_id: String::new(),
            oper_id: -1,
            oper_list_ids: Vec::new(),
            transform: Transform::identity(),
            is_identity: true,
            suffix: String::new(),
        }
    }
}

pub type Operator = UnitOperator;

#[derive(Clone, Debug)]
pub struct StructureUnit {
    pub id: usize,
    pub invariant_id: usize,
    pub chain_group_id: usize,
    pub kind: AtomicUnitKind,
    pub traits: UnitTraits,
    pub model_index: usize,
    pub chain_index: usize,
    pub chain_indices: Vec<usize>,
    pub elements: Vec<usize>,
    pub atom_indices: Vec<usize>,
    pub residue_indices: Vec<usize>,
    pub residue_index_by_element: Vec<usize>,
    pub chain_index_by_element: Vec<usize>,
    pub props: UnitProps,
    pub operator: Operator,
}

impl StructureUnit {
    pub fn are_same_chain_operator_group(&self, other: &StructureUnit) -> bool {
        self.chain_group_id == other.chain_group_id && self.operator.name == other.operator.name
    }

    fn chain_index_at(&self, element_index: usize) -> Option<usize> {
        match self.kind {
            UnitKind::Atomic => self
                .elements
                .get(element_index)
                .and_then(|element| self.chain_index_by_element.get(*element))
                .copied(),
            UnitKind::Spheres | UnitKind::Gaussians => {
                self.chain_index_by_element.get(element_index).copied()
            }
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct UnitProps {
    pub residue_count: usize,
    pub protein_elements: Vec<usize>,
    pub nucleotide_elements: Vec<usize>,
    pub polymer_elements: Vec<usize>,
    pub gap_elements: Vec<usize>,
    pub water_elements: Vec<usize>,
    pub boundary: Boundary,
    pub principal_axes: PrincipalAxes,
    pub lookup3d: UnitLookup3D,
    pub intra_unit_bonds: IntraUnitBonds,
    pub intra_unit_bond_count: usize,
    pub inter_unit_bond_count: usize,
}

impl UnitProps {
    fn from_elements(
        hierarchy: &AtomicHierarchy,
        ranges: &AtomicRanges,
        elements: &[usize],
        residue_indices: &[usize],
        operator: &UnitOperator,
    ) -> Self {
        let mut protein_elements = Vec::new();
        let mut nucleotide_elements = Vec::new();
        let mut water_elements = Vec::new();
        let mut positions = Vec::new();
        let mut invariant_positions = Vec::new();
        let mut radii = Vec::new();
        for &element in elements {
            if let Some(atom) = hierarchy.atoms.get(element) {
                invariant_positions.push(atom.position);
                positions.push(operator.transform.apply(atom.position));
                radii.push(molstar_vdw_radius(&atom.type_symbol));
            }
        }
        let boundary = Boundary::from_positions_and_radii(&positions, &radii);
        let principal_axes = PrincipalAxes::of_positions(&invariant_positions);
        let lookup3d = UnitLookup3D::new(positions, boundary.clone());
        let polymer_elements = atomic_polymer_elements(hierarchy, ranges, elements);
        let gap_elements = atomic_gap_elements(hierarchy, ranges, elements);
        for &residue_index in residue_indices {
            let Some(residue) = hierarchy.residues.get(residue_index) else {
                continue;
            };
            if !residue_intersects_elements(residue, elements) {
                continue;
            }
            let element = residue_trace_or_first_element(hierarchy, residue_index);
            let molecule_type = hierarchy
                .derived
                .residue
                .molecule_type
                .get(residue_index)
                .copied()
                .unwrap_or_default();
            if matches!(molecule_type, MoleculeType::Protein) {
                protein_elements.push(element);
            }
            if matches!(
                molecule_type,
                MoleculeType::Rna | MoleculeType::Dna | MoleculeType::Pna
            ) {
                nucleotide_elements.push(element);
            }
        }
        for &element in elements {
            if hierarchy
                .derived
                .atom
                .is_water
                .get(element)
                .copied()
                .unwrap_or(false)
            {
                water_elements.push(element);
            }
        }
        UnitProps {
            residue_count: residue_indices.len(),
            protein_elements,
            nucleotide_elements,
            polymer_elements,
            gap_elements,
            water_elements,
            boundary,
            principal_axes,
            lookup3d,
            intra_unit_bonds: IntraUnitBonds::default(),
            intra_unit_bond_count: 0,
            inter_unit_bond_count: 0,
        }
    }

    fn from_coarse_positions(
        positions: Vec<Vec3>,
        invariant_positions: Vec<Vec3>,
        radii: Vec<f32>,
        elements: &[usize],
        polymer_ranges: &[coarse::CoarseRange],
        gap_ranges: &[coarse::CoarseRange],
    ) -> Self {
        let boundary = Boundary::from_positions_and_radii(&positions, &radii);
        let principal_axes = PrincipalAxes::of_positions(&invariant_positions);
        let residue_count = positions.len();
        let lookup3d = UnitLookup3D::new(positions, boundary.clone());
        UnitProps {
            residue_count,
            protein_elements: Vec::new(),
            nucleotide_elements: Vec::new(),
            polymer_elements: coarse_polymer_elements(elements, polymer_ranges),
            gap_elements: coarse_gap_elements(elements, gap_ranges),
            water_elements: Vec::new(),
            boundary,
            principal_axes,
            lookup3d,
            intra_unit_bonds: IntraUnitBonds::default(),
            intra_unit_bond_count: 0,
            inter_unit_bond_count: 0,
        }
    }
}

fn residue_trace_or_first_element(hierarchy: &AtomicHierarchy, residue_index: usize) -> usize {
    hierarchy
        .derived
        .residue
        .trace_element_index
        .get(residue_index)
        .and_then(|index| *index)
        .or_else(|| {
            hierarchy
                .residues
                .get(residue_index)
                .map(|residue| residue.start_atom)
        })
        .unwrap_or(0)
}

fn residue_intersects_elements(residue: &AtomicResidue, elements: &[usize]) -> bool {
    elements
        .iter()
        .any(|element| (residue.start_atom..residue.end_atom).contains(element))
}

fn atomic_range_pair_intersects_elements(
    start: usize,
    end_inclusive: usize,
    elements: &[usize],
) -> bool {
    elements
        .iter()
        .any(|element| (start..=end_inclusive).contains(element))
}

fn atomic_polymer_elements(
    hierarchy: &AtomicHierarchy,
    ranges: &AtomicRanges,
    elements: &[usize],
) -> Vec<usize> {
    let mut indices = Vec::new();
    for pair in ranges.polymer_ranges.chunks_exact(2) {
        let start_element = pair[0];
        let end_element = pair[1];
        if !atomic_range_pair_intersects_elements(start_element, end_element, elements) {
            continue;
        }
        let start_residue = hierarchy
            .residue_atom_segments
            .index
            .get(start_element)
            .copied()
            .unwrap_or(0);
        let end_residue = hierarchy
            .residue_atom_segments
            .index
            .get(end_element)
            .copied()
            .unwrap_or(start_residue);
        for residue_index in start_residue..=end_residue {
            let Some(residue) = hierarchy.residues.get(residue_index) else {
                continue;
            };
            if residue_intersects_elements(residue, elements) {
                indices.push(residue_trace_or_first_element(hierarchy, residue_index));
            }
        }
    }
    indices
}

fn atomic_gap_elements(
    hierarchy: &AtomicHierarchy,
    ranges: &AtomicRanges,
    elements: &[usize],
) -> Vec<usize> {
    let mut indices = Vec::new();
    for pair in ranges.gap_ranges.chunks_exact(2) {
        let start_element = pair[0];
        let end_element = pair[1];
        if !atomic_range_pair_intersects_elements(start_element, end_element, elements) {
            continue;
        }
        let Some(first_element) = elements
            .iter()
            .copied()
            .find(|element| (start_element..=end_element).contains(element))
        else {
            continue;
        };
        let Some(last_element) = elements
            .iter()
            .rev()
            .copied()
            .find(|element| (start_element..=end_element).contains(element))
        else {
            continue;
        };
        let start_residue = hierarchy
            .residue_atom_segments
            .index
            .get(first_element)
            .copied()
            .unwrap_or(0);
        let end_residue = hierarchy
            .residue_atom_segments
            .index
            .get(last_element)
            .copied()
            .unwrap_or(start_residue);
        if hierarchy.residues.get(start_residue).is_some()
            && hierarchy.residues.get(end_residue).is_some()
        {
            indices.push(residue_trace_or_first_element(hierarchy, start_residue));
            indices.push(residue_trace_or_first_element(hierarchy, end_residue));
        }
    }
    indices
}

fn coarse_polymer_elements(elements: &[usize], ranges: &[coarse::CoarseRange]) -> Vec<usize> {
    let mut indices = Vec::new();
    for range in ranges {
        for element in range.start_element..=range.end_element {
            if elements.contains(&element) {
                indices.push(element);
            }
        }
    }
    indices
}

fn coarse_gap_elements(elements: &[usize], ranges: &[coarse::CoarseRange]) -> Vec<usize> {
    let mut indices = Vec::new();
    for range in ranges {
        if elements.contains(&range.start_element) && elements.contains(&range.end_element) {
            indices.push(range.start_element);
            indices.push(range.end_element);
        }
    }
    indices
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Boundary {
    pub box_min: Vec3,
    pub box_max: Vec3,
    pub sphere: BoundingSphere,
}

impl Boundary {
    pub fn from_positions(positions: &[Vec3]) -> Self {
        Boundary::from_positions_and_optional_radii(positions, None)
    }

    pub(crate) fn from_positions_and_radii(positions: &[Vec3], radii: &[f32]) -> Self {
        Boundary::from_positions_and_optional_radii(positions, Some(radii))
    }

    fn from_positions_and_optional_radii(positions: &[Vec3], radii: Option<&[f32]>) -> Self {
        if positions.is_empty() {
            return Boundary::default();
        }

        if positions.len() > 250_000 {
            return Boundary::fast_from_positions(positions);
        }

        let mut helper = BoundaryHelper::new(if positions.len() > 10_000 {
            EposQuality::Coarse
        } else {
            EposQuality::Fine
        });
        for (i, &position) in positions.iter().enumerate() {
            helper.include_position_radius(
                position,
                radii.and_then(|r| r.get(i)).copied().unwrap_or(0.0),
            );
        }
        helper.finished_include_step();
        for (i, &position) in positions.iter().enumerate() {
            helper.radius_position_radius(
                position,
                radii.and_then(|r| r.get(i)).copied().unwrap_or(0.0),
            );
        }
        let mut boundary = helper.boundary();

        if radii.is_none() && positions.len() <= boundary.sphere.extrema.len() {
            boundary.sphere.extrema = positions.to_vec();
            boundary.sphere.extrema64 = positions
                .iter()
                .map(|point| [point.x as f64, point.y as f64, point.z as f64])
                .collect();
        }
        boundary
    }

    fn fast_from_positions(positions: &[Vec3]) -> Self {
        let mut min = positions[0];
        let mut max = positions[0];
        for &position in &positions[1..] {
            min = min.min(position);
            max = max.max(position);
        }
        let center = (min + max) * 0.5;
        let radius = center.distance(max);
        let extrema = box_corners(min, max);
        Boundary {
            box_min: min,
            box_max: max,
            sphere: BoundingSphere {
                center,
                radius,
                extrema,
                center64: Some([center.x as f64, center.y as f64, center.z as f64]),
                radius64: Some(radius as f64),
                extrema64: box_corners(min, max)
                    .into_iter()
                    .map(|point| [point.x as f64, point.y as f64, point.z as f64])
                    .collect(),
            },
        }
    }

    #[cfg(test)]
    pub(crate) fn union(self, other: Boundary) -> Boundary {
        Boundary::from_spheres(&[self.sphere, other.sphere], EposQuality::Fine)
    }

    pub(crate) fn from_bounding_spheres(spheres: &[BoundingSphere]) -> Boundary {
        Boundary::from_spheres(spheres, EposQuality::Fine)
    }

    fn from_spheres(spheres: &[BoundingSphere], quality: EposQuality) -> Boundary {
        if spheres.is_empty() {
            return Boundary::default();
        }
        let mut helper = BoundaryHelper::new(quality);
        for sphere in spheres {
            helper.include_sphere(sphere);
        }
        helper.finished_include_step();
        for sphere in spheres {
            helper.radius_sphere(sphere);
        }
        helper.boundary()
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BoundingSphere {
    pub center: Vec3,
    pub radius: f32,
    pub extrema: Vec<Vec3>,
    pub(crate) center64: Option<[f64; 3]>,
    pub(crate) radius64: Option<f64>,
    pub(crate) extrema64: Vec<[f64; 3]>,
}

impl BoundingSphere {
    pub(crate) fn center64(&self) -> [f64; 3] {
        self.center64.unwrap_or([
            self.center.x as f64,
            self.center.y as f64,
            self.center.z as f64,
        ])
    }

    pub(crate) fn radius64(&self) -> f64 {
        self.radius64.unwrap_or(self.radius as f64)
    }
}

#[derive(Clone, Copy)]
enum EposQuality {
    Coarse,
    Fine,
}

struct BoundaryHelper {
    directions: Vec<Vec3d>,
    min_dist: Vec<f64>,
    max_dist: Vec<f64>,
    extrema: Vec<Vec3d>,
    center: Vec3d,
    radius_sq: f64,
    count: usize,
}

impl BoundaryHelper {
    fn new(quality: EposQuality) -> Self {
        let directions = epos_directions(quality);
        let len = directions.len();
        BoundaryHelper {
            directions,
            min_dist: vec![f64::INFINITY; len],
            max_dist: vec![f64::NEG_INFINITY; len],
            extrema: vec![Vec3d::default(); len * 2],
            center: Vec3d::default(),
            radius_sq: 0.0,
            count: 0,
        }
    }

    fn include_sphere(&mut self, sphere: &BoundingSphere) {
        if sphere.extrema64.len() > 1 {
            for &point in &sphere.extrema64 {
                self.include_position64(point);
            }
        } else if sphere.extrema.len() > 1 {
            for &point in &sphere.extrema {
                self.include_position(point);
            }
        } else {
            self.include_position_radius(sphere.center, sphere.radius);
        }
    }

    fn include_position(&mut self, point: Vec3) {
        let point = Vec3d::from(point);
        self.include_position_vec3d(point);
    }

    fn include_position64(&mut self, point: [f64; 3]) {
        self.include_position_vec3d(Vec3d::from(point));
    }

    fn include_position_vec3d(&mut self, point: Vec3d) {
        for i in 0..self.directions.len() {
            self.compute_extrema(i, point);
        }
    }

    fn include_position_radius(&mut self, center: Vec3, radius: f32) {
        let center = Vec3d::from(center);
        let radius = radius as f64;
        for i in 0..self.directions.len() {
            self.compute_sphere_extrema(i, center, radius);
        }
    }

    fn compute_extrema(&mut self, i: usize, point: Vec3d) {
        let distance = self.directions[i].dot(point);
        if distance < self.min_dist[i] {
            self.min_dist[i] = distance;
            self.extrema[i * 2] = point;
        }
        if distance > self.max_dist[i] {
            self.max_dist[i] = distance;
            self.extrema[i * 2 + 1] = point;
        }
    }

    fn compute_sphere_extrema(&mut self, i: usize, center: Vec3d, radius: f64) {
        let direction = self.directions[i];
        let distance = direction.dot(center);
        if distance - radius < self.min_dist[i] {
            self.min_dist[i] = distance - radius;
            self.extrema[i * 2] = center - direction * radius;
        }
        if distance + radius > self.max_dist[i] {
            self.max_dist[i] = distance + radius;
            self.extrema[i * 2 + 1] = center + direction * radius;
        }
    }

    fn finished_include_step(&mut self) {
        for point in &self.extrema {
            self.center = self.center + *point;
            self.count += 1;
        }
        if self.count > 0 {
            self.center = self.center / self.count as f64;
        }
    }

    fn radius_sphere(&mut self, sphere: &BoundingSphere) {
        if sphere.extrema64.len() > 1 {
            for &point in &sphere.extrema64 {
                self.radius_position64(point);
            }
        } else if sphere.extrema.len() > 1 {
            for &point in &sphere.extrema {
                self.radius_position(point);
            }
        } else {
            self.radius_position_radius(sphere.center, sphere.radius);
        }
    }

    fn radius_position(&mut self, point: Vec3) {
        self.radius_position_vec3d(Vec3d::from(point));
    }

    fn radius_position64(&mut self, point: [f64; 3]) {
        self.radius_position_vec3d(Vec3d::from(point));
    }

    fn radius_position_vec3d(&mut self, point: Vec3d) {
        let distance = point.squared_distance(self.center);
        if distance > self.radius_sq {
            self.radius_sq = distance;
        }
    }

    fn radius_position_radius(&mut self, center: Vec3, radius: f32) {
        let distance = Vec3d::from(center).distance(self.center) + radius as f64;
        let distance_sq = distance * distance;
        if distance_sq > self.radius_sq {
            self.radius_sq = distance_sq;
        }
    }

    fn boundary(self) -> Boundary {
        let extrema64 = self
            .extrema
            .iter()
            .map(|point| [point.x, point.y, point.z])
            .collect::<Vec<_>>();
        let extrema = self
            .extrema
            .iter()
            .map(|point| point.to_vec3())
            .collect::<Vec<_>>();
        let (box_min, box_max) = box_from_points(&extrema);
        Boundary {
            box_min,
            box_max,
            sphere: BoundingSphere {
                center: self.center.to_vec3(),
                radius: self.radius_sq.sqrt() as f32,
                extrema,
                center64: Some([self.center.x, self.center.y, self.center.z]),
                radius64: Some(self.radius_sq.sqrt()),
                extrema64,
            },
        }
    }
}

#[derive(Clone, Copy, Debug, Default)]
struct Vec3d {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3d {
    fn dot(self, other: Vec3d) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn length(self) -> f64 {
        self.dot(self).sqrt()
    }

    fn normalized(self) -> Vec3d {
        let length = self.length();
        if length > 0.0 {
            self / length
        } else {
            self
        }
    }

    fn squared_distance(self, other: Vec3d) -> f64 {
        let d = self - other;
        d.dot(d)
    }

    fn distance(self, other: Vec3d) -> f64 {
        self.squared_distance(other).sqrt()
    }

    fn to_vec3(self) -> Vec3 {
        Vec3::new(self.x as f32, self.y as f32, self.z as f32)
    }
}

impl From<Vec3> for Vec3d {
    fn from(value: Vec3) -> Self {
        Vec3d {
            x: value.x as f64,
            y: value.y as f64,
            z: value.z as f64,
        }
    }
}

impl From<[f64; 3]> for Vec3d {
    fn from(value: [f64; 3]) -> Self {
        Vec3d {
            x: value[0],
            y: value[1],
            z: value[2],
        }
    }
}

impl std::ops::Add for Vec3d {
    type Output = Vec3d;
    fn add(self, rhs: Vec3d) -> Vec3d {
        Vec3d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl std::ops::Sub for Vec3d {
    type Output = Vec3d;
    fn sub(self, rhs: Vec3d) -> Vec3d {
        Vec3d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl std::ops::Mul<f64> for Vec3d {
    type Output = Vec3d;
    fn mul(self, rhs: f64) -> Vec3d {
        Vec3d {
            x: self.x * rhs,
            y: self.y * rhs,
            z: self.z * rhs,
        }
    }
}

impl std::ops::Div<f64> for Vec3d {
    type Output = Vec3d;
    fn div(self, rhs: f64) -> Vec3d {
        Vec3d {
            x: self.x / rhs,
            y: self.y / rhs,
            z: self.z / rhs,
        }
    }
}

fn box_from_points(points: &[Vec3]) -> (Vec3, Vec3) {
    let Some(first) = points.first().copied() else {
        return (Vec3::default(), Vec3::default());
    };
    let mut min = first;
    let mut max = first;
    for &point in &points[1..] {
        min = min.min(point);
        max = max.max(point);
    }
    (min, max)
}

fn box_corners(min: Vec3, max: Vec3) -> Vec<Vec3> {
    vec![
        Vec3::new(min.x, min.y, min.z),
        Vec3::new(max.x, max.y, max.z),
        Vec3::new(max.x, min.y, min.z),
        Vec3::new(min.x, max.y, max.z),
        Vec3::new(min.x, min.y, max.z),
        Vec3::new(max.x, min.y, max.z),
        Vec3::new(max.x, max.y, min.z),
        Vec3::new(min.x, max.y, min.z),
    ]
}

fn epos_directions(quality: EposQuality) -> Vec<Vec3d> {
    let mut directions = Vec::new();
    directions.extend(TYPE_001);
    directions.extend(TYPE_111);
    if matches!(quality, EposQuality::Fine) {
        directions.extend(TYPE_011);
        directions.extend(TYPE_012);
        directions.extend(TYPE_112);
        directions.extend(TYPE_122);
    }
    directions
        .into_iter()
        .map(|[x, y, z]| Vec3d { x, y, z }.normalized())
        .collect()
}

const TYPE_001: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

const TYPE_111: [[f64; 3]; 4] = [
    [1.0, 1.0, 1.0],
    [-1.0, 1.0, 1.0],
    [-1.0, -1.0, 1.0],
    [1.0, -1.0, 1.0],
];

const TYPE_011: [[f64; 3]; 6] = [
    [1.0, 1.0, 0.0],
    [1.0, -1.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 0.0, -1.0],
    [0.0, 1.0, 1.0],
    [0.0, 1.0, -1.0],
];

const TYPE_012: [[f64; 3]; 12] = [
    [0.0, 1.0, 2.0],
    [0.0, 2.0, 1.0],
    [1.0, 0.0, 2.0],
    [2.0, 0.0, 1.0],
    [1.0, 2.0, 0.0],
    [2.0, 1.0, 0.0],
    [0.0, 1.0, -2.0],
    [0.0, 2.0, -1.0],
    [1.0, 0.0, -2.0],
    [2.0, 0.0, -1.0],
    [1.0, -2.0, 0.0],
    [2.0, -1.0, 0.0],
];

const TYPE_112: [[f64; 3]; 12] = [
    [1.0, 1.0, 2.0],
    [2.0, 1.0, 1.0],
    [1.0, 2.0, 1.0],
    [1.0, -1.0, 2.0],
    [1.0, 1.0, -2.0],
    [1.0, -1.0, -2.0],
    [2.0, -1.0, 1.0],
    [2.0, 1.0, -1.0],
    [2.0, -1.0, -1.0],
    [1.0, -2.0, 1.0],
    [1.0, 2.0, -1.0],
    [1.0, -2.0, -1.0],
];

const TYPE_122: [[f64; 3]; 12] = [
    [2.0, 2.0, 1.0],
    [1.0, 2.0, 2.0],
    [2.0, 1.0, 2.0],
    [2.0, -2.0, 1.0],
    [2.0, 2.0, -1.0],
    [2.0, -2.0, -1.0],
    [1.0, -2.0, 2.0],
    [1.0, 2.0, -2.0],
    [1.0, -2.0, -2.0],
    [2.0, -1.0, 2.0],
    [2.0, 1.0, -2.0],
    [2.0, -1.0, -2.0],
];

fn molstar_vdw_radius(element: &str) -> f32 {
    match crate::chemistry::atomic_number(element) {
        1 => 1.10,
        2 => 1.40,
        3 => 1.81,
        4 => 1.53,
        5 => 1.92,
        6 => 1.70,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        10 => 1.54,
        11 => 2.27,
        12 => 1.73,
        13 => 1.84,
        14 => 2.10,
        15 => 1.80,
        16 => 1.80,
        17 => 1.75,
        18 => 1.88,
        19 => 2.75,
        20 => 2.31,
        21 => 2.30,
        22 => 2.15,
        23 => 2.05,
        24 => 2.05,
        25 => 2.05,
        26 => 2.05,
        27 => 2.00,
        28 => 2.00,
        29 => 2.00,
        30 => 2.10,
        31 => 1.87,
        32 => 2.11,
        33 => 1.85,
        34 => 1.90,
        35 => 1.83,
        36 => 2.02,
        37 => 3.03,
        38 => 2.49,
        39 => 2.40,
        40 => 2.30,
        41 => 2.15,
        42 => 2.10,
        43 => 2.05,
        44 => 2.05,
        45 => 2.00,
        46 => 2.05,
        47 => 2.10,
        48 => 2.20,
        49 => 2.20,
        50 => 1.93,
        51 => 2.17,
        52 => 2.06,
        53 => 1.98,
        54 => 2.16,
        55 => 3.43,
        56 => 2.68,
        57 => 2.50,
        58 => 2.48,
        59 => 2.47,
        60 => 2.45,
        61 => 2.43,
        62 => 2.42,
        63 => 2.40,
        64 => 2.38,
        65 => 2.37,
        66 => 2.35,
        67 => 2.33,
        68 => 2.32,
        69 => 2.30,
        70 => 2.28,
        71 => 2.27,
        72 => 2.25,
        73 => 2.20,
        74 => 2.10,
        75 => 2.05,
        76 => 2.00,
        77 => 2.00,
        78 => 2.05,
        79 => 2.10,
        80 => 2.05,
        81 => 1.96,
        82 => 2.02,
        83 => 2.07,
        84 => 1.97,
        85 => 2.02,
        86 => 2.20,
        87 => 3.48,
        88 => 2.83,
        89..=109 => 2.00,
        _ => 1.70,
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LookupHit {
    pub index: usize,
    pub squared_distance: f32,
}

#[derive(Clone, Debug, Default)]
pub struct UnitLookup3D {
    positions: Vec<Vec3>,
    pub boundary: Boundary,
    grid: UnitLookupGrid,
}

impl UnitLookup3D {
    pub(crate) fn new(positions: Vec<Vec3>, boundary: Boundary) -> Self {
        let grid = UnitLookupGrid::new(&positions, boundary.clone());
        UnitLookup3D {
            positions,
            boundary,
            grid,
        }
    }

    pub fn find(&self, point: Vec3, radius: f32) -> Vec<LookupHit> {
        self.grid.find(&self.positions, point, radius, false).0
    }

    pub fn check(&self, point: Vec3, radius: f32) -> bool {
        self.grid.find(&self.positions, point, radius, true).1
    }

    pub fn nearest(&self, point: Vec3, k: usize) -> Vec<LookupHit> {
        let mut hits = self
            .positions
            .iter()
            .enumerate()
            .map(|(index, position)| LookupHit {
                index,
                squared_distance: position.squared_distance(point),
            })
            .collect::<Vec<_>>();
        hits.sort_by(|a, b| a.squared_distance.total_cmp(&b.squared_distance));
        hits.truncate(k);
        hits
    }

    pub fn len(&self) -> usize {
        self.positions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    #[cfg(test)]
    pub(crate) fn grid_size(&self) -> [usize; 3] {
        self.grid.size
    }

    #[cfg(test)]
    pub(crate) fn grid_delta(&self) -> Vec3 {
        self.grid.delta
    }

    #[cfg(test)]
    pub(crate) fn bucket_offsets(&self) -> &[usize] {
        &self.grid.bucket_offsets
    }

    #[cfg(test)]
    pub(crate) fn bucket_counts(&self) -> &[usize] {
        &self.grid.bucket_counts
    }

    #[cfg(test)]
    pub(crate) fn bucket_array(&self) -> &[usize] {
        &self.grid.bucket_array
    }
}

#[derive(Clone, Debug, Default)]
struct UnitLookupGrid {
    size: [usize; 3],
    min: Vec3,
    delta: Vec3,
    grid: Vec<u32>,
    bucket_offsets: Vec<usize>,
    bucket_counts: Vec<usize>,
    bucket_array: Vec<usize>,
}

impl UnitLookupGrid {
    const DEFAULT_CELL_COUNT: f32 = 32.0;
    const MAX_VOLUME: usize = 1 << 24;

    fn new(positions: &[Vec3], boundary: Boundary) -> Self {
        let expanded_min = boundary.box_min - Vec3::new(0.5, 0.5, 0.5);
        let expanded_max = boundary.box_max + Vec3::new(0.5, 0.5, 0.5);
        let size_vec = expanded_max - expanded_min;
        let element_count = positions.len();

        let (size, delta) = if element_count > 0 {
            let required_volume = (element_count as f32 / Self::DEFAULT_CELL_COUNT).ceil();
            let box_volume = size_vec.x * size_vec.y * size_vec.z;
            let factor = (required_volume / box_volume).powf(1.0 / 3.0);
            let mut size = [
                ceil_positive(size_vec.x * factor),
                ceil_positive(size_vec.y * factor),
                ceil_positive(size_vec.z * factor),
            ];
            let mut delta = Vec3::new(
                size_vec.x / size[0] as f32,
                size_vec.y / size[1] as f32,
                size_vec.z / size[2] as f32,
            );
            let volume = size[0] * size[1] * size[2];
            if volume > Self::MAX_VOLUME {
                let factor = (volume as f32 / Self::MAX_VOLUME as f32).cbrt();
                size = [
                    ceil_positive(size[0] as f32 / factor),
                    ceil_positive(size[1] as f32 / factor),
                    ceil_positive(size[2] as f32 / factor),
                ];
                delta = Vec3::new(
                    size_vec.x / size[0] as f32,
                    size_vec.y / size[1] as f32,
                    size_vec.z / size[2] as f32,
                );
            }
            (size, delta)
        } else {
            ([1, 1, 1], size_vec)
        };

        let grid_volume = size[0] * size[1] * size[2];
        let mut grid = vec![0u32; grid_volume];
        let mut bucket_index = vec![0usize; element_count];
        let mut bucket_count = 0usize;

        for (position_index, position) in positions.iter().enumerate() {
            let x = grid_axis(position.x, expanded_min.x, delta.x, size[0]);
            let y = grid_axis(position.y, expanded_min.y, delta.y, size[1]);
            let z = grid_axis(position.z, expanded_min.z, delta.z, size[2]);
            let index = grid_index(x, y, z, size);
            grid[index] += 1;
            if grid[index] == 1 {
                bucket_count += 1;
            }
            bucket_index[position_index] = index;
        }

        let mut bucket_counts = vec![0usize; bucket_count];
        let mut next_bucket = 0usize;
        for cell in &mut grid {
            let count = *cell;
            if count > 0 {
                *cell = (next_bucket + 1) as u32;
                bucket_counts[next_bucket] = count as usize;
                next_bucket += 1;
            }
        }

        let mut bucket_offsets = vec![0usize; bucket_count];
        for i in 1..bucket_count {
            bucket_offsets[i] = bucket_offsets[i - 1] + bucket_counts[i - 1];
        }

        let mut bucket_fill = vec![0usize; bucket_count];
        let mut bucket_array = vec![0usize; element_count];
        for (position_index, index) in bucket_index.iter().copied().enumerate() {
            let bucket = grid[index];
            if bucket == 0 {
                continue;
            }
            let bucket = bucket as usize - 1;
            bucket_array[bucket_offsets[bucket] + bucket_fill[bucket]] = position_index;
            bucket_fill[bucket] += 1;
        }

        UnitLookupGrid {
            size,
            min: expanded_min,
            delta,
            grid,
            bucket_offsets,
            bucket_counts,
            bucket_array,
        }
    }

    fn find(
        &self,
        positions: &[Vec3],
        point: Vec3,
        radius: f32,
        is_check: bool,
    ) -> (Vec<LookupHit>, bool) {
        let mut hits = Vec::new();
        if positions.is_empty() {
            return (hits, false);
        }
        let radius_sq = radius * radius;

        let lo_x = query_min_axis(point.x, radius, self.min.x, self.delta.x, self.size[0]);
        let lo_y = query_min_axis(point.y, radius, self.min.y, self.delta.y, self.size[1]);
        let lo_z = query_min_axis(point.z, radius, self.min.z, self.delta.z, self.size[2]);
        let hi_x = query_max_axis(point.x, radius, self.min.x, self.delta.x, self.size[0]);
        let hi_y = query_max_axis(point.y, radius, self.min.y, self.delta.y, self.size[1]);
        let hi_z = query_max_axis(point.z, radius, self.min.z, self.delta.z, self.size[2]);

        let (Some(lo_x), Some(lo_y), Some(lo_z), Some(hi_x), Some(hi_y), Some(hi_z)) =
            (lo_x, lo_y, lo_z, hi_x, hi_y, hi_z)
        else {
            return (hits, false);
        };

        if lo_x > hi_x || lo_y > hi_y || lo_z > hi_z {
            return (hits, false);
        }

        for x in lo_x..=hi_x {
            for y in lo_y..=hi_y {
                for z in lo_z..=hi_z {
                    let bucket = self.grid[grid_index(x, y, z, self.size)];
                    if bucket == 0 {
                        continue;
                    }
                    let bucket = bucket as usize - 1;
                    let offset = self.bucket_offsets[bucket];
                    let end = offset + self.bucket_counts[bucket];
                    for i in offset..end {
                        let index = self.bucket_array[i];
                        let squared_distance = positions[index].squared_distance(point);
                        if squared_distance <= radius_sq {
                            if is_check {
                                return (hits, true);
                            }
                            hits.push(LookupHit {
                                index,
                                squared_distance,
                            });
                        }
                    }
                }
            }
        }

        let found = !hits.is_empty();
        (hits, found)
    }

    fn find_with_radii(
        &self,
        positions: &[Vec3],
        radii: &[f32],
        max_radius: f32,
        point: Vec3,
        radius: f32,
        is_check: bool,
    ) -> (Vec<LookupHit>, bool) {
        let mut hits = Vec::new();
        if positions.is_empty() {
            return (hits, false);
        }
        let query_radius = radius + max_radius;
        let radius_sq = query_radius * query_radius;

        let lo_x = query_min_axis(
            point.x,
            query_radius,
            self.min.x,
            self.delta.x,
            self.size[0],
        );
        let lo_y = query_min_axis(
            point.y,
            query_radius,
            self.min.y,
            self.delta.y,
            self.size[1],
        );
        let lo_z = query_min_axis(
            point.z,
            query_radius,
            self.min.z,
            self.delta.z,
            self.size[2],
        );
        let hi_x = query_max_axis(
            point.x,
            query_radius,
            self.min.x,
            self.delta.x,
            self.size[0],
        );
        let hi_y = query_max_axis(
            point.y,
            query_radius,
            self.min.y,
            self.delta.y,
            self.size[1],
        );
        let hi_z = query_max_axis(
            point.z,
            query_radius,
            self.min.z,
            self.delta.z,
            self.size[2],
        );

        let (Some(lo_x), Some(lo_y), Some(lo_z), Some(hi_x), Some(hi_y), Some(hi_z)) =
            (lo_x, lo_y, lo_z, hi_x, hi_y, hi_z)
        else {
            return (hits, false);
        };

        if lo_x > hi_x || lo_y > hi_y || lo_z > hi_z {
            return (hits, false);
        }

        for x in lo_x..=hi_x {
            for y in lo_y..=hi_y {
                for z in lo_z..=hi_z {
                    let bucket = self.grid[grid_index(x, y, z, self.size)];
                    if bucket == 0 {
                        continue;
                    }
                    let bucket = bucket as usize - 1;
                    let offset = self.bucket_offsets[bucket];
                    let end = offset + self.bucket_counts[bucket];
                    for i in offset..end {
                        let index = self.bucket_array[i];
                        let squared_distance = positions[index].squared_distance(point);
                        if squared_distance <= radius_sq {
                            if max_radius > 0.0 && squared_distance.sqrt() - radii[index] > radius {
                                continue;
                            }
                            if is_check {
                                return (hits, true);
                            }
                            hits.push(LookupHit {
                                index,
                                squared_distance,
                            });
                        }
                    }
                }
            }
        }

        let found = !hits.is_empty();
        (hits, found)
    }
}

fn ceil_positive(value: f32) -> usize {
    value.ceil().max(1.0) as usize
}

fn grid_axis(value: f32, min: f32, delta: f32, size: usize) -> usize {
    (((value - min) / delta).floor() as isize).clamp(0, size.saturating_sub(1) as isize) as usize
}

fn query_min_axis(value: f32, radius: f32, min: f32, delta: f32, size: usize) -> Option<usize> {
    let axis = ((value - radius - min) / delta).floor() as isize;
    let axis = axis.max(0);
    (axis < size as isize).then_some(axis as usize)
}

fn query_max_axis(value: f32, radius: f32, min: f32, delta: f32, size: usize) -> Option<usize> {
    let axis = ((value + radius - min) / delta).floor() as isize;
    let axis = axis.min(size as isize - 1);
    (axis >= 0).then_some(axis as usize)
}

fn grid_index(x: usize, y: usize, z: usize, size: [usize; 3]) -> usize {
    ((x * size[1]) + y) * size[2] + z
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StructureLookupHit {
    pub unit_id: usize,
    pub index: usize,
    pub squared_distance: f32,
}

#[derive(Clone, Debug, Default)]
pub struct StructureLookup3D {
    unit_lookups: Vec<(usize, UnitLookup3D)>,
    unit_centers: Vec<Vec3>,
    unit_radii: Vec<f32>,
    max_unit_radius: f32,
    unit_grid: UnitLookupGrid,
    pub boundary: Boundary,
}

impl StructureLookup3D {
    fn from_units(units: &[StructureUnit]) -> Self {
        let boundary = structure_boundary(units);
        let mut unit_centers = Vec::with_capacity(units.len());
        let mut unit_radii = Vec::with_capacity(units.len());
        let mut max_unit_radius = 0.0f32;
        for unit in units {
            let sphere = &unit.props.boundary.sphere;
            unit_centers.push(sphere.center);
            unit_radii.push(sphere.radius);
            max_unit_radius = max_unit_radius.max(sphere.radius);
        }
        let unit_grid = UnitLookupGrid::new(&unit_centers, boundary.clone());
        let unit_lookups = units
            .iter()
            .map(|unit| (unit.id, unit.props.lookup3d.clone()))
            .collect();
        StructureLookup3D {
            unit_lookups,
            unit_centers,
            unit_radii,
            max_unit_radius,
            unit_grid,
            boundary,
        }
    }

    pub fn find(&self, point: Vec3, radius: f32) -> Vec<StructureLookupHit> {
        let mut hits = Vec::new();
        let close_units = self.find_unit_indices(point, radius, false).0;
        for close_unit in close_units {
            let Some((unit_id, lookup)) = self.unit_lookups.get(close_unit.index) else {
                continue;
            };
            hits.extend(
                lookup
                    .find(point, radius)
                    .into_iter()
                    .map(|hit| StructureLookupHit {
                        unit_id: *unit_id,
                        index: hit.index,
                        squared_distance: hit.squared_distance,
                    }),
            );
        }
        hits
    }

    pub fn check(&self, point: Vec3, radius: f32) -> bool {
        let close_units = self.find_unit_indices(point, radius, false).0;
        close_units.into_iter().any(|close_unit| {
            self.unit_lookups
                .get(close_unit.index)
                .is_some_and(|(_, lookup)| lookup.check(point, radius))
        })
    }

    pub fn nearest(&self, point: Vec3, k: usize) -> Vec<StructureLookupHit> {
        let mut hits = self
            .unit_lookups
            .iter()
            .flat_map(|(unit_id, lookup)| {
                lookup
                    .nearest(point, lookup.len())
                    .into_iter()
                    .map(|hit| StructureLookupHit {
                        unit_id: *unit_id,
                        index: hit.index,
                        squared_distance: hit.squared_distance,
                    })
            })
            .collect::<Vec<_>>();
        hits.sort_by(|a, b| a.squared_distance.total_cmp(&b.squared_distance));
        hits.truncate(k);
        hits
    }

    fn find_unit_indices(
        &self,
        point: Vec3,
        radius: f32,
        is_check: bool,
    ) -> (Vec<LookupHit>, bool) {
        self.unit_grid.find_with_radii(
            &self.unit_centers,
            &self.unit_radii,
            self.max_unit_radius,
            point,
            radius,
            is_check,
        )
    }

    #[cfg(test)]
    pub(crate) fn unit_grid_size(&self) -> [usize; 3] {
        self.unit_grid.size
    }

    #[cfg(test)]
    pub(crate) fn unit_bucket_offsets(&self) -> &[usize] {
        &self.unit_grid.bucket_offsets
    }

    #[cfg(test)]
    pub(crate) fn unit_bucket_counts(&self) -> &[usize] {
        &self.unit_grid.bucket_counts
    }

    #[cfg(test)]
    pub(crate) fn unit_bucket_array(&self) -> &[usize] {
        &self.unit_grid.bucket_array
    }

    #[cfg(test)]
    pub(crate) fn close_unit_indices(&self, point: Vec3, radius: f32) -> Vec<usize> {
        self.find_unit_indices(point, radius, false)
            .0
            .into_iter()
            .map(|hit| hit.index)
            .collect()
    }
}

fn structure_boundary(units: &[StructureUnit]) -> Boundary {
    let spheres = units
        .iter()
        .map(|unit| unit.props.boundary.sphere.clone())
        .collect::<Vec<_>>();
    Boundary::from_spheres(
        &spheres,
        if units.len() > 500 {
            EposQuality::Coarse
        } else {
            EposQuality::Fine
        },
    )
}

fn structure_principal_axes(units: &[StructureUnit]) -> PrincipalAxes {
    let positions = units
        .iter()
        .flat_map(|unit| unit.props.lookup3d.positions.iter().copied())
        .collect::<Vec<_>>();
    PrincipalAxes::of_positions(&positions)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InterUnitBond {
    pub unit_a: usize,
    pub index_a: usize,
    pub unit_b: usize,
    pub index_b: usize,
    pub source_bond: usize,
    pub order: i8,
    pub flags: BondFlags,
    pub key: i32,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct InterUnitBonds {
    pub edge_count: usize,
    pub edges: Vec<InterUnitBondEdge>,
    unit_pairs: Vec<(usize, Vec<InterUnitBondUnitPairEdges>)>,
    edge_key_index: BTreeMap<u64, BTreeMap<u64, usize>>,
    vertex_key_index: BTreeMap<u64, Vec<usize>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InterUnitBondUnitPairEdges {
    pub unit_a: usize,
    pub unit_b: usize,
    pub edge_count: usize,
    pub connected_indices: Vec<usize>,
    edge_map: Vec<(usize, Vec<InterUnitBondInfo>)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InterUnitBondInfo {
    pub index_b: usize,
    pub props: InterUnitBondProps,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InterUnitBondEdge {
    pub unit_a: usize,
    pub unit_b: usize,
    pub index_a: usize,
    pub index_b: usize,
    pub props: InterUnitBondProps,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct InterUnitBondProps {
    pub order: i8,
    pub flag: BondFlags,
    pub key: i32,
}

impl InterUnitBonds {
    fn from_bonds(bonds: &[InterUnitBond]) -> Self {
        let mut builder = InterUnitBondGraphBuilder::default();
        let mut canonical = bonds.iter().collect::<Vec<_>>();
        canonical.sort_by_key(|bond| inter_unit_bond_sort_key(bond));
        for bond in canonical {
            let props = InterUnitBondProps {
                order: bond.order,
                flag: bond.flags,
                key: bond.key,
            };
            if bond.unit_a <= bond.unit_b {
                builder.add(bond.unit_a, bond.unit_b, bond.index_a, bond.index_b, props);
            } else {
                builder.add(bond.unit_b, bond.unit_a, bond.index_b, bond.index_a, props);
            }
        }
        builder.finish()
    }

    pub fn get_connected_units(&self, unit: usize) -> &[InterUnitBondUnitPairEdges] {
        self.unit_pairs
            .iter()
            .find_map(|(unit_id, pairs)| (*unit_id == unit).then_some(pairs.as_slice()))
            .unwrap_or(&[])
    }

    pub fn get_edge_index(
        &self,
        index_a: usize,
        unit_a: usize,
        index_b: usize,
        unit_b: usize,
    ) -> Option<usize> {
        self.edge_key_index
            .get(&cantor_pairing(unit_a, unit_b))
            .and_then(|indices| indices.get(&cantor_pairing(index_a, index_b)))
            .copied()
    }

    pub fn has_edge(&self, index_a: usize, unit_a: usize, index_b: usize, unit_b: usize) -> bool {
        self.get_edge_index(index_a, unit_a, index_b, unit_b)
            .is_some()
    }

    pub fn get_edge(
        &self,
        index_a: usize,
        unit_a: usize,
        index_b: usize,
        unit_b: usize,
    ) -> Option<&InterUnitBondEdge> {
        self.get_edge_index(index_a, unit_a, index_b, unit_b)
            .and_then(|index| self.edges.get(index))
    }

    pub fn get_edge_indices(&self, index: usize, unit: usize) -> &[usize] {
        self.vertex_key_index
            .get(&cantor_pairing(index, unit))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }
}

fn inter_unit_bond_sort_key(bond: &InterUnitBond) -> (usize, usize, usize, usize, i32, usize) {
    if bond.unit_a <= bond.unit_b {
        (
            bond.unit_a,
            bond.unit_b,
            bond.index_a,
            bond.index_b,
            bond.key,
            bond.source_bond,
        )
    } else {
        (
            bond.unit_b,
            bond.unit_a,
            bond.index_b,
            bond.index_a,
            bond.key,
            bond.source_bond,
        )
    }
}

impl InterUnitBondUnitPairEdges {
    pub fn has_edges(&self, index_a: usize) -> bool {
        self.edge_map.iter().any(|(index, _)| *index == index_a)
    }

    pub fn get_edges(&self, index_a: usize) -> &[InterUnitBondInfo] {
        self.edge_map
            .iter()
            .find_map(|(index, edges)| (*index == index_a).then_some(edges.as_slice()))
            .unwrap_or(&[])
    }

    pub fn are_units_ordered(&self) -> bool {
        self.unit_a < self.unit_b
    }
}

#[derive(Default)]
struct InterUnitBondGraphBuilder {
    pairs: Vec<InterUnitBondPairBuilder>,
}

struct InterUnitBondPairBuilder {
    unit_a: usize,
    unit_b: usize,
    map_ab: Vec<(usize, Vec<InterUnitBondInfo>)>,
    map_ba: Vec<(usize, Vec<InterUnitBondInfo>)>,
    linked_a: Vec<usize>,
    linked_b: Vec<usize>,
    link_count: usize,
}

impl InterUnitBondGraphBuilder {
    fn add(
        &mut self,
        unit_a: usize,
        unit_b: usize,
        index_a: usize,
        index_b: usize,
        props: InterUnitBondProps,
    ) {
        let pair_index = self
            .pairs
            .iter()
            .position(|pair| pair.unit_a == unit_a && pair.unit_b == unit_b)
            .unwrap_or_else(|| {
                self.pairs.push(InterUnitBondPairBuilder {
                    unit_a,
                    unit_b,
                    map_ab: Vec::new(),
                    map_ba: Vec::new(),
                    linked_a: Vec::new(),
                    linked_b: Vec::new(),
                    link_count: 0,
                });
                self.pairs.len() - 1
            });
        let pair = &mut self.pairs[pair_index];
        add_inter_unit_edge_map_entry(
            &mut pair.map_ab,
            index_a,
            InterUnitBondInfo { index_b, props },
        );
        add_inter_unit_edge_map_entry(
            &mut pair.map_ba,
            index_b,
            InterUnitBondInfo {
                index_b: index_a,
                props,
            },
        );
        add_unique_index(&mut pair.linked_a, index_a);
        add_unique_index(&mut pair.linked_b, index_b);
        pair.link_count += 1;
    }

    fn finish(self) -> InterUnitBonds {
        let mut unit_pairs = Vec::<(usize, Vec<InterUnitBondUnitPairEdges>)>::new();
        for pair in self.pairs {
            if pair.link_count == 0 {
                continue;
            }
            add_unit_pair_edges(
                &mut unit_pairs,
                pair.unit_a,
                InterUnitBondUnitPairEdges {
                    unit_a: pair.unit_a,
                    unit_b: pair.unit_b,
                    edge_count: pair.link_count,
                    connected_indices: pair.linked_a,
                    edge_map: pair.map_ab,
                },
            );
            add_unit_pair_edges(
                &mut unit_pairs,
                pair.unit_b,
                InterUnitBondUnitPairEdges {
                    unit_a: pair.unit_b,
                    unit_b: pair.unit_a,
                    edge_count: pair.link_count,
                    connected_indices: pair.linked_b,
                    edge_map: pair.map_ba,
                },
            );
        }
        InterUnitBonds::from_unit_pairs(unit_pairs)
    }
}

impl InterUnitBonds {
    fn from_unit_pairs(unit_pairs: Vec<(usize, Vec<InterUnitBondUnitPairEdges>)>) -> Self {
        let mut edge_count = 0usize;
        let mut edges = Vec::new();
        let mut edge_key_index = BTreeMap::<u64, BTreeMap<u64, usize>>::new();
        let mut vertex_key_index = BTreeMap::<u64, Vec<usize>>::new();
        for (_, pair_edges_array) in &unit_pairs {
            for pair_edges in pair_edges_array {
                edge_count += pair_edges.edge_count;
                for &index_a in &pair_edges.connected_indices {
                    for edge_info in pair_edges.get_edges(index_a) {
                        let unit_a = pair_edges.unit_a;
                        let unit_b = pair_edges.unit_b;
                        let edge_unit_key = cantor_pairing(unit_a, unit_b);
                        let edge_index_key = cantor_pairing(index_a, edge_info.index_b);
                        edge_key_index
                            .entry(edge_unit_key)
                            .or_default()
                            .insert(edge_index_key, edges.len());
                        vertex_key_index
                            .entry(cantor_pairing(index_a, unit_a))
                            .or_default()
                            .push(edges.len());
                        edges.push(InterUnitBondEdge {
                            unit_a,
                            unit_b,
                            index_a,
                            index_b: edge_info.index_b,
                            props: edge_info.props,
                        });
                    }
                }
            }
        }
        InterUnitBonds {
            edge_count,
            edges,
            unit_pairs,
            edge_key_index,
            vertex_key_index,
        }
    }
}

fn add_inter_unit_edge_map_entry(
    map: &mut Vec<(usize, Vec<InterUnitBondInfo>)>,
    index: usize,
    edge: InterUnitBondInfo,
) {
    if let Some((_, edges)) = map
        .iter_mut()
        .find(|(entry_index, _)| *entry_index == index)
    {
        edges.push(edge);
    } else {
        map.push((index, vec![edge]));
    }
}

fn add_unit_pair_edges(
    unit_pairs: &mut Vec<(usize, Vec<InterUnitBondUnitPairEdges>)>,
    unit: usize,
    pair_edges: InterUnitBondUnitPairEdges,
) {
    if let Some((_, pairs)) = unit_pairs
        .iter_mut()
        .find(|(entry_unit, _)| *entry_unit == unit)
    {
        pairs.push(pair_edges);
    } else {
        unit_pairs.push((unit, vec![pair_edges]));
    }
}

fn add_unique_index(indices: &mut Vec<usize>, index: usize) {
    if !indices.contains(&index) {
        indices.push(index);
    }
}

fn cantor_pairing(a: usize, b: usize) -> u64 {
    let a = a as u64;
    let b = b as u64;
    (a + b) * (a + b + 1) / 2 + b
}

#[derive(Clone, Debug)]
pub struct UnitSymmetryGroup {
    pub kind: UnitKind,
    pub model_id: usize,
    pub invariant_id: usize,
    pub elements: Vec<usize>,
    pub unit_ids: Vec<usize>,
    pub operator_names: Vec<String>,
    pub operator_instance_ids: Vec<String>,
    pub unit_index_map: Vec<(usize, usize)>,
    pub hash_code: u32,
    pub transform_hash: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ChainKey {
    id: String,
    entity_id: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ResidueKey {
    chain_index: usize,
    label_seq_id: String,
    auth_seq_id: String,
    insertion_code: String,
}

fn structure_units(
    molecule: &Molecule,
    model: &AtomicModel,
    coarse: &CoarseModel,
    ranges: &AtomicRanges,
) -> Vec<StructureUnit> {
    let asymmetric = asymmetric_structure_unit_build(model, ranges);
    let mut units = if let Some(assembly) = &molecule.selected_assembly {
        assembly_structure_units(assembly, model, ranges, &asymmetric.units)
    } else {
        asymmetric.units
    };
    append_coarse_units(
        &mut units,
        molecule.selected_assembly.as_ref(),
        coarse,
        asymmetric.next_invariant_id,
        asymmetric.next_chain_group_id,
    );
    sort_structure_units_molstar(&mut units);
    units
}

pub(crate) fn sort_structure_units_molstar(units: &mut [StructureUnit]) {
    if !units_are_molstar_sorted(units) {
        units.sort_by_key(|unit| unit.id);
    }
}

fn units_are_molstar_sorted(units: &[StructureUnit]) -> bool {
    units.windows(2).all(|pair| pair[0].id <= pair[1].id)
}

struct StructureUnitBuild {
    units: Vec<StructureUnit>,
    next_invariant_id: usize,
    next_chain_group_id: usize,
}

fn asymmetric_structure_unit_build(
    model: &AtomicModel,
    ranges: &AtomicRanges,
) -> StructureUnitBuild {
    let hierarchy = &model.hierarchy;
    let mut builder = StructureUnitBuilder::new(model, ranges);
    let is_coarse_grained = model_is_coarse_grained(model);
    let mut chain_index = 0usize;
    while chain_index < hierarchy.chains.len() {
        let start_chain = chain_index;
        let mut end_chain = chain_index;
        let mut is_multi_chain = false;
        let is_water = is_water_chain(hierarchy, chain_index);

        if is_water {
            while end_chain + 1 < hierarchy.chains.len()
                && is_water_chain(hierarchy, end_chain + 1)
                && chains_have_same_operator(hierarchy, end_chain, end_chain + 1)
            {
                is_multi_chain = true;
                end_chain += 1;
            }
        } else {
            while end_chain + 1 < hierarchy.chains.len()
                && chain_atom_count(hierarchy, end_chain) == 1
                && chain_atom_count(hierarchy, end_chain + 1) == 1
                && chains_have_same_entity_and_auth_asym(hierarchy, end_chain, end_chain + 1)
                && chains_have_same_operator(hierarchy, end_chain, end_chain + 1)
            {
                is_multi_chain = true;
                end_chain += 1;
            }
        }

        let mut traits = UnitTraits::NONE;
        if is_water {
            traits = traits.union(UnitTraits::WATER);
        }
        if is_multi_chain {
            traits = traits.union(UnitTraits::MULTI_CHAIN);
        }
        if is_coarse_grained {
            traits = traits.union(UnitTraits::COARSE_GRAINED);
        }
        let operator_name = first_chain_operator(hierarchy, start_chain)
            .unwrap_or("1_555")
            .to_string();
        builder.add_unit(UnitBuildInput {
            chain_indices: (start_chain..=end_chain).collect(),
            operator_name: Some(&operator_name),
            kind: UnitKind::Atomic,
            operator: asymmetric_operator(&operator_name),
            traits,
            invariant_id: None,
            chain_group_id: None,
        });
        chain_index = end_chain + 1;
    }
    StructureUnitBuild {
        units: builder.units,
        next_invariant_id: builder.next_invariant_id,
        next_chain_group_id: builder.next_chain_group_id,
    }
}

fn model_is_coarse_grained(model: &AtomicModel) -> bool {
    let hierarchy = &model.hierarchy;
    let mut polymer_residue_count = 0usize;
    let mut polymer_direction_count = 0usize;
    for (polymer_type, direction_to) in hierarchy
        .derived
        .residue
        .polymer_type
        .iter()
        .zip(&hierarchy.derived.residue.direction_to_element_index)
    {
        if *polymer_type != PolymerType::None {
            polymer_residue_count += 1;
            if direction_to.is_some() {
                polymer_direction_count += 1;
            }
        }
    }

    let mut has_bb = false;
    let mut has_sc1 = false;
    let mut has_trace = false;
    for atom in &hierarchy.atoms {
        match atom.name.as_str() {
            "BB" => has_bb = true,
            "SC1" => has_sc1 = true,
            "CA" | "CA1" | "BAS" | "O3'" | "O3*" | "N4'" | "N4*" => has_trace = true,
            _ => {}
        }
        if has_bb && has_sc1 && has_trace {
            break;
        }
    }

    !hierarchy.atoms.is_empty()
        && polymer_residue_count > 0
        && ((has_bb && has_sc1)
            || (hierarchy.atoms.len() as f32 / polymer_residue_count as f32) < 3.0
            || (polymer_direction_count == 0 && has_trace))
}

fn assembly_structure_units(
    assembly: &Assembly,
    model: &AtomicModel,
    ranges: &AtomicRanges,
    base_units: &[StructureUnit],
) -> Vec<StructureUnit> {
    let hierarchy = &model.hierarchy;
    let generators = if assembly.generators.is_empty() {
        vec![AssemblyGenerator::from_transforms(
            &assembly.id,
            assembly.asym_ids.clone(),
            0,
            assembly.transforms.clone(),
            vec![Vec::new(); assembly.transforms.len()],
        )]
    } else {
        assembly.generators.clone()
    };
    let mut builder = StructureUnitBuilder::new(model, ranges);
    let mut operator_offset = 0usize;
    for generator in &generators {
        let operators = generator.operators_for_assembly(&assembly.id, operator_offset);
        operator_offset += operators.len();
        for operator in &operators {
            for unit in base_units {
                let chain_indices =
                    unit_chain_indices_matching_asym_ids(unit, hierarchy, &generator.asym_ids);
                if chain_indices.is_empty() {
                    continue;
                }
                builder.add_unit(UnitBuildInput {
                    chain_indices,
                    operator_name: Some(&unit.operator.name),
                    kind: unit.kind,
                    operator: UnitOperator {
                        name: operator.name.clone(),
                        instance_id: operator.instance_id.clone(),
                        assembly_id: operator.assembly_id.clone(),
                        oper_id: operator.oper_id as i32,
                        oper_list_ids: operator.oper_list_ids.clone(),
                        transform: operator.transform,
                        is_identity: operator.transform.is_identity(),
                        suffix: if operator.transform.is_identity() {
                            String::new()
                        } else {
                            format!("_{}", operator.oper_id)
                        },
                    },
                    traits: unit.traits,
                    invariant_id: Some(unit.invariant_id),
                    chain_group_id: Some(unit.chain_group_id),
                });
            }
        }
    }
    builder.units
}

fn append_coarse_units(
    units: &mut Vec<StructureUnit>,
    assembly: Option<&Assembly>,
    coarse: &CoarseModel,
    next_invariant_id: usize,
    next_chain_group_id: usize,
) {
    let mut next_invariant_id = next_invariant_id;
    let mut next_chain_group_id = next_chain_group_id;
    append_coarse_kind_units(
        units,
        assembly,
        UnitKind::Spheres,
        &coarse.hierarchy.spheres,
        &(0..coarse.conformation.spheres.len())
            .filter_map(|index| {
                let position = coarse.conformation.spheres.position(index)?;
                let radius = *coarse.conformation.spheres.radius.get(index)?;
                Some((position, radius))
            })
            .collect::<Vec<_>>(),
        &mut next_invariant_id,
        &mut next_chain_group_id,
    );
    append_coarse_kind_units(
        units,
        assembly,
        UnitKind::Gaussians,
        &coarse.hierarchy.gaussians,
        &(0..coarse.conformation.gaussians.len())
            .filter_map(|index| {
                let position = coarse.conformation.gaussians.position(index)?;
                Some((position, 0.0))
            })
            .collect::<Vec<_>>(),
        &mut next_invariant_id,
        &mut next_chain_group_id,
    );
}

fn append_coarse_kind_units(
    units: &mut Vec<StructureUnit>,
    assembly: Option<&Assembly>,
    kind: UnitKind,
    elements: &CoarseElements,
    positions_and_radii: &[(Vec3, f32)],
    next_invariant_id: &mut usize,
    next_chain_group_id: &mut usize,
) {
    if elements.elements.is_empty() {
        return;
    }
    let operator_groups = coarse_operator_groups(assembly);
    let invariant_offset = *next_invariant_id;
    let chain_group_offset = *next_chain_group_id;
    *next_invariant_id += elements.chain_element_segments.count;
    *next_chain_group_id += elements.chain_element_segments.count;
    for chain_index in 0..elements.chain_element_segments.count {
        let Some(&start) = elements.chain_element_segments.offsets.get(chain_index) else {
            continue;
        };
        let Some(&end) = elements.chain_element_segments.offsets.get(chain_index + 1) else {
            continue;
        };
        let chain_elements = &elements.elements[start..end];
        let Some(first) = chain_elements.first() else {
            continue;
        };
        for operator in &operator_groups {
            if !coarse_operator_matches_asym(operator, &first.asym_id) {
                continue;
            }
            let element_indices = (start..end).collect::<Vec<_>>();
            let source_indices = chain_elements
                .iter()
                .map(|element| element.source_index)
                .collect::<Vec<_>>();
            let mut transformed_positions = Vec::new();
            let mut invariant_positions = Vec::new();
            let mut radii = Vec::new();
            for index in &element_indices {
                if let Some((position, radius)) = positions_and_radii.get(*index) {
                    invariant_positions.push(*position);
                    transformed_positions.push(operator.operator.transform.apply(*position));
                    radii.push(*radius);
                }
            }
            let props = UnitProps::from_coarse_positions(
                transformed_positions,
                invariant_positions,
                radii,
                &element_indices,
                &elements.polymer_ranges,
                &elements.gap_ranges,
            );
            let unit_id = units.len();
            units.push(StructureUnit {
                id: unit_id,
                invariant_id: invariant_offset + chain_index,
                chain_group_id: chain_group_offset + chain_index,
                kind,
                traits: UnitTraits::NONE,
                model_index: 0,
                chain_index,
                chain_indices: vec![chain_index],
                elements: element_indices.clone(),
                atom_indices: source_indices,
                residue_indices: element_indices.clone(),
                residue_index_by_element: element_indices,
                chain_index_by_element: vec![chain_index; end - start],
                props,
                operator: operator.operator.clone(),
            });
        }
    }
}

#[derive(Clone, Debug)]
struct CoarseOperatorGroup {
    asym_ids: Vec<String>,
    operator: UnitOperator,
}

fn coarse_operator_groups(assembly: Option<&Assembly>) -> Vec<CoarseOperatorGroup> {
    let Some(assembly) = assembly else {
        return vec![CoarseOperatorGroup {
            asym_ids: Vec::new(),
            operator: UnitOperator::default(),
        }];
    };
    let generators = if assembly.generators.is_empty() {
        vec![AssemblyGenerator::from_transforms(
            &assembly.id,
            assembly.asym_ids.clone(),
            0,
            assembly.transforms.clone(),
            vec![Vec::new(); assembly.transforms.len()],
        )]
    } else {
        assembly.generators.clone()
    };
    let mut groups = Vec::new();
    let mut operator_offset = 0usize;
    for generator in generators {
        let operators = generator.operators_for_assembly(&assembly.id, operator_offset);
        operator_offset += operators.len();
        for operator in operators {
            groups.push(CoarseOperatorGroup {
                asym_ids: generator.asym_ids.clone(),
                operator: UnitOperator {
                    name: operator.name,
                    instance_id: operator.instance_id,
                    assembly_id: operator.assembly_id,
                    oper_id: operator.oper_id as i32,
                    oper_list_ids: operator.oper_list_ids,
                    transform: operator.transform,
                    is_identity: operator.transform.is_identity(),
                    suffix: if operator.transform.is_identity() {
                        String::new()
                    } else {
                        format!("_{}", operator.oper_id)
                    },
                },
            });
        }
    }
    groups
}

fn coarse_operator_matches_asym(operator: &CoarseOperatorGroup, asym_id: &str) -> bool {
    operator.asym_ids.is_empty() || operator.asym_ids.iter().any(|id| id == asym_id)
}

struct UnitBuildInput<'a> {
    chain_indices: Vec<usize>,
    operator_name: Option<&'a str>,
    kind: UnitKind,
    operator: UnitOperator,
    traits: UnitTraits,
    invariant_id: Option<usize>,
    chain_group_id: Option<usize>,
}

struct StructureUnitBuilder<'a> {
    model: &'a AtomicModel,
    ranges: &'a AtomicRanges,
    units: Vec<StructureUnit>,
    next_invariant_id: usize,
    next_chain_group_id: usize,
    single_element_keys: Vec<(usize, u32, u32, u32)>,
}

impl<'a> StructureUnitBuilder<'a> {
    fn new(model: &'a AtomicModel, ranges: &'a AtomicRanges) -> Self {
        Self {
            model,
            ranges,
            units: Vec::new(),
            next_invariant_id: 0,
            next_chain_group_id: 0,
            single_element_keys: Vec::new(),
        }
    }

    fn add_unit(&mut self, input: UnitBuildInput<'_>) {
        let hierarchy = &self.model.hierarchy;
        let UnitBuildInput {
            chain_indices,
            operator_name,
            kind,
            operator,
            traits,
            invariant_id,
            chain_group_id,
        } = input;
        let Some(&chain_index) = chain_indices.first() else {
            return;
        };
        let elements = chain_indices
            .iter()
            .flat_map(|chain_index| {
                hierarchy
                    .atoms
                    .iter()
                    .enumerate()
                    .filter(move |(_, atom)| atom.chain_index == *chain_index)
                    .filter(move |(_, atom)| match operator_name {
                        Some(operator_name) => atom.operator_name == operator_name,
                        None => true,
                    })
                    .map(|(element_index, _)| element_index)
            })
            .collect::<Vec<_>>();
        if elements.is_empty() {
            return;
        }
        let invariant_id = invariant_id.unwrap_or_else(|| {
            let id = self.next_invariant_id;
            self.next_invariant_id += 1;
            id
        });
        let chain_group_id = chain_group_id.unwrap_or_else(|| {
            let id = self.next_chain_group_id;
            self.next_chain_group_id += 1;
            id
        });
        if elements.len() == 1 {
            let position = operator
                .transform
                .apply(self.model.conformation.positions[elements[0]]);
            let key = (
                invariant_id,
                canonical_single_element_key_coordinate(position.x).to_bits(),
                canonical_single_element_key_coordinate(position.y).to_bits(),
                canonical_single_element_key_coordinate(position.z).to_bits(),
            );
            if self.single_element_keys.contains(&key) {
                return;
            }
            self.single_element_keys.push(key);
        }
        let atom_indices = elements
            .iter()
            .map(|element_index| hierarchy.atoms[*element_index].source_index)
            .collect::<Vec<_>>();
        let residue_indices = chain_indices
            .iter()
            .flat_map(|chain_index| {
                let chain = &hierarchy.chains[*chain_index];
                (chain.start_residue..chain.end_residue).filter(move |residue_index| {
                    hierarchy
                        .residues
                        .get(*residue_index)
                        .is_some_and(|residue| residue.chain_index == *chain_index)
                })
            })
            .collect::<Vec<_>>();
        let residue_index_by_element = hierarchy.residue_atom_segments.index.clone();
        let chain_index_by_element = hierarchy.chain_atom_segments.index.clone();
        let props = UnitProps::from_elements(
            hierarchy,
            self.ranges,
            &elements,
            &residue_indices,
            &operator,
        );
        let unit_id = self.units.len();
        self.units.push(StructureUnit {
            id: unit_id,
            invariant_id,
            chain_group_id,
            kind,
            traits,
            model_index: 0,
            chain_index,
            chain_indices,
            elements,
            atom_indices,
            residue_indices,
            residue_index_by_element,
            chain_index_by_element,
            props,
            operator,
        });
    }
}

fn canonical_single_element_key_coordinate(value: f32) -> f32 {
    if value == 0.0 {
        0.0
    } else {
        value
    }
}

fn chain_atom_count(hierarchy: &AtomicHierarchy, chain_index: usize) -> usize {
    hierarchy
        .chain_atom_segments
        .offsets
        .get(chain_index + 1)
        .zip(hierarchy.chain_atom_segments.offsets.get(chain_index))
        .map(|(end, start)| end.saturating_sub(*start))
        .unwrap_or(0)
}

fn is_water_chain(hierarchy: &AtomicHierarchy, chain_index: usize) -> bool {
    hierarchy.index.entity_type_from_chain(chain_index) == Some("water")
}

fn chains_have_same_entity_and_auth_asym(hierarchy: &AtomicHierarchy, a: usize, b: usize) -> bool {
    let Some(chain_a) = hierarchy.chains.get(a) else {
        return false;
    };
    let Some(chain_b) = hierarchy.chains.get(b) else {
        return false;
    };
    chain_a.entity_id == chain_b.entity_id && chain_a.auth_id == chain_b.auth_id
}

fn chains_have_same_operator(hierarchy: &AtomicHierarchy, a: usize, b: usize) -> bool {
    first_chain_operator(hierarchy, a) == first_chain_operator(hierarchy, b)
}

fn first_chain_operator(hierarchy: &AtomicHierarchy, chain_index: usize) -> Option<&str> {
    hierarchy
        .atoms
        .iter()
        .find(|atom| atom.chain_index == chain_index)
        .map(|atom| atom.operator_name.as_str())
}

fn asymmetric_operator(operator_name: &str) -> UnitOperator {
    UnitOperator {
        name: operator_name.to_string(),
        instance_id: operator_name.to_string(),
        assembly_id: String::new(),
        oper_id: -1,
        oper_list_ids: if operator_name == "1_555" {
            Vec::new()
        } else {
            vec![operator_name.to_string()]
        },
        transform: Transform::identity(),
        is_identity: operator_name == "1_555",
        suffix: String::new(),
    }
}

fn unit_chain_indices_matching_asym_ids(
    unit: &StructureUnit,
    hierarchy: &AtomicHierarchy,
    asym_ids: &[String],
) -> Vec<usize> {
    if asym_ids.is_empty() {
        return unit.chain_indices.clone();
    }
    unit.chain_indices
        .iter()
        .copied()
        .filter(|chain_index| {
            hierarchy
                .chains
                .get(*chain_index)
                .is_some_and(|chain| asym_ids.iter().any(|asym_id| asym_id == &chain.id))
        })
        .collect()
}

#[derive(Default)]
struct IntraUnitBondBuild {
    atom_a: Vec<usize>,
    atom_b: Vec<usize>,
    metadata: Vec<BondMetadata>,
}

fn assign_unit_bonds(
    units: &mut [StructureUnit],
    bonds: &[Bond],
    metadata: &[BondMetadata],
    index_pair_bonds: Option<&IndexPairBonds>,
) -> (usize, Vec<InterUnitBond>, InterUnitBonds) {
    let mut intra_count = 0usize;
    let mut inter = Vec::new();
    let mut intra_builders = (0..units.len())
        .map(|_| IntraUnitBondBuild::default())
        .collect::<Vec<_>>();
    for (source_bond, bond) in bonds.iter().enumerate() {
        let mut endpoints_a = Vec::<(usize, usize, usize)>::new();
        let mut endpoints_b = Vec::<(usize, usize, usize)>::new();
        for (unit_index, unit) in units.iter().enumerate() {
            if let Some(index) = unit
                .atom_indices
                .iter()
                .position(|source_index| *source_index == bond.a)
            {
                endpoints_a.push((unit.id, unit_index, index));
            }
            if let Some(index) = unit
                .atom_indices
                .iter()
                .position(|source_index| *source_index == bond.b)
            {
                endpoints_b.push((unit.id, unit_index, index));
            }
        }
        let bond_metadata = metadata.get(source_bond).cloned().unwrap_or_default();
        for (unit_a_id, unit_a_index, index_a) in &endpoints_a {
            for (unit_b_id, unit_b_index, index_b) in &endpoints_b {
                if unit_a_id == unit_b_id {
                    if !metadata_allows_intra_bond(&bond_metadata, &units[*unit_a_index])
                        || !metadata_allows_bond_distance(
                            &bond_metadata,
                            source_bond,
                            index_pair_bonds,
                            &units[*unit_a_index],
                            *index_a,
                            &units[*unit_b_index],
                            *index_b,
                        )
                    {
                        continue;
                    }
                    intra_count += 1;
                    if let Some(builder) = intra_builders.get_mut(*unit_a_index) {
                        builder.atom_a.push(*index_a);
                        builder.atom_b.push(*index_b);
                        builder.metadata.push(bond_metadata.clone());
                    }
                    if let Some(unit) = units.get_mut(*unit_a_index) {
                        unit.props.intra_unit_bond_count += 1;
                    }
                } else if metadata_allows_inter_unit_bond(
                    &bond_metadata,
                    &units[*unit_a_index],
                    &units[*unit_b_index],
                ) && metadata_allows_bond_distance(
                    &bond_metadata,
                    source_bond,
                    index_pair_bonds,
                    &units[*unit_a_index],
                    *index_a,
                    &units[*unit_b_index],
                    *index_b,
                ) {
                    inter.push(InterUnitBond {
                        unit_a: *unit_a_id,
                        index_a: *index_a,
                        unit_b: *unit_b_id,
                        index_b: *index_b,
                        source_bond,
                        order: bond_metadata.order,
                        flags: bond_metadata.flags,
                        key: bond_metadata.key,
                    });
                    if let Some(unit) = units.get_mut(*unit_a_index) {
                        unit.props.inter_unit_bond_count += 1;
                    }
                    if let Some(unit) = units.get_mut(*unit_b_index) {
                        unit.props.inter_unit_bond_count += 1;
                    }
                }
            }
        }
    }
    for (unit, builder) in units.iter_mut().zip(intra_builders) {
        unit.props.intra_unit_bonds = if unit.kind == UnitKind::Atomic
            && (!builder.atom_a.is_empty() || unit.elements.len() > 1)
        {
            IntraUnitBonds::from_edges(
                &builder.atom_a,
                &builder.atom_b,
                &builder.metadata,
                unit.elements.len(),
                false,
                index_pair_bonds.is_some_and(|index_pairs| index_pairs.cacheable),
            )
        } else {
            IntraUnitBonds::default()
        };
    }
    let inter_graph = InterUnitBonds::from_bonds(&inter);
    (intra_count, inter, inter_graph)
}

fn metadata_allows_intra_bond(metadata: &BondMetadata, unit: &StructureUnit) -> bool {
    if let Some(struct_conn) = &metadata.struct_conn {
        return normalized_symmetry(&struct_conn.partner_a_symmetry)
            == normalized_symmetry(&struct_conn.partner_b_symmetry);
    }
    if metadata.source == BondSource::IndexPair {
        let op_key = unit.operator.oper_id;
        return (metadata.operator_a < 0 || metadata.operator_a == op_key)
            && (metadata.operator_b < 0 || metadata.operator_b == op_key);
    }
    true
}

fn metadata_allows_inter_unit_bond(
    metadata: &BondMetadata,
    unit_a: &StructureUnit,
    unit_b: &StructureUnit,
) -> bool {
    if metadata.struct_conn.is_some() {
        return true;
    }
    if metadata.source == BondSource::IndexPair {
        if metadata.operator_a >= 0 && metadata.operator_b >= 0 {
            return metadata.operator_a != metadata.operator_b
                && operator_matches_key(&unit_a.operator, metadata.operator_a)
                && operator_matches_key(&unit_b.operator, metadata.operator_b);
        }
        return true;
    }
    if metadata.operator_a >= 0 || metadata.operator_b >= 0 {
        return operator_matches_key(&unit_a.operator, metadata.operator_a)
            && operator_matches_key(&unit_b.operator, metadata.operator_b);
    }
    unit_a.operator.instance_id == unit_b.operator.instance_id
}

fn metadata_allows_bond_distance(
    metadata: &BondMetadata,
    source_bond: usize,
    index_pair_bonds: Option<&IndexPairBonds>,
    unit_a: &StructureUnit,
    index_a: usize,
    unit_b: &StructureUnit,
    index_b: usize,
) -> bool {
    if metadata.struct_conn.is_some() && unit_a.id != unit_b.id {
        let Some(position_a) = unit_position(unit_a, index_a) else {
            return false;
        };
        let Some(position_b) = unit_position(unit_b, index_b) else {
            return false;
        };
        return position_a.distance(position_b) <= 4.0;
    }
    if metadata.source != BondSource::IndexPair {
        return true;
    }
    let Some(position_a) = unit_position(unit_a, index_a) else {
        return false;
    };
    let Some(position_b) = unit_position(unit_b, index_b) else {
        return false;
    };
    let distance = position_a.distance(position_b);
    if let Some(expected) = metadata.distance {
        return (distance - expected).abs() <= 0.3;
    }
    if let Some(index_pair_bonds) = index_pair_bonds {
        if index_pair_bonds.contains_bond(source_bond) && index_pair_bonds.max_distance >= 0.0 {
            return distance < index_pair_bonds.max_distance;
        }
    }
    true
}

fn unit_position(unit: &StructureUnit, index: usize) -> Option<Vec3> {
    unit.props.lookup3d.positions.get(index).copied()
}

fn normalized_symmetry(symmetry: &str) -> String {
    match symmetry.trim() {
        "" | "." | "?" => String::new(),
        value => value.to_string(),
    }
}

fn operator_matches_key(operator: &UnitOperator, key: i32) -> bool {
    key < 0 || operator.oper_id == key
}

fn coordinate_system_operator(assembly: Option<&Assembly>) -> UnitOperator {
    let Some(assembly) = assembly else {
        return UnitOperator {
            name: "1_555".to_string(),
            instance_id: "1_555".to_string(),
            assembly_id: String::new(),
            oper_id: -1,
            oper_list_ids: Vec::new(),
            transform: Transform::identity(),
            is_identity: true,
            suffix: String::new(),
        };
    };
    UnitOperator {
        name: assembly.id.clone(),
        instance_id: assembly.id.clone(),
        assembly_id: assembly.id.clone(),
        oper_id: 0,
        oper_list_ids: Vec::new(),
        transform: Transform::identity(),
        is_identity: true,
        suffix: String::new(),
    }
}

fn atomic_models(atoms: &[Atom], molecule: &Molecule) -> Vec<AtomicModel> {
    if !molecule.ihm_model_list.is_empty() {
        let atom_windows = model_windows_by(atoms, |atom| atom.model_num);
        return molecule
            .ihm_model_list
            .iter()
            .enumerate()
            .map(|(model_index, ihm_model)| {
                let atoms = atom_windows
                    .get(&ihm_model.model_id)
                    .map(|(start, end)| &atoms[*start..*end])
                    .unwrap_or(&[]);
                AtomicModel::from_atoms(model_index, ihm_model.model_id, atoms, molecule)
            })
            .collect();
    }

    if atoms.is_empty() {
        return vec![AtomicModel::from_atoms(0, 1, &[], molecule)];
    }

    let mut models = Vec::new();
    let mut start = 0;
    while start < atoms.len() {
        let model_num = atoms[start].model_num;
        let mut end = start + 1;
        while end < atoms.len() && atoms[end].model_num == model_num {
            end += 1;
        }
        models.push(AtomicModel::from_atoms(
            models.len(),
            model_num,
            &atoms[start..end],
            molecule,
        ));
        start = end;
    }
    models
}

fn coarse_models(molecule: &Molecule) -> Vec<CoarseModel> {
    if molecule.ihm_model_list.is_empty() {
        return vec![CoarseModel::from_molecule(molecule)];
    }

    let sphere_windows = model_windows_by(&molecule.coarse_spheres, |sphere| sphere.model_num);
    let gaussian_windows =
        model_windows_by(&molecule.coarse_gaussians, |gaussian| gaussian.model_num);
    molecule
        .ihm_model_list
        .iter()
        .map(|ihm_model| {
            let spheres = sphere_windows
                .get(&ihm_model.model_id)
                .map(|(start, end)| &molecule.coarse_spheres[*start..*end])
                .unwrap_or(&[]);
            let gaussians = gaussian_windows
                .get(&ihm_model.model_id)
                .map(|(start, end)| &molecule.coarse_gaussians[*start..*end])
                .unwrap_or(&[]);
            CoarseModel::from_parts(spheres, gaussians, &molecule.entity_index)
        })
        .collect()
}

fn model_windows_by<T, F>(items: &[T], model_num: F) -> BTreeMap<i32, (usize, usize)>
where
    F: Fn(&T) -> i32,
{
    let mut windows = BTreeMap::new();
    let mut start = 0;
    while start < items.len() {
        let id = model_num(&items[start]);
        let mut end = start + 1;
        while end < items.len() && model_num(&items[end]) == id {
            end += 1;
        }
        windows.insert(id, (start, end));
        start = end;
    }
    windows
}

fn atom_operator_name(atom: &Atom) -> String {
    if atom.operator_name.is_empty() {
        "1_555".to_string()
    } else {
        atom.operator_name.clone()
    }
}

fn unit_symmetry_groups(units: &[StructureUnit]) -> Vec<UnitSymmetryGroup> {
    let mut groups = Vec::<UnitSymmetryGroup>::new();
    for unit in units {
        if let Some(group) = groups.iter_mut().find(|group| {
            group.model_id == unit.model_index
                && group.invariant_id == unit.invariant_id
                && group.elements == unit.elements
        }) {
            group.unit_ids.push(unit.id);
            group.operator_names.push(unit.operator.name.clone());
            group
                .operator_instance_ids
                .push(unit.operator.instance_id.clone());
            group
                .unit_index_map
                .push((unit.id, group.unit_ids.len() - 1));
            group.transform_hash = molstar_fnv32a(&group.unit_ids);
        } else {
            let unit_ids = vec![unit.id];
            let hash_code = molstar_hash_unit(unit.invariant_id, &unit.elements);
            groups.push(UnitSymmetryGroup {
                kind: unit.kind,
                model_id: unit.model_index,
                invariant_id: unit.invariant_id,
                elements: unit.elements.clone(),
                unit_ids: unit_ids.clone(),
                operator_names: vec![unit.operator.name.clone()],
                operator_instance_ids: vec![unit.operator.instance_id.clone()],
                unit_index_map: vec![(unit.id, 0)],
                hash_code,
                transform_hash: molstar_fnv32a(&unit_ids),
            });
        }
    }
    groups
}

fn molstar_hash_unit(invariant_id: usize, elements: &[usize]) -> u32 {
    molstar_hash2(
        usize_to_i32_bits(invariant_id),
        molstar_sorted_array_hash(elements) as i32,
    )
}

fn molstar_sorted_array_hash(values: &[usize]) -> u32 {
    let size = values.len();
    if size == 0 {
        return 0;
    }
    let first = usize_to_i32_bits(values[0]);
    let last = usize_to_i32_bits(values[size - 1]);
    if size > 2 {
        return molstar_hash4(
            size as i32,
            first,
            last,
            usize_to_i32_bits(values[size >> 1]),
        );
    }
    molstar_hash3(size as i32, first, last)
}

fn molstar_hash2(i: i32, j: i32) -> u32 {
    let mut a = 23i32;
    a = a.wrapping_mul(31).wrapping_add(i);
    a = a.wrapping_mul(31).wrapping_add(j);
    molstar_finish_hash(a)
}

fn molstar_hash3(i: i32, j: i32, k: i32) -> u32 {
    let mut a = 23i32;
    a = a.wrapping_mul(31).wrapping_add(i);
    a = a.wrapping_mul(31).wrapping_add(j);
    a = a.wrapping_mul(31).wrapping_add(k);
    molstar_finish_hash(a)
}

fn molstar_hash4(i: i32, j: i32, k: i32, l: i32) -> u32 {
    let mut a = 23i32;
    a = a.wrapping_mul(31).wrapping_add(i);
    a = a.wrapping_mul(31).wrapping_add(j);
    a = a.wrapping_mul(31).wrapping_add(k);
    a = a.wrapping_mul(31).wrapping_add(l);
    molstar_finish_hash(a)
}

fn molstar_finish_hash(mut value: i32) -> u32 {
    value ^= value >> 4;
    value = (value ^ 0xdead_beefu32 as i32).wrapping_add(value.wrapping_shl(5));
    value ^= value >> 11;
    value as u32
}

fn molstar_fnv32a(values: &[usize]) -> u32 {
    let mut hash = 0x811c_9dc5u32;
    for &value in values {
        hash ^= value as u32;
        hash = hash
            .wrapping_add(hash.wrapping_shl(1))
            .wrapping_add(hash.wrapping_shl(4))
            .wrapping_add(hash.wrapping_shl(7))
            .wrapping_add(hash.wrapping_shl(8))
            .wrapping_add(hash.wrapping_shl(24));
    }
    hash
}

fn usize_to_i32_bits(value: usize) -> i32 {
    value as u32 as i32
}
