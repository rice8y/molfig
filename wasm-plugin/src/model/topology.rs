use std::collections::{BTreeMap, HashSet};

use super::derived::{is_dna_residue, is_protein_residue, is_rna_residue};
use super::{Atom, Bond, BondFlags, BondMetadata, BondSource, Molecule, Vec3};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Ring {
    pub atom_indices: Vec<usize>,
    pub aromatic: bool,
    pub fingerprint: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Resonance {
    pub ring_count: usize,
    pub aromatic_ring_count: usize,
    pub delocalized_bond_count: usize,
    pub delocalized_triplets: Vec<[usize; 3]>,
    pub delocalized_triplet_lookup: DelocalizedTriplets,
    pub element_ring_indices: Vec<Vec<usize>>,
    pub element_aromatic_ring_indices: Vec<Vec<usize>>,
    pub ring_component_index: Vec<usize>,
    pub ring_components: Vec<Vec<usize>>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DelocalizedTriplets {
    pub triplets: Vec<[usize; 3]>,
    third_element_by_pair: BTreeMap<usize, usize>,
    triplet_indices_by_element: BTreeMap<usize, Vec<usize>>,
}

impl DelocalizedTriplets {
    pub fn get_third_element(&self, a: usize, b: usize) -> Option<usize> {
        self.third_element_by_pair
            .get(&sorted_cantor_pairing(a, b))
            .copied()
    }

    pub fn get_triplet_indices(&self, element: usize) -> Option<&[usize]> {
        self.triplet_indices_by_element
            .get(&element)
            .map(Vec::as_slice)
    }

    fn add(&mut self, center: usize, neighbor_a: usize, neighbor_b: usize) {
        let index = self.triplets.len();
        self.triplets
            .push(sorted_triplet(center, neighbor_a, neighbor_b));
        self.third_element_by_pair
            .insert(sorted_cantor_pairing(center, neighbor_a), neighbor_b);
        self.triplet_indices_by_element
            .entry(center)
            .or_default()
            .push(index);
    }
}

pub(super) fn apply_chemical_component_bonds(molecule: &mut Molecule) {
    if molecule.chemical_component_bonds.is_empty() || molecule.atoms.is_empty() {
        return;
    }
    let mut existing = molecule
        .bonds
        .iter()
        .map(|bond| normalized_bond_pair(bond.a, bond.b))
        .collect::<HashSet<_>>();
    let mut residue_atoms = Vec::<(i32, String, String, String, String, Vec<usize>)>::new();
    for (atom_index, atom) in molecule.atoms.iter().enumerate() {
        let key = (
            atom.model_num,
            atom.chain.clone(),
            atom.residue.clone(),
            atom.residue_seq.clone(),
            atom.insertion_code.clone(),
        );
        if let Some((_, _, _, _, _, atoms)) = residue_atoms.iter_mut().find(|entry| {
            entry.0 == key.0
                && entry.1 == key.1
                && entry.2 == key.2
                && entry.3 == key.3
                && entry.4 == key.4
        }) {
            atoms.push(atom_index);
        } else {
            residue_atoms.push((key.0, key.1, key.2, key.3, key.4, vec![atom_index]));
        }
    }
    for (_, _, comp_id, _, _, atom_indices) in residue_atoms {
        for chem_bond in molecule
            .chemical_component_bonds
            .iter()
            .filter(|bond| bond.comp_id.eq_ignore_ascii_case(&comp_id))
        {
            let a_candidates = chemical_component_bond_atom_candidates(
                &molecule.atoms,
                &atom_indices,
                &comp_id,
                &chem_bond.atom_id_1,
            );
            let b_candidates = chemical_component_bond_atom_candidates(
                &molecule.atoms,
                &atom_indices,
                &comp_id,
                &chem_bond.atom_id_2,
            );
            for &a in &a_candidates {
                for &b in &b_candidates {
                    if a == b
                        || !chemical_component_bond_alt_locs_are_compatible(
                            &molecule.atoms[a],
                            &molecule.atoms[b],
                        )
                    {
                        continue;
                    }
                    let pair = normalized_bond_pair(a, b);
                    if existing.insert(pair) {
                        molecule.bonds.push(Bond {
                            a: pair.0,
                            b: pair.1,
                        });
                        molecule.bond_metadata.push(BondMetadata {
                            source: BondSource::ChemComp,
                            order: chem_bond.order,
                            flags: BondFlags::COVALENT.union(chem_bond.flags),
                            key: chem_bond.ordinal.unwrap_or(0),
                            distance: None,
                            operator_a: -1,
                            operator_b: -1,
                            struct_conn: None,
                        });
                    }
                }
            }
        }
    }
}

fn chemical_component_bond_atom_candidates(
    atoms: &[Atom],
    atom_indices: &[usize],
    comp_id: &str,
    chem_atom_id: &str,
) -> Vec<usize> {
    atom_indices
        .iter()
        .copied()
        .filter(|&index| {
            atoms.get(index).is_some_and(|atom| {
                atom_matches_chemical_component_bond(comp_id, atom, chem_atom_id)
            })
        })
        .collect()
}

fn chemical_component_bond_alt_locs_are_compatible(a: &Atom, b: &Atom) -> bool {
    a.alt_id.is_empty() || b.alt_id.is_empty() || a.alt_id == b.alt_id
}

fn atom_matches_chemical_component_bond(comp_id: &str, atom: &Atom, chem_atom_id: &str) -> bool {
    if atom.name.eq_ignore_ascii_case(chem_atom_id) {
        return true;
    }
    if comp_id.eq_ignore_ascii_case("DOD") || !is_hydrogen_like_atom(atom) {
        return false;
    }
    let Some(rest) = atom.name.strip_prefix('D') else {
        return false;
    };
    let normalized = format!("H{rest}");
    normalized.eq_ignore_ascii_case(chem_atom_id)
}

fn is_hydrogen_like_atom(atom: &Atom) -> bool {
    atom.element.eq_ignore_ascii_case("H") || atom.element.eq_ignore_ascii_case("D")
}

fn normalized_bond_pair(a: usize, b: usize) -> (usize, usize) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

pub(super) fn assign_intra_bond_orders(molecule: &mut Molecule) {
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        let (Some(a), Some(b)) = (molecule.atoms.get(bond.a), molecule.atoms.get(bond.b)) else {
            continue;
        };
        if a.chain != b.chain
            || a.residue != b.residue
            || a.residue_seq != b.residue_seq
            || a.insertion_code != b.insertion_code
        {
            continue;
        }
        let order = intra_bond_order_from_table(&a.residue, &a.name, &b.name);
        if let Some(metadata) = molecule.bond_metadata.get_mut(bond_index) {
            metadata.order = metadata.order.max(order);
        }
    }
}

pub(super) fn detect_rings(molecule: &Molecule) -> Vec<Ring> {
    let graph = RingBondGraph::new(molecule);
    let residue_ranges = molecule_residue_ranges(molecule);
    let largest_residue = residue_ranges
        .iter()
        .map(|range| range.end - range.start)
        .max()
        .unwrap_or(0);
    let mut state = RingSearchState::new(largest_residue, &graph, molecule);

    for range in residue_ranges {
        state.process_residue(range.start, range.end);
    }

    state
        .rings
        .into_iter()
        .map(|atom_indices| {
            let aromatic = ring_is_aromatic(molecule, &graph, &atom_indices);
            Ring {
                aromatic,
                fingerprint: ring_fingerprint(molecule, &atom_indices),
                atom_indices,
            }
        })
        .collect()
}

fn ring_fingerprint(molecule: &Molecule, ring_order: &[usize]) -> String {
    let elements = ring_order
        .iter()
        .map(|index| {
            molecule.atoms.get(*index).map_or_else(String::new, |atom| {
                if atom.type_symbol.is_empty() {
                    atom.element.clone()
                } else {
                    atom.type_symbol.clone()
                }
            })
        })
        .collect::<Vec<_>>();
    molstar_ring_fingerprint(&elements)
}

fn molstar_ring_fingerprint(elements: &[String]) -> String {
    let len = elements.len();
    if len == 0 {
        return String::new();
    }
    let reversed = elements.iter().rev().cloned().collect::<Vec<_>>();
    let rot_normal = minimal_rotation(elements);
    let rot_reversed = minimal_rotation(&reversed);
    let mut normal_smaller = false;
    for i in 0..len {
        let normal = &elements[(i + rot_normal) % len];
        let reverse = &reversed[(i + rot_reversed) % len];
        if normal != reverse {
            normal_smaller = normal < reverse;
            break;
        }
    }
    if normal_smaller {
        build_ring_fingerprint(elements, rot_normal)
    } else {
        build_ring_fingerprint(&reversed, rot_reversed)
    }
}

fn minimal_rotation(elements: &[String]) -> usize {
    let len = elements.len();
    if len <= 1 {
        return 0;
    }

    let mut failure = vec![-1isize; len * 2];
    let mut start = 0usize;

    for cursor in 1..failure.len() {
        let mut matched = failure[cursor - start - 1];
        while matched != -1 {
            let a = &elements[cursor % len];
            let b = &elements[(start + matched as usize + 1) % len];
            if a == b {
                break;
            }
            if a < b {
                start = cursor - matched as usize - 1;
            }
            matched = failure[matched as usize];
        }

        let slot = cursor - start;
        if matched == -1 {
            let a = &elements[cursor % len];
            let b = &elements[start % len];
            if a != b {
                if a < b {
                    start = cursor;
                }
                failure[cursor - start] = -1;
            } else {
                failure[slot] = 0;
            }
        } else {
            failure[slot] = matched + 1;
        }
    }

    start
}

fn build_ring_fingerprint(elements: &[String], offset: usize) -> String {
    (0..elements.len())
        .map(|i| elements[(i + offset) % elements.len()].as_str())
        .collect::<Vec<_>>()
        .join("-")
}

const MOLSTAR_RING_MAX_DEPTH: usize = 5;

#[derive(Clone, Copy)]
struct RingResidueRange {
    start: usize,
    end: usize,
}

fn molecule_residue_ranges(molecule: &Molecule) -> Vec<RingResidueRange> {
    let mut ranges = Vec::new();
    let mut start = 0usize;
    while start < molecule.atoms.len() {
        let mut end = start + 1;
        while end < molecule.atoms.len()
            && atoms_are_same_molstar_residue(&molecule.atoms[start], &molecule.atoms[end])
        {
            end += 1;
        }
        ranges.push(RingResidueRange { start, end });
        start = end;
    }
    ranges
}

fn atoms_are_same_molstar_residue(a: &Atom, b: &Atom) -> bool {
    a.model_num == b.model_num
        && a.chain == b.chain
        && a.entity_id == b.entity_id
        && a.residue_seq == b.residue_seq
        && a.auth_residue_seq == b.auth_residue_seq
        && a.insertion_code == b.insertion_code
}

struct RingBondGraph {
    offset: Vec<usize>,
    b: Vec<usize>,
    flags: Vec<BondFlags>,
}

impl RingBondGraph {
    fn new(molecule: &Molecule) -> Self {
        let atom_count = molecule.atoms.len();
        let mut bucket_sizes = vec![0usize; atom_count];
        let valid_edges = molecule
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(bond_index, bond)| {
                (bond.a < atom_count && bond.b < atom_count).then_some((bond_index, bond))
            })
            .collect::<Vec<_>>();
        for (_, bond) in &valid_edges {
            bucket_sizes[bond.a] += 1;
            bucket_sizes[bond.b] += 1;
        }

        let mut offset = vec![0usize; atom_count + 1];
        let mut cursor = 0usize;
        for (index, size) in bucket_sizes.iter().enumerate() {
            offset[index] = cursor;
            cursor += *size;
        }
        offset[atom_count] = cursor;

        let mut graph = RingBondGraph {
            offset,
            b: vec![0; cursor],
            flags: vec![BondFlags::NONE; cursor],
        };
        let mut bucket_fill = vec![0usize; atom_count];
        for (bond_index, bond) in valid_edges {
            let metadata = molecule
                .bond_metadata
                .get(bond_index)
                .cloned()
                .unwrap_or_default();
            let slot_ab = graph.offset[bond.a] + bucket_fill[bond.a];
            bucket_fill[bond.a] += 1;
            graph.b[slot_ab] = bond.b;
            graph.flags[slot_ab] = metadata.flags;

            let slot_ba = graph.offset[bond.b] + bucket_fill[bond.b];
            bucket_fill[bond.b] += 1;
            graph.b[slot_ba] = bond.a;
            graph.flags[slot_ba] = metadata.flags;
        }
        graph
    }

    fn bond_count(&self, atom_index: usize) -> usize {
        self.offset
            .get(atom_index + 1)
            .zip(self.offset.get(atom_index))
            .map(|(end, start)| end - start)
            .unwrap_or(0)
    }
}

struct RingSearchState<'a> {
    start_vertex: usize,
    end_vertex: usize,
    count: usize,
    is_ring_atom: Vec<i32>,
    marked: Vec<i32>,
    queue: Vec<usize>,
    color: Vec<i32>,
    pred: Vec<isize>,
    depth: Vec<usize>,
    left: [usize; MOLSTAR_RING_MAX_DEPTH],
    right: [usize; MOLSTAR_RING_MAX_DEPTH],
    current_color: i32,
    current_alt_loc: String,
    has_alt_loc: bool,
    rings: Vec<Vec<usize>>,
    current_rings: Vec<Vec<usize>>,
    graph: &'a RingBondGraph,
    molecule: &'a Molecule,
}

impl<'a> RingSearchState<'a> {
    fn new(capacity: usize, graph: &'a RingBondGraph, molecule: &'a Molecule) -> Self {
        RingSearchState {
            start_vertex: 0,
            end_vertex: 0,
            count: 0,
            is_ring_atom: vec![0; capacity],
            marked: vec![0; capacity],
            queue: vec![0; capacity],
            color: vec![0; capacity],
            pred: vec![-1; capacity],
            depth: vec![0; capacity],
            left: [0; MOLSTAR_RING_MAX_DEPTH],
            right: [0; MOLSTAR_RING_MAX_DEPTH],
            current_color: 0,
            current_alt_loc: String::new(),
            has_alt_loc: false,
            rings: Vec::new(),
            current_rings: Vec::new(),
            graph,
            molecule,
        }
    }

    fn process_residue(&mut self, start: usize, end: usize) {
        self.start_vertex = start;
        self.end_vertex = end;
        if end.saturating_sub(start) < 3 {
            return;
        }

        self.current_rings.clear();
        let mut alt_locs = Vec::<String>::new();
        for atom_index in start..end {
            let alt_loc = &self.molecule.atoms[atom_index].alt_id;
            if !alt_loc.is_empty() && !alt_locs.iter().any(|entry| entry == alt_loc) {
                alt_locs.push(alt_loc.clone());
            }
        }

        let mut mark = 1i32;
        if alt_locs.is_empty() {
            self.reset_state();
            for i in 0..self.count {
                if !self.is_start_index(i) {
                    continue;
                }
                self.reset_depth();
                mark = self.find_rings(i, mark);
            }
        } else {
            for alt_loc in alt_locs {
                self.reset_state();
                self.has_alt_loc = true;
                self.current_alt_loc = alt_loc;
                for i in 0..self.count {
                    if !self.is_start_index(i) {
                        continue;
                    }
                    let atom_alt_loc = &self.molecule.atoms[self.start_vertex + i].alt_id;
                    if !atom_alt_loc.is_empty() && atom_alt_loc != &self.current_alt_loc {
                        continue;
                    }
                    self.reset_depth();
                    mark = self.find_rings(i, mark);
                }
            }
        }

        self.rings.extend(self.current_rings.iter().cloned());
    }

    fn reset_state(&mut self) {
        self.count = self.end_vertex - self.start_vertex;
        self.is_ring_atom[..self.count].fill(0);
        self.marked[..self.count].fill(-1);
        self.color[..self.count].fill(0);
        self.pred[..self.count].fill(-1);
        self.depth[..self.count].fill(0);
        self.current_color = 0;
        self.current_alt_loc.clear();
        self.has_alt_loc = false;
    }

    fn reset_depth(&mut self) {
        self.depth[..self.count].fill(self.count + 1);
    }

    fn is_start_index(&self, i: usize) -> bool {
        let atom_index = self.start_vertex + i;
        let bond_count = self.graph.bond_count(atom_index);
        if bond_count <= 1 || (self.is_ring_atom[i] != 0 && bond_count == 2) {
            return false;
        }
        true
    }

    fn find_rings(&mut self, from: usize, mark: i32) -> i32 {
        self.marked[from] = mark;
        self.depth[from] = 0;
        self.queue[0] = from;
        let mut head = 0usize;
        let mut size = 1usize;

        while head < size {
            let top = self.queue[head];
            head += 1;
            let d = self.depth[top];
            let atom_index = self.start_vertex + top;
            let start = self.graph.offset[atom_index];
            let end = self.graph.offset[atom_index + 1];

            for slot in start..end {
                let neighbor = self.graph.b[slot];
                if neighbor < self.start_vertex
                    || neighbor >= self.end_vertex
                    || !self.graph.flags[slot].contains(BondFlags::COVALENT)
                {
                    continue;
                }

                if self.has_alt_loc {
                    let alt_loc = &self.molecule.atoms[neighbor].alt_id;
                    if !alt_loc.is_empty() && alt_loc != &self.current_alt_loc {
                        continue;
                    }
                }

                let other = neighbor - self.start_vertex;
                if self.marked[other] == mark {
                    if self.pred[other] != top as isize
                        && self.pred[top] != other as isize
                        && self.add_ring(top, other)
                    {
                        return mark + 1;
                    }
                    continue;
                }

                let new_depth = self.depth[other].min(d + 1);
                if new_depth > MOLSTAR_RING_MAX_DEPTH {
                    continue;
                }
                self.depth[other] = new_depth;
                self.marked[other] = mark;
                self.queue[size] = other;
                size += 1;
                self.pred[other] = top as isize;
            }
        }
        mark + 1
    }

    fn add_ring(&mut self, a: usize, b: usize) -> bool {
        if b < a {
            return false;
        }

        self.current_color += 1;
        let color = self.current_color;
        let mut current = a as isize;
        for _ in 0..MOLSTAR_RING_MAX_DEPTH {
            if current < 0 {
                break;
            }
            self.color[current as usize] = color;
            current = self.pred[current as usize];
        }

        let mut left_offset = 0usize;
        let mut right_offset = 0usize;
        let mut found = false;
        let mut target = 0usize;
        current = b as isize;
        for _ in 0..MOLSTAR_RING_MAX_DEPTH {
            if current < 0 {
                break;
            }
            let current_index = current as usize;
            if self.color[current_index] == color {
                target = current_index;
                found = true;
                break;
            }
            self.right[right_offset] = current_index;
            right_offset += 1;
            current = self.pred[current_index];
        }
        if !found {
            return false;
        }

        current = a as isize;
        for _ in 0..MOLSTAR_RING_MAX_DEPTH {
            if current < 0 {
                break;
            }
            let current_index = current as usize;
            self.left[left_offset] = current_index;
            left_offset += 1;
            if target == current_index {
                break;
            }
            current = self.pred[current_index];
        }

        let len = left_offset + right_offset;
        if len < 3 {
            return false;
        }

        let mut ring = Vec::with_capacity(len);
        for index in 0..left_offset {
            let atom_index = self.start_vertex + self.left[index];
            ring.push(atom_index);
            self.is_ring_atom[self.left[index]] = 1;
        }
        for index in (0..right_offset).rev() {
            let atom_index = self.start_vertex + self.right[index];
            ring.push(atom_index);
            self.is_ring_atom[self.right[index]] = 1;
        }
        ring.sort_unstable();

        for existing in &self.current_rings {
            if ring.len() == existing.len() {
                if ring == *existing {
                    return false;
                }
            } else if ring.len() > existing.len() && sorted_slice_contains_subset(&ring, existing) {
                return false;
            }
        }

        self.current_rings.push(ring);
        true
    }
}

fn sorted_slice_contains_subset(a: &[usize], b: &[usize]) -> bool {
    let mut i = 0usize;
    let mut j = 0usize;
    let mut equal = 0usize;
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            i += 1;
        } else if a[i] > b[j] {
            j += 1;
        } else {
            i += 1;
            j += 1;
            equal += 1;
        }
    }
    equal == b.len()
}

fn ring_is_aromatic(molecule: &Molecule, graph: &RingBondGraph, atom_indices: &[usize]) -> bool {
    let Some(first) = atom_indices
        .first()
        .and_then(|index| molecule.atoms.get(*index))
    else {
        return false;
    };
    let residue = first.residue.to_ascii_uppercase();
    if residue == "PRO" {
        return false;
    }

    let mut aromatic_bond_count = 0usize;
    let mut has_aromatic_ring_element = false;
    for &atom_index in atom_indices {
        if !has_aromatic_ring_element
            && molecule
                .atoms
                .get(atom_index)
                .is_some_and(atom_can_be_aromatic)
        {
            has_aromatic_ring_element = true;
        }
        let Some(start) = graph.offset.get(atom_index).copied() else {
            continue;
        };
        let Some(end) = graph.offset.get(atom_index + 1).copied() else {
            continue;
        };
        for slot in start..end {
            if graph.flags[slot].contains(BondFlags::AROMATIC)
                && atom_indices.binary_search(&graph.b[slot]).is_ok()
            {
                aromatic_bond_count += 1;
            }
        }
    }

    if aromatic_bond_count == 2 * atom_indices.len() {
        return true;
    }
    if !has_aromatic_ring_element {
        return false;
    }
    if atom_indices.len() < 5 {
        return false;
    }
    if aromatic_bond_count > 0 {
        return false;
    }
    ring_planarity_deviation(molecule, atom_indices) < 0.05
}

fn atom_can_be_aromatic(atom: &Atom) -> bool {
    matches!(
        aromatic_element_symbol(atom).as_str(),
        "B" | "C" | "N" | "O" | "SI" | "P" | "S" | "GE" | "AS" | "SN" | "SB" | "BI"
    )
}

fn aromatic_element_symbol(atom: &Atom) -> String {
    if atom.type_symbol.is_empty() {
        atom.element.to_ascii_uppercase()
    } else {
        atom.type_symbol.to_ascii_uppercase()
    }
}

fn ring_planarity_deviation(molecule: &Molecule, atom_indices: &[usize]) -> f32 {
    if atom_indices.len() <= 1 {
        return 0.0;
    }
    let positions = atom_indices
        .iter()
        .filter_map(|index| molecule.atoms.get(*index).map(|atom| atom.position))
        .collect::<Vec<_>>();
    if positions.len() != atom_indices.len() {
        return f32::INFINITY;
    }

    let center = positions
        .iter()
        .copied()
        .fold(Vec3::default(), |sum, position| sum + position)
        / positions.len() as f32;
    let mut xx = 0.0;
    let mut xy = 0.0;
    let mut xz = 0.0;
    let mut yy = 0.0;
    let mut yz = 0.0;
    let mut zz = 0.0;
    for position in positions {
        let delta = position - center;
        xx += delta.x * delta.x;
        xy += delta.x * delta.y;
        xz += delta.x * delta.z;
        yy += delta.y * delta.y;
        yz += delta.y * delta.z;
        zz += delta.z * delta.z;
    }
    let lambda = smallest_symmetric_3x3_eigenvalue(xx, xy, xz, yy, yz, zz).max(0.0);
    (lambda / atom_indices.len() as f32).sqrt()
}

fn smallest_symmetric_3x3_eigenvalue(xx: f32, xy: f32, xz: f32, yy: f32, yz: f32, zz: f32) -> f32 {
    let p1 = xy * xy + xz * xz + yz * yz;
    if p1 <= f32::EPSILON {
        return xx.min(yy).min(zz);
    }

    let trace = xx + yy + zz;
    let q = trace / 3.0;
    let axx = xx - q;
    let ayy = yy - q;
    let azz = zz - q;
    let p2 = axx * axx + ayy * ayy + azz * azz + 2.0 * p1;
    let p = (p2 / 6.0).sqrt();
    if p <= f32::EPSILON {
        return q;
    }

    let bxx = axx / p;
    let bxy = xy / p;
    let bxz = xz / p;
    let byy = ayy / p;
    let byz = yz / p;
    let bzz = azz / p;
    let det_b = bxx * (byy * bzz - byz * byz) - bxy * (bxy * bzz - byz * bxz)
        + bxz * (bxy * byz - byy * bxz);
    let r = (det_b / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let eig1 = q + 2.0 * p * phi.cos();
    let eig3 = q + 2.0 * p * (phi + 2.0 * std::f32::consts::PI / 3.0).cos();
    let eig2 = 3.0 * q - eig1 - eig3;
    eig1.min(eig2).min(eig3)
}

pub(super) fn assign_ring_resonance(molecule: &mut Molecule) {
    let membership = RingMembership::new(&molecule.rings, molecule.atoms.len());
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        if membership.atoms_share_aromatic_ring(bond.a, bond.b) {
            if let Some(metadata) = molecule.bond_metadata.get_mut(bond_index) {
                if !metadata.flags.contains(BondFlags::RESONANCE) {
                    molecule.derived_resonance_bonds.insert(bond_index);
                    metadata.flags = metadata.flags.union(BondFlags::RESONANCE);
                }
                if !metadata.flags.contains(BondFlags::AROMATIC) {
                    molecule.derived_aromatic_bonds.insert(bond_index);
                    metadata.flags = metadata.flags.union(BondFlags::AROMATIC);
                }
            }
        }
    }
}

pub(super) fn build_resonance(molecule: &Molecule) -> Resonance {
    let membership = RingMembership::new(&molecule.rings, molecule.atoms.len());
    let delocalized_triplet_lookup = build_delocalized_triplets(molecule, &membership);

    let mut parent = (0..molecule.rings.len()).collect::<Vec<_>>();
    for indices in &membership.element_ring_indices {
        for window in indices.windows(2) {
            union_components(&mut parent, window[0], window[1]);
        }
    }
    let mut component_roots = Vec::<usize>::new();
    let mut ring_component_index = Vec::with_capacity(molecule.rings.len());
    for ring_index in 0..molecule.rings.len() {
        let root = find_component(&mut parent, ring_index);
        let component_index = component_roots
            .iter()
            .position(|existing| *existing == root)
            .unwrap_or_else(|| {
                component_roots.push(root);
                component_roots.len() - 1
            });
        ring_component_index.push(component_index);
    }
    let mut ring_components = vec![Vec::new(); component_roots.len()];
    for (ring_index, component_index) in ring_component_index.iter().enumerate() {
        ring_components[*component_index].push(ring_index);
    }

    Resonance {
        ring_count: molecule.rings.len(),
        aromatic_ring_count: molecule.rings.iter().filter(|ring| ring.aromatic).count(),
        delocalized_bond_count: molecule
            .bond_metadata
            .iter()
            .filter(|metadata| {
                metadata.flags.contains(BondFlags::AROMATIC)
                    || metadata.flags.contains(BondFlags::RESONANCE)
            })
            .count(),
        delocalized_triplets: delocalized_triplet_lookup.triplets.clone(),
        delocalized_triplet_lookup,
        element_ring_indices: membership.element_ring_indices,
        element_aromatic_ring_indices: membership.element_aromatic_ring_indices,
        ring_component_index,
        ring_components,
    }
}

fn build_delocalized_triplets(
    molecule: &Molecule,
    membership: &RingMembership,
) -> DelocalizedTriplets {
    let mut aromatic_neighbors = vec![Vec::<usize>::new(); molecule.atoms.len()];
    for (bond_index, bond) in molecule.bonds.iter().enumerate() {
        let Some(metadata) = molecule.bond_metadata.get(bond_index) else {
            continue;
        };
        if !metadata.flags.contains(BondFlags::AROMATIC) {
            continue;
        }
        if bond.a < molecule.atoms.len() && bond.b < molecule.atoms.len() {
            aromatic_neighbors[bond.a].push(bond.b);
            aromatic_neighbors[bond.b].push(bond.a);
        }
    }

    let mut triplets = DelocalizedTriplets::default();
    for (center, neighbors) in aromatic_neighbors.iter().enumerate() {
        if membership
            .element_aromatic_ring_indices
            .get(center)
            .is_some_and(|rings| !rings.is_empty())
        {
            continue;
        }
        if neighbors.len() >= 2 {
            triplets.add(center, neighbors[0], neighbors[1]);
            for index in 1..neighbors.len() {
                triplets.add(center, neighbors[index], neighbors[0]);
            }
        }
    }
    triplets
}

fn sorted_triplet(a: usize, b: usize, c: usize) -> [usize; 3] {
    let mut triplet = [a, b, c];
    triplet.sort_unstable();
    triplet
}

fn sorted_cantor_pairing(a: usize, b: usize) -> usize {
    if a < b {
        cantor_pairing(a, b)
    } else {
        cantor_pairing(b, a)
    }
}

fn cantor_pairing(a: usize, b: usize) -> usize {
    (a + b) * (a + b + 1) / 2 + b
}

#[derive(Debug, Default)]
struct RingMembership {
    element_ring_indices: Vec<Vec<usize>>,
    element_aromatic_ring_indices: Vec<Vec<usize>>,
}

impl RingMembership {
    fn new(rings: &[Ring], atom_count: usize) -> Self {
        let mut element_ring_indices = vec![Vec::new(); atom_count];
        let mut element_aromatic_ring_indices = vec![Vec::new(); atom_count];

        for (ring_index, ring) in rings.iter().enumerate() {
            for &atom_index in &ring.atom_indices {
                if let Some(indices) = element_ring_indices.get_mut(atom_index) {
                    indices.push(ring_index);
                }
                if ring.aromatic {
                    if let Some(indices) = element_aromatic_ring_indices.get_mut(atom_index) {
                        indices.push(ring_index);
                    }
                }
            }
        }

        Self {
            element_ring_indices,
            element_aromatic_ring_indices,
        }
    }

    fn atoms_share_aromatic_ring(&self, a: usize, b: usize) -> bool {
        let (Some(a_rings), Some(b_rings)) = (
            self.element_aromatic_ring_indices.get(a),
            self.element_aromatic_ring_indices.get(b),
        ) else {
            return false;
        };
        a_rings
            .iter()
            .any(|ring_index| b_rings.contains(ring_index))
    }
}

fn find_component(parent: &mut [usize], index: usize) -> usize {
    if parent[index] != index {
        parent[index] = find_component(parent, parent[index]);
    }
    parent[index]
}

fn union_components(parent: &mut [usize], a: usize, b: usize) {
    let root_a = find_component(parent, a);
    let root_b = find_component(parent, b);
    if root_a != root_b {
        parent[root_b] = root_a;
    }
}

pub(crate) fn intra_bond_order_from_table(comp_id: &str, atom_id_1: &str, atom_id_2: &str) -> i8 {
    let comp_id = comp_id.to_ascii_uppercase();
    let mut atom_id_1 = atom_id_1.to_ascii_uppercase();
    let mut atom_id_2 = atom_id_2.to_ascii_uppercase();
    if atom_id_1 > atom_id_2 {
        std::mem::swap(&mut atom_id_1, &mut atom_id_2);
    }
    if is_protein_residue(&comp_id) && atom_id_1 == "C" && atom_id_2 == "O" {
        return 2;
    }
    if (is_rna_residue(&comp_id) || is_dna_residue(&comp_id))
        && atom_id_1 == "OP1"
        && atom_id_2 == "P"
    {
        return 2;
    }
    match (comp_id.as_str(), atom_id_1.as_str(), atom_id_2.as_str()) {
        ("HIS", "CD2", "CG") | ("HIS", "CE1", "ND1") | ("ARG", "CZ", "NH2") => 2,
        ("PHE", "CD1", "CG") | ("PHE", "CD2", "CE2") | ("PHE", "CE1", "CZ") => 2,
        ("TRP", "CD1", "CG")
        | ("TRP", "CD2", "CE2")
        | ("TRP", "CE3", "CZ3")
        | ("TRP", "CH2", "CZ2") => 2,
        ("ASN", "CG", "OD1") | ("GLN", "CD", "OE1") => 2,
        ("TYR", "CD1", "CG") | ("TYR", "CD2", "CE2") | ("TYR", "CE1", "CZ") => 2,
        ("ASP", "CG", "OD1") | ("GLU", "CD", "OE1") => 2,
        ("G" | "DG", "C2", "N3")
        | ("G" | "DG", "C4", "C5")
        | ("G" | "DG", "C6", "O6")
        | ("G" | "DG", "C8", "N7")
        | ("C" | "DC", "C2", "O2")
        | ("C" | "DC", "C4", "N3")
        | ("C" | "DC", "C5", "C6")
        | ("A" | "DA", "C2", "N3")
        | ("A" | "DA", "C4", "C5")
        | ("A" | "DA", "C6", "N1")
        | ("A" | "DA", "C8", "N7")
        | ("U" | "DT", "C2", "O2")
        | ("U" | "DT", "C4", "O4")
        | ("U" | "DT", "C5", "C6") => 2,
        _ => 1,
    }
}
