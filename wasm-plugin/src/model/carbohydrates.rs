use std::collections::{BTreeMap, BTreeSet};

use super::{
    saccharide_component, AtomicStructure, BondFlags, Molecule, SaccharideComponent, StructureUnit,
    UnitKind, Vec3,
};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CarbohydrateLink {
    pub carbohydrate_index_a: usize,
    pub carbohydrate_index_b: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CarbohydrateTerminalLink {
    pub carbohydrate_index: usize,
    pub element_unit_id: usize,
    pub element_index: usize,
    pub from_carbohydrate: bool,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CarbohydrateSymbolGeometry {
    pub center: Vec3,
    pub normal: Vec3,
    pub direction: Vec3,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CarbohydrateElement {
    pub geometry: CarbohydrateSymbolGeometry,
    pub unit_id: usize,
    pub residue_index: usize,
    pub component: SaccharideComponent,
    pub ring_index: usize,
    pub ring_element_indices: Vec<usize>,
    pub alt_id: String,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PartialCarbohydrateElement {
    pub unit_id: usize,
    pub residue_index: usize,
    pub component: SaccharideComponent,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Carbohydrates {
    pub links: Vec<CarbohydrateLink>,
    pub terminal_links: Vec<CarbohydrateTerminalLink>,
    pub elements: Vec<CarbohydrateElement>,
    pub partial_elements: Vec<PartialCarbohydrateElement>,
    element_indices: BTreeMap<(usize, usize), Vec<usize>>,
    link_indices: BTreeMap<(usize, usize), Vec<usize>>,
    terminal_link_indices: BTreeMap<(usize, usize), Vec<usize>>,
}

impl Carbohydrates {
    pub(crate) fn from_structure(molecule: &Molecule, structure: &AtomicStructure) -> Self {
        if !molecule
            .atoms
            .iter()
            .any(|atom| saccharide_component(&atom.residue).is_some())
        {
            return Carbohydrates::default();
        }

        let mut elements = Vec::<CarbohydrateElement>::new();
        let mut partial_elements = Vec::<PartialCarbohydrateElement>::new();
        let mut links = Vec::<CarbohydrateLink>::new();
        let mut terminal_links = Vec::<CarbohydrateTerminalLink>::new();
        let mut elements_with_ring_map = BTreeMap::<(usize, usize, String), Vec<usize>>::new();

        for unit in &structure.units {
            if unit.kind != UnitKind::Atomic {
                continue;
            }
            for &residue_index in &unit.residue_indices {
                let Some(residue) = structure.model.hierarchy.residues.get(residue_index) else {
                    continue;
                };
                let Some(component) = saccharide_component(&residue.comp_id) else {
                    continue;
                };
                let sugar_rings = filter_fused_rings(
                    molecule,
                    unit,
                    sugar_ring_indices_for_residue(molecule, unit, residue_index),
                );
                if sugar_rings.is_empty() {
                    partial_elements.push(PartialCarbohydrateElement {
                        unit_id: unit.id,
                        residue_index,
                        component,
                    });
                    continue;
                }

                let mut ring_elements = Vec::new();
                for ring_index in sugar_rings {
                    let Some(ring_element_indices) = ring_unit_indices(molecule, unit, ring_index)
                    else {
                        continue;
                    };
                    let alt_id = ring_alt_id(molecule, unit, &ring_element_indices);
                    let element_index = elements.len();
                    ring_elements.push((ring_index, element_index, ring_element_indices.clone()));
                    add_ring_element(
                        &mut elements_with_ring_map,
                        residue_index,
                        unit.id,
                        &alt_id,
                        element_index,
                    );
                    if !alt_id.is_empty() {
                        add_ring_element(
                            &mut elements_with_ring_map,
                            residue_index,
                            unit.id,
                            "",
                            element_index,
                        );
                    }
                    let geometry = carbohydrate_symbol_geometry(
                        structure,
                        molecule,
                        unit,
                        &ring_element_indices,
                    );
                    elements.push(CarbohydrateElement {
                        geometry,
                        unit_id: unit.id,
                        residue_index,
                        component: component.clone(),
                        ring_index,
                        ring_element_indices,
                        alt_id,
                    });
                }

                for a in 0..ring_elements.len() {
                    for b in (a + 1)..ring_elements.len() {
                        let (_, element_a, ring_a) = &ring_elements[a];
                        let (_, element_b, ring_b) = &ring_elements[b];
                        if elements[*element_a].alt_id == elements[*element_b].alt_id
                            && are_vertex_sets_connected(
                                unit,
                                ring_a,
                                ring_b,
                                MOLSTAR_CARBOHYDRATE_LINK_MAX_DISTANCE,
                            )
                        {
                            fix_link_direction(&mut elements, *element_a, *element_b);
                            fix_link_direction(&mut elements, *element_b, *element_a);
                            links.push(CarbohydrateLink {
                                carbohydrate_index_a: *element_a,
                                carbohydrate_index_b: *element_b,
                            });
                            links.push(CarbohydrateLink {
                                carbohydrate_index_a: *element_b,
                                carbohydrate_index_b: *element_a,
                            });
                        }
                    }
                }
            }
        }

        for index_a in 0..elements.len() {
            let unit_a_id = elements[index_a].unit_id;
            let Some(unit) = structure.unit_by_id(unit_a_id) else {
                continue;
            };
            for index_b in (index_a + 1)..elements.len() {
                if elements[index_b].unit_id != unit_a_id
                    || elements[index_a].residue_index == elements[index_b].residue_index
                {
                    continue;
                }
                if are_vertex_sets_connected(
                    unit,
                    &elements[index_a].ring_element_indices,
                    &elements[index_b].ring_element_indices,
                    MOLSTAR_CARBOHYDRATE_LINK_MAX_DISTANCE,
                ) {
                    fix_link_direction(&mut elements, index_a, index_b);
                    fix_link_direction(&mut elements, index_b, index_a);
                    links.push(CarbohydrateLink {
                        carbohydrate_index_a: index_a,
                        carbohydrate_index_b: index_b,
                    });
                    links.push(CarbohydrateLink {
                        carbohydrate_index_a: index_b,
                        carbohydrate_index_b: index_a,
                    });
                }
            }
        }

        for unit in &structure.units {
            if unit.kind != UnitKind::Atomic {
                continue;
            }
            for pair_bonds in structure.inter_unit_bond_graph.get_connected_units(unit.id) {
                for &index_a in &pair_bonds.connected_indices {
                    for edge in pair_bonds.get_edges(index_a) {
                        if !edge.props.flag.contains(BondFlags::COVALENT) {
                            continue;
                        }
                        let Some(unit_a) = structure.unit_by_id(pair_bonds.unit_a) else {
                            continue;
                        };
                        let Some(unit_b) = structure.unit_by_id(pair_bonds.unit_b) else {
                            continue;
                        };
                        let index_b = edge.index_b;
                        let ring_element_indices_a = ring_element_indices_for_unit_index(
                            &elements_with_ring_map,
                            molecule,
                            unit_a,
                            index_a,
                        );
                        let ring_element_indices_b = ring_element_indices_for_unit_index(
                            &elements_with_ring_map,
                            molecule,
                            unit_b,
                            index_b,
                        );
                        if !ring_element_indices_a.is_empty() && !ring_element_indices_b.is_empty()
                        {
                            let len_a = ring_element_indices_a.len();
                            let len_b = ring_element_indices_b.len();
                            let atom_id_a = atom_name_for_unit_index(molecule, unit_a, index_a);
                            for index in 0..len_a.max(len_b) {
                                let carbohydrate_index_a =
                                    ring_element_indices_a[index.min(len_a - 1)];
                                let carbohydrate_index_b =
                                    ring_element_indices_b[index.min(len_b - 1)];
                                if atom_id_a.starts_with("O1") || atom_id_a.starts_with("C1") {
                                    fix_link_direction(
                                        &mut elements,
                                        carbohydrate_index_a,
                                        carbohydrate_index_b,
                                    );
                                }
                                links.push(CarbohydrateLink {
                                    carbohydrate_index_a,
                                    carbohydrate_index_b,
                                });
                            }
                        } else if ring_element_indices_b.is_empty() {
                            let atom_id_a = atom_name_for_unit_index(molecule, unit_a, index_a);
                            for carbohydrate_index in ring_element_indices_a {
                                if atom_id_a.starts_with("O1") || atom_id_a.starts_with("C1") {
                                    fix_terminal_link_direction(
                                        structure,
                                        &mut elements,
                                        carbohydrate_index,
                                        unit_b.id,
                                        index_b,
                                    );
                                }
                                terminal_links.push(CarbohydrateTerminalLink {
                                    carbohydrate_index,
                                    element_unit_id: unit_b.id,
                                    element_index: index_b,
                                    from_carbohydrate: true,
                                });
                            }
                        } else if ring_element_indices_a.is_empty() {
                            for carbohydrate_index in ring_element_indices_b {
                                terminal_links.push(CarbohydrateTerminalLink {
                                    carbohydrate_index,
                                    element_unit_id: unit_a.id,
                                    element_index: index_a,
                                    from_carbohydrate: false,
                                });
                            }
                        }
                    }
                }
            }
        }

        let mut carbohydrates = Carbohydrates {
            links,
            terminal_links,
            elements,
            partial_elements,
            ..Carbohydrates::default()
        };
        carbohydrates.build_lookups(structure);
        carbohydrates
    }

    pub fn get_element_indices(&self, unit_id: usize, source_atom_index: usize) -> &[usize] {
        self.element_indices
            .get(&(unit_id, source_atom_index))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn get_link_indices(&self, unit_id: usize, source_atom_index: usize) -> &[usize] {
        self.link_indices
            .get(&(unit_id, source_atom_index))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn get_terminal_link_indices(&self, unit_id: usize, source_atom_index: usize) -> &[usize] {
        self.terminal_link_indices
            .get(&(unit_id, source_atom_index))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    fn build_lookups(&mut self, structure: &AtomicStructure) {
        for (element_index, element) in self.elements.iter().enumerate() {
            for source_atom in ring_source_atoms(structure, element) {
                add_unique_map_index(
                    &mut self.element_indices,
                    element.unit_id,
                    source_atom,
                    element_index,
                );
            }
        }

        for (link_index, link) in self.links.iter().enumerate() {
            let Some(element) = self.elements.get(link.carbohydrate_index_a) else {
                continue;
            };
            for source_atom in ring_source_atoms(structure, element) {
                add_unique_map_index(
                    &mut self.link_indices,
                    element.unit_id,
                    source_atom,
                    link_index,
                );
            }
        }

        for (terminal_index, link) in self.terminal_links.iter().enumerate() {
            if link.from_carbohydrate {
                let Some(element) = self.elements.get(link.carbohydrate_index) else {
                    continue;
                };
                for source_atom in ring_source_atoms(structure, element) {
                    add_unique_map_index(
                        &mut self.terminal_link_indices,
                        element.unit_id,
                        source_atom,
                        terminal_index,
                    );
                }
            } else if let Some(source_atom) = structure
                .units
                .get(link.element_unit_id)
                .and_then(|unit| unit.atom_indices.get(link.element_index))
                .copied()
            {
                add_unique_map_index(
                    &mut self.terminal_link_indices,
                    link.element_unit_id,
                    source_atom,
                    terminal_index,
                );
            }
        }
    }
}

const MOLSTAR_CARBOHYDRATE_LINK_MAX_DISTANCE: usize = 3;
const SUGAR_RING_FINGERPRINTS: &[&str] = &["C-C-C-O", "C-C-C-C-O", "C-C-C-C-C-O", "C-C-C-C-C-C-O"];

fn carbohydrate_symbol_geometry(
    structure: &AtomicStructure,
    molecule: &Molecule,
    unit: &StructureUnit,
    ring: &[usize],
) -> CarbohydrateSymbolGeometry {
    let positions = ring_positions(structure, molecule, unit, ring);
    let center = average_position(&positions);
    let normal = ring_normal(&positions);
    let anomeric_position = anomeric_carbon_unit_index(molecule, unit, ring)
        .and_then(|unit_index| unit_position(structure, molecule, unit, unit_index))
        .unwrap_or_else(|| positions.first().copied().unwrap_or(center));
    let direction = orthogonalize(normal, (center - anomeric_position).normalized());

    CarbohydrateSymbolGeometry {
        center,
        normal,
        direction,
    }
}

fn ring_positions(
    structure: &AtomicStructure,
    molecule: &Molecule,
    unit: &StructureUnit,
    ring: &[usize],
) -> Vec<Vec3> {
    ring.iter()
        .filter_map(|&unit_index| unit_position(structure, molecule, unit, unit_index))
        .collect()
}

fn unit_position(
    structure: &AtomicStructure,
    molecule: &Molecule,
    unit: &StructureUnit,
    unit_index: usize,
) -> Option<Vec3> {
    structure.position(unit.id, unit_index).or_else(|| {
        unit.atom_indices
            .get(unit_index)
            .and_then(|&source_atom| molecule.atoms.get(source_atom))
            .map(|atom| atom.position)
    })
}

fn average_position(positions: &[Vec3]) -> Vec3 {
    if positions.is_empty() {
        return Vec3::default();
    }
    positions
        .iter()
        .copied()
        .fold(Vec3::default(), |sum, position| sum + position)
        / positions.len() as f32
}

fn ring_normal(positions: &[Vec3]) -> Vec3 {
    let mut normal = Vec3::default();
    if positions.len() >= 3 {
        for i in 0..positions.len() {
            let a = positions[i];
            let b = positions[(i + 1) % positions.len()];
            normal.x += (a.y - b.y) * (a.z + b.z);
            normal.y += (a.z - b.z) * (a.x + b.x);
            normal.z += (a.x - b.x) * (a.y + b.y);
        }
    }
    if normal.length() <= 0.000_001 && positions.len() >= 3 {
        let center = average_position(positions);
        for pair in positions.windows(2) {
            normal = normal + (pair[0] - center).cross(pair[1] - center);
        }
    }
    if normal.length() <= 0.000_001 {
        Vec3::new(0.0, 0.0, 1.0)
    } else {
        normal.normalized()
    }
}

fn orthogonalize(normal: Vec3, direction: Vec3) -> Vec3 {
    let mut out = normal.cross(direction).cross(normal).normalized();
    if out.length() > 0.000_001 {
        return out;
    }
    out = normal
        .cross(Vec3::new(1.0, 0.0, 0.0))
        .cross(normal)
        .normalized();
    if out.length() > 0.000_001 {
        return out;
    }
    out = normal
        .cross(Vec3::new(0.0, 1.0, 0.0))
        .cross(normal)
        .normalized();
    if out.length() > 0.000_001 {
        return out;
    }
    out = direction.normalized();
    if out.length() > 0.000_001 {
        out
    } else {
        Vec3::new(1.0, 0.0, 0.0)
    }
}

fn anomeric_carbon_unit_index(
    molecule: &Molecule,
    unit: &StructureUnit,
    ring: &[usize],
) -> Option<usize> {
    let mut index_has_oxygen_and_carbon = None;
    let mut index_has_c1_name = None;
    let mut index_is_carbon = None;

    for &unit_index in ring {
        let Some(&source_atom) = unit.atom_indices.get(unit_index) else {
            continue;
        };
        let Some(atom) = molecule.atoms.get(source_atom) else {
            continue;
        };
        if !atom_is_element(atom, "C") {
            continue;
        }

        let (linked_oxygen_count, linked_carbon_count) =
            linked_element_counts(molecule, source_atom);
        if linked_oxygen_count == 2 {
            return Some(unit_index);
        } else if linked_oxygen_count == 1 && linked_carbon_count == 1 {
            index_has_oxygen_and_carbon = Some(unit_index);
        } else if atom.name.starts_with("C1") {
            index_has_c1_name = Some(unit_index);
        } else if index_is_carbon.is_none() {
            index_is_carbon = Some(unit_index);
        }
    }

    index_has_oxygen_and_carbon
        .or(index_has_c1_name)
        .or(index_is_carbon)
        .or_else(|| ring.first().copied())
}

fn linked_element_counts(molecule: &Molecule, source_atom: usize) -> (usize, usize) {
    let mut oxygen = 0;
    let mut carbon = 0;
    for bond in &molecule.bonds {
        let neighbor = if bond.a == source_atom {
            Some(bond.b)
        } else if bond.b == source_atom {
            Some(bond.a)
        } else {
            None
        };
        let Some(neighbor) = neighbor else {
            continue;
        };
        let Some(atom) = molecule.atoms.get(neighbor) else {
            continue;
        };
        if atom_is_element(atom, "O") {
            oxygen += 1;
        } else if atom_is_element(atom, "C") {
            carbon += 1;
        }
    }
    (oxygen, carbon)
}

fn atom_is_element(atom: &super::Atom, element: &str) -> bool {
    atom.type_symbol.eq_ignore_ascii_case(element)
        || atom.element.eq_ignore_ascii_case(element)
        || atom.name.starts_with(element)
}

fn atom_name_for_unit_index<'a>(
    molecule: &'a Molecule,
    unit: &StructureUnit,
    unit_index: usize,
) -> &'a str {
    unit.atom_indices
        .get(unit_index)
        .and_then(|&source_atom| molecule.atoms.get(source_atom))
        .map(|atom| atom.name.as_str())
        .unwrap_or_default()
}

fn fix_link_direction(elements: &mut [CarbohydrateElement], index_a: usize, index_b: usize) {
    let Some(center_b) = elements.get(index_b).map(|element| element.geometry.center) else {
        return;
    };
    let Some(element_a) = elements.get_mut(index_a) else {
        return;
    };
    element_a.geometry.direction = (center_b - element_a.geometry.center).normalized();
}

fn fix_terminal_link_direction(
    structure: &AtomicStructure,
    elements: &mut [CarbohydrateElement],
    carbohydrate_index: usize,
    unit_id: usize,
    element_index: usize,
) {
    let Some(position) = structure.position(unit_id, element_index) else {
        return;
    };
    let Some(element) = elements.get_mut(carbohydrate_index) else {
        return;
    };
    element.geometry.direction = (position - element.geometry.center).normalized();
}

fn add_ring_element(
    map: &mut BTreeMap<(usize, usize, String), Vec<usize>>,
    residue_index: usize,
    unit_id: usize,
    alt_id: &str,
    element_index: usize,
) {
    map.entry((residue_index, unit_id, alt_id.to_string()))
        .or_default()
        .push(element_index);
}

fn sugar_ring_indices_for_residue(
    molecule: &Molecule,
    unit: &StructureUnit,
    residue_index: usize,
) -> Vec<usize> {
    molecule
        .rings
        .iter()
        .enumerate()
        .filter_map(|(ring_index, ring)| {
            SUGAR_RING_FINGERPRINTS
                .contains(&ring.fingerprint.as_str())
                .then_some(())?;
            let ring_elements = ring_unit_indices(molecule, unit, ring_index)?;
            ring_elements
                .iter()
                .all(|&index| unit.residue_index_by_element.get(index) == Some(&residue_index))
                .then_some(ring_index)
        })
        .collect()
}

fn filter_fused_rings(
    molecule: &Molecule,
    unit: &StructureUnit,
    ring_indices: Vec<usize>,
) -> Vec<usize> {
    if ring_indices.len() < 2 {
        return ring_indices;
    }
    let mut fused = BTreeSet::new();
    for a in 0..ring_indices.len() {
        for b in (a + 1)..ring_indices.len() {
            let Some(ring_a) = ring_unit_indices(molecule, unit, ring_indices[a]) else {
                continue;
            };
            let Some(ring_b) = ring_unit_indices(molecule, unit, ring_indices[b]) else {
                continue;
            };
            if sorted_intersects(&ring_a, &ring_b)
                && ring_alt_id(molecule, unit, &ring_a) == ring_alt_id(molecule, unit, &ring_b)
            {
                fused.insert(ring_indices[a]);
                fused.insert(ring_indices[b]);
            }
        }
    }
    if fused.is_empty() {
        ring_indices
    } else {
        ring_indices
            .into_iter()
            .filter(|ring_index| !fused.contains(ring_index))
            .collect()
    }
}

fn sorted_intersects(a: &[usize], b: &[usize]) -> bool {
    let mut ia = 0;
    let mut ib = 0;
    while ia < a.len() && ib < b.len() {
        if a[ia] == b[ib] {
            return true;
        }
        if a[ia] < b[ib] {
            ia += 1;
        } else {
            ib += 1;
        }
    }
    false
}

fn ring_unit_indices(
    molecule: &Molecule,
    unit: &StructureUnit,
    ring_index: usize,
) -> Option<Vec<usize>> {
    let ring = molecule.rings.get(ring_index)?;
    let mut indices = Vec::with_capacity(ring.atom_indices.len());
    for &source_atom in &ring.atom_indices {
        indices.push(
            unit.atom_indices
                .iter()
                .position(|&unit_source| unit_source == source_atom)?,
        );
    }
    indices.sort_unstable();
    Some(indices)
}

fn ring_alt_id(molecule: &Molecule, unit: &StructureUnit, ring: &[usize]) -> String {
    ring.iter()
        .filter_map(|&unit_index| unit.atom_indices.get(unit_index))
        .filter_map(|&source_atom| molecule.atoms.get(source_atom))
        .find_map(|atom| (!atom.alt_id.is_empty()).then(|| atom.alt_id.clone()))
        .unwrap_or_default()
}

fn ring_element_indices_for_unit_index(
    elements_with_ring_map: &BTreeMap<(usize, usize, String), Vec<usize>>,
    molecule: &Molecule,
    unit: &StructureUnit,
    index: usize,
) -> Vec<usize> {
    let Some(&source_atom) = unit.atom_indices.get(index) else {
        return Vec::new();
    };
    let Some(atom) = molecule.atoms.get(source_atom) else {
        return Vec::new();
    };
    let residue_index = unit
        .residue_index_by_element
        .get(index)
        .copied()
        .unwrap_or(0);
    elements_with_ring_map
        .get(&(residue_index, unit.id, atom.alt_id.clone()))
        .cloned()
        .unwrap_or_default()
}

fn are_vertex_sets_connected(
    unit: &StructureUnit,
    vertices_a: &[usize],
    vertices_b: &[usize],
    max_distance: usize,
) -> bool {
    let target = vertices_b.iter().copied().collect::<BTreeSet<_>>();
    if vertices_a.iter().any(|vertex| target.contains(vertex)) {
        return true;
    }
    let mut visited = BTreeSet::new();
    let mut frontier = vertices_a.to_vec();
    for &vertex in &frontier {
        visited.insert(vertex);
    }

    for _ in 0..max_distance {
        let mut next_frontier = Vec::new();
        for vertex in frontier {
            let Some(start) = unit.props.intra_unit_bonds.offset.get(vertex).copied() else {
                continue;
            };
            let Some(end) = unit.props.intra_unit_bonds.offset.get(vertex + 1).copied() else {
                continue;
            };
            for slot in start..end {
                let neighbor = unit.props.intra_unit_bonds.b[slot];
                if target.contains(&neighbor) {
                    return true;
                }
                if visited.insert(neighbor) {
                    next_frontier.push(neighbor);
                }
            }
        }
        if next_frontier.is_empty() {
            break;
        }
        frontier = next_frontier;
    }
    false
}

fn ring_source_atoms(structure: &AtomicStructure, element: &CarbohydrateElement) -> Vec<usize> {
    let Some(unit) = structure.unit_by_id(element.unit_id) else {
        return Vec::new();
    };
    element
        .ring_element_indices
        .iter()
        .filter_map(|&index| unit.atom_indices.get(index).copied())
        .collect()
}

fn add_unique_map_index(
    map: &mut BTreeMap<(usize, usize), Vec<usize>>,
    unit_id: usize,
    source_atom: usize,
    index: usize,
) {
    let entry = map.entry((unit_id, source_atom)).or_default();
    if !entry.contains(&index) {
        entry.push(index);
    }
}
