use std::collections::BTreeMap;
use std::sync::atomic::{AtomicU64, Ordering};

use super::{EntityIndexMap, Molecule, Vec3};

#[derive(Clone, Debug, PartialEq)]
pub struct CoarseSphere {
    pub id: usize,
    pub model_num: i32,
    pub entity_id: String,
    pub asym_id: String,
    pub seq_id_begin: i32,
    pub seq_id_end: i32,
    pub position: Vec3,
    pub radius: f32,
    pub rmsf: f32,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CoarseGaussian {
    pub id: usize,
    pub model_num: i32,
    pub entity_id: String,
    pub asym_id: String,
    pub seq_id_begin: i32,
    pub seq_id_end: i32,
    pub position: Vec3,
    pub weight: f32,
    pub covariance: [[f32; 3]; 3],
}

#[derive(Clone, Debug, Default)]
pub struct CoarseModel {
    pub hierarchy: CoarseHierarchy,
    pub conformation: CoarseConformation,
}

impl CoarseModel {
    pub(super) fn from_molecule(molecule: &Molecule) -> Self {
        Self::from_parts(
            &molecule.coarse_spheres,
            &molecule.coarse_gaussians,
            &molecule.entity_index,
        )
    }

    pub(super) fn from_parts(
        sphere_rows: &[CoarseSphere],
        gaussian_rows: &[CoarseGaussian],
        entity_index: &EntityIndexMap,
    ) -> Self {
        let spheres = CoarseElements::from_spheres(sphere_rows, entity_index);
        let gaussians = CoarseElements::from_gaussians(gaussian_rows, entity_index);
        let index = CoarseIndex::new(&spheres, &gaussians);
        CoarseModel {
            hierarchy: CoarseHierarchy {
                is_defined: !sphere_rows.is_empty() || !gaussian_rows.is_empty(),
                spheres,
                gaussians,
                index,
            },
            conformation: CoarseConformation {
                id: next_coarse_conformation_id(),
                spheres: CoarseSphereConformation::from_spheres(sphere_rows),
                gaussians: CoarseGaussianConformation::from_gaussians(gaussian_rows),
            },
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct CoarseHierarchy {
    pub is_defined: bool,
    pub spheres: CoarseElements,
    pub gaussians: CoarseElements,
    pub index: CoarseIndex,
}

#[derive(Clone, Debug, Default)]
pub struct CoarseElements {
    pub count: usize,
    pub elements: Vec<CoarseElement>,
    pub entity_id: Vec<String>,
    pub asym_id: Vec<String>,
    pub seq_id_begin: Vec<i32>,
    pub seq_id_end: Vec<i32>,
    pub chain_element_segments: CoarseSegmentation,
    pub chain_key: Vec<Option<usize>>,
    pub entity_key: Vec<Option<usize>>,
    pub chain_to_entity: Vec<Option<usize>>,
    pub polymer_ranges: Vec<CoarseRange>,
    pub gap_ranges: Vec<CoarseRange>,
    sequence_maps: BTreeMap<usize, Vec<CoarseSequenceMapEntry>>,
    chain_lookup: BTreeMap<(usize, String), usize>,
}

impl CoarseElements {
    fn from_spheres(spheres: &[CoarseSphere], entity_index: &EntityIndexMap) -> Self {
        let elements = spheres
            .iter()
            .enumerate()
            .map(|(source_index, sphere)| CoarseElement {
                source_index,
                id: sphere.id,
                model_num: sphere.model_num,
                entity_id: sphere.entity_id.clone(),
                asym_id: sphere.asym_id.clone(),
                seq_id_begin: sphere.seq_id_begin,
                seq_id_end: sphere.seq_id_end,
            })
            .collect::<Vec<_>>();
        CoarseElements::from_elements(elements, entity_index)
    }

    fn from_gaussians(gaussians: &[CoarseGaussian], entity_index: &EntityIndexMap) -> Self {
        let elements = gaussians
            .iter()
            .enumerate()
            .map(|(source_index, gaussian)| CoarseElement {
                source_index,
                id: gaussian.id,
                model_num: gaussian.model_num,
                entity_id: gaussian.entity_id.clone(),
                asym_id: gaussian.asym_id.clone(),
                seq_id_begin: gaussian.seq_id_begin,
                seq_id_end: gaussian.seq_id_end,
            })
            .collect::<Vec<_>>();
        CoarseElements::from_elements(elements, entity_index)
    }

    fn from_elements(elements: Vec<CoarseElement>, entity_index: &EntityIndexMap) -> Self {
        let chain_element_segments = CoarseSegmentation::from_elements(&elements);
        let keys =
            CoarseElementKeys::from_elements(&elements, &chain_element_segments, entity_index);
        let (polymer_ranges, gap_ranges) = coarse_ranges(&elements, &chain_element_segments);
        let count = elements.len();
        let entity_id = elements
            .iter()
            .map(|element| element.entity_id.clone())
            .collect();
        let asym_id = elements
            .iter()
            .map(|element| element.asym_id.clone())
            .collect();
        let seq_id_begin = elements
            .iter()
            .map(|element| element.seq_id_begin)
            .collect();
        let seq_id_end = elements.iter().map(|element| element.seq_id_end).collect();
        CoarseElements {
            count,
            elements,
            entity_id,
            asym_id,
            seq_id_begin,
            seq_id_end,
            chain_element_segments,
            chain_key: keys.chain_key,
            entity_key: keys.entity_key,
            chain_to_entity: keys.chain_to_entity,
            polymer_ranges,
            gap_ranges,
            sequence_maps: keys.sequence_maps,
            chain_lookup: keys.chain_lookup,
        }
    }

    pub fn find_chain_key(&self, entity_id: &str, asym_id: &str) -> Option<usize> {
        self.entity_key_for_id(entity_id).and_then(|entity_key| {
            self.chain_lookup
                .get(&(entity_key, asym_id.to_string()))
                .copied()
        })
    }

    pub fn find_sequence_key(&self, entity_id: &str, asym_id: &str, seq_id: i32) -> Option<usize> {
        let chain_key = self.find_chain_key(entity_id, asym_id)?;
        self.sequence_maps
            .get(&chain_key)?
            .iter()
            .find(|entry| entry.seq_id_begin <= seq_id && seq_id <= entry.seq_id_end)
            .map(|entry| entry.element_index)
    }

    pub fn get_entity_from_chain(&self, chain_index: usize) -> Option<usize> {
        self.chain_to_entity.get(chain_index).copied().flatten()
    }

    fn entity_key_for_id(&self, entity_id: &str) -> Option<usize> {
        self.elements
            .iter()
            .zip(self.entity_key.iter())
            .find_map(|(element, key)| (element.entity_id == entity_id).then_some(*key).flatten())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoarseElement {
    pub source_index: usize,
    pub id: usize,
    pub model_num: i32,
    pub entity_id: String,
    pub asym_id: String,
    pub seq_id_begin: i32,
    pub seq_id_end: i32,
}

#[derive(Clone, Debug, Default)]
pub struct CoarseConformation {
    pub id: String,
    pub spheres: CoarseSphereConformation,
    pub gaussians: CoarseGaussianConformation,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CoarseSphereConformation {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub radius: Vec<f32>,
    pub rmsf: Vec<f32>,
}

impl CoarseSphereConformation {
    fn from_spheres(spheres: &[CoarseSphere]) -> Self {
        CoarseSphereConformation {
            x: spheres.iter().map(|sphere| sphere.position.x).collect(),
            y: spheres.iter().map(|sphere| sphere.position.y).collect(),
            z: spheres.iter().map(|sphere| sphere.position.z).collect(),
            radius: spheres.iter().map(|sphere| sphere.radius).collect(),
            rmsf: spheres.iter().map(|sphere| sphere.rmsf).collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    pub fn position(&self, index: usize) -> Option<Vec3> {
        Some(Vec3::new(
            *self.x.get(index)?,
            *self.y.get(index)?,
            *self.z.get(index)?,
        ))
    }
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct CoarseGaussianConformation {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub z: Vec<f32>,
    pub weight: Vec<f32>,
    pub covariance_matrix: Vec<[[f32; 3]; 3]>,
}

impl CoarseGaussianConformation {
    fn from_gaussians(gaussians: &[CoarseGaussian]) -> Self {
        CoarseGaussianConformation {
            x: gaussians
                .iter()
                .map(|gaussian| gaussian.position.x)
                .collect(),
            y: gaussians
                .iter()
                .map(|gaussian| gaussian.position.y)
                .collect(),
            z: gaussians
                .iter()
                .map(|gaussian| gaussian.position.z)
                .collect(),
            weight: gaussians.iter().map(|gaussian| gaussian.weight).collect(),
            covariance_matrix: gaussians
                .iter()
                .map(|gaussian| gaussian.covariance)
                .collect(),
        }
    }

    pub fn len(&self) -> usize {
        self.x.len()
    }

    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }

    pub fn position(&self, index: usize) -> Option<Vec3> {
        Some(Vec3::new(
            *self.x.get(index)?,
            *self.y.get(index)?,
            *self.z.get(index)?,
        ))
    }
}

static COARSE_CONFORMATION_ID_COUNTER: AtomicU64 = AtomicU64::new(1);

fn next_coarse_conformation_id() -> String {
    let id = COARSE_CONFORMATION_ID_COUNTER.fetch_add(1, Ordering::Relaxed);
    format!("{id:022}")
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CoarseSegmentation {
    pub offsets: Vec<usize>,
    pub index: Vec<usize>,
    pub count: usize,
}

impl CoarseSegmentation {
    fn from_elements(elements: &[CoarseElement]) -> Self {
        let mut offsets = Vec::new();
        let mut index = vec![0; elements.len()];
        if elements.is_empty() {
            return CoarseSegmentation {
                offsets: vec![0],
                index,
                count: 0,
            };
        }
        offsets.push(0);
        let mut chain_count = 0usize;
        let mut current = &elements[0].asym_id;
        for (i, element) in elements.iter().enumerate() {
            if &element.asym_id != current {
                chain_count += 1;
                offsets.push(i);
                current = &element.asym_id;
            }
            index[i] = chain_count;
        }
        offsets.push(elements.len());
        CoarseSegmentation {
            offsets,
            index,
            count: chain_count + 1,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoarseRange {
    pub start_element: usize,
    /// Inclusive max element index, matching Mol* SortedRanges storage.
    pub end_element: usize,
}

#[derive(Clone, Debug, Default)]
pub struct CoarseIndex {
    sphere_mapping: CoarseElementMapping,
    gaussian_mapping: CoarseElementMapping,
}

impl CoarseIndex {
    fn new(spheres: &CoarseElements, gaussians: &CoarseElements) -> Self {
        CoarseIndex {
            sphere_mapping: CoarseElementMapping::from_elements(&spheres.elements),
            gaussian_mapping: CoarseElementMapping::from_elements(&gaussians.elements),
        }
    }

    pub fn find_sphere_element(&self, key: &CoarseElementKey) -> Option<usize> {
        self.sphere_mapping.find(key)
    }

    pub fn find_gaussian_element(&self, key: &CoarseElementKey) -> Option<usize> {
        self.gaussian_mapping.find(key)
    }

    pub fn find_element(&self, key: &CoarseElementKey) -> Option<CoarseElementReference> {
        self.find_sphere_element(key)
            .map(|index| CoarseElementReference {
                kind: CoarseElementKind::Spheres,
                index,
            })
            .or_else(|| {
                self.find_gaussian_element(key)
                    .map(|index| CoarseElementReference {
                        kind: CoarseElementKind::Gaussians,
                        index,
                    })
            })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoarseElementKey {
    pub label_entity_id: String,
    pub label_asym_id: String,
    pub label_seq_id: i32,
}

impl CoarseElementKey {
    pub fn new(entity_id: impl Into<String>, asym_id: impl Into<String>, seq_id: i32) -> Self {
        CoarseElementKey {
            label_entity_id: entity_id.into(),
            label_asym_id: asym_id.into(),
            label_seq_id: seq_id,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoarseElementKind {
    Spheres,
    Gaussians,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoarseElementReference {
    pub kind: CoarseElementKind,
    pub index: usize,
}

#[derive(Clone, Debug, Default)]
struct CoarseElementKeys {
    chain_key: Vec<Option<usize>>,
    entity_key: Vec<Option<usize>>,
    chain_to_entity: Vec<Option<usize>>,
    sequence_maps: BTreeMap<usize, Vec<CoarseSequenceMapEntry>>,
    chain_lookup: BTreeMap<(usize, String), usize>,
}

impl CoarseElementKeys {
    fn from_elements(
        elements: &[CoarseElement],
        segmentation: &CoarseSegmentation,
        entity_index: &EntityIndexMap,
    ) -> Self {
        let entity_key = elements
            .iter()
            .map(|element| entity_index.get_entity_index(&element.entity_id))
            .collect::<Vec<_>>();
        let mut chain_key = vec![None; elements.len()];
        let mut chain_to_entity = vec![None; segmentation.count];
        let mut sequence_maps = BTreeMap::new();
        let mut chain_lookup = BTreeMap::new();
        let mut chain_counter = 0usize;

        for (chain_index, chain_entity) in chain_to_entity.iter_mut().enumerate() {
            let start = segmentation.offsets[chain_index];
            let end = segmentation.offsets[chain_index + 1];
            if start >= end {
                continue;
            }
            let Some(entity_key_for_chain) = entity_key[start] else {
                continue;
            };
            *chain_entity = Some(entity_key_for_chain);
            let lookup_key = (entity_key_for_chain, elements[start].asym_id.clone());
            let element_chain_key = if let Some(key) = chain_lookup.get(&lookup_key) {
                *key
            } else {
                let key = chain_counter;
                chain_counter += 1;
                chain_lookup.insert(lookup_key, key);
                key
            };
            for key in chain_key.iter_mut().take(end).skip(start) {
                *key = Some(element_chain_key);
            }
            let sequence_map = elements[start..end]
                .iter()
                .enumerate()
                .map(|(offset, element)| CoarseSequenceMapEntry {
                    element_index: start + offset,
                    seq_id_begin: element.seq_id_begin,
                    seq_id_end: element.seq_id_end,
                })
                .collect();
            sequence_maps.insert(element_chain_key, sequence_map);
        }

        CoarseElementKeys {
            chain_key,
            entity_key,
            chain_to_entity,
            sequence_maps,
            chain_lookup,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct CoarseSequenceMapEntry {
    element_index: usize,
    seq_id_begin: i32,
    seq_id_end: i32,
}

#[derive(Clone, Debug, Default)]
struct CoarseElementMapping {
    ranges: Vec<CoarseElementMappingRange>,
}

impl CoarseElementMapping {
    fn from_elements(elements: &[CoarseElement]) -> Self {
        CoarseElementMapping {
            ranges: elements
                .iter()
                .enumerate()
                .map(|(element_index, element)| CoarseElementMappingRange {
                    entity_id: element.entity_id.clone(),
                    asym_id: element.asym_id.clone(),
                    seq_id_begin: element.seq_id_begin,
                    seq_id_end: element.seq_id_end,
                    element_index,
                })
                .collect(),
        }
    }

    fn find(&self, key: &CoarseElementKey) -> Option<usize> {
        self.ranges
            .iter()
            .rev()
            .find(|range| {
                range.entity_id == key.label_entity_id
                    && range.asym_id == key.label_asym_id
                    && range.seq_id_begin <= key.label_seq_id
                    && key.label_seq_id <= range.seq_id_end
            })
            .map(|range| range.element_index)
    }
}

#[derive(Clone, Debug)]
struct CoarseElementMappingRange {
    entity_id: String,
    asym_id: String,
    seq_id_begin: i32,
    seq_id_end: i32,
    element_index: usize,
}

fn coarse_ranges(
    elements: &[CoarseElement],
    segmentation: &CoarseSegmentation,
) -> (Vec<CoarseRange>, Vec<CoarseRange>) {
    let mut polymer_ranges = Vec::new();
    let mut gap_ranges = Vec::new();
    for chain_index in 0..segmentation.count {
        let Some(&start) = segmentation.offsets.get(chain_index) else {
            continue;
        };
        let Some(&end) = segmentation.offsets.get(chain_index + 1) else {
            continue;
        };
        if start >= end {
            continue;
        }
        let mut range_start = start;
        let mut prev_seq_end = elements[start].seq_id_end;
        for (i, element) in elements.iter().enumerate().take(end).skip(start + 1) {
            if element.seq_id_begin - prev_seq_end > 1 {
                polymer_ranges.push(CoarseRange {
                    start_element: range_start,
                    end_element: i - 1,
                });
                gap_ranges.push(CoarseRange {
                    start_element: i - 1,
                    end_element: i,
                });
                range_start = i;
            }
            prev_seq_end = element.seq_id_end;
        }
        polymer_ranges.push(CoarseRange {
            start_element: range_start,
            end_element: end - 1,
        });
    }
    (polymer_ranges, gap_ranges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Entity;

    fn entity_index() -> EntityIndexMap {
        EntityIndexMap::from_entities(
            &[
                Entity {
                    id: "1".to_string(),
                    type_name: "polymer".to_string(),
                    description: String::new(),
                },
                Entity {
                    id: "2".to_string(),
                    type_name: "polymer".to_string(),
                    description: String::new(),
                },
            ],
            &[],
            &[],
        )
    }

    fn sphere(id: usize, entity_id: &str, asym_id: &str, begin: i32, end: i32) -> CoarseSphere {
        CoarseSphere {
            id,
            model_num: 1,
            entity_id: entity_id.to_string(),
            asym_id: asym_id.to_string(),
            seq_id_begin: begin,
            seq_id_end: end,
            position: Vec3::default(),
            radius: 1.0,
            rmsf: 0.0,
        }
    }

    fn gaussian(id: usize, entity_id: &str, asym_id: &str, begin: i32, end: i32) -> CoarseGaussian {
        CoarseGaussian {
            id,
            model_num: 1,
            entity_id: entity_id.to_string(),
            asym_id: asym_id.to_string(),
            seq_id_begin: begin,
            seq_id_end: end,
            position: Vec3::default(),
            weight: 1.0,
            covariance: [[0.0; 3]; 3],
        }
    }

    #[test]
    fn coarse_elements_store_molstar_style_fields_keys_and_segments() {
        let entity_index = entity_index();
        let spheres = vec![
            sphere(1, "1", "A", 1, 10),
            sphere(2, "1", "A", 11, 20),
            sphere(3, "1", "B", 1, 5),
            sphere(4, "2", "A", 1, 5),
        ];
        let coarse = CoarseModel::from_parts(&spheres, &[], &entity_index);
        let elements = &coarse.hierarchy.spheres;

        assert_eq!(elements.count, 4);
        assert_eq!(elements.entity_id, vec!["1", "1", "1", "2"]);
        assert_eq!(elements.asym_id, vec!["A", "A", "B", "A"]);
        assert_eq!(elements.seq_id_begin, vec![1, 11, 1, 1]);
        assert_eq!(elements.seq_id_end, vec![10, 20, 5, 5]);
        assert_eq!(elements.chain_element_segments.offsets, vec![0, 2, 3, 4]);
        assert_eq!(elements.chain_element_segments.index, vec![0, 0, 1, 2]);
        assert_eq!(
            elements.entity_key,
            vec![Some(0), Some(0), Some(0), Some(1)]
        );
        assert_eq!(elements.chain_key, vec![Some(0), Some(0), Some(1), Some(2)]);
        assert_eq!(elements.find_chain_key("1", "A"), Some(0));
        assert_eq!(elements.find_chain_key("1", "B"), Some(1));
        assert_eq!(elements.find_chain_key("2", "A"), Some(2));
        assert_eq!(elements.get_entity_from_chain(2), Some(1));
        assert_eq!(elements.find_sequence_key("1", "A", 11), Some(1));
        assert_eq!(elements.find_sequence_key("1", "A", 21), None);
    }

    #[test]
    fn coarse_ranges_match_molstar_inclusive_sorted_range_pairs() {
        let entity_index = entity_index();
        let spheres = vec![
            sphere(1, "1", "A", 1, 10),
            sphere(2, "1", "A", 11, 20),
            sphere(3, "1", "A", 25, 30),
            sphere(4, "1", "B", 1, 5),
        ];
        let gaussians = vec![gaussian(1, "1", "A", 1, 10)];
        let coarse = CoarseModel::from_parts(&spheres, &gaussians, &entity_index);
        let elements = &coarse.hierarchy.spheres;

        assert_eq!(
            elements.polymer_ranges,
            vec![
                CoarseRange {
                    start_element: 0,
                    end_element: 1,
                },
                CoarseRange {
                    start_element: 2,
                    end_element: 2,
                },
                CoarseRange {
                    start_element: 3,
                    end_element: 3,
                },
            ]
        );
        assert_eq!(
            elements.gap_ranges,
            vec![CoarseRange {
                start_element: 1,
                end_element: 2,
            }]
        );
        assert_eq!(
            coarse.hierarchy.gaussians.polymer_ranges,
            vec![CoarseRange {
                start_element: 0,
                end_element: 0,
            }]
        );
        assert!(coarse.hierarchy.gaussians.gap_ranges.is_empty());
    }

    #[test]
    fn coarse_model_empty_hierarchy_uses_molstar_empty_shape() {
        let entity_index = entity_index();
        let coarse = CoarseModel::from_parts(&[], &[], &entity_index);

        assert!(!coarse.hierarchy.is_defined);
        assert_eq!(coarse.hierarchy.spheres.count, 0);
        assert_eq!(coarse.hierarchy.gaussians.count, 0);
        assert_eq!(
            coarse.hierarchy.spheres.chain_element_segments.offsets,
            vec![0]
        );
        assert_eq!(
            coarse.hierarchy.gaussians.chain_element_segments.offsets,
            vec![0]
        );
        assert!(coarse.conformation.spheres.is_empty());
        assert!(coarse.conformation.gaussians.is_empty());
        assert_eq!(
            coarse
                .hierarchy
                .index
                .find_element(&CoarseElementKey::new("1", "A", 1)),
            None
        );
    }

    #[test]
    fn coarse_conformation_stores_molstar_columnar_arrays_and_ids() {
        let entity_index = entity_index();
        let spheres = vec![
            CoarseSphere {
                position: Vec3::new(1.0, 2.0, 3.0),
                radius: 4.0,
                rmsf: 5.0,
                ..sphere(1, "1", "A", 1, 10)
            },
            CoarseSphere {
                position: Vec3::new(6.0, 7.0, 8.0),
                radius: 9.0,
                rmsf: 10.0,
                ..sphere(2, "1", "A", 11, 20)
            },
        ];
        let gaussians = vec![CoarseGaussian {
            position: Vec3::new(11.0, 12.0, 13.0),
            weight: 14.0,
            covariance: [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            ..gaussian(1, "1", "A", 21, 30)
        }];

        let first = CoarseModel::from_parts(&spheres, &gaussians, &entity_index);
        let second = CoarseModel::from_parts(&spheres, &gaussians, &entity_index);

        assert_eq!(first.conformation.id.len(), 22);
        assert_eq!(second.conformation.id.len(), 22);
        assert_ne!(first.conformation.id, second.conformation.id);
        assert_eq!(first.conformation.spheres.x, vec![1.0, 6.0]);
        assert_eq!(first.conformation.spheres.y, vec![2.0, 7.0]);
        assert_eq!(first.conformation.spheres.z, vec![3.0, 8.0]);
        assert_eq!(first.conformation.spheres.radius, vec![4.0, 9.0]);
        assert_eq!(first.conformation.spheres.rmsf, vec![5.0, 10.0]);
        assert_eq!(
            first.conformation.spheres.position(1),
            Some(Vec3::new(6.0, 7.0, 8.0))
        );
        assert_eq!(first.conformation.gaussians.x, vec![11.0]);
        assert_eq!(first.conformation.gaussians.y, vec![12.0]);
        assert_eq!(first.conformation.gaussians.z, vec![13.0]);
        assert_eq!(first.conformation.gaussians.weight, vec![14.0]);
        assert_eq!(
            first.conformation.gaussians.covariance_matrix,
            vec![[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]]
        );
        assert_eq!(
            first.conformation.gaussians.position(0),
            Some(Vec3::new(11.0, 12.0, 13.0))
        );
    }

    #[test]
    fn coarse_index_finds_spheres_before_gaussians_and_uses_last_overlapping_range() {
        let entity_index = entity_index();
        let spheres = vec![sphere(1, "1", "A", 1, 10), sphere(2, "1", "A", 5, 8)];
        let gaussians = vec![gaussian(1, "1", "A", 5, 8)];
        let coarse = CoarseModel::from_parts(&spheres, &gaussians, &entity_index);
        let key = CoarseElementKey::new("1", "A", 6);

        assert_eq!(
            coarse.hierarchy.spheres.find_sequence_key("1", "A", 6),
            Some(0)
        );
        assert_eq!(coarse.hierarchy.index.find_sphere_element(&key), Some(1));
        assert_eq!(coarse.hierarchy.index.find_gaussian_element(&key), Some(0));
        assert_eq!(
            coarse.hierarchy.index.find_element(&key),
            Some(CoarseElementReference {
                kind: CoarseElementKind::Spheres,
                index: 1,
            })
        );
        assert_eq!(
            coarse
                .hierarchy
                .index
                .find_element(&CoarseElementKey::new("1", "A", 11)),
            None
        );
    }
}
