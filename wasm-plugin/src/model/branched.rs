use std::collections::BTreeMap;

use super::{Molecule, PdbxBranchScheme, PdbxEntityBranchLink};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BranchedSequenceMap {
    pub entries: Vec<PdbxBranchScheme>,
    pub by_entity_id: BTreeMap<String, Vec<usize>>,
    pub by_asym_id: BTreeMap<String, Vec<usize>>,
    pub by_pdb_asym_id: BTreeMap<String, Vec<usize>>,
    pub by_auth_asym_id: BTreeMap<String, Vec<usize>>,
    pub by_entity_num: BTreeMap<(String, i32), Vec<usize>>,
    pub by_asym_num: BTreeMap<(String, i32), usize>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BranchedEntityLinkMap {
    pub links: Vec<PdbxEntityBranchLink>,
    pub placements: Vec<BranchedEntityLinkPlacement>,
    pub by_entity_id: BTreeMap<String, Vec<usize>>,
    pub by_entity_num: BTreeMap<(String, i32), Vec<usize>>,
    pub by_entity_link_id: BTreeMap<(String, i32), usize>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BranchedEntityLinkPlacement {
    pub link_index: usize,
    pub entry_index_1: usize,
    pub entry_index_2: usize,
}

impl BranchedSequenceMap {
    pub fn from_branch_scheme(entries: &[PdbxBranchScheme]) -> Self {
        let mut map = BranchedSequenceMap {
            entries: entries.to_vec(),
            ..BranchedSequenceMap::default()
        };
        for (index, entry) in map.entries.iter().enumerate() {
            push_index(&mut map.by_entity_id, &entry.entity_id, index);
            push_index(&mut map.by_asym_id, &entry.asym_id, index);
            push_index(&mut map.by_pdb_asym_id, &entry.pdb_asym_id, index);
            push_index(&mut map.by_auth_asym_id, &entry.auth_asym_id, index);
            map.by_entity_num
                .entry((entry.entity_id.clone(), entry.num))
                .or_default()
                .push(index);
            map.by_asym_num
                .entry((entry.asym_id.clone(), entry.num))
                .or_insert(index);
        }
        map
    }

    pub fn entries_for_entity(&self, entity_id: &str) -> &[usize] {
        self.by_entity_id
            .get(entity_id)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn entries_for_asym(&self, asym_id: &str) -> &[usize] {
        self.by_asym_id
            .get(asym_id)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn entries_for_entity_num(&self, entity_id: &str, num: i32) -> &[usize] {
        self.by_entity_num
            .get(&(entity_id.to_string(), num))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn entry_for_asym_num(&self, asym_id: &str, num: i32) -> Option<&PdbxBranchScheme> {
        let index = self.by_asym_num.get(&(asym_id.to_string(), num)).copied()?;
        self.entries.get(index)
    }
}

impl BranchedEntityLinkMap {
    pub fn from_links(links: &[PdbxEntityBranchLink], sequence_map: &BranchedSequenceMap) -> Self {
        let mut map = BranchedEntityLinkMap {
            links: links.to_vec(),
            ..BranchedEntityLinkMap::default()
        };
        for (index, link) in map.links.iter().enumerate() {
            push_index(&mut map.by_entity_id, &link.entity_id, index);
            map.by_entity_link_id
                .entry((link.entity_id.clone(), link.link_id))
                .or_insert(index);
            push_unique_index(
                &mut map.by_entity_num,
                (link.entity_id.clone(), link.entity_branch_list_num_1),
                index,
            );
            push_unique_index(
                &mut map.by_entity_num,
                (link.entity_id.clone(), link.entity_branch_list_num_2),
                index,
            );

            let entries_1 =
                sequence_map.entries_for_entity_num(&link.entity_id, link.entity_branch_list_num_1);
            let entries_2 =
                sequence_map.entries_for_entity_num(&link.entity_id, link.entity_branch_list_num_2);
            if entries_1.is_empty() || entries_2.is_empty() {
                continue;
            }

            let placement_start = map.placements.len();
            for &entry_index_1 in entries_1 {
                let Some(entry_1) = sequence_map.entries.get(entry_index_1) else {
                    continue;
                };
                for &entry_index_2 in entries_2 {
                    let Some(entry_2) = sequence_map.entries.get(entry_index_2) else {
                        continue;
                    };
                    if entry_1.asym_id == entry_2.asym_id {
                        map.placements.push(BranchedEntityLinkPlacement {
                            link_index: index,
                            entry_index_1,
                            entry_index_2,
                        });
                    }
                }
            }
            if map.placements.len() != placement_start {
                continue;
            }

            for &entry_index_1 in entries_1 {
                for &entry_index_2 in entries_2 {
                    map.placements.push(BranchedEntityLinkPlacement {
                        link_index: index,
                        entry_index_1,
                        entry_index_2,
                    });
                }
            }
        }
        map
    }

    pub fn links_for_entity(&self, entity_id: &str) -> &[usize] {
        self.by_entity_id
            .get(entity_id)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn links_for_entity_num(&self, entity_id: &str, num: i32) -> &[usize] {
        self.by_entity_num
            .get(&(entity_id.to_string(), num))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub fn link_for_entity_link_id(
        &self,
        entity_id: &str,
        link_id: i32,
    ) -> Option<&PdbxEntityBranchLink> {
        let index = self
            .by_entity_link_id
            .get(&(entity_id.to_string(), link_id))
            .copied()?;
        self.links.get(index)
    }
}

impl Molecule {
    pub fn branched_sequence_map(&self) -> BranchedSequenceMap {
        BranchedSequenceMap::from_branch_scheme(&self.pdbx_branch_scheme)
    }

    pub fn branched_entity_link_map(&self) -> BranchedEntityLinkMap {
        BranchedEntityLinkMap::from_links(
            &self.pdbx_entity_branch_links,
            &self.branched_sequence_map(),
        )
    }
}

fn push_index(map: &mut BTreeMap<String, Vec<usize>>, key: &str, index: usize) {
    if key.is_empty() {
        return;
    }
    map.entry(key.to_string()).or_default().push(index);
}

fn push_unique_index<K: Ord>(map: &mut BTreeMap<K, Vec<usize>>, key: K, index: usize) {
    let indices = map.entry(key).or_default();
    if !indices.contains(&index) {
        indices.push(index);
    }
}
