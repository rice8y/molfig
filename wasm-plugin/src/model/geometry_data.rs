use super::Vec3;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Face {
    pub(crate) a: usize,
    pub(crate) b: usize,
    pub(crate) c: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct MeshMaterial {
    pub(crate) color: u32,
    pub(crate) alpha_tenths: u8,
}

impl MeshMaterial {
    pub(crate) fn opaque(color: u32) -> Self {
        Self {
            color: color & 0x00ff_ffff,
            alpha_tenths: 10,
        }
    }

    pub(crate) fn with_alpha_tenths(color: u32, alpha_tenths: u8) -> Self {
        Self {
            color: color & 0x00ff_ffff,
            alpha_tenths: alpha_tenths.min(10),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MeshSection {
    pub(crate) key: String,
    pub(crate) vertex_start: usize,
    pub(crate) vertex_end: usize,
    pub(crate) face_start: usize,
    pub(crate) face_end: usize,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Mesh {
    pub(crate) vertices: Vec<Vec3>,
    pub(crate) normals: Vec<Vec3>,
    pub(crate) faces: Vec<Face>,
    pub(crate) vertex_groups: Vec<usize>,
    pub(crate) face_groups: Vec<usize>,
    pub(crate) face_materials: Vec<MeshMaterial>,
    pub(crate) sections: Vec<MeshSection>,
    pub(crate) group_count: usize,
}

impl Mesh {
    pub(crate) fn face_group(&self, face_index: usize) -> usize {
        self.face_groups.get(face_index).copied().unwrap_or(0)
    }

    pub(crate) fn face_material(&self, face_index: usize) -> Option<MeshMaterial> {
        self.face_materials.get(face_index).copied()
    }

    pub(crate) fn effective_group_count(&self) -> usize {
        self.face_groups
            .iter()
            .copied()
            .max()
            .map(|group| self.group_count.max(group + 1))
            .unwrap_or_else(|| {
                if self.faces.is_empty() {
                    self.group_count
                } else {
                    self.group_count.max(1)
                }
            })
    }
}

#[derive(Clone, Debug)]
pub(crate) struct TraceResidue {
    pub(crate) chain: String,
    #[allow(dead_code)]
    pub(crate) residue: String,
    pub(crate) seq: i32,
    pub(crate) insertion_code: String,
    pub(crate) position: Vec3,
    pub(crate) direction: Option<Vec3>,
    pub(crate) initial: bool,
    pub(crate) final_residue: bool,
    pub(crate) sec_struc_first: bool,
    pub(crate) sec_struc_last: bool,
    pub(crate) is_nucleotide: bool,
    pub(crate) base_center: Option<Vec3>,
    pub(crate) base_normal: Option<Vec3>,
    pub(crate) nucleotide_atoms: Option<NucleotideAtoms>,
    pub(crate) nucleotide_base_kind: Option<NucleotideBaseKind>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum NucleotideBaseKind {
    Purine,
    Pyrimidine,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct NucleotideAtoms {
    pub(crate) trace: Option<Vec3>,
    pub(crate) n1: Option<Vec3>,
    pub(crate) c1: Option<Vec3>,
    pub(crate) c2: Option<Vec3>,
    pub(crate) n3: Option<Vec3>,
    pub(crate) c4: Option<Vec3>,
    pub(crate) c5: Option<Vec3>,
    pub(crate) n5: Option<Vec3>,
    pub(crate) c6: Option<Vec3>,
    pub(crate) n7: Option<Vec3>,
    pub(crate) c7: Option<Vec3>,
    pub(crate) c8: Option<Vec3>,
    pub(crate) n9: Option<Vec3>,
    pub(crate) c1_prime: Option<Vec3>,
    pub(crate) c2_prime: Option<Vec3>,
    pub(crate) c3_prime: Option<Vec3>,
    pub(crate) c4_prime: Option<Vec3>,
    pub(crate) o4_prime: Option<Vec3>,
}

impl NucleotideAtoms {
    pub(crate) fn record_atom(&mut self, name: &str, position: Vec3) {
        match name {
            "N1" => set_first(&mut self.n1, position),
            "C1" => set_first(&mut self.c1, position),
            "C2" => set_first(&mut self.c2, position),
            "N3" => set_first(&mut self.n3, position),
            "C4" => set_first(&mut self.c4, position),
            "C5" => set_first(&mut self.c5, position),
            "N5" => set_first(&mut self.n5, position),
            "C6" => set_first(&mut self.c6, position),
            "N7" => set_first(&mut self.n7, position),
            "C7" => set_first(&mut self.c7, position),
            "C8" => set_first(&mut self.c8, position),
            "N9" => set_first(&mut self.n9, position),
            "C1'" | "C1*" => set_first(&mut self.c1_prime, position),
            "C2'" | "C2*" => set_first(&mut self.c2_prime, position),
            "C3'" | "C3*" => set_first(&mut self.c3_prime, position),
            "C4'" | "C4*" => set_first(&mut self.c4_prime, position),
            "O4'" | "O4*" => set_first(&mut self.o4_prime, position),
            _ => {}
        }
    }

    pub(crate) fn set_trace(&mut self, position: Vec3) {
        self.trace = Some(position);
    }
}

fn set_first(slot: &mut Option<Vec3>, position: Vec3) {
    if slot.is_none() {
        *slot = Some(position);
    }
}
