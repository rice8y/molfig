use crate::model::{Atom, Bond};

pub(crate) fn infer_element_from_name(name: &str) -> String {
    let letters: String = name
        .chars()
        .filter(|c| c.is_ascii_alphabetic())
        .take(2)
        .collect();
    if letters.is_empty() {
        "C".to_string()
    } else {
        letters
    }
}

pub(crate) fn normalize_element(value: String) -> String {
    let mut chars = value.trim().chars().filter(|c| c.is_ascii_alphabetic());
    let Some(first) = chars.next() else {
        return "C".to_string();
    };
    let mut out = String::new();
    out.push(first.to_ascii_uppercase());
    if let Some(second) = chars.next() {
        out.push(second.to_ascii_lowercase());
    }
    out
}

pub(crate) fn atomic_number(element: &str) -> u8 {
    match normalize_element(element.to_string()).as_str() {
        "H" | "D" | "T" => 1,
        "He" => 2,
        "Li" => 3,
        "Be" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Ne" => 10,
        "Na" => 11,
        "Mg" => 12,
        "Al" => 13,
        "Si" => 14,
        "P" => 15,
        "S" => 16,
        "Cl" => 17,
        "Ar" => 18,
        "K" => 19,
        "Ca" => 20,
        "Sc" => 21,
        "Ti" => 22,
        "V" => 23,
        "Cr" => 24,
        "Mn" => 25,
        "Fe" => 26,
        "Co" => 27,
        "Ni" => 28,
        "Cu" => 29,
        "Zn" => 30,
        "Ga" => 31,
        "Ge" => 32,
        "As" => 33,
        "Se" => 34,
        "Br" => 35,
        "Kr" => 36,
        "Rb" => 37,
        "Sr" => 38,
        "Y" => 39,
        "Zr" => 40,
        "Nb" => 41,
        "Mo" => 42,
        "Tc" => 43,
        "Ru" => 44,
        "Rh" => 45,
        "Pd" => 46,
        "Ag" => 47,
        "Cd" => 48,
        "In" => 49,
        "Sn" => 50,
        "Sb" => 51,
        "Te" => 52,
        "I" => 53,
        "Xe" => 54,
        "Cs" => 55,
        "Ba" => 56,
        "La" => 57,
        "Ce" => 58,
        "Pr" => 59,
        "Nd" => 60,
        "Pm" => 61,
        "Sm" => 62,
        "Eu" => 63,
        "Gd" => 64,
        "Tb" => 65,
        "Dy" => 66,
        "Ho" => 67,
        "Er" => 68,
        "Tm" => 69,
        "Yb" => 70,
        "Lu" => 71,
        "Hf" => 72,
        "Ta" => 73,
        "W" => 74,
        "Re" => 75,
        "Os" => 76,
        "Ir" => 77,
        "Pt" => 78,
        "Au" => 79,
        "Hg" => 80,
        "Tl" => 81,
        "Pb" => 82,
        "Bi" => 83,
        "Po" => 84,
        "At" => 85,
        "Rn" => 86,
        "Fr" => 87,
        "Ra" => 88,
        "Ac" => 89,
        "Th" => 90,
        "Pa" => 91,
        "U" => 92,
        "Np" => 93,
        "Pu" => 94,
        "Am" => 95,
        "Cm" => 96,
        "Bk" => 97,
        "Cf" => 98,
        "Es" => 99,
        "Fm" => 100,
        "Md" => 101,
        "No" => 102,
        "Lr" => 103,
        "Rf" => 104,
        "Db" => 105,
        "Sg" => 106,
        "Bh" => 107,
        "Hs" => 108,
        "Mt" => 109,
        _ => 0,
    }
}

#[allow(dead_code)]
fn covalent_radius(element: &str) -> f32 {
    match element {
        "H" => 0.31,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "P" => 1.07,
        "S" => 1.05,
        "F" => 0.57,
        "Cl" => 1.02,
        "Br" => 1.20,
        "I" => 1.39,
        "Fe" => 1.24,
        "Mg" => 1.30,
        "Zn" => 1.22,
        "Ca" => 1.74,
        _ => 0.77,
    }
}

pub(crate) fn vdw_radius(element: &str) -> f32 {
    vdw_radius64(element) as f32
}

pub(crate) fn vdw_radius64(element: &str) -> f64 {
    match atomic_number(element) {
        1 => 1.1,
        2 => 1.4,
        3 => 1.81,
        4 => 1.53,
        5 => 1.92,
        6 => 1.7,
        7 => 1.55,
        8 => 1.52,
        9 => 1.47,
        10 => 1.54,
        11 => 2.27,
        12 => 1.73,
        13 => 1.84,
        14 => 2.1,
        15 => 1.8,
        16 => 1.8,
        17 => 1.75,
        18 => 1.88,
        19 => 2.75,
        20 => 2.31,
        21 => 2.3,
        22 => 2.15,
        23 => 2.05,
        24 => 2.05,
        25 => 2.05,
        26 => 2.05,
        27 => 2.0,
        28 => 2.0,
        29 => 2.0,
        30 => 2.1,
        31 => 1.87,
        32 => 2.11,
        33 => 1.85,
        34 => 1.9,
        35 => 1.83,
        36 => 2.02,
        37 => 3.03,
        38 => 2.49,
        39 => 2.4,
        40 => 2.3,
        41 => 2.15,
        42 => 2.1,
        43 => 2.05,
        44 => 2.05,
        45 => 2.0,
        46 => 2.05,
        47 => 2.1,
        48 => 2.2,
        49 => 2.2,
        50 => 1.93,
        51 => 2.17,
        52 => 2.06,
        53 => 1.98,
        54 => 2.16,
        55 => 3.43,
        56 => 2.68,
        57 => 2.5,
        58 => 2.48,
        59 => 2.47,
        60 => 2.45,
        61 => 2.43,
        62 => 2.42,
        63 => 2.4,
        64 => 2.38,
        65 => 2.37,
        66 => 2.35,
        67 => 2.33,
        68 => 2.32,
        69 => 2.3,
        70 => 2.28,
        71 => 2.27,
        72 => 2.25,
        73 => 2.2,
        74 => 2.1,
        75 => 2.05,
        76 => 2.0,
        77 => 2.0,
        78 => 2.05,
        79 => 2.1,
        80 => 2.05,
        81 => 1.96,
        82 => 2.02,
        83 => 2.07,
        84 => 1.97,
        85 => 2.02,
        86 => 2.2,
        87 => 3.48,
        88 => 2.83,
        89..=109 => 2.0,
        _ => 1.7,
    }
}

#[allow(dead_code)]
pub(crate) fn infer_bonds(atoms: &[Atom]) -> Vec<Bond> {
    let mut bonds = Vec::new();
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            if atoms[i].model_num != atoms[j].model_num {
                continue;
            }
            let d = atoms[i].position.distance(atoms[j].position);
            let cutoff =
                covalent_radius(&atoms[i].element) + covalent_radius(&atoms[j].element) + 0.45;
            if d > 0.35 && d <= cutoff {
                bonds.push(Bond { a: i, b: j });
            }
        }
    }
    bonds
}

pub(crate) fn infer_xyz_bonds_molstar(atoms: &[Atom]) -> Vec<Bond> {
    let mut bonds = Vec::new();
    for i in 0..atoms.len() {
        let element_i = molstar_bond_element_index(&atoms[i].element);
        for j in (i + 1)..atoms.len() {
            if atoms[i].model_num != atoms[j].model_num {
                continue;
            }
            let element_j = molstar_bond_element_index(&atoms[j].element);
            if element_i == 0 && element_j == 0 {
                continue;
            }

            let delta = atoms[i].position - atoms[j].position;
            let distance = ((delta.x as f64) * (delta.x as f64)
                + (delta.y as f64) * (delta.y as f64)
                + (delta.z as f64) * (delta.z as f64))
                .sqrt();
            if distance == 0.0 {
                continue;
            }
            let threshold = molstar_pairing_threshold(element_i, element_j);
            if distance <= threshold {
                bonds.push(Bond { a: i, b: j });
            }
        }
    }
    bonds
}

fn molstar_bond_element_index(element: &str) -> i16 {
    match atomic_number(element) {
        0 => -1,
        1 => 0,
        number => number as i16,
    }
}

fn molstar_pairing_threshold(a: i16, b: i16) -> f64 {
    if let Some(threshold) = molstar_element_pair_threshold(a, b) {
        return threshold;
    }
    let threshold_a = molstar_element_bond_threshold(a);
    if b < 0 {
        threshold_a
    } else {
        (threshold_a + molstar_element_bond_threshold(b)) / 1.95
    }
}

fn molstar_element_bond_threshold(element: i16) -> f64 {
    match element {
        0 => 1.42,
        3 | 4 | 11..=13 | 19..=31 | 37..=50 | 55..=83 | 87..=108 => 2.7,
        6 => 1.75,
        7 => 1.6,
        8 => 1.52,
        14 => 1.9,
        15 => 2.0,
        16 => 1.9,
        17 => 1.8,
        33 => 2.68,
        109 => 2.88,
        _ => 2.001,
    }
}

fn molstar_element_pair_threshold(a: i16, b: i16) -> Option<f64> {
    if a < 0 || b < 0 {
        return None;
    }
    let (min, max) = if a < b { (a, b) } else { (b, a) };
    let key = (min + max) * (min + max + 1) / 2 + max;
    Some(match key {
        0 => 0.8,
        20 => 1.31,
        27 => 1.2,
        35 => 1.15,
        44 => 1.1,
        54 => 1.0,
        60 => 1.84,
        72 => 1.88,
        84 => 1.75,
        85 => 1.56,
        86 => 1.76,
        98 => 1.6,
        99 => 1.68,
        100 => 1.63,
        112 => 1.6,
        113 => 1.59,
        114 => 1.36,
        129 => 1.45,
        135 => 1.47,
        144 => 1.6,
        152 => 1.45,
        170 => 1.4,
        180 => 1.55,
        202 => 2.4,
        222 => 2.24,
        224 => 1.91,
        225 => 1.98,
        243 => 2.02,
        269 => 2.0,
        293 => 1.9,
        316 => 1.8,
        420 => 2.37,
        480 => 2.3,
        512 => 2.3,
        544 => 2.3,
        612 => 2.1,
        629 => 1.54,
        665 => 1.0,
        813 => 2.6,
        851 => 2.65,
        854 => 2.27,
        894 => 1.93,
        896 => 2.1,
        937 => 2.05,
        938 => 2.06,
        981 => 1.62,
        1258 => 2.68,
        1309 => 2.33,
        1484 => 1.0,
        1763 => 2.14,
        1823 => 2.48,
        1882 => 2.1,
        1944 => 1.72,
        2063 => 2.72,
        2380 => 2.34,
        3132 => 2.6,
        3367 => 2.44,
        3733 => 2.11,
        3819 => 2.6,
        3821 => 2.36,
        4736 => 2.75,
        5724 => 2.73,
        5959 => 2.63,
        6519 => 2.84,
        6750 => 2.87,
        8991 => 2.81,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(element: &str, x: f32, model_num: i32) -> Atom {
        Atom {
            id: 0,
            source_index: 0,
            model_num,
            name: element.to_string(),
            type_symbol: element.to_string(),
            element: element.to_string(),
            chain: "A".to_string(),
            auth_chain: "A".to_string(),
            entity_id: "1".to_string(),
            residue: "MOL".to_string(),
            auth_residue: "MOL".to_string(),
            group_pdb: String::new(),
            residue_seq: "1".to_string(),
            auth_residue_seq: "1".to_string(),
            insertion_code: String::new(),
            alt_id: String::new(),
            auth_name: element.to_string(),
            occupancy: 1.0,
            b_iso: 0.0,
            formal_charge: 0,
            position: crate::model::Vec3::new(x, 0.0, 0.0),
            position64: [x as f64, 0.0, 0.0],
            het: false,
            operator_name: String::new(),
        }
    }

    #[test]
    fn xyz_bond_inference_uses_molstar_pair_thresholds() {
        assert_eq!(
            infer_xyz_bonds_molstar(&[atom("C", 0.0, 0), atom("H", 1.19, 0)]).len(),
            1
        );
        assert!(infer_xyz_bonds_molstar(&[atom("C", 0.0, 0), atom("H", 1.21, 0)]).is_empty());
        assert!(infer_xyz_bonds_molstar(&[atom("H", 0.0, 0), atom("H", 0.741, 0)]).is_empty());
    }

    #[test]
    fn xyz_bond_inference_never_connects_models() {
        assert!(infer_xyz_bonds_molstar(&[atom("C", 0.0, 0), atom("H", 1.0, 1)]).is_empty());
    }
}
