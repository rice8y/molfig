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

pub(crate) fn infer_bonds(atoms: &[Atom]) -> Vec<Bond> {
    let mut bonds = Vec::new();
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
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
