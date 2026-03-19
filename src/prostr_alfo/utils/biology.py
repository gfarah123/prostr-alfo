"""Biochemical constants used across the project."""

from __future__ import annotations

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
EXTENDED_AMINO_ACIDS = AMINO_ACIDS | {"X", "U", "O"}

THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "SEC": "U",
    "PYL": "O",
}

ONE_TO_THREE = {value: key for key, value in THREE_TO_ONE.items() if len(value) == 1}

RESIDUE_ATOMS = {
    "A": {"N", "CA", "C", "O", "CB"},
    "R": {"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "N": {"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"},
    "D": {"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"},
    "C": {"N", "CA", "C", "O", "CB", "SG"},
    "Q": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"},
    "E": {"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"},
    "G": {"N", "CA", "C", "O"},
    "H": {"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "I": {"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"},
    "L": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"},
    "K": {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"},
    "M": {"N", "CA", "C", "O", "CB", "CG", "SD", "CE"},
    "F": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "P": {"N", "CA", "C", "O", "CB", "CG", "CD"},
    "S": {"N", "CA", "C", "O", "CB", "OG"},
    "T": {"N", "CA", "C", "O", "CB", "OG1", "CG2"},
    "W": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "Y": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "V": {"N", "CA", "C", "O", "CB", "CG1", "CG2"},
}

HYDROPHOBICITY_KD = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

RESIDUE_CHARGE = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1,
    "H": 0.5,
}

RESIDUE_VOLUME = {
    "A": 88.6,
    "R": 173.4,
    "N": 114.1,
    "D": 111.1,
    "C": 108.5,
    "Q": 143.8,
    "E": 138.4,
    "G": 60.1,
    "H": 153.2,
    "I": 166.7,
    "L": 166.7,
    "K": 168.6,
    "M": 162.9,
    "F": 189.9,
    "P": 112.7,
    "S": 89.0,
    "T": 116.1,
    "W": 227.8,
    "Y": 193.6,
    "V": 140.0,
}

SPECIAL_RESIDUES = {
    "G": "Glycine contributes backbone flexibility.",
    "P": "Proline constrains backbone geometry.",
    "C": "Cysteine can participate in disulfide bonds or redox-sensitive chemistry.",
}
