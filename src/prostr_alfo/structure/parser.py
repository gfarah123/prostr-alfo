"""PDB parsing helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio.PDB import PDBParser

from prostr_alfo.utils.biology import THREE_TO_ONE


@dataclass(slots=True)
class ParsedResidue:
    structure_index: int
    chain_id: str
    pdb_resseq: int
    residue_name: str
    sequence_code: str
    plddt: float | None
    ca_coord: tuple[float, float, float] | None
    sg_coord: tuple[float, float, float] | None


def parse_structure(path: Path) -> list[ParsedResidue]:
    """Parse a PDB file into residue records."""

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("alphafold", path)
    model = next(structure.get_models())

    residues: list[ParsedResidue] = []
    structure_index = 0
    for chain in model:
        for residue in chain:
            hetflag, resseq, _icode = residue.id
            if hetflag.strip():
                continue

            residue_name = residue.resname.strip().upper()
            sequence_code = THREE_TO_ONE.get(residue_name)
            if sequence_code is None:
                continue

            structure_index += 1
            ca_atom = residue["CA"] if "CA" in residue else None
            sg_atom = residue["SG"] if "SG" in residue else None
            b_factors = [atom.bfactor for atom in residue.get_atoms()]
            residues.append(
                ParsedResidue(
                    structure_index=structure_index,
                    chain_id=chain.id.strip() or "A",
                    pdb_resseq=int(resseq),
                    residue_name=sequence_code,
                    sequence_code=sequence_code,
                    plddt=(sum(b_factors) / len(b_factors)) if b_factors else None,
                    ca_coord=tuple(float(value) for value in ca_atom.coord) if ca_atom is not None else None,
                    sg_coord=tuple(float(value) for value in sg_atom.coord) if sg_atom is not None else None,
                )
            )
    if not residues:
        raise ValueError(f"No protein residues could be parsed from {path}.")
    return residues
