"""Mutant proxy structure generation."""

from __future__ import annotations

from pathlib import Path

from prostr_alfo.models.schemas import Mutation, StructureAnalysis
from prostr_alfo.structure.parser import parse_structure
from prostr_alfo.utils.biology import ONE_TO_THREE, RESIDUE_ATOMS


def create_mutant_proxy_structure(
    *,
    structure_analysis: StructureAnalysis,
    mutations: list[Mutation],
    output_path: Path,
) -> Path | None:
    """Create a lightweight mutant proxy PDB by renaming residues and pruning incompatible atoms.

    The coordinates remain WT-derived. No side-chain repacking or relaxation is performed.
    """

    if not mutations:
        return None

    parsed_residues = parse_structure(structure_analysis.structure_path)
    parsed_by_index = {residue.structure_index: residue for residue in parsed_residues}
    target_residues: dict[tuple[str, int], Mutation] = {}

    for mutation in mutations:
        structure_index = structure_analysis.sequence_to_structure.get(mutation.position)
        if structure_index is None:
            continue
        parsed_residue = parsed_by_index.get(structure_index)
        if parsed_residue is None:
            continue
        target_residues[(parsed_residue.chain_id, parsed_residue.pdb_resseq)] = mutation

    if not target_residues:
        return None

    mutated_lines: list[str] = []
    for line in structure_analysis.structure_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            chain_id = (line[21] or "A").strip() or "A"
            try:
                pdb_resseq = int(line[22:26].strip())
            except ValueError:
                mutated_lines.append(line)
                continue

            mutation = target_residues.get((chain_id, pdb_resseq))
            if mutation is None:
                mutated_lines.append(line)
                continue

            atom_name = line[12:16].strip()
            allowed_atoms = RESIDUE_ATOMS.get(mutation.mutant, {"N", "CA", "C", "O"})
            if atom_name not in allowed_atoms:
                continue

            new_resname = ONE_TO_THREE[mutation.mutant]
            mutated_lines.append(f"{line[:17]}{new_resname:>3}{line[20:]}")
            continue

        mutated_lines.append(line)

    output_path.write_text("\n".join(mutated_lines) + "\n", encoding="utf-8")
    return output_path
