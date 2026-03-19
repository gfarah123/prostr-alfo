"""PyMOL script generation."""

from __future__ import annotations

from pathlib import Path

from prostr_alfo.models.schemas import Mutation, StructureAnalysis


def write_pymol_script(
    *,
    structure_analysis: StructureAnalysis,
    mutations: list[Mutation],
    output_path: Path,
) -> Path:
    """Create a PyMOL script with pLDDT coloring and mutation highlighting."""

    structure_path = structure_analysis.structure_path.resolve()
    mutation_positions = "+".join(str(mutation.position) for mutation in mutations)
    low_conf_positions = "+".join(
        str(residue.sequence_position)
        for residue in structure_analysis.residues
        if residue.plddt is not None and residue.plddt < 50
    )
    disulfide_positions = "+".join(
        str(position)
        for bond in structure_analysis.disulfide_bonds
        for position in (bond.residue_a, bond.residue_b)
    )

    commands = [
        f'load "{structure_path.as_posix()}", protein',
        "hide everything, protein",
        "show cartoon, protein",
        "spectrum b, red_yellow_green_cyan_blue, protein, minimum=0, maximum=100",
        "set cartoon_transparency, 0.1",
    ]
    if mutation_positions:
        commands.extend(
            [
                f"select mutation_sites, resi {mutation_positions}",
                "show sticks, mutation_sites",
                "color magenta, mutation_sites",
            ]
        )
    if low_conf_positions:
        commands.extend(
            [
                f"select low_confidence, resi {low_conf_positions}",
                "show sticks, low_confidence",
                "color orange, low_confidence",
            ]
        )
    if disulfide_positions:
        commands.extend(
            [
                f"select disulfides, resi {disulfide_positions}",
                "show sticks, disulfides",
                "color yellow, disulfides",
            ]
        )
    commands.extend(["bg_color white", "zoom protein"])
    script = "\n".join(commands)
    output_path.write_text(script + "\n", encoding="utf-8")
    return output_path
