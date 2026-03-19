"""PyMOL script generation."""

from __future__ import annotations

from pathlib import Path

from prostr_alfo.models.schemas import Mutation, StructureAnalysis


def write_pymol_script(
    *,
    structure_analysis: StructureAnalysis,
    mutations: list[Mutation],
    mutant_structure_path: Path | None,
    output_path: Path,
) -> Path:
    """Create a PyMOL script with WT/mutant views and mutation highlighting."""

    structure_path = structure_analysis.structure_path.resolve()
    mutant_path = mutant_structure_path.resolve() if mutant_structure_path is not None else None
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
        f'load "{structure_path.as_posix()}", wt',
        "hide everything, wt",
        "show cartoon, wt",
        "spectrum b, red_yellow_green_cyan_blue, wt, minimum=0, maximum=100",
        "set cartoon_transparency, 0.1",
    ]
    if mutant_path is not None:
        commands.extend(
            [
                f'load "{mutant_path.as_posix()}", mutant',
                "hide everything, mutant",
                "show cartoon, mutant",
                "color lightblue, mutant",
                "translate [35, 0, 0], object=mutant",
            ]
        )
    if mutation_positions:
        commands.extend(
            [
                f"select wt_mutation_sites, wt and resi {mutation_positions}",
                "show sticks, wt_mutation_sites",
                "color magenta, wt_mutation_sites",
            ]
        )
        if mutant_path is not None:
            commands.extend(
                [
                    f"select mutant_mutation_sites, mutant and resi {mutation_positions}",
                    "show sticks, mutant_mutation_sites",
                    "color red, mutant_mutation_sites",
                ]
            )
    if low_conf_positions:
        commands.extend(
            [
                f"select low_confidence, wt and resi {low_conf_positions}",
                "show sticks, low_confidence",
                "color orange, low_confidence",
            ]
        )
    if disulfide_positions:
        commands.extend(
            [
                f"select disulfides, wt and resi {disulfide_positions}",
                "show sticks, disulfides",
                "color yellow, disulfides",
            ]
        )
    commands.extend(["bg_color white", "orient wt"])
    script = "\n".join(commands)
    output_path.write_text(script + "\n", encoding="utf-8")
    return output_path
