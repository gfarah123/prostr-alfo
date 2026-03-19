"""Plot generation for reports and the frontend."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt

from prostr_alfo.models.schemas import Mutation, StructureAnalysis


def plot_plddt_profile(
    *,
    structure_analysis: StructureAnalysis,
    mutations: list[Mutation],
    output_path: Path,
) -> Path:
    """Generate a pLDDT profile plot."""

    positions = [residue.sequence_position for residue in structure_analysis.residues if residue.plddt is not None]
    plddt_values = [residue.plddt for residue in structure_analysis.residues if residue.plddt is not None]

    plt.figure(figsize=(11, 4.5))
    plt.plot(positions, plddt_values, color="#1d4ed8", linewidth=1.6, label="pLDDT")
    plt.axhline(50, color="#ef4444", linestyle="--", linewidth=1, label="Low confidence")
    plt.axhline(70, color="#f59e0b", linestyle="--", linewidth=1, label="Confident")
    plt.axhline(90, color="#16a34a", linestyle="--", linewidth=1, label="Very high confidence")

    mutation_positions = [mutation.position for mutation in mutations]
    residue_lookup = {residue.sequence_position: residue.plddt for residue in structure_analysis.residues}
    valid_points = [(position, residue_lookup.get(position)) for position in mutation_positions if residue_lookup.get(position) is not None]
    if valid_points:
        plt.scatter(
            [position for position, _ in valid_points],
            [value for _, value in valid_points],
            color="#7c3aed",
            s=60,
            zorder=3,
            label="Mutation sites",
        )

    plt.xlabel("Residue position")
    plt.ylabel("pLDDT")
    plt.title("AlphaFold pLDDT Profile")
    plt.ylim(0, 100)
    plt.xlim(min(positions, default=1), max(positions, default=1))
    plt.legend(loc="lower right")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    plt.close()
    return output_path
