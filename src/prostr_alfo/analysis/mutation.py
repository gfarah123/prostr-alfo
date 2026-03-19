"""Mutation-centric structural and biochemical analysis."""

from __future__ import annotations

import math

from prostr_alfo.analysis.scoring import classify_impact, compute_mutation_impact_score
from prostr_alfo.models.schemas import Mutation, MutationAssessment, StructureAnalysis
from prostr_alfo.utils.biology import HYDROPHOBICITY_KD, RESIDUE_CHARGE, RESIDUE_VOLUME, SPECIAL_RESIDUES


def _distance(
    point_a: tuple[float, float, float] | None,
    point_b: tuple[float, float, float] | None,
) -> float | None:
    if point_a is None or point_b is None:
        return None
    return math.sqrt(sum((coord_a - coord_b) ** 2 for coord_a, coord_b in zip(point_a, point_b, strict=False)))


def assess_mutation(
    mutation: Mutation,
    structure_analysis: StructureAnalysis | None,
) -> MutationAssessment:
    """Assess one mutation against the wild-type AlphaFold structure."""

    if structure_analysis is None:
        return MutationAssessment(
            mutation=mutation,
            mapped=False,
            message="No AlphaFold structure was available, so structure-based mutation analysis was skipped.",
        )

    structure_index = structure_analysis.sequence_to_structure.get(mutation.position)
    if structure_index is None:
        return MutationAssessment(
            mutation=mutation,
            mapped=False,
            message=f"Sequence position {mutation.position} could not be mapped onto the available structure model.",
        )

    residues_by_structure = {
        residue.structure_index: residue for residue in structure_analysis.residues if residue.structure_index is not None
    }
    target = residues_by_structure.get(structure_index)
    if target is None:
        return MutationAssessment(
            mutation=mutation,
            mapped=False,
            message=f"Mapped structure residue {structure_index} is missing from the parsed coordinate set.",
        )

    target_coord = (target.x, target.y, target.z) if None not in (target.x, target.y, target.z) else None
    neighboring_residues: list[str] = []
    local_contact_count = 0
    for residue in structure_analysis.residues:
        if residue.structure_index == target.structure_index:
            continue
        residue_coord = (residue.x, residue.y, residue.z) if None not in (residue.x, residue.y, residue.z) else None
        distance = _distance(target_coord, residue_coord)
        if distance is not None and distance <= 8.0:
            local_contact_count += 1
            neighboring_residues.append(f"{residue.residue_name}{residue.sequence_position}")

    hydrophobicity_delta = HYDROPHOBICITY_KD[mutation.mutant] - HYDROPHOBICITY_KD[mutation.wild_type]
    charge_delta = RESIDUE_CHARGE.get(mutation.mutant, 0.0) - RESIDUE_CHARGE.get(mutation.wild_type, 0.0)
    volume_delta = RESIDUE_VOLUME[mutation.mutant] - RESIDUE_VOLUME[mutation.wild_type]
    burial_score = target.burial_score or 0.0
    contact_change_proxy = round((volume_delta / 60.0) * (1.0 + burial_score + (local_contact_count / 12.0)), 2)

    special_flags: list[str] = []
    for residue_code in (mutation.wild_type, mutation.mutant):
        if residue_code in SPECIAL_RESIDUES:
            description = SPECIAL_RESIDUES[residue_code]
            if description not in special_flags:
                special_flags.append(description)
    for bond in structure_analysis.disulfide_bonds:
        if mutation.position in (bond.residue_a, bond.residue_b):
            note = "This position participates in a putative disulfide bond in the AlphaFold model."
            if note not in special_flags:
                special_flags.append(note)

    impact_score = compute_mutation_impact_score(
        residue_plddt=target.plddt,
        burial_score=burial_score,
        local_contact_count=local_contact_count,
        contact_change_proxy=contact_change_proxy,
        charge_delta=charge_delta,
        hydrophobicity_delta=hydrophobicity_delta,
        special_flag_count=len(special_flags),
    )
    impact_level = classify_impact(impact_score)

    interpretation_parts = [
        f"{mutation.code} maps to a residue with pLDDT {target.plddt:.1f}." if target.plddt is not None else None,
        f"Local contact count is {local_contact_count} with burial score {burial_score:.2f}.",
        f"Hydrophobicity delta is {hydrophobicity_delta:+.2f} and charge delta is {charge_delta:+.1f}.",
        f"Predicted heuristic impact is {impact_level.lower()} ({impact_score:.2f}/100)." if impact_level else None,
    ]
    if special_flags:
        interpretation_parts.append("Special residue considerations apply.")

    return MutationAssessment(
        mutation=mutation,
        mapped=True,
        message="Structure-based mutation assessment completed.",
        residue_plddt=target.plddt,
        burial_score=burial_score,
        local_contact_count=local_contact_count,
        contact_change_proxy=contact_change_proxy,
        hydrophobicity_delta=round(hydrophobicity_delta, 2),
        charge_delta=round(charge_delta, 2),
        special_residue_flags=special_flags,
        neighboring_residues=neighboring_residues[:15],
        impact_score=impact_score,
        impact_level=impact_level,
        interpretation=" ".join(part for part in interpretation_parts if part),
    )


def assess_mutations(
    mutations: list[Mutation],
    structure_analysis: StructureAnalysis | None,
) -> list[MutationAssessment]:
    """Assess all requested mutations."""

    return [assess_mutation(mutation, structure_analysis) for mutation in mutations]
