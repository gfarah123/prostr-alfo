"""Heuristic mutation impact scoring."""

from __future__ import annotations


def classify_impact(score: float | None) -> str | None:
    """Convert a numeric score into a qualitative impact class."""

    if score is None:
        return None
    if score >= 67:
        return "High"
    if score >= 34:
        return "Moderate"
    return "Low"


def compute_mutation_impact_score(
    *,
    residue_plddt: float | None,
    burial_score: float | None,
    local_contact_count: int | None,
    contact_change_proxy: float | None,
    charge_delta: float | None,
    hydrophobicity_delta: float | None,
    special_flag_count: int,
) -> float:
    """Compute a heuristic mutation impact score in the 0-100 range."""

    confidence_component = ((residue_plddt or 0.0) / 100.0) * 18.0
    burial_component = min(1.0, burial_score or 0.0) * 22.0
    contact_component = min(1.0, (local_contact_count or 0) / 12.0) * 18.0
    contact_change_component = min(1.0, abs(contact_change_proxy or 0.0) / 3.0) * 17.0
    charge_component = min(1.0, abs(charge_delta or 0.0) / 2.0) * 13.0
    hydrophobicity_component = min(1.0, abs(hydrophobicity_delta or 0.0) / 4.5) * 8.0
    special_component = min(1.0, special_flag_count / 2.0) * 4.0
    return round(
        confidence_component
        + burial_component
        + contact_component
        + contact_change_component
        + charge_component
        + hydrophobicity_component
        + special_component,
        2,
    )
