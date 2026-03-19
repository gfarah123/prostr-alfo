"""Structure feature extraction."""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

from Bio.PDB import DSSP, PDBParser

from prostr_alfo.models.schemas import DisulfideBond, LowConfidenceRegion, ResidueFeature, StructureAnalysis
from prostr_alfo.structure.mapping import build_sequence_to_structure_map
from prostr_alfo.structure.parser import ParsedResidue, parse_structure


def _distance(point_a: tuple[float, float, float], point_b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((coord_a - coord_b) ** 2 for coord_a, coord_b in zip(point_a, point_b, strict=False)))


def compute_contact_density(parsed_residues: list[ParsedResidue], cutoff: float = 10.0) -> dict[int, int]:
    """Compute a CA-based local density proxy."""

    contacts: dict[int, int] = {}
    for residue in parsed_residues:
        if residue.ca_coord is None:
            contacts[residue.structure_index] = 0
            continue
        count = 0
        for other in parsed_residues:
            if other.structure_index == residue.structure_index or other.ca_coord is None:
                continue
            if _distance(residue.ca_coord, other.ca_coord) <= cutoff:
                count += 1
        contacts[residue.structure_index] = count
    return contacts


def compute_disulfide_bonds(
    parsed_residues: list[ParsedResidue],
    sequence_to_structure: dict[int, int | None],
    max_distance: float = 2.5,
) -> list[DisulfideBond]:
    """Identify plausible disulfide bonds from SG atom distances."""

    inverse_map = {structure_index: sequence_index for sequence_index, structure_index in sequence_to_structure.items() if structure_index}
    cysteines = [residue for residue in parsed_residues if residue.residue_name == "C" and residue.sg_coord is not None]
    bonds: list[DisulfideBond] = []
    for index, residue in enumerate(cysteines):
        for partner in cysteines[index + 1 :]:
            distance = _distance(residue.sg_coord, partner.sg_coord)
            if distance <= max_distance:
                bonds.append(
                    DisulfideBond(
                        residue_a=inverse_map.get(residue.structure_index, residue.structure_index),
                        residue_b=inverse_map.get(partner.structure_index, partner.structure_index),
                        distance=round(distance, 3),
                    )
                )
    return bonds


def summarize_low_confidence_regions(residues: list[ResidueFeature], threshold: float = 50.0) -> list[LowConfidenceRegion]:
    """Build contiguous low-confidence segments."""

    regions: list[LowConfidenceRegion] = []
    current: list[ResidueFeature] = []
    for residue in sorted(residues, key=lambda item: item.sequence_position):
        if residue.plddt is not None and residue.plddt < threshold:
            if not current or residue.sequence_position == current[-1].sequence_position + 1:
                current.append(residue)
            else:
                regions.append(
                    LowConfidenceRegion(
                        start=current[0].sequence_position,
                        end=current[-1].sequence_position,
                        mean_plddt=round(sum(item.plddt or 0.0 for item in current) / len(current), 2),
                    )
                )
                current = [residue]
        elif current:
            regions.append(
                LowConfidenceRegion(
                    start=current[0].sequence_position,
                    end=current[-1].sequence_position,
                    mean_plddt=round(sum(item.plddt or 0.0 for item in current) / len(current), 2),
                )
            )
            current = []
    if current:
        regions.append(
            LowConfidenceRegion(
                start=current[0].sequence_position,
                end=current[-1].sequence_position,
                mean_plddt=round(sum(item.plddt or 0.0 for item in current) / len(current), 2),
            )
        )
    return regions


def assign_secondary_structure(structure_path: Path, parsed_residues: list[ParsedResidue]) -> tuple[dict[tuple[str, int], str], list[str]]:
    """Assign DSSP secondary structure when the external executable is available."""

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("alphafold", structure_path)
    model = next(structure.get_models())
    assignments: dict[tuple[str, int], str] = {}
    notes: list[str] = []

    try:
        dssp = DSSP(model, str(structure_path))
    except Exception:
        notes.append("DSSP was not available, so secondary structure annotations were skipped.")
        return assignments, notes

    for residue in parsed_residues:
        key = (residue.chain_id, (" ", residue.pdb_resseq, " "))
        if key in dssp:
            annotation = dssp[key][2]
            assignments[(residue.chain_id, residue.pdb_resseq)] = annotation if annotation.strip() else "C"
    return assignments, notes


def analyze_structure(structure_path: Path, sequence: str) -> StructureAnalysis:
    """Extract AlphaFold-derived structure features for one sequence."""

    parsed_residues = parse_structure(structure_path)
    structure_sequence = "".join(residue.sequence_code for residue in parsed_residues)
    sequence_to_structure = build_sequence_to_structure_map(sequence, structure_sequence)
    inverse_map = {structure_index: sequence_index for sequence_index, structure_index in sequence_to_structure.items() if structure_index}

    contact_density = compute_contact_density(parsed_residues)
    secondary_assignments, notes = assign_secondary_structure(structure_path, parsed_residues)
    residues: list[ResidueFeature] = []
    ss_counter: Counter[str] = Counter()
    plddt_values: list[float] = []
    surface_count = 0

    for residue in parsed_residues:
        density = contact_density[residue.structure_index]
        burial_score = min(1.0, density / 12.0)
        is_surface = density < 8
        if is_surface:
            surface_count += 1

        secondary = secondary_assignments.get((residue.chain_id, residue.pdb_resseq))
        if secondary:
            ss_counter[secondary] += 1
        if residue.plddt is not None:
            plddt_values.append(residue.plddt)

        mapped_position = inverse_map.get(residue.structure_index, residue.structure_index)
        x, y, z = residue.ca_coord if residue.ca_coord is not None else (None, None, None)
        residues.append(
            ResidueFeature(
                sequence_position=mapped_position,
                structure_index=residue.structure_index,
                chain_id=residue.chain_id,
                residue_name=residue.residue_name,
                plddt=round(residue.plddt, 2) if residue.plddt is not None else None,
                secondary_structure=secondary,
                contact_density=density,
                burial_score=round(burial_score, 3),
                is_surface=is_surface,
                x=x,
                y=y,
                z=z,
            )
        )

    mean_plddt = round(sum(plddt_values) / len(plddt_values), 2) if plddt_values else None
    percent_above_70 = round((sum(value > 70 for value in plddt_values) / len(plddt_values)) * 100, 2) if plddt_values else None
    percent_above_90 = round((sum(value > 90 for value in plddt_values) / len(plddt_values)) * 100, 2) if plddt_values else None

    low_confidence_regions = summarize_low_confidence_regions(residues)
    disulfide_bonds = compute_disulfide_bonds(parsed_residues, sequence_to_structure)

    return StructureAnalysis(
        structure_path=structure_path,
        chain_ids=sorted({residue.chain_id for residue in parsed_residues}),
        structure_sequence=structure_sequence,
        residues=residues,
        mean_plddt=mean_plddt,
        percent_above_70=percent_above_70,
        percent_above_90=percent_above_90,
        low_confidence_regions=low_confidence_regions,
        disulfide_bonds=disulfide_bonds,
        secondary_structure_counts=dict(ss_counter),
        surface_fraction=round((surface_count / len(residues)) * 100, 2) if residues else None,
        notes=notes,
        sequence_to_structure=sequence_to_structure,
    )
