"""Point mutation parsing and application."""

from __future__ import annotations

import json
import re
from pathlib import Path

from prostr_alfo.io.fasta import write_fasta
from prostr_alfo.models.schemas import Mutation, VariantRecord
from prostr_alfo.utils.biology import AMINO_ACIDS

MUTATION_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def parse_mutations(expression: str | None) -> list[Mutation]:
    """Parse a mutation expression like A23V;C105S."""

    if expression is None or not expression.strip():
        return []

    mutations: list[Mutation] = []
    seen_positions: set[int] = set()
    for chunk in expression.split(";"):
        token = chunk.strip().upper()
        if not token:
            continue
        match = MUTATION_PATTERN.match(token)
        if match is None:
            raise ValueError(f"Invalid mutation token: {token}. Expected format like A23V.")
        wild_type, position_text, mutant = match.groups()
        position = int(position_text)
        if wild_type not in AMINO_ACIDS or mutant not in AMINO_ACIDS:
            raise ValueError(f"Mutation {token} uses unsupported amino-acid codes.")
        if position in seen_positions:
            raise ValueError(f"Mutation position {position} was provided more than once.")
        seen_positions.add(position)
        mutations.append(Mutation(wild_type=wild_type, position=position, mutant=mutant))
    return mutations


def validate_mutations(sequence: str, mutations: list[Mutation]) -> None:
    """Validate that mutations match the wild-type sequence."""

    for mutation in mutations:
        if mutation.position < 1 or mutation.position > len(sequence):
            raise ValueError(f"Mutation {mutation.code} is outside the sequence length {len(sequence)}.")
        observed = sequence[mutation.position - 1]
        if observed != mutation.wild_type:
            raise ValueError(
                f"Mutation {mutation.code} expects wild-type residue {mutation.wild_type} "
                f"at position {mutation.position}, but the sequence contains {observed}."
            )


def apply_mutations(sequence: str, mutations: list[Mutation]) -> str:
    """Apply one or more validated mutations to a sequence."""

    validate_mutations(sequence, mutations)
    residues = list(sequence)
    for mutation in mutations:
        residues[mutation.position - 1] = mutation.mutant
    return "".join(residues)


def variant_name_for_mutations(mutations: list[Mutation]) -> str:
    """Create a deterministic variant name."""

    return "WT" if not mutations else "_".join(mutation.code for mutation in mutations)


def create_variant_registry(
    *,
    sequence: str,
    mutations: list[Mutation],
    output_dir: Path,
    prefix: str,
) -> tuple[list[VariantRecord], Path]:
    """Create FASTA files for WT and cumulative variants, plus a registry JSON."""

    output_dir.mkdir(parents=True, exist_ok=True)
    prefix_clean = prefix or "protein"
    variants: list[VariantRecord] = []

    cumulative: list[Mutation] = []
    for mutation_set in [[]] + [[mutation] for mutation in mutations]:
        if mutation_set:
            cumulative.append(mutation_set[0])
        current_mutations = cumulative.copy()
        variant_name = variant_name_for_mutations(current_mutations)
        variant_sequence = apply_mutations(sequence, current_mutations) if current_mutations else sequence
        fasta_path = output_dir / f"{prefix_clean}_{variant_name}.fasta"
        write_fasta(fasta_path, f"{prefix_clean}|{variant_name}", variant_sequence)
        variants.append(
            VariantRecord(
                name=variant_name,
                mutations=current_mutations,
                sequence=variant_sequence,
                fasta_path=fasta_path,
            )
        )

    registry_path = output_dir / f"{prefix_clean}_variant_registry.json"
    payload = [
        {
            "name": variant.name,
            "mutations": [mutation.code for mutation in variant.mutations],
            "sequence_length": len(variant.sequence),
            "fasta_path": str(variant.fasta_path),
        }
        for variant in variants
    ]
    registry_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return variants, registry_path
