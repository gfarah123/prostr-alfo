"""Sequence to structure mapping helpers."""

from __future__ import annotations

from Bio.Align import PairwiseAligner


def build_sequence_to_structure_map(sequence: str, structure_sequence: str) -> dict[int, int | None]:
    """Map 1-based sequence positions to 1-based structure residue indices."""

    mapping = {index: None for index in range(1, len(sequence) + 1)}
    if not sequence or not structure_sequence:
        return mapping

    if sequence == structure_sequence:
        return {index: index for index in range(1, len(sequence) + 1)}

    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(sequence, structure_sequence)[0]

    for sequence_block, structure_block in zip(alignment.aligned[0], alignment.aligned[1], strict=False):
        seq_start, seq_end = sequence_block
        struct_start, struct_end = structure_block
        block_length = min(seq_end - seq_start, struct_end - struct_start)
        for offset in range(block_length):
            sequence_index = seq_start + offset
            structure_index = struct_start + offset
            if sequence[sequence_index] == structure_sequence[structure_index]:
                mapping[sequence_index + 1] = structure_index + 1
    return mapping
