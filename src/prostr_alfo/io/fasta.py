"""FASTA parsing and writing utilities."""

from __future__ import annotations

from io import StringIO
from pathlib import Path

from Bio import SeqIO

from prostr_alfo.utils.biology import EXTENDED_AMINO_ACIDS


def normalize_sequence(sequence: str) -> str:
    """Normalize and validate an amino-acid sequence."""

    normalized = "".join(sequence.split()).upper()
    invalid = sorted(set(normalized) - EXTENDED_AMINO_ACIDS)
    if not normalized:
        raise ValueError("Sequence input is empty.")
    if invalid:
        invalid_text = ", ".join(invalid)
        raise ValueError(f"Sequence contains unsupported residues: {invalid_text}.")
    return normalized


def read_fasta(path: Path) -> tuple[str, str]:
    """Read the first FASTA record from disk."""

    with path.open("r", encoding="utf-8") as handle:
        record = next(SeqIO.parse(handle, "fasta"), None)
    if record is None:
        raise ValueError(f"No FASTA records were found in {path}.")
    return record.description, normalize_sequence(str(record.seq))


def parse_fasta_text(text: str) -> tuple[str, str]:
    """Read the first FASTA record from raw text."""

    handle = StringIO(text)
    record = next(SeqIO.parse(handle, "fasta"), None)
    if record is None:
        raise ValueError("Uploaded FASTA content did not contain any records.")
    return record.description, normalize_sequence(str(record.seq))


def write_fasta(path: Path, header: str, sequence: str) -> None:
    """Write a sequence as FASTA using 80-character wrapping."""

    path.parent.mkdir(parents=True, exist_ok=True)
    wrapped = "\n".join(sequence[i : i + 80] for i in range(0, len(sequence), 80))
    content = f">{header}\n{wrapped}\n"
    path.write_text(content, encoding="utf-8")
