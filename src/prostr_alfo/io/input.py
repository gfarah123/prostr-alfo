"""Input resolution for UniProt, raw sequence, and FASTA sources."""

from __future__ import annotations

import re
from pathlib import Path

import requests

from prostr_alfo.io.fasta import normalize_sequence, read_fasta
from prostr_alfo.models.schemas import ProteinInput

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"


def infer_uniprot_id_from_header(header: str | None) -> str | None:
    """Infer a UniProt accession from a FASTA header when possible."""

    if not header:
        return None

    swissprot_match = re.search(r"(?:sp|tr)\|([A-Z0-9]{6,10})\|", header)
    if swissprot_match:
        return swissprot_match.group(1)

    token_match = re.match(r"([A-Z0-9]{6,10})\b", header)
    if token_match:
        return token_match.group(1)
    return None


def fetch_uniprot_sequence(uniprot_id: str, session: requests.Session | None = None) -> tuple[str, str]:
    """Fetch a FASTA sequence from UniProt."""

    client = session or requests.Session()
    url = UNIPROT_FASTA_URL.format(uniprot_id=uniprot_id)
    response = client.get(url, timeout=30)
    if response.status_code == 404:
        raise ValueError(f"UniProt accession {uniprot_id} was not found.")
    response.raise_for_status()

    lines = response.text.strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"UniProt accession {uniprot_id} did not return FASTA content.")
    header = lines[0][1:].strip()
    sequence = normalize_sequence("".join(lines[1:]))
    return header, sequence


def resolve_protein_input(
    *,
    uniprot_id: str | None = None,
    sequence: str | None = None,
    fasta_path: Path | None = None,
    session: requests.Session | None = None,
) -> ProteinInput:
    """Resolve one supported input source into a normalized ProteinInput."""

    provided_sources = [value is not None for value in (uniprot_id, sequence, fasta_path)]
    if sum(provided_sources) == 0:
        raise ValueError("Provide one of: uniprot_id, sequence, or fasta_path.")

    if sequence is not None:
        normalized = normalize_sequence(sequence)
        return ProteinInput(
            source_type="sequence",
            source_value=normalized,
            sequence=normalized,
            uniprot_id=uniprot_id.strip().upper() if uniprot_id else None,
        )

    if fasta_path is not None:
        header, fasta_sequence = read_fasta(Path(fasta_path))
        inferred = infer_uniprot_id_from_header(header)
        return ProteinInput(
            source_type="fasta",
            source_value=str(fasta_path),
            sequence=fasta_sequence,
            uniprot_id=(uniprot_id or inferred or "").strip().upper() or None,
            fasta_header=header,
        )

    normalized_uniprot = uniprot_id.strip().upper()
    header, uniprot_sequence = fetch_uniprot_sequence(normalized_uniprot, session=session)
    return ProteinInput(
        source_type="uniprot",
        source_value=normalized_uniprot,
        sequence=uniprot_sequence,
        uniprot_id=normalized_uniprot,
        fasta_header=header,
    )
