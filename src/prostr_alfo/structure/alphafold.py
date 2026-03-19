"""AlphaFold DB client with local caching."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests

from prostr_alfo.config import Settings, get_settings

ALPHAFOLD_METADATA_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"


class StructureNotFoundError(RuntimeError):
    """Raised when no AlphaFold structure is available."""


@dataclass(slots=True)
class AlphaFoldClient:
    settings: Settings
    session: requests.Session | None = None

    def __post_init__(self) -> None:
        self.settings.ensure_directories()
        self.session = self.session or requests.Session()

    def get_structure_path(self, uniprot_id: str, refresh: bool = False) -> Path:
        """Return a cached PDB path or download it from AlphaFold DB."""

        normalized = uniprot_id.strip().upper()
        cached_path = self.settings.structures_dir / f"AF-{normalized}.pdb"
        if cached_path.exists() and not refresh:
            return cached_path

        metadata = self.fetch_prediction_metadata(normalized)
        pdb_url = self.extract_pdb_url(metadata)
        self.download_file(pdb_url, cached_path)
        return cached_path

    def fetch_prediction_metadata(self, uniprot_id: str) -> Any:
        """Fetch prediction metadata from AlphaFold DB."""

        url = ALPHAFOLD_METADATA_URL.format(uniprot_id=uniprot_id)
        response = self.session.get(url, timeout=30)
        if response.status_code == 404:
            raise StructureNotFoundError(f"No AlphaFold prediction was found for {uniprot_id}.")
        response.raise_for_status()
        metadata = response.json()
        if not metadata:
            raise StructureNotFoundError(f"AlphaFold DB returned no metadata for {uniprot_id}.")
        return metadata

    def extract_pdb_url(self, metadata: Any) -> str:
        """Extract a PDB download URL from AlphaFold metadata."""

        records = metadata if isinstance(metadata, list) else [metadata]
        for record in records:
            direct_url = record.get("pdbUrl") or record.get("pdb_url")
            if direct_url:
                return direct_url

            files = record.get("files")
            if isinstance(files, list):
                for file_item in files:
                    if str(file_item.get("fileType", "")).lower() == "pdb" and file_item.get("url"):
                        return file_item["url"]
        raise StructureNotFoundError("AlphaFold metadata did not contain a usable PDB URL.")

    def download_file(self, url: str, destination: Path) -> None:
        """Download a file to disk."""

        response = self.session.get(url, timeout=60)
        response.raise_for_status()
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(response.content)


def get_alphafold_client(session: requests.Session | None = None) -> AlphaFoldClient:
    """Return an AlphaFold client with default settings."""

    return AlphaFoldClient(settings=get_settings(), session=session)
