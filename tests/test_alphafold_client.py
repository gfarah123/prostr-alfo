from pathlib import Path

from prostr_alfo.config import Settings
from prostr_alfo.structure.alphafold import AlphaFoldClient


class DummyResponse:
    def __init__(self, *, status_code: int = 200, json_data=None, content: bytes = b"", text: str = "") -> None:
        self.status_code = status_code
        self._json_data = json_data
        self.content = content
        self.text = text

    def json(self):
        return self._json_data

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP error {self.status_code}")


class RecordingSession:
    def __init__(self, responses):
        self._responses = list(responses)
        self.calls = []

    def get(self, url: str, timeout: int):
        self.calls.append(url)
        if not self._responses:
            raise AssertionError("Unexpected network request")
        return self._responses.pop(0)


def test_cached_structure_avoids_network(tmp_path: Path) -> None:
    settings = Settings(project_root=tmp_path)
    settings.ensure_directories()
    cached_path = settings.structures_dir / "AF-P12345.pdb"
    cached_path.write_text("MODEL", encoding="utf-8")

    session = RecordingSession([])
    client = AlphaFoldClient(settings=settings, session=session)

    resolved = client.get_structure_path("P12345")
    assert resolved == cached_path
    assert session.calls == []


def test_download_structure_uses_metadata_then_pdb_url(tmp_path: Path) -> None:
    settings = Settings(project_root=tmp_path)
    session = RecordingSession(
        [
            DummyResponse(json_data=[{"pdbUrl": "https://example.org/AF-P12345.pdb"}]),
            DummyResponse(content=b"ATOM  TEST\n"),
        ]
    )
    client = AlphaFoldClient(settings=settings, session=session)

    resolved = client.get_structure_path("P12345")
    assert resolved.exists()
    assert resolved.read_text(encoding="utf-8") == "ATOM  TEST\n"
    assert session.calls == [
        "https://alphafold.ebi.ac.uk/api/prediction/P12345",
        "https://example.org/AF-P12345.pdb",
    ]
