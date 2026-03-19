"""Project configuration helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class Settings:
    """Filesystem locations used by the project."""

    project_root: Path = Path(__file__).resolve().parents[2]

    @property
    def data_dir(self) -> Path:
        return self.project_root / "data"

    @property
    def cache_dir(self) -> Path:
        return self.data_dir / "cache"

    @property
    def structures_dir(self) -> Path:
        return self.data_dir / "structures"

    @property
    def variants_dir(self) -> Path:
        return self.data_dir / "variants"

    @property
    def reports_dir(self) -> Path:
        return self.data_dir / "reports"

    @property
    def references_path(self) -> Path:
        return self.project_root / "src" / "prostr_alfo" / "reporting" / "references.json"

    def ensure_directories(self) -> None:
        for path in (
            self.data_dir,
            self.cache_dir,
            self.structures_dir,
            self.variants_dir,
            self.reports_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)


def get_settings() -> Settings:
    """Return default settings."""

    settings = Settings()
    settings.ensure_directories()
    return settings
