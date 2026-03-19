"""Path and naming helpers."""

from __future__ import annotations

import re


def slugify(value: str) -> str:
    """Create a filesystem-safe slug."""

    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return slug.strip("_") or "analysis"
