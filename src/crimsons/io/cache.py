from __future__ import annotations

from pathlib import Path

from ..config import RunConfig


def cache_path(cache_dir, config: RunConfig) -> Path:
    return Path(cache_dir) / f"{config.run_id()}.h5"


def try_load_cache(cache_dir, config: RunConfig):
    """Return a cached EnrichmentResult if one exists for this exact
    config, else None. Import of EnrichmentResult is local to avoid a
    circular import (results.py itself uses this module)."""
    from ..results import EnrichmentResult

    path = cache_path(cache_dir, config)
    if path.exists():
        return EnrichmentResult.load(path)
    return None