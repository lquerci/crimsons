from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class RunConfig:
    """Everything needed to reproduce a Simulation run.

    This is deliberately small and JSON-serializable (not the actual IMF
    or Channel objects) so it can be hashed into a stable run_id used for
    on-disk caching, and stored as metadata alongside results.
    """

    imf_name: str
    imf_mass_min: float
    imf_mass_max: float
    n_bins: int
    mass_formed: float
    metallicity: float
    n_realizations: int
    seed: int | None
    channel_names: tuple

    def __post_init__(self):
        # normalize to a tuple regardless of what was passed in (list,
        # generator, etc.) so hashing is stable
        object.__setattr__(self, "channel_names", tuple(self.channel_names))

    def run_id(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(payload.encode()).hexdigest()[:16]