from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import RunConfig
from .io.hdf5 import load_result, save_result


@dataclass
class EnrichmentResult:
    """Output of a Simulation run: N realizations of sampled stellar
    populations (as log-mass bins, not individual stars -- see
    IMF.sample_binned) and time-resolved cumulative element enrichment.
    """

    config: RunConfig
    time: np.ndarray
    elements: list
    enrichment: np.ndarray  # (n_realizations, n_time, n_elements)
    fates: list  # list of (n_bins_i,) object arrays, one per realization
    masses: list  # list of (n_bins_i,) float arrays -- mean mass per populated bin
    counts: list  # list of (n_bins_i,) float arrays -- stars represented by each bin

    def mean(self) -> np.ndarray:
        """Mean enrichment across realizations, shape (n_time, n_elements)."""
        return self.enrichment.mean(axis=0)

    def std(self) -> np.ndarray:
        """Std dev across realizations, shape (n_time, n_elements)."""
        return self.enrichment.std(axis=0)

    def save(self, path):
        save_result(self, Path(path))

    @classmethod
    def load(cls, path) -> EnrichmentResult:
        return load_result(cls, Path(path))