from __future__ import annotations

import numpy as np

from ..chemistry import HDF_COLUMNS, ZSUN, check_metallicity
from ..config import RunConfig
from ..enrichment.engine import bin_enrichment, run_realization
from ..imf.defaults import default_imf
from ..io.cache import cache_path, try_load_cache
from ..results import EnrichmentResult
from ..stars.lifetimes import StellarLifetime
from ..yields.channels import default_channels


class Simulation:
    """Ties an IMF, a lifetime function, and a set of channels together
    into N Monte Carlo realizations of a stellar population's chemical
    enrichment.

    Parameters
    ----------
    mass_formed : float, total stellar mass formed (Msun) per realization
    metallicity : float, metallicity of the formed stars. If negative is log Z/Z_sun, if positive is absolute Z
    imf : crimsons.imf.base.IMF, optional -- defaults to a population-
        appropriate crimsons.imf.defaults.default_imf(metallicity) if
        not given, so `Simulation(mass_formed=..., metallicity=...)` on
        its own already produces something physically reasonable
    lifetime_fn : crimsons.stars.lifetimes.LifetimeFunction, optional --
        defaults to StellarLifetime() if not given
    channels : list of crimsons.yields.base.Channel, optional -- defaults
        to crimsons.yields.channels.default_channels(metallicity) if not
        given (also population-appropriate)
    n_realizations : int
    seed : int or None -- a master seed; each realization gets its own
        independently-spawned child seed (via numpy's SeedSequence), so
        realizations are reproducible individually and as a set
    time_grid : array of times (Myr) at which enrichment is reported
    """

    def __init__(
        self,
        mass_formed: float,
        metallicity: float,
        imf=None,  # Optional: defaults to default_imf(metallicity) if None
        lifetime_fn=None,  # Optional: defaults to StellarLifetime() if None
        channels=None,  # Optional: defaults to default_channels(metallicity) if None
        n_realizations: int = 10,
        seed: int | None = None,
        time_grid=None,
        n_bins: int = 1000,
        sample_chunk_size: int = 1_000_000,
        verbose : bool = False,
    ):
        self.mass_formed = mass_formed
        self.n_realizations = n_realizations
        self.seed = seed
        self.time_grid = (
            np.asarray(time_grid) if time_grid is not None else np.logspace(0, 4, 200)
        )
        self.n_bins = n_bins
        self.sample_chunk_size = sample_chunk_size
        self.verbose = verbose

        # Fallback to standard defaults if not explicitly provided
        self.lifetime_fn = (
            lifetime_fn if lifetime_fn is not None else StellarLifetime()
        )

        # handling metallicity -- resolved before imf/channels below,
        # since both of their metallicity-adaptive defaults need it
        self.zsun = ZSUN

        if check_metallicity(metallicity=metallicity):
            self.metallicity = metallicity

        self.imf = imf if imf is not None else default_imf(self.metallicity)
        self.channels = (
            channels if channels is not None else default_channels(self.metallicity)
        )

        # define imf interval in the metallicity
        self.imf.bind_metallicity(self.metallicity)

        

    def config(self) -> RunConfig:
        
        return RunConfig(
            imf_name=type(self.imf).__name__,
            imf_mass_min=float(self.imf.m_min),
            imf_mass_max=float(self.imf.m_max),
            mass_formed=float(self.mass_formed),
            metallicity=float(self.metallicity),
            n_realizations=int(self.n_realizations),
            seed=self.seed,
            channel_names=tuple(c.name for c in self.channels),
            n_bins= int(self.n_bins)
        )

    def run(self, cache_dir=None, overwrite: bool = False) -> EnrichmentResult:
        cfg = self.config()

        if cache_dir is not None and not overwrite:
            cached = try_load_cache(cache_dir, cfg)
            if cached is not None:
                return cached

        seed_seq = np.random.SeedSequence(self.seed)
        child_seeds = seed_seq.spawn(self.n_realizations)

        all_masses, all_counts, all_fates, all_events = [], [], [], []
        enrichment = np.zeros((self.n_realizations, len(self.time_grid), len(HDF_COLUMNS)))

        # managing of the progress bar
        realization_iterable = self._get_progress_bar(
                range(len(child_seeds)), verbose=self.verbose
            )

        for i in realization_iterable:
            child_seed = child_seeds[i]
            rng = np.random.default_rng(child_seed)

            bin_masses, bin_counts, fates, events = run_realization(
                self.imf,
                self.mass_formed,
                self.metallicity,
                self.channels,
                self.lifetime_fn,
                self.time_grid,
                rng,
                n_bins=self.n_bins,
                chunk_size=self.sample_chunk_size,
            )
            all_masses.append(bin_masses)
            all_counts.append(bin_counts)
            all_fates.append(fates)
            enrichment[i] = bin_enrichment(events, self.time_grid, len(HDF_COLUMNS))

            # remove the yields to save memory, keeping only (name, times, counts)
            memory_safe_events = [(name, times, counts) for name, times, yields, counts in events]
            all_events.append(memory_safe_events)
        result = EnrichmentResult(
            config=cfg,
            time=self.time_grid,
            elements=list(HDF_COLUMNS),
            enrichment=enrichment,
            fates=all_fates,
            masses=all_masses,
            counts=all_counts,
            events_history=all_events,
        )

        if cache_dir is not None:
            result.save(cache_path(cache_dir, cfg))

        return result

    def _get_progress_bar(self, iterable, verbose: bool):
        """Helper to safely load tqdm or raise a clear error if missing."""
        if not verbose:
            return iterable

        try:
            from tqdm import tqdm

            return tqdm(
                iterable,
                desc="[CRIMSONS] Simulating",
                unit="realization",
                leave=True,
            )
        except ImportError:
            raise ImportError(
                "The 'tqdm' package is required for progress display (verbose=True).\n"
                "Please install it via 'pip install tqdm' or set 'verbose=False'."
            )