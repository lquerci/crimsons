from __future__ import annotations

import os
import pickle
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

from ..chemistry import ELEMENTS, HDF_COLUMNS, ZSUN, check_metallicity
from ..config import RunConfig
from ..enrichment.engine import run_and_bin_realization
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
        realizations are reproducible individually and as a set --
        *regardless* of n_jobs (see n_jobs below)
    time_grid : array of times (Myr) at which enrichment is reported
    n_jobs : int or None, optional -- how many realizations to run at
        once, in separate worker processes:

        - None (default): auto. Realizations run in parallel only if
          there's more than one CPU available *and* at least
          `Simulation.AUTO_PARALLEL_MIN_REALIZATIONS` (4) of them to run
          -- below that, process-pool startup tends to cost more than
          serial execution saves. This is a simple realization-count
          heuristic, not a real profile of your workload: if a handful
          of realizations is still worth parallelizing for you (e.g. a
          very large mass_formed and/or n_bins making each one
          individually expensive), set n_jobs explicitly instead of
          relying on auto.
        - 1: force serial (a plain for loop, no multiprocessing at all
          -- the previous, only, behavior).
        - -1: use every available CPU.
        - a positive int: use that many worker processes (capped at
          n_realizations only -- if you ask for more than the number of
          CPUs available you'll get a warning, but it's still honored,
          matching joblib/scikit-learn's n_jobs convention).

        Every realization produces exactly the same result regardless
        of n_jobs -- it only changes how the work is scheduled, not the
        physics, so it isn't part of the run's cache identity
        (`config()`/`RunConfig`): switching n_jobs between runs of an
        otherwise-identical Simulation still hits the same cache entry.

        Caveats: parallelizing needs `imf`, `channels`, and `lifetime_fn`
        to be picklable -- true for everything built into this library,
        but a custom yield distribution, IMF shape function, or channel
        that uses a lambda or other local closure won't pickle (use a
        module-level function or a callable class instead). This is
        checked up front (before spawning any worker), so you get a
        clear error naming the problem rather than a run half-started.
        On Windows specifically, using n_jobs != 1 from a script (not a
        notebook) requires guarding the entry point with
        `if __name__ == "__main__":`, a standard Python multiprocessing
        requirement.
    verbose : bool -- show a progress bar (requires tqdm). Reflects
        actual completions either way: realizations finishing as they
        finish when n_jobs > 1, not submission order.
    """

    # see n_jobs above -- the auto heuristic's realization-count cutoff
    AUTO_PARALLEL_MIN_REALIZATIONS = 4

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
        n_jobs: int | None = None,
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
        self.n_jobs = n_jobs
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

    def _resolve_n_workers(self) -> int:
        """How many worker processes `run()` should actually use --
        1 means "just use a plain for loop", including whenever the
        n_jobs=None auto heuristic decides against parallelizing. See
        the n_jobs parameter docs on the class for what each setting
        means.
        """
        cpu_count = os.cpu_count() or 1

        if self.n_jobs is None:
            if cpu_count > 1 and self.n_realizations >= self.AUTO_PARALLEL_MIN_REALIZATIONS:
                return min(self.n_realizations, cpu_count)
            return 1

        if self.n_jobs == -1:
            return min(self.n_realizations, cpu_count)

        if self.n_jobs < 1:
            raise ValueError(f"n_jobs must be a positive int, -1, or None, got {self.n_jobs!r}")

        if self.n_jobs > cpu_count:
            warnings.warn(
                f"n_jobs={self.n_jobs} is more than the {cpu_count} CPU(s) "
                "available on this machine -- using that many worker "
                "processes anyway, but you likely won't see speedup past "
                "cpu_count workers",
                stacklevel=3,
            )
        return min(self.n_jobs, self.n_realizations)

    def _check_picklable_for_parallel_run(self):
        """Fail fast with a clear error if imf/channels/lifetime_fn can't
        be pickled, rather than a confusing traceback from deep inside
        concurrent.futures after a worker process has already been
        spawned."""
        try:
            pickle.dumps((self.imf, self.channels, self.lifetime_fn))
        except (pickle.PicklingError, AttributeError, TypeError) as exc:
            raise pickle.PicklingError(
                "Simulation.imf/channels/lifetime_fn must be picklable to use "
                "n_jobs != 1 (parallel realizations run in separate worker "
                "processes). A custom yield distribution, IMF shape function, "
                "or channel that uses a lambda or other local closure won't "
                "pickle -- use a module-level function or a callable class "
                f"instead. Set n_jobs=1 to run serially without this "
                f"restriction. Original error: {exc}"
            ) from exc

    def run(self, cache_dir=None, overwrite: bool = False) -> EnrichmentResult:
        cfg = self.config()

        if cache_dir is not None and not overwrite:
            cached = try_load_cache(cache_dir, cfg)
            if cached is not None:
                return cached

        seed_seq = np.random.SeedSequence(self.seed)
        child_seeds = seed_seq.spawn(self.n_realizations)
        n_columns = len(HDF_COLUMNS)

        n_workers = self._resolve_n_workers()
        if n_workers > 1:
            self._check_picklable_for_parallel_run()
            all_masses, all_counts, all_fates, all_events, enrichment = self._run_parallel(
                child_seeds, n_columns, n_workers
            )
        else:
            all_masses, all_counts, all_fates, all_events, enrichment = self._run_serial(
                child_seeds, n_columns
            )

        # split the explosion-energy column out of the raw (n_realizations,
        # n_time, n_elements+1) array bin_enrichment produced -- everything
        # downstream (EnrichmentResult.elements/enrichment) is chemistry
        # only, energy is tracked separately (see EnrichmentResult).
        energy_idx = HDF_COLUMNS.index("Energy")
        energy = enrichment[:, :, energy_idx]
        chem_enrichment = np.delete(enrichment, energy_idx, axis=2)

        result = EnrichmentResult(
            config=cfg,
            time=self.time_grid,
            elements=list(ELEMENTS),
            enrichment=chem_enrichment,
            energy=energy,
            fates=all_fates,
            masses=all_masses,
            counts=all_counts,
            events_history=all_events,
        )

        if cache_dir is not None:
            result.save(cache_path(cache_dir, cfg))

        return result

    def _run_serial(self, child_seeds, n_columns):
        """Plain for-loop realization path -- used when n_jobs resolves
        to 1 (explicitly, or via the n_jobs=None auto heuristic)."""
        all_masses, all_counts, all_fates, all_events = [], [], [], []
        enrichment = np.zeros((self.n_realizations, len(self.time_grid), n_columns))

        realization_iterable = self._get_progress_bar(
            range(self.n_realizations), verbose=self.verbose
        )

        for i in realization_iterable:
            bin_masses, bin_counts, fates, memory_safe_events, enrichment_row = (
                run_and_bin_realization(
                    self.imf,
                    self.mass_formed,
                    self.metallicity,
                    self.channels,
                    self.lifetime_fn,
                    self.time_grid,
                    child_seeds[i],
                    n_bins=self.n_bins,
                    chunk_size=self.sample_chunk_size,
                    n_columns=n_columns,
                )
            )
            all_masses.append(bin_masses)
            all_counts.append(bin_counts)
            all_fates.append(fates)
            all_events.append(memory_safe_events)
            enrichment[i] = enrichment_row

        return all_masses, all_counts, all_fates, all_events, enrichment

    def _run_parallel(self, child_seeds, n_columns, n_workers):
        """Realizations dispatched to `n_workers` worker processes, each
        running run_and_bin_realization for one realization -- see the
        n_jobs parameter docs on the class. Results are collected as
        each realization actually finishes (not submission order), so
        the progress bar (verbose=True) reflects real completions."""
        all_masses = [None] * self.n_realizations
        all_counts = [None] * self.n_realizations
        all_fates = [None] * self.n_realizations
        all_events = [None] * self.n_realizations
        enrichment = np.zeros((self.n_realizations, len(self.time_grid), n_columns))

        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_index = {
                executor.submit(
                    run_and_bin_realization,
                    self.imf,
                    self.mass_formed,
                    self.metallicity,
                    self.channels,
                    self.lifetime_fn,
                    self.time_grid,
                    child_seeds[i],
                    self.n_bins,
                    self.sample_chunk_size,
                    n_columns,
                ): i
                for i in range(self.n_realizations)
            }

            completed = self._get_progress_bar(
                as_completed(future_to_index),
                verbose=self.verbose,
                total=self.n_realizations,
            )

            for future in completed:
                i = future_to_index[future]
                bin_masses, bin_counts, fates, memory_safe_events, enrichment_row = (
                    future.result()
                )
                all_masses[i] = bin_masses
                all_counts[i] = bin_counts
                all_fates[i] = fates
                all_events[i] = memory_safe_events
                enrichment[i] = enrichment_row

        return all_masses, all_counts, all_fates, all_events, enrichment

    def _get_progress_bar(self, iterable, verbose: bool, total: int | None = None):
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
                total=total,
            )
        except ImportError:
            raise ImportError(
                "The 'tqdm' package is required for progress display (verbose=True).\n"
                "Please install it via 'pip install tqdm' or set 'verbose=False'."
            )