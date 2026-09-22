"""Tests for Simulation's n_jobs parallel-realizations feature: the
n_jobs=None auto heuristic (_resolve_n_workers), explicit n_jobs values,
the picklability pre-check, and that parallel execution produces exactly
the same results as serial (same seeds -> same physics, regardless of
how the work is scheduled).
"""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from crimsons import SNII, ZSUN, Salpeter1955, Simulation


def _make_simulation(**overrides):
    kwargs = {
        "imf": Salpeter1955(),
        "mass_formed": 5000.0,
        "metallicity": ZSUN,
        "n_realizations": 6,
        "seed": 123,
        "time_grid": np.linspace(0, 13.8, 30),
    }
    kwargs.update(overrides)
    return Simulation(**kwargs)


# ---------------------------------------------------------------------------
# _resolve_n_workers: the n_jobs=None auto heuristic and explicit settings
# ---------------------------------------------------------------------------


class TestResolveNWorkers:
    def test_n_jobs_1_is_always_serial(self):
        sim = _make_simulation(n_jobs=1, n_realizations=50)
        with patch("os.cpu_count", return_value=8):
            assert sim._resolve_n_workers() == 1

    def test_auto_stays_serial_below_the_realization_threshold(self):
        sim = _make_simulation(
            n_jobs=None, n_realizations=Simulation.AUTO_PARALLEL_MIN_REALIZATIONS - 1
        )
        with patch("os.cpu_count", return_value=8):
            assert sim._resolve_n_workers() == 1

    def test_auto_parallelizes_at_the_realization_threshold(self):
        sim = _make_simulation(
            n_jobs=None, n_realizations=Simulation.AUTO_PARALLEL_MIN_REALIZATIONS
        )
        with patch("os.cpu_count", return_value=8):
            assert sim._resolve_n_workers() > 1

    def test_auto_stays_serial_on_a_single_cpu_machine(self):
        sim = _make_simulation(n_jobs=None, n_realizations=50)
        with patch("os.cpu_count", return_value=1):
            assert sim._resolve_n_workers() == 1

    def test_auto_never_exceeds_available_cpus(self):
        sim = _make_simulation(n_jobs=None, n_realizations=50)
        with patch("os.cpu_count", return_value=4):
            assert sim._resolve_n_workers() == 4

    def test_auto_never_exceeds_n_realizations(self):
        sim = _make_simulation(n_jobs=None, n_realizations=5)
        with patch("os.cpu_count", return_value=32):
            assert sim._resolve_n_workers() == 5

    def test_n_jobs_minus_one_uses_all_cpus(self):
        sim = _make_simulation(n_jobs=-1, n_realizations=50)
        with patch("os.cpu_count", return_value=6):
            assert sim._resolve_n_workers() == 6

    def test_n_jobs_minus_one_caps_at_n_realizations(self):
        sim = _make_simulation(n_jobs=-1, n_realizations=3)
        with patch("os.cpu_count", return_value=32):
            assert sim._resolve_n_workers() == 3

    def test_explicit_n_jobs_is_honored_even_above_cpu_count(self):
        # matches joblib/scikit-learn's n_jobs convention: an explicit
        # positive value is respected (with a warning), not silently
        # capped to cpu_count
        sim = _make_simulation(n_jobs=4, n_realizations=50)
        with (
            patch("os.cpu_count", return_value=1),
            pytest.warns(UserWarning, match="n_jobs=4"),
        ):
            assert sim._resolve_n_workers() == 4

    def test_explicit_n_jobs_still_caps_at_n_realizations(self):
        sim = _make_simulation(n_jobs=10, n_realizations=3)
        with patch("os.cpu_count", return_value=32):
            assert sim._resolve_n_workers() == 3

    def test_n_jobs_zero_is_invalid(self):
        sim = _make_simulation(n_jobs=0)
        with pytest.raises(ValueError, match="n_jobs"):
            sim._resolve_n_workers()

    def test_n_jobs_negative_below_minus_one_is_invalid(self):
        sim = _make_simulation(n_jobs=-2)
        with pytest.raises(ValueError, match="n_jobs"):
            sim._resolve_n_workers()


# ---------------------------------------------------------------------------
# n_jobs doesn't change the physics, or the cache identity
# ---------------------------------------------------------------------------


class TestNJobsDoesNotAffectConfigOrCaching:
    def test_config_does_not_depend_on_n_jobs(self):
        cfg_serial = _make_simulation(n_jobs=1).config()
        cfg_parallel = _make_simulation(n_jobs=4).config()
        assert cfg_serial.run_id() == cfg_parallel.run_id()

    def test_cache_written_by_one_n_jobs_setting_is_used_by_another(self, tmp_path):
        result_serial = _make_simulation(n_jobs=1, seed=7).run(cache_dir=tmp_path)
        # a differently-parallelized, but config-identical, Simulation
        # should hit the exact same cache entry rather than recomputing
        result_parallel = _make_simulation(n_jobs=2, seed=7).run(cache_dir=tmp_path)
        np.testing.assert_allclose(result_serial.enrichment, result_parallel.enrichment)
        assert result_serial.config.run_id() == result_parallel.config.run_id()


# ---------------------------------------------------------------------------
# Parallel execution produces identical results to serial
# ---------------------------------------------------------------------------


class TestParallelMatchesSerial:
    @pytest.fixture
    def serial_and_parallel_results(self):
        kwargs = {
            "mass_formed": 100_000.0, 
            "metallicity": ZSUN, 
            "n_realizations": 6, 
            "seed": 42
        }
        result_serial = Simulation(n_jobs=1, **kwargs).run()
        result_parallel = Simulation(n_jobs=2, **kwargs).run()
        return result_serial, result_parallel

    def test_enrichment_matches(self, serial_and_parallel_results):
        result_serial, result_parallel = serial_and_parallel_results
        np.testing.assert_allclose(result_serial.enrichment, result_parallel.enrichment)

    def test_energy_matches(self, serial_and_parallel_results):
        result_serial, result_parallel = serial_and_parallel_results
        np.testing.assert_allclose(result_serial.energy, result_parallel.energy)

    def test_masses_and_counts_match_per_realization(self, serial_and_parallel_results):
        result_serial, result_parallel = serial_and_parallel_results
        for a, b in zip(result_serial.masses, result_parallel.masses):
            np.testing.assert_allclose(a, b)
        for a, b in zip(result_serial.counts, result_parallel.counts):
            np.testing.assert_allclose(a, b)

    def test_fates_match_per_realization(self, serial_and_parallel_results):
        result_serial, result_parallel = serial_and_parallel_results
        for a, b in zip(result_serial.fates, result_parallel.fates):
            assert np.array_equal(a, b)

    def test_results_are_not_just_all_zero(self, serial_and_parallel_results):
        # a sanity check that this test would actually catch a real
        # divergence, not just two empty results agreeing trivially
        result_serial, _ = serial_and_parallel_results
        assert result_serial.enrichment.sum() > 0

    def test_a_realistic_50_realization_run(self):
        """The scenario from the feature request: many realizations,
        n_jobs handed a small worker count -- must complete and agree
        with the serial baseline."""
        kwargs = {
            "mass_formed":20_000.0, 
            "metallicity":ZSUN, 
            "n_realizations":50, 
            "seed":99, 
            "n_bins":200,
        }
        result_serial = Simulation(n_jobs=1, **kwargs).run()
        result_parallel = Simulation(n_jobs=2, **kwargs).run()
        assert result_parallel.enrichment.shape[0] == 50
        np.testing.assert_allclose(result_serial.enrichment, result_parallel.enrichment)


# ---------------------------------------------------------------------------
# Picklability pre-check
# ---------------------------------------------------------------------------


class TestPicklabilityGuard:
    def test_lambda_based_model_params_fails_fast_with_a_clear_error(self):
        bad_snii = SNII(
            model="Heger10",
            model_params={"energy": lambda rng, n: rng.uniform(1, 10, n), "mixing": 63.1},
            metallicity=1e-7,
        )
        sim = Simulation(
            mass_formed=1000.0,
            metallicity=1e-7,
            channels=[bad_snii],
            n_realizations=10,
            n_jobs=2,
        )
        with pytest.raises(Exception, match="picklable"):
            sim.run()

    def test_serial_path_does_not_require_picklability(self):
        # the same unpicklable channel must still work fine at n_jobs=1
        # -- the picklability requirement is specific to parallel execution
        bad_snii = SNII(
            model="Heger10",
            model_params={"energy": lambda rng, n: rng.uniform(1, 10, n), "mixing": 63.1},
            metallicity=1e-7,
        )
        sim = Simulation(
            mass_formed=1000.0,
            metallicity=1e-7,
            channels=[bad_snii],
            n_realizations=3,
            n_jobs=1,
        )
        result = sim.run()
        assert result.enrichment.shape[0] == 3

    def test_picklable_channels_pass_the_check(self):
        sim = _make_simulation(n_jobs=2)
        sim._check_picklable_for_parallel_run()  # must not raise


# ---------------------------------------------------------------------------
# Progress bar total (as_completed doesn't have __len__, so needs one)
# ---------------------------------------------------------------------------


class TestProgressBarTotal:
    def test_total_is_passed_through_to_tqdm(self):
        sim = _make_simulation()
        bar = sim._get_progress_bar(iter(range(5)), verbose=True, total=7)
        assert bar.total == 7
        bar.close()

    def test_serial_path_infers_total_from_range(self):
        sim = _make_simulation()
        bar = sim._get_progress_bar(range(9), verbose=True)
        assert bar.total == 9
        bar.close()

    def test_verbose_false_ignores_total(self):
        sim = _make_simulation()
        result = sim._get_progress_bar(range(5), verbose=False, total=999)
        assert result == range(5)