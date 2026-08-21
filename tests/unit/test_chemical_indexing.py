"""Tests for EnrichmentResult's __getitem__ (element and [X/Y] abundance
ratio access) and the supporting crimsons.chemistry solar-abundance
machinery.

Assumes a standard src-layout project with `pip install -e .` done, so
`import crimsons` resolves to your local checkout. Move this file to
wherever your test suite actually lives if that differs (e.g. under
`tests/unit/`).

Run with: pytest tests/test_results_indexing.py -v
Skip the (slower, real-simulation) integration test with:
    pytest tests/test_results_indexing.py -v -m "not slow"
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

import crimsons.chemistry as chem
from crimsons import ELEMENTS, EnrichmentResult, SolarAbundances
from crimsons.chemistry import (
    ATOMIC_WEIGHTS,
    default_solar_abundances,
    load_solar_abundances,
    set_default_solar_abundances,
)
from crimsons.config import RunConfig

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_config(n_realizations: int) -> RunConfig:
    return RunConfig(
        imf_name="Kroupa2001",
        imf_mass_min=0.1,
        imf_mass_max=100.0,
        n_bins=10,
        mass_formed=1.0,
        metallicity=0.02,
        n_realizations=n_realizations,
        seed=0,
        channel_names=("SNII",),
    )


def _make_result(enrichment, elements, time=None) -> EnrichmentResult:
    """Build a minimal EnrichmentResult around a hand-picked enrichment
    array, so __getitem__ can be tested without running a full
    Simulation. fates/masses/counts/events_history aren't exercised by
    indexing, so they're left as empty per-realization placeholders.
    """
    enrichment = np.asarray(enrichment, dtype=float)
    n_real, n_time, n_elem = enrichment.shape
    assert n_elem == len(elements)
    if time is None:
        time = np.linspace(1.0, 10.0, n_time)
    return EnrichmentResult(
        config=_make_config(n_real),
        time=time,
        elements=list(elements),
        enrichment=enrichment,
        fates=[np.array([]) for _ in range(n_real)],
        masses=[np.array([]) for _ in range(n_real)],
        counts=[np.array([]) for _ in range(n_real)],
        events_history=[[] for _ in range(n_real)],
    )


@pytest.fixture
def simple_result() -> EnrichmentResult:
    """2 realizations x 4 time steps x 5 elements (H, He, C, O, Fe), all
    strictly positive from the first time step -- no NaN edge cases, good
    for exact-value checks.
    """
    elements = ["H", "He", "C", "O", "Fe"]
    enrichment = np.zeros((2, 4, len(elements)))
    enrichment[0, :, elements.index("H")] = [10.0, 20.0, 30.0, 40.0]
    enrichment[0, :, elements.index("He")] = [1.0, 2.0, 3.0, 4.0]
    enrichment[0, :, elements.index("C")] = [0.10, 0.20, 0.30, 0.40]
    enrichment[0, :, elements.index("O")] = [0.50, 0.60, 0.70, 0.80]
    enrichment[0, :, elements.index("Fe")] = [0.010, 0.020, 0.030, 0.040]
    enrichment[1, :, elements.index("H")] = [12.0, 22.0, 32.0, 42.0]
    enrichment[1, :, elements.index("He")] = [1.2, 2.2, 3.2, 4.2]
    enrichment[1, :, elements.index("C")] = [0.12, 0.22, 0.32, 0.42]
    enrichment[1, :, elements.index("O")] = [0.55, 0.65, 0.75, 0.85]
    enrichment[1, :, elements.index("Fe")] = [0.012, 0.022, 0.032, 0.042]
    return _make_result(enrichment, elements)


@pytest.fixture
def result_with_delayed_iron() -> EnrichmentResult:
    """Fe is exactly zero (not yet produced) for the first 2 of 5 time
    steps in *every* realization, and zero in only *one* of two
    realizations at the 3rd step -- exercises both the "every realization
    undefined" and "some realizations undefined" NaN cases for [X/Fe].
    """
    elements = ["H", "C", "Fe"]
    enrichment = np.zeros((2, 5, len(elements)))
    enrichment[:, :, elements.index("H")] = 10.0  # never zero
    enrichment[:, :, elements.index("C")] = 1.0  # never zero
    enrichment[:, :, elements.index("Fe")] = [
        [0.00, 0.00, 0.00, 0.05, 0.06],  # realization 0: Fe appears at t=3
        [0.00, 0.00, 0.02, 0.03, 0.04],  # realization 1: Fe appears at t=2
    ]
    return _make_result(enrichment, elements)


@pytest.fixture
def toy_solar() -> SolarAbundances:
    """A simple, hand-computable solar reference table (mass-fraction
    style, so the numbers below are the *actual* ratios used -- no
    atomic-weight conversion involved) for checking ratio math by hand
    rather than trusting the bundled Asplund+09 numbers.
    """
    proxy = {el: 1.0 for el in ELEMENTS}  # harmless filler for unused elements
    proxy.update({"H": 1.0, "He": 0.1, "C": 0.01, "O": 0.02, "Fe": 0.001})
    return SolarAbundances.from_mass_fractions(proxy)


@pytest.fixture
def restore_default_solar():
    """Any test that calls set_default_solar_abundances mutates
    module-level state in crimsons.chemistry -- reach in and restore it
    afterwards so tests stay order-independent.
    """
    original = chem._default_solar_abundances
    yield
    chem._default_solar_abundances = original


# ---------------------------------------------------------------------------
# Element access: result['C']
# ---------------------------------------------------------------------------


class TestElementAccess:
    def test_mean_matches_manual_indexing(self, simple_result):
        idx = simple_result.elements.index("Fe")
        manual = simple_result.mean()[:, idx]
        np.testing.assert_allclose(simple_result["Fe"].mean(), manual)

    def test_std_matches_manual_indexing(self, simple_result):
        idx = simple_result.elements.index("C")
        manual = simple_result.std()[:, idx]
        np.testing.assert_allclose(simple_result["C"].std(), manual)

    def test_every_tracked_element_is_reachable(self, simple_result):
        for el in simple_result.elements:
            idx = simple_result.elements.index(el)
            np.testing.assert_allclose(
                simple_result[el].mean(), simple_result.mean()[:, idx]
            )

    def test_values_shape_and_len(self, simple_result):
        series = simple_result["H"]
        assert series.values.shape == (2, 4)
        assert len(series) == 2

    def test_label_is_the_lookup_key(self, simple_result):
        assert simple_result["Fe"].label == "Fe"

    def test_time_is_shared_not_copied(self, simple_result):
        assert simple_result["Fe"].time is simple_result.time

    def test_unknown_element_raises_keyerror_listing_available(self, simple_result):
        with pytest.raises(KeyError, match="Unknown element 'Xx'"):
            simple_result["Xx"]

    def test_non_string_key_raises_typeerror(self, simple_result):
        with pytest.raises(TypeError):
            simple_result[0]

    def test_slash_without_brackets_is_not_treated_as_a_ratio(self, simple_result):
        # "C/Fe" (no brackets) must NOT be silently interpreted as a
        # ratio -- it's an invalid (nonexistent) element symbol.
        with pytest.raises(KeyError):
            simple_result["C/Fe"]

    def test_contains(self, simple_result):
        assert "Fe" in simple_result
        assert "Xx" not in simple_result


# ---------------------------------------------------------------------------
# Abundance ratios: result['[X/Y]'] -- including denominators other than Fe
# ---------------------------------------------------------------------------


class TestAbundanceRatios:
    def test_ratio_matches_hand_computed_value(
        self, simple_result, toy_solar, restore_default_solar
    ):
        set_default_solar_abundances(toy_solar)
        i_c = simple_result.elements.index("C")
        i_fe = simple_result.elements.index("Fe")
        m_c = simple_result.enrichment[:, :, i_c]
        m_fe = simple_result.enrichment[:, :, i_fe]
        expected = np.log10(m_c / m_fe) - np.log10(toy_solar.ratio("C", "Fe"))
        np.testing.assert_allclose(simple_result["[C/Fe]"].values, expected)

    @pytest.mark.parametrize(
        "numerator,denominator",
        [
            ("C", "Fe"),
            ("O", "Fe"),
            ("C", "O"),   # <- answers "does [X/O] work?" -- yes
            ("Fe", "O"),  # <- any element can be the denominator
            ("O", "H"),
            ("He", "C"),
        ],
    )
    def test_ratio_is_not_hardcoded_to_any_one_denominator(
        self, simple_result, toy_solar, restore_default_solar, numerator, denominator
    ):
        """The [X/Y] machinery doesn't special-case Fe (or any other
        element) as the denominator -- any two tracked elements work.
        """
        set_default_solar_abundances(toy_solar)
        key = f"[{numerator}/{denominator}]"
        i_num = simple_result.elements.index(numerator)
        i_den = simple_result.elements.index(denominator)
        m_num = simple_result.enrichment[:, :, i_num]
        m_den = simple_result.enrichment[:, :, i_den]
        expected = np.log10(m_num / m_den) - np.log10(
            toy_solar.ratio(numerator, denominator)
        )
        np.testing.assert_allclose(simple_result[key].values, expected)

    def test_ratio_antisymmetry(self, simple_result, toy_solar, restore_default_solar):
        set_default_solar_abundances(toy_solar)
        cfe = simple_result["[C/Fe]"].values
        fec = simple_result["[Fe/C]"].values
        # atol matters here: some entries land at/near exactly zero, where
        # rtol alone is meaningless (relative diff blows up on ~1e-16 noise)
        np.testing.assert_allclose(cfe, -fec, atol=1e-12)

    def test_unknown_element_inside_brackets_raises_keyerror(self, simple_result):
        with pytest.raises(KeyError, match="Unknown element 'Xx'"):
            simple_result["[Xx/Fe]"]
        with pytest.raises(KeyError, match="Unknown element 'Xx'"):
            simple_result["[Fe/Xx]"]

    def test_contains_for_ratios(self, simple_result):
        assert "[C/Fe]" in simple_result
        assert "[C/O]" in simple_result
        assert "[Xx/Fe]" not in simple_result
        assert "[C/Xx]" not in simple_result

    def test_malformed_brackets_are_not_a_ratio(self, simple_result):
        # missing closing bracket, missing slash, extra element -- none
        # of these should be misparsed as a valid ratio; they should
        # fail as an (invalid) plain element lookup instead.
        for bad_key in ["[C/Fe", "C/Fe]", "[C]", "[C/Fe/H]"]:
            with pytest.raises(KeyError):
                simple_result[bad_key]


class TestNaNHandling:
    def test_undefined_before_any_yield_is_nan(self, result_with_delayed_iron):
        series = result_with_delayed_iron["[C/Fe]"]
        assert np.isnan(series.values[:, 0]).all()  # t=0: both realizations Fe=0
        assert np.isnan(series.values[:, 1]).all()  # t=1: both realizations Fe=0

    def test_all_nan_timestep_gives_nan_mean_without_warning(
        self, result_with_delayed_iron
    ):
        series = result_with_delayed_iron["[C/Fe]"]
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # any warning becomes a failure
            mean = series.mean()
            std = series.std()
        assert np.isnan(mean[0])
        assert np.isnan(std[0])

    def test_partially_defined_timestep_excludes_nan_realizations(
        self, result_with_delayed_iron
    ):
        series = result_with_delayed_iron["[C/Fe]"]
        # at t=2: realization 0 has Fe=0 (nan), realization 1 has Fe=0.02 (valid)
        assert np.isnan(series.values[0, 2])
        assert not np.isnan(series.values[1, 2])
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            mean_t2 = series.mean()[2]
        np.testing.assert_allclose(mean_t2, series.values[1, 2])

    def test_plain_element_access_never_produces_nan(self, result_with_delayed_iron):
        # a channel legitimately returning zero mass (not yet produced)
        # is a real value for a single element -- only *ratios* should
        # turn that into NaN.
        assert not np.isnan(result_with_delayed_iron["Fe"].values).any()


# ---------------------------------------------------------------------------
# crimsons.chemistry: SolarAbundances, load_solar_abundances
# ---------------------------------------------------------------------------


class TestSolarAbundances:
    def test_missing_elements_raises_value_error(self):
        with pytest.raises(ValueError, match="missing elements"):
            SolarAbundances({"H": 1.0, "He": 1.0})

    def test_ratio_is_reciprocal(self, toy_solar):
        assert toy_solar.ratio("C", "Fe") == pytest.approx(
            1.0 / toy_solar.ratio("Fe", "C")
        )

    def test_unknown_element_mass_ratio_raises_keyerror(self, toy_solar):
        with pytest.raises(KeyError):
            toy_solar.mass_ratio("Xx")

    def test_from_log_eps_uses_atomic_weights(self):
        # A(C) - A(Fe) = 1.0 dex of *number* ratio; converting to a mass
        # ratio must scale by the atomic-weight ratio.
        log_eps = {el: 5.0 for el in ELEMENTS}
        log_eps["C"] = 8.50
        log_eps["Fe"] = 7.50
        table = SolarAbundances.from_log_eps(log_eps)
        expected = 10 ** (8.50 - 7.50) * (ATOMIC_WEIGHTS["C"] / ATOMIC_WEIGHTS["Fe"])
        assert table.ratio("C", "Fe") == pytest.approx(expected, rel=1e-9)

    def test_from_mass_fractions_ignores_atomic_weights(self):
        # a *mass*-based reference is used as-is; the ratio must come out
        # exactly equal to the input ratio, with no atomic-weight factor.
        mf = {el: 1.0 for el in ELEMENTS}
        mf["C"] = 4.0
        mf["Fe"] = 2.0
        table = SolarAbundances.from_mass_fractions(mf)
        assert table.ratio("C", "Fe") == pytest.approx(2.0)


class TestLoadSolarAbundances:
    def test_log_eps_format_auto_detected(self, tmp_path):
        rows = {el: 5.0 for el in ELEMENTS}
        rows.update({"H": 12.00, "He": 10.93, "C": 8.43, "O": 8.69, "Fe": 7.50})
        path = tmp_path / "solar_logeps.csv"
        path.write_text(
            "element,A\n" + "\n".join(f"{el},{val}" for el, val in rows.items())
        )
        table = load_solar_abundances(path)
        expected = 10 ** (8.43 - 7.50) * (ATOMIC_WEIGHTS["C"] / ATOMIC_WEIGHTS["Fe"])
        assert table.ratio("C", "Fe") == pytest.approx(expected, rel=1e-9)

    def test_mass_fraction_format_auto_detected(self, tmp_path):
        rows = {el: 1e-6 for el in ELEMENTS}
        rows.update({"H": 0.70, "He": 0.28, "O": 6.0e-3, "Fe": 1.0e-3})
        path = tmp_path / "solar_mass.csv"
        path.write_text(
            "element,mass_fraction\n"
            + "\n".join(f"{el},{val}" for el, val in rows.items())
        )
        table = load_solar_abundances(path)
        assert table.ratio("O", "Fe") == pytest.approx(6.0e-3 / 1.0e-3)

    def test_explicit_column_overrides_autodetect(self, tmp_path):
        rows = {el: 1.0 for el in ELEMENTS}
        rows["C"] = 2.0
        path = tmp_path / "solar_two_cols.csv"
        lines = ["element,mass_fraction,notes"]
        lines += [f"{el},{val},ok" for el, val in rows.items()]
        path.write_text("\n".join(lines))
        table = load_solar_abundances(path, column="mass_fraction")
        assert table.ratio("C", "H") == pytest.approx(2.0)

    def test_comment_lines_before_header_are_skipped(self, tmp_path):
        # Regression test: an earlier version choked on '#'-prefixed
        # citation/comment lines before the header, misreading one of
        # them as the header row itself.
        rows = {el: 1.0 for el in ELEMENTS}
        path = tmp_path / "solar_with_comments.csv"
        path.write_text(
            "# Some Paper et al. (2009)\n"
            "# second comment line\n"
            "element,mass_fraction\n"
            + "\n".join(f"{el},{val}" for el, val in rows.items())
        )
        table = load_solar_abundances(path)  # must not raise
        assert table.ratio("H", "He") == pytest.approx(1.0)

    def test_no_data_rows_raises_value_error(self, tmp_path):
        path = tmp_path / "empty.csv"
        path.write_text("element,mass_fraction\n")
        with pytest.raises(ValueError):
            load_solar_abundances(path)


class TestBundledDefault:
    def test_covers_every_tracked_element(self):
        table = default_solar_abundances()
        for el in ELEMENTS:
            ratio = table.mass_ratio(el)
            assert np.isfinite(ratio) and ratio > 0

    def test_iron_is_a_trace_element_relative_to_hydrogen(self):
        # sanity check on the bundled Asplund+09 numbers: Fe is a trace
        # species by mass relative to H in the sun.
        table = default_solar_abundances()
        assert table.ratio("Fe", "H") < 0.01


class TestSetDefaultSolarAbundances:
    def test_override_changes_ratio_lookup(
        self, simple_result, toy_solar, restore_default_solar
    ):
        before = simple_result["[C/Fe]"].values.copy()
        set_default_solar_abundances(toy_solar)
        after = simple_result["[C/Fe]"].values
        assert not np.allclose(before, after)

    def test_accepts_a_path_directly(self, tmp_path, restore_default_solar):
        rows = {el: 1.0 for el in ELEMENTS}
        rows["C"] = 3.0
        path = tmp_path / "solar.csv"
        path.write_text(
            "element,mass_fraction\n"
            + "\n".join(f"{el},{val}" for el, val in rows.items())
        )
        set_default_solar_abundances(path)
        assert default_solar_abundances().ratio("C", "H") == pytest.approx(3.0)


# ---------------------------------------------------------------------------
# Integration: a real (small) Simulation, not a hand-built EnrichmentResult
# ---------------------------------------------------------------------------


@pytest.mark.slow
class TestIntegrationWithRealSimulation:
    def test_element_and_ratio_access_on_a_real_run(self):
        from crimsons import Kroupa2001, Simulation

        sim = Simulation(
            imf=Kroupa2001(),
            mass_formed=2e5,
            metallicity=-1.0,
            n_realizations=4,
            seed=123,
        )
        result = sim.run()

        # every tracked element is reachable and matches manual indexing
        for el in result.elements:
            idx = result.elements.index(el)
            np.testing.assert_allclose(
                result[el].mean(), result.mean()[:, idx]
            )

        # ratios beyond Fe-as-denominator all come back finite by the
        # final time step, with the expected shape
        for key in ["[C/Fe]", "[O/Fe]", "[Mg/Fe]", "[C/O]", "[Fe/O]"]:
            series = result[key]
            assert series.values.shape == result.enrichment.shape[:2]
            assert np.isfinite(series.mean()[-1])