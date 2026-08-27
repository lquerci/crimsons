from __future__ import annotations

from importlib.resources import files

import numpy as np

from ..chemistry import POPIII_THRESHOLD, check_metallicity
from .base import MassRangeChannel, PopulationChannel
from .io import describe_model, list_models, load_yield_table_hdf5

_DATA_DIR = files("crimsons.yields") / "data"
_STELLAR_YIELDS_H5 = _DATA_DIR / "stellar_yields.h5"

# I generated a synthetic stellar_yields.h5 to build and test this
# against, since I don't have the real one -- these default model names
# match what's in that synthetic file (see
# scripts/generate_stellar_yields_h5.py). Swap them for whatever you
# actually want as defaults once this is pointed at your real file, and
# tell me if your model names or axis layout differ from what
# yields/io.py assumes -- I'll adjust the loader.
#
# SNII defaults to NK rather than LC deliberately: LC needs a rotation
# choice with no principled silent default (unlike an axis with only one
# grid value, "which rotation velocity" is a real physics choice this
# library shouldn't make quietly) -- pass model="LC",
# model_params={"rotation": ...} yourself when you want it.


def _regime_for(metallicity: float | None) -> str | None:
    """'Population III', 'Population II/I', or None if `metallicity`
    isn't known -- mirrors the split every bundled YieldTable already
    makes internally (see yields.base.POPIII_THRESHOLD)."""
    if metallicity is None:
        return None
    return "Population III" if metallicity < POPIII_THRESHOLD else "Population II/I"


def _hw_energy_distribution(rng, n):
    """SNII/HW's tabulated explosion-energy grid (1e51 erg units), drawn
    uniformly per star -- the default Population III
    model_params['energy'] (see _REGIME_DEFAULTS)."""
    energies = [3.0, 6.0, 12.0, 15.0, 18.0, 30.0, 50.0, 100.0]
    return rng.choice(energies, size=n)


def _hw_mixing_distribution(rng, n):
    """SNII/HW's tabulated mixing-parameter grid, drawn uniformly per
    star -- the default Population III model_params['mixing'] (see
    _REGIME_DEFAULTS)."""
    mixing = [39.8, 63.1, 100.0, 158.5]
    return rng.choice(mixing, size=n)


# Population-appropriate (model, model_params) per channel -- the single
# source of truth both default_channels(metallicity) and
# _resolve_model_defaults draw from, so a channel constructed directly
# with metallicity= (but no model=) always agrees with what
# default_channels(that same metallicity) would have picked.
_REGIME_DEFAULTS = {
    "SNII": {
        "Population III": (
            "HW",
            {"energy": _hw_energy_distribution, "mixing": _hw_mixing_distribution},
        ),
        "Population II/I": ("LC", {"rotation": 0}),
    },
    "AGB": {
        "Population III": ("MM", {}),
        "Population II/I": ("VAN", {}),
    },
    "PISN": {
        # PISN is physically a Population III / extremely metal-poor
        # phenomenon either way -- HW is used regardless of regime.
        "Population III": ("HW", {}),
        "Population II/I": ("HW", {}),
    },
    "SNIa": {
        "Population III": ("Iwamoto", {"model": "W70"}),
        "Population II/I": ("Iwamoto", {"model": "W7"}),
    },
}


def _resolve_model_defaults(channel: str, model, model_params, metallicity):
    """Fill in (model, model_params) for a channel constructor when
    `model` wasn't given explicitly.

    - `model` given: returned as-is (model_params just copied, untouched).
    - `model` not given, `metallicity` given: the population-appropriate
      default from `_REGIME_DEFAULTS`, with any caller-supplied
      model_params layered on top (caller values win).
    - neither given: assumes Population II/I and returns the appropriate
      '_REGIME_DEFAULTS'.
    """
    model_params = dict(model_params or {})
    if model is not None:
        return model, model_params

    regime = _regime_for(metallicity)

    # fall back to Population II/I if metallicity is not specified
    if regime is None: regime = "Population II/I"

    default_model, default_params = _REGIME_DEFAULTS[channel][regime]
    return default_model, {**default_params, **model_params}


def _default_yield_table(channel, model, model_params, h5_path):
    h5_path = h5_path or _STELLAR_YIELDS_H5
    return load_yield_table_hdf5(h5_path, channel, model, model_params=model_params)


def list_available_models(channel: str, h5_path=None) -> list:
    """List the models bundled for `channel` ('SNII', 'AGB', or 'PISN')
    in stellar_yields.h5, or in your own file via h5_path=."""
    return list_models(h5_path or _STELLAR_YIELDS_H5, channel)


def describe_available_model(channel: str, model: str, h5_path=None) -> dict:
    """Introspect a bundled model's axes and grid values -- e.g. to find
    out what model_params it needs. See yields.io.describe_model."""
    return describe_model(h5_path or _STELLAR_YIELDS_H5, channel, model)


class SNII(MassRangeChannel):
    """Core-collapse supernovae.

    `model` selects among the sources bundled in stellar_yields.h5 (see
    list_available_models('SNII')). `model_params` fixes any extra
    parameter axes the chosen model has (e.g. LC's rotation velocity) --
    required if that axis has more than one value; see
    describe_available_model('SNII', model).

    If `model` isn't given, pass `metallicity=` to get a population-
    appropriate default (HW for Population III, LC for Population II/I)
    with its required model_params already filled in, instead of this
    library's single metallicity-agnostic default (NK) -- see
    `default_channels`, which always does this.
    """

    def __init__(
        self,
        mass_min: float = 8.0,
        mass_max: float = 40.0,
        model: str | None = None,
        model_params: dict | None = None,
        metallicity: float | None = None,
        yield_table=None,
        h5_path=None,
    ):
        if metallicity: check_metallicity(metallicity=metallicity)
        self.model, model_params = _resolve_model_defaults("SNII", model, model_params, metallicity)
        table = yield_table or _default_yield_table("SNII", self.model, model_params, h5_path)
        super().__init__("SNII", mass_min, mass_max, table)


class AGB(MassRangeChannel):
    """Low/intermediate-mass stars enriching via stellar winds. See SNII
    for how `model`/`model_params`/`metallicity` work (default models:
    MM for Population III, VAN for Population II/I)."""

    def __init__(
        self,
        mass_min: float = 0.8,
        mass_max: float = 8.0,
        model: str | None = None,
        model_params: dict | None = None,
        metallicity: float | None = None,
        yield_table=None,
        h5_path=None,
    ):
        if metallicity: check_metallicity(metallicity=metallicity)
        self.model, model_params = _resolve_model_defaults("AGB", model, model_params, metallicity)
        table = yield_table or _default_yield_table("AGB", self.model, model_params, h5_path)
        super().__init__("AGB", mass_min, mass_max, table)


class PISN(MassRangeChannel):
    """Pair-instability supernovae: complete disruption (no compact
    remnant) of very massive, metal-free/extremely metal-poor stars in a
    narrow mass window -- effectively Population III only. Mass range is
    illustrative; adjust to your models. See SNII for how
    `model`/`model_params`/`metallicity` work.
    """

    def __init__(
        self,
        mass_min: float = 140.0,
        mass_max: float = 260.0,
        model: str | None = None,
        model_params: dict | None = None,
        metallicity: float | None = None,
        yield_table=None,
        h5_path=None,
    ):
        if metallicity: check_metallicity(metallicity=metallicity)
        self.model, model_params = _resolve_model_defaults("PISN", model, model_params, metallicity)
        table = yield_table or _default_yield_table("PISN", self.model, model_params, h5_path)
        super().__init__("PISN", mass_min, mass_max, table)


def _dtd_mannucci(tau, tau_min, tau_max):
    """Mannucci+06 / Matteucci+06 double-Gaussian-in-log-time DTD shape
    (unnormalized -- population_events normalizes it). tau, tau_min,
    tau_max in Myr; ported directly from the Fortran dtd_mannucci
    function, including its internal Myr -> yr conversion and the
    prompt/delayed split at t0 = 10^7.93 yr."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    out = np.zeros_like(tau)
    in_range = (tau >= tau_min) & (tau <= tau_max)
    if not np.any(in_range):
        return out

    t_yr = tau * 1.0e6
    log_t = np.log10(np.clip(t_yr, 1e-30, None))
    t0 = 10.0**7.93

    a1, b1, c1 = 1.4, -50.0, -7.7
    a2, b2, c2 = -0.8, -0.9, -8.7

    prompt = in_range & (t_yr <= t0)
    delayed = in_range & (t_yr > t0)
    out[prompt] = 10.0 ** (a1 + b1 * (log_t[prompt] + c1) ** 2)
    out[delayed] = 10.0 ** (a2 + b2 * (log_t[delayed] + c2) ** 2)
    return out


def _dtd_maoz(tau, min_time, slope: float = 1.12):
    """Maoz+12 power-law DTD shape (unnormalized). tau, min_time in Myr;
    ported directly from the Fortran dtd_maoz function."""
    tau = np.atleast_1d(np.asarray(tau, dtype=float))
    out = np.zeros_like(tau)
    mask = tau > min_time
    out[mask] = (tau[mask] / min_time) ** (-slope)
    return out


_DTD_SHAPES = {"mannucci": _dtd_mannucci, "maoz": _dtd_maoz}
# literature SNIa-per-Msun-formed rates, one per DTD shape -- matching
# the Fortran's SNIa_number_per_Msun for flag_SNIa_mode 2 and 3
_DEFAULT_RATE_PER_MSUN = {"mannucci": 0.0025, "maoz": 0.0013}


def _integerize_with_carry(expected_per_bin):
    """Convert a (generally fractional) expected-count-per-bin array
    into integer counts via deterministic remainder carry-over: each
    bin's leftover fraction is added to the next bin's expected count
    rather than rounded away, so the running total is conserved exactly.
    Ported directly from the Fortran's SNIa_remainder bookkeeping in
    SNaIa_chem_evolution. Not random -- same input always gives the same
    output, unlike a Poisson draw.
    """
    counts = np.zeros_like(expected_per_bin)
    remainder = 0.0
    for i, expected in enumerate(expected_per_bin):
        total = expected + remainder
        n = np.floor(total)
        remainder = total - n
        counts[i] = n
    return counts


class SNIa(PopulationChannel):
    """SN Ia enrichment from the population as a whole -- ported from a
    Fortran supernovae_Ia module that offered exactly two modes
    (flag_SNIa_mode):

    mode="dtd" (default): explosions follow a delay-time distribution
        SHAPE (dtd_shape="mannucci" (Mannucci+06/Matteucci+06, default)
        or "maoz" (Maoz+12)), normalized to integrate to 1 over the
        run's time_grid, then scaled by rate_per_msun (SNIa per Msun of
        stars formed -- defaults to the literature value for the chosen
        shape: 0.0025 for mannucci, 0.0013 for maoz) to get an absolute
        expected number of explosions per time bin. The DTD's support is
        bounded below/above by the progenitor mass range's (default
        0.8-8 Msun) lifetimes, via whatever lifetime_fn the Simulation
        uses.

    mode="single_burst": every eligible SNIa (rate_per_msun * mass_formed
        of them) explodes at one fixed delay after formation
        (burst_delay_myr) -- the discretized limit of a delta-function
        DTD, i.e. "a fixed fraction of the mass formed goes off as SN Ia
        at a specific time". rate_per_msun has no literature default in
        this mode (it's an inherently simplified treatment, calibrated
        per-model rather than from a population-integrated rate) and
        must be given explicitly.

    Either way, the (generally fractional) expected number of explosions
    per time bin is converted to an integer count via deterministic
    remainder carry-over, not a random draw -- see
    _integerize_with_carry, ported from the Fortran's SNIa_remainder.

    One simplification from the Fortran: explosion counts here scale
    directly with mass_formed * rate_per_msun, without the extra
    correction the original applied for how ITS specific fixed-mass-bin
    IMF discretization diverged from the analytic IMF integral over the
    progenitor range (number_stars_08_8 / number_stars / (P(...)-P(...))
    in initialize_DTD). This library's bins already track actual sampled
    star counts directly rather than a separately-tabulated analytic CDF,
    so that particular correction doesn't have an equivalent role here --
    say if you want it added back for closer numerical parity with the
    Fortran.

    Yields are a single fixed per-explosion composition (doesn't vary
    with progenitor mass, matching the Fortran's fixed SNIa_ele array)
    -- pass yield_table= for your own, or it loads the bundled
    snia_placeholder.csv.
    """

    def __init__(
        self,
        mode: str = "dtd",
        dtd_shape: str = "maoz",
        rate_per_msun: float | None = None,
        burst_delay_myr: float | None = None,
        yield_table=None,
        mass_min: float = 3.0,
        mass_max: float = 8.0,
        model: str | None = None,
        model_params: dict | None = None,
        metallicity: float | None = None,
        h5_path=None,
    ):
        
        if mode not in ("dtd", "single_burst"):
            raise ValueError(f"mode must be 'dtd' or 'single_burst', got {mode!r}")
        if mode == "dtd" and dtd_shape not in _DTD_SHAPES:
            raise ValueError(f"dtd_shape must be one of {list(_DTD_SHAPES)}, got {dtd_shape!r}")
        if mode == "single_burst" and burst_delay_myr is None:
            raise ValueError("mode='single_burst' needs burst_delay_myr")
        if mode == "single_burst" and rate_per_msun is None:
            raise ValueError(
                "mode='single_burst' has no literature-default rate_per_msun "
                "(unlike the dtd modes) -- pass the SNIa-per-Msun-formed rate you want"
            )

        if metallicity: check_metallicity(metallicity=metallicity)


        # fall back on the default model -- population-aware if
        # metallicity is given (see _resolve_model_defaults), else this
        # library's single metallicity-agnostic default
        model, model_params = _resolve_model_defaults("SNIa", model, model_params, metallicity)

        if model == "Iwamoto" and "model" not in model_params:
            # sub-choice among Iwamoto's W7/W70/... deflagration models;
            # W70 is the Population III-appropriate one (see
            # _REGIME_DEFAULTS), W7 the metallicity-agnostic default.
            model_params["model"] = "W70" if _regime_for(metallicity) == "Population III" else "W7"

        self.name = "SNIa"
        self.mode = mode
        self.model = model
        self.dtd_shape = dtd_shape
        self.rate_per_msun = (
            rate_per_msun if rate_per_msun is not None else _DEFAULT_RATE_PER_MSUN[dtd_shape]
        )
        self.progenitor_mass_range = (mass_min, mass_max)
        self.burst_delay_myr = burst_delay_myr
        #self._yield_table = yield_table or _default_snia_table()
        self._yield_table = yield_table or _default_yield_table("SNIa", self.model, model_params, h5_path)

    def yield_table(self):
        return self._yield_table

    def population_events(self, mass_formed, metallicity, lifetime_fn, time_grid, rng):
        time_grid = np.asarray(time_grid, dtype=float)
        dt = np.diff(time_grid)
        dt = np.append(dt, dt[-1]) if len(dt) else np.ones_like(time_grid)

        if self.mode == "single_burst":
            idx = min(int(np.searchsorted(time_grid, self.burst_delay_myr)), len(time_grid) - 1)
            expected = np.zeros_like(time_grid)
            expected[idx] = self.rate_per_msun * mass_formed
        else:
            mass_lo, mass_hi = self.progenitor_mass_range
            tau_min = float(lifetime_fn(np.array([mass_hi]), metallicity)[0])  # shorter-lived, more massive
            tau_max = float(lifetime_fn(np.array([mass_lo]), metallicity)[0])  # longer-lived, less massive

            shape_fn = _DTD_SHAPES[self.dtd_shape]
            shape = (
                shape_fn(time_grid, tau_min, tau_max)
                if self.dtd_shape == "mannucci"
                else shape_fn(time_grid, tau_min)
            )
            norm = np.sum(shape * dt)
            if norm <= 0:
                return None
            expected = mass_formed * self.rate_per_msun * (shape / norm) * dt

        counts = _integerize_with_carry(expected)
        mask = counts > 0
        if not np.any(mask):
            return None

        times = time_grid[mask]
        n_events = counts[mask]
        # yields don't depend on progenitor mass -- pass the midpoint of
        # the progenitor range as a formality the yield table still needs
        representative_mass = np.array([sum(self.progenitor_mass_range) / 2.0])
        one_event_yield = self.yield_table()(representative_mass, metallicity, rng=rng)
        return times, one_event_yield * n_events[:, None], n_events


def default_channels(metallicity):
    """Population-appropriate SNII + AGB + SNIa (+ PISN for Population
    III), for `metallicity` (absolute Z).

    Which regime `metallicity` falls in (see
    `crimsons.chemistry.POPIII_THRESHOLD`) decides both which models are
    used and, where needed, their model_params -- e.g. SNII's Population
    III default (HW) needs per-star 'energy'/'mixing' draws that its
    Population II/I default (LC's 'rotation') doesn't. See
    `_REGIME_DEFAULTS` for exactly what each regime uses; passing
    `metallicity=` to SNII/AGB/PISN/SNIa yourself (instead of `model=`)
    resolves the same way, so a channel list you build by hand can stay
    consistent with what this function would have picked.

    SN Ia here uses mode='single_burst' with rate_per_msun=0 in both
    regimes -- i.e. no SN Ia events at all -- as a placeholder until a
    real rate is calibrated; pass your own SNIa(...) in a hand-built
    `channels=` list if you want SN Ia actually contributing.
    """
    if _regime_for(metallicity) == "Population III":
        return [
            PISN(metallicity=metallicity),
            SNII(mass_min=10.0, mass_max=140.0, metallicity=metallicity),
            AGB(metallicity=metallicity),
            SNIa(
                mode="single_burst",
                burst_delay_myr=50.0,
                rate_per_msun=0.0,
                metallicity=metallicity,
            ),
        ]
    return [
        SNII(metallicity=metallicity),
        AGB(metallicity=metallicity),
        SNIa(
            mode="dtd",
            metallicity=metallicity,
        ),
    ]