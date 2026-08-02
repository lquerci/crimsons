from __future__ import annotations

from importlib.resources import files

import numpy as np

from .base import Channel, MassRangeChannel
from .io import describe_model, list_models, load_yield_table, load_yield_table_hdf5

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
DEFAULT_MODELS = {"SNII": "NK", "AGB": "VAN", "PISN": "HW"}


def _default_yield_table(channel, model, model_params, h5_path):
    h5_path = h5_path or _STELLAR_YIELDS_H5
    return load_yield_table_hdf5(h5_path, channel, model, model_params=model_params)


def _default_snia_table():
    return load_yield_table(_DATA_DIR / "snia_placeholder.csv")


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
    """

    def __init__(
        self,
        mass_min: float = 8.0,
        mass_max: float = 40.0,
        model: str | None = None,
        model_params: dict | None = None,
        yield_table=None,
        h5_path=None,
    ):
        self.model = model or DEFAULT_MODELS["SNII"]
        table = yield_table or _default_yield_table("SNII", self.model, model_params, h5_path)
        super().__init__("SNII", mass_min, mass_max, table)


class AGB(MassRangeChannel):
    """Low/intermediate-mass stars enriching via stellar winds. See SNII
    for how `model`/`model_params` work."""

    def __init__(
        self,
        mass_min: float = 0.8,
        mass_max: float = 8.0,
        model: str | None = None,
        model_params: dict | None = None,
        yield_table=None,
        h5_path=None,
    ):
        self.model = model or DEFAULT_MODELS["AGB"]
        table = yield_table or _default_yield_table("AGB", self.model, model_params, h5_path)
        super().__init__("AGB", mass_min, mass_max, table)


class PISN(MassRangeChannel):
    """Pair-instability supernovae: complete disruption (no compact
    remnant) of very massive, metal-free/extremely metal-poor stars in a
    narrow mass window -- effectively Population III only. Mass range is
    illustrative; adjust to your models. See SNII for how
    `model`/`model_params` work.
    """

    def __init__(
        self,
        mass_min: float = 140.0,
        mass_max: float = 260.0,
        model: str | None = None,
        model_params: dict | None = None,
        yield_table=None,
        h5_path=None,
    ):
        self.model = model or DEFAULT_MODELS["PISN"]
        table = yield_table or _default_yield_table("PISN", self.model, model_params, h5_path)
        super().__init__("PISN", mass_min, mass_max, table)


class SNIa(Channel):
    """SN Ia progenitors via a binary fraction and a delay-time
    distribution (DTD), following the standard approach in stochastic
    chemical evolution models (e.g. Matteucci & Greggio 1986; Greggio
    2005) -- since SN Ia progenitors are white dwarfs in binaries, not
    single stars, this is NOT a MassRangeChannel like SNII/AGB/PISN.

    The progenitor mass range, binary fraction, and DTD shape/normalization
    here are illustrative placeholders, as is the bundled CSV yield table.
    This channel's whole approach -- including how it interacts with the
    log-mass-binned sampling in the engine -- is a known TODO to rework
    separately; it's unaffected by the SNII/AGB/PISN changes here.
    """

    def __init__(
        self,
        progenitor_mass_range=(3.0, 8.0),
        binary_fraction: float = 0.05,
        dtd_power: float = 1.0,
        dtd_min_delay_gyr: float = 0.04,
        dtd_max_delay_gyr: float = 13.0,
        yield_table=None,
    ):
        self.name = "SNIa"
        self.progenitor_mass_range = progenitor_mass_range
        self.binary_fraction = binary_fraction
        self.dtd_power = dtd_power
        self.dtd_min_delay_gyr = dtd_min_delay_gyr
        self.dtd_max_delay_gyr = dtd_max_delay_gyr
        self._yield_table = yield_table or _default_snia_table()

    def contributes(self, mass, metallicity, rng):
        lo, hi = self.progenitor_mass_range
        in_range = (mass >= lo) & (mass <= hi)
        is_binary = rng.random(mass.shape) < self.binary_fraction
        return in_range & is_binary

    def delay_time(self, mass, metallicity, lifetime, rng):
        # power-law DTD: dN/dt ~ t^-dtd_power, sampled via inverse CDF, then
        # added on top of the progenitor's own lifetime
        u = rng.random(mass.shape)
        tmin, tmax = self.dtd_min_delay_gyr, self.dtd_max_delay_gyr
        n = 1.0 - self.dtd_power
        if np.isclose(n, 0.0):
            extra = tmin * (tmax / tmin) ** u
        else:
            extra = (u * (tmax**n - tmin**n) + tmin**n) ** (1.0 / n)
        return lifetime + extra

    def yield_table(self):
        return self._yield_table


def default_channels():
    """SNII + AGB + SNIa with default models/placeholder SNIa table. Add
    PISN yourself for Population III / extremely metal-poor runs, e.g.
    default_channels() + [PISN()]."""
    return [SNII(), AGB()]