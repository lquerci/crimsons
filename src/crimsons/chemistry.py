"""The 30 elements tracked through the enrichment pipeline, their atomic
weights, and solar reference abundances for [X/Y] ratios.

Default choice of elements: the first 30 elements of the periodic table
(H through Zn). If your yield tables track a different set of 30 species
(e.g. with isotopes, or skipping the noble gases), change `ELEMENTS` --
it is read by the bundled yield tables and by everything downstream, so
this is the one place to edit.
"""

from __future__ import annotations

import csv
from pathlib import Path

ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
]

HDF_COLUMNS = ["Energy"] + ELEMENTS

assert len(ELEMENTS) == 30
assert len(HDF_COLUMNS) == 31

ZSUN = 0.0142

POPIII_THRESHOLD = 3e-7  # roughly log(Z/Z_sun) < -4.5


# --- Atomic weights --------------------------------------------------------

ATOMIC_WEIGHTS = {
    # Standard atomic weights (IUPAC), g/mol. Only used to convert a
    # log number-abundance solar reference (see `SolarAbundances.from_log_eps`)
    # into mass units -- yields themselves are always in Msun, never
    # converted through this table.
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
}

assert set(ATOMIC_WEIGHTS) == set(ELEMENTS)


# --- Solar reference abundances, for [X/Y] ratios ---------------------------

_DATA_DIR = Path(__file__).parent / "data"
_DEFAULT_SOLAR_ABUNDANCES_FILE = _DATA_DIR / "solar_abundances_asplund2009.csv"


class SolarAbundances:
    """A solar reference abundance table, for computing ``[X/Y]`` ratios.

    Internally this stores, per element, a value *proportional to* that
    element's solar abundance by mass -- proportional, not absolute,
    because :meth:`ratio` (and therefore every ``[X/Y]`` computed from it,
    see `EnrichmentResult.__getitem__`) only ever uses the ratio between
    two elements. A common normalization constant cancels in that ratio
    and is never tracked.

    Most users won't construct this directly -- load a table with
    :func:`load_solar_abundances`, or use the bundled default via
    :func:`default_solar_abundances`.

    Parameters
    ----------
    mass_proxy : dict[str, float]
        Mapping of element symbol to a mass-proportional quantity. Must
        cover every element in `ELEMENTS`.
    """

    def __init__(self, mass_proxy: dict):
        missing = set(ELEMENTS) - set(mass_proxy)
        if missing:
            raise ValueError(
                f"Solar abundance table is missing elements: {sorted(missing)}"
            )
        self._mass_proxy = dict(mass_proxy)

    @classmethod
    def from_log_eps(cls, log_eps: dict) -> SolarAbundances:
        """Build from log number-abundances, the standard astronomical
        convention :math:`A(X) = \\log_{10}(N_X / N_H) + 12` (e.g. Table 1
        of Asplund et al. 2009).

        Parameters
        ----------
        log_eps : dict[str, float]
            Mapping of element symbol to ``A(X)``.

        Returns
        -------
        SolarAbundances
        """
        mass_proxy = {el: 10.0**a * ATOMIC_WEIGHTS[el] for el, a in log_eps.items()}
        return cls(mass_proxy)

    @classmethod
    def from_mass_fractions(cls, mass_fractions: dict) -> SolarAbundances:
        """Build directly from a mass-proportional quantity per element --
        a true mass fraction, a mass ratio to hydrogen, or anything else
        proportional to solar mass abundance. The overall normalization
        doesn't matter, only the ratios between elements.

        Parameters
        ----------
        mass_fractions : dict[str, float]
            Mapping of element symbol to a mass-proportional abundance.

        Returns
        -------
        SolarAbundances
        """
        return cls(dict(mass_fractions))

    def mass_ratio(self, element: str) -> float:
        """A value proportional to `element`'s solar mass abundance.

        Parameters
        ----------
        element : str

        Returns
        -------
        float
        """
        try:
            return self._mass_proxy[element]
        except KeyError:
            raise KeyError(
                f"No solar abundance for {element!r}. "
                f"Available: {sorted(self._mass_proxy)}"
            ) from None

    def ratio(self, numerator: str, denominator: str) -> float:
        """The solar mass ratio ``M_numerator / M_denominator``.

        Parameters
        ----------
        numerator : str
        denominator : str

        Returns
        -------
        float
        """
        return self.mass_ratio(numerator) / self.mass_ratio(denominator)


def load_solar_abundances(path, column: str | None = None) -> SolarAbundances:
    """Load a solar reference abundance table from a CSV file.

    The file must have an ``element`` column and one abundance column.
    The abundance column is auto-detected by name (case-insensitive)
    unless `column` is given explicitly:

    - a name containing ``"mass"`` (e.g. ``mass_fraction``) is read as a
      mass-proportional quantity directly, via
      `SolarAbundances.from_mass_fractions`
    - anything else (e.g. ``A``, ``logeps``, ``photosphere``) is read as
      :math:`A(X) = \\log_{10}(N_X/N_H) + 12`, via
      `SolarAbundances.from_log_eps`

    Parameters
    ----------
    path : str or pathlib.Path
        CSV file to read.
    column : str, optional
        Force a specific column name instead of auto-detecting one.

    Returns
    -------
    SolarAbundances

    Examples
    --------
    >>> table = load_solar_abundances("my_solar_abundances.csv")
    >>> table.ratio("C", "Fe")  # doctest: +SKIP
    """
    path = Path(path)
    with open(path, newline="") as f:
        # `#`-prefixed comment/citation lines are allowed before the
        # header -- csv.DictReader has no concept of comments, so strip
        # them ourselves rather than mis-reading one as the header row.
        data_lines = (line for line in f if not line.lstrip().startswith("#"))
        rows = list(csv.DictReader(data_lines))
    if not rows:
        raise ValueError(f"{path} has no data rows")

    if column is None:
        candidates = [c for c in rows[0] if c.lower() != "element"]
        if not candidates:
            raise ValueError(f"{path} has no abundance column besides 'element'")
        mass_like = [c for c in candidates if "mass" in c.lower()]
        column = mass_like[0] if mass_like else candidates[0]

    values = {row["element"]: float(row[column]) for row in rows}
    if "mass" in column.lower():
        return SolarAbundances.from_mass_fractions(values)
    return SolarAbundances.from_log_eps(values)


_default_solar_abundances: SolarAbundances | None = None


def default_solar_abundances() -> SolarAbundances:
    """The bundled default `SolarAbundances` (Asplund et al. 2009, Table 1
    photospheric composition), loaded once and cached for the rest of the
    process.

    Used automatically by ``result['[X/Y]']`` lookups (see
    `crimsons.results.EnrichmentResult.__getitem__`) unless you've called
    `set_default_solar_abundances`.

    Returns
    -------
    SolarAbundances
    """
    global _default_solar_abundances
    if _default_solar_abundances is None:
        _default_solar_abundances = load_solar_abundances(
            _DEFAULT_SOLAR_ABUNDANCES_FILE
        )
    return _default_solar_abundances


def set_default_solar_abundances(table) -> None:
    """Replace the default `SolarAbundances` used by ``result['[X/Y]']``
    lookups for the rest of the process.

    Parameters
    ----------
    table : SolarAbundances or str or pathlib.Path
        A ready-made table, or a path to load one from (via
        `load_solar_abundances`).

    Examples
    --------
    >>> import crimsons
    >>> crimsons.set_default_solar_abundances("my_solar_abundances.csv")  # doctest: +SKIP
    """
    global _default_solar_abundances
    if not isinstance(table, SolarAbundances):
        table = load_solar_abundances(table)
    _default_solar_abundances = table