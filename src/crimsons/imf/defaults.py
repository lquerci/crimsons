from __future__ import annotations

from ..chemistry import POPIII_THRESHOLD
from .standard import Chabrier2003

# Characteristic (turnover) mass used by default_imf, Msun -- Population
# III star formation is commonly modeled as much more top-heavy than
# solar-neighborhood star formation, hence the much larger value below
# POPIII_THRESHOLD.
POPIII_M_CH = 10.0
STANDARD_M_CH = 0.35


def default_imf(metallicity, **kwargs) -> Chabrier2003:
    """A population-appropriate default IMF: `Chabrier2003`, with a much
    larger characteristic mass (`m_ch=10`) below
    `crimsons.chemistry.POPIII_THRESHOLD` than above it (`m_ch=0.35`).

    Used automatically by `Simulation` when `imf=` isn't given -- call
    this directly if you want the same metallicity-appropriate choice
    standalone. Any keyword `Chabrier2003` accepts (`m_min`, `m_max`,
    `alpha`, ...) can be overridden via `kwargs`; passing `m_ch=`
    yourself just skips the population-based choice entirely.

    Parameters
    ----------
    metallicity : float or None
        Absolute metallicity Z. `None` is treated as "unknown", which
        falls back to the standard (non-Population-III) `m_ch`.

    Returns
    -------
    Chabrier2003
    """
    is_popiii = metallicity is not None and metallicity < POPIII_THRESHOLD
    kwargs.setdefault("m_ch", POPIII_M_CH if is_popiii else STANDARD_M_CH)
    return Chabrier2003(metallicity=metallicity, **kwargs)