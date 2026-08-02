"""Default stellar mass sampling range, adaptive to metallicity.

Deliberately separate from any specific IMF shape: it answers "what mass
range should I sample over by default", which is a property of the
stellar population being formed (how metal-poor it is), not of the IMF's
functional form. Every IMF in this package uses this by default when
m_min/m_max aren't given explicitly.

Override the policy itself by passing `mass_range_fn=` to any IMF
constructor (must have the same signature as `default_mass_range`), or
bypass it entirely by passing m_min/m_max directly -- explicit bounds
always win over metallicity, at every level including inside a
Simulation.
"""

from ..chemistry import POPIII_THRESHOLD

STANDARD_MASS_RANGE = (0.1, 100.0)
EXTREMELY_METAL_POOR_MASS_RANGE = (0.8, 1000.0)


def default_mass_range(metallicity):
    """(m_min, m_max) in Msun for a given metallicity.

    Below PopIII_THRESHOLD, the IMF is often assumed to
    be more top-heavy (e.g. Population-III-like), so the default range
    is extended and raised at the low-mass end. `metallicity=None` (its
    value isn't known/relevant yet) falls back to the standard range.
    """
    if metallicity is not None and metallicity < POPIII_THRESHOLD:
        return EXTREMELY_METAL_POOR_MASS_RANGE
    return STANDARD_MASS_RANGE