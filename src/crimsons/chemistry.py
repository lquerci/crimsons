"""The 30 elements tracked through the enrichment pipeline.

Default choice: the first 30 elements of the periodic table (H through
Zn). If your yield tables track a different set of 30 species (e.g. with
isotopes, or skipping the noble gases), change this list -- it is read
by the bundled yield tables and by everything downstream, so this is the
one place to edit.
"""

ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
]

assert len(ELEMENTS) == 30

ZSUN = 0.0142

POPIII_THRESHOLD = 3e-7  # roughly log(Z/Z_sun) < -4.5
