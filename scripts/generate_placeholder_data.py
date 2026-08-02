"""Generate small placeholder yield tables so the pipeline runs end to end.

These numbers are NOT real nucleosynthesis yields -- they are simple
synthetic placeholders so you can see the full pipeline work before
plugging in your Fortran code's actual tables (e.g. converted from its
existing yield files). This script also doubles as a template for writing
your own converter: read the CSV format expected by
crimsons.yields.io.load_yield_table and reuse `write_table` as-is.

Run with: python scripts/generate_placeholder_data.py
"""

import csv
from pathlib import Path

import numpy as np

ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
]

OUT_DIR = Path(__file__).resolve().parents[0] / "src" / "crimsons" / "yields" / "data"


def write_table(path, masses, metallicities, yield_fn):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["mass", "metallicity", *ELEMENTS])
        for m in masses:
            for z in metallicities:
                row = yield_fn(m, z)
                writer.writerow([m, z, *[f"{v:.6e}" for v in row]])


def snii_yield(mass, z):
    idx = np.arange(1, len(ELEMENTS) + 1, dtype=float)
    return 1e-3 * mass * idx * (1.0 + 5.0 * z)


def agb_yield(mass, z):
    idx = np.arange(1, len(ELEMENTS) + 1, dtype=float)
    return 2e-4 * mass * idx * (1.0 + 5.0 * z)


def snia_yield(mass, z):
    # roughly flat with mass/metallicity, with a bump on Fe to be vaguely
    # reminiscent of a real deflagration model (~0.6 Msun of Fe) -- still
    # just a placeholder, not a citation
    idx = np.arange(1, len(ELEMENTS) + 1, dtype=float)
    row = 5e-4 * idx
    row[ELEMENTS.index("Fe")] = 0.6
    return row


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_table(OUT_DIR / "snii_placeholder.csv", [8, 15, 25, 40], [0.0, 0.02], snii_yield)
    write_table(OUT_DIR / "agb_placeholder.csv", [0.8, 2, 4, 6, 8], [0.0, 0.02], agb_yield)
    write_table(OUT_DIR / "snia_placeholder.csv", [3, 8], [0.0, 0.02], snia_yield)
    print("wrote placeholder yield tables to", OUT_DIR)