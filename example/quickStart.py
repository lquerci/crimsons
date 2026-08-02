"""Minimal end-to-end example with IMF sampling and Enrichment visualization.

Run with: python examples/quickstart.py
"""
import matplotlib.pyplot as plt
import numpy as np

from crimsons import Simulation, StellarLifetime, default_channels
from crimsons.imf.functional import FunctionalIMF
from crimsons.imf.standard import *
from crimsons.yields.channels import *

def lognormal_func(m):
    mc = 0.079  # characteristic mass
    sigma = 0.69
    return (1.0 / m) * np.exp(-((np.log10(m) - np.log10(mc)) ** 2) / (2 * sigma**2))

def chabrier(m, m_ch = None):
    m_ch = m_ch if m_ch is not None else 10  # characteristic mass
    alpha = 1.35
    return m ** -(1+alpha)* np.exp(- m_ch / (m))


# 1. Setup and Run Simulation
sim = Simulation(
    #imf=FunctionalIMF(chabrier, imf_params={"m_ch":100 }),
    imf=Salpeter1955(),
    mass_formed=1e4,  # Msun formed in this stellar population
    metallicity=-5,
    n_realizations=50,
    verbose= True,
    time_grid=np.logspace(1,4,50),
    channels=[AGB()]
)


result = sim.run()  # Reruns with same config load from cache

mean_enrichment = result.mean()  # Shape: (n_time, n_elements)
std_enrichment = result.std()    # Shape: (n_time, n_elements)

# ==============================================================================
# 2. Plotting IMF Sampling & Chemical Enrichment
# ==============================================================================
fig, (ax_imf, ax_chem) = plt.subplots(1, 2, figsize=(14, 5.5))

# ------------------------------------------------------------------------------
# PANEL 1: IMF vs. Sampled Masses (Realization 0)
# ------------------------------------------------------------------------------
n_plot_bins = 30

# 1. Define global macro bin edges across the entire IMF range
macro_edges = np.logspace(
    np.log10(sim.imf.m_min), np.log10(sim.imf.m_max), n_plot_bins + 1
)
macro_centers = np.sqrt(macro_edges[:-1] * macro_edges[1:])
macro_dm = np.diff(macro_edges)

m_grid = np.logspace(
    np.log10(sim.imf.m_min), np.log10(sim.imf.m_max), 500
)
pdf_theoretical = sim.imf.pdf(m_grid)
ax_imf.plot(m_grid, pdf_theoretical, label="Theoretical IMF", color="black", lw=2, zorder=4)

for i in range(sim.n_realizations):
    fine_centers = result.masses[i]
    fine_counts = result.counts[i]
    total_stars = fine_counts.sum()

    if total_stars == 0:
        continue

    # Re-aggregate fine counts into the 30 macro-bins using weights
    macro_counts, _ = np.histogram(
        fine_centers, bins=macro_edges, weights=fine_counts
    )

    # Convert counts to probability density: N_macro / (N_total * dm_macro)
    pdf_macro = macro_counts.astype(float) / (total_stars * macro_dm)

    # Mask zero-count bins with NaN so log-scale doesn't attempt to plot log(0)
    #pdf_macro[macro_counts == 0] = np.nan

    ax_imf.step(
        macro_centers,
        pdf_macro,
        where="mid",
        alpha=0.7,
        linewidth=1.8,
        label=f"Realization {i+1}"
        if sim.n_realizations <= 5
        else ("Realizations" if i == 0 else None),
    )

    
ax_imf.set_xscale("log")
ax_imf.set_yscale("log")
ax_imf.set_ylim(1e-6,10)
ax_imf.set_xlabel(r"Stellar Mass [$M_\odot$]", fontsize=12)
ax_imf.set_ylabel(r"PDF $dN/dM$ [$\mathrm{M}_\odot^{-1}$]", fontsize=12)
ax_imf.set_title("IMF vs. Sampled Stellar Population", fontsize=13)
ax_imf.grid(True, which="both", linestyle="--", alpha=0.5)
ax_imf.legend(fontsize=10)

# ------------------------------------------------------------------------------
# PANEL 2: Time-Resolved Chemical Enrichment
# ------------------------------------------------------------------------------
# Choose a subset of key elements to display
elements_to_plot = ["Fe", "O", "C"]
colors = {"Fe": "crimson", "O": "dodgerblue", "C": "forestgreen"}

for elem in elements_to_plot:
    if elem in result.elements:
        idx = result.elements.index(elem)
        mean_curve = mean_enrichment[:, idx]
        std_curve = std_enrichment[:, idx]

        # Plot mean line
        ax_chem.plot(
            result.time,
            mean_curve,
            label=f"{elem} (Mean)",
            color=colors.get(elem, None),
            lw=2,
        )
        # Plot 1-sigma spread across realizations
        ax_chem.fill_between(
            result.time,
            np.maximum(0, mean_curve - std_curve),
            mean_curve + std_curve,
            color=colors.get(elem, None),
            alpha=0.2,
        )

ax_chem.set_xlabel("Time [Myr]", fontsize=12)
ax_chem.set_ylabel(r"Cumulative Returned Mass [$M_\odot$]", fontsize=12)
ax_chem.set_title("Chemical Enrichment Over Time (20 Realizations)", fontsize=13)
ax_chem.grid(True, linestyle="--", alpha=0.5)
ax_chem.legend(fontsize=10)
ax_chem.set_yscale('log')
ax_chem.set_xscale('log')

plt.tight_layout()
plt.show()