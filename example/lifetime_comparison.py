import matplotlib

matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
import numpy as np

# Load our LifetimeFunction models defined previously
# (Assuming FortranLifetime and PowerLawLifetime are available in your code)
from crimsons import StellarLifetime

# 1. Instantiate the Lifetime models
# Standard stellar mass range (0.08 to 500 M_sun to capture Pop III regimes)
m_min, m_max = 0.1, 500.0

models = {
    "Raiteri+96 (Pop II/I, Z = Zsun)": StellarLifetime(),
    "Pop III (Raiteri+96 / Scharer+02, Z = 1e-5)": StellarLifetime(),
}

# Define corresponding metallicities [log10(Z/Z_sun)]
metallicities = {
    "Raiteri+96 (Pop II/I, Z = Zsun)": 0.0142,       # Solar metallicity
    "Pop III (Raiteri+96 / Scharer+02, Z = 1e-5)": 10**(-7.0), # Sub-critical (Pop III)
}

# 2. Generate logarithmically spaced mass points
mass_grid = np.logspace(np.log10(m_min), np.log10(m_max), 500)

# 3. Plot Stellar Lifetime vs Mass
plt.figure(figsize=(8, 6))

for label, model in models.items():
    z_val = metallicities[label]
    lifetime_gyr = model(mass_grid, metallicity=z_val)
    
    plt.plot(mass_grid, lifetime_gyr*1000, label=label, linewidth=2)

plt.scatter(list[:,0], list[:,1]/1e6)

# Set log-log scaling
plt.xscale("log")
plt.yscale("log")

# Annotate axes and title
plt.xlabel(r"Stellar Mass [$M_\odot$]", fontsize=12)
plt.ylabel(r"Stellar Lifetime $\tau$ [Myr]", fontsize=12)
plt.title("Stellar Lifetime Comparison", fontsize=14)
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend(fontsize=11)

plt.tight_layout()
plt.show()
