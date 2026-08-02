import matplotlib
matplotlib.use('QtAgg')
import matplotlib.pyplot as plt
import numpy as np

# Import your IMFs from the package
from crimsons.imf.standard import FlatIMF, Salpeter1955, Kroupa2001
#from crimsons.imf.numerical import NumericalIMF
from crimsons.imf.functional import FunctionalIMF

def lognormal_func(m):
    mc = 0.079  # characteristic mass
    sigma = 0.69
    return (1.0 / m) * np.exp(-((np.log10(m) - np.log10(mc)) ** 2) / (2 * sigma**2))

def chabrier(m):
    m_ch = 10  # characteristic mass
    alpha = 1.35
    return m ** -(1+alpha)* np.exp(- m_ch / (m))


# 1. Instantiate the IMFs with a standard mass range (e.g., 0.1 to 120 M_sun)
m_min, m_max = 0.8, 100.0

imfs = {
    "Salpeter (1955)": Salpeter1955(m_min=m_min, m_max=m_max),
    "Flat IMF": FlatIMF(m_min=m_min, m_max=m_max),
    "Kroupa (2001)": Kroupa2001(m_min=m_min, m_max=m_max),
    "LogNormal": FunctionalIMF(lognormal_func,m_min=m_min, m_max=m_max),
    "Chabrier": FunctionalIMF(chabrier,m_min=m_min, m_max=m_max),
    
}

# check integration
from scipy.integrate import quad

for label, imf in imfs.items(): 
    area, _ = quad(lambda m: imf.pdf(m)[0], m_min, m_max)

    print(f"{label} Area Integral: {area:.6f}") # Should print: 1.000000



# 2. Generate logarithmically spaced mass points for smooth evaluation
mass_grid = np.logspace(np.log10(m_min), np.log10(m_max), 500)

# 3. Plot dN/dM vs Mass
plt.figure(figsize=(8, 6))

for label, imf in imfs.items():
    pdf_values = imf.pdf(mass_grid)
    plt.plot(mass_grid, pdf_values, label=label, linewidth=2)

# Set log-log scaling
plt.xscale("log")
plt.yscale("log")

# Annotate axes and title
plt.xlabel(r"Stellar Mass [$M_\odot$]", fontsize=12)
plt.ylabel(r"PDF $dN/dM$ [$\mathrm{M}_\odot^{-1}$]", fontsize=12)
plt.title("Initial Mass Function (IMF) Comparison", fontsize=14)
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend(fontsize=11)

plt.tight_layout()
plt.show()