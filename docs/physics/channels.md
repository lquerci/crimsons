# Yield Channels & Physical Assumptions

`crimsons` models nucleosynthetic enrichment by partitioning stellar mass ranges into channels:

$$
M_{\text{channel}} = 
\begin{cases} 
\text{AGB} & M < 8.1 \, \mathrm{M_\odot} \\
\text{SNII} & 8.1 \, \mathrm{M_\odot} \le M \le 120 \, \mathrm{M_\odot} \\
\text{PISN} & M > 120 \, \mathrm{M_\odot}
\end{cases}
$$

## Core Assumptions

1. **Pop III Metallicity Thresholding**:
   For metallicities $Z \le Z_{\text{PopIII}}$ (default $10^{-4}$), tables do not interpolate across $Z$. Instead, they map directly to zero-metallicity Pop III yields:
   $$ E_i(M, Z) = E_i(M, Z_{\text{PopIII}}) \quad \text{for } Z \le Z_{\text{PopIII}} $$

2. **Stochastic Sampling**:
   Individual stellar masses $m_j$ are drawn from the IMF inverse cumulative distribution function $F^{-1}(u)$ until the total target mass is reached:
   $$ M_{\text{formed}} = \sum_{j=1}^{N} m_j $$