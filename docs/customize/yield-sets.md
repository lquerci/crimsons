# Yield Sets

This page lists every nucleosynthesis yield set included in CRIMSONS. The sets are grouped by enrichment channel, reporting any assumptions or modification to the published version.


## Summary
 
| Channel | Model name | Reference | Mass range (M<sub>⊙</sub>) | Metallicity (Z, absolute) | Explosion energy ($10^{51}$ erg) |
|---|---|---|---|---|---|
| AGB/SNII/PISN | `Nomoto13`[^nomoto] | [Nomoto et al. (2013)](https://scixplorer.org/abs/2013ARA&A..51..457N/abstract){: target="_blank" } | 1-300[^nomoto] | $10^{-7} - 0.05$  | 1 |
| SNII/PISN (hypernovae) | `Nomoto13_Hypernovae`[^nomoto] | [Nomoto et al. (2013)](https://scixplorer.org/abs/2013ARA&A..51..457N/abstract){: target="_blank" } | 20-140[^nomoto] | $10^{-7} - 0.05$ | $10 - 60$[^nomoto]  |
| SNII | `Heger10` | [Heger & Woosley (2010)](https://scixplorer.org/abs/2010ApJ...724..341H/abstract){: target="_blank" } | 10–100 | $10^{-7}$ | 0.3-10 |
| SNII | `Woosley95` | [Woosley & Weaver (1995)](https://scixplorer.org/abs/1995ApJS..101..181W/abstract){: target="_blank" } | 12-40 | $10^{-7} - 0.0142$ | 1 |
| SNII| `Iwamoto05` | [Iwamoto et al,2005](https://scixplorer.org/abs/2005Sci...309..451I/abstract){: target="_blank" } | 25 | $10^{-7}$ | 0.79
| PISN | `Heger02` | [Heger & Woosley (2002)](https://scixplorer.org/abs/2002ApJ...567..532H/abstract){: target="_blank" } | 140-260 | $10^{-7}$ | $5-87$ |
| SNIa | `Iwamoto99` | [Iwamoto et al. (1999)](https://scixplorer.org/abs/1999ApJS..125..439I/abstract){: target="_blank" } | 1.2 | $10^{-7} - 0.0142$ | 1.2 |
| AGB | `VanDenHoek97` | [van den Hoek & Groenewegen (1997)](https://scixplorer.org/abs/1997A&AS..123..305V/abstract){: target="_blank" } | 0.8–8 | $10^{-7} - 0.0142$ | – |
| AGB | `Meynet02` | [Meynet & Maeder (2002)](https://scixplorer.org/abs/2002A&A...390..561M/abstract){: target="_blank" } | 2 - 8 | $10^{-7}$ | – |


Tables provided in solar metallicity have been converted assuming a solar metallicity of 0.0142 (Asplund et al,2009). In the following all the yiled sets are discussed divided in the enrichment channels they contribute to. Tables with metallicity $10^-7$ relate to metal-free stars. 
 
## Core-collapse supernovae (SNII)
 
### Nomoto13
 
[Nomoto, Kobayashi & Tominaga (2013)](https://scixplorer.org/abs/2013ARA&A..51..457N/abstract) (`Nomoto13`) provide the yields of core-collapse supernovae. The mass grid differs between the PopIII and PopII case resulting in two different model names `Nomoto13_PopIII` and `Nomoto13_PopIII`, respectively. The two different mass grids are 11-100 M<sub>⊙</sub> in the PopIII case and 13-40 M<sub>⊙</sub>  in the PopII case.

### Nomoto13_Hypernovae
 
[Nomoto, Kobayashi & Tominaga (2013)](https://scixplorer.org/abs/2013ARA&A..51..457N/abstract) (`Nomoto13`) provide also the yields of hypernovae. Similar to the core-collapse, the mass grid differs between the PopIII and PopII case resulting in two different model names `Nomoto13_Hypernovae_PopIII` and `Nomoto13_Hypernovae_PopIII`, respectively. The two different mass grids are 20-100 M<sub>⊙</sub> in the PopIII case and 20-40 M<sub>⊙</sub>  in the PopII case. Explosion energies increase with the stellar mass with stars of 20 (40) M<sub>⊙</sub> releaseing $10 (30) \times 10^{51}$ erg, while stars of 100  M<sub>⊙</sub> release $60 \times 10^{51}$ erg. 
 
### Heger10
 
[Heger & Woosley (2010)](https://scixplorer.org/abs/2010ApJ...724..341H/abstract) (`Heger10`) provide the yields of non-rotating, metal-free stars of 10–100 M<sub>⊙</sub>, evolved through core collapse and exploded with a piston at the base of the oxygen shell. They provide diffent sets of explosion energies and mixing, resulting in two additional axis that can be selected (see [extra axis](../examples/yield-models.md#extra-axis)). Explosion enegies span the range 0.3-10 foe, while the mixing parameter is multiplied by 100, resulting in a range spanned of 0-251. Because the grid is metal-free it applies only to PopIII stars.
 
### Woosley95
 
[Woosley & Weaver (1995)](https://scixplorer.org/abs/1995ApJS..101..181W/abstract) (`Woosley95`) provide the yields of Type II supernovae for progenitors on a grid of 10 masses and five metallicities. The adoped model is the A model from the paper and we consider isotopes decays in the yiled table.

### Iwamoto05
 
[Iwamoto et al,2005)](https://scixplorer.org/abs/2005Sci...309..451I/abstract) (`Iwamoto05`) provide the abundance patter between Carbon and Zinc of the enrichment for a metal-free star of 25 <sub>⊙</sub>. The explosion energy of such star is $0.79 \times 10^{51}$ erg.
 
## Pair-instability supernovae (PISN)
 
### Heger02
 
[Heger & Woosley (2002)](https://scixplorer.org/abs/2002ApJ...567..532H/abstract) (`PISN/Heger02`) provide the yields of metal-free pair-instability supernovae. The explosion energy increases with stellar mass and the entire star is disrupted, leaving no remnant.


### Nomoto13
 
[Nomoto, Kobayashi & Tominaga (2013)](https://scixplorer.org/abs/2013ARA&A..51..457N/abstract) (`Nomoto13`) provide also the yields of PISN. These energetic explosion are expected only for PopIII stars and, for consistency with the core-collapse case, the name is `Nomoto13_PopIII`. The mass explosion energy increases with stellar mass from $10^{51}$ erg for the 140 M<sub>⊙</sub> star, up to $5 \times 10^{52}$ erg for the more massive star of 300 M<sub>⊙</sub>.

### Nomoto13_Hypernovae
 
[Nomoto, Kobayashi & Tominaga (2013)](https://scixplorer.org/abs/2013ARA&A..51..457N/abstract) (`Nomoto13`) provide the yields of hypernovae in the PISN regime. Specifically they provide the yields for a star of 140 M<sub>⊙</sub>  and explosion energy $7 \times 10^{52}$ erg. Similar to the core-collapse, the model name is  `Nomoto13_Hypernovae_PopIII`. 
 
## Type Ia supernovae (SNIa)
 
### Iwamoto99
 
[Iwamoto et al. (1999)](https://scixplorer.org/abs/1999ApJS..125..439I/abstract) (`SNIa/Iwamoto99`) provide the yields of Chandrasekhar-mass white-dwarf explosions, including the W7 deflagration model, the delayed-detonation models (WDD1–3, CDD1–2) and their zero-metallicity variants (W70, WDD…). Seven models are included in the tables and can be selected as extra axis (see [extra axis](../examples/yield-models.md#extra-axis)).
 
## AGB stars
 
### VanDenHoek97
 
[van den Hoek & Groenewegen (1997)](https://scixplorer.org/abs/1997A&AS..123..305V/abstract) (`AGB/VanDenHoek97`) provide yields from synthetic AGB evolution models for initial masses of 0.8–8 M<sub>⊙</sub> for both PopIII and PopII stars. They include mass loss, third dredge-up and hot-bottom burning in the most massive stars, and the explosion energy is not applicable. The elements included in the enrichment are  C, N, and O.
 
### Meynet02
 
[Meynet & Maeder (2002)](https://scixplorer.org/abs/2002A&A...390..561M/abstract) (`AGB/Meynet02`) provide the wind yields of rotating massive stars (300 km/s) at very low metallicity (PopIII). The chemical enrichment evolves the chemical elements up to sulfur, with most of the enrichment happening for CNO elements.

### Nomoto13
 
[Nomoto, Kobayashi & Tominaga (2013)](https://scixplorer.org/abs/2013ARA&A..51..457N/abstract) (`Nomoto13`) provide the for AGB stars up to sulfur. The mass range is different between PopIII and PopII stars, and therefore we adopted the same distintion as the core-collapse case, naming the two models `Nomoto13_PopIII` and `Nomoto13_PopII` respectivelu. In the former the mass range is 0.8-3 M<sub>⊙</sub>, in the latter is 1-6 M<sub>⊙</sub>.


[^nomoto]: might change with metallicity, check the entry in the list