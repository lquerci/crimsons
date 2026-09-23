---
title: CRIMSONS
description: Chemical Enrichment & IMF Stochastic Engine
---

# CRIMSONS

Welcome to the CRIMSONS documentation! 


!!! warning "Current status: under development"
    This Python package is currently under development and things might be unstable and untested. A more stable version of CRIMSONS can be found as web-based version of the [CRIMSONS tool](https://martina-rossi.it/crimsons.html){: target="_blank" }. Thanks for your patience. 


CRIMSONS is a tool for modelling the chemical enrichment and energy feedback from a single stellar population of Population III or Population II/I stars. It performs a stochastic sampling of the Initial Mass Function (IMF) and evolves the individual stars with four chemical enrichment channels: Supernovae type II (SNII) and Ia (SNIa), Asymptotic Giant Branch (AGB), and pair instability supernovae (PISN). The tool is distributed in two forms: an online tool version  that can be found at [CRIMSONS tool](https://martina-rossi.it/crimsons.html){: target="_blank" } and a Python package. This is the documentation page for the Python package.



## Before starting

CRIMSONS enables a high level of customization of the chemical enrichment, allowing modifications in the IMF form and shape, the enrichment channels, the stellar lifetime prescription, and the return yields. At the same time, it's also possible to run CRIMSONS with as few as three parameters. In the latter case, the code adopts a fiducial set of parameters based on the metallicity, described in Rossi et al. 2026. See [Customize](customize/index.md) for a description of the fiducial parameters and every option available to change them, and [Examples](examples/basic-simulation.md) for worked customizations.

Head to [Getting Started](getting-started/installation.md) for the full
walkthrough, or straight to [Examples](examples/basic-simulation.md) for more
worked scenarios.
 
## Where to go next

<div class="grid cards" markdown>

- :material-download:{ .lg .middle } **Getting Started**

    ---

    Install CRIMSONS and run your first simulation in a few minutes.

    [:octicons-arrow-right-24: Installation](getting-started/installation.md)
    [:octicons-arrow-right-24: Quick start](getting-started/quickstart.md)

- :material-flask-outline:{ .lg .middle } **Examples**

    ---

    Worked examples: custom IMF shapes, yield-model selection, reading
    results, caching.

    [:octicons-arrow-right-24: Basic simulation](examples/basic-simulation.md)
    [:octicons-arrow-right-24: Custom IMF](examples/custom-imf.md)
    [:octicons-arrow-right-24: Yield models](examples/yield-models.md)

- :material-atom:{ .lg .middle } **Customize**

    ---

    The IMF, lifetime, channel, and metallicity models behind the code --
    and every parameter you can change on each.

    [:octicons-arrow-right-24: Customize overview](customize/index.md)

- :material-api:{ .lg .middle } **API Reference**

    ---

    Full, auto-generated reference for every public class and function.

    [:octicons-arrow-right-24: API overview](api/index.md)

</div>

## Citing

If CRIMSONS contributes to a publication, please cite  the literature sources for whichever yield tables and IMF/lifetime prescriptions you actually used and  consider citing the 
presentation paper 
> Rossi, M., et al. (2026). *CRIMSONS: An Online Tool for Modeling Chemical Enrichment with Stochastic IMF Sampling*. arXiv preprint [arXiv:2609.25200](https://arxiv.org/abs/2609.25200).

**BibTeX:**
```bibtex
@misc{rossi2026crimsonsonlinetoolmodeling,
      title={CRIMSONS: An Online Tool for Modeling Chemical Enrichment with Stochastic IMF Sampling}, 
      author={Martina Rossi and Lapo Querci and Stefano Ciabattini and Stefania Salvadori and Irene Vanni and Davide Massari and Edoardo Ceccarelli and Donatella Romano and Viola Gelli and Raffaele Pascale and Alice Mori and Elka Rusta and Ioanna Koutsouridou and Ása Skúladóttir and Laura Magrini and Riano Giribaldi and Jose Schiappacasse-Ulloa},
      year={2026},
      eprint={2609.25200},
      archivePrefix={arXiv},
      primaryClass={astro-ph.GA},
      url={[https://arxiv.org/abs/2609.25200](https://arxiv.org/abs/2609.25200)}, 
}
```