# CRIMSONS

Chemical evolution with the Random sampling of the Initial Mass function: Studying the Origin of Nucleosynthetic Stellar products

This python package performs the stochastic chemical-enrichment simulations of stellar populations: sample
an initial mass function, evolve it through core-collapse supernovae, AGB
winds, Type Ia supernovae and (optionally) pair-instability supernovae,
and get back the time-resolved abundances of 30 elements released to the
interstellar medium.

> **Status:** first public version

## Install

```bash
pip install crimsons 
```

or from source

```bash
git clone https://github.com/lquerci/crimsons.git
cd crimsons
pip install .
```


## Quickstart

```python
from crimsons import Kroupa2001, Simulation

sim = Simulation(
    imf=Kroupa2001(),
    mass_formed=1e6,
    metallicity=0.0142,
    n_realizations=20,
    seed=42,
)
result = sim.run()
mean_fe = result["Fe"].mean()
```

## Documentation

Full documentation -- installation, examples, the physics behind each
piece, and the complete API reference -- lives at
`https://crimsons.readthedocs.io`.

## License

MIT -- see [LICENSE](LICENSE).