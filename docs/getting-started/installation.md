# Installation

## Requirements

CRIMSONS targets Python 3.10+ and depends on:

| Package | Used for |
|---|---|
| `numpy` | array math throughout |
| `scipy` | grid interpolation of yield tables (`RegularGridInterpolator`) |
| `h5py` | reading the bundled yield tables and reading/writing cached results |
| `tqdm` | progress bar when running `Setup(..., verbose=True)` |

## Install from PyPI

You can install CRIMSONS via `pip`

```bash
pip install crimsons
```

!!! note
    The project hasn't been published to PyPI yet-- update this page once it has.

## Install from source

It is also possible to install from source

```bash
git clone https://github.com/lquerci/crimsons.git
cd crimsons
pip install .
```

## Verify the install

to verfy the installation run 
```bash
python -c "import crimsons; print(crimsons.__version__)"
```

The current version is 0.1.0

Next: [Quickstart](quickstart.md).
