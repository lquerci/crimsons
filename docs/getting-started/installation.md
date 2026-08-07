# Installation

## Requirements

CRIMSONS targets Python 3.10+ and depends on:

| Package | Used for |
|---|---|
| `numpy` | array math throughout |
| `scipy` | grid interpolation of yield tables (`RegularGridInterpolator`) |
| `h5py` | reading the bundled yield tables and reading/writing cached results |

Optional:

| Package | Used for |
|---|---|
| `tqdm` | progress bar when running `Simulation(..., verbose=True)` |

## Install from PyPI

```bash
pip install crimsons
```

!!! note
    If the project hasn't been published to PyPI yet, install from source
    instead (below) -- update this page once it has.

## Install from source

```bash
git clone https://github.com/your-org/crimsons.git
cd crimsons
pip install -e .
```

The `-e` (editable) install is convenient during development: changes to
`src/crimsons` take effect immediately without reinstalling.

## Verify the install

```bash
python -c "import crimsons; print(crimsons.__version__)"
```

## Building the documentation locally

This site is built with [MkDocs](https://www.mkdocs.org/) and
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/), with
[mkdocstrings](https://mkdocstrings.github.io/) generating the
[API Reference](../api/index.md) directly from docstrings in `src/crimsons`.

```bash
pip install -e .                      # crimsons itself, so mkdocstrings can import it
pip install -r docs/requirements.txt  # the doc toolchain
mkdocs serve                          # http://127.0.0.1:8000, live-reloading
```

To build a static site (what Read the Docs / CI produce):

```bash
mkdocs build --strict
```

Next: [Quickstart](quickstart.md).
