# Contributing

!!! info "Template"
    This page is a reasonable starting point, not a description of an
    existing process -- adjust the commands below to match your actual
    test runner, linter, and branching model before publishing.

## Development setup

```bash
git clone https://github.com/your-org/crimsons.git
cd crimsons
pip install -e ".[dev]"   # or: pip install -e . && pip install pytest ruff
```

## Running the test suite

```bash
pytest
```

## Docstring conventions

Docstrings in `src/crimsons` follow **NumPy style** (`Parameters`,
`Returns`, `Raises`, `Examples`, ... section headers) -- `mkdocs.yml` is
configured with `docstring_style: numpy` to match. Keep new docstrings
consistent with this, and prefer the plural `Examples` header over
`Example`, since that's what the NumPy docstring convention (and
mkdocstrings' parser) actually recognizes as a section.

## Working on the documentation

```bash
pip install -e .
pip install -r docs/requirements.txt
mkdocs serve
```

This live-reloads at `http://127.0.0.1:8000`. `mkdocs build --strict`
(what CI / Read the Docs run) will fail on broken internal links or
missing pages -- run it before opening a PR that touches `docs/`.

When adding a new public class or function, add (or extend) a `:::`
mkdocstrings block for it under `docs/api/` so it shows up in the
[API Reference](api/index.md) -- nothing there is discovered
automatically from `nav:` alone.

## Pull requests

1. Fork and branch from `main`.
2. Keep changes focused; add/update tests and docs alongside code changes.
3. Make sure `pytest` and `mkdocs build --strict` both pass locally.
4. Open a PR describing what changed and why.
