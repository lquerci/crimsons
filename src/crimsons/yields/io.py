from __future__ import annotations

import csv
from pathlib import Path

import h5py
import numpy as np

from .base import StochasticYieldTable, YieldTable, _is_distribution


def load_yield_table(csv_path) -> YieldTable:
    """Load a yield table from a CSV with columns: mass, metallicity, then
    one column per element (see yields/data/snia_placeholder.csv for the
    expected shape). SNII/AGB/PISN load from stellar_yields.h5 instead
    (see load_yield_table_hdf5 below) -- this CSV path is still here for
    SN Ia and for anyone who wants a simple custom table without going
    through HDF5.
    """
    csv_path = Path(csv_path)
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        elements = header[2:]
        rows = [list(map(float, row)) for row in reader]

    rows = np.array(rows)
    masses = np.unique(rows[:, 0])
    metallicities = np.unique(rows[:, 1])

    yields = np.zeros((len(masses), len(metallicities), len(elements)))
    m_index = {m: i for i, m in enumerate(masses)}
    z_index = {z: i for i, z in enumerate(metallicities)}
    for row in rows:
        mi = m_index[row[0]]
        zi = z_index[row[1]]
        yields[mi, zi, :] = row[2:]

    return YieldTable(masses=masses, metallicities=metallicities, elements=elements, yields=yields)


# --- stellar_yields.h5: multi-model, multi-parameter yield tables ---------
#
# Schema (see also scripts/generate_stellar_yields_h5.py, which builds a
# synthetic file matching this exactly):
#
#   /<channel>/<model>/          e.g. /SNII/LC, /AGB/VAN, /PISN/HW
#       yields          N-D dataset; axes as listed in attrs["axes"]
#       metallicity     1-D dataset (always present)
#       mass            1-D dataset (always present)
#       elements        1-D string dataset (always present, one axis)
#       <extra axes>    1-D dataset per extra physical parameter, e.g.
#                       "rotation" (LC), "energy"/"mixing" (HW)
#       attrs["axes"]   axis order of `yields`, e.g.
#                       ["metallicity", "mass", "rotation", "elements"]
#
# Every model must have "metallicity", "mass", and "elements" axes. Any
# other axis is a discrete model parameter that gets resolved to a single
# slice at load time (nearest grid value, not interpolated) via
# model_params -- see load_yield_table_hdf5.


def list_models(h5_path, channel: str) -> list:
    """List the model names available for `channel` (e.g. 'SNII') in a
    stellar-yields HDF5 file, by scanning /<channel>/* groups."""
    h5_path = Path(h5_path)
    with h5py.File(h5_path, "r") as f:
        if channel not in f:
            raise KeyError(f"no channel group {channel!r} in {h5_path.name} (available: {list(f.keys())})")
        return sorted(f[channel].keys())


def describe_model(h5_path, channel: str, model: str) -> dict:
    """Introspect one model's grid: its axis order (from the group's
    "axes" attribute) and, for every axis, its available values -- useful
    for finding out what model_params a model needs before loading it.
    """
    h5_path = Path(h5_path)
    with h5py.File(h5_path, "r") as f:
        group = _get_model_group(f, channel, model)
        axes = list(group.attrs["axes"])
        info = {"axes": axes, "yields_shape": tuple(group["yields"].shape)}
        for axis in axes:
            info[axis] = _read_axis(group, axis)
        return info


def load_yield_table_hdf5(h5_path, channel: str, model: str, model_params: dict | None = None):
    """Load one model's yield table from a stellar-yields HDF5 file,
    resolving any extra parameter axes (e.g. HW's energy/mixing, LC's
    rotation).

    Every model has "metallicity", "mass", and "elements" axes. Any
    additional axes are discrete model variants. For each one, a
    model_params entry can be either:

      - a fixed value (e.g. {"rotation": 300}) -- resolved once, the
        same for every star, snapped to the nearest tabulated grid
        point; or
      - a distribution to draw from per star (e.g.
        {"rotation": scipy.stats.norm(loc=200, scale=50)}, or any plain
        callable(rng, n) -> array of n draws) -- each star gets its own
        independent draw, again snapped to the nearest grid point. This
        returns a StochasticYieldTable in that case, which needs an rng
        at call time (Channel.yields provides one automatically).

    An axis with only one grid value is picked automatically either way
    -- no model_params entry needed. A ValueError names the axis and its
    available values if a required (multi-valued, unspecified) axis is
    missing from model_params.

    Snapping to the nearest value (rather than interpolating) applies
    here regardless of fixed-or-distribution: these are separate model
    calculations (different literature sources' physics), not a smoothly
    varying physical dimension.
    """
    model_params = dict(model_params or {})
    h5_path = Path(h5_path)

    with h5py.File(h5_path, "r") as f:
        group = _get_model_group(f, channel, model)
        axes = list(group.attrs["axes"])
        yields_nd = group["yields"][:]
        grids = {axis: _read_axis(group, axis) for axis in axes}

    for required in ("metallicity", "mass", "elements"):
        if required not in axes:
            raise ValueError(f"{channel}/{model} is missing a required '{required}' axis; has {axes}")

    extra_axes = [a for a in axes if a not in ("metallicity", "mass", "elements")]
    pdf_axes = [a for a in extra_axes if _is_distribution(model_params.get(a))]
    fixed_axes = [a for a in extra_axes if a not in pdf_axes]

    # resolve FIXED axes first: collapse them out of yields_nd immediately
    slicer = [slice(None)] * yields_nd.ndim
    chosen_fixed = {}
    for axis in fixed_axes:
        values = np.asarray(grids[axis], dtype=float)
        if len(values) == 1:
            idx = 0
        elif axis in model_params:
            idx = int(np.argmin(np.abs(values - model_params[axis])))
        else:
            raise ValueError(
                f"{channel}/{model} has {len(values)} possible '{axis}' values "
                f"{values.tolist()} -- pass model_params={{'{axis}': <value or distribution>}} to pick one"
            )
        slicer[axes.index(axis)] = idx
        chosen_fixed[axis] = float(values[idx])

    collapsed = yields_nd[tuple(slicer)]
    remaining_axes = [a for a in axes if a not in fixed_axes]  # metallicity, mass, [pdf axes...], elements

    def build_table(pdf_idx_by_axis: dict) -> YieldTable:
        s = [slice(None)] * collapsed.ndim
        for axis, idx in pdf_idx_by_axis.items():
            s[remaining_axes.index(axis)] = idx
        sliced = collapsed[tuple(s)]
        surviving = [a for a in remaining_axes if a not in pdf_idx_by_axis]
        perm = [surviving.index("mass"), surviving.index("metallicity"), surviving.index("elements")]
        sliced = np.transpose(sliced, axes=perm)
        return YieldTable(
            masses=np.asarray(grids["mass"], dtype=float),
            metallicities=np.asarray(grids["metallicity"], dtype=float),
            elements=list(grids["elements"]),
            yields=sliced,
        )

    if not pdf_axes:
        table = build_table({})
        table.model = model
        table.model_params = chosen_fixed
        return table

    tables_by_index = {}
    grid_lengths = [len(grids[a]) for a in pdf_axes]
    for idx_tuple in np.ndindex(*grid_lengths):
        tables_by_index[idx_tuple] = build_table(dict(zip(pdf_axes, idx_tuple)))

    stochastic = StochasticYieldTable(
        tables_by_index=tables_by_index,
        axis_names=pdf_axes,
        axis_grids={a: np.asarray(grids[a], dtype=float) for a in pdf_axes},
        axis_distributions={a: model_params[a] for a in pdf_axes},
    )
    stochastic.model = model
    stochastic.model_params = {**chosen_fixed, **{a: model_params[a] for a in pdf_axes}}
    return stochastic


def _get_model_group(f, channel, model):
    if channel not in f:
        raise KeyError(f"no channel group {channel!r} (available: {list(f.keys())})")
    if model not in f[channel]:
        raise KeyError(
            f"no model {model!r} under channel {channel!r} (available: {sorted(f[channel].keys())})"
        )
    return f[channel][model]


def _read_axis(group, axis):
    values = group[axis][:]
    if axis == "elements":
        return [v.decode() if isinstance(v, bytes) else v for v in values]
    return values.tolist()