import csv
from pathlib import Path

import h5py
import numpy as np

from .base import YieldTable


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


def load_yield_table_hdf5(h5_path, channel: str, model: str, model_params: dict | None = None) -> YieldTable:
    """Load one model's yield table from a stellar-yields HDF5 file,
    resolving any extra parameter axes (e.g. HW's energy/mixing, LC's
    rotation) to a fixed slice.

    Every model has "metallicity", "mass", and "elements" axes. Any
    additional axes are discrete model variants -- if an axis has more
    than one grid value, you must say which one to use via model_params
    (e.g. {"rotation": 300} or {"energy": 1.2, "mixing": 0.02}); axes
    with only one grid value are picked automatically, no need to specify
    them. The chosen value is snapped to the nearest available grid
    point, not interpolated across -- these are separate model
    calculations, not a smoothly varying physical dimension the engine
    has any way to assign per star.

    Raises a ValueError naming the axis and its available values if a
    required model_params entry is missing.
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
    slicer = [slice(None)] * yields_nd.ndim
    chosen = {}
    for axis in extra_axes:
        values = np.asarray(grids[axis], dtype=float)
        if len(values) == 1:
            idx = 0
        elif axis in model_params:
            idx = int(np.argmin(np.abs(values - model_params[axis])))
        else:
            raise ValueError(
                f"{channel}/{model} has {len(values)} possible '{axis}' values "
                f"{values.tolist()} -- pass model_params={{'{axis}': <value>}} to pick one"
            )
        slicer[axes.index(axis)] = idx
        chosen[axis] = float(values[idx])

    sliced = yields_nd[tuple(slicer)]  # now 3D, axes in their original relative order
    surviving_axes = [a for a in axes if a not in extra_axes]
    perm = [surviving_axes.index("mass"), surviving_axes.index("metallicity"), surviving_axes.index("elements")]
    sliced = np.transpose(sliced, axes=perm)

    table = YieldTable(
        masses=np.asarray(grids["mass"], dtype=float),
        metallicities=np.asarray(grids["metallicity"], dtype=float),
        elements=list(grids["elements"]),
        yields=sliced,
    )
    table.model = model  # lightweight provenance, not used by the interpolator
    table.model_params = chosen
    return table


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