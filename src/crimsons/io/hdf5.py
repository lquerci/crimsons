from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import h5py
import numpy as np


def save_result(result, path: Path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.attrs["config_json"] = json.dumps(asdict(result.config))
        f.create_dataset("time", data=result.time)
        f.create_dataset("elements", data=np.array(result.elements, dtype="S"))
        f.create_dataset("enrichment", data=result.enrichment, compression="gzip")

        real_group = f.create_group("realizations")
        for i, (m, c, fate) in enumerate(zip(result.masses, result.counts, result.fates)):
            g = real_group.create_group(str(i))
            g.create_dataset("masses", data=m, compression="gzip")
            g.create_dataset("counts", data=c, compression="gzip")
            g.create_dataset("fates", data=np.asarray(fate, dtype="S"), compression="gzip")


def load_result(cls, path: Path):
    from ..config import RunConfig

    with h5py.File(path, "r") as f:
        cfg_dict = json.loads(f.attrs["config_json"])
        config = RunConfig(**cfg_dict)
        time = f["time"][:]
        elements = [e.decode() for e in f["elements"][:]]
        enrichment = f["enrichment"][:]

        n_real = len(f["realizations"])
        masses, counts, fates = [], [], []
        for i in range(n_real):
            g = f["realizations"][str(i)]
            masses.append(g["masses"][:])
            counts.append(g["counts"][:])
            fates.append(np.array([s.decode() for s in g["fates"][:]], dtype=object))

    return cls(
        config=config,
        time=time,
        elements=elements,
        enrichment=enrichment,
        fates=fates,
        masses=masses,
        counts=counts,
    )