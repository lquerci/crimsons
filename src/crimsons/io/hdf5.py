from __future__ import annotations

import json
import warnings
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
        f.create_dataset("energy", data=result.energy, compression="gzip")

        real_group = f.create_group("realizations")

        for i, (m, c, fate, events) in enumerate(zip(result.masses, result.counts, result.fates, result.events_history)):
            g = real_group.create_group(str(i))
            g.create_dataset("masses", data=m, compression="gzip")
            g.create_dataset("counts", data=c, compression="gzip")
            g.create_dataset("fates", data=np.asarray(fate, dtype="S"), compression="gzip")
            
            # ---  Save events history ---
            events_grp = g.create_group("events")
            for j, (name, times, counts) in enumerate(events):
                ev_g = events_grp.create_group(str(j))
                ev_g.attrs["name"] = name
                ev_g.create_dataset("times", data=times, compression="gzip")
                ev_g.create_dataset("counts", data=counts, compression="gzip")


def load_result(cls, path: Path):
    from crimsons.config import RunConfig  # Adjust import based on your actual path

    with h5py.File(path, "r") as f:
        cfg_dict = json.loads(f.attrs["config_json"])
        config = RunConfig(**cfg_dict)
        time = f["time"][:]
        elements = [e.decode() for e in f["elements"][:]]
        enrichment = f["enrichment"][:]

        if "energy" in f:
            energy = f["energy"][:]
        else:
            warnings.warn(
                f"{path} was saved before explosion-energy tracking was split "
                "out of EnrichmentResult.elements -- filling `energy` with "
                "zeros. If this file's `elements` still has 31 entries "
                "(includes 'Energy'), delete and regenerate it instead of "
                "relying on this fallback.",
                UserWarning,
                stacklevel=2,
            )
            energy = np.zeros(enrichment.shape[:2])

        n_real = len(f["realizations"])
        masses, counts, fates, events_history = [], [], [], []
        
        for i in range(n_real):
            g = f["realizations"][str(i)]
            masses.append(g["masses"][:])
            counts.append(g["counts"][:])
            fates.append(np.array([s.decode() for s in g["fates"][:]], dtype=object))
            
            # --- Load events history ---
            real_events = []
            if "events" in g:
                events_grp = g["events"]
                # Iterate in exact numerical order
                n_events = len(events_grp)
                for j in range(n_events):
                    ev_g = events_grp[str(j)]
                    name = ev_g.attrs["name"]
                    times = ev_g["times"][:]
                    counts_arr = ev_g["counts"][:]
                    real_events.append((name, times, counts_arr))
            events_history.append(real_events)

    return cls(
        config=config,
        time=time,
        elements=elements,
        enrichment=enrichment,
        energy=energy,
        fates=fates,
        masses=masses,
        counts=counts,
        events_history=events_history,  # Pass it back into the constructor
    )