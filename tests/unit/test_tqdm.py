from unittest.mock import patch

import numpy as np
import pytest

from crimsons import (
    Salpeter1955,
    Simulation,
)


def _make_simulation(**overrides):
    kwargs = {
        "imf": Salpeter1955(),
        "mass_formed": 2000.0,
        "metallicity": 0.02,
        "n_realizations": 3,
        "seed": 123,
        "time_grid": np.linspace(0, 13.8, 30),
        "verbose": False
    }
    kwargs.update(overrides)
    return Simulation(**kwargs)

def test_tqdm_missing_raises_error():
    """Test that setting verbose=True without tqdm raises an informative ImportError."""
    # Simulate tqdm not being installed
    with patch.dict("sys.modules", {"tqdm": None}):
        with pytest.raises(ImportError) as exc_info:
            _make_simulation(verbose=True).run()

        assert "pip install tqdm" in str(exc_info.value)