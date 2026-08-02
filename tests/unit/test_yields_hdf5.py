from importlib.resources import files

import h5py
import numpy as np
import pytest

from crimsons.chemistry import ELEMENTS, ZSUN
from crimsons.yields.channels import (
    AGB,
    PISN,
    SNII,
    describe_available_model,
    list_available_models,
)
from crimsons.yields.io import describe_model, list_models, load_yield_table_hdf5

H5_PATH = files("crimsons.yields") / "data" / "stellar_yields.h5"


def test_list_models_scans_channel_groups():
    assert set(list_models(H5_PATH, "SNII")) == {"LC", "NK", "NK_HN", "HW"}
    assert set(list_models(H5_PATH, "AGB")) == {"VAN", "MM", "NK"}
    assert set(list_models(H5_PATH, "PISN")) == {"HW", "NK", "NK_HN"}


def test_list_models_unknown_channel_raises_with_available_listed():
    with pytest.raises(KeyError, match="AGB"):
        list_models(H5_PATH, "NOT_A_CHANNEL")


def test_describe_model_reports_axes_and_grids():
    info = describe_model(H5_PATH, "SNII", "LC")
    assert info["axes"] == ["metallicity", "rotation", "mass", "elements"]
    assert info["metallicity"] == pytest.approx([1.420e-5, 1.420e-4, 1.420e-3,1.420e-2])
    assert info["rotation"] == [0.0, 150.0, 300.0]
    assert info["mass"] == [13,15,20,25,30,40,60,80,120]
    assert info["elements"] == list(ELEMENTS)
    assert info["yields_shape"] == (4, 3, 9, 30)


def test_plain_3d_model_loads_and_matches_grid_point_exactly():
    """NK has no extra parameters -- the simple case."""
    table = load_yield_table_hdf5(H5_PATH, "SNII", "NK")
    NK_mass_list = [10,11,13,15,18,20,25,30,40,100]
    NK_met_list = [1.000e-7,0.001,0.004,0.008,0.02,0.05]
    assert table.masses.tolist() == NK_mass_list
    assert table.metallicities.tolist() == NK_met_list

    # synth_yield(mass=20, z=0.0, elem_idx=0) = 1e-3*20*1*(1+0)*1.0
    y = table(mass=20.0, metallicity=0.0)
    assert y[0, 0] == pytest.approx(8.7742)


def test_model_with_one_extra_axis_requires_model_params_when_ambiguous():
    with pytest.raises(ValueError, match="rotation"):
        load_yield_table_hdf5(H5_PATH, "SNII", "LC")


def test_model_with_one_extra_axis_resolves_to_exact_slice():
    """LC has a 'rotation' axis with 3 values -- picking one should
    reproduce exactly what the generator wrote for that slice."""
    table = load_yield_table_hdf5(H5_PATH, "SNII", "LC", model_params={"rotation": 300})
    assert table.model_params == {"rotation": 300.0}

    # synth_yield(mass=25, z=0.0, elem_idx=0, factor=1+300/1000=1.3)
    y = table(mass=25.0, metallicity=ZSUN)
    assert y[0, 0] == pytest.approx(6.4478000005477645)


def test_nearest_rotation_value_is_snapped_not_interpolated():
    exact = load_yield_table_hdf5(H5_PATH, "SNII", "LC", model_params={"rotation": 150})
    nearby = load_yield_table_hdf5(H5_PATH, "SNII", "LC", model_params={"rotation": 200})
    assert nearby.model_params["rotation"] == 150.0  # 200 is closer to 150 than to 300
    np.testing.assert_allclose(nearby.yields, exact.yields)


def test_model_with_two_extra_axes_resolves_to_exact_slice():
    """HW has both 'energy' and 'mixing' axes."""
    table = load_yield_table_hdf5(
        H5_PATH, "SNII", "HW", model_params={"energy": 3, "mixing": 0.0}
    )
    assert table.model_params == {"energy": 3, "mixing": 0.0}

    c_index = ELEMENTS.index("C")
    # synth_yield(mass=20, z=0.0, elem_idx=c_index, factor=(1.2/1.2)*(1+0.01*2)=1.02)
    y = table(mass=20.0, metallicity=1e-7)
    expected = 1.869e-10
    assert y[0, c_index] == pytest.approx(expected)


def test_missing_one_of_two_required_params_names_the_axis():
    with pytest.raises(ValueError, match="mixing"):
        load_yield_table_hdf5(H5_PATH, "SNII", "HW", model_params={"energy": 1.2})


def test_single_metallicity_model_uses_1d_fallback():
    """PISN/HW has only Z=0 -- this must not crash trying to build a 2D
    interpolator, and metallicity should be irrelevant to the result."""
    table = load_yield_table_hdf5(H5_PATH, "PISN", "HW")
    assert table.metallicities.tolist() == [1e-7]

    y_z0 = table(mass=200.0, metallicity=1e-7)
    y_other_z = table(mass=200.0, metallicity=0.02)
    np.testing.assert_allclose(y_z0, y_other_z)


def test_extra_axis_with_a_single_value_is_auto_selected(tmp_path):
    """If an extra axis only has one grid value, no model_params entry
    should be required for it."""
    path = tmp_path / "mini.h5"
    with h5py.File(path, "w") as f:
        g = f.create_group("SNII/ONEVAL")
        g.create_dataset("yields", data=np.arange(2 * 2 * 1 * 3, dtype=float).reshape(2, 2, 1, 3))
        g.create_dataset("metallicity", data=[0.0, 0.02])
        g.create_dataset("mass", data=[10.0, 20.0])
        g.create_dataset("extra", data=[5.0])
        g.create_dataset("elements", data=np.array(["H", "He", "C"], dtype=h5py.string_dtype()))
        g.attrs["axes"] = ["metallicity", "mass", "extra", "elements"]

    table = load_yield_table_hdf5(path, "SNII", "ONEVAL")  # no model_params needed
    assert table.masses.tolist() == [10.0, 20.0]
    assert table.model_params == {"extra": 5.0}


def test_unknown_model_raises_with_available_listed():
    with pytest.raises(KeyError, match="NK"):
        load_yield_table_hdf5(H5_PATH, "SNII", "NOT_A_MODEL")


def test_list_and_describe_available_models_convenience_wrappers():
    assert "VAN" in list_available_models("AGB")
    info = describe_available_model("AGB", "VAN")
    assert info["axes"] == ["metallicity", "mass", "elements"]


def test_channels_construct_with_default_and_explicit_models():
    snii_default = SNII()  # default model is NK, needs no model_params
    snii_rotating = SNII(model="LC", model_params={"rotation": 0})
    agb = AGB(model="MM")
    pisn = PISN()

    for channel, mass, mass_range in [
        (snii_default, 25.0, (8.0, 40.0)),
        (snii_rotating, 25.0, (8.0, 40.0)),
        (agb, 3.0, (0.8, 8.0)),
        (pisn, 200.0, (140.0, 260.0)),
    ]:
        mask = channel.contributes(np.array([mass]), metallicity=0.0, rng=None)
        assert mask[0]
        y = channel.yields(np.array([mass]), metallicity=0.0)
        assert y.shape == (1, len(ELEMENTS))
        assert np.all(y >= 0)
        assert channel.mass_min, channel.mass_max == mass_range


def test_snii_with_no_args_does_not_require_model_params():
    """The zero-config path (as used by default_channels()) must not
    force a rotation/energy/mixing choice on the user."""
    channel = SNII()
    assert channel.model == "NK"