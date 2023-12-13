import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "fsi"

# Fields to test
fields = ["Displacement", "Pressure", "Velocity"]


@pytest.mark.parametrize("name_inp", ["svFSI.xml", "svFSI_petsc.xml"])
def test_pipe_3d(name_inp, n_proc):
    test_folder = "pipe_3d"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max, name_inp=name_inp)
