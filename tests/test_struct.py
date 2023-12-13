import os
import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "struct"

# Fields to test
fields = [
    "Displacement",
    "Velocity",
    "Jacobian",
    "Stress",
    "Strain",
    "Caucy_stress",
    "Def_grad",
    "VonMises_stress",
]

# Fields to test
fields = [
    "Displacement",
    "Velocity",
    "Jacobian",
    "Stress",
    "Strain",
    "Caucy_stress",
    "Def_grad",
    "VonMises_stress",
]

# Fields to test
fields = [
    "Displacement",
    "Velocity",
    "Jacobian",
    "Stress",
    "Strain",
    "Caucy_stress",
    "Def_grad",
    "VonMises_stress",
]


def test_LV_Guccione_passive(n_proc):
    test_folder = "LV_Guccione_passive"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_block_compression(n_proc):
    test_folder = "block_compression"
    run_with_reference(base_folder, test_folder, fields, n_proc)


def test_robin(n_proc):
    test_folder = "robin"
    run_with_reference(base_folder, test_folder, fields, n_proc)


@pytest.mark.parametrize("n_proc", [1])
def test_gr_equilibrated(n_proc):
    folder = "gr_equilibrated"
    t_max = 11
    run_with_reference(base_folder, folder, fields, n_proc, t_max)
