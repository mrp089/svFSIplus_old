import os
import pytest

import pandas as pd

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = os.path.join("cases", "struct")


def test_LV_Guccione_passive(n_proc):
    folder = os.path.join(base_folder, "LV_Guccione_passive")
    fields = ["Displacement", "Velocity", "Jacobian"]
    run_with_reference(folder, fields, n_proc)


@pytest.mark.parametrize(
    "material",
    ["nHK"],
)
def test_block_compression(n_proc, material):
    folder = os.path.join(base_folder, "block_compression")
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
    t_max = 1
    name_inp = "svFSI_" + material + ".xml"
    name_ref = "result_" + material + "_" + str(t_max).zfill(3) + ".vtu"
    run_with_reference(folder, fields, n_proc, t_max, name_ref, name_inp)

def test_gr_equilibrated(n_proc):
    folder = os.path.join(base_folder, "gr_equilibrated")
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
    t_max = 11
    run_with_reference(folder, fields, n_proc, t_max)
