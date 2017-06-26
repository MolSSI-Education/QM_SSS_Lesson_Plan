"""
Unit and regression test for the MP2 classes.
"""

# Import our module and shorten the name

import psi4
import quantum_python as qp
import pytest

def test_mp2_water():
    mol = psi4.geometry("""
    O
    H 1 1.1
    H 1 1.1 2 104.5
    """)

    rhf_object = qp.RHF(mol, "sto-3g")
    rhf_energy = rhf_object.compute_energy()
    assert rhf_energy == pytest.approx(-74.9418022615317909, 1.e-5)

    mp2_object = qp.MP2(rhf_object)
    mp2_energy = mp2_object.compute_energy()
    assert mp2_energy == pytest.approx(-74.9908766157, 1.e-5)
