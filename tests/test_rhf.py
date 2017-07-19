"""
Unit and regression test for the numeric file.
"""

# Import our module and shorten the name
import quantum_python as qp
import psi4
import pytest


# Test python_template/text/levenstein
mol = psi4.geometry("""
O
H 1 1.1
H 1 1.1 2 104.5
""")

def test_rhf_water():
    
    rhf_object = qp.RHF(mol, "sto-3g")
    rhf_energy = rhf_object.compute_energy()

    assert rhf_energy == pytest.approx(-74.9418022615317909, 1.e-5)

def test_rhf_c_water():
    
    rhf_object = qp.RHF(mol, "sto-3g", use_c=True)
    rhf_energy = rhf_object.compute_energy()

    assert rhf_energy == pytest.approx(-74.9418022615317909, 1.e-5)

def test_df_rhf_water():
    
    rhf_object = qp.RHF(mol, "sto-3g", scf_type="DF")
    rhf_energy = rhf_object.compute_energy()

    assert rhf_energy == pytest.approx(-74.9418983834437853, 1.e-5)


    

