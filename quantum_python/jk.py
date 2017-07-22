"""
A file for JK builders.
"""

import time
import numpy as np
import psi4.core as pc
from . import core


def build_JK(molecule, basis_name, jk_type, use_c=True):
    """
    Construct a JK object of various types with the given specifications.

    Parameters
    ----------
    molecule : psi4.core.Molecule
        Psi4 Molecule object
    basis_name : str
        The requested orbital basis name
    jk_type : str (PK, DF)
        The type of JK object to construct

    Returns
    -------
    jk_object : JK
        A initialized JK object

    Notes
    -----
    For DF the basis is automatically constructed as the complementary JK object.

    Examples
    --------

    mol = psi4.geometry("He")
    jk = build_JK(mol, "cc-pVDZ", "DF")

    jk.compute_JK(C_left)
    ...

    """
    orbital = pc.BasisSet.build(molecule, "ORBITAL", basis_name)
    auxiliary = pc.BasisSet.build(molecule, None, "", "JKFIT", basis_name)

    if jk_type == "PK":
        return PKJK(orbital, use_c)
    elif jk_type == "DF":
        return DFJK(orbital, auxiliary, use_c)
    else:
        raise KeyError("build_JK: Unknown JK type '%s'" % jk_type)


class PKJK(object):
    """
    Constructs a "PK" JK object. This effectively holds two supermatrices which the inner product is
    then taken for speed.

    J_pq[D_rs] = I_pqrs D_rds
    K_pq[D_rs] = I.swapaxes(1,2)_pqrs D_rds


    """

    def __init__(self, orbital, use_c):
        self.orbital = orbital
        self.nbf = self.orbital.nbf()
        self.use_c = use_c

        self.mints = pc.MintsHelper(orbital)
        self.eri = np.array(self.mints.ao_eri())

    def compute_JK(self, C_left, C_right=None, do_J=True, do_K=True):

        if C_right is None:
            C_right = C_left

        D_list = []
        for Cl, Cr in zip(C_left, C_right):
            D_list.append(np.dot(Cl, Cr.T))

        if self.use_c:
            return self._compute_JK_c(D_list, do_J, do_K)
        else:
            return self._compute_JK_np(D_list, do_J, do_K)

    def _compute_JK_c(self, D_list, do_J, do_K):

        J_list = []
        K_list = []
        for D in D_list:
            J = np.zeros((self.nbf, self.nbf))
            K = np.zeros((self.nbf, self.nbf))

            # Call our C++ function
            core.compute_PKJK(self.eri, D, J, K)

            J_list.append(J)
            K_list.append(K)

        if not do_J:
            J_list = None
        if not do_K:
            K_list = None

        return (J_list, K_list)

    def _compute_JK_np(self, D_list, do_J, do_K):

        # Compute J
        if do_J:
            J = []
            for D in D_list:
                J.append(np.einsum('pqrs,rs->pq', self.eri, D))
        else:
            J = None

        # Compute K
        if do_K:
            K = []
            for D in D_list:
                K.append(np.einsum('prqs,rs->pq', self.eri, D))
        else:
            K = None

        return (J, K)


class DFJK(object):
    """
    Constructs a "DF" JK object. This exploits density-fitting to reduce the total number of
    integrals and compress the overall information in the ERI tensor.

    I_pqrs \approx g_pqQ g_Qrs
    Where g_Qrs = \phi_Q (r_1) (1 / r_{12} ) \phi_r (r_2) \phi_s (r_2)

    Uppercase letters (Q) denote a auxiliary basis.


    """

    def __init__(self, orbital, auxiliary, use_c):
        self.orbital = orbital
        self.auxiliary = auxiliary
        self.zero_basis = pc.BasisSet.zero_ao_basis_set()
        self.nbf = self.orbital.nbf()
        self.naux = self.auxiliary.nbf()
        self.use_c = use_c

        # Build required integrals
        self.mints = pc.MintsHelper(orbital)

        # Form (P|Q) ^ (-0.5)
        metric = self.mints.ao_eri(self.auxiliary, self.zero_basis, self.auxiliary, self.zero_basis)
        metric.power(-0.5, 1.e-14)
        self.metric = np.array(metric).squeeze()

        # Form AO Integrals
        Qpq = np.asarray(self.mints.ao_eri(self.auxiliary, self.zero_basis, self.orbital, self.orbital)).squeeze()

        # Apply metric contraction
        self.Ppq = np.dot(self.metric, Qpq.reshape(self.naux, -1)).reshape(self.naux, self.nbf, self.nbf)

    def compute_JK(self, C_left, C_right=None, do_J=True, do_K=True):

        if C_right is None:
            C_right = C_left

        if do_J:
            J = self.compute_J(C_left, C_right)
        else:
            J = None

        if do_K:
            K = self.compute_K(C_left, C_right)
        else:
            K = None

        return (J, K)

    def compute_J(self, C_left, C_right):

        D_list = []
        for Cl, Cr in zip(C_left, C_right):
            D_list.append(np.dot(Cl, Cr.T))

        # Form Ppq, pq -> P
        P_list = []
        for D in D_list:
            tmp = np.dot(self.Ppq.reshape(self.naux, -1), D.ravel())
            P_list.append(tmp)

        # Form rsP, P -> rs
        J_list = []
        for P in P_list:
            tmp = np.dot(P, self.Ppq.reshape(self.naux, -1)).reshape(self.nbf, self.nbf)

            # Make sure its symmetrized
            tmp = (tmp + tmp.T) / 2
            J_list.append(tmp)

        return J_list

    def compute_K(self, C_left, C_right):

        K_list = []
        for Cl, Cr in zip(C_left, C_right):

            # Qpq Cqi -> Qpi
            tmp_l = np.dot(self.Ppq.reshape(-1, self.nbf), Cl).reshape(self.naux, self.nbf, -1)
            if id(Cl) == id(Cr):
                tmp_r = tmp_l
            else:
                tmp_r = np.dot(self.Ppq.reshape(-1, self.nbf), Cr).reshape(self.naux, self.nbf, -1)

            # Qpi Qqi -> pq
            K = np.einsum('Ppi,Pqi->pq', tmp_l, tmp_r)

            K_list.append(K)

        return K_list
