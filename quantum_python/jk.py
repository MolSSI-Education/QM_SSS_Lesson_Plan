"""
A file for JK builders.
"""


import time
import numpy as np
import psi4.core as pc


def DFJK(object):
    def __init__(self, orbital, auxiliary):
        self.orbital = orbital
        self.auxiliary = auxiliary
        self.zero_basis = pc.BasisSet.zero_ao_basis_set()
        self.nbf = self.orbital.nbf()
        self.naux = self.auxiliary.nbf()

        # Build required integrals
        self.mints = pc.MintsHelper(orbital)

        # Form (P|Q) ^ (-0.5)
        metric = self.mints.ao_eri(self.auxiliary, self.zero_basis, self.auxiliary, self.zero_basis)
        metric.power(-0.5, 1.e-14)
        self.metric = np.array(metric)

        # Form AO Integrals
        Qpq = self.mints.ao_eri(self.auxiliary, self.zero_basis, self.orbital, self.orbital)

        # Apply metric contraction
        self.Ppq = np.dot(self.metric, Qpq.reshape(self.naux, -1)).reshape(self.naux, self.nbf, self.nbf)


    def compute_JK(C_left, C_right=None, do_J=True, do_K=True):

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


    def compute_J(C_left, C_right):

        D_list = []
        for Cl, Cr in zip(C_left, C_right):
            D_list.append(np.dot(Cl, Cr.T))

        # Form Ppq, pq -> P
        P_list = []
        for D in D_list:
            tmp = np.dot(P.reshape(self.naux, -1), D.ravel())
            P_list.append(tmp)

        # Form rsP, P -> rs
        J_list = []
        for P in P_list:
            tmp = np.dot(P, Ppq.reshape(self.naux, -1)).reshape(self.nbf, self.nbf)

            # Make sure its symmetrized
            tmp = (tmp + tmp.T) / 2
            J_list.append(tmp)

        return J_list


    def compute_K(C_left, C_right):

        K_list = []
        for Cl, Cr in zip(C_left, C_right):

            tmp_l = np.dot(self.Ppq.reshape(-1, self.nbf), Cl)
            if id(Cl) == id(Cr):
                tmp_r = tmp_l
            else:
                tmp_r = np.dot(self.Ppq.reshape(-1, self.nbf), Cr)

            K = np.einsum('Pqi,Ppi->qp', tmp_l, tmp_r)

            K_list.append(K)

        return K_list

