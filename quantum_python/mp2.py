import time
import numpy as np
import psi4.core as pc


class MP2(object):
    def __init__(self, wfn):

        # Copy over old wavefunctions data
        self.basis = wfn.basis
        self.mints = wfn.mints
        self.molecule = wfn.molecule
        self.nbf = wfn.nbf
        self.nel = wfn.nel
        self.ndocc = wfn.ndocc
        self.eps = wfn.eps
        self.C = wfn.C
        self.SCF_E = wfn.SCF_E

    def compute_energy(self):

        print('Building MO integrals.')
        # Integral generation from Psi4's MintsHelper
        t = time.time()
        mints = pc.MintsHelper(self.basis)
        Co = pc.Matrix.from_array(self.C[:, :self.ndocc])
        Cv = pc.Matrix.from_array(self.C[:, self.ndocc:])
        MO = np.asarray(mints.mo_eri(Co, Cv, Co, Cv))

        Eocc = self.eps[:self.ndocc]
        Evirt = self.eps[self.ndocc:]

        print('Shape of MO integrals: %s' % str(MO.shape))
        print('\n...finished SCF and integral build in %.3f seconds.\n' % (time.time() - t))

        print('Computing MP2 energy...')
        t = time.time()
        e_denom = 1 / (Eocc.reshape(-1, 1, 1, 1) - Evirt.reshape(-1, 1, 1) + Eocc.reshape(-1, 1) - Evirt)

        # Get the two spin cases
        MP2corr_OS = np.einsum('iajb,iajb,iajb->', MO, MO, e_denom)
        MP2corr_SS = np.einsum('iajb,iajb,iajb->', MO - MO.swapaxes(1, 3), MO, e_denom)
        print('...MP2 energy computed in %.3f seconds.\n' % (time.time() - t))

        MP2corr_E = MP2corr_SS + MP2corr_OS
        MP2_E = self.SCF_E + MP2corr_E

        SCS_MP2corr_E = MP2corr_SS / 3 + MP2corr_OS * (6. / 5)
        SCS_MP2_E = self.SCF_E + SCS_MP2corr_E

        print('MP2 SS correlation energy:         %16.10f' % MP2corr_SS)
        print('MP2 OS correlation energy:         %16.10f' % MP2corr_OS)

        print('\nMP2 correlation energy:            %16.10f' % MP2corr_E)
        print('MP2 total energy:                  %16.10f' % MP2_E)

        print('\nSCS-MP2 correlation energy:        %16.10f' % MP2corr_SS)
        print('SCS-MP2 total energy:              %16.10f' % SCS_MP2_E)

        return MP2_E
