import time
import numpy as np
import psi4.core as pc


class RHF(object):
    def __init__(self, molecule, basis, numpy_memory=2.e9):

        # Set defaults
        maxiter = 40
        E_conv = 1.0E-6
        D_conv = 1.0E-3

        # Integral generation from Psi4's MintsHelper
        start_time = time.time()
        self.molecule = molecule
        self.basis = pc.BasisSet.build(self.molecule, "ORBITAL", basis)
        self.mints = pc.MintsHelper(self.basis)
        self.S = np.asarray(self.mints.ao_overlap())

        # Get nbf and ndocc for closed shell molecules
        self.nbf = self.S.shape[0]
        self.nel = sum(self.molecule.Z(n) for n in range(self.molecule.natom()))
        self.nel -= self.molecule.molecular_charge()

        if not (self.nel / 2.0).is_integer():
            raise ValueError("RHF: Molecule did not have an even number of electrons!")

        self.ndocc = int(self.nel / 2.0)

        print('\nNumber of occupied orbitals: %d' % self.ndocc)
        print('Number of basis functions: %d' % self.nbf)

        # Run a quick check to make sure everything will fit into memory
        I_Size = (self.nbf**4) * 8.e-9
        print("\nSize of the ERI tensor will be %4.2f GB." % I_Size)

        # Estimate memory usage
        memory_footprint = I_Size * 1.5
        if I_Size > numpy_memory:
            pc.clean()
            raise Exception("Estimated memory utilization (%4.2f GB) exceeds numpy_memory \
                            limit of %4.2f GB." % (memory_footprint, numpy_memory))

        # Compute required quantities for SCF
        self.V = np.asarray(self.mints.ao_potential())
        self.T = np.asarray(self.mints.ao_kinetic())
        self.I = np.asarray(self.mints.ao_eri())

        self.Enuc = self.molecule.nuclear_repulsion_energy()

        print('\nTotal time taken for integrals: %.3f seconds.' % (time.time() - start_time))

        t = time.time()

        # Build H_core
        self.H = self.T + self.V

        # Orthogonalizer A = S^(-1/2) using Psi4's matrix power.
        A = self.mints.ao_overlap()
        A.power(-0.5, 1.e-16)
        self.A = np.asarray(A)

        print('\nTotal time taken for setup: %.3f seconds' % (time.time() - start_time))

    def form_orbitals(self, F):
        # Calculate orthogonal orbitals
        Fp = (self.A).dot(F).dot(self.A)
        e, C2 = np.linalg.eigh(Fp)

        # Set various quantities
        self.eps = e
        self.C = (self.A).dot(C2)
        self.Cocc = self.C[:, :self.ndocc]
        self.D = np.dot(self.Cocc, self.Cocc.T)

    def build_fock(self, D):
        # Build Fock matrix
        J = np.einsum('pqrs,rs->pq', self.I, D)
        K = np.einsum('prqs,rs->pq', self.I, D)
        self.F = self.H + J * 2 - K
        return self.F

    def compute_energy(self, maxiter=12, E_conv=1.e-6, D_conv=1.e-4):
        print('\nStarting SCF iterations:\n')

        t = time.time()
        E = 0.0
        Eold = 0.0
        Dold = np.zeros_like(self.H)

        self.form_orbitals(self.H)

        for SCF_ITER in range(1, maxiter + 1):

            # Build Fock matrix
            self.build_fock(self.D)

            # Build DIIS error vector
            diis_e = np.dot(self.F, self.D).dot(self.S) - np.dot(self.S, self.D).dot(self.F)

            # Make sure that error is normalized!
            diis_e = (self.A).dot(diis_e).dot(self.A)

            # SCF energy and update
            SCF_E = np.einsum('pq,pq->', self.F + self.H, self.D) + self.Enuc
            dRMS = np.mean(diis_e**2)**0.5

            print('SCF Iteration %3d: Energy = %4.16f   dE = % 1.5E   dRMS = %1.5E' % (SCF_ITER, SCF_E,
                                                                                       (SCF_E - Eold), dRMS))
            if (abs(SCF_E - Eold) < E_conv) and (dRMS < D_conv):
                break

            Eold = SCF_E
            Dold = self.D.copy()

            self.form_orbitals(self.F)

            if SCF_ITER == maxiter:
                pc.clean()
                raise Exception("Maximum number of SCF cycles exceeded.")

        print('Total time for SCF iterations: %.3f seconds \n' % (time.time() - t))

        print('Final SCF energy: %.8f hartree' % SCF_E)

        return SCF_E
