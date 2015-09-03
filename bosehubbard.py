# File: bosehubbard.py
#
# Efficient but simplistic approach to generating the
# many-body Hamiltonian for Bose-Hubbard models.
# Can be used for exact diagonalization.
#
# Code released under a MIT license, by
# Kim G. L. Pedersen, 2015
# (unless otherwise noted in function description)
#
# Waiver: No guarantees given. Please use completely at own risk.

import numpy as np
from scipy.special import binom
import scipy.sparse as sparse


class Model:
    """
    Defines some Bose-Hubbard model by specifying the onsite energies (omegas),
    the links between sites on the network, and
    the onsite interaction strength (U)

    Parameters
    ----------
    omegas : list
        Onsite energies
    links : list of lists
        Each hopping is of the form [from, to, amplitude].
    U : float
        Interaction strength
    """
    def __init__(self, omegas, links, U):
        self.omegas = np.array(omegas) - U/2
        self.links = links
        self.U = U
        self.n = len(omegas)

    @property
    def hopping(self):
        """
        The single particle hopping Hamiltonian
        """
        H0 = np.zeros([self.n]*2)

        for link in self.links:
            H0[link[0], link[1]] = link[2] if len(link) > 2 else -1

        return H0 + H0.T

    def chargestate(self, nb):
        """
        Returns a specific charge state object based on this model.
        """
        return ChargeState(self.n, nb, model=self)


class ChargeState:
    """
    Defines a specific charge state of a given Bose-Hubbard model.
    """
    def __init__(self, N, nb, model=None):
        self.N = N
        self.nb = nb
        self.basis = ChargeStateBasis(N, nb)
        if model is not None:
            self.model = model

    def hamiltonian(self):
        """
        Generates the (sparse) Hamiltonian
        """
        m = self.model
        nbas = self.basis.len

        HDi = np.arange(nbas)
        HD = self.onsite_hamiltonian(m.omegas, self.basis.vs) \
            + self.interaction_hamiltonian(m.U, self.basis.vs)
        Hki, Hkj, Hkv = self.hopping_hamiltonian(self.basis, m.hopping, self.basis.vs)

        return sparse.coo_matrix((Hkv, (Hki, Hkj)), shape=(nbas, nbas)).tocsr() \
            + sparse.coo_matrix((HD, (HDi, HDi)), shape=(nbas, nbas)).tocsr()

    @staticmethod
    def onsite_hamiltonian(omegas, states):
        """
        Onsite hamiltonian
        """
        return states.dot(omegas)

    @staticmethod
    def hopping_hamiltonian(basis, H0, states):
        """
        Hopping Hamiltonian expressed in the many-particle basis
        """
        H1s, H2s, Hvs = [], [], []

        for i in xrange(H0.shape[0]):
            js = np.nonzero(states[:, i])[0]  # affected states
            nj = len(js)

            ls = np.nonzero(H0[i, :])[0]  # relevant hoppings
            nl = len(ls)

            ks = np.zeros((nj*nl,))  # storing result states
            vs = np.zeros((nj*nl,))  # storing result states

            for k, l in enumerate(ls):
                nstates = states[js, :]
                nstates[:, i] -= 1  # remove one element
                nstates[:, l] += 1  # add it here

                ks[k*nj:(k+1)*nj] = basis.index(nstates)  # the new states
                vs[k*nj:(k+1)*nj] = H0[i, l]*np.sqrt(states[js, i]*(states[js, l] + 1))

            H1s += np.tile(js, nl).tolist()
            H2s += ks.tolist()
            Hvs += vs.tolist()

        return H1s, H2s, Hvs

    @staticmethod
    def interaction_hamiltonian(U, states):
        return np.sum(U/2*states**2, axis=1)


class ChargeStateBasis:
    """
    Many-body basis of specific <nb> charge state on a lattice with <N> sites.
    """
    def __init__(self, N, nb):
        """
        Parameters
        ----------
        N : int
            Number of sites
        nb : int
            Number of bosons
        """
        self.N = N  # number of sites
        self.nb = nb  # number of bosons

        self.len = ChargeStateBasis.size(N, nb)
        self.vs = ChargeStateBasis.generate(N, nb)
        self.hashes = ChargeStateBasis.hash(self.vs)
        self.sorter = ChargeStateBasis.argsort(self.hashes)

    def index(self, state):
        """
        Find the index of the state in self.basis.

        Parameters
        ----------
        state : ndarray
            One or more states.
        """
        return ChargeStateBasis.stateindex(state, self.hashes, self.sorter)

    @staticmethod
    def size(N, nb):
        """
        Return the size of the boson many-body basis. 

        Parameters
        ----------
        N : int
            Number of sites
        nb : int
            Number of bosons
        """
        return int(binom(nb+N-1, nb))

    @staticmethod
    def stateindex(state, hashes, sorter):
        """
        Converts state to hash and searches for the hash among hashes,
        which are sorted by the sorter list.
        
        Parameters
        ----------
        state : ndarray
            An array of one or more states
        hashes : ndarray
            List of hashes so search among
        sorter : ndarray
            Sorting indicies which sorts hashes
            (generated from ChargeStateBasis.argsort).
        """
        key = ChargeStateBasis.hash(state)
        return sorter[np.searchsorted(hashes, key, sorter=sorter)]

    @staticmethod
    def generate(N, nb):
        """
        Generate basis incrementally based on the method described in e.g.
        http://iopscience.iop.org/article/10.1088/0143-0807/31/3/016

        Parameters
        ----------
        N : int
            Number of sites
        nb : int
            Number of bosons
        """
        states = np.zeros((ChargeStateBasis.size(N, nb), N), dtype=int)

        states[0, 0] = nb
        ni = 0  # init
        for i in xrange(1, states.shape[0]):

            states[i, :N-1] = states[i-1, :N-1]
            states[i, ni] -= 1
            states[i, ni+1] += 1 + states[i-1, N-1]

            if ni >= N-2:
                if np.any(states[i, :N-1]):
                    ni = np.nonzero(states[i, :N-1])[0][-1]
            else:
                ni += 1

        return states

    @staticmethod
    def hash(states):
        """
        Hash function as given in:
        http://iopscience.iop.org/article/10.1088/0143-0807/31/3/016

        Parameters
        ----------
        states : ndarray
            List of states to hash
        """
        n = states.shape[1] if len(states.shape) > 1 else len(states)
        ps = np.sqrt(lowest_primes(n))
        return states.dot(ps)

    @staticmethod
    def argsort(hashes):
        """
        Argsort our hashes for searching, using e.g. quicksort.
        """
        return np.argsort(hashes, 0, 'quicksort')


def lowest_primes(n):
    """
    Return the lowest n primes

    Parameters
    ----------
    n : int
        Number of primes to return
    """
    return primes(n**2)[:n]


def primes(upto):
    """
    Prime sieve below an <upto> value.
    Copied from http://rebrained.com/?p=458

    Parameters
    ----------
    upto : int
        Find all primes leq this limit.
    """
    primes = np.arange(3, upto+1, 2)

    isprime = np.ones((upto-1)/2, dtype=bool)

    for factor in primes[:int(np.sqrt(upto))]:

        if isprime[(factor-2)/2]:
            isprime[(factor*3-2)/2::factor] = 0

    return np.insert(primes[isprime], 0, 2)
