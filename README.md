# bosehubbard

Lite python library to generate the basis and many-body Hamiltonian for the Bose-Hubbard model.

Requires numpy and scipy.sparse installed on your system.

## How to use it:

As an example consider a ring with 8 sites. 

Construct the model from a list of onsite-energies `omegas`, a list of links and hopping strengths, and the interaction energy `U`.

  import bosehubbard as bohu

  # The model
  omegas = [0]*8
  links = [[i, (i+1) % 8, 1] for i in xrange(8)]
  U = 2

  m = bohu.Model(omegas, links, U)

  # Construct a specific chargestate with two bosons:
  m2 = m.chargestate(2)

  # Construct the many-body Hamiltonian (in sparse format)
  H = m2.hamiltonian()


One can access the many-body basis from the `Basis` object,

  vs = m2.basis.vs
  print(vs)

Produces,
  [[2 0 0 0 0 0 0 0]
   [1 1 0 0 0 0 0 0]
   [1 0 1 0 0 0 0 0]
   [1 0 0 1 0 0 0 0]
   [1 0 0 0 1 0 0 0]
   [1 0 0 0 0 1 0 0]
   [1 0 0 0 0 0 1 0]
   [1 0 0 0 0 0 0 1]
   [0 2 0 0 0 0 0 0]
          ...
   [0 0 0 0 0 0 1 1]
   [0 0 0 0 0 0 0 2]]
