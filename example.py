import bosehubbard

# The model parameters
omegas = [0]*8
links = [[i, (i+1) % 8, .1] for i in xrange(8)]
U = 2

# Construc the model
m = bosehubbard.Model(omegas, links, U)

# The single-particle hopping becomes,
print('Hopping Hamiltonian:')
print(m.hopping)

# Construct a specific chargestate with two bosons:
m2 = m.numbersector(2)

# Construct the many-body Hamiltonian (in sparse format)
H = m2.hamiltonian
print('Many-Body Hamiltonian:')
print(H)

# One can access the many-body basis from the `Basis` object,
vs = m2.basis.vs
print('Many-Body Basis:')
print(vs)