"""
Test file for the functions on Circuit Block.

author: Eduardo Villasenor
created-on: 20/11/17
"""

import qutip as qt
import numpy as np
import stabilizer

# Determine parameters
ps = 0.006
pm = 0.006
pg = 0.006




# Initialize  objects
stab = stabilizer.Stabilizer(ps, pm, pg)
stab_ideal = stabilizer.Stabilizer(0, 0, 0)

print("------------------TEST SWAP PAIR--------------------------")
pair = qt.bell_state('00') * qt.bell_state('00').dag()
pair = qt.tensor(pair, qt.basis(2,0) * qt.basis(2,0).dag())

pair_swapped = stab._swap_pair(pair, [0, 1])
print(pair)
print(pair_swapped)
print(qt.fidelity(pair, pair_swapped)**2)


print("------------------TEST SWAP GHZ--------------------------")
ghz = qt.ghz_state(4) * qt.ghz_state(4).dag()

ghz_swapped = stab._swap_ghz(ghz)
# print(ghz_swapped)
print(qt.fidelity(ghz, ghz_swapped)**2)
