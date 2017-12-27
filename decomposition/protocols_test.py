"""
File to test with different protocols to purify a Bell state.

author: Eduardo Villasenor
created-on: 21/11/17
"""

import numpy as np
import qutip as qt
import collections
import circuit_block
import circuit
import error_models

# Determine parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 83.33
a1 = 1/3.
eta = (0.1)*(0.03)*(0.8)
theta = .86

iterations = 1

# Initialize objects
cb = circuit_block.Blocks(ps, pm, pg, eta, a0, a1, theta)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

# Lists to save results
SIMPLE_fidelity = []
SIMPLE_steps = []
MEDIUM_fidelity = []
MEDIUM_steps = []
COMPLEX_fidelity = []
COMPLEX_steps = []

'''PROTOCOLS START HERE'''
print("-------------------PROTOCOL TEST------------------")
test = circuit.Circuit(a0=a0, a1=a1,
                                     circuit_block=cb.start_epl)
test.add_circuit(circuit_block=cb.single_selection,
                               operation_qubits=[0, 1],
                               sigma="X")
p, c, rho = test.run(None)
fidelity = qt.fidelity(rho, rho_ref)
print(fidelity)
