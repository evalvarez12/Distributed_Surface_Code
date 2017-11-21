import qutip as qt
import circuit_block
import error_models as errs
import operations as ops
import numpy as np

# Determine parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.1

# Initialize  objects
cb = circuits_blocks.Blocks(ps, pm, pg, pn)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()


print("------------------NOISY BELL GENERATOR-------------------")
n, rho = cb.generate_bell_pair()
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------SINGLE SELECTION-------------------")
n, rho = cb.generate_bell_pair()
p, n, rho = cb.single_selection(rho, [0, 1], "X")
print("p_success: ", p)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------DOUBLE SELECTION-------------------")
n, rho = cb.generate_bell_pair()
p, n, rho = cb.double_selection(rho, [0, 1], "X")
print("p_success: ", p)
print("n steps: ", n)
print("F: ", qt.fidelity(rho, rho_ref))
