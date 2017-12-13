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
a0 = 20
a1 = 1/3.
pd = 1/1000.
theta = np.pi/4.
# Initialize  objects
cb = circuit_block.Blocks(ps, pm, pg, pn, pd, a0, a1, theta)
cb_ideal = circuit_block.Blocks(0, 0, 0, 0, 1, 0, 0, np.pi/4.)
rho_ref = qt.bell_state('00') * qt.bell_state('00').dag()


print("------------------NOISY BELL GENERATOR-------------------")
p, check, rho = cb.start_bell_pair()
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------SINGLE SELECTION-------------------")
_, _, rho = cb.start_bell_pair()
p, check, rho = cb.single_selection(rho, [0, 1], "X")
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

print("------------------DOUBLE SELECTION-------------------")
_, _, rho = cb.start_bell_pair()
p, check, rho = cb.double_selection(rho, [0, 1], "X")
print("p_success: ", p)
print("check: ", check)
print("F: ", qt.fidelity(rho, rho_ref))

# print("------------------TWO DOTS-------------------")
# n, rho = cb.generate_bell_pair()
# p, n, rho = cb.two_dots(rho, [0, 1], "X")
# print("p_success: ", p)
# print("n steps: ", n)
# print("F: ", qt.fidelity(rho, rho_ref))
