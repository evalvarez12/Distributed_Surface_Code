import qutip as qt
import numpy as np
import error_models as errs
import tools.operations as ops

# Fidelity of initial states
p = 0.1
def p_ref(f): return (f**2 +2*f*(1-f)/3 + 5*((1-f)/3)**2)
def single_selection(F1, F2):
    F = F1*F2 + (1 - F1)*(1 - F2)/9
    F = F/(F1*F2 + F1*(1 - F2)/3 + F2*(1 - F1)/3 + 5*(1 - F1)*(1 - F2)/9)
    return F

### Purification protocols for psi
print("--------------------PSI--------------------")
ref = qt.bell_state('10')
ref2 = qt.bell_state('11')
rho_i = errs.bell_pair_phi(p)
print("Initial F: ", qt.fidelity(ref, rho_i)**2)

rho_i = qt.tensor(rho_i, rho_i)

# Apply two_qubit_gates
CNOT1 = qt.cnot(4, 2, 0)
CNOT2 = qt.cnot(4, 3, 1)
GATE = CNOT1 * CNOT2
# CPHASE1 = qt.cphase(np.pi, 4, 0, 2)
# CPHASE2 = qt.cphase(np.pi, 4, 1, 3)
# GATE = CPHASE1 * CPHASE2
rho = GATE * rho_i * GATE.dag()

p1, rho = ops.forced_measure_single_Xbasis(rho, 4, 2, 0, True)
print("p1: ", p1)
p2, rho = ops.forced_measure_single_Xbasis(rho, 3, 2, 0, True)
print("p2: ", p2, p_ref(1-p))
print("F: ", qt.fidelity(ref, rho)**2, single_selection(1-p, 1-p))

print("--------------------")
print(qt.fidelity(ref, rho)**2)
H = qt.snot(2,0) * qt.snot(2,1)
rho = H * rho
rho2i = qt.tensor(rho,  errs.bell_pair_phi(p))
# rho2i = qt.tensor(rho, rho)
# CPHASE1 = qt.cphase(np.pi, 4, 2, 0)
# CPHASE2 = qt.cphase(np.pi, 4, 3, 1)
# GATE = CPHASE1 * CPHASE2
CNOT1 = qt.cnot(4, 2, 0)
CNOT2 = qt.cnot(4, 3, 1)
GATE = CNOT1 * CNOT2
rho2 = GATE * rho2i * GATE.dag()

p1, rho2 = ops.forced_measure_single_Xbasis(rho2, 4, 2, 0, True)
print("p1: ", p1)
p2, rho2 = ops.forced_measure_single_Xbasis(rho2, 3, 2, 0, True)
print("p2: ", p2, p_ref(0.9263959390862941))
print("F: ", qt.fidelity(ref, rho2)**2, single_selection(0.9, 0.9263959390862941))


## Purification protocols for phi+
print("--------------------PHI--------------------")
ref = qt.bell_state('00')
rho_i = errs.bell_pair(p)
print("Initial F: ", qt.fidelity(ref, rho_i)**2)

rho_i = qt.tensor(rho_i, rho_i)

# Apply two_qubit_gates
CNOT1 = qt.cnot(4, 2, 0)
CNOT2 = qt.cnot(4, 3, 1)
GATE = CNOT1 * CNOT2
rho = GATE * rho_i * GATE.dag()

p1, rho = ops.forced_measure_single_Xbasis(rho, 4, 2, 0, True)
print("p1: ", p1)
p2, rho = ops.forced_measure_single_Xbasis(rho, 3, 2, 0, True)
print("p2: ", p2, p_ref(1-p))
print("F: ", qt.fidelity(ref, rho)**2, single_selection(1-p, 1-p))


print("----------2----------")
print(qt.fidelity(ref, rho)**2)
H = qt.snot(2,0) * qt.snot(2,1)
rho = H * rho
rho2i = qt.tensor(rho,  errs.bell_pair(p))
# rho2i = qt.tensor(rho, rho)
# CPHASE1 = qt.cphase(np.pi, 4, 2, 0)
# CPHASE2 = qt.cphase(np.pi, 4, 3, 1)
# GATE = CPHASE1 * CPHASE2
CNOT1 = qt.cnot(4, 2, 0)
CNOT2 = qt.cnot(4, 3, 1)
GATE = CNOT1 * CNOT2
rho2 = GATE * rho2i * GATE.dag()

p1, rho2 = ops.forced_measure_single_Xbasis(rho2, 4, 2, 0, True)
print("p1: ", p1)
p2, rho2 = ops.forced_measure_single_Xbasis(rho2, 3, 2, 0, True)
print("p2: ", p2, p_ref(0.9263959390862941))
print("F: ", qt.fidelity(ref, rho2)**2, single_selection(0.9, 0.9263959390862941))


print("---------3-----------")
print(qt.fidelity(ref, rho2)**2)
H = qt.snot(2,0) * qt.snot(2,1)
# rho2 = H * rho2
rho2i = qt.tensor(rho2,  errs.bell_pair(p))
# rho2i = qt.tensor(rho, rho)
# CPHASE1 = qt.cphase(np.pi, 4, 2, 0)
# CPHASE2 = qt.cphase(np.pi, 4, 3, 1)
# GATE = CPHASE1 * CPHASE2
CNOT1 = qt.cnot(4, 2, 0)
CNOT2 = qt.cnot(4, 3, 1)
GATE = CNOT1 * CNOT2
rho2 = GATE * rho2i * GATE.dag()

p1, rho2 = ops.forced_measure_single_Xbasis(rho2, 4, 2, 0, True)
print("p1: ", p1)
p2, rho2 = ops.forced_measure_single_Xbasis(rho2, 3, 2, 0, True)
print("p2: ", p2, p_ref(0.9263959390862941))
print("F: ", qt.fidelity(ref, rho2)**2, single_selection(0.9, 0.9263959390862941))
