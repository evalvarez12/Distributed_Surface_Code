import qutip as qt
import protocols_nn as pr
import error_models as errs
import tools.operations as ops
import tools.p_success as ps

rho = errs.bell_pair(.4)
rho_ideal = qt.bell_state("00")
# print(de.fidelity(rho, rho_ideal))

prot = pr.Protocols(0.0, 0.0075, 0.0075, .1)
prot_perf = pr.Protocols(0., 0., 0., 0.,)

print("------------------SINGLE SELECTION-------------------")
single = prot.single_selection(rho, [0, 1], "Z")
# print(single)
print(qt.fidelity(single, rho_ideal))

print("------------------ONE DOT-------------------")
one_dot = prot.one_dot(rho, [0, 1], "Z")
# print(one_dot)
print(qt.fidelity(one_dot, rho_ideal))

print("------------------TWO DOTS-------------------")
two_dot = prot.two_dots(rho, [0, 1], "X")
# print(two_dot)
print(qt.fidelity(two_dot, rho_ideal))

print("-------------------DOUBLE SELECTION-----------")
double = prot.double_selection(rho, [0, 1], "Z")
# print(double)
print(qt.fidelity(double, rho_ideal))

print("-------------------EXPEDIENT------------------")
state = qt.rand_ket(16, dims=[[2]*4, [1]*4])
state = state * state.dag()
probs, rhos = prot.expedient(state, list(range(4)), "X")
print(sum(probs))

print("-------------------STRINGENT------------------")
state = qt.rand_ket(16, dims=[[2]*4, [1]*4])
state = state * state.dag()
probs, rhos = prot.expedient(state, list(range(4)), "X")
print(sum(probs))

# print("-----------------LOCAL STABILIZER------------------")
# rho_initial = qt.snot() * qt.basis(2, 0)
# for i in range(3):
#     rho_initial = qt.tensor(rho_initial, qt.snot() * qt.basis(2, 1))
# rho_initial = rho_initial * rho_initial.dag()
# m, rho = prot_perf.local_stabilizer(rho_initial, [0, 1, 2, 3], "X")
# # print(rho)
# print(m)

# print("--------------STRINGENT/EXPEDIENT")
# rho_initial = qt.snot() * qt.basis(2, 0)
# for i in range(2):
#     rho_initial = qt.tensor(rho_initial, qt.snot() * qt.basis(2, 1))
# rho_initial = rho_initial * rho_initial.dag()
# a = prot.stringent(rho_initial, [0, 1, 2], "X")
# b = prot.expedient(rho_initial, [0, 1, 2], "X")


print("----GHZ FIDELITY--------------")
ghz = prot.make_ghz_expedient(4)
ghz_ideal = prot_perf.make_ghz_expedient(4)
print((ghz_ideal * ghz * ghz_ideal).tr())

print("-----------P SUCCESS SINGLE SELECTION----------")
rho = prot.generate_bell_pair()
print("Fidelity: ", qt.fidelity(rho, qt.bell_state('00')))
# Generate raw bell pair
rho = qt.tensor(rho, prot.generate_bell_pair())
N = len(rho.dims[0])

# Apply two qubit gates
CNOT = qt.cnot(N, 2, 0) * qt.cnot(N, 3, 1)
rho = CNOT * rho * CNOT.dag()

p = ps.single_sel(rho, N, [2, 3])
def p_ref(f): return (f**2 +2*f*(1-f)/3 + 5*((1-f)/3)**2)
print("p: ", p)
print("p_ref: ", p_ref(.9))

print("-----------P SUCCESS DOUBLE SELECTION----------")
rho = prot.generate_bell_pair()
# Generate raw bell pair
rho = qt.tensor(rho, prot.generate_bell_pair())
rho = qt.tensor(rho, prot.generate_bell_pair())
N = len(rho.dims[0])

# Apply first set of two qubit gates
CNOT = qt.cnot(N, 2, 0) * qt.cnot(N, 3, 1)
rho = CNOT * rho * CNOT.dag()
# Second set
CNOT = qt.cnot(N, 4, 2) * qt.cnot(N, 5, 3)
rho = CNOT * rho * CNOT.dag()

p = ps.double_sel(rho, N, [2, 3], [4,5])
print(p)
