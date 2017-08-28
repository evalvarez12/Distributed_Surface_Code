import qutip as qt
import protocols
import error_models as errs
import decomposition

rho = errs.bell_pair(.4)
rho_ideal = qt.bell_state("00")
print(decomposition.fidelity(rho, rho_ideal))

prot = protocols.Protocols(0.0, 0.0, 0., .3, 2)

print("------------------SINGLE SELECTION-------------------")
success, single = prot.single_selection(rho, [0, 1], "Z")
print(single)
if success:
    print(decomposition.fidelity(single, rho_ideal))

print("------------------ONE DOT-------------------")
success, one_dot = prot.one_dot(rho, [0, 1], "Z")
print(one_dot)
if success:
    print(decomposition.fidelity(one_dot, rho_ideal))

print("------------------TWO DOTS-------------------")
success, two_dot = prot.two_dots(rho, [0, 1], "X")
print(two_dot)
if success:
    print(decomposition.fidelity(two_dot, rho_ideal))

print("-------------------DOUBLE SELECTION-----------")
success, double = prot.double_selection(rho, [0, 1], "Z")
print(double)
if success:
    print(decomposition.fidelity(double, rho_ideal))


# print("-----------------EXPEDIENT------------------")
# m, rho = prot.expedient("Z")
# print(rho)
# print(m)

print("-----------------MONOLITHIC------------------")
rho_initial = qt.snot() * qt.basis(2, 0)
for i in range(3):
    rho_initial = qt.tensor(rho_initial, qt.snot()*qt.basis(2, 0))
rho_initial = rho_initial * rho_initial.dag()
m, rho = prot.monolitic(rho_initial, "X")
print(rho)
print(m)
