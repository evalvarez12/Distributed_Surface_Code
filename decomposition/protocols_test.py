import qutip as qt
import protocols
import error_models as errs
import numpy as np

rho = errs.bell_pair(.4)
rho_ideal = qt.bell_state("00")
# print(de.fidelity(rho, rho_ideal))

prot = protocols.Protocols(0.0, 0.0075, 0.0075, .1)
prot_perf = protocols.Protocols(0., 0., 0., 0.,)

print("------------------SINGLE SELECTION-------------------")
p, single = prot.single_selection(rho, [0, 1], "Z")
print("p = ", p)
# print(single)
print(qt.fidelity(single, rho_ideal))

print("------------------ONE DOT-------------------")
p, one_dot = prot.one_dot(rho, [0, 1], "Z")
print("p = ", p)
# print(one_dot)
print(qt.fidelity(one_dot, rho_ideal))

print("------------------TWO DOTS-------------------")
p, two_dot = prot.two_dots(rho, [0, 1], "X")
print("p = ", p)
# print(two_dot)
print(qt.fidelity(two_dot, rho_ideal))

print("-------------------DOUBLE SELECTION-----------")
p, double = prot.double_selection(rho, [0, 1], "Z")
print("p = ", p)
# print(double)
print(qt.fidelity(double, rho_ideal))


print("-----------------LOCAL STABILIZER------------------")
rho_initial = qt.snot() * qt.basis(2, 0)
for i in range(3):
    rho_initial = qt.tensor(rho_initial, qt.snot() * qt.basis(2, 1))
rho_initial = rho_initial * rho_initial.dag()
m, rho = prot_perf.local_stabilizer(rho_initial, [0, 1, 2, 3], "X")
print(rho)
print(m)

print("--------------STRINGENT/EXPEDIENT")
rho_initial = qt.snot() * qt.basis(2, 0)
for i in range(2):
    rho_initial = qt.tensor(rho_initial, qt.snot() * qt.basis(2, 1))
rho_initial = rho_initial * rho_initial.dag()
a = prot.stringent(rho_initial, [0, 1, 2], "X")
b = prot.expedient(rho_initial, [0, 1, 2], "X")


print("----GHZ Fidelity--------------")
ghz = prot.make_ghz_expedient(4)
ghz_ideal = prot_perf.make_ghz_expedient(4)
print((ghz_ideal * ghz * ghz_ideal).tr())
