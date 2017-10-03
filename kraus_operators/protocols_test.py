import qutip as qt
import protocols
import error_models as errs

rho = errs.bell_pair(.4)
rho_ideal = qt.bell_state("00")
# print(de.fidelity(rho, rho_ideal))

prot = protocols.Protocols(0.1, 0.1, 0.1, .1)
prot_ideal = protocols.Protocols(0, 0, 0, 0)
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



print("-----------------MONOLITHIC------------------")
rho_initial = qt.snot() * qt.basis(2, 1)
for i in range(3):
    rho_initial = qt.tensor(rho_initial, qt.snot()*qt.basis(2, 1))
rho_initial = rho_initial * rho_initial.dag()
m, rho = prot.monolithic(rho_initial, [0, 1, 2, 3], "X")
print(rho)
print(m)


# print("----------------PARITY MEASUREMENT-------------")
# psi = qt.bell_state('00')
# targets = [0, 2]
# psi_collapsed = prot.parity_projection_ket(psi, targets, 1, "Z")
# print(psi_collapsed)
#
#
# print("--------------PERFECT COMPARATION-------------")
# # Compare parity measurements with perfect case
# psi = qt.tensor(qt.bell_state('00'), qt.bell_state('00'))
# psi = qt.tensor(psi, psi)
# targets = [0, 2, 4, 6]
# psi_ref0 = prot_perf.parity_projection_ket(psi, targets, 0, "Z")
# rho_ref0 = psi_ref0 * psi_ref0.dag()
# psi_ref1 = prot_perf.parity_projection_ket(psi, targets, 1, "Z")
# rho_ref1 = psi_ref1 * psi_ref1.dag()

#
# m, rho_perf = prot_perf.monolithic(psi * psi.dag(), targets, "Z")
# print(m)
# print(rho_ref0 == rho_perf)
# print(rho_ref1 == rho_perf)
