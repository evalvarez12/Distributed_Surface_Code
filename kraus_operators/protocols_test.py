import qutip as qt
import protocols
import error_models as errs
import decomposition

rho = errs.bell_pair(.4)
rho_ideal = qt.bell_state("00")
print(decomposition.fidelity(rho, rho_ideal))

prot = protocols.Protocols(0.0, 0.0, 0., .4, 2)
succes, single = prot.single_selection(rho, [0, 1], "Z")
print(succes)
print(single)
print(decomposition.fidelity(single, rho_ideal))


succes, double = prot.double_selection(rho, [0, 1], "Z")
print(succes)
print(double)
print(decomposition.fidelity(double, rho_ideal))
