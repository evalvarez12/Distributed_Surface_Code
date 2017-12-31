"""
File to generate a analysis of how the error types depend on the
Fidelity of the GHZ state used into making the stabilizer measurements
in the distributed surface code.
"""
import numpy as np
import qutip as qt
from os.path import dirname, realpath
import matplotlib.pyplot as plt
import protocols_det
import noise_modeling
import error_models as errs
import pickle



def generate_name(ps, pm, pg, eta, a0, a1, theta, stab_size, parity, protocol):
    param_names = ["ps=" + str(ps), "pm=" + str(pm),
                   "pg=" + str(pg), "eta=" + str(eta),
                   "a0=" + str(a0), "a1=" + str(a1),
                   "theta=" + str(theta)]

    param_names = "_".join(param_names)
    file_name = [protocol, parity, str(stab_size)]
    file_name = "_".join(file_name)
    # The address of the parent directory
    script_path = dirname(dirname(realpath(__file__)))
    file_name = (script_path + "/data/" + file_name
                 + "_" + param_names + ".dict")
    return file_name


# Initial parameters
ps = 0.003
pm = 0.003
pg = 0.003
system_size = 4

parity = "Z"

eta = 1/100.
theta = .24
a1 = 1/80.

# Initialize objects
model = noise_modeling.NoiseModel(system_size, parity)
model.separate_basis_parity()
prot = protocols_det.ProtocolsDeterministic(ps=ps, pm=pm, pg=pg, pn=0)

# Choi state for noise noise modeling
choi = model._choi_state_ket(system_size)
choi = choi * choi.dag()
targets = list(range(system_size))

for a0 in [12., 10., 8., 6., 4., 2.]:
    prot.change_parameters(ps=ps, pm=pm, pg=pg, pn=0)

    ghz_file = "ghz_4_a0_" + str(round(a0, 3))
    ghz = qt.qload(ghz_file)

    p_res, rhos = prot.measure_ghz_stabilizer(choi, ghz, targets, parity)
    model.set_rho(rhos, p_res)
    model.make_chi_matrix()

    print("a0: ", a0)
    print("Total sum check: ", model.check_total_sum())

    file_name = generate_name(ps, pm, pg, eta, a0, a1, theta,
                              system_size, parity, protocol="GHZ")

    pickle_out = open(file_name, "wb")
    pickle.dump(model.chi, pickle_out, protocol=2)
    pickle_out.close()


    model.reset_chi()
