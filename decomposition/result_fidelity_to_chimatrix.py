"""
File to generate a analysis of how the error types depend on the
Fidelity of the GHZ state used into making the stabilizer measurements
in the distributed surface code.
"""
import matplotlib.pyplot as plt
import stabilizer
import noise_modeling
import error_models as errs
import tools.names as names
import protocols_nn
import pickle
import numpy as np

# Dummy function to find errors in the chi dictionary
def get_errors_dict(dict, ok):
    e = 0
    if "XIII_" + ok in dict :
        return dict["XIII_" + ok]
    if "IXII_" + ok in dict :
        return dict["IXII_" + ok]
    if "IIXI_" + ok in dict :
        return dict["IIXI_" + ok]
    if "IIIX_" + ok in dict :
        return dict["IIIX_" + ok]



# Initial parameters
ps = 0.0075
pm = 0.0075
pg = 0.0075
system_size = 4
parity = "X"

eta = 0.0
a0 = 0.0
a1 = 0.0
theta = 0.0
protocol = "LOCAL"

# Initialize objects
model = noise_modeling.NoiseModel(system_size, parity)
model.separate_basis_parity()
stab = stabilizer.Stabilizer(ps=ps, pm=pm, pg=pg)

# # Choi state for noise noise modeling
choi = model._choi_state_ket(system_size)
choi = choi * choi.dag()
targets = list(range(system_size))

I_OK_full = []
I_NOK_full = []
E_full = []

pgs = np.arange(0.0030, 0.0101, 0.0001)
# pgs = [0.003, 0.006, 0.009, 0.012]
# pgs = [0.0075]
fs = np.linspace(.5, 1, 50)

for parity in ["X", "Z"]:
    for pg in pgs:
        ps = pg
        pm = pg
        pg = pg
        pn = 0.1

        stab.change_parameters(ps=ps, pm=pm, pg=pg)
        # stab = protocols_nn.Protocols(ps, pm, pg, pn)
        I_OK = []
        I_NOK = []
        E = []


        for f in [0]:
            # ghz = errs.generate_noisy_ghz(f, system_size)
            # probs, rhos = stab.measure_ghz_stabilizer(choi, ghz, targets, parity)
            # probs, rhos = stab.local_stabilizer(choi, targets, parity)
            # model.set_rho(rhos, probs)

            # Define function and apply superoperator
            superoperator_function = stab.local_stabilizer
            # superoperator_function = stab.stringent
            model.apply_superoperator(superoperator_function)
            model.make_chi_matrix()

            print("Total sum check: ", model.check_total_sum())
            print("LEN: ", len(model.chi))
            I_OK += [model.chi["IIII_OK"]]
            I_NOK += [model.chi["IIII_NOK"]]
            # The sum of all physical errors
            E += [(1 - model.chi["IIII_OK"] - model.chi["IIII_NOK"])/4.]

            print("Errors:", pg)
            print("OK: ", I_OK)
            print("NOK: ", I_NOK)
            print("E: ", E)
            print("-------------")

            file_name = names.chi(ps, pm, pg, eta, a0, a1, theta,
                                  system_size, parity, protocol)

            pickle_out = open(file_name, "wb")
            pickle.dump(model.chi, pickle_out, protocol=2)
            pickle_out.close()

            model.reset_chi()
            # print(model.chi)
            # chi = model.chi

            # Remove the 3 qubit errors part
            # for i in list(chi.keys()):
            #     if i.count('I') < 2:
            #         chi[i] = 0
            #
            # rest = 1 - sum(chi.values())
            # chi['IIII_OK'] += rest

            # pickle_out = open(file_name, "wb")
            # pickle.dump(chi, pickle_out, protocol=2)
            # pickle_out.close()

    I_OK_full += [I_OK]
    I_NOK_full += [I_NOK]
    E_full += [E]

# np.save("I_OK_full.npy", I_OK_full)
# np.save("I_NOK_full.npy", I_NOK_full)
# np.save("E_full.npy", E_full)
############################ GRAPH STUFF
# I_OK_full = np.load("I_OK_full.npy")
# I_NOK_full = np.load("I_NOK_full.npy")
# E_full = np.load("E_full.npy")
#
#
# fig, ax = plt.subplots(nrows=3, ncols=1)
# titles = iter([r"NO ERROR", r"MEASUREMENT ERROR", r"PHYSICAL ERROR"])
# it = iter(pgs)
# subplot_i = iter([1, 2, 3])
#
# flag = True
# for i in [I_OK_full, I_NOK_full, E_full]:
#     plt.subplot(3, 1, next(subplot_i))
#     for j in i:
#         if flag:
#             plt.plot(fs, j, label=r"$p_g=$" +str(next(it)))
#             plt.legend()
#         else:
#             plt.plot(fs, j)
#     plt.ylabel(next(titles), fontsize=17)
#     # plt.xlabel(r"Fidelity", fontsize=17)
#     plt.xticks(fontsize=15)
#     plt.yticks(fontsize=15)
#     plt.legend(fontsize=17)
#     plt.tight_layout()
#     flag = False
#
# plt.xlabel(r"Fidelity", fontsize=17)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.legend(fontsize=17)
# plt.tight_layout()
# plt.savefig('ftoerrors.pdf', format='pdf', dpi=300)
# plt.show()
