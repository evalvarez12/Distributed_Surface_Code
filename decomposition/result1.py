"""
File to generate a analysis of how the error types depend on the
Fidelity of the GHZ state used into making the stabilizer measurements
in the distributed surface code.
"""
import numpy as np
import matplotlib.pyplot as plt
import protocols_det
import noise_modeling
import error_models as errs

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
ps = 0.0
pm = 0.015
pg = 0.015
pn = 0.1
system_size = 4
parity = "X"

# Initialize objects
model = noise_modeling.NoiseModel(system_size, parity)
model.separate_basis_parity()
prot = protocols_det.ProtocolsDeterministic(ps, pm, pg, pm)

# Choi state for noise noise modeling
choi = model._choi_state_ket(system_size)
choi = choi * choi.dag()
targets = list(range(system_size))

I_OK_full = []
I_NOK_full = []
E_full = []

pgs = [0.006, 0.0075, 0.009, 0.0105]
fs = np.linspace(.5, 1, 50)
for pg in pgs:
    prot.change_parameters(ps, pm, pg, pn)
    I_OK = []
    I_NOK = []
    E = []


    for f in fs:
        ghz = errs.generate_noisy_ghz(f, system_size)
        ps, rhos = prot.measure_ghz_stabilizer(choi, ghz, targets, parity)
        model.set_rho(rhos, ps)
        model.make_chi_matrix()

        print("Total sum check: ", model.check_total_sum())
        I_OK += [model.chi["IIII_OK"]]
        I_NOK += [model.chi["IIII_NOK"]]
        # The sum of all physical errors
        E += [1 - model.chi["IIII_OK"] - model.chi["IIII_NOK"]]


        model.reset_chi()

    I_OK_full += [I_OK]
    I_NOK_full += [I_NOK]
    E_full += [E]

fig, ax = plt.subplots(nrows=3, ncols=1)
titles = iter([r"NO ERROR", r"MEASUREMENT ERROR", r"PHYSICAL ERROR"])
it = iter(pgs)
subplot_i = iter([1, 2, 3])

flag = True
for i in [I_OK_full, I_NOK_full, E_full]:
    plt.subplot(3, 1, next(subplot_i))
    for j in i:
        if flag:
            plt.plot(fs, j, label=r"$p_g=$" +str(next(it)))
            plt.legend()
        else:
            plt.plot(fs, j)
    plt.title(next(titles))
    flag = False
# for i in range(4):
#     for j in [I_OK_full[i], I_NOK_full[i], E_OK_full[i], E_NOK_full[i]]:
#         plt.plot(f, j)
#     plt.figure()

plt.show()
