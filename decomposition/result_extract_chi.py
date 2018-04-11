"""
Routines to generate a analysis of how the error types depend on the
Fidelity of the GHZ state used into making the stabilizer measurements
in the distributed surface code.

author: Eduardo Villasenor
created-on: 03/01/18
"""
import qutip as qt
import stabilizer
import noise_modeling
import pickle
import tools.names as names
# Improved parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 30.0
a1 = 1/30.
eta = 1/100.
theta = .63

# GHZ info
ghz_size = 4
stab_size = 4
protocol = "thres_eta"



for eta in [0.0100, 0.0095, 0.0090, 0.0085, 0.0080, 0.0075, 0.0070, 0.0065, 0.0060, 0.0055, 0.0050, 0.0045, 0.0040, 0.0035, 0.0030]:
# for a0 in [6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0, 10500.0, 11000.0, 11500.0, 12000.0, 12500.0, 13000.0]:
# for nothing in [0]:
    for parity in ["X", "Z"]:
        # Initialize objects
        model = noise_modeling.NoiseModel(stab_size, parity)
        model.separate_basis_parity()
        stab = stabilizer.Stabilizer(ps=ps, pm=pm, pg=pg)

        # Choi state for noise noise modeling
        choi = model._choi_state_ket(stab_size)
        choi = choi * choi.dag()
        targets = list(range(stab_size))


        # Load GHZ state
        ghz_file = names.ghz(ps, pm, pg, eta, a0, a1, theta,
                             ghz_size, protocol)
        ghz = qt.qload(ghz_file)
        if len(ghz.dims[0]) == 3:
            protocol += "_3on4"

        p_res, rhos = stab.measure_ghz_stabilizer(choi, ghz, targets, parity)

        # Set channel output and make chi matrix
        model.set_rho(rhos, p_res)
        model.make_chi_matrix()

        print("a0: ", a0)
        print("Total sum check: ", model.check_total_sum())

        I_OK = model.chi["IIII_OK"]
        I_NOK = model.chi["IIII_NOK"]
        # The sum of all physical errors
        E = (1 - model.chi["IIII_OK"] - model.chi["IIII_NOK"])/4.
        print("Errors:")
        print("OK: ", I_OK)
        print("NOK: ", I_NOK)
        print("E: ", E)
        print("-------------")

        file_name = names.chi(ps, pm, pg, eta, a0, a1, theta,
                              stab_size, parity, protocol)

        print(file_name)
        pickle_out = open(file_name, "wb")
        pickle.dump(model.chi, pickle_out, protocol=2)
        pickle_out.close()

        # print(model.chi)
        model.reset_chi()
