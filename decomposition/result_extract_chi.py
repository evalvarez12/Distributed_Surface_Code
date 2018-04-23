"""
Routines to generate a analysis of how the error types depend on the
Fidelity of the GHZ state used into making the stabilizer measurements
in the distributed surface code.

author: Eduardo Villasenor
created-on: 03/01/18
"""
import qutip as qt
import numpy as np
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
protocol = "thres_a0"

extra = True
ignore_percent = 5

TIME = []

# for eta in [0.0100, 0.0095, 0.0090, 0.0085, 0.0080, 0.0075, 0.0070, 0.0065, 0.0060, 0.0055, 0.0050, 0.0045, 0.0040, 0.0035, 0.0030]:
# for a0 in [1000.0, 1500.0, 2000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7000.0]:
for a0 in [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]:
# for pg in [0.00325]:
    ps = pg
    pm = pg
    # Load GHZ state files
    ghz_file = names.ghz(ps, pm, pg, eta, a0, a1, theta,
                         ghz_size, protocol)
    times_file = names.ghz_times(ps, pm, pg, eta, a0, a1, theta,
                                 ghz_size, protocol)

    ghzs = qt.qload(ghz_file)
    times = np.load(times_file)
    if extra:
        ghz_file2 = names.ghz(ps, pm, pg, eta, a0, a1, theta,
                              ghz_size, protocol+"2")
        times_file2 = names.ghz_times(ps, pm, pg, eta, a0, a1, theta,
                                      ghz_size, protocol+"2")

        ghzs = np.append(ghzs, qt.qload(ghz_file2))
        times = np.append(times, np.load(times_file2))

    ###################################################################
    # CHOSE PERCENTILE
    N = len(times)
    ignore_number = int(N*ignore_percent/100)

    indices_sorted = np.argsort(times)
    t_sorted = times[indices_sorted]

    if ignore_number != 0:
        t_max = t_sorted[:-ignore_number][-1]
        tavg = np.average(times[indices_sorted][:-ignore_number])
        tstd = np.std(times[indices_sorted][:-ignore_number])
        ghz = np.sum(ghzs[indices_sorted][:-ignore_number])
    else:
        t_max = t_sorted[-1]
        tavg = np.average(times)
        tstd = np.std(times)
        ghz = np.sum(ghzs)

    ghz = ghz/(N - ignore_number)
    TIME += [t_max]
    # print("F: ", qt.fidelity(rho, rho_ref))
    print("_____________________________________________________")
    print("N: ", N)
    print("T_avg: ", tavg, tstd)
    print("TIME_MAX:", t_max)
    ################################################
    # CHOI decomposition
    for parity in ["X", "Z"]:
        # Initialize objects
        model = noise_modeling.NoiseModel(stab_size, parity)
        model.separate_basis_parity()
        stab = stabilizer.Stabilizer(ps=ps, pm=pm, pg=pg)

        # Choi state for noise noise modeling
        choi = model._choi_state_ket(stab_size)
        choi = choi * choi.dag()
        targets = list(range(stab_size))

        if len(ghz.dims[0]) == 3:
            protocol += "_3on4"
        elif len(ghz.dims[0]) == 2:
            protocol += "_2on4"

        ghz = stab.twirl_ghz(ghz)
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

TIME = np.array(TIME)
