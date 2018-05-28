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
import tools.fgh as fgh
# Improved parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 5.0
a1 = 1/30.
eta = 1/100.
theta = .63

# GHZ info
ghz_size = 3
stab_size = 4
protocol = "thres_tau_paired"
protocol_chi = "thres_tau_paired"
mode = 3
extra = False
ignore_percent = 5

TIME = []

# for eta in [0.01, 0.0095, 0.0090, 0.0085, 0.0080, 0.0075, 0.0070, 0.0065, 0.0060, 0.0055, 0.0050]:
# for eta in [0.0005, 0.00075, 0.0010, 0.00125, 0.0015, 0.00175, 0.002]:
# for eta in [0.00055]:
# for a0 in [2000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0]:
# for a0 in [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]:
# for pg in [0.0031, 0.0032, 0.0033, 0.0034, 0.0035, 0.0036, 0.0037, 0.0038, 0.0039, 0.0040, 0.0041]:
# for pg in [0.0032, 0.0033, 0.0034, 0.0035, 0.0036, 0.0037, 0.0038, 0.0039, 0.0040, 0.0041, 0.0042, 0.0043, 0.0044, 0.0045]:
for tau in range(10):
    f, g, h = fgh.fgh(mode)
    a0 = f(tau)
    eta = g(tau)
    pg = h(tau)

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

        ghz = stab.twirl_ghz(ghz)

        # Choi state for noise noise modeling
        choi = model._choi_state_ket(stab_size)
        choi = choi * choi.dag()
        targets = list(range(stab_size))

        extra_name = ""
        if len(ghz.dims[0]) == 2:
            extra_name = "_extra"
        protocol_chi2 =  protocol_chi + extra_name

        p_res, rhos = stab.measure_ghz_stabilizer(choi, ghz, targets, parity)

        # Set channel output and make chi matrix
        model.set_rho(rhos, p_res)
        model.make_chi_matrix()

        print("PARAMS: ", a0, eta, pg)
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
                              stab_size, parity, protocol_chi2)

        print(file_name)
        pickle_out = open(file_name, "wb")
        pickle.dump(model.chi, pickle_out, protocol=2)
        pickle_out.close()

        # print(model.chi)
        model.reset_chi()

TIME = np.array(TIME)


def env_error_rate(t, a):
    # Function to calculate the error to the enviroment for step of stabilizers
    # measurements
    x = a * t
    p_env = (1 - np.exp(-x))/4.
    return p_env
