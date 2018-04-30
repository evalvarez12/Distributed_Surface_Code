import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import tools.names as names
import stabilizer
import matplotlib.patches as mpatches


# Improved parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 5.0
a1 = 1/30.
eta = 1/100.
theta = .63

stab = stabilizer.Stabilizer(ps=ps, pm=pm, pg=pg)

# GHZ info
ghz_size = 4
protocol = "thres_pg_bk"
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

ignore_percent = 10
extra = False

def extract_ft(ignore_percent, ghzs, times):
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
    ghz = stab.twirl_ghz(ghz)
    fidelity = qt.fidelity(ghz, ghz_ref)

    print("_____________________________________________________")
    print("N: ", N)
    print("T_avg: ", tavg, tstd)
    print("TIME_MAX:", t_max)
    print("F: ", fidelity)


    return fidelity, t_max


it = iter(["o", "v", "s"])
# for eta in [0.01, 0.0095, 0.0090, 0.0085, 0.0080, 0.0075, 0.0070, 0.0065, 0.0060, 0.0055, 0.0050, 0.0045]:
# for a0 in [2000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0]:
# for a0 in [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]:
for protocol in ["thres_pg", "thres_pg_singlesel", "thres_pg_doublesel"]:
    TIME = []
    FIDELITY = []
    for pg in [0.0031, 0.0032, 0.0033, 0.0034, 0.0035, 0.0036, 0.0037, 0.0038, 0.0039, 0.0040, 0.0041]:
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
        f, t_max = extract_ft(ignore_percent, ghzs, times)
        FIDELITY += [f]
        TIME += [t_max]

        ################################################


    time = np.array(TIME)
    fidelity = np.array(FIDELITY)

    time = np.sort(time)
    # fidelity = np.sort(fidelity)

    plt.plot(time, fidelity, "r" + next(it) + "-")

# Improved parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 5.0
a1 = 1/30.
eta = 1/100.
theta = .63

it = iter(["o", "v", "s"])

for protocol in ["thres_eta", "thres_eta_singlesel", "thres_eta_doublesel"]:
    TIME = []
    FIDELITY = []
    for eta in [0.01, 0.0095, 0.0090, 0.0085, 0.0080, 0.0075, 0.0070, 0.0065, 0.0060, 0.0055, 0.0050]:
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
        f, t_max = extract_ft(ignore_percent, ghzs, times)
        FIDELITY += [f]
        TIME += [t_max]

        ################################################


    time = np.array(TIME)
    fidelity = np.array(FIDELITY)

    time = np.sort(time)
    # fidelity = np.sort(fidelity)

    plt.plot(time, fidelity, "b" + next(it) + "-")


# Improved parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 5.0
a1 = 1/30.
eta = 1/100.
theta = .63

it = iter(["o", "v", "s"])

for protocol in ["thres_a0", "thres_a0_singlesel", "thres_a0_doublesel"]:
    TIME = []
    FIDELITY = []
    for a0 in [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]:
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
        f, t_max = extract_ft(ignore_percent, ghzs, times)
        FIDELITY += [f]
        TIME += [t_max]

        ################################################


    time = np.array(TIME)
    fidelity = np.array(FIDELITY)

    time = np.sort(time)
    # fidelity = np.sort(fidelity)

    plt.plot(time, fidelity, "g" + next(it) + "-")



plt.ylabel(r"Fidelity", fontsize=17)
plt.xlabel(r"Time (seg)", fontsize=17)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

blue_patch = mpatches.Patch(color='blue', label=r'$\eta$')
red_patch = mpatches.Patch(color='red', label=r'$p_g$')
green_patch = mpatches.Patch(color='green', label=r'$a_0$')

plt.legend(fontsize=14, handles=[blue_patch, red_patch, green_patch])

plt.axis([0.059, 0.25, .942, 0.981])

plt.tight_layout()
plt.savefig('FvsT.pdf', format='pdf', dpi=300)
plt.show()
