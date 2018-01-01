"""
Simple simulation to test the surface code simulation is working
without looking at plots.

created-on: 09/12/17
@author: eduardo
"""
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from os.path import dirname, realpath
from mpi4py import MPI
import surface_code
import layers
import matching


def lambda_env(t, a0, a1):
    a = (a0 + a1)*t
    lamb = (1 + np.exp(-a * t))/2.
    return 1 - lamb


# Start the comm for mpi4py
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Get the arguments
pythonfile = sys.argv[0]
args = dict(arg.split('=') for arg in sys.argv[1:])

# Parameters for noisy measurement
ps = 0.003
pm = 0.003
pg = 0.003
eta = 0.01
a0 = 2.0
a1 = 1/80.
protocol = "GHZ"
theta = .24
NOISY_MEASUREMENT = True

# Set parameters
distance = int(args["distance"])
topology = args["topology"]
iterations = int(args["iterations"])
p = float(args["p"])
q = float(args["q"])
if q != 0:
    cycles = int(args["cycles"])
else:
    cycles = 0

"""
Surface code simulations for one set of parameters
"""
# Initialize fail rate
fail_rate = 0

# Initialize objects
sc = surface_code.SurfaceCode(distance, topology)
lc = layers.Layers(sc)
sc.init_error_obj(topology, ps, pm, pg, eta, a0, a1, theta, protocol)

# Set time for each GHZ generation
t = 0.30347
lamb = lambda_env(t, a0, a1)

# Perform measurements
for i in range(iterations):

    # Errors and measurements
    # if q != 0:
    #     for t in range(cycles):
    #         sc.apply_qubit_error(p, 0)
    #         sc.measure_all_stablizers()
    #         sc.apply_measurement_error(q)
    #         lc.add()
    # else:
    #     sc.apply_qubit_error(p, 0)
    #     sc.measure_all_stablizers()
    #     lc.add()

    # Noisy measurements
    for t in range(cycles):
        # sc.noisy_measurement("star")
        # sc.noisy_measurement("plaq")
        sc.noisy_measurement_cycle(lamb)
        lc.add()


    # Get anyons
    anyons_star, anyons_plaq = lc.find_anyons_all()

    # Decode
    match_star = matching.match(distance, anyons_star, topology,
                                "star", time=cycles, weights=[1, 1])
    match_plaq = matching.match(distance, anyons_plaq, topology,
                                "plaq", time=cycles, weights=[1, 1])

    # Apply corrections
    sc.correct_error("star", match_star)
    sc.correct_error("plaq", match_plaq)

    # Round of perfect detection to eliminate stray errors
    if q!= 0 or NOISY_MEASUREMENT:
        lc.reset()
        sc.measure_all_stablizers()
        lc.add()
        anyons_star, anyons_plaq = lc.find_anyons_all()
        match_star = matching.match(distance, anyons_star, topology,
                                    "star", time=0, weights=[1, 1])
        match_plaq = matching.match(distance, anyons_plaq, topology,
                                    "plaq", time=0, weights=[1, 1])
        sc.correct_error("star", match_star, cycles)
        sc.correct_error("plaq", match_plaq, cycles)

    # # Check for errors in decoding and correcting
    # sc.measure_stabilizer_type("star")
    # sc.measure_stabilizer_type("plaq")
    # if (sc.qubits[:, sc.tags != "Q"] == -1).any():
    #     print("FAILURE CORRECTING")
    #     fail_rate = -9999

    # Measure logical qubit
    logical = sc.measure_logical()

    if -1 in logical[0] or -1 in logical[1]:
        fail_rate += 1

    lc.reset()
    sc.reset()

fail_rate = fail_rate/float(iterations)

# Initializing variables. mpi4py requires that we pass numpy objects.
f_rate = np.zeros(1)
total = np.zeros(1)

f_rate[0] = fail_rate

# Communication: root node receives results with a collective "reduce"
comm.Reduce(f_rate, total, op=MPI.SUM, root=0)

# Root process saves the results
if comm.rank == 0:
        total = total/float(size)
        print("size: ", size)
        print("id: ", rank)
        print("TOTAL: ", total)
        args_str = json.dumps(args)
        script_path = dirname(realpath(__file__))
        file_name = (script_path + "/data/" + args_str)
        # np.save(file_name, total)
