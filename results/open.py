"""
Simple simulation to test the surface code simulation is working
without looking at plots.

created-on: 09/12/17
@author: eduardo
"""
import sys
import json
import numpy as np
from os.path import dirname, realpath


def get_file_name(params):
    protocol = "protocol=" + params["protocol"]
    topology = "topology=" + params["topology"]
    distance = "distance=" + str(params["distance"])
    iterations = "iterations=" + str(params["iterations"])
    cycles = "cycles=" + str(params["cycles"])
    eta = "eta=" + str(params["eta"])
    a0 = "a0=" + str(params["a0"])
    a1 = "a1=" + str(params["a1"])
    time = "time=" + str(params["time"])

    param_names = [protocol, topology, distance, iterations, cycles,
                   eta, a0, a1, time]
    file_name = "_".join(param_names)
    return file_name + ".npy"


# Set parameters
parameters = {"topology": "toric",
              "cycles": 50,
              "iterations": 500,
              "protocol": "GHZ",
              "eta": 0.01,
              "a1": 0.0125,
              "time": 0.30347}

for a0 in [2.0]:
    for distance in [8]:

        parameters["distance"] = distance
        parameters["a0"] = a0

        args_str = json.dumps(parameters)
        script_path = dirname(realpath(__file__))
        file_name = get_file_name(parameters)
        fail_rate = np.load(file_name)

        print("Distance: ", distance)
        print("a0: ", a0)
        print("Fail rate:", fail_rate)
        print("--------------------")
