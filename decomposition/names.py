"""
Functions to generate the names of the data to be used.

author: Eduardo Villasenor
created-on: 20/01/18
"""

from os.path import dirname, realpath


def chi(ps, pm, pg, eta, a0, a1, theta, stab_size, parity, protocol):
    """Name for CHI matrix file."""
    param_names = ["ps=" + str(ps), "pm=" + str(pm),
                   "pg=" + str(pg), "eta=" + str(round(eta, 4)),
                   "a0=" + str(round(a0, 4)), "a1=" + str(round(a1, 4)),
                   "theta=" + str(theta)]

    param_names = "_".join(param_names)
    file_name = ["CHI", protocol, parity, str(stab_size)]
    file_name = "_".join(file_name)
    # The address of the parent directory
    script_path = dirname(dirname(realpath(__file__)))
    file_name = (script_path + "/data/" + file_name
                 + "_" + param_names + ".dict")
    return file_name


def ghz(ps, pm, pg, eta, a0, a1, theta, size, protocol):
    """Name for GHZ data file."""
    param_names = ["ps=" + str(ps), "pm=" + str(pm),
                   "pg=" + str(pg), "eta=" + str(round(eta, 4)),
                   "a0=" + str(round(a0, 4)), "a1=" + str(round(a1, 4)),
                   "theta=" + str(theta)]

    param_names = "_".join(param_names)
    file_name = ["GHZ", protocol, str(size)]
    file_name = "_".join(file_name)
    # The address of the parent directory
    script_path = dirname(dirname(realpath(__file__)))
    file_name = (script_path + "/data/" + file_name
                 + "_" + param_names + ".dict")
    return file_name
