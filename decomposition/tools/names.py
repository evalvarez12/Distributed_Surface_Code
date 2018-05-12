"""
Functions to generate the names of the data to be used.
The purpose of using this functions is to set a standard name system.

author: Eduardo Villasenor
created-on: 20/01/18
"""

from os.path import dirname, realpath


def chi(ps, pm, pg, eta, a0, a1, theta, stab_size, parity, protocol):
    """Name generator for CHI matrix file outputs."""
    param_names = ["ps=" + str(round(ps, 5)), "pm=" + str(round(pm, 5)),
                   "pg=" + str(round(pg, 5)), "eta=" + str(round(eta, 5)),
                   "a0=" + str(round(a0, 5)), "a1=" + str(round(a1, 5)),
                   "theta=" + str(round(theta, 5))]

    param_names = "_".join(param_names)
    file_name = ["CHI", protocol, parity, str(stab_size)]
    file_name = "_".join(file_name)
    # The address of the parent parent parent (3 levels) directory
    script_path = dirname(dirname(dirname(realpath(__file__))))
    file_name = (script_path + "/data/" + file_name
                 + "_" + param_names + ".dict")
    return file_name


def ghz(ps, pm, pg, eta, a0, a1, theta, size, protocol):
    """Name generator for saved GHZ states."""
    param_names = ["ps=" + str(round(ps, 5)), "pm=" + str(round(pm, 5)),
                   "pg=" + str(round(pg, 5)), "eta=" + str(round(eta, 5)),
                   "a0=" + str(round(a0, 5)), "a1=" + str(round(a1, 5)),
                   "theta=" + str(round(theta, 5))]

    param_names = "_".join(param_names)
    file_name = ["GHZ", protocol, str(size)]
    file_name = "_".join(file_name)
    # The address of the parent parent directory
    script_path = dirname(dirname(realpath(__file__)))
    file_name = (script_path + "/data/" + file_name
                 + "_" + param_names)
    return file_name

def ghz_times(ps, pm, pg, eta, a0, a1, theta, size, protocol):
    """Name generator for saved GHZ states."""
    param_names = ["ps=" + str(round(ps, 5)), "pm=" + str(round(pm, 5)),
                   "pg=" + str(round(pg, 5)), "eta=" + str(round(eta, 5)),
                   "a0=" + str(round(a0, 5)), "a1=" + str(round(a1, 5)),
                   "theta=" + str(round(theta, 5))]

    param_names = "_".join(param_names)
    file_name = ["TIMES", protocol, str(size)]
    file_name = "_".join(file_name)
    # The address of the parent parent directory
    script_path = dirname(dirname(realpath(__file__)))
    file_name = (script_path + "/data/" + file_name
                 + "_" + param_names)
    return file_name + ".npy"
