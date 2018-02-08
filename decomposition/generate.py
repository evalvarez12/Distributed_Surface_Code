""""
Wrapper object to feed error models into the surface code from the noisy circuits
making the GHZ states.

NOTE: Due to changes in the way the GHZ states are calculated this routines are no longer
utilizable.

author: Eduardo Villasenor
created-on: 028/08/17
"""

import pickle
from os.path import dirname, realpath
# from . import noise_modeling
# from . import protocols
import noise_modeling
import protocols_tools

class Generator:
    def __init__(self):
        # Initialize obeject protocols
        self.prot = protocols_det.ProtocolsDeterministic(0, 0, 0, 0)

        # Initialize protocol names strings dictionary
        self.funcs = {"STRINGENT": self.prot.stringent,
                      "EXPEDIENT": self.prot.expedient,
                      "LOCAL": self.prot.local_stabilizer}

    def generate_model(self, ps, pm, pg, pn, stab_size, parity, protocol):
        if protocol not in self.funcs:
            raise ValueError("protocol not in protocols.py")

        self.prot.change_parameters(ps, pm, pg, pn)

        # Model takes the superoperator as a parameter
        noise_model = noise_modeling.NoiseModel(stab_size, parity)

        # Calculate chi matrix
        noise_model.apply_superoperator(self.funcs[protocol])
        noise_model.separate_basis_parity()
        noise_model.make_chi_matrix()

        file_name = self.generate_name(ps, pm, pg, pn,
                                       stab_size, parity, protocol)

        pickle_out = open(file_name, "wb")
        pickle.dump(noise_model.chi, pickle_out, protocol=2)
        pickle_out.close()

    def generate_name(self, ps, pm, pg, pn, stab_size, parity, protocol):
        param_names = ["ps=" + str(ps), "pm=" + str(pm),
                       "pg=" + str(pg), "pn=" + str(pn)]
        param_names = "_".join(param_names)
        file_name = [protocol, parity, str(stab_size)]
        file_name = "_".join(file_name)
        # The address of the parent directory
        script_path = dirname(dirname(realpath(__file__)))
        file_name = (script_path + "/data/" + file_name
                     + "_" + param_names + ".dict")
        return file_name

    def load_model(self, ps, pm, pg, pn, stab_size, parity, protocol):
        file_name = self.generate_name(ps, pm, pg, pn,
                                       stab_size, parity, protocol)
        try:
            pickle_in = open(file_name, "rb")
            return pickle.load(pickle_in)
        except:
            return None

    def ask_model(self, ps, pm, pg, pn, stab_size, parity, protocol):
        model = self.load_model(ps, pm, pg, pn,
                                stab_size, parity, protocol)
        if model:
            return model
        else:
            self.generate_model(ps, pm, pg, pn,
                                stab_size, parity, protocol)
            model = self.load_model(ps, pm, pg, pn,
                                    stab_size, parity, protocol)
            return model
