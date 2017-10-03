import noise_modeling
import protocols
import pickle

# # Set error values
# ps = 0.0
# pm = 0.0075
# pg = 0.0075
# pn = 0.1
#
# # Initial parameters
# system_size = 3
# parity = "X"
#
# # Initialize objects
# prot = protocols.Protocols(ps, pm, pg, pn)
# # Model takes the superoperator as a parameter
# model = noise_modeling.NoiseModel(system_size, prot.stringent, parity)
#
# # Calculate chi matrix
# model.separate_basis_parity()
# model.make_chi_matrix()

class ErrorWrapper:
    def __init__(self):
        # Initialize function names strings
        self.prot = protocols.Protocols(0, 0, 0, 0)

        # Dictionary of functions on the protocols
        self.funcs = {"STRINGENT": self.prot.stringent,
                      "EXPEDIENT": self.prot.expedient,
                      "LOCAL": self.prot.local_stabilizer}

    def generate_model(self, ps, pm, pg, pn, stab_size, parity, function):
        if function not in self.funcs:
            raise ValueError("Function not in protocols.py")

        self.prot.change_parameters(ps, pm, pg, pn)

        # Model takes the superoperator as a parameter
        model = noise_modeling.NoiseModel(stab_size, self.funcs[function], parity)

        # Calculate chi matrix
        model.separate_basis_parity()
        model.make_chi_matrix()

        file_name = self.generate_name(ps, pm, pg, pn, stab_size, parity, function)

        pickle_out = open(file_name, "wb")
        pickle.dump(model.chi, pickle_out)
        pickle_out.close()

    def generate_name(self, ps, pm, pg, pn, stab_size, parity, function):
        param_names = ["ps=" + str(ps), "pm=" + str(pm), "pg=" + str(pg), "pn=" + str(pn)]
        param_names = "_".join(param_names)
        file_name = [function, parity, str(stab_size)]
        file_name = "_".join(file_name)
        file_name = "data/" + file_name + "_" + param_names + ".dict"
        return file_name

    def load_model(self, ps, pm, pg, pn, stab_size, parity, function):
        file_name = self.generate_name(ps, pm, pg, pn, stab_size, parity, function)
        pickle_in = open(file_name, "rb")
        return pickle.load(pickle_in)
