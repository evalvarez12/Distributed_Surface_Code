"""
Internal circuits involved in the different steps of purification protocols.

author: Eduardo Villasenor
created-on: 20/11/17
"""
import error_models as errs
import numpy as np

class Circuit:
    """
    Class to hold a block in the overall purification protocol.
    Subcircuit object for the circuit that goes inside this block.
    Circuits are assembled using recursive objects and circuit blocks.
    """
    def  __init__(self, circuit_block, **kwargs):
        # Save circuit block of this level, a function to execute of circuit block.
        self.circuit = circuit_block
        self.circuit_kwargs = kwargs

        # Subcircuit for this level
        self.subcircuit = None


    def run(self, rho, p_parent=0, steps_parent=0):
        # First run self circuit
        p_success, steps, rho = self.circuit(rho, **self.circuit_kwargs)

        # If this circuit dependends on the parent, take their probability of
        # success
        # First check if this circuit is linked to its parent
        if p_parent:
            p_success *= p_parent
            steps += steps_parent

        # Now check if it has a subcircuit
        if self.subcircuit:
                _, steps, rho = self.subcircuit.run(rho, p_success, steps)
                return 1, steps, rho
        else:
            # If this is the end of the dependency calculate success event
            # starts from 0, where 0 means success on the first try
            n_extra_attempts = self.success_number_of_tries(p_success)

            if n_extra_attempts != 0:
                steps += n_extra_attempts * (steps + steps_parent)
                if "operation_qubits" in self.circuit_kwargs and p_parent:
                    except_q = self.circuit_kwargs["operation_qubits"]
                    N = len(rho.dims[0])
                    qs = np.arange(N)
                    qs = np.delete(qs, except_q)
                    rho = errs.env_dephasing(rho, steps, True, N, qs)
                else:
                    rho = errs.env_dephasing_all(rho, steps, True)

        return 1, steps, rho



    def add_circuit(self, circuit_block, **kwargs):
        if not self.subcircuit:
            self.subcircuit = Circuit(circuit_block, **kwargs)
        else:
            self.subcircuit.add_circuit(circuit_block, **kwargs)

    def success_number_of_tries(self, p_success):
        # Up to 20 tries for success
        i = np.arange(100)
        d = self.distribution(p_success, i)
        return np.random.choice(i, 1, p=d)[0]

    def distribution(self, p, n):
        # Distribution for the probability in the number of tries of
        # each event
        return p*(1-p)**n
