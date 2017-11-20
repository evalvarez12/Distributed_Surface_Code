"""
Internal circuits involved in the different steps of purification protocols.

author: Eduardo Villasenor
created-on: 20/11/17
"""
import circuit_blocks as cb
import error_models as errs

class Circuit:
    """
    Class to hold a block in the overall purification protocol.

    Subcircuit object for the circuit that goes inside this block.
    """
    def  __init__(self, circuit_block, **kwargs):
        # Save circuit block of this level, a function to execute of circuit block.
        self.circuit = circuit_block
        self.circuit_kwargs = kwargs

        # Subcircuit for this level
        self.subcircuit = None


    def run(self, rho, p_parent, steps_parent):
        # First run self circuit
        p_success, steps, rho = self.circuit(rho, **self.circuit_kwargs)

        # If this circuit dependends on the parent, take their probability of
        # success
        if p_parent:
            p_success *= p_parent

        # # Count the number of steps so far
        # steps += steps_parent

        # If this is the end of the dependency calculate success event
        if not self.subcircuit_dependant:
            n_extra_attempts = self.success_number_of_tries(p_success) - 1
            if n_extra_attempts != 0:
                steps = n_extra_attempts * steps
                rho = errs.env_dephasing_all(rho, n_steps, True)
                p_success = 0

        # If there is another circuit run it else return
        if self.subcircuit:
            steps_sub, rho = self.subcircuit.run(rho, p_success, steps)
            # Count the steps so far
            steps = steps_sub + steps
            return steps, rho
        else:
            return steps, rho

    def add_circuit(self, dependant, circuit_block, **kwargs):
        if self.subcircuit:
            self.subcircuit.add_circuit(dependant, circuit_block, **kwargs)
        else:
            self.subciruit = Circuit(circuit_block, **kwargs)
            self.subcircuit_dependendant = dependant

    def success_number_of_tries(self, p_success):
        # Up to 20 tries for success
        i = np.arange(100)
        d = self.distribution(p_success, i)
        return np.random.choice(i, 1, p=d)

    def distribution(p, n):
        # Distribution for the probability in the number of tries of
        # each event
        return p*(1-p)**n
