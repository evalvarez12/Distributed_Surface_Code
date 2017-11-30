"""
Internal circuits involved in the different steps of purification protocols.

author: Eduardo Villasenor
created-on: 20/11/17
"""
import error_models as errs
import numpy as np
import qutip as qt
import collections

class Circuit:
    """
    Class to hold a block in the overall purification protocol.
    Subcircuit object for the circuit that goes inside this block.
    Circuits are assembled using recursive objects and circuit blocks.
    """
    def  __init__(self, p_env, circuit_block, **kwargs):
        # Save circuit block of this level, a function to execute of circuit block.
        self.circuit = circuit_block
        self.circuit_kwargs = kwargs

        # Subcircuit for this level
        self.subcircuit = None

        self.p_env = p_env

    def run(self, rho, p_parent=1, check_parent=collections.Counter({})):
        # First run self circuit
        p_success, check, rho = self.circuit(rho, **self.circuit_kwargs)

        # If this circuit dependends on the parent, take their probability of
        # success
        p_success *= p_parent
        check += check_parent

        # Now check if it has a subcircuit
        if self.subcircuit:
                _, check, rho = self.subcircuit.run(rho, p_success, check)
        else:
            # If this is the end of the dependency calculate success event
            # starts from 0, where 0 means success on the first try
            n_extra_attempts = self.success_number_of_tries(p_success)

            if n_extra_attempts != 0:
                for k in check:
                    check[k] *= n_extra_attempts

                # steps += n_extra_attempts * (steps + steps_parent)

                # TODO remove dephasing here. it circuit append and run parallel should be used instead
                # if "operation_qubits" in self.circuit_kwargs:
                #     except_q = self.circuit_kwargs["operation_qubits"]
                #     N = len(rho.dims[0])
                #     qs = np.arange(N)
                #     qs = np.delete(qs, except_q)
                #
                #     rho = errs.env_dephasing(rho, self.p_env, steps, True, N, qs)
                # else:
                #     rho = errs.env_dephasing_all(rho, self.p_env, steps, True)
            # print("End of chain")
            # print("p_success:", p_success)
            # print("steps: ", steps)

        return 1, check, rho

    def run_parallel(self, rho=None):
        """
        Run circuit two times in parallel, tensoring the resulting states,
        and dephasing the one that was generated first accordingly.
        Cicuits must be self contained events to be able to run in parallel.
        """
        _, check1, rho1 = self.run(rho)
        _, check2, rho2 = self.run(rho)

        # Only take the check with the longest time
        diff_time = np.abs(check1["time"] - check2["time"])
        if check1["time"] > check2["time"]:
            rho2 = errs.env_dephasing_all(rho2, self.p_env, diff_time, True)
            check = check1
        else:
            rho1 = errs.env_dephasing_all(rho1, self.p_env, diff_time, True)
            check = check2
        rho = qt.tensor(rho1, rho2)
        return 1, check, rho

    def append_circuit(self, rho):
        """
        Appended circuit, and dephase accordingly.
        Must be self contained event
        """
        _, check, rho_app = self.run(None)
        time = check["time"]
        rho = errs.env_dephasing_all(rho, self.p_env, time, True)
        rho = qt.tensor(rho, rho_app)
        return 1, check, rho

    def add_circuit(self, circuit_block, **kwargs):
        """
        Add a circuit to the current chain.
        """
        if not self.subcircuit:
            self.subcircuit = Circuit(self.p_env, circuit_block, **kwargs)
        else:
            self.subcircuit.add_circuit(circuit_block, **kwargs)

    def success_number_of_tries(self, p_success):
        """
        Draw a number of attempts according to the Distribution
        """
        # Up to 20 tries for success
        i = np.arange(100)
        d = self._distribution(p_success, i)
        # return np.random.choice(i, 1, p=d)[0]
        return 0

    def _distribution(self, p, n):
        # Distribution for the probability in the number of tries of
        # each event
        return p*(1-p)**n
