"""
Wrapper functions to put together circuits involved
in the process of making GHZ states out of entangled Bell pairs.
Must be used with 'circuit_block.py'

author: Eduardo Villasenor
created-on: 20/11/17
"""
import error_models as errs
import numpy as np
import qutip as qt
import collections


class Circuit:
    """
    Class to hold a block in the overall GHZ state formation protocol.
    Subcircuit object for the circuit that goes inside this block.

    Circuits are assembled using recursive objects and circuit blocks.
    After each level with a success probability a number of attempts is calculated
    and used to evaluate the time it required unit success. Similar to an event
    oriented simulation.
    """

    def __init__(self, a0, a1, circuit_block, **kwargs):
        """
        Init function.

        Parameters
        -----------
        a0 : (scalar)  environmental error parameter
        a1 : (scalar) environmental error parameter
        cicuit_block : (Blocks) function from circuit_block.py to be evaluated
                       in a step of the circuit
        **kwargs : arguments for evaluating the cicuit_block
        """
        # Save circuit block function of this level
        self.circuit = circuit_block
        self.circuit_kwargs = kwargs

        # Subcircuit for this level
        self.subcircuit = None

        # Environmental dephasing parameters
        self.a0 = a0
        self.a1 = a1

    def run(self, rho):
        """
        Main function to run circuits recurively.

        Parameters
        -----------
        rho : (densmat) density matrix involved in the circuit, can be None depending
              on the circuit block
        """
        return self._run(rho)

    def _run(self, rho, p_parent=1, check_parent=collections.Counter({})):
        # First run self circuit
        p_success, check, rho = self.circuit(rho, **self.circuit_kwargs)

        # If this circuit dependends on the parent, take their probability of
        # success
        p_success *= p_parent
        check += check_parent

        # Now check if it has a subcircuit
        if self.subcircuit:
                _, check, rho = self.subcircuit._run(rho, p_success, check)
        else:
            # If this is the end of the dependency calculate success event
            # starts from 0, where 0 means success on the first try
            n_extra_attempts = self._number_of_attempts(p_success)

            if n_extra_attempts != 0:
                for k in check:
                    check[k] *= n_extra_attempts

                # steps += n_extra_attempts * (steps + steps_parent)

        # NOTE: No dephasing here. it circuit append and run parallel
        # should be used instead
        return 1, check, rho

    def run_parallel(self, rho=None):
        """
        Run circuit two times in parallel, tensoring the resulting states,
        and dephasing the one that was generated first accordingly.
        Cicuits must be self contained events to be able to run in parallel.

        Parameters
        -----------
        rho : (densmat) density matrix involved in the circuit, can be None depending
              on the circuit block
        """
        _, check1, rho1 = self.run(rho)
        _, check2, rho2 = self.run(rho)

        # Only take the check with the longest time
        diff_time = np.abs(check1["time"] - check2["time"])
        if check1["time"] > check2["time"]:
            rho2 = errs.env_error_all(rho2, 0,
                                      self.a1, diff_time)
            check = check1
        else:
            rho1 = errs.env_error_all(rho1, 0,
                                      self.a1, diff_time)
            check = check2
        rho = qt.tensor(rho1, rho2)
        return 1, check, rho

    def append_circuit(self, rho):
        """
        Appended circuit, and dephase accordingly.
        Must be self contained event

        Parameters
        -----------
        rho : (densmat) density matrix involved in the circuit, can be None depending
              on the circuit block
        """
        _, check, rho_app = self.run(None)
        time0 = check["time0"]
        rho = errs.env_error_all(rho, self.a0,
                                 self.a1, time0)
        time1 = check["time1"]
        rho = errs.env_error_all(rho, 0,
                                 self.a1, time1)

        rho = qt.tensor(rho, rho_app)
        return 1, check, rho

    def append_circuit_diff_node(self, rho):
        """
        Appended circuit, and dephase considering the appended ciruit is
        executed entirely in different nodes.
        Must be self contained event.

        Parameters
        -----------
        rho : (densmat) density matrix involved in the circuit, can be None depending
              on the circuit block
        """
        _, check, rho_app = self.run(None)
        time = check["time"]
        rho = errs.env_error_all(rho, 0,
                                 self.a1, time)

        rho = qt.tensor(rho, rho_app)
        return 1, check, rho

    def add_circuit(self, circuit_block, **kwargs):
        """
        Add a circuit to the current circuit chain.

        Parameters
        -----------
        cicuit_block : (Blocks) function from circuit_block.py to be evaluated
                       in a step of the circuit
        **kwargs : arguments for evaluating the cicuit_block
        """
        if not self.subcircuit:
            self.subcircuit = Circuit(self.a0, self.a1, circuit_block, **kwargs)
        else:
            self.subcircuit.add_circuit(circuit_block, **kwargs)

    def _number_of_attempts(self, p_success):
        # Draw a number of attempts according to the Distribution

        # Up to 1000 tries for success
        i = np.arange(1000)
        d = self._distribution(p_success, i)
        return np.random.choice(i, 1, p=d)[0]
        # return 0

    def _distribution(self, p, n):
        # Distribution for the probability in the number of tries of
        # each event
        return p*(1-p)**n
