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

    def run(self):
        """
        Main function to run circuits recurively.

        Parameters
        -----------
        rho : (densmat) density matrix involved in the circuit, can be None depending
              on the circuit block
        """
        check, rho = self._run()
        return rho, check

    def _run(self):
        check = self._empty_check()
        rho = None
        p_success = 0
        r = np.random.rand()
        while  r > p_success:
            # Check if it has a subcircuit
            if self.subcircuit:
                c, rho = self.subcircuit._run()
                check += c

            # Run self circuit keeping track of the time
            p_success, c, rho = self.circuit(rho, **self.circuit_kwargs)
            check += c
            # print(check)
            r = np.random.rand()

        return check, rho

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
        check1, rho1 = self._run()
        check2, rho2 = self._run()

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

    def append_circuit(self, rho=None):
        """
        Appended circuit, and depolarize accordingly.
        Must be self contained event

        Parameters
        -----------
        rho : (densmat) density matrix involved in the circuit, can be None depending
              on the circuit block
        """
        check, rho_app = self._run()
        time0 = check["time0"]
        rho = errs.env_error_all(rho, self.a0,
                                 self.a1, time0)
        time1 = check["time1"]
        rho = errs.env_error_all(rho, 0,
                                 self.a1, time1)

        rho = qt.tensor(rho, rho_app)
        return 1, check, rho

    def start(self, rho=None):
        """
        Star the first element of the circuit.
        Must be self contained event

        Parameters
        -----------
        rho : (densmat) density matrix involved in the circuit, can be None depending
              on the circuit block
        """
        check, rho = self._run()
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

    def _empty_check(self):
        return collections.Counter({"bell_pair": 0,
                                    "two_qubit_gate": 0,
                                    "single_qubit_gate": 0,
                                    "measurement": 0,
                                    "time0": 0,
                                    "time1": 0,
                                    "time": 0})
