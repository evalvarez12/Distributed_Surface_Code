## Disitruted implementation.

Routines into generating GHZ states in a quantum network and making stabilizer measuements
to obtain error models.


### Files:
* circuit_block.py : Individual circuits involved the remote generation of entanglement
and in different steps of purification protocols. Must be used together
with 'circuit.py' to assemble complete circuits.
* circuit.py : Wrapper functions to put together circuits involved
in the process of making GHZ states out of entangled Bell pairs.
Must be used with 'circuit_block.py'
* error_models.py : Error models used to simulate the decoherence effects due to imperfect operations.
* noise_modeling.py : Noise modeling using the Choi-Jamiolkowski isomorphism
to find the chi matrix that caracterizes the process
* protocols.py : Functions with the different protocols for generating a GHZ state.
They rely on circuit.py and circuit_block.py to enssemble the circuit that
generates the GHZ.
* stabilizer.py : Routines for measuring stabilizers in a local system or in a disitributed
one by means of GHZ states.
* result_* : Code for obtaining results/testing the simulations.
