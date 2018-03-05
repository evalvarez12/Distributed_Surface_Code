# Fault tolerant surface code simulations.

Surface code simulations for a distributed implementation using considering NV centers
for the quantum nodes.

The simulations depend first on the creation of GHZ states from which a error model
can be obtained see *decomposition/*.
The simulations rely on Kolmogorov's implementation of the Blossom algorithm : Blossom V
for the decoder.


### Files:
* surface_code.py : Main class for the surface code. Contains all relevant functions
required for measuring noisy stabilizers and applying the corrections given by the decoder.
* layers.py : Class for embedding a surface code object, that allows for recording several
continuous stabilizer measurements as required when doing imperfect measurements.
* matching.py : Functions for interacting with the decoder. This included transforming the measurements into the graphs required.
* result_* : Code for obtaining results/testing the simulations.
