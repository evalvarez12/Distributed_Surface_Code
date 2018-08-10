
Brief code overview as interpreted by Tim

Last updated: 2018_08_10


The class `SurfaceCode` (in `surface_code.py`) represents the set of qubits involved in an implementation of the surface code and has methods such as `_measure_stabilizer_bulk` and `apply_qubit_error`.

A `Layers` object (in `layers.py`) is a collection of "sheets of measurements", each of which is performed at a different instance in time. Dealing with anyons and decoding happens in this class, by importing functions from `matching.py`, which contains functions that basically wrap around the implementation of Edmond's blossom algorithm is the `blossom5` folder.

The class `Generator` (in `errors.py`) reads an error model in, stored in the folder `config_files`.

**Folder `decomposition`**

The class `Blocks` (in `decomposition/circuit_block.py`) represents compositions of the (distributed) circuits for creating a (noisy) Bell pair. Possible entanglement-creation-protocols are single-click, Barrett-Kok and EPL. It also contains a few functions that correspond to operations on these Bell pairs, such as performing a "single selection round" or "collapsing the ancillas in the nodes to create a GHZ state".

The class `Circuit` (in `decomposition/circuit.py`) represents compositions of general circuits (so not only for Bell pair generation and operations on it). It can hold either `Block` instances or `Circuit` instances itself. This `Circuit` object can be `run`, which results in a final state as output of the circuit.

The file `decomposition/error_models.py` contains a whole bunch of functions which take a state :math:`\rho` and returns :math:`\mathcal{N}(\rho)`, where :math:`\mathcal{N}` is a particular (error) map. Also includes measurements.

The file `decomposition/protocols.py` contains a number of functions, each of which outputs a `Circuit object`. These functions are `pair_EPL`, `pair_single_sel`, `ghz4_epl`, etc. The Circuit object they output is the circuit that generates an EPL state, a state outputted by the single-selection-distillation-protocol, etc. Each function only takes error parameters as input.

The class `Stabilizer` (in `decomposition/stabilizer.py`) contains functions used for the production of particular states without having to go through the `Circuit` objects, as in `decomposition/protocols.py`. It also has function for processing of the states, such as applying gates and performing measurements. I presume these functions are mainly meant for operations on GHZ-states, which are used for the stabilizer measurements; hence the class's name.

Several standard parameters values are defined in the directory `decomposition/tools`, but also some operations such as measurements and gates (e.g. `operations.py`).
 
**Folder `blossom5`**

This folder contains Vladimir Kolmogorov's implementation of Edmonds blossom algorithm that is used in the minimum-weight-perfect-matching decoding of the surface code.


**Folder `config_files`**

Contains a large number noisy channels, represented as Chi-matrices, each for a different set of noise parameters.

.. figure:: ../../faulttolerance/pics/dependency_tree.jpg
   :width: 1000

   Dependency tree of the folder `faulttolerance/decomposition`. `For better resolution, click here <./../../faulttolerance/dependency_tree.pdf>`_


Recommendations for refactoring
-------------------------------

**file blocks and circuit**

  * make a class for Blocks/Circuits that output a bell pair, with as initialization parameter `generation_method`, which takes a string from ["single_click", "Barrett-Kok", "EPL"]
  * all the operations in the `Block` class are quite vague; it looks like a lot of refactoring, subclassing et cetera could definitely help here. Also, there are many of them, which makes testing a whole lot of work.



**file `error_models.py`**

  * the functions `bell_pair` (which outputs a noisy Phi+) and `bell_pair_psi` (which outputs a noisy Psi+) can be put together in a single function
  * (optional) create a class `State` (inherit from `qutip.Quobj`?) which has functions that, when called, apply the noise to `self`


**file `protocols.py`**

  * each of these functions outputs a Circuit object. Each of these functions could be rewritten to a subclass of `Circuit`.
  * it seems like `protocols.py` and `stabilizer.py` have functions that perform a similar function. Can the two be merged?


**file `stabilizer.py`**:

  *  The name "Stabilizer" is maybe not the most intuitive name for this class; maybe `GHZoperations` or something like that better covers it purpose?


**file `surface_code.py`**:

  * the current code is a bit hard to test
  * the class has almost 40 methods, which is quite a lot
  * it would be a lot simpler if there were a separate base class `Lattice`, which has functions such as `get_neighbours(pos)` where `pos` is a tuple (x,y)-coordinates. The class `SurfaceCode` then either inherits from `Lattice` or holds one as an attribute. This prevents that you need, like now, functions such as `_measure_stabilizer_side`, `_measure_stabilizer_qubits_bulk`, et cetera. A function like `_stabilizer_qubits_boundary` would then move to the `Lattice` class
    * rather than using `c=0` for "stars" and `c=1` for "plaquettes", readability would be improved when one makes a class variable `_STARS = 0` and `_PLAQS = 1`.


**file `layers.py`**:

  * what about making a class `Layer` or `MeasurementSheet`, which contains a single layer of measurements? The current class `Layers` (which we would then rename to `LayerList` or something like it) would be a collection of such `Layer` object, possibly with some additional methods in the way it is now.
