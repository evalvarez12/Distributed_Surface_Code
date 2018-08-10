
Brief code overview as interpreted by Tim

Last updated: 2018_08_09


The class `Blocks` (in `decomposition/circuit_block.py`) represents compositions of the (distributed) circuits for creating a (noisy) Bell pair. Possible entanglement-creation-protocols are single-click, Barrett-Kok and EPL. It also contains a few functions that correspond to operations on these Bell pairs, such as performing a "single selection round" or "collapsing the ancillas in the nodes to create a GHZ state".

The class `Circuit` (in `decomposition/circuit.py`) represents compositions of general circuits (so not only for Bell pair generation and operations on it). It can hold either `Block` instances or `Circuit` instances itself. This `Circuit` object can be `run`, which results in a final state as output of the circuit.

The file `decomposition/error_models.py` contains a whole bunch of functions which take a state :math:`\rho` and returns :math:`\mathcal{N}(\rho)`, where :math:`\mathcal{N}` is a particular (error) map. Also includes measurements.

The file `decomposition/protocols.py` contains a number of functions, each of which outputs a `Circuit object`. These functions are `pair_EPL`, `pair_single_sel`, `ghz4_epl`, etc. The Circuit object they output is the circuit that generates an EPL state, a state outputted by the single-selection-distillation-protocol, etc. Each function only takes error parameters as input.

The class `Stabilizer` (in `decomposition/stabilizer.py`) contains functions used for the production of particular states without having to go through the `Circuit` objects, as in `decomposition/protocols.py`. It also has function for processing of the states, such as applying gates and performing measurements. I presume these functions are mainly meant for operations on GHZ-states, which are used for the stabilizer measurements; hence the class's name.

Several standard parameters values are defined in the directory `decomposition/tools`, but also some operations such as measurements and gates (e.g. `operations.py`).
 

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





