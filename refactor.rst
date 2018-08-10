2018_08_10

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
