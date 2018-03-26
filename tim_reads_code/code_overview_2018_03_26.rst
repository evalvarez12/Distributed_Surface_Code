2018_03_02

Note that many classes take functions as input parameters.

Code overview of folder `decomposition`
---------------------------------------

.. admonition:: README of the folder:

   ## Disitruted implementation.
   
   Routines into generating GHZ states in a quantum network and making stabilizer measuements
   to obtain error models.
   
   
   ### Files:
   
  * circuit_block.py : Individual circuits involved the remote generation of entanglement    and in different steps of purification protocols. Must be used together    with 'circuit.py' to assemble complete circuits.
   
  * circuit.py : Wrapper functions to put together circuits involved    in the process of making GHZ states out of entangled Bell pairs.    Must be used with 'circuit_block.py'
   
  * error_models.py : Error models used to simulate the decoherence effects due to imperfect operations.
   
  * noise_modeling.py : Noise modeling using the Choi-Jamiolkowski isomorphism    to find the chi matrix that caracterizes the process
   
  * protocols.py : Functions with the different protocols for generating a GHZ state.    They rely on circuit.py and circuit_block.py to enssemble the circuit that    generates the GHZ.
   
  * stabilizer.py : Routines for measuring stabilizers in a local system or in a disitributed    one by means of GHZ states.
   
  * result_* : Code for obtaining results/testing the simulations.

.. figure:: dependency_tree.jpg
   :width: 1000

   Dependency tree of the folder `fault_tolerance/decomposition`. `For better resolution, click here <./dependency_tree.pdf>`_




**class `Circuit`**

(circuit.py)

  * parameters:
        * "environmental" error parameters :math:`a_0, a_1`
        * circuit block: this is a method from the class `Block`
        * more arguments, which are then used for evaluating the circuit block
  * contains bunch of useful functions to:
        * take tensor product of circuit `self` with circuit `other`
        * append a circuit
        * concatenate circuit `self` with circuit `other`
  * depends only on `error_models.py`

**class `Block`**

(circuit_block.py)
 
  * parameters: a bunch of error parameters
  * contains a large bunch of private functions, named like `_generate_bell_single_click` that return the tuple (time, state), where `state` is e.g. a Bell state and `time` is the time it takes to generate
  * there are a lot of parameters in here which are set to a specific value. Maybe worth checking that these are correct.
  * depends on `error_models.py` and some functions from `tools`


**file `error_models.py`**

  * contains functions that takes a state :math:`\rho` and returns :math:`\mathcal{N}(\rho)`, where :math:`\mathcal{N}` is a particular (error) map. Also includes measurements.
  * depends only on `tools.operations`

**class Generator**

(generate.py)
  
  * DEPRECATED?!
  * depends on `noise_modeling` and `protocols_tools`

**class NoiseModel**

(noise_modeling.py)

  * contains the :math:`\chi`-matrix that characterizes a noise map (=quantum channel)
  * parameters:
      * system size
      * parity (i.e. "X or Z parity of the stabilizer measurement")
  * depends on some functions from `tools`
  * TODO: FINISH CHECKING THIS CLASS

**file `protocols.py`**
  * the docstring:"""Functions with the different protocols for generating a GHZ state. They rely on circuit.py and circuit_block.py to enssemble the circuit that generates the GHZ. Each function returns a circuit object that must be executed to obtain a GHZ."""
  * each function only takes error parameters as input
  * methods in this file are: `EPL_4`, `BK_4`, `BK_3`, et cetera.
  * depends on `circuit_block` and `circuit`


**class `Stabilizer`**
  * contains functions used for the processing of the GHZ states  (the name "Stabilizer" is not such a good name)
  * depends on `error_models` and `tools`
  * TO DO: FINISH READING

**directory `tools`**
  * STILL TO READ

**_A list of files for running the simulations and generating the statistics:_**

**file `result_entanglement_generation.py**
  * depends on `circuit` and `circuit_block`

**file `result_extract_chi.py**
  * depends on `noise_modeling`, `stabilizer` and `tools.names`

**file `result_fidelity_to_chimatrix.py**
  * depends on `noise_modeling`, `error_models`, `stabilizer` and `tools.names`

**file `results_protocols.py**
  * depends on `stabilizer`, `protocols`, `tools.names`




Questions ("towards finding bugs")
----------------------------------
circuit.py: line 109, 105, 130: why is the error parameter :math:`a_0` here set to zero?

circuit_block.py: maybe check parameter values

noise_modeling.py:
        * the example in `_chi_reduce_permutations` is not a valid equality? What should this be?


Code questions/comments
-----------------------
  * what is the difference between `append_circuit` and `run_parallel`?
  * the method `run_parallel` still talks about dephasing in the docstring....
  * what does "depolarize accordingly" mean in the docstring of method `append_circuit`?
  * class `Generator` does not have a docstring (although the file itself does)
  * noise_modeling.py: 
        * line 35: does the symbol * perform matrix multiplication?


Code overview of home folder `fault_tolerance`
----------------------------------------------

.. admonition:: README of the folder:

		# Fault tolerant surface code simulations.
		
		Surface code simulations for a distributed implementation using considering NV centers
		for the quantum nodes.
		
		The simulations depend first on the creation of GHZ states from which a error model can be obtained see *decomposition/*.
		The simulations rely on Kolmogorov's implementation of the Blossom algorithm : Blossom V for the decoder.
		
		
		### Files:
		  * surface_code.py : Main class for the surface code. Contains all relevant functions required for measuring noisy stabilizers and applying the corrections given by the decoder.
		  * layers.py : Class for embedding a surface code object, that allows for recording several continuous stabilizer measurements as required when doing imperfect measurements.
		  * matching.py : Functions for interacting with the decoder. This included transforming the measurements into the graphs required.
		  * `result_`  * : Code for obtaining results/testing the simulations.




**class `SurfaceCode`**

(surface_code.py)

  * parameters:
      * *distance*: (int) ???
      * *surface*: string "toric" or "planar"
  * depends on: *errors*
  * in the *init* function: 
        * if we have the toric code (with *d* the distance):
            * the number of data qubits equals :math:`2 d^2`
            * the number of stabilizers is :math:`d^2`
            * the *side* equals :math:`2d`
        * if we have the planar code:
            * the number of data qubits equals :math:`d^2+ (d-1)^2`
            * the number of stabilizers is :math:`(d-1)*d`
            * the *side* equals :math:`2d-1`
                
     where I presume that *side* refers to 

        * the variable *ind1* (*ind2*) contains all odd (even) numbers between 1 and *side*.
        * the variable *qubits* is a tuple of two matrices of size *size x size*. The first contains the Z errors, the latter contains the X errors. Both are initialized with 1's in every entry.
        * for the remainder of the variables, see the jupyter notebook `Disecting_file_surface_code`

   * method `init_error_obj`: takes a bunch of error parameters, whether planar or toric surface code is used, and which protocol is applied, and stores the result of `errors.Generator` (these parameters).
     
   * method `_stabilizer_qubits_bulk` obtains, given a position on the grid, the four positions adjacent to this position.

   * method `_measure_stabilizer_bulk` then computes, given a position on the grid, the product of the measurement outcomes of the four adjacent nodes, where the measurement outcome is given (I think) in `self.qubits`. Subsequently, set the value of self.qubits[0] to this product.

   * IN THIS FILE, I FIND IT VERY UNCLEAR WHICH MATRIX REFERS TO WHAT, AND WHAT IS MEASURED PRECISELY WHEN!!!




**class `Layers`**

(layers.py)

  * looks like a class that only has been created to ease administration
  * depends on: *matching*
  * docstring:*"""Layers to store the syndrome measurements in the surface codes. Works as a 3d surface code where time is the third dimension."""*
  * only single parameter: *surface_code* is a SurfaceCode object
  * Attributes:
      * *number_stabs*: number of stabilizers, obtained from *surface_code*
      * past_syndrome_star and past_syndrome_plaq: this contains the **previous** initialized to a vector of ones of length *number_stabs*
      * *syndromes_star*/*syndromes_plaq*: lists in which the different "layers" in time will be stored 
      * *stars_pos*/*plaqs_pos*

   * selection of methods:
      * *decode* TO DO
      * *add*: check the previous layer (*past_syndrome_star*) for all the star positions (plaquette positions)
      * *find_anyons*: 

Bugs?
-----
  * in surface_code.py:
      * line 89: after this command, the variable `self.qubits` is a tuple of two :math:`size\times size` matrices, which take the form

		>>>  # this is supposed to "mark the Z errors"
		>>>  [[[1. 1. 1. 1. 1. 1. 1.]
		>>>  [1. 1. 1. 1. 1. 1. 1.]
		>>>  [1. 1. 1. 1. 1. 1. 1.]
		>>>  [1. 1. 1. 1. 1. 1. 1.]
		>>>  [1. 1. 1. 1. 1. 1. 1.]
		>>>  [1. 1. 1. 1. 1. 1. 1.]
		>>>  [1. 1. 1. 1. 1. 1. 1.]]

		>>>  # this is supposed to "mark the X errors"
		>>> [[1. 0. 1. 0. 1. 0. 1.]
		>>>  [0. 1. 0. 1. 0. 1. 0.]
		>>>  [1. 0. 1. 0. 1. 0. 1.]
		>>>  [0. 1. 0. 1. 0. 1. 0.]
		>>>  [1. 0. 1. 0. 1. 0. 1.]
		>>>  [0. 1. 0. 1. 0. 1. 0.]
		>>>  [1. 0. 1. 0. 1. 0. 1.]]]

        What does this mean? That there are already a lot of X-errors but only few Z-errors?

        * line 93-100: the variable self.plane becomes

		>>> [['-' 't' '-' 't' '-' 't' '-']
		>>> ['l' '-' 'o' '-' 'o' '-' 'r']
		>>> ['-' 'o' '-' 'o' '-' 'o' '-']
		>>> ['l' '-' 'o' '-' 'o' '-' 'r']
		>>> ['-' 'o' '-' 'o' '-' 'o' '-']
		>>> ['l' '-' 'o' '-' 'o' '-' 'r']
		>>> ['-' 'b' '-' 'b' '-' 'b' '-']]

           I presume t=top, l=left, r=right, b=bottom. What do the 'o' and the '-' refer to?

       * line 172-3: now what happens if *p_not_complete=0*? And in the case of :math:`p_not_complete\neq 0`, then is another perfect measurement performed after the incomplete measurement? Which commands are here executed in which order?

       * line 187, 204: the matrix self.qubits contains zeroes and ones. Does this mean that *vals* =1 if and only if there are just ones in the list and no zeroes?

       * line 189, 208: don't you need a separate list to store the new values? Now, in order to compute some of the new values, you're already using the new values of adjacent nodes if they were updated earlier!

       READ ON FROM `_stabilizer_qubits_boundary`

  * in layers.py:
      * line 115: does the index [0,0] now always also end up as an error?

Question on the code
--------------------
  * in general: when documenting a method, provide: which parameters does it take, what does it return, and what does it do in the meantime?
  * in surface_code.py: 
      * lines 77-78: why first y and x second?
      * line 88: where is the '9' in the code then?
      * when defining *self.qubits*: for readability, maybe better to make a class out of this, with attributes *xerrors* and *zerrors*. Will they only contain zeroes and ones? What do these entries mean?
      * line 140: maybe best to change name of `stabilizer` into *stabilizer_type* or *XorZ* or something like this.
      * line 173, 251: typo: measuerement
      * line 206: why is the parameter axis needed here and what does it do in this context?
      * line 191 and many other lines: by just looking at it, it is not clear to the reader what the variable `c` corresponds to. Maybe better to have a globally defined list `[index_planar, index_surface] = [0,1]` and call `c` `surface_code_type` instead.
  * layers.py: 
      * generally well documented!
      * maybe rename *past_syndrome_star* to *previous_syndrome_star*
      * line 74, method *add*: what do you mean by *current measurement status*?
      * line 97: best not to call "plaq" or "star" a *stabilizer*
      * line 77: what do you mean by "physical errors are only carried once" and how is this implemented in the code?





To write up in thesis: MWPM, anyons-stuff, decoding procedure
