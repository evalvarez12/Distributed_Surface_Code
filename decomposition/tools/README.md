## Tools
Different extra operations required.


### Files:
..* ghz_correction.py : Functions to apply the required corrections when generating GHZ states under
the different numbers of Bell pairs and GHZ sizes.
..* names.py : Functions to generate the names of the data to be used.
The purpose of using this functions is to set a standard name system.
..* operations.py : Routines to perform basic operations not supplied by qutip already.
..* p_success.py : Functions to calculate the success probability of certain protocols.
Relies in 'projectors.py' and 'operations.py' to generate the required projectors.
..* pauli_basis.py : Routines to assemble the pauli basis based on constructing the
operators based on symbols ("IIXI", "YZZX", ... )
..* projectors.py : Some tools to generate projectors based on even/odd or binary numbers.
