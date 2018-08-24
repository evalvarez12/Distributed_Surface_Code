"""
This module contains functions that we can use to inspect the rest of
the codebase. It's not a test module, so we don't need to import other
modules here, but we can import it in notebooks or interpreters to
learn more about what's under the hood.

    import sys
    sys.path.append('../')
    from analysis import *

"""


#-------Variables that we'll probably want to use more than once------#

from qutip import sigmax, sigmay, sigmaz
X, Y, Z = sigmax(), sigmay(), sigmaz()
I = X * X

#---------------------------------------------------------------------#

#---------------------Shortcuts to External Functions-----------------#
"""
I don't want to mess up the namespace by importing * from qutip, for
example, so I'll import * from analysis instead.
"""
from qutip import tensor

#---------------------------------------------------------------------#

#------------------------Analysis Functions---------------------------#
def properties(obj):
    """
    Returns a set of names which correspond to properties of an object:
    these names can be used to access parts of the object, but not
    methods, class methods, etc.

    For example, if you used this function on a Point object with
    internal `x`, `y`, and `z`, you'd get back `{'x', 'y', 'z'}`, but 
    nothing about `__add__`.
    """
    return set(dir(obj)) - set(dir(type(obj)))

def internals(obj):
    """
    For convenience, gets the values associated with 'properties' and 
    returns them in a dict.
    """
    return {_: obj.__dict__[_] for _ in properties(obj)}
#---------------------------------------------------------------------#