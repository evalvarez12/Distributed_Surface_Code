"""
This module contains functions that we can use to inspect the rest of
the codebase. It's not a test module, so we don't need to import other
modules here, but we can import it in notebooks or interpreters to
learn more about what's under the hood.
"""

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