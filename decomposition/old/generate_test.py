"""
Test routines for generate.py

author: Eduardo Villasenor
created-on: 028/08/17
"""

import generate

generator = generate.Generator()

# Initial parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.0
system_size = 4
parity = "Z"
function = "LOCAL"

# Calculate chi matrix
chi = generator.ask_model(ps=ps, pm=pm, pg=pg, pn=pn, stab_size=system_size,
                          parity=parity, protocol=function)
