import decompose


generator = decompose.Generator()

# Initial parameters
ps = 0.0
pm = 0.009
pg = 0.009
pn = 0.0
system_size = 4
parity = "X"
function = "LOCAL"

# Calculate chi matrix
model = decompose.Generator()
chi = generator.ask_model(ps=ps, pm=pm, pg=pg, pn=pn, stab_size=system_size, parity="X", protocol=function)
