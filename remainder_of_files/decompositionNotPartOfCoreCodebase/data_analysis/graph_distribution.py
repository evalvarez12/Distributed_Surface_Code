from scipy import stats
import numpy as np
import matplotlib.pylab as plt
import tools.names as names
import qutip as qt

# Parameters
ps = 0.003
pm = 0.003
pg = 0.003
a0 = 5.0
a1 = 1/30.
eta = 1/100.
theta = .63

ghz_size = 4
protocol = "thres_eta"
ghz_ref = qt.ghz_state(4) * qt.ghz_state(4).dag()

ghz_file = names.ghz(ps, pm, pg, eta, a0, a1, theta,
                     ghz_size, protocol)
times_file = names.ghz_times(ps, pm, pg, eta, a0, a1, theta,
                             ghz_size, protocol)



ghzs = qt.qload(ghz_file)
times = np.load(times_file)
print("N: ", len(times))

fidelities = []
for i in ghzs:
    fidelities += [qt.fidelity(i, ghz_ref)]
fidelities = np.array(fidelities)

indices_sorted = np.argsort(times)
times = times[indices_sorted][:-1]
fidelities = fidelities[indices_sorted][:-1]

ignore_number = int(len(times)/100*5.)
t_max = times[:-ignore_number][-1]

fig, ax1 = plt.subplots()
left, bottom, width, height = [0.52, 0.50, 0.43, 0.43]
ax2 = fig.add_axes([left, bottom, width, height])


ax2.plot(times, fidelities, 'k.')

vline = np.linspace(min(fidelities), max(fidelities))
o = np.ones_like(vline)*(t_max)
ax2.plot(o, vline, 'r--')

# plot normed histogram
ax1.hist(times, normed=True, bins=100)
plt.xlabel(r"$T$", fontsize=17)

vline = np.linspace(0, 10.2)
o = np.ones_like(vline)*(t_max)
ax1.plot(o, vline, 'r--')

# find minimum and maximum of xticks, so we know
# where we should compute theoretical distribution
xt = plt.xticks()[0]
xmin, xmax = min(xt), max(xt)
lnspc = np.linspace(xmin, xmax, 10000)

# lets try the normal distribution first
m, s, sk = stats.lognorm.fit(times) # get mean and standard deviation
pdf_g = stats.lognorm.pdf(lnspc, m, s, sk) # now get theoretical values in our interval
ax1.plot(lnspc, pdf_g, label="Norm") # plot it


ax1.set_xlabel(r"$T$", fontsize=17)
ax1.set_xlim(0.006, .155)



ax2.set_xlabel(r"$T$", fontsize=17)
ax2.set_ylabel(r"$F$", fontsize=17)

plt.tight_layout()
plt.savefig('histT.pdf', format='pdf', dpi=300)
plt.show()
