from scipy import stats
import numpy as np
import matplotlib.pylab as plt

ignore = [0.0, 2.5, 5, 7.5, 10]
thresholds = [0.00693, 0.00688, 0.00685, 0.0067, 0.00658]
errs = [0.00005, 0.00004, 0.00005, 0.00004, 0.00005]

plt.errorbar(ignore, thresholds, yerr=errs, fmt='ro')
plt.ylabel(r"threshold", fontsize=17)
# plt.xlabel(r"Error rate", fontsize=17)
plt.xlabel(r"% clasical erasure", fontsize=17)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
# plt.legend(fontsize=17)
plt.tight_layout()
plt.savefig('erasure_thresholds.pdf', format='pdf', dpi=300)

plt.show()
