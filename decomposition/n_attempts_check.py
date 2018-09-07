from circuit_block import Blocks
# from scipy.stats import expon
from scipy.stats import planck
from numpy.random import exponential
import matplotlib.pyplot as plt
import numpy as np

def count_plot_tpl(nums):
    """
    input: nums (list of integers)
    Returns a pair (ns, freqs) with the unique numbers in the list in
    ns, and the counts in freqs. 
    """
    ns = list(sorted(set(nums)))
    return (ns, [nums.count(_) for _ in ns])

# size: number of trials
sz = 10**3

# parameters to be passed to model
ps, pm, pg, eta, a0, a1, theta = \
(0.003, 0.003, 0.003, 0.01, 1.5, 0.0125, 0.24)

test_block = Blocks(ps, pm, pg, eta, a0, a1, theta)

# probability of trial success
p_succ = 0.1

n_attempts_lst = [test_block._success_number_of_attempts(p_succ) 
                    for _ in range(sz)]

# expon_rvs = list(map(lambda _: int(np.floor(_)),
#                     expon.rvs(scale=p_succ**-1, size=sz)))
planck_lambda = -np.log(1. - p_succ)
expon_rvs = list(planck.rvs(planck_lambda, size=sz))
# expon_rvs = list(map(lambda _: int(np.floor(_)),
                    # exponential(scale=p_succ**-1, size=sz)))

plt.plot(*count_plot_tpl(n_attempts_lst), 'k.')
plt.plot(*count_plot_tpl(expon_rvs), 'b+')
plt.show()