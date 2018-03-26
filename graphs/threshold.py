import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# PQ TEST toric topology
pq = [0.026, 0.0265, 0.027, 0.0275, 0.028, 0.0285, 0.029, 0.0295, 0.03, 0.0305,
      0.031, 0.0315, 0.032]
d8 = [0.0112125, 0.012475, 0.0138781, 0.0149594, 0.0162938, 0.0177094,
      0.0192875, 0.0204875, 0.021975, 0.0239344, 0.0255, 0.0280031, 0.029025]

d10 = [0.0096625, 0.0113425, 0.01227, 0.013765, 0.014925, 0.0172475, 0.018215,
       0.020365, 0.0225, 0.0241825, 0.025835, 0.027935, 0.029765]

d12 = [0.0092104, 0.0099708, 0.0118313, 0.0132479, 0.0147583, 0.0170083,
       0.0188437, 0.0207437, 0.0226458, 0.0248146, 0.0266188, 0.02855, 0.0305646]

d14 = [0.0086411, 0.0099393, 0.0113268, 0.0133018, 0.0153054, 0.0171875,
       0.0189214, 0.0214036, 0.0233786, 0.0255839, 0.0279054, 0.0298821,
       0.0320571]

d16 = [0.0084203, 0.0097656, 0.0106703, 0.0128609, 0.0144922,
       0.0169984, 0.0189656, 0.0220359, 0.0241891, 0.0264484, 0.0285422,
       0.0306531, 0.0326547]


# PQ TEST planar topology
# pq = [0.026, 0.0265, 0.027, 0.0275, 0.028, 0.0285, 0.029, 0.0295, 0.03, 0.0305,
#       0.031, 0.0315, 0.032]
#
# d8 = [0.0092458, 0.0100453, 0.0106766, 0.0116203, 0.0126844, 0.0137844,
#       0.0141406, 0.0153156, 0.0164438, 0.0170594, 0.0183859, 0.0195109,
#       0.0202813]
#
# d10 = [0.00814, 0.008925, 0.0098638, 0.0106175, 0.0117562, 0.0125588,
#        0.0136525, 0.0148325, 0.015835, 0.016705, 0.0178, 0.0189713, 0.0205463]
#
# d12 = [0.0073944, 0.0083042, 0.0091708, 0.0103167, 0.0113833, 0.0121271, 0.0131896,
#        0.0145542, 0.0154104, 0.0164187, 0.018025, 0.0192167, 0.0204208]
#
# d14 = [0.0068333, 0.0076982, 0.0087714, 0.0098929, 0.0110732, 0.0121679,
#        0.0134946, 0.014675, 0.0160054, 0.016625, 0.0184357, 0.0196179, 0.0207125]
#
# d16 = [0.0066875, 0.0076734, 0.0087562, 0.0097859, 0.0108766, 0.0121094,
#        0.0134984, 0.0149484, 0.0160281, 0.0171891, 0.0187953, 0.0199812, 0.0210406]


# Transform to arrays
pq = np.array(pq)
# d6 = np.array(d6)
d8 = np.array(d8)
d10 = np.array(d10)
d12 = np.array(d12)
d14 = np.array(d14)
d16 = np.array(d16)


# # plt.plot(pq, 1/d14, 'o-', label=r"$d=14$")
# plt.plot(pq*100, 1/d8, 'o-', label=r"$d=8$")
# plt.plot(pq*100, 1/d12, 'o-', label=r"$d=12$")
# # plt.plot(pq, 1/d10, 'o-', label=r"$d=10$")
# plt.plot(pq*100, 1/d16, 'o-', label=r"$d=16$")





def threshold(X, A, B, C, pth, v):
    p, d = X
    pl = A + B*(p - pth)*(d**(1/v)) + C*((p - pth)*(d**(1/v)))**2
    return pl

o = np.ones_like(d8)
pl = np.concatenate((d8, d10,  d12, d14, d16))
pqs = np.concatenate((pq, pq, pq, pq, pq))
ds = np.concatenate((o*8, o*10, o*12, o*14, o*16))

vals, pconv = curve_fit(threshold, (pqs, ds), pl, (0.5, 0.5, 0.5, 0.1, .2))
perr = np.sqrt(np.diag(pconv))
# plt.plot(pq, 1-d6, 'ro', label=r"$d=6$")
plt.plot(pq, 1-d8, 'gv', label=r"$d=8$")
plt.plot(pq, 1-d10, 'b*', label=r"$d=10$")
plt.plot(pq, 1-d12, 'c>', label=r"$d=12$")
plt.plot(pq, 1-d14, 'yo', label=r"$d=14$")
plt.plot(pq, 1-d16, 'rs', label=r"$d=16$")

# plt.plot(pq, 1- threshold((pq, 6), vals[0], vals[1], vals[2], vals[3], vals[4]))
plt.plot(pq, 1- threshold((pq, 8), vals[0], vals[1], vals[2], vals[3], vals[4]), 'g--')
plt.plot(pq, 1- threshold((pq, 10), vals[0], vals[1], vals[2], vals[3], vals[4]), 'b--')
plt.plot(pq, 1- threshold((pq, 12), vals[0], vals[1], vals[2], vals[3], vals[4]), 'c--')
plt.plot(pq, 1-threshold((pq, 14), vals[0], vals[1], vals[2], vals[3], vals[4]), 'y--')
plt.plot(pq, 1-threshold((pq, 16), vals[0], vals[1], vals[2], vals[3], vals[4]), 'r--')




# plt.ylim([0.4, 1])
# plt.title("EXPEDIENT")
plt.ylabel(r"Success rate", fontsize=17)
plt.xlabel(r"Error rate", fontsize=17)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=17)
plt.tight_layout()
plt.savefig('sc.pdf', format='pdf', dpi=300)
plt.show()
