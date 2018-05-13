import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from pyexcel_ods import get_data
from matplotlib.ticker import FormatStrFormatter

# PQ TEST toric topology
# pq = [0.026, 0.0265, 0.027, 0.0275, 0.028, 0.0285, 0.029, 0.0295, 0.03, 0.0305,
#       0.031, 0.0315, 0.032]
# d8 = [0.0112125, 0.012475, 0.0138781, 0.0149594, 0.0162938, 0.0177094,
#       0.0192875, 0.0204875, 0.021975, 0.0239344, 0.0255, 0.0280031, 0.029025]
#
# d10 = [0.0096625, 0.0113425, 0.01227, 0.013765, 0.014925, 0.0172475, 0.018215,
#        0.020365, 0.0225, 0.0241825, 0.025835, 0.027935, 0.029765]
#
# d12 = [0.0092104, 0.0099708, 0.0118313, 0.0132479, 0.0147583, 0.0170083,
#        0.0188437, 0.0207437, 0.0226458, 0.0248146, 0.0266188, 0.02855, 0.0305646]
#
# d14 = [0.0086411, 0.0099393, 0.0113268, 0.0133018, 0.0153054, 0.0171875,
#        0.0189214, 0.0214036, 0.0233786, 0.0255839, 0.0279054, 0.0298821,
#        0.0320571]
#
# d16 = [0.0084203, 0.0097656, 0.0106703, 0.0128609, 0.0144922,
#        0.0169984, 0.0189656, 0.0220359, 0.0241891, 0.0264484, 0.0285422,
#        0.0306531, 0.0326547]


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


# Threshols over a0
# pq = [20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]
# d8 = [0.0131875, 0.014625, 0.0158125, 0.0166563, 0.01825, 0.0218125, 0.02475, 0.02575, 0.0256875, 0.0305625, 0.0349687, 0.0381563, 0.0391875, 0.0423125]
# d10= [0.011425,	0.01375, 0.014725, 0.0163, 0.0177, 0.020375, 0.024375, 0.0267, 0.02645, 0.031075, 0.035175, 0.039925, 0.040275, 0.0447]
# d12= [0.0106875, 0.0121042, 0.0140208, 0.015625, 0.0175625, 0.020375, 0.0249375, 0.0271042, 0.0266875, 0.0318333, 0.0369375, 0.0423542, 0.0430833, 0.0463958]
# d14= [0.0101607, 0.0123571, 0.0136964, 0.0152321, 0.0179286, 0.0212143, 0.0261429, 0.02775, 0.0275, 0.0333036, 0.0365536, 0.0429821, 0.0430893, 0.0468571]
# d16= [0.0102344, 0.0117031, 0.0139219, 0.0157344, 0.017875, 0.0228594, 0.0282656, 0.0294375, 0.0292969, 0.0349609, 0.0388594, 0.0438906, 0.0443906, 0.0473906]


# Threshols over eta
# pq = [0.01, 0.0095, 0.009, 0.0085, 0.008, 0.0075, 0.007, 0.0065, 0.006, 0.0055, 0.005]
# d8 = [0.0139687, 0.0155312, 0.0179062, 0.018, 0.019875, 0.0222812, 0.0246563, 0.0256875, 0.0296875, 0.0335, 0.0423125]
# d10= [0.012875, 0.012275, 0.01535, 0.016575, 0.017275, 0.019375, 0.02385, 0.023275, 0.0286, 0.032875, 0.04215]
# d12= [0.0121042, 0.0125833, 0.0152917, 0.0157708, 0.0172292, 0.0199792, 0.0237083, 0.0244167, 0.0297083, 0.0332917, 0.045875]
# d14= [0.0108214, 0.0122143, 0.015125, 0.0156786, 0.0171964, 0.0204821, 0.0243393, 0.0251429, 0.0311429, 0.03575, 0.0473214]
# d16= [0.0116406, 0.0127187, 0.0152031, 0.0164367, 0.0182187, 0.0210625, 0.025875, 0.0269219, 0.0327656, 0.0371094, 0.0457344]


data = get_data("data.ods")
data = np.array(data["Sheet3"]).transpose()
# Transform to arrays
pq = data[0]
# d6 = np.array(d6)
# d8 = data[1]
# d10 = data[1]
# d12 = data[2]
# d14 = data[3]
# d16 = data[4]
# d20 = data[4]/20

# d6 = data[1]
# d8 = data[2]
# d10 = data[1]
# d12 = data[2]
# d14 = data[3]
# d16 = data[4]

d9 = data[1]
d12 = data[2]
d15 = data[3]
d18 = data[4]

# # plt.plot(pq, 1/d14, 'o-', label=r"$d=14$")
# plt.plot(pq*100, 1/d8, 'o-', label=r"$d=8$")
# plt.plot(pq*100, 1/d12, 'o-', label=r"$d=12$")
# # plt.plot(pq, 1/d10, 'o-', label=r"$d=10$")
# plt.plot(pq*100, 1/d16, 'o-', label=r"$d=16$")





def threshold(X, A, B, C, pth, v):
    p, d = X
    pl = A + B*(p - pth)*(d**(1/v)) + C*((p - pth)*(d**(1/v)))**2
    return pl



# plt.plot(pq, d6, 'm<', label=r"$d=6$")
# plt.plot(pq, d8, 'gv', label=r"$d=8$")
# plt.plot(pq, d9, 'b*', label=r"$d=9$")
plt.plot(pq, d12, 'c>', label=r"$d=12$")
plt.plot(pq, d15, 'yo', label=r"$d=15$")
plt.plot(pq, d18, 'rs', label=r"$d=18$")
# plt.plot(pq, 1-d20, 'mo-', label=r"$d=20$")


################################################################
################################################################
##### Curve fit
o = np.ones_like(d12)
pl = np.concatenate((d12, d15, d18))
pqs = np.concatenate((pq, pq, pq))
ds = np.concatenate((o*12, o*15, o*18))

vals, pconv = curve_fit(threshold, (pqs, ds), pl, (1., 1., 1., 1., 1.))
perr = np.sqrt(np.diag(pconv))
# # plt.plot(pq, 1- threshold((pq, 6), vals[0], vals[1], vals[2], vals[3], vals[4]))
# plt.plot(pq, threshold((pq, 8), vals[0], vals[1], vals[2], vals[3], vals[4]), 'g--')
# plt.plot(pq, threshold((pq, 9), vals[0], vals[1], vals[2], vals[3], vals[4]), 'b--')
plt.plot(pq, threshold((pq, 12), vals[0], vals[1], vals[2], vals[3], vals[4]), 'c--')
plt.plot(pq, threshold((pq, 15), vals[0], vals[1], vals[2], vals[3], vals[4]), 'y--')
plt.plot(pq, threshold((pq, 18), vals[0], vals[1], vals[2], vals[3], vals[4]), 'r--')

print("VALUES")
print(vals)

print("ERROR")
print(perr)
##################################################################
##################################################################

# plt.ylim([0.4, 1])
# plt.title("EXPEDIENT")
plt.ylabel(r"$p_{logical}$", fontsize=17)
# plt.xlabel(r"Error rate", fontsize=17)
plt.xlabel(r"$p$", fontsize=17)


plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%g'))
plt.tight_layout()
# plt.xticks(pq, pq, fontsize=15)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=17)
plt.savefig('threshold_pg_hybrid.pdf', format='pdf', dpi=300)
plt.show()
