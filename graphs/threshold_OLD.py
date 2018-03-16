import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Over a0
# a0 = [2, 4, 6, 8]
# GHZ_8 = [.0607, .3914, .4369, .7445]
# GHZ_12 = [.0186, .1997, .3227, .7609]
# GHZ_16 = [.0071, .1976, .3172, .7904]
#
#
# a03 = [2, 4, 6, 8, 12]
# GHZ3_8 = [.3745, .7415, .8347, .91, .9143]
# GHZ3_12 = [.231, .7215, .7444, .8817, .9374]
# GHZ3_16 = [.1945, .7569, .8136, .9104, .9206]


# Over eta
# eta = [1/30, 1/40, 1/50, 1/60, 1/70, 1/80]
# GHZ_8 = [0.0144, 0.0449, 0.8002, 0.934, 0.9362, 0.9372]
# GHZ_12 = [0.0014, 0.0093, 0.8289, 0.9401, 0.9425, 0.9493]
# GHZ_16 = [0.0009, 0.003, 0.8117, 0.944, 0.949, 0.953]


# PQ test cycles=50
# pq = [0.0060, 0.0080, 0.0100, 0.0120, 0.0140, 0.0160, 0.0180, 0.0200, 0.0220, 0.0240, 0.0260, 0.0280, 0.0300, 0.0320, 0.0340]
# d4 = [0.0479, 0.0827, 0.1380, 0.2012, 0.2735, 0.3440, 0.4209, 0.4877, 0.5516, 0.6133, 0.6542, 0.6831, 0.7151, 0.7169, 0.7339]
# d8 = [0.0006, 0.0019, 0.0033, 0.0099, 0.0275, 0.0535, 0.1130, 0.1919, 0.2902, 0.4146, 0.5455, 0.6370, 0.6969, 0.7189, 0.7392]
# d12 =[0.0000, 0.0001, 0.0001, 0.0011, 0.0041, 0.0148, 0.0403, 0.0866, 0.1951, 0.3364, 0.5012, 0.6349, 0.7074, 0.7361, 0.7521]
# d16 =[0.0000, 0.0000, 0.0000, 0.0004, 0.0008, 0.0049, 0.0167, 0.0513, 0.1463, 0.3073, 0.5068, 0.6589, 0.7261, 0.7427, 0.7476]

# With cycles=distance
# d4 = [0.0062, 0.0087, 0.0143, 0.0213, 0.0281, 0.0393, 0.0527, 0.0668, 0.0810, 0.1067, 0.1265, 0.1473, 0.1722, 0.1936, 0.2183]
# d8 = [0.0001, 0.0004, 0.0019, 0.0040, 0.0079, 0.0132, 0.0270, 0.0467, 0.0745, 0.1043, 0.1487, 0.2052, 0.2644, 0.3373, 0.3930]
# d12= [0.0000, 0.0001, 0.0004, 0.0011, 0.0027, 0.0068, 0.0170, 0.0335, 0.0755, 0.1336, 0.2086, 0.2916, 0.3980, 0.5007, 0.5706]
# d16= [0.0000, 0.0000, 0.0001, 0.0001, 0.0009, 0.0031, 0.0095, 0.0301, 0.0685, 0.1462, 0.2700, 0.3878, 0.5233, 0.6310, 0.6825]


# d4 = np.array(d4)
# d8 = np.array(d8)
# d12 = np.array(d12)
# d16 = np.array(d16)
#
# plt.plot(pq, 1-d4, 'o-', label=r"$d=4$")
# plt.plot(pq, 1-d8, 'o-', label=r"$d=8$")
# plt.plot(pq, 1-d12, 'o-', label=r"$d=12$")
# plt.plot(pq, 1-d16, 'o-', label=r"$d=16$")


# cycles=2distance
# d4 = [0.0078, 0.0152, 0.0240, 0.0404, 0.0554, 0.0742, 0.0927, 0.1320, 0.1611, 0.1879, 0.2183, 0.2623, 0.2971, 0.3371, 0.3698]
# d8 = [0.0000, 0.0005, 0.0033, 0.0052, 0.0111, 0.0198, 0.0429, 0.0798, 0.1245, 0.1864, 0.2585, 0.3472, 0.4274, 0.4986, 0.5801]
# d12= [0.0000, 0.0000, 0.0000, 0.0004, 0.0026, 0.0101, 0.0245, 0.0518, 0.1105, 0.1995, 0.3313, 0.4606, 0.5856, 0.6654, 0.7151]
# d16= [0.0001, 0.0000, 0.0001, 0.0003, 0.0014, 0.0024, 0.0115, 0.0386, 0.1097, 0.2165, 0.4007, 0.5613, 0.6832, 0.7327, 0.7467]

# cycles=3distance
# d4 = [0.0126, 0.0247, 0.0387, 0.0558, 0.0759, 0.1067, 0.1368, 0.1795, 0.2090, 0.2605, 0.3042, 0.3500, 0.3791, 0.4332, 0.4809]
# d8 = [0.0001, 0.0004, 0.0024, 0.0053, 0.0123, 0.0313, 0.0619, 0.1021, 0.1708, 0.2503, 0.3429, 0.4416, 0.5356, 0.6090, 0.6716]
# d12= [0.0000, 0.0000, 0.0002, 0.0010, 0.0043, 0.0097, 0.0287, 0.0728, 0.1538, 0.2717, 0.4189, 0.5622, 0.6632, 0.7217, 0.7492]
# d16= [0.0000, 0.0000, 0.0000, 0.0000, 0.0014, 0.0041, 0.0155, 0.0464, 0.1366, 0.2967, 0.4857, 0.6524, 0.7243, 0.7529, 0.7456]

###################3 Using the cheat decoder, c=d
# pq = [0.0060, 0.0080, 0.0100, 0.0120, 0.0140, 0.0160, 0.0180, 0.0200, 0.0220, 0.0240, 0.0260, 0.0280, 0.0300, 0.0320, 0.0340]
# d4 = [0.0038, 0.0073, 0.0130, 0.0148, 0.0242, 0.0283, 0.0412, 0.0519, 0.0714, 0.0832, 0.1012, 0.1151, 0.1360, 0.1600, 0.1785]
# d8 = [0.0000, 0.0001, 0.0003, 0.0008, 0.0018, 0.0057, 0.0133, 0.0214, 0.0361, 0.0598, 0.0903, 0.1337, 0.1775, 0.2261, 0.2876]
# d12= [0.0000, 0.0000, 0.0000, 0.0003, 0.0007, 0.0016, 0.0044, 0.0126, 0.0262, 0.0572, 0.1107, 0.1765, 0.2641, 0.3674, 0.4745]
# d16= [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0005, 0.0013, 0.0084, 0.0238, 0.0589, 0.1381, 0.2498, 0.3843, 0.5233, 0.6177]
#
# Cheat decoder c=100
# d4 = [0.0890, 0.1594, 0.2392, 0.3359, 0.4301, 0.5176, 0.5943, 0.6570, 0.6954, 0.7172, 0.7348, 0.7458, 0.7511, 0.7564, 0.7511]
# d8 = [0.0004, 0.0010, 0.0039, 0.0155, 0.0424, 0.0943, 0.1785, 0.3018, 0.4455, 0.5802, 0.6709, 0.7307, 0.7431, 0.7374, 0.7533]
# d12= [0.0000, 0.0001, 0.0002, 0.0015, 0.0034, 0.0148, 0.0490, 0.1318, 0.2991, 0.4861, 0.6413, 0.7249, 0.7543, 0.7495, 0.7467]



# P test Q=0
# pq = [0.095, 0.096, 0.097, 0.098, 0.100, 0.101, 0.102, 0.103, 0.104, 0.105]
# d4 =[0.2474, 0.2644, 0.2641, 0.2699, 0.2892, 0.2786, 0.2824, 0.2786, 0.3005, 0.2982]
# d8 =[0.2221, 0.2295, 0.2331, 0.2402, 0.2664, 0.2568, 0.2814, 0.2812, 0.3008, 0.2970]
# d12=[0.2028, 0.2094, 0.2242, 0.2277, 0.2403, 0.2596, 0.2709, 0.2754, 0.2964, 0.3041]
# d16=[0.1903, 0.1952, 0.2086, 0.2204, 0.2446, 0.2553, 0.2812, 0.2867, 0.2896, 0.3008]

# plt.plot(a03, GHZ3_8, 'o-', label=r"$d=8$")
# plt.plot(a03, GHZ3_12, 'o-', label=r"$d=12$")
# plt.plot(a03, GHZ3_16, 'o-', label=r"$d=16$")
#

# plt.plot(eta, GHZ_8, 'o-', label=r"$d=8$")
# plt.plot(eta, GHZ_12, 'o-', label=r"$d=12$")
# plt.plot(eta, GHZ_16, 'o-', label=r"$d=16$")


# Recalculation of the PQ test ------------------------------------------------------------------------------
#### Equal X Y Z error p/3
#### cycles = d
pq = [0.0220, 0.0240, 0.0260, 0.0280, 0.0300, 0.0320, 0.0340]
# d8 = [0.0067, 0.0130, 0.0230, 0.0327, 0.0520, 0.0770, 0.0970]
# d12= [0.0018, 0.0028, 0.0087, 0.0212, 0.0335, 0.0627, 0.0985]
# d16= [0.0005, 0.0013, 0.0042, 0.0093, 0.0239, 0.0520, 0.1062]
# #
# # #### cycles = 3d
# d8 = [0.0323, 0.0541, 0.0918, 0.1326, 0.1946, 0.2634, 0.3446]
# d12= [0.0067, 0.0174, 0.0353, 0.0752, 0.1362, 0.2366, 0.3664]
# d16= [0.0017, 0.0061, 0.0142, 0.0451, 0.1037, 0.2208, 0.3981]
# # #

# # ### Only Z err and q stab lie in stars
# # #### cycles = d
d8 = [0.0361, 0.0620, 0.0869, 0.1315, 0.1798, 0.2300, 0.2899]
d12= [0.0291, 0.0591, 0.1117, 0.1777, 0.2718, 0.3722, 0.4676]
d16= [0.0206, 0.0604, 0.1320, 0.2499, 0.3927, 0.5176, 0.6273]

#### GOOD cycles = 16
# d8 = [0.0918, 0.1361, 0.2048, 0.2855, 0.3773, 0.4559, 0.5421]
# d12= [0.0438, 0.0881, 0.1639, 0.2541, 0.3709, 0.4834, 0.5845]
# d16= [0.0215, 0.0582, 0.1326, 0.2521, 0.3909, 0.5263, 0.6187]

## GOOD one cycles = d already divided
# pq=[0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033]
# d8= [0.01136, 0.01354, 0.01627, 0.019, 0.02303, 0.02551, 0.02909, 0.03301]
# d10= [0.01015, 0.01207, 0.01499, 0.0185, 0.02216, 0.02574, 0.02975, 0.03366]
# d12= [0.00894, 0.01163, 0.01546, 0.0187, 0.02298, 0.02725, 0.03089, 0.03509]
# d14= [0.00875, 0.01178, 0.01508, 0.0191, 0.02391, 0.02781, 0.03203, 0.03622]





# Circuit error rate p test --------------------------------------------
#### cycles = d
# pq = [0.0040, 0.0045, 0.0050, 0.0055, 0.0060, 0.0065, 0.0070, 0.0075, 0.0080, 0.0085, 0.0090]

# WRONG IMPLEMENTATION
# d8= [0.00045, 0.00030, 0.00105, 0.00175, 0.00195, 0.00420, 0.00790, 0.00720, 0.01290, 0.01965, 0.01585]
# d12=[0.00000, 0.00005, 0.00015, 0.00045, 0.00075, 0.00035, 0.00150, 0.00160, 0.00835, 0.00350, 0.01025]
# d16=[0.00000, 0.00000, 0.00005, 0.00000, 0.00005, 0.00035, 0.00120, 0.00080, 0.00380, 0.00315, 0.00790]

# STILL WRONG
# d8 = [0.00280, 0.00715, 0.00975, 0.01420, 0.02070, 0.02575, 0.03380, 0.05875, 0.07865, 0.08710, 0.11765]
# d12= [0.00045, 0.00130, 0.00260, 0.00335, 0.00790, 0.01670, 0.02275, 0.04380, 0.05690, 0.08790, 0.12240]
# d16= [0.00005, 0.00025, 0.00055, 0.00180, 0.00420, 0.00840, 0.01770, 0.02985, 0.05855, 0.08950, 0.19810]

# STILL WRONG
# d8 = [0.00216667, 0.00296667, 0.00456667, 0.00671667, 0.01050000, 0.01550000, 0.02168333, 0.02840000, 0.03746667, 0.04950000, 0.06291667]
# d12= [0.00010000, 0.00041667, 0.00080000, 0.00141667, 0.00331667, 0.00550000, 0.00961667, 0.01551667, 0.02301667, 0.03625000, 0.05433333]
# d16= [0.00003333, 0.00010000, 0.00020000, 0.00043333, 0.00078333, 0.00216667, 0.00465000, 0.00861667, 0.01568333, 0.02828333, 0.05267773]

# cycles = d WRONG!!!!!
# pq = [0.0040, 0.0045, 0.0050, 0.0055, 0.0060, 0.0065, 0.0070, 0.0075, 0.0080, 0.0085, 0.0090, 0.0095, 0.0100, 0.0105, 0.0110]
# d8 = [0.00100, 0.00250, 0.00420, 0.00705, 0.01085, 0.01630, 0.02080, 0.03025, 0.04260, 0.05880, 0.07310, 0.09485, 0.12130, 0.14930, 0.18585]
# d12= [0.00000, 0.00030, 0.00065, 0.00170, 0.00265, 0.00550, 0.00975, 0.01710, 0.02435, 0.04440, 0.06490, 0.10475, 0.14100, 0.18810, 0.24800]
# d16= [0.00000, 0.00000, 0.00010, 0.00025, 0.00065, 0.00215, 0.00485, 0.00865, 0.01945, 0.03895, 0.06370, 0.11195, 0.16875, 0.23805, 0.33380]



# NN monolithic cycles = 100
# pq = [0.006, 0.0065, 0.007, 0.0075, 0.008, 0.0085, 0.009, 0.0095, 0.01, 0.0105,
# 0.011, 0.0115, 0.012, 0.0125, 0.013, 0.0135, 0.014]
#
# d8 = [0.00689, 0.00773, .00831, 0.00874, .00901, 0.00916, 0.00926, 0.00934,
#       0.00932, 0.00936, .00936, 0.00938, .00937, 0.00939, .00938, 0.00938,
#       0.00937]
# d10= [0.00527, 0.00650, 0.00759, 0.00833, 0.00881, 0.00913, 0.00926, 0.00933, 0.00935,
#       0.00937, 0.00937,
#       0.00937, 0.00939,
#       0.00936, 0.00936,
#       0.00938, 0.00936]
#
# d12= [0.004, 0.00553, 0.00696, 0.00801, 0.00872, 0.00911, 0.00927, 0.00934,
#       0.00933, 0.00938, 0.00938, 0.00937, 0.00937, 0.00939, 0.0094, 0.00937,
#       0.00939]
#
# d14= [0.00292, 0.00461, 0.00638, 0.00774, 0.00865, 0.0091, 0.0093, 0.00935,
#       0.00939, 0.00938, 0.00936, 0.00938, 0.00936, 0.00939, 0.00938, 0.00938,
#       0.00938]


## NN nat comm comparation
## EXPEDIENT
# pq = [0.00600, 0.00625, 0.00650]
# d8 = [0.00839, 0.00874, 0.00898]
# d12= [0.00801, 0.00864, 0.00903]
# d16= [0.00793, 0.00868, 0.00909]

# # STRINGENT
# pq = [0.00775, 0.00800, 0.00825]
# d8 = [0.00884, 0.00905, 0.00918]
# d12= [0.00882, 0.00906, 0.00921]
# d16= [0.00876, 0.00909, 0.00923]




#### Understatiding WTF
# #### d = 8 cycles = [d, 2d, 3d, 4d, 5d] resuls non divided by cycles
# cycles = [8, 16, 24, 32, 40, 48]
# ### Below threshold pq = 0.015
# b = [0.005, 0.007, 0.01375, 0.02, 0.028, 0.0265]
# ### Close to threshold pq = 0.027
# c = [0.109, 0.25050, 0.35925, 0.42375, 0.5115, 0.56425]
# ### Beyond threshold pq = 0.04
# a = [0.478, 0.679, 0.72925, 0.7295, 0.755, 0.7445]
#
# ### Beyond threshold pq = 0.10
# a = [0.7585, 0.746]
#
#
#





# fact=.61
pq = np.array(pq)
# d4 = np.array(d4)
d8 = np.array(d8)*10/8.
# d10 = np.array(d10)
d12 = np.array(d12)*10/12.
# d14 = np.array(d14)
d16 = np.array(d16)*10/16.
# d8 = np.array(d82)/16
# d12 = np.array(d122)/16
# d16 = np.array(d162)/16


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
pl = np.concatenate((d8, d12, d16))
pqs = np.concatenate((pq, pq, pq))
ds = np.concatenate((o*8, o*12, o*16))

vals, errors = curve_fit(threshold, (pqs, ds), pl, (0.5, 0.5, 0.5, 0.1, .2))

plt.plot(pq, 1-d8, 'o', label=r"$d=8$")
plt.plot(pq, 1-d12, 'o', label=r"$d=12$")
# plt.plot(pq, 1-d14, 'o-', label=r"$d=14$")
plt.plot(pq, 1-d16, 'o', label=r"$d=16$")

plt.plot(pq, 1- threshold((pq, 8), vals[0], vals[1], vals[2], vals[3], vals[4]))
plt.plot(pq, 1- threshold((pq, 12), vals[0], vals[1], vals[2], vals[3], vals[4]))
plt.plot(pq, 1-threshold((pq, 16), vals[0], vals[1], vals[2], vals[3], vals[4]))




# plt.ylim([0.4, 1])
# plt.title("EXPEDIENT")
plt.ylabel(r"Average time until failure", fontsize=13)
plt.xlabel(r"Error rate on all operations", fontsize=13)
plt.legend(fontsize=13)
plt.show()




# ('p=', 0.006, ' : ', 0.00839)
# ('p=', 0.00625, ' : ', 0.00874)
# ('p=', 0.0065, ' : ', 0.00898)
# ('p=', 0.006, ' : ', 0.00801)
# ('p=', 0.00625, ' : ', 0.00864)
# ('p=', 0.0065, ' : ', 0.00903)
# ('p=', 0.006, ' : ', 0.00793)
# ('p=', 0.00625, ' : ', 0.00868)
# ('p=', 0.0065, ' : ', 0.00909)
# ('p=', 0.00775, ' : ', 0.00884)
# ('p=', 0.008, ' : ', 0.00905)
# ('p=', 0.00825, ' : ', 0.00916)
# ('p=', 0.00775, ' : ', 0.00876)
# ('p=', 0.008, ' : ', 0.00906)
# ('p=', 0.00825, ' : ', 0.00921)
# ('p=', 0.00775, ' : ', 0.00882)
# ('p=', 0.008, ' : ', 0.00911)
# ('p=', 0.00825, ' : ', 0.00927)