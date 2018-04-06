import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

def lambda_env(t, a0, a1):
    a = (a0 + a1)*t
    lamb = (1 + np.exp(-a * t))/2.
    return 1 - lamb

# Data
T4epl = [0.1233188]
T4single = [0.1347472]
T4double= [0.2303306]

T3epl = [0.0570291]
T3single = [0.0831743]
T3double= [0.1318343]

T2epl = [0.0398776]
T2single = [0.0629081]
T2double= [0.1057237]

F4epl = [0.979317452327, 0.00195496048824]
F4single = [0.966544977088, 0.00210544022824]
F4double = [0.974064936268, 0.00702994574387]

F3epl = [0.985404250679, 0.000731628944008]
F3single = [0.976730448884, 0.000936538399082]
F3double = [0.984126242197, 0.000878146598433]

F2epl = [0.991980711676, 4.49439761249e-06]
F2single = [0.996317975745, 1.31525620288e-05]
F2double = [0.996929904009, 6.90980714079e-05]



plt.errorbar(T4epl, F4epl[0], F4epl[1], fmt='bs', label="EPL", markersize=8)
plt.errorbar(T3epl, F3epl[0], F3epl[1], fmt='b^', label="EPL", markersize=8)
plt.errorbar(T2epl, F2epl[0], F2epl[1], fmt='bx', label="EPL", markersize=8)
plt.errorbar(T4single, F4single[0], F4single[1], fmt='rs', label="Purification1", markersize=8)
plt.errorbar(T3single, F3single[0], F3single[1], fmt='r^', label="Purification1", markersize=8)
plt.errorbar(T2single, F2single[0], F2single[1], fmt='rx', label="Purification1", markersize=8)
plt.errorbar(T4double, F4double[0], F4double[1], fmt='gs', label="Purification2", markersize=8)
plt.errorbar(T3double, F3double[0], F3double[1], fmt='g^', label="Purification2", markersize=8)
plt.errorbar(T2double, F2double[0], F2double[1], fmt='gx', label="Purification2", markersize=8)



# Horizontals and verticals
ver = np.linspace(0.5, 1)
hor = np.linspace(0, 2)
ones = np.ones(50)
# plt.plot(0.031* ones, ver, 'y--')
# plt.plot(0.8* ones, ver, 'r--')
# plt.plot(hor, 0.95* ones, 'g--')



#
# plt.ylim([0.7, 1])
# plt.xlim([0., 1.4])
plt.ylabel(r"Fidelity", fontsize=17)
plt.xlabel(r"Time (seg)", fontsize=17)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

blue_patch = mpatches.Patch(color='blue', label='EPL')
red_patch = mpatches.Patch(color='red', label='Purification 1')
green_patch = mpatches.Patch(color='green', label='Purification 2')

plt.legend(fontsize=14, handles=[blue_patch, red_patch, green_patch])
plt.show()
