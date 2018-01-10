import matplotlib.pyplot as plt
import numpy as np

# Over a0
a0 = [2, 4, 6, 8]
GHZ_8 = [.0607, .3914, .4369, .7445]
GHZ_12 = [.0186, .1997, .3227, .7609]
GHZ_16 = [.0071, .1976, .3172, .7904]


a03 = [2, 4, 6, 8, 12]
GHZ3_8 = [.3745, .7415, .8347, .91, .9143]
GHZ3_12 = [.231, .7215, .7444, .8817, .9374]
GHZ3_16 = [.1945, .7569, .8136, .9104, .9206]


# Over eta
# eta = [1/30, 1/40, 1/50, 1/60, 1/70, 1/80]
# GHZ_8 = [0.0144, 0.0449, 0.8002, 0.934, 0.9362, 0.9372]
# GHZ_12 = [0.0014, 0.0093, 0.8289, 0.9401, 0.9425, 0.9493]
# GHZ_16 = [0.0009, 0.003, 0.8117, 0.944, 0.949, 0.953]


plt.plot(a03, GHZ3_8, 'o-', label=r"$d=8$")
plt.plot(a03, GHZ3_12, 'o-', label=r"$d=12$")
plt.plot(a03, GHZ3_16, 'o-', label=r"$d=16$")


# plt.plot(eta, GHZ_8, 'o-', label=r"$d=8$")
# plt.plot(eta, GHZ_12, 'o-', label=r"$d=12$")
# plt.plot(eta, GHZ_16, 'o-', label=r"$d=16$")


# plt.ylim([0.4, 1])
plt.ylabel(r"Fail rate", fontsize=13)
plt.xlabel(r"$a0$", fontsize=13)
plt.legend(fontsize=13)
plt.show()
