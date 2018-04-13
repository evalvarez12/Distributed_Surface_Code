from scipy import stats
import numpy as np
import matplotlib.pylab as plt

# create some normal random noisy data
# ser = 50*np.random.rand() * np.random.normal(10, 10, 100) + 20
ser = np.load("../decomposition/TIME_EPL_2000it.npy")

# plot normed histogram
plt.hist(ser, normed=True, bins=100)

# find minimum and maximum of xticks, so we know
# where we should compute theoretical distribution
xt = plt.xticks()[0]
xmin, xmax = min(xt), max(xt)
lnspc = np.linspace(xmin, xmax, 10000)

# lets try the normal distribution first
m, s, sk = stats.lognorm.fit(ser) # get mean and standard deviation
pdf_g = stats.lognorm.pdf(lnspc, m, s, sk) # now get theoretical values in our interval
plt.plot(lnspc, pdf_g, label="Norm") # plot it

# exactly same as above
# al, bl = stats.levy.fit(ser)
# pdf_levy = stats.levy.pdf(lnspc, al, bl)
# plt.plot(lnspc, pdf_levy, label="Levy")

# Exponential distribution
ae,be = stats.expon.fit(ser)
pdf_expon = stats.expon.pdf(lnspc, ae, be)
plt.plot(lnspc, pdf_expon, label="Exponential")

# guess what :)
# ab,bb,cb,db = stats.beta.fit(ser)
# pdf_beta = stats.beta.pdf(lnspc, ab, bb,cb, db)
# plt.plot(lnspc, pdf_beta, label="Beta")

plt.legend()

plt.plot(lnspc, pdf_expon > 1, 'k,')
tmax = np.argmax(pdf_expon[1000:] < 1)
print("TMAX: ", tmax)
# plt.figure()
# plt.plot(np.sort(ser), 'b.')

plt.show()
