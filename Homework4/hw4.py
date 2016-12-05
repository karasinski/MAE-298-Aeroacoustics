from scipy.special import jv, yv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
sns.set(style="ticks", palette="muted", color_codes=True)


# B = 18
# V = 1
# m = n * B - k * V

# M = 0.525
# RPM = 8326.3042
# c = 13503.937009  # in/s
# rho = 1.4988E-5  # slug/in^3

# dom_noise = RPM * (2 * np.pi / 60) * B

###############################################################################
# Problem 1 ###################################################################
###############################################################################

def Jm_(m, val):
    return 0.5 * (jv(m - 1, val) - jv(m + 1, val))


def Ym_(m, val):
    return 0.5 * (yv(m - 1, val) - yv(m + 1, val))


def func(mu, m, Ri, Ro):
    '''Equation taken from "Turbomachinery Noise" notes, page 10.'''
    return Jm_(m, mu * Ro) * Ym_(m, mu * Ri) - Jm_(m, mu * Ri) * Ym_(m, mu * Ro)


Ri, Ro = 3, 13  # inches

eigenvalues = []
for m in np.arange(15, 19):
    for initial_guess in np.linspace(0.1, 3, 300):
        try:
            mu = scipy.optimize.broyden1(lambda guess: func(guess, m, Ri, Ro), initial_guess)
            eigenvalues.append([m, mu])
        except Exception:
            pass
eigenvalues = pd.DataFrame(eigenvalues)
eigenvalues.columns = ['m', 'mu']

eigenvalues.mu = eigenvalues.mu.astype(float).round(4)
eigenvalues.drop_duplicates(inplace=True)
eigenvalues.sort(['m', 'mu'], inplace=True)

eigenvalues = eigenvalues.groupby('m').head()
eigenvalues['n'] = 4 * list(np.arange(0, 5))
print(eigenvalues.pivot('n', 'm'))

###############################################################################
# Problem 2 ###################################################################
###############################################################################


def eigenfunction(m, mu, r, Ri):
    '''Equation taken from "Turbomachinery Noise" notes, page 10.'''
    return jv(m, mu * r) - (Jm_(m, mu * Ri) / Ym_(m, mu * Ri)) * yv(m, mu * r)

f, ax = plt.subplots(nrows=5, sharex=True, sharey=True, squeeze=True)
r = np.linspace(Ri, Ro, 1000)
for n, group in eigenvalues.groupby('n'):
    for _, g in group.iterrows():
        ax[n].plot(r, eigenfunction(g.m, g.mu, r, Ri), color=sns.color_palette()[n])
        # ax[n].legend([n], loc='upper left')
        ax[n].set_ylabel("n={}".format(n))
        ax[n].plot([Ri, Ro], [0, 0], '--', alpha=0.5, color='k')

plt.xlim(Ri, Ro)
plt.ylim(-.5, .5)
plt.xlabel('r')
# plt.show()
plt.savefig('tex/figs/problem2.pdf')

###############################################################################
# Problem 3 ###################################################################
###############################################################################

M = 0.525
RPM = 8326.3042  # rotations per minute
B = 18
w = RPM * (2 * np.pi / 60) * B  # angular velocity
c = 13503.937009  # in/s
K = w / c

res = []
for i, (m, mu, n) in eigenvalues.iterrows():
    Kz = K * (- M + np.emath.sqrt(1 - (1 - M**2) * (mu / K)**2)) / (1 - M**2)
    Kzr = float(np.real(Kz))
    Kzi = float(np.imag(Kz))
    res.append([m, n, mu, Kz, Kzr, Kzi])
res = pd.DataFrame(res, columns=['m', 'n', 'mu', 'Kz', 'Kzr', 'Kzi'])

###############################################################################
# Problem 4 ###################################################################
###############################################################################

p = pd.DataFrame.from_csv('pressure_input.dat', sep='\t', header=None, index_col=None)


