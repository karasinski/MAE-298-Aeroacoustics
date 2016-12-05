import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import jv, yv
from scipy.optimize import fsolve
import scipy


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


def eigenfunction(mu, m, Ri, Ro):
    '''Equation taken from "Turbomachinery Noise" notes, page 10.'''
    return Jm_(m, mu * Ro) * Ym_(m, mu * Ri) - Jm_(m, mu * Ri) * Ym_(m, mu * Ro)


Ri, Ro = 3, 13  # inches

eigenvalues = []
for m in np.arange(15, 19):
    for initial_guess in np.linspace(0.1, 3, 300):
        try:
            mu = scipy.optimize.broyden1(lambda guess: eigenfunction(guess, m, Ri, Ro), initial_guess)
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

# ###############################################################################
# # Problem 2 ############################################################
# ###############################################################################


# # def eigenfunction(m, mu, r, Ri):
# #     '''Equation taken from "Turbomachinery Noise" notes, page 10.'''
# #     return jv(m, mu * r) - (Jm_(m, mu * Ri) / Ym_(m, mu * Ri)) * yv(m, mu * r)

# r = np.linspace(Ri, Ro, 1000)
# for m, mu in eigenvalues:
#     plt.plot(r, eigenfunction(r, m, Ri, Ro))
#     # plt.plot(r, eigenfunction(m, mu, r, Ri))
#     plt.ylim(-1, 1)
# plt.legend(eigenvalues[:, 0].astype(int))
# plt.xlim(Ri, Ro)
# plt.show()


# ###############################################################################
# # Problem 3 ############################################################
# ###############################################################################

# M = 0.525
# RPM = 8326.3042  # rotations per minute
# w = RPM * (2 * np.pi / 60)  # angular velocity
# c = 13503.937009  # in/s
# K = w / c

# res = []
# for [m, mu] in eigenvalues:
#     Kz = K * (- M + np.emath.sqrt(1 - (1 - M**2) * (mu / K)**2)) / (1 - M**2)
#     res.append([m, mu, Kz])
