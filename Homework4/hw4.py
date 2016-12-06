from scipy.special import jv, yv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
sns.set(style="ticks", palette="muted", color_codes=True)


###############################################################################
# Problem 1 ###################################################################
###############################################################################

def Jm_(m, val):
    return 0.5 * (jv(m - 1, val) - jv(m + 1, val))

def Ym_(m, val):
    return 0.5 * (yv(m - 1, val) - yv(m + 1, val))

def func(mu, m, Ri, Ro):
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

eigenvalues.mu = eigenvalues.mu.astype(float).round(8)
eigenvalues.drop_duplicates(inplace=True)
eigenvalues.sort(['m', 'mu'], inplace=True)

eigenvalues = eigenvalues.groupby('m').head()
eigenvalues['n'] = 4 * list(np.arange(0, 5))
print(eigenvalues.pivot('n', 'm'))

###############################################################################
# Problem 2 ###################################################################
###############################################################################

def eigenfunction(m, mu, r, Ri):
    return jv(m, mu * r) - (Jm_(m, mu * Ri) / Ym_(m, mu * Ri)) * yv(m, mu * r)

f, ax = plt.subplots(nrows=5, sharex=True, sharey=True, squeeze=True)
r = np.linspace(Ri, Ro, 1000)
for n, group in eigenvalues.groupby('n'):
    for i, g in group.reset_index(drop=True).iterrows():
        ax[n].plot(r, eigenfunction(g.m, g.mu, r, Ri),
                   color=sns.color_palette()[i], alpha=0.75)
        ax[n].set_ylabel("n={}".format(n))

for ax_ in ax:
    ax_.plot([Ri, Ro], [0, 0], '--', alpha=0.25, color='k')

legend = ax[0].legend([15, 16, 17, 18], title='m',
                      loc='center left', shadow=True, bbox_to_anchor=(1, 0))
legend.get_frame().set_facecolor('#333333')

plt.xlim(Ri, Ro)
plt.ylim(-.5, .5)
plt.xlabel('Radius [in]')
plt.savefig('tex/figs/problem2.pdf')
plt.close()

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
print(res)

###############################################################################
# Problem 4 ###################################################################
###############################################################################

rho = 1.4988E-5  # slug/in^3
p = pd.DataFrame.from_csv('pressure_input.dat', sep='\t', header=None, index_col=None)
p.columns = ['r', 'Pr', 'Pi']
p['P'] = p.Pr + p.Pi * 1j

PWLs = []
m = 18
for _, (n, mu, Kz) in res.query('m == @m and n < 3')[['n', 'mu', 'Kz']].iterrows():
    mu = float(mu)
    gamma = (+ 0.5 * (Ro ** 2 - m ** 2 / mu ** 2) * eigenfunction(m, mu, Ro, Ri)**2
             - 0.5 * (Ri ** 2 - m ** 2 / mu ** 2) * eigenfunction(m, mu, Ri, Ri)**2)

    p['psi'] = eigenfunction(m, mu, p.r, Ri)
    A = (1 / gamma) * np.trapz(y=p.P * p.psi * p.r, x=p.r)
    W1 = (np.pi / (rho * c)) * gamma * A * np.conjugate(A)
    W2 = ((1 + M ** 2) * np.real(Kz / (w / c - Kz * M))
          + M * (1 + abs(Kz / (w / c - Kz * M))**2))
    Wmn = W1 * W2
    PWL = 10 * np.log10(abs(Wmn)) - 10 * np.log10(7.3756E-13)
    # print(int(n), PWL)
    PWLs.append([int(n), PWL, Wmn, Kz, mu, gamma, A])
PWLs = pd.DataFrame(PWLs)
PWLs.columns = ['n', 'PWL', 'Wmn', 'Kz', 'mu', 'gamma', 'A']
print(PWLs[['n', 'PWL']])
