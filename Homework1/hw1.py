import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal
import soundfile as sf

###############################################################################
# Problem 1
###############################################################################

# Load the audio file
audio, fs = sf.read('Boom_F1B2_6.wav')

df = pd.DataFrame(audio)
df.columns = ['Voltage']
df['Time'] = np.array(df.index) / fs

# Convert to pascals
df['Pascals'] = df.Voltage * -116.

# Find the peak pressure
maxPeak = df.ix[df.Pascals.idxmax()]
minPeak = df.ix[df.Pascals.idxmin()]


def plot1a():
    f, ax = plt.subplots()
    df.plot(x='Time', y='Pascals', legend=False, ax=ax)

    plt.xlim(xmin=0.2)
    plt.ylim(-100, 100)
    plt.ylabel('Pressure (Pascals)')

    arrowprops = dict(arrowstyle="-",
                      shrinkA=5, shrinkB=5,
                      connectionstyle="arc3")

    ax.annotate('Peak Pressure, {:2.2f} Pa'.format(maxPeak.Pascals),
                xy=(maxPeak.Time, maxPeak.Pascals),
                xytext=(25, 25), textcoords='offset points',
                arrowprops=arrowprops,
                horizontalalignment='left', verticalalignment='bottom')

    ax.annotate('',
                xy=(minPeak.Time, minPeak.Pascals),
                xytext=(maxPeak.Time, minPeak.Pascals),
                arrowprops=arrowprops)

    plt.text((maxPeak.Time + minPeak.Time) / 2, 1.1 * minPeak.Pascals,
             '{:2.2f} sec'.format(minPeak.Time - maxPeak.Time),
             horizontalalignment='center')

    plt.vlines(x=maxPeak.Time, ymax=0, ymin=minPeak.Pascals, color='k', linestyles='--')

    plt.tight_layout()
    plt.savefig('tex/figs/Pascals_vs_Time.pdf')
    plt.show()
#plot1a()

################################################################################
# Calculate the single sided power spectral density function
f, Gxx = scipy.signal.periodogram(df.Pascals, fs=fs, return_onesided=True)


def plot1b():
    plt.plot(f, Gxx, marker='o', color='k')

    plt.xlim(0, 50)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density')

    plt.tight_layout()
    plt.savefig('tex/figs/G_xx.pdf')
    plt.show()
#plot1b()

###############################################################################
# Calculate and plot standard narrowband sound pressure level
Pref = 20E-6

T = len(df) / fs
SPLf = 10 * np.log10((Gxx / T) / Pref ** 2)
df = pd.DataFrame([f, SPLf]).T
df.drop(0, inplace=True)
df.columns = ['Frequency', 'SPL']


def plot1c():
    plt.semilogx(f, SPLf, 'k')

    plt.xlim(f[1], f[-1])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Narrowband Sound Pressure Level')

    plt.tight_layout()
    plt.savefig('tex/figs/narrowband.pdf')
    plt.show()
#plot1c()

###############################################################################
# Problem 2
###############################################################################


def octaves(third=False, upper_frequency=10E3):
    dx = 1.
    if third:
        dx = 1 / 3

    fc30 = 1000
    m = np.arange(1, 80)

    center = fc30 * 2 ** (-10 + (m * dx))
    upper = center * 2 ** (dx / 2)
    lower = center / 2 ** (dx / 2)

    freqs = pd.DataFrame([lower, center, upper]).T
    freqs.columns = ['lower', 'center', 'upper']

    return freqs.query('upper < @upper_frequency')

o = octaves(third=True)
