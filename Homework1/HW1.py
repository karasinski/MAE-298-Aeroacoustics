import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal
import seaborn as sns
import soundfile as sf

sns.reset_orig()
colors = sns.color_palette('deep', 10)

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
    df.plot(x='Time', y='Pascals', legend=False, color=colors[0], ax=ax)

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

    plt.vlines(x=maxPeak.Time, ymax=0, ymin=minPeak.Pascals, linestyles='--')

    plt.tight_layout()
    plt.savefig('tex/figs/Pascals_vs_Time.pdf')
    plt.clf()
#plot1a()

################################################################################
# Calculate the single sided power spectral density function
f, Gxx = scipy.signal.periodogram(df.Pascals, fs=fs, return_onesided=True)


def plot1b():
    plt.plot(f, Gxx, marker='o', color=colors[0])

    plt.xlim(0, 50)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Power Spectral Density')

    plt.tight_layout()
    plt.savefig('tex/figs/G_xx.pdf')
    plt.clf()
#plot1b()

###############################################################################
# Calculate and plot standard narrowband sound pressure level
Pref = 20E-6

T = len(df) / fs
SPL = 10 * np.log10((Gxx / T) / Pref ** 2)
df = pd.DataFrame([f, SPL]).T
df.drop(0, inplace=True)
df.columns = ['Frequency', 'SPL']


def plot1c():
    df.plot(x='Frequency', y='SPL', logx=True, color=colors[0])

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Narrowband Sound Pressure Level (dB)')

    plt.tight_layout()
    plt.savefig('tex/figs/narrowband.pdf')
    plt.clf()
#plot1c()

###############################################################################
# Problem 2
###############################################################################


def generate_octaves(octave=1., upper_frequency=10E3):
    fc30 = 1000
    m = np.arange(1, 80)

    center = fc30 * 2 ** (-10 + (m * octave))
    upper = center * 2 ** (octave / 2)
    lower = center / 2 ** (octave / 2)

    freqs = pd.DataFrame([lower, center, upper]).T
    freqs.columns = ['lower', 'center', 'upper']

    return freqs.query('upper < @upper_frequency')


def sum_octaves(octaves, octave=1.):
    octave = generate_octaves(octave, upper_frequency=df.Frequency.max())
    res = []
    for i, band in octave.iterrows():
        b = octaves.query('@band.lower < Frequency < @band.upper')
        Lp = 10 * np.log10((10 ** (b.SPL / 10)).sum())
        if Lp > 0:
            res.append([band.center, Lp])

    res = pd.DataFrame(res)
    res.columns = ['Frequency', 'SPL']
    return res


third = sum_octaves(df, octave=1 / 3)
full = sum_octaves(third, octave=1.)
overall = 10 * np.log10((10 ** (full.SPL / 10)).sum())


def plot2():
    f, ax = plt.subplots()
    df.plot(x='Frequency', y='SPL', logx=True,
            color=colors[0], legend=False, ax=ax)
    third.plot(x='Frequency', y='SPL', logx=True,
               marker='o', color=colors[1], legend=False, ax=ax)
    full.plot(x='Frequency', y='SPL', logx=True,
              marker='o', color=colors[2], legend=False, ax=ax)

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Sound Pressure Level (dB)')
    plt.legend(['Narrowband', 'One Third', 'Octave'])
    plt.title('Overall SPL: {:2.2f} dB'.format(overall))

    plt.tight_layout()
    plt.savefig('tex/figs/octaves.pdf')
    plt.clf()
#plot2()
