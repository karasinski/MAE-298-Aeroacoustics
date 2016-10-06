import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io.wavfile
import scipy.signal

###############################################################################
# Problem 1
###############################################################################

# Load the audio file
fs, audio = scipy.io.wavfile.read('Boom_F1B2_6.wav')
df = pd.DataFrame(audio)
df.columns = ['Voltage']
df['Time'] = np.array(df.index) / fs

# Convert to pascals
df['Pascals'] = df.Voltage * -116.

# Plot the pressure as a function of time
df.plot(x='Time', y='Pascals')
plt.show()

# Find the peak pressure
print(df.Pascals.max())

###############################################################################
# Calculate the single sided power spectral density function
f, Gxx = scipy.signal.periodogram(df.Pascals, fs=fs, return_onesided=True)

# Plot G_xx
plt.plot(f, Gxx)
plt.xlim(0, 50)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power Spectral Density')
plt.show()

###############################################################################
# Calculate and plot standard narrowband sound pressure level
df['dB'] = 10 * np.log10(0.5 * df.Pascals**2 / (2E-5)**2)
df.plot(x='Time', y='dB')
plt.show()

###############################################################################
# Problem 2
###############################################################################
