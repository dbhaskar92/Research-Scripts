#!/usr/bin/env python

#
# Compute frequency of a signal using fourier transform
# Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
# Last modified: 28 June 2017
#

import numpy as np
import matplotlib.pyplot as plt

# sampling rate: 0.02
# sampling duration: 10
t = np.arange(0, 10, 0.02);

a = np.empty(500);
b = np.empty(500);

# signal with frequency 0.5Hz
a = np.sin(2*np.pi*0.5*t) + np.random.normal(0, 0.5, 500);
# signal with frequency 2Hz
b = np.sin(2*np.pi*2*t) + np.random.normal(0, 0.5, 500);
# Specify test signal
sig = b

# plot test data
plt.plot(t, sig)
plt.xlabel('Time (s)')
plt.title('Signal')
plt.show()

fft_res = np.fft.rfft(sig);

# length of result is 251 (corresponding to time indices 0-250)
x_range = len(fft_res);

# frequency plot
frequencies = np.empty(x_range)
for i in range(x_range):
	# to convert fft indices to frequency values:
	# divide by number of indices (500) and multiply by # of samples per second (freq)
	frequencies[i] = (float(i)/500)*(1/0.02);
		
plt.plot(frequencies, np.absolute(fft_res))
# The following are equivalent: np.sqrt(np.multiply(fft_res, np.conj(fft_res))) and np.absolute(fft_res)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('FFT of Signal')
plt.show()

# find index corresponding to dominant frequency (i.e. frequency with largest amplitude)
indx = np.argmax(np.absolute(fft_res));

# Print result
freq = (float(indx)/500)*(1/0.02);
print freq
