#!/usr/bin/env python

#
# Compute phase difference between 2 signals using FFT
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
# signal with frequency 0.5Hz at phase offset of 1.57 (\pi/2)
b = np.sin(2*np.pi*0.5*t - 1.57) + np.random.normal(0, 0.5, 500);

# plot test data
plt.plot(t, a, t, b)
plt.xlabel('Time (s)')
plt.title('Signals')
plt.show()

fft_a = np.fft.rfft(a);
fft_b = np.fft.rfft(b);

# find dominant frequency of both signals and record corresponding index
indx_a = np.argmax(np.absolute(fft_a));
indx_b = np.argmax(np.absolute(fft_b));

# calculate phase corresponding to dominant frequency
phi_a = np.angle(fft_a)[indx_a];
phi_b = np.angle(fft_b)[indx_b];

# plot phases on a circle
vector_a = np.exp(1j * phi_a);
vector_b = np.exp(1j * phi_b);
# make a circle, by plotting e^{i\theta} for \theta in \[0,2\pi\]
pts = np.exp(1j * np.arange(0, 2*np.pi, 0.001)); 
plt.plot(pts.real, pts.imag, '-')
plt.plot(vector_a.real, vector_a.imag, '*')
plt.plot(vector_b.real, vector_b.imag, '*')
plt.xlim(-1.2, 1.2)
plt.ylim(-1.2, 1.2)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

# print result
if np.absolute(phi_a - phi_b) > np.pi:
	print (np.absolute(phi_a - phi_b))/np.pi;
else:
 print np.absolute(phi_a - phi_b)
