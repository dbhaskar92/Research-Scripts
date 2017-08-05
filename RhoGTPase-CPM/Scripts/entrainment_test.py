#!/usr/bin/env python

#
# Quantify entrainment over time between 2 signals using FFT
# Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
# Last modified: 28 June 2017
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# sampling rate: 0.02
# sampling duration: 10
t = np.arange(0, 10, 0.02);

a = np.empty(500);
b = np.empty(500);

# phase change function
# zero phase difference (entrained) at t=5, maximum phase difference at t=0 and t=10
delta_phi = np.divide(np.multiply((5 - t), 1.57), 5);

# signal with frequency 0.5Hz
a = np.sin(2*np.pi*0.5*t) + np.random.normal(0, 0.5, 500);
# signal with frequency 0.5Hz at variable phase offset
b = np.sin(2*np.pi*0.5*t - delta_phi) + np.random.normal(0, 0.5, 500);

# plot test data
plt.plot(t, a, t, b)
plt.xlabel('Time (s)')
plt.title('Signals')
plt.show()

# perform windowed FFT to compute phase difference in signals over time
windowsize = 100;
# 1 second = (1/10)*500 time steps
num_shifts = 500 - windowsize;
# number of possible shifts of the window to the right, given the window size 
phi_a = np.empty(num_shifts);
phi_b = np.empty(num_shifts);

for shift in np.arange(0, num_shifts, 1):

	a_window = a[shift:shift+windowsize];
	b_window = b[shift:shift+windowsize];
	
	fft_a = np.fft.rfft(a_window);
	fft_b = np.fft.rfft(b_window);

	# find dominant frequency of both signals and record corresponding index
	indx_a = np.argmax(np.absolute(fft_a));
	indx_b = np.argmax(np.absolute(fft_b));

	# calculate phase corresponding to dominant frequency
	phi_a[shift] = np.angle(fft_a)[indx_a];
	phi_b[shift] = np.angle(fft_b)[indx_b];

vector_a = np.exp(1j * phi_a);
vector_b = np.exp(1j * phi_b);
phase_diffs = np.absolute(vector_b - vector_a);
plt.plot(np.arange(0, num_shifts, 1), phase_diffs, '-')
plt.xlabel('Window Shift')
plt.ylabel('Phase Difference')
plt.show()

# Kuramoto order parameter
R = np.empty(num_shifts)
for shift in np.arange(0, num_shifts, 1):
	total = vector_a[shift] + vector_b[shift]
	z = total / 2
	R[shift] = np.absolute(z)	
plt.plot(np.arange(0, num_shifts, 1), R, '-')
plt.xlabel('Window Shift')
plt.ylabel('R')
plt.title('Kuramoto order parameter')
plt.show()
