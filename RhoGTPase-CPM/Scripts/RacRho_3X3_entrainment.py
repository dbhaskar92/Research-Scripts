#!/usr/bin/env python

#
# Quantify entrainment in Rho GTPase oscillators over time
# Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
# Last modified: 28 June 2017
#

import sys
import matplotlib
import numpy as np
from matplotlib import rc
from lxml import etree as et
import matplotlib.pyplot as plt

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
matplotlib.rcParams.update({'font.size': 16})

if len(sys.argv) == 1:
    print '\nError: No XML file path specified.\n'
    exit()

tree = et.parse(sys.argv[1])
root = tree.getroot()
timeframe_elements = root.findall(".//time")

# Building data array
# rows = timeframes, columns = cells
data = 0  
timeframe_num = len(timeframe_elements)
cell_num = len(timeframe_elements[0].findall(".//cell"))

for timeframe in timeframe_elements:

    cell_elements = timeframe.findall(".//cell")
    timeframe_data = np.empty(cell_num)

    for cell in cell_elements:
        timeframe_data[cell_elements.index(cell)] = cell.get('gtpase')
    if timeframe_elements.index(timeframe) == 0:
        data = [timeframe_data]
    else:
        data = np.concatenate((data, [timeframe_data]))

# Plot
plt.subplot(211)
plt.plot(data)
plt.xlabel('Time')
plt.ylabel(r'Rho GTPase, $G$')
plt.grid(True)

# rows = cells, columns = timeframes
data = np.transpose(data)  
(num_cells, timelimit) = np.shape(data)

# Get rid of DC component
for cell in np.arange(0, num_cells, 1):
	data[cell] = data[cell] - np.mean(data[cell])

# Perform windowed FFT to compute phase difference in signals over time

# Calculate number of possible shifts of the window to the right, given the window size
windowsize = 120
offset = 0
num_shifts = timelimit - windowsize - offset
 
phase = np.empty([num_cells, num_shifts])

for cell in np.arange(0, num_cells, 1):

	cell_phi = np.empty(num_shifts)

	for shift in np.arange(0, num_shifts, 1):
		
		gtpase_window = data[cell, offset+shift:offset+shift+windowsize]
	
		gtpase_fft = np.fft.rfft(gtpase_window)

		# Find dominant frequency and record corresponding index
		indx = np.argmax(np.absolute(gtpase_fft))

		# Calculate phase corresponding to dominant frequency
		cell_phi[shift] = np.angle(gtpase_fft)[indx]
		
	phase[cell] = cell_phi

# Kuramoto order parameter 
R = np.empty(num_shifts)

# Variance of phases
V = np.empty(num_shifts)

for shift in np.arange(0, num_shifts, 1):
	cell_phases = np.empty(num_cells, dtype='complex128')
	for cell in np.arange(0, num_cells, 1):
		cell_phases[cell] = np.exp(1j * phase[cell, shift])
	z = np.divide(np.sum(cell_phases), num_cells)
	R[shift] = np.absolute(z)
	V[shift] = np.var(cell_phases)

plt.subplot(223)	
plt.plot(np.arange(0, num_shifts, 1), R, '-')
plt.xlabel('Time')
plt.ylabel(r'Kuramoto Order Parameter, $R$')
plt.ylim([0.2,1.0])

vector_a = np.exp(1j * phase[0]);
vector_b = np.exp(1j * phase[1]);
phase_diffs = np.absolute(vector_b - vector_a);

plt.subplot(224)
plt.plot(np.arange(0, num_shifts, 1), V, '-')
plt.xlabel('Time')
plt.ylabel(r'Variance($\phi_{i}$)', fontsize=18)
plt.ylim([0.0,0.8])

plt.show()
