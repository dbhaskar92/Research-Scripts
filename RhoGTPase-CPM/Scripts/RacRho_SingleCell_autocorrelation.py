#!/usr/bin/env python

#
# Compute frequency of Rho GTPase oscillators
# Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
# Last modified: 29 June 2017
#

import sys
import numpy as np
from lxml import etree as et
import matplotlib.pyplot as plt

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
plt.plot(data)
plt.xlabel('Time')
plt.ylabel('GTPase')
plt.grid(True)
plt.show()

# rows = cells, columns = timeframes
data = np.transpose(data)  
(num_cells, timelimit) = np.shape(data)

# Get rid of DC component
for cell in np.arange(0, num_cells, 1):
	data[cell] = data[cell] - np.mean(data[cell])

fft_res = np.fft.rfft(data[0])

x_range = len(fft_res);

# Frequency plot
frequencies = np.empty(x_range)

for i in range(x_range):

	# to convert fft indices to frequency values:
	# divide by number of indices and multiply by # of samples per MCS (freq)
	# ODE solver params: 'stepSize':0.001, 'steps':10, 'stiff':True
	# SBML time scale param: alpha = 2000
	frequencies[i] = float(float(i)/timelimit)*(1/(2000*0.001));
		
plt.plot(frequencies, np.absolute(fft_res))
# The following are equivalent: np.sqrt(np.multiply(fft_res, np.conj(fft_res))) and np.absolute(fft_res)

plt.xlabel('Frequency (Hz)')
plt.ylabel('Amplitude')
plt.title('FFT of Rho GTPase Activity')
plt.xlim([0.0,0.2])
plt.grid(True)
plt.show()

# Find index corresponding to dominant frequency (i.e. frequency with largest amplitude)
indx = np.argmax(np.absolute(fft_res));
print indx

# Print result
freq = float(float(indx)/timelimit)*(1/(2000*0.001));
print freq
