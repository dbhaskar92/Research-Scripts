#!/usr/bin/env python

#
# Last modified: 3 Jun 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
#

import re
import sys
import collections
import numpy as NP
import matplotlib.pyplot as PLT

if (len(sys.argv) != 2):
	print 'Usage: /usr/bin/python2 plot_MATLAB_statistics.py /PATH/TO/CSV/FILE'
	sys.exit()

data = NP.genfromtxt(sys.argv[1], delimiter=',', skip_header=1, skip_footer=0, usecols=(1,0,4,5), names=['CellID','Frame','Area','Perimeter'])

# Calculate area statistics

f_area = zip(data['Frame'],data['Area'])
areaMap = collections.defaultdict(list)
for (frame,area) in f_area:
	areaMap[int(frame)].append(int(area))
	
frames = areaMap.keys()
frames.sort()
areaStat = collections.defaultdict(list)
caption = "Number of data points: "+`len(frames)`
for frame in frames:
	areaStat[frame] = [0.64*NP.percentile(areaMap[frame],50),0.64*NP.percentile(areaMap[frame],90),0.64*NP.percentile(areaMap[frame],10)]

# Make area plot

fig, ax1 = PLT.subplots()
line1, = ax1.plot(frames, [areaStat[t][0] for t in frames], color='b', linestyle='-')
ax1.set_xlabel('Frame/Image Number')
ax1.set_ylabel(r'Area (${\mu m}^2$)', color='b')
ax1.set_ylim([50,450])
for tl in ax1.get_yticklabels():
    tl.set_color('b')

# Calculate perimeter statistics

f_perimeter = zip(data['Frame'],data['Perimeter'])
pMap = collections.defaultdict(list)
for (frame,perimeter) in f_perimeter:
	pMap[int(frame)].append(float(perimeter))
	
frames = pMap.keys()
pStat = collections.defaultdict(list)
frames.sort()
for frame in frames:
	pStat[frame] = [0.8*NP.percentile(pMap[frame],50),0.8*NP.percentile(pMap[frame],90),0.8*NP.percentile(pMap[frame],10)]

# Make perimeter plot

ax2 = ax1.twinx()
line2, = ax2.plot(frames, [pStat[t][0] for t in frames], color='r', linestyle='-')
ax2.set_ylabel(r'Perimeter ($\mu m$)', color='r')
ax2.set_ylim([0,100])
ax2.set_xlim([0,480])
for tl in ax2.get_yticklabels():
    tl.set_color('r')
#PLT.legend([line1, line2], ['Median Area', 'Median Perimeter'])
#PLT.text(5, 5, caption)
PLT.savefig('MATLAB_Cells_plot_area_perim.png', dpi=120)

# Plot cell count

f_cell = zip(data['Frame'],data['CellID'])
cellcntMap = collections.defaultdict(list)
for (frame,cid) in f_cell:
	if frame not in cellcntMap.keys():
		cellcntMap[int(frame)] = 0
	else:	
		cellcntMap[int(frame)] += 1

frames = cellcntMap.keys()
PLT.figure()
PLT.title('Cell Count')
PLT.plot(frames, [cellcntMap[t] for t in frames], color='b', linestyle='-')
PLT.xlabel('Frame/Image Number')
PLT.ylabel('Number of Cells')
PLT.xlim([0,480])
PLT.ylim([100,700])
PLT.savefig('MATLAB_cell_count.png')
