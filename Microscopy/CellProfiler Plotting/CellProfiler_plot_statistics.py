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

if (len(sys.argv) != 3):
	print 'Usage: /usr/bin/python2 CellProfiler_plot_statistics.py /path/to/AllMyExpt_Image.csv /path/to/AllMyExpt_MyCells.csv'
	sys.exit()
	
# Parse AllMyExpt_Image.csv

data = NP.genfromtxt(sys.argv[1], delimiter=',', skip_header=1, skip_footer=0, usecols=(1,24,48,49,50,51), names=['cell_count','f_num','lost_obj','merged_obj','new_obj','split_obj'])

totalframes = len(data['f_num'])

# Calculate number of frames with zero detected cells
 
lostframes = 0
for num_cells in data['cell_count']:
	if num_cells == 0:
		lostframes = lostframes + 1
remaining_frames = totalframes - lostframes

caption = "Total frames: " + `totalframes` + " Lost frames: " + `lostframes` + " Remaining frames: " + `remaining_frames`

# Plot cell count

PLT.figure(1)
PLT.title('Cell Count')
PLT.plot(data['f_num'], data['cell_count'], color='b', linestyle='-')
PLT.xlabel('Frame/Image Number')
PLT.ylabel('Number of Cells')
#PLT.text(5, 10, caption)
axes = PLT.gca()
axes.set_ylim([100, 700])
axes.set_xlim([0, 480])
PLT.savefig('CellProfiler_cell_count.png')

# Plot number of splits, merges and lost object count

PLT.figure(2)
PLT.title('CellProfiler Statistics')
line1, = PLT.plot(data['f_num'], data['lost_obj'], color='r', linestyle='-', label='Lost Objects')
line2, = PLT.plot(data['f_num'], data['merged_obj'], color='g', linestyle='-', label='Merged Objects')
line3, = PLT.plot(data['f_num'], data['split_obj'], color='b', linestyle='-', label='Split Objects')
PLT.legend([line1, line2, line3], ['Lost Objects', 'Merged Objects', 'Split Objects'])
PLT.xlabel('Frame/Image Number')
PLT.ylabel('Object Count')
axes = PLT.gca()
axes.set_xlim([0, 480])
axes.set_ylim([-5, 300])
PLT.savefig('CellProfiler_obj_detection_stats.png')

# Parse AllMyExpt_MyCells.csv 

data = NP.genfromtxt(sys.argv[2], delimiter=',', skip_header=1, skip_footer=0, usecols=(1,4,8,24), names=['obj_num','f_num','area','perimeter'])

# Calculate area statistics

f_area = zip(data['f_num'],data['area'])
areaMap = collections.defaultdict(list)
for (frame,area) in f_area:
	areaMap[int(frame)].append(int(area))
	
frames = areaMap.keys()
frames.sort()
areaStat = collections.defaultdict(list)
caption = "Number of data points: " + `len(frames)`
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

f_perimeter = zip(data['f_num'],data['perimeter'])
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
PLT.savefig('CellProfiler_Cells_plot_area_perim.png', dpi=120)
