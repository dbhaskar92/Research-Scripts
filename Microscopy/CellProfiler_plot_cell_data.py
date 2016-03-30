#!/usr/bin/env python

#
# Last modified: 2 Feb 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# 
import re
import sys
import collections
import numpy as NP
import matplotlib.pyplot as PLT

if (len(sys.argv) != 2):
	print 'Usage: python CellProfiler_plot_statistics.py /path/to/AllMyExpt_MyCells.csv'
	sys.exit()

data = NP.genfromtxt(sys.argv[1], delimiter=',', skip_header=1, skip_footer=0, usecols=(1,4,8,24), names=['obj_num','f_num','area','perimeter'])

f_area = zip(data['f_num'],data['area'])
areaMap = collections.defaultdict(list)
for (frame,area) in f_area:
	areaMap[int(frame)].append(int(area))
	
frames = areaMap.keys()
frames.sort()
areaStat = collections.defaultdict(list)
caption = "Number of data points: "+`len(frames)`
for frame in frames:
	areaStat[frame] = [0.64*NP.percentile(areaMap[frame],50),0.64*NP.percentile(areaMap[frame],90),0.64*NP.percentile(areaMap[frame],10)]
	
PLT.figure(1)
PLT.title('Cell Area Statistics')
PLT.xlabel('Frame/Image Number')
PLT.ylabel(r'Area (${\mu m}^2$)')
line1, = PLT.plot(frames, [areaStat[t][0] for t in frames], color='b', marker='.', linestyle='o', label='Median Area')
line2, = PLT.plot(frames, [areaStat[t][1] for t in frames], color='g', marker='.', linestyle='o', label='Q90 Area')
line3, = PLT.plot(frames, [areaStat[t][2] for t in frames], color='r', marker='.', linestyle='o', label='Q10 Area')
PLT.legend([line1, line2, line3], ['Q50', 'Q90', 'Q10'])
PLT.text(5,100,caption)
PLT.savefig('CellProfiler_Cells_plot1.png')
del frames
del f_area
del areaStat
del areaMap

f_perimeter = zip(data['f_num'],data['perimeter'])
pMap = collections.defaultdict(list)
for (frame,perimeter) in f_perimeter:
	pMap[int(frame)].append(float(perimeter))
	
frames = pMap.keys()
pStat = collections.defaultdict(list)
frames.sort()
caption = "Number of data points: "+`len(frames)`
for frame in frames:
	pStat[frame] = [0.8*NP.percentile(pMap[frame],50), 0.8*NP.percentile(pMap[frame],90), 0.8*NP.percentile(pMap[frame],10)]


PLT.figure(2)
PLT.title('Cell Perimeter Statistics')
PLT.xlabel('Frame/Image Number')
PLT.ylabel(r'Perimeter ($\mu m$)')
line1, = PLT.plot(frames, [pStat[t][0] for t in frames], color='b', marker='.', linestyle='o', label='Median Perimeter')
line2, = PLT.plot(frames, [pStat[t][1] for t in frames], color='g', marker='.', linestyle='o', label='Q90 Perimeter')
line3, = PLT.plot(frames, [pStat[t][2] for t in frames], color='r', marker='.', linestyle='o', label='Q10 Perimeter')
PLT.legend([line1, line2, line3], ['Q50', 'Q90', 'Q10'])
PLT.text(5,25,caption)
PLT.tight_layout()
PLT.savefig('CellProfiler_Cells_plot2.png')
del frames
del f_perimeter
del pStat
del pMap
