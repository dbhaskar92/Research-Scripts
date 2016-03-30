#!/usr/bin/env python

#
# Last modified: 31 Jan 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# 
import re
import sys
import numpy as NP
import matplotlib.pyplot as PLT

if (len(sys.argv) != 2):
	print 'Usage: python CellProfiler_plot_statistics.py /path/to/AllMyExpt_Image.csv'
	sys.exit()

data = NP.genfromtxt(sys.argv[1], delimiter=',', skip_header=1, skip_footer=0, usecols=(1,24,48,49,50,51), names=['cell_count','f_num','lost_obj','merged_obj','new_obj','split_obj'])

totalframes = len(data['f_num'])
# Lost frames due to lamp 
lostframes = 0
for num_cells in data['cell_count']:
	if num_cells == 0:
		lostframes = lostframes + 1
remaining_frames = totalframes - lostframes
caption = "Total frames: "+`totalframes`+" Lost frames: "+`lostframes`+" Remaining frames: "+`remaining_frames`

PLT.figure(1)
PLT.title('Identified Cell Count')
PLT.plot(data['f_num'], data['cell_count'], color='b', marker='.', linestyle='o')
PLT.xlabel('Frame/Image Number')
PLT.ylabel('Cell Count')
PLT.text(5,100,caption)
axes = PLT.gca()
axes.set_ylim([-50, 1200])
PLT.savefig('CellProfiler_plot1.png')

PLT.figure(2)
PLT.title('Cell Profiler Statistics')
line1, = PLT.plot(data['f_num'], data['lost_obj'], color='r', marker='.', linestyle='o', label='Lost Objects')
line2, = PLT.plot(data['f_num'], data['merged_obj'], color='g', marker='.', linestyle='o', label='Merged Objects')
line3, = PLT.plot(data['f_num'], data['split_obj'], color='b', marker='.', linestyle='o', label='Split Objects')
PLT.legend([line1, line2, line3], ['Lost Objects', 'Merged Objects', 'Split Objects'])
PLT.xlabel('Frame/Image Number')
PLT.ylabel('Cell Count')
axes = PLT.gca()
axes.set_xlim([0, 1500])
axes.set_ylim([-50, 1200])
PLT.savefig('CellProfiler_plot2.png')
