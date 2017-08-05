#!/usr/bin/env python

#
# Last modified: 13 Feb 2017
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# 

import sys
import re
import math
import collections
import numpy as NP
import matplotlib.pyplot as PLT
import lxml.etree as ET

if (len(sys.argv) != 2):
	print 'Usage: /usr/bin/python2 single_cell_gtpase_stats.py /path/to/sim_sml_gtpase.xml'
	sys.exit()
	
sim = ET.parse(open(sys.argv[1],'r'))

root = sim.getroot()
cell_type_list = re.split(',',root.get('type_list'))
num_cell_types = len(cell_type_list)

areaStat = collections.defaultdict(list)
perimeterStat = collections.defaultdict(list)
tareaStat = collections.defaultdict(list)
pressureStat = collections.defaultdict(list)
gtpaseStat = collections.defaultdict(list)

for time in sim.getiterator('time'):

	numcells = len(time)
	
	if numcells != 1:
		print "Error: More than one cell."
	
	mcs = int(time.get('t'))
	
	
	for cell in time.getchildren():
	
		if int(cell.get('cell_id')) != 1:
			print "Error: Cell ID not equal to 1."
			exit
	
		areaStat[mcs] = float(cell.get('area'))
		perimeterStat[mcs] = float(cell.get('perimeter'))
		tareaStat[mcs] = float(cell.get('target_volume'))
		pressureStat[mcs] = float(cell.get('pressure'))
		gtpaseStat[mcs] = float(cell.get('gtpase'))

# Plot Figures

numFigures = 0
timeseries = areaStat.keys()
timeseries.sort()

numFigures += 1
PLT.figure()
PLT.subplot(311)
PLT.title('Cell Area and Target Area')
line1, = PLT.plot(timeseries, [areaStat[t] for t in timeseries], color='b', marker='.', linestyle='-', label='Cell Area')
line2, = PLT.plot(timeseries, [tareaStat[t] for t in timeseries], color='g', marker='.', linestyle='-', label='Target Area')
PLT.legend([line1, line2], ['Area', 'Target Area'])

PLT.subplot(312)
PLT.title('GTPase Activity')
PLT.plot(timeseries, [gtpaseStat[t] for t in timeseries], color='k', marker='.', linestyle='-', label='GTPase')

PLT.subplot(313)
PLT.title('Pressure (Target Area - Area)')
PLT.plot(timeseries, [pressureStat[t] for t in timeseries], color='b', marker='.', linestyle='-', label='Median Pressure')
PLT.xlabel('Monte-Carlo Steps (MCS)')
PLT.tight_layout()
PLT.savefig('Plot1.png')	
# PLT.show()

numFigures += 1
PLT.figure(numFigures)
PLT.title('Cell Perimeter')
PLT.ylabel('# pixels')
PLT.xlabel('Monte-Carlo Steps (MCS)')
PLT.plot(timeseries, [perimeterStat[t] for t in timeseries], color='k', marker='.', linestyle='-')
PLT.savefig('Plot2.png')	
# PLT.show()
