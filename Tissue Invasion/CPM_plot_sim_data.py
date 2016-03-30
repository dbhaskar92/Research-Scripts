#!/usr/bin/env python

#
# Last modified: 15 Aug 2015
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# 

import sys
import re
import math
import collections
import matplotlib.pyplot as PLT
import lxml.etree as ET

if (len(sys.argv) != 2):
	print 'Usage: python plot_sim_data.py /path/to/simulation_results.xml'
	sys.exit()
	
sim = ET.parse(open(sys.argv[1],'r'))

root = sim.getroot()
cell_type_list = re.split(',',root.get('type_list'))
num_cell_types = len(cell_type_list)

cellNumStat = collections.defaultdict(list)
areaStat = collections.defaultdict(list)
perimeterStat = collections.defaultdict(list)
pressureStat = collections.defaultdict(list)
neighborStat = collections.defaultdict(list)

cellTracking = collections.defaultdict(list)

for time in sim.getiterator('time'):
	numcells = len(time)
	mcs = int(time.get('t'))
	
	for cell_type in cell_type_list:
		cellNumStat[mcs].append(0)
	
	cellDict = collections.defaultdict(list)
	for cell in time.getchildren():
		cellDict[cell.get('cell_id')].append(int(cell.get('type')))
		cellDict[cell.get('cell_id')].append(float(cell.get('x')))
		cellDict[cell.get('cell_id')].append(float(cell.get('y')))
		cellDict[cell.get('cell_id')].append(float(cell.get('area')))
		cellDict[cell.get('cell_id')].append(float(cell.get('perimeter')))
		cellDict[cell.get('cell_id')].append(float(cell.get('pressure')))
		cellDict[cell.get('cell_id')].append(float(cell.get('contact_perimeter')))
		cellDict[cell.get('cell_id')].append(len(re.split('\s+',cell.get('neighbors'))))
		
		case = int(cell.get('type'))
		if case == 1:
			cellTracking[cell.get('cell_id')].append(['g',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][0] += 1 
		elif case == 2:
			cellTracking[cell.get('cell_id')].append(['r',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][1] += 1
		elif case == 3:
			cellTracking[cell.get('cell_id')].append(['b',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][2] += 1
		else:
			cellTracking[cell.get('cell_id')].append(['k',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][3] += 1
			
	avg_area = sum ( cellDict[cell][3] for cell in cellDict ) / len(cellDict)
	avg_perimeter = sum ( cellDict[cell][4] for cell in cellDict ) / len(cellDict)
	avg_pressure = sum ( cellDict[cell][5] for cell in cellDict ) / len(cellDict)
	avg_contact_perimeter = sum ( cellDict[cell][6] for cell in cellDict ) / len(cellDict)
	avg_neighbours = sum ( cellDict[cell][7] for cell in cellDict ) / len(cellDict) 
	
	max_area = max ( cellDict[cell][3] for cell in cellDict )
	max_perimeter = max ( cellDict[cell][4] for cell in cellDict )
	max_pressure = max ( cellDict[cell][5] for cell in cellDict )
	max_contact_perimeter = max ( cellDict[cell][6] for cell in cellDict )
	max_neighbours = max ( cellDict[cell][7] for cell in cellDict )  
	
	min_area = min ( cellDict[cell][3] for cell in cellDict )
	min_perimeter = min ( cellDict[cell][4] for cell in cellDict )
	min_pressure = min ( cellDict[cell][5] for cell in cellDict )
	min_contact_perimeter = min ( cellDict[cell][6] for cell in cellDict )
	min_neighbours = min ( cellDict[cell][7] for cell in cellDict )    
	
	areaStat[mcs] = [avg_area, max_area, min_area]
	perimeterStat[mcs] = [avg_perimeter, max_perimeter, min_perimeter]
	pressureStat[mcs] = [avg_pressure, max_pressure, min_pressure]
	neighborStat[mcs] = [avg_neighbours, max_neighbours, min_neighbours]


# Plot figures
numFigures = 0
timeseries = areaStat.keys()
timeseries.sort()


# Tissue Statistics Plot
numFigures += 1
PLT.figure()
PLT.subplot(311)
PLT.title('Area Statistics')
line1, = PLT.plot(timeseries, [areaStat[t][0] for t in timeseries], color='b', marker='.', linestyle='-', label='Avg. Area')
line2, = PLT.plot(timeseries, [areaStat[t][1] for t in timeseries], color='g', marker='.', linestyle='-', label='Max. Area')
line3, = PLT.plot(timeseries, [areaStat[t][2] for t in timeseries], color='r', marker='.', linestyle='-', label='Min. Area')
PLT.legend([line1, line2, line3], ['Avg.', 'Max.', 'Min.'])

PLT.subplot(312)
PLT.title('Perimeter Statistics')
PLT.plot(timeseries, [perimeterStat[t][0] for t in timeseries], color='b', marker='.', linestyle='-', label='Avg. Perimeter')
PLT.plot(timeseries, [perimeterStat[t][1] for t in timeseries], color='g', marker='.', linestyle='-', label='Max. Perimeter')
PLT.plot(timeseries, [perimeterStat[t][2] for t in timeseries], color='r', marker='.', linestyle='-', label='Min. Perimeter')

PLT.subplot(313)
PLT.title('Cell Pressure Statistics')
PLT.plot(timeseries, [pressureStat[t][0] for t in timeseries], color='b', marker='.', linestyle='-', label='Avg. Pressure')
PLT.plot(timeseries, [pressureStat[t][1] for t in timeseries], color='g', marker='.', linestyle='-', label='Max. Pressure')
PLT.plot(timeseries, [pressureStat[t][2] for t in timeseries], color='r', marker='.', linestyle='-', label='Min. Pressure')
PLT.xlabel('Monte-Carlo Steps (MCS)')
PLT.tight_layout()
PLT.savefig('cpm_plot1.png');	#PLT.show()
del areaStat
del perimeterStat
del pressureStat

# Neighbour Statistics Plot
numFigures += 1
PLT.figure(numFigures)
PLT.subplot(211)
PLT.title('Neighbour Statistics')
PLT.ylabel('Number of Neighbours')
line1, = PLT.plot(timeseries, [neighborStat[t][0] for t in timeseries], color='b', marker='.', linestyle='steps', label='Average')
line2, = PLT.plot(timeseries, [neighborStat[t][1] for t in timeseries], color='g', marker='.', linestyle='steps', label='Maximum')
line3, = PLT.plot(timeseries, [neighborStat[t][2] for t in timeseries], color='r', marker='.', linestyle='steps', label='Minimum')
[xmin, xmax, ymin, ymax] = PLT.axis()
PLT.axis([0,xmax,ymin-2,ymax+2])
PLT.legend([line1, line2, line3], ['Average', 'Maximum', 'Minimum'])

PLT.subplot(212)
if num_cell_types >= 2:
	heterotypicboundary = collections.defaultdict(list)
	for time in sim.getiterator('time'):
		mcs = int(time.get('t'))
		heterotypicboundary[mcs] = float(time.get('heterotypic_boundary_length'))

	PLT.title('Heterotypic Boundary Length Between Cell Types: ' + cell_type_list[0] + ' and ' + cell_type_list[1])
	PLT.ylabel('Boundary Length (# Lattice sites)')
	PLT.plot(timeseries, [heterotypicboundary[t] for t in timeseries], color='k', marker='.', linestyle='-')
	PLT.xlabel('Monte-Carlo Steps (MCS)')
PLT.tight_layout()
PLT.savefig('cpm_plot2.png')	#PLT.show()
del neighborStat
del heterotypicboundary


# Cell COM Tracking Plot
numFigures += 1	    
PLT.figure(numFigures)
PLT.title('Cell Centroid Tracking')
xdata = []
ydata = []
xmax, ymax = 0, 0
for cell_id in cellTracking:
	clr = ''
	for com in cellTracking[cell_id]:
		if (clr != '' and clr != com[0]):
			# check for cell type changes
			PLT.plot(xdata, ydata, color=clr, marker='.', linestyle='o')
			xdata = [com[1]]
			ydata = [com[2]]
			xmax = max(xmax, com[1])
			ymax = max(ymax, com[2])
			clr = ''
		else:
			clr = com[0]	
			xdata.append(com[1])
			ydata.append(com[2])
			xmax = max(xmax, com[1])
			ymax = max(ymax, com[2])
	PLT.plot(xdata, ydata, color=clr, marker='.', linestyle='o')
	xdata = []
	ydata = []
PLT.xlabel('X Lattice Sites')
PLT.ylabel('Y Lattice Sites')
PLT.axis([0,xmax+5,0,ymax+5])
PLT.grid(True)
PLT.savefig('cpm_plot3.png')	#PLT.show()


# Cell Type Statistics
numFigures += 1
PLT.figure(numFigures)
PLT.subplot(211)
PLT.title('Cell Count Statistics')
PLT.ylabel('Number of Cells')
PLT.xlabel('Monte-Carlo Steps (MCS)')
line1, = PLT.plot(timeseries, [cellNumStat[t][0] for t in timeseries], color='g', marker='.', linestyle='-', label=cell_type_list[0])
line2, = PLT.plot(timeseries, [cellNumStat[t][1] for t in timeseries], color='r', marker='.', linestyle='steps', label=cell_type_list[1])
line3, = PLT.plot(timeseries, [sum(cellNumStat[t]) for t in timeseries], color='k', marker='.', linestyle='steps', label='Total')
[xmin, xmax, ymin, ymax] = PLT.axis()
PLT.axis([0,xmax,ymin-10,ymax+10])
PLT.legend([line1, line2, line3], [cell_type_list[0], cell_type_list[1], 'Total'])

cellType1Displacement = []
cellType2Displacement = []
tol = 50
for cell_id in cellTracking:
	xl, yl, xh, yh, ctype = -1, -1, -1, -1, -1
	displacement = 0
	for com in cellTracking[cell_id]:
		if xl == -1 and yl == -1:
			xl,xh = com[1], com[1]
			yl,yh = com[2], com[2]
			if com[0] == 'g':
				ctype = 1
			elif com[0] == 'r':
				ctype = 2
		if abs(com[1] - xh) > tol or abs(com[2] - yh) > tol:
			# check for periodic boundary conditions
			displacement += math.sqrt((xh - xl)**2 + (yh - yl)**2)
			xl,xh = com[1], com[1]
			yl,yh = com[2], com[2]
		else:
			xh = com[1]
			yh = com[2]
	if xh != -1 and yh != -1:
		displacement += math.sqrt((xh - xl)**2 + (yh - yl)**2)			
	if ctype == 1:
		cellType1Displacement.append(displacement)
	elif ctype == 2:
		cellType2Displacement.append(displacement)

PLT.subplot(212)
PLT.title('Cell Centroid Displacement Histogram')
PLT.xlabel('Displacement (length units)')
PLT.ylabel('Number of Cells')
PLT.hist([cellType1Displacement, cellType2Displacement], min(25,len(cellTracking)), normed=False, cumulative=False, color=['g','r'], label=[cell_type_list[0], cell_type_list[1]])
PLT.legend()
PLT.tight_layout()
PLT.savefig('cpm_plot4.png')	#PLT.show()
