#!/usr/bin/env python

#
# Last modified: 12 Feb 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# 

import sys
import re
import math
import collections
import numpy as NP
import lxml.etree as ET
import matplotlib.pyplot as PLT

if (len(sys.argv) != 6):
	print 'Usage: python plot_sim_data.py /path/to/simulation_results.xml /path/to/heterotypic.dat ...'
	print '/path/to/T1SwapLocations.dat /path/to/T2SwapLocations.dat /path/to/T3SwapLocations.dat'
	sys.exit()
	
def itertodt(numiterlist):
	xdata = []
	for element in numiterlist:
		xdata.append((0.02/10)*float(element))
	xdata.sort()	
	return xdata
	
sim = ET.parse(open(sys.argv[1],'r'))

root = sim.getroot()
cell_type_list = re.split(',',root.get('type_list'))
num_cell_types = len(cell_type_list)

cellNumStat = collections.defaultdict(list)
areaStat = collections.defaultdict(list)
perimeterStat = collections.defaultdict(list)
pressureStat = collections.defaultdict(list)
neighborStat = collections.defaultdict(list)
hetneighborStat = collections.defaultdict(list)
homneighborStat = collections.defaultdict(list)

cellTracking = collections.defaultdict(list)

for time in sim.getiterator('time'):
	numcells = len(time)
	mcs = int(time.get('tau'))
	tstep = float(time.get('t'))
	
	for cell_type in cell_type_list:
		cellNumStat[mcs].append(0)
	
	cellDict = collections.defaultdict(list)
	neighborDict = collections.defaultdict(list)
	cellHetNeighborTracking = collections.defaultdict(list)
	cellHomNeighborTracking = collections.defaultdict(list)
	
	for cell in time.getchildren():
		cell_id = int(cell.get('cell_id'))
		cellDict[cell_id].append(int(cell.get('type')))
		cellDict[cell_id].append(float(cell.get('x')))
		cellDict[cell_id].append(float(cell.get('y')))
		cellDict[cell_id].append(float(cell.get('area')))
		cellDict[cell_id].append(float(cell.get('perimeter')))
		cellDict[cell_id].append(float(cell.get('pressure')))
		cellDict[cell_id].append(float(cell.get('contact_perimeter')))
		cellDict[cell_id].append(len(str.split(cell.get('neighbors'))))
		
		case = int(cell.get('type'))
		if case == 1:
			cellTracking[cell_id].append(['g',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][0] += 1 
		elif case == 2:
			cellTracking[cell_id].append(['r',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][1] += 1
		elif case == 3:
			cellTracking[cell_id].append(['b',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][2] += 1
		else:
			cellTracking[cell_id].append(['k',float(cell.get('x')),float(cell.get('y'))])
			cellNumStat[mcs][3] += 1
		
		neighborlist = str.split(cell.get('neighbors'))
		for neighborId in neighborlist:
			neighborDict[cell_id].append(int(neighborId))
			
					
	for cid in cellDict.keys():
		for nid in neighborDict[cid]:
			try:
				ctype = cellDict[cid][0]
				ntype = cellDict[nid][0]
				if ctype != ntype:
					cellHetNeighborTracking[cid].append(nid)	
				else:
					cellHomNeighborTracking[cid].append(nid)
			except IndexError:
				print "Possible T3 Swap or Mesh Error. Iter: ",mcs," Time step: ",tstep," Cell id: ",cid," Neighbour id: ",nid
				continue
	
	area_list = []
	perimeter_list = []
	pressure_list = []
	contact_perimeter_list = []
	het_neighbour_count = []
	hom_neighbour_count = []
	num_neighbours = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	
	for cell in cellDict.keys():
		try:
			area_list.append(cellDict[cell][3])
			perimeter_list.append(cellDict[cell][4])
			pressure_list.append(cellDict[cell][5])
			contact_perimeter_list.append(cellDict[cell][6])
			num_het = len(cellHetNeighborTracking[cell])
			num_hom = len(cellHomNeighborTracking[cell])
			het_neighbour_count.append(num_het)
			hom_neighbour_count.append(num_hom)
			cell_neighbour_count = cellDict[cell][7]
			if cell_neighbour_count < 12:
				num_neighbours[cell_neighbour_count] += 1
			if num_het + num_hom != cell_neighbour_count:
				print "Warning Het count: ",num_het," Hom count: ",num_hom," Total count: ",cell_neighbour_count
		except IndexError:
			print "Warning Cell ID: ",cell," not in dictionary"
			continue	
	
	areaStat[mcs] = [NP.percentile(area_list,50), NP.percentile(area_list,90), NP.percentile(area_list,10)]
	perimeterStat[mcs] = [NP.percentile(perimeter_list,50), NP.percentile(perimeter_list,90), NP.percentile(perimeter_list,10)]
	pressureStat[mcs] = [NP.percentile(pressure_list,50), NP.percentile(pressure_list,90), NP.percentile(pressure_list,10)]
	neighborStat[mcs] = num_neighbours
	hetneighborStat[mcs] = sum(het_neighbour_count)/2
	homneighborStat[mcs] = sum(hom_neighbour_count)/2
	
del area_list
del perimeter_list
del contact_perimeter_list
del num_neighbours
del het_neighbour_count
del hom_neighbour_count

# Plot figures
numFigures = 0
timeseries = areaStat.keys()
timeseries.sort()

# Tissue Statistics Plot

numFigures += 1
PLT.figure()
PLT.subplot(311)
PLT.title('Area Statistics')
line1, = PLT.plot(timeseries, [areaStat[t][0] for t in timeseries], color='b', marker='.', linestyle='-', label='Median Area')
line2, = PLT.plot(timeseries, [areaStat[t][1] for t in timeseries], color='g', marker='.', linestyle='-', label='Q90 Area')
line3, = PLT.plot(timeseries, [areaStat[t][2] for t in timeseries], color='r', marker='.', linestyle='-', label='Q10 Area')
PLT.legend([line1, line2, line3], ['Q50', 'Q90', 'Q10'])

PLT.subplot(312)
PLT.title('Perimeter Statistics')
PLT.plot(timeseries, [perimeterStat[t][0] for t in timeseries], color='b', marker='.', linestyle='-', label='Median Perimeter')
PLT.plot(timeseries, [perimeterStat[t][1] for t in timeseries], color='g', marker='.', linestyle='-', label='Q90 Perimeter')
PLT.plot(timeseries, [perimeterStat[t][2] for t in timeseries], color='r', marker='.', linestyle='-', label='Q10 Perimeter')

PLT.subplot(313)
PLT.title('Cell Pressure Statistics')
PLT.plot(timeseries, [pressureStat[t][0] for t in timeseries], color='b', marker='.', linestyle='-', label='Median Pressure')
PLT.plot(timeseries, [pressureStat[t][1] for t in timeseries], color='g', marker='.', linestyle='-', label='Q90 Pressure')
PLT.plot(timeseries, [pressureStat[t][2] for t in timeseries], color='r', marker='.', linestyle='-', label='Q10 Pressure')
PLT.xlabel('Iteration #')
PLT.tight_layout()
PLT.savefig('CHASTE_plot1.png')
del areaStat
del perimeterStat
del pressureStat

# Neighbour Statistics Plot

numFigures += 1
PLT.figure()
PLT.subplot(111)
PLT.title('Neighbour Statistics')
PLT.ylabel('Number of Neighbours')
PLT.xlabel('Log(Time Step)')
ax = PLT.subplot(111)
for i in xrange(12):
	line, = ax.semilogx(itertodt(timeseries), [neighborStat[t][i] for t in timeseries], marker='.', linestyle='-', label='$P%i$'%i)

lgd = ax.legend(bbox_to_anchor=(0.0, 1.05, 1.0, 1.05), loc=3, ncol=4, mode="expand", borderaxespad=0.1, fancybox=True, shadow=True)
PLT.savefig('CHASTE_plot2.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
del neighborStat

# Heterotypic Boundary Measures

if num_cell_types >= 2:

	t_step_series = []
	heterotypicboundary = collections.defaultdict(list)
	
	with open(sys.argv[2],'r') as f:
		for line in f:
			values = re.split(r'\t+', line)
			t_step = float(values[0])
			t_step_series.append(t_step)
			heterotypicboundary[t_step] = float(values[2])		
	t_step_series.sort()
	
	numFigures += 1
	PLT.figure(numFigures)
	PLT.title('Heterotypic Boundary Between Cell Types: ' + cell_type_list[0] + ' and ' + cell_type_list[1])
	PLT.ylabel('Total Boundary Length')
	PLT.xlabel('Time Step (dt)')
	PLT.plot(t_step_series, [heterotypicboundary[t] for t in t_step_series], color='k', marker='.', linestyle='-')
	PLT.savefig('CHASTE_plot_Heterotypic.png')

	numFigures += 1
	PLT.figure(numFigures)
	PLT.title('Heterotypic Boundary Between Cell Types: ' + cell_type_list[0] + ' and ' + cell_type_list[1])
	PLT.ylabel('Total Boundary Length')
	PLT.xlabel('Log(Time Step)')
	PLT.semilogx(t_step_series, [heterotypicboundary[t] for t in t_step_series], color='k', marker='.', linestyle='-')
	PLT.savefig('CHASTE_plot_Heterotypic_semilogx.png')
	
	numFigures += 1
	PLT.figure(numFigures)
	PLT.title('Heterotypic Neighbour Count Between Cell Types: ' + cell_type_list[0] + ' and ' + cell_type_list[1])
	PLT.ylabel('Total Neighbour Count')
	PLT.xlabel('Time Step (dt)')
	PLT.plot(itertodt(timeseries), [hetneighborStat[t] for t in timeseries], color='k', marker='.', linestyle='-')
	PLT.savefig('CHASTE_plot_Heterotypic_neighbours.png')
	
	numFigures += 1
	PLT.figure(numFigures)
	PLT.title('Homotypic Neighbour Count Between Cell Types: ' + cell_type_list[0] + ' and ' + cell_type_list[1])
	PLT.ylabel('Total Neighbour Count')
	PLT.xlabel('Time Step (dt)')
	PLT.plot(itertodt(timeseries), [homneighborStat[t] for t in timeseries], color='k', marker='.', linestyle='-')
	PLT.savefig('CHASTE_plot_Homotypic_neighbours.png')
	
del heterotypicboundary

# Swap Statistics

numFigures += 1
PLT.figure(numFigures)
t_step_series = []
numswaps = collections.defaultdict(list)
	
with open(sys.argv[3],'r') as f:
	for line in f:
		values = re.split(r'\t+', line)
		t_step = float(values[0])
		t_step_series.append(t_step)
		numswaps[t_step].append(float(values[1]))
				
t_step_series.sort()

with open(sys.argv[4],'r') as f:
	for line in f:
		values = re.split(r'\t+', line)
		t_step = float(values[0])
		numswaps[t_step].append(float(values[1]))
		
with open(sys.argv[5],'r') as f:
	for line in f:
		values = re.split(r'\t+', line)
		t_step = float(values[0])
		numswaps[t_step].append(float(values[1]))	

PLT.subplot(311)
PLT.title('T1 Swap Statistics')
PLT.plot(t_step_series, [numswaps[t][0] for t in t_step_series], color='b', marker='.', linestyle='.')
[xmin, xmax, ymin, ymax] = PLT.axis()
PLT.axis([0,xmax,ymin-0.5,ymax+0.5])

PLT.subplot(312)
PLT.title('T2 Swap Statistics')
PLT.plot(t_step_series, [numswaps[t][1] for t in t_step_series], color='g', marker='.', linestyle='.')
[xmin, xmax, ymin, ymax] = PLT.axis()
PLT.axis([0,xmax,ymin-0.5,ymax+0.5])

PLT.subplot(313)
PLT.title('T3 Swap Statistics')
PLT.plot(t_step_series, [numswaps[t][2] for t in t_step_series], color='r', marker='.', linestyle='.')
[xmin, xmax, ymin, ymax] = PLT.axis()
PLT.axis([0,xmax,ymin-0.5,ymax+0.5])

PLT.xlabel('Time Step (dt)')
PLT.tight_layout()
PLT.savefig('CHASTE_plot3.png')

# Cell Type Statistics

numFigures += 1
PLT.figure(numFigures)
PLT.subplot(211)
PLT.title('Cell Count Statistics')
PLT.ylabel('Number of Cells')
PLT.xlabel('Time Step (dt)')
line1, = PLT.plot(itertodt(timeseries), [cellNumStat[t][0] for t in timeseries], color='g', marker='.', label=cell_type_list[0])
line2, = PLT.plot(itertodt(timeseries), [cellNumStat[t][1] for t in timeseries], color='r', marker='.', label=cell_type_list[1])
line3, = PLT.plot(itertodt(timeseries), [sum(cellNumStat[t]) for t in timeseries], color='k', marker='.', label='Total')
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
PLT.savefig('CHASTE_plot4.png')
