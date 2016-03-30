#!/usr/bin/env python

#
# Last modified: 20 Aug 2015
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# 

import os
import re
import sys
import math
import signal
import datetime
import collections
import numpy as NP
import lxml.etree as ET
import matplotlib.pyplot as PLT
from optparse import OptionParser
from matplotlib.widgets import Slider, Button, RadioButtons

def signal_handler(signal, frame):
	print('\nCtrl+C Interrupt: Terminating')
	sys.exit()
signal.signal(signal.SIGINT, signal_handler)

xmlFile = ''
drawscale = 0.1
animate = False
toggle = 1

parser = OptionParser()

parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to XML FILE", metavar="FILE")
parser.add_option("-d", "--draw", action="store", type="float", dest="scale", help="drawing SCALE (between 0.1 and 1)", metavar="SCALE")
parser.add_option("-a", "--animate", action="store_true", dest="animate", default=False, help="output images for animation")
parser.add_option("-f", "--force", action="store_true", dest="plotForce", default=True, help="plot force vector field")
parser.add_option("-v", "--velocity", action="store_false", dest="plotForce", default=True, help="plot velocity vector field")

(options, args) = parser.parse_args()

if options.inputfile:
	xmlFile = options.inputfile
if options.plotForce:
	toggle = 0
	drawscale = 1
if options.scale:
	drawscale = float(options.scale)
if options.animate:
	animate = options.animate


if (os.path.isfile(xmlFile) == False):
	print '\nPlease provide a valid XML file. See help: python widget_sim_data.py --help'
	sys.exit()
	
sim = ET.parse(open(xmlFile,'r'))

root = sim.getroot()
timeseries = []
cellTracking = collections.defaultdict(list)

for time in sim.getiterator('time'):
	
	mcs = int(time.get('t'))
	timeseries.append(mcs)
	
	for cell in time.getchildren():
	
		p_vec = [float(cell.get('x')), float(cell.get('y'))]
		f_vec = re.split('\s+',cell.get('ext_force'))
		
		case = int(cell.get('type'))
		if case == 1:
			cellTracking[cell.get('cell_id')].append(['green',mcs,float(p_vec[0]),float(p_vec[1]),float(f_vec[0]),float(f_vec[1])])
		elif case == 2:
			cellTracking[cell.get('cell_id')].append(['red',mcs,float(p_vec[0]),float(p_vec[1]),float(f_vec[0]),float(f_vec[1])])
		elif case == 3:
			cellTracking[cell.get('cell_id')].append(['blue',mcs,float(p_vec[0]),float(p_vec[1]),float(f_vec[0]),float(f_vec[1])])
		else:
			cellTracking[cell.get('cell_id')].append(['black',mcs,float(p_vec[0]),float(p_vec[1]),float(f_vec[0]),float(f_vec[1])])

		
axcolor = 'lightgoldenrodyellow'

fig, ax = PLT.subplots()
PLT.subplots_adjust(bottom=0.25)
mcs = NP.median(NP.array(timeseries))

# Plot Force/Velocity Vector Field

def plot_figure(tg, ts, ax):
	p_idata = []
	p_fdata = []
	fdata = []
	tol = 50
	xmax, ymax = 0, 0
	for cell_id in cellTracking:
		xl, yl = -1, -1
		for com in cellTracking[cell_id]:
			if (com[1] == ts):
				xl, yl = com[2], com[3]
				p_idata.append([com[0], com[2], com[3]])
				fdata.append([com[0], -1*com[4], -1*com[5]])
			elif (com[1] == ts + 1 and xl != -1 and yl != -1):
				# check if the cell is newborn and check for periodic boundary conditions
				xh, yh = -1, -1
				if com[2] < xl - tol:
					xh = xl + 5
				elif com[2] > xl + tol:
					xh = 0
				if com[3] < yl - tol:
					yh = yl + 5
				elif com[3] > yl + tol:
					yh = 0
				if xh == -1:
					xh = com[2]
				if yh == -1:
					yh = com[3]
				p_fdata.append([com[0], xh, yh])
			xmax = max(xmax, com[2])
			ymax = max(ymax, com[3])
	X = []
	Y = []
	U = []
	V = []
	C = []
			
	for i in range(len(p_idata)):
	
		C.append(p_idata[i][0])
		X.append(p_idata[i][1])
		Y.append(p_idata[i][2])
		
		if tg == 0:
			U.append(fdata[i][1])
			V.append(fdata[i][2])
		elif tg == 1:
			U.append(p_fdata[i][1] - p_idata[i][1])
			V.append(p_fdata[i][2] - p_idata[i][2])
	
	ax.quiver(X, Y, U, V, color=C, angles='xy', scale_units='xy', scale=drawscale)
	ax.set_xlim([0,xmax+5])
	ax.set_ylim([0,ymax+5])

ax = PLT.gca()

# Generate PNG files to animate

if (animate == True):
	totalMCS = len(timeseries)
	totalDigits = len(str(totalMCS-1))
	path = os.getcwd() + os.path.sep
	if toggle == 0:
		path += "Force_cc3d_"+datetime.datetime.now().strftime('%m_%d_%y_%H_%M_%S')
	elif toggle == 1:
		path += "Velocity_cc3d_"+datetime.datetime.now().strftime('%m_%d_%y_%H_%M_%S')
		totalMCS = max(timeseries) - 1
	os.makedirs(path)
	print "\nWriting Files to Directory: "+path
	for j in range(totalMCS):
		plot_figure(toggle,j,ax)
		filename = ""
		if toggle == 0:
			filename += "Force"
		elif toggle == 1:
			filename += "Velocity"
		filename += "_cc3d_"+str(j).zfill(totalDigits)+".png"
		percent = "%d" % ((j*100.0)/totalMCS)
		print "Generating PNG Files ( "+percent+"% Complete )             \r",
		filename = os.path.join(path, filename)
		PLT.savefig(filename, bbox_inches='tight')
		ax.cla()
	print ""
	sys.exit()

plot_figure(toggle,mcs,ax)

# Interactive Mode - Widget Controls 

axMCS = PLT.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
sliderMCS = Slider(axMCS, 'MCS', 1, max(timeseries) - 1, valinit=mcs, valfmt='%0.0f')

def update(val):
    mcs = sliderMCS.val
    ax.cla()
    plot_figure(toggle,int(round(mcs)),ax)
    if toggle == 0:
    	PLT.title('Force Vector Field at MCS = ' + str(int(round(mcs))))
    elif toggle == 1:
    	PLT.title('Velocity Vector Field at MCS = ' + str(int(round(mcs))))
    fig.subplots_adjust()
    PLT.grid(True)
    PLT.draw()

sliderMCS.on_changed(update)

update(mcs)
PLT.grid(True)
PLT.show()
