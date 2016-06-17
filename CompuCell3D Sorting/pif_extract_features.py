#!/usr/bin/env python

#
# Last modified: 19 May 2016
# Author: Darrick Lee <y.l.darrick@gmail.com>, Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Extract cell shape features from PIF file
# 

import os
import sys
import math
import shlex
import string
import threading
import subprocess
import extractcellfeatures

import numpy as NP
import lxml.etree as ET
import matplotlib.pyplot as PLT

from PIL import Image
from scipy import interpolate
from optparse import OptionParser
from scipy import ndimage as NDI
from scipy.special import expit
from matplotlib.patches import Ellipse, Polygon

imgDPI = 200
radius_thresh = 50 # Threshold value for circle radius (also used for ellipse)

pifFile = None
xmlFile = None
outFolder = './'
l_height, l_width = -1, -1
xmlTime = -1
testCell = -1
plotZoom = -1
plotfits, comparefits, separate = None, None, None
splineSmooth = 10

parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-x", "--xml", action="store", type="string", dest="xmlfile", help="path to XML file", metavar="XML")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")
parser.add_option("-p","--plot", action="store_true", dest="plotfits", help="plot fitted geometry", default=False)
parser.add_option("-c","--compare", action="store_true", dest="comparefits", help="compare fitted geometry", default=False)
parser.add_option("-t","--time", action="store", type="int", dest="time", help="time in XML file to compare", metavar="TIME")
parser.add_option("-T","--test", action="store", type="int", dest="testcell", help="runs the test code on specified cell")
parser.add_option("-o","--output", action="store", type="string", dest="outputfolder", help="the folder to store output plots", metavar="OUTPUT")
parser.add_option("-s","--separate", action="store_true", dest="separate", help="separate subfolders for plots", default=False, metavar="SEPARATE")
parser.add_option("-m","--multithread", action="store_true", dest="multithread", help="use multithreading", default=False)
parser.add_option("-z","--zoom", action="store", type="int", dest="zoom", help="units to zoom in each direction", metavar="ZOOM")
parser.add_option("-S","--smooth", action="store", type="float", dest="smooth", help="amount of smoothing for spline", metavar="SMOOTH")
parser.add_option("-C","--curvSpline", action="store_true", dest="curvSpline", default=False, help="create plot for curvature-dependent boundary for spline", metavar="CURVSPLINE")

# Options parsing
(options, args) = parser.parse_args()
if options.inputfile:
	pifFile = options.inputfile
if options.xmlfile:
	xmlFile = options.xmlfile
if options.height:
	l_height = options.height
if options.width:
	l_width = options.width
if options.plotfits:
	plotfits = True
else:
	plotfits = False
if options.comparefits:
	comparefits = True
else:
	comparefits = False
if options.time is not None:
	xmlTime = options.time
if options.testcell:
	testCell = options.testcell
if options.outputfolder:
	outFolder = options.outputfolder
if options.separate:
	separate = True
else:
	separate = False
if options.multithread:
	multithreading = True
else:
	multithreading = False
if options.zoom:
	plotZoom = options.zoom
if options.smooth:
	splineSmooth = options.smooth
if options.curvSpline:
	curvSpline = True
else:
	curvSpline = False

plotXML = xmlFile is not None and xmlTime != -1

pifFileName = ''

if os.path.isfile(pifFile):
	pifFileName_withPath = os.path.splitext(pifFile)[0]
	pifFileName = pifFileName_withPath.split('/')[-1]
else:
	print("Error: PIF file does not exist.\n")
	sys.exit(1)

if not os.path.isdir(outFolder):
	print("Error: Output directory does not exist.\n")
	sys.exit(1)

# Check if separate folders have been made and initialize output folders
if separate:
	circleOut = outFolder + 'CircleFit/'
	ellipseOut = outFolder + 'EllipseFit/'
	splineOut = outFolder + 'SplineFit/'
	splineBdyOut = outFolder + 'SplineBdyFit/'
	polyOut = outFolder + 'PolyFit/'
	boundaryOut = outFolder + 'BoundaryFit/'
	vectorizeOut = outFolder + 'vectorizeFit/'

	if (not os.path.isdir(circleOut) or not os.path.isdir(ellipseOut)
		or not os.path.isdir(splineOut) or not os.path.isdir(polyOut)
		or not os.path.isdir(boundaryOut) or not os.path.isdir(vectorizeOut)
		or not os.path.isdir(splineBdyOut)):
		print("Error: One of the output folders do not exist.\n")
		sys.exit(1)
		
else:
	circleOut = ellipseOut = splineOut = polyOut = boundaryOut = vectorizeOut = outFolder

# Parse PIF file	
lattice = open(pifFile).read().split('\n')[1:-1]
cellDict = dict()
cellTypeDict = dict()
lattice_data = NP.zeros((l_width, l_height))
lattice_matrix = NP.empty((l_width, l_height), NP.uint32)
lattice_matrix.fill(0xFFFFFFFF)

for pix in lattice:

	fields = pix.split()

	cell_id = int(fields[0])
	c_type = str(fields[1])
	pix_x = int(fields[3])
	pix_y = int(fields[5])
            
	if cell_id not in cellDict.keys():
    
		cellDict[cell_id] = [[pix_x, pix_y]]

	cellDict[cell_id].append([pix_x, pix_y])
	cellTypeDict[cell_id] = c_type
	lattice_data[pix_x, pix_y] = cell_id
	
	if c_type == 'CellU':
		lattice_matrix[pix_x, pix_y] = 0x8000EE00	# alpha, blue, green, red
	if c_type == 'CellV':
		lattice_matrix[pix_x, pix_y] = 0x800000EE


# Check if lattice contains isolated cells

def contains_isolated_cells():

	'''
	Description: This returns true if lattice_data contains more than one connected component
	and false otherwise. This currently uses 1-connectivity.

	Reference: http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label
	'''

	global lattice_data
	
	clipped_lattice_data = NP.clip(lattice_data,0,1)

	[labeled_img, num_labels] = NDI.measurements.label(clipped_lattice_data)

	if num_labels > 1:
		return True
	else:
		return False

# Conversion functions
def conv_distance(x):
	conv_factor = 0.8
	units = 'um'

	return x*conv_factor, units

def conv_time(x):
	conv_factor = 0.2
	units = 'min'

	return x*conv_factor, units

def conv_area(x):
	conv_factor = 0.64
	units = 'um^2'

	return x*conv_factor, units

# Plot fits for a single cell (for testing and illustration purposes)

def TestSingleCellPlot(extractor):

	extractor.basic_props(splineSmooth)
	extractor.ellipse_props()
	extractor.cell_centre_fit()

	fixed_perim = NP.transpose(extractor.perim_img)

	perim_img_ind = NP.where(fixed_perim == 1)

	xlim_min = min(perim_img_ind[1])-2
	xlim_max = max(perim_img_ind[1])+2

	ylim_min = min(perim_img_ind[0])-2
	ylim_max = max(perim_img_ind[0])+2

	U = extractor.spl_u
	OUT = interpolate.splev(U, extractor.spl_poly)

	# Create Circle/Ellipse plot
	fig = PLT.figure(1)
	ax = fig.add_subplot(111, aspect='equal')

	# Create cirlce plot
	c = PLT.Circle((extractor.ccm_fvector[0],
			extractor.ccm_fvector[1]),
			extractor.ccm_fvector[2])

	# Create ellipse plot
	e = Ellipse(xy=NP.array([extractor.ellipse_fvector[0], extractor.ellipse_fvector[1]]),
			width = extractor.ellipse_fvector[6],
			height = extractor.ellipse_fvector[5],
			angle = extractor.ellipse_fvector[7]/(2*NP.pi)*360)

	PLT.imshow(fixed_perim, interpolation='nearest', cmap='Greys')
	PLT.plot(extractor.perim_coord_poly[:,0], extractor.perim_coord_poly[:,1], label='Poly', linewidth=3, color='g')
	
	ax.add_artist(c)
	c.set_alpha(1)
	c.set_facecolor('none')
	c.set_edgecolor('blue')
	c.set_linewidth(3)
	c.set_label('Circle')

	ax.add_artist(e)
	e.set_alpha(1)
	e.set_facecolor('none')
	e.set_edgecolor('orange')
	e.set_linewidth(3)
	e.set_label('Ellipse')

	PLT.plot(0,0,color='blue',label='Circle',lw=3)
	PLT.plot(0,0,color='orange',label='Ellipse',lw=3)

	lgd = PLT.legend(bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)

	ax.set_xlim([xlim_min, xlim_max])
	ax.set_ylim([ylim_min, ylim_max])
	PLT.savefig(pifFileName + '_Fits.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 400)

	# Create spline plot with boundary color based on magnitude and parity of curvature
	fig = PLT.figure(2)
	ax = fig.add_subplot(111, aspect='equal')
	knorm = expit(extractor.spl_k/max(abs(extractor.spl_k))*10)
	pcolor = PLT.cm.bwr(knorm)

	PLT.imshow(fixed_perim, interpolation='nearest', cmap='Greys')
	
	for i in range(len(U)):
		PLT.plot(OUT[0][i:i+2],OUT[1][i:i+2],color=pcolor[i],linewidth=3)

	ax.set_xlim([xlim_min, xlim_max])
	ax.set_ylim([ylim_min, ylim_max])

	PLT.savefig(pifFileName + '_SplineCurvature.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 400)

	# Create spline plot with a binary boundary color scheme based on parity of curvature
	fig = PLT.figure(3)
	ax = fig.add_subplot(111, aspect='equal')
	knorm = NP.sign(extractor.spl_k)/2 + 0.5
	pcolor = PLT.cm.bwr(knorm)

	PLT.imshow(fixed_perim, interpolation='nearest', cmap='Greys')
	
	for i in range(len(U)):
		PLT.plot(OUT[0][i:i+2],OUT[1][i:i+2],color=pcolor[i],linewidth=3)

	ax.set_xlim([xlim_min, xlim_max])
	ax.set_ylim([ylim_min, ylim_max])

	PLT.plot(0,0,color='blue',label='Negative Curvature',lw=3)
	PLT.plot(0,0,color='red',label='Positive Curvature',lw=3)

	lgd = PLT.legend(bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=2, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)

	PLT.savefig(pifFileName + '_SplineCurvatureBin.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 400)

		
# Compute features for all cells

featureDict = dict()
polyPtDict = dict()
splinePtDict = dict()
curvatureDict = dict()

for cell_id in cellDict.keys():
	
	if testCell != -1:
		if cell_id == testCell:
			extractor = extractcellfeatures.ExtractFeatures(cellDict[cell_id])
			TestSingleCellPlot(extractor)
			sys.exit(0)
	else:

		extractor = extractcellfeatures.ExtractFeatures(cellDict[cell_id])

		if multithreading:
			thread_list = []

			thread_list.append(threading.Thread(target=extractor.basic_props(splineSmooth), args=(), kwargs={}))
			thread_list.append(threading.Thread(target=extractor.ellipse_props(), args=(), kwargs={}))
			thread_list.append(threading.Thread(target=extractor.cell_centre_fit(), args=(), kwargs={}))
		
			for thread in thread_list:
				thread.start()
		
			for thread in thread_list:
				thread.join()
		else:
			extractor.basic_props(splineSmooth)
			extractor.ellipse_props()
			extractor.cell_centre_fit()
	
		featureDict[cell_id] = [extractor.area_cell, extractor.perim_1sqrt2, extractor.equiv_diameter]
		featureDict[cell_id] = featureDict[cell_id] +  [extractor.perim_3pv]
		featureDict[cell_id] = featureDict[cell_id] + [extractor.perim_poly, extractor.area_poly]
		featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_fvector + extractor.ccm_fvector
		featureDict[cell_id] = featureDict[cell_id] + [extractor.perim_eroded,  extractor.area_eroded]
		featureDict[cell_id] = featureDict[cell_id] + extractor.bdy_fvector

		polyPtDict[cell_id] = extractor.perim_coord_poly

		# Put spline points and curvature data to be used in plots
		splinePtDict[cell_id] = NP.transpose(interpolate.splev(extractor.spl_u, extractor.spl_poly))
		curvatureDict[cell_id] = extractor.spl_k

# Construct featIndexDict

featIndexDict = dict(BASIC=None, ELLIPSE=None, CCM=None, TPV=None, ERODED=None, POLY=None, BDY=None)
BASIC_numfeat = 3
TPV_numfeat = 1
POLY_numfeat = 2
ELLIPSE_numfeat = 11
CCM_numfeat = 5
ERODED_numfeat = 2
BDY_numfeat = 6

TPV_start = BASIC_numfeat
POLY_start = TPV_start + TPV_numfeat
ELLIPSE_start = POLY_start + POLY_numfeat
CCM_start = ELLIPSE_start + ELLIPSE_numfeat
ERODED_start = CCM_start + CCM_numfeat
BDY_start = ERODED_start + ERODED_numfeat

featIndexDict['BASIC'] = dict(
	area = 0,
	perimeter =1,
	equiv_diameter=2)
	
featIndexDict['TPV'] = dict(
	perimeter=TPV_start)
	
featIndexDict['POLY'] = dict(
	perimeter=POLY_start,
	area=POLY_start+1)
	
featIndexDict['ELLIPSE'] = dict(
	centroid_x=ELLIPSE_start, 
	centroid_y=ELLIPSE_start+1,
	eccentricity=ELLIPSE_start+2,
	euler_number=ELLIPSE_start+3,
	extent=ELLIPSE_start+4,
	major_axis_length=ELLIPSE_start+5,
	minor_axis_length=ELLIPSE_start+6,
	orientation=ELLIPSE_start+7,
	solidity=ELLIPSE_start+8,
	area=ELLIPSE_start+9,
	perimeter=ELLIPSE_start+10)
	
featIndexDict['CCM'] = dict(
	centroid_x=CCM_start,
	centroid_y=CCM_start+1,
	radius=CCM_start+2,
	perimeter=CCM_start+3,
	area=CCM_start+4)
	
featIndexDict['ERODED'] = dict(
	perimeter=ERODED_start,
	area=ERODED_start+1)

featIndexDict['BDY'] = dict(
	maxpos=BDY_start,
	minpos=BDY_start+1,
	maxneg=BDY_start+2,
	minneg=BDY_start+3,
	num_extrema=BDY_start+4,
	num_signflip=BDY_start+5)

# Plot ellipsoid fit, cell-centre model fit, fitted polygons
numFigs = 0;

if plotfits:

	# Plot original cells
	lattice_matrix = NP.ascontiguousarray(NP.flipud(NP.transpose(lattice_matrix)))
	pilImage = Image.frombuffer('RGBA', (l_width, l_height), lattice_matrix, 'raw', 'RGBA', 0, 1)
	pilImage = pilImage.convert('RGB')

	pilImage.save(boundaryOut + pifFileName + '_Boundary.png')
	
	# # Plot polygonized lattice
	# sp = None
	# if not contains_isolated_cells():
	# 	dirname = os.path.dirname(os.path.abspath(__file__))
	# 	if os.path.isfile(os.path.join(dirname, 'vectorize.py')):
	# 		cmd = "python3 vectorize.py --file " + pifFile + " --size " + str(l_width) + "," + str(l_height) + " --output " + vectorizeOut
	# 		args = shlex.split(cmd)
	# 		sp = subprocess.Popen(args)

	# Plot ellipse fits
	numFigs += 1
	fig = PLT.figure(numFigs)
	ax = fig.add_subplot(111, aspect='equal')

	ells = [[Ellipse(xy=NP.array([featureDict[cell_id][featIndexDict['ELLIPSE']['centroid_x']]
		,featureDict[cell_id][featIndexDict['ELLIPSE']['centroid_y']]]),
		width = featureDict[cell_id][featIndexDict['ELLIPSE']['minor_axis_length']],
		height = featureDict[cell_id][featIndexDict['ELLIPSE']['major_axis_length']],
		angle = featureDict[cell_id][featIndexDict['ELLIPSE']['orientation']]/(2*NP.pi)*360),
		cellTypeDict[cell_id]] for cell_id in cellDict.keys()]
	
	for el in ells:
		e = el[0]
		ctype = el[1]

		# Don't plot the ellipse if it is too big
		if (e.height > radius_thresh) or (e.width > radius_thresh):
			continue

		ax.add_artist(e)
		e.set_clip_box(ax.bbox)
		e.set_alpha(0.3)
		if ctype == 'CellU':
			e.set_facecolor([0,1,0])
		elif ctype == 'CellV':
			e.set_facecolor([1,0,0])
		else:
			e.set_facecolor([0,0,1])
	
	if plotZoom != -1:
		ax.set_xlim([plotZoom,l_width-plotZoom])
		ax.set_ylim([plotZoom,l_height-plotZoom])
	else:
		ax.set_xlim([0,l_width])
		ax.set_ylim([0,l_height])
	PLT.savefig(ellipseOut + pifFileName + '_EllipseFit.png', bbox_inches='tight', dpi = imgDPI)
	
	# Plot cell-centre model (disk fit)
	numFigs += 1
	fig = PLT.figure(numFigs)
	ax = fig.add_subplot(111, aspect='equal')

	circles = [[PLT.Circle((featureDict[cell_id][featIndexDict['CCM']['centroid_x']],
		featureDict[cell_id][featIndexDict['CCM']['centroid_y']]),
		featureDict[cell_id][featIndexDict['CCM']['radius']]),
		cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

	for circle in circles:
		c = circle[0]
		ctype = circle[1]

		# Don't plot the circle if it is too big
		if c.radius > radius_thresh:
			continue

		ax.add_artist(c)
		c.set_alpha(0.3)
		if ctype == 'CellU':
			c.set_facecolor([0,1,0])
		elif ctype == 'CellV':
			c.set_facecolor([1,0,0])
		else:
			c.set_facecolor([0,0,1])

	if plotZoom != -1:
		ax.set_xlim([plotZoom,l_width-plotZoom])
		ax.set_ylim([plotZoom,l_height-plotZoom])
	else:
		ax.set_xlim([0,l_width])
		ax.set_ylim([0,l_height])
	PLT.savefig(circleOut + pifFileName + '_CircleFit.png', bbox_inches='tight', dpi = imgDPI)

	# Plot 3pv-polygon model
	numFigs += 1
	fig = PLT.figure(numFigs)
	ax = fig.add_subplot(111, aspect='equal')

	polys = [[Polygon(NP.array(polyPtDict[cell_id])),
		cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

	for poly in polys:
		p = poly[0]
		ctype = poly[1]
		ax.add_artist(p)
		p.set_alpha(0.3)
		if ctype == 'CellU':
			p.set_facecolor([0,1,0])
		elif ctype == 'CellV':
			p.set_facecolor([1,0,0])
		else:
			p.set_facecolor([0,0,1])

	if plotZoom != -1:
		ax.set_xlim([plotZoom,l_width-plotZoom])
		ax.set_ylim([plotZoom,l_height-plotZoom])
	else:
		ax.set_xlim([0,l_width])
		ax.set_ylim([0,l_height])
	PLT.savefig(polyOut + pifFileName + '_PolyFit.png', bbox_inches='tight', dpi = imgDPI)
	
	# Plot spline model
	numFigs += 1
	fig = PLT.figure(numFigs)
	ax = fig.add_subplot(111, aspect='equal')

	polys = [[Polygon(NP.array(splinePtDict[cell_id])),
		cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

	for poly in polys:
		p = poly[0]
		ctype = poly[1]
		ax.add_artist(p)
		p.set_alpha(0.3)
		if ctype == 'CellU':
			p.set_facecolor([0,1,0])
		elif ctype == 'CellV':
			p.set_facecolor([1,0,0])
		else:
			p.set_facecolor([0,0,1])

	if plotZoom != -1:
		ax.set_xlim([plotZoom,l_width-plotZoom])
		ax.set_ylim([plotZoom,l_height-plotZoom])
	else:
		ax.set_xlim([0,l_width])
		ax.set_ylim([0,l_height])
	PLT.savefig(splineOut + pifFileName + '_SplineFit.png', bbox_inches='tight', dpi = imgDPI)

	# Plot spline model with curvature dependent boundary
	if curvSpline:
		numFigs += 1
		fig = PLT.figure(numFigs)
		ax = fig.add_subplot(111, aspect='equal')

		polys = [[Polygon(NP.array(splinePtDict[cell_id])),
			cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

		for poly in polys:
			p = poly[0]
			ctype = poly[1]
			ax.add_artist(p)
			p.set_alpha(0.3)
			p.set_edgecolor('none')
			if ctype == 'CellU':
				p.set_facecolor([0,1,0])
			elif ctype == 'CellV':
				p.set_facecolor([1,0,0])
			else:
				p.set_facecolor([0,0,1])

		for cell_id in cellDict.keys():
			k = curvatureDict[cell_id]
			OUT = NP.transpose(splinePtDict[cell_id])

			knorm = NP.sign(k)/2 + 0.5
			pcolor = PLT.cm.bwr(knorm)

			for i in range(len(k)):
				PLT.plot(OUT[0][i:i+2],OUT[1][i:i+2],color=pcolor[i],linewidth=0.75)

		if plotZoom != -1:
			ax.set_xlim([plotZoom,l_width-plotZoom])
			ax.set_ylim([plotZoom,l_height-plotZoom])
		else:
			ax.set_xlim([0,l_width])
			ax.set_ylim([0,l_height])
		PLT.savefig(splineBdyOut + pifFileName + '_SplineBdyFit.png', bbox_inches='tight', dpi = imgDPI)
	
	# if sp is not None:
	# 	sp.wait()


# Compare features
 
if comparefits:

	# Initialize data vectors
	xData = []

	perim_1sqrt2 = []
	perim_ellipse = []
	perim_circle = []
	perim_eroded = []
	perim_3pv = []
	perim_poly = []
	perim_xml = []

	area_ellipse = []
	area_circle = []
	area_basic = []
	area_eroded = []
	area_poly = []
	area_xml = []

	# Prepare XML file
	if plotXML:
		infile = open(xmlFile,'r')
		xml1 = ET.parse(infile)
		root = xml1.getroot()
		times = root.getchildren()
		cells = times[xmlTime].getchildren()

	for cell_id in cellDict.keys():
		xData.append(cell_id)

	# Populate the data vectors
	for cell_id in xData:
		perim_1sqrt2.append(featureDict[cell_id][featIndexDict['BASIC']['perimeter']])
		perim_ellipse.append(featureDict[cell_id][featIndexDict['ELLIPSE']['perimeter']])
		perim_circle.append(featureDict[cell_id][featIndexDict['CCM']['perimeter']])
		perim_eroded.append(featureDict[cell_id][featIndexDict['ERODED']['perimeter']])
		perim_3pv.append(featureDict[cell_id][featIndexDict['TPV']['perimeter']])
		perim_poly.append(featureDict[cell_id][featIndexDict['POLY']['perimeter']])

		area_ellipse.append(featureDict[cell_id][featIndexDict['ELLIPSE']['area']])
		area_circle.append(featureDict[cell_id][featIndexDict['CCM']['area']])
		area_eroded.append(featureDict[cell_id][featIndexDict['ERODED']['area']])
		area_basic.append(featureDict[cell_id][featIndexDict['BASIC']['area']])
		area_poly.append(featureDict[cell_id][featIndexDict['POLY']['area']])

		if plotXML:
			perim_xml.append(cells[cell_id].get('perimeter'))
			area_xml.append(cells[cell_id].get('area'))

	# Reorder the vectors so that xml parameters are sorted from smallest to largest
	# Note: We switch to NP.array rather than list so we can input a list to select elements
	perim_reorder_ind = NP.argsort(perim_1sqrt2)
	perim_1sqrt2, perim_unit = conv_distance(NP.array(perim_1sqrt2)[perim_reorder_ind])
	perim_ellispe, _ = conv_distance(NP.array(perim_ellipse)[perim_reorder_ind])
	perim_circle, _ = conv_distance(NP.array(perim_circle)[perim_reorder_ind])
	perim_poly, _ = conv_distance(NP.array(perim_poly)[perim_reorder_ind])
	perim_eroded, _ = conv_distance(NP.array(perim_eroded)[perim_reorder_ind])
	perim_3pv, _ = conv_distance(NP.array(perim_3pv)[perim_reorder_ind])

	area_reorder_ind = NP.argsort(area_basic)
	area_ellipse, area_unit = conv_area(NP.array(area_ellipse)[area_reorder_ind])
	area_circle, _ = conv_area(NP.array(area_circle)[area_reorder_ind])
	area_poly, _ = conv_area(NP.array(area_poly)[area_reorder_ind])
	area_eroded, _ = conv_area(NP.array(area_eroded)[area_reorder_ind])
	area_basic, _ = conv_area(NP.array(area_circle)[area_reorder_ind])

	if plotXML:
		perim_xml = NP.array(perim_xml)[perim_reorder_ind]
		area_xml = NP.array(area_xml)[area_reorder_ind]

	cell_range = range(len(xData))

	# Plot the figures
	numFigs += 1
	PLT.figure(numFigs)
	if plotXML:
		PLT.plot(cell_range, perim_xml)
	PLT.plot(cell_range, perim_ellipse)
	PLT.plot(cell_range, perim_circle)
	PLT.plot(cell_range, perim_poly)
	PLT.plot(cell_range, perim_eroded)
	PLT.plot(cell_range, perim_1sqrt2)
	PLT.plot(cell_range, perim_3pv)

	PLT.xlabel('Cell (arbitrary)')
	PLT.ylabel('Perimeter (' + perim_unit + ')')
	PLT.title('Comparison of Perimeter vs. Cell')
	lgd = PLT.legend(["Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "1, sqrt(2) method", "3pv method"], 
	bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)
	
	if plotXML:
		lgd = PLT.legend(["xml", "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "1, sqrt(2) method", "3pv method"],
		bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)
	
	PLT.savefig(outFolder + pifFileName + '_PerimCompare.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = imgDPI)

	numFigs += 1
	PLT.figure(numFigs)
	if plotXML:
		PLT.plot(cell_range, area_xml)
	PLT.plot(cell_range, area_ellipse)
	PLT.plot(cell_range, area_circle)
	PLT.plot(cell_range, area_poly)
	PLT.plot(cell_range, area_eroded)
	PLT.plot(cell_range, area_basic)

	PLT.xlabel('Cell (arbitrary)')
	PLT.ylabel('Area (' + area_unit + ')')
	PLT.title('Comparison of Area vs. Cell')
	lgd = PLT.legend([ "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "Pixel counting"],
	bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)
	
	if plotXML:
		PLT.legend([ "xml", "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "Pixel counting"],
		bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)
				
	PLT.savefig(outFolder + pifFileName + '_AreaCompare.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = imgDPI)
