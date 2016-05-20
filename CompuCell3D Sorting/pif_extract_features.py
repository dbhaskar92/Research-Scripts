#!/usr/bin/env python

#
# Last modified: 11 May 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Extract cell shape features from PIF file
# 

import os
import sys
import math
import shlex
import string
import threading
import subprocess
import scipy.special
import skimage.measure
import skimage.morphology
import numpy as NP
import lxml.etree as ET
import matplotlib.pyplot as PLT
from PIL import Image
from operator import itemgetter
from optparse import OptionParser
from matplotlib.patches import Ellipse, Polygon
from scipy import ndimage as NDI
from scipy import interpolate
from perimeter_3pvm import perimeter_3pvm

###################################################################

def testFunction(extractor):
	# This function is used to test new features

	extractor.basic_props()

	perim_img_ind = NP.where(extractor.perim_img == 1)

	xlim_min = min(perim_img_ind[1])-2
	xlim_max = max(perim_img_ind[1])+2

	ylim_min = min(perim_img_ind[0])-2
	ylim_max = max(perim_img_ind[0])+2

	fig = PLT.figure(1)
	ax = fig.add_subplot(111, aspect='equal')

	U = NP.linspace(0,1.01,100)
	OUT = interpolate.splev(U, extractor.spl_poly)
	OUT = NP.flipud(OUT)

	PLT.imshow(extractor.perim_img, interpolation='nearest', cmap='Greys')
	PLT.plot(extractor.perim_coord_poly[:,1], extractor.perim_coord_poly[:,0], linewidth=3, color='g')
	PLT.plot(OUT[0], OUT[1], linewidth=3, color='r')

	# PLT.plot(extractor.perim_coord_eroded[:,1], extractor.perim_coord_eroded[:,0], linewidth=3, color='r')
	# PLT.plot(extractor.perim_coord_dp[:,1], extractor.perim_coord_dp[:,0], linewidth=3, color='m')

	PLT.legend(['Poly','Spline','D-P'], loc=2)

	ax.set_xlim([xlim_min, xlim_max])
	ax.set_ylim([ylim_min, ylim_max])
	PLT.savefig(pifFileName + '_SplineFits.png', bbox_inches='tight', dpi = 400)

###################################################################

pifFile = None
xmlFile = None
l_height, l_width = -1, -1
xmlTime = -1
testCell = -1
plotfits, comparefits = None, None

parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-x", "--xml", action="store", type="string", dest="xmlfile", help="path to XML file", metavar="XML")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")
parser.add_option("-p","--plot", action="store_true", dest="plotfits", help="plot fitted geometry", default=False)
parser.add_option("-c","--compare", action="store_true", dest="comparefits", help="compare fitted geometry", default=False)
parser.add_option("-t","--time", action="store", type="int", dest="time", help="time in XML file to compare", metavar="TIME")
parser.add_option("-T","--test", action="store", type="int", dest="testcell", help="runs the test code on specified cell")

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

plotXML = xmlFile is not None and xmlTime != -1

pifFileName = ''

if os.path.isfile(pifFile):
	pifFileName = os.path.splitext(pifFile)[0]
else:
	print("Error: PIF file does not exist.\n")
	exit()

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
		
# Extract features from a list of pixels representing a cell

class ExtractFeatures:

	def __init__(self, cell_pixel_list):
	
		self.pix_list = cell_pixel_list

		self.cell_img = None # Binary image of full cell
		self.perim_img = None # Binary image of cell perimeter
		self.eroded_img = None # Binary image of eroded cell perimeter

		self.perim_coord = None # Coordinates of the cell perimeter pixels
		self.perim_coord_dp = None # Coordinates of approximate perimeter (Douglas-Peucker)
		self.perim_coord_poly = None # Coordinates of polygon (derived from 3pv)
		self.perim_coord_eroded = None # Coordinates of eroded polygon

		self.area_cell = None # Area of cell by pixel counting
		self.area_poly = None # Area of polygon derived from 3pv
		self.area_eroded = None # Area of eroded polygon

		self.perim_3pv = None # Perimeter of cell using 3pv method
		self.perim_1sqrt2 = None # 1 sqrt 2 method
		self.perim_poly = None # Perimeter of polygon derived from 3pv
		self.perim_eroded = None # Perimeter of polygon of eroded cell

		self.equiv_diameter = None # The equivalent diameter of a circle with same area as cell
		self.shape_factor = None

		self.ellipse_fvector = None
		self.ccm_fvector = None

		self.spl_poly = None # Spline tck variables approximating 3pv-polygon

		self.cell_to_image()

	def cell_to_image(self):
	
		# Find x, y coordinate bounds
		x_res = max(self.pix_list, key=itemgetter(0))[0]
		y_res = max(self.pix_list, key=itemgetter(1))[1]

		# Creating labeled_img
		self.cell_img = NP.zeros([x_res+2, y_res+2], dtype=NP.int_)
		
		for (x_pix, y_pix) in self.pix_list:
			self.cell_img[x_pix-1, y_pix-1] = 1

		# Find the pixels that make up the perimeter
		eroded_image = NDI.binary_erosion(self.cell_img)
		eroded_image2 = NDI.binary_erosion(eroded_image)

		# self.perim_img = self.cell_img - eroded_image
		self.eroded_img = eroded_image - eroded_image2
		self.perim_img = self.cell_img - eroded_image

		# Create a list of the coordinates of the pixels (use the center of the pixels)
		perim_image_ind = NP.where(self.perim_img == 1)
		perim_image_coord = NP.array([perim_image_ind[0], perim_image_ind[1]])
		self.perim_coord = NP.transpose(perim_image_coord)

		return

	def basic_props(self):
		'''
		Description: Calculates the perimeter and area using basic methods. For perimeter,
		we use the 3pv, 1 sqrt2 method, and look at the 3pv-polygon perimeter. For area,
		we use pixel counting, and look at the 3pv-polygon area.

		For 3pv perimeter: Use three-pixel vector method to compute perimeter and shape factor
		Reference: http://www.sciencedirect.com/science/article/pii/0308912687902458
		'''
		# Perimeter: 3pv and polygon perimeter (polygon from 3pv)
		self.perim_3pv, self.perim_poly, self.perim_coord_poly = perimeter_3pvm(self.perim_img)
		_, self.perim_eroded, self.perim_coord_eroded = perimeter_3pvm(self.eroded_img)

		# Perimeter: Approximate polygon using Douglas-Peucker algorithm
		self.perim_coord_dp = skimage.measure.approximate_polygon(NP.array(self.perim_coord_poly), 0.75)

		# Create cubic spline
		self.spl_poly, _ = interpolate.splprep(NP.transpose(self.perim_coord_poly), per=1)

		# Perimeter: 1 sqrt2 (function from regionprops)
		props = skimage.measure.regionprops(self.cell_img)
		self.perim_1sqrt2 = props[0].perimeter

		# Area: Pixel Counting
		self.area_cell = len(self.pix_list)	

		# Area: Polygon area (from 3pv polygon)
		# Extract x and y coordinates
		# We subtract 0.5 because PLT.imshow() shows coordinates as the centers of pixels
		# Using the shoelace formula: https://en.wikipedia.org/wiki/Shoelace_formula
		YY = self.perim_coord_poly[:,0]
		XX = self.perim_coord_poly[:,1]

		self.area_poly = 0.5*NP.abs(NP.dot(XX,NP.roll(YY,1))-NP.dot(YY,NP.roll(XX,1)))

		# Area: Eroded polygon area
		YY = self.perim_coord_eroded[:,0]
		XX = self.perim_coord_eroded[:,1]

		self.area_eroded = 0.5*NP.abs(NP.dot(XX,NP.roll(YY,1))-NP.dot(YY,NP.roll(XX,1)))

		# Equivalent Diameter: The diameter of the circle with the same area (pixel counted) as the cell
		self.equiv_diameter = NP.sqrt(4*self.area_cell/NP.pi)
		
		return
		
	def shape_factor(self):
		# TO DO: Check what shape factor gives us
	
		if self.perim_3pv is None:
			perimeter(self)
		if self.area_cell is None:
			area(self)
			
		self.shape_factor = self.perim_3pv/(4*NP.pi*self.area_cell)
		return 

	def ellipse_props(self):
		'''
		Description: Returns list of properties derived from fitting ellipse (in the following order)
		centroid_x, centroid_y, eccentricity, eulerNumber, extent, majorAxisLength,
		minorAxisLength, orientation, perimeter, solidity, ellipseArea and ellipsePerimeter.

		This uses regionprops() fom skimage.measure. The ellipse fit is done by
		fitting an ellipse with the same second central moment as the image. By looking
		at the code, this is done by calculating the inertia tensor of the matrix,
		finding the eigenvalues (the second central moments using the principal axes),
		and matching those with the equations for second central moment of an ellipse.

		Reference: https://en.wikipedia.org/wiki/Image_moment

		Ellipse perimeter: Equation given in https://en.wikipedia.org/wiki/Ellipse#Circumference
		The elliptic integral of the second kind implemented in scipy: 
		http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe
		Note that the scipy definition of the integral differs slightly than wiki, so we take E(e^2) rather than E(e).
		'''

		props = skimage.measure.regionprops(self.cell_img)

		centroid = props[0].centroid

		ellipse_prop_list = [centroid[0]]
		ellipse_prop_list.append(centroid[1])
		ellipse_prop_list.append(props[0].eccentricity)
		ellipse_prop_list.append(props[0].euler_number)
		ellipse_prop_list.append(props[0].extent) # Ratio of pixels in the region to pixels in the total bounding box
		ellipse_prop_list.append(props[0].major_axis_length)
		ellipse_prop_list.append(props[0].minor_axis_length)
		ellipse_prop_list.append(props[0].orientation) # In degrees starting from the x-axis
		ellipse_prop_list.append(props[0].solidity) # Ratio of pixels in the region to pixels of the convex hull image
		ellipse_prop_list.append(NP.pi*ellipse_prop_list[5]*ellipse_prop_list[6]/4.0) # Ellipse area
		ellipse_prop_list.append(2.0*ellipse_prop_list[5]*scipy.special.ellipe(ellipse_prop_list[2]**2)) # Ellipse perimeter

		self.ellipse_fvector = ellipse_prop_list
		return
		
	def cell_centre_fit(self):
		'''
		Description: Returns a list of features derived from fitting a circle (in the following order):
		centroid_x, centroid_y, radius, perimeter, area.

		This uses a least-squares estimator for the circle, using the points on the boundary of the cell.
		These points are chosen to be at the center of the boundary pixels.
		'''

		c_model = skimage.measure.CircleModel()
		c_model.estimate(self.perim_coord)

		if skimage.__version__ == '0.9.3':
			(xc, yc, r) = c_model._params	
		else:									# For newer versions
			(xc, yc, r) = c_model.params

		cell_centre_features = [xc]
		cell_centre_features.append(yc)
		cell_centre_features.append(r)
		cell_centre_features.append(2*NP.pi*r)
		cell_centre_features.append(NP.pi*r**2)
	
		self.ccm_fvector = cell_centre_features
		return

# Check if lattice contains isolated cells
def contains_isolated_cells():
	'''
	Description: This returns true if lattice_data contains more than one connected component
	and false otherwise. This currently uses 1-connectivity.

	Reference: http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label
	'''

	global lattice_data
	global lattice_matrix
	
	clipped_lattice_data = NP.clip(lattice_data,0,1)

	[labeled_img, num_labels] = NDI.measurements.label(clipped_lattice_data)

	if num_labels > 1:
		return True
	else:
		return False
	
# Compute featues for all cells
featureDict = dict()
polyPtDict = dict()
splinePtDict = dict()

spl_u = NP.linspace(0,1, 100) # Spline Parameter

for cell_id in cellDict.keys():
	'''
	The following is code to test implementation of features.
	Only calculate features for the desired test cell.
	'''
	if testCell != -1 and cell_id != testCell:
		continue
	elif testCell != -1 and cell_id == testCell:
		extractor = ExtractFeatures(cellDict[cell_id])
		testFunction(extractor)

	extractor = ExtractFeatures(cellDict[cell_id])

	thread_list = []

	thread_list.append(threading.Thread(target=extractor.basic_props(), args=(), kwargs={}))
	thread_list.append(threading.Thread(target=extractor.ellipse_props(), args=(), kwargs={}))
	thread_list.append(threading.Thread(target=extractor.cell_centre_fit(), args=(), kwargs={}))
	
	for thread in thread_list:
		thread.start()
	
	for thread in thread_list:
		thread.join()
	
	featureDict[cell_id] = [extractor.area_cell, extractor.perim_1sqrt2, extractor.equiv_diameter]
	featureDict[cell_id] = featureDict[cell_id] +  [extractor.perim_3pv]
	featureDict[cell_id] = featureDict[cell_id] + [extractor.perim_poly, extractor.area_poly]
	featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_fvector + extractor.ccm_fvector
	featureDict[cell_id] = featureDict[cell_id] + [extractor.perim_eroded,  extractor.area_eroded]

	polyPtDict[cell_id] = extractor.perim_coord_poly

	# Need flipud() to get x and y coordinates back in (x,y) form
	splinePtDict[cell_id] = NP.transpose(NP.flipud(interpolate.splev(spl_u, extractor.spl_poly)))

# End program here if we're just testing features
if testCell != -1:
	sys.exit()

# Construct featIndexDict
featIndexDict = dict(BASIC=None, ELLIPSE=None, CCM=None, TPV=None, ERODED=None, POLY=None)
BASIC_numfeat = 3
TPV_numfeat = 1
POLY_numfeat = 2
ELLIPSE_numfeat = 11
CCM_numfeat = 5
ERODED_numfeat = 2

TPV_start = BASIC_numfeat
POLY_start = TPV_start + TPV_numfeat
ELLIPSE_start = POLY_start + POLY_numfeat
CCM_start = ELLIPSE_start + ELLIPSE_numfeat
ERODED_start = CCM_start + CCM_numfeat

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

# Plot ellipsoid fit, cell-centre spherical fit, minimum perimeter polygon (MPP) fit
if plotfits:

	# Plot original cells
	lattice_matrix = NP.ascontiguousarray(NP.flipud(NP.transpose(lattice_matrix)))
	pilImage = Image.frombuffer('RGBA', (l_width, l_height), lattice_matrix, 'raw', 'RGBA', 0, 1)
	pilImage = pilImage.convert('RGB')
	pilImage.save(pifFileName + '_Boundary.png')
	
	# Plot polygonized lattice
	sp = None
	if not contains_isolated_cells():
		dirname = os.path.dirname(os.path.abspath(__file__))
		cmd = "python3 vectorize.py --file " + pifFile + " --size " + str(l_width) + "," + str(l_height) + " --output " + dirname
		args = shlex.split(cmd)
		sp = subprocess.Popen(args)

	# Plot ellipse fits
	fig = PLT.figure(1)
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
		ax.add_artist(e)
		e.set_clip_box(ax.bbox)
		e.set_alpha(0.3)
		if ctype == 'CellU':
			e.set_facecolor([0,1,0])
		elif ctype == 'CellV':
			e.set_facecolor([1,0,0])
		else:
			e.set_facecolor([0,0,1])
	
	ax.set_xlim([0,l_width])
	ax.set_ylim([0,l_height])
	PLT.savefig(pifFileName + '_EllipseFit.png', bbox_inches='tight', dpi = 400)
	
	# Plot cell-centre model (disk fit)
	fig = PLT.figure(2)
	ax = fig.add_subplot(111, aspect='equal')

	circles = [[PLT.Circle((featureDict[cell_id][featIndexDict['CCM']['centroid_x']],
		featureDict[cell_id][featIndexDict['CCM']['centroid_y']]),
		featureDict[cell_id][featIndexDict['CCM']['radius']]),
		cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

	for circle in circles:
		c = circle[0]
		ctype = circle[1]
		ax.add_artist(c)
		c.set_alpha(0.3)
		if ctype == 'CellU':
			c.set_facecolor([0,1,0])
		elif ctype == 'CellV':
			c.set_facecolor([1,0,0])
		else:
			c.set_facecolor([0,0,1])

	ax.set_xlim([0,l_width])
	ax.set_ylim([0,l_height])
	PLT.savefig(pifFileName + '_CircleFit.png', bbox_inches='tight', dpi = 400)

	# Plot 3pv-polygon model
	fig = PLT.figure(3)
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

	ax.set_xlim([0,l_width])
	ax.set_ylim([0,l_height])
	PLT.savefig(pifFileName + '_3pvPolyFit.png', bbox_inches='tight', dpi = 400)
	
	# Plot spline model
	fig = PLT.figure(4)
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

	ax.set_xlim([0,l_width])
	ax.set_ylim([0,l_height])
	PLT.savefig(pifFileName + '_SplineFit.png', bbox_inches='tight', dpi = 400)
	
	if sp is not None:
		sp.wait()
	
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
	perim_1sqrt2 = NP.array(perim_1sqrt2)[perim_reorder_ind]
	perim_ellispe = NP.array(perim_ellipse)[perim_reorder_ind]
	perim_circle = NP.array(perim_circle)[perim_reorder_ind]
	perim_poly = NP.array(perim_poly)[perim_reorder_ind]
	perim_eroded = NP.array(perim_eroded)[perim_reorder_ind]
	perim_3pv = NP.array(perim_3pv)[perim_reorder_ind]

	area_reorder_ind = NP.argsort(area_basic)
	area_ellipse = NP.array(area_ellipse)[area_reorder_ind]
	area_circle = NP.array(area_circle)[area_reorder_ind]
	area_poly = NP.array(area_poly)[area_reorder_ind]
	area_eroded = NP.array(area_eroded)[area_reorder_ind]
	area_basic = NP.array(area_circle)[area_reorder_ind]

	if plotXML:
		perim_xml = NP.array(perim_xml)[perim_reorder_ind]
		area_xml = NP.array(area_xml)[area_reorder_ind]

	cell_range = range(len(xData))

	# Plot the figures
	PLT.figure(1)
	if plotXML:
		PLT.plot(cell_range, perim_xml)
	PLT.plot(cell_range, perim_ellipse)
	PLT.plot(cell_range, perim_circle)
	PLT.plot(cell_range, perim_poly)
	PLT.plot(cell_range, perim_eroded)
	PLT.plot(cell_range, perim_1sqrt2)
	PLT.plot(cell_range, perim_3pv)

	PLT.xlabel('Cell (arbitrary)')
	PLT.ylabel('Perimeter')
	PLT.title('Comparison of Perimeter vs. Cell')
	if plotXML:
		PLT.legend(["xml", "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "1, sqrt(2) method", "3pv method"], loc=2)
	else:
		PLT.legend(["Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "1, sqrt(2) method", "3pv method"], loc=2)	
	PLT.savefig(pifFileName + '_PerimCompare.png', bbox_inches='tight', dpi = 400)

	PLT.figure(2)
	if plotXML:
		PLT.plot(cell_range, area_xml)
	PLT.plot(cell_range, area_ellipse)
	PLT.plot(cell_range, area_circle)
	PLT.plot(cell_range, area_poly)
	PLT.plot(cell_range, area_eroded)
	PLT.plot(cell_range, area_basic)

	PLT.xlabel('Cell (arbitrary)')
	PLT.ylabel('Area')
	PLT.title('Comparison of Area vs. Cell')
	if plotXML:
		PLT.legend([ "xml", "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "Pixel counting"], loc=2)
	else:
		PLT.legend([ "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "Pixel counting"], loc=2)		
	PLT.savefig(pifFileName + '_AreaCompare.png', bbox_inches='tight', dpi = 400)

