#!/usr/bin/env python

#
# Last modified: 11 May 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Extract cell shape features from PIF file
# 

import os
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
from matplotlib.patches import Ellipse
from scipy import ndimage as NDI
from perimeter_3pvm import perimeter_3pvm

pifFile = None
xmlFile = None
l_height, l_width = -1, -1
xmlTime = -1
plotfits, comparefits = None, None

parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-x", "--xml", action="store", type="string", dest="xmlfile", help="path to XML file", metavar="XML")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")
parser.add_option("-p","--plot", action="store_true", dest="plotfits", help="plot fitted geometry", default=False)
parser.add_option("-c","--compare", action="store_true", dest="comparefits", help="compare fitted geometry", default=False)
parser.add_option("-t","--time", action="store", type="int", dest="time", help="time in XML file to compare", metavar="TIME")

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
if options.time is not None:
	xmlTime = options.time
else:
	plotfits = False
if options.comparefits:
	comparefits = True
else:
	comparefits = False

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
		self.bin_img = None
		self.perim_img = None
		self.cell_area = None
		self.cell_perimeter = None
		self.shape_factor = None
		self.ellipse_fvector = None
		self.mpp_fvector = None
		self.ccm_fvector = None

		self.cell_to_image()

	def cell_to_image(self):
	
		# Find x, y coordinate bounds
		x_res = max(self.pix_list, key=itemgetter(0))[0]
		y_res = max(self.pix_list, key=itemgetter(1))[1]

		# Creating labeled_img
		self.bin_img = NP.zeros([x_res+2, y_res+2], dtype=NP.int_)
		
		for (x_pix, y_pix) in self.pix_list:
			self.bin_img[x_pix-1, y_pix-1] = 1

		# Find the pixels that make up the perimeter
		eroded_image = NDI.binary_erosion(self.bin_img)
		self.perim_img = self.bin_img - eroded_image

		return

	def area(self):
	
		self.cell_area = len(self.pix_list)
		return

	def perimeter(self):
		'''
		Description: Use three-pixel vector method to compute perimeter and shape factor
		Reference: http://www.sciencedirect.com/science/article/pii/0308912687902458
		'''
		self.cell_perimeter = perimeter_3pvm(self.perim_img)
		
		return
		
	def shape_factor(self):
	
		if self.cell_perimeter is None:
			perimeter(self)
		if self.cell_area is None:
			area(self)
			
		self.shape_factor = self.cell_perimeter/(4*NP.pi*self.cell_area)
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

		props = skimage.measure.regionprops(self.bin_img)

		centroid = props[0].centroid

		ellipse_prop_list = [centroid[0]]
		ellipse_prop_list.append(centroid[1])
		ellipse_prop_list.append(props[0].eccentricity)
		ellipse_prop_list.append(props[0].euler_number)
		ellipse_prop_list.append(props[0].extent) #Ratio of pixels in the region to pixels in the total bounding box
		ellipse_prop_list.append(props[0].major_axis_length)
		ellipse_prop_list.append(props[0].minor_axis_length)
		ellipse_prop_list.append(props[0].orientation) # In degrees starting from the x-axis
		ellipse_prop_list.append(props[0].perimeter) # Uses 1, sqrt(2) method using center of boundary pixels
		ellipse_prop_list.append(props[0].solidity) # Ratio of pixels in the region to pixels of the convex hull image
		ellipse_prop_list.append(NP.pi*ellipse_prop_list[5]*ellipse_prop_list[6]/4.0) # Ellipse area
		ellipse_prop_list.append(2.0*ellipse_prop_list[5]*scipy.special.ellipe(ellipse_prop_list[2]**2)) # Ellipse perimeter

		self.ellipse_fvector = ellipse_prop_list
		return
		
	def minimum_perimeter_polygon(self):
	
		MPP_feature_list = []
		
		# TODO: Compute minimum perimeter polygon
	
		self.mpp_fvector = MPP_feature_list
		return
		
	def cell_centre_fit(self):
		'''
		Description: Returns a list of features derived from fitting a circle (in the following order):
		centroid_x, centroid_y, radius, perimeter, area.

		This uses a least-squares estimator for the circle, using the points on the boundary of the cell.
		These points are chosen to be at the center of the boundary pixels.
		'''

		# Create a list of the coordinates of the pixels (use the center of the pixels)
		perim_image_ind = NP.where(self.perim_img == 1)
		perim_image_coord = NP.array([perim_image_ind[0], perim_image_ind[1]])
		perim_image_coord = NP.transpose(perim_image_coord)
		perim_image_coord = perim_image_coord + 0.5

		c_model = skimage.measure.CircleModel()
		c_model.estimate(perim_image_coord)

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

for cell_id in cellDict.keys():

	extractor = ExtractFeatures(cellDict[cell_id])

	thread_list = []

	thread_list.append(threading.Thread(target=extractor.area(), args=(), kwargs={}))
	thread_list.append(threading.Thread(target=extractor.ellipse_props(), args=(), kwargs={}))
	thread_list.append(threading.Thread(target=extractor.perimeter(), args=(), kwargs={}))
	thread_list.append(threading.Thread(target=extractor.minimum_perimeter_polygon(), args=(), kwargs={}))
	thread_list.append(threading.Thread(target=extractor.cell_centre_fit(), args=(), kwargs={}))
	
	for thread in thread_list:
		thread.start()
	
	for thread in thread_list:
		thread.join()
	
	featureDict[cell_id] = [extractor.cell_area, extractor.cell_perimeter, extractor.shape_factor]
	featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_fvector + extractor.ccm_fvector

# Construct featIndexDict
featIndexDict = dict(BASIC=None, ELLIPSE=None, CCM=None, MPP=None)
BASIC_numfeat = 3
ELLIPSE_numfeat = 12
CCM_numfeat = 5

ELLIPSE_start = BASIC_numfeat
CCM_start = BASIC_numfeat + ELLIPSE_numfeat

featIndexDict['BASIC'] = dict(area = 0,
	perimeter=1,
	shape_factor=2)
featIndexDict['ELLIPSE'] = dict(centroid_x=ELLIPSE_start, 
	centroid_y=ELLIPSE_start+1,
	eccentricity=ELLIPSE_start+2,
	euler_number=ELLIPSE_start+3,
	extent=ELLIPSE_start+4,
	major_axis_length=ELLIPSE_start+5,
	minor_axis_length=ELLIPSE_start+6,
	orientation=ELLIPSE_start+7,
	perimeter=ELLIPSE_start+8,
	solidity=ELLIPSE_start+9,
	area=ELLIPSE_start+10,
	ellipse_perim=ELLIPSE_start+11)
featIndexDict['CCM'] = dict(centroid_x=CCM_start,
	centroid_y=CCM_start+1,
	radius=CCM_start+2,
	perimeter=CCM_start+3,
	area=CCM_start+4)

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

	ells = [[Ellipse(xy=NP.array([featureDict[cell_id][1],featureDict[cell_id][2]]),
		width = featureDict[cell_id][7], height = featureDict[cell_id][6],
		angle = featureDict[cell_id][8]/(2*NP.pi)*360), cellTypeDict[cell_id]] for cell_id in cellDict.keys()]
	
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
	PLT.savefig(pifFileName + '_EllipseFit.png')
	
	# Plot cell-centre model (disk fit)
	fig = PLT.figure(2)
	ax = fig.add_subplot(111, aspect='equal')

	circles = [[PLT.Circle((featureDict[cell_id][13],featureDict[cell_id][14]),
		featureDict[cell_id][15]), cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

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
	PLT.savefig(pifFileName + '_CircleFit.png')
	PLT.show()
	
	# Plot MPP fit
	
	# TODO
	
	if sp is not None:
		sp.wait()
	

# Compare features 
if comparefits:
	# Prepare XML file
	infile = open(xmlFile,'r')
	xml1 = ET.parse(infile)
	root = xml1.getroot()
	times = root.getchildren()
	cells = times[xmlTime].getchildren()

	# Initialize data vectors
	xData = []

	perim_1sqrt2 = []
	perim_ellipse = []
	perim_circle = []
	perim_3pv = []
	perim_xml = []

	area_ellipse = []
	area_circle = []
	area_basic = []
	area_xml = []

	for cell_id in cellDict.keys():
		xData.append(cell_id)

	# Populate the data vectors
	for cell_id in xData:
		perim_1sqrt2.append(featureDict[cell_id][featIndexDict['ELLIPSE']['perimeter']])
		perim_ellipse.append(featureDict[cell_id][featIndexDict['ELLIPSE']['ellipse_perim']])
		perim_circle.append(featureDict[cell_id][featIndexDict['CCM']['perimeter']])
		perim_3pv.append(featureDict[cell_id][featIndexDict['BASIC']['perimeter']])

		perim_xml.append(cells[cell_id].get('perimeter'))

		area_ellipse.append(featureDict[cell_id][featIndexDict['ELLIPSE']['area']])
		area_circle.append(featureDict[cell_id][featIndexDict['CCM']['area']])
		area_basic.append(featureDict[cell_id][featIndexDict['BASIC']['area']])

		area_xml.append(cells[cell_id].get('area'))

	# Reorder the vectors so that xml parameters are sorted from smallest to largest
	# Note: We switch to NP.array rather than list so we can input a list to select elements
	perim_reorder_ind = NP.argsort(perim_xml)
	perim_xml = NP.array(perim_xml)[perim_reorder_ind]
	perim_1sqrt2 = NP.array(perim_1sqrt2)[perim_reorder_ind]
	perim_ellispe = NP.array(perim_ellipse)[perim_reorder_ind]
	perim_circle = NP.array(perim_circle)[perim_reorder_ind]
	perim_3pv = NP.array(perim_3pv)[perim_reorder_ind]

	area_reorder_ind = NP.argsort(area_xml)
	area_xml = NP.array(area_xml)[area_reorder_ind]
	area_ellipse = NP.array(area_ellipse)[area_reorder_ind]
	area_circle = NP.array(area_circle)[area_reorder_ind]

	cell_range = range(len(xData))

	# Plot the figures
	PLT.figure(1)
	PLT.plot(cell_range, perim_xml)
	PLT.plot(cell_range, perim_ellipse)
	PLT.plot(cell_range, perim_circle)
	PLT.plot(cell_range, perim_1sqrt2)
	PLT.plot(cell_range, perim_3pv)

	PLT.xlabel('Cell ID')
	PLT.ylabel('Perimeter')
	PLT.title('Comparison of Perimeter vs. Cell ID')
	PLT.legend(["xml", "Ellipse fit", "CCM fit", "1, sqrt(2) method", "3pv method"], loc=2)
	PLT.savefig('perim_compare.png', bbox_inches='tight', dpi = 400)

	PLT.figure(2)
	PLT.plot(cell_range, area_xml)
	PLT.plot(cell_range, area_ellipse)
	PLT.plot(cell_range, area_circle)
	PLT.plot(cell_range, area_basic)

	PLT.xlabel('Cell ID')
	PLT.ylabel('Area')
	PLT.title('Comparison of Area vs. Cell ID')
	PLT.legend([ "xml", "Ellipse fit", "CCM fit", "Pixel counting"], loc=2)
	PLT.savefig('area_compare.png', bbox_inches='tight', dpi = 400)

	# Plot centroid location
	
	# Plot perimeter estimates
	
	# Plot area estimates
