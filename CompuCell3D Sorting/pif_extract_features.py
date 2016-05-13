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
import matplotlib.pyplot as PLT
from PIL import Image
from operator import itemgetter
from optparse import OptionParser
from matplotlib.patches import Ellipse
from scipy import ndimage as NDI

pifFile = None
xmlFile = None
l_height, l_width = -1, -1
plotfits, comparefits = None, None

parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-x", "--xml", action="store", type="string", dest="xmlfile", help="path to XML file", metavar="XML")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")
parser.add_option("-p","--plot", action="store_true", dest="plotfits", help="plot fitted geometry", default=False)
parser.add_option("-c","--compare", action="store_true", dest="comparefits", help="compare fitted geometry", default=False)

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
		self.labeled_img = None
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
		self.labeled_img = NP.zeros([x_res, y_res], dtype=NP.int_)
		
		for (x_pix, y_pix) in self.pix_list:
			self.labeled_img[x_pix-1, y_pix-1] = 1
		return

		    
	def area(self):
	
		self.cell_area = len(self.pix_list)
		return

	def perimeter(self):
	
		'''
		Description: Use three-pixel vector method to compute perimeter and shape factor
		Reference: http://www.sciencedirect.com/science/article/pii/0308912687902458
		'''
		
		# TODO: Implement method described in paper
		
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

		props = skimage.measure.regionprops(self.labeled_img)

		centroid = props[0].centroid

		ellipse_prop_list = [centroid[0]]
		ellipse_prop_list.append(centroid[1])
		ellipse_prop_list.append(props[0].eccentricity)
		ellipse_prop_list.append(props[0].euler_number)
		ellipse_prop_list.append(props[0].extent)
		ellipse_prop_list.append(props[0].major_axis_length)
		ellipse_prop_list.append(props[0].minor_axis_length)
		ellipse_prop_list.append(props[0].orientation)
		ellipse_prop_list.append(props[0].perimeter)
		ellipse_prop_list.append(props[0].solidity)
		ellipse_prop_list.append(NP.pi*ellipse_prop_list[5]*ellipse_prop_list[6]/4.0)
		ellipse_prop_list.append(2.0*ellipse_prop_list[5]*scipy.special.ellipe(ellipse_prop_list[2]**2))

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
		
		# First find the pixels that make up the perimeter
		eroded_image = NDI.binary_erosion(self.labeled_img)
		perim_image = self.labeled_img - eroded_image

		# Create a list of the coordinates of the pixels (use the center of the pixels)
		perim_image_ind = NP.where(perim_image == 1)
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
		
		# TODO: Fit a disk, compute its radius, perimeter and goodness of fit
	
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
	
	featureDict[cell_id] = [extractor.cell_area]
	featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_fvector + extractor.ccm_fvector

# Construct featIndexDict
featIndexDict = dict(BASIC=None, ELLIPSE=None, CCM=None, MPP=None)
BASIC_numfeat = 1
ELLIPSE_numfeat = 12
CCM_numfeat = 5

ELLIPSE_start = BASIC_numfeat
CCM_start = BASIC_numfeat + ELLIPSE_numfeat

featIndexDict['BASIC'] = dict(area = 0)
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
	area=CCM_start+3,
	perimeter=CCM_start+4)

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
# if comparefits:

	# Plot centroid location
	
	# Plot perimeter estimates
	
	# Plot area estimates
