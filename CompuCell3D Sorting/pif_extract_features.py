#!/usr/bin/env python

#
# Last modified: 11 May 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Extract cell shape features from PIF file
# 

import os
import math
import string
import threading
import scipy.special
import skimage.measure
import numpy as NP
from PIL import Image
from operator import itemgetter
from optparse import OptionParser
import matplotlib.pyplot as PLT
from matplotlib.patches import Ellipse

pifFile = ''
l_height = -1
l_width = -1

parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")
parser.add_option("-p","--plot", action="store_true", dest="plotfits", help="plot fitted geometry", default=False)

(options, args) = parser.parse_args()
if options.inputfile:
	pifFile = options.inputfile
if options.height:
	l_height = options.height
if options.width:
	l_width = options.width
if options.plotfits:
	plotfits = True
else:
	plotfits = False

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
		self.cell_area = None
		self.cell_perimeter = None
		self.shape_factor = None
		self.ellipse_fvector = None
		self.mpp_fvector = None
		self.cc_fvector = None
		    
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

		# Find x, y coordinate bounds
		x_res = max(self.pix_list, key=itemgetter(0))[0]
		y_res = max(self.pix_list, key=itemgetter(1))[1]

		# Creating labeled_img
		labeled_img = NP.zeros([x_res, y_res])
		
		for (x_pix, y_pix) in self.pix_list:
			labeled_img[x_pix-1, y_pix-1] = 1

		props = skimage.measure.regionprops(NP.int_(labeled_img))

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
	
		cell_centre_features = []
		
		# TODO: Fit a disk, compute its radius, perimeter and goodness of fit
		
		self.cc_fvector = cell_centre_features
		return

# Check if lattice contains isolated cells

def contains_isolated_cells():
	
	global lattice_data
	
	# TODO: return true iff lattice_matrix contains isolated cells
	
	return True
	
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
	
	for thred in thread_list:
		thread.join()
	
	featureDict[cell_id] = [extractor.cell_area]
	featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_fvector

# Plot ellipsoid fit, cell-centre spherical fit, minimum perimeter polygon (MPP) fit

if plotfits:

	# Plot original cells
	
	lattice_matrix = NP.ascontiguousarray(NP.flipud(NP.transpose(lattice_matrix)))
	pilImage = Image.frombuffer('RGBA', (l_width, l_height), lattice_matrix, 'raw', 'RGBA', 0, 1)
	pilImage = pilImage.convert('RGB')
	pilImage.save('PIFimage.png')
	
	# Plot polygonized lattice
	
	# TODO

	# Plot ellipse fits
	fig = PLT.figure(1)
	ax = fig.add_subplot(111, aspect='equal')
	ellipse_objs = dict()

	# Plot ellipse fits
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
	PLT.savefig('EllipseFit.png')
	
	# Plot cell-centre model (disk fit)
	
	# TODO
	
	# Plot MPP fit
	
	# TODO
	

# Compare features 

# Plot centroid location
