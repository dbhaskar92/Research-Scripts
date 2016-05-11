#!/usr/bin/env python

#
# Last modified: 5 May 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Extract cell shape features from PIF file
# 

import math
import string
import numpy as NP
import skimage.measure
import scipy.special
import matplotlib.pyplot as PLT
from matplotlib.patches import Ellipse
from operator import itemgetter
from optparse import OptionParser

pifFile = ''
l_height = -1
l_width = -1

parser = OptionParser()

parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")
parser.add_option("-p","--plot", action="store_true", dest="plotfits", default=False, help="plot the various fits")

(options, args) = parser.parse_args()

if options.inputfile:
	pifFile = options.inputfile
if options.height:
	l_height = options.height
if options.width:
	l_width = options.width

if options.plotfits:
	plotfits = 1
else:
	plotfits = 0
	
lattice = open(pifFile).read().split('\n')[1:-1]

cellDict = dict()

for pix in lattice:

	fields = pix.split()

	cell_id = int(fields[0])
	pix_x = int(fields[3])
	pix_y = int(fields[5])
            
	if cell_id not in cellDict.keys():
    
		cellDict[cell_id] = [[pix_x,pix_y]]

	cellDict[cell_id].append([pix_x,pix_y])

# Extract features from a list of pixels representing a cell

class ExtractFeatures:

	def __init__(self, cell_pixel_list):
	
		self.pix_list = cell_pixel_list
		    
	def area(self):
	
		return len(self.pix_list) 

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
		The elliptic integral of the second kind implemented in scipy: http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe
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

		return ellipse_prop_list

# Compute feature vector of each cell for further classification
# Order of features in the feature vector:
# 	0-area
# 	1-centroid_x
# 	2-centroid_y
# 	3-eccentricity
# 	4-eulerNumber
# 	5-extent
# 	6-majorAxisLength
# 	7-minorAxisLength
# 	8-orientation
# 	9-perimeter
# 	10-solidity
#	11-ellipseArea (the area of the ellipse calculating using major/minor axes)
#	12-ellipsePerim 

featureDict = dict()

for cell_id in cellDict.keys():

	extractor = ExtractFeatures(cellDict[cell_id])

	featureDict[cell_id] = [extractor.area()]
	featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_props()
	
print featureDict[1]

if plotfits:
	# Plot original cells
	cell_img = NP.zeros([l_width, l_height])

	for cell_id in cellDict.keys():
		for(x_pix, y_pix) in cellDict[cell_id]:
			# Note: I think imshow() plots the axes backwards, so thats what x and y are flipped
			cell_img[y_pix-1, x_pix-1] = 1.0

	fig = PLT.figure(0)
	ax = fig.add_subplot(111, aspect='equal')
	ax.imshow(cell_img)

	# Plot ellipse fits
	ells = [Ellipse(xy=NP.array([featureDict[cell_id][1],featureDict[cell_id][2]]),
		width = featureDict[cell_id][7], height = featureDict[cell_id][6],
		angle = featureDict[cell_id][8]/(2*NP.pi)*360) for cell_id in cellDict.keys()]

	# ax = fig.add_subplot(111, aspect='equal')
	for e in ells:
	    ax.add_artist(e)
	    e.set_clip_box(ax.bbox)
	    e.set_alpha(0.3)
	    e.set_facecolor([0,1,0])

	PLT.show()


