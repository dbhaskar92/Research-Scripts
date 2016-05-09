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
from optparse import OptionParser

pifFile = ''
l_height = -1
l_width = -1

parser = OptionParser()

parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")

(options, args) = parser.parse_args()

if options.inputfile:
	pifFile = options.inputfile
if options.height:
	l_height = options.height
if options.width:
	l_width = options.width
	
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

# Find the maximum x and y pixel location (used to create labeled img for regionprops)
all_xpix = [int(pix.split()[3]) for pix in lattice]
all_ypix = [int(pix.split()[5]) for pix in lattice]

xres = max(all_xpix)
yres = max(all_ypix)
	
# Extract features from a list of pixels representing a cell

class ExtractFeatures:

	def __init__(self, cell_pixel_list):
	
		self.pix_list = cell_pixel_list
		    
	def area(self):
	
		return len(self.pix_list) 

	def ellipse_props(self):
		# Returns list of properties derived from fitting ellipse (in the following order)
		# centroid_x, centroid_y, eccentricity, eulerNumber, extent, majorAxisLength,
		# minorAxisLength, orientation, perimeter, solidity

		# Start by creating labeled_img 
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

		return ellipse_prop_list

	
# Compute feature vector of each cell for further classification
# Order of features in the feature vector:
# 	area
# 	centroid_x
# 	centroid_y
# 	eccentricity
# 	eulerNumber
# 	extent
# 	majorAxisLength
# 	minorAxisLength
# 	orientation
# 	perimeter
# 	solidity

featureDict = dict()

for cell_id in cellDict.keys():

	extractor = ExtractFeatures(cellDict[cell_id])

	featureDict[cell_id] = [extractor.area()]
	featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_props()
	
print featureDict[1]
