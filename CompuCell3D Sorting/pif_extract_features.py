#!/usr/bin/env python

#
# Last modified: 5 May 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Extract cell shape features from PIF file
# 

import math
import string
import numpy as NP
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
	
	
# Extract features from a list of pixels representing a cell

class ExtractFeatures:

	def __init__(self, cell_pixel_list):
	
		self.pix_list = cell_pixel_list
		    
	def area(self):
	
		return len(self.pix_list)
	
# Compute feature vector of each cell for further classification

featureDict = dict()

for cell_id in cellDict.keys():

	extractor = ExtractFeatures(cellDict[cell_id])

	featureDict[cell_id] = [extractor.area()]
	
print featureDict[1]
