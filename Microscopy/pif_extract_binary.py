#!/usr/bin/env python

#
# Last modified: March 2017
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Extract binary cell images from multiple PIF files
# 

import os
import math
import argparse
import numpy as NP
from PIL import Image
from operator import itemgetter

parser = argparse.ArgumentParser(description='Extract binary cell images from PIF files')
parser.add_argument('pifs', metavar='PIF FILE', type=file, nargs='+', help='list of paths to PIF files')

args = parser.parse_args()

parseImgVector = []
cellIdentifierVector = []
res_size = 0

# Parse PIF files
for pifFile in args.pifs:

	lattice = pifFile.read().split('\n')[1:-1]
	imgName = os.path.splitext(os.path.basename(pifFile.name))[0]
	
	cellDict = dict()

	for pix in lattice:

		fields = pix.split()

		cell_id = int(fields[0])
		c_type = str(fields[1])
		pix_x = int(fields[3])
		pix_y = int(fields[5])
				
		if cell_id not in cellDict.keys():
			cellDict[cell_id] = [[pix_x, pix_y]]

		cellDict[cell_id].append([pix_x, pix_y])
	
	
	for cell_id in cellDict.keys():
		
		x_max = max(cellDict[cell_id], key=itemgetter(0))[0]
		y_max = max(cellDict[cell_id], key=itemgetter(1))[1]
		
		x_min = min(cellDict[cell_id], key=itemgetter(0))[0]
		y_min = min(cellDict[cell_id], key=itemgetter(1))[1]
		
		width = x_max - x_min + 1
		height = y_max - y_min + 1
		
		cell_img = NP.zeros([width, height], dtype=NP.int_)
		res_size = max(width, height, res_size)
			
		for (x_pix, y_pix) in cellDict[cell_id]:
			cell_img[x_pix-x_min, y_pix-y_min] = 1
			
		parseImgVector.append(cell_img)
		cellIdentifierVector.append([imgName, cell_id])


# Check vector lengths
if len(cellIdentifierVector) != len(parseImgVector):
	print "Error: Vector length mismatch."	
		
# Convert images to same size
if res_size % 2 != 0:
	res_size = res_size + 1

# Save images
for cnt in range(len(parseImgVector)):

	img = parseImgVector[cnt]
	folderName = cellIdentifierVector[cnt][0]
	cid = str(cellIdentifierVector[cnt][1])

	r, c = img.shape
	padded = NP.zeros([res_size, res_size], dtype=NP.int_)
	
	r_start = res_size/2 - math.floor(r/2)
	c_start = res_size/2 - math.floor(c/2)
	
	padded[r_start:r_start+r, c_start:c_start+c] = img
	
	im = Image.fromarray(padded.astype('uint8')*255)
	
	if not os.path.exists(folderName):
		os.makedirs(folderName)
	
	im.save(folderName + '/' + cid + '.png')
