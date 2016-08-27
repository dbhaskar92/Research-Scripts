#!/usr/bin/env python

#
# Last modified: 08 August 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Build a vector of labelled cell images for convolutional neural network from multiple PIF files
# 

import os
import math
import pickle
import argparse
import numpy as NP
from PIL import Image
from operator import itemgetter

parser = argparse.ArgumentParser(description='Labelled cell images for CNN classifier')
parser.add_argument('pifs', metavar='PIF FILE', type=file, nargs='+', help='list of paths to PIF files')

args = parser.parse_args()

parseImgVector = []
cellLabelVector = []
res_size = 0

for pifFile in args.pifs:

	lattice = pifFile.read().split('\n')[1:-1]
	cellDict = dict()
	cellTypeDict = dict()

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
	
	
	for cell_id in cellDict.keys():

		if cellTypeDict[cell_id] != 'CellU' and cellTypeDict[cell_id] != 'CellV':
			continue
		
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
		if cellTypeDict[cell_id] == 'CellU':
			cellLabelVector.append(0)
		elif cellTypeDict[cell_id] == 'CellV':
			cellLabelVector.append(1)

# Check vector lengths
if len(cellLabelVector) != len(parseImgVector):
	print "Error: Vector length mismatch."	
		
# Convert images to same size (print for debugging)
if res_size % 2 != 0:
	res_size = res_size + 1
	
cnt = 0
if not os.path.exists('CellImgs'):
    os.makedirs('CellImgs')

cellImgVector = []
for img in parseImgVector:
	r, c = img.shape
	padded = NP.zeros([res_size, res_size], dtype=NP.int_)
	r_start = res_size/2 - math.floor(r/2)
	c_start = res_size/2 - math.floor(c/2)
	padded[r_start:r_start+r, c_start:c_start+c] = img
	cellImgVector.append(padded)
	im = Image.fromarray(padded.astype('uint8')*255)
	im.save('CellImgs/' + str(cnt).zfill(4) + '.png')
	cnt = cnt + 1

# Serialize
output = open('labelled_data.pkl', 'wb')
pickle.dump(cellImgVector, output)
pickle.dump(cellLabelVector, output)
output.close()
