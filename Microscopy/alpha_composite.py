#!/usr/bin/env python

#
# Last modified: 09 Mar 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Alpha blend two series of images
#

import Image
import glob
import cv2
import os

path1 = "../Output-Green-Enhanced/OutlineCells/*.png"
path2 = "../Output-Red-Enhanced/OutlineCells/*.png"

video = None
vInit = False

frameStart = 0
numFrames = 100

for f1name in sorted(glob.glob(path1))[frameStart:frameStart+numFrames]:
	for f2name in sorted(glob.glob(path2))[frameStart:frameStart+numFrames]:
		
		file1 = os.path.basename(f1name)
		file2 = os.path.basename(f2name)

		if file1 == file2:
			
			background = Image.open(f1name)
			overlay = Image.open(f2name)
			
			# For JPEG inputs
			# background = background.convert("RGBA")	
			# overlay = overlay.convert("RGBA")
			new_img = Image.blend(background, overlay, 0.5)
			new_img.save(file1,"PNG")
			
			if vInit == False:
				img = cv2.imread(file1)
				height, width, layers = img.shape

				fourcc = cv2.cv.CV_FOURCC('m', 'p', '4', 'v')
				video = cv2.VideoWriter('overlay.avi', fourcc, 25, (width, height), True)

				video.write(img)
				vInit = True
			
			if vInit == True:
				video.write(cv2.imread(file1))
				
if video is not None:	
	cv2.destroyAllWindows()
	video.release()
