#!/usr/bin/env python

#
# Last modified: 26 May 2016
# Author: Darrick Lee <y.l.darrick@gmail.com>, Dhananjay Bhaskar <dbhaskar92@gmail.com>
# A wrapper for pif_extract_features.py to execute for a set of PIF files.
# 
# Requriements: vectorize.py


import os
import glob
import time
import threading
import subprocess
import numpy as NP
from optparse import OptionParser


pifFolder = ''
startTime, endTime = None, None
threads = 1
pifFolder = None

parser = OptionParser()
parser.add_option("-p", "--path", action="store", type="string", dest="pifFolder", help="path to folder with PIF files", metavar="PIF")
parser.add_option("--start", action="store", type="int", dest="start", help="first time to process", metavar="START")
parser.add_option("--end", action="store", type="int", dest="end", help="last time to process", metavar="END")
parser.add_option('-t', "--threads", action="store", type="int", dest="threads", help="number of threads to run", metavar="THREAD")
parser.add_option("-o","--output", action="store", type="string", dest="outputfolder", help="the folder to store output plots", metavar="OUTPUT")
parser.add_option('-r', "--rewrite", action="store_true", dest="rewrite", help="rewrite old files", default=False)

# Options parsing
(options, args) = parser.parse_args()
if options.pifFolder:
	pifFolder = options.pifFolder
if options.start is not None:
	startTime = options.start
else:
	startTime = 0
if options.end is not None:
	endTime = options.end
else:
	endTime = 9999999
if options.threads:
	threads = options.threads
if options.outputfolder:
	outFolder = options.outputfolder
else:
	outFolder = pifFolder
if options.rewrite:
	rewrite = True
else:
	rewrite = False

boundaryFolder = outFolder + "BoundaryFit/"

def run_pif_extract_features(pifList):
	# Function that executes pif_extract_features.py for all given pif files
	for pif in pifList:
		pifName = os.path.splitext(pif)[0]
		pifName_nopath = pifName.split('/')[-1]
		
		boundaryFile = boundaryFolder + pifName_nopath + '_Boundary.png'

		# Check if the file already exists, only run the function if the file doesn't exist yet
		if not os.path.isfile(boundaryFile) or rewrite:
			subprocess.call(['python', 'pif_extract_features.py',
				'-i', pif,
				'-o', outFolder,
				'-l', '800',
				'-w', '800',
				'-p', '-s'])

t0 = time.time()

# Create a list of all pif file names
if pifFolder:
	if os.path.isdir(pifFolder):
		allPif = glob.glob(pifFolder + '*.pif')
	else:
		print("Error: Directory does not exist.\n")
		exit()
else:
	allPif = glob.glob('*.pif')

# Create equal partitions of allPif for each thread
numPif_total = len(allPif)
if endTime > numPif_total:
	endTime = numPif_total

numPif = endTime - startTime

allPif = NP.array(allPif)
pifPartition = dict()
allPifNum = NP.array([int(pif.split('/')[-1].split('.')[0]) for pif in allPif])

for i in range(threads):
	part_startTime = startTime + i*numPif/(threads)
	part_endTime = startTime + (i+1)*numPif/(threads)

	inRange = (allPifNum >= part_startTime) * (allPifNum < part_endTime)
	pifPartition[i] = list(allPif[inRange])

# Ensure output directory exist; if not create it
if not os.path.isdir(outFolder):
	os.mkdir(outFolder)

# Ensure output folders exist; if not, create them
circleOut = outFolder + 'CircleFit/'
ellipseOut = outFolder + 'EllipseFit/'
splineOut = outFolder + 'SplineFit/'
polyOut = outFolder + 'PolyFit/'
boundaryOut = outFolder + 'BoundaryFit/'
vectorizeOut = outFolder + 'vectorizeFit/'
splineBdyOut = outFolder + 'SplineBdyFit/'

if not os.path.isdir(circleOut):
	os.mkdir(circleOut)
if not os.path.isdir(ellipseOut):
	os.mkdir(ellipseOut)
if not os.path.isdir(splineOut):
	os.mkdir(splineOut)
if not os.path.isdir(polyOut):
	os.mkdir(polyOut)
if not os.path.isdir(boundaryOut):
	os.mkdir(boundaryOut)
if not os.path.isdir(vectorizeOut):
	os.mkdir(vectorizeOut)
if not os.path.isdir(splineBdyOut):
	os.mkdir(splineBdyOut)

# Split the job up into the given number of threads
thread_list = []

for partNum in pifPartition:
	thread_list.append(threading.Thread(target=run_pif_extract_features, args=(pifPartition[partNum],)))

for thread in thread_list:
	thread.start()

for thread in thread_list:
	thread.join()

t1 = time.time()

print("Total runtime is {0} seconds.").format(t1-t0)


