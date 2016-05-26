#!/usr/bin/env python

#
# Last modified: 26 May 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Requires package: pip install sortedcontainers
#

import sys
import csv
import collections
import numpy as NP
from sortedcontainers import SortedSet

if (len(sys.argv) != 2):
	print 'Usage: python CellProfiler_plot_statistics.py /path/to/AllMyExpt_MyCells.csv'
	sys.exit()

# params
get_time = 5
interval_length = 10

# hold in memory
frame_path_map = {}


# parse file
with open(sys.argv[1]) as csvfile:

	reader = csv.DictReader(csvfile)
	
	for row in reader:
		frame_path_map[int(row['Metadata_FrameNumber'])] = row['Metadata_FileLocation']
		

# track all cells identified at the given time (frame)
def track_cells(time):

	frames = SortedSet()							# List of image frames

	c_area_map = collections.defaultdict(list)		# Cell Area
	c_aXY_map = collections.defaultdict(list)		# AreaShape_Center
	c_cmi_map = collections.defaultdict(list)		# Location_CenterMassIntensity
	c_cnt_map = collections.defaultdict(list)		# Location_Center

	lines = open(sys.argv[1]).readlines()
	linereader = csv.DictReader(lines)
	
	for row in linereader:
	
		if int(row['Metadata_FrameNumber']) == time:
		
			frames.add(time)
		
			c_num = int(row['ObjectNumber'])
			
			c_area_map[c_num].append([int(row['AreaShape_Area'])])
			c_aXY_map[c_num].append([float(row['AreaShape_Center_X']), float(row['AreaShape_Center_Y'])])
			
			cmi_X = float(row['Location_CenterMassIntensity_X_Outlines'])
			cmi_Y = float(row['Location_CenterMassIntensity_Y_Outlines'])
			c_cmi_map[c_num].append([cmi_X, cmi_Y])
			
			c_cnt_map[c_num].append([float(row['Location_Center_X']), float(row['Location_Center_Y'])])
			
			# peer into the future
			max_frame = get_cell_future(c_num, time, c_area_map, c_aXY_map, c_cmi_map, c_cnt_map)
			if max_frame > time:
				frames.add(max_frame)

			print("DEBUG Area Centroid Track:")
			for item in c_aXY_map[c_num]:
				print " ".join(map(str, item[:]))
			print "\n"
				
			# Generate Centroid Tracks
			# TODO

			# Generate Displacement Histogram
			# TODO
			
	print "DEBUG Largest track for time t=" + str(time) + " ends at frame number: " + str(frames[len(frames)-1]) + "\n"
							
						
def get_cell_future(obj_num, time, area_future, aXY_future, cmi_future, cnt_future):
		
	with open(sys.argv[1]) as csvh:
		
		hreader = csv.DictReader(csvh)
	
		curr_lifetime = -1
		curr_frame = -1

		for data in hreader:

			if int(data['Metadata_FrameNumber']) < time:
				continue
			
			elif int(data['Metadata_FrameNumber']) == time and int(data['ObjectNumber']) == obj_num:
				curr_lifetime = int(data['TrackObjects_Lifetime_30'])
				curr_frame = time
		
			elif int(data['ObjectNumber']) == obj_num:
		
				if curr_lifetime == -1 or curr_frame == -1:
					print "ERROR Current lifetime and frame not set."
					break
		
				t = int(data['Metadata_FrameNumber'])
				parentid = int(data['TrackObjects_ParentObjectNumber_30'])
				lifetime = int(data['TrackObjects_Lifetime_30'])
			
				if parentid == obj_num and lifetime == curr_lifetime + 1 and t == curr_frame + 1:
				
					area_future[obj_num].append([int(data['AreaShape_Area'])])
				
					aXY_future[obj_num].append([float(data['AreaShape_Center_X']), float(data['AreaShape_Center_Y'])])
			
					cmi_X = float(data['Location_CenterMassIntensity_X_Outlines'])
					cmi_Y = float(data['Location_CenterMassIntensity_Y_Outlines'])
					cmi_future[obj_num].append([cmi_X, cmi_Y])
			
					cnt_future[obj_num].append([float(row['Location_Center_X']), float(row['Location_Center_Y'])])
				
					curr_lifetime = lifetime
					curr_frame = t

				else:
					break
				
			else:
				continue
			
	print "DEBUG Cell id: " + str(obj_num) + " Frames found: " + str(curr_frame - time + 1) + " Max tracking frame: " + str(curr_frame)
	return curr_frame
	
print "RUNNING TESTS\n"
print "Tracking time t=10\n"
track_cells(10)
