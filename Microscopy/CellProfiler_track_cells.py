#!/usr/bin/env python

#
# Last modified: 31 May 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Requires package: pip install sortedcontainers
#

import sys
import csv
import collections
import numpy as NP
import matplotlib.pyplot as PLT
from scipy.misc import imread
from sortedcontainers import SortedSet

if (len(sys.argv) != 2):
	print 'Usage: python CellProfiler_plot_statistics.py /path/to/AllMyExpt_MyCells.csv'
	sys.exit()

# hold in memory
frame_path_map = {}

# parse file
with open(sys.argv[1]) as csvfile:

	reader = csv.DictReader(csvfile)
	
	for row in reader:
		frame_path_map[int(row['Metadata_FrameNumber'])] = row['Metadata_FileLocation']
		

# track all cells identified at the given time (frame)
def track_cells(time):

	frames = SortedSet()						# List of image frames

	c_area_map = collections.defaultdict(list)		# Cell Area
	c_aXY_map = collections.defaultdict(list)		# AreaShape_Center
	c_cmi_map = collections.defaultdict(list)		# Location_CenterMassIntensity
	c_cnt_map = collections.defaultdict(list)		# Location_Center

	fh = open(sys.argv[1], "rb", 10485760)
	lines = fh.readlines()
	linereader = csv.DictReader(lines)
	
	centroid_X = []
	c_color_map = {}
	
	for row in linereader:
	
		if int(row['Metadata_FrameNumber']) == time:
		
			centroid_X.append(float(row['AreaShape_Center_X']))
	
	arr = NP.asarray(centroid_X)
	bins = []
	
	for i in range(10, 100, 10):
		bins.append(NP.percentile(arr, i))
	
	linereader = csv.DictReader(lines)
	for row in linereader:
	
		if int(row['Metadata_FrameNumber']) == time:
		
			cent_X = float(row['AreaShape_Center_X'])
			c_num = int(row['ObjectNumber'])
			if cent_X < bins[0]:
				c_color_map[c_num] = 'royalblue'
			elif cent_X < bins[1]:
				c_color_map[c_num] = 'orangered'
			elif cent_X < bins[2]:
				c_color_map[c_num] = 'lime'
			elif cent_X < bins[3]:
				c_color_map[c_num] = 'fuchsia'
			elif cent_X < bins[4]:
				c_color_map[c_num] = 'crimson'
			elif cent_X < bins[5]:
				c_color_map[c_num] = 'dodgerblue'
			elif cent_X < bins[6]:
				c_color_map[c_num] = 'darksalmon'
			elif cent_X < bins[7]:
				c_color_map[c_num] = 'gold'
			elif cent_X < bins[8]:
				c_color_map[c_num] = 'aqua'
			else:
				c_color_map[c_num] = 'magenta'
	
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
	plotFrame = frames[len(frames)-1]
	imgFile = frame_path_map[plotFrame]
	print "DEBUG Largest track for time t = " + str(time) + " ends at frame: " + str(plotFrame) + " image file: " + imgFile + "\n"
	
	cnt = 0
	
	for i in range(time, plotFrame+1):
		
		img = "./OutlineCells/OutlineCells" + "{0:0>3}".format(i) + ".png"
		fig = PLT.figure(cnt)
		axes = PLT.gca()
		axes.set_xlim([0,1600])
		axes.set_ylim([0,1200])
		axes.invert_yaxis()
		axes.xaxis.tick_top()
		axes.yaxis.tick_left()
		bg = imread(img)
		PLT.imshow(bg, zorder=0)
		
		for cid in c_aXY_map.keys():
		
			xdata = []
			ydata = []
		
			if cnt == 0:
				[x, y] = c_aXY_map[cid][0]
				PLT.scatter(x, y, s=4, zorder=1)
			else:
				for k in range(0, cnt+1):
					try:
						xdata.append(c_aXY_map[cid][k][0])
						ydata.append(c_aXY_map[cid][k][1])
					except IndexError:
						break
				if len(xdata) == cnt+1:
					lines = PLT.plot(xdata, ydata, zorder=1)
					PLT.setp(lines, 'color', c_color_map[cid], 'linewidth', 1.0)
		
		PLT.savefig("Track_" + "{0:0>3}".format(i) + ".png", bbox_inches='tight', dpi=200)
		PLT.close(fig)
			
		cnt = cnt + 1
	

	# Generate Displacement and Velocity Histogram
	# TODO
			
	
							
						
def get_cell_future(obj_num, time, area_future, aXY_future, cmi_future, cnt_future):

	with open(sys.argv[1]) as csvh:
		
		hreader = csv.DictReader(csvh)
	
		curr_lifetime = -1
		curr_frame = -1
		curr_obj_num = -1
		num_daughters = 0

		for data in hreader:

			if int(data['Metadata_FrameNumber']) < time:
				continue
		
			elif int(data['Metadata_FrameNumber']) == time and int(data['ObjectNumber']) == obj_num:
				curr_lifetime = int(data['TrackObjects_Lifetime_30'])
				curr_frame = time
				curr_obj_num = obj_num
	
			else:
		
				c_num = int(data['ObjectNumber'])
	
				if curr_lifetime == -1 or curr_frame == -1 or curr_obj_num == -1:
					continue
	
				t = int(data['Metadata_FrameNumber'])
				parentid = int(data['TrackObjects_ParentObjectNumber_30'])
				lifetime = int(data['TrackObjects_Lifetime_30'])
			
				if t > curr_frame + 1:
					break
		
				elif parentid == curr_obj_num and lifetime == curr_lifetime + 1 and t == curr_frame + 1:
			
					area_future[obj_num].append([int(data['AreaShape_Area'])])
			
					aXY_future[obj_num].append([float(data['AreaShape_Center_X']), float(data['AreaShape_Center_Y'])])
		
					cmi_X = float(data['Location_CenterMassIntensity_X_Outlines'])
					cmi_Y = float(data['Location_CenterMassIntensity_Y_Outlines'])
					cmi_future[obj_num].append([cmi_X, cmi_Y])
		
					cnt_future[obj_num].append([float(row['Location_Center_X']), float(row['Location_Center_Y'])])
			
					curr_lifetime = lifetime
					curr_frame = t
					curr_obj_num = c_num
				
				elif parentid == curr_obj_num and lifetime == 1 and t == curr_frame + 1:
			
					num_daughters = num_daughters + 1

				else:
					continue
			
	print "DEBUG Cell id: " + str(obj_num) + " Frames found: " + str(curr_frame - time + 1) + " Num daughters: " + str(num_daughters)
	return curr_frame

	
print "RUNNING TESTS\n"
print "Tracking time t=10\n"
track_cells(10)
