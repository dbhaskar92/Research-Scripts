#!/usr/bin/env python

#
# Last modified: 8 June 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Requires package: pip install sortedcontainers
#

import sys
import csv
import math
import collections
from scipy.misc import imread
from sortedcontainers import SortedSet
from matplotlib import collections as MC
import numpy as NP
import matplotlib.cm as cmx
import matplotlib.pyplot as PLT
import matplotlib.colors as colors

if (len(sys.argv) != 2):
	print 'Usage: python CellProfiler_plot_statistics.py /path/to/AllMyExpt_MyCells.csv'
	sys.exit()


# map frames to image path
frame_path_map = {}

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
	
	for i in range(20, 100, 20):
		bins.append(NP.percentile(arr, i))
	
	linereader = csv.DictReader(lines)
	for row in linereader:
	
		if int(row['Metadata_FrameNumber']) == time:
		
			cent_X = float(row['AreaShape_Center_X'])
			c_num = int(row['ObjectNumber'])
			if cent_X < bins[0]:
				c_color_map[c_num] = 'lime'
			elif cent_X < bins[1]:
				c_color_map[c_num] = 'dodgerblue'
			elif cent_X < bins[2]:
				c_color_map[c_num] = 'gold'
			elif cent_X < bins[3]:
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
	
					
	# Print debug information
	plotFrame = frames[len(frames)-1]
	imgFile = frame_path_map[plotFrame]
	print "DEBUG Largest track for time t = " + str(time) + " ends at frame: " + str(plotFrame) + " image file: " + imgFile + "\n"
	
	
	figcnt = 0
	
	# Tracks over segmented cell image
	cnt = 0
	lineclr = 'aqua'						# aqua for red channel, gold for green channel 
	for i in range(time, plotFrame+1):
		
		img = "./OutlineCells/OutlineCells" + "{0:0>3}".format(i) + ".png"
		fig = PLT.figure(figcnt)
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
				PLT.scatter(x, y, color=lineclr, s=4, zorder=1)
			else:
				for k in range(0, cnt+1):
					try:
						xdata.append(c_aXY_map[cid][k][0])
						ydata.append(c_aXY_map[cid][k][1])
					except IndexError:
						break
				if len(xdata) == cnt+1:
					lines = PLT.plot(xdata, ydata, zorder=1)
					PLT.setp(lines, 'color', lineclr, 'linewidth', 1.0)
					PLT.scatter(xdata[-1], ydata[-1], color=lineclr, s=4, zorder=2)
		
		PLT.savefig("Track_" + "{0:0>3}".format(i) + ".png", bbox_inches='tight', dpi=200)
		PLT.close(fig)
		figcnt = figcnt + 1	
		cnt = cnt + 1
	
	
	# Tracks color coded by horizontal position
	cnt = 0
	for i in range(time, plotFrame+1):
		
		fig = PLT.figure(figcnt)
		axes = PLT.gca()
		axes.set_xlim([0,1600])
		axes.set_ylim([0,1200])
		axes.invert_yaxis()
		axes.xaxis.tick_top()
		axes.yaxis.tick_left()
		
		for cid in c_aXY_map.keys():
		
			xdata = []
			ydata = []
		
			if cnt == 0:
				[x, y] = c_aXY_map[cid][0]
				PLT.scatter(x, y, color=c_color_map[cid], s=4)
			else:
				for k in range(0, cnt+1):
					try:
						xdata.append(c_aXY_map[cid][k][0])
						ydata.append(c_aXY_map[cid][k][1])
					except IndexError:
						break
				lines = PLT.plot(xdata, ydata)
				PLT.setp(lines, 'color', c_color_map[cid], 'linewidth', 1.0)
				PLT.scatter(xdata[-1], ydata[-1], color=c_color_map[cid], s=4)
		
		PLT.savefig("Color_" + "{0:0>3}".format(i) + ".png", bbox_inches='tight', dpi=200)
		PLT.close(fig)
		figcnt = figcnt + 1	
		cnt = cnt + 1


	# Calculate displacement, distance, avg. and instantaneous velocity
	distance_map = {}
	displacement_map = {}
	avg_velocity_map = {}
	inst_velocity_map = collections.defaultdict(list)
	out_csv_file = "Velocity_Frame" + "{0:0>3}".format(time) + ".csv"
	
	time_conv_factor = 0.2		# 5 minutes time interval
	dist_conv_factor = 0.8		# 0.8 microns distance
	
	with open(out_csv_file, 'w') as csvfile:
		
		fieldnames = ['ObjectNumber', 'Velocity_X', 'Velocity_Y']
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()
	
		for cid in c_aXY_map.keys():
		
			prev_x = None
			prev_y = None
			pp_x = None
			pp_y = None
			first_x = None
			first_y = None
		
			num_velocity_vecs = 0
		
			for [cur_x, cur_y] in c_aXY_map[cid]:
		
				if first_x is None and first_y is None:
					first_x = cur_x
					first_y = cur_y
				elif prev_x is None and prev_y is None:
					distance_map[cid] = 0
					displacement_map[cid] = 0
					prev_x = cur_x
					prev_y = cur_y
				elif pp_x is None and pp_y is None:
					distance_map[cid] = (math.hypot(cur_x - first_x, cur_y - first_y))*dist_conv_factor
					displacement_map[cid] = displacement_map[cid] + (math.hypot(cur_x - prev_x, cur_y - prev_y))*dist_conv_factor
					avg_velocity_map[cid] = 0
					
					pp_x = prev_x
					pp_y = prev_y
					prev_x = cur_x
					prev_y = cur_y
				else:
					distance_map[cid] = (math.hypot(cur_x - first_x, cur_y - first_y))*dist_conv_factor
					displacement_map[cid] = displacement_map[cid] + (math.hypot(cur_x - prev_x, cur_y - prev_y))*dist_conv_factor
					new_velocity_vec = (math.hypot(cur_x - pp_x, cur_y - pp_y))*dist_conv_factor*time_conv_factor*0.5
					inst_velocity_map[cid].append([cur_x - pp_x, cur_y - pp_y])
					if num_velocity_vecs == 0:
						writer.writerow({'ObjectNumber': cid, 'Velocity_X': cur_x - pp_x, 'Velocity_Y': cur_y - pp_y})
					avg_velocity_map[cid] = ((avg_velocity_map[cid]*num_velocity_vecs) + new_velocity_vec)/(num_velocity_vecs + 1)
					pp_x = prev_x
					pp_y = prev_y
					prev_x = cur_x
					prev_y = cur_y
					num_velocity_vecs = num_velocity_vecs + 1
	
	
	# Tracks color coded by speed
	cnt = 0
	drawscale = 0.1
	for i in range(time, plotFrame+1):
		
		fig = PLT.figure(figcnt)
		axes = PLT.gca()
		axes.set_xlim([0,1600])
		axes.set_ylim([0,1200])
		axes.invert_yaxis()
		axes.xaxis.tick_top()
		axes.yaxis.tick_left()
		
		speed = []
		for cid in c_aXY_map.keys():
			for v_vec in inst_velocity_map[cid]:
				speed.append((math.hypot(v_vec[0], v_vec[1]))*dist_conv_factor*time_conv_factor*0.5)
		
		speed.sort()
		
		jet = cm = PLT.get_cmap('jet') 
		cNorm  = colors.Normalize(vmin=0, vmax=speed[-1])
		scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
		
		c_f = 0.2						# scaling factor for velocity vector
		
		for cid in c_aXY_map.keys():
		
			xd = []
			yd = []
			ud = []
			vd = []			
			color_list = []

			if cnt == 0:
				try:
					[x, y] = c_aXY_map[cid][0]
					[u, v] = inst_velocity_map[cid][0]
					PLT.quiver(x, y, u*c_f, v*c_f, color='limegreen', angles='xy', scale_units='xy', scale=drawscale)
				except IndexError:
					continue
			else:
				for k in range(0, cnt+1):
					try:
						xd.append(c_aXY_map[cid][k][0])
						yd.append(c_aXY_map[cid][k][1])
						vx = inst_velocity_map[cid][k][0]
						vy = inst_velocity_map[cid][k][1]
						ud.append(vx)
						vd.append(vy)
						color_list.append(scalarMap.to_rgba((math.hypot(vx, vy))*dist_conv_factor*time_conv_factor*0.5))
					except IndexError:
						break
				try:
					PLT.quiver(xd[-1], yd[-1], ud[-1]*c_f, vd[-1]*c_f, color='limegreen', angles='xy', scale_units='xy', scale=drawscale)
					points = NP.array([xd, yd]).T.reshape(-1, 1, 2)
					segments = NP.concatenate([points[:-1], points[1:]], axis=1)
					lc = MC.LineCollection(segments, colors=color_list, linewidths=1)
					PLT.gca().add_collection(lc)
				except IndexError:
					continue
		
		PLT.savefig("Velocity_" + "{0:0>3}".format(i) + ".png", bbox_inches='tight', dpi=200)
		PLT.close(fig)
		figcnt = figcnt + 1	
		cnt = cnt + 1
	
				
	# Plot histograms
	
	fig = PLT.figure(figcnt)
	dist_hist_data = collections.defaultdict(list)
	max_dist = 0
	for cid in distance_map.keys():
		dist_hist_data[c_color_map[cid]].append(distance_map[cid])
		if distance_map[cid] > max_dist:
			max_dist = distance_map[cid]
	bins = NP.linspace(0, max_dist, num=math.ceil(math.sqrt(len(distance_map.keys()))))
	f, axarr = PLT.subplots(len(dist_hist_data.keys()), sharex=True)
	ind = 0
	for clr in dist_hist_data.keys():
		axarr[ind].hist(dist_hist_data[clr], bins, normed=False, cumulative=False, color=clr)
		axarr[0].set_title('Cell Distance Travelled Histogram')
		ind = ind + 1
		figcnt = figcnt + 1
	axarr[ind-1].set_xlabel(r'Distance ($\mu m$)')
	PLT.tight_layout()
	PLT.savefig('DistanceHist.png')
	
	fig = PLT.figure(figcnt)
	PLT.title('Cell Displacement Histogram')
	PLT.xlabel(r'Displacement ($\mu m$)')
	PLT.ylabel('Number of Cells')
	disp_hist_data = collections.defaultdict(list)
	max_disp = 0
	for cid in displacement_map.keys():
		disp_hist_data[c_color_map[cid]].append(displacement_map[cid])
		if displacement_map[cid] > max_disp:
			max_disp = displacement_map[cid]
	bins = NP.linspace(0, max_disp, num=math.ceil(math.sqrt(len(displacement_map.keys()))))
	for clr in disp_hist_data.keys():
		PLT.hist(disp_hist_data[clr], bins, normed=False, cumulative=False, color=clr)
	PLT.tight_layout()
	PLT.savefig('DisplacementHist.png')
	figcnt = figcnt + 1
	
	fig = PLT.figure(figcnt)
	PLT.title('Cell Avg. Speed Histogram')
	PLT.xlabel(r'Speed ($\mu m$ per min)')
	PLT.ylabel('Number of Cells')
	speed_hist_data = collections.defaultdict(list)
	max_speed = 0
	for cid in avg_velocity_map.keys():
		speed_hist_data[c_color_map[cid]].append(avg_velocity_map[cid])
		if avg_velocity_map[cid] > max_speed:
			max_speed = avg_velocity_map[cid]
	bins = NP.linspace(0, max_speed, num=math.ceil(math.sqrt(len(avg_velocity_map.keys()))))
	for clr in speed_hist_data.keys():
		PLT.hist(speed_hist_data[clr], bins, normed=False, cumulative=False, color=clr)
	PLT.tight_layout()
	PLT.savefig('SpeedHist.png')
	figcnt = figcnt + 1	

							
# Helper function						
def get_cell_future(obj_num, time, area_future, aXY_future, cmi_future, cnt_future):

	with open(sys.argv[1]) as csvh:
		
		hreader = csv.DictReader(csvh)
	
		curr_lifetime = -1
		curr_frame = -1
		curr_obj_num = -1
		num_daughters = 0
		
		# For testing
		cutoff = 500

		for data in hreader:

			if int(data['Metadata_FrameNumber']) < time:
				continue
				
			elif int(data['Metadata_FrameNumber']) > time + cutoff:
				break
		
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


# Main	
print "RUNNING TESTS\n"
print "Tracking time t=10\n"
track_cells(10)
