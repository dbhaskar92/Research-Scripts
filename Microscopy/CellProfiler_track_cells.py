#!/usr/bin/env python

#
# Last modified: 13 June 2016
# Author: Dhananjay Bhaskar <dbhaskar92@gmail.com>
#

import sys
import csv
import math
import collections
from scipy.misc import imread
from scipy.signal import savgol_filter
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
				c_color_map[c_num] = 'forestgreen'
			elif cent_X < bins[1]:
				c_color_map[c_num] = 'purple'
			elif cent_X < bins[2]:
				c_color_map[c_num] = 'darkorange'
			elif cent_X < bins[3]:
				c_color_map[c_num] = 'mediumblue'
			else:
				c_color_map[c_num] = 'magenta'
	
	linereader = csv.DictReader(lines)
	division_events = {}
	plotFrame = track_all_cells(time, linereader, c_area_map, c_aXY_map, c_cmi_map, c_cnt_map, division_events)
					
	# Print debug information
	imgFile = frame_path_map[plotFrame]
	print "DEBUG Largest track for time t = " + str(time) + " ends at frame: " + str(plotFrame) + " image file: " + imgFile + "\n"
	
	
	figcnt = 0
	
	# Tracks over segmented cell image
	cnt = 0
	lineclr = 'dodgerblue'						# dodgerblue for green channel, goldenrod for red channel 
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
	avg_speed_map = {}
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
			position_x_vec = []
			position_y_vec = []
		
			num_velocity_vecs = 0
		
			for [cur_x, cur_y] in c_aXY_map[cid]:
			
				position_x_vec.append(cur_x)
				position_y_vec.append(cur_y)
		
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
					pp_x = prev_x
					pp_y = prev_y
					prev_x = cur_x
					prev_y = cur_y
				else:
					distance_map[cid] = (math.hypot(cur_x - first_x, cur_y - first_y))*dist_conv_factor
					displacement_map[cid] = displacement_map[cid] + (math.hypot(cur_x - prev_x, cur_y - prev_y))*dist_conv_factor
					if num_velocity_vecs == 0:
						writer.writerow({'ObjectNumber': cid, 'Velocity_X': cur_x - pp_x, 'Velocity_Y': cur_y - pp_y})
					pp_x = prev_x
					pp_y = prev_y
					prev_x = cur_x
					prev_y = cur_y
					num_velocity_vecs = num_velocity_vecs + 1
					
			if len(position_x_vec) >= 5 and len(position_y_vec) >= 5:
				velocity_x_vec = savgol_filter(position_x_vec, 5, 3, deriv=1, mode='nearest')
				velocity_y_vec = savgol_filter(position_y_vec, 5, 3, deriv=1, mode='nearest')
				speed = []
				for [v_x, v_y] in zip(velocity_x_vec, velocity_y_vec):
					inst_velocity_map[cid].append([v_x, v_y])
					speed.append((math.hypot(v_x, v_y))*dist_conv_factor*time_conv_factor*0.5)
				avg_speed_map[cid] = NP.mean(speed)
	
	
	# Tracks color coded by speed
	cnt = 0
	dscale = 0.1		# drawing scale
	c_f = 0.1			# scaling factor for velocity vector
	
	speed = []
	for cid in c_aXY_map.keys():
		for v_vec in inst_velocity_map[cid]:
			speed.append((math.hypot(v_vec[0], v_vec[1]))*dist_conv_factor*time_conv_factor*0.5)
	
	speed.sort()
	
	jet = cm = PLT.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=speed[0], vmax=speed[-1])
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	
	lineclr = 'dodgerblue'		# dodgerblue for green channel, goldenrod for red channel 
	
	for i in range(time, plotFrame+1):
		
		fig = PLT.figure(figcnt)
		axes = PLT.gca()
		axes.set_xlim([0,1600])
		axes.set_ylim([0,1200])
		axes.invert_yaxis()
		axes.xaxis.tick_top()
		axes.yaxis.tick_left()
		
		for cid in c_aXY_map.keys():
		
			xd = []
			yd = []
			ud = []
			vd = []			
			color_list = []

			if cnt == 0:
				try:
					[x, y] = c_aXY_map[cid][cnt]
					[u, v] = inst_velocity_map[cid][cnt-2]
					PLT.scatter(x, y, color=lineclr, s=4)
				except IndexError:
					continue
			else:
				for k in range(0, cnt+1):
					try:
						vx = inst_velocity_map[cid][k][0]
						vy = inst_velocity_map[cid][k][1]
						xd.append(c_aXY_map[cid][k][0])
						yd.append(c_aXY_map[cid][k][1])
						ud.append(vx)
						vd.append(vy)
						color_list.append(scalarMap.to_rgba((math.hypot(vx, vy))*dist_conv_factor*time_conv_factor*0.5))
					except IndexError:
						break
				try:
					if len(xd) == cnt+1:
						PLT.quiver(xd[-1], yd[-1], ud[-1]*c_f, vd[-1]*c_f, color=lineclr, angles='xy', scale_units='xy', scale=dscale)
					else:
						PLT.scatter(xd[-1], yd[-1], color=lineclr, s=4)	
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
	for cid in avg_speed_map.keys():
		speed_hist_data[c_color_map[cid]].append(avg_speed_map[cid])
		if avg_speed_map[cid] > max_speed:
			max_speed = avg_speed_map[cid]
	bins = NP.linspace(0, max_speed, num=math.ceil(math.sqrt(len(avg_speed_map.keys()))))
	for clr in speed_hist_data.keys():
		PLT.hist(speed_hist_data[clr], bins, normed=False, cumulative=False, color=clr)
	PLT.tight_layout()
	PLT.savefig('SpeedHist.png')
	figcnt = figcnt + 1

							
# Helper function						
def track_all_cells(time, linereader, c_area_map, c_aXY_map, c_cmi_map, c_cnt_map, division_events):

	max_frame = 0
	cutoff = 20
	
	area_map = collections.defaultdict(list)
	aXY_map = collections.defaultdict(list)
	cmi_map = collections.defaultdict(list)
	cnt_map = collections.defaultdict(list)
	
	curr_lifetime = {}
	curr_frame = {}
	obj_parent_num = {}
	parent_last_accessed = {}
	
	for row in linereader:
	
		c_num = int(row['ObjectNumber'])
		t = int(row['Metadata_FrameNumber'])

		if t < time:
			continue
			
		elif t > time + cutoff:
			break
	
		elif t == time:
			
			area_map[(c_num, t)].append([int(row['AreaShape_Area'])])
			aXY_map[(c_num, t)].append([float(row['AreaShape_Center_X']), float(row['AreaShape_Center_Y'])])
			
			cmi_X = float(row['Location_CenterMassIntensity_X_Outlines'])
			cmi_Y = float(row['Location_CenterMassIntensity_Y_Outlines'])
			
			cmi_map[(c_num, t)].append([cmi_X, cmi_Y])
			cnt_map[(c_num, t)].append([float(row['Location_Center_X']), float(row['Location_Center_Y'])])
		
		
			curr_lifetime[(int(row['ObjectNumber']), time)] = int(row['TrackObjects_Lifetime_30'])
			curr_frame[(int(row['ObjectNumber']), time)] = time
			obj_parent_num[(int(row['ObjectNumber']), time)] = -1
			parent_last_accessed[(int(row['ObjectNumber']), time)] = time
			max_frame = time
			
			c_area_map[c_num] = area_map[(c_num, t)]
			c_aXY_map[c_num] = aXY_map[(c_num, t)]
			c_cmi_map[c_num] = cmi_map[(c_num, t)]
			c_cnt_map[c_num] = cnt_map[(c_num, t)]

		else:
	
			parentid = int(row['TrackObjects_ParentObjectNumber_30'])
			lifetime = int(row['TrackObjects_Lifetime_30'])
			
			if (parentid, t-1) not in obj_parent_num:
				continue
	        
			elif lifetime == curr_lifetime[(parentid, t-1)] + 1 and t == curr_frame[(parentid, t-1)] + 1:
				
				# Note: Cannot detect entity merges. One cell is retained. 
				# Merge is indistinguishable from losing one object and continuing to recognize the other
				# In future, we will disable merges in Cell Profiler 
			
				# If cell division occurs
				if parent_last_accessed[(parentid, t-1)] == t:
					division_events[t-1] = [parentid, aXY_map[(parentid, t-1)][-1]]
					# TODO:
					# create new original id
					# set new origin as parent
					# pad upto time t-1 with (-1, -1)
					# divided cell has no entry in colour map, pick a separate color to visualize
					continue
				
				area_map[(c_num, t)] = area_map[(parentid, t-1)]
				area_map[(c_num, t)].append([int(row['AreaShape_Area'])])
				aXY_map[(c_num, t)] = aXY_map[(parentid, t-1)]
				aXY_map[(c_num, t)].append([float(row['AreaShape_Center_X']), float(row['AreaShape_Center_Y'])])
			
				cmi_X = float(row['Location_CenterMassIntensity_X_Outlines'])
				cmi_Y = float(row['Location_CenterMassIntensity_Y_Outlines'])
			
				cmi_map[(c_num, t)] = cmi_map[(parentid, t-1)]
				cmi_map[(c_num, t)].append([cmi_X, cmi_Y])
				cnt_map[(c_num, t)] = cnt_map[(parentid, t-1)]
				cnt_map[(c_num, t)].append([float(row['Location_Center_X']), float(row['Location_Center_Y'])])
					
				curr_lifetime[(c_num, t)] = lifetime
				curr_frame[(c_num, t)] = t
				obj_parent_num[(c_num, t)] = parentid
				parent_last_accessed[(parentid, t-1)] = t
				parent_last_accessed[(c_num, t)] = t
				
				if t > max_frame:
					max_frame = t
					
				orig_parent = parentid
				orig_time = t-1
				while obj_parent_num[(orig_parent, orig_time)] != -1:
					orig_parent = obj_parent_num[(orig_parent, orig_time)]
					orig_time -= 1
				
				c_area_map[orig_parent] = area_map[(c_num, t)]
				c_aXY_map[orig_parent] = aXY_map[(c_num, t)]
				c_cmi_map[orig_parent] = cmi_map[(c_num, t)]
				c_cnt_map[orig_parent] = cnt_map[(c_num, t)]

			else:
				continue
		
	return max_frame


# Main	
print "RUNNING TESTS\n"
print "Tracking time t=10\n"
track_cells(10)
