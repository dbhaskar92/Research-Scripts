#!/usr/bin/env python
#
# Last modified: 16 May 2016
# Author: Darrick Lee <y.l.darrick@gmail.com>
# Description: This function uses the 3-pixel vector method described in the paper
# 	below to calculate the perimeter of binary image.
# Reference: http://www.sciencedirect.com/science/article/pii/0308912687902458
#
# NOTE: If the total number of pixels in perim_img is odd, the perimeter MAY be off by +/- 1
#	This still needs to be checked more thoroughly, but this function should be accurate
#	for the large majority of images.

import numpy as NP
from scipy import ndimage as NDI

def perimeter_3pvm(perim_img):
	# Calculate (clockwise) chain code (convention in paper)
	# Initialize potential movement directions (corresponding to chain code convention in paper)
	dir_y = NP.array([0, -1, -1, -1, 0, 1, 1, 1])
	dir_x = NP.array([1, 1, 0, -1, -1, -1, 0, 1])

	# First, find the leftmost, uppermost pixel
	# NOTE: When we consider an image as a picture, index 0 corresponds to the
	# y-axis whereas index 1 corresponds to the x-axis
	perim_img_ind = NP.where(perim_img == 1)
	leftmost_ind = min(perim_img_ind[1])
	leftmost_pix = NP.where(perim_img_ind[1] == leftmost_ind)[0]
	uppermost_ind = max(perim_img_ind[0][leftmost_pix])

	cur_y = first_y = uppermost_ind
	cur_x = first_x = leftmost_ind

	bordering_values = perim_img[cur_y + dir_y, cur_x + dir_x]
	possible_directions = NP.where(bordering_values==1)[0]

	next_dir = NP.sort(NP.mod(possible_directions-6,8))[0] + 6
	chain_code = [NP.mod(next_dir,8)]

	cur_y = cur_y + dir_y[chain_code[-1]]
	cur_x = cur_x + dir_x[chain_code[-1]]

	# Continue to append to chain code until we reach the beginning
	while [cur_y, cur_x] != [first_y, first_x]:
		bordering_values = perim_img[cur_y + dir_y, cur_x + dir_x]
		possible_directions = NP.where(bordering_values==1)[0]

		'''
		The following line picks the next direction to go in the chain code.
		The first thing to note is that it the next direction must be going in the
		FORWARD direction.
		The second thing to note is that the chain code is defined to be going in the
		CLOCKWISE direction. This means that the first direction to check is mod(chain_code[-1]-2,8),
		which corresponds to going left relative to the previous direction. If this direction
		is not in possible_directions, then we check the subsequent numbers in ascending order.

		This is done by sorting and choosing the minimum element.
		'''
		next_dir = NP.sort(NP.mod(possible_directions-chain_code[-1]+2,8))[0] + chain_code[-1] -2
		chain_code.append(NP.mod(next_dir,8))

		cur_y = cur_y + dir_y[chain_code[-1]]
		cur_x = cur_x + dir_x[chain_code[-1]]

	# Define l_b, D, and delta_l matrices
	def D(i, n):
		D1 = {
			1: 4*(4-n),
			2: 4*(4-n),
			3: NP.mod(18-4*n,16),
			4: NP.mod(14-4*n,16),
			5: NP.mod(17-4*n,16),
			6: NP.mod(17-4*n,16),
			7: NP.mod(17-4*n,16),
			8: NP.mod(17-4*n,16),
			9: 15-4*n,
			10: 15-4*n,
			11: 15-4*n,
			12: 15-4*n,
			13: 15-4*n
		}.get(i)

		D2 = {
			1: 4*(4-n),
			2: 4*(4-n),
			3: NP.mod(18-4*n,16),
			4: NP.mod(14-4*n,16),
			5: NP.mod(17-4*n,16),
			6: 13-4*n,
			7: NP.mod(27-4*n,16),
			8: NP.mod(25-4*n,16),
			9: NP.mod(19-4*n,16),
			10: 15-4*n,
			11: NP.mod(27-4*n,16),
			12: NP.mod(25-4*n,16),
			13: NP.mod(23-4*n,16)
		}.get(i)

		return (D1, D2)

	l_b = NP.array(
		[NP.sqrt(5), NP.sqrt(5), NP.sqrt(5), NP.sqrt(5), 2, 3, 3+NP.sqrt(2), 4,
		NP.sqrt(2), 2*NP.sqrt(2), 1+2*NP.sqrt(2), 2+NP.sqrt(2), 2+2*NP.sqrt(2)])

	delta_l = NP.array([
		[0,0,0,0,NP.sqrt(2)-NP.sqrt(5),1-NP.sqrt(2),2-NP.sqrt(5),0,2,2,2,2,1,1,1,1],
		[0,0,NP.sqrt(2)-NP.sqrt(5),NP.sqrt(2)-NP.sqrt(5),2-5,0,0,2,2,2,1,1,1,1,0,0],
		[0,-1,-1,2-NP.sqrt(5),0,0,0,0,2,1,1,1,1,0,0,0],
		[0,0,0,NP.sqrt(2)-NP.sqrt(5),NP.sqrt(5)-2*NP.sqrt(2),2-NP.sqrt(5),0,2,2,2,2,1,1,1,1,0]])

	# Calculate perimeter
	num_chain = len(chain_code)
	isodd = NP.mod(num_chain,2) == 1

	if isodd:
		chain_pair = zip(chain_code[0:-1:2],chain_code[1:-1:2])
	else:
		chain_pair = zip(chain_code[0::2],chain_code[1::2])

	perimeter = 0

	for c1, c2 in chain_pair:
		n = c1//2

		# p_type = pattern type
		# g_type = group type
		if NP.mod(c1,2)==0:
			# Even c1 cases
			if c2 == 2*n+1:
				p_type = 1
				g_type = 1
			elif c2 == NP.mod(2*n+7,8):
				p_type = 3
				g_type = 2
			elif c2 == 2*n:
				p_type = 5
				g_type = 3
			elif c2 == NP.mod(2*n+2,8):
				p_type = 6
				g_type = 3
			elif c2 == NP.mod(2*n+3,8):
				p_type = 7
				g_type = 3
			elif c2 == NP.mod(2*n+4,8):
				p_type = 8
				g_type  = 3
		else:
			# Odd c1 cases
			if c2 == 2*n:
				p_type = 2
				g_type = 1
			elif c2 == NP.mod(2*n+2,8):
				p_type = 4
				g_type = 2
			elif c2 == NP.mod(2*n+7,8):
				p_type = 9
				g_type = 4
			elif c2 == 2*n+1:
				p_type = 10
				g_type = 4
			elif c2 == NP.mod(2*n+3,8):
				p_type = 11
				g_type = 4
			elif c2 == NP.mod(2*n+4,8):
				p_type = 12
				g_type = 4
			elif c2 == NP.mod(2*n+5,8):
				p_type = 13
				g_type = 4

		if perimeter == 0:
			perimeter = l_b[p_type-1]

			D1_last, D2_prev = D(p_type, n)

			# This is used to calculate the final delta_l term 
			D1_first = D1_last
			g_type_first = g_type
		else:
			D1_cur, D2_cur = D(p_type, n)
			Delta_D =NP.mod(D1_cur - D2_prev + 16, 16) + 1

			perimeter = perimeter + l_b[p_type-1] + delta_l[g_type-1][Delta_D-1]
			D2_prev = D2_cur

	# D2_prev and g_type is from the iteration of for loop
	Delta_D = NP.mod(D1_first - D2_prev + 16, 16) + 1
	perimeter = perimeter + delta_l[g_type_first-1][Delta_D-1]

	if isodd:
		if(NP.mod(chain_code[-1],2) == 0):
			perimeter = perimeter + 1
		else:
			perimeter = perimeter + NP.sqrt(2)
	else:
		# D2_prev and g_type is from the iteration of for loop
		Delta_D =NP.mod(D1_first - D2_prev + 16, 16) + 1
		perimeter = perimeter + delta_l[g_type-1][Delta_D-1]

	return perimeter