#!/usr/bin/env python
#
# Last modified: 18 May 2016
# Author: Darrick Lee <y.l.darrick@gmail.com>
# Description: This function uses the 3-pixel vector method described in the paper
# 	below to calculate the perimeter of binary image.
#
#	In addition, it returns an polygon that is approximately the one used for the 3pv-method
#	along with the perimeter of this approximate polygon.
#
# Reference: http://www.sciencedirect.com/science/article/pii/0308912687902458
#
# NOTE: If the total number of pixels in perim_img is odd, the perimeter MAY be off by +/- 1
#	This still needs to be checked more thoroughly, but this function should be accurate
#	for the large majority of images.

import numpy as NP
from scipy import ndimage as NDI

def chain_pair_props(cur_pt, c1, c2):
	'''
	Description: This function calculates and returns properties of the chain pair.
	Specifically, it returns the pattern type (p_type), the group type (g_type),
	and the points to add to the overall polygon. The coordinate of the current pixel
	is assumed to be at the bottom left corner.

	NOTE: As with the rest of the code here, coordinates are labelled (y,x).

	Output: next_pt, p_type, g_type, n, poly_pts
	'''
	# The next two lists are to index the corners of pixels
	# 0 - bottom left
	# 1 - upper left
	# 2 - upper right
	# 3 - bottom right
	c = NP.array([[0,0],[1,0],[1,1],[0,1]])

	# Directions corresponding to chain codes
	d = NP.array([[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1],[1,0],[1,1]])

	n = c1//2
	poly_pts = []

	# p_type = pattern type
	# g_type = group type
	if NP.mod(c1,2)==0:
		# Even c1 cases
		if c2 == 2*n+1:
			p_type = 1
			g_type = 1
			
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c1] + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])

		elif c2 == NP.mod(2*n+7,8):
			p_type = 3
			g_type = 2

			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])
			cur_pt = cur_pt + d[c1] + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])

		elif c2 == 2*n:
			p_type = 5
			g_type = 3

			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])
			cur_pt = cur_pt + d[c1] + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])

		elif c2 == NP.mod(2*n+2,8):
			p_type = 6
			g_type = 3

			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])
			cur_pt = cur_pt + d[c1]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])

		elif c2 == NP.mod(2*n+3,8):
			p_type = 7
			g_type = 3

			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])
			cur_pt = cur_pt + d[c1]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])
			cur_pt = cur_pt + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])

		elif c2 == NP.mod(2*n+4,8):
			p_type = 8
			g_type = 3

			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])
			cur_pt = cur_pt + d[c1]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])
			cur_pt = cur_pt + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])

	else:
		# Odd c1 cases
		if c2 == 2*n:
			p_type = 2
			g_type = 1

			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c1] + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])

		elif c2 == NP.mod(2*n+2,8):
			p_type = 4
			g_type = 2

			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c1] + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])

		elif c2 == NP.mod(2*n+7,8):
			p_type = 9
			g_type = 4

			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)] + 0.5*d[c1])
			cur_pt = cur_pt + d[c1] + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(1+n,4)])

		elif c2 == 2*n+1:
			p_type = 10
			g_type = 4

			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c1] + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])

		elif c2 == NP.mod(2*n+3,8):
			p_type = 11
			g_type = 4

			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c1]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])
			cur_pt = cur_pt + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])

		elif c2 == NP.mod(2*n+4,8):
			p_type = 12
			g_type = 4

			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c1]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])
			cur_pt = cur_pt + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])

		elif c2 == NP.mod(2*n+5,8):
			p_type = 13
			g_type = 4

			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			cur_pt = cur_pt + d[c1]
			poly_pts.append(cur_pt + c[NP.mod(2+n,4)])
			poly_pts.append(cur_pt + c[NP.mod(3+n,4)])
			poly_pts.append(cur_pt + c[NP.mod(0+n,4)])
			cur_pt = cur_pt + d[c2]
			poly_pts.append(cur_pt + c[NP.mod(0+n,4)])

	return cur_pt, p_type, g_type, n, poly_pts

def poly_pts_refine(poly_pts):
	'''
	Description: This function refines the polygon points. Directly using the points
	derived from the 3pv method will result in intersections. This function goes through
	the edge vectors in poly_pts, and finds adjacent edge vectors that form an angle larger
	than 90 degrees, and hence a negative dot product. It removes the point common to these
	two vectors.
	'''
	# Delete repeated points
	for i, pt in reversed(list(enumerate(poly_pts))):
		if i == len(poly_pts)-1:
			continue
		elif NP.all(pt == poly_pts[i+1]):
			del poly_pts[i+1]

	np_poly_pts = NP.array(poly_pts)
	num_pts = len(poly_pts)

	poly_direc = np_poly_pts[1:] - np_poly_pts[:-1]

	poly_dotprod = [NP.dot(poly_direc[n+1],poly_direc[n]) for n in range(num_pts-2)]

	poly_negdot = list(NP.array(poly_dotprod)<0)
	prevneg = False # Usually negative dot products will come in pairs, we only need to consider one

	for i, negdot in reversed(list(enumerate(poly_negdot))):
		if negdot and not prevneg:
			del poly_pts[i]
			prevneg = True
		else:
			prevneg = False

	return poly_pts

def calculate_chain_code(perim_img):
	'''
	Description: This function calculates the chain code of the binary perimeter image.
	This function begins on the left-most, upper-most pixel and moves clockwise. The numbering
	follows the convention given in the 3pv paper.
	'''
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

	return chain_code, first_y, first_x


def perimeter_3pvm(perim_img):
	
	chain_code, first_y, first_x = calculate_chain_code(perim_img)

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
	cur_pt = NP.array([first_y, first_x])
	poly_pts = []
	delta_l_seq = []

	for c1, c2 in chain_pair:

		cur_pt, p_type, g_type, n, cur_poly_pts = chain_pair_props(cur_pt, c1, c2)
		poly_pts = poly_pts + cur_poly_pts

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

			delta_l_seq.append(delta_l[g_type-1][Delta_D-1])
			D2_prev = D2_cur

	# D2_prev and g_type is from the iteration of for loop
	Delta_D = NP.mod(D1_first - D2_prev + 16, 16) + 1
	perimeter = perimeter + delta_l[g_type_first-1][Delta_D-1]

	if isodd:
		if(NP.mod(chain_code[-1],2) == 0):
			perimeter = perimeter + 1
		else:
			perimeter = perimeter + NP.sqrt(2)

	# Refine the poly_pts
	poly_pts.append(poly_pts[0]) # Add the first point back into the array
	poly_pts = poly_pts_refine(poly_pts)

	# Calculate the perimeter to check if it's the same as the other perimeter
	poly_pts = NP.array(poly_pts) - 0.5 # Minus 0.5 because pixels are centered on integers
	poly_direc = poly_pts[1:] - poly_pts[:-1]
	poly_perimeter = 0

	for vec in poly_direc:
		poly_perimeter = poly_perimeter + NP.linalg.norm(vec)

	return perimeter, poly_perimeter, poly_pts
