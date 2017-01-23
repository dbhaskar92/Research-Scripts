#
# Last modified: 19 May 2016
# Authors: Darrick Lee <y.l.darrick@gmail.com>, Dhananjay Bhaskar <dbhaskar92@gmail.com>
# Description: Compute features from list of pixels representing a Cellular Potts Model (CPM) cell
#

import numpy as NP
import scipy.special
import skimage.measure
import skimage.morphology
import perimeter_3pvm
import matplotlib.pyplot as PLT

from operator import itemgetter
from scipy import interpolate
from scipy import ndimage as NDI
from scipy.signal import argrelextrema

class ExtractFeatures:

	def __init__(self, cell_pixel_list):

		self.pix_list = cell_pixel_list

		self.cell_img = None # Binary image of full cell
		self.perim_img = None # Binary image of cell perimeter
		self.eroded_img = None # Binary image of eroded cell perimeter

		self.perim_coord = None # Coordinates of the cell perimeter pixels
		self.perim_coord_dp = None # Coordinates of approximate perimeter (Douglas-Peucker)
		self.perim_coord_poly = None # Coordinates of polygon (derived from 3pv)
		self.perim_coord_eroded = None # Coordinates of eroded polygon

		self.area_cell = None # Area of cell by pixel counting
		self.area_poly = None # Area of polygon derived from 3pv
		self.area_eroded = None # Area of eroded polygon

		self.perim_3pv = None # Perimeter of cell using 3pv method
		self.perim_1sqrt2 = None # 1 sqrt 2 method
		self.perim_poly = None # Perimeter of polygon derived from 3pv
		self.perim_eroded = None # Perimeter of polygon of eroded cell

		self.equiv_diameter = None # The equivalent diameter of a circle with same area as cell
		self.shape_factor = None

		self.ellipse_fvector = None # Features from ellipse fit
		self.shape_fvector = None # Shape factors
		self.ccm_fvector = None # Features from CCM fit
		self.bdy_fvector = None # Features based on boundary

		self.spl_poly = None # Spline tck variables approximating 3pv-polygon
		self.spl_u = None # Spline parameter
		self.spl_k = None # Spline curvature

		self.cell_to_image()

		# Check how many connected components there are
		s = [[1,1,1],[1,1,1],[1,1,1]] # Used to allow diagonal connections
		_, self.connectedComp = NDI.measurements.label(self.cell_img,structure=s)

	def cell_to_image(self):

		# Find x, y coordinate bounds
		x_res = max(self.pix_list, key=itemgetter(0))[0]
		y_res = max(self.pix_list, key=itemgetter(1))[1]

		# Creating labeled_img
		self.cell_img = NP.zeros([x_res+2, y_res+2], dtype=NP.int_)

		for (x_pix, y_pix) in self.pix_list:
			self.cell_img[x_pix-1, y_pix-1] = 1

		# Find the pixels that make up the perimeter
		eroded_image = NDI.binary_erosion(self.cell_img)

		eroded_image_open = NDI.binary_opening(eroded_image, structure=NP.ones((3,3)))
		eroded_image_open2 = NDI.binary_erosion(eroded_image_open)

		# self.perim_img = self.cell_img - eroded_image
		self.eroded_img = eroded_image_open - eroded_image_open2
		self.perim_img = self.cell_img - eroded_image

		# Create a list of the coordinates of the pixels (use the center of the pixels)
		perim_image_ind = NP.where(self.perim_img == 1)
		perim_image_coord = NP.array([perim_image_ind[0], perim_image_ind[1]])
		self.perim_coord = NP.transpose(perim_image_coord)

		return


	def basic_props(self, splineSmooth=10):

		'''
		Description: Calculates the perimeter and area using basic methods. For perimeter,
		we use the 3pv, 1 sqrt2 method, and look at the 3pv-polygon perimeter. For area,
		we use pixel counting, and look at the 3pv-polygon area.

		For 3pv perimeter: Use three-pixel vector method to compute perimeter and shape factor
		Reference: http://www.sciencedirect.com/science/article/pii/0308912687902458

		For cubic spline: Use the built-in function from scipy. Note about smoothing parameter in reference.
		Reference: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html
		'''

		# Perimeter: 3pv and polygon perimeter (polygon from 3pv)
		self.perim_3pv, self.perim_poly, self.perim_coord_poly = perimeter_3pvm.perimeter_3pvm(self.perim_img)
		_, self.perim_eroded, self.perim_coord_eroded = perimeter_3pvm.perimeter_3pvm(self.eroded_img)

		# Perimeter: Approximate polygon using Douglas-Peucker algorithm
		self.perim_coord_dp = skimage.measure.approximate_polygon(NP.array(self.perim_coord_poly), 0.75)

		# Create cubic spline
		self.spl_poly, _ = interpolate.splprep(NP.transpose(self.perim_coord_poly), per=1,s=splineSmooth)
		self.spl_u = NP.linspace(0,1.0,100)

		# Calculate spline curvature
		D1 = interpolate.splev(self.spl_u, self.spl_poly, der=1)
		D2 = interpolate.splev(self.spl_u, self.spl_poly, der=2)
		self.spl_k = (D1[0]*D2[1] - D1[1]*D2[0])/((D1[0]**2+D1[1]**2)**(3./2))

		# Calculate boundary features
		ksign = NP.sign(self.spl_k)
		signchange = ((NP.roll(ksign, 1) - ksign) != 0).astype(int)
		num_signchange = sum(signchange)

		localmax = argrelextrema(self.spl_k, NP.greater)[0]
		localmin = argrelextrema(self.spl_k, NP.less)[0]
		num_extrema = len(localmax) + len(localmin)

		posk = [n for n in self.spl_k if n > 0]
		negk = [n for n in self.spl_k if n < 0]

		# Calculate the max/min positive and negative curvatures
		if posk:
			maxpos = max([n for n in self.spl_k if n > 0])
			minpos = min([n for n in self.spl_k if n > 0])
		else:
			maxpos = 0
			minpos = 0

		if negk:
			maxneg = max([n for n in self.spl_k if n < 0])
			minneg = min([n for n in self.spl_k if n < 0])
		else:
			maxneg = 0
			minneg = 0

		self.bdy_fvector = [maxpos, minpos, maxneg, minneg, num_extrema, num_signchange]

		# Perimeter: 1 sqrt2 (function from regionprops)
		props = skimage.measure.regionprops(self.cell_img)
		self.perim_1sqrt2 = props[0].perimeter

		# Area: Pixel Counting
		self.area_cell = len(self.pix_list)

		# Area: Polygon area (from 3pv polygon)
		# Extract x and y coordinates
		# We subtract 0.5 because PLT.imshow() shows coordinates as the centers of pixels
		# Using the shoelace formula: https://en.wikipedia.org/wiki/Shoelace_formula
		YY = self.perim_coord_poly[:,0]
		XX = self.perim_coord_poly[:,1]

		self.area_poly = 0.5*NP.abs(NP.dot(XX,NP.roll(YY,1))-NP.dot(YY,NP.roll(XX,1)))

		# Area: Eroded polygon area
		YY = self.perim_coord_eroded[:,0]
		XX = self.perim_coord_eroded[:,1]

		self.area_eroded = 0.5*NP.abs(NP.dot(XX,NP.roll(YY,1))-NP.dot(YY,NP.roll(XX,1)))

		# Equivalent Diameter: The diameter of the circle with the same area (pixel counted) as the cell
		self.equiv_diameter = NP.sqrt(4*self.area_cell/NP.pi)

		return


	def shape_props(self):

		'''
		Description: Returns list of properties derived from fitting ellipse (in the following order)
		centroid_x, centroid_y, eccentricity, majorAxisLength, minorAxisLength, orientation,
		area and perimeter.

		This uses regionprops() fom skimage.measure. The ellipse fit is done by
		fitting an ellipse with the same second central moment as the image. By looking
		at the code, this is done by calculating the inertia tensor of the matrix,
		finding the eigenvalues (the second central moments using the principal axes),
		and matching those with the equations for second central moment of an ellipse.

		In addition, this method returns a set of shape factors such as extent, euler number,
		solidity, compactness, elongation, convexity, and circularity.

		NOTE: We assume that basic_props() has already been called.

		Reference: https://en.wikipedia.org/wiki/Image_moment

		Ellipse perimeter: Equation given in https://en.wikipedia.org/wiki/Ellipse#Circumference
		The elliptic integral of the second kind implemented in scipy:
		http://docs.scipy.org/doc/scipy/reference/generated/scipy.special.ellipe.html#scipy.special.ellipe
		Note that the scipy definition of the integral differs slightly than wiki, so we take E(e^2) rather than E(e).
		'''

		props = skimage.measure.regionprops(self.cell_img)

		centroid = props[0].centroid

		# Calculate ellipse variance
		perim_coord = NP.transpose(self.perim_coord)
		numPt = len(perim_coord[0])
		cov = NP.mat(NP.cov(perim_coord))
		V = NP.array([perim_coord[0] - centroid[0], perim_coord[1] - centroid[1]])
		cV = NP.array(NP.linalg.inv(cov)*NP.mat(V))
		d = NP.sqrt(V[0]*cV[0] + V[1]*cV[1])
		mu = NP.sum(d)/numPt
		sigma = NP.sqrt(NP.sum((d-mu)**2)/numPt)

		ellipse_prop_list = [centroid[0]]
		ellipse_prop_list.append(centroid[1])
		ellipse_prop_list.append(props[0].eccentricity)
		ellipse_prop_list.append(props[0].major_axis_length)
		ellipse_prop_list.append(props[0].minor_axis_length)
		ellipse_prop_list.append(props[0].orientation) # In degrees starting from the x-axis
		ellipse_prop_list.append(NP.pi*ellipse_prop_list[4]*ellipse_prop_list[5]/4.0) # Ellipse area
		ellipse_prop_list.append(2.0*ellipse_prop_list[4]*scipy.special.ellipe(props[0].eccentricity**2)) # Ellipse perimeter, parameter to special.ellipe should be eccentricity**2
		ellipse_prop_list.append(sigma/mu) # Ellipse variance

		self.ellipse_fvector = ellipse_prop_list

		# NOTE: For shape factors, we use the perim_poly for perimeter, and pixel counting for area

		# Calculate values needed for shape factors
		inertia_ev = props[0].inertia_tensor_eigvals
		area = self.area_cell
		perim = self.perim_poly

		# Calculate convex hull perimeter
		cvx_img = props[0].convex_image # Find the pixels that make up the perimeter
		eroded_cvx_img = NDI.binary_erosion(cvx_img)
		cvx_perim_img = cvx_img - eroded_cvx_img
		cvx_perim_img = NP.lib.pad(cvx_perim_img,(1,1),'constant') # Pad with 0's for perimeter code to work properly
		_, cvx_perim, _ = perimeter_3pvm.perimeter_3pvm(cvx_perim_img)


		# Calculate shape factors
		# compactness = area**2/(NP.pi*2*NP.sqrt(inertia_ev[0]**2 + inertia_ev[1]**2))
		compactness = NP.sqrt((4*area)/NP.pi)/(props[0].major_axis_length) # sqrt(4(area)/pi)/(MaximumDiameter of Ellipse Fit)
		elongation = 1-(props[0].minor_axis_length/props[0].major_axis_length)
		convexity = cvx_perim/perim
		circularity = 4*NP.pi*area/(perim**2)
		extension = area/(NP.pi*ellipse_prop_list[4]*ellipse_prop_list[5]/4.0) #Extension = Area of cell/area of ellipse fit?
		dispersion = 0 #set to 0 for now

		# Create shape feature vector
		self.shape_fvector = []
		self.shape_fvector.append(props[0].extent) # Ratio of pixels in the region to pixel in bounding box (from 0 to 1)
		self.shape_fvector.append(props[0].euler_number) # Euler number
		self.shape_fvector.append(props[0].solidity) # Ratio of pixels in the region to pixels of the convex hull image (from 0 to 1)
		self.shape_fvector.append(compactness) # Ratio of squared area to magnitude of second moments (circle = 1, I-shape << 1)
		self.shape_fvector.append(elongation) # Square root of ratio of two second moments (smaller over larger) (from 0 to 1)
		self.shape_fvector.append(convexity) # Ratio of convex hull perimeter to perimeter (from 0 to 1)
		self.shape_fvector.append(circularity) # Ratio of area to perimeter squared (circle = 1, starfish << 1)
		self.shape_fvector.append(extension)
		self.shape_fvector.append(dispersion)

		return


	def cell_centre_fit(self):

		'''
		Description: Returns a list of features derived from fitting a circle (in the following order):
		centroid_x, centroid_y, radius, perimeter, area.

		This uses a least-squares estimator for the circle, using the points on the boundary of the cell.
		These points are chosen to be at the center of the boundary pixels.

		Circle variance is a goodness of fit measure for the circle fit and is defined in this reference:
		http://www.math.uci.edu/icamp/summer/research_11/park/shape_descriptors_survey.pdf
		'''

		c_model = skimage.measure.CircleModel()
		c_model.estimate(self.perim_coord)

		if skimage.__version__ == '0.9.3':
			(xc, yc, r) = c_model._params
		else:									# For newer versions
			(xc, yc, r) = c_model.params

		# Calculate the circle variance
		perim_coord = NP.transpose(self.perim_coord)
		numPt = len(perim_coord[0])
		d = NP.sqrt((perim_coord[0]-xc)**2 + (perim_coord[1]-yc)**2)
		mu = NP.sum(d)/numPt
		sigma = NP.sqrt(NP.sum((d-mu)**2)/numPt)

		cell_centre_features = [xc]
		cell_centre_features.append(yc)
		cell_centre_features.append(r)
		cell_centre_features.append(2*NP.pi*r)
		cell_centre_features.append(NP.pi*r**2)
		cell_centre_features.append(sigma/mu) # Circle variance (lower is better)

		self.ccm_fvector = cell_centre_features
		return
