#!/usr/bin/env python

#
# Last modified: 21 March 2018
# Author: Darrick Lee <y.l.darrick@gmail.com>, Dhananjay Bhaskar <dbhaskar92@gmail.com>, MoHan Zhang <mohan_z@hotmail.com>
# Extract cell shape features from PIF file
# 

from __future__ import division

import os
import sys
import csv
import time
import math
import shlex
import string
import datetime
import threading
import subprocess
import extractcellfeatures

import numpy as NP
import lxml.etree as ET

import matplotlib
import matplotlib.pyplot as PLT
import matplotlib.patches as patches

from PIL import Image
from scipy import interpolate
from optparse import OptionParser
from scipy import ndimage as NDI
from scipy.special import expit
from matplotlib.patches import Ellipse, Polygon
from matplotlib.path import Path

imgDPI = 200

pifFile = None
xmlFile = None
outFolder = './'
l_height, l_width = -1, -1
xmlTime = -1
testCell = -1
plotZoom = -1
plotfits, comparefits, separate = None, None, None
splineSmooth = 10 # Reference: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.interpolate.splprep.html

## OPTIONS PARSING ############################################################

parser = OptionParser()
parser.add_option("-i", "--input", action="store", type="string", dest="inputfile", help="path to PIF file", metavar="PIF")
parser.add_option("-x", "--xml", action="store", type="string", dest="xmlfile", help="path to XML file", metavar="XML")
parser.add_option("-l", "--height", action="store", type="int", dest="height", help="lattice HEIGHT", metavar="HEIGHT")
parser.add_option("-w", "--width", action="store", type="int", dest="width", help="lattice WIDTH", metavar="WIDTH")
parser.add_option("-p","--plot", action="store_true", dest="plotfits", help="plot fitted geometry", default=False)
parser.add_option("-c","--compare", action="store_true", dest="comparefits", help="compare fitted geometry", default=False)
parser.add_option("-t","--time", action="store", type="int", dest="time", help="time in XML file to compare", metavar="TIME")
parser.add_option("-T","--test", action="store", type="int", dest="testcell", help="runs the test code on specified cell")
parser.add_option("-o","--output", action="store", type="string", dest="outputfolder", help="the folder to store output plots", metavar="OUTPUT")
parser.add_option("-s","--separate", action="store_true", dest="separate", help="separate subfolders for plots", default=False, metavar="SEPARATE")
parser.add_option("-m","--multithread", action="store_true", dest="multithread", help="use multithreading", default=False)
parser.add_option("-z","--zoom", action="store", type="int", dest="zoom", help="units to zoom in each direction", metavar="ZOOM")
parser.add_option("-S","--smooth", action="store", type="float", dest="smooth", help="amount of smoothing for spline", metavar="SMOOTH")
parser.add_option("-C","--curvSpline", action="store_true", dest="curvSpline", default=False, help="create plot for curvature-dependent boundary for spline", metavar="CURVSPLINE")
parser.add_option("-v","--csv", action="store_true", dest="createcsv", default=False, help="create feature list csv")

# Options parsing
(options, args) = parser.parse_args()
if options.inputfile:
    pifFile = options.inputfile
if options.xmlfile:
    xmlFile = options.xmlfile
if options.height:
    l_height = options.height
if options.width:
    l_width = options.width
if options.plotfits:
    plotfits = True
else:
    plotfits = False
if options.comparefits:
    comparefits = True
else:
    comparefits = False
if options.time is not None:
    xmlTime = options.time
if options.testcell:
    testCell = options.testcell
if options.outputfolder:
    outFolder = options.outputfolder
if options.separate:
    separate = True
else:
    separate = False
if options.multithread:
    multithreading = True
else:
    multithreading = False
if options.zoom:
    plotZoom = options.zoom
if options.smooth:
    splineSmooth = options.smooth
if options.curvSpline:
    curvSpline = True
else:
    curvSpline = False
if options.createcsv:
    createcsv = True
else:
    createcsv = False

## FOLDER INITIALIZATION ######################################################

plotXML = xmlFile is not None and xmlTime != -1
pifFileName = ''

if os.path.isfile(pifFile):
    pifFileName_withPath = os.path.splitext(pifFile)[0]
    pifFileName = pifFileName_withPath.split('/')[-1]
else:
    print("Error: PIF file does not exist.\n")
    sys.exit(1)

if not os.path.isdir(outFolder):
    print("Error: Output directory does not exist.\n")
    sys.exit(1)

# Check if separate folders have been made and initialize output folders
if separate:
    circleOut = outFolder + 'CircleFit/'
    ellipseOut = outFolder + 'EllipseFit/'
    splineOut = outFolder + 'SplineFit/'
    splineBdyOut = outFolder + 'SplineBdyFit/'
    polyOut = outFolder + 'PolyFit/'
    boundaryOut = outFolder + 'BoundaryFit/'
    vectorizeOut = outFolder + 'vectorizeFit/'

    if (not os.path.isdir(circleOut) or not os.path.isdir(ellipseOut)
        or not os.path.isdir(splineOut) or not os.path.isdir(polyOut)
        or not os.path.isdir(boundaryOut) or not os.path.isdir(vectorizeOut)
        or not os.path.isdir(splineBdyOut)):
        print("Error: One of the output folders do not exist.\n")
        sys.exit(1)

else:
    circleOut = ellipseOut = splineOut = splineBdyOut = polyOut = boundaryOut = vectorizeOut = outFolder

## PARSE PIF FILE #############################################################

lattice = open(pifFile).read().split('\n')[1:-1]
cellDict = dict()
cellTypeDict = dict()
lattice_data = NP.zeros((l_height, l_width))
lattice_matrix = NP.empty((l_height, l_width), NP.uint32)
lattice_matrix.fill(0xFFFFFFFF)

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
    lattice_data[pix_x, pix_y] = cell_id

    if c_type == 'CellU':
        lattice_matrix[pix_x, pix_y] = 0x8000EE00	# alpha, blue, green, red
    if c_type == 'CellV':
        lattice_matrix[pix_x, pix_y] = 0x800000EE


## HELPER FUNCTIONS ###########################################################

def contains_isolated_cells():

    '''
    Description: This returns true if lattice_data contains more than one connected component
    and false otherwise. This currently uses 1-connectivity.

    Reference: http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.label
    '''

    global lattice_data

    clipped_lattice_data = NP.clip(lattice_data, 0, 1)

    [labeled_img, num_labels] = NDI.measurements.label(clipped_lattice_data)

    if num_labels > 1:
        return True
    else:
        return False

## SCALING ####################################################################

# for MIAPaCa (Pancreatic Cancer) Phase Contrast Microscopy

def conv_distance(x):
    conv_factor = 0.4
    units = 'um'
    return x*conv_factor, units

def conv_time(x):
    conv_factor = 0.2
    units = 'min'
    return x*conv_factor, units

def conv_area(x):
    conv_factor = 0.16
    units = 'um_2'
    return x*conv_factor, units

## SINGLE CELL PLOT ###########################################################

def TestSingleCellPlot(extractor):

    extractor.basic_props(splineSmooth)
    extractor.shape_props()
    extractor.cell_centre_fit()
    
    padding = 10

    fixed_perim = NP.transpose(extractor.perim_img)

    perim_img_ind = NP.where(fixed_perim == 1)

    xlim_min = min(perim_img_ind[1])
    xlim_max = max(perim_img_ind[1])

    ylim_min = min(perim_img_ind[0])
    ylim_max = max(perim_img_ind[0])

    U = extractor.spl_u
    OUT = interpolate.splev(U, extractor.spl_poly)
    BDY_FEATS = extractor.bdy_fvector

    # Create Circle/Ellipse plot
    fig = PLT.figure(1)
    ax = fig.add_subplot(111, aspect='equal')

    # Create circle plot
    c = PLT.Circle((extractor.ccm_fvector[0],
            extractor.ccm_fvector[1]),
            extractor.ccm_fvector[2])
            
    xlim_min = min(xlim_min, extractor.ccm_fvector[0] - extractor.ccm_fvector[2])
    xlim_max = max(xlim_max, extractor.ccm_fvector[0] + extractor.ccm_fvector[2])
    ylim_min = min(ylim_min, extractor.ccm_fvector[1] - extractor.ccm_fvector[2])
    ylim_max = max(ylim_max, extractor.ccm_fvector[1] + extractor.ccm_fvector[2])

    # Create ellipse plot
    e = Ellipse(xy=NP.array([extractor.ellipse_fvector[0], extractor.ellipse_fvector[1]]),
            width = extractor.ellipse_fvector[4],
            height = extractor.ellipse_fvector[3],
            angle = extractor.ellipse_fvector[5]/(2*NP.pi)*360)

    # Compute horizontal width and height of the oriented ellipse
    a = extractor.ellipse_fvector[4]/2.0
    b = extractor.ellipse_fvector[3]/2.0
    alpha = extractor.ellipse_fvector[5]
    x_c = extractor.ellipse_fvector[0]
    y_c = extractor.ellipse_fvector[1]
    x_b = 0
    y_b = 0
    if a > b:
        x_b = abs(a * math.cos(alpha))
        y_b = abs(a * math.cos(alpha))
    else:
        x_b = abs(b * math.sin(alpha))
        y_b = abs(b * math.cos(alpha))
        
    xlim_min = min(xlim_min, x_c - x_b)
    xlim_max = max(xlim_max, x_c + x_b)
    ylim_min = min(ylim_min, y_c - y_b)
    ylim_max = max(ylim_max, y_c + y_b)

    PLT.imshow(fixed_perim, interpolation='nearest', cmap='Greys')
    perimeter = conv_distance(extractor.perim_3pv)[0]
    PLT.xticks([])
    PLT.yticks([])
    PLT.plot(extractor.perim_coord_poly[:,0], extractor.perim_coord_poly[:,1], label='Polygon Fit (Perimeter = %.2f um)'%perimeter, color='g', lw=2)

    ax.add_artist(c)
    c.set_alpha(1)
    c.set_facecolor('none')
    c.set_edgecolor('blue')
    c.set_linewidth(3)
    c.set_label('Circle')

    ax.add_artist(e)
    e.set_alpha(1)
    e.set_facecolor('none')
    e.set_edgecolor('orange')
    e.set_linewidth(3)
    e.set_label('Ellipse')
    
    PLT.plot(0, 0, color='blue', label='Circle Fit (Variance = %.2f)'%extractor.ccm_fvector[5], lw=2)
    PLT.plot(0, 0, color='orange', label='Ellipse Fit (Variance = %.2f)'%extractor.ellipse_fvector[8], lw=2)

    lgd = PLT.legend(bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=1, mode=None, fontsize="small", borderaxespad=0.2, fancybox=True, shadow=True)
    
    ax.set_xlim([xlim_min - padding, xlim_max + padding])
    ax.set_ylim([ylim_min - padding, ylim_max + padding])
    
    PLT.xticks([])
    PLT.yticks([])
    
    PLT.savefig(outFolder + pifFileName + '_Fits.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 400)

    # Create spline plot with boundary color based on magnitude and parity of curvature
    fig, (ax1, ax2) = PLT.subplots(1, 2, gridspec_kw = {'width_ratios':[1,10]})
    knorm = expit(extractor.spl_k/max(abs(extractor.spl_k))*10)
    
    norm = matplotlib.colors.Normalize(vmin=NP.min(extractor.spl_k), vmax=NP.max(extractor.spl_k))
    cb = matplotlib.colorbar.ColorbarBase(ax1, cmap=matplotlib.cm.coolwarm, norm=norm, orientation='vertical')
    cb.set_label('Magnitude of Curvature', labelpad=-100)
    
    pcolor = PLT.cm.coolwarm(knorm)

    for i in range(len(U)):
        ax2.plot(OUT[0][i:i+2], OUT[1][i:i+2], color=pcolor[i], linewidth=2)

    xlim_min = min(perim_img_ind[1])
    xlim_max = max(perim_img_ind[1])
    ylim_min = min(perim_img_ind[0])
    ylim_max = max(perim_img_ind[0])
    
    PLT.xticks([])
    PLT.yticks([])
    ax2.set_aspect(1)
    ax2.set_xlim([xlim_min - padding, xlim_max + padding])
    ax2.set_ylim([ylim_min - padding, ylim_max + padding])

    k_protrusions = BDY_FEATS[2]
    k_indentations = BDY_FEATS[3]
    PLT.gcf().text(0.01, 1.02, "Number of Protrusions: %d"%k_protrusions, fontsize=11)
    PLT.gcf().text(0.01, 0.99, "Number of Indentations: %d"%k_indentations, fontsize=11)
    PLT.savefig(outFolder + pifFileName + '_SplineCurvature.png', bbox_inches='tight', dpi = 400)
    
    # Plot curvature function
    min_idx = NP.argmin(NP.absolute(extractor.spl_k))
    spl_k_shifted = NP.roll(extractor.spl_k, min_idx)
    x_data = range(0, len(extractor.spl_k))
    tick_data = NP.roll(x_data, min_idx) 
    
    fig, (ax1, ax2) = PLT.subplots(1, 2, gridspec_kw = {'width_ratios':[6,4]})
     
    ax1.set_xticks(x_data[0:-1:50])
    ax1.set_xticklabels(tick_data[0:-1:50])
    ax1.plot(x_data, spl_k_shifted)
    ax1.set_xlabel('Spline Parameterization Index')
    ax1.set_ylabel('Curvature')
    
    ax2.boxplot(extractor.spl_k)
    ax2.set_xticks([])
    
    PLT.savefig(outFolder + pifFileName + '_CurvatureFcn.png', bbox_inches='tight', dpi = 400)

    # Create spline plot with a binary boundary color scheme based on parity of curvature
    fig = PLT.figure(4)
    ax = fig.add_subplot(111, aspect='equal')
    knorm = NP.sign(extractor.spl_k)/2 + 0.5
    pcolor = PLT.cm.bwr(knorm)

    PLT.imshow(fixed_perim, interpolation='nearest', cmap='Greys')

    for i in range(len(U)):
        PLT.xticks([])
        PLT.yticks([])
        PLT.plot(OUT[0][i:i+2], OUT[1][i:i+2], color=pcolor[i], linewidth=2)

    ax.set_xlim([xlim_min - padding, xlim_max + padding])
    ax.set_ylim([ylim_min - padding, ylim_max + padding])

    PLT.xticks([])
    PLT.yticks([])
    
    PLT.plot(0, 0, color='blue', label='Negative Curvature', lw=2)
    PLT.plot(0, 0, color='red', label='Positive Curvature', lw=2)

    lgd = PLT.legend(bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=2, mode=None, fontsize="small", borderaxespad=0.2, fancybox=True, shadow=True)
    PLT.savefig(outFolder + pifFileName + '_SplineCurvatureBin.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=400)

    # Create oriented bounding rectangle plot
    fig = PLT.figure(5)
    ax = fig.add_subplot(111, aspect='equal')
    PLT.imshow(fixed_perim, interpolation='nearest', cmap='Greys')
    verts = extractor.ret_pvector
    codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]
    path = Path(verts,codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor='red', lw=2)
    ax.add_patch(patch)
    length = conv_distance(extractor.ferret_max)[0]
    width = conv_distance(extractor.ferret_min)[0]
    PLT.xticks([])
    PLT.yticks([])
    PLT.plot(0, 0, color='red', label='Rectangle Fit (Length = %.2f um, Width = %.2f um)'%(length,width), lw=2)
    
    xlim_rec_min = min(verts[0][0], verts[3][0], verts[1][0], verts[2][0]) - padding
    xlim_rec_max = max(verts[1][0], verts[2][0], verts[0][0], verts[3][0]) + padding 
    ylim_rec_min = min(verts[3][1], verts[2][1], verts[0][1], verts[1][1]) - padding 
    ylim_rec_max = max(verts[0][1], verts[1][1], verts[3][1], verts[2][1]) + padding
    
    ax.set_xlim([xlim_rec_min, xlim_rec_max])
    ax.set_ylim([ylim_rec_min, ylim_rec_max])
    
    lgd = PLT.legend(bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=1, mode=None, fontsize="small", borderaxespad=0.2, fancybox=True, shadow=True)
    PLT.savefig(outFolder + pifFileName + '_RectangularFit.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi=400)


## COMPUTE FEATURES FOR ALL CELLS #############################################

featureDict = dict()
polyPtDict = dict()
splinePtDict = dict()
curvatureDict = dict()
rectPtDict = dict()

for cell_id in cellDict.keys():	
    if testCell != -1:
        if cell_id == testCell:
            extractor = extractcellfeatures.ExtractFeatures(cellDict[cell_id])
            TestSingleCellPlot(extractor)
            sys.exit(0)
    else:

        extractor = extractcellfeatures.ExtractFeatures(cellDict[cell_id])

        # If there is more than one connected component in the cell, ignore this cell
        # since we are not able to perform circle or ellipse fits on it
        if extractor.connectedComp > 1:
            del cellDict[cell_id]
            continue

        if multithreading:
            thread_list = []

            thread_list.append(threading.Thread(target=extractor.basic_props(splineSmooth), args=(), kwargs={}))
            thread_list.append(threading.Thread(target=extractor.shape_props(), args=(), kwargs={}))
            thread_list.append(threading.Thread(target=extractor.cell_centre_fit(), args=(), kwargs={}))
            thread_list.append(threading.Thread(target=extractor.moments(), args=(), kwargs={}))

            for thread in thread_list:
                thread.start()

            for thread in thread_list:
                thread.join()
        else:
            extractor.basic_props(splineSmooth)
            extractor.shape_props()
            extractor.cell_centre_fit()
            extractor.moments()

        featureDict[cell_id] = [extractor.area_cell, extractor.perim_1sqrt2, extractor.equiv_diameter]
        featureDict[cell_id] = featureDict[cell_id] + [extractor.perim_3pv]
        featureDict[cell_id] = featureDict[cell_id] + [extractor.perim_poly, extractor.area_poly]
        featureDict[cell_id] = featureDict[cell_id] + extractor.ellipse_fvector + extractor.ccm_fvector
        featureDict[cell_id] = featureDict[cell_id] + [extractor.perim_eroded,  extractor.area_eroded]
        featureDict[cell_id] = featureDict[cell_id] + extractor.bdy_fvector
        featureDict[cell_id] = featureDict[cell_id] + extractor.shape_fvector
        featureDict[cell_id] = featureDict[cell_id] + extractor.ret_fvector
        featureDict[cell_id] = featureDict[cell_id] + extractor.hu_moments

        polyPtDict[cell_id] = extractor.perim_coord_poly
        rectPtDict[cell_id] = extractor.ret_pvector

        # Generate warning if euler number is less than 1:
        if(extractor.shape_fvector[1] < 1):
            print('Warning: Cell ' + str(int(cell_id)) + ' has Euler number less than 1.')

        # Put spline points and curvature data to be used in plots
        splinePtDict[cell_id] = NP.transpose(interpolate.splev(extractor.spl_u, extractor.spl_poly))
        curvatureDict[cell_id] = extractor.spl_k

## CONSTRUCT FEATINDEXDICT ####################################################

featIndexDict = dict(BASIC=None, ELLIPSE=None, CCM=None, TPV=None, ERODED=None, POLY=None, BDY=None, SHAPE=None, RECT=None, MOMENTS=None)

BASIC_numfeat = 3
TPV_numfeat = 1
POLY_numfeat = 2
ELLIPSE_numfeat = 9
CCM_numfeat = 6
ERODED_numfeat = 2
BDY_numfeat = 6
SHAPE_numfeat = 7
RECT_numfeat = 5
MOMENTS_numfeat = 7

TPV_start = BASIC_numfeat
POLY_start = TPV_start + TPV_numfeat
ELLIPSE_start = POLY_start + POLY_numfeat
CCM_start = ELLIPSE_start + ELLIPSE_numfeat
ERODED_start = CCM_start + CCM_numfeat
BDY_start = ERODED_start + ERODED_numfeat
SHAPE_start = BDY_start + BDY_numfeat
RECT_start = SHAPE_start + SHAPE_numfeat
MOMENTS_start = RECT_start + RECT_numfeat

featIndexDict['BASIC'] = dict(
    area=0,
    perimeter=1,
    equiv_diameter=2,
)

featIndexDict['TPV'] = dict(
    perimeter=TPV_start
)

featIndexDict['POLY'] = dict(
    perimeter=POLY_start,
    area=POLY_start+1
)

featIndexDict['ELLIPSE'] = dict(
    centroid_x=ELLIPSE_start,
    centroid_y=ELLIPSE_start+1,
    eccentricity=ELLIPSE_start+2,
    major_axis_length=ELLIPSE_start+3,
    minor_axis_length=ELLIPSE_start+4,
    orientation=ELLIPSE_start+5,
    area=ELLIPSE_start+6,
    perimeter=ELLIPSE_start+7,
    variance=ELLIPSE_start+8
)

featIndexDict['CCM'] = dict(
    centroid_x=CCM_start,
    centroid_y=CCM_start+1,
    radius=CCM_start+2,
    perimeter=CCM_start+3,
    area=CCM_start+4,
    variance=CCM_start+5
)
    
featIndexDict['ERODED'] = dict(
    perimeter=ERODED_start,
    area=ERODED_start+1
)
    
featIndexDict['BDY'] = dict(
    mean=BDY_start,
    std_dev=BDY_start+1,
    protrusions=BDY_start+2,
    indentations=BDY_start+3,
    global_max=BDY_start+4,
    global_min=BDY_start+5
)
    
featIndexDict['SHAPE'] = dict(
    extent=SHAPE_start,
    euler_number=SHAPE_start+1,
    solidity=SHAPE_start+2,
    compactness=SHAPE_start+3,
    elongation=SHAPE_start+4,
    convexity=SHAPE_start+5,
    circularity=SHAPE_start+6,
)

featIndexDict['RECT'] = dict(
    centroid_x=RECT_start,
    centroid_y=RECT_start+1,
    orientation=RECT_start+2,
    ferret_max=RECT_start+3,
    ferret_min=RECT_start+4
)

featIndexDict['MOMENTS'] = dict(
    hu_moment_one=MOMENTS_start,
    hu_moment_two=MOMENTS_start+1,
    hu_moment_three=MOMENTS_start+2,
    hu_moment_four=MOMENTS_start+3,
    hu_moment_five=MOMENTS_start+4,
    hu_moment_six=MOMENTS_start+5,
    hu_moment_seven=MOMENTS_start+6
)

## PRODUCE LATTICE PLOTS WITH FITS ############################################

numFigs = 0;
if plotfits:

    # Plot original cells
    lattice_matrix = NP.ascontiguousarray(lattice_matrix)
    pilImage = Image.frombuffer('RGBA', (l_width, l_height), lattice_matrix, 'raw', 'RGBA', 0, 1)
    pilImage = pilImage.convert('RGB')

    pilImage.save(boundaryOut + pifFileName + '_Boundary.png')

    # Plot polygonized lattice
    sp = None
    if not contains_isolated_cells():
        dirname = os.path.dirname(os.path.abspath(__file__))
        if os.path.isfile(os.path.join(dirname, 'vectorize.py')):
            cmd = "/usr/bin/python3 vectorize.py --file " + pifFile + " --size " + \
            str(l_width) + "," + str(l_height) + " --output " + vectorizeOut
            args = shlex.split(cmd)
            sp = subprocess.Popen(args)

    # Plot ellipse fits
    numFigs += 1
    fig = PLT.figure(numFigs)
    ax = fig.add_subplot(111, aspect='equal')

    ells = [[Ellipse(xy=NP.array([featureDict[cell_id][featIndexDict['ELLIPSE']['centroid_x']]
        ,featureDict[cell_id][featIndexDict['ELLIPSE']['centroid_y']]]),
        width = featureDict[cell_id][featIndexDict['ELLIPSE']['minor_axis_length']],
        height = featureDict[cell_id][featIndexDict['ELLIPSE']['major_axis_length']],
        angle = featureDict[cell_id][featIndexDict['ELLIPSE']['orientation']]/(2*NP.pi)*360),
        cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

    for el in ells:
        e = el[0]
        ctype = el[1]

        ax.add_artist(e)
        e.set_clip_box(ax.bbox)
        e.set_alpha(0.3)
        if ctype == 'CellU':
            e.set_facecolor([0,1,0])
        elif ctype == 'CellV':
            e.set_facecolor([1,0,0])
        else:
            e.set_facecolor([0,0,1])

    if plotZoom != -1:
        ax.set_xlim([plotZoom,l_height-plotZoom])
        ax.set_ylim([plotZoom,l_width-plotZoom])
    else:
        ax.set_xlim([0,l_height])
        ax.set_ylim([0,l_width])
    
    PLT.xticks(rotation=90)
    PLT.yticks(rotation=90)
    PLT.savefig(ellipseOut + pifFileName + '_EllipseFit.png', bbox_inches='tight', dpi = imgDPI)

    # Plot cell-centre model (disk fit)
    numFigs += 1
    fig = PLT.figure(numFigs)
    ax = fig.add_subplot(111, aspect='equal')

    circles = [[PLT.Circle((featureDict[cell_id][featIndexDict['CCM']['centroid_x']],
        featureDict[cell_id][featIndexDict['CCM']['centroid_y']]),
        featureDict[cell_id][featIndexDict['CCM']['radius']]),
        cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

    for circle in circles:
        c = circle[0]
        ctype = circle[1]

        ax.add_artist(c)
        c.set_alpha(0.3)
        if ctype == 'CellU':
            c.set_facecolor([0,1,0])
        elif ctype == 'CellV':
            c.set_facecolor([1,0,0])
        else:
            c.set_facecolor([0,0,1])

    if plotZoom != -1:
        ax.set_xlim([plotZoom,l_height-plotZoom])
        ax.set_ylim([plotZoom,l_width-plotZoom])
    else:
        ax.set_xlim([0,l_height])
        ax.set_ylim([0,l_width])
        
    PLT.xticks(rotation=90)
    PLT.yticks(rotation=90)
    PLT.savefig(circleOut + pifFileName + '_CircleFit.png', bbox_inches='tight', dpi = imgDPI)

    # Plot 3pv-polygon model
    numFigs += 1
    fig = PLT.figure(numFigs)
    ax = fig.add_subplot(111, aspect='equal')

    polys = [[Polygon(NP.array(polyPtDict[cell_id])),
        cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

    for poly in polys:
        p = poly[0]
        ctype = poly[1]
        ax.add_artist(p)
        p.set_alpha(0.3)
        if ctype == 'CellU':
            p.set_facecolor([0,1,0])
        elif ctype == 'CellV':
            p.set_facecolor([1,0,0])
        else:
            p.set_facecolor([0,0,1])

    if plotZoom != -1:
        ax.set_xlim([plotZoom,l_height-plotZoom])
        ax.set_ylim([plotZoom,l_width-plotZoom])
    else:
        ax.set_xlim([0,l_height])
        ax.set_ylim([0,l_width])
        
    PLT.xticks(rotation=90)
    PLT.yticks(rotation=90)
    PLT.savefig(polyOut + pifFileName + '_PolyFit.png', bbox_inches='tight', dpi = imgDPI)

    # Plot spline model
    numFigs += 1
    fig = PLT.figure(numFigs)
    ax = fig.add_subplot(111, aspect='equal')

    polys = [[Polygon(NP.array(splinePtDict[cell_id])),
        cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

    for poly in polys:
        p = poly[0]
        ctype = poly[1]
        ax.add_artist(p)
        p.set_alpha(0.3)
        if ctype == 'CellU':
            p.set_facecolor([0,1,0])
        elif ctype == 'CellV':
            p.set_facecolor([1,0,0])
        else:
            p.set_facecolor([0,0,1])

    if plotZoom != -1:
        ax.set_xlim([plotZoom,l_height-plotZoom])
        ax.set_ylim([plotZoom,l_width-plotZoom])
    else:
        ax.set_xlim([0,l_height])
        ax.set_ylim([0,l_width])
        
    PLT.xticks(rotation=90)
    PLT.yticks(rotation=90)
    PLT.savefig(splineOut + pifFileName + '_SplineFit.png', bbox_inches='tight', dpi = imgDPI)

    # Plot Rectangular Fits
    numFigs += 1
    fig = PLT.figure(numFigs)
    ax = fig.add_subplot(111, aspect='equal')

    rects = [[rectPtDict[cell_id],
        cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

    for re in rects:
        verts = re[0]
        ctype = re[1]
        codes = [Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         ]
        path = Path(verts,codes)
        if ctype == 'CellU':
            patch = patches.PathPatch(path, facecolor=[0,1,0], edgecolor='green', lw=1, alpha=0.3)
        elif ctype == 'CellV':
            patch = patches.PathPatch(path, facecolor=[1,0,0], edgecolor='green', lw=1, alpha=0.3)
        else:
            patch = patches.PathPatch(path, facecolor=[0,0,1], edgecolor='green', lw=1, alpha=0.3)
        ax.add_patch(patch)
        PLT.plot(0,0,color='green',label='Rectangle',lw=1)

    if plotZoom != -1:
        ax.set_xlim([plotZoom,l_height-plotZoom])
        ax.set_ylim([plotZoom,l_width-plotZoom])
    else:
        ax.set_xlim([0,l_height])
        ax.set_ylim([0,l_width])
        
    PLT.xticks(rotation=90)
    PLT.yticks(rotation=90)
    PLT.savefig(splineOut + pifFileName + '_RectangularFit.png', bbox_inches='tight', dpi = imgDPI)

    # Plot spline model with curvature dependent boundary
    if curvSpline:
        numFigs += 1
        fig = PLT.figure(numFigs)
        ax = fig.add_subplot(111, aspect='equal')

        polys = [[Polygon(NP.array(splinePtDict[cell_id])),
            cellTypeDict[cell_id]] for cell_id in cellDict.keys()]

        for poly in polys:
            p = poly[0]
            ctype = poly[1]
            ax.add_artist(p)
            p.set_alpha(0.3)
            p.set_edgecolor('none')
            if ctype == 'CellU':
                p.set_facecolor([0,1,0])
            elif ctype == 'CellV':
                p.set_facecolor([1,0,0])
            else:
                p.set_facecolor([0,0,1])

        for cell_id in cellDict.keys():
            k = curvatureDict[cell_id]
            OUT = NP.transpose(splinePtDict[cell_id])

            knorm = NP.sign(k)/2 + 0.5
            pcolor = PLT.cm.bwr(knorm)

            for i in range(len(k)):
                PLT.plot(OUT[0][i:i+2],OUT[1][i:i+2],color=pcolor[i],linewidth=0.75)

        if plotZoom != -1:
            ax.set_xlim([plotZoom,l_height-plotZoom])
            ax.set_ylim([plotZoom,l_width-plotZoom])
        else:
            ax.set_xlim([0,l_height])
            ax.set_ylim([0,l_width])
            
        PLT.xticks(rotation=90)
        PLT.yticks(rotation=90)
        PLT.savefig(splineBdyOut + pifFileName + '_SplineBdyFit.png', bbox_inches='tight', dpi = imgDPI)

    if sp is not None:
        sp.wait()

## PRODUCE COMPARISON PLOTS ###################################################
 
if comparefits:

    # Initialize data vectors
    xData = []

    perim_1sqrt2 = []
    perim_ellipse = []
    perim_circle = []
    perim_eroded = []
    perim_3pv = []
    perim_poly = []
    perim_xml = []

    area_ellipse = []
    area_circle = []
    area_basic = []
    area_eroded = []
    area_poly = []
    area_xml = []

    # Prepare XML file
    if plotXML:
        infile = open(xmlFile,'r')
        xml1 = ET.parse(infile)
        root = xml1.getroot()
        times = root.getchildren()
        cells = times[xmlTime].getchildren()

    for cell_id in cellDict.keys():
        xData.append(cell_id)

    # Populate the data vectors
    for cell_id in xData:
    
        perim_1sqrt2.append(featureDict[cell_id][featIndexDict['BASIC']['perimeter']])
        perim_ellipse.append(featureDict[cell_id][featIndexDict['ELLIPSE']['perimeter']])
        perim_circle.append(featureDict[cell_id][featIndexDict['CCM']['perimeter']])
        perim_eroded.append(featureDict[cell_id][featIndexDict['ERODED']['perimeter']])
        perim_3pv.append(featureDict[cell_id][featIndexDict['TPV']['perimeter']])
        perim_poly.append(featureDict[cell_id][featIndexDict['POLY']['perimeter']])

        area_ellipse.append(featureDict[cell_id][featIndexDict['ELLIPSE']['area']])
        area_circle.append(featureDict[cell_id][featIndexDict['CCM']['area']])
        area_eroded.append(featureDict[cell_id][featIndexDict['ERODED']['area']])
        area_basic.append(featureDict[cell_id][featIndexDict['BASIC']['area']])
        area_poly.append(featureDict[cell_id][featIndexDict['POLY']['area']])

        if plotXML:
            perim_xml.append(cells[cell_id].get('perimeter'))
            area_xml.append(cells[cell_id].get('area'))

    # Reorder the vectors so that xml parameters are sorted from smallest to largest
    # Note: We switch to NP.array rather than list so we can input a list to select elements
    perim_reorder_ind = NP.argsort(perim_1sqrt2)
    perim_1sqrt2, perim_unit = conv_distance(NP.array(perim_1sqrt2)[perim_reorder_ind])
    perim_ellipse, _ = conv_distance(NP.array(perim_ellipse)[perim_reorder_ind])
    perim_circle, _ = conv_distance(NP.array(perim_circle)[perim_reorder_ind])
    perim_poly, _ = conv_distance(NP.array(perim_poly)[perim_reorder_ind])
    perim_eroded, _ = conv_distance(NP.array(perim_eroded)[perim_reorder_ind])
    perim_3pv, _ = conv_distance(NP.array(perim_3pv)[perim_reorder_ind])

    area_reorder_ind = NP.argsort(area_basic)
    area_ellipse, area_unit = conv_area(NP.array(area_ellipse)[area_reorder_ind])
    area_circle, _ = conv_area(NP.array(area_circle)[area_reorder_ind])
    area_poly, _ = conv_area(NP.array(area_poly)[area_reorder_ind])
    area_eroded, _ = conv_area(NP.array(area_eroded)[area_reorder_ind])
    area_basic, _ = conv_area(NP.array(area_circle)[area_reorder_ind])

    if plotXML:
        perim_xml = NP.array(perim_xml)[perim_reorder_ind]
        area_xml = NP.array(area_xml)[area_reorder_ind]

    cell_range = range(len(xData))

    # Plot the figures
    numFigs += 1
    PLT.figure(numFigs)
    if plotXML:
        PLT.plot(cell_range, perim_xml)
    PLT.plot(cell_range, perim_ellipse)
    PLT.plot(cell_range, perim_circle)
    PLT.plot(cell_range, perim_poly)
    PLT.plot(cell_range, perim_eroded)
    PLT.plot(cell_range, perim_1sqrt2)
    PLT.plot(cell_range, perim_3pv)

    PLT.xlabel('Cell (arbitrary)')
    PLT.ylabel('Perimeter (' + perim_unit + ')')
    PLT.title('Comparison of Perimeter vs. Cell')
    lgd = PLT.legend(["Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "1, sqrt(2) method", "3pv method"],
    bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)

    if plotXML:
        lgd = PLT.legend(["xml", "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "1, sqrt(2) method", "3pv method"],
        bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)

    PLT.savefig(outFolder + pifFileName + '_PerimCompare.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = imgDPI)

    numFigs += 1
    PLT.figure(numFigs)
    if plotXML:
        PLT.plot(cell_range, area_xml)
    PLT.plot(cell_range, area_ellipse)
    PLT.plot(cell_range, area_circle)
    PLT.plot(cell_range, area_poly)
    PLT.plot(cell_range, area_eroded)
    PLT.plot(cell_range, area_basic)

    PLT.xlabel('Cell (arbitrary)')
    PLT.ylabel('Area (' + area_unit + ')')
    PLT.title('Comparison of Area vs. Cell')
    lgd = PLT.legend([ "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "Pixel counting"],
    bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)

    if plotXML:
        PLT.legend([ "xml", "Ellipse fit", "CCM fit", "Poly fit", "Eroded Poly", "Pixel counting"],
        bbox_to_anchor=(0.0, 1.1, 1.0, 1.5), loc=3, ncol=3, mode="expand", borderaxespad=0.2, fancybox=True, shadow=True)

    PLT.savefig(outFolder + pifFileName + '_AreaCompare.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = imgDPI)

 ## CREATE CSV FILE ###########################################################

if createcsv:
    ## BUILD FEATURE LIST #########################################################

    # Build location list - features not used for clustering
    # Build feature list
    locationList = []
    nondimList = []
    momentList = []
    featListOrig = []

    for cell_id in cellDict.keys():
    
        # location features
        cellLoc = []
        cellLoc.append(cell_id)
        cellLoc.append(featureDict[cell_id][featIndexDict['ELLIPSE']['centroid_x']])
        cellLoc.append(featureDict[cell_id][featIndexDict['ELLIPSE']['centroid_y']])
        cellLoc.append(featureDict[cell_id][featIndexDict['CCM']['centroid_x']])
        cellLoc.append(featureDict[cell_id][featIndexDict['CCM']['centroid_y']])
        cellLoc.append(featureDict[cell_id][featIndexDict['RECT']['centroid_x']])
        cellLoc.append(featureDict[cell_id][featIndexDict['RECT']['centroid_y']])
        locationList.append(list(cellLoc))

        # non-dimensional features
        cellNondim = []
        cellNondim.append(featureDict[cell_id][featIndexDict['SHAPE']['extent']])
        cellNondim.append(featureDict[cell_id][featIndexDict['SHAPE']['euler_number']])
        cellNondim.append(featureDict[cell_id][featIndexDict['SHAPE']['solidity']])
        cellNondim.append(featureDict[cell_id][featIndexDict['SHAPE']['compactness']])
        cellNondim.append(featureDict[cell_id][featIndexDict['SHAPE']['elongation']])
        cellNondim.append(featureDict[cell_id][featIndexDict['SHAPE']['convexity']])
        cellNondim.append(featureDict[cell_id][featIndexDict['SHAPE']['circularity']])
        nondimList.append(list(cellNondim))
        
        cellHu = []
        cellHu.append(featureDict[cell_id][featIndexDict['MOMENTS']['hu_moment_one']])
        cellHu.append(featureDict[cell_id][featIndexDict['MOMENTS']['hu_moment_two']])
        cellHu.append(featureDict[cell_id][featIndexDict['MOMENTS']['hu_moment_three']])
        cellHu.append(featureDict[cell_id][featIndexDict['MOMENTS']['hu_moment_four']])
        cellHu.append(featureDict[cell_id][featIndexDict['MOMENTS']['hu_moment_five']])
        cellHu.append(featureDict[cell_id][featIndexDict['MOMENTS']['hu_moment_six']])
        cellHu.append(featureDict[cell_id][featIndexDict['MOMENTS']['hu_moment_seven']])
        momentList.append(list(cellHu))

        # Create feature vector for the cell, but remove centroids as they are in the location list
        # Also remove nondimensional shape factors and hu moments
        
        cellFeat = list(featureDict[cell_id])
        
        cellFeat.pop(featIndexDict['MOMENTS']['hu_moment_seven'])
        cellFeat.pop(featIndexDict['MOMENTS']['hu_moment_six'])
        cellFeat.pop(featIndexDict['MOMENTS']['hu_moment_five'])
        cellFeat.pop(featIndexDict['MOMENTS']['hu_moment_four'])
        cellFeat.pop(featIndexDict['MOMENTS']['hu_moment_three'])
        cellFeat.pop(featIndexDict['MOMENTS']['hu_moment_two'])
        cellFeat.pop(featIndexDict['MOMENTS']['hu_moment_one'])
        
        cellFeat.pop(featIndexDict['RECT']['centroid_y'])
        cellFeat.pop(featIndexDict['RECT']['centroid_x'])
        
        cellFeat.pop(featIndexDict['SHAPE']['circularity'])
        cellFeat.pop(featIndexDict['SHAPE']['convexity'])
        cellFeat.pop(featIndexDict['SHAPE']['elongation'])
        cellFeat.pop(featIndexDict['SHAPE']['compactness'])
        cellFeat.pop(featIndexDict['SHAPE']['solidity'])
        cellFeat.pop(featIndexDict['SHAPE']['euler_number'])
        cellFeat.pop(featIndexDict['SHAPE']['extent'])
        
        cellFeat.pop(featIndexDict['CCM']['centroid_y'])
        cellFeat.pop(featIndexDict['CCM']['centroid_x'])
        cellFeat.pop(featIndexDict['ELLIPSE']['centroid_y'])
        cellFeat.pop(featIndexDict['ELLIPSE']['centroid_x'])

        featListOrig.append(cellFeat)

    # Create location name list
    locationNames = ['Cell_ID', 'ELLIPSE_centroid_x', 'ELLIPSE_centroid_y', 'CCM_centroid_x', 'CCM_centroid_y', 'RECT_centroid_x', 'RECT_centroid_y']

    # Create nondimensional factors name list
    nondimNames = ['SHAPE_extent', 'SHAPE_euler_number', 'SHAPE_solidity', 'SHAPE_compactness', 'SHAPE_elongation', 'SHAPE_convexity', 'SHAPE_circularity']
    
    # Create hu moment name list
    huNames = ['Hu_moment_1', 'Hu_moment_2', 'Hu_moment_3', 'Hu_moment_4', 'Hu_moment_5', 'Hu_moment_6','Hu_moment_7']

    # Create feature name list
    featNamesOrig = []
    featNamesOrig = featNamesOrig + ['BASIC_area', 'BASIC_perimeter', 'BASIC_equiv_diameter']
    featNamesOrig = featNamesOrig + ['TPV_perimeter']
    featNamesOrig = featNamesOrig + ['POLY_perimeter', 'POLY_area']
    featNamesOrig = featNamesOrig + ['ELLIPSE_eccentricity', 'ELLIPSE_major_axis_length', 'ELLIPSE_minor_axis_length', 'ELLIPSE_orientation', 'ELLIPSE_area', 'ELLIPSE_perimeter', 'ELLIPSE_variance']
    featNamesOrig = featNamesOrig + ['CCM_radius', 'CCM_perimeter', 'CCM_area', 'CCM_variance']
    featNamesOrig = featNamesOrig + ['ERODED_perimeter', 'ERODED_area']
    featNamesOrig = featNamesOrig + ['BDY_mean', 'BDY_std_dev', 'BDY_protrusions', 'BDY_indentations', 'BDY_global_max', 'BDY_global_min']
    featNamesOrig = featNamesOrig + ['RECT_orientation', 'RECT_major_axis_length', 'RECT_minor_axis_length']

    numFeat = len(featNamesOrig)    

    # Transpose lists and create standardized lists
    featListOrig = list(NP.transpose(featListOrig))
    locationList = list(NP.transpose(locationList))
    nondimList = list(NP.transpose(nondimList))
    huList = list(NP.transpose(momentList))
    featListStd = list(featListOrig)
    featNamesStd = [fname + '_std' for fname in featNamesOrig]

    # Select the index of elements to scale
    dist_index = [1,2,3,4,7,8,11,13,14,17,26,27]
    loc_index = [1,2,3,4,5,6]
    area_index = [0,5,10,15,18]

    for i in dist_index:
        featListOrig[i][:], dist_unit = conv_distance(featListOrig[i][:])
        featNamesOrig[i] = featNamesOrig[i] + '_' + dist_unit

    for i in loc_index:
        locationList[i][:], dist_unit = conv_distance(locationList[i][:])
        locationNames[i] = locationNames[i] + '_' + dist_unit

    for i in area_index:
        featListOrig[i][:], area_unit = conv_area(featListOrig[i][:])
        featNamesOrig[i] = featNamesOrig[i] + '_' + area_unit

    ## STANDARDIZE FEATURES #######################################################

    featListStd = NP.array(featListStd)
    # Don't standardize the centroids
    for k in range(numFeat):
        featListStd[k] = (featListStd[k] - NP.mean(featListStd[k]))/NP.sqrt(NP.var(featListStd[k]))
    featListStd = list(featListStd)

    ## CREATE CSV FILE ############################################################

    # Put all features into a single list
    featList = NP.transpose(locationList + nondimList + huList + featListStd + featListOrig)
    featNamesOrig = locationNames + nondimNames + huNames + featNamesStd + featNamesOrig
    numCells = len(cellDict.keys())
    outFile = pifFileName + '.csv'

    with open(outFile, 'wb') as csvfile:
        featWriter = csv.writer(csvfile)

        # Write the comments
        ts = time.time()
        timestamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

        featWriter.writerow(['# Simulation Data Feature List'])
        featWriter.writerow(['# Date Created: ' + timestamp])
        featWriter.writerow(['# Original PIF file: ' + pifFile])

        featWriter.writerow(['# Identifier (Cell ID) Column: 0'])
        featWriter.writerow(['# Location (Ellipse Centroid) Columns: 1-2'])
        featWriter.writerow(['# Location (Circle Centroid) Columns: 3-4'])
        featWriter.writerow(['# Location (Rect Centroid) Columns: 5-6'])

        featWriter.writerow(['# Nondimensional (Shape Factor) Columns: 7-13'])
        featWriter.writerow(['# Hu Moments Columns: 14-20'])
        featWriter.writerow(['# Standardized Feature Columns: 21-' + str(20+numFeat)])
        featWriter.writerow(['# Original Feature Columns: ' + str(21+numFeat) + '-' + str(20+2*numFeat)])
        featWriter.writerow(['# Number of Cells: ' + str(numCells)])

        featWriter.writerow([])
        featWriter.writerow(featNamesOrig)

        for row in featList:
            featWriter.writerow(row)
