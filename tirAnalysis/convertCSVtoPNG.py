#!/usr/bin/python 

"""
convertCSVtoGreyscalePNG.py : Reads a floating-point csv file and returns an 
    integer-valued, greyscale png.  
    For use with thermal images, the maximum and minimum values (23 and 80 degC 
    by default) need to be given.
"""

from __future__ import print_function
from numpy import genfromtxt

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.misc

fileIn = sys.argv[1]  # input file name

# Read the file to a numpy array 
data = genfromtxt(fileIn, delimiter=';')

# rescale data 
range = np.array([23.0, 80.0])  # bounds of tempeature range
rescaled = (data - range.min()) / range.max() * 255.0

# clip out-of-bounds values to 0-255
for i, col in enumerate(rescaled):
    for j, elem in enumerate(col):
        if elem < 0:
            rescaled[i][j] =   0.0
        if elem > 255:
            rescaled[i][j] = 255.0

fileOut = fileIn[::-1][4:][::-1] + '.png'  # reverse the string, cuts the last 4
                                           # chars, reverse this + add '.png'

cmap = plt.cm.CMRmap

# Save to file, conserving the max and min colourmap range
plt.imsave(fileOut, rescaled, cmap=cmap, vmin=0.0, vmax=255.0)
#scipy.misc.toimage(rescaled.round(), cmin=0.0, cmax=255.0).save(fileOut)

print("For file %s, data range is %5.1f, %5.1f --> %s " % (fileIn, 
                                                           data.min(), 
                                                           data.max(), 
                                                           fileOut))
