#!/usr/bin/python3

"""
convertCSVtoGreyscalePNG.py
"""

#from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import sys
#import scipy.misc

def convertCSVtoGreyscalePNG(fileIn, Trange=None, delimiter=';', postfix=''):
    """
    Reads a floating-point csv file and returns an integer-valued, 
    greyscale (i.e. 1 "colour" channel + alpha layer) png.  
    Typically for use with thermal images.

    Parameters:
    -----------
    fileIn : string
        Absolute or relative path of the file to be converted
    Trange : list, tuple or 1d array
        Range of temperatures that the images should be clipped to.  This
        is useful when a sequence of images are to be converted so as to keep
        the temperature scale between them equal.
        Default is for Trange to be of type None, in which case the upper and 
        lower limits of the data range will be used.
    delimiter : string
        Separator or delimiting character between data fields.  Default is
        a semicolon, ';'
    postfix : string
        A (short) sequence of characters that will be appended to the output 
        filename (before the file extension) in order to distinguish the 
        output images from other images.

    TO DO:
    - Introduce tifffile.imsave for saving 16-bit depth images    

    """
    # Read the file to a numpy array 
    data = np.genfromtxt(fileIn, delimiter=delimiter)

    # Deal with Trange input
    if Trange is None:
        Trange = [data.min(), data.max()]
    for ind, T in enumerate(Trange):
        if T is None:
            if ind == 0:
                Trange[ind] = data.min()
            elif ind == 1:
                Trange[ind] = data.max()
    
    # Add a check for the length of Trange
    if not len(Trange) == 2:
        print('Error: Trange must be of length 2')
        return
    Trange.sort()  # make sure that the minimum value is first
    # Rescale the data so be within the range 0-255.
    # Note that the upper bound needs to be made flexible if higher bit-depth
    # images are to be output.
    rescaled = (data - Trange[0]) / (Trange[1] - Trange[0]) * 255.0

    # clip out-of-bounds values to 0-255
    # upper bound needs to be adjusted for 16-bit depth images
    # # This is slow and unpythonic.  Use np.where instead?
    # for i, col in enumerate(rescaled):
    #     for j, elem in enumerate(col):
    #         if elem < 0:
    #             rescaled[i][j] =   0.0
    #         if elem > 255:
    #             rescaled[i][j] = 255.0
    rescaled[rescaled < 0.]   =   0.
    rescaled[rescaled > 255.] = 255.

    # cut last 4 chars, add PF + extn.
    fileOut = fileIn[:-4] + postfix + '.png'  

    cmap = plt.cm.gray

    # Save to file, conserving the max and min colourmap range
    plt.imsave(fileOut, rescaled, cmap=cmap, vmin=0.0, vmax=255.0)
    #scipy.misc.toimage(rescaled.round(), cmin=0.0, cmax=255.0).save(fileOut)

    print("For file %s, data range is %5.1f, %5.1f --> %s " % (fileIn, 
                                                               data.min(), 
                                                               data.max(), 
                                                               fileOut))


if __name__ == '__main__':
    fileIn   = sys.argv[1]      # input file name
    Trange   = [23., 123.]      # Temperature range
    postfix  = ''               # '_HR' #or ''

    convertCSVtoGreyscalePNG(fileIn, Trange, ';', postfix)
