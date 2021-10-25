#!/usr/bin/python3

""" irbisCSV2file.py
Provides a module for reading a csv file produced by IRBIS 3 (IR camera 
software) and outputting images.  There is an option to produce either 
standard 8-bit images or 16-bit "bigTiffs".

This file was originally called readCSVfiles.py

If called as a stand-alone function, use multiprocessing module on all csv 
file in the current directory.

To do:
------
- write exifdata to image file

"""

import numpy as np
from matplotlib import cm as cm
#import piexif


def irbisCSV2data(fname, sep=','):
    """
    Returns the data and metadata (header info) contained within "fname"

    Parameters:
    -----------
    fname : str
        the name of the file to open
    sep : str, optional
        the separator of the input file.  Default is "comma"

    Returns:
    --------
    data : np.array
        the data 
    metadata : dict
        the header information

    TO DO
    -----
    - Include support for (g)zipped files
    """
    f = open(fname, 'r', encoding='utf-8', errors='ignore')

    found = False
    data  = []
    nread = 0
    metadata = {}
    new_metadata = {}
    mdread = 0

    ##--------------------READ METADATA--------------------##
    key = None
    while key != '[Data]':
        # Remove leading and trailing whitespace, as well as any newline
        # characters.  Split the remaining output.
        key, *params = map(str.strip, f.readline().split('='))
        if len(params) > 0: 
            params = params[0].split(';')
        if mdread > 0 and key != '':
            metadata[key] = params
        mdread += 1

    try:
        # Loop once over dict to remove [*] keys
        for key in metadata.keys():
            if '[' not in key:
                new_metadata[key] = metadata[key]
        metadata = new_metadata
        # Loop again to clean up some other stuff.  Can these two loops
        # be combined?  Also, when type(vals) == list and len(vals) == 1,
        # can vals be converted to another data type, such as float or int?
        for key in metadata.keys():
            if key != 'TempUnit':
                vals = metadata[key]
                floatVals = []
                if type(vals) == list and len(vals) > 1:
                    for val in vals:
                        floatVals.append(float(val))
                else:
                    floatVals = float(vals[0])
                vals = floatVals
            else:
                vals = metadata[key][0]  # 'TempUnit'
            metadata[key] = vals
    except ValueError:
        pass
        
    ##--------------------READ THE DATA--------------------##
    f.seek(0)
    for line in f:
        if found:
            line = line.strip().split(sep)
            # Sometimes an additional, empty element is found at the 
            # end of the line.  Cut it off.
            if line[-1] == '':
                line = line[:-1]
            data.append(line)
            nread += 1  # should equal 480 if the code completes successfully
        else:
            if line.strip() == '[Data]': found = True

    f.close()

    data = np.float64(data)

    return data, metadata


def irbisData2File(savename, data, limits=(None, None),
                   filetype='png', metadata=None, bigtiff=False, 
                   focal_length=30, verbose=False, cmap=cm.gray):
    ''' 
    Takes a data array and saves it as an image. 
    
    Parameters:
    -----------
    savename : str
        the root of the filename under which the file should be saved.  Note 
        that the postfix '.{png,jpg,tif}' (or '_bt.tif' if bigtiff is True) 
        will be added to this.
    data :  array_like
        the data array to save

    Optional:
    ---------
    limits : tuple or list
        lower and upper bounds of the data range to be entered as a 
        tuple (or list).  Defaults to (None, None)
    filetype : str
        Post-fix or extention type of the image file
    metadata : dict
        exif data to be written to a metadata header
    bigtiff : boolean
        option to save as a 16-bit tiff
    focal_length : float
        focal length of the lens used.  Will be added to metadata (TODO!!)
    verbose : bool
        Verbose output of filenames (useful for batch runs)
    cmap : matplotlib.cm.colormap
        colourmap for plotting
    
    TO DO:
        - EXIF write function
        - Add focal length info to MD
    '''

    from matplotlib.cm import gray

    # Clip the data to the range of "limits", if given
    TLimits = [data.min(), data.max()]
    if any(limits) is not None:  # This if statement probably unnecessary
        if limits[0] is not None:
            TLimits[0] = limits[0]
        if limits[1] is not None:
            TLimits[1] = limits[1]
        TRange = TLimits[1] - TLimits[0]
        data[data < TLimits[0]] = TLimits[0]
        data[data > TLimits[1]] = TLimits[1]
    TRange = np.diff(TLimits)
        
    # Options according to the required image bit depth
    if bigtiff:
        from tifffile import imsave
        savename += '_bt.tif'
        bitDepth = 2**16 - 1    # N.B. range is 0--65536
        imsave(savename, data)
    else:
        # imsave is depreciated in SciPy v1.0.0 and will be replaced in 
        # v1.2.0 by imageio.imwrite (or use {plt,tifffile}.imsave?)
        from matplotlib.pyplot import imsave
        # from imageio import imwrite as imsave
        savename += '.' + filetype
        bitDepth = 2**8 - 1     # N.B. range is 0--255
        # Create rescaled image
        rescaledData = np.round((data - TLimits[0]) / TRange * bitDepth)
        if verbose: print(rescaledData.min(), rescaledData.max())
        # Save the image to file whilst keeping the grey levels constant
        imsave(savename, rescaledData, cmap=cmap, vmin=0, vmax=bitDepth)

        
    if metadata:
        # Create exif data from metadata
        zeroth_ifd = {}
        # exif_ifd = {piexif.ExifIFD.FNumber :         1.,
        #             piexif.ExifIFD.FocalLength :     u'%d mm' % focal_length,
        #             piexif.ExifIFD.LensMake :        u'Janoptik',
        #             piexif.ExifIFD.SubjectDistance : u'4.59 m',
        # }

    if verbose:
        print('Saving file %s...  Data limits: [%3.1f, %3.1f]' %
              (savename, TLimits[0], TLimits[1]))


def aveStdImages(root, start=1, stop=3000, sep=','):
    '''
    Returns mean and std dev. images from a sequence of thermal images

    Parameters
    ----------
    root : str
        base part of the file name, up to underscore separating the integer 
        sequence.  Assumes filenames of the form './AA080104_%04d.csv'
    start : int
        starting image in the sequence to be used.  Default is first image.
    stop : int
        last image in sequence to be used.  Default is 3000th image.
    sep : str
        separator in the CSV files.  Default is a comma ','

    '''
    T_min, T_max = 200., 0.  # Minimum and maximum temperatures
    data, md = irbisCSV2data(root + '_%04d.csv' % start, sep=sep)
    imSum = np.zeros_like(data)
    counter = 0
    for i in range(start, stop):
        filename = root + '_%04d.csv' % i
        try:
            data, md = irbisCSV2data(filename, sep=sep)
            imSum   += data
            counter += 1
            if data.min() < T_min:
                T_min = data.min()
            if data.max() > T_max:
                T_max = data.max()
        except:
            pass
    imAve = imSum / counter

    imSum = np.zeros_like(data)
    counter = 0
    for i in range(start, stop):
        filename = root + '_%04d.csv' % i
        try:
            data, md = irbisCSV2data(filename, sep=sep)
            imSum   += (data - imAve)**2
            counter += 1
        except:
            pass
    imStd = np.sqrt(imSum / (counter - 1))

    return imAve, imStd, T_min, T_max


if __name__ == '__main__':
    '''
    If called as a stand-alone function, use multiprocessing module on all 
    csv files in the current directory.
    '''
    from multiprocessing import Pool
    from glob import glob

    def multiprocess_allcsv(filename):
        data, md = irbisCSV2data(filename, sep=',')
        irbisData2File(filename[:-4], data, limits=(17, 36))

        
    filenames = glob('*csv')
    p = Pool(6)
    p.map(multiprocess_allcsv, filenames)
