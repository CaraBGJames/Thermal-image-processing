#!/usr/bin/python3

""" readConvertAsterData.py
Reads ASTER Level-1B *.hdf files and outputs the required data in georeferenced 
form.  The routines below depends heavily on the "pyhdf" module that can be 
downloaded and automatically installed by doing 

[sudo] pip install pyhdf

OR downloading and manually installing the tar file:  

wget -c http://hdfeos.org/software/pyhdf/pyhdf-0.9.0.tar.gz

A good reference for information (and one that I've shamelessly pilaged) is the 
aster_user_guide which can be found by googling it.

L1B files format typically contain the following Scientific Datasets (SDS):
0 ImageData10                   Data for band 10
1 ImageData11                   Data for band 11
2 ImageData12                   Data for band 12
3 ImageData13                   Data for band 13
4 ImageData14                   Data for band 14
5 Latitude                      Latitude matrix
6 Longitude                     Longitude matrix
7 TIR_Supplement_Chopper
8 TIR_Supplement_Encoder
9 TIR_Supplement_Temp

They also have metadata such as:
0 HDFEOSVersion                 HDF EOS version (typically v2.17)
1 StructMetadata.0              In strange format
2 coremetadata.0                Contains granule capture information such as 
                                date and time, bounding rectangle coordinates,
                                and map projection information.
3 productmetadata.0             ASTER generic metadata
4 productmetadata.1             GSD generic metadata
5 productmetadata.t             Product specific metadata for TIR: band-by-band
                                image data and statistics.


FUNCTION DESCRIPTIONS:
----------------------
- readMetadata: reads metadata stored in the hdf file and prints a list of the 
  attributes (metadata names).  If an attribute name is given, the data 
  contained within it is returned.
- metadataToDict: returns a dictionary of an attribute's metadata.
- DNToRadiance: converts scaled radiance (arbitrary units) into spectral 
  radiance in units of W/(m**2 sr micron).
- geocentricToGeodetic: converts geocentric latitude to geodetit latitude based
  on the WGS 84 datum.  Note that longitude in these two coordinate systems is 
  identical
- print_sds: returns a list of the datasets available within the hdf file.
- find_max_in_scene: returns the location (image indices) and value of the 
  hottest pixel in a scene.

TO DO:
------
- Apply the latitude and longitude
- Write functions to print and extract the metadata
- Convert metadata to a dictionary, extracting only important information.
- Make this into a class to that the ASTER file is stored as an object whose
properties can be accessed through various methods?

CHANGE LOG:
-----------
2018-06-14      DEJ     Started working on getting metadataToDict to actually 
                        output a dictionary.
"""

from __future__ import print_function
from pyhdf.SD import SD, SDC

import numpy as np
import matplotlib.pyplot as plt
import re


def readMetadata(hdf_file, attr=None, as_list=True):
    """ Reads metadata stored in the hdf file and prints a list of the 
    attributes (metadata names).  If an attribute name is given, the data 
    contained within it is returned.
    """
    metadata = hdf_file.attributes()  # reads metadata into a dictionary
    # If no attribute is stated, print out a list of the available types
    if not attr:
        # If no attribute given, enumerate tuples returned by metadata items 
        for ind, key_val_tpl in enumerate(sorted(metadata.items())):
            print('% 3d %s' % (ind, key_val_tpl[0]))  # key is first in tuple
        return 
    else:
        if as_list == True:
            return metadata[attr].split('\n') 
        else:
            return metadata[attr]  # or metadata.get(attr)


def metadataToDict(hdf_file, attr, as_list=True):
    """ Return a dictionary of an attribute's metadata.
    
    Currently the metadata is requested from 'readMetaData' in a single string 
    format, which is then split into a list of strings at the start of the for
    loop.  It may be more worthwhile to keep the single string which can then
    be searched through.
    """

    md_value = readMetadata(hdf_file, attr)
    md_dict  = {}  # An empty dictionary

    # print function to tidy up the following code
    def print_depths(line_no, obj_depth, grp_depth):
        print('match on line %4d, group depth: % 2d, object depth % 2d' % 
              (line_no, obj_depth, grp_depth))

    # These counters will record the depth of the objects/groups
    obj_depth = 0
    grp_depth = 0

    # The for loop below assumes that md_values is a list as, by default
    # readMetadata (and this function) is such.
    if as_list is not True:
        md_value = md_value.split('\n')

    for line_no, line in enumerate(md_value): 
        # Remove excess whitespace
        line = re.sub(r' ', '', line)

        # Regex pattern for searching 
        pattern = r'(END_)?(GROUP|OBJECT)'
        match = re.search(pattern, line)

        # Go through matches and, according to the depth, add keys and values
        # to md_dict
        if match:
            group = match.group()
            if group == 'GROUP':
                grp_depth += 1
                print_depths(line_no, grp_depth, obj_depth)
            elif group == 'END_GROUP':
                grp_depth -= 1
                print_depths(line_no, grp_depth, obj_depth)
            elif group == 'OBJECT':
                obj_depth += 1
                key, val = line.split(r'=')
                md_dict[key] = val
                print_depths(line_no, grp_depth, obj_depth)
            elif group == 'END_OBJECT':
                obj_depth -= 1            
                print_depths(line_no, grp_depth, obj_depth)
            else:
                continue
    return md_dict


def DNToRadiance(data, band):
    """ p. 25-26 The ASTER Level-1B data are offered in terms of scaled 
    radiance.  To convert from DN to radiance at the sensor, the unit 
    conversion coefficients (defined as radiance per 1 DN) are used.  Radiance 
    (spectral radiance) is expressed in unit of W/(m2*sr*micron).  The relation 
    between DN values and radiances is given by
    
        Radiance = (DN value - 1) x Unit conversion coefficient.
    
    NOTE:
    (i) a DN value of zero is allocated to dummy pixels
    (ii) a DN value of 1 is allocated to zero radiance
    (iii) a DN value of 254 is allocated to the maximum radiance for VNIR and 
        SWIR bands
    (iv) a DN value of 4094 is allocated to the maximum radiance for TIR bands
    (v) a DN value of 255 is allocated to saturated pixels for VNIR and SWIR 
        bands
    (vi) a DN value of 4095 is allocated to saturated pixels for TIR bands    
    """

    # Conversion factors for the TIR bands (10--14) are taken from 
    # Table 5 (p.26 of aster_user_guide_V2.pdf)
    if band == 10:
        conversionFactor = 6.822e-3
    elif band == 11:
        conversionFactor = 6.780e-3
    elif band == 12:
        conversionFactor = 6.590e-3
    elif band == 13:
        conversionFactor = 5.693e-3
    elif band == 14:
        conversionFactor = 5.225e-3
    else:
        raise ValueError('TIR bands are 10--14')
    return (data - 1.0) * conversionFactor


def geocentricToGeodetic(Latitude):
    """ Convert swath-based geocentric latitude to geodetic latitude.
    The equation used is for the WGS 84 datum only.  Geocentric longitudes do 
    not need to be converted to geodetic since they both are the same and they 
    also share the same reference meridian and axis. 

    For more information see pp. 61 & 79 of the ASTER user guide
    """
    return np.arctan((np.tan(Latitude)) / 0.99330562)


def print_sds(hdf_file):
    """ Returns a list of the datasets available within the current hdf file """
    datasets_dic = hdf_file.datasets()
    # A sorted list of the datasets
    datasets_srt = [x[0] for x in sorted(datasets_dic.items())]
    for idx, sds in enumerate(datasets_srt):
        print(idx, sds)


def find_max_in_scene(scene):
    """ Returns the location (image indices) and value of the hottest pixel
    in a scene.

    Input 'scene' must be a 2D array, i.e. intensity only
    """
    if len(scene.shape) > 2:
        print('Warning: scene must be a 2D array')
        return

    maxValue = np.nanmax(scene)
    
    M, N = scene.shape
    for i in range(M):
        for j in range(N):
            if scene[i, j] == maxValue:
                return i, j, maxValue
    # If the above code terminates and maxValue hasn't been found, print an 
    # error message and return None
    print('Warning: scene maximum not found')
    return


if __name__ == "__main__":
    from scipy.ndimage.interpolation import rotate

    import sys

    
    # open hdf file as an SD (scientific database) interface
    N_inputs = len(sys.argv)
 
    if N_inputs == 2:
        file_name = sys.argv[1]
    else:
        file_name = './Vulcano/L1B_20160915210341.hdf'
    print('Number of input args: % 2d\nInput file: %s' % (N_inputs, file_name))
    hdf_file  = SD(file_name, SDC.READ)
    
    # Longitude and latitude are given for the corners of 10 sub-squares
    Longitude = hdf_file.select('Longitude').get()
    Latitude  = hdf_file.select('Latitude').get()

    """ From aster_user_guide_v2.pdf : p.23 
    The L1B latitude and longitude geolocation arrays are two 11 x 11 matrices 
    of geocentric latitude and geodetic longitude in units of degrees. The 
    block size of the geolocation array is 420 lines by 498 samples for the 
    VNIR bands; 210 lines by 249 samples for the SWIR bands; and 70 lines by 83 
    samples for the TIR bands.
    ----------------------------------------------------------------------
    Latitude and longitude appear to vary linearly along rows and columns of 
    these matrices, so we only really care about the lat/long values at the 
    four corners for the georeferencing of a scene.  You'll notice that the 
    useful data in the scenes represent some kind of skewed and sheared 
    rectangle.  Check out LongLatTest.pdf to see what I mean here.  
    """
    
    # UPPERLEFT, UPPERRIGHT, LOWERLEFT, LOWERRIGHT corners
    SCENEFOURCORNERS = [[Longitude[0][0], Latitude[0][0]],
                        [Longitude[-1][0], Latitude[-1][0]],
                        [Longitude[0][-1], Latitude[0][-1]],
                        [Longitude[-1][-1], Latitude[-1][-1]]]

    # REPRESENTATION OF THE LOCATIONS OF THE LONGITUDE AND LATITUDE MATRICES
    plt.figure(figsize=(8, 6))
    xs = [x[0] for x in SCENEFOURCORNERS]
    ys = [x[1] for x in SCENEFOURCORNERS]
    plt.plot(xs, ys, 'ko', label='scene corners')
    for lng, lat in zip(Longitude, Latitude):
        plt.plot(lng, lat, 'r.', label='submatrix corners')
    
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Longitude and latitude of the scene boxes')
    #plt.legend(loc='best')

    plt.savefig('LongLatTest.pdf', bbox_inches='tight')
    plt.close()

    # Dictionary of conversion factors for the scaled radiance data.  This is
    # the same information as in DNtoRadiance but in a dictionary form.  Thus, 
    # to access the CF for band 10, say, do: 
    # >>> cf = conversionFactors[10]
    conversionFactors = {10: 6.822e-3, 
        11: 6.780e-3, 
        12: 6.590e-3, 
        13: 5.693e-3, 
        14: 5.225e-3}

    # Plot the data to make it easier to understand
    plt.figure(figsize=(8.27, 11.69))
    for band in range(10, 15):
        data = hdf_file.select('ImageData' + str(band)).get()

        # Note that this angle (in degrees) is for the Soufriere-Guadeloupe
        # L1B_20160926023116.hdf granule only
        angle = -8.62486636728127

        data = rotate(data, angle)
        # Convert the data to radiance [W/(m**2 sr um)]
        convertedData = DNToRadiance(data, band)

        # Radiance values shouldn't be negative - convert to nans
        convertedData[convertedData < 0] = np.nan

        # Plot each band in a separate subplot
        plt.subplot(3, 2, band - 9)
        im = plt.imshow(convertedData)
        im.set_clim((6, 12))
        cbar = plt.colorbar()
        cbar.set_label('W/(m$^2$sr$\mu$m)')
        plt.title('TIR Band ' + str(band))

    plt.savefig('test.pdf', bbox_inches='tight', dpi=300)
    plt.close('all')
        
        
        
        
