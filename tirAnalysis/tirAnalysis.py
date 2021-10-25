#!/usr/bin/python

"""
tirAnalysis.py : this code defines functions and returns values for the 
heat flux calculated from TIR images using image mosaics.  Uses afine 
transform parameters from QGIS georeferencing.  

Several layers are imported:
Tlayer   : thermal image mosaic - a greyscale image georeferenced to the DEM
DEMlayer : DEM

CHANGE LOG:
2018-07-13      Converting functions to be python3 compatible.  This involves
                removing all the qgis dependencies, which only support 
                python2.7
2018-07-14      Wrote a plotting function for the flux maps so that the same
                code can be applied to the radiative, convective and total 
                surface flux cases.
                Note that there is a discrepancy between calculations of 
                radiative flux involving this new function

TO DO:
- Set cbar upper lim to be equal to the last cbar label
- Is pixelAreaInMeters actually in meters?
- Account for atmospheric transmissivity
- Account for ground emissivity (partially done) - make emissivity map
- Formally calculate the heat transfer coefficient from base data.
- Calculate Sekioka + Yuhara geothermal flux
"""

#from qgis.utils import iface.activeLayer
#from qgis.core import QgsRasterLayer, QgsMapLayerRegistry
#from PyQt4.QtCore import QFileInfo # Not being used.
# try:
#     from osgeo import gdal #, osr
# except ImportError:
#     import gdal #, osr
from skimage.io import imread
#from scipy.misc import imsave
from colormaps.colormaps import cmaps  # new matplotlib colormaps
from cartography.cartography import latLonElevToUTM
from CoolProp.CoolProp import PropsSI
from thermo_fluids.heat_mass_transfer import heat_transfer_coeff
from scipy.constants import sigma  # sigma = 5.67e-8 Wm**-2K**-4

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import tifffile as tiff
import pandas
import utm


def openRaster(fileName):
    """
    Opens a raster layer from a file name using gdal utilities.  This 
    function is not currently being used, and will be removed in the future.
    """
    try:
        # fileInfo = QFileInfo(fileName)
        # path     = fileInfo.filePath()
        # baseName = fileInfo.baseName()
        #layer    = QgsRasterLayer(path, baseName)
        layerDS  = gdal.Open(filename)
        layer    = layerDS.GetRasterBand(1).ReadAsArray()
    except:
        print("Failed to open layer")
    # QgsMapLayerRegistry.instance().addMapLayer(layer)

    # if layer.isValid():
    #     print("Successful")
    #     return layer
    # else:
    #     print("Unsuccessful")
    #     return
    return layer


def pixel2coord(col, row, transform):
    """Returns global coordinates to pixel center using base-0 raster index"""
    c, a, b, f, d, e = transform
    xp = a * (col + .5) + b * (row + .5) + c
    yp = d * (col + .5) + e * (row + .5) + f
    return np.array((xp, yp)).T


# def val_raster(point, raster):
#     """
#     Returns the values of a raster layer at a location point.  Depends upon
#     gdal and thus will be removed
#     """
#     ident = raster.dataProvider().identify(point, 
#                                            QgsRaster.IdentifyFormatValue)
#     return ident.results().values()


def atmosphericTransmissivity(dist, atmConds=None, CO2=400, CH4=1.7):
    """
    Returns the transmissivity of an atmosphere over a given distance.

    Further to this, some assumptions for the concentration of certain
    gas species and scattering due to aerosols must be made.  These are 
    summarised as:
    SPECIES     CONC.   UNITS
    CO2         400     ppm
    CH4         1.7     ppm
    aerosols    

    N.B. The current version ignores all these details and instead returns
    a fit to the data of Fig. 7.7, Harris (2013), pp. 430

    Parameters:
    -----------
    dist : float
        The path length (in metres).
    atmConds : pandas.DataFrame
        Atmospheric conditions: 
                TEMP/[deg C]
            PRESSURE/[mbar]
            HUMIDITY/[%]
    CO2 : float
        Atmospheric concentration of CO2
    CH4 : float
        Atmospheric concentration of CH4

    Returns:
    --------
    Transmissivity : float
    """
    transPer500m = .88    # From Harris (2013)
    transPer1m   = np.power(transPer500m, 1/500)

    # This function derived by fitting data from Harris (2013), pp. ???
    return .9289976938 * np.exp(-1.235159e-4 * dist)

#    return 1.   # transPer1m ** dist


def unpackLatLonElev(latLonElevTuple):
    """ 
    Unpack lat, lon, elevation tuple.  
    If no elevation data is available, set this to zero.

    Parameters:
    -----------
    latLonElevTuple : tuple or array-like
        Tuple containng latitude, longitude, elevation 

    Returns:
    --------
    lat : float
        latitude from latLonElevTuple
    lon : float
        longitude from latLonElevTuple        
    elev : float
        elevation from latLonElevTuple
    """
    if len(latLonElevTuple) == 3:
        lat, lon, elev = latLonElevTuple
    elif len(latLonElevTuple) == 2:
        lat, lon = latLonElevTuple
        elev     = 0.
    return lat, lon, elev


def distanceToTarget(latLonElevObserver, latLonElevTarget, isUTM=True):
    """
    Returns the approximate distance from the observer to the target (in 
    meters) based on a pair of latitude-longitude-elevation tuples.  Assumes
    flat-world geometry.

    Parameters:
    -----------
    latLonElevObserver : tuple or array-like
        A list, array or tuple containing the latitude, longitude and 
        elevation of the observer station.
    latLonElevTarget : tuple or array-like
        A list, array or tuple containing the latitude, longitude and 
        elevation of the target object.
    isUTM : boolean
        Coordinates in (local) UTM format?  (as opposed to global, degrees)

    Returns:
    --------
    dist : float
        Distance from observer to target.
    """
    R = 6371e3  # Mean radius of earth in meters (from Wikipedia)
    if ~isUTM:
        # Observer
        obsEasting, obsNorthing, obsElev = latLonElevToUTM(latLonElevObserver)
        obsVec = np.array([obsEasting, obsNorthing, obsElev])
        print("Observer UTM: %.2f E, %.2f N, %d m a.s.l." % tuple(obsVec))

        # Target
        tarEasting, tarNorthing, tarElev = latLonElevToUTM(latLonElevTarget)
        tarVec = np.array([tarEasting, tarNorthing, tarElev])
        print("Target   UTM: %.2f E, %.2f N, %d m a.s.l." % tuple(tarVec))
    else:
        obsVec = np.array(unpackLatLonElev(latLonElevObserver))
        tarVec = np.array(unpackLatLonElev(latLonElevTarget))
    # Calculate distance using vector L2 norm (shorter than sqrt method)
    dist = np.linalg.norm(obsVec - tarVec)
    return dist #, obsVec, tarVec


def UTM2tick(eastings, northings, transform):
    """ A utility for converting UTM loactions to tick marks using the 
        GDAL-defined afine transform """
    return ((eastings - transform[0])/transform[1],
            (northings - transform[3])/transform[5])


def tick2UTM(xticks, yticks, transform):
    """ A utility for converting tick marks to UTM loactions using the 
        GDAL-defined afine transform """
    return (transform[0] + xticks * transform[1], 
            transform[3] + yticks * transform[5])


def saveAsGeotiff(destFilename, data, transform, mask=None):
    """ 
    Saves input "image" data, and georeferencing transform in geotif format 
    
    Parameters:
    -----------
    destFilename : str
        The file to which data should be written
    data : array-like
        Input data to be saved to file
    transform : list or array-like
        The geotiff transform array consisting of 6 elements: 
        pixelWidth, shearXY, shearYX, pixelHeight, originX, originY
        N.B. this order is different from the GDAL format (X0, W, XY, Y0, YX, H).
        To convert from GDAL do: transform[np.array([1,2,4,5,0,3])]
    mask : array-like boolean
        
    """

    H, W = data.shape

    # Check that the output filename extension is for "tiff" format
    ext = destFilename.split('.')[-1]
    if ext[:3] != 'tif':        # Also covers .tiff case 
        destFilename += '.tif'
    
    if mask is not None:
        try:
            data[~mask] = 0
        except IndexError as e:
            print(e)
        finally:
            tiff.imsave(destFilename, data)
            
    worldCoordsFile = destFilename + 'w'
    with open(worldCoordsFile, 'w') as f:
        for val in transform:
            f.write(str(val))
            f.write('\n')
            

def plotFluxes(flux, fluxType, figAx, extent=None, clim=None):
    """
    A plotting utilty to ensure that the same treatment is applied to 
    all of the fluxes.
    
    Parameters
    ----------
    flux : array-like
       Geotransformed flux image
    fluxStyle : str
        'Radiative', 'Convective' or 'Surface'
    figAx : tuple
        matplotlib fig and axes
    extent : list or array-like
        extents of the image (left, right, bottom, top)
    CLim : tuple (of ints/floats)
        Upper and lower intensity limits for the image

    Returns
    -------
    QTot : array-like
    
    """
    labelFontSize = 8
    fig, ax = figAx
    CLim = clim
    #eastings, northings = EN

    im = ax.imshow(flux, extent=extent,
                   cmap=cmaps.get('inferno')) #, clim=CLim)# magma
    if CLim is not None:
        im.set_clim(CLim)
        
    print("%s CLIMs : " % fluxType, im.get_clim())

    # colourbar properties
    cbar = fig.colorbar(im, shrink=.65)
    cbar.set_label(label='%s heat flux/[W/m$^2$]' % fluxType.capitalize(), 
                   rotation=270, labelpad=15)
    labels = cbar.ax.get_yticklabels()
    for label in labels:
        label.set_fontsize(labelFontSize)
    # # Set last cbar label to be '> 150' (for example)
    # labels[-1].set_text('> ' + labels[-1].get_text())
    # cbar.ax.set_yticklabels(labels)
    
    # labels + title
    plt.xlabel('UTM Easting/[metres]')
    plt.ylabel('UTM Northing/[metres]')
    
    # Calculate total fluxes and convert to MW for easier typing later
    QTot = np.nansum(flux) * 1e-6 * pixelAreaInMeters
    
    plt.title('Total %s flux: %.2f MW' % (fluxType.lower(), QTot))
    
    # The following code stops the axis tick labels from being
    # cut off or shortened in any way
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useOffset=False))
    
    plt.draw()
    fig.tight_layout()
    if fluxType.lower() == 'radiative':
        keepChars = 3
    else:
        keepChars = 4
        outname = 'Q%s' % fluxType.lower()[:keepChars]
        plt.savefig(workingDir + outname + '.pdf', bbox_inches='tight')
        plt.savefig(workingDir + outname + '.png', bbox_inches='tight')
        plt.close()
    return QTot


if __name__ == "__main__":
    # Import layers
    # Thermal layer
    SITE       = '11'
    workingDir = ('/home/david/FieldWork/Vulcano/2016-09-Vulcano/' + 
                  '2016-09-17_nighttimeFLIR/')
    fileName    = (workingDir + './SITE' + SITE +
                   '/georeferenced/SITE' + SITE + '_georef.tif')
    TlayerDS    = gdal.Open(fileName) 
    # Temp data is in first layer (or layers 1-3), transparency is last layer
    alphaLayer = TlayerDS.RasterCount;
    if alphaLayer == 2:
        try:
            Tlayer = TlayerDS.GetRasterBand(1).ReadAsArray()
        except AttributeError:
            pass
    if alphaLayer == 4:
        # Read each of the layers as one of the RGB bands
        R = TlayerDS.GetRasterBand(1).ReadAsArray()
        G = TlayerDS.GetRasterBand(2).ReadAsArray()
        B = TlayerDS.GetRasterBand(3).ReadAsArray()
        # Depth-wise stack the bands, and dot them with the "MATLAB" grey-scale
        # conversion vector
        Tlayer = np.dstack((R, G, B)).dot([.2989, .5870, .1140])

    # Thermal image contains two layers, greyscale and alpha [0 255].  If alpha
    # layer has value = 255, that corresponding pixel is transparent.  Create a 
    # mask to exploit this and only extract temperature if mask[i,j] == 1
    mask = TlayerDS.GetRasterBand(alphaLayer).ReadAsArray() == 255

    # DEM layer
    DEMfileName = (workingDir +
                   '../../DEM/06_final_20cm_model_ortho/dem_20cm_hshd.bmp')
    DEMlayerDS  = gdal.Open(DEMfileName)

    # Unravel GDAL affine transform parameters
    transform = TlayerDS.GetGeoTransform()

    # Load log file for TIR measurements and drop invalid columns
    names = [u'POS', u'GPS E', u'GPS N', u'TIME', u'TEMP', u'HUMIDITY',
             u'PRESSURE', u'CLOUD', u'NOTES']
    df = pandas.read_excel(workingDir + '2016-09-17_nighttimeFLIR-log.xls', 
                           sheet_name='Sheet1', names=names,
                           header=2, skiprows=0).drop(['NOTES'], axis=1)
    # Cloud cover column has non-numeric values - convert these
    df.loc[df.CLOUD == "< 5", 'CLOUD'] = 1  # suppose that 1 % < 5 % ...
    U, L = 5, 5
    T, p = df['TEMP'].mean() + 273.16, df['PRESSURE'].mean() / 10
    
    # Get width and height of raster layer (np.array -> shape).  
    H, W = Tlayer.shape

    # Ordinarily we would use the vectors joining nodes on two of the sizes of 
    # a pixel and take the cross product to define the area of the pixel.
    # However, all pixels are the same size so we can get away with being lazy 
    # here simplifying.  Is this truly in meters???
    northingPixelSize =  transform[1]
    eastingPixelSize  = -transform[5]
    pixelAreaInMeters =  northingPixelSize * eastingPixelSize
    # How many pixels actually correspond to thermally-active zones
    areaCount = np.count_nonzero(mask)

    # Calculate path length from each pixel to the measurement station
    easting  = df[df['POS'] == int(SITE)]['GPS E'].values[0]
    northing = df[df['POS'] == int(SITE)]['GPS N'].values[0]
    obsVec   = unpackLatLonElev((easting, northing)) 
    x = np.arange(W) * transform[1] + transform[0]
    y = np.arange(H) * transform[5] + transform[3]
    X, Y = np.meshgrid(x, y)
    #Y    = np.flipud(Y)
    pathLength = np.zeros_like(X)
    for i, row in enumerate(mask):  
        for j, val in enumerate(row):
            if val:     # Only do the operations if the mask value is True
                tarVec = np.array(unpackLatLonElev((X[i,j], Y[i,j])))
                # dist[i,j] = distanceToTarget((easting, northing),
                #                              (X[i,j], Y[i,j]), isUTM=True)
                pathLength[i,j] = np.linalg.norm(obsVec - tarVec)

    # Define emissivities (8-14 micron band)
    epsSurface = 0.975          # See table 2.3 of Harris (2013, pp. 81)
    tauAtmosph = atmosphericTransmissivity(pathLength)  # Fill in function
    epsAtmosph = 1.0 - tauAtmosph
    # Heat transfer coefficient.  # hc = .4
    hc, *foo = heat_transfer_coeff(U, L, T, p)  


    # Define dict of temperature ranges, and select that for the current site
    TrangeDict = {}
    TrangeDict['01'] = np.array([15.6,  60.6]) + 273.14
    TrangeDict['02'] = np.array([18.6, 101. ]) + 273.14
    TrangeDict['05'] = np.array([19.8,  99.3]) + 273.14
    TrangeDict['07'] = np.array([20.4,  63.4]) + 273.14
    TrangeDict['10'] = np.array([19.9,  99.7]) + 273.14 
    TrangeDict['11'] = np.array([12.7, 107. ]) + 273.14 

    Trange = TrangeDict[SITE]

    # Convert pixel values to real (absolute) temperatures
    TempField  = Tlayer / 255.0               # Normalise values
    TempField *= Trange.max() - Trange.min()  # Scale by temperature range
    TempField += Trange.min()                 # Add offset

    # Ambient temperature for measurements taken at SITE05
    Tamb        = df[df.POS == 5].TEMP.values[0] + 273.14
    TambPow4    = Tamb ** 4  
    epsTambPow4 = epsAtmosph * TambPow4    

    # Convert from brightness to kinetic temperature
    Tstar = np.copy(TempField)
    T     = np.power((Tstar**4 - epsTambPow4)
                     / (epsSurface * tauAtmosph), .25)
    
    Qrad  = (epsSurface * sigma * (T**4 - TambPow4))
    Qconv = hc * (T - Tamb)
    Qsurf = Qrad + Qconv

    ## PLOT AND SAVE FLUX MAPS

    # Maxima and minima for the plots in UTM
    yMax, xMin, foo, bar = utm.from_latlon(transform[0]/1e5, transform[4]/1e5)
    yMin, xMax, foo, bar = utm.from_latlon((transform[0]+H*transform[1])/1e5,
                                           (transform[3] + W*transform[5])/1e5)

    # Define UTM eastings and northings that will be used on the plots
    eastings  =  496000 + np.arange(500, 1000, 100)
    northings = 4250000 + np.arange(600, 950, 50)
    newXticks, newYticks = UTM2tick(eastings, northings, transform)
    extent = [transform[0], transform[0] + H * transform[1],
              transform[3] + W * transform[5], transform[3]]

    # RADIATIVE FLUXES
    fig, ax = plt.subplots(figsize=(8,6), dpi=300)
    QradOrig = Qrad.copy()
    Qrad[T <= Trange[0]] = None  # So that the "empty" pixels are "transparent"
    QradTot = plotFluxes(Qrad, 'Radiative', (fig, ax),
                         extent=extent) #, clim=(0., 140.), )
    
    # CONVECTIVE FLUXES
    fig, ax = plt.subplots(figsize=(8,6), dpi=300)
    QconvOrig = Qconv.copy()
    Qconv[T <= Trange[0]] = None #
    QconvTot = plotFluxes(Qconv, 'Convective', (fig, ax),
                          extent=extent) #clim=(0., 10.), 

    # TOTAL SURFACE FLUXES (= QradTot + QconvTot)
    fig, ax = plt.subplots(figsize=(8,6), dpi=300)
    QsurfOrig = Qsurf.copy()
    Qsurf[T <= Trange[0]] = None #
    QsurfTot = plotFluxes(Qsurf, 'Surface', (fig, ax),
                          extent=extent) #clim=(0., 140.), 

    # # Calculate total fluxes and convert to MW for easier typing later
    # QradTot  = np.nansum(Qrad)  * 1e-6 * pixelAreaInMeters 
    # QconvTot = np.nansum(Qconv) * 1e-6 * pixelAreaInMeters 
    # QsurfTot = np.nansum(Qsurf) * 1e-6 * pixelAreaInMeters 

    # Total radiative flux, not accounting for atmospheric and ground 
    # emissivities
    print('')  # Clear line
    print('Total radiative  flux : %.2f MW' % (QradTot))
    print('Total convective flux : %.2f MW' % (QconvTot))
    print('Total surface flux    : %.2f MW' % (QsurfTot))
    print('#-----Sanity check: %.2f = %.2f + %.2f = %.2f-----#\n' % (
            QsurfTot, 
            QradTot, 
            QconvTot, 
            QradTot + QconvTot))

    # SAVE A COPY OF THE BASIC IMAGE
    fileName = workingDir + './SITE' + SITE + '/stitched/SITE' + SITE + '.tif'
    im = imread(fileName, as_grey=True, flatten=False)
    # Convert to grey scale as per MATLAB
    imGrey = np.dot(im[...,:3], [.2989, .5870, .1140])
    # Normalise to [0, 1]
    imGrey /= 255.0
    # Rescale to temperature range
    imRescaled = imGrey * (Trange.max() - Trange.min()) + Trange.min()
    # Set transparent pixels in original image to be beyond the T range
    imRescaled[im[...,3] == 0.0] = 0.0  # Absolute zero!
    

    ## PLOT THE IMAGE AND SAVE
    fig, ax = plt.subplots(figsize=(8,6), dpi=300)
    im = ax.imshow(imRescaled, cmap=cmaps.get('inferno'), 
                   clim=[Trange.min(), 325]) #Trange.max()])
    tickPosns = np.arange(300, 330, 5) 
    cbar = fig.colorbar(im, shrink=.4, ticks=tickPosns)
    cbar.set_label('Apparent temp./[K]', rotation=270, labelpad=15)

    # Set cbar tick labels
    tickPosnLabels = list(map(str, tickPosns))
#    tickPosnLabels[0] = '< ' + tickPosnLabels[0]
    tickPosnLabels[-1] = '> ' + tickPosnLabels[-1]
    cbar.ax.set_yticklabels(tickPosnLabels)

    fig.tight_layout()
    plt.savefig(workingDir + './SITE' + SITE + '/stitched/SITE'
                + SITE + '-cbar.pdf', 
                bbox_inches='tight')
    plt.close()

    # Save the emissivity and transmissivity corrected georeferenced
    # temperature map as a geotif
    filename = (workingDir + './SITE' + SITE + '/georeferenced/SITE'
                + SITE + '_T_georef.tif')
    transform = np.array(transform)[np.array([1,2,4,5,0,3])]
    saveAsGeotiff(filename, T, transform, mask=mask)


def measure_hotspot_area(image, thresh, pixel_area=None):
    '''
    Return the area of hotspots in a thermal image, where these are
    defined as having a temperature above a given value.

    Uses the watershed segmentation method from scikit-image

    Parameters
    ----------
    image : image or array-like
        The (2D) image to be processed
    thesh : float
        A threshold value for segmentation.  
        See skimage.filters.threshold_* for useful implementations of many
        common threshold filters.
    pixel_area : float (optional)
        A physical conversion factor (e.g. m^2/pixel^2) that will be 
        applied to all returned values, if supplied.

    Code pilfered from 
        https://scipy-lectures.org/packages/scikit-image/index.html
    '''

    from skimage.feature import peak_local_max
    from skimage import morphology
    from skimage import measure
    from scipy import ndimage


    # Create mask of all pixels above threshold value
    blobs = image > thresh

    # Now we want to separate the objects in the image
    # Generate the markers as local maxima of the distance
    # to the background.
    distance  = ndimage.distance_transform_edt(blobs)
    local_max = peak_local_max(distance, indices=False,
                               footprint=np.ones((3, 3)), labels=blobs)
    # Create markers based on the local maxima and use these to
    # segment the image.
    markers = morphology.label(local_max)
    labels  = morphology.watershed(-distance, markers, mask=blobs)

    # Compute size, perimeter etc of the segmented regions.
    properties = measure.regionprops(labels)
    areas      = [prop.area for prop in properties]

    if pixel_area is not None:
        areas = list(map(lambda x: x * pixel_area, areas))

    return areas, labels
