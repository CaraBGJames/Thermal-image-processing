try:
    from osgeo import gdal #, osr
except ImportError:
    import gdal #, osr
from skimage.io import imread
#from scipy.misc import imsave
from colormap.colormaps import cmaps  # new matplotlib colormaps


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import tifffile as tiff
import pandas
import utm


# Import layers
# Thermal layer
SITE       = 'FTY01'
workingDir = ('/home/david/FieldWork/Soufriere-Guadeloupe/Thermography/'
              + '2018-06-27_FTY/FTY01/TIR/tir/')
fileName    = (workingDir + '/georeferenced/20180627_FTY01_georef.tif')
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

# Unravel GDAL affine transform parameters
transform = TlayerDS.GetGeoTransform()

# # Load log file for TIR measurements and drop invalid columns
# names = [u'POS', u'GPS E', u'GPS N', u'TIME', u'TEMP', u'HUMIDITY',
#          u'PRESSURE', u'CLOUD', u'NOTES']
# df = pandas.read_excel(workingDir + '2016-09-17_nighttimeFLIR-log.xls', 
#                        sheet_name='Sheet1', names=names,
#                        header=2, skiprows=0).drop(['NOTES'], axis=1)
# # Cloud cover column has non-numeric values - convert these
# df.loc[df.CLOUD == "< 5", 'CLOUD'] = 1  # suppose that 1 % < 5 % ...

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

# # Calculate path length from each pixel to the measurement station
# easting  = df[df['POS'] == int(SITE)]['GPS E'].values[0]
# northing = df[df['POS'] == int(SITE)]['GPS N'].values[0]
elev, easting, northing = 1164.236,643155.051,1773691.924
obsVec   = unpackLatLonElev((easting, northing)) 

x    = np.arange(W) * transform[1] + transform[0]
y    = np.arange(H) * transform[5] + transform[3]
X, Y = np.meshgrid(x, y)

pathLength = np.zeros_like(X)
for i, row in enumerate(mask):  
    for j, val in enumerate(row):
        if val:     # Only do the operations if the mask value is True
            tarVec = np.array(unpackLatLonElev((X[i,j], Y[i,j])))
            # dist[i,j] = distanceToTarget((easting, northing),
            #                              (X[i,j], Y[i,j]), isUTM=True)
            pathLength[i,j] = np.linalg.norm(obsVec - tarVec)
            
# Define Stefan-Boltzmann constant and emissivities (8-14 micron band)
stefanBoltzmann = 5.67e-8   # W m**-2 K**-4
epsSurface = 0.975          # See table 2.3 of Harris (2013, pp. 81)
tauAtmosph = atmosphericTransmissivity(pathLength)  # Fill in function
epsAtmosph = 1.0 - tauAtmosph
# Heat transfer coefficient.  This needs to be formally calculated
hc = .4

# Define dict of temperature ranges, and select that for the current site
TrangeDict = {}
TrangeDict['01'] = np.array([15.6, 60.6]) + 273.14  # Kelvin
TrangeDict['02'] = np.array([18.6, 101.]) + 273.14  # Kelvin
TrangeDict['05'] = np.array([19.8, 99.3]) + 273.14  # Kelvin
TrangeDict['07'] = np.array([20.4, 63.4]) + 273.14  # Kelvin
TrangeDict['10'] = np.array([19.9, 99.7]) + 273.14  # Kelvin 
TrangeDict['11'] = np.array([12.7, 107.]) + 273.14  # Kelvin 

Trange = TrangeDict[SITE]

# Convert pixel values to real (absolute) temperatures
TempField  = Tlayer / 255.0              # Normalise values 
TempField *= Trange.max() - Trange.min() # Scale by temperature range
TempField += Trange.min()                # Add offset

# Ambient temperature for measurements taken at SITE05
Tamb        = df[df.POS == 5].TEMP.values[0] + 273.14  # Kelvin
TambPow4    = Tamb ** 4  # To avoid excessive operations
epsTambPow4 = epsAtmosph * TambPow4    

# Pull all the row and col calculations out of the loop to save time
Tstar = np.copy(TempField)
T     = np.power((Tstar**4 - epsTambPow4) 
                 / (epsSurface * tauAtmosph), .25)

Qrad  = (epsSurface * stefanBoltzmann* (T**4 - TambPow4))
Qconv = hc * (T - Tamb)
Qsurf = Qrad + Qconv

## PLOT AND SAVE FLUX MAPS

# Maxima and minima for the plots in UTM
yMax, xMin, foo, bar = utm.from_latlon(transform[0]/1e5, transform[4]/1e5)
yMin, xMax, foo, bar = utm.from_latlon((transform[0]+H*transform[1])/1e5,
                                       (transform[3] + W*transform[5])/1e5)

# Define UTM eastings and northings that will be used on the plots
eastings  = 496000  + np.arange(500, 1000, 100)
northings = 4250000 + np.arange(600, 950, 50)
newXticks, newYticks = UTM2tick(eastings, northings, transform)
extent = [transform[0], transform[0] + H * transform[1],
          transform[3] + W * transform[5], transform[3]]

# RADIATIVE FLUXES
QradOrig = Qrad.copy()
Qrad[T <= Trange[0]] = None  # So that the "empty" pixels are "transparent"
fig, ax = plt.subplots(figsize=(8,6), dpi=300)
QradTot = plotFluxes(Qrad, 'Radiative', (fig, ax),
                     extent=extent) #, clim=(0., 140.), )
    
# CONVECTIVE FLUXES
QconvOrig = Qconv.copy()
Qconv[T <= Trange[0]] = None #
fig, ax = plt.subplots(figsize=(8,6), dpi=300)
QconvTot = plotFluxes(Qconv, 'Convective', (fig, ax),
                      extent=extent) #clim=(0., 10.), 

# TOTAL SURFACE FLUXES (= QradTot + QconvTot)
QsurfOrig = Qsurf.copy()
Qsurf[T <= Trange[0]] = None #
fig, ax = plt.subplots(figsize=(8,6), dpi=300)
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

# Plot the image and save
fig, ax = plt.subplots(figsize=(8,6), dpi=300)
im = ax.imshow(imRescaled, cmap=cmaps.get('inferno'), 
               clim=[Trange.min(), 325]) #Trange.max()])

# Set up colorbar properties
# cbar tick marks
tickPosns = np.arange(300, 330, 5) 

# Define colorbar
cbar = fig.colorbar(im, shrink=.4, ticks=tickPosns)

# Set cbar label properties
cbar.set_label('Apparent temp./[K]', rotation=270, labelpad=15)

# Set cbar tick labels
tickPosnLabels = list(map(str, tickPosns))
#    tickPosnLabels[0] = '< ' + tickPosnLabels[0]
tickPosnLabels[-1] = '> ' + tickPosnLabels[-1]
cbar.ax.set_yticklabels(tickPosnLabels)

#plt.show()
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

