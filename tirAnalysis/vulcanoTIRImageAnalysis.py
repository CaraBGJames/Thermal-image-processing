""" vulcanoTIRImageAnalysis.py

Contents of
~/FieldWork/Vulcano/2016-09-Vulcano/2016-09-17_nighttimeFLIR/allGeorefScenes/emissivityCorrected/

mergeScenes.gdal
SITE01_T_georef.tif
SITE01_T_georef.tifw
SITE02_T_georef.tif
SITE02_T_georef.tifw
SITE05_T_georef.tif
SITE05_T_georef.tifw
SITE07_T_georef.tif
SITE07_T_georef.tifw
SITE10_T_georef.tif
SITE10_T_georef.tifw
SITE11_T_georef.tif
SITE11_T_georef.tifw

Each of the *tif files above is a geotif (containing georeferencing data)
"""

# from rasterio.errors import RasterioIOError
# from rasterio.plot import show
from rasterio import features
from scipy.stats import gamma
from random import sample

import numpy as np
import rasterio as rio
import geopandas as gpd
# import pandas as pd
# import matplotlib.pyplot as plt
import fiona
import os


def affine_transform(filename):
    """
    """
    from rasterio.transform import Affine
    
    with open(filename, 'r') as f: 
        transform = np.float64(f.read().split())
    transform[2], transform[4] = transform[4], transform[2]
    return Affine(*transform)


def gamma_dist_fit(data, axis=None):
    """
    Returns MLE of scale and location parameters for a gamma distribution.  
    If optional parameter "axis" is not None, the distribution is plotted on 
    the specified figure axis.  If the supplied data is too large (>1e6 data
    points), it is subsampled.

    Parameters
    ----------
    data: array_like
        "one"-dimensional array of the random variable samples
    axis: matplotlib axis

    Returns
    -------
    data: array_like
    popt: tuple
        MLE of the distribution parameters (
    """
    try:
        data = data[data.mask]
    except AttributeError:
        pass

    if len(data) > 1e6:
        data = sample(data, 1e6)
        
    popt = gamma.fit(data)

    return data, popt


# plt.close('all')

# fig0, ax0 = plt.subplots(figsize=(8, 8))
# fig1, ax1 = plt.subplots(figsize=(8, 4))
# fig2, ax2 = plt.subplots(figsize=(8, 4))

BASE = '/home/david/FieldWork/Vulcano/2016-09-Vulcano/' + \
       '2016-09-17_nighttimeFLIR/'

# Open orthophotos and DEM
# '/home/david/FieldWork/Vulcano/DEM/' +
OrthoImLayer = rio.open(
    BASE + '../../DEM/06_final_20cm_model_ortho/rgb_20cm_hillshaded.jpg', 'r')
Oim = OrthoImLayer.read()

DEM = rio.open(
    BASE + '../../DEM/06_final_20cm_model_ortho/dem_20cm_ortho.flt', 'r')
DEM_mask       = DEM.read_masks()
DEM_transform  = DEM.transform
DEM_elevations = DEM.read()[0, :, :]
# np.ma.masked_array(DEM_mask, DEM.read()[0, :, :])

# GeoDataFrame of the survey locations and conditions at each site
gdf = gpd.read_file(BASE + 'vulcanoShapefiles/SITES.shp')
# gdf.crs = {'init': 'epsg:32633'}  

# Open shapefiles for the vulcano features, converting to UTM zone 33N
craterFloor = gpd.read_file(
    BASE + 'vulcanoShapefiles/craterFloor.shp').to_crs(epsg=32633)
lowerRidge  = gpd.read_file(
    BASE + 'vulcanoShapefiles/lowerRidge.shp').to_crs(epsg=32633) 
middleRidge = gpd.read_file(
    BASE + 'vulcanoShapefiles/middleRidge.shp').to_crs(epsg=32633)
rimRidge = gpd.read_file(
    BASE + 'vulcanoShapefiles/rimRidge.shp').to_crs(epsg=32633) 


# Plot orthophoto or dem
#show(Oim, transform=OrthoImLayer.transform, ax=ax0, zorder=1)
# show(DEM_elevations, transform=DEM_transform, contours=True, ax=ax0, 
#      cmap=plt.cm.gray, zorder=1)
#ax0.get_images()[-1].set_clim(190, 390)

# A dict of the layers used, plus positional argument for zorder
layer_orders = {'1' : 4,
                '2' : 1,
                '5' : 5,
                '7' : 3,
                '10': 2,
                '11': 1}

histInfo = []
bins = np.linspace(295, 395, 251)  # np.append(0, )

for i in layer_orders:
    i = int(i)
    PATH = BASE + 'SITE%02d/georeferenced/' % i
    filename  = PATH + 'SITE%02d_T_georef.tif' % i
    transform_filename = filename + 'w'

    image     = rio.open(filename, crs={'init' : 'epsg:32633'})
    #image.crs = {'init' : 'epsg:32633'}
    T         = np.float32(image.read(1))  # float64 too large for shapely
    mask      = T != 0  # Interior of scenes
    Tm        = np.ma.masked_where(T == 0, T)
    transform = image.transform  # affine_transform(transform_filename)
    with rio.open(BASE + 'allGeorefScenes/rasterio_test.tif', 'w',
                  driver='GTiff',
                  height=T.shape[0], width=T.shape[1],
                  count=1, dtype=T.dtype,
                  crs={'init' : 'epsg:32633'},
                  transform=transform) as dst:
        dst.write(T, 1)

    print(filename, '\n', image.meta, '\n')
    # Extract boundaries of visible raster images to shapefiles, with
    # a check to see if file already exists
    SHAPEFILE = PATH + 'SITE%02d.shp' % i
    if not os.path.isfile(SHAPEFILE):
        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for i, (s, v) in enumerate(
                    features.shapes(mask.astype('int16'), mask=mask,
                                    transform=transform)))
        print('Writing shapefile for SITE%02d...' % i)
        with fiona.open(SHAPEFILE, 'w',
                        driver='Shapefile', crs={'init' : 'epsg:32633'},
                        schema={'properties': [('raster_val', 'int')],
                                'geometry': 'Polygon'}) as dst:
            dst.writerecords(results)
        print('...done!')

    geom = list(
        {'properties': {'raster_val': v}, 'geometry': s}
        for i, (s, v) in enumerate(
                features.shapes(Tm, mask=Tm.mask, transform=transform)))
    polygon = gpd.GeoDataFrame.from_features(geom)
    polygon.to_file('SITE%02d/georeferenced/SITE%02d.shp' % (i, i))
    #polygon.plot(ax=ax0, linestyle='-')

    # popt = gamma.fit(T[T > 295].ravel())  # alpha, loc, scale
    # binCentres = .5 * (bins[1:] + bins[:-1])
    # histInfo.append(np.histogram(Tm.ravel(), bins=bins, density=True))
    #                 #alpha=.75, label='SITE%02d' % i))
    # np.save(PATH + 'T', T)
    # np.save(PATH + 'mask', mask)
    # np.save(PATH + 'popt', popt)


# Add features such as the measurement sites, crater rim etc
# gdf.plot(ax=ax0, marker='v', markersize=80, c='C1', 
#          edgecolor='k', linewidth=.5, zorder=7, label='SITES')
# rimRidge.plot(
#     ax=ax0, edgecolor='C3', facecolor='', zorder=6, label='Rim ridge')
# craterFloor.plot(
#     ax=ax0, edgecolor='C0', facecolor='', zorder=6, label='Crater floor')
# craterFloor.centroid.plot(
#     ax=ax0, marker='.', markersize=80, color='r', zorder=6)

# ax0.ticklabel_format(useOffset=False)
# ax0.set_ylim(OrthoImLayer.bounds[1:4:2])        

# for im in ax0.get_images(): 
#     im.set_clim(294, 365)

# ax0.set_title('WGS84 UTM Zone 33N projection')
# ax0.set_xlabel('Eastings/[m]')
# ax0.set_ylabel('Northings/[m]')
    
# ax1.set_xlabel('(Uncorrected) pixel-averaged temperature/[K]')
# ax1.set_ylabel('Pixel counts/[-]')
# ax1.set_yscale('log')
# ax1.legend()

# fig0.tight_layout()
# fig0.savefig('/home/david/Dropbox/articles/vulcanoThermalSurvey2016/images/'
#              + 'vulcanoThermalSurvey.png', dpi=300)

# fig1.tight_layout()
# fig1.savefig(('/home/david/Dropbox/articles/vulcanoThermalSurvey2016/images/'
#               + 'vulcanoThermalSurvey_siteHistograms.png'), dpi=300)

# plt.show()
