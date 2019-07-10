# -*- coding: utf-8 -*-

import time
import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.geometry import Polygon

def dfToGrid(pdData, sEPSG, iCellSize, verbose):
    epsg = 'epsg:'+sEPSG
    r = iCellSize/2
    if verbose:
        print('\t-------------------------------------------------')
        print('\tstart of making polygons')
        print('\t-------------------------------------------------')
        debut = time.time()
    pdData['geometry'] = [Polygon([(x-r, y+r),
                                 (x+r, y+r),
                                 (x+r, y-r),
                                 (x-r, y-r),
                                 (x-r, y+r)])
                              for x,y in zip(pdData.x,pdData.y)]
    if verbose:
        duree = time.time() - debut
        print('\t-------------------------------------------------')
        print(f'\tend of making polygons - duration in seconds = {duree}')
        print('\t-------------------------------------------------')
        print('\tstart of saving geodataframe')
        print('\t-------------------------------------------------')
        debut = time.time()
    geoResultat = gpd.GeoDataFrame(pdData, crs = {'init' :epsg}, geometry = 'geometry')
    if verbose:
        duree = time.time() - debut
        print('\t-------------------------------------------------')
        print(f'\tend of saving geodataframe - duration in seconds = {duree}')
        print('\t-------------------------------------------------')
    return(geoResultat)