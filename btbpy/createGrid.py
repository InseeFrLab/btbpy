# -*- coding: utf-8 -*-

import time
import numpy as np
import pandas as pd
import geopandas as gpd

from shapely.geometry import Polygon

def dfToGrid(pdData, sEPSG, iCellSize, verbose):
    """
    
    Compute a grid from a pandas.DataFrame
    Génère une grille à partir d'un pandas.DataFrame
    
    Parameters
    ----------
    pdObservations: pandas.DataFrame
        pandas.DataFrame with (x,y) coordinates of the centroids of the squares to draw.
        pandas.DataFrame comportant les coordonnées (x,y) des centroïdes des cellules à dessiner.
    sEpsg: str
        EPSG code of the grid returned by the function. For example, the RGF93 / Lambert 93
        projection has "2154" epsg code. The WGS 84 projection has the "4356" epsg code.
        Code EPSG du système de projection de la grille retournée par la fonction. Une projection en Lambert 93 
        a le code epsg "2154" et le système de projection WGS 84 a le code epsg "4356".
    iCellSize: int
        Size of grid cells (side's length of squares in the same unit as (x,y) coordinates) - 
        Taille des carreaux de la grille, exprimée dans l'unité des coordonnées géographiques.
    verbose: boolean, optional
        Print (or not) information about sequences and time during the processing (default to True) - 
        Affiche (ou non) les informations sur les différentes étapes de calcul et les temps d'exécution (par défaut : True)
    
    Return
    ------
    A geopandas.DataFrame made of polygons (the squares of the grid).
    Un geopandas.DataFrame composé de polygônes (les carrés de la grille).
    
    Example
    -------
    
    >>> df = pd.DataFrame({'x': pd.Series([100,100,300,300,500]),
                       'y': pd.Series([100,300,100,300,100])})
    >>> dfResult = dfToGrid(df, '2154', 200, True)
    """
    
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