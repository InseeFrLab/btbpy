# -*- coding: utf-8 -*-

import sys
import copy
import time
import functools
import multiprocessing as mp
import math as m
import numpy as np
import pandas as pd
import geopandas as gpd

from btbpy.algo_smooth import lissage_loop
from btbpy.createGrid import dfToGrid

from shapely.geometry import Polygon

def smooth_slot(allArgs):
    """
    Call the smooth algorithm in case of parallelization.
    """
    npObsPart, nomsColonnes, iNeighbor, iCellSize, iBandwidth, listVar, nbColCentr, nbRowCentr, xminCentr, yminCentr, mIcentroides, iNbCentroidesUniques = allArgs
    mVariablesLissees = lissage_loop(
        npObsPart[:,nomsColonnes.index('x')]
        , npObsPart[:,nomsColonnes.index('y')]
        , npObsPart[:,nomsColonnes.index('row')].astype(int) + iNeighbor
        , npObsPart[:,nomsColonnes.index('col')].astype(int) + iNeighbor
        , iCellSize
        , iBandwidth
        , iNeighbor
        , npObsPart[:,[nomsColonnes.index(var) for var in listVar]]
        , nbColCentr
        , nbRowCentr
        , xminCentr
        , yminCentr
        , mIcentroides
        , iNbCentroidesUniques
    )
    return(mVariablesLissees)

def kernelSmoothing(pdObservations, sEPSG, iCellSize, iBandwidth, pdCentroids = None, iNbCore = 1, verbose = True):
    """
    Smoothing function with a bisquare kernel - Fonction de lissage à partir d'un noyau bisquare
    
    Parameters
    ----------
    pdObservations: pandas.DataFrame
        Geolocated data with (x,y) cartesian geographical coordinates and variables to smooth
        Tableau de données localisées par des coordonnées géographiques cartésiennes (x,y)  et 
        comportant les variables à lisser
    sEpsg: str
        EPSG code of the projection of geographical data. For example, the RGF93 / Lambert 93
        projection has "2154" epsg code. The WGS 84 projection has the "4356" epsg code.
        Code EPSG du système de projection géographique des coordonnées. Une projection en Lambert 93 
        a le code epsg "2154" et le système de projection WGS 84 a le code epsg "4356".
    iCellSize: int
        Size of grid cells (side's length of squares in the same unit as (x,y) coordinates) - 
        Taille des carreaux de la grille, exprimée dans l'unité des coordonnées géographiques.
    iBandwidth: int
        Radius of the Kernel Density Estimator. This bandwitdh acts as a smoothing parameter, controlling the balance between bias and variance.
        A large bandwidth leads to a very smooth (i.e. high-bias) density distribution. 
        A small bandwith leads to an unsmooth (i.e. high-variance) density distribution.
        It must be the same unit as (x,y) coordinates and iCellSize parameter.
        Rayon de lissage de l'estimation d'intensité par noyau. Cette bande-passante se comporte comme un paramètre de lissage, contrôlant l'équilibre entre biais et variance.
        Un rayon élevé conduit à une densité très lissée, avec un biais élevé? Un rayon faible génère une densité peu lissée avec une forte variance. 
        L'unité est identique à celle des coordonnées géographiques et au paramètre iCellSize.
    pdCentroids: pandas.DataFrame, optional
        Grid with (x,y) coordinates of the centroids to use (Default = None : the centroids are computed by the function).
        The projection has to be the same as the pdObservations' projection.
        DataFrame comportant les coordonnées x et y des centroïdes à utiliser pour le lissage.
        Le système de projection doit être le même que celui des observations (pdObservations).
    iNbCore: int, optional
        Number of cores to use (if > 1, the smoothing algorithm is parallelized - default not parallelized : iNbCore = 1) - 
        Nombre de coeurs du processeur utilisés (si >1, l'algorithme de lissage est parallélisé - par défaut iNbCore = 1).
    verbose: boolean, optional
        print (or not) information about sequences and time during the processing (default to True) - 
        Affiche (ou non) les informations sur les différentes étapes de calcul et les temps d'exécution (par défaut : True)
    
    Return
    ------
    A geopandas.DataFrame with (x,y) coordinates of the centroids, the smoothed variables and a geometry.
    
    Examples
    --------
    
    >>> carburants = pd.read_csv('data/dfPrix_SP95_2016.csv')
    >>> fr_1km = pd.read_csv('data/fr_metro_grid1km.csv') # Grid of 1000 meters side's squares of Metropolitan France
    >>> result = kernelSmoothing(carburants, '2154', 1000, 50000, fr_1km, 4)
        
    """
    iCellSize = int(iCellSize)
    iBandwidth = int(iBandwidth)
    dRayonMinimum = iCellSize * m.sqrt(2) / 2
    
    nomColonnes = list(pdObservations.columns)
    nomColonnes.append('col')
    nomColonnes.append('row')
    iXindex = nomColonnes.index('x')
    iYindex = nomColonnes.index('y')
    iColindex = nomColonnes.index('col')
    iRowindex = nomColonnes.index('row')
    listVar = list(set(nomColonnes) - set(['x','y','col','row']))
    
    npObservations = copy.deepcopy(pdObservations).to_numpy()
    del pdObservations
    
    if (pdCentroids is None):
        iNeighbor = max(0, m.ceil(iBandwidth / iCellSize / 2) - 1)
    else :
        iNeighbor = 0
            
    if(iBandwidth < dRayonMinimum):
        sys.exit("iBandwidth must be greater than iCellSize * sqrt(2) / 2")
    
    if(iCellSize <= 0):
        sys.exit("iCellSize must be greater than 0")
        
    nanData = list(np.sum(np.isnan(npObservations), axis = 1))
    indVar = [ i for i in range(len(nomColonnes)) if i not in [iXindex,iYindex]]
    if(nanData[iXindex] > 0 or nanData[iYindex] > 0):
        sys.exit("NA coordinates are not allowed")
    elif(sum([nanData[i] for i in indVar]) > 0):
        print("Be careful! NA values detected in your observations\n")
    
    # vérifier la regularite des centroides fournis par l'utilisateur
    #l'offset est l'éventuel décalage de la grille en fonction de la taille de carreaux
    if (pdCentroids is not None):
        nomColsCentr = list(pdCentroids.columns)
        iXindexCentr = nomColsCentr.index('x')
        iYindexCentr = nomColsCentr.index('y')
        npCentroids = pdCentroids.copy().to_numpy()
        del pdCentroids
        
        xOffset = (npCentroids[:,iXindexCentr] + iCellSize / 2) % iCellSize
        yOffset = (npCentroids[:,iYindexCentr] + iCellSize / 2) % iCellSize
        if ((not all(xOffset == xOffset[1])) | (not all(yOffset == yOffset[0]))):
            print("Centroids are not regular")
            sys.exit()
        
        obsEtCentroides = np.concatenate((npCentroids[:,[iXindexCentr,iYindexCentr]], npObservations[:,[iXindex,iYindex]]), axis=0)
        #obsEtCentroides = pd.concat([pdCentroids, pdObservations[['x','y']]], axis=0, ignore_index=True, sort=False)
        indiceMinX = np.floor(min(obsEtCentroides[:,0] / iCellSize))
        indiceMinY = np.floor(min(obsEtCentroides[:,1] / iCellSize))

        #col et row correspondent à la coordonnée du centre du carreau dans lequel se trouve l'observation de coordonnées x,y
        colObs = np.floor((npObservations[:,iXindex] - xOffset[0]) / iCellSize - indiceMinX).astype(int)
        rowObs = np.floor((npObservations[:,iYindex] - yOffset[0]) / iCellSize - indiceMinY).astype(int)
        npObservations = np.append(npObservations, colObs.reshape(npObservations.shape[0],1), axis = 1)
        npObservations = np.append(npObservations, rowObs.reshape(npObservations.shape[0],1), axis = 1)
        
        # calcul de l'indice des centroides
        colCentr = np.floor(npCentroids[:,iXindexCentr] / iCellSize - indiceMinX).astype(int)
        rowCentr = np.floor(npCentroids[:,iYindexCentr] / iCellSize - indiceMinY).astype(int)
        
        nbColCentr, nbRowCentr, xminCentr, yminCentr = int(max(colCentr)+1), int(max(rowCentr)+1), int(min(npCentroids[:,iXindexCentr])), int(min(npCentroids[:,iYindexCentr]))
        
        del obsEtCentroides
    else :
        xOffset = 0
        yOffset = 0
        # calcul de l'indice des observations - on prend le rectangle englobant et on positionne le debut de la numérotation sur la première observation
        colObs = (np.floor((npObservations[:,iXindex] - xOffset) / iCellSize) - np.floor(min(npObservations[:,iXindex]/ iCellSize))).astype(int)
        rowObs = (np.floor((npObservations[:,iYindex] - yOffset) / iCellSize) - np.floor(min(npObservations[:,iYindex]/ iCellSize))).astype(int)
        npObservations = np.append(npObservations, colObs.reshape(npObservations.shape[0],1), axis = 1)
        npObservations = np.append(npObservations, rowObs.reshape(npObservations.shape[0],1), axis = 1)
        # calcul des centroides
        npCentroids = np.zeros((rowObs.shape[0], 2))
        npCentroids[:,0] = (np.floor(npObservations[:,iXindex] / iCellSize) * iCellSize + (iCellSize / 2)).astype(int)
        npCentroids[:,1] = (np.floor(npObservations[:,iYindex] / iCellSize) * iCellSize + (iCellSize / 2)).astype(int)
        # les observations sont positionnées sur une matrice. mIndices[col, row] == 1 indique qu'il y a au moins 1 observation pour le carreau (col, row)
        mIndices = np.zeros((max(colObs+1), max(rowObs+1)), dtype = int)
        mIndices[colObs, rowObs] = 1
        # construction d'une matrice des indices des centroides étendue au voisinage
        mIndicesEtendus = np.zeros((max(colObs) + 1 + 2 * iNeighbor, max(rowObs) + 1 + 2 * iNeighbor), dtype = 'int')
        
        for voisin_x in range(-iNeighbor,iNeighbor+1):
            for voisin_y in range(-iNeighbor,iNeighbor+1):
                indiceRows = np.array([iNeighbor + i + voisin_x for i in range(0,mIndices.shape[0])],dtype=np.intp)
                indiceCols = np.array([iNeighbor + j + voisin_y for j in range(0,mIndices.shape[1])],dtype=np.intp)
                mIndicesEtendus[indiceRows[:,np.newaxis], indiceCols] = mIndicesEtendus[indiceRows[:,np.newaxis], indiceCols] + mIndices

        
        indicesNonNuls = np.where(mIndicesEtendus > 0 )
        vIndicesEtendus = np.array([indicesNonNuls[0], indicesNonNuls[1]]).transpose()
        del mIndicesEtendus, mIndices, indicesNonNuls
                  
        #coordonnées             
        npX = (np.round(min(npCentroids[:,0]) + (vIndicesEtendus[:,0] - iNeighbor)*iCellSize)).astype(int)
        npY = (np.round(min(npCentroids[:,1]) + (vIndicesEtendus[:,1] - iNeighbor)*iCellSize)).astype(int)
        npCentroids = np.array([npX, npY]).transpose()
        
        # calcul de l'indice des centroides
        colCentr = vIndicesEtendus[:,0].astype(int)
        rowCentr = vIndicesEtendus[:,1].astype(int)
        nbColCentr, nbRowCentr, xminCentr, yminCentr = int(max(colCentr)+1), int(max(rowCentr)+1), int(min(npCentroids[:,0])), int(min(npCentroids[:,1]))
        del vIndicesEtendus, npX, npY
    
    iNbCentroidesUniques = npCentroids.shape[0]
    index = [i for i in range(npCentroids.shape[0])]    
    mIcentroides = np.full((max(rowCentr)+1, max(colCentr)+1), -1)
    mIcentroides[rowCentr,colCentr] = index
        
    iNbCore = mp.cpu_count() if iNbCore is None else iNbCore
    if verbose :
        print('-------------------------------------------------')
        print(f'number of asked cores : {iNbCore}')
        print('-------------------------------------------------')
        print('start of smoothing loop')
        print('-------------------------------------------------')
        debut = time.time()
    if(iNbCore == 1):
        mVariablesLissees = smooth_slot((npObservations, nomColonnes, iNeighbor, iCellSize, iBandwidth, listVar, nbColCentr, nbRowCentr, xminCentr, yminCentr, mIcentroides, iNbCentroidesUniques))
        del npObservations, mIcentroides
    else:        
        lNpObservations = np.array_split(npObservations, iNbCore, axis = 0)
        del npObservations
        def parallelize(nbproc, func, observ):
            with mp.Pool(iNbCore) as pool:       
                resultat = pool.map(func, [(part, nomColonnes, iNeighbor, iCellSize, iBandwidth, listVar, nbColCentr, nbRowCentr, xminCentr, yminCentr, mIcentroides, iNbCentroidesUniques) for part in observ])
            return(resultat)
        lMatVariablesLissees = parallelize(iNbCore, smooth_slot, lNpObservations)
        mVariablesLissees = np.array(functools.reduce(lambda a,b : a+b, lMatVariablesLissees))
        del lMatVariablesLissees, lNpObservations, mIcentroides
    
    if verbose:
        duree = time.time() - debut
        print('-------------------------------------------------')
        print( f'end of smoothing loop - duration time in seconds = {duree}')
        print('-------------------------------------------------')
    
    mResultat = np.concatenate((npCentroids, mVariablesLissees), axis = 1)
    
    nomsColonnes = ['x','y']
    nomsColonnes.extend(listVar)
    pdResultat = pd.DataFrame(mResultat, columns = nomsColonnes)    
    del mVariablesLissees, mResultat
    if verbose:
        print('-------------------------------------------------')
        print('start of the final grid map cooking')
        print('-------------------------------------------------')
        debut = time.time()
    gpdResultat = dfToGrid(pdResultat, sEPSG = sEPSG, iCellSize = iCellSize, verbose = verbose)
    if verbose:
        duree = time.time() - debut
        print('-------------------------------------------------')
        print( f'end of the final grid map cooking - duration = {duree}')
        print('=================================================')
        print('end of program - enjoy')
        print('=================================================') 
    return(gpdResultat)