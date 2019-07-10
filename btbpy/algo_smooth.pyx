# distutils: language=c++

import cython
from libc.math cimport ceil
import numpy as np
cimport numpy as np
import math as m

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing
@cython.cdivision(True) # Deactivate division by zero checking
cpdef lissage_loop(
    double [:] vXObservation
    , double [:] vYObservation
    , long [:] vLigneObservation
    , long [:] vColonneObservation
    , int iPas
    , int iRayon
    , int iNeighbor
    , double [:,:] mVariables
    , int iNumberCols
    , int iNumberRows
    , int iMinXCentroide
    , int iMinYCentroide
    , long [:,:] mIcentroide
    , int iNbCentroides
):
    """
    Traduction/Adaptation la plus fidèle possible du programme rcppLissage.cpp 
    écrit par Arlindo Dos Santos pour le package sous R.
    """
    cdef double [:] vXObservations = vXObservation
    cdef double [:] vYObservations = vYObservation
    cdef long [:] vLigneObs = vLigneObservation
    cdef long [:] vColonneObs = vColonneObservation
    cdef double [:,:] mVar = mVariables
    cdef long [:,:] mIcentroides = mIcentroide

    cdef int iNbVoisins = int(ceil(iRayon / iPas - 0.5))
    cdef unsigned int iNbCentroidesAbscisse = iNumberCols
    cdef unsigned int iNbCentroidesOrdonnee = iNumberRows
    cdef unsigned int iNbVars = mVar.shape[1]
    cdef unsigned int iNbObs = vXObservations.shape[0]
    cdef unsigned int iRayonCarre = iRayon ** 2    

    cdef double [:,:] mPonderation = np.zeros((iNbCentroidesOrdonnee, iNbCentroidesAbscisse), dtype = 'float64')
    cdef double [:,:] mVariablesLissees = np.zeros((iNbCentroides, iNbVars), dtype = 'float64')
    
    cdef double dSommePonderation
    
    cdef unsigned int iIndiceObsCourante = 0
    cdef unsigned int iCol
    cdef unsigned int iLigne
    cdef unsigned int iVarCourante
    cdef unsigned int iColMin
    cdef unsigned int iColMax
    cdef unsigned int iLigneMin
    cdef unsigned int iLigneMax
    
    cdef double dDistanceCarre

    for iIndiceObsCourante in range(iNbObs):
        dSommePonderation = 0
        iColMin = max(vColonneObs[iIndiceObsCourante] - iNbVoisins, 0)
        iColMax = min(vColonneObs[iIndiceObsCourante] + iNbVoisins, iNbCentroidesAbscisse-1)
        iLigneMin = max(vLigneObs[iIndiceObsCourante] - iNbVoisins, 0)
        iLigneMax = min(vLigneObs[iIndiceObsCourante] + iNbVoisins, iNbCentroidesOrdonnee-1)

        for iCol in range(iColMin,iColMax+1):
            for iLigne  in range(iLigneMin,iLigneMax+1):
                if mIcentroides[iLigne, iCol] == -1:
                    continue
                dDistanceCarre = (vXObservations[iIndiceObsCourante] - (iMinXCentroide + iCol * iPas)) ** 2 + (vYObservations[iIndiceObsCourante] - (iMinYCentroide + iLigne * iPas)) ** 2

                if dDistanceCarre < iRayonCarre:
                    mPonderation[iLigne, iCol] = (1 - (dDistanceCarre / iRayonCarre)) ** 2
                    dSommePonderation += mPonderation[iLigne, iCol]
                else:
                    mPonderation[iLigne, iCol] = 0

        if dSommePonderation > 0:
            for iCol in range(iColMin,iColMax+1):
                for iLigne  in range(iLigneMin,iLigneMax+1):
                    if mIcentroides[iLigne, iCol] == -1:
                        continue

                    for iVarCourante in range(iNbVars):
                        mVariablesLissees[mIcentroides[iLigne, iCol], iVarCourante] += mVar[iIndiceObsCourante, iVarCourante] * mPonderation[iLigne, iCol] / dSommePonderation

                    mPonderation[iLigne, iCol] = 0
    return(np.asarray(mVariablesLissees))