#!/usr/bin/python

"""
fixdata
=========
    Provides functions for removing erroneous data from arrays as grids

    Contains:
        findbadindices
        fillbaddata


"""

import numpy as np

from . import iotools

#==============================================================================
def findbadindices (Z,badvalue=inf,removenans=True):
    """
        remove nans or infs from an array
    """

    badind=[]
    for i in range(Z.shape[0]) :
        for j in range(Z.shape[1]):

            if Z[i,j]==inf:
                badind.append((i,j))
            elif Z[i,j]==badvalue:
                badind.append((i,j))

            if Z[i,j]!=Z[i,j] and removenans:
                badind.append((i,j))


    return badind

#===============================================================================
def fillbaddata (Z,badinds):

    """
    fill data in array Z, at indice tuples in list badinds
    by averaging surrounding good data.
    return new array.
    """
    
    m=Z.shape[0]
    n=Z.shape[1]

    for ind in badinds :
        i=ind[0]
        j=ind[1]
        r=0
        indbad=True
        while indbad and r < np.max(m,n):
            r=r+1 #radius of ball around badinds in inf-norm (square)
            irange=list(np.range(np.max(0,i-r),np.min(i+r+1,m)))
            jrange=list(np.range(np.max(0,j-r),np.min(j+r+1,n)))
            summands=0
            sum=0.
            for ii in irange:
                for jj in jrange:
                    ballind=(ii,jj)
                    if not ballind in badinds:
                        sum = sum + Z[ballind[0],ballind[1]]
                        summands=summands+1
            if summands >0 : 
                Z[ind[0],ind[1]] = sum/summands
                indbad=False

    return Z


#=====================================================================================
def filterdata (Z,filterinds,radius=1):
    """
    filter data in array z, at indice tuples in list filterinds
    by averaging surrounding data, ball with radius=radius in inf-norm
    acts as a low-band pass filter and removes oscillatory data
    """

    m=Z.shape[0]
    n=Z.shape[1]

    for ind in filterinds :
        i=ind[0]
        j=ind[1]
        r=radius

        irange=list(np.range(np.max(0,i-r),np.min(i+r+1,m)))
        jrange=list(np.range(np.max(0,j-r),np.min(j+r+1,n)))
        summands=0
        sum=0.
        for ii in irange:
            for jj in jrange:
                ballind=(ii,jj)
                sum = sum + Z[ballind[0],ballind[1]]
                summands=summands+1
        if summands >0 : 
            Z[ind[0],ind[1]] = sum/summands

    return Z
