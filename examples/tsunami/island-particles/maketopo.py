
"""
Module to create topo and qinit data files for this example.
"""

import clawpack.geoclaw.topotools as topotools
import numpy as np

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 201
    nypoints = 241
    xlower = 0.e0
    xupper = 100.e0
    ylower = 0.e0
    yupper = 50.e0
    outfile= "island.tt3"     

    topography = topotools.Topography(topo_func=topo)
    topography.x = np.linspace(xlower,xupper,nxpoints)
    topography.y = np.linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=3, Z_format="%22.15e")

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 101
    nypoints = 101
    xlower = -50.e0
    xupper = 50.e0
    yupper = 50.e0
    ylower = -50.e0
    outfile= "qinit.xyz"     

    topography = topotools.Topography(topo_func=qinit)
    topography.x = np.linspace(xlower,xupper,nxpoints)
    topography.y = np.linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=1)

def topo(x,y):
    """
    Island
    """
    ze = -((x-40.)**2 + (y-35.)**2)/20.
    
    #z_island = where(ze>-10., 100.*exp(ze), 0.)
    z_island = np.where(ze>-10., 150.*np.exp(ze), 0.)
    z = -50 + z_island
    return z


def qinit(x,y):
    """
    Dam break
    """
    return np.where(x<10, 40., 0.)

if __name__=='__main__':
    maketopo()
    makeqinit()
