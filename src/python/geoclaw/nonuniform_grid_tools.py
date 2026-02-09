"""
Some tools for making 1D grids.

make_mapc2p - based on a file celledges.data with cell edges, create a
    mapc2p file making 0 <= xc <= 1 to physical cell edges.

make_celledges_cfl - create a celledges.data file with celledges chosen so
    that the Courant number is nearly 1 in each cell in the ocean.

make_pwlin_topo_fcn - take a list of (x,z) pairs and create a function that
    is the piecewise linear interpolant through these points.

"""

import numpy as np
import scipy.interpolate
import os


def make_mapc2p(fname_celledges='celledges.data'):
    """
    Create a mapc2p function that maps computational cell edges xc
    with 0 <= xc <= 1 to the physical cell edges.  The physical
    cell edges should be in the file fname_celledges, 
    starting in the second row (following mx_edge, the number of edges).

    Returns the mapc2p function and also mx_edge, xp_edge, which may be
    needed for other purposes.
    """

    path = os.path.abspath(fname_celledges)
    d = np.loadtxt(path, skiprows=1) 
    mx_edge = d.shape[0]

    print('make_mapc2p: Using %i cell edge values from %s' % (d.shape[0], path))

    xc_edge = np.linspace(0,1,mx_edge)  # assumes xlower=0, xupper=1 in setrun
    xp_edge = d[:,0]
    mapc2p = scipy.interpolate.interp1d(xc_edge, xp_edge, kind='linear')

    return mapc2p, mx_edge, xp_edge


def make_pwlin_topo_fcn(xzpairs):
    """
    Input: xzpairs should be a list of tuples (xi,zi).
    Output: a function that is the piecwise linear interpolant.
    """

    xi = np.array([xz[0] for xz in xzpairs])
    zi = np.array([xz[1] for xz in xzpairs])
    z_fcn = scipy.interpolate.interp1d(xi, zi, kind='linear', bounds_error=False, 
                     fill_value='extrapolate')
    return z_fcn


def make_celledges_cfl(xlower, xupper, mx, topo_fcn, hmin,
                       fname='celledges.data', plot_topo=False):

    grav = 9.81

    cmin = np.sqrt(grav*hmin)

    def c(x):
        z = topo_fcn(x)
        h = np.where(-z > hmin, -z, hmin)
        c = np.sqrt(grav*h)
        return c

    xunif = np.linspace(xlower, xupper, 2*mx)
    cunif = c(xunif)
    csum = np.cumsum(1./cunif)
    csum = csum - csum[0]

    csum = csum / csum[-1]
    cinv = scipy.interpolate.interp1d(csum, xunif)

    xc = np.linspace(0, 1, mx+1)   # computational grid
    xp = cinv(xc)
    z = topo_fcn(xp)
    dxp = np.diff(xp)
    
    if plot_topo:
        import matplotlib.pyplot as plt
        plt.figure(97, figsize=(6,8))
        plt.clf()
        plt.subplot(311)
        #plot(csum, xunif, 'b')
        plt.plot(xunif, csum, 'b')
        plt.ylabel('computational coordinate xc')
        plt.grid(True)
        plt.axis([xlower,xupper,0,1])
        plt.title('inverse of mapc2p function')

        plt.subplot(312)
        xcell = 0.5*(xp[1:] + xp[:-1])
        plt.plot(xcell, dxp, 'b')
        #xlabel('physical coordinate xp')
        plt.xlabel('physical coordinate xp')
        plt.grid(True)
        plt.xlim(xlower,xupper)
        plt.title('Mesh width')
        
        plt.subplot(313)
        #xcell = 0.5*(xp[1:] + xp[:-1])
        dxratio = dxp[1:]/dxp[:-1]
        print('dx ratio between adjacent cells varies between %.4f and %.4f' \
            % (dxratio.min(), dxratio.max()))
        plt.plot(xcell[1:], dxratio, 'b')
        plt.xlabel('physical coordinate xp')
        plt.ylabel('delta x ratio')
        plt.grid(True)
        plt.xlim(xlower,xupper)
        plt.title('Mesh width ratio between adjacent cells')
        plt.tight_layout()
        
        png_fname = 'cellmap.png'
        plt.savefig(png_fname)
        print("Created ",png_fname)

    with open(fname,'w') as f:
        f.write('%i   # number of cell edges\n' % (mx+1))

        for i in range(mx+1):
            f.write('%15.8f %15.8f\n' % (xp[i],z[i]))
        f.close()

    print("Created %s, containing %i cell edges" % (fname,mx+1))
    print("Min dx = %g, Max dx = %g" % (dxp.min(),dxp.max()))

    if plot_topo:
        import matplotlib.pyplot as plt
        plt.figure(99, figsize=(8,4))
        plt.clf()
        plt.fill_between(xp,np.where(z<0,z,np.nan),0.,color=[.5,.5,1])
        plt.plot(xp,z,'g')
        plt.xlim(xlower,xupper)
        zmax = np.max(z.max(), 0)
        zmargin = 0.1*(zmax-z.min())
        plt.ylim(z.min()-zmargin,zmax+zmargin)
        plt.grid(True)
        plt.title('Topography')
        png_fname = 'topo.png'
        plt.savefig(png_fname)
        print("Created ",png_fname)

    return xp,z
