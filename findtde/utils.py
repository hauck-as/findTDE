#!/usr/bin/env python3
"""Python module of utility functions for findTDE calculations."""
from math import gcd
from fractions import Fraction
import numpy as np


# math functions
def unit_vector(vector):
    # https://stackoverflow.com/a/13849249
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    # https://stackoverflow.com/a/13849249
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))*(180/np.pi)


def cart2sph(x, y, z):
    dxy = np.sqrt(x**2 + y**2)
    r = np.sqrt(dxy**2 + z**2)
    theta = np.arctan2(y, x)
    phi = np.arctan2(dxy, z)
    theta, phi = np.rad2deg([theta, phi])
    return r, theta % 360, phi


def sph2cart(theta, phi, r=0.5):
    theta, phi = np.deg2rad([theta, phi])
    z = r * np.cos(phi)
    rsinphi = r * np.sin(phi)
    x = rsinphi * np.cos(theta)
    y = rsinphi * np.sin(theta)
    return x, y, z


# crystallography functions
# https://ssd.phys.strath.ac.uk/resources/crystallography/crystallographic-direction-calculator/
def lattice_iconv(indices, ind_type='UVTW'):
    """
    Function converting between rectangular and hexagonal lattice indices. Given a
    tuple of indices and the initial index type as a keyword argument. Index type
    is either [uvw] for rectangular or [UVTW] for hexagonal. Returns a tuple of
    the opposite type.
    """
    print(indices)
    
    if ind_type == 'UVTW':
        U, V, T, W = indices
        u = 2*U + V
        v = 2*V + U
        w = W
        n = gcd(u, gcd(v, w))
        
        u /= n
        v /= n
        w /= n
        
        new_indices = (int(u), int(v), int(w))
    
    elif ind_type == 'uvw':
        u, v, w = indices
        U = Fraction(((2*u) - v), 3)
        V = Fraction(((2*v) - u), 3)
        T = -1*(U + V)
        W = w
        
        denom_gcd = gcd(U.denominator, gcd(V.denominator, T.denominator))
        U *= denom_gcd
        V *= denom_gcd
        T *= denom_gcd
        W *= denom_gcd
        
        new_indices = (int(U), int(V), int(T), int(W))
    
    return new_indices


def lat2sph(uvw, ai):
    """
    Function converting from lattice directions to spherical coordinates. Given two 2D
    arrays, one for lattice directions and one for lattice vectors, returns a 2D array
    of the spherical coordinates corresponding to the lattice directions.
    """
    xyz, rpt = np.zeros(uvw.shape), np.zeros(uvw.shape)
    
    for i in range(uvw.shape[0]):
        xyz[i, :] = uvw[i, :]@ai
        rpt[i, :] = np.around(cart2sph(xyz[i, 0], xyz[i, 1], xyz[i, 2]), decimals=2)
        rpt[i, 0] = 1.00
    
    return rpt


def sph2lat(rpt, ai):
    """
    Function converting from lattice directions to spherical coordinates. Given two 2D
    arrays, one for lattice directions and one for lattice vectors, returns a 2D array
    of the spherical coordinates corresponding to the lattice directions.
    """
    xyz, uvw = np.zeros(rpt.shape), np.zeros(rpt.shape)
    
    for i in range(rpt.shape[0]):
        xyz[i, :] = sph2cart(rpt[i, 1], rpt[i, 2], r=rpt[i, 0])
        uvw[i, :] = xyz[i, :]@np.linalg.inv(ai)
        n = gcd(round(uvw[i, 0]), gcd(round(uvw[i, 1]), round(uvw[i, 2])))
        uvw[i, :] /= n
    
    return uvw


def scale_C3z_symmetry(rpt):
    """
    Function to scale all given spherical coordinates to a single zone based on threefold
    rotation symmetry about the c-axis. Given a 2D array of spherical coordinates, returns
    another 2D array of spherical coordinates with polar angle scaled by 120 deg increments
    to the desired range.
    """
    for i in range(rpt.shape[0]):
        if rpt[i, 1] >= 30. and rpt[i, 1] <= 150.:
            pass
        elif rpt[i, 1] >= 0. and rpt[i, 1] < 30.:
            rpt[i, 1] += 120.
        elif rpt[i, 1] > 150. and rpt[i, 1] <= 270.:
            rpt[i, 1] -= 120.
        elif rpt[i, 1] > 270. and rpt[i, 1] <= 360.:
            rpt[i, 1] -= 240.
        else:
            raise Exception('Angle outside of 0-360 deg')
    return rpt


def scale_C6z_symmetry(rpt):
    """
    Function to scale all given spherical coordinates to a single zone based on sixfold
    rotation symmetry about the c-axis. Given a 2D array of spherical coordinates, returns
    another 2D array of spherical coordinates with polar angle scaled by 60 deg increments
    to the desired range.
    """
    rpt_C3z = scale_C3z_symmetry(rpt)
    
    for i in range(rpt_C3z.shape[0]):
        if rpt_C3z[i, 1] >= 30. and rpt_C3z[i, 1] <= 90.:
            pass
        elif rpt_C3z[i, 1] > 90. and rpt_C3z[i, 1] <= 150.:
            rpt[i, 1] = 180. - rpt[i, 1]
        else:
            raise Exception('Angle outside of 30-150 deg')
    return rpt


# interpolation functions
def idw(samples, tx, ty, P=5):
    """
    Function to compute a single IDW interpolation of sample data.
    Change P value for different interpolation (P>2 is recommended).
    """
    def dist(a, b):
        return ((a[0]-b[0])**2 + (a[1]-b[1])**2)**(1./2.)

    num = 0.
    den = 0.
    for i in range(0, len(samples)):
        d = (dist(samples[i], [tx, ty]))**P
        if(d < 1.e-5):
            return samples[i][2]

        w = 1/d
        num += w*samples[i][2]
        den += w

    return num/den


def idw_heatmap(inputdata, RES=360, P=5):
    """
    Perform full IDW interpolation on a (nsamples x 3) array (x, y, f(x, y))
    and produce data able to be used in a heatmap. From Victor.
    Change P value for different interpolation (P>2 is recommended).
    RES is the resolution of the final image.
    """
    # from Victor, plot heatmap
    # Setup z as a function of interpolated x, y
    minx, maxx = min([d[0] for d in inputdata]), max([d[0] for d in inputdata])
    miny, maxy = min([d[1] for d in inputdata]), max([d[1] for d in inputdata])
    minz, maxz = min([d[2] for d in inputdata]), max([d[2] for d in inputdata])    # useful to scale 0 < z < 1
    dx, dy = (maxx - minx)/(RES-1), (maxy - miny)/(RES-1)
    xs = [minx + i*dx for i in range(0, RES)]
    ys = [miny + i*dy for i in range(0, RES)]
    zs = [[None for i in range(0, RES)] for j in range(0, RES)]
    for i in range(0, RES):
        for j in range(0, RES):
            zs[i][j] = idw(samples=inputdata, tx=xs[i], ty=ys[j], P=P)
    # print(xs, ys, zs)
    return (xs, ys, np.transpose(zs))
