#!/usr/bin/env python3
"""Python module of utility functions for findTDE calculations."""
from pathlib import Path
from typing import TYPE_CHECKING, Optional, List, Optional, Tuple, Dict, Any
from collections.abc import Callable

from math import gcd
from fractions import Fraction
import random as rand
import numpy as np

import pymatgen.core as mg
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.util.typing import PathLike
from ase import Atoms
from ase.io.lammpsdata import read_lammps_data, write_lammps_data


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


def binary_to_random_ternary(
    struc_path: PathLike,
    new_struc_path: PathLike,
    to_replace: str = 'Ga',
    replace_with: str = 'Al',
    x_comp: float = 0.5,
    struc_format: str = 'LAMMPS',
    specific_idx: List[int] | None = None,
    rand_seed: int | None = None,
    sort_key: Callable | None = lambda x: x.species.elements[0].symbol.lower()
) -> Tuple[Structure | Atoms, List[int]]:
    """
    Creates a ternary structure by randomly replacing a percentage of a
    given atom type in a binary structure with a third atom type.

    Args
    ---------
        struc_path (PathLike):
            Path to the original binary structure file.
        new_struc_path (PathLike):
            Path at which the new structure file will be created.
        to_replace (str):
            String of an atom symbol in the original structure of
            the type to be replaced. Defaults to 'Ga'.
        replace_with (str):
            String of an atom symbol that will replace atoms in
            the original structure. Defaults to 'Al'.
        x_comp (float):
            Float between 0 and 1 corresponding to the percentage
            of 'replace_with' introduced into the structure. 0 is
            no 'replace_with' in the structure, and 1 replaces all
            of the 'to_replace' atoms with 'replace_with'. Defaults
            to 0.5 (half the atoms are replaced).
        struc_format (str):
            String corresponding to the format of the structure file
            used. Currently supports either LAMMPS data files (specify
            using string starting with 'l' or 'd') or VASP POSCARs
            (specify using string starting with 'v' or 'p'). Defaults
            to 'LAMMPS'.
        specific_idx (list(int) or None):
            List of integers corresponding to specific atom indices
            to be replaced. Defaults to None (random indices).
        rand_seed (int or None):
            Integer to specify the seed for generating a list of
            random indices. Defaults to None (chosen by random).
        sort_key (Callable | None):
            Key for the sorting function used by Python lists/pymatgen.
            Defaults to `lambda x: x.species.elements[0].symbol.lower()`
            (sorts by atomic symbols alphabetically).

    Returns
    ---------
        Tuple of the structure information (either a pymatgen Structure
        or ASE Atoms object) and list of indices replaced.
    """
    if rand_seed is not None:
        rand.seed(rand_seed)
    
    if struc_format.lower()[0] == 'l' or struc_format.lower()[0] == 'd':
        # LAMMPS or Data file
        lmp_data = read_lammps_data(struc_path, atom_style='atomic', units='metal', sort_by_id=True)

        # gather info from structure
        symbols = lmp_data.get_chemical_symbols()
        masses = lmp_data.get_masses()
        atomicnum = lmp_data.get_atomic_numbers()
    elif struc_format.lower()[0] == 'v' or struc_format.lower()[0] == 'p':
        # VASP or POSCAR
        vasp_data = Poscar.from_file(struc_path)
        data_struc = vasp_data.structure

        # gather info from structure
        elements = data_struc.species
        symbols = [i.symbol for i in elements]
        masses = [i.atomic_mass for i in elements]
        atomicnum = [i.Z for i in elements]
    else:
        raise ValueError('Please use a valid structure format (LAMMPS/Data or VASP/POSCAR).')

    other_symbol = [i for i in list(set(symbols)) if i != to_replace][0]
    
    try:
        first_idx_to_replace = symbols.index(to_replace)
    except ValueError:
        raise ValueError(f'No {to_replace} found')

    num_to_replace = symbols.count(to_replace)

    idx_to_replace = []
    if specific_idx is not None:
        idx_to_replace = specific_idx
        for idx in idx_to_replace:
            if symbols[idx] != to_replace:
                raise ValueError(f'Index {idx} is not a {to_replace} atom (found {symbols[idx]})')

    else:
        num_x = int(num_to_replace * x_comp)
        replaced_idx = [i for i, sym in enumerate(symbols) if sym == to_replace]
        idx_to_replace = rand.sample(replaced_idx, num_x)

    idx_to_replace.sort()

    if struc_format.lower()[0] == 'l' or struc_format.lower()[0] == 'd':
        for idx in idx_to_replace:
            symbols[idx] = replace_with
            masses[idx] = mg.Element(replace_with).atomic_mass
            atomicnum[idx] = mg.Element(replace_with).Z
            
        # replace info in structure
        lmp_data.set_chemical_symbols(symbols)
        lmp_data.set_masses(masses)
        lmp_data.set_atomic_numbers(atomicnum)

        for j, val in enumerate(idx_to_replace):
            lmp_data.append(lmp_data[val-j])
            lmp_data.pop(val-j)

        new_symbols = lmp_data.get_chemical_symbols()
    elif struc_format.lower()[0] == 'v' or struc_format.lower()[0] == 'p':
        # replace info in structure
        for idx in idx_to_replace:
            data_struc.replace(idx, replace_with)

        data_struc.sort(key=sort_key)
        
        new_elements = data_struc.species
        new_symbols = [i.symbol for i in new_elements]
    else:
        raise ValueError('Please use a valid structure format (LAMMPS/Data or VASP/POSCAR).')

    print(f'Replaced {len(idx_to_replace)} {to_replace} atoms with {replace_with}')
    print(f'Total atoms: {len(new_symbols)}')
    print(f'{to_replace}: {new_symbols.count(to_replace)}')
    print(f'{replace_with}: {new_symbols.count(replace_with)}')
    print(f'{other_symbol}:  {new_symbols.count(other_symbol)}')

    if struc_format.lower()[0] == 'l' or struc_format.lower()[0] == 'd':
        new_struc = lmp_data
        write_lammps_data(
            new_struc_path,
            new_struc,
            specorder=[to_replace, other_symbol, replace_with],
            masses=True,
            velocities=True,
            units='metal',
            atom_style='atomic'
        )
    elif struc_format.lower()[0] == 'v' or struc_format.lower()[0] == 'p':
        new_struc = Poscar(data_struc)
        new_struc.write_file(new_struc_path)
    else:
        raise ValueError('Please use a valid structure format (LAMMPS/Data or VASP/POSCAR).')
    

    return new_struc, idx_to_replace


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
