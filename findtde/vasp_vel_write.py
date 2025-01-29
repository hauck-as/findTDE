#!/usr/bin/env python3
"""Python module used to read lattice direction and KE text files and generate velocity vectors for VASP."""
import os
import subprocess
import sys
import re
import math as m
from pymatgen.io.vasp import Poscar, Potcar
from pymatgen.io.lammps.inputs import LammpsInputFile
import numpy as np

atom_type, atom_num, lat_dir_pseudo, sim_prog = str(sys.argv[1]), int(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4])

# set directory variables, using PWD as the base directory
base_path = os.getcwd()
bin_path, inp_path, perfect_path = os.path.relpath('bin', base_path), os.path.relpath('inp', base_path), os.path.relpath('perfect', base_path)
if sim_prog == 'vasp':
    tde_run_path = os.path.relpath(lat_dir_pseudo + '_' + atom_type + str(atom_num), base_path)
    # set atomic weights from POTCAR
    mass_lines = subprocess.run(['grep', 'POMASS', os.path.join(inp_path, 'POTCAR')], capture_output = True, text = True)
    mass_lines_list = mass_lines.stdout.split('\n ')
    atomic_masses = []
    for i in range(len(mass_lines_list)):
        atomic_masses.append(float(re.search(r'\d+.\d+', mass_lines_list[i]).group()))
elif sim_prog == 'lammps':
    tde_run_path = os.path.relpath(lat_dir_pseudo + '_' + atom_type + str(atom_num) + '_lmp', base_path)
    # set atomic weights from LAMMPS input file
    lmp_inp_f = LammpsInputFile.from_file(os.path.join(inp_path, 'input.tde'))
    mass_list, group_list = lmp_inp_f.get_args('mass'), lmp_inp_f.get_args('group')
    atom_mass_dict = {}
    for i in group_list:
        if 'type' in i:
            atom_id_list = i.split(' type ')
            for j in mass_list:
                id_mass_list = j.split(' ')
                if atom_id_list[1] == id_mass_list[0]:
                    atom_mass_dict.update({atom_id_list[0]: float(id_mass_list[1])})
    atomic_masses = list(atom_mass_dict.values())
tde_run_outfile_path = os.path.relpath(os.path.join(tde_run_path, lat_dir_pseudo + '_' + atom_type + str(atom_num) + '_out.txt'), base_path)


# utility functions
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


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


def vasp_inp_match_check(poscar_filepath, potcar_filepath):
    """
    Function to compare the order of elements in the POSCAR file to the order of
    elements in the POTCAR file and ensure the order matches.
    """
    # read in POSCAR and POTCAR using pymatgen
    pos_file, pot_file = Poscar.from_file(poscar_filepath), Potcar.from_file(potcar_filepath)

    # get list of species in POSCAR and POTCAR files
    pos_spec_list, pot_spec_list = pos_file.site_symbols, pot_file.symbols

    # compare lists
    for ele_idx in range(len(pos_spec_list)):
        pos_ele, pot_ele = pos_spec_list[ele_idx], pot_spec_list[ele_idx]
        pos_ele = pos_ele.lower()
        pot_ele = pot_ele.lower().split('_')[0]
        if pos_ele == pot_ele:
            pos_pot_ele_match = True
        elif pos_ele != pot_ele:
            pos_pot_ele_match = False
            break
        else:
            print('Error regarding POSCAR/POTCAR species')
    
    return pos_pot_ele_match


# velocity vector calculation function
def lat_vel_calc(M, E, dir):
    """
    Function to calculate velocity vector of an atom in the lattice.
    Given an atomic mass (u), energy (eV), and crystallographic direction,
    returns the velocity vector for the lattice atom (Ang./fs).
    The direction and velocity are 1D numpy arrays.
    """
    M *= 931.49432*1e6    # 1/u * MeV/c^2 * eV/MeV
    
    v_mag = np.sqrt(np.divide(2*E, M))
    v_mag *= 299792458*1e10*1e-15    # m/s * Angstrom/m * s/fs
    
    v = np.zeros((v_mag.shape[0], dir.shape[0]))
    
    dir_mag = np.linalg.norm(dir)
    dir_unit_vec = dir/dir_mag
    
    for i in range(v_mag.shape[0]):
        v[i, :] = np.multiply(v_mag[i], dir_unit_vec)
    return v


# writing formatted velocity vector to POSCAR file
def vel_to_POSCAR(vel_vec, atom_type, atom_no, filepath):
    """
    Function to write velocity vector of an atom in the lattice to the POSCAR file.
    Given a velocity vector (Ang./fs), atom type, and the number of the atom
    (from the POSCAR file). The velocity vector is a 1D numpy array, the atom
    type is a pymatgen Element, and the atom number is an integer.
    """
    # read in POSCAR with pymatgen
    pos_file = Poscar.from_file(filepath)
    spec_list = pos_file.site_symbols
    no_list = pos_file.natoms
    total_ions = sum(no_list)

    # set velocity list to array of zeros for all atoms
    vel_arr_list = [np.zeros(3) for i in range(total_ions)]
    
    # replace velocity line for desired atom with converted velocity vector
    for i in range(len(spec_list)):
        if atom_type.lower() in spec_list[i].lower():
            vel_arr_list[atom_no+sum(no_list[:i])-1] = vel_vec

    # set velocities in POSCAR to velocity array
    pos_file.velocities = vel_arr_list
    
    pos_file.write_file(filepath)


# read lattice directions and kinetic energies from text files
# finds lattice direction from pseudo passed from bash
lat_dir_list_file = os.path.relpath('latt_dirs_to_calc.csv', base_path)
with open(lat_dir_list_file) as ld_f:
    for line in ld_f:
        if line[:len(lat_dir_pseudo)] == lat_dir_pseudo:
            lat_dir_list = re.split(r'[;,\s]\s*', line.strip('\n'))  # [u v w] or (rho, phi, theta)
            break

kinE_list_file = os.path.relpath(os.path.join(tde_run_path, 'KE_calcs_list.txt'), base_path)
with open(kinE_list_file) as ke_f:
    for line in ke_f:
        pass
    kinE_line = line  # eV

kinE = np.array([float(kinE_line)])
print('KE:', kinE, 'eV')

tde_kinE_path = os.path.relpath(os.path.join(tde_run_path, kinE_line.strip('\n') + 'eV'), base_path)
pos_file = os.path.join(tde_kinE_path, 'POSCAR')
pot_file = os.path.join(tde_kinE_path, 'POTCAR')

# check if POTCAR and POSCAR are congruent
if vasp_inp_match_check(pos_file, pot_file) == False:
    raise ValueError

# read in POSCAR file using pymatgen
pos = Poscar.from_file(pos_file)
lat_vecs = pos.structure.lattice.matrix
pos_spec_list = pos.site_symbols

if lat_dir_pseudo[-1] == 'L':
    lat_dir = np.array(lat_dir_list[5:], dtype=int)
    vel_dir = lat_dir@lat_vecs
    print('lattice direction:', str(lat_dir))
elif lat_dir_pseudo[-1] == 'S':
    lat_dir = np.array(lat_dir_list[5:], dtype=float)
    # rho can be any value since magnitude is scaled, choose 1 for simplicity
    # first angle (phi) needs to be the angle given in the heatmap - 60 deg since heatmap adds 60 deg
    cart_coords = np.array(sph2cart(float(lat_dir_list[5:][1]), float(lat_dir_list[5:][2]), r=1.))
    vel_dir = cart_coords
    print('spherical direction:', str(lat_dir))

# calculate velocity vector and write to POSCAR file
for i in range(len(pos_spec_list)):
    if atom_type.lower() == pos_spec_list[i].lower():
        vel_vecs = lat_vel_calc(atomic_masses[i], kinE, vel_dir)
vel_to_POSCAR(vel_vecs[0], atom_type, atom_num, pos_file)
print('velocity vectors:', vel_vecs)

out_file = open(tde_run_outfile_path, 'a')
out_file.write('direction:  '+str(lat_dir)+'\nvelocity vectors:  '+str(vel_vecs)+'\n')
out_file.close()
