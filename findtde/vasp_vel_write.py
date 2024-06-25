#!/usr/bin/env python3
"""Python module used to read lattice direction and KE text files and generate velocity vectors for VASP."""
import os
import subprocess
import sys
import re
import math as m
import fortranformat as ff
import numpy as np

atom_type, atom_num, lat_dir_pseudo, sim_prog = str(sys.argv[1]), int(sys.argv[2]), str(sys.argv[3]), str(sys.argv[4])

# set directory variables, using PWD as the base directory
base_path = os.getcwd()
bin_path, inp_path, perfect_path = os.path.relpath('bin', base_path), os.path.relpath('inp', base_path), os.path.relpath('perfect', base_path)
if sim_prog == 'vasp':
    tde_run_path = os.path.relpath(lat_dir_pseudo + '_' + atom_type + str(atom_num), base_path)
elif sim_prog == 'lammps':
    tde_run_path = os.path.relpath(lat_dir_pseudo + '_' + atom_type + str(atom_num) + '_lmp', base_path)
tde_run_outfile_path = os.path.relpath(os.path.join(tde_run_path, lat_dir_pseudo + '_' + atom_type + str(atom_num) + '_out.txt'), base_path)

# set atomic weights from POTCAR
mass_lines = subprocess.run(['grep', 'POMASS', os.path.join(inp_path, 'POTCAR')], capture_output = True, text = True)
mass_lines_list = mass_lines.stdout.split('\n ')
atomic_masses = []
for i in range(len(mass_lines_list)):
    atomic_masses.append(float(re.search(r'\d+.\d+', mass_lines_list[i]).group()))


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
    # open POSCAR file and read lines into a list
    pos_f = open(poscar_filepath, 'r')
    pos_lines = pos_f.readlines()
    pos_f.close()
    
    # convert species names line from POSCAR file from Fortran to a list 
    pos_spec_line = ff.FortranRecordReader('(10A5)')  # arbitrarily choose max species number of 10
    pos_species = pos_lines[5]
    pos_spec_list = pos_spec_line.read(pos_species)
    for i in range(len(pos_spec_list)):
        pos_spec_list[i] = pos_spec_list[i].strip()
    pos_spec_list = list(filter(None, pos_spec_list)) # remove empty values from list
    
    # open POTCAR file and read lines into a list
    pot_f = open(potcar_filepath, 'r')
    pot_lines = pot_f.readlines()
    pot_f.close()
    
    # find POTCAR lines with species names
    pot_spec_lines = subprocess.run(['grep', 'TITEL', potcar_filepath], capture_output = True, text = True)
    pot_spec_lines_list = pot_spec_lines.stdout.split('\n ')
    pot_spec_list = []
    for i in range(len(pot_spec_lines_list)):
        pot_spec_list.append(re.split(r'[;,\s]\s*', pot_spec_lines_list[i])[4])
    
    if pos_spec_list == pot_spec_list:
        return True
    else:
        return False


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
def vel_to_POSCAR(vel, atom_type, atom_no, filepath):
    """
    Function to write velocity vector of an atom in the lattice to the POSCAR file.
    Given a velocity vector (Ang./fs), atom type, and the number of the atom
    (from the POSCAR file). The velocity vector is a 1D numpy array, the atom
    type is a pymatgen Element, and the atom number is an integer.
    """
    # open POSCAR file and read lines into a list
    f = open(filepath, 'r')
    f_lines = f.readlines()
    f.close()
    
    # convert velocity vector into Fortran format
    vel_line = ff.FortranRecordWriter('(3E16.8)')
    vel_vector = vel_line.write(vel)+'\n'
    
    # convert species names line from POSCAR file from Fortran to a list 
    spec_line = ff.FortranRecordReader('(10A5)')  # arbitrarily choose max species number of 10
    species = f_lines[5]
    spec_list = spec_line.read(species)
    for i in range(len(spec_list)):
        spec_list[i] = spec_list[i].strip()
    spec_list = list(filter(None, spec_list)) # remove empty values from list
    
    # convert ions per species line from POSCAR file from Fortran to a list
    no_line = ff.FortranRecordReader('(10I6)')  # arbitrarily choose max species number of 10
    ion_nos = f_lines[6]
    no_list = no_line.read(ion_nos)
    no_list = list(filter(None, no_list)) # remove empty values from list
    total_ions = sum(no_list)
    
    # check if empty velocity lines exist, if not, add lines
    if len(f_lines) < 2*total_ions:
        if f_lines[-1].isspace() == False:
            f_lines += ['\n']
        f_lines += ['\n'] + [vel_line.write(np.array([0, 0, 0]))+'\n' for i in range(total_ions)]
    
    # replace velocity line for desired atom with converted velocity vector
    for i in range(len(spec_list)):           
        if atom_type.lower() in spec_list[i].lower():
            f_lines[8+total_ions+atom_no+sum(no_list[:i])] = vel_vector
    
    f = open(filepath, 'w+')
    f.writelines(f_lines)
    f.close()
    
    return


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

# open POSCAR file and read lines into a list
pos_f = open(pos_file, 'r')
pos_lines = pos_f.readlines()
pos_f.close()

# open POSCAR file and read lattice vector lines into a list
ai = [pos_lines[2][1:-1], pos_lines[3][1:-1], pos_lines[4][1:-1]]

# convert lattice vectors line from POSCAR file from Fortran to a list
lat_vec_line = ff.FortranRecordReader('(3F22.16)')
lat_vecs = np.array([lat_vec_line.read(ai[j]) for j in range(len(ai))])

# convert species names line from POSCAR file from Fortran to a list
pos_spec_line = ff.FortranRecordReader('(10A5)')
pos_species = pos_lines[5]
pos_spec_list = pos_spec_line.read(pos_species)
for i in range(len(pos_spec_list)):
    pos_spec_list[i] = pos_spec_list[i].strip()

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
