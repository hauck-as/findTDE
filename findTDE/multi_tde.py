#!/usr/bin/env python3
"""Python module used to submit several find_tde jobs at once."""
import os
import re
import subprocess
import math as m
import numpy as np

base_path = os.getcwd()
bin_path, inp_path, perfect_path = os.path.relpath('bin', base_path), os.path.relpath('inp', base_path), os.path.relpath('perfect', base_path)


def write_find_tde_calcs(direction, atom_type, atom_number, ke_i=25, ke_cut=40, mode='L', directions_filepath=os.path.relpath('latt_dirs_to_calc.csv', os.getcwd())):
    """
    Writes the text file used by find_tde to perform calculations. Requires a direction list (L or S)
    as well as the atom type (string) and number (int), uses 25 eV initial KE and lattice direction
    mode by default (can use any integer KE or spherical direction mode (S)).
    """
    lat_dir_list_header_text = '########################################\n' + \
    '# format of text file\n' + \
    '# nL    atom_type    atom_number    ke_i    ke_cut    u    v    w\n' + \
    '# n+1S    atom_type    atom_number    ke_i    ke_cut    r    p    t\n' + \
    '########################################'
    
    # if file does exist, append directions to it starting from previous number
    if os.path.isfile(directions_filepath):
        ld_f = open(directions_filepath, 'r+')
        print('Appending to file')
        n, n_all = 0, []
        ld_f_lines = ld_f.readlines()
        for line in ld_f_lines:
            line = line.strip('\n')
            if line[0] != '#':
                n_all.append(re.split(r'[;,\s]\s*', line)[0][:-1])
                dir_n = re.split(r'[;,\s]\s*', line)[5:]
                if np.all(np.array(dir_n, dtype=float) == direction.astype(float)):
                    print('Direction matches previous pseudo')
                    n = re.split(r'[;,\s]\s*', line)[0][:-1]
        n_prev_max = int(np.amax(np.array((n_all), dtype=int)))
        if n != 0:
            ld_f.write(', '.join(['\n'+str(n)+mode, atom_type, str(atom_number), str(ke_i), str(ke_cut), str(direction[0]), str(direction[1]), str(direction[2])]))
        else:
            ld_f.write(', '.join(['\n'+str(1+n_prev_max)+mode, atom_type, str(atom_number), str(ke_i), str(ke_cut), str(direction[0]), str(direction[1]), str(direction[2])]))
        ld_f.close()
    
    # if file does not exist, create file with header and add direction to it
    else:
        ld_f = open(directions_filepath, 'w+')
        print('Creating new file')
        ld_f.write(lat_dir_list_header_text)
        ld_f.write(', '.join(['\n'+str(1)+mode, atom_type, str(atom_number), str(ke_i), str(ke_cut), str(direction[0]), str(direction[1]), str(direction[2])]))
        ld_f.close()
    
    return


def find_multiple_tde(directions, atom_type, atom_number, ke_i=25, ke_cut=40, mode='L', conv='standard', prog='vasp', lmp_ff='AlGaN.sw', screen_num=0):
    """
    Function used to write and submit multiple instances for find_tde.
    """
    for i in range(directions.shape[0]):
        write_find_tde_calcs(directions[i, :], atom_type, atom_number, ke_i=ke_i, ke_cut=ke_cut, mode=mode)
        # os.system('screen -S TDE'+str(i)+' find_tde -c '+str(conv))
        if prog == 'vasp':
            subprocess.run(['screen', '-S', 'TDE'+str(screen_num+i)+atom_type.lower()+str(atom_number), 'find_tde', '-c', str(conv), '-p', str(prog)], capture_output = True, text = True)
        elif prog == 'lammps':
            subprocess.run(['screen', '-S', 'TDE'+str(screen_num+i)+atom_type.lower()+str(atom_number), 'find_tde', '-c', str(conv), '-p', str(prog), '-f', str(lmp_ff)], capture_output = True, text = True)
    return


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
        rpt[i, 0] = 1.
    
    return rpt


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


def write_sph_dirs(rpt, directions_filepath=os.path.relpath('sph_directions.csv', os.getcwd())):
    """
    Function to write spherical directions to a text/csv file for multi_tde.py use.
    Given a 2D array of spherical coordinates, writes array values to the file
    represented by the filepath string.
    """
    sph_f_header = '########################################\n' + \
    'Spherical directions for multi_tde.py calculations\n' + \
    'r, p, t\n' + \
    '########################################\n'
    
    sph_f = open(directions_filepath, 'a+')
    if os.stat(directions_filepath).st_size == 0:
        sph_f.write(sph_f_header)
    sph_f.write('########################################\n')
    for i in range(rpt.shape[0]):
        sph_line = ', '.join([str(round(j, 2)) for j in rpt[i, :]])
        sph_f.write(sph_line+'\n')
    sph_f.close()
    return


# generate directions for calculations
"""l_directions = np.array([
    [-1, 4, 1],
    [-1, 7, 3],
    [-1, 8, 4],
    [-2, -1, 0],
    [-4, 11, 0],
    [-4, 4, 3],
    [0, 0, -1],
    [0, 0, 1],
    [0, 1, -5],
    [0, 1, 0],
    [0, 2, 1],
    [0, 4, 1],
    [0, 5, 2],
    [1, -1, 1],
    [1, 1, 1],
    [1, 2, -5],
    [1, 2, 0],
    [1, 7, 0],
    [2, -2, 1],
    [2, 0, -5],
    [2, 1, 0],
    [2, 4, -5],
    [4, -3, 0],
    [4, -5, 4],
    [4, 1, 3],
    [4, 3, 1],
    [5, 2, 1],
    [5, 3, 0],
    [7, 0, 3],
    [7, 4, 5],
    [8, 4, 3],
    [1, 1, 0],
    [-1, -1, 0],
    [-1, 1, 0]
])

s_directions = np.array([
    # [1., 30., 105.],
    # [1., 90., 60.]
    [1., 90., 97.5]
])"""

ga_s_gan_dirs = np.array([
    [1., 45., 112.5],
    [1., 30., 105.],
    [1., 90., 60.]
])

ga_l_gan_dirs = np.array([
    [5, 2, 1],
    [-4, 4, 3],
    [-2, -1, 0]
])

n_s_gan_dirs = np.array([
    [1., 90., 112.5],
    [1., 90., 97.5]
])

n_l_gan_dirs = np.array([
    [-1, 1, 0]
])

al_s_aln_dirs = np.array([
    [1., 82.5, 52.5],
    [1., 82.5, 60.],
    [1., 82.5, 67.5],
    [1., 82.5, 45.]
])

n_s_aln_dirs = np.array([
    [1., 90., 45.],
    # [1., 90., 112.5]
])

# write_sph_dirs(scale_C3z_symmetry(lat2sph(all_tde_data[1], gan_lattice_vecs)), filepath=os.path.relpath('sph_directions.csv', base_path))
# sph_directions = np.genfromtxt(os.path.relpath('sph_directions.csv', base_path), delimiter=',', comments='#', skip_header=4)
# sph_directions_ext_s1 = np.genfromtxt(os.path.relpath('sph_directions_ext_set1.csv', base_path), delimiter=',', comments='#', skip_header=4)
# sph_directions_ext_s2 = np.genfromtxt(os.path.relpath('sph_directions_ext_set2.csv', base_path), delimiter=',', comments='#', skip_header=4)
# sph_directions_ext_s3 = np.genfromtxt(os.path.relpath('sph_directions_ext_set3.csv', base_path), delimiter=',', comments='#', skip_header=4)

# write directions and run find_tde
# ga 34, n 35
# find_multiple_tde(s_directions, 'ga', 34, ke_i=25, ke_cut=45, mode='S', conv='midpoint', prog='vasp')

# calculate lowest energy directions from GaN/AlN calcs, both Ga/Al & N knockouts, in VASP
find_multiple_tde(ga_s_gan_dirs, 'al', 34, ke_i=25, ke_cut=45, mode='S', conv='midpoint', prog='vasp', screen_num=4410)
find_multiple_tde(ga_l_gan_dirs, 'al', 34, ke_i=25, ke_cut=45, mode='L', conv='midpoint', prog='vasp', screen_num=4420)
find_multiple_tde(n_s_gan_dirs, 'n', 35, ke_i=25, ke_cut=45, mode='S', conv='midpoint', prog='vasp', screen_num=4430)
find_multiple_tde(n_l_gan_dirs, 'n', 35, ke_i=25, ke_cut=45, mode='L', conv='midpoint', prog='vasp', screen_num=4440)
find_multiple_tde(al_s_aln_dirs, 'al', 34, ke_i=25, ke_cut=45, mode='S', conv='midpoint', prog='vasp', screen_num=4450)
find_multiple_tde(n_s_aln_dirs, 'n', 35, ke_i=25, ke_cut=45, mode='S', conv='midpoint', prog='vasp', screen_num=4460)

# calculate extended spherical direction mesh in LAMMPS
# find_multiple_tde(sph_directions_ext_s3, 'ga', 34, ke_i=10, ke_cut=100, mode='S', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=1100)
# find_multiple_tde(sph_directions_ext_s3, 'n', 35, ke_i=10, ke_cut=100, mode='S', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=1400)

# calculate lowest energy directions from GaN/AlN calcs, both Ga/Al & N knockouts, in LAMMPS
# find_multiple_tde(ga_s_gan_dirs, 'al', 34, ke_i=10, ke_cut=100, mode='S', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=7010)
# find_multiple_tde(ga_l_gan_dirs, 'al', 34, ke_i=10, ke_cut=100, mode='L', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=7020)
# find_multiple_tde(n_s_gan_dirs, 'n', 35, ke_i=10, ke_cut=100, mode='S', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=7030)
# find_multiple_tde(n_l_gan_dirs, 'n', 35, ke_i=10, ke_cut=100, mode='L', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=7040)
# find_multiple_tde(al_s_aln_dirs, 'al', 34, ke_i=10, ke_cut=100, mode='S', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=7050)
# find_multiple_tde(n_s_aln_dirs, 'n', 35, ke_i=10, ke_cut=100, mode='S', conv='midpoint', prog='lammps', lmp_ff='AlGaN.sw', screen_num=7060)