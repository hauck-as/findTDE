#!/usr/bin/env python3
"""Python module used to analyze find_tde calculations."""
import os
import glob
from pathlib import Path

import subprocess
import re
import fortranformat as ff
import pprint

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
# import seaborn as sns

from findtde.utils import *
from findtde.plotting import *
from findtde.io import *

base_path = Path.cwd()
bin_path, inp_path, perfect_path = base_path / 'bin', base_path / 'inp', base_path / 'perfect'


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
    pos_spec_line = ff.FortranRecordReader('(2A5)')
    pos_species = pos_lines[5]
    pos_spec_list = pos_spec_line.read(pos_species)
    for i in range(len(pos_spec_list)):
        pos_spec_list[i] = pos_spec_list[i].strip()
    
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


# function to analyze find_tde runs and determine if there are any issues
def check_find_tde_runs(tde_calc_dir=base_path, program='vasp', ke_tol=1, temp_tol=0.6):
    """
    Function used to check if TDE calcs have proper velocity given to the correct atom, correct KE, if the final temperature is proper relative to initial temperature.
    Input is a directory path string containing find_tde calculations, and function returns a dictionary of runs with errors.
    """
    # empty dicts for final checked simulations and found errors
    all_find_tde_checks_dict, all_find_tde_errs_dict = {}, {}
    
    # find each overall findTDE run directory
    data_filepaths = np.array(glob.glob(os.path.join(tde_calc_dir, '*', '*_data.csv'), recursive=True))
    
    # determine values for a set of findTDE runs corresponding to a given data file
    for i in range(data_filepaths.shape[0]):
        pseudo_path = data_filepaths[i].rsplit('/', 1)[0]    # relative path for findTDE pseudo
        # vasp runs use pseudo_atominfo, lammps runs have an additional _lmp after
        if program == 'vasp':
            pseudo_atom = pseudo_path.rsplit('/', 1)[-1]    # string for pseudo and atom type/number
        elif program == 'lammps':
            pseudo_atom = pseudo_path.rsplit('/', 1)[-1].rsplit('_', 1)[0]    # string for pseudo and atom type/number
        pseudo, atom_info = pseudo_atom.rsplit('_', 1)    # separate pseudo/atom info
        atom_type_ideal, atom_num_ideal = re.search(r'\D+', atom_info).group(), int(re.search(r'\d+', atom_info).group())    # specified atom type and number for findTDE run
        find_tde_checks_dict = {
            'tde': [],
            'knockout': [],
            'ke': [],
            'temp': []
                             }
        
        # check if TDE was found from last line of outfile
        out_filepath = glob.glob(os.path.join(pseudo_path, '*_out.txt'))[0]
        with open(out_filepath) as of:
            for line in of:
                pass
            last_line = line
        
        if last_line.startswith('Loop complete, TDE found:'):
            find_tde_checks_dict['tde'].append(True)
        elif last_line.startswith('Loop not complete, TDE above cutoff:'):
            find_tde_checks_dict['tde'].append('Cutoff')
        else:
            find_tde_checks_dict['tde'].append(False)
            
        # checking values for a VASP run
        if program == 'vasp':
            for subdir, dirs, files in os.walk(pseudo_path):
                for k in range(len(dirs)):
                    if dirs[k] == 'cgm':    # ignore post-AIMD CGM runs
                        continue
                    else:
                        ke_calc_path = os.path.join(subdir, dirs[k])

                        # compare specified KE from directory to the actual initial value from the OUTCAR
                        # ke_tol gives tolerance in eV, necessary due to momentum correction in VASP
                        ke_i_ideal = float(re.search(r'\d+', dirs[k]).group())
                        ke_i_actual = float(re.search(r'-?\d*\.{0,1}\d+', subprocess.run(['grep', 'EKIN', os.path.join(ke_calc_path, 'OUTCAR')], capture_output = True, text = True).stdout.strip('\n').split('\n')[0]).group())
                        if (ke_i_actual >= ke_i_ideal - ke_tol) and (ke_i_actual <= ke_i_ideal + ke_tol):
                            find_tde_checks_dict['ke'].append(True)
                        else:
                            find_tde_checks_dict['ke'].append(False)

                        # determine final to initial lattice temperature ratio, should be roughly half or less due to equipartition theorem of thermodynamics
                        # temp_tol gives tolerance value for ratio, halving of temperature corresponds to 0.5
                        temps = subprocess.run(['grep', 'EKIN_LAT', os.path.join(ke_calc_path, 'OUTCAR')], capture_output = True, text = True).stdout.strip('\n').split('\n')
                        temp_i = float(re.search(r'\d*\.{0,1}\d+', temps[0].split('(')[-1]).group())
                        temp_f = float(re.search(r'\d*\.{0,1}\d+', temps[-1].split('(')[-1]).group())
                        tfti = temp_f/temp_i
                        if tfti < temp_tol:
                            find_tde_checks_dict['temp'].append(True)
                        else:
                            find_tde_checks_dict['temp'].append(False)

                        # determine from initial POSCAR if velocity vector is given to correct atom
                        pos_f = open(os.path.join(ke_calc_path, 'POSCAR'), 'r')
                        pos_lines = pos_f.readlines()
                        pos_f.close()

                        # convert species names line from POSCAR file from Fortran to a list 
                        spec_line = ff.FortranRecordReader('(2A5)')
                        species = pos_lines[5]
                        spec_list = spec_line.read(species)
                        for i in range(len(spec_list)):
                            spec_list[i] = spec_list[i].strip()

                        # convert ions per species line from POSCAR file from Fortran to a list
                        no_line = ff.FortranRecordReader('(2I6)')
                        ion_nos = pos_lines[6]
                        no_list = no_line.read(ion_nos)
                        total_ions = sum(no_list)

                        # convert velocity vector line from POSCAR file from Fortran format to a list
                        vel_line = ff.FortranRecordWriter('(3E16.8)')
                        for i in range(len(spec_list)):           
                            if atom_type_ideal.lower() in spec_list[i].lower():
                                vel_vector = pos_lines[8+total_ions+atom_num_ideal+sum(no_list[:i])].strip('\n')

                        # check if velocity vector for the desired knockout atom is nonzero
                        if vel_vector == vel_line.write(np.array([0, 0, 0])):
                            find_tde_checks_dict['knockout'].append(False)
                        elif vel_vector != vel_line.write(np.array([0, 0, 0])):
                            find_tde_checks_dict['knockout'].append(True)

        elif program == 'lammps':
            pass
        
        # create summary of errors found and detailed checks done
        all_find_tde_errs_dict[pseudo_atom] = []
        if False in find_tde_checks_dict['tde']:
            all_find_tde_errs_dict[pseudo_atom].append('tde')
        if False in find_tde_checks_dict['knockout']:
            all_find_tde_errs_dict[pseudo_atom].append('knockout')
        if False in find_tde_checks_dict['ke']:
            all_find_tde_errs_dict[pseudo_atom].append('ke')
        if False in find_tde_checks_dict['temp']:
            all_find_tde_errs_dict[pseudo_atom].append('temp')
        if not all_find_tde_errs_dict[pseudo_atom]:
            all_find_tde_errs_dict.pop(pseudo_atom)
        
        all_find_tde_checks_dict[pseudo_atom] = find_tde_checks_dict
    
    return (all_find_tde_errs_dict, all_find_tde_checks_dict)


# create dataframes for final energies and energy differences
def find_tde_analysis(atom_types, atom_nums, datafile='all_tde_data.csv', keyfile=os.path.join(base_path, 'latt_dirs_to_calc.csv')):
    """
    Function to analyze gathered find_tde data and format appropriately. Inputs are
    two lists for each atom type (strings) and atom numbers (integers) for which data
    is desired. Returns a tuple with two dictionaries. The first dictionary contains
    tuples with each dictionary key being the atom type and number concatenated and the
    tuples containing two pandas DataFrames. The first DataFrame is the final energy for
    each calculation direction at each calculated kinetic energy, and the second DataFrame
    is the same for final energy difference from the perfect crystal energy. The second
    dictionary contains keys for the lattice direction pseudos.
    """
    find_tde_data = np.genfromtxt(datafile, delimiter=', ', skip_header=1, dtype=str)
    find_tde_key = np.genfromtxt(keyfile, delimiter=', ', skip_header=4, dtype=str)
    
    ##### need to return the cutoff energy as well as the lattice directions #####
    
    # create overall dictionary for each instance of atom types and numbers
    all_df_dict, pseudo_keys = {}, {}
    
    # find minimum cutoff energy to use for analysis
    ke_cutoff = np.amin(find_tde_key[:, 4].astype('int'))
    
    # sort data by atom types, atom numbers, and calculation directions
    for i in range(len(atom_types)):
        for j in range(len(atom_nums)):
            atom_indices = np.where((find_tde_data[:, 1] == atom_types[i].lower()) & (find_tde_data[:, 2] == str(atom_nums[j])))[0]
            toten_dict, toten_dif_dict = {}, {}
            for c, k in enumerate(atom_indices):
                # determine directions from pseudos and group data by pseudo/direction
                if find_tde_data[k, 0] == find_tde_data[atom_indices[c-1], 0]:
                    toten_dict[find_tde_data[k, 0]].update({int(find_tde_data[k, 3]):float(find_tde_data[k, 4])})
                    toten_dif_dict[find_tde_data[k, 0]].update({int(find_tde_data[k, 3]):float(find_tde_data[k, 5])})
                elif find_tde_data[k, 0] != find_tde_data[atom_indices[c-1], 0]:
                    if find_tde_data[k, 0][-1] == 'L':
                        latdir = find_tde_key[np.where(find_tde_data[k, 0] == find_tde_key[:, 0])[0][0], :][-3:]
                    elif find_tde_data[k, 0][-1] == 'S':
                        latdir = find_tde_key[np.where(find_tde_data[k, 0] == find_tde_key[:, 0])[0][0], :][-2:]
                    pseudo_keys[find_tde_data[k, 0]] = latdir
                    toten_dict[find_tde_data[k, 0]], toten_dif_dict[find_tde_data[k, 0]] = {int(find_tde_data[k, 3]):float(find_tde_data[k, 4])}, {int(find_tde_data[k, 3]):float(find_tde_data[k, 5])}

            # add to overall dictionary if data has been found
            if toten_dict and toten_dif_dict:
                # dataframe from dictionary with final TOTEN values
                toten_df = pd.DataFrame(data=toten_dict)
                toten_df.sort_index(inplace=True)

                # dataframe from dictionary with calculated differences from perfect crystal TOTEN
                toten_dif_df = pd.DataFrame(data=toten_dif_dict)
                toten_dif_df.sort_index(inplace=True)
            
                # add dataframes 
                all_df_dict[atom_types[i].lower()+str(atom_nums[j])] = (toten_df, toten_dif_df)
    
    return (all_df_dict, pseudo_keys)


def find_tde_sph_analysis(find_tde_dE, pseudo_keys, lattice_vecs, e_tol=1.0, ke_cut=45):
    """
    Given a dataframe with lattice direction pseudo columns, initial KE indices, and delta E
    values, returns a dictionary listing this info and the knockout's spherical coordinates.
    """
    find_tde_sph_dict = {}
    
    for i in range(find_tde_dE.columns.shape[0]):
        tde_info_dict = {}
        for j in range(find_tde_dE.index.shape[0]):
            if np.isnan(find_tde_dE[find_tde_dE.columns[i]].loc[find_tde_dE.index[j]]) == True:
                continue
            elif find_tde_dE[find_tde_dE.columns[i]].loc[find_tde_dE.index[j]] > e_tol:
                if find_tde_dE.index[j] <= ke_cut:
                    tde_info_dict['TDE'] = find_tde_dE.index[j]
                    tde_info_dict['dE'] = find_tde_dE[find_tde_dE.columns[i]].loc[find_tde_dE.index[j]]
                elif find_tde_dE.index[j] > ke_cut:
                    tde_info_dict['TDE'] = ke_cut
                break
            else:
                tde_info_dict['TDE'] = ke_cut
            
            if find_tde_dE.columns[i][-1] == 'S':
                tde_info_dict['phi'], tde_info_dict['theta'] = pseudo_keys[find_tde_dE.columns[i]].astype(float)
            elif find_tde_dE.columns[i][-1] == 'L':
                tde_info_dict['latt_dir'] = pseudo_keys[find_tde_dE.columns[i]]
                rho, tde_info_dict['phi'], tde_info_dict['theta'] = lat2sph(np.reshape(pseudo_keys[find_tde_dE.columns[i]], (1, 3)).astype(int), lattice_vecs)[0]
        
        find_tde_sph_dict[find_tde_dE.columns[i]] = tde_info_dict
    
    return find_tde_sph_dict


def generate_tde_sph_arr(tde_data_df, pseudo_keys, lattice_vecs, e_tol=1.0, ke_cut=45, polar_offset=60):
    """
    Function to reorganize dataframe into a Nx3 array with spherical coordinates.
    """
    find_tde_sph_dict = find_tde_sph_analysis(tde_data_df, pseudo_keys=pseudo_keys, lattice_vecs=lattice_vecs, e_tol=e_tol, ke_cut=ke_cut)
    find_tde_pseudos = find_tde_sph_dict.keys()
    find_tde_sph_arr = np.zeros((len(find_tde_pseudos), 3))
    for k in range(len(find_tde_pseudos)):
        key = list(find_tde_pseudos)[k]
        try:
            rpt = np.array([[1.0, find_tde_sph_dict[key]['phi'], find_tde_sph_dict[key]['theta']]])
            rpt_scaled = scale_C6z_symmetry(scale_C3z_symmetry(rpt))    # scale_C3z_symmetry(rpt)
            find_tde_sph_arr[k, :] = np.array([
                # 60 DEG ADDED TO POLAR ANGLE SO THAT 0 deg ALIGNS TO THE RIGHT IN VISUALIZATION
                (round(rpt_scaled[0, 1], 2)+polar_offset),
                round(rpt_scaled[0, 2], 2),
                round(find_tde_sph_dict[key]['TDE'], 2)
            ])
        except:
            pass
    find_tde_sph_arr = find_tde_sph_arr[~np.all(find_tde_sph_arr == 0, axis=1)]
    
    return (find_tde_sph_arr, find_tde_pseudos)


# function to gather data from find_tde lammps runs, recalculated after thermalizing
def tde_thrm_data_gather(ofile='all_tde_data_lmp_thrm.csv', tde_calc_dir=base_path, e_tol=1.0):
    """
    Gathers all TDE data formatted as find_tde data.
    """
    tde_data_header = 'Lattice Direction Pseudo, Atom Type, Atom Number, KE, Final Energy, Delta Energy\n'
    thrm_data_header = 'Lattice Direction Pseudo_Atom TypeAtom Number, KE, Calcs with Defects, Total Calcs\n'
    out_filepaths = np.array(glob.glob(os.path.join(tde_calc_dir, '*/*_out.txt'), recursive=True))
    
    for i in range(out_filepaths.shape[0]):
        pseudo_path = out_filepaths[i].rsplit('/', 1)[0]
        pseudo_atom = pseudo_path.rsplit('/', 1)[-1].rsplit('_', 1)[0]
        pseudo_datafile = os.path.join(pseudo_path, pseudo_atom + '_data.csv')
        pf = open(pseudo_datafile, 'w+')
        pf.write(tde_data_header)
        
        for subdir, dirs, files in os.walk(pseudo_path):
            for k in range(len(dirs)):
                ke_calc_path = os.path.join(subdir, dirs[k])
                ke = re.search(r'\d+', dirs[k]).group()
                # Ei from log_thermal.lammps may be different depending on size
                # Ei = re.search(r'-?\d*\.{0,1}\d+', subprocess.run(['grep', 'Initial energy of atoms', os.path.join(ke_calc_path, 'log_thermal.lammps')], capture_output = True, text = True).stdout.strip('\n').split('\n')[-1]).group()
                log_filepaths = np.array(glob.glob(os.path.join(ke_calc_path, 'log_*.lammps'), recursive=True))
                
                for j in range(log_filepaths.shape[0]):
                    if log_filepaths[j].rsplit('/', 1)[-1] == 'log_thermal.lammps':
                        pass
                    else:
                        pseudo, atom_info = pseudo_atom.rsplit('_', 1)
                        atom_type, atom_num = re.search(r'\D+', atom_info).group(), re.search(r'\d+', atom_info).group()
                        Ei = re.search(r'-?\d*\.{0,1}\d+', subprocess.run(['grep', 'Initial energy of atoms', log_filepaths[j]], capture_output = True, text = True).stdout.strip('\n').split('\n')[-1]).group()
                        Ef = re.search(r'-?\d*\.{0,1}\d+', subprocess.run(['grep', 'Final energy of atoms', log_filepaths[j]], capture_output = True, text = True).stdout.strip('\n').split('\n')[-1]).group()
                        dE = float(Ef) - float(Ei)
                        data_line = pseudo + ', ' + atom_type + ', ' + atom_num + ', ' + ke + ', ' + Ef + ', ' + str(dE) + '\n'
                        pf.write(data_line)
        pf.close()
    
    
    data_filepaths = np.array(glob.glob(os.path.join(tde_calc_dir, '*/*_data.csv'), recursive=True))
    tde_thrm_dict = {}
    
    for file in data_filepaths:
        pseudo_data = np.genfromtxt(file, delimiter=', ', skip_header=1,
                                    dtype=[('f0', '<U4'), ('f1', '<U4'), ('f2', '<i4'),
                                           ('f3', '<i4'), ('f4', '<f8'), ('f5', '<f8')]
                                   )
        pseudo_atom = pseudo_data[0][0] + '_' + pseudo_data[0][1] + str(pseudo_data[0][2])
        tde_thrm_dict[pseudo_atom] = np.array([[pseudo_data[0][3], 0, 0]], dtype=int)
        
        for i in range(pseudo_data.shape[0]):
            if pseudo_data[i][3] not in tde_thrm_dict[pseudo_atom][:, 0]:
                tde_thrm_dict[pseudo_atom] = np.append(tde_thrm_dict[pseudo_atom], np.array([[pseudo_data[i][3], 0, 0]]), axis=0)
            for j in range(tde_thrm_dict[pseudo_atom].shape[0]):
                if pseudo_data[i][3] == tde_thrm_dict[pseudo_atom][j, 0]:
                    tde_thrm_dict[pseudo_atom][j, 2] += 1
                    if pseudo_data[i][-1] >= e_tol:
                        tde_thrm_dict[pseudo_atom][j, 1] += 1
                    elif pseudo_data[i][-1] < e_tol:
                        pass
                elif pseudo_data[i][3] != tde_thrm_dict[pseudo_atom][j, 0]:
                    pass
    
    of = open(tde_calc_dir+ofile, 'w+')
    of.write(thrm_data_header)
    of_lines = []
    
    for psatan in tde_thrm_dict.keys():
        for n in range(tde_thrm_dict[psatan].shape[0]):
            of_lines.append(psatan + ', ' + str(tde_thrm_dict[psatan][n, 0]) + ', ' + str(tde_thrm_dict[psatan][n, 1]) + ', ' + str(tde_thrm_dict[psatan][n, 2]) + '\n')
    for line in of_lines:
        of.write(line)
    
    of.close()
    
    return tde_thrm_dict


if __name__ == "__main__":
    FILES_LOC = base_path
    
    run_type, knockout_atom = 'vasp', 'ga34'
    
    # read in VASP and LAMMPS perfect structure information
    if run_type == 'vasp':
        vasp_ref_file = os.path.join(FILES_LOC, 'inp', 'POSCAR')
    elif run_type == 'lammps':
        vasp_ref_file, lmp_ref_file = os.path.join(FILES_LOC, 'inp', 'POSCAR'), os.path.join(FILES_LOC, 'inp', 'read_data_perfect.lmp')
    pos_f = open(vasp_ref_file, 'r')
    pos_lines = pos_f.readlines()
    pos_f.close()

    # open POSCAR file and read lattice vector lines into a list
    ai = [pos_lines[2][1:-1], pos_lines[3][1:-1], pos_lines[4][1:-1]]

    # convert lattice vectors line from POSCAR file from Fortran to a list
    lat_vec_line = ff.FortranRecordReader('(3F22.16)')
    vasp_lattice_vecs = np.array([lat_vec_line.read(ai[j]) for j in range(len(ai))])
    vasp_params = np.array([np.linalg.norm(vasp_lattice_vecs[0, :]), np.linalg.norm(vasp_lattice_vecs[1, :]), np.linalg.norm(vasp_lattice_vecs[2, :])])

    # check if runs have issues
    if run_type == 'vasp':
        print(check_find_tde_runs(tde_calc_dir=FILES_LOC, program='vasp')[0])
    elif run_type == 'lammps':
        print(check_find_tde_runs(tde_calc_dir=FILES_LOC, program='lammps')[0])
    
    # gather TDE data from calculations 
    tde_data_gather(ofile='all_tde_data.csv', tde_calc_dir=FILES_LOC)
    
    ######## TDE data ########
    if run_type == 'vasp':
        # organize data from csv file into dataframes
        all_find_tde_dfs, pseudo_keys = find_tde_analysis(['Ga', 'N'], [34, 35], datafile='all_tde_data.csv', keyfile=os.path.join(FILES_LOC, 'latt_dirs_to_calc.csv'))
        ga_find_tde_df, n_find_tde_df = all_find_tde_dfs['ga34'][1], all_find_tde_dfs['n35'][1]
        
        # reorganize data into array
        ga_tde_sph_arr, ga_tde_pseudos = generate_tde_sph_arr(ga_find_tde_df, pseudo_keys, lattice_vecs=vasp_lattice_vecs, e_tol=1.0, ke_cut=45, polar_offset=angle_between([1., 0., 0.], vasp_lattice_vecs[0]))
        n_tde_sph_arr, n_tde_pseudos = generate_tde_sph_arr(n_find_tde_df, pseudo_keys, lattice_vecs=vasp_lattice_vecs, e_tol=1.0, ke_cut=45, polar_offset=angle_between([1., 0., 0.], vasp_lattice_vecs[0]))
        
        # from Victor, read data into a (nsamples x 3) array (x, y, f(x, y)), interpolate and plot heatmap data
        # ps_ga_tde, ts_ga_tde, es_ga_tde = idw_heatmap(ga_tde_sph_arr, RES=ga_tde_sph_arr.shape[0], P=5)
        # ps_n_tde, ts_n_tde, es_n_tde = idw_heatmap(n_tde_sph_arr, RES=n_tde_sph_arr.shape[0], P=5)
        
    elif run_type == 'lammps':
        # organize data from csv file into dataframes
        all_find_tde_lmp_dfs, lmp_pseudo_keys = find_tde_analysis(['Ga', 'N'], [34, 35], datafile='all_tde_data_lmp.csv', keyfile=os.path.join(FILES_LOC, 'latt_dirs_to_calc.csv'))
        ga_find_tde_lmp_df, n_find_tde_lmp_df = all_find_tde_lmp_dfs['ga34'][1], all_find_tde_lmp_dfs['n35'][1]
        
        # reorganize data into array
        ga_tde_lmp_sph_arr, ga_tde_lmp_pseudos = generate_tde_sph_arr(ga_find_tde_lmp_df, lmp_pseudo_keys, lattice_vecs=vasp_lattice_vecs, e_tol=1.0, ke_cut=100, polar_offset=angle_between([1., 0., 0.], vasp_lattice_vecs[0]))
        n_tde_lmp_sph_arr, n_tde_lmp_pseudos = generate_tde_sph_arr(n_find_tde_lmp_df, lmp_pseudo_keys, lattice_vecs=vasp_lattice_vecs, e_tol=1.0, ke_cut=100, polar_offset=angle_between([1., 0., 0.], vasp_lattice_vecs[0]))
        
        # from Victor, read data into a (nsamples x 3) array (x, y, f(x, y)), interpolate and plot heatmap data
        # ps_ga_lmp_tde, ts_ga_lmp_tde, es_ga_lmp_tde = idw_heatmap(ga_tde_lmp_sph_arr, RES=ga_tde_lmp_sph_arr.shape[0], P=5)
        # ps_n_lmp_tde, ts_n_lmp_tde, es_n_lmp_tde = idw_heatmap(n_tde_lmp_sph_arr, RES=n_tde_lmp_sph_arr.shape[0], P=5)
    
    generate_tde_line_plot(ga_find_tde_df, im_write=True, im_name='find_tde_lineplot.png')
    generate_tde_scatter_plot(ga_tde_sph_arr, ga_tde_pseudos, txt_show=False, im_write=True, im_name='gan_ga_tde.png')
    # generate_tde_scatter_plot(n_tde_sph_arr, n_tde_pseudos, txt_show=False, im_write=True, im_name='gan_n_tde.png')
