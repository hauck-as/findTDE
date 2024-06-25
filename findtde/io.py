#!/usr/bin/env python3
"""Python module of file management functions for findTDE calculations."""
import os
import glob
from pathlib import Path

import numpy as np
import pandas as pd

base_path = Path.cwd()
bin_path, inp_path, perfect_path = base_path / 'bin', base_path / 'inp', base_path / 'perfect'


def find_contcars(tde_calc_dir=base_path, pseudo='*'):
    """
    Function to create an array of post-CGM CONTCAR files for defect analysis.
    """
    contcar_files = np.array(glob.glob(os.path.join(tde_calc_dir, pseudo, '*', 'cgm', 'CONTCAR'), recursive=True))
    return contcar_files


def find_dumpfiles(tde_calc_dir=base_path, pseudo='*'):
    """
    Function to create an array of post-MD dump files for defect analysis.
    """
    dump_files = np.array(glob.glob(os.path.join(tde_calc_dir, pseudo, '*', 'dump.final'), recursive=True))
    return dump_files

def pseudo_keys_to_file(pseudo_keys, pseudo_keyfile=os.path.join(base_path, 'pseudo_keys.csv')):
    """
    Function to convert a pseudo keys dictionary to a CSV.
    """
    pk_df = pd.DataFrame.from_dict(pseudo_keys, orient='index')
    pk_df.to_csv(pseudo_keyfile)
    return


def pseudo_keys_from_file(pseudo_keyfile=(base_path / 'pseudo_keys.csv')):
    """
    Function to convert a CSV file of pseudo keys to a dictionary.
    """
    pk_df = pd.read_csv(pseudo_keyfile, index_col=0)
    pks_dict = pk_df.to_dict('split')
    pks_dict = dict(zip(pks_dict['index'], pks_dict['data']))
    return pks_dict


def write_sph_dirs(rpt, filepath=os.path.join(base_path, 'sph_directions.csv')):
    """
    Function to write spherical directions to a text/csv file for multi_tde.py use.
    Given a 2D array of spherical coordinates, writes array values to the file
    represented by the filepath string.
    """
    sph_f_header = '########################################\n' + \
    'Spherical directions for multi_tde.py calculations\n' + \
    'r, p, t\n' + \
    '########################################\n'
    
    sph_f = open(filepath, 'a+')
    if os.stat(filepath).st_size == 0:
        sph_f.write(sph_f_header)
    sph_f.write('########################################\n')
    for i in range(rpt.shape[0]):
        sph_line = ', '.join([str(round(j, 2)) for j in rpt[i, :]])
        sph_f.write(sph_line+'\n')
    sph_f.close()
    return


# function to gather data from find_tde runs
def tde_data_gather(ofile='all_tde_data.csv', tde_calc_dir=base_path):
    """
    Gathers all TDE data formatted as find_tde data.
    """
    tde_data_header = 'Lattice Direction Pseudo, Atom Type, Atom Number, KE, Final Energy, Delta Energy\n'
    
    data_filepaths = glob.glob(os.path.join(tde_calc_dir, '*', '*_data.csv'), recursive=True)
    of = open(ofile, 'a+')
    of.write(tde_data_header)
    
    for file in data_filepaths:
        data = np.genfromtxt(file, delimiter=', ', skip_header=1,
                             dtype=[('f0', '<U4'), ('f1', '<U4'), ('f2', '<i4'),
                                    ('f3', '<i4'), ('f4', '<f8'), ('f5', '<f8')]
                            )
        np.savetxt(of, data, delimiter=', ', newline='\n',
                   fmt=['%s', '%s', '%d', '%d', '%f', '%f']
                  )
    
    of.close()


# function to determine TDE from a single datafile
def get_tde_from_datafile(data_filepath='*_data.csv', e_tol=1.0, ke_cut=45):
    """
    Scans a findTDE *_data.csv file and returns the TDE value for the specified knockout type.
    """
    # read in tde data into structured array
    tde_data = np.genfromtxt(data_filepath, delimiter=', ', skip_header=1,
                             dtype=[('pseudo', '<U4'), ('atom_type', '<U4'), ('atom_num', '<i4'),
                                ('ke', '<i4'), ('e_fin', '<f8'), ('dE', '<f8')]
                            )
    
    # sort array from lowest KE to highest KE
    tde_data = np.sort(tde_data, order='ke')
    
    # find TDE value from energy difference
    # first KE corresponding to a dE > e_tol is the TDE in the sorted array
    tde_value = 0
    for i in range(tde_data.shape[0]):
        if np.isnan(tde_data[i]['dE']) == True:
            continue
        elif tde_data[i]['dE'] > e_tol:
            if tde_data[i]['ke'] <= ke_cut:
                tde_value = tde_data[i]['ke']
            elif tde_data[i]['ke'] > ke_cut:
                tde_value = ke_cut
            break
        else:
            tde_value = ke_cut
    
    return tde_value
