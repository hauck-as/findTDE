#!/usr/bin/env python3
"""Python module used to characterize defects resulting from TDE calculations."""
import os
import glob
import sys

import subprocess
import re
import fortranformat as ff

import math
from fractions import Fraction
import random as rand

import numpy as np
import pandas as pd

import pymatgen
import pymatgen.core as mg
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Locpot
import pymatgen.analysis.defects as pad
from pymatgen.analysis.defects import core, ccd, finder

import ovito as ov
from ovito.io import import_file, export_file
from ovito.pipeline import FileSource
from ovito.vis import Viewport
from ovito.modifiers import *

base_path = os.getcwd()
poscar_path, defects_path = os.path.join(base_path, 'inp/POSCAR'), os.path.join(base_path, 'defect_analysis')


def find_contcars(base_path, pseudo='*'):
    """
    Function to create an array of post-CGM CONTCAR files for defect analysis.
    """
    contcar_files = np.array(glob.glob(os.path.join(base_path, pseudo+'/*/cgm/CONTCAR'), recursive=True))
    return contcar_files


contcar_paths = np.concatenate((find_contcars(base_path, pseudo='17L_ga34'), find_contcars(base_path, pseudo='57S_ga34')))

# open POSCAR file and read lines into a list
pos_f = open(poscar_path, 'r')
pos_lines = pos_f.readlines()
pos_f.close()

# open POSCAR file and read lattice vector lines into a list
ai = [pos_lines[2][1:-1], pos_lines[3][1:-1], pos_lines[4][1:-1]]

# convert lattice vectors line from POSCAR file from Fortran to a list
lat_vec_line = ff.FortranRecordReader('(3F22.16)')
gan_lattice_vecs = np.array([lat_vec_line.read(ai[j]) for j in range(len(ai))])

# convert species names line from POSCAR file from Fortran to a list
pos_spec_line = ff.FortranRecordReader('(2A5)')
pos_species = pos_lines[5]
pos_spec_list = pos_spec_line.read(pos_species)
for i in range(len(pos_spec_list)):
    pos_spec_list[i] = pos_spec_list[i].strip()


###########################################################################
######################## PYMATGEN DEFECTS ANALYSIS ########################
###########################################################################
def find_defects_multi_files(base_pos_path, cont_paths):
    """
    Given an initial POSCAR file and any number of post-AIMD, post-CGM CONTCAR files, returns the
    defect positions by comparing final atomic positions to initial atomic positions. File input
    is a string for the POSCAR filepath and an array of strings for the CONTCAR filepaths.
    """
    base_struct = Structure.from_file(base_pos_path)

    dsf = finder.DefectSiteFinder()
    defect_pos = np.zeros((cont_paths.shape[0], 3))
    
    for q in range(cont_paths.shape[0]):
        defect_struct = Structure.from_file(cont_paths[q])
        fpos = dsf.get_defect_fpos(defect_structure=defect_struct, base_structure=base_struct)
        fpos -= np.round(fpos)
        defect_pos[q, :] = fpos
    
    return defect_pos


def fp_sep_dist(base_pos_path, cont_paths, atom_type='ga', atom_no=34, lattice_vecs=gan_lattice_vecs):
    """
    Given an initial POSCAR file and any number of post-AIMD, post-CGM CONTCAR files, returns the
    Frenkel pair separation distance after determining defect positions in the structure. File
    input is a string for the POSCAR filepath and an array of strings for the CONTCAR filepaths.
    """
    # this method assumes the vacancy is at the initial atomic position
    defect_pos = find_defects_multi_files(base_pos_path, cont_paths)@lattice_vecs
    v_pos = np.zeros((defect_pos.shape[0], 3))
    fp_sep = np.zeros((defect_pos.shape[0], 1))

    # open POSCAR file and read lines into a list
    pos_f = open(base_pos_path, 'r')
    pos_lines = pos_f.readlines()
    pos_f.close()

    # convert species names line from POSCAR file from Fortran to a list
    pos_spec_line = ff.FortranRecordReader('(2A5)')
    pos_species = pos_lines[5]
    pos_spec_list = pos_spec_line.read(pos_species)
    for i in range(len(pos_spec_list)):
        pos_spec_list[i] = pos_spec_list[i].strip()
    
    # convert ions per species line from POSCAR file from Fortran to a list
    pos_no_line = ff.FortranRecordReader('(2I6)')
    ion_nos = pos_lines[6]
    no_list = pos_no_line.read(ion_nos)
    total_ions = sum(no_list)
    
    # convert position line from POSCAR file from Fortran to a list
    pos_line = ff.FortranRecordReader('(3F20.16)')
    
    for i in range(len(pos_spec_list)):
        if atom_type.lower() in pos_spec_list[i].lower():
            atom_pos_line = 8+atom_no+sum(no_list[:i])
    
    for n in range(defect_pos.shape[0]):
        v_pos[n, :] = pos_line.read(pos_lines[atom_pos_line].strip('\n'))
        fp_sep[n] = math.dist(defect_pos[n, :], np.array([0, 0, 0]))
    
    return fp_sep

print('######################## PYMATGEN DEFECTS ANALYSIS ########################')
print(find_defects_multi_files(poscar_path, contcar_paths)@gan_lattice_vecs)
print(fp_sep_dist(poscar_path, contcar_paths, atom_type='ga', atom_no=34, lattice_vecs=gan_lattice_vecs))

###########################################################################
############################ OVITO WS ANALYSIS ############################
###########################################################################
"""
how this needs to work
can probably just use a single function for everything, don't really have much repetition

WS analysis: output_displaced = False
- returns the vacant site location and the other atom closest to final interstitial
- will not know which is the vacant site yet
- can get the vacancy and interstitial count info here

WS analysis: output_displaced = True
- returns interstitial site location and other atom closest to final interstitial
- site that matches previous WS analysis can be thrown away, other site from first is v and other site from second is i
- need to choose a cutoff such that the additional site is either thrown away or considered as a dumbbell
- need for returning proper displacement vectors, so should be done second
- also need site type/index/identifier properties from this

Displacements modifier
- find displacement vectors, should be able to use the site index from second WS to determine final position of interstitial
- use this and site type/identifier to determine what type of interstitial and what position

final info of desire
v & i count (get directly from first WS)
v & i type: Ga/N (get directly from second WS?)
i site: t/o/db (disp mod and something else idk yet)
if v & i type are same (they should be, i think only case where they aren't is if there are multiple defects): separation distance (both WS and calculate distance formula)


main issues currently:
- problem regarding cutoff separating atoms inhabiting the lattice site closest to the interstitial being regarded as so vs. split dumbbell interstitial
- could be a problem with defining the array masks to delete rows (should they be the same? or could they be in different orders)
- this trickles down into the FP separation distance calculation
- need to use the displacement magnitude/vectors to determine whether interstitials are tetrahedral or octahedral
- also potential issues with runs like the [120] direction (17L?) where ovito recognizes 2 FPs (1 Ga, 1 N) as 1 FP
- need to bring in other info such as direction associated with pseudo (should be easy, use key_dict again)
"""


def determine_defect_properties(base_pos_path, cont_paths, r_cut=0.5):
    """
    Given an initial POSCAR filepath and an array of CONTCAR filepaths corresponding to structures
    with defects, yields information on the defects including the vacancy and interstitial count,
    atom type corresponding to each defect, what type of interstitial site(s) are occupied (if
    applicable), and the separation distance for Frenkel pairs. POSCAR filepath should be a
    string and the CONTCAR filepaths should be a 1D array of strings. A keyword argument (float)
    is defined in Angstroms for the cutoff distance to differentiate near-interstitial atoms as
    remaining on-site or forming dumbbell configurations.
    """
    defects_dict = {}

    for i in range(cont_paths.shape[0]):
        pipeline_dict = {}
        tde_dir = str(re.split(r'/', cont_paths[i])[-4])
        tde_energy = int(re.split(r'/', cont_paths[i])[-3][:-2])
        pipeline_dict['E_{d}'] = tde_energy

        # Import file into pipeline
        pipeline = import_file(cont_paths[i])

        # Determine displacement of interstitial atom
        disp = CalculateDisplacementsModifier(
            affine_mapping = ov.pipeline.ReferenceConfigurationModifier.AffineMapping.ToReference)
        disp.reference = FileSource()
        disp.reference.load(poscar_path)
        pipeline.modifiers.append(disp)

        data = pipeline.compute()
        disp_vecs, disp_mags = data.particles['Displacement'][...], data.particles['Displacement Magnitude'][...]
        # print(disp_mags)

        # Define a modifier function that identifies vacancies and interstitials.
        def get_defects_modifier(frame, data):
            """
            Modifier function for Ovito to identify vacancies and interstitials.
            """

            # Retrieve the two-dimensional array with the site occupancy numbers.
            # Use [...] to cast it to a Numpy array.
            occupancies = data.particles['Occupancy']

            # Calculate total occupancy of every site:
            total_occupancy = np.sum(occupancies, axis=1)

            # Set up a particle selection by creating the Selection property:
            selection = data.particles_.create_property('Selection')

            # Select A-sites occupied by exactly one B-atom (the second entry of the Occupancy
            # array must be 1, and all others 0). Note that the Occupancy array uses 0-based
            # indexing, while atom type IDs are typically 1-based.
            selection[...] = (total_occupancy == 1)  # (occupancies[:,0] == 1) & (occupancies[:,1] == 1) & 

            # Remove all other atoms.
            pipeline.modifiers.append(DeleteSelectedModifier())


        # Perform Wigner-Seitz analysis with output mode "Sites"
        ws_sites = WignerSeitzAnalysisModifier(
            output_displaced = False,
            per_type_occupancies = True,
            affine_mapping = ov.pipeline.ReferenceConfigurationModifier.AffineMapping.ToReference,
            minimum_image_convention = True)
        ws_sites.reference = FileSource()
        ws_sites.reference.load(base_pos_path)
        pipeline.modifiers.append(ws_sites)

        data = pipeline.compute()
        v_cnt, i_cnt = data.attributes['WignerSeitz.vacancy_count'], data.attributes['WignerSeitz.interstitial_count']
        pipeline_dict['v_cnt'], pipeline_dict['i_cnt'] = v_cnt, i_cnt
        
        print(i, v_cnt, i_cnt)
        
        if v_cnt == 0 and i_cnt == 0:
            continue
        elif tde_dir in defects_dict.keys():
            continue

        # Insert Python modifier into the data pipeline to remove non-defect atoms/sites
        pipeline.modifiers.append(get_defects_modifier)
        data = pipeline.compute()

        # Define output file path variables
        sites_xyz_path, atoms_xyz_path = os.path.join(defects_path, f'{tde_dir}_sites.xyz'), os.path.join(defects_path, f'{tde_dir}_atoms.xyz')

        # Export the XYZ coordinates of the defect sites
        export_file(pipeline, sites_xyz_path, 'xyz',
            columns = ['Position.X', 'Position.Y', 'Position.Z'],
            multiple_frames = False)

        # Disable Winger-Seitz analysis "Sites" output
        # pipeline.modifiers[0].enabled = False

        # Perform Wigner-Seitz analysis with output mode "Atoms"
        ws_atoms = WignerSeitzAnalysisModifier(
            output_displaced = True,
            per_type_occupancies = True,
            affine_mapping = ov.pipeline.ReferenceConfigurationModifier.AffineMapping.ToReference,
            minimum_image_convention = True)
        ws_atoms.reference = FileSource()
        ws_atoms.reference.load(base_pos_path)
        pipeline.modifiers[1] = ws_atoms

        # Export the XYZ coordinates of the displaced atom positions
        export_file(pipeline, atoms_xyz_path, 'xyz',
            columns = ['Position.X', 'Position.Y', 'Position.Z'],
            multiple_frames = False)

        # Gather defect positions for sites and displaced atoms
        sites_xyz, atoms_xyz = np.genfromtxt(sites_xyz_path, delimiter=' ', skip_header=2), np.genfromtxt(atoms_xyz_path, delimiter=' ', skip_header=2)
        pipeline_dict['sites_xyz'], pipeline_dict['atoms_xyz'] = sites_xyz, atoms_xyz
        
        print('sites and atoms '+str(i), sites_xyz, atoms_xyz, sep='\n')
        
        # Identify extra atom associated with site and site left vacant
        for n in range(sites_xyz.shape[0]):
            for m in range(atoms_xyz.shape[0]):
                print(sites_xyz[n, :], atoms_xyz[m, :])
                r_nm = math.dist(sites_xyz[n, :], atoms_xyz[m, :])
                print(r_nm)
                if r_nm <= r_cut:
                    sites_mask, atoms_mask = np.ones(sites_xyz.shape, dtype=bool), np.ones(atoms_xyz.shape, dtype=bool)
                    sites_mask[n, :], atoms_mask[m, :] = False, False
                    i_pos, v_pos = sites_xyz[sites_mask,...], atoms_xyz[atoms_mask,...]
                elif r_nm > r_cut:
                    i_pos, v_pos = sites_xyz, atoms_xyz

        if len(i_pos.shape) == 1 and len(v_pos.shape) == 1:
            i_pos, v_pos = np.reshape(i_pos, (1, i_pos.shape[0])), np.reshape(v_pos, (1, v_pos.shape[0]))
        elif len(i_pos.shape) == 1 and len(v_pos.shape) >= 1:
            i_pos = np.reshape(i_pos, (1, i_pos.shape[0]))
        elif len(i_pos.shape) >= 1 and len(v_pos.shape) == 1:
            v_pos = np.reshape(v_pos, (1, v_pos.shape[0]))
        # potential issues: should use same mask for both arrays? unsure if order could change --> can resolve by writing the rest and then comparing to hand calculation in ovito
        # not sure how to treat the i and v sites if there isn't something found --> what to do if dumbbell and what to do if multiple FPs

        # Calculate distance for Frenkel pair separation
        d_sep = np.zeros((i_pos.shape[0], 1))

        for t in range(i_pos.shape[0]):
            d_sep[t] = math.dist(i_pos[t, :], v_pos[t, :])
        pipeline_dict['d_sep'] = d_sep

        # Determine defect site types and the displacement vector/magnitude
        data = pipeline.compute()
        site_type, site_index = data.particles['Site Type'][...], data.particles['Site Index'][...]
        pipeline_dict['site_type'], pipeline_dict['site_index'] = site_type, site_index
        pipeline_dict['disp_vecs'], pipeline_dict['disp_mags'] = disp_vecs[site_index[0]], disp_mags[site_index[0]]

        defects_dict[tde_dir] = pipeline_dict
    
    return defects_dict

print('############################ OVITO WS ANALYSIS ############################')
tde_defects_dict = determine_defect_properties(poscar_path, contcar_paths, r_cut=0.5)
print(tde_defects_dict)

defect_ovito_file_path = os.path.join(defects_path, 'defect_ovito.txt')
defect_file = open(defect_ovito_file_path, 'w+')
defect_file.write(pd.DataFrame(data=tde_defects_dict).to_string())
defect_file.close()

def parse_defect_properties(defects_dict, atom_types=['Ga', 'N']):
    """
    Function to analyze the defect properties for each calculation contained
    in a dictionary.
    """
    # for each contcar
    # if defects are present (first check)
    # organize type of defect and number of defects
    # also give frenkel pair separation
    # maybe make like a table type thing? not sure if that's necessary but it could be good ig

    # create initial DataFrame to emulate table
    defects_df = pd.DataFrame(data=defects_dict)
    
    # remove unnecessary info for the purposes of this analysis
    try:
        defects_df = defects_df.drop(['sites_xyz', 'atoms_xyz'])
    except:
        print('No sites_xyz or atoms_xyz in df')
    
    # add empty rows to Dataframe for additional data
    defects_df.loc['Defect Type'] = pd.Series([None for n in range(defects_df.columns.shape[0])])

    # iterate through columns (calculation directions)
    for i in range(defects_df.columns.shape[0]):
        try:
            site_type_arr = defects_df[defects_df.columns[i]].loc['site_type']
            if np.all(site_type_arr == site_type_arr[0]):
                defects_df[defects_df.columns[i]].loc['site_type'] = site_type_arr[0]
        except:
            print('No site_type in df')

        try:
            site_index_arr = defects_df[defects_df.columns[i]].loc['site_index']
            if np.all(site_index_arr == site_index_arr[0]):
                defects_df[defects_df.columns[i]].loc['site_index'] = site_index_arr[0]
        except:
            print('No site_index in df')

        try:
            d_sep_arr = defects_df[defects_df.columns[i]].loc['d_sep']
            if d_sep_arr.shape == (1, 1):
                defects_df[defects_df.columns[i]].loc['d_sep'] = d_sep_arr[0, 0]
        except:
            print('No d_sep in df')
        
        vacancies_list = ['V'+'_{'+atom_types[defects_df[defects_df.columns[i]].loc['site_type']-1]+'}' for v in range(defects_df[defects_df.columns[i]].loc['v_cnt'])]
        interstitials_list = [atom_types[defects_df[defects_df.columns[i]].loc['site_type']-1]+'_{i}' for t in range(defects_df[defects_df.columns[i]].loc['i_cnt'])]
        defects_df[defects_df.columns[i]].loc['Defect Type'] = ' + '.join(vacancies_list+interstitials_list)
        # could determine O or T by having a list of possible O and T sites, determine distance from each via displacement vector, choose lowest distance

    # remove unnecessary info used in analysis for the purposes of displaying
    try:
        defects_df = defects_df.drop(['i_cnt', 'v_cnt'])
    except:
        print('No i_cnt or v_cnt in df')
    
    try:
        defects_df = defects_df.drop('site_type')
    except:
        print('No site_type in df')
    
    try:
        defects_df = defects_df.drop('site_index')
    except:
        print('No site_index in df')
    
    # sort DataFrame
    defects_df.sort_index(inplace=True)

    # would be good to include information such as the knockout atom and the Ed value

    return defects_df


print(parse_defect_properties(tde_defects_dict).T)

defect_table_file_path = os.path.join(defects_path, 'defect_table.txt')
defect_file = open(defect_table_file_path, 'w+')
defect_file.write(parse_defect_properties(tde_defects_dict).T.to_string())
defect_file.close()
# parse_defect_properties(tde_defects_dict).T.to_csv(defect_info_file_path, sep=',')
