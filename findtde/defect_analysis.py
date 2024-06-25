#!/usr/bin/env python3
"""Python module used to characterize defects resulting from TDE calculations."""
import os
import sys
import glob
import pprint
import fortranformat as ff

import math
import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')

import ovito._extensions.pyscript

from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
from ovito.pipeline import *


def find_contcars(base_path, pseudo='*'):
    """
    Function to create an array of post-CGM CONTCAR files for defect analysis.
    """
    contcar_files = np.array(glob.glob(os.path.join(base_path, pseudo, '*', 'cgm', 'CONTCAR'), recursive=True))
    return contcar_files


def find_dumpfiles(base_path, pseudo='*'):
    """
    Function to create an array of post-MD dump files for defect analysis.
    """
    dump_files = np.array(glob.glob(os.path.join(base_path, pseudo, '*', 'dump.final'), recursive=True))
    return dump_files


def pseudo_keys_from_file(pseudo_keyfile='pseudo_keys.csv'):
    """
    Function to convert a CSV file of pseudo keys to a dictionary.
    """
    pk_df = pd.read_csv(pseudo_keyfile, index_col=0)
    pks_dict = pk_df.to_dict('split')
    pks_dict = dict(zip(pks_dict['index'], pks_dict['data']))
    return pks_dict


def PBC_distance(x0, x1, dimensions):
    """
    Source: https://stackoverflow.com/a/11109336
    Function that takes 3 numpy arrays: an initial position array, a final position array,
    and the dimensions of the boundary. Returns an array for the distance, taking into
    account Periodic Boundary Conditions (PBC).
    """
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


def cart2sph(x, y, z):
    dxy = np.sqrt(x**2 + y**2)
    r = np.sqrt(dxy**2 + z**2)
    theta = np.arctan2(y, x)
    phi = np.arctan2(dxy, z)
    theta, phi = np.rad2deg([theta, phi])
    return r, theta % 360, phi


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


def fp_sep_dist(base_pos_path, cont_paths, lattice_vecs, atom_type='ga', atom_no=34):
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

# print('######################## PYMATGEN DEFECTS ANALYSIS ########################')
# print(find_defects_multi_files(poscar_path, contcar_paths)@gan_lattice_vecs)
# print(fp_sep_dist(poscar_path, contcar_paths, gan_lattice_vecs, atom_type='ga', atom_no=34))

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


main issues currently:    <---- this was for the previous version of this file, may not be relevant anymore
- problem regarding cutoff separating atoms inhabiting the lattice site closest to the interstitial being regarded as so vs. split dumbbell interstitial
- could be a problem with defining the array masks to delete rows (should they be the same? or could they be in different orders)
- this trickles down into the FP separation distance calculation
- need to use the displacement magnitude/vectors to determine whether interstitials are tetrahedral or octahedral
- also potential issues with runs like the [120] direction (17L?) where ovito recognizes 2 FPs (1 Ga, 1 N) as 1 FP
- need to bring in other info such as direction associated with pseudo (should be easy, use key_dict again)
"""
def defect_analysis(defect_file, ref_file, verbose=False):
    """
    Perform WS analyses
    """
    ws_comb_dict = {
        'sites': {},
        'atoms': {}
    }

    pipeline = import_file(defect_file)
    # displace atom config
    data = pipeline.compute()
    disp_type = data.particles['Particle Type'] 
    disp_pos = data.particles['Position'][...]

    # Perform Wigner-Seitz analysis (sites mode):
    ws = WignerSeitzAnalysisModifier(
        output_displaced = False,
        per_type_occupancies = True,
        affine_mapping = ReferenceConfigurationModifier.AffineMapping.ToReference,
        minimum_image_convention = True
        )
    
    # Load the reference config from a separate input file.
    ws.reference = FileSource() # Note: You may have to import FileSource from the ovito.pipeline module. 
    ws.reference.load(ref_file)
    pipeline.modifiers.append(ws)
    output_sites = pipeline.compute()


    def get_defects_from_WS(data):
        """
        Parse data gathered from Wigner-Seitz Analysis
        """
        ws_dict = {}

        occupancies = data.particles['Occupancy']
        site_pos = data.particles['Position'][...]

        # Get the site types as additional input:
        site_type = data.particles['Particle Type']
        # Calculate total occupancy of every site:
        total_occupancy = np.sum(occupancies, axis=1)

        ######## get vacancy information ########
        nAvac = np.sum((site_type == 1) & (total_occupancy == 0))
        nBvac = np.sum((site_type == 2) & (total_occupancy == 0))
        nvac = data.attributes['WignerSeitz.vacancy_count']
        print('\nIdenfity %d vacancies, including %d type1 and %d type2'%(nvac,nAvac,nBvac)) if verbose else None
        print('type1 vacancy site positions:') if verbose else None
        Avac_pos = site_pos[(site_type == 1) & (total_occupancy == 0)]
        print(Avac_pos) if verbose else None
        print('type2 vacancy site positions:') if verbose else None
        Bvac_pos = site_pos[(site_type == 2) & (total_occupancy == 0)]
        print(Bvac_pos) if verbose else None
        ws_dict['nAvac'], ws_dict['nBvac'], ws_dict['Avac_pos'], ws_dict['Bvac_pos'] = nAvac, nBvac, Avac_pos, Bvac_pos

        ######## get antisite information ########
        B2A = np.sum((site_type == 1) & (occupancies[:,1] == 1) & (total_occupancy == 1))
        A2B = np.sum((site_type == 2) & (occupancies[:,0] == 1) & (total_occupancy == 1))
        print('\nIdenfity %d antisites, including %d type1occupy2 and %d type2occupy1'%(B2A+A2B,A2B,B2A)) if verbose else None
        print('type2occupy1 site positions:') if verbose else None
        B2A_pos = site_pos[(site_type == 1) & (occupancies[:,1] == 1) & (total_occupancy == 1)]
        print(B2A_pos) if verbose else None
        print('type1occupy2 site positions:') if verbose else None
        A2B_pos = site_pos[(site_type == 2) & (occupancies[:,0] == 1) & (total_occupancy == 1)]
        print(A2B_pos) if verbose else None
        ws_dict['nB2A'], ws_dict['nA2B'], ws_dict['B2A_pos'], ws_dict['A2B_pos'] = B2A, A2B, B2A_pos, A2B_pos

        ######## get interstitial information ########
        int_sites = site_pos[total_occupancy > 1]
        int_site_types = site_type[total_occupancy > 1]
        finder = NearestNeighborFinder(np.max(total_occupancy), data)
        ws_dict['Aint_pos'], ws_dict['Bint_pos'] = [np.empty((0, 3), dtype=np.float64) for n in range(0, 2)]
        ws_dict['nAint'], ws_dict['nBint'] = 0, 0  # np.array([[]])
        for i,pos in enumerate(int_sites):
            print('\nIntersitial at site',pos,'with site type=',int_site_types[i]) if verbose else None
            if int_site_types[i] == 1 and pos not in ws_dict['Aint_pos']:
                ws_dict['Aint_pos'] = np.append(ws_dict['Aint_pos'], np.array([pos]), axis=0)
                ws_dict['nAint'] += 1
            elif int_site_types[i] == 2 and pos not in ws_dict['Bint_pos']:
                ws_dict['Bint_pos'] = np.append(ws_dict['Bint_pos'], np.array([pos]), axis=0)
                ws_dict['nBint'] += 1
            
            """ # determine nearest neighbors for each atom
            print('Atoms at this site')
            for neigh in finder.find_at(pos):
                ind = neigh.index
                print(disp_type[ind],disp_pos[ind])"""

        return ws_dict

    print('######## Wigner-Seitz Analysis (Sites) ########') if verbose else None
    ws_comb_dict['sites'] = get_defects_from_WS(output_sites)

    # Update Wigner-Seitz analysis (atoms mode):
    pipeline.modifiers[0].output_displaced = True
    output_atoms = pipeline.compute()

    print('######## Wigner-Seitz Analysis (Atoms) ########') if verbose else None
    ws_comb_dict['atoms'] = get_defects_from_WS(output_atoms)

    pprint.pprint(ws_comb_dict) if verbose else None

    # Compare both Wigner-Seitz analyses, return cross-analyzed data
    defect_dict = {
        'nAvac': int(ws_comb_dict['sites']['nAvac']),
        'nBvac': int(ws_comb_dict['sites']['nBvac']),
        'nA2B': int(ws_comb_dict['sites']['nA2B']),
        'nB2A': int(ws_comb_dict['sites']['nB2A']),
        'nAint': int(ws_comb_dict['atoms']['nAint']) - int(ws_comb_dict['sites']['nAint']),
        'nBint': int(ws_comb_dict['atoms']['nBint']) - int(ws_comb_dict['sites']['nBint']),
        'Avac_pos': ws_comb_dict['sites']['Avac_pos'],
        'Bvac_pos': ws_comb_dict['sites']['Bvac_pos'],
        'A2B_pos': ws_comb_dict['sites']['A2B_pos'],
        'B2A_pos': ws_comb_dict['sites']['B2A_pos'],
        'Aint_pos': ws_comb_dict['atoms']['Aint_pos'],  # np.empty((0, 3), dtype=np.float64),
        'Bint_pos': ws_comb_dict['atoms']['Bint_pos']  # np.empty((0, 3), dtype=np.float64)
    }

    # something is going wrong, here's some pseudocode
    """
    maybe set original lattice site either way and add to dict

    for each interstitial (use counts); need to
    1. determine distance from interstitial and displaced atom to original lattice site
        interstitial and displaced atom are both in defect_dict
        can both be in Aint_pos/Bint_pos, or have one in each
        original lattice site is in ws_comb_dict['sites']
        will be in same Aint_pos/Bint_pos based on location of displaced atom
    2a. if int&disp are in the same array, will need to determine which is interstitial based on which one is closer to original lattice site
    2b. if int&disp are in different arrays, maybe just treat the same?
    3. after determining int&disp, int will remain in interstitial list and disp will be categorized as...
        either nothing (and will be removed), if distance from site is < ~0.5-0.75 Ang
        or dumbbell config, if site is roughly shared, distance from site is ~ 0.75-1.75 Ang
    4. remove displaced atom either way (can therefore treat this the same, just adding dumbbell to dictionary either way)
        need to remove disp either way and add dumbbell to dict either way
        will just add a value if a dumbbell exists, will add empty if not
    5. use parsing function to add dumbbell into config if exists
    """

    r_i, r_ij = [], []
    og_A_latt_pos, disp_A_fin_pos, db_A_fin_pos = [np.empty((0, 3), dtype=np.float64) for n in range(0, 3)]
    og_B_latt_pos, disp_B_fin_pos, db_B_fin_pos = [np.empty((0, 3), dtype=np.float64) for n in range(0, 3)]
    for iA in range(defect_dict['nAint']):
        #og_A_latt_pos = og_A_latt_pos.append(ws_comb_dict['sites']['Aint_pos'][iA])
        pass
    for iB in range(defect_dict['nBint']):
        pass

    """r_i, r_ij = [], []
    disp_atom_fin_pos, db_atom_fin_pos = [], []
    if defect_dict['nAint'] < defect_dict['Aint_pos'].shape[0]:
        if defect_dict['nAint'] == 0:
            # Aint_pos contains only displaced atom, Bint_pos contains interstitial
            for i in range(defect_dict['Aint_pos'].shape[0]):
                r_i.append(defect_dict['Aint_pos'][i])
                for j in range(ws_comb_dict['sites']['Aint_pos'].shape[0]):
                    r_ij.append(PBC_distance(defect_dict['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j], gan_params))
            # need to find index of minimum r_ij value corresponding to first set of r_i values
            # aka need to slice r_ij list based on [(i*ws_comb_dict['sites']['Aint_pos'].shape[0]):(i+1)*ws_comb_dict['sites']['Aint_pos'].shape[0])]
            for i in range(len(r_i)):
                j_start = i*ws_comb_dict['sites']['Aint_pos'].shape[0]
                j_end = (i + 1)*ws_comb_dict['sites']['Aint_pos'].shape[0]
                j_min = np.argmin(r_ij[j_start:j_end]) + j_start
                if r_ij[j_min] <= 0.75:  # small displacement
                    disp_atom_fin_pos.append(r_i[i])
                elif r_ij[j_min] > 0.75 and r_ij[j_min] <= 1.75:  # dumbbell config
                    db_atom_fin_pos.append(r_i[i])
                else:
                    pass
                defect_dict['Aint_pos'] = np.delete(defect_dict['Aint_pos'], r_i[i])
        else:
            # Aint_pos contains both displaced atom and interstitial
            # need to break in case of interstitial
            for i in range(defect_dict['Aint_pos'].shape[0]):
                r_i.append(defect_dict['Aint_pos'][i])
                for j in range(ws_comb_dict['sites']['Aint_pos'].shape[0]):
                    r_ij.append(PBC_distance(defect_dict['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j], gan_params))
            for i in range(len(r_i)):
                j_start = i*ws_comb_dict['sites']['Aint_pos'].shape[0]
                j_end = (i + 1)*ws_comb_dict['sites']['Aint_pos'].shape[0]
                j_min = np.argmin(r_ij[j_start:j_end]) + j_start
                if r_ij[j_min] <= 0.75:  # small displacement
                    disp_atom_fin_pos.append(r_i[i])
                elif r_ij[j_min] > 0.75 and r_ij[j_min] <= 1.75:  # dumbbell config
                    db_atom_fin_pos.append(r_i[i])
                else:
                    pass
                defect_dict['Aint_pos'] = np.delete(defect_dict['Aint_pos'], r_i[i])
    elif defect_dict['nBint'] < defect_dict['Bint_pos'].shape[0]:
        if defect_dict['nBint'] == 0:
            # Bint_pos contains only displaced atom, Aint_pos contains interstitial
            for i in range(defect_dict['Bint_pos'].shape[0]):
                r_i.append(defect_dict['Bint_pos'][i])
                for j in range(ws_comb_dict['sites']['Bint_pos'].shape[0]):
                    r_ij.append(PBC_distance(defect_dict['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j], gan_params))
            for i in range(len(r_i)):
                j_start = i*ws_comb_dict['sites']['Bint_pos'].shape[0]
                j_end = (i + 1)*ws_comb_dict['sites']['Bint_pos'].shape[0]
                j_min = np.argmin(r_ij[j_start:j_end]) + j_start
                if r_ij[j_min] <= 0.75:  # small displacement
                    disp_atom_fin_pos.append(r_i[i])
                elif r_ij[j_min] > 0.75 and r_ij[j_min] <= 1.75:  # dumbbell config
                    db_atom_fin_pos.append(r_i[i])
                else:
                    pass
                defect_dict['Bint_pos'] = np.delete(defect_dict['Bint_pos'], r_i[i])
        else:
            # Bint_pos contains both displaced atom and interstitial
            for i in range(defect_dict['Bint_pos'].shape[0]):
                r_i.append(defect_dict['Bint_pos'][i])
                for j in range(ws_comb_dict['sites']['Bint_pos'].shape[0]):
                    r_ij.append(PBC_distance(defect_dict['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j], gan_params))
            for i in range(len(r_i)):
                j_start = i*ws_comb_dict['sites']['Bint_pos'].shape[0]
                j_end = (i + 1)*ws_comb_dict['sites']['Bint_pos'].shape[0]
                j_min = np.argmin(r_ij[j_start:j_end]) + j_start
                if r_ij[j_min] <= 0.75:  # small displacement
                    disp_atom_fin_pos.append(r_i[i])
                elif r_ij[j_min] > 0.75 and r_ij[j_min] <= 1.75:  # dumbbell config
                    db_atom_fin_pos.append(r_i[i])
                else:
                    pass
                defect_dict['Bint_pos'] = np.delete(defect_dict['Bint_pos'], r_i[i])"""

    # can use now accurate Aint & Bint values
    # if either Aint or Bint array is longer than the nAint or nBint, from array that is longer, determine which value is closest to original site position
    # set this closer position to a different variable, remove from int_pos array
    # if distance is > 0.5ish, add dumbbell position to dictionary
    """lattice_site_org_pos, disp_atom_fin_pos, db_atom_fin_pos = [], [], []
    for i in range(ws_comb_dict['atoms']['Aint_pos'].shape[0]):
        for j in range(ws_comb_dict['sites']['Aint_pos'].shape[0]):
            print(ws_comb_dict['atoms']['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j])
            print(PBC_distance(ws_comb_dict['atoms']['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j], gan_params))

            if (ws_comb_dict['atoms']['Aint_pos'][i] == ws_comb_dict['sites']['Aint_pos'][j]).all():
                lattice_site_org_pos.append(ws_comb_dict['atoms']['Aint_pos'][i])
            elif PBC_distance(ws_comb_dict['atoms']['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j], gan_params) <= 0.75 and PBC_distance(ws_comb_dict['atoms']['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j], gan_params) > 0.:
                disp_atom_fin_pos.append(ws_comb_dict['sites']['Aint_pos'][i])
            elif PBC_distance(ws_comb_dict['atoms']['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j], gan_params) <= 1.0 and PBC_distance(ws_comb_dict['atoms']['Aint_pos'][i], ws_comb_dict['sites']['Aint_pos'][j], gan_params) > 0.75:
                db_atom_fin_pos.append(ws_comb_dict['sites']['Aint_pos'][i])
            else:
                defect_dict['Aint_pos'] = ws_comb_dict['atoms']['Aint_pos']
    
    for i in range(ws_comb_dict['atoms']['Bint_pos'].shape[0]):
        for j in range(ws_comb_dict['sites']['Bint_pos'].shape[0]):
            print(ws_comb_dict['atoms']['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j])
            print(PBC_distance(ws_comb_dict['atoms']['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j], gan_params))

            if (ws_comb_dict['atoms']['Bint_pos'][i] == ws_comb_dict['sites']['Bint_pos'][j]).all():
                lattice_site_org_pos.append(ws_comb_dict['atoms']['Bint_pos'][i])
            elif PBC_distance(ws_comb_dict['atoms']['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j], gan_params) <= 0.75 and PBC_distance(ws_comb_dict['atoms']['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j], gan_params) > 0.:
                disp_atom_fin_pos.append(ws_comb_dict['sites']['Bint_pos'][i])
            elif PBC_distance(ws_comb_dict['atoms']['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j], gan_params) <= 1.0 and PBC_distance(ws_comb_dict['atoms']['Bint_pos'][i], ws_comb_dict['sites']['Bint_pos'][j], gan_params) > 0.75:
                db_atom_fin_pos.append(ws_comb_dict['sites']['Bint_pos'][i])
            else:
                defect_dict['Bint_pos'] = ws_comb_dict['atoms']['Bint_pos']
    
    defect_dict['nAint'], defect_dict['nBint'] = defect_dict['Aint_pos'].shape[0], defect_dict['Bint_pos'].shape[0]"""

    #print(og_atom_latt_pos, disp_atom_fin_pos, db_atom_fin_pos)

    return defect_dict


######################################
# function to determine TDE from a single datafile
######################################
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


######################################
# function to average TDE values for each atom type for a given direction range
# function currently assumes rad-induced defect config is Frenkel pair
######################################
def average_tde_vals(defects_dict, pseudo_keys, lattice_vecs, dir_range=((0, 360), (0, 180))):
    """
    Function to average the TDE values for each atom type over the given range of directions.
    """
    tde_averages = {
        'avg': {
            'E_d_avg_A': [],
            'E_d_avg_B': []
        },
        'std': {
            'E_d_std_A': [],
            'E_d_std_B': []
        },
        'weights': {
            'E_d_wi_A': [],
            'E_d_wi_B': []
        }
    }
    
    cnt_A_FP, cnt_B_FP, cnt_complex, cnt_oor, cnt_cutoff = 0, 0, 0, 0, 0
    # need variable dtheta, dphi
    dtheta, dphi = np.deg2rad(7.5), np.deg2rad(7.5)

    for key in defects_dict.keys():
        # check the direction and pass if outside of range
        if key[-1] == 'S':
            phi, theta = pseudo_keys[key][0], pseudo_keys[key][1]
        elif key[-1] == 'L':
            phi, theta = lat2sph(np.array([pseudo_keys[key]]), lattice_vecs)[0][1], lat2sph(np.array([pseudo_keys[key]]), lattice_vecs)[0][2]
        else:
            raise ValueError('Incorrect direction pseudo.')
        
        if (phi < dir_range[0][0]) or (phi > dir_range[0][1]) or (theta < dir_range[1][0]) or (theta > dir_range[1][1]):
            cnt_oor += 1
            continue

        # check the resultant defect type
        if (defects_dict[key]['nAint'] > 0) and (defects_dict[key]['nAvac'] > 0):
            defect_type = 'A'
            cnt_A_FP += 1
        elif (defects_dict[key]['nBint'] > 0) and (defects_dict[key]['nBvac'] > 0):
            defect_type = 'B'
            cnt_B_FP += 1
        else:
            #print('{} --> Other type: nAint = {}, nAvac = {}, nBint = {}, nBvac = {}'.format(key, defects_dict[key]['nAint'], defects_dict[key]['nAvac'], defects_dict[key]['nBint'], defects_dict[key]['nBvac']))
            #print(defects_dict[key])
            cnt_complex += 1
            continue
        
        # categorize TDE values by the defect type
        tde_averages['avg'][f'E_d_avg_{defect_type}'].append(defects_dict[key]['E_d'])
        tde_averages['weights'][f'E_d_wi_{defect_type}'].append(np.sin(np.deg2rad(theta))*dtheta*dphi)
    
    # find standard deviation of the TDE values for each defect type
    for key in tde_averages['std'].keys():
        tde_averages['std'][key] = np.std(tde_averages['avg'][key.replace('std', 'avg')], dtype=np.float64)
    
    # average the TDE values for each defect type
    for key in tde_averages['avg'].keys():
        # standard average
        # tde_averages['avg'][key] = np.mean(tde_averages['avg'][key], dtype=np.float64)
        
        # weighted average
        tde_averages['avg'][key] = np.average(tde_averages['avg'][key], weights=tde_averages['weights'][key.replace('avg', 'wi')])
    del tde_averages['weights']
    
    print('A FPs:', cnt_A_FP, '\tB FPs:', cnt_B_FP, '\tOther:', cnt_complex, '\tOut of range:', cnt_oor)

    return tde_averages


######################################
# can then take the code i started to make previously to analyze this info
# maybe need to change how it is stored first in the original function
######################################
def parse_defect_properties(defects_dict, atom_types=['Ga', 'N']):
    """
    Function to analyze the defect properties for each calculation contained
    in a dictionary.
    """

    # create initial DataFrame to emulate table
    defects_df = pd.DataFrame(data=defects_dict)
    
    # add empty rows to Dataframe for additional data
    defects_df.loc['Defect Type'] = pd.Series([None for n in range(defects_df.columns.shape[0])])
    defects_df.loc['d_sep'] = pd.Series([None for n in range(defects_df.columns.shape[0])])

    # iterate through columns (calculation directions)
    for i in range(defects_df.columns.shape[0]):
        vacancies_list = ['V_{}'.format(atom_types[0]) for v in range(defects_df[defects_df.columns[i]].loc['nAvac'])] + \
            ['V_{}'.format(atom_types[1]) for v in range(defects_df[defects_df.columns[i]].loc['nBvac'])]
        antisites_list = ['{}_{}'.format(atom_types) for a in range(defects_df[defects_df.columns[i]].loc['nA2B'])] + \
            ['{}_{}'.format(atom_types[::-1]) for a in range(defects_df[defects_df.columns[i]].loc['nB2A'])]
        interstitials_list = ['{}_i'.format(atom_types[0]) for t in range(defects_df[defects_df.columns[i]].loc['nAint'])] + \
            ['{}_i'.format(atom_types[1]) for t in range(defects_df[defects_df.columns[i]].loc['nBint'])]
        defects_df[defects_df.columns[i]].loc['Defect Type'] = ' + '.join(vacancies_list+antisites_list+interstitials_list)
        # could determine O or T by having a list of possible O and T sites, determine distance from each via displacement vector, choose lowest distance

        if defects_df[defects_df.columns[i]].loc['nAvac'] == 1 and defects_df[defects_df.columns[i]].loc['nAint'] == 1:
            defects_df[defects_df.columns[i]].loc['d_sep'] = PBC_distance(defects_df[defects_df.columns[i]].loc['Aint_pos'], defects_df[defects_df.columns[i]].loc['Avac_pos'], gan_params)
        elif defects_df[defects_df.columns[i]].loc['nBvac'] == 1 and defects_df[defects_df.columns[i]].loc['nBint'] == 1:
            defects_df[defects_df.columns[i]].loc['d_sep'] = PBC_distance(defects_df[defects_df.columns[i]].loc['Bint_pos'], defects_df[defects_df.columns[i]].loc['Bvac_pos'], gan_params)

    # remove unnecessary info used in analysis for the purposes of displaying
    defects_df = defects_df.drop(['nAvac', 'nBvac', 'nA2B', 'nB2A', 'nAint', 'nBint'])
    
    # sort DataFrame
    defects_df.sort_index(inplace=True)

    # would be good to include information such as the knockout atom and the Ed value

    return defects_df


#####command: ovitos xx.py
if __name__ == "__main__":
    FILES_LOC = os.getcwd()    # os.path.dirname(__file__)
    
    run_type, knockout_atom = 'vasp', 'ga34'

    if run_type == 'vasp':
        defect_files = find_contcars(FILES_LOC, pseudo=f'*{knockout_atom}*')
    elif run_type == 'lammps':
        defect_files = find_dumpfiles(FILES_LOC, pseudo=f'*{knockout_atom}*')    # test with 10*S_ga34*
    # pprint.pprint(defect_files)
    
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

    # add analysis information for each direction to a dict
    defect_configs_dict = {}

    for i in range(defect_files.shape[0]):
        if run_type == 'vasp':
            calc_desc, calc_tde = defect_files[i].split('/')[-4], int(defect_files[i].split('/')[-3][:-2])
            calc_pseudo, calc_atom_info = calc_desc.split('_')[0], calc_desc.split('_')[1]
            
            # get TDE value for direction from data file
            pseudo_datafile = os.path.abspath(os.path.join(defect_files[i], os.pardir, os.pardir, os.pardir, calc_pseudo+'_'+calc_atom_info+'_data.csv'))
            pseudo_tde_value = get_tde_from_datafile(pseudo_datafile, e_tol=1.0, ke_cut=45)
        elif run_type == 'lammps':
            calc_desc, calc_tde = defect_files[i].split('/')[-3], int(defect_files[i].split('/')[-2][:-2])
            calc_pseudo, calc_atom_info = calc_desc.split('_')[0], calc_desc.split('_')[1]
            
            # get TDE value for direction from data file
            pseudo_datafile = os.path.abspath(os.path.join(defect_files[i], os.pardir, os.pardir, calc_pseudo+'_'+calc_atom_info+'_data.csv'))
            pseudo_tde_value = get_tde_from_datafile(pseudo_datafile, e_tol=1.0, ke_cut=100)
        
        # print(calc_desc, ': ', calc_tde, ' eV\t (TDE: ', pseudo_tde_value, ' eV)' sep='')
        
        # perform defect analysis
        if calc_tde == pseudo_tde_value:
            defect_configs_dict[calc_pseudo] = defect_analysis(defect_files[i], vasp_ref_file)
            defect_configs_dict[calc_pseudo]['E_d'] = pseudo_tde_value
        elif calc_tde != pseudo_tde_value:
            # potentially need to add cutoff handling here
            if pseudo_tde_value >= 45:
                print('cutoff needed')
            continue
        else:
            continue

    print('############################ OVITO WS ANALYSIS SUMMARY ############################')
    # pprint.pprint(defect_configs_dict)
    if run_type == 'vasp':
        vasp_pseudo_keys = pseudo_keys_from_file(pseudo_keyfile=os.path.join(FILES_LOC, 'gan_vasp_pseudo_keys.csv'))
    elif run_type == 'lammps':
        lmp_pseudo_keys = pseudo_keys_from_file(pseudo_keyfile=os.path.join(FILES_LOC, 'gan_lmp_pseudo_keys.csv'))
    pprint.pprint(average_tde_vals(defect_configs_dict, vasp_pseudo_keys, vasp_lattice_vecs, dir_range=((30, 150), (0, 180))))

    # defect_configs_df = parse_defect_properties(defect_configs_dict, atom_types=['Ga', 'N']).T
    # pprint.pprint(defect_configs_df)
