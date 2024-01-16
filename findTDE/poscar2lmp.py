#!/usr/bin/env python3
### Convert POSCAR to lammps read_data file####
import numpy as np
import os
import sys
import itertools


def pbc_distance(r1,r2,dimensions):
    # periodic distance given two points; dimensions is a numpy array, 
    delta = np.abs(r1 - r2)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))


#https://docs.lammps.org/Howto_triclinic.html
def convert_basis(A,B,C):
    ax = np.linalg.norm(A)
    bx = np.dot(B,A/ax)
    by = np.linalg.norm(np.cross(A/ax,B))
    cx = np.dot(C,A/ax)
    cy = (np.dot(B,C)-bx*cx)/by
    cz = np.sqrt(np.dot(C,C)-cx**2-cy**2)
    return np.array([[ax,bx,cx],[0,by,cy],[0,0,cz]])


class poscar2lmp:
    def __init__(self,posfile,outfile):
        self.read_pos(posfile)
        self.write_data(outfile)

    def read_pos(self,posfile):
        f = open(posfile,'r')
        headerlines = f.readlines()

        elemsline = 5
        dimlines = [2,3,4]
        atomsline = elemsline + 1
        headerline = 8
        if all(map(str.isdigit,headerlines[elemsline].strip().split()[0])): # if no element names exist
            atomsline = elemsline 
            headerline = 7

        scaling = float(headerlines[1])
        self.dims = []
        self.los = np.zeros(3)
        self.his = np.zeros(3)
        for i,item in enumerate(dimlines):
            lin = headerlines[item].strip().split()
            self.dims.append([float(k) for k in lin])
        self.dims = np.array(self.dims)
        atom_counts = [int(i) for i in headerlines[atomsline].strip().split()]
        self.natom = np.sum(atom_counts)
        f.close()

        data = [] #position data 
        k = headerline
        for line in headerlines[headerline:]:
             lin = line.strip().split()
             k += 1
             if len(lin) == 0:
                 continue
             data.append([float(i) for i in lin[:3]])
             if len(data) == self.natom:
                 break
        data = np.array(data)
        if headerlines[headerline-1].strip().split()[0] in ['direct','Direct','D']:
            data = np.dot(data,self.dims)*scaling
        data_v = [] #velocity data
        for line in headerlines[k:]:
             lin = line.strip().split()
             if len(lin) == 0:
                 continue
             data_v.append([float(i) for i in lin[:3]])
             if len(data_v) == self.natom:
                 break
        data_v = np.array(data_v)*1000  # scale to Ang/ps
        
        A,B,C = self.dims 
        Vol = np.linalg.norm(np.dot(np.cross(A,B),C))
        abc = convert_basis(A,B,C)
        factor = 1/Vol*np.dot(abc,np.array([np.cross(B,C),np.cross(C,A),np.cross(A,B)]))
        pos = np.dot(data,factor.transpose())
        vel = np.dot(data_v,factor.transpose())
        self.his[0] = abc[0,0] 
        self.his[1] = abc[1,1]
        self.his[2] = abc[2,2]
        if abs(abc[0,1]) <= abs((self.his[0] - self.los[0])/2.):
            self.xy_xz_yz = [abc[0,1],abc[0,2],abc[1,2]]
        elif abs(abc[0,1]) > abs((self.his[0] - self.los[0])/2.) and abc[0,1] < 0.:
            abc[0,1] += abs(self.his[0] - self.los[0])
            self.xy_xz_yz = [abc[0,1],abc[0,2],abc[1,2]]
        elif abs(abc[0,1]) > abs((self.his[0] - self.los[0])/2.) and abc[0,1] > 0.:
            abc[0,1] -= abs(self.his[0] - self.los[0])
            self.xy_xz_yz = [abc[0,1],abc[0,2],abc[1,2]]


        atomids = np.arange(1,self.natom+1)
        atomtypes = np.concatenate([np.ones(item)*(i+1) for i,item in enumerate(atom_counts)])
        self.atom_types = len(atom_counts)
        self.data = np.hstack((atomids.reshape(-1,1),atomtypes.reshape(-1,1),pos))
        self.data_v = np.hstack((atomids.reshape(-1,1),vel))


    def write_data(self,outfile):
        # write to readdata for lammps
            f = open(outfile,'w')
            f.write('Generate lammps data file from poscar\n\n')
            f.write('%d  atoms\n' %(self.natom))
            f.write('%d  %s\n\n' %(self.atom_types,'atom types'))
            f.write('%.8e  %.8e  %s\n' %(self.los[0],self.his[0],'xlo xhi'))
            f.write('%.8e  %.8e  %s\n' %(self.los[1],self.his[1],'ylo yhi'))
            f.write('%.8e  %.8e  %s\n' %(self.los[2],self.his[2],'zlo zhi'))
            f.write('%.8e  %.8e  %.8e xy xz yz\n\n' %(tuple(self.xy_xz_yz)))
            f.write('%s\n\n' %('Atoms'))
            for d in self.data:
                f.write('%d %d %.8e %.8e %.8e\n' %(tuple(d)))

            f.write('\nVelocities\n\n')
            for d in self.data_v:
                f.write('%d %.8e %.8e %.8e\n' %(tuple(d)))
            f.close()


###### IF XY IS TOO LARGE, ADD XHI TO IT. SIMILIAR TREATMENT FOR XZYZ I THINK
if __name__ == '__main__':
    posfile = os.path.join(os.getcwd(), 'POSCAR')
    outfile = os.path.join(os.getcwd(), 'read_data.lmp')
    poscar2lmp(posfile,outfile)
    
