import sys, os, math
import MDAnalysis as mda
from MDAnalysis.core.topologyobjects import Bond
from MDAnalysis.core.topologyobjects import Dihedral
import numpy as np
import pandas as pd
from sklearn.metrics import euclidean_distances
from matplotlib import pyplot as plt, colors

def readin(fname):
    f = open(fname, 'r')
    content = f.read()
    return content

def writeout(fname, content):
    f = open(fname, 'w')
    f.write(content)
    f.close()

def CheckDihedrals():

    dcd       = 'NPT.part0001fix.trr'
    psf       = 'NPT.part0001.gro'
    traj      = mda.Universe(psf, dcd)
    selection = 'resname MOA0 and name N or name CB or name CG or name OD'
    specific  = np.split(traj.select_atoms(selection).ix,2)

    for j, each in enumerate(specific):
        complete = []
        counter = 0
        for ts in traj.trajectory:
            dihed = Dihedral(each, traj).dihedral(pbc=True)
            print(dihed,counter,each)
            counter +=1
            complete.append(dihed)
        if j == 0:
            both = np.array(complete)
        else:
            both = np.vstack((both,complete))

    ## TEST #############################
    #sub = both[1,:] - np.array(complete)
    #print(sub,both,both.shape)
    #####################################

    for i, k in enumerate(both):
        plt.hist(k,label=f'residue {i+1}',histtype='step',linewidth=1.75)
    plt.xlim(-180,180)
    plt.xticks(np.arange(-180, 181, step=30), fontsize=8)
    plt.legend()
    plt.ylabel('Population')
    plt.xlabel('Angle ($^\circ$)')
    plt.savefig('SidechainTorsions.png', dpi=100)
    plt.show()
    plt.close()


CheckDihedrals()
