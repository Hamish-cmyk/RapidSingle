import sys, os, math
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import dihedral
from MDAnalysis.core.topologyobjects import Dihedral
from matplotlib import pyplot as plt, colors
import numpy as np
import pandas as pd
from ase.io import read
from copy import copy
from MDAnalysis.core.topologyobjects import Bond
from MDAnalysis.core.topologyobjects import Dihedral
from ase import Atoms

def readin(fname):
    f = open(fname, 'r')
    content = f.read()
    return content

def writeout(fname, content):
    f = open(fname, 'w')
    f.write(content)
    f.close()

def testdir(path):
    kill_test = False
    if os.path.exists(path):
        kill_test = True
    return kill_test

def SpecialIndex():
    ##################################### make original index
    cmd1 = f"echo q | gmx_mpi make_ndx -f i{name}.gro -o index.ndx"
    os.system(cmd1)
    original = readin("index.ndx")
    ##################################### select relevant IDs using MDA
    u = mda.Universe(f'i{name}.pdb')
    solvent = u.select_atoms('resname TIP3 or resname CLA')
    protein = u.select_atoms(f'resname KE00 or resname FES0 or resname NN or resname NQ or resname NFE or resname NKE or resname NF or resname NA')
    ##################################### write atom IDs to a new list as a loop
    Solute = []
    Solvent = []
    for atom in protein:
        Solute.append(atom.ix + 1)
    for mol in solvent:
        Solvent.append(mol.ix + 1)
    print("Peptoid:", len(Solute))
    print("Solvent:", len(Solvent))
    ##################################### write special IDs to a new index file
    index = open("SpecialIndex.ndx", 'w')
    index.write(original)
    index.write("[ Peptoid ]\n")
    i = 0
    for ndx in Solute:
        index.write(str(ndx) + " ")
        i += 1
        if i % 15 == 0:
            index.write("\n")
    index.write("\n")
    ############################
    index.write("[ Solvent ]\n")
    for ndx in Solvent:
        index.write(str(ndx) + " ")
        i += 1
        if i % 15 == 0:
            index.write("\n")
    index.write("\n\n")
    index.close()

def FixYouanSH():
    reads = readin('/users/fjb15167/LSAtomistic/run_0.sh')
    with open(f'{each}/run.sh','w') as sh:
        sh.write(reads.replace('YYYY',f'{each}'))

#####################################################################################

test =  True
root = '/users/fjb15167/LSAtomistic'

if test == True:
    sequences = ['Nf-Nk-Nf', 'Nf-Nke-Nf', 'Nk-Nf-Nf', 'Nke-Nf-Nf']
    ions      = [2,2,2,2]
    ##########################
    # sequences = ['Nf-Na-Nf']
    # ions      = [1]
    # states    = ['CC', 'CT', 'TC', 'TT']
    # sequences = ['Nfe-Nke-Nfe']
    # ions      = [2]
    # states    = ['CC', 'CT', 'TC', 'TT']
else:
    pass
    # sequences = ['Nf-Nn-Nf', 'Nfe-Nke-Nfe', 'Nfes-Nke-Nfes', 'Nfes-Nn-Nfes', 'Nfes-Nq-Nfes']
    # ions      = [1,2,2,1,1]
    # states    = ['CC', 'CT', 'TC', 'TT']
    # sequences = ['Nfe-Nke-Nfe', 'Nfes-Nke-Nfes', 'Nfes-Nn-Nfes', 'Nfes-Nq-Nfes']
    # sequences = ['Nfes-Nke-Nfes', 'Nfes-Nn-Nfes', 'Nfes-Nq-Nfes']
    # ions      = [2,1,1]
    # states    = ['CC', 'CT', 'TC', 'TT']

for i, each in enumerate(sequences):
    FixYouanSH()
    for j, those in enumerate(states):
        ######################################
        subject = f'{each}/{those}'
        name    = each
        ######################################
        if testdir(subject) != True:
            os.system(f'mkdir {subject}')
        ######################################
        os.chdir(subject)
        ######################################
        box  = '3 3 3'  # !! if you change this then update SpecialIndex function!!
        cmd1 = f'gmx_mpi insert-molecules -ci ../{those}.pdb -nmol 1 -box {box} -o box.pdb'
        cmd2 = f'gmx_mpi insert-molecules -f box.pdb -ci ../../chloride.pdb -nmol {ions[i]} -box {box} -o box.pdb'
        cmd3 = f'gmx_mpi solvate -cp box.pdb -box {box} -o box.pdb'
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        ######################################
        box      = readin('box.pdb')
        boxwrite = box.replace('SOL ', 'TIP3').replace('HW1', 'H1 ').replace('HW2', 'H2 ').replace('OW ', 'OH2')
        with open('box.pdb','w') as box:
            box.write(boxwrite)
        ##################################################################
        os.system('cat box.pdb | grep -v CLA | grep -v TIP3 > peptoid.pdb')
        os.system('cat box.pdb | grep CLA    > ions.pdb')
        os.system('cat box.pdb | grep TIP3   > water.pdb')                      ## pull water...
        psfgen = readin('../../patch_temp.pgn')
        ########################################
        if 'Nfes-' in name:
            patch = 'NTER2'
        elif 'Nf-' in name:
            patch = 'NTER'
        elif 'Nfe-' in name:
            patch = 'NTER3'
        else:
            print('Patch Undefined!')
            sys.exit()
        ########################################
        with open('patch.pgn','w') as pgn:
            pgn.write(psfgen.replace('XXX',name).replace('ACTION',patch))
        os.system('vmd -dispdev text -e patch.pgn > psf.log')
        ########################################
        reader = readin('../../generate_topology.tcl')
        with open('gen.tcl','w') as tcl:
            tcl.write(reader.replace('XXX',name))
        os.system('vmd -dispdev text -e gen.tcl')
        ########################################
        cmd1 = f'gmx_mpi editconf -f i{name}.pdb -o i{name}.gro -box 3 3 3'
        os.system(cmd1)
        SpecialIndex()
        #########################################################################
        os.system(f'gmx_mpi grompp -f ../../min.mdp -c i{name}.gro -p fiber.top -o min.tpr -n SpecialIndex.ndx > grompp_min.log')
        os.system('mpirun -np 4 gmx_mpi mdrun -s min.tpr -v -deffnm min -ntomp 1')
        os.system('gmx_mpi grompp -f ../../NPT.mdp -c min.gro -p fiber.top -o NPT.tpr -n SpecialIndex.ndx > grompp_NPT.log')
        #########################################################################
        os.system(f"rm *#*")
        os.chdir(root)

sys.exit()
