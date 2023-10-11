import sys, os, math
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import angle
from MDAnalysis.lib.mdamath import dihedral
from MDAnalysis.core.topologyobjects import Bond
from MDAnalysis.core.topologyobjects import Dihedral
import numpy as np
from matplotlib import pyplot as plt, colors
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import wasserstein_distance

def readin(fname):
    f = open(fname, 'r')
    content = f.read()
    return content

def writeout(fname, content):
    f = open(fname, 'w')
    f.write(content)
    f.close()

def GiveCOM():
    print(f'Parsing {path}')
    file1    = path + 'NPTm.part0001.gro'
    file2    = path + 'NPTm.part0001fix.xtc'
    traj     = mda.Universe(file1,file2)
    duration = len(traj.trajectory)
    empty    = np.ones((duration))
    for i in range(0,duration):
        print(f'reading frame {i}')
        traj.trajectory[i]
        backbone      = traj.select_atoms('name CA or name NP or name C or name O or name N')
        # backbone      = traj.select_atoms('name CA')                                                                  ## testing (see Nf-Na-Nf/CC)
        residue_ids   = np.unique(backbone.resids)
        points        = []
        for each in residue_ids:
            selection = traj.select_atoms(f'name CA or name NP or name C or name O or name N and resid {each}')
            # selection = traj.select_atoms(f'name CA and resid {each}')                                                ## testing (see Nf-Na-Nf/CC)
            com       = selection.center_of_mass()
            points.append(com)
        vec1          = points[0] - points[1]
        vec2          = points[2] - points[1]
        total_angle   = round(math.degrees(angle(vec1, vec2)), 3)
        empty[i]      = total_angle
    if state == 'CC':
        np.savetxt(path+'AA_total_angle.csv', empty, delimiter=',', fmt='%s')
        sys.exit()
    return empty

def SpecialIndex(species):
    residues    = np.unique(species.split('-'))
    resname_str = ' or '.join(f'resname {name}' for name in residues).replace('Na', 'Nx')
    #####################################
    cmd = "echo q | gmx_mpi make_ndx -f min.gro -o index.ndx"
    os.system(cmd)
    original = readin("index.ndx")
    #####################################
    u = mda.Universe("min.gro")
    solvent = u.select_atoms('resname W or resname ION or resname WF')
    protein = u.select_atoms(resname_str)
    #####################################
    Solute = []
    Solvent = []
    for atom in protein:
        Solute.append(atom.ix + 1)
    for mol in solvent:
        Solvent.append(mol.ix + 1)
    print("Peptoid:", len(Solute))
    print("Solvent:", len(Solvent))
    #####################################
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
    index.write("[ Solvent ]\n")
    for ndx in Solvent:
        index.write(str(ndx) + " ")
        i += 1
        if i % 15 == 0:
            index.write("\n")
    index.write("\n\n")
    index.close()

def GoToMartinoid(each):
    os.chdir('/users/fjb15167/Martinoid-main/lsaStudy')
    os.system(f'python ../Martinoid2.py {each} Linear')
    os.chdir(root)
    pass

def FixYouaTOP(martinoid,waters):
    reads = readin('/users/fjb15167/LSAtomistic/temp.top')
    with open(f'sys.top','w') as topology:
        topology.write(reads.replace('SPECIFIC',f'{martinoid}.itp').replace('HOWMANY',f'{waters}'))

def FixYouanSH(martinoid):
    reads = readin('/users/fjb15167/LSAtomistic/CGrun.sh')
    with open(f'runCG.sh','w') as SH:
        SH.write(reads.replace('XYZ',f'{martinoid}'))

def CGSetUP():
    dir = subject+'/sysCG'
    if not os.path.exists(dir):
        os.system(f'mkdir {dir}')
    # if not os.path.exists(dir+'/eq1.gro'):
    #################################################################
    martinoid   = subject.replace('Nfes','Nfex').replace('Nfer','Nfex')
    water     = '/users/fjb15167/ResSASA/water_4A.gro'
    cgmdp     = '/users/fjb15167/ResSASA/cg_min.mdp'
    GoToMartinoid(martinoid)
    #################################################################
    os.chdir(dir)
    os.system(f'gmx_mpi insert-molecules  -f {water} -ci /users/fjb15167/Martinoid-main/lsaStudy/{martinoid}.pdb -nmol 1 -box 4 4 4 -o solvate.gro -replace W')
    read   = mda.Universe('solvate.gro')
    waters = read.select_atoms('resname W').positions
    waters = str(waters.shape[0])
    #########################
    FixYouaTOP(martinoid,waters)
    #########################
    cmd1 = f'gmx_mpi grompp -f {cgmdp} -c solvate.gro -p sys.top -o min.tpr -maxwarn 1'
    cmd2 = f'echo W | gmx_mpi genion -s min.tpr -p sys.top -nname CL -pname NA -neutral -o solvate.gro'
    cmd3 = 'mpirun -np 4 gmx_mpi mdrun -tableb ../../angle5_a0.xvg -v -deffnm min -ntomp 1'
    ########################
    os.system(cmd1)  # make tpr
    os.system(cmd2)  # add ions
    os.system(cmd1)  # regen tpr
    os.system(cmd3)  # minimize
    ########################
    SpecialIndex(martinoid)
    os.system(f"rm *#*")
    ########################
    if martinoid in longer:
        name = 'eq1more.mdp'
    else:
        name = 'eq1.mdp'
    cmd4 = f'gmx_mpi grompp -f ../../{name} -c min.gro -p sys.top -o eq1.tpr -n SpecialIndex.ndx'
    os.system(cmd4)
    ########################
    FixYouanSH(martinoid)
    os.system('sbatch runCG.sh')
    os.chdir(root)

def dihedral_iterator(add):
    i, j, k, l = 0+add, 1+add, 2+add, 3+add
    dihedrals = []
    for each in range(0,len(CoreGroup)):
        if l < len(CoreGroup):
            v1 = CoreGroup[i] - CoreGroup[j]
            v2 = CoreGroup[j] - CoreGroup[k]
            v3 = CoreGroup[k] - CoreGroup[l]
            ######################################
            dihed = math.degrees(dihedral(v1, v2, v3))
            dihedrals.append(dihed)
            #print(i, j, k, l,dihed)
            i += 3
            j += 3
            k += 3
            l += 3
        else:
            break
    return dihedrals

def CGAngle():
    dir = subject+'/sysCG'
    os.chdir(dir)
    ##################################################################################
    cmd = 'echo System | gmx_mpi trjconv -f eq1.xtc -s eq1.tpr -pbc mol -o eq1fix.xtc'
    os.system(cmd)
    os.system(f"rm *#*")
    ##################################################################################
    file1    = 'eq1.gro'
    file2    = 'eq1fix.xtc'
    traj     = mda.Universe(file1, file2)
    duration = len(traj.trajectory)
    empty    = np.ones((duration))
    for i in range(0, duration):
        print(f'reading frame {i}')
        traj.trajectory[i]
        selection     = traj.select_atoms('name BB').positions
        vec1          = selection[0,:] - selection[1,:]
        vec2          = selection[2,:] - selection[1,:]
        total_angle   = round(math.degrees(angle(vec1, vec2)), 3)
        empty[i]      = total_angle
    #np.savetxt('total_angle.csv',empty,delimiter=',', fmt='%s')
    os.chdir(root)
    return empty


###################################################################################################################################
root      = os.getcwd()
subject   = 'Nfes-Nke-Nfes'
states    = ['CC','CT','TC','TT']
col       = {'CC': 'orangered', 'CT': 'deepskyblue', 'TC': 'deeppink', 'TT': 'darkorange'}

## Controls ##
CGRUN        = False
Measure      = True
TorsionCheck = False
##############

subjects = ['Nf-Na-Nf','Nf-Nk-Nf','Nf-Nke-Nf','Nf-Nn-Nf','Nfe-Nke-Nfe','Nfes-Nke-Nfes','Nfes-Nn-Nfes','Nfes-Nq-Nfes','Nk-Nf-Nf','Nke-Nf-Nf']
longer   = ['Nf-Nk-Nf','Nf-Nke-Nf','Nk-Nf-Nf','Nke-Nf-Nf']
#######################
subjects = ['Nf-Na-Nf']
#######################

if TorsionCheck == True:
    for k, subject in enumerate(subjects):
        fig = plt.figure()
        for l, state in enumerate(states):
            file1 = f'{subject}/{state}/NPTm.part0001.gro'
            file2 = f'{subject}/{state}/NPTm.part0001fix.xtc'
            traj = mda.Universe(file1, file2)
            duration = len(traj.trajectory)
            print(duration)
            for each in range(0, len(traj.trajectory)):
                print(f'Parsing Frame: {each}')
                traj.trajectory[each]
                CoreGroup = traj.select_atoms('name N or name CA or name C or name NP').positions
                omega     = dihedral_iterator(1)
                if each == 0:
                    omega_set = np.array(omega)
                else:
                    omega_set = np.vstack((omega_set, omega))
            print(omega_set,omega_set.shape)
            colour = col.get(state)
            ax     = fig.add_subplot(2,2,l+1)
            ax.set_title(f'{state}', loc="left", pad=-14)
            ax.hist(omega_set[:,0], histtype='step', bins=30, label=f'res1', edgecolor=f'r')
            ax.hist(omega_set[:,1], histtype='step', bins=30, label=f'res2', edgecolor=f'b')
            ax.set_xlim(-180,180)
            ax.set_ylim(0,2000)
            ax.set_yticks(np.arange(0, 1601, step=400))
            ###############################
            legend = ax.legend(fontsize=9,loc='upper center')
            frame  = legend.get_frame()
            frame.set_linewidth(0)
            ###############################
        fig.suptitle(f'{subject}')
        fig.supylabel('Population')
        fig.supxlabel('Angle ($^\circ$)')
        fig.tight_layout()
        plt.savefig(f"{subject}/{subject}_omega.jpg")
        plt.close()

for k,subject in enumerate(subjects):
    if CGRUN == True:
        CGSetUP()
if Measure == True:
    ######################################################################### CG-MD
    for k, subject in enumerate(subjects):
        CGangles = CGAngle()
    ######################################################################### AA-MD
    for k, subject in enumerate(subjects):
        for j,state in enumerate(states):
            path = subject + f'/{state}/'
            arr  = GiveCOM()
            if j == 0:
                combined = arr
            else:
                combined = np.vstack((combined,arr))
        ##################################################################### plot AA-MD
        fig = plt.figure()
        for j,state in enumerate(states):
            colour = col.get(state)
            ax     = fig.add_subplot(2,2, j + 1)
            ax.set_xlim(0, 180)
            ax.set_ylim(0, 800)
            ax.set_title(f'{state}', y=1.0, pad=-14)
            ax.hist(CGangles, histtype='step', bins=30, label=f'CG-MD', edgecolor=f'b')
            ax.hist(combined[j,:], histtype='step',bins=30,label=f'{state}',edgecolor=f'{colour}')
            ax.set_xticks(np.arange(0,181,step=30))
        fig.supylabel('Population')
        fig.supxlabel('Angle ($^\circ$)')
        fig.tight_layout()
        plt.savefig(f"{subject}/{subject}_backbone.jpg")
        plt.close()










