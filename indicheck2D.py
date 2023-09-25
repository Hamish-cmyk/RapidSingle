import sys, os, math
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import dihedral
from MDAnalysis.core.topologyobjects import Bond
from MDAnalysis.core.topologyobjects import Dihedral
import numpy as np
import pandas as pd
from sklearn.metrics import euclidean_distances
from matplotlib import pyplot as plt, colors
from scipy.stats import wasserstein_distance

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

def beauty2d(empty_diheds,empty_rhos):
    fig, ax = plt.subplots()
    ##############
    minv  = 1.0
    maxv  = 100
    nbins = 60
    ##############
    if 'N' in interest:
        hh = ax.hist2d(empty_diheds, empty_rhos, bins=nbins, norm=colors.LogNorm(vmin=minv, vmax=maxv), cmap=plt.cm.binary)
    else:
        hh = ax.hist2d(empty_diheds, empty_rhos, bins=nbins, norm=colors.LogNorm(vmin=minv, vmax=maxv), cmap=plt.cm.Blues)
    ax.set_facecolor('antiquewhite')
    fig.colorbar(hh[3], ax=ax)
    ###########################################################
    plt.ylim(0.4, 2.2)
    plt.xlim(-180,180)
    plt.yticks(np.arange(0.4, 2.3,  step=0.2),fontsize=8)
    plt.xticks(np.arange(-180, 181,  step=30),fontsize=8)
    plt.savefig(root + f'/{subjects[i]}/{subjects[i]}_{states[j]}.png', dpi=100)

def testplot2():
    if not os.path.exists(interest+f'diheds.csv'):
        #############################################
        if testdir(f'{interest}NPTm.part0001fix.xtc') != True:
            cmd1 = f'echo System | gmx_mpi trjconv -f {interest}NPTm.part0001.xtc -pbc mol -o {interest}NPTm.part0001fix.xtc -s {interest}NPTm.tpr'
            os.system(cmd1)
        #############################################
        xtc = f'{interest}NPTm.part0001fix.xtc'
        gro = f'{interest}NPTm.part0001.gro'
        ###################################################################
        traj = mda.Universe(gro,xtc)                                        ## names error means xtg,gro order is wrong
        # #################################################################
        #empty_diheds, empty_rhos = np.ones((4001, 1)), np.ones((4001, 1))
        empty_diheds, empty_rhos = np.ones((5001, 1)), np.ones((5001, 1))
        for ts in range(0, len(traj.trajectory)):
            print(f'parsing frame {ts}')
            traj.trajectory[ts]
            reduced2a = traj.select_atoms('resname F000 or resname NF or resname FES0 or resname NFE or resname PHE and name CG or name CZ')
            count = 0
            for mols in range(0, 1):
                id3   = np.array((reduced2a.indices[1 + count], reduced2a.indices[0 + count], reduced2a.indices[2 + count], reduced2a.indices[3 + count]))
                ###################################################
                dihed = Dihedral(id3,traj).dihedral(pbc=True)
                empty_diheds[ts, mols] = dihed
                ###################################################
                id1 = np.array((reduced2a.indices[0 + count], reduced2a.indices[2 + count]))
                id2 = np.array((reduced2a.indices[1 + count], reduced2a.indices[3 + count]))
                ipso = Bond(id1, traj).length(pbc=True)
                para = Bond(id2, traj).length(pbc=True)
                rho = ipso / para
                ###################################################
                empty_rhos[ts, mols] = rho
                count += 4
        np.savetxt(interest+f'diheds.csv', empty_diheds, delimiter=",", fmt='%s')
        np.savetxt(interest+f'rhos.csv', empty_rhos, delimiter=",", fmt='%s')
        #####################################
        empty_diheds = empty_diheds.flatten()
        empty_rhos = empty_rhos.flatten()
        beauty2d(empty_diheds, empty_rhos)
    else:
        pass

def WeCompare():
    for i in range(0, len(subjects)):
        interest = root + f'/{subjects[i]}/CC/'
        cc = np.genfromtxt(interest+f'empty_diheds.csv', delimiter=',')
        for j in range(0, len(states)):
            interest  = root + f'/{subjects[i]}/{states[j]}/'
            case      = np.genfromtxt(interest + f'empty_diheds.csv', delimiter=',')
            distance1 = wasserstein_distance(cc[:,3], case[:,3])
            print(distance1)

########################################################################################################################
root = '/users/fjb15167/LSAtomistic'
# subjects = ['Nf-Nn-Nf','Nfes-Nke-Nfes','Nfes-Nn-Nfes','Nfes-Nq-Nfes','Nfe-Nke-Nfe']
subjects = ['Ac-Nf-Nf','Nf-Nk-Nf','Nf-Nke-Nf','Nk-Nf-Nf','Nke-Nf-Nf']
states   = ['CC','CT','TC','TT']

test     = False
explore  = True                                 ## measure frequency of 30 degree lambda swings from trajectory
EMD      = False                                ## measure WS metric of dihedral array of sequence from FF reference
stats    = False                                ## extract occupancy of alpha, beta and gamma domains

if test == True:
    subjects = ['Nf-Nk-Nf']
    states = ['CC', 'CT', 'TC', 'TT']
    # subjects = ['Nf-Na-Nf', 'Nfe-Nke-Nfe', 'Nf-Nn-Nf', 'Nfes-Nke-Nfes', 'Nfes-Nn-Nfes', 'Nfes-Nq-Nfes']
    # states = ['CC', 'CT', 'TC', 'TT']
else:
    pass
    # subjects  = ['KFF','DFF','FKF','FDF','FF']
    # states    = ['TT']
    # subjects = ['Nf-Nke-Nf', 'Nf-Nk-Nf', 'Nke-Nf-Nf', 'Nk-Nf-Nf']
    # states = ['CC', 'CT', 'TC', 'TT']

SwitchResults = pd.DataFrame(columns=states)
for i in range(0,len(subjects)):
    print('_____________________________________________________________________________________________________________\n')
    print(subjects[i])
    if stats == True:
        f = open(subjects[i]+'/summary.txt','w')
        f.write(subjects[i]+'\n\n')
    for j in range(0,len(states)):
        interest = root + f'/{subjects[i]}/{states[j]}/'
        testplot2()
        ##########################
        empty_diheds = np.genfromtxt(interest+f'diheds.csv', delimiter=",").flatten()
        empty_rhos   = np.genfromtxt(interest+f'rhos.csv', delimiter=",").flatten()
        beauty2d(empty_diheds, empty_rhos)
        plt.close()
        ##########################
        alpha   = 0
        beta    = 0
        gamma   = 0
        ##########################
        rho   = np.genfromtxt(interest + f'rhos.csv', delimiter=',').flatten()
        lambd = np.genfromtxt(interest + f'diheds.csv', delimiter=',').flatten()
        if EMD == True:
            diff_rho = wasserstein_distance(ref_rho,rho)
            diff_lam = wasserstein_distance(ref_lambd,lambd)
            print(diff_rho,diff_lam,interest)
        ##########################
        frames    = 4000
        molecules = 1
        plt.title(interest)
        ##########################
        plt.hist(lambd[:],bins=12, histtype='step')
        plt.savefig(root + f'/{subjects[i]}/self_{subjects[i]}_{states[j]}.png', dpi=100)
        plt.close()
        ##########################
        delta  = []
        switch = 0
        if explore == True:
            lambd = np.genfromtxt(interest + f'diheds.csv', delimiter=',')
            print(lambd.shape)
            for kol in range(0,4000):
                next = kol + 1
                if next <= 4000:
                    ##########################################
                    if lambd[kol] < 0:
                        xi = lambd[kol] + 360
                    else:
                        xi = lambd[kol]
                    ##########################################
                    if lambd[next] < 0:
                        xi1 = lambd[next] + 360
                    else:
                        xi1 = lambd[next]
                    ##########################################
                    distance = xi - xi1
                    ##########################################
                    if (distance**2)**0.5 >= 20:
                        switch += 1
                else:
                    pass
                delta.append(switch)
            if j == 0:
                combi = np.array(delta)
            else:
                combi = np.vstack((combi,delta))

        if explore == True:
            print('Update @ '+root + f'/{subjects[i]}/switches360.csv')
            np.savetxt(root + f'/{subjects[i]}/switches360.csv',combi, delimiter=",")

    if explore == True:
        lams       = np.genfromtxt(root + f'/{subjects[i]}/switches360.csv', delimiter=",")
        state_data = []
        for each in range(0,4):
            state_data.append(np.mean(lams[each,:]))
            print(state_data)
        state_data    = np.array(state_data).reshape(1,4)
        lamSwitch     = pd.DataFrame(state_data, columns=states, index=[f'{subjects[i]}'])
        SwitchResults = pd.concat([SwitchResults, lamSwitch])

if stats == True:
    f.close()

if explore == True:
    SwitchResults.to_csv(root + '/SwitchResults.csv')

sys.exit()