import re
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis import *
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.lib.distances import self_distance_array
from MDAnalysis.core.topologyobjects import Angle, Dihedral
from tqdm import tqdm
import seaborn as sns
from matplotlib.colors import to_rgba
from MDAnalysis.analysis.rdf import InterRDF
from tqdm import tqdm
import pandas as pd
import seaborn as sns

def read_neighbour_lists(sys_name="c20_0p_s0_iter2", tf=-1):

    # nlist is a dataframe of water neighbour lists
    # column id = index of water closest to solute (defined as O)
    # content of the column = indices of 8 nearest waters to O, sorted from closest to farthest
    # note this list is for a single time frame

    nlist = pd.read_csv(f"{sys_name}.o1o2{tf}.csv", index_col=0)
    print(nlist.info())
    print(nlist.head())

    return nlist


def compute_order(nlist, TPR="../c20_0p_s0_iter2.md_PME.tpr", XTC="../water-mediated/c20_0p_s0_iter2.md_PME.dt-1000.xtc",
                   grp1="name OW", grp2="name OW",
                   ti=-2, tf=-1, skip=1, excl_block=(1,1)):

    # IMPORTANT: use the same timeframe for this analysis as the neighbour list

    u = Universe(TPR, XTC)

    q = []

    for ts in u.trajectory[ti:tf:skip]:


        for label, content in nlist.items():
            # iterate each column
            # original waters/waters that touch the solute
            oi = u.select_atoms(f"index {label}")
            temp_o1 = " ".join(map(str, content[0:4]))
            o1 = u.select_atoms(f"index {temp_o1}")

            # v1 = o1.positions[0,:] - oi.positions[0,:]
            # v2 = o1.positions[1,:] - oi.positions[0,:]
            # dp = np.dot(v1/np.linalg.norm(v1), v2/np.linalg.norm(v2))
            # angle = np.arccos(dp)
            # print(angle)

            cos_psi_jk = []
            j=0
            while j<3:
                k=j+1
                while k<4:
                    #v1 = o1.positions[j,:] - oi.positions[0,:]
                    #v2 = o1.positions[k,:] - oi.positions[0,:]
                    #dp = np.dot(v1/np.linalg.norm(v1), v2/np.linalg.norm(v2))
                    #cos_psi_jk.append(dp)

                    x = Angle([o1.indices[j], oi.indices[0], o1.indices[k]], u)
                    #print(x.angle(pbc=True))
                    cos_psi_jk.append( np.cos( x.angle(pbc=True) * np.pi/180. ) )
                    k+=1
                j+=1
            
            #print(cos_psi_jk)
            cos_terms = [ (c + 1/3. )**2 for c in cos_psi_jk ]

            qi = 1 - (3./8.)*np.sum(cos_terms)

            q.append(qi)

    return q
        
# for sys_name in ['spc'][:]:

#     q_dict = {}
#     for ti, tf in zip(range(-14,-2), range(-13,-1)):
#         nlist = read_neighbour_lists(sys_name=f"{sys_name}", tf=tf)
#         q = compute_order(nlist,  TPR=f"../{sys_name}.md.tpr", XTC=f"../{sys_name}.md.dt-1000.xtc",
#                                 ti=ti, tf=tf)
#         q_dict[tf] = q
#     q_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in q_dict.items() ])) # This way uneven columns don't raise an error
#     q_df.to_csv(f"{sys_name}.q.csv", index=False, header=False)


# for sys_name in ['c20_0p_s0_iter2', 'c20_100p_s512_iter2', 'c20_0p_s512']: # , 'spc']:#[:]:

#     q_dict = {}
#     for ti, tf in zip(range(-14,-2), range(-13,-1)):
#         nlist = read_neighbour_lists(sys_name=f"{sys_name}", tf=tf)
#         q = compute_order(nlist, TPR=f"../{sys_name}.md_PME.tpr", XTC=f"../water-mediated/{sys_name}.md_PME.dt-1000.xtc",
#                                 ti=ti, tf=tf)
#         q_dict[tf] = q
#     q_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in q_dict.items() ])) # This way uneven columns don't raise an error
#     q_df.to_csv(f"{sys_name}.q.csv", index=False, header=False)


# for sys_name in ['c20_0p_s512']: # , 'spc']:#[:]:

#     q_dict = {}
#     for ti, tf in zip(range(-6,-2), range(-5,-1)):
#         nlist = read_neighbour_lists(sys_name=f"{sys_name}", tf=tf)
#         q = compute_order(nlist, TPR=f"../{sys_name}.md_PME.tpr", XTC=f"../water-mediated/{sys_name}.md_PME.dt-1000.xtc",
#                                 ti=ti, tf=tf)
#         q_dict[tf] = q
#     q_df = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in q_dict.items() ])) # This way uneven columns don't raise an error
#     q_df.to_csv(f"{sys_name}.q.csv", index=False, header=False)



##### to plot only

sns.set_theme(style="ticks")
palette = sns.color_palette("colorblind")
#darkpalette = sns.color_palette("dark")
print(palette)
sns.set_palette("colorblind")
TINY_SIZE   =   7
SMALL_SIZE  =   7
MED_SIZE    =   10
LARGE_SIZE  =   7

plt.rc('font', size=MED_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MED_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title


fig, ax = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=False, dpi=300, figsize=(7,7/3))
colors = sns.color_palette('colorblind')
means=[]
sem=[]
m2 = []
for i, sys_name in enumerate(['c20_0p_s0_iter2', 'c20_0p_s512', 'c20_100p_s512_iter2', 'spc']):

    q = pd.read_csv(f"{sys_name}.q.csv", header=None)
    #print(q.info())
    q = pd.melt(q.iloc[:,:15])
    #names+=[sys_name]*len(q['value'].tolist())
    #qq+=q['value'].tolist()

    m = q['value'].mean()
    st = q['value'].sem()

    m2.append(q['value'])
    #print(q['value'])
    print(m,st)
    means.append(m)
    sem.append(st)
    #print(q)
    #h, be = np.histogram(q['value'], range=(-3,1), density=True, bins=200)
    
    # plt.plot([(be[i+1]+be[i])/2. for i in range(len(h))],
    #             h, label=f"{sys_name},{m:.2f}", color=colors[i])
    # sns.histplot(q['value'], label=f"{sys_name},{m:.2f}", stat='density', 
    #                 binrange=(-3,1), binwidth=0.05, color=colors[i], element='step', fill=False)
    sns.kdeplot(q['value'], label=f"{sys_name},{m:.2f}", ax=ax[0])   

#q_mega = pd.DataFrame.from_dict({'name': names, 'q': qq})

#sns.kdeplot(data=q_mega, x='q', hue='name', common_grid=False, label=f"{sys_name}", stats="density")

ax[0].set_xlim((-1.5,1))
ax[0].grid()

ax[0].spines['top'].set_visible(False)
ax[0].spines['right'].set_visible(False)


ax[0].legend(['NB-0', 'NB-180', 'NA-180', 'Bulk water'])
ax[0].set_ylabel(r"$P(q)$")
ax[0].set_xlabel(r"$q$")

mean_q = {'name': ['NB-0', 'NB-180', 'NA-180', 'Bulk water'],
            'meanq' : means}

mean_q2 = pd.DataFrame.from_dict( {'NB-0': m2[0],
            'NB-180': m2[1],
            'NA-180': m2[2],
            'Bulk water': m2[3]} )

#mean_q2 = pd.melt(mean_q2)

#q2 = pd.concat(m2)

#sns.barplot(data=mean_q, x='name', y='meanq', ci="sd", capsize=5, ax=ax[1])
#sns.barplot(data=mean_q2, x='variable', y='value', ci=95, capsize=0.5, ax=ax[1])
sns.barplot(data=mean_q2, ci=100, ax=ax[1])
ax[1].set_ylim((0,0.5))
ax[1].grid()

ax[1].spines['top'].set_visible(False)
ax[1].spines['right'].set_visible(False)

ax[1].set_ylabel(r"$\langle q\rangle$")

plt.tight_layout()

#plt.savefig("water-order.svg")
#plt.savefig("water-order.png")
plt.show()
plt.close()
# for sys_name in ['c20_0p_s0_iter2', 'c20_100p_s512_iter2', 'c20_0p_s512', 'spc']:

#     q = pd.read_csv(f"{sys_name}.q.csv", header=None)
#     #print(q.info())
#     #q = pd.melt(q)
#     #m = q['value'].mean()
#     #print(q)
#     #sns.kdeplot(q['value'], label=f"{sys_name},{m:.2f}")
#     sns.kdeplot(data=q)
#     plt.legend([-5, -4, -3, -2, -1])
#     plt.show()
#     plt.close()
# #read_neighbour_lists("c20_0p_s0_iter2.o1o2.csv")
