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

def rdf_solute_ow(TPR="../c20_0p_s0_iter2.md_PME.tpr", XTC="../water-mediated/c20_0p_s0_iter2.md_PME.dt-1000.xtc",
                   grp1="resname SDS NGLU NGUC NGUO NGLP NGPO NGPC", grp2="name OW",
                   ti=-2, tf=-1, skip=1, excl_block=(1,1)):
    """
    Plot RDF for waters
    """
    u = Universe(TPR, XTC)

    grp1 = u.select_atoms(grp1)
    grp2 = u.select_atoms(grp2)

    rdf = InterRDF(grp2, grp2, nbins=150, range=[2, 15], excl_block=(1,1))
    rdf_result = rdf.run(start=ti, stop=tf, step=skip, verbose=True)
    
    result = np.transpose(np.vstack((rdf.bins, rdf.rdf)))

    return result

def ow_indices(u, ts, solute_ow_cut=3.5):


    """
    Returns list of OW indices for waters that are within cutoff from solute
    """
    solute = u.select_atoms("resname SDS NGLU NGUC NGUO NGLP NGPO NGPC and name S1 O11 O12 O13 O1 O4 N2 O3 O5 O6")
    waters = u.select_atoms(f"name OW and around 12 resname SDS NGLU NGUC NGUO NGLP NGPO NGPC and name S1 O11 O12 O13 O1 O4 N2 O3 O5 O6")

    dist = distance_array(solute.positions, waters.positions, box=ts.dimensions) 

    # convert to dataframe to keep track of atom indices
    dist = pd.DataFrame(dist, columns=[a.index for a in waters], index = [a.index for a in solute])

    
    dist = dist.where(dist<solute_ow_cut, 0)    # set far away waters as 0
    dist = dist.where(dist==0, 1)               # set waters within cutoff as 1
    ow_indices = list(dist.any()[dist.any() == 1].index) # returns column names that contain `1`, aka, water ids in contact with solute

    return ow_indices

def ow_neighbours(u, ts, ow_indices, ow_ow_cutoff = 3):

    """
    Created sorted neighbor list of waters in contact with solute, let's call them O
    For each O, 8 closest waters are identifeid and stored
    First 4 nearest neighbours are called O1
    5th-8th nearest neighbours are called O2
    Output is a dataframe of indices of O, with sorted list of their 8 nearest neighbours (O1O2)
    """

    temp = " ".join(map(str,ow_indices))
    o = u.select_atoms(f"index {temp}")
    waters = u.select_atoms("name OW and around 12 resname SDS NGLU NGUC NGUO NGLP NGPO NGPC and name S1 O11 O12 O13 O1 O4 N2 O3 O5 O6")

    dist = distance_array(waters.positions, o.positions, box=ts.dimensions)
    dist = pd.DataFrame(dist, index=[a.index for a in waters], columns = [a.index for a in o])
    dist = dist.where(dist<ow_ow_cutoff)       # set far away waters as NaN
    
    o1o2_dict = {}
    # parse columns, sort and store 8 nearest neighbor
    for label, content in dist.items():
        # parse columns
        neighbours = tuple(content.sort_values()[1:9].index) # omit self by skipping the closest atom
        o1o2_dict[label] = neighbours

    o1o2 = pd.DataFrame.from_dict(o1o2_dict, orient='columns')

    return o1o2

def run_neighbour_search(TPR="../c20_0p_s0_iter2.md_PME.tpr", XTC="../water-mediated/c20_0p_s0_iter2.md_PME.dt-1000.xtc",
                   grp1="name OW", grp2="name OW",
                   ti=-2, tf=-1, skip=1, excl_block=(1,1)):
    """
    Plot RDF for waters
    """
    u = Universe(TPR, XTC)
    for ts in u.trajectory[ti:tf:skip]:

        o = ow_indices(u, ts) # get shell waters
        #print(o)
        o1o2 = ow_neighbours(u,ts, o)
    
    return o1o2

def run_neighbour_search_BULK_WATER(TPR="../spc.md.tpr", XTC="../spc.md.dt-1000.xtc",
                   grp1="name OW", grp2="name OW",
                   ti=-2, tf=-1, skip=1, excl_block=(1,1)):
    """
    Plot RDF for waters
    """
    u = Universe(TPR, XTC)
    for ts in u.trajectory[ti:tf:skip]:

        o = u.select_atoms("name OW").indices
        #o = ow_indices(u, ts) # get shell waters
        #print(o)
        o1o2 = ow_neighbours(u,ts, o)

    return o1o2

tis = list(range(-14, -6))
tfs = [ti+1 for ti in tis]


# for ti,tf in zip(tis, tfs): 
#     x = run_neighbour_search_BULK_WATER(ti=ti, tf=tf)
#     x.to_csv(f'spc.o1o2{tf}.csv')

# for ti,tf in zip(tis, tfs): 
    
#     for sys_name in ['c20_0p_s0_iter2', 'c20_100p_s512_iter2' ,'c20_0p_s512']:

#         x = run_neighbour_search(TPR=f"../{sys_name}.md_PME.tpr", XTC=f"../water-mediated/{sys_name}.md_PME.dt-1000.xtc",
#                             ti=ti, tf=tf)

#         x.to_csv(f'{sys_name}.o1o2{tf}.csv')


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


fig, ax = plt.subplots(nrows=1, ncols=2, sharex=True, sharey=True, dpi=300, figsize=(7,7/3))
colors = sns.color_palette('colorblind')

for sys_name in ['c20_0p_s0_iter2', 'c20_0p_s512', 'c20_100p_s512_iter2'][:]:

    x = rdf_solute_ow(TPR=f"../{sys_name}.md_PME.tpr", XTC=f"../water-mediated/{sys_name}.md_PME.dt-1000.xtc",
                    grp1="resname SDS NGLU NGUC NGUO NGLP NGPO NGPC and name S1 O11 O12 O13 O1 O4 N2 O3 O5 O6",
                    grp2="name OW",
                    ti=-2, tf=-1)

    ax[0].plot(x[:,0]/10., x[:,1])

    y = rdf_solute_ow(TPR=f"../{sys_name}.md_PME.tpr", XTC=f"../water-mediated/{sys_name}.md_PME.dt-1000.xtc",
                    grp1="resname SDS NGLU NGUC NGUO NGLP NGPO NGPC and name S1 O11 O12 O13 O1 O4 N2 O3 O5 O6",
                    grp2="name HW1",
                    ti=-2, tf=-1)

    ax[1].plot(y[:,0]/10., y[:,1])

    #np.savetxt(f"{sys_name}.rdf.txt", x)
ax[0].legend(['NB-0', 'NB-180', 'NA-180'])
#plt.title("SO4-HW rdf")
ax[0].set_ylabel(r"$g(r)_{OW}$")
ax[1].set_ylabel(r"$g(r)_{HW}$")
ax[0].set_xlabel(r"$r$ (nm)")
ax[1].set_xlabel(r"$r$ (nm)")
ax[0].grid()
ax[1].grid()
plt.tight_layout()
plt.savefig("solute-water.rdf.png")
plt.show()