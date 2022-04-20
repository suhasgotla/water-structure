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

"""
Sript for analysizing the structure of water in different hydrogel network conditions
"""

def rdf_ow(TPR="../c20_0p_s0_iter2.md_PME.tpr", XTC="../water-mediated/c20_0p_s0_iter2.md_PME.dt-1000.xtc",
                   grp1="name OW", grp2="name OW",
                   ti=-2, tf=-1, skip=1, excl_block=(1,1)):
    """
    Plot RDF for waters
    """
    u = Universe(TPR, XTC)

    grp1 = u.select_atoms(grp1)
    grp2 = u.select_atoms(grp2)

    rdf = InterRDF(grp2, grp2, nbins=150, range=[2, 15], excl_block=(1,1))
    rdf_result = rdf.run(start=ti, stop=tf, step=skip, verbose=True)

    plt.plot(rdf.bins, rdf.rdf)
    
    result = np.transpose(np.vstack((rdf.bins, rdf.rdf)))

    return result

for sys_name in ['c20_0p_s0_iter2', 'c20_100p_s512_iter2', 'c20_0p_s512']:

    #x = rdf_ow(TPR=f"../{sys_name}.md_PME.tpr", XTC=f"../water-mediated/{sys_name}.md_PME.dt-1000.xtc",
    #                ti=-5, tf=-1)

    #np.savetxt(f"{sys_name}.rdf.txt", x)
    """
    Use below to plot
    """
    x = np.loadtxt(f"{sys_name}.rdf.txt")

    plt.plot(x[:,0]/10., x[:,1], label=sys_name)
    
plt.xlabel("r (nm)")
plt.ylabel(r"$rdf_{ow-ow}$")
plt.legend()
plt.savefig("ow-ow.png")
#plt.show()
