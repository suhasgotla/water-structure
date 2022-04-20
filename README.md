# water-structure
Orientational order parameter analysis for the solvation shells of chitosan-SDS hydrogel networks/


RDF functions in water-rdf.py and water-neighbor_search.py are used to determine cutoffs for solavation shell detrmination.

water-neighbour_search.py parses trajectories for solvation waters, and generates lists of nearest neighbors for each.
> Each timeframe generates one .csv file of the 4 closest water neighbors (called O1s), and the next 4 nearest neighbors (called O2s)
> Files are named <sys_nmae>.o1o2<time>.csv. By default, <time> is  ordered as 0, -1, -2, -3 ... and so on, ordered in reverse, where 0 is the last frame of the simulation.
  
water-order.py uses water neighbor lists, along with **corresponding** trajectory timeframes todetermine the orientational order parameter.
> Output is a single csv file containing the ensemble of values of q for each solvation shell water molecule, across all timeframes.
> Files are named <sys_name>.q.csv --> the final result
  
See Errington & Debenedetti, Nature (2001) for the definition of orientational order parameter.
