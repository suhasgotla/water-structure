# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 10:40:17 2019

@author: abhik
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
from numpy import (array, dot, arccos, clip)
from numpy.linalg import norm
from MDAnalysis import *
import MDAnalysis
import MDAnalysis.lib.distances
from MDAnalysis.lib.distances import distance_array
import numpy.linalg
import scipy.stats
import matplotlib.pyplot as plt
import math
from MDAnalysis.core.topologyobjects import Angle
from itertools import combinations
import pandas as pd
from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def smooth(y, box_pts):
	box = np.ones(box_pts)/box_pts
	y_smooth = np.convolve(y, box, mode='same')
	return y_smooth

def makehist(data, datalabel, colormark):
    y,binEdges = np.histogram(data,bins=100)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    bincenters = bincenters/10.0
    plt.plot(bincenters, y/np.trapz(y, bincenters), colormark, label=datalabel, lw =3)
    
    
def watervec_anglewnormal(u, molnum):
    Z = np.array([0, 0, 1])
    oh2 = u.select_atoms("molnum " + str(molnum) + " and name OH2").positions
    h1 = u.select_atoms("molnum " + str(molnum) + " and name H1").positions
    h2 = u.select_atoms("molnum " + str(molnum) + " and name H2").positions
    vec = np.add((oh2-h2), (oh2-h1))[0]
    return dot(vec,Z)/norm(vec)/norm(Z)
      
#u = Universe("../prod_cr.tpr", "../prod_cr.xtc")
#prod_noca.xtc
#u = Universe("../prod.tpr", "../../Backmap/traj.xtc")
u = Universe("../prod_noca.tpr", "../prod_noca.xtc")
def costhetawater(u):
	S_large = np.array([])
	Q = np.array([])
	Z = np.array([])
	DENS = []
	All_int = []
	atomgroups= []
	cross_interactors = np.zeros((7, 7))
	side = 'up'
	DPC = []
	DPS = []
	for ts in u.trajectory[-1000:-1:1]:
		print ts
		cut = 100.0
		cutoff = 3.5
		wat = u.select_atoms("resname TIP3 and name OH2")
		water = u.select_atoms("resname TIP3")
		water_fixed = water.positions.reshape((wat.n_atoms, 3, 3))
		water_vecs = np.array([(i[0]-i[1]) + (i[0]-i[2]) for i in water_fixed])
		
		pc_choline = u.select_atoms("resname POPC and name N")
		pc_phosphate = u.select_atoms("resname POPC and name P")
		pc_o21 = u.select_atoms("resname POPC and name O21")
		pc_o22 = u.select_atoms("resname POPC and name O22")
		pc_o31 = u.select_atoms("resname POPC and name O31")
		pc_o32 = u.select_atoms("resname POPC and name O32")
		pc_o = pc_o21 + pc_o22 + pc_o31 + pc_o32 
		
		ps_serine_1 = u.select_atoms("resname POPS and name N")
		ps_serine_2 = u.select_atoms("resname POPS and name O13A")
		ps_serine_3 = u.select_atoms("resname POPS and name O13B")
		ps_ser = ps_serine_1 + ps_serine_2 + ps_serine_3
		ps_phosphate = u.select_atoms("resname POPS and name P")
		ps_o21 = u.select_atoms("resname POPS and name O21")
		ps_o22 = u.select_atoms("resname POPS and name O22")
		ps_o31 = u.select_atoms("resname POPS and name O31")
		ps_o32 = u.select_atoms("resname POPS and name O32")
		ps_o = ps_o21 + ps_o22 + ps_o31 + ps_o32
		
		po4 = u.select_atoms("name P")
		ca = u.select_atoms("name CAL")
		
		wat_pos_z = wat.positions[:,2]
		bilcenter_z = u.select_atoms("name P").center_of_mass()[2]
		dist_pc_N = np.abs(pc_choline.positions[:,2] - bilcenter_z)
		dist_pc_P = np.abs(pc_phosphate.positions[:,2] - bilcenter_z) 
		dist_pc_O = np.abs(pc_o.positions[:,2] - bilcenter_z) 
		
		dist_ps_N = np.abs(ps_ser.positions[:,2] - bilcenter_z) 
		dist_ps_P = np.abs(ps_phosphate.positions[:,2] - bilcenter_z) 
		dist_ps_O = np.abs(ps_o.positions[:,2] - bilcenter_z)  
		
		DPC.append([dist_pc_N , dist_pc_P , dist_pc_O]) 
		DPS.append([dist_ps_N , dist_ps_P , dist_ps_O]) 
		
		if side == 'up':
		    Z = np.array([0, 0, 1.0])
		    water_angle = np.array([dot(vec,Z)/norm(vec)/norm(Z) for vec in water_vecs])
		    wat_up_oh2 = wat[(wat.positions[:,2] > bilcenter_z)]
		    water_up = water_fixed[(wat.positions[:,2] > bilcenter_z)]
		    water_vecs_up = np.array([(i[0]-i[1]) + (i[0]-i[2]) for i in water_up])
		    water_angle_up = np.array([dot(vec,Z)/norm(vec)/norm(Z) for vec in water_vecs_up])
		    water_up_molnum = wat_up_oh2.molnums
		    dict_up = dict(zip(water_up_molnum, water_angle_up))
		    
		    #PC-UP
		    D_pc_chol_up = distance_array(wat_up_oh2.positions, pc_choline[(pc_choline.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_phos_up = distance_array(wat_up_oh2.positions, pc_phosphate[(pc_phosphate.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o21_up = distance_array(wat_up_oh2.positions, pc_o21[(pc_o21.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o22_up = distance_array(wat_up_oh2.positions, pc_o22[(pc_o22.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o31_up = distance_array(wat_up_oh2.positions, pc_o31[(pc_o31.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o32_up = distance_array(wat_up_oh2.positions, pc_o32[(pc_o32.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		
		    #PS-UP
		    D_ps_ser1_up = distance_array(wat_up_oh2.positions, ps_serine_1[(ps_serine_1.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_ser2_up = distance_array(wat_up_oh2.positions, ps_serine_2[(ps_serine_2.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_ser3_up = distance_array(wat_up_oh2.positions, ps_serine_3[(ps_serine_3.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_phos_up = distance_array(wat_up_oh2.positions, ps_phosphate[(ps_phosphate.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o21_up = distance_array(wat_up_oh2.positions, ps_o21[(ps_o21.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o22_up = distance_array(wat_up_oh2.positions, ps_o22[(ps_o22.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o31_up = distance_array(wat_up_oh2.positions, ps_o31[(ps_o31.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o32_up = distance_array(wat_up_oh2.positions, ps_o32[(ps_o32.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		    
		    d_caup = distance_array(wat_up_oh2.positions, ca[(ca.positions[:,2]>bilcenter_z)].positions, np.copy(ts.dimensions))
		else:
		    Z = np.array([0, 0, 1.0])
		    water_angle = np.array([dot(vec,Z)/norm(vec)/norm(Z) for vec in water_vecs])
		    wat_up_oh2 = wat[(wat.positions[:,2] < bilcenter_z)]
		    water_up = water_fixed[(wat.positions[:,2] < bilcenter_z)]
		    water_vecs_up = np.array([(i[0]-i[1]) + (i[0]-i[2]) for i in water_up])
		    water_angle_up = -1.0*np.array([dot(vec,Z)/norm(vec)/norm(Z) for vec in water_vecs_up])
		    water_up_molnum = wat_up_oh2.molnums
		    dict_up = dict(zip(water_up_molnum, water_angle_up))
		
		    #PC-UP
		    D_pc_chol_up = distance_array(wat_up_oh2.positions, pc_choline[(pc_choline.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_phos_up = distance_array(wat_up_oh2.positions, pc_phosphate[(pc_phosphate.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o21_up = distance_array(wat_up_oh2.positions, pc_o21[(pc_o21.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o22_up = distance_array(wat_up_oh2.positions, pc_o22[(pc_o22.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o31_up = distance_array(wat_up_oh2.positions, pc_o31[(pc_o31.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_pc_o32_up = distance_array(wat_up_oh2.positions, pc_o32[(pc_o32.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		
		    #PS-UP
		    D_ps_ser1_up = distance_array(wat_up_oh2.positions, ps_serine_1[(ps_serine_1.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_ser2_up = distance_array(wat_up_oh2.positions, ps_serine_2[(ps_serine_2.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_ser3_up = distance_array(wat_up_oh2.positions, ps_serine_3[(ps_serine_3.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_phos_up = distance_array(wat_up_oh2.positions, ps_phosphate[(ps_phosphate.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o21_up = distance_array(wat_up_oh2.positions, ps_o21[(ps_o21.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o22_up = distance_array(wat_up_oh2.positions, ps_o22[(ps_o22.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o31_up = distance_array(wat_up_oh2.positions, ps_o31[(ps_o31.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    D_ps_o32_up = distance_array(wat_up_oh2.positions, ps_o32[(ps_o32.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		    
		    d_caup = distance_array(wat_up_oh2.positions, ca[(ca.positions[:,2]<bilcenter_z)].positions, np.copy(ts.dimensions))
		
		pc_up_N = wat_up_oh2[np.where(D_pc_chol_up<5.85)[0]]
		pc_up_P = wat_up_oh2[np.where(D_pc_phos_up<4.55)[0]]
		pc_up_O21 = wat_up_oh2[np.where(D_pc_o21_up<3.65)[0]]
		pc_up_O22 = wat_up_oh2[np.where(D_pc_o22_up<3.25)[0]]
		pc_up_O31 = wat_up_oh2[np.where(D_pc_o31_up<3.65)[0]]
		pc_up_O32 = wat_up_oh2[np.where(D_pc_o32_up<3.25)[0]]
		pc_up_O = (pc_up_O21 + pc_up_O22 + pc_up_O31 + pc_up_O32).unique
		
		ps_up_N1 = wat_up_oh2[np.where(D_ps_ser1_up<3.45)[0]]
		ps_up_N2 = wat_up_oh2[np.where(D_ps_ser2_up<3.45)[0]]
		ps_up_N3 = wat_up_oh2[np.where(D_ps_ser3_up<3.45)[0]]
		ps_up_N = ps_up_N1.union(ps_up_N2).union(ps_up_N3)
		ps_up_P = wat_up_oh2[np.where(D_ps_phos_up<4.55)[0]]
		ps_up_O21 = wat_up_oh2[np.where(D_ps_o21_up<3.65)[0]]
		ps_up_O22 = wat_up_oh2[np.where(D_ps_o22_up<3.25)[0]]
		ps_up_O31 = wat_up_oh2[np.where(D_ps_o31_up<3.65)[0]]
		ps_up_O32 = wat_up_oh2[np.where(D_ps_o32_up<3.25)[0]]
		ps_up_O = (ps_up_O21 + ps_up_O22 + ps_up_O31 + ps_up_O32).unique
		
		ca_up = wat_up_oh2[np.where(D_pc_chol_up<3.00)[0]]
		
		all_interactors = pc_up_N + pc_up_P + pc_up_O + ps_up_N + ps_up_P + ps_up_O
		other_interactor = wat_up_oh2 - all_interactors.unique
		pc_up_NP = pc_up_N.intersection(pc_up_P)
		pc_up_NO = pc_up_N.intersection(pc_up_O)
		pc_up_PO = pc_up_P.intersection(pc_up_O)
		pc_up_NPO = pc_up_NP.intersection(pc_up_O)
		
		ps_up_NP = ps_up_N.intersection(ps_up_P)
		ps_up_NO = ps_up_N.intersection(ps_up_O)
		ps_up_PO = ps_up_P.intersection(ps_up_O)
		ps_up_NPO = ps_up_NP.intersection(ps_up_O)
		
		Superdict = {}
		
		AllPC = [pc_up_N, pc_up_P, pc_up_O, pc_up_NP, pc_up_NO, pc_up_PO, pc_up_NPO]
		pc_names = ['pcN', 'pcP', 'pcO', 'pcNP', 'pcNO', 'pcPO', 'pcNPO']
		
		a=0
		for partition in AllPC:
		    anglelist = np.array([dict_up[i] for i in partition.molnums])
		    z = np.abs(partition.positions[:,2]-bilcenter_z)
		    fin = np.array(zip(anglelist, z))
		    Superdict.update( {pc_names[a] : fin} )
		    a+=1
		
		AllPS = [ps_up_N, ps_up_P, ps_up_O, ps_up_NP, ps_up_NO, ps_up_PO, ps_up_NPO]
		ps_names = ['psN', 'psP', 'psO', 'psNP', 'psNO', 'psPO', 'psNPO']
		
		a=0
		for partition in AllPS:
		    anglelist = np.array([dict_up[i] for i in partition.molnums])
		    z = np.abs(partition.positions[:,2]-bilcenter_z)
		    fin = np.array(zip(anglelist, z))
		    Superdict.update( {ps_names[a] : fin} )
		    a+=1
		        
		PCPS = []
		a=0
		for partitioni in AllPC:
		    b=0
		    for partitionj in AllPS:
		        newgroup = partitioni.intersection(partitionj)
		        PCPS.append(newgroup)
		        newname = pc_names[a]+ps_names[b]
		        anglelist = np.array([dict_up[i] for i in newgroup.molnums])
		        z = np.abs(newgroup.positions[:,2] - bilcenter_z)
		        fin = np.array(zip(anglelist, z))
		        Superdict.update( {newname : fin} )
		        cross_interactors[a,b]+=newgroup.n_atoms
		        b+=1
		    a+=1
		
		othergroupname = 'Other'
		otheranglelist = np.array([dict_up[i] for i in other_interactor.molnums])
		otherz = np.abs(other_interactor.positions[:,2] - bilcenter_z)
		fin = np.array(zip(otheranglelist, otherz))
		Superdict.update( {othergroupname : fin} )
		All_int.append(Superdict)
		atomgroups.append([AllPC, AllPS])
	return np.array(All_int)
    

u = Universe("../prod_noca.tpr", "../prod_noca.xtc")
all_noca = costhetawater(u) 

u = Universe("../prod_cr.tpr", "../prod_cr.xtc")
all_ca = costhetawater(u) 
   
#atomgroups = np.array(atomgroups)    


'''
makehist(np.hstack(DPC[:,0]), "pcN", "b")
makehist(np.hstack(DPC[:,1]), "pcP", "g")
makehist(np.hstack(DPC[:,2]), "pcO", "r")

makehist(np.hstack(DPS[:,0]), "pSN", "b--")
makehist(np.hstack(DPS[:,1]), "pSP", "g--")
makehist(np.hstack(DPS[:,2]), "psO", "r--")

plt.legend() 
'''


'''

plt.clf()
counter = 0
for labeli in range(len(pc_names)):
    for labelj in range(len(ps_names)):
        try:
            counter+=1
            a = np.hstack(np.array([i[pc_names[labeli] + ps_names[labelj]][:,0] for i in All_int]))   
            b = np.hstack(np.array([i[pc_names[labeli] + ps_names[labelj]][:,1] for i in All_int])) 
            n, bins = numpy.histogram(b, bins=100, range=(0,100), weights=a)
            n1, bins1 = numpy.histogram(b, bins=100, range=(0,100))
            Xaxis = np.linspace(min(b), max(b), 100)  
            plt.plot(Xaxis, np.nan_to_num(n), label=pc_names[labeli]+ps_names[labelj], color=2) 
        except IndexError:
            pass
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n<cos($\Theta$)>')
plt.savefig("./Repeat/ncostheta_cross" + pc_names[labeli]+ps_names[labelj] + ".png", dpi=300)
            
'''    


plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorpc[counter], label=pc_names[labeli]) 
    
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorps[counter], label=ps_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n<cos($\Theta$)>')
plt.savefig("./Repeat/ca_ncostheta.svg", dpi=300)

plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorpc[counter], label=pc_names[labeli]) 
    
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorps[counter], label=ps_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n<cos($\Theta$)>')
plt.savefig("./Repeat/noca_ncostheta.svg", dpi=300)


plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n1), colorpc[counter], label=pc_names[labeli]) 
    
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n1), colorps[counter], label=ps_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n')
plt.savefig("./Repeat/ca_n.svg", dpi=300)

plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n1), colorpc[counter], label=pc_names[labeli]) 
    
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n1), colorps[counter], label=ps_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n')
plt.savefig("./Repeat/noca_n.svg", dpi=300)








#PC
plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorpc[counter], label=pc_names[labeli]) 
    
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorps[counter], label=ps_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n<cos($\Theta$)>')
plt.savefig("./Repeat/noca_n.svg", dpi=300)


#PS
plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n1), colorpc[counter], label="noca_"+pc_names[labeli]) 
    
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorps[counter], label="ca_"+pc_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n<cos($\Theta$)>')
plt.savefig("./Repeat/pc_canoca.svg", dpi=300)

plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n1), colorpc[counter], label="noca_"+ps_names[labeli]) 
    
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a)
    n1, bins1 = numpy.histogram(b, bins=100)
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorps[counter], label="ca_"+ps_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n<cos($\Theta$)>')
plt.savefig("./Repeat/ps_canoca.svg", dpi=300)





plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorpc[counter], label="noca_"+pc_names[labeli]) 
    
    a = np.hstack(np.array([i[pc_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[pc_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorps[counter], label="ca_"+pc_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n')
plt.savefig("./Repeat/pc_canoca.svg", dpi=300)

plt.clf()
counter = 0
colorpc = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
colorps = ['b--', 'g--', 'r--', 'c--', 'm--', 'y--', 'k--']
for labeli in range(len(pc_names)):
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_noca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_noca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorpc[counter], label="noca_"+ps_names[labeli]) 
    
    a = np.hstack(np.array([i[ps_names[labeli]][:,0] for i in all_ca]))   
    b = np.hstack(np.array([i[ps_names[labeli]][:,1] for i in all_ca])) 
    n, bins = numpy.histogram(b, bins=100, weights=a, range=(1,50))
    n1, bins1 = numpy.histogram(b, bins=100, range=(1,50))
    Xaxis = np.linspace(min(b), max(b), 100)  
    plt.plot(Xaxis, np.nan_to_num(n), colorps[counter], label="ca_"+ps_names[labeli]) 
    counter+=1
plt.legend()
plt.xlabel('Z (in angstroms)')
plt.ylabel('n')
plt.savefig("./Repeat/ps_canoca.svg", dpi=300)
