import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import os 
import sys
import numpy
import gc

numpy.set_printoptions(threshold=sys.maxsize)

def grouping_residue(res_id):
    group=0
    if  403<=res_id<=415 or 500<=res_id<=507:
        group=1
    elif 416<=res_id<=430 or 446<=res_id<=462 or 491<=res_id<=498:
        group=2
    elif 470<=res_id<=490:
        group=3
    return group
 
def tracking_distance_mutation(repname,sysname,neighbors_list,cv_list,mutation_name,frame):
    mutation_list = neighbors_list.residues
    res_count = len(mutation_list)
    res_name = mutation_list.residues.resnames
    res_id = mutation_list.residues.resids
    res_chain = mutation_list.residues.segids
    if res_count ==0:
        print(repname,sysname,"nan",frame,*cv_list['pca1'],"nan","nan","nan","nan","nan")
        output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(repname,sysname,"nan",frame,*cv_list['pca1'],"nan","nan","nan","nan","nan"))
        output.flush()
        gc.collect()

    elif res_count==1:
        group = grouping_residue(res_id)
        if group==0:
            print(repname,sysname,mutation_name,frame,*cv_list['pca1'],res_count,*res_chain,*res_name,*res_id,group)
        else:
            output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(repname,sysname,mutation_name,frame,*cv_list['pca1'],res_count,*res_chain,*res_name,*res_id,group))
            output.flush()
            gc.collect()
    elif res_count >1:
        for i in range(res_count):
            group = grouping_residue(res_id[i])
            if group==0:
                print(repname,sysname,mutation_name,frame,*cv_list['pca1'],res_count,res_chain[i],res_name[i],res_id[i],group)
            else:
                output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%(repname,sysname,mutation_name,frame,*cv_list['pca1'],res_count,res_chain[i],res_name[i],res_id[i],group))
                output.flush()
                gc.collect()

### I/O ###
path=''
cv_all = pd.read_csv('cv_files/cv_WT_mwMETAD_Milosz_all',delim_whitespace=True)

# print(cv_all.head())
output = open(path+'WT_mwMETAD_contact_count.dat', 'w')
output.write("rep\tsystem\tmutant\tframe\tpca1\tres_count\tres_chain\tres_name\tres_id\tgroup\n")

for rep in range(1,19):
    u = mda.Universe(path+'../rechained_wt.pdb',path+'../nojump-'+str(rep)+'-wt.xtc')

    neighbors_S371 = u.select_atoms('(around 5 segid A and resid 371) and not segid A and protein',updating=True)
    neighbors_S373 = u.select_atoms('(around 5 segid A and resid 373) and not segid A and protein',updating=True)
    neighbors_S375 = u.select_atoms('(around 5 segid A and resid 375) and not segid A and protein',updating=True)
    # print(neighbors_S371.resids)
    # print(rep)
    for ts in u.trajectory[1::1]:
        cv_list = cv_all.query('system=="WT" and rep==%s and frame==%s'%(rep,ts.frame))
        tracking_distance_mutation(rep,"WT",neighbors_S371,cv_list,"S371",ts.frame)
        tracking_distance_mutation(rep,"WT",neighbors_S373,cv_list,"S373",ts.frame)
        tracking_distance_mutation(rep,"WT",neighbors_S375,cv_list,"S375",ts.frame)
    del u, neighbors_S371, neighbors_S373, neighbors_S375
output.close()
