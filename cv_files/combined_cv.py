import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 
import numpy,sys,os
numpy.set_printoptions(threshold=sys.maxsize)


cv_all_list = pd.DataFrame(columns=["frame","pca1","rep","system"])

for i in range(1,19):
	cv_file = pd.read_csv('cv_'+str(i)+'_PBC',comment='#',names=['frame','pca1'],delim_whitespace=True)
	cv_file['rep']=i
	cv_file['system']="WT"
	cv_all_list = pd.concat([cv_all_list,cv_file],ignore_index='True')

cv_all_list.to_csv('cv_WT_mwMETAD_Milosz_all', sep='\t')



