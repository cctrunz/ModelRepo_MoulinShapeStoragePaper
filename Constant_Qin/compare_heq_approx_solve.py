# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 11:59:56 2022

@author: truc3007
"""

# COMPARISON        heq approx and heq ...


import FR_mainCode_plus as frplus
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle #to unpickle

H = np.arange(100,1000,100)

reslist_cylinder = [] 
for H_fix in H:
    resdic = frplus.calc_sim(shape = 'linear', H_fix=H_fix)
    reslist_cylinder.append(resdic) # add dictionnary to the bottom of the list

rhoi = 910 #ice density $, kg/m3
rhow = 1000 #water density, kg/m3
g =  9.8 #Gravity
#H = 1000


#%%

heq_approx = []
heq_solve = []

for i in np.arange(len(reslist_cylinder)):
    Pi = rhoi*g*H[i]
    hfl = Pi/rhow/g
    heq_approx.append(reslist_cylinder[i]['h_eq_approx']*hfl)
    heq_solve.append(reslist_cylinder[i]['h_eq_nd']*hfl)
    
#%%


plt.figure()
plt.plot([0,600],[0,600], color='black', linestyle='--')
plt.plot(heq_solve,heq_approx, 'o')  

plt.title('simulations for ice thickness from 100 -- 1000 m')
plt.xlabel('$h_{eq}$ solver') 
plt.ylabel('$h_{eq}$ approx') 
