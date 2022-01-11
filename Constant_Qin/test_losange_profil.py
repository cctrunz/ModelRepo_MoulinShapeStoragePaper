import FR_mainCode as fr
import FR_mainCode_plus as frplus
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle #to unpickle

profile = np.arange(1000,61000,1000)
H0 = 1000/np.sqrt(30000) #assures ice thickness of 1500 m at 60 km from edge
H = H0*np.sqrt(profile[2]) # ice thickness for a square root glacier.
#%%
reslist_cylinder = [] # Create empty list
for r in [5,7.5,10,12.5,15]:
    resdic = frplus.calc_sim(shape = 'linear', r_fix=r, H_fix=H)
    reslist_cylinder.append(resdic) # add dictionnary to the bottom of the list
frplus.plot_FR_JGR(reslist_cylinder,position='top')    


reslist_H2 = [] # Create empty list
for m in [-0.02,
          -0.01,
          0,
          0.01,
          0.02]:
    resdic = frplus.calc_sim(shape = 'linear', r_fix=10,m=m,z_fix='H2', H_fix=H)
    reslist_H2.append(resdic) # add dictionnary to the bottom of the list
frplus.plot_FR_JGR(reslist_H2,letter=['d','e','f'],position='middle')  


reslist_heq = [] # Create empty list
for m in [-0.06,
          -0.03,
          0,
          0.03,
          0.06]:
    resdic = frplus.calc_sim(shape = 'linear', r_fix=10, m=m,z_fix='heq', H_fix=H)
    reslist_heq.append(resdic) # add dictionnary to the bottom of the list
frplus.plot_FR_JGR(reslist_heq,letter=['g','h','i'],position='middle')


#%%

text1 = ['Cylinder r=5m','Cylinder r=7.5m','Cylinder r=10m','Cylinder r=12.5m','Cylinder r=15m']
text2 = ['Cone H/2 m=-0.02','Cone H/2 m=-0.01','Cone H/2 m=0','Cone H/2 m=0.01','Cone H/2 m=0.02']
text3 = ['Cone heq m=-0.06','Cone heq m=-0.03','Cone heq m=0','Cone heq m=0.03','Cone heq m=0.06']

plt.figure(figsize=(7,8),dpi=300)
res1 = reslist_cylinder
res2 = reslist_H2
res3 = reslist_heq
index = [1,4,7,10]
for i in [0,1,2,3]:
    ii = index[i]
  
    ax1 = plt.subplot(5,3,ii)
    plt.plot(res1[i]['td'],res1[i]['hd'],label='h')
    plt.plot(res1[i]['td'],res1[i]['hd_fit'],label='h_fit',linestyle='--')
    plt.ylabel('h')
    #plt.xlim( -1,30)
    #plt.ylim(670,820)
    plt.xticks([])
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    #ax1.spines['left'].set_bounds(700, 800)
    ax1.tick_params(direction='in')
    
    ax2 = plt.subplot(5,3,ii+1)
    plt.plot(res2[i]['td'],res2[i]['hd'],label='h')
    plt.plot(res2[i]['td'],res2[i]['hd_fit'],label='h_fit',linestyle='--')
    #plt.xlim(-1,30)
    #plt.ylim(670,820)
    plt.xticks([])
    plt.yticks([])
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    #ax2.spines['left'].set_bounds(700, 800)
    
    
    ax3 = plt.subplot(5,3,ii+2)
    plt.plot(res3[i]['td'],res3[i]['hd'],label='h')
    plt.plot(res3[i]['td'],res3[i]['hd_fit'],label='h_fit',linestyle='--')  
    #plt.xlim(-1,30)
    #plt.ylim(670,820)
    plt.xticks([])
    plt.yticks([])
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    #ax3.spines['left'].set_bounds(700, 800)
    ax3.tick_params(direction='in')
    
    plt.subplots_adjust(hspace=0.1)
    plt.subplots_adjust(wspace=0.1)
    
    
    
    ax1.text(0.9, 0.9, text1[i], ha='right', va='top', transform=ax1.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    ax2.text(0.9, 0.9, text2[i], ha='right', va='top', transform=ax2.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    ax3.text(0.9, 0.9, text3[i], ha='right', va='top', transform=ax3.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))

for i in [4]:
    ii = 13
  
    ax1 = plt.subplot(5,3,ii)
    plt.plot(res1[i]['td'],res1[i]['hd'],label='h')
    plt.plot(res1[i]['td'],res1[i]['hd_fit'],label='h_fit',linestyle='--')
    plt.ylabel('h')
    plt.xlabel('t')
    #plt.xlim( -1,30)
    #plt.ylim(670,820)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    #ax1.spines['left'].set_bounds(700, 800)
    ax1.spines['bottom'].set_bounds(0, 30)
    ax1.tick_params(direction='in')
    
    ax2 = plt.subplot(5,3,ii+1)
    plt.plot(res2[i]['td'],res2[i]['hd'],label='h')
    plt.plot(res2[i]['td'],res2[i]['hd_fit'],label='h_fit',linestyle='--')
    plt.xlabel('t')
    #plt.xlim(-1,30)
    #plt.ylim(670,820)
    plt.yticks([])
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.spines['left'].set_bounds(700, 800)
    ax2.tick_params(direction='in')
    ax2.spines['bottom'].set_bounds(0, 30)
    
    ax3 = plt.subplot(5,3,ii+2)
    plt.plot(res3[i]['td'],res3[i]['hd'],label='h')
    plt.plot(res3[i]['td'],res3[i]['hd_fit'],label='h_fit',linestyle='--') 
    plt.xlabel('t')
    #plt.xlim(-1,30)
    #plt.ylim(670,820)
    plt.yticks([])
    plt.setp(ax3.get_yticklabels(), visible=False)

    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    #ax3.spines['left'].set_bounds(700, 800)
    ax3.tick_params(direction='in')
    ax3.spines['bottom'].set_bounds(0, 30)
    
    ax1.text(0.9, 0.9, text1[i], ha='right', va='top', transform=ax1.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    ax2.text(0.9, 0.9, text2[i], ha='right', va='top', transform=ax2.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    ax3.text(0.9, 0.9, text3[i], ha='right', va='top', transform=ax3.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    
    
            
plt.savefig('Figures/FR_compare_Fit.pdf')        
    