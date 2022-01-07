import FR_mainCode as fr
import FR_mainCode_plus as frplus
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import pickle #to unpickle
import seaborn as sns

#%%
profile = np.arange(1000,61000,1000)
radius_cylinder = [5,7.5,10,12.5,15]
#radius_losange = [1,2,5,10,19]
slope_h2 = [-0.02,
          -0.01,
          0,
          0.01,
          0.02]
slope_heq = [-0.06,
          -0.03,
          0,
          0.03,
          0.06]
slope_losange = [-0.06,
          -0.03,
          0,
          0.03,
          0.06]
#%%

profile = np.arange(1000,61000,1000)


damping_cylinder = []
oscillation_cylinder = []
r_vector_cylinder = []
z_vector_cylinder = []
for r in radius_cylinder:
    col_damp = []
    col_oscillation = []
    for L in profile:
        resdic = frplus.calc_sim(r_fix=r,L=L,profile=True)
        col_damp.append(resdic['damping'])
        col_oscillation.append(resdic['oscillation'])
    damping_cylinder.append(col_damp) 
    oscillation_cylinder.append(col_oscillation)
    r_vector_cylinder.append(resdic['r_vector'])
    z_vector_cylinder.append(resdic['z_vector'])
#%%     
damping_H2 = []
oscillation_H2 = []
r_vector_H2 = []
z_vector_H2 = []
for m in slope_h2:
    col_damp = []
    col_oscillation = []
    for L in profile:
        resdic = frplus.calc_sim(r_fix=10,m=m,L=L,profile=True,z_fix='H2')
        col_damp.append(resdic['damping'])
        col_oscillation.append(resdic['oscillation'])
    damping_H2.append(col_damp) 
    oscillation_H2.append(col_oscillation)
    r_vector_H2.append(resdic['r_vector'])
    z_vector_H2.append(resdic['z_vector'])
    
damping_heq = []
oscillation_heq = []
r_vector_heq = []
z_vector_heq = []
for m in slope_heq:
    col_damp = []
    col_oscillation = []
    for L in profile:
        resdic = frplus.calc_sim(r_fix=10,m=m,L=L,profile=True,z_fix='heq')
        col_damp.append(resdic['damping'])
        col_oscillation.append(resdic['oscillation'])   
    damping_heq.append(col_damp) 
    oscillation_heq.append(col_oscillation)
    r_vector_heq.append(resdic['r_vector'])
    z_vector_heq.append(resdic['z_vector'])

#%%       
damping_losange = []
oscillation_losange = []
r_vector_losange = []
z_vector_losange = []
for m in slope_losange:
    col_damp = []
    col_oscillation = []
    for L in profile:
        resdic = frplus.calc_sim(shape = 'nodes',dh_below = 300,
                                 m=m)
        col_damp.append(resdic['damping'])
        col_oscillation.append(resdic['oscillation'])  
    damping_losange.append(col_damp) 
    oscillation_losange.append(col_oscillation)
    r_vector_losange.append(resdic['r_vector'])
    z_vector_losange.append(resdic['z_vector'])

#%%    
def plot_moulin(ax,r_vector,z_vector):
    ax.plot(r_vector[i],z_vector[i],color=colors[i],lw=2)
    ax.plot(-np.array(r_vector[i]),z_vector[i],color=colors[i],lw=2)
    ax.set_xlim([-60,30]) 
    ax.set_ylim([0,1500])


#colors:
Red = '#ED2224'
Orange = '#FBB263'
#Green = '#A2D39E'
Blue = '#15B4E9'
Purple = '#6B52A2'

colors = (Red, Orange ,(0,0,0), Blue, Purple) #5colors
linestyles = ('solid','solid',(0, (1, 0.5)),'solid','solid')

fig, ax = plt.subplots(figsize=(7.48,3),dpi=300)
gs1 = gridspec.GridSpec(3,4)

ax0a = plt.subplot(gs1[0])
ax0b = plt.subplot(gs1[1])
ax0c = plt.subplot(gs1[2])
ax0d = plt.subplot(gs1[3])
ax1 = plt.subplot(gs1[4])
ax2 = plt.subplot(gs1[8])
ax3 = plt.subplot(gs1[5])
ax4 = plt.subplot(gs1[9])
ax5 = plt.subplot(gs1[6])
ax6 = plt.subplot(gs1[10])
ax7 = plt.subplot(gs1[7])
ax8 = plt.subplot(gs1[11])


#frameon=False
offset=(5,5)
sns.despine(ax=ax0a, offset=offset, bottom=True, left=True)
sns.despine(ax=ax0b, offset=offset, bottom=True, left=True)
sns.despine(ax=ax0c, offset=offset, bottom=True, left=True)
sns.despine(ax=ax0d, offset=offset, bottom=True, left=True)
sns.despine(ax=ax1, offset=offset,  bottom=True)
sns.despine(ax=ax2, offset=offset)
sns.despine(ax=ax3, offset=offset, bottom=True, left=True)
sns.despine(ax=ax4, offset=offset, left=True)
sns.despine(ax=ax5, offset=offset, bottom=True, left=True)
sns.despine(ax=ax6, offset=offset, left=True)
sns.despine(ax=ax7, offset=offset, bottom=True, left=True)
sns.despine(ax=ax8, offset=offset, left=True)

ax0a.set_xticks([])
ax0b.set_xticks([])
ax0c.set_xticks([])
ax0d.set_xticks([])
ax1.set_xticks([])
ax3.set_xticks([])
ax5.set_xticks([])
ax7.set_xticks([])

ax1.set_xlim([0,60])
ax2.set_xlim([0,60])
ax3.set_xlim([0,60])
ax4.set_xlim([0,60])
ax5.set_xlim([0,60])
ax6.set_xlim([0,60])
ax7.set_xlim([0,60])
ax8.set_xlim([0,60])

ylim1=[0,8]
ylim2=[0,40]
ax1.set_ylim(ylim1)
ax2.set_ylim(ylim2)
ax3.set_ylim(ylim1)
ax4.set_ylim(ylim2)
ax5.set_ylim(ylim1)
ax6.set_ylim(ylim2)
ax7.set_ylim(ylim1)
ax8.set_ylim(ylim2)

ax1.set_yticks([0,2,4,6,8])
ax2.set_yticks([0,10,20,30,40]) 
ax0a.set_yticks([])
ax0b.set_yticks([])
ax0c.set_yticks([])
ax0d.set_yticks([])
ax3.set_yticks([])
ax4.set_yticks([])
ax5.set_yticks([])
ax6.set_yticks([])
ax7.set_yticks([])
ax8.set_yticks([])

ax1.set_ylabel('$\\tau_{osc}$ (days)')
ax2.set_ylabel('$\\tau_{damp}$ (days)')
ax4.set_xlabel('Distance from margin (km)')  

elements_ax0a = [Line2D([0], [0],  color=colors[0], lw=2, label=radius_cylinder[0]),
                 Line2D([0], [0],  color=colors[1], lw=2, label=radius_cylinder[1]),
                 Line2D([0], [0],  color=colors[2], lw=2, label=radius_cylinder[2]),
                 Line2D([0], [0],  color=colors[3], lw=2, label=radius_cylinder[3]),
                 Line2D([0], [0],   color=colors[4], lw=2, label=radius_cylinder[4])]

elements_ax0b = [Line2D([0], [0],  color=colors[0], lw=2, label=slope_h2[0]),
                 Line2D([0], [0],  color=colors[1], lw=2, label=slope_h2[1]),
                 Line2D([0], [0],  color=colors[2], lw=2, label=slope_h2[2]),
                 Line2D([0], [0],  color=colors[3], lw=2, label=slope_h2[3]),
                 Line2D([0], [0],   color=colors[4], lw=2, label=slope_h2[4])]

elements_ax0c = [Line2D([0], [0],  color=colors[0], lw=2, label=slope_heq[0]),
                 Line2D([0], [0],  color=colors[1], lw=2, label=slope_heq[1]),
                 Line2D([0], [0],  color=colors[2], lw=2, label=slope_heq[2]),
                 Line2D([0], [0],  color=colors[3], lw=2, label=slope_heq[3]),
                 Line2D([0], [0],   color=colors[4], lw=2, label=slope_heq[4])]

elements_ax0d = [Line2D([0], [0],  color=colors[0], lw=2, label=round((slope_losange[0]-5)/300,2)),
                 Line2D([0], [0],  color=colors[1], lw=2, label=round((slope_losange[1]-5)/300,2)),
                 Line2D([0], [0],  color=colors[2], lw=2, label=round((slope_losange[2]-5)/300,2)),
                 Line2D([0], [0],  color=colors[3], lw=2, label=round((slope_losange[3]-5)/300,2)),
                 Line2D([0], [0],   color=colors[4], lw=2, label=round((slope_losange[4]-5)/300,2))]

ax0a.legend(handles=elements_ax0a, loc=1,  labelspacing=0, handletextpad=1,  
       prop={'size': 6}, title='radius',title_fontsize=8, bbox_to_anchor=(0.45, 1))
ax0b.legend(handles=elements_ax0b, loc=1,  labelspacing=0, handletextpad=1,  
       prop={'size': 6}, title='slope',title_fontsize=8, bbox_to_anchor=(0.45, 1))
ax0c.legend(handles=elements_ax0c, loc=1,  labelspacing=0, handletextpad=1,  
       prop={'size': 6}, title='slope',title_fontsize=8, bbox_to_anchor=(0.45, 1))
ax0d.legend(handles=elements_ax0d, loc=1,  labelspacing=0, handletextpad=1,  
       prop={'size': 6}, title='slope',title_fontsize=8, bbox_to_anchor=(0.45, 1))
for i in np.arange(5): 

    #plot moulin shapes
    plot_moulin(ax0a,r_vector_cylinder,z_vector_cylinder)
    plot_moulin(ax0b,r_vector_H2,z_vector_H2)
    plot_moulin(ax0c,r_vector_heq,z_vector_heq)
    plot_moulin(ax0d,r_vector_losange,z_vector_losange)
    
    
    
    '''Change radius with fixed slope'''  
    ax1.plot(profile/1000,oscillation_cylinder[i],color=colors[i], lw=2, linestyle=linestyles[i])
    ax2.plot(profile/1000,damping_cylinder[i],color=colors[i], lw=2, linestyle=linestyles[i])
    
    '''Change slope with h_middle fixed'''
    ax3.plot(profile/1000,oscillation_H2[i],color=colors[i], lw=2, linestyle=linestyles[i])
    ax4.plot(profile/1000,damping_H2[i],color=colors[i], lw=2, linestyle=linestyles[i])
    
    '''Change slope with h_eq fixed'''
    ax5.plot(profile/1000,oscillation_heq[i],color=colors[i], lw=2, linestyle=linestyles[i])
    ax6.plot(profile/1000,damping_heq[i],color=colors[i], lw=2, linestyle=linestyles[i])

    '''Change radius with h_eq fixed for losange'''
    ax7.plot(profile/1000,oscillation_losange[i],color=colors[i], lw=2, linestyle=linestyles[i])
    ax8.plot(profile/1000,damping_losange[i],color=colors[i], lw=2, linestyle=linestyles[i])


#plt.subplots_adjust(left=0.08, right=0.99, top=0.96, bottom=0.14)


#plt.savefig('Figures/Timescales_osc_damp.pdf')






