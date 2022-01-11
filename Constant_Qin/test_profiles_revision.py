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
alpha_fit_cylinder = []
beta_fit_cylinder = []
pcov_cylinder = []
for r in radius_cylinder:
    col_damp = []
    col_oscillation = []
    col_alpha_fit = []
    col_beta_fit = []
    col_pcov = []
    for L in profile:
        resdic = frplus.calc_sim(r_fix=r,L=L,profile=True)
        col_damp.append(resdic['damping_fit'])
        col_oscillation.append(resdic['oscillation_fit'])
        col_alpha_fit.append(resdic['alpha_fit'])
        col_beta_fit.append(resdic['beta_fit'])
        col_pcov.append(resdic['pcov'])
    damping_cylinder.append(col_damp) 
    oscillation_cylinder.append(col_oscillation)
    r_vector_cylinder.append(resdic['r_vector'])
    z_vector_cylinder.append(resdic['z_vector'])
    alpha_fit_cylinder.append(col_alpha_fit)
    beta_fit_cylinder.append(col_beta_fit)
    pcov_cylinder.append(col_pcov)
#%%     
damping_H2 = []
oscillation_H2 = []
r_vector_H2 = []
z_vector_H2 = []
alpha_fit_H2 = []
beta_fit_H2 = []
pcov_H2 = []
for m in slope_h2:
    col_damp = []
    col_oscillation = []
    col_alpha_fit = []
    col_beta_fit = []
    col_pcov = []
    for L in profile:
        resdic = frplus.calc_sim(r_fix=10,m=m,L=L,profile=True,z_fix='H2')
        col_damp.append(resdic['damping_fit'])
        col_oscillation.append(resdic['oscillation_fit'])
        col_alpha_fit.append(resdic['alpha_fit'])
        col_beta_fit.append(resdic['beta_fit'])
        col_pcov.append(resdic['pcov'])
    damping_H2.append(col_damp) 
    oscillation_H2.append(col_oscillation)
    r_vector_H2.append(resdic['r_vector'])
    z_vector_H2.append(resdic['z_vector'])
    alpha_fit_H2.append(col_alpha_fit)
    beta_fit_H2.append(col_beta_fit)
    pcov_H2.append(col_pcov)

    
damping_heq = []
oscillation_heq = []
r_vector_heq = []
z_vector_heq = []
alpha_fit_heq = []
beta_fit_heq = []
pcov_heq = []
for m in slope_heq:
    col_damp = []
    col_oscillation = []
    col_alpha_fit = []
    col_beta_fit = []
    col_pcov = []
    for L in profile:
        resdic = frplus.calc_sim(r_fix=10,m=m,L=L,profile=True,z_fix='heq')
        col_damp.append(resdic['damping_fit'])
        col_oscillation.append(resdic['oscillation_fit'])   
        col_alpha_fit.append(resdic['alpha_fit'])
        col_beta_fit.append(resdic['beta_fit'])
        col_pcov.append(resdic['pcov'])
    damping_heq.append(col_damp) 
    oscillation_heq.append(col_oscillation)
    r_vector_heq.append(resdic['r_vector'])
    z_vector_heq.append(resdic['z_vector'])
    alpha_fit_heq.append(col_alpha_fit)
    beta_fit_heq.append(col_beta_fit)
    pcov_heq.append(col_pcov)
    
#%%       
damping_losange = []
oscillation_losange = []
r_vector_losange = []
z_vector_losange = []
alpha_fit_losange = []
beta_fit_losange = []
pcov_losange = []
for m in slope_losange:
    col_damp = []
    col_oscillation = []
    col_alpha_fit = []
    col_beta_fit = []
    col_pcov = []
    for L in profile:
        resdic = frplus.calc_sim(shape = 'nodes',#,dh_below = 300,
                                 r_heq = 10,
                                 profile = True,
                                 L=L,
                                 slope = True,
                                 m=m)
        col_damp.append(resdic['damping_fit'])
        col_oscillation.append(resdic['oscillation_fit']) 
        col_alpha_fit.append(resdic['alpha_fit'])
        col_beta_fit.append(resdic['beta_fit'])
        col_pcov.append(resdic['pcov'])
    damping_losange.append(col_damp) 
    oscillation_losange.append(col_oscillation)
    r_vector_losange.append(resdic['r_vector'])
    z_vector_losange.append(resdic['z_vector'])
    alpha_fit_losange.append(col_alpha_fit)
    beta_fit_losange.append(col_beta_fit)
    pcov_losange.append(col_pcov)
    
#%%    
def plot_moulin(ax,r_vector,z_vector):
    ax.plot(r_vector[i],z_vector[i],color=colors[i],lw=lw)
    ax.plot(-np.array(r_vector[i]),z_vector[i],color=colors[i],lw=lw)
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

fig, ax = plt.subplots(figsize=(7.48,4),dpi=300)
gs1 = gridspec.GridSpec(3,4,height_ratios=[2, 3,3])

lw = 1.5

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
sns.despine(ax=ax0a, offset=offset, bottom=True)
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

ax0a.set_xticks([-20,0,20])
ax0a.set_xticklabels(['-r', '0', 'r'])
ax0b.set_xticks([-20,0,20])
ax0b.set_xticklabels(['-r', '0', 'r'])
ax0c.set_xticks([-20,0,20])
ax0c.set_xticklabels(['-r', '0', 'r'])
ax0d.set_xticks([-20,0,20])
ax0d.set_xticklabels(['-r', '0', 'r'])
ax1.set_xticks([])
ax3.set_xticks([])
ax5.set_xticks([])
ax7.set_xticks([])

ax1.set_xlim([0,40])
ax2.set_xlim([0,40])
ax3.set_xlim([0,40])
ax4.set_xlim([0,40])
ax5.set_xlim([0,40])
ax6.set_xlim([0,40])
ax7.set_xlim([0,40])
ax8.set_xlim([0,40])

ylim1=[0,8]
ylim2=[0,20]
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
ax0a.set_yticks([0,1500/2,1200,1500])
ax0a.set_yticklabels(['0','H/2','h_eq','H'])
ax0b.set_yticks([])
ax0c.set_yticks([])
ax0d.set_yticks([])
ax3.set_yticks([])
ax4.set_yticks([])
ax5.set_yticks([])
ax6.set_yticks([])
ax7.set_yticks([])
ax8.set_yticks([])

ax0a.set_ylabel('z')
ax1.set_ylabel('$\\tau_{osc}$ (days)')
ax2.set_ylabel('$\\tau_{damp}$ (days)')
ax4.set_xlabel('Distance from margin (km)')  

elements_ax0a = [Line2D([0], [0],  color=colors[0], lw=lw, label=radius_cylinder[0]),
                 Line2D([0], [0],  color=colors[1], lw=lw, label=radius_cylinder[1]),
                 Line2D([0], [0],  color=colors[2], lw=lw, label=radius_cylinder[2]),
                 Line2D([0], [0],  color=colors[3], lw=lw, label=radius_cylinder[3]),
                 Line2D([0], [0],   color=colors[4], lw=lw, label=radius_cylinder[4])]

elements_ax0b = [Line2D([0], [0],  color=colors[0], lw=lw, label=slope_h2[0]),
                 Line2D([0], [0],  color=colors[1], lw=lw, label=slope_h2[1]),
                 Line2D([0], [0],  color=colors[2], lw=lw, label=slope_h2[2]),
                 Line2D([0], [0],  color=colors[3], lw=lw, label=slope_h2[3]),
                 Line2D([0], [0],   color=colors[4], lw=lw, label=slope_h2[4])]

elements_ax0c = [Line2D([0], [0],  color=colors[0], lw=lw, label=slope_heq[0]),
                 Line2D([0], [0],  color=colors[1], lw=lw, label=slope_heq[1]),
                 Line2D([0], [0],  color=colors[2], lw=lw, label=slope_heq[2]),
                 Line2D([0], [0],  color=colors[3], lw=lw, label=slope_heq[3]),
                 Line2D([0], [0],   color=colors[4], lw=lw, label=slope_heq[4])]

elements_ax0d = [Line2D([0], [0],  color=colors[0], lw=lw, label=slope_losange[0]),
                 Line2D([0], [0],  color=colors[1], lw=lw, label=slope_losange[1]),
                 Line2D([0], [0],  color=colors[2], lw=lw, label=slope_losange[2]),
                 Line2D([0], [0],  color=colors[3], lw=lw, label=slope_losange[3]),
                 Line2D([0], [0],   color=colors[4], lw=lw, label=slope_losange[4])]

ax0a.legend(handles=elements_ax0a, loc=1,  labelspacing=0, handletextpad=0.5,  
       prop={'size': 5}, title='radius',title_fontsize=8, bbox_to_anchor=(0.35, 0.8))
ax0b.legend(handles=elements_ax0b, loc=1,  labelspacing=0, handletextpad=0.5,  
       prop={'size': 5}, title='slope',title_fontsize=8, bbox_to_anchor=(0.35, 0.8))
ax0c.legend(handles=elements_ax0c, loc=1,  labelspacing=0, handletextpad=0.5,  
       prop={'size': 5}, title='slope',title_fontsize=8, bbox_to_anchor=(0.35, 0.8))
ax0d.legend(handles=elements_ax0d, loc=1,  labelspacing=0, handletextpad=0.5,  
       prop={'size': 5}, title='slope',title_fontsize=8, bbox_to_anchor=(0.35, 0.8))

pos2 = 0.85

ax0a.text(0.1, 1, 'a', ha='right', va='top', transform=ax0a.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax0b.text(0.1,1,'b', ha='right', va='top', transform=ax0b.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax0c.text(0.1, 1,'c', ha='right', va='top', transform=ax0c.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax0d.text(0.1, 1, 'd', ha='right', va='top', transform=ax0d.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax1.text(0.1, pos2,'e', ha='right', va='top', transform=ax1.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax2.text(0.1, pos2,'f', ha='right', va='top', transform=ax2.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax3.text(0.1, pos2, 'g', ha='right', va='top', transform=ax3.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax4.text(0.1, pos2,'h', ha='right', va='top', transform=ax4.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax5.text(0.1, pos2,'i', ha='right', va='top', transform=ax5.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax6.text(0.1, pos2,'j', ha='right', va='top', transform=ax6.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax7.text(0.1, pos2,'k', ha='right', va='top', transform=ax7.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
ax8.text(0.1, pos2,'l', ha='right', va='top', transform=ax8.transAxes,fontsize=8, 
          bbox=dict(facecolor='white', edgecolor='none', pad=1.0))


ax0a.set_title('Cylinder', fontsize=10)
ax0b.set_title('Cone (common \n radii at H/2)', fontsize=10)
ax0c.set_title('Cone \n(equilibrium head)', fontsize=10)
ax0d.set_title('Diamond-Hourglass \n(equilibrium head)', fontsize=10)

for i in np.arange(5): 

    #plot moulin shapes
    plot_moulin(ax0a,r_vector_cylinder,z_vector_cylinder)
    plot_moulin(ax0b,r_vector_H2,z_vector_H2)
    plot_moulin(ax0c,r_vector_heq,z_vector_heq)
    plot_moulin(ax0d,r_vector_losange,z_vector_losange)
    
    cond_cylinder = np.array(alpha_fit_cylinder[i]) != 9999
    cond_h2 = np.array(alpha_fit_H2[i]) != 9999
    cond_heq = np.array(alpha_fit_heq[i]) != 9999
    cond_losange = np.array(alpha_fit_losange[i]) != 9999
    
    '''Change radius with fixed slope'''
    #ax1.plot(profile/1000,oscillation_cylinder[i],color=colors[i], lw=lw, linestyle=linestyles[i])
    ax1.plot(profile[cond_cylinder]/1000,np.array(oscillation_cylinder[i])[cond_cylinder],color=colors[i], lw=lw, linestyle=linestyles[i])
    #ax1.plot(profile[cond_cylinder]/1000,np.array(oscillation_cylinder[i])[cond_cylinder],'.', color=colors[i])
    ax2.plot(profile[cond_cylinder]/1000,np.array(damping_cylinder[i])[cond_cylinder],color=colors[i], lw=lw, linestyle=linestyles[i])
    
    '''Change slope with h_middle fixed'''
    ax3.plot(profile[cond_h2]/1000,np.array(oscillation_H2[i])[cond_h2],color=colors[i], lw=lw, linestyle=linestyles[i])
    ax4.plot(profile[cond_h2]/1000,np.array(damping_H2[i])[cond_h2],color=colors[i], lw=lw, linestyle=linestyles[i])
    
    '''Change slope with h_eq fixed'''
    ax5.plot(profile[cond_heq]/1000,np.array(oscillation_heq[i])[cond_heq],color=colors[i], lw=lw, linestyle=linestyles[i])
    ax6.plot(profile[cond_heq]/1000,np.array(damping_heq[i])[cond_heq],color=colors[i], lw=lw, linestyle=linestyles[i])

    '''Change radius with h_eq fixed for losange'''
    ax7.plot(profile[cond_losange]/1000,np.array(oscillation_losange[i])[cond_losange],color=colors[i], lw=lw, linestyle=linestyles[i])
    ax8.plot(profile[cond_losange]/1000,np.array(damping_losange[i])[cond_losange],color=colors[i], lw=lw, linestyle=linestyles[i])


#plt.subplots_adjust(left=0.08, right=0.99, top=0.96, bottom=0.14)


plt.savefig('timescales_across_icesheet.pdf')

#%%


fig,axs = plt.subplots(5,4,figsize=(7.48,6),dpi=300)
for i in np.arange(5):
    cond_cylinder = np.array(alpha_fit_cylinder[i]) != 9999
    cond_h2 = np.array(alpha_fit_H2[i]) != 9999
    cond_heq = np.array(alpha_fit_heq[i]) != 9999
    cond_losange = np.array(alpha_fit_losange[i]) != 9999
    
    axs[0, 0].plot(profile[cond_cylinder]/1000,np.array(damping_cylinder[i])[cond_cylinder],'.')
    axs[1, 0].plot(profile[cond_cylinder]/1000,np.array(oscillation_cylinder[i])[cond_cylinder],'.')
    axs[2, 0].plot(profile[cond_cylinder]/1000,np.array(alpha_fit_cylinder[i])[cond_cylinder],'.')
    axs[3, 0].plot(profile[cond_cylinder]/1000,np.array(beta_fit_cylinder[i])[cond_cylinder],'.')
    axs[4, 0].plot(profile[cond_cylinder]/1000,np.array(alpha_fit_cylinder[i])[cond_cylinder]/np.array(beta_fit_cylinder[i])[cond_cylinder],'.')

    axs[0, 1].plot(profile[cond_h2]/1000,np.array(damping_H2[i])[cond_h2],'.')
    axs[1, 1].plot(profile[cond_h2]/1000,np.array(oscillation_H2[i])[cond_h2],'.')
    axs[2, 1].plot(profile[cond_h2]/1000,np.array(alpha_fit_H2[i])[cond_h2],'.')
    axs[3, 1].plot(profile[cond_h2]/1000,np.array(beta_fit_H2[i])[cond_h2],'.')
    axs[4, 1].plot(profile[cond_h2]/1000,np.array(alpha_fit_H2[i])[cond_h2]/np.array(beta_fit_H2[i])[cond_h2],'.')
    
    axs[0, 2].plot(profile[cond_heq]/1000,np.array(damping_heq[i])[cond_heq],'.')
    axs[1, 2].plot(profile[cond_heq]/1000,np.array(oscillation_heq[i])[cond_heq],'.')    
    axs[2, 2].plot(profile[cond_heq]/1000,np.array(alpha_fit_heq[i])[cond_heq],'.')
    axs[3, 2].plot(profile[cond_heq]/1000,np.array(beta_fit_heq[i])[cond_heq],'.')
    axs[4, 2].plot(profile[cond_heq]/1000,np.array(alpha_fit_heq[i])[cond_heq]/np.array(beta_fit_heq[i])[cond_heq],'.')

    axs[0, 3].plot(profile[cond_losange]/1000,np.array(damping_losange[i])[cond_losange],'.')
    axs[1, 3].plot(profile[cond_losange]/1000,np.array(oscillation_losange[i])[cond_losange],'.')
    axs[2, 3].plot(profile[cond_losange]/1000,np.array(alpha_fit_losange[i])[cond_losange],'.')
    axs[3, 3].plot(profile[cond_losange]/1000,np.array(beta_fit_losange[i])[cond_losange],'.')
    axs[4, 3].plot(profile[cond_losange]/1000,np.array(alpha_fit_losange[i])[cond_losange]/np.array(beta_fit_losange[i])[cond_losange],'.')

axs[0,0].set_ylabel('damping')
axs[1,0].set_ylabel('oscillation')
axs[2,0].set_ylabel('alpha')
axs[3,0].set_ylabel('beta')
axs[4,0].set_ylabel('alpha/beta')

axs[2,1].set_xlabel('distance from margin (km)')

axs[0,0].set_title('Cylinder')
axs[0,1].set_title('H/2')
axs[0,2].set_title('heq')
axs[0,3].set_title('diamond-hourglass')

plt.savefig('comparison_alpha_beta_test')
#%%
for r in [5,10,15]:
    plt.figure()
    for i in np.arange(20):
        resdic = frplus.calc_sim(r_fix=r,L=profile[i],profile=True)
        plt.plot(resdic['td'],resdic['hd'], label=profile[i]/1000)
        if resdic['alpha_fit']!=9999:
            plt.plot(resdic['td'],resdic['hd_fit'], color='black', linestyle='--')
        plt.ylim([0,600])
    plt.legend(title='km from margin', fontsize=6)
    plt.ylabel('head (m)')
    plt.xlabel('time (days)')
    plt.title('Cylinder, r=%s'%r)
    plt.savefig('head_cylinder_r%s.pdf'%r)
#%%
for m in [-0.02,0.02]:
    plt.figure()
    for i in np.arange(20):   
        resdic = frplus.calc_sim(r_fix=10,m=m,L=profile[i],profile=True,z_fix='H2')
        plt.plot(resdic['td'],resdic['hd'], label=profile[i]/1000)
        if resdic['alpha_fit']!=9999:
            plt.plot(resdic['td'],resdic['hd_fit'], color='black', linestyle='--')
        plt.ylim([0,600])
    plt.legend(title='km from margin', fontsize=6)
    plt.ylabel('head (m)')
    plt.xlabel('time (days)')
    plt.title('H2, m=%s'%m)
    plt.savefig('head_H2_m%s.pdf'%m)

for m in [-0.06,0.06]:
    plt.figure()
    for i in np.arange(20):   
        resdic = frplus.calc_sim(r_fix=10,m=m,L=profile[i],profile=True,z_fix='heq')
        plt.plot(resdic['td'],resdic['hd'], label=profile[i]/1000)
        if resdic['alpha_fit']!=9999:
            plt.plot(resdic['td'],resdic['hd_fit'], color='black', linestyle='--')
        plt.ylim([0,600])
    plt.legend(title='km from margin', fontsize=6)
    plt.ylabel('head (m)')
    plt.xlabel('time (days)')
    plt.title('heq, m=%s'%m)
    plt.savefig('head_heq_m%s.pdf'%m)

# plt.figure()
# resdic = frplus.calc_sim(r_fix=r,L=profile[59],profile=True, initial_ratio=1.1)
# plt.plot(resdic['t']/3600,resdic['h'])

# #%%
# plt.figure()
# for i in np.arange(5):
#     cond_cylinder = np.array(alpha_fit[i]) != 9999
#     plt.plot(profile[cond_cylinder]/1000,np.array(alpha_fit[i])[cond_cylinder],'.',color=colors[i])#, lw=lw, linestyle=linestyles[i])
#     #ax1.plot(profile[cond_cylinder]/1000,np.array(oscillation_cylinder[i])[cond_cylinder],color=colors[i], lw=lw, linestyle=linestyles[i])
#     #plt.plot(profile[cond_cylinder]/1000,np.array(beta_fit[i])[cond_cylinder],'o',color=colors[i])#, lw=lw, linestyle=linestyles[i])

# plt.figure()
# for i in np.arange(5):
#     cond_cylinder = np.array(beta_fit[i]) != 9999
#     plt.plot(profile[cond_cylinder]/1000,np.array(beta_fit[i])[cond_cylinder],'o',color=colors[i])#, lw=lw, linestyle=linestyles[i])
    
# plt.figure()
# for i in np.arange(5):
#     cond_cylinder = np.array(beta_fit[i]) != 9999
#     plt.plot(profile[cond_cylinder]/1000,np.array(beta_fit[i])[cond_cylinder]/np.array(alpha_fit[i])[cond_cylinder],'o',color=colors[i])#, lw=lw, linestyle=linestyles[i])

