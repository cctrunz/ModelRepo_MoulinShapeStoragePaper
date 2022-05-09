import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import seaborn as sns
#Hamemade package
import FR_mainCode_plus as frplus
import math

#%%
profile = np.arange(1000,61000,1000)
radius_cylinder = [5,7.5,10,12.5,15]
#radius_losange = [1,2,5,10,19]
# slope_h2 = [-0.02,
#           -0.01,
#           0,
#           0.01,
#           0.02]
# slope_heq = [-0.06,
#           -0.03,
#           0,
#           0.03,
#           0.06]
# slope_losange = [-0.06,
#           -0.03,
#           0,
#           0.03,
#           0.06]
#%%

# #Calculate simulation results for cylinders
# damping_cylinder = []
# oscillation_cylinder = []
# r_vector_cylinder = []
# z_vector_cylinder = []
# alpha_fit_cylinder = []
# for r in radius_cylinder:
#     col_damp = []
#     col_oscillation = []
#     col_alpha_fit = []
#     for L in profile:
#         resdic = frplus.calc_sim(r_fix=r,L=L,profile=True)
#         col_damp.append(resdic['damping_fit'])
#         col_oscillation.append(resdic['oscillation_fit'])
#         col_alpha_fit.append(resdic['alpha_fit'])
#     damping_cylinder.append(col_damp) 
#     oscillation_cylinder.append(col_oscillation)
#     r_vector_cylinder.append(resdic['r_vector'])
#     z_vector_cylinder.append(resdic['z_vector'])
#     alpha_fit_cylinder.append(col_alpha_fit)

# #%%

# plt.figure()
# for i in np.arange(5): 
#     cond_cylinder = np.array(alpha_fit_cylinder[i]) != 9999
#     plt.plot(profile[cond_cylinder]/1000,np.array(oscillation_cylinder[i])[cond_cylinder])
    
    
#%%
profile = np.arange(1000,40000,2000)
sns.set_palette(sns.color_palette("icefire",n_colors=len(profile)))


for r in  [5,7.5,10,12.5,15]:
    fig,ax = plt.subplots(1,2,figsize=(6,6))
    ax[-1].axis('off')
    for L in profile:#[5*1000]:
        res = frplus.calc_sim(r_fix=r,L=L, profile=True)
         
        # RMSE
        t = res['td'][res['td']<3*res['oscillation_fit']]
        y_actual = res['hd'][res['td']<3*res['oscillation_fit']]
        y_predicted = res['hd_fit'][res['td']<3*res['oscillation_fit']]
         
        # mse = np.square(np.subtract(y_actual,y_predicted)).mean() 
        # rmse = math.sqrt(mse)
        # print('Radius = ', r, 'm')
        # # print("RMSPE = ", rmse)
        
        rmspe = (np.sqrt(np.mean(np.square((y_actual - y_predicted) / y_actual)))) * 100
        #rmspe = np.sqrt(np.mean(np.square(((y_actual - y_predicted) / y_actual)), axis=0))
        print('RMSPE = ',np.round(rmspe,2),'%')
        if rmspe<1:
            # plt.figure()
            ax[0].plot(t, y_actual, linewidth=2, label='L=%skm, RMSPE=%spct'%(L/1000,np.round(rmspe,2)))
            ax[0].plot(t, y_predicted, linestyle='--', color='grey', linewidth=2) 
            #plt.text(1,max(y_actual),'RMSPE:%s percent'%(np.round(rmspe,2)))
    ax[0].set_xlabel('Days')
    ax[0].set_ylabel('Head (m)')
    ax[0].set_xlim([0,4])
    ax[0].set_ylim([0,1000])
    ax[0].legend(loc=1, bbox_to_anchor=(2.4,1))
    plt.title('Cylinder, r=%sm'%r)
    plt.savefig('Oscillation_RMSPE_r%sm.pdf'%r)
            #plt.tight_layout()
            #plt.xlim([0,3*res['damping_fit']])
    
#%%

#Calculate simulation results for cylinders
damping_cylinder = []
oscillation_cylinder = []
r_vector_cylinder = []
z_vector_cylinder = []
alpha_fit_cylinder = []
for r in radius_cylinder:
    col_damp = []
    col_oscillation = []
    col_alpha_fit = []
    for L in profile:
        resdic = frplus.calc_sim(r_fix=r,L=L,profile=True)
        t = resdic['td'][resdic['td']<3*resdic['damping_fit']]
        y_actual = resdic['hd'][resdic['td']<3*resdic['damping_fit']]
        y_predicted = resdic['hd_fit'][resdic['td']<3*resdic['damping_fit']]
        rmspe = (np.sqrt(np.mean(np.square((y_actual - y_predicted) / y_actual)))) * 100
        if rmspe < 1:
            col_damp.append(resdic['damping_fit'])
            col_oscillation.append(resdic['oscillation_fit'])
            col_alpha_fit.append(resdic['alpha_fit'])
        else:
            col_damp.append(np.nan)
            col_oscillation.append(np.nan)
            col_alpha_fit.append(np.nan)
    damping_cylinder.append(col_damp) 
    oscillation_cylinder.append(col_oscillation)
    r_vector_cylinder.append(resdic['r_vector'])
    z_vector_cylinder.append(resdic['z_vector'])
    alpha_fit_cylinder.append(col_alpha_fit)
    
    