

## Compare Moush (0D), fixed cylinder (0D), and fixed cylinder (1D)
import numpy as np
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import seaborn as sns
from scipy.integrate import solve_ivp


#for the 1D simulation. This requires to have Matt's modules installed
#from conduits1D_landlab_matt_02112021 import run1Dsim,plot_3panels_singlestep, plot_3panels_movie, plot_2panel_overtime_multiposition #for the 1D simulation. This requires to have Matt's modules installed
#import Parameters as param
import Constant_JGRpaper as cst

from onedim_conduit_model_matt_20220505 import one_dim_sim


secinday = 24*3600
ZERO_KELVIN = 273.15

#Glacier parameters
moulin_radii = 1. #m (based out of the highest radius when running moush for JEME)
ice_thickness = 1000 #m
initial_head = 900 #m
initial_subglacial_area = 1 #m
channel_length = 30e3 #m 

Qin_mean = 3
dQ = 0.3

#%% Simulation with Matt new 1D subglacial channel model

h_1d_moulin = []
h_1d_margin = []
s_1d_moulin = []
s_1d_margin = []
t_1d = []

params = {'R_func': 'sin', 
          'R_mean': Qin_mean, 
          'R_min': Qin_mean - dQ, 
          'A_R': np.pi * moulin_radii**2, 
          'L': channel_length, 
          'h0_init': initial_head,
          "S0": initial_subglacial_area,  # Conduit area in m^2
          "dt": 500,  # time in seconds
          }
sim = one_dim_sim(params = params)
nsteps = 200000
for i in np.arange(nsteps):
     sim.run_one_step(dt=50)
     h_1d_moulin.append(sim.h_node[0])
     h_1d_margin.append(sim.h_node[-2])
     s_1d_moulin.append(sim.S[0])
     s_1d_margin.append(sim.S[-2])
     t_1d.append(sim.t)
     
t_1d = np.array(t_1d)
     
#%% Simulation with 0D
t0= 0 #starting time
tf= 30 *secinday #(s) indicate number of days until solver stops

def dy_dt(t,y):
    # print(f'y: {y}')
    h = y[0] #(–) S*: Non dimentionalized channel cross-section area
    S = y[1] #(–) h*: Non dimentionalized moulin head
    A_R = np.pi*moulin_radii**2
    L=channel_length
    H=ice_thickness
    Pi = cst.rhoi*cst.g*H 
    
    #R(t)
    Rt =   dQ * np.sin(2.*np.pi*t/secinday) + Qin_mean
    #Head partial differential equation
    dhdt = 1/A_R * ( Rt - cst.C3*S**(5/4)*np.sqrt((cst.rhow*cst.g*h)/L) )
    #Channel cross section area partial differential equation
    dSdt = cst.C1 * cst.C3 * S**(5/4) * ((cst.rhow*cst.g*h)/L)**(3/2) - cst.C2 * ( Pi - cst.rhow*cst.g*h )**cst.n * S
    return [dhdt, dSdt]

sol = solve_ivp(dy_dt,(t0, tf),(initial_head,initial_subglacial_area),method = 'LSODA',max_step = 10*60)
h_0d = sol.y[0]
S_0d = sol.y[1]
t_0d = sol.t



#%% Plot head for all simulations
#xlim = [195,200]
#xlim = [180,250]
Qin = dQ * np.sin(2.*np.pi*t_0d/secinday) + Qin_mean
Qin_1d = dQ * np.sin(2.*np.pi*t_1d/secinday) + Qin_mean

fig,ax = plt.subplots(3,1)#, sharex=True)
#plt.figure()

#Plot Qin
ax[0].plot(t_0d/secinday,Qin, color='grey')
ax[0].set_ylabel('Q$_{in}$ (m$^3$/s)')
ax[0].tick_params(labelbottom=False)
#Oax[0].set_xticks([])
ax[0].set_ylim([2.4,3.6])
ax[0].set_yticks([2.4,2.6,2.8,3,3.2,3.4,3.6])
ax[0].set_xlim([4,10])
ax[0].set_xticks([])

#Plot head
ax[1].plot(t_0d/secinday,h_0d, linestyle='-',label='0D (moulin)')
ax[1].plot((t_1d-300)/secinday,h_1d_moulin,linestyle='--',label='1D (moulin)')
ax[1].set_ylim([500,1000])
ax[1].set_yticks([500,600,700,800,900,1000])
ax[1].set_ylabel('h (m)')
ax[1].legend(loc=4, prop={'size': 6}, bbox_to_anchor=(0.2,-0.2))
ax[1].tick_params(labelbottom=False)
ax[1].set_xlim([4,10])
ax[1].set_xticks([])

#Plot subglacial channel
ax[2].plot(t_0d/secinday, S_0d,linestyle='-',label='0D (moulin)')
ax[2].plot(t_1d/secinday,s_1d_moulin,linestyle='--',label='1D (moulin)')
ax[2].set_ylim([1,1.6])
ax[2].set_xlim([4,10])
ax[2].set_yticks([1,1.2,1.4,1.6])
ax[2].set_xticks([4,5,6,7,8,9,10])
ax[2].set_ylabel('S (m$^2$)')
ax[2].set_xlabel('Days')

ax[0].text(0.3,3.8,'(a)',fontsize=8)
ax[1].text(0.3,1100,'(b)',fontsize=8)
ax[2].text(0.3,1.3,'(c)',fontsize=8)

sns.despine(ax=ax[0],bottom=True)#trim=True)
sns.despine(ax=ax[1],bottom=True)
sns.despine(ax=ax[2],bottom=False, offset=(0,10))

plt.savefig('compare_0D1D.pdf')
plt.savefig('compare_0D1D.png')

# add subplot for subglacial radius?
# add subplot for moulin radius at h?

#%%

# #plt.figure()
# df_fix = {'t_fix':time/secinday, 'head_fix': head_fix}
# df_disc = {'t_disc':t_1d2/secinday, 'head_disc': h_1d_moulin2}

# df_fix = pd.DataFrame(df_fix)
# df_disc = pd.DataFrame(df_disc)

# # fig,ax = plt.subplots()
# # ax.plot(t_1d2/secinday,head_fix-h_1d_moulin2)#,label='1D (moulin)')
# # ax.set_ylabel('h (m)')













