import FR_mainCode_plus as frplus
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle #to unpickle
import seaborn as sns
import itertools

#Compare fit for different initial values

# CYLINDER
reslist_fit = [] 
initial_ratio_list = np.arange(0.5,1.5,0.1)
r_list = np.arange(1,6,1)
tau_damp = []
tau_osc = []
alpha = []
beta = []

#palette = itertools.cycle(sns.color_palette("Spectral"),n_colors=len(initial_ratio_list))
sns.set_palette(sns.color_palette("icefire",n_colors=len(initial_ratio_list)))

#plt.figure()
for r in r_list:
    col_damp = []
    col_osc =[]
    col_alpha = []
    col_beta = []
    plt.figure()
    plt.ylabel('head')
    plt.title('Radius = %s'%r)
    plt.xlim([0,4])

    for initial_ratio in [1]: #initial_ratio_list:   
        resdic = frplus.calc_sim(shape='linear',
                                r_fix=r,
                                initial_ratio=initial_ratio)
        #reslist_fit.append(resdic) # add dictionnary to the bottom of the list
        #print(resdic['damping_fit'])
        col_damp.append(resdic['damping_fit'])
        col_osc.append(resdic['oscillation_fit'])
        col_alpha.append(resdic['alpha_fit'])
        col_beta.append(resdic['beta_fit'])
        plt.plot(resdic['td'],resdic['hd'],label='initial_ratio =%s' %np.round(initial_ratio,2))#, color=next(palette))
    plt.legend()
    plt.savefig('initial_ratio_radius_%sm.pdf'%r)
    plt.savefig('initial_ratio_radius_%sm.png'%r)
        
    tau_damp.append(col_damp)
    tau_osc.append(col_osc)
    alpha.append(col_alpha)
    beta.append(col_beta)
#plt.plot(tau_damp)

#transform to  a numpy array to be able to transpose
tau_damp_t = np.asarray(tau_damp)
tau_osc_t = np.asarray(tau_osc)
alpha = np.asarray(alpha)
beta = np.asarray(beta)

#%%
sns.set_palette(sns.color_palette("ch:s=-.2,r=.6",n_colors=len(r_list)))

plt.figure(figsize=(6,4))
plt.plot(initial_ratio_list, tau_damp_t.T,'o', linestyle='--')
plt.xlabel('initial ratio')
plt.ylabel('Tau damp')
plt.legend(r_list, title='Moulin radius (m)')
plt.show()
plt.savefig('Comparison_Initial_ratio_TAUdamp.pdf')
plt.savefig('Comparison_Initial_ratio_TAUdamp.png')

plt.figure()
plt.plot(initial_ratio_list, tau_osc_t.T,'o', linestyle='--')
plt.xlabel('initial ratio')
plt.ylabel('Tau osc')
plt.legend(r_list)
plt.legend(r_list, title='Moulin radius (m)')
plt.show()
plt.savefig('Comparison_Initial_ratio_TAUosc.pdf')
plt.savefig('Comparison_Initial_ratio_TAUosc.png')


plt.figure()
plt.plot(initial_ratio_list, alpha.T,'o', linestyle='--')
plt.xlabel('initial ratio')
plt.ylabel('alpha')
plt.legend(r_list)
plt.legend(r_list, title='Moulin radius (m)')
plt.show()
plt.savefig('Comparison_Initial_ratio_alpha.pdf')
plt.savefig('Comparison_Initial_ratio_alpha.png')

plt.figure()
plt.plot(initial_ratio_list, beta.T,'o', linestyle='--')
plt.xlabel('initial ratio')
plt.ylabel('beta')
plt.legend(r_list)
plt.legend(r_list, title='Moulin radius (m)')
plt.show()
plt.savefig('Comparison_Initial_ratio_beta.pdf')
plt.savefig('Comparison_Initial_ratio_beta.png')