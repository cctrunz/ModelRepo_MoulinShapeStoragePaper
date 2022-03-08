""" Extract timescales of oscillation and decay (alpha and beta) for CONE shaped moulin

Use the non-dimentionalised version of Schoof 2010 model and the the approximation.
Both are derived by Matt Covington. This code is not intended to be used long term,
but it is intended to serve as an easier platform to implement the changes from a
constant moulin cross-section area to one depending on the head

Code based on Matt Covington code. Modified by Celia Trunz. Fall 2019, Buffalo NY

"""
#%%

#Import packages
from scipy.integrate import solve_ivp #odeint
from scipy.optimize import root#,curve_fit
#from mpl_toolkits.axes import Axes

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import cm
#import matplotlib.tri as tri
import pickle
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
#from cycler import cycler
#from matplotlib.cm import get_cmap
#from mpltools import color


#import os #so that pickle create the file in the right folder


"""
Fixed Parameters
"""
#Physical constants
rhow = 1000 #water density, kg/m3
rhoi = 910 #ice density $, kg/m3
g =  9.8 #Gravity
f =  0.1 #Darcy weissbach friction factor #different in Schoof(2010)
Lf = 3.32e5 #Latent heat of fusion (from schoof2010) J/kg  #(CT) 3.35e5 in schoof(2010)
A =  6e-24 #Pa^-3 s^-1
n =  3 #flow law parameter
#Calculated parameters
C1 = 1/(rhoi*Lf)
C2 = 2*A*n**-n
C3 = ( 2.**(5./4.)/np.pi**(1./4.) * np.sqrt(np.pi/(np.pi + 2.)) )/ np.sqrt(rhow*f)#(2**(1/4) * (np.pi+2)**(1/2)) / (np.pi**(1/4) * (rhow*f)**(1/2))


"""
Definition of non-dimentionalized equation
"""

#Modified from Matt Covington "Derivations of subglacial conduit timescales and approximate solutions", unpublished
def dy_dt(t,y):
    h = y[0] #(–) h*: Non dimentionalized moulin head
    S = y[1] #(–) S*: Non dimentionalized channel cross-section area

    #define AR* specifically for a cone
    AR_h_cone = np.pi * (m * h* hfl + r_base)**2 #(m) AR(h) for a cone
    AR = AR_h_cone / AR_heq #(–) 
#    #define AR* for a general function reproducing a realistic moulin shape
#    TBD

    #define R* for any function
    Rt = R_func(t) #(m3/s)3
    R = Rt/R_mean #(–)

    # Definition of the non-dimentional equation:
    dh_dt = (R - S**(5/4)*h**0.5)/AR #Moulin head oscillations
#    dh_dt = (R - S**(5/4)*h**0.5)/AR #Moulin head oscillations
    dS_dt = T1*S**(5/4)*h**(3/2) - T2*S*(1-h)**3 #Channel creep and closure
    return (dh_dt, dS_dt)

#Definition of no time equation
def dy_dt_no_time(y):
    return dy_dt(0,y)

#Definition of Timescales equations
def calcTauRes():
    return AR_heq*hfl/R_mean #(s)
def calcTauCreep():
    return 1/C2/Pi**3 #(1/s)
def calcTauMelt():
    return (L/Pi)**(7/5) * 1/C1/C3**(4/5)/R_mean**(1/5)# #(1/s)


#Definition of R(t) equation
def R_func(t):
    #R_constant
    return R_mean
    #R_sin
#    return (R_mean - R_min) * np.sin(2.*np.pi*t*TauRes/R_period) + R_mean
    #R_daily
    #return p['R_mean'] * (1. + np.sin(2.*np.pi*t/p['R_period'])/2.) #this come from one of matt's code
def calcTauFrequency(t,h):
    #AR = AR_func(h)
    #R = R_func(t)
    R = R_mean
    AR_h_cone = np.pi * (m * h* hfl + r_base)**2 #(m) AR(h) for a cone
    AR = AR_h_cone / AR_heq #(–) 
    return AR * Pi/rhow/g/R


#Function for sinusoidal recharge
#def R_sin(t):
#    (R_mean - R_min) * np.sin(2.*np.pi*t*TauRes/R_period) + R_mean
#def R_real(t): 
    


"""
Define elements
"""

reslist = [] # Create empty list
sec_in_day = 24*60*60 
# Define time vector
t0 = 0
tf = 100 #(non dim time) 
initial_ratio = 1.1
#Define profile vector for a realistic glacier
H0 = 1500./np.sqrt(60000.) #assures ice thickness of 1500 m at 60 km from edge
L_profile = np.arange(0,61000,1000)
#L_profile = [30000]#Define recharge (should we define recharge relative to the moulin size maybe?)
#R_vec = np.arange(1,300,10) #(m3/s)  
#Define recharge in function of t
R_vec = [3] #m^3/s, mean discharge into moulin
#R_min = 1 #Minimum daily recharge (used for sine func)
#R_period = sec_in_day #Period for oscillatory Recharge
#radius at equilibrium (m)
#r_heqvec = np.arange(1,31,1) 
#r_heqvec = [20] #(m)
#mean radius at the middle of the moulin (m)
#r_hmiddlevec = np.arange(1,31,1) 
r_hmiddlevec = [20] #(m)
#radius at equilibrium (m)
#Moulin slope (with r is y and z(r) )
#m_vec = np.array([-0.2, -0.15, -0.1, -0.08, -0.05,  
#                  0,
#                  0.05, 0.08,  0.1,  0.12])## Option 2 #Negative slope means base is bigger than top. 
#
#m_vec = np.linspace(-0.2,0.12,50)
m_vec = np.array([-1.15, -0.57, 0, 0.57, 1.15]) # cone H/2
#m_vec = np.linspace(-0.05,0.05,50)
#m_vec = np.array([-0.05,-0.01,-0.001, 
#                  0,
#                  0.001, 0.01, 0.05 ])## Option 2 #Negative slope means base is bigger than top. 
#m_vec = [0]

# res_filename = 'results_TauProfile_slope_hmiddle_AGU19' #change file name each time and describe what it contains in overleaf
# para_filename = 'param_TauProfile_slope_hmiddle_AGU19' #change file name each time and describe what it contains in overleaf



"""
CALCULATE OUTPUTS
"""
#fig, ax = plt.subplots(figsize=(10,4))
# or
#fig, (ax1, ax2) = plt.subplots(1, 2,  figsize=(15,5),sharey=True, gridspec_kw={'width_ratios': [3, 1], 'wspace': 0}) #, 'wspace': 0 , sharey=True,
## Define the colors to be used using rainbow map (or any other map)
#colors = [plt.cm.rainbow(i) for i in np.linspace(0, 1, len(m_vec))]
# or
fig = plt.figure(figsize=(15,5))
grid = plt.GridSpec(3, 3, wspace=0, hspace=0)
ax1 = fig.add_subplot(grid[0:, :2])#,xticklabels=[]
ax2 = fig.add_subplot(grid[0:, 2], sharey=ax1)# , yticklabels=[])
ax3 = fig.add_subplot(grid[2, :2], sharex=ax1)
#ax4 = fig.add_subplot(grid[3, :2], sharex=ax1)
#ax5 = fig.add_subplot(grid[0:, 3])

colors = [plt.cm.rainbow(i) for i in np.linspace(0, 1, len(m_vec))] 

#plt.rc('font', size=14)          # controls default text sizes
#plt.rc('axes', titlesize=14)     # fontsize of the axes title
#plt.rc('axes', labelsize=14)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
#plt.rc('legend', fontsize=10)    # legend fontsize
##plt.rc('figure', titlesize=14)  # fontsize of the figure title


#for r_heq in r_heqvec:
for r_hmiddle in r_hmiddlevec:    
    #print('r_heq =', r_heq)
    for R_mean in R_vec:
#        for L in L_profile: #for timeseries
        for i,m in zip(range(len(m_vec)),m_vec): #for profile
#            print('m = ', m)
#            for i,m in zip(range(len(m_vec)),m_vec): #for timeseries
            for L in L_profile: #for profile
                #Calculate variable depending on L and z       
                #z = round(H0*np.sqrt(L),-1) # ice thickness for a square root glacier.
                z = H0*np.sqrt(L) # ice thickness for a square root glacier.
                Pi = rhoi*g*z #Units?? Ice pressure   
                hfl = Pi/rhow/g #(m) head at overburden pressure. to dimentionalize h        
            
                #Define initial AR_heq, T1 and T2. h_eq and S_eq do not depend on the moulin radius or shape.
                AR_heq = 1 #this in only to enable to solve for h_eq and S_eq. this value is replace later.
                r_base = 1 # option 4. This in only to enable to solve for h_eq and S_eq. this value is replace later.
                T1 = calcTauRes() / calcTauMelt() # 6e-8 – 6 #(–)
                T2 = calcTauRes() / calcTauCreep() # 6e-6 – 41200 #(–)

                """
                Solve equations for h*, S*
                """   
                #Calculate approximation of h_eq* and S_eq* equlibrium
                h_eq_approx = 1/((T1/T2)**(5/7) + 15/7) #(–)
                S_eq_approx = h_eq_approx**(-2/5) #(–)
                
                ##Calcule de h* et S* a l'équilibre
                eq_sol = root(dy_dt_no_time, (h_eq_approx,S_eq_approx))
                h_eq = eq_sol.x[0] #(–)
                S_eq = eq_sol.x[1] #(–)
                
                ##print(h_eq)
                ##Calculate r_base and r_top based on r_heq.
                #r_base = r_heq - m*h_eq*hfl #(m) hfl is used to dimensionalized h_eq
                #r_top = m*z + r_base #(m)
               
                #Calculate r_base and r_top based on r_hmiddle.
                r_base = r_hmiddle - m*(z/2)#(m)
                r_top = m*z + r_base #(m) #(m)
                r_heq = m * h_eq *hfl + r_base #(m) hfl is used to dimensionalized h_eq
                
                #Find the cross-section area at h_eq
                AR_heq =  np.pi * r_heq**2 #(m^2)
                
                #Find T1 and T2 with real h_eq
                T1 = calcTauRes() / calcTauMelt() #(–)
                T2 = calcTauRes() / calcTauCreep() #(–)
                
                ##Calculate h* and S*   
                h0 = h_eq*initial_ratio
                S0 = S_eq*initial_ratio
                
                sol = solve_ivp(dy_dt, 
                                (t0, tf), #initial time and end time. !! in non-dim!                                (param['t0'], param['tf']), #initial time and end time. !! in non-dim!
                                (h0,S0), #initial head and channel cross-section area. !!in non dim!
                                method = 'LSODA', #this method is the same 
                #                atol = 1e-6, #tolerance. can make it faster
                #                rtol = 1e-3,
                                max_step = 0.01 #!! non-dim time step!
                                )

                hnd = sol.y[0]  #(non-dim) Head
                Snd = sol.y[1] #(non-dim) Channel cross-section
                tnd = sol.t # (non-dim) Time vector 
                
                ##Calculate alpha* and beta*
                a = 5/4 *T1 *S_eq**0.25 *h_eq**1.5 - T2 *(1-h_eq)**3
                b = 3/2 *T1 *S_eq**(5/4) *h_eq**0.5 + 3. *T2 *S_eq *(1-h_eq)**2
                c = -5/4 *h_eq**0.5 *S_eq**0.25
                d = -1/2 *h_eq**0.5 *S_eq**(5/4)
                
                p = a+d
                q = a*d - b*c
                
                alpha = p/2.
                beta = np.imag(np.sqrt(np.complex(p**2.-4*q)))/2.
                
                
                """ 
                Recovering dimensions and calculate timescales
                """
             
                #Calculate TauRes, TauMelt, TauCreep
                TauRes = calcTauRes()/sec_in_day #(days)
                TauCreep = calcTauCreep()/sec_in_day #(1/days)
                TauMelt = calcTauMelt()/sec_in_day # #(1/days)
                #Calculate t, h and S with dimension
                hd = hnd * hfl #(m) 
                Sd = Snd * (R_mean**2*L/C3**2/Pi)**(2/5) #(m) Sstar * S0
                td = tnd * TauRes # (days) to transform in days and set t0 to zero
                h_eq_d = h_eq*hfl
                #Calculate damping and oscillation with dimension
                damping = np.abs(1/alpha)  * TauRes #(day)
                oscillation = 2*np.pi/beta * TauRes #(day)
                #Calculate Q
                Q = C3*Sd**(5/4)*np.sqrt(rhow*g*hd/L)
                
                TauFreq = calcTauFrequency(td,hd)
                
                #    #Recovering AR
                #    hdummy = np.arange(1,z,10)
                #    AR_h_cone = np.pi * (m*hdummy + r_base)**2 #(m) AR(h) for a cone
                #    AR_test = AR_h_cone / AR_heq #(–) should use h_eq in futur
                
                #save in a dictionnary  #'r_mean':r_mean,
                resdic = {  'T1':T1,'T2':T2,\
                            'z':z, 'r_top':r_top, 'r_base':r_base, 'r_heq':r_heq, 'm':m, 'R':R_mean, 'L':L,\
                            'S0':S0, 'h0':h0,\
                            'S_eq':S_eq, 'h_eq':h_eq, 'S_eq_approx':S_eq_approx, 'h_eq_approx': h_eq_approx,\
                            'S':Sd, 'h':hd, 't':td, 'h_eq_d':h_eq_d,\
                            'Q':Q,\
                            'alpha':alpha, 'beta':beta,\
                            'damping':damping, 'oscillation':oscillation,\
                            'Sd':Sd, 'hd':hd, 'td':td,\
                            'TauCreep':TauCreep, 'TauRes':TauRes, 'TauMelt':TauMelt
                            } 
                #save to a list of dictionnary                   }
                #to call an element from the list: results_list['name of the variable'][position of value of vector if needed]
                #to extract all the values from a single key: [result['KeyName'] for result in reslist] -- change KeyName with the needed key
                reslist.append(resdic) # add dictionnary to the bottom of the list
                                 
#%% To plot head timeseries for various slope and displaying the moulin shape
                
#                #ax1.plot(td,hd,label='$m=$%s, $r_{top}=$%s, $r_{base}=$%s, $r_{heq}=$%s'%(m,round(r_top),round(r_base),round(r_heq)),color=colors[i]) 
#                ax1.plot(td,hd,label='$m=$%s'%(m),color=colors[i]) 
#                #ax1.set_ylim([650,880])
#                ax1.set_ylim([0,z])
#                ax1.set_xlim([-1,20])
#                #ax1.legend(loc='lower right')
#                ax1.set_xlabel('Days')
#                ax1.set_ylabel('Head (m)')
#                #ax1.set_title('Moulin head')
#                
#                z_plot = np.linspace(0,z)
#                r_plot = m*z_plot + r_base
#                ax2.plot(r_plot,z_plot,color=colors[i])
#                ax2.plot(-r_plot,z_plot,color=colors[i])
#                ax2.set_xlim([-100,100]) 
#                #ax2.set_ylim([0,z])
#                ax2.yaxis.tick_right()
#                ax2.set_xlabel('Radius from center (m)')
#                ax2.set_title('Moulin shape')
#                
#                ax3.plot(td,Sd,color=colors[i])
#                #ax3.set_ylim([1.1,1.3])
#                #ax3.set_ylabel('S (m)') 
#                
##                ax4.plot(td,TauFreq,color=colors[i])
##                #ax4.set_ylabel('TauRes/24h (–)')
##                ax4.set_xlabel('Days')
##                #plt.ylim([650,870])
#                
#                #textbox:
#                textstr = '\n'.join((
#                    r'$R=%.f m^3/s$ ' %R_mean,
#                    r'$Z=%.f m$' %z,
#                    r'$L=%.f km$' %(L/1000)#,
#                    #r'$r_{h_middle}=%.f m$' %r_hmiddle
#                    ))
#                
#                # these are matplotlib.patch.Patch properties
#                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#                
#                # place a text box in upper left in axes coords
#                plt.text(0.8, 1, textstr, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props)
#                plt.show()
#                
#%% 
                #plt.plot(resdic['L'],resdic['damping'],'o')               
res_filename = 'results_TauProfile_H2' #change file name each time and describe what it contains in overleaf
para_filename = 'param_TauProfile_H2' #change file name each time and describe what it contains in overleaf               
# Save list of dictionary. !!Change name each time and move to Result_lists folder,

outfile = open(res_filename,'wb') #This means that the data will be written in the form of byte objects.
pickle.dump(reslist, outfile) #takes two arguments: the object you want to pickle and the file to which the object has to be saved
outfile.close() #close the file we just opened

# Save parameters 
#parameters = {'R':R_vec, 'r_hmiddle':r_hmiddlevec, #'r_heq':r_heqvec, 
#              'm':m_vec, 'L_profile': L_profile,  't0':t0, 'tf':tf, 'H0':H0, 'initial_ratio':initial_ratio}
#
#outfile = open(para_filename,'wb') #This means that the data will be written in the form of byte objects.
#pickle.dump(parameters, outfile) #takes two arguments: the object you want to pickle and the file to which the object has to be saved
#outfile.close() #close the file we just opened            
 #%%
#                #Figure to compare slopes              
#                plt.plot(td,hd,
#                         label='$m=$%s, $r_{top}=$%s, $r_{base}=$%s, $r_{heq}=$%s'%(m,round(r_top),round(r_base),round(r_heq)),
#                         color=colors[i])
#              
#                plt.legend()
#                plt.xlabel('Days')
#                plt.ylabel('Head (m)')
#                plt.ylim([0,z])
#                #textbox:
#                textstr = '\n'.join((
#                    r'$R=%.f m^3/s$ ' %R_mean,
#                    r'$Z=%.f m$' %z,
#                    r'$L=%.f km$' %(L/1000)#,
#                    #r'$r_{h_middle}=%.f m$' %r_hmiddle
#                    ))
#                
#                # these are matplotlib.patch.Patch properties
#                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#                
#                # place a text box in upper left in axes coords
#                plt.text(0.4, 0.95, textstr, transform=ax.transAxes, fontsize=12,
#                        verticalalignment='top', bbox=props)
#                plt.show()

#%%
#                #Figures to visualize moulin shape and oscillations
#                fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(15,5))    
#                
#                ax1.plot(td,hd,label='$m=$%s, $r_{top}=$%s, $r_{base}=$%s, $r_{heq}=$%s'%(m,round(r_top),round(r_base),round(r_heq)))             
#                ax1.legend()
#                ax1.set_xlabel('Days')
#                ax1.set_ylabel('Head (m)')
#                
#                z_plot = np.linspace(0,z)
#                r_plot = m*z_plot + r_base
#                ax2.plot(r_plot,z_plot)
#                ax2.plot(-r_plot,z_plot)
#                ax2.set_xlim([-50,50]) 
#                
#                plt.ylim([650,875])
#                
#                #textbox:
#                textstr = '\n'.join((
#                    r'$R=%.f m^3/s$ ' %R_mean,
#                    r'$Z=%.f m$' %z,
#                    r'$L=%.f km$' %(L/1000)#,
#                    #r'$r_{h_middle}=%.f m$' %r_hmiddle
#                    ))
#                
#                # these are matplotlib.patch.Patch properties
#                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#                
#                # place a text box in upper left in axes coords
#                plt.text(0.8, 0.20, textstr, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props)
#                plt.show()
                
                                #Figures to visualize moulin shape and oscillations
                                


#%% To plot Sd 
                
                #ax1.plot(td,hd,label='$m=$%s, $r_{top}=$%s, $r_{base}=$%s, $r_{heq}=$%s'%(m,round(r_top),round(r_base),round(r_heq)),color=colors[i]) 
#                ax1.plot(td,Sd,label='$m=$%s'%(m),color=colors[i])        
#                ax1.legend(loc='lower right')
#                ax1.set_xlabel('Days')
#                ax1.set_ylabel('Area ($m^2$)')
#                ax1.set_title('Channel cross-section')
                
#                z_plot = np.linspace(0,z)
#                r_plot = m*z_plot + r_base
#                ax2.plot(r_plot,z_plot,color=colors[i])
#                ax2.plot(-r_plot,z_plot,color=colors[i])
#                ax2.set_xlim([-50,50]) 
#                ax2.set_xlabel('Distance center (m)')
#                ax2.set_title('Moulin diameter')
#                
                
                #plt.ylim([650,870])
                
#                #textbox:
#                textstr = '\n'.join((
#                    r'$R=%.f m^3/s$ ' %R_mean,
#                    r'$Z=%.f m$' %z,
#                    r'$L=%.f km$' %(L/1000)#,
#                    #r'$r_{h_middle}=%.f m$' %r_hmiddle
#                    ))
#                
#                # these are matplotlib.patch.Patch properties
#                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#                
#                # place a text box in upper left in axes coords
#                plt.text(0.8, 0.9, textstr, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props)
#                plt.show()                
