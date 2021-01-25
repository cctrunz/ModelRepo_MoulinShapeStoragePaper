#%%
""" Extract timescales of oscillation and decay (alpha and beta) for CONE shaped moulin

Use the non-dimentionalised version of Schoof 2010 model and the the approximation.
Both are derived by Matt Covington. This code is not intended to be used long term,
but it is intended to serve as an easier platform to implement the changes from a
constant moulin cross-section area to one depending on the head

Code based on Matt Covington code. Modified by Celia Trunz. Fall 2019, Buffalo NY

"""
#%%
#Test editing by matt

#Import packages
from scipy.integrate import solve_ivp #odeint
from scipy.optimize import root,curve_fit
import pandas as pd
#from mpl_toolkits.axes import Axes

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from matplotlib import cm
#import matplotlib.tri as tri
import pickle

import matplotlib
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rcParams['axes.labelsize'] = 8
matplotlib.rcParams['axes.labelweight'] = 'bold'
matplotlib.rc('xtick', labelsize=8)     
matplotlib.rc('ytick', labelsize=8)

#colors:
Red = '#ED2224'
Orange = '#FBB263'
Green = '#A2D39E'
Blue = '#15B4E9'
Purple = '#6B52A2'


#%%



"""
DEFINE FUNCTIONS
"""
def calc_zr_rectangle(h_eq_nd,hfl,dh_above,dh_below,r1, r2):
    ''' This function produces a vector of moulin radius and a position vector 
    z, upward from the base of the moulin.
    '''   
    h_max = h_eq_nd*hfl + dh_above #dh are in m
    h_min = h_eq_nd*hfl - dh_below #dh are in m
    
    z = np.array([0, h_min-0.001, h_min, h_max])#, h_max+0.01, param.H])
    r = np.array([r1, r1,  r2, r2])#,  r_else,r_else])
    return [z,r]



#def interp_r(zi,z,r):
#    return np.interp(zi,z_moulin,r_moulin)    

def calc_r_linear(zi,m,r_fix,z_fix):
    r_base = r_fix-m*z_fix
    return m*zi+r_base

def calc_AR(r):
    return np.pi*r**2

#Modified from Matt Covington "Derivations of subglacial conduit timescales and approximate solutions", unpublished
def dy_dt(t,y,hfl,z,r,AR_heq,T1,T2): #r_fix,z_fix,
    h = y[0] #(–) h*: Non dimentionalized moulin head
    S = y[1] #(–) S*: Non dimentionalized channel cross-section area
    hd = h*hfl
    #print('hd=',hd,'z=',z,'r=',r)
    #print(type(z))
    #z=[1,2,3,4,5,6,7]      
    #r=[1,2,3,4,5,6,7]
    r_h = np.interp(hd,z,r) 
    AR_h= calc_AR(r_h)
    AR = AR_h/ AR_heq #(–) 
    dh_dt = (1 - S**(5/4)*h**0.5)/AR #Moulin head oscillations
    dS_dt = T1*S**(5/4)*h**(3/2) - T2*S*(1-h)**3 #Channel creep and closure
    
    return (dh_dt, dS_dt)



#Definition of no time equation
def dy_dt_no_time(y,T1_approx,T2_approx):
    return dy_dt(0,y,1,[1,1],[1,1],1,T1_approx,T2_approx)

#Definition of Timescales equations
def calcTauRes(AR_heq,hfl,R):
    return AR_heq*hfl/R #(s)
def calcTauCreep(C2,Pi):
    return 1/C2/Pi**3 #(1/s)
def calcTauMelt(L,Pi,C1,C3,R):
    return (L/Pi)**(7/5) * 1/C1/C3**(4/5)/R**(1/5)# #(1/s)




def calc_sim(   R = 3, #m^3/s, mean discharge into moulin
                r_fix = 5,#np.linspace(0.1,10,5)##
                r1 = 5,
                r2 = 5,
                z_fix = 'H2',
                shape = 'linear', #linear for cone, nodes for bottle and goblet
                m = 0,
                dh_below = 10, # distance below heq in meters!!  #dh2 =  , # distance above heq
                t0 = 0,
                tf = 100, #(non dim time) 
                initial_ratio = 1.1,
                profile = False,
                H_fix = 1000,
                L = 30000): #Define recharge (should we define recharge relative to the moulin size maybe?)
            
    
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
    sec_in_day = 24*60*60

    #Calculate variable depending on L and z       
    #z = round(H0*np.sqrt(L),-1) # ice thickness for a square root glacier.
    if profile == True:
        H0 = 1500/np.sqrt(60000) #assures ice thickness of 1500 m at 60 km from edge
        H = H0*np.sqrt(L) # ice thickness for a square root glacier.
    else:
        H = H_fix
    Pi = rhoi*g*H #Units?? Ice pressure   
    hfl = Pi/rhow/g #(m) head at overburden pressure. to dimentionalize h  
   
    
    #Define initial AR_heq, T1 and T2. h_eq and S_eq do not depend on the moulin radius or shape.
    T1_approx = calcTauRes(1,hfl,R) / calcTauMelt(L,Pi,C1,C3,R) # 6e-8 – 6 #(–)
    T2_approx = calcTauRes(1,hfl,R) / calcTauCreep(C2,Pi) # 6e-6 – 41200 #(–)

    
    """
    Solve equations for h*, S*
    """   
    #Calculate approximation of h_eq* and S_eq* equlibrium
    h_eq_approx = 1/((T1_approx/T2_approx)**(5/7) + 15/7) #(–)
    S_eq_approx = h_eq_approx**(-2/5) #(–)
    
    ##Calcule de h* et S* a l'équilibre
    eq_sol = root(dy_dt_no_time, (h_eq_approx,S_eq_approx), args=(T1_approx,T2_approx))
    h_eq_nd = eq_sol.x[0] #(–)
    S_eq_nd = eq_sol.x[1] #(–)   

    if shape == 'linear':
        #Find the cross-section area at h_eq
        if z_fix == 'H2':
            z_fix = H/2 #h_eq*hfl
        if z_fix == 'heq':
            z_fix = h_eq_nd*hfl
        r_top = calc_r_linear(H,m,r_fix,z_fix)
        r_base = calc_r_linear(0,m,r_fix,z_fix) 
        z_vector = [0,H]
        r_vector = [r_base,r_top]

    if shape == 'nodes':
        dh_above = H-h_eq_nd*hfl # distance from heq to top of the moulin, in meters
        z_vector,r_vector = calc_zr_rectangle(h_eq_nd,hfl,dh_above,dh_below,r1,r2)
        
    
    r_heq = np.interp(h_eq_nd*hfl,z_vector,r_vector) 
    AR_heq =  calc_AR(r_heq)
    
    #Find T1 and T2 with real h_eq
    T1 = calcTauRes(AR_heq,hfl,R) / calcTauMelt(L,Pi,C1,C3,R) #(–)
    T2 = calcTauRes(AR_heq,hfl,R) / calcTauCreep(C2,Pi) #(–)
    #print(T1,T2)
    
    ##Calculate h* and S*   
    h0 = h_eq_nd*initial_ratio
    S0 = S_eq_nd*initial_ratio
    #print("Z_vector=",z_vector, '  r_vector=',r_vector)
    sol = solve_ivp(dy_dt, 
                    (t0, tf), #initial time and end time. !! in non-dim!                                (param['t0'], param['tf']), #initial time and end time. !! in non-dim!
                    (h0,S0), #initial head and channel cross-section area. !!in non dim!
                    method = 'LSODA', #this method is the same 
                    args = (hfl,z_vector,r_vector,AR_heq,T1,T2), #r_fix,z_fix
    #                atol = 1e-6, #tolerance. can make it faster
    #                rtol = 1e-3,
                    max_step = 0.01 #!! non-dim time step!
                    )
    
    hnd = sol.y[0]  #(non-dim) Head
    Snd = sol.y[1] #(non-dim) Channel cross-section
    tnd = sol.t # (non-dim) Time vector 

    
    ##Calculate alpha* and beta*
    a = 5/4 *T1 *S_eq_nd**0.25 *h_eq_nd**1.5 - T2 *(1-h_eq_nd)**3
    b = 3/2 *T1 *S_eq_nd**(5/4) *h_eq_nd**0.5 + 3. *T2 *S_eq_nd *(1-h_eq_nd)**2
    c = -5/4 *h_eq_nd**0.5 *S_eq_nd**0.25
    d = -1/2 *h_eq_nd**0.5 *S_eq_nd**(5/4)
    
    p = a+d
    q = a*d - b*c
    
    alpha = p/2.
    beta = np.imag(np.sqrt(np.complex(p**2.-4*q)))/2.
     
    
    ''' Calculate alpha and beta with curvefit'''
    
    def MouTimFit(t, alpha, beta, Cst, phi):
        return np.exp(alpha*t) * Cst * np.sin(beta*t + phi) + h_eq_nd  
    
    p0 = [alpha,beta,0,np.pi] 
    try:
        pop, pcov = curve_fit(MouTimFit, tnd , hnd , p0=p0)#, bounds=(lowerbound,upperbound))   
    except RuntimeError:
        alpha_fit = 9999
        beta_fit = 9999
        Cst_fit = 9999
        phi_fit = 9999
    else:           
        alpha_fit = pop[0]
        beta_fit = pop[1]
        Cst_fit = pop[2]
        phi_fit = pop[3]
        
    #print("alpha=",alpha_fit, "beta=",beta_fit, 'C=',Cst_fit, 'phi=',phi_fit)
    hnd_fit = MouTimFit(tnd, alpha_fit, beta_fit, Cst_fit,phi_fit)#Cst_fit * np.exp(alpha_fit*tnd) * np.sin(beta_fit*tnd + phi) + h_eq 
    hnd_approx = MouTimFit(tnd, alpha, beta, Cst_fit,phi_fit) 
    hd_fit = hnd_fit *Pi/rhow/g
    hd_approx = hnd_approx *Pi/rhow/g                
    
    """ 
    Recovering dimensions and calculate timescales
    """
     
    #Calculate TauRes, TauMelt, TauCreep
    TauRes = calcTauRes(AR_heq,hfl,R)/sec_in_day #(day) from seconds to days
    TauCreep = calcTauCreep(C2,Pi)/sec_in_day #(1/days)
    TauMelt = calcTauMelt(L,Pi,C1,C3,R)/sec_in_day # #(1/days)
    #Calculate t, h and S with dimension
    hd = hnd * hfl #(m) 
    Sd = Snd * (R**2*L/C3**2/Pi)**(2/5) #(m) Sstar * S0
    td = tnd * TauRes # (days) to transform in days and set t0 to zero
    h_eq_d = h_eq_nd*hfl
    #Calculate damping and oscillation with dimension 
    damping = np.abs(1/alpha) * TauRes #(day)
    damping_fit = np.abs(1/alpha_fit)  * TauRes
    oscillation = 2*np.pi/beta * TauRes #(day)
    oscillation_fit = 2*np.pi/beta_fit * TauRes #(day)
    #Calculate Q
    Q = C3*Sd**(5/4)*np.sqrt(rhow*g*hd/L)


 
    #save in a dictionnary
    if shape == 'linear':
        resdic = {  'T1':T1,'T2':T2,\
                    'H':H, 'r_top':r_top, 'r_vector':r_vector, 'z_vector':z_vector,'r_base':r_base, 'r_heq':r_heq, 'm':m, 'R':R, 'L':L,\
                    'S0':S0, 'h0':h0,\
                    'tnd':tnd, 'hnd':hnd, 'Snd':Snd, 'hnd_fit':hnd_fit, 'hnd_approx':hnd_approx, 'hd_approx':hd_approx,\
                    'S_eq_nd':S_eq_nd, 'h_eq_nd':h_eq_nd, 'S_eq_approx':S_eq_approx, 'h_eq_approx': h_eq_approx,\
                    'S':Sd, 'h':hd, 't':td, 'h_eq_d':h_eq_d,\
                    'Q':Q,\
                    'alpha':alpha, 'beta':beta,\
                    'alpha_fit':alpha_fit, 'beta_fit':beta_fit, 'Cst_fit':Cst_fit, 'phi_fit':phi_fit, \
                    'damping':damping, 'damping_fit':damping_fit, 'oscillation_fit':oscillation_fit, 'oscillation':oscillation,\
                    'Sd':Sd, 'hd':hd, 'td':td, 'hd_fit':hd_fit,\
                    'TauCreep':TauCreep, 'TauRes':TauRes, 'TauMelt':TauMelt
                    } 

    if shape == 'nodes':
        resdic = {  'T1':T1,'T2':T2,\
            'H':H, 'r_vector':r_vector, 'z_vector':z_vector, 'r_heq':r_heq, 'R':R, 'L':L,\
            'S0':S0, 'h0':h0,\
            'r1':r1, 'r2':r2,\
            'tnd':tnd, 'hnd':hnd, 'Snd':Snd, 'hnd_fit':hnd_fit, 'hnd_approx':hnd_approx, 'hd_approx':hd_approx,\
            'S_eq_nd':S_eq_nd, 'h_eq_nd':h_eq_nd, 'S_eq_approx':S_eq_approx, 'h_eq_approx': h_eq_approx,\
            'S':Sd, 'h':hd, 't':td, 'h_eq_d':h_eq_d,\
            'Q':Q,\
            'alpha':alpha, 'beta':beta,\
            'alpha_fit':alpha_fit, 'beta_fit':beta_fit, 'Cst_fit':Cst_fit, 'phi_fit':phi_fit, \
            'damping':damping, 'damping_fit':damping_fit, 'oscillation_fit':oscillation_fit, 'oscillation':oscillation,\
            'Sd':Sd, 'hd':hd, 'td':td, 'hd_fit':hd_fit,\
            'TauCreep':TauCreep, 'TauRes':TauRes, 'TauMelt':TauMelt
            } 

    return resdic   


def plot_FR_JGR(results_dictionnary,variable='r',letter=['a','b','c']):
    
    res = results_dictionnary
    
    fig = plt.figure(figsize=(7,2.7),dpi=300)
    grid = plt.GridSpec(2, 30, wspace=0, hspace=0)
    ax1 = fig.add_subplot(grid[0:1, :19])
    ax2 = fig.add_subplot(grid[0:, 20:])
    ax3 = fig.add_subplot(grid[1, :19], sharex=ax1)

    colors = (Red, Orange,(0,0,0), Blue, Purple) #5colors

    
    for i in np.arange(len(res)):
        
                    
        ax1.plot(res[i]['td'],res[i]['hd'],color=colors[i],lw=2) 
        ax1.set_ylim([500,res[i]['H']])
        ax1.set_ylabel('h (m)')
        ax1.set_yticks([500,700,900]) 
        ax1.set_yticklabels([500,700,900])                
        ax1.get_yaxis().set_label_coords(-0.08,0.4)
        ax1.axes.xaxis.set_visible(False)
        

        ax2.plot(res[i]['r_vector'],res[i]['z_vector'],color=colors[i],lw=2)
        ax2.plot(-res[i]['r_vector'],res[i]['z_vector'],color=colors[i],lw=2)
        ax2.set_xlim([-30,30]) 
        ax2.set_ylim([0,1000])
        ax2.yaxis.tick_right()
        ax2.set_yticks([100,300,500,700,900]) 
        ax2.set_yticklabels([100,300,500,700,900])
        ax2.set_xlabel('r (m)')
        ax2.set_ylabel('z (m)')
        ax2.set_title('Moulin profile',fontsize=10)
        ax2.yaxis.set_label_position("right")
        
        ax3.plot(res[i]['td'],res[i]['Sd'],color=colors[i],lw=2)
        ax3.set_ylim([0.9,2])
        ax3.set_yticks([1,1.2,1.4,1.6]) 
        ax3.set_yticklabels([1,1.2,1.4,1.6])
        ax3.set_xticks([0,2,4,6,8,10,12,14,16,18,20]) 
        ax3.set_xticklabels([0,2,4,6,8,10,12,14,16,18,20])
        ax3.set_xlim([-1,20])
        ax3.set_ylabel('S (m)')
        ax3.set_xlabel('Days')
        ax3.get_yaxis().set_label_coords(-0.08,0.4)
        
        ax1.spines['top'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        
        ax1.spines['left'].set_bounds(500,900)
        ax2.spines['right'].set_bounds(100,900)
        ax2.spines['bottom'].set_bounds(-20,20)
        ax3.spines['left'].set_bounds(1,1.6)
        ax3.spines['bottom'].set_bounds(0,20)   
        
        ax1.tick_params(direction='in')
        ax2.tick_params(direction='in')
        ax3.tick_params(direction='in')
                        
    # zz = np.linspace(0,res[0]['h_eq_d'])
    # x1 =  res[0]['m']*zz + res[0]['r_base']
    # x2 =  res[-1]['m']*zz + res[-1]['r_base']
        
    # ax2.fill([-x1[0], -x1[-1], x1[-1], x1[0]],
    #     [zz[0],zz[-1],zz[-1],zz[0]],
    #     facecolor='#3990AC', alpha=0.2)
    # ax2.fill([-x2[0], -x2[-1], x2[-1], x2[0]],
    #     [zz[0],zz[-1],zz[-1],zz[0]],
    #     facecolor='#3990AC', alpha=0.2)
    
    
    ax1.text(0.07, 0.9, letter[0], ha='right', va='top', transform=ax1.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    ax2.text(0.1, 0.98, letter[1], ha='right', va='top', transform=ax2.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    ax3.text(0.07, 0.8, letter[2], ha='right', va='top', transform=ax3.transAxes,fontsize=8, 
              bbox=dict(facecolor='white', edgecolor='none', pad=1.0))
    

             
def timescales_subplot(reslist_cylinder,reslist_H2,reslist_heq):
    # matplotlib.rcParams.update({'font.size': 24})
    # figsize_1page=(7.48,3)
    
    # import matplotlib
    # matplotlib.rcParams['font.sans-serif'] = "Arial"
    # matplotlib.rcParams['font.family'] = "sans-serif"
    # matplotlib.rcParams['axes.labelsize'] = 8
    # matplotlib.rcParams['axes.labelweight'] = 'bold'
    # matplotlib.rc('xtick', labelsize=8)     
    # matplotlib.rc('ytick', labelsize=8)
    
    #colors:
    Red = '#ED2224'
    Orange = '#FBB263'
    #Green = '#A2D39E'
    Blue = '#15B4E9'
    Purple = '#6B52A2'
    
    colors = (Red, Orange ,(0,0,0), Blue, Purple) #5colors
            #[plt.cm.Greys(i) for i in np.linspace(0.4,1,5)]#[plt.cm.cividis(i) for i in np.linspace(1,0.3,5)]
    linestyles = ('solid','solid',(0, (1, 0.5)),'solid','solid')
    profile = reslist_cylinder['L']/1000
    osc_max = 5
    osc_min = 0
    dam_max = 10
    dam_min = -2
    
    fig, ax = plt.subplots(figsize=(7.48,3),dpi=300)
    #n_lines = 5
    gs1 = gridspec.GridSpec(2,3)
    gs1.update(wspace=0, hspace=0) # set the spacing between axes. 
    #frameon=False
    
    '''Change radius with fixed slope'''   
    ax1 = plt.subplot(gs1[0])
    for param, i in zip(reslist_cylinder['r'],range(5)):    
        timescale = np.array ([ result['oscillation'] 
                               for result in reslist_cylinder 
                               if result['r'] == param])  
        lines = plt.plot(profile,timescale,color=colors[i], lw=2, linestyle=linestyles[i], label='$=$%s'%param)
        plt.ylim([osc_min,osc_max])
        plt.yticks([1,2,3,4])
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.ylabel('$T_{\rm osc}$ (days)')
        ax1.text(2, osc_max, 'a',fontsize=12)
    
        
    ax4 = plt.subplot(gs1[3])#,sharex=ax1
    for param, i in zip(reslist_cylinder['r'],range(5)):    
        timescale = np.array ([ result['damping'] 
                               for result in reslist_cylinder 
                               if result['r'] == param])  
        lines = plt.plot(profile,timescale,color=colors[i], lw=2, linestyle=linestyles[i], label='$=$%s'%param)
        plt.ylim([dam_min,dam_max])
        plt.yticks([0,2,4,6,8])
        plt.ylabel('$T_{\rm damp}$ (days)')
        ax4.text(2, dam_max, 'd',fontsize=12)
        
        
        
    '''Change slope with h_middle fixed'''
    ax2 = plt.subplot(gs1[1])#, sharey=ax1
    for param, i in zip(reslist_H2['m'],range(5)):    
        timescale = np.array ([ result['oscillation'] 
                               for result in reslist_H2 
                               if result['m'] == param])  
        lines = plt.plot(profile,timescale,color=colors[i], lw=2, linestyle=linestyles[i], label='$=$%s'%param)
        plt.ylim([osc_min,osc_max])
        plt.setp(ax2.get_xticklabels(), visible=False)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.yticks([])
        ax2.text(2, osc_max, 'b',fontsize=12)
    
        
    ax5 = plt.subplot(gs1[4])#,sharex=ax2, sharey=ax4
    for param, i in zip(reslist_H2['m'],range(5)):    
        timescale = np.array ([ result['damping'] 
                               for result in reslist_H2 
                               if result['m'] == param])  
        lines = plt.plot(profile,timescale,color=colors[i], lw=2, linestyle=linestyles[i], label='$=$%s'%param)
        plt.ylim([dam_min,dam_max])
        plt.setp(ax5.get_yticklabels(), visible=False)
        plt.xlabel('Distance from margin (km)')
        plt.yticks([])
        ax5.text(2, dam_max, 'e',fontsize=12)
        
        
        
    '''Change slope with h_eq fixed'''
    ax3 = plt.subplot(gs1[2])#, sharey=ax1
    for param, i in zip(reslist_heq['m'],range(5)):    
        timescale = np.array ([ result['oscillation'] 
                               for result in reslist_heq 
                               if result['m'] == param])  
        lines = plt.plot(profile,timescale,color=colors[i], lw=2, linestyle=linestyles[i], label='$=$%s'%param)
        plt.ylim([osc_min,osc_max])
        plt.setp(ax3.get_xticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.yticks([])
        ax3.text(2, osc_max, 'c',fontsize=12)
    
        
    ax6 = plt.subplot(gs1[5])#,sharex=ax3, sharey=ax4
    for param, i in zip(reslist_heq['m'],range(5)):    
        timescale = np.array ([ result['damping'] 
                               for result in reslist_heq 
                               if result['m'] == param])  
        lines = plt.plot(profile,timescale,color=colors[i], lw=2, linestyle=linestyles[i], label='$=$%s'%param)
        plt.ylim([dam_min,dam_max])
        plt.setp(ax6.get_yticklabels(), visible=False)
        plt.yticks([])
        ax6.text(2, dam_max, 'f',fontsize=12)
    
    plt.subplots_adjust(left=0.08, right=0.99, top=0.96, bottom=0.14)


    plt.savefig('Timescales_osc_damp.pdf')
    
    
#resdic = calc_sim(r_fix = 5)
        
# %%


def plot_movie_JGR(results_dictionnary,idx):
    
    
    res = results_dictionnary
    Qout = np.array(res['Q'])
    H = res['H']
    td = np.array(res['td'])
    hd = np.array(res['hd'])
    Sd = np.sqrt(np.array(res['Sd']) *2 /np.pi)
    r_vec = np.array(res['r_vector'])
    z_vec = np.array(res['z_vector'])
    Qin = np.ones(len(td))*3
    
    
    #colors:
    # Red = '#ED2224'
    # Orange = '#FBB263'
    # Green = '#A2D39E'
    # Blue = '#15B4E9'
    # Purple = '#6B52A2'
    
    color_background = 'white'
    
    #colors = (1,1,1)# (Red, Orange,(0,0,0), Blue, Purple) #5colors
    
    fig = plt.figure(figsize=(7,2.7),dpi=300,facecolor=color_background)
    grid = plt.GridSpec(5, 30, wspace=0, hspace=0)
    

    ax2 = fig.add_subplot(grid[0:, :19]) 
    ax1 = fig.add_subplot(grid[3, :19], sharex=ax2)
    ax3 = fig.add_subplot(grid[4, :19], sharex=ax2)
    ax4 = fig.add_subplot(grid[0:, 20:])
    
    '''Qin'''
    
    ax1.plot(td,Qin,'--',color='grey',lw=1,label='Qin')
    ax1.plot(td[0:idx],Qin[0:idx],'-',color='black',lw=1.5,label='Qin')
    #ax1.plot(td,Qout,color='black',lw=1.5,label='Qout')
    ax1.set_ylabel('($m^3/s$)')
    #ax1.legend(loc="upper right")
    ax1.set_ylim([2,4])  
    ax1.set_yticks([2.5,3,3.5]) 
    ax1.set_yticklabels([2.5,3,3.5])                
    ax1.get_yaxis().set_label_coords(-0.08,0.4)
    ax1.axes.xaxis.set_visible(False)

    
    
    '''Head'''  
    ax2.plot(td,hd,'--',color='grey',lw=1)           
    ax2.plot(td[0:idx],hd[0:idx],color='black',lw=1.5) 
    #ax1.axhspan(500, 900, facecolor ='white', alpha = 1,zorder=1)
    ax2.set_ylim([0,H])
    ax2.set_ylabel('($m$)')
    ax2.set_yticks([500,700,900]) 
    ax2.set_yticklabels([500,700,900])            
    ax2.get_yaxis().set_label_coords(-0.08,0.7)
    ax2.axes.xaxis.set_visible(False)
    
    '''Moulin'''
    ax4.plot(r_vec,z_vec,color='black',lw=2)
    ax4.plot(-r_vec,z_vec,color='black',lw=2)
    ax4.set_xlim([-30,30]) 
    ax4.set_ylim([0,1000])
    ax4.yaxis.tick_right()
    ax4.set_yticks([100,300,500,700,900]) 
    ax4.set_yticklabels([100,300,500,700,900])
    ax4.set_xlabel('r ($m$)')
    ax4.set_ylabel('z ($m$)')
    ax4.set_title('Moulin profile',fontsize=10)
    ax4.yaxis.set_label_position("right")
    ax4.axhspan(0, hd[idx], facecolor ='lightblue', alpha = 1,zorder=1)
    ax4.axhspan(-20, 0, facecolor =color_background, alpha = 1,zorder=1)
    ax4.fill_betweenx(z_vec,-40,-r_vec,color='aliceblue',zorder=2)
    ax4.fill_betweenx(z_vec,r_vec,40,color='aliceblue',zorder=2)  
    
    '''Subglacial channel'''
    ax3.plot(td,Sd,'--',color='grey',lw=1)   
    ax3.plot(td[0:idx],Sd[0:idx],color='black',lw=1.5)
    #ax3.axhspan(1, 1.6, facecolor ='white', alpha = 1,zorder=1)
    ax3.set_ylim([0.8,1])
    ax3.set_yticks([0.85,0.9,0.95]) 
    ax3.set_yticklabels([0.85,0.9,0.95])
    ax3.set_xticks([0,2,4,6,8,10,12,14,16,18,20]) 
    ax3.set_xticklabels([0,2,4,6,8,10,12,14,16,18,20])
    ax3.set_xlim([-1,20])
    ax3.set_ylabel('($m$)')
    ax3.set_xlabel('Days')
    ax3.get_yaxis().set_label_coords(-0.08,0.4)
    
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.spines['left'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    ax1.spines['left'].set_bounds(2.5,3.5)
    ax2.spines['left'].set_bounds(500,900)
    ax4.spines['right'].set_bounds(100,900)
    ax4.spines['bottom'].set_bounds(-20,20)
    ax3.spines['left'].set_bounds(0.85,0.95)
    ax3.spines['bottom'].set_bounds(0,20)   
    
    ax1.tick_params(direction='in')
    ax2.tick_params(direction='in')
    ax4.tick_params(direction='in')
    ax3.tick_params(direction='in')
    
    ax1.patch.set_color(color_background)
    ax2.patch.set_color(color_background)
    ax4.patch.set_color(color_background)
    ax3.patch.set_color(color_background)
    
    color_text = 'black'
    #change color
    ax2.spines['left'].set_color(color_text)
    ax2.yaxis.label.set_color(color_text)
    ax2.tick_params(axis='y', colors=color_text)
    
    ax4.spines['left'].set_color(color_text)
    ax4.spines['bottom'].set_color(color_text)
    ax4.yaxis.label.set_color(color_text)
    ax4.tick_params(axis='y', colors=color_text)
    ax4.xaxis.label.set_color(color_text)
    ax4.tick_params(axis='x', colors=color_text)
    
    ax3.spines['right'].set_color(color_text)
    ax3.spines['bottom'].set_color(color_text)
    ax3.yaxis.label.set_color(color_text)
    ax3.tick_params(axis='y', colors=color_text)
    ax3.xaxis.label.set_color(color_text)
    ax3.tick_params(axis='x', colors=color_text)
                        
    
    ax1.text(0.1, 1, 'Discharge', ha='left', va='top', color='black',transform=ax1.transAxes,fontsize=8)
    ax2.text(0.1, 3.5, 'Head', ha='left', va='top', color='black',transform=ax1.transAxes,fontsize=8)
    ax3.text(0.1, 1, 'Subglacial channel radius', ha='left', va='top', transform=ax3.transAxes,fontsize=8)
    
#%%
#import FR_mainCode_plus as fr
# r=5
# resdic = calc_sim(r_fix=r)
# for idx in np.arange(0,len(resdic['td']),10):
#     plot_movie_JGR(resdic,idx,variable='r')
#     plt.pause(0.01)
